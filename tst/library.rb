module Library

  def lib_build(buildspec,lock,activejobs)
    cmd="cd #{buildspec['buildsrc']} && #{buildspec['cmd']}"
    ext(cmd,{:msg=>"Build failed"})
  end

  def lib_build_post(buildspec,output=nil)
    buildspec['buildrun']
  end

  def lib_build_prep(buildspec)
    def copy_files(desc,name,key,rsync_args,topdir,dstdir,buildspec)
      src=valid_dir(File.join(topdir,name))
      dst=File.join(dstdir,name)
      FileUtils.mkdir_p(dst) unless Dir.exist?(dst)
      cmd="rsync #{rsync_args} #{src}/* #{dst}"
      ext(cmd)
      logd "Copied #{src} into #{dst}"
      d=valid_dir(dst)
      buildspec[key]=d
      logd "Set build #{desc} directory: #{d}"
    end
    dstdir=File.join(buildspec['buildroot'],buildspec['build'])
    FileUtils.mkdir_p(dstdir) unless Dir.exist?(dstdir)
    logd "Made directory: #{dstdir}"
    topdir=File.expand_path("..")
    copy_files('source','src','buildsrc','-av',topdir,dstdir,buildspec)
    copy_files('run','run','buildrun','-av --no-recursive',topdir,dstdir,buildspec)
    ['bin'].each do |e|
      d=File.join(dstdir,e)
      FileUtils.mkdir_p(d) unless Dir.exist?(d)
      logd "Made directory: #{d}"
    end
  end

  def lib_dataspecs
    f="ipedata.tgz"
#nm20140206    path="/scratch2/portfolios/BMC/acb/IPEdata/ipedata.tgz"
    path="/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/IPEdata/ipedata.tgz"
    cmd="cp --force #{path} data.tgz"
#nm20140206    md5='8598ce1c1eee6f9fce573335ab3f3714'
#nm20140415    md5='d2c71804a4c097db143e6b4745cac9d2'
    md5='f944709c93f2daf6a62ce10ed8d93006'
    [cmd,md5]
  end

  def lib_outfiles(path)
    restrs=['(.*/)(plasma\d\d)','(.*/)(fort\.20\d\d)']
    res=restrs.map { |e| Regexp.new(e) }
    outs=[]
    Find.find(path) do |e|
      res.each do |re|
        m=re.match(e)
        unless m.nil?
          outs << [m[1],m[2]]
          break
        end
      end
    end
    outs
  end

  def lib_prep_job(rundir,runspec)
    buildrun=runspec['buildrun']
    logd "Copying #{buildrun} -> #{rundir}"
    FileUtils.cp_r(buildrun,rundir)
    bindir=File.expand_path(File.join(buildrun,"..","bin"))
    logd "Linking #{bindir} -> #{rundir}"
    FileUtils.ln_s(bindir,rundir)
    rundir=File.join(rundir,File.basename(buildrun))
    inpsrc=valid_file(File.join(rundir,runspec['inpfile']))
    inpdst=File.join(rundir,"IPE.inp")
    logd "Copying #{inpsrc} -> #{inpdst}"
    FileUtils.rm_f(inpdst)
    FileUtils.cp(inpsrc,inpdst)
    nlfile=valid_file(File.join(rundir,'SMSnamelist'))
    mod_namelist_file(nlfile,runspec['namelists'])
    rundir
  end

  def lib_queue_del_cmd
    'qdel'
  end

  def lib_re_str_job_id
    'The job (\d+).* has been submitted.'
  end

  def lib_re_str_run_dir
    'Created (.*)'
  end

  def lib_re_str_success
    'IPE completed successfully'
  end

  def lib_run_job(rundir,runspec,lock,activejobs)
    jobid=nil
    subdir=nil
    re1=Regexp.new(lib_re_str_job_id)
    re2=Regexp.new(lib_re_str_run_dir)
    datadir=valid_dir(File.join(FileUtils.pwd,"data"))
    ipedata="IPEDATA=#{datadir}"
    ipequeue="IPEQUEUE=debug"
    cmd="cd #{rundir} && #{ipedata} #{ipequeue} ./#{runspec['qsubcmd']} #{runspec['tasks']}"
    logd "Submitting job with command: #{cmd}"
    output,status=ext(cmd,{:msg=>"ERROR: Job submission failed"})
    output.each do |e|
      e.chomp!
      logd e unless e=~/^\s*$/
      jobid=e.gsub(re1,'\1') if re1.match(e)
      subdir=e.gsub(re2,'\1') if re2.match(e)
    end
    if jobid.nil?
      logi "ERROR: Job ID not found in queue-submission output"
      return nil
    end
    lock.synchronize { activejobs[jobid]=self }
    if subdir.nil?
      logi "ERROR: Run directory not found in queue-submission output"
      return nil
    end
    runspec['subdir']=subdir
    rundir.replace(File.join(rundir,subdir))
    qs="Queued with job ID #{jobid}"
    logi qs
    lib_wait_for_job(jobid)
    lock.synchronize { activejobs.delete(jobid) }
    valid_file(File.join(rundir,'output'))
  end

  def lib_wait_for_job(jobid)
    ok=%w[E H Q R T W S]
    # 'tolerance' is the number of seconds the batch system retains information
    # about completed jobs. If this interval passes without a non-error response
    # to queries, we may never receive confirmation that the job completed.
    # Consider this a batch-system failure and abort.
    tolerance=600
    batch_failure=false
    last_response=Time.now
    begin
      sleep 10
      tolog=[]
      cmd="qstat -f #{jobid}"
      output,status=ext(cmd,{:die=>false,:out=>false})
      if status==0
        live=false
        last_response=Time.now
        output.each do |e|
          tolog.push(e)
          live=true if ok.include?(e.chomp.sub(/^ *job_state = (.)$/,'\1'))
        end
      else
        live=true
        now=Time.now
        logd "#{cmd} set error status #{status} at #{now}"
        if now-last_response > tolerance
          batch_failure=true
          live=false
        end
      end
    end while live
    logd "* Final batch info for job ID #{jobid}:"
    logd "--"
    tolog.each { |e| logd e }
    logd "--"
    die "* Batch system unresponsive for #{tolerance} seconds" if batch_failure
    re=Regexp.new('^ *exit_status = (\d+)')
    tolog.each do |e|
      m=re.match(e)
      return m[1].to_i if m
    end
    'unknown'
  end

end
