module Library

  # REQUIRED METHODS (CALLED BY DRIVER)

  def lib_build(env,prepkit)
    cmd="cd #{prepkit[:srcdir]} && #{env.build.cmd}"
    Thread.exclusive { ext(cmd,{:msg=>"Build failed"}) }
    prepkit
  end

  def lib_build_post(env,buildkit)
    buildkit
  end

  def lib_build_prep(env)
    srcdir=File.join(env.build.ddts_root,"src")
    cmd="rsync -a #{File.join($DDTSHOME,"..","src")}/ #{srcdir}"
    Thread.exclusive { ext(cmd,{:msg=>"Error copying 'src' to '#{srcdir}'"}) }
    rundir=File.join(env.build.ddts_root,"run")
    cmd="rsync -a --no-recursive #{File.join($DDTSHOME,"..","run")}/* #{rundir}"
    Thread.exclusive {ext(cmd,{:msg=>"Error copying 'run' to '#{rundir}'"}) }
    bindir=File.join(env.build.ddts_root,"bin")
    FileUtils.mkdir_p(bindir)
    {:bindir=>bindir,:rundir=>rundir,:srcdir=>srcdir}
  end

  def lib_data_zeus(env)
    f="/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/IPEdata/ipedata.tgz"
    cmd="cp --force #{f} data.tgz"
    md5='f944709c93f2daf6a62ce10ed8d93006'
    [cmd,md5]
  end

  # def lib_outfiles(path)
  def lib_outfiles(env,path)
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

  def lib_queue_del_cmd(env)
    'qdel'
  end

  #### HERE ####

  # def lib_run_job(rundir,runspec,lock,activejobs)
  def lib_run(env,prepkit)
    jobid=nil
    subdir=nil
    re1=Regexp.new(re_str_job_id)
    re2=Regexp.new(re_str_run_dir)
    datadir=valid_dir(File.join(FileUtils.pwd,"data"))
    ipedata="IPEDATA=#{datadir}"
    ipequeue="IPEQUEUE=batch"
    cmd="cd #{rundir} && #{ipedata} #{ipequeue} ./#{runspec['qsubcmd']} #{runspec['tasks']}"
    logd "Submitting job with command: #{cmd}"
    Thread.exclusive { output,status=ext(cmd,{:msg=>"ERROR: Job submission failed"}) }
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

  def lib_run_check(env,postkit)
  end

  def lib_run_prep(env)
    uniq=env.run.ddts_root
    run=env.build.ddts_result[:rundir]
    logd "Copying #{run} -> #{uniq}"
    FileUtils.cp_r(run,uniq)
    bin=env.build.ddts_result[:bindir]
    logd "Linking #{bin} -> #{uniq}"
    FileUtils.ln_s(bin,uniq)
    rundir=File.join(uniq,File.basename(run))
    inpsrc=valid_file(File.join(rundir,env.run.inpfile))
    inpdst=File.join(rundir,"IPE.inp")
    logd "Copying #{inpsrc} -> #{inpdst}"
    FileUtils.rm_f(inpdst)
    FileUtils.cp(inpsrc,inpdst)
    nlfile=valid_file(File.join(rundir,'SMSnamelist'))
    mod_namelist_file(nlfile,env.run.namelists)
    rundir
  end

  def lib_run_post(env,runkit)
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

  # CUSTOM METHODS (NOT CALLED BY DRIVER)

# def lib_re_str_success
#   'IPE completed successfully'
# end

  def mod_namelist_file(nlfile,nlenv)
    h=convert_o2h(nlenv)
    sets=h.reduce([]) do |m0,(n,kv)|
      inner=kv.reduce([]) do |m1,(k,v)|
        v="\"#{quote_string(v)}\""
        logd "Set namelist #{n}:#{k}=#{v}"
        m1.push("-s #{n}:#{k}=#{v}")
      end
      m0.concat(inner)
    end
    nml=valid_file(File.expand_path(File.join($DDTSHOME,"nml")))
    cmd="#{nml} -i #{nlfile} -o #{nlfile} #{sets.join(" ")}"
    Thread.exclusive { ext(cmd,{:msg=>"Failed to edit #{nlfile}"}) }
  end

  def re_str_job_id
    'The job (\d+).* has been submitted.'
  end

  def re_str_run_dir
    'Created (.*)'
  end

end
