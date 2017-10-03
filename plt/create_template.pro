;date: 20170215
;purpose: create/save a template to read ascii akebono data
pro create_template

  WORKDIR='/home/naomi.maruyama/scrub/tmp/tmp20170215akebono/'
  myfile=WORKDIR+'ne-20130316.txt'
  mytemplate = ascii_template(myfile)
  save,mytemplate,FILENAME=WORKDIR+'mytemplateAkebono.sav'

  print,'pro create_template finished'
end;pro create_template
