###remember:module_plasma.3d.dbg.f90 is used for trasn dbg purpose!!!
#! DATE: 08 September, 2011
#!********************************************
#!***      Copyright 2011 NAOMI MARUYAMA   ***
#!***      ALL RIGHTS RESERVED             ***
#!********************************************
#! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
#! DEVELOPER: Dr. Naomi Maruyama
#! CONTACT INFORMATION:
#! E-MAIL : Naomi.Maruyama@noaa.gov
#! PHONE  : 303-497-4857
#! ADDRESS: 325 Broadway, Boulder, CO 80305
#!--------------------------------------------  
# Makefile
## you have to setenv/export
#(1) COMPILER=pgf95, lf648, ifort, or xlf90
#(2) machine=$(machine)_$(COMPILER)
include  Make.$(machine)
TEST=v18
PACKAGE=../bin/ipe_${TEST}.$(COMPILER).exe
# DR0: modules
DR0=${HOME_dev}/main
FFLAG0=${FFLAGS90} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
# DR1: neutral
DR1=${HOME_dev}/neutral
FFLAG1=${FFLAGS90}  ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR3}
# DR2: FLIP
DR2=${HOME_dev}/flip
FFLAG2=${FFLAGS}      ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR1} ${INCLUDE_PATH}${DR3} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
# DR3: plasma
DR3=${HOME_dev}/plasma
FFLAG3=${FFLAGS90}    ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR6} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
# DR6: eldyn
DR6=${HOME_dev}/eldyn
FFLAG6=${FFLAGS_e}    ${INCLUDE_PATH}${DR0}
# DR4: main driver
DR4=${HOME_dev}/main
FFLAG4=${FFLAGS90} ${INCLUDE_PATH}${DR1} ${INCLUDE_PATH}${DR3} ${INCLUDE_PATH}${DR6} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}

SRC0S	=\
	${DR0}/module_precision.f90 \
	${DR0}/module_IPE_dimension.f90 \
	${DR0}/module_physical_constants.f90 \
	${DR0}/module_unit_conversion.f90 \
	${DR0}/module_io.f90 \
	${DR0}/module_open_file.f90 \
	${DR0}/module_output.f90 \
	${DR0}/module_close_files.f90
SRC0P	=\
	${DR0}/module_decomp.f90 \
	${DR0}/module_input_parameters.f90 \
	${DR0}/module_open_output_files.f90 \
	${DR0}/module_field_line_grid.f90 \
	${DR0}/allocate_arrays.f90 \
	${DR0}/get_flip_grid.f90 \
	${DR0}/module_read_plasma_grid_global.f90 \
	${DR0}/module_init_plasma_grid.f90 \
	${DR0}/output_plasma_grid.f90 \
	${DR0}/stop.f90
# DR1: neutrals
SRC1S	=\
	${DR1}/nrlmsise00.f90 \
	${DR1}/hwm93.f90 \
	${DR1}/thermosphere_1d_v2.f90 
SRC1P	=\
	${DR1}/module_neutral.f90
# DR2: FLIP
SRC2	=\
	${DR2}/CTIPE-int.f \
	${DR2}/initialize_module_parameters.f \
	${DR2}/CMINOR.f \
	${DR2}/ELECXS.f \
	${DR2}/INIT-PROFILES.f \
	${DR2}/KEMPRN.f \
	${DR2}/MINORA.f \
	${DR2}/Neut_Heating.f \
	${DR2}/Photoel-Freqs.f \
	${DR2}/RSDENA_EXB.f \
	${DR2}/RSJACA.f \
	${DR2}/RSLPSD.f \
	${DR2}/RSLPST.f \
	${DR2}/RSPE2B.f \
	${DR2}/RSPRIM.f \
	${DR2}/RSTEMC_EXB.f \
	${DR2}/Rates.f \
	${DR2}/FLIP_GRID.f
# DR3: plasma
SRC3	=\
	${DR3}/module_plasma.f90 \
	${DR3}/module_sub_plasma.f90 \
	${DR3}/get_pvalue_dipole.f90 \
	${DR3}/get_sza.f90 \
	${DR3}/get_sinim.f90 \
	${DR3}/module_heating_rate.f90 \
	${DR3}/flux_tube_solver.f90 \
	${DR3}/perpendicular_transport.f90 \
	${DR3}/stepback_mag.f90 \
	${DR3}/interpolate_flux_tube_${TEST}.f90 \
	${DR3}/io_plasma_bin.f90 \
	${DR3}/read_vexb.f90


#dbg20120509	${DR3}/activate_perp_transport.f90 \
# DR6: eldyn
SRC6S	=\
	${DR6}/efield.f \
	${DR6}/module_ff.f \
	${DR6}/module_pnm.f \
	${DR6}/module_prep_pnm.f \
	${DR6}/module_index_quiet.f \
	${DR6}/module_adj_S_a.f \
	${DR6}/module_constants.f \
	${DR6}/module_prep_fk.f \
	${DR6}/module_set_fkflfs.f \
	${DR6}/module_efield_mid.f \
	${DR6}/module_FSVal.f \
	${DR6}/module_SetModel.f \
	${DR6}/module_ADJUST.f \
	${DR6}/module_JULDAY.f \
	${DR6}/module_GET_TILT.f \
	${DR6}/module_prep_weimer.f \
	${DR6}/module_pot_latsmo.f \
	${DR6}/module_pot_latsmo2.f \
	${DR6}/module_pot_lonsmo.f \
	${DR6}/module_efield_sunloc.f \
	${DR6}/module_fun_MLT.f \
	${DR6}/module_fun_MagLong.f \
	${DR6}/module_highlat_adjust.f \
	${DR6}/module_interp_poten.f \
	${DR6}/module_DerivPotential.f \
	${DR6}/module_LEGENDRE.f \
	${DR6}/module_EpotVal.f \
	${DR6}/module_GECMP.f \
	${DR6}/module_ReadCoef.f \
	${DR6}/module_highlat_getbnd.f \
	${DR6}/module_svbksb.f \
	${DR6}/module_svdcmp.f \
	${DR6}/module_bnd_sinus.f \
	${DR6}/module_GlobalElPotential.f \
	${DR6}/module_get_efield.f \
	${DR6}/module_magfield.f \
	${DR6}/module_sunloc.f \
	${DR6}/apex_sunloc.f \
	${DR6}/supim_exb.f
SRC6P	=\
	${DR6}/module_read_acoef.f \
	${DR6}/module_efield_init.f \
	${DR6}/module_eldyn.f      \
	${DR6}/module_init_eldyn.f \
	${DR6}/module_sub_eldyn.f  \
	${DR6}/get_efield90km.f 
# DR4: main driver
SRC4P	=\
	${DR4}/driver_ipe.f90

#HEADS	= $(PACKAGE).h
OBJ0S	= $(SRC0S:.f90=.o)
OBJ0P	= $(SRC0P:.f90=.o)
OBJ1S	= $(SRC1S:.f90=.o)
OBJ1P	= $(SRC1P:.f90=.o)
OBJ2	= $(SRC2:.f=.o)
OBJ3	= $(SRC3:.f90=.o)
OBJ4S	= 
OBJ4P	= $(SRC4P:.f90=.o)
OBJ6S	= $(SRC6S:.f=.o)
OBJ6P	= $(SRC6P:.f=.o)
#
FILES	= README Makefile $(HEADS) $(SRC0S) $(SRC0P) $(SRC1S) $(SRC1P) $(SRC2) $(SRC3) $(SRC6S) $(SRC6p) $(SRC4P)
VER	= `date +%Y%m%d`


### to avoid m2c... !!!didn't work!!!
###%.o: %.mod

all: $(PACKAGE)

FINAL_OBJS=\
	$(OBJ0S) \
	$(OBJ0P) \
	$(OBJ1S) \
	$(OBJ1P) \
	$(OBJ6S) \
	$(OBJ6P) \
	$(OBJ3) \
	$(OBJ2) \
	$(OBJ4P)

$(PACKAGE): D0_OBJ  D1_OBJ  D6_OBJ  D3_OBJ  D2_OBJ  D4_OBJ

ifeq ($(PAR),sms)
	which mpif90
	$(FCparallel) $(LDFLAGS) $(FINAL_OBJS) -o $@ $(LDLIBS) $(GPTL_LDFLAGS) $(SMS_LDFLAGS)
else
	$(FC)         $(LDFLAGS) $(FINAL_OBJS) -o $@ $(LDLIBS) $(GPTL_LDFLAGS) $(SMS_LDFLAGS)
endif

$(FINAL_OBJS): $(HEADS) Makefile

 D0_OBJ :
	cd ${DR0} ; ${MAKE} all \
	"OBJS=${OBJ0S}" \
	"OBJP=${OBJ0P}" \
	"FC=${FC}" \
	"FCparallel=${FCparallel}" \
	"FFLAGS=${FFLAG0}" \
	"LFLAGS=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}" \
	"PAR=${PAR}" \
	"CPP=${CPP}" \
	"CPP_FLAGS=${CPP_FLAGS}" \
	"PPP=${PPP}" \
	"PPP_FLAGS=${PPP_FLAGS}" \
	"PPP_FLAGS1=${PPP_FLAGS1}" \
	"SMSFLAGS=${SMSFLAGS}" \
	"PINCLUDES=${PINCLUDES}" 
 D1_OBJ :
	cd ${DR1} ; ${MAKE} all \
	"OBJS=${OBJ1S}" \
	"OBJP=${OBJ1P}" \
	"FC=${FC}" \
	"FCparallel=${FCparallel}" \
	"FFLAGS=${FFLAG1}" \
	"LFLAGS=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}" \
	"PAR=${PAR}" \
	"CPP=${CPP}" \
	"CPP_FLAGS=${CPP_FLAGS}" \
	"PPP=${PPP}" \
	"PPP_FLAGS=${PPP_FLAGS}" \
	"PPP_FLAGS1=${PPP_FLAGS1}" \
	"SMSFLAGS=${SMSFLAGS}" \
	"PINCLUDES=${PINCLUDES}" 

 D2_OBJ :
	cd ${DR2} ; ${MAKE} all \
	"OBJ2=${OBJ2}" \
	"FC=${FC}" \
	"FFLAG2=${FFLAG2}" \
	"LFLAG2=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}"
 D3_OBJ :
	cd ${DR3} ; ${MAKE} all \
	"OBJS=${OBJ3}" \
	"FC=${FC}" \
	"FCparallel=${FCparallel}" \
	"FFLAGS=${FFLAG3}" \
	"LFLAGS=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}" \
	"PAR=${PAR}" \
	"CPP=${CPP}" \
	"CPP_FLAGS=${CPP_FLAGS}" \
	"PPP=${PPP}" \
	"PPP_FLAGS=${PPP_FLAGS}" \
	"PPP_FLAGS1=${PPP_FLAGS1}" \
	"SMSFLAGS=${SMSFLAGS}" \
	"PINCLUDES=${PINCLUDES}" 
 D6_OBJ :
	cd ${DR6} ; ${MAKE} all \
	"OBJS=${OBJ6S}" \
	"OBJP=${OBJ6P}" \
	"FC=${FC}" \
	"FCparallel=${FCparallel}" \
	"FFLAGS=${FFLAG6}" \
	"LFLAGS=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}" \
	"PAR=${PAR}" \
	"CPP=${CPP}" \
	"CPP_FLAGS=${CPP_FLAGS}" \
	"PPP=${PPP}" \
	"PPP_FLAGS=${PPP_FLAGS}" \
	"PPP_FLAGS1=${PPP_FLAGS1}" \
	"SMSFLAGS=${SMSFLAGS}" \
	"PINCLUDES=${PINCLUDES}" 
 D4_OBJ :
	cd ${DR4} ; ${MAKE} all \
	"OBJP=${OBJ4S}" \
	"OBJP=${OBJ4P}" \
	"FC=${FC}" \
	"FCparallel=${FC}" \
	"FFLAGS=${FFLAG4}" \
	"LFLAGS=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}" \
	"PAR=${PAR}" \
	"CPP=${CPP}" \
	"CPP_FLAGS=${CPP_FLAGS}" \
	"PPP=${PPP}" \
	"PPP_FLAGS=${PPP_FLAGS}" \
	"PPP_FLAGS1=${PPP_FLAGS1}" \
	"SMSFLAGS=${SMSFLAGS}" \
	"PINCLUDES=${PINCLUDES}" 

### useful commands ###

cln:
	$(RM) $(PACKAGE) $(FINAL_OBJS)
	$(RM) ${DR0}/*.mod ${DR0}/*.tmp* ${DR0}/*.f ${DR0}/*~
	$(RM) ${DR0}/*_sms.f90 ${DR0}/*ppptmp*
	$(RM) ${DR1}/*.mod ${DR2}/*.mod ${DR3}/*.mod ${DR4}/*.mod ${DR6}/*.mod 
	$(RM) ${DR1}/*.mod ${DR1}/*.tmp* ${DR1}/*.f ${DR1}/*~
	$(RM) ${DR1}/*_sms.f90 ${DR1}/*ppptmp*
	$(RM) ${DR3}/*.mod ${DR3}/*.tmp* ${DR3}/*.f ${DR3}/*~
	$(RM) ${DR3}/*_sms.f90 ${DR3}/*ppptmp*
	$(RM) ${DR6}/*.mod ${DR6}/*.tmp* ${DR6}/*~
	$(RM) ${DR6}/*_sms.f90 ${DR6}/*ppptmp*

tar:
	@echo $(PACKAGE)-$(VER) > .package
	@$(RM) -r `cat .package`
	@mkdir `cat .package`
	@ln $(FILES) `cat .package`
	tar cvf - `cat .package` | gzip -9 > `cat .package`.tar.gz
	@$(RM) -r `cat .package` .package

zip:
	zip -9 $(PACKAGE)-$(VER).zip $(FILES)


prof: run
	$(PROF) $(PACKAGE) | less

run: all
	./$(PACKAGE) < sample-data | less
