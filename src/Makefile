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
PACKAGE=../bin/ipe_${TEST}.exe
# DR0: modules
DR0=${HOME_dev}/main
FFLAG0=${FFLAGS90} ${INCLUDE_PATH}${DR3}
# DR1: neutral
DR1=${HOME_dev}/neutral
FFLAG1=${FFLAGS90}  ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR3} ${SMSFLAGS}
# DR2: FLIP
DR2=${HOME_dev}/flip
FFLAG2=${FFLAGS}      ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR1} ${INCLUDE_PATH}${DR3} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS} ${SMSFLAGS}
# DR3: plasma
DR3=${HOME_dev}/plasma
FFLAG3=${FFLAGS90}    ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR6} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS} ${SMSFLAGS}
# DR6: eldyn
DR6=${HOME_dev}/eldyn
FFLAG6=${FFLAGS_e}    ${INCLUDE_PATH}${DR0} ${SMSFLAGS}
# DR4: main driver
DR4=${HOME_dev}/main
FFLAG4=${FFLAGS90} ${INCLUDE_PATH}${DR1} ${INCLUDE_PATH}${DR3} ${INCLUDE_PATH}${DR6} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
SRC0	=\
	${DR0}/module_decomp.f90 \
	${DR0}/module_precision.f90 \
	${DR0}/module_IPE_dimension.f90 \
	${DR0}/module_physical_constants.f90 \
	${DR0}/module_unit_conversion.f90 \
	${DR0}/module_input_parameters.f90 \
	${DR0}/module_io.f90 \
	${DR0}/module_sub_io.f90 \
	${DR0}/module_field_line_grid.f90 \
	${DR0}/module_sub_field_line_grid.f90
# DR1: neutrals
SRC1	=\
	${DR1}/nrlmsise00.f90 \
	${DR1}/hwm93.f90 \
	${DR1}/thermosphere_1d_v2.f90 \
	${DR1}/module_neutral.f90 \
	${DR1}/module_heating_rate.f90
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
	${DR3}/flux_tube_solver.f90 \
	${DR3}/perpendicular_transport.f90 \
	${DR3}/stepback_mag.f90 \
	${DR3}/interpolate_flux_tube_${TEST}.f90 \
	${DR3}/io_plasma_bin.f90 \
	${DR3}/read_vexb.f90


#dbg20120509	${DR3}/activate_perp_transport.f90 \
# DR6: eldyn
SRC6	=\
	${DR6}/efield.f \
	${DR6}/module_eldyn.f \
	${DR6}/module_sub_eldyn.f \
	${DR6}/module_magfield.f \
	${DR6}/module_sunloc.f \
	${DR6}/get_efield90km.f \
	${DR6}/apex_sunloc.f \
	${DR6}/supim_exb.f
# DR4: main driver
SRC4	=\
	${DR4}/get_flip_grid.f90 \
	${DR4}/allocate_arrays.f90 \
	${DR4}/output_plasma_grid.f90 \
	${DR4}/stop.f90 \
	${DR4}/driver_ipe.f90
#
#HEADS	= $(PACKAGE).h
OBJ0	= $(SRC0:.f90=.o)
OBJ1	= $(SRC1:.f90=.o)
OBJ2	= $(SRC2:.f=.o)
OBJ3	= $(SRC3:.f90=.o)
OBJ4	= $(SRC4:.f90=.o)
OBJ6	= $(SRC6:.f=.o)
#
FILES	= README Makefile $(HEADS) $(SRC0) $(SRC1) $(SRC2) $(SRC3) $(SRC6) $(SRC4)
VER	= `date +%Y%m%d`


### to avoid m2c... !!!didn't work!!!
###%.o: %.mod

all: $(PACKAGE)

FINAL_OBJS=\
	$(OBJ0) \
	$(OBJ1) \
	$(OBJ6) \
	$(OBJ3) \
	$(OBJ2) \
	$(OBJ4)

$(PACKAGE): D0_OBJ  D1_OBJ  D6_OBJ  D3_OBJ  D2_OBJ  D4_OBJ

ifeq ($(PAR),sms)
	$(FCparallel) $(LDFLAGS) $(FINAL_OBJS) -o $@ $(LDLIBS) $(GPTL_LDFLAGS) $(SMS_LDFLAGS)
else
	$(FC) $(LDFLAGS) $(FINAL_OBJS) -o $@ $(LDLIBS) $(GPTL_LDFLAGS) $(SMS_LDFLAGS)
endif

$(FINAL_OBJS): $(HEADS) Makefile

 D0_OBJ :
	cd ${DR0} ; ${MAKE} all \
	"OBJS=${OBJ0}" \
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
	"OBJ1=${OBJ1}" \
	"FC=${FC}" \
	"FFLAG1=${FFLAG1}" \
	"LFLAG1=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}"
 D2_OBJ :
	cd ${DR2} ; ${MAKE} all \
	"OBJ2=${OBJ2}" \
	"FC=${FC}" \
	"FFLAG2=${FFLAG2}" \
	"LFLAG2=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}"
 D3_OBJ :
	cd ${DR3} ; ${MAKE} all \
	"OBJ3=${OBJ3}" \
	"FC=${FC}" \
	"FFLAG3=${FFLAG3}" \
	"LFLAG3=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}"
 D6_OBJ :
	cd ${DR6} ; ${MAKE} all \
	"OBJ6=${OBJ6}" \
	"FC=${FC}" \
	"FFLAG6=${FFLAG6}" \
	"LFLAG6=${LDFLAGS}" \
	"HOME_dev=${HOME_dev}"
 D4_OBJ :
	cd ${DR4} ; ${MAKE} all \
	"OBJS=${OBJ4}" \
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
	$(RM) ${DR0}/*.mod ${DR0}/*.tmp ${DR0}/*.f ${DR0}/*~
	$(RM) ${DR0}/*_sms.f90 ${DR0}/*ppptmp*
	$(RM) ${DR1}/*.mod ${DR2}/*.mod ${DR3}/*.mod ${DR4}/*.mod ${DR6}/*.mod 

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
