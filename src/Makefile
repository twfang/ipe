# IPE Makefile

# Static/Unconditional/Default Variables

CPP=/lib/cpp
DEBUG=no
GPTL_AUTOINST=no
GPTL_HAVE_PAPI=no
GPTL_CPPFLAGS=-I../dummygptl
GPTL_LDFLAGS=-L../../dummygptl -lgptl
GPTL_USE_GPTL=no
GPTL_USE_PMPI=no
LDLIBS=-lm
PPP_FLAGS=--Fmodule=module_decomp --Free --comment --HaloSize=10
PROF=gprof
TARGETS=zeus_intel_parallel zeus_intel_serial zeus_lahey_serial clean

# Dynamic/Conditional Variables

CPP_FLAGS=-traditional -P $(GPTL_CPPFLAGS)
FFLAGS=$(OPT_FLAG) $(DEBUG_FLAGS) $(GPTL_CPPFLAGS) $(GPTL_LDFLAGS)
LDFLAGS=$(OPT_FLAG) $(DEBUG_FLAGS)
PPP=$(SMS)/bin/ppp
SMS_LDFLAGS=-L$(SMS)/lib -lsms
SMSFLAGS=-I $(SMS)/include

TOP=$(CURDIR)

ifeq ($(MACHINE),zeus)
  ifeq ($(COMPILER),intel)
    ENVSET=module load intel mpt
    FCP=mpif90
    FCS=ifort
    GPTL_LDFLAGS=-L/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/lib -lgptl_pmpi
    OPT_FLAG=-O3 -fp-model precise
    SMS=/home/Jacques.Middlecoff/SMS/SMS_r395_intel_mpt/
    ifeq ($(DEBUG),yes)
      DEBUG_FLAGS=-g -O0 -ftz -traceback -check bounds -fpe0
    endif
    ifeq ($(PARALLELISM),parallel)
      ifeq ($(GPTL_USE_GPTL),yes)
        GPTL_CPPFLAGS=-I/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_mpi/include
        ifeq ($(GPTL_USE_PMPI),no)
          GPTL_LDFLAGS=-L/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_mpi/lib -lgptl
        endif
      endif
    endif
    ifeq ($(PARALLELISM),serial)
      ifeq ($(GPTL_USE_GPTL),yes)
        GPTL_CPPFLAGS=-I/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_nompi/include
        ifeq ($(GPTL_USE_PMPI),no)
          GPTL_LDFLAGS=-L/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_nompi/lib -lgptl
        endif
      endif
    endif
  endif
  ifeq ($(COMPILER),lahey)
    ENVSET=module load lahey mpt
    FCP=mpif90
    FCS=lf95
    SMS=/home/Jacques.Middlecoff/SMS/SMS_r395_lahey_mpt/
    ifeq ($(PARALLELISM),parallel)
    endif
    ifeq ($(PARALLELISM),serial)
    endif
  endif
endif

ifeq ($(GPTL_AUTOINST),yes)
  GPTL_LDFLAGS+=-finstrument-functions
endif

ifeq ($(GPTL_HAVE_PAPI),yes)
  GPTL_LDFLAGS+=-L/path/to/papi/lib -lpapi -lpfm
endif

VARS="ENVSET="$(ENVSET)" FCP="$(FCP)" FCS="$(FCS)" FFLAGS="$(FFLAGS)" PARALLELISM="$(PARALLELISM)" PPP="$(PPP)" PPP_FLAGS="$(PPP_FLAGS)" SMSFLAGS="$(SMSFLAGS)""

# Functions

USAGE=$(error Valid targets: $(TARGETS))

# Targets

all:
	$(call $(USAGE))

zeus_intel_parallel:
	$(MAKE) ipe MACHINE=zeus COMPILER=intel PARALLELISM=parallel

zeus_intel_serial:
	$(MAKE) ipe MACHINE=zeus COMPILER=intel PARALLELISM=serial

zeus_lahey_serial:
	$(MAKE) ipe MACHINE=zeus COMPILER=lahey PARALLELISM=serial

ipe:
	$(if $(filter 0,$(MAKELEVEL)),$(call USAGE))
	$(MAKE) -C main "$(VARS)" EXTRADEP="-d module_decomp.o"
	$(MAKE) -C neutral "$(VARS)"
	$(MAKE) -C flip "$(VARS)"

clean:
	$(MAKE) -C main clean
	$(MAKE) -C neutral clean
	$(MAKE) -C flip clean

#TEST=v18
#PACKAGE=../bin/ipe_${TEST}.$(COMPILER)$(PAR).exe
## DR0: modules
#DR0=${TOP}/main
#FFLAG0=${FFLAGS} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
## DR1: neutral
#DR1=${TOP}/neutral
#FFLAG1=${FFLAGS}  ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR3}
## DR2: FLIP
#DR2=${TOP}/flip
#FFLAG2=${FFLAGS}      ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR1} ${INCLUDE_PATH}${DR3} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
## DR3: plasma
#DR3=${TOP}/plasma
#FFLAG3=${FFLAGS}    ${INCLUDE_PATH}${DR0} ${INCLUDE_PATH}${DR6} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
## DR6: eldyn
#DR6=${TOP}/eldyn
#FFLAG6=${FFLAGS}    ${INCLUDE_PATH}${DR0}
## DR4: main driver
#DR4=${TOP}/main
#FFLAG4=${FFLAGS} ${INCLUDE_PATH}${DR1} ${INCLUDE_PATH}${DR3} ${INCLUDE_PATH}${DR6} ${GPTL_CPPFLAGS} ${GPTL_LDFLAGS}
##DR7: dummy gptl
#DR7=${TOP}/dummygptl
#FFLAGS7=
#
#SRC0S	=\
#	${DR0}/module_precision.f90 \
#	${DR0}/module_IPE_dimension.f90 \
#	${DR0}/module_physical_constants.f90 \
#	${DR0}/module_unit_conversion.f90 \
#	${DR0}/module_io.f90 \
#	${DR0}/module_open_file.f90 \
#	${DR0}/module_output.f90 \
#	${DR0}/module_close_files.f90
#SRC0P	=\
#	${DR0}/module_decomp.f90 \
#	${DR0}/module_input_parameters.f90 \
#	${DR0}/module_open_output_files.f90 \
#	${DR0}/module_field_line_grid.f90 \
#	${DR0}/allocate_arrays.f90 \
#	${DR0}/get_flip_grid.f90 \
#	${DR0}/module_read_plasma_grid_global.f90 \
#	${DR0}/module_init_plasma_grid.f90 \
#	${DR0}/output_plasma_grid.f90 \
#	${DR0}/stop.f90
## DR1: neutrals
#SRC1S	=\
#	${DR1}/nrlmsise00.f90 \
#	${DR1}/hwm93.f90 \
#	${DR1}/thermosphere_1d_v2.f90
#SRC1P	=\
#	${DR1}/module_neutral.f90
## DR2: FLIP
#SRC2	=\
#	${DR2}/CTIPE-int.f \
#	${DR2}/initialize_module_parameters.f \
#	${DR2}/CMINOR.f \
#	${DR2}/ELECXS.f \
#	${DR2}/INIT-PROFILES.f \
#	${DR2}/KEMPRN.f \
#	${DR2}/MINORA.f \
#	${DR2}/Neut_Heating.f \
#	${DR2}/Photoel-Freqs.f \
#	${DR2}/RSDENA_EXB.f \
#	${DR2}/RSJACA.f \
#	${DR2}/RSLPSD.f \
#	${DR2}/RSLPST.f \
#	${DR2}/RSPE2B.f \
#	${DR2}/RSPRIM.f \
#	${DR2}/RSTEMC_EXB.f \
#	${DR2}/Rates.f \
#	${DR2}/FLIP_GRID.f
## DR3: plasma
#SRC3	=\
#	${DR3}/module_plasma.f90 \
#	${DR3}/module_stepback_mag_R.f90 \
#	${DR3}/module_stepback_mag_th.f90 \
#	${DR3}/module_find_neighbor_grid_R.f90 \
#	${DR3}/module_find_neighbor_grid_th.f90 \
#	${DR3}/module_perpendicular_transport.f90 \
#	${DR3}/module_sub_plasma.f90 \
#	${DR3}/get_pvalue_dipole.f90 \
#	${DR3}/get_sza.f90 \
#	${DR3}/get_sinim.f90 \
#	${DR3}/module_heating_rate.f90 \
#	${DR3}/flux_tube_solver.f90 \
#	${DR3}/module_Qinterpolation.f90 \
#	${DR3}/interpolate_flux_tube_${TEST}.f90 \
#	${DR3}/io_plasma_bin.f90 \
#	${DR3}/read_vexb.f90
#
#
##dbg20120509	${DR3}/activate_perp_transport.f90 \
## DR6: eldyn
#SRC6S	=\
#	${DR6}/efield.f \
#	${DR6}/module_ff.f \
#	${DR6}/module_pnm.f \
#	${DR6}/module_prep_pnm.f \
#	${DR6}/module_index_quiet.f \
#	${DR6}/module_adj_S_a.f \
#	${DR6}/module_constants.f \
#	${DR6}/module_prep_fk.f \
#	${DR6}/module_set_fkflfs.f \
#	${DR6}/module_efield_mid.f \
#	${DR6}/module_FSVal.f \
#	${DR6}/module_SetModel.f \
#	${DR6}/module_ADJUST.f \
#	${DR6}/module_JULDAY.f \
#	${DR6}/module_GET_TILT.f \
#	${DR6}/module_prep_weimer.f \
#	${DR6}/module_pot_latsmo.f \
#	${DR6}/module_pot_latsmo2.f \
#	${DR6}/module_pot_lonsmo.f \
#	${DR6}/module_efield_sunloc.f \
#	${DR6}/module_fun_MLT.f \
#	${DR6}/module_fun_MagLong.f \
#	${DR6}/module_highlat_adjust.f \
#	${DR6}/module_interp_poten.f \
#	${DR6}/module_DerivPotential.f \
#	${DR6}/module_LEGENDRE.f \
#	${DR6}/module_EpotVal.f \
#	${DR6}/module_GECMP.f \
#	${DR6}/module_ReadCoef.f \
#	${DR6}/module_highlat_getbnd.f \
#	${DR6}/module_svbksb.f \
#	${DR6}/module_svdcmp.f \
#	${DR6}/module_bnd_sinus.f \
#	${DR6}/module_GlobalElPotential.f \
#	${DR6}/module_get_efield.f \
#	${DR6}/module_magfield.f \
#	${DR6}/module_sunloc.f \
#	${DR6}/apex_sunloc.f \
#	${DR6}/supim_exb.f
#SRC6P	=\
#	${DR6}/module_read_acoef.f \
#	${DR6}/module_efield_init.f \
#	${DR6}/module_eldyn.f      \
#	${DR6}/module_init_eldyn.f \
#       ${DR6}/module_cal_monthday.f \
#	${DR6}/module_sub_eldyn.f  \
#	${DR6}/get_efield90km.f
## DR4: main driver
#SRC4P	=\
#	${DR4}/driver_ipe.f90
#
##HEADS	= $(PACKAGE).h
#OBJ0S	= $(SRC0S:.f90=.o)
#OBJ0P	= $(SRC0P:.f90=.o)
#OBJ1S	= $(SRC1S:.f90=.o)
#OBJ1P	= $(SRC1P:.f90=.o)
#OBJ2	= $(SRC2:.f=.o)
#OBJ3	= $(SRC3:.f90=.o)
#OBJ4S	=
#OBJ4P	= $(SRC4P:.f90=.o)
#OBJ6S	= $(SRC6S:.f=.o)
#OBJ6P	= $(SRC6P:.f=.o)
##
#FILES	= README Makefile $(HEADS) $(SRC0S) $(SRC0P) $(SRC1S) $(SRC1P) $(SRC2) $(SRC3) $(SRC6S) $(SRC6p) $(SRC4P)
#VER	= `date +%Y%m%d`
#
#
#### to avoid m2c... !!!didn't work!!!
####%.o: %.mod
#
##all: $(PACKAGE)
##
##FINAL_OBJS=\
##	$(OBJ0S) \
##	$(OBJ0P) \
##	$(OBJ1S) \
##	$(OBJ1P) \
##	$(OBJ6S) \
##	$(OBJ6P) \
##	$(OBJ3) \
##	$(OBJ2) \
##	$(OBJ4P)
#
#$(PACKAGE):  D7_OBJ D0_OBJ D1_OBJ D6_OBJ D3_OBJ D2_OBJ D4_OBJ
#
#ifeq ($(PAR),sms)
#	which mpif90
#	$(FCparallel) $(LDFLAGS) $(FINAL_OBJS) -o $@ $(LDLIBS) $(GPTL_LDFLAGS) $(SMS_LDFLAGS)
#else
#	$(FC)         $(LDFLAGS) $(FINAL_OBJS) -o $@ $(LDLIBS) $(GPTL_LDFLAGS)
#endif
#
#$(FINAL_OBJS): $(HEADS) Makefile
#
#D0_OBJ :
#	cd ${DR0} ; ${MAKE} all \
#	"OBJS=${OBJ0S}" \
#	"OBJP=${OBJ0P}" \
#	"FC=${FC}" \
#	"FCparallel=${FCparallel}" \
#	"FFLAGS=${FFLAG0}" \
#	"LDFLAGS=${LDFLAGS}" \
#	"TOP=${TOP}" \
#	"PAR=${PAR}" \
#	"CPP=${CPP}" \
#	"CPP_FLAGS=${CPP_FLAGS}" \
#	"PPP=${PPP}" \
#	"PPP_FLAGS=${PPP_FLAGS}" \
#	"SMSFLAGS=${SMSFLAGS}" \
#D1_OBJ :
#	cd ${DR1} ; ${MAKE} all \
#	"OBJS=${OBJ1S}" \
#	"OBJP=${OBJ1P}" \
#	"FC=${FC}" \
#	"FCparallel=${FCparallel}" \
#	"FFLAGS=${FFLAG1}" \
#	"LDFLAGS=${LDFLAGS}" \
#	"TOP=${TOP}" \
#	"PAR=${PAR}" \
#	"CPP=${CPP}" \
#	"CPP_FLAGS=${CPP_FLAGS}" \
#	"PPP=${PPP}" \
#	"PPP_FLAGS=${PPP_FLAGS}" \
#	"SMSFLAGS=${SMSFLAGS}" \
#
#D2_OBJ :
#	cd ${DR2} ; ${MAKE} all \
#	"OBJ2=${OBJ2}" \
#	"FC=${FC}" \
#	"FFLAG2=${FFLAG2}" \
#	"LFLAG2=${LDFLAGS}" \
#	"TOP=${TOP}"
#D3_OBJ :
#	cd ${DR3} ; ${MAKE} all \
#	"OBJS=${OBJ3}" \
#	"FC=${FC}" \
#	"FCparallel=${FCparallel}" \
#	"FFLAGS=${FFLAG3}" \
#	"LDFLAGS=${LDFLAGS}" \
#	"TOP=${TOP}" \
#	"PAR=${PAR}" \
#	"CPP=${CPP}" \
#	"CPP_FLAGS=${CPP_FLAGS}" \
#	"PPP=${PPP}" \
#	"PPP_FLAGS=${PPP_FLAGS}" \
#	"SMSFLAGS=${SMSFLAGS}" \
#D6_OBJ :
#	cd ${DR6} ; ${MAKE} all \
#	"OBJS=${OBJ6S}" \
#	"OBJP=${OBJ6P}" \
#	"FC=${FC}" \
#	"FCparallel=${FCparallel}" \
#	"FFLAGS=${FFLAG6}" \
#	"LDFLAGS=${LDFLAGS}" \
#	"TOP=${TOP}" \
#	"PAR=${PAR}" \
#	"CPP=${CPP}" \
#	"CPP_FLAGS=${CPP_FLAGS}" \
#	"PPP=${PPP}" \
#	"PPP_FLAGS=${PPP_FLAGS}" \
#	"SMSFLAGS=${SMSFLAGS}" \
#D4_OBJ :
#	cd ${DR4} ; ${MAKE} all \
#	"OBJP=${OBJ4S}" \
#	"OBJP=${OBJ4P}" \
#	"FC=${FC}" \
#	"FCparallel=${FC}" \
#	"FFLAGS=${FFLAG4}" \
#	"LDFLAGS=${LDFLAGS}" \
#	"TOP=${TOP}" \
#	"PAR=${PAR}" \
#	"CPP=${CPP}" \
#	"CPP_FLAGS=${CPP_FLAGS}" \
#	"PPP=${PPP}" \
#	"PPP_FLAGS=${PPP_FLAGS}" \
#	"SMSFLAGS=${SMSFLAGS}" \
#D7_OBJ :
#	cd ${DR7} ; ${MAKE} "FC=${FC}" \
#
#
#### useful commands ###
#
#cln:
#	$(RM) $(FINAL_OBJS)
#	$(RM) ${DR0}/*.mod ${DR0}/*.tmp* ${DR0}/*.f ${DR0}/*~
#	$(RM) ${DR0}/*_sms.f90 ${DR0}/*ppptmp*
#	$(RM) ${DR1}/*.mod ${DR2}/*.mod ${DR3}/*.mod ${DR4}/*.mod ${DR6}/*.mod
#	$(RM) ${DR1}/*.mod ${DR1}/*.tmp* ${DR1}/*.f ${DR1}/*~
#	$(RM) ${DR1}/*_sms.f90 ${DR1}/*ppptmp*
#	$(RM) ${DR3}/*.mod ${DR3}/*.tmp* ${DR3}/*.f ${DR3}/*~ ${DR3}/module_decomp.* ${DR3}/module_input_parameters.*
#	$(RM) ${DR3}/*_sms.f90 ${DR3}/*ppptmp*
#	$(RM) ${DR6}/*.mod ${DR6}/*.tmp* ${DR6}/*~
#	$(RM) ${DR6}/*_sms.f90 ${DR6}/*ppptmp*
#	cd ${DR7} ; ${MAKE} clean
#tar:
#	@echo $(PACKAGE)-$(VER) > .package
#	@$(RM) -r `cat .package`
#	@mkdir `cat .package`
#	@ln $(FILES) `cat .package`
#	tar cvf - `cat .package` | gzip -9 > `cat .package`.tar.gz
#	@$(RM) -r `cat .package` .package
#
#zip:
#	zip -9 $(PACKAGE)-$(VER).zip $(FILES)
#
#
#prof: run
#	$(PROF) $(PACKAGE) | less
#
#run: all
#	./$(PACKAGE) < sample-data | less
#

