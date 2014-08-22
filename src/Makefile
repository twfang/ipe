# IPE Makefile

# Static/Unconditional/Default Variables

DEBUG         = no
DIRS          = main neutral flip eldyn plasma dummygptl
GPTL_INCFLAGS = -I../dummygptl
GPTL_LDFLAGS  = -L../dummygptl -lgptl
GPTL_USE_GPTL = no
MAKETHREADS   = 8
PPP_FLAGS     = --Fmodule=module_decomp --Free --comment --HaloSize=10
PROF          = gprof
TARGETS       = zeus_intel_parallel zeus_intel_serial zeus_lahey_serial clean
TEST          = v18

# Dynamic/Conditional Variables

EXE           = ipe_${TEST}_$(MACHINE)_$(COMPILER)_$(PARALLELISM).exe
FFLAGS        = $(OPT_FLAG) $(DEBUG_FLAGS) $(GPTL_INCFLAGS) $(GPTL_LDFLAGS)
PARALLEL_MAKE = $(MAKE) -j $(MAKETHREADS)
PPP           = $(SMS)/bin/ppp
SMS_INCFLAGS  = -I$(SMS)/include
SMS_LDFLAGS   = -L$(SMS)/lib -lsms

ifeq ($(MACHINE),zeus)
  ifeq ($(COMPILER),intel)
    ENVSET   = module load intel mpt
    FCP      = mpif90
    FCS      = ifort
    OPT_FLAG = -O3 -fp-model precise
    SMS      = /contrib/acs/SMS/legacy_r616/zeus-intel-mpt
    ifeq ($(DEBUG),yes)
      DEBUG_FLAGS = -g -O0 -ftz -traceback -check bounds -fpe0
    endif
    ifeq ($(PARALLELISM),parallel)
#     ifeq ($(GPTL_USE_GPTL),yes)
#     endif
    endif
    ifeq ($(PARALLELISM),serial)
#     ifeq ($(GPTL_USE_GPTL),yes)
#     endif
    endif
  endif
  ifeq ($(COMPILER),lahey)
    ENVSET = module load lahey mpt
    FCP    = mpif90
    FCS    = lf95
    SMS    = /contrib/acs/SMS/legacy_r616/zeus-lahey-mpt
    ifeq ($(PARALLELISM),parallel)
    endif
    ifeq ($(PARALLELISM),serial)
    endif
  endif
endif

VARS="DIRS="$(DIRS)" ENVSET="$(ENVSET)" FCP="$(FCP)" FCS="$(FCS)" FFLAGS="$(FFLAGS)" PARALLELISM="$(PARALLELISM)" PPP="$(PPP)" PPP_FLAGS="$(PPP_FLAGS)" SMS_INCFLAGS="$(SMS_INCFLAGS)""

# Functions

USAGE=$(error Valid targets: $(TARGETS))

# Targets

.PHONY: clean

all:
	$(call $(USAGE))

zeus_intel_parallel:
	$(PARALLEL_MAKE) ipe MACHINE=zeus COMPILER=intel PARALLELISM=parallel

zeus_intel_serial:
	$(PARALLEL_MAKE) ipe MACHINE=zeus COMPILER=intel PARALLELISM=serial

zeus_lahey_serial:
	$(PARALLEL_MAKE) ipe MACHINE=zeus COMPILER=lahey PARALLELISM=serial

ipe:
	$(if $(filter 0,$(MAKELEVEL)),$(call USAGE))
	$(MAKE) -C main      "$(VARS)" EXTRADEP="-d module_decomp.o"
	$(MAKE) -C neutral   "$(VARS)"
	$(MAKE) -C flip      "$(VARS)"
	$(MAKE) -C eldyn     "$(VARS)"
	$(MAKE) -C plasma    "$(VARS)" TEST="$(TEST)"
	$(MAKE) -C dummygptl "$(VARS)"
	$(MAKE) -C main bin  "$(VARS)" EXE="$(EXE)" SMS_LDFLAGS="$(SMS_LDFLAGS)"

clean:
	$(foreach x,$(DIRS),$(MAKE) -C $(x) clean;)

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

