# IPE Makefile

# Static / Unconditional / Default variables

DIRS          = main neutral flip eldyn plasma dummygptl driver
GPTL_INCFLAGS = -I../dummygptl
GPTL_LDFLAGS  = -L../dummygptl -lgptl
GPTL_USE_GPTL = yes
MAKEJOBS      = 8
PPP_FLAGS     = --Fmodule=module_decomp --Free --comment --HaloSize=10
TARGETS       = clean zeus_intel_parallel zeus_intel_serial zeus_lahey_serial
TEST          = v18

# Dynamic / Conditional variables

EXE           = $(CURDIR)/../bin/ipe_${TEST}_$(MACHINE)_$(COMPILER)_$(PARALLELISM).exe
FFLAGS        = $(strip $(OPT_FLAG) $(DEBUG_FLAGS) $(GPTL_INCFLAGS) $(GPTL_LDFLAGS))
PARALLEL_MAKE = $(MAKE) -j $(MAKEJOBS)
PPP           = $(SMS)/bin/ppp
SMS_INCFLAGS  = -I$(SMS)/include
SMS_LDFLAGS   = -L$(SMS)/lib -lsms

ifndef DEBUG
  DEBUG = no
endif

ifeq ($(MACHINE),zeus)
  ifeq ($(COMPILER),intel)
    ENVSET   = module load intel mpt
    FCP      = mpif90
    FCS      = ifort
    OPT_FLAG = -O3 -fp-model precise
    SMS      = /contrib/acs/SMS/legacy_r617/zeus-intel-mpt
    ifeq ($(DEBUG),yes)
      DEBUG_FLAGS = -g -O0 -ftz -traceback -check bounds -fpe0
    endif
    ifeq ($(PARALLELISM),parallel)
      ifeq ($(GPTL_USE_GPTL),yes)
        GPTL_INCFLAGS = -I/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_mpi/include
        GPTL_LDFLAGS  = -L/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_mpi/lib -lgptl
      endif
    endif
    ifeq ($(PARALLELISM),serial)
      ifeq ($(GPTL_USE_GPTL),yes)
        GPTL_INCFLAGS = -I/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_nompi/include
        GPTL_LDFLAGS  = -L/scratch1/portfolios/BMC/fim/GPTL/gptl-v5.0/intel_nompi/lib -lgptl
      endif
    endif
  endif
  ifeq ($(COMPILER),lahey)
    ENVSET = module load lahey mpt
    FCP    = mpif90
    FCS    = lf95
    SMS    = /contrib/acs/SMS/legacy_r617/zeus-lahey-mpt
    ifeq ($(DEBUG),yes)
      DEBUG_FLAGS = --chk -g --trap --trace --info
    endif
  endif
endif

VARS =                                \
       DIRS="$(DIRS)"                 \
       ENVSET="$(ENVSET)"             \
       EXE="$(EXE)"                   \
       FCP="$(FCP)"                   \
       FCS="$(FCS)"                   \
       FFLAGS="$(FFLAGS)"             \
       PARALLELISM="$(PARALLELISM)"   \
       PPP="$(PPP)"                   \
       PPP_FLAGS="$(PPP_FLAGS)"       \
       SMS_INCFLAGS="$(SMS_INCFLAGS)" \
       SMS_LDFLAGS="$(SMS_LDFLAGS)"   \
       TEST="$(TEST)"

# Functions

USAGE=$(error Valid targets: $(TARGETS))

# Targets

.PHONY: $(TARGETS) ipe

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
	$(foreach x,$(DIRS),$(PARALLEL_MAKE) -C $(x) $(VARS) &&) true

clean:
	$(foreach x,$(DIRS),$(MAKE) -C $(x) clean &&) true
