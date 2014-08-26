# IPE Makefile

# Static / Unconditional / Default variables

DIRS          = main neutral flip eldyn plasma dummygptl driver
GPTL_INCFLAGS = -I../dummygptl
GPTL_LDFLAGS  = -L../dummygptl -lgptl
PPP_FLAGS     = --Fmodule=module_decomp --Free --comment --HaloSize=10
TEST          = v18

#TARGETS = clean                  \
#         trillian_cray_parallel \
#         trillian_cray_serial   \
#         trillian_pgi_parallel  \
#         trillian_pgi_serial    \
#         zeus_intel_parallel    \
#         zeus_intel_serial      \
#         zeus_lahey_serial
TARGETS = clean                  \
          trillian_pgi_parallel  \
          trillian_pgi_serial    \
          zeus_intel_parallel    \
          zeus_intel_serial      \
          zeus_lahey_serial

# Dynamic / Conditional variables

EXE           = $(CURDIR)/../bin/ipe_${TEST}_$(MACHINE)_$(COMPILER)_$(PARALLELISM).exe
FFLAGS        = $(strip $(OPT_FLAGS) $(DEBUG_FLAGS) $(GPTL_INCFLAGS) $(GPTL_LDFLAGS))
PARALLEL_MAKE = $(MAKE) -j $(MAKEJOBS)
PPP           = $(SMS)/bin/ppp
SMS_INCFLAGS  = -I$(SMS)/include
SMS_LDFLAGS   = -L$(SMS)/lib -lsms

ifndef DEBUG
  DEBUG = no
endif

# trilian stanza

ifeq ($(MACHINE),trillian)
  GPTL_USE_GPTL = no
  MAKEJOBS      = 6
  ifeq ($(GPTL_USE_GPTL),yes)
    $(error GPTL support currently unavailable on trillian)
  endif
# ifeq ($(COMPILER),cray)
#   MODCMD    = $(shell ./modcmd trillian cray)
#   FCP       = ftn
#   FCS       = ftn
#   OPT_FLAGS = -O3 -em -ef -eZ -ra -h noomp
#   SMS       = /mnt/lustre/lus0/space/madden/SMS/legacy_r618/trillian-cray
#   ifeq ($(DEBUG),yes)
#     DEBUG_FLAGS = -g -O0
#   endif
# endif
  ifeq ($(COMPILER),pgi)
    MODCMD    = $(shell ./modcmd trillian pgi)
    FCP       = ftn
    FCS       = ftn
    OPT_FLAGS = -O2
    SMS       = /mnt/lustre/lus0/space/madden/SMS/legacy_r618/trillian-pgi
    ifeq ($(DEBUG),yes)
      DEBUG_FLAGS = -g -O0 -Mbounds
    endif
  endif
endif

# zeus stanza

ifeq ($(MACHINE),zeus)
  GPTL_USE_GPTL = yes
  MAKEJOBS      = 8
  ifeq ($(COMPILER),intel)
    MODCMD    = $(shell ./modcmd zeus intel)
    FCP       = mpif90
    FCS       = ifort
    OPT_FLAGS = -O3 -fp-model precise
    SMS       = /contrib/acs/SMS/legacy_r617/zeus-intel-mpt
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
    MODCMD = $(shell ./modcmd zeus lahey)
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
       EXE="$(EXE)"                   \
       FCP="$(FCP)"                   \
       FCS="$(FCS)"                   \
       FFLAGS="$(FFLAGS)"             \
       MODCMD="$(MODCMD)"             \
       PARALLELISM="$(PARALLELISM)"   \
       PPP="$(PPP)"                   \
       PPP_FLAGS="$(PPP_FLAGS)"       \
       SMS_INCFLAGS="$(SMS_INCFLAGS)" \
       SMS_LDFLAGS="$(SMS_LDFLAGS)"   \
       TEST="$(TEST)"

# Targets

.PHONY: $(TARGETS) ipe

all:
	$(info Valid targets are:)
	$(foreach x,$(TARGETS),$(info - $(x)))
	$(error Invalid target specified)

#trillian_cray_parallel:
#	$(MAKE) ipe MACHINE=trillian COMPILER=cray PARALLELISM=parallel
#
#trillian_cray_serial:
#	$(MAKE) ipe MACHINE=trillian COMPILER=cray PARALLELISM=serial

trillian_pgi_parallel:
	$(MAKE) ipe MACHINE=trillian COMPILER=pgi PARALLELISM=parallel

trillian_pgi_serial:
	$(MAKE) ipe MACHINE=trillian COMPILER=pgi PARALLELISM=serial

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
