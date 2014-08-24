# driver Makefile

TOP = ..

PSRCS = driver_ipe.f90

OBJS = $(wildcard $(addprefix ../,$(addsuffix /*.o,$(DIRS))))

.PHONY: exe

exe: $(EXE)

include $(TOP)/Makefile.common

$(EXE): driver_ipe.o $(OBJS)
ifeq ($(PARALLELISM),parallel)
	$(ENVSET) && $(FCP) $(FFLAGS) $(INCS) -o $@ $^ $(SMS_INCFLAGS) $(SMS_LDFLAGS)
else
	$(ENVSET) && $(FCS) $(FFLAGS) $(INCS) -o $@ $^
endif
