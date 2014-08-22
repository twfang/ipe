# main Makefile

TOP=..

PSRCS=allocate_arrays.f90                \
      get_flip_grid.f90                  \
      module_decomp.f90                  \
      module_field_line_grid.f90         \
      module_init_plasma_grid.f90        \
      module_input_parameters.f90        \
      module_open_output_files.f90       \
      module_read_plasma_grid_global.f90 \
      output_plasma_grid.f90             \
      stop.f90

SSRCS=module_close_files.f90             \
      module_io.f90                      \
      module_IPE_dimension.f90           \
      module_open_file.f90               \
      module_output.f90                  \
      module_physical_constants.f90      \
      module_precision.f90               \
      module_unit_conversion.f90

include $(TOP)/common

BIN  = $(TOP)/../bin/$(EXE)
OBJS = $(wildcard $(addprefix ../,$(addsuffix /*.o,$(DIRS))))

bin: $(BIN)

$(BIN): $(OBJS) driver_ipe.f90
ifeq ($(PARALLELISM),parallel)
	$(ENVSET) && $(FCP) $(FFLAGS) $(INCS) -o $@ $^ $(SMS_INCFLAGS) $(SMS_LDFLAGS)
else
	$(ENVSET) && $(FCS) $(FFLAGS) $(INCS) -o $@ $^
endif

module_decomp.o: module_decomp.f90
ifeq ($(PARALLELISM),parallel)
	$(PPP) $(filter-out --Fmodule=module_decomp,$(PPP_FLAGS)) -o $*_sms.f90 $<
	$(ENVSET) && $(FCP) -c $(FFLAGS) $(SMS_INCFLAGS) -o $@ $*_sms.f90
else
	$(ENVSET) && $(FCS) -c $(FFLAGS) -o $@ $<
endif
