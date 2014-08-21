# main Makefile

TOP=..

PBASES=allocate_arrays                \
       get_flip_grid                  \
       module_decomp                  \
       module_field_line_grid         \
       module_init_plasma_grid        \
       module_input_parameters        \
       module_open_output_files       \
       module_read_plasma_grid_global \
       output_plasma_grid stop

SBASES=module_close_files             \
       module_io                      \
       module_IPE_dimension           \
       module_open_file               \
       module_output                  \
       module_physical_constants      \
       module_precision               \
       module_unit_conversion

include $(TOP)/common

module_decomp.o: module_decomp.f90
ifeq ($(PARALLELISM),parallel)
	$(PPP) $(filter-out --Fmodule=module_decomp,$(PPP_FLAGS)) -o $*_sms.f90 $<
	$(ENVSET) && $(FCP) -c $(FFLAGS) $(SMSFLAGS) -o $@ $*_sms.f90
else
	$(ENVSET) && $(FCS) -c $(FFLAGS) -o $@ $<
endif
