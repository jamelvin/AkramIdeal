# AkramIdeal

# copy this code at the bottom of the Makefile after setup
# note you may need to fix the spacing for the indent lines
# by deleting the spaces and putting in a tab

ifeq ($(FLASHBINARY),true)
iso_c_binding.mod :
        touch $@
#gcc version 4.9.1 results in MPI communication errors
#unless we compile with -O0
mpi_amr_1blk_guardcell.o : %.o : %.F90
        $(FCOMP) $(FFLAGS) -O0 $(F90FLAGS) $(FDEFINES)  $<
endif
