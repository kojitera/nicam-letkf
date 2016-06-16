include ../../Mkinclude

MODS = common_nicam.o \
       common_mpi_nicam.o \
       common_obs_nicam.o \
       common_tvs_nicam.o 

all:
	make allmod

allmod: $(MODS)

#
common_nicam.o   : common_nicam.f90 ../../common/common.o
common_mpi_nicam.o   : common_mpi_nicam.f90 common_nicam.o ../../common/common_mpi.o 
common_obs_nicam.o   : common_obs_nicam.f90 common_nicam.o ../../common/common.o 
common_tvs_nicam.o   : common_tvs_nicam.f90 common_nicam.o ../../common/common.o 

.SUFFIXES: .o .f90 .c .f
.f.o:
	$(FC)  $(FFLAGS) -o $@ -c $<
.f90.o:
	$(FC)  $(F90OPT) -o $@ -c $<
.c.o:
	$(CC)  $(CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY: clean

clean:
	rm -f *.mod *.o *~ *.lst
