include ../Mkinclude

#$FC = mpiifort
#F90OPT= -xHOST -O3 -convert big_endian -assume byterecl
#CC = mpicc
#LD = $(FC)
#LFLAGS = $(F90OPT)
#LFLAGS_OMP = $(F90OPT) -openmp

MODS = common.o \
       common_mpi.o \
       common_mtx.o \
       common_letkf.o \
       SFMT.o   \
       netlib2.o \
       netlibblas.o \
       netlibrs.o \
       minimizelib.o \
       lbfgs.o

all:
	make allmod

allmod: $(MODS)

#
common.o       : common.f90
common_mpi.o   : common_mpi.f90
common_mtx.o   : common_mtx.f90
common_letkf.o : common_letkf.f90 common.o common_mtx.o
SFMT.o         : SFMT.f90
netlib2.o      : netlib2.f
netlibblas.o   : netlibblas.f
netlibrs.o     : netlibrs.f
mt95.o         : mt95.f90
mt19937ar.o    : mt19937ar.f90
minimizelib.o  : minimizelib.f90
mbfgs.o        : lbfgs.f

.SUFFIXES: .o .f90 .c .f
.f.o:
	$(FC)  $(F90OPT) -o $@ -c $<
.f90.o:
	$(FC)  $(F90OPT) -o $@ -c $<
.c.o:
	$(CC)  $(CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY: clean

clean:
	rm -f *.mod *.o *~ *.lst ./mod_oprt.f90 ./mod_oprt3d.f90
