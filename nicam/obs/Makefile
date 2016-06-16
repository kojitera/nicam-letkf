#include ../../Mkinclude

INSTALL=install
INSTALLDIR=../../bin/

FC=mpif90
F90OPT=-O3 -assume byterecl -convert big_endian

MODS = mod_calendar.o prg_read_amsua_thinning.o prg_dec_prepbufr.o
PRJS = read_amsua_thinning read_prepbufr_thinning dec_prepbufr

#all: allmod
all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

read_amsua_thinning: prg_read_amsua_thinning.o
	$(FC) $(F90OPT) -o read_amsua_thinning prg_read_amsua_thinning.o mod_calendar.o -L/volume63/data/ra000015/koji/lib/bufrlib_ppb -lbufrlib

read_prepbufr_thinning: prg_read_prepbufr_thinning.o
	$(FC) $(F90OPT) -o read_prepbufr_thinning prg_read_prepbufr_thinning.o mod_calendar.o

dec_prepbufr: prg_dec_prepbufr.f90
	$(FC) $(F90OPT) -o dec_prepbufr prg_dec_prepbufr.f90 -L/volume63/data/ra000015/koji/lib/bufrlib_ppb -lbufrlib

#
prg_read_amsua_thinning.o: prg_read_amsua_thinning.f90 mod_calendar.o

.SUFFIXES:
.SUFFIXES: .o .f90 .f
.f.o:
	$(FC)  $(F90OPT) -o $@ -c $<
.f90.o:
	$(FC)  $(F90OPT) -o $@ -c $<

#%.o: %.mod

.PHONY: clean

clean:
	rm -f *.mod *.o *~ *.lst $(PRJS)
