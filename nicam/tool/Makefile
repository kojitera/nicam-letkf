include ../../Mkinclude

MODS = prg_icormse.o prg_increment.o z2p.o prg_llrmse.o prg_make_gph.f90
PRJS = icormse increment z2p llrmse make_gph

all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

icormse: prg_icormse.o
	$(FC) $(F90OPT) -o icormse prg_icormse.o ../../common/common.o ../../common/SFMT.o
increment: prg_increment.o
	$(FC) $(F90OPT) -o increment prg_increment.o
z2p: z2p.o
	$(FC) $(F90OPT) -o z2p z2p.o
llrmse: prg_llrmse.o
	$(FC) $(F90OPT) -o llrmse prg_llrmse.o ../../common/common.o ../../common/SFMT.o

make_gph: prg_make_gph.o
	$(FC) $(F90OPT) -o make_gph prg_make_gph.o
#
prg_icormse.o: prg_icormse.f90 ../../common/common.o ../../common/SFMT.o
prg_increment.o: prg_increment.f90
z2p.o: z2p.f90
prg_llrmse.o: prg_llrmse.f90 ../../common/common.o ../../common/SFMT.o
prg_make_gph.o: prg_make_gph.f90

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
