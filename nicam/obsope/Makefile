include ../../Mkinclude

MODS = amsua_fwd.o atms_fwd.o iasi_fwd.o mod_obsope_common.o mod_read_history.o mod_read_firstguess.o mod_scanbias.o mod_vbc.o mod_obsope_prepbufr.o mod_obsope_amsua.o mod_obsope_atms.o mod_obsope_gsmap.o mod_obsope_iasi.o mod_obsope_cris.o mod_obsope_airs.o mod_obsope_driver.o prg_obsope.o prg_vbc.o
PRJS = amsua_fwd.o atms_fwd.o iasi_fwd.o obsope  vbc

all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

amsua_fwd.o: amsua_fwd.F90
	$(FC) -c amsua_fwd.F90 $(F90OPT)

atms_fwd.o: atms_fwd.F90
	$(FC) -c atms_fwd.F90 $(F90OPT)

iasi_fwd.o: iasi_fwd.F90
	$(FC) -c iasi_fwd.F90 $(F90OPT)

obsope: prg_obsope.o
	$(FC) -o obsope prg_obsope.o mod_*.o amsua_fwd.o atms_fwd.o iasi_fwd.o ../../common/common_mtx.o ../../common/netlib2.o $(NICAMDIR)/lib*.a $(F90OPT)

vbc: prg_vbc.o
	$(FC) -o vbc prg_vbc.o mod_*.o amsua_fwd.o atms_fwd.o iasi_fwd.o ../../common/common_mtx.o ../../common/netlib2.o $(NICAMDIR)/lib*.a $(F90OPT)

mod_obsope_common.o: mod_obsope_common.f90
mod_read_firstguess.o: mod_read_firstguess.f90
mod_read_history.o: mod_read_history.f90
mod_scanbias.o: mod_scanbias.f90
mod_vbc.o: mod_vbc.f90
mod_obsope_prepbufr.o: mod_obsope_prepbufr.f90
mod_obsope_amsua.o: mod_obsope_amsua.f90 amsua_fwd.F90
mod_obsope_atms.o: mod_obsope_atms.f90 atms_fwd.F90
mod_obsope_gsmap.o: mod_obsope_gsmap.f90
mod_obsope_iasi.o: mod_obsope_iasi.f90 iasi_fwd.F90
mod_obsope_cris.o: mod_obsope_cris.f90
mod_obsope_airs.o: mod_obsope_airs.f90
mod_obsope_driver.o: mod_obsope_driver.f90 mod_obsope_common.o mod_obsope_prepbufr.o mod_obsope_amsua.o mod_obsope_atms.o mod_obsope_gsmap.o mod_obsope_iasi.o mod_obsope_cris.o mod_obsope_airs.o mod_scanbias.o mod_vbc.o 
prg_obsope.o: prg_obsope.f90 mod_obsope_driver.o mod_read_firstguess.o mod_scanbias.o 
prg_vbc.o: prg_obsope.f90 mod_obsope_driver.o mod_read_history.o mod_scanbias.o
#

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .f
.f.o:
	$(FC)  $(F90OPT) -o $@ -c $<
.f90.o:
	$(FC)  $(F90OPT) -o $@ -c $<

#%.o: %.mod

.PHONY: clean

clean:
	rm -f *.mod *.o *~ *.lst $(PRJS)
