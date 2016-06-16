include ../../Mkinclude

MODS = prg_obsope_prepbufr.o prg_obsope_amsua.o prg_obsope_amsua_mean.o prg_vbc.o prg_makemean_landocn.o prg_makeproduct.o prg_obsope_amsua_one.o prg_estimateR.o prg_estimateR_mpi.o prg_estimateR_interchannel.o prg_icolatlon.o
PRJS = obsope_prepbufr obsope_amsua obsope_amsua_mean vbc  makemean_landocn makeproduct obsope_amsua_one estimateR estimateR_mpi estimateR_interchannel icolatlon

all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

compute_uv: prg_compute_uv.o
	$(FC) $(F90OPT) -o compute_uv prg_compute_uv.o $(NICAMDIR)/lib*.a

obsope_prepbufr: prg_obsope_prepbufr.o
	$(FC) $(F90OPT) -o obsope_prepbufr prg_obsope_prepbufr.o $(NICAMDIR)/lib*.a

makeproduct: prg_makeproduct.o
	$(FC) $(F90OPT) -o makeproduct prg_makeproduct.o $(NICAMDIR)/lib*.a

histgram: prg_histgram.o
	$(FC) $(F90OPT) -o histgram prg_histgram.o $(NICAMDIR)/lib*.a

obsope_amsua: prg_obsope_amsua.o
	$(FC) -o obsope_amsua prg_obsope_amsua.o ../common/common_tvs_nicam.o ../../common/common_mpi.o $(NICAMDIR)/lib*.a amsua_fwd.F90 $(F90OPT)

obsope_amsua_mean: prg_obsope_amsua_mean.o
	$(FC) -o obsope_amsua_mean prg_obsope_amsua_mean.o ../common/common_tvs_nicam.o ../../common/common_mpi.o $(NICAMDIR)/lib*.a amsua_fwd.F90 $(F90OPT)

obsope_amsua_one: prg_obsope_amsua_one.o
	$(FC) -o obsope_amsua_one prg_obsope_amsua_one.o ../common/common_tvs_nicam.o ../../common/common_mpi.o $(NICAMDIR)/lib*.a amsua_fwd.F90 $(F90OPT)

vbc: prg_vbc.o
	$(FC) -o vbc prg_vbc.o ../common/common_tvs_nicam.o ../../common/common_mpi.o ../../common/common_mtx.o ../../common/netlib2.o $(NICAMDIR)/lib*.a amsua_fwd.F90 $(F90OPT)

estimateR: prg_estimateR.o
	$(FC) -o estimateR prg_estimateR.o ../common/common_tvs_nicam.o ../../common/common_mpi.o ../common/common_nicam.o $(NICAMDIR)/lib*.a $(F90OPT)

estimateR_mpi: prg_estimateR_mpi.o
	$(FC) -o estimateR_mpi prg_estimateR_mpi.o ../common/common_tvs_nicam.o ../../common/common_mpi.o ../common/common_nicam.o $(NICAMDIR)/lib*.a $(F90OPT)

estimateR_interchannel: prg_estimateR_interchannel.o
	$(FC) -o estimateR_interchannel prg_estimateR_interchannel.o ../common/common_tvs_nicam.o ../../common/common_mpi.o ../common/common_nicam.o $(NICAMDIR)/lib*.a $(F90OPT)

makemean_landocn: prg_makemean_landocn.o
	$(FC) $(F90OPT) -o makemean_landocn prg_makemean_landocn.o $(NICAMDIR)/lib*.a

icolatlon: prg_icolatlon.o
	$(FC) $(F90OPT) -o icolatlon prg_icolatlon.o $(NICAMDIR)/lib*.a

#
prg_compute_uv.o: prg_compute_uv.f90
prg_obsope_amsua.o: prg_obsope_amsua.f90 ../common/common_tvs_nicam.f90 ../../common/common_mpi.f90 amsua_fwd.F90
prg_vbc.o: prg_vbc.f90 ../common/common_tvs_nicam.f90 ../../common/common_mpi.f90 ../../common/common_mtx.f90 amsua_fwd.F90 

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
