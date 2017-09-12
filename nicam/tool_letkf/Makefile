include ../../Mkinclude

MODS = prg_makeproduct_forecast.o prg_spread_uv.o prg_makemean_landocn.o prg_makeproduct.o prg_estimateR.o prg_estimateR_mpi.o prg_estimateR_interchannel.o prg_icolatlon.o
PRJS = makeproduct_forecast spread_uv makemean_landocn makeproduct estimateR estimateR_mpi estimateR_interchannel icolatlon

all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

spread_uv: prg_spread_uv.o
	$(FC) $(F90OPT) -o spread_uv prg_spread_uv.o $(NICAMDIR)/lib*.a

compute_uv: prg_compute_uv.o
	$(FC) $(F90OPT) -o compute_uv prg_compute_uv.o $(NICAMDIR)/lib*.a

makeproduct: prg_makeproduct.o
	$(FC) $(F90OPT) -o makeproduct prg_makeproduct.o $(NICAMDIR)/lib*.a

makeproduct_forecast: prg_makeproduct_forecast.o
	$(FC) $(F90OPT) -o makeproduct_forecast prg_makeproduct_forecast.o $(NICAMDIR)/lib*.a

histgram: prg_histgram.o
	$(FC) $(F90OPT) -o histgram prg_histgram.o $(NICAMDIR)/lib*.a

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

cdfmake: prg_cdfmake.o
	$(FC) $(F90OPT) -o cdfmake    prg_cdfmake.o    fnc_transgauss.o $(NICAMDIR)/lib*.a
dstat: prg_dstat.o
	$(FC) $(F90OPT) -o dstat      prg_dstat.o      fnc_transgauss.o $(NICAMDIR)/lib*.a
rain2trans: prg_rain2trans.o
	$(FC) $(F90OPT) -o rain2trans prg_rain2trans.o fnc_transgauss.o $(NICAMDIR)/lib*.a ../../common/SFMT.o
trans2rain: prg_trans2rain.o
	$(FC) $(F90OPT) -o trans2rain prg_trans2rain.o fnc_transgauss.o $(NICAMDIR)/lib*.a ../../common/SFMT.o
#
prg_compute_uv.o: prg_compute_uv.f90

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
