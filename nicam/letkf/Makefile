include ../../Mkinclude

MODS = letkf_obs.o letkf_tools.o letkf.o prg_icolatlon.o
PRJS = letkf icolatlon

#RTTOV_DIR = /scratch/ra000015/koji/rttov_O3/11.1
#FFLAGS_RTTOV = $(RTTOV_DIR)/src/test/amsua_fwd.F90 -L$(RTTOV_DIR)/lib -lrttov11.1.0_test -lrttov11.1.0_other -lrttov11.1.0_coef_io -lrttov11.1.0_parallel -lrttov11.1.0_main  -I$(RTTOV_DIR)/mod/ -I$(RTTOV_DIR)/include

all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

letkf: letkf.o
	$(FC) $(F90OPT) -o letkf letkf.o letkf_obs.o letkf_tools.o ../../common/SFMT.o ../../common/common.o ../../common/common_mpi.o ../../common/common_mtx.o ../../common/netlib2.o ../../common/common_letkf.o ../common/common_nicam.o ../common/common_obs_nicam.o ../common/common_mpi_nicam.o ../common/common_tvs_nicam.o $(NICAMDIR)/lib*.a
icolatlon: prg_icolatlon.o
	$(FC) $(F90OPT) -o icolatlon prg_icolatlon.o $(NICAMDIR)/lib*.a

	

#
letkf_obs.o   : letkf_obs.f90 ../../common/*.o ../common/*.o
letkf_tools.o : letkf_tools.f90 letkf_obs.o ../../common/*.o ../common/*.o
letkf.o       : letkf.f90 letkf_obs.o letkf_tools.o ../../common/*.o ../common/*.o
prg_icolatlon.o: prg_icolatlon.f90 

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
