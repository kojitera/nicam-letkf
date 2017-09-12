include ../../Mkinclude

MODS = letkf_obs.o letkf_tools.o letkf.o 
PRJS = letkf

all: allmod allprg
allmod: $(MODS)
allprg: $(PRJS)
	$(INSTALL) $(PRJS) $(INSTALLDIR)

letkf: letkf.o
	$(FC) $(F90OPT) -o letkf letkf.o letkf_obs.o letkf_tools.o ../../common/SFMT.o ../../common/common.o ../../common/common_mpi.o ../../common/common_mtx.o ../../common/netlib2.o ../../common/common_letkf.o ../common/common_nicam.o ../common/common_obs_nicam.o ../common/common_mpi_nicam.o ../common/common_tvs_nicam.o $(NICAMDIR)/lib*.a

#
letkf_obs.o   : letkf_obs.f90 ../../common/*.o ../common/*.o
letkf_tools.o : letkf_tools.f90 letkf_obs.o ../../common/*.o ../common/*.o
letkf.o       : letkf.f90 letkf_obs.o letkf_tools.o ../../common/*.o ../common/*.o

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
