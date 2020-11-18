#

.SUFFIXES: .f .f90

FC=/opt/intel/bin/ifort
#FC_ALL = 
#FC = ifort

#FFLAGS = -g -traceback
FFLAGS = -check all -g -traceback
#FFLAGS = -O2 -c -parallel  
#FLAGS = -O2 -fast -qopenmp -c -parallel #-ipo -par-threshold10 
# O2 paraleliza corretamente

LIB_BLAS   = -L/opt/intel/mkl/lib/intel64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -qopt-matmul
LIB_LAPACK = -L/opt/intel/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -qopt-matmul

INCS_MKL   = -I/opt/intel/mkl/include/intel64/lp64



LIB  = $(LIB_BLAS) $(LIB_LAPACK)
INCS = $(INCS_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------
SOURCE = types.o\
       constants.o\
     parameters.o\
     functions.o\
     overlap.o\
     system_hamiltonian.o\
     verlet.o\
     rdftensor.o\
          ode.o\
	 gnugraphs.o\
     time_evolution.o\
       main.o\




a: $(SOURCE)
	-rm -f a
	$(FC) $(FC_ALL) $(INCS) -o a $(SOURCE) $(LIB)
	-rm -f *.log
.f.o:
	$(FC) -fpp -free $(FC_ALL) $(FFLAGS) $(INCS) -c $<
.f90.o:
	$(FC) -fpp -free $(FC_ALL) $(FFLAGS) $(INCS) -c $<
clean:
	-rm -f *.o *.mod; touch *.f *.f90;
#safe: FC_ALL += -check all -traceback -fstack-protector -assume protect_parens
#safe: a 	
