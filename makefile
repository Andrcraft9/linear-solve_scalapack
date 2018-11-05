LIB_COMPILE = -I/opt/ibmmath/essl/6.1/include
LIB_LINK = -L/opt/ibmmath/pessl/5.3/lib64 -L/opt/ibmmath/essl/6.1/lib64
FLAGS = -lmpi_ibm -lmpi_mpifh

FORTRAN_LIB = -L/opt/ibm/xlf/15.1.6/lib/ -L/opt/ibm/xlsmp/4.1.6/lib/ -lxlf90_r -lxl -lxlsmp -lxlfmath -lm 
SCALAPACK_LIB = -lscalapack -llapack
BLAS_BLACS_LIB = -L/opt/ibmmath/essl/6.1/lib64/ -lesslsmp -L/opt/ibmmath/pessl/5.3/lib64 -lblacssmp 
#ALL_LIB = $(SCALAPACK_LIB) $(BLAS_BLACS_LIB) $(FORTRAN_LIB) 
ALL_LIB = $(SCALAPACK_LIB)

ibm:
	mpifort -c tools.f pmatgeninc.f pdmatgen.f $(LIB_COMPILE)
	mpifort -c linearsolve.f90 $(LIB_COMPILE)
	mpifort $(FLAGS) -o prog linearsolve.o pmatgeninc.o pdmatgen.o tools.o $(LIB_LINK) -lesslsmp -lblacssmp -lpesslsmp 

openmpi:
	mpifort -fopenmp -c tools.f pmatgeninc.f pdmatgen.f
	mpifort -fopenmp -c linearsolve.f90
	mpifort -fopenmp -o prog linearsolve.o pmatgeninc.o pdmatgen.o tools.o $(ALL_LIB)
clean: 
	rm *.o prog
