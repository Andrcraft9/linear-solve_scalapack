LIB_COMPILE = -I/opt/ibmmath/essl/6.1/include
LIB_LINK = -L/opt/ibmmath/pessl/5.3/lib64 -L/opt/ibmmath/essl/6.1/lib64
FLAGS = -qsmp=omp 

make:
	mpifort -c tools.f pmatgeninc.f pdmatgen.f $(LIB_COMPILE)
	mpifort -c linearsolve.f90 $(LIB_COMPILE)
	mpifort $(FLAGS) -o prog linearsolve.o pmatgeninc.o pdmatgen.o tools.o $(LIB_LINK) -lessl6464 -lblacssmp -lblas -lpesslsmp 

clean: 
	rm *.o prog
