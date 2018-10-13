LIB_DIR = /home/andr/lib/

make:
	mpif90 -c pmatgeninc.f pdmatgen.f
	mpif90 -c linearsolve.f90
	mpif90 -o prog linearsolve.o pmatgeninc.o pdmatgen.o -L $(LIB_DIR) -lscalapack -llapack -lblas

clean: 
	rm *.o prog
