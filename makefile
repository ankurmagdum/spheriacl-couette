# set the fortran compiler
FC = gfortran

PGM = main

MAINF = $(wildcard src/*.f90)
OBJS = $(patsubst %.f90, %.o, $(MAINF))

# set the path to the libraries
LDFLAGS = -L/usr/local/lib -lshtns -lfftw3 -L/usr/lib -llapack

main: $(OBJS) main.o
	$(FC) -o $(PGM) main.o $(OBJS) $(LDFLAGS)

src/%.o: src/%.f90
	cd src; $(FC) -c $*.f90

main.o: main.f90
	$(FC) -I./src -c main.f90

clean:
	rm -r *.o *.dat							
