main: para.o main_sub.o func.o sub.o
	gfortran -g -O3 -fdefault-real-8 -o main para.o main_sub.o func.o sub.o -L/usr/local/lib -lshtns -lfftw3 -L/usr/lib -llapack
para.mod: para.o para.f90
	gfortran -c para.f90
para.o: para.f90
	gfortran -c para.f90	
func.mod: para.mod func.o func.f90
	gfortran -c func.f90
func.o: para.mod func.f90
	gfortran -c func.f90
sub.mod: para.mod sub.o sub.f90
	gfortran -c sub.f90
sub.o: para.mod sub.f90
	gfortran -c sub.f90
main_sub.o: para.mod func.mod sub.mod main_sub.f90
	gfortran -c main_sub.f90
clean:
	rm para.mod func.mod sub.mod para.o func.o sub.o main_sub.o							
