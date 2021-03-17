# Spherical Gap Flow

This program simulates incompressible flow in a spherica gap based on a pseudo-spectral technique developed by G. Dumas. The technique relies on the use of special  divergence-free test functions that vanish on the boundary to eliminate pressure from the set of unknown variables. The time-evolution is performed in spectral space which has fewer degrees of freedom than the physical space. Nevertheless, every time-step requires transforming back to physical space but this cab be performed efficiently by the Fourier Transform routine FFTW3.

The program is implemented with FORTRAN in the file __main.f90__ and relies on the following libraries:
* FFTW3 (<http://www.fftw.org/>)
* SHTns (<https://bitbucket.org/nschaeff/shtns/src/master/>)
* lapack (<https://www.netlib.org/lapack/>)

The program is compiled with the command `make` and generates the executable __main__ which can be run with the command `./main`. Although, before compliling, one must ensure that the path to the libraries in the makefile are set correctly.

A report of the theoretical basis and the implemntation can be found in the file __report.pdf__.
