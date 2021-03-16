      module parameters 
      
      implicit none
      
      integer :: plan, nlmd, layout
      integer, parameter       :: n=32, l=36, m=0, md=1.50d0*m+1, mres=1
      integer, parameter       :: nd=1.5d0*(n+4) + 1, ld=1.5d0*(l+1) + 2
      real*8, parameter        :: pi=3.1415926d0, re=50.d0
      real*8, parameter        :: reinv=1.d0/re, radi=1.d0, rado=2.d0
      real*8, parameter        :: km=.5d0*(radi+rado) 
      complex*8, parameter     :: iota = (0.d0,1.d0)
      
      real*8, dimension(0:nd)    :: rad, c1
      real*8, dimension(-7:n+5)  :: c, d
      real*8, dimension(-5:n+5)  :: e, f, g, p, q, r
      real*8, dimension(0:n,-5:5):: b1, b13, b20, b21, b23, b24, b25, &
     &                              b26, b27, b28, b32, b33, b34, b35,&
     &                              b36, b37, b16, b22 
      real*8, dimension(-5:n+5,-2:2)   :: auxb14, auxb24, auxb25
      
      complex*8, dimension(0:n,0:m)    :: vbc
      complex*8, dimension(1:l,0:n,0:n):: t1m, t1p, t2m, t2p
      complex*8, dimension(0:n,0:n)    :: tband
      real*8, dimension(0:n,0:n)       :: p1m, p1p, p2m, p2p
      real*8, dimension(1:l,0:n,0:nd)  :: nlx, nlv, nlw                    !non-linear term convrsion matrices
      real*8, dimension(1:l,0:n+5,0:n) :: cux, cuv, cuw, csx, csv, csw     !velocity and vorticity conversion matrices
      
      end module
