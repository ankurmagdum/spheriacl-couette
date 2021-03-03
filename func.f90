      module func
      
      implicit none
      contains
      
      function mult1(a,b)
      
      use para
      
      integer :: i, j, k
      real*8, dimension(0:n,-5:5) :: a, b
      real*8, dimension(0:n,0:n) :: mult1
      
      include "cheb_cons.f"
      
      mult1(:,:) = 0.d0
                  
      do i=0,n
      do j=0,n
      do k=j-5,j+5
      if((k-i).ge.(-5).and.(k-i).le.(5)) then
      mult1(i,j) = mult1(i,j) + b(j,k-j)*a(i,k-i)*.5d0*pi*c(k)
      end if
      end do
      end do
      end do
      
      return
      end function
      
      function mult2(a,b)
      
      use para
      
      integer :: i, j, k
      real*8, dimension(0:n,-5:5) :: a,b
      real*8, dimension(0:n,-7:7) :: auxb
      real*8, dimension(0:n,0:n) :: mult2
      !real*16, dimension(-5:n+5,-2:2) :: auxb25
      
      include "cheb_cons.f"
      include "cheb_arrays.f"
      
      mult2(:,:) = 0.d0
      
      do i=0,n
      do j=i-7,i+7
      do k=-2,2
      if ((j-k-i).ge.(-5).and.(j-k-i).le.(5)) then
      auxb(i,j-i) = b(i,j-k-i)*auxb25(j-i,k)
      end if
      end do
      end do
      end do
      
      mult2 = mult1(a,auxb)
      
      end function
      
      function mult3(a,b)
      
      use para
      
      integer :: i, j, k
      real*8, dimension(0:n,-5:5) :: a,b
      real*8, dimension(0:n,-7:7) :: auxa, auxb
      real*8, dimension(0:n,0:n) :: mult3
      !real*16, dimension(-5:n+5,-2:2) :: auxb14, auxb25
      
      include "cheb_cons.f"
      include "cheb_arrays.f"
      
      mult3(:,:) = 0.d0
      
      do i=0,n
      do j=i-7,i+7
      do k=-2,2
      if ((j-k-i).ge.(-5).and.(j-k-i).le.(5)) then
      auxa(i,j-i) = b(i,j-k-i)*auxb14(j-i,k)
      end if
      end do
      end do
      end do
            
      do i=0,n
      do j=i-7,i+7
      do k=-2,2
      if ((j-k-i).ge.(-5).and.(j-k-i).le.(5)) then
      auxb(i,j-i) = b(i,j-k-i)*auxb25(j-i,k)
      end if
      end do
      end do
      end do
      
      do i=0,n
      do j=0,n
      do k=j-7,j+7
      if((k-i).ge.(-7).and.(k-i).le.(7)) then
      mult3(i,j) = mult3(i,j) + b(j,k-j)*a(i,k-i)*.5d0*pi*c(k)
      end if
      end do
      end do
      end do
      
      end function
      
      function mult4(a,b)
      
      use para
      
      integer :: i, j, k
      real*8, dimension(0:n,-5:5) :: a,b
      real*8, dimension(0:n,-7:7) :: auxa, auxb
      real*8, dimension(0:n,0:n) :: mult4
      !real*16, dimension(-5:n+5,-2:2) :: auxb24, auxb25
      
      include "cheb_cons.f"
      include "cheb_arrays.f"
      
      mult4(:,:) = 0.d0
      
      do i=0,n
      do j=i-7,i+7
      do k=-2,2
      if ((j-k-i).ge.(-5).and.(j-k-i).le.(5)) then
      auxa(i,j-i) = b(i,j-k-i)*auxb24(j-i,k)
      end if
      end do
      end do
      end do
            
      do i=0,n
      do j=i-7,i+7
      do k=-2,2
      if ((j-k-i).ge.(-5).and.(j-k-i).le.(5)) then
      auxb(i,j-i) = b(i,j-k-i)*auxb25(j-i,k)
      end if
      end do
      end do
      end do
      
      do i=0,n
      do j=0,n
      do k=j-7,j+7
      if((k-i).ge.(-7).and.(k-i).le.(7)) then
      mult4(i,j) = mult4(i,j) + b(j,k-j)*a(i,k-i)*.5d0*pi*c(k)
      end if
      end do
      end do
      end do
      
      end function
      
      end module
