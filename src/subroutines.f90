      module subroutines
      
      use parameters
            
      implicit none
      
      include "../fftw3.f"
      include "../shtns.f"
      
      contains
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine for_vsh_trans (a1,a2,a3)                 !working
                  
      integer :: i0, j, k, x, lm
      complex*8, dimension(0:ld,-md-1:md) :: a1, a2, a3
      real*8, dimension(-md-1:md,0:ld)    :: h1, h2, h3
      complex*16, allocatable             :: qlm(:), slm(:), tlm(:)
      
      allocate(qlm(1:nlmd),slm(1:nlmd),tlm(1:nlmd))
      
      do k=0,md
      call shtns_lmidx(lm,k,k)
      do j=lm,lm+ld-k
      x = j-lm+k
      qlm(j) = -sqrt((dble(x) + 1.d0)/(2.d0*dble(x) + 1.d0))*a2(x,k) &
     &         + sqrt(dble(x)/(2.d0*dble(x) + 1.d0))*a3(x,k)
      slm(j) = sqrt(1.d0/((dble(x) + 1.d0)*(2.d0*dble(x) + 1.d0)))   &
     &         *a2(x,k) + sqrt(1.d0/(dble(x)*(2.d0*dble(x) + 1.d0))) &
     &         *a3(x,k) 
      tlm(j) = iota*a1(x,k)/sqrt(dble(x)*(dble(x) + 1.d0))
      end do
      end do
      
      slm(1) = 0.d0
      tlm(1) = 0.d0
      
      call shtns_qst_to_spat (qlm,slm,tlm,h1,h2,h3)
      
      do j=-md-1,md
      do i0=0,ld
      a1(i0,j) = h1(j,i0)
      a2(i0,j) = h2(j,i0)
      a3(i0,j) = h3(j,i0)
      end do
      end do
      
      end subroutine
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine back_vsh_trans (a1,a2,a3)                !working
      
      integer :: i0, j, k, x, lm 
      complex*8, dimension(0:ld,-md-1:md) :: a1, a2, a3
      real*8, dimension(-md-1:md,0:ld)    :: h1, h2, h3
      complex*16, allocatable             :: qlm(:), slm(:), tlm(:)
      
      allocate(qlm(1:nlmd),slm(1:nlmd),tlm(1:nlmd))
      
      do j=-md-1,md
      do i0=0,ld
      h1(j,i0) = real(a1(i0,j))
      h2(j,i0) = real(a2(i0,j))
      h3(j,i0) = real(a3(i0,j))
      end do
      end do
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      
      call shtns_spat_to_qst (h1,h2,h3,qlm,slm,tlm)
      
      do k=0,md
      call shtns_lmidx(lm,k,k)
      do j=lm,lm+ld-k
      x = j-lm+k
      a1(x,k) = -iota*sqrt(dble(x)*(dble(x) + 1.d0))*tlm(j)
      a2(x,k) = sqrt(dble(x+1)/(2.d0*dble(x) + 1.d0))*(qlm(j) - &
     &          dble(x)*slm(j))
      a3(x,k) = sqrt(dble(x)/(2.d0*dble(x) + 1.d0))*(qlm(j) +   &
     &          dble(x+1)*slm(j))
      end do
      end do
      
      end subroutine
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine for_cheb_trans (a1)                         !working
      
      integer :: i0
      complex*8, dimension(0:nd) :: a1
      real*8, dimension(0:nd)    :: z1
      real*8, dimension(0:nd)    :: z2
      
      c1(:) = 2.d0
      c1(0) = 1.d0
      c1(nd) = 1.d0
            
      do i0=0,nd
      z1(i0) = real(a1(i0))/c1(i0)
      z2(i0) = aimag(a1(i0))/c1(i0)
      end do
                 
      call dfftw_plan_r2r_1d (plan, nd+1, z1, z1, &
     &     fftw_redft00, fftw_estimate)
      call dfftw_execute_r2r (plan, z1, z1)
      call dfftw_destroy_plan (plan)
      
      call dfftw_plan_r2r_1d (plan, nd+1, z2, z2, &
     &     fftw_redft00, fftw_estimate)
      call dfftw_execute_r2r (plan, z2, z2)
      call dfftw_destroy_plan (plan)
      
      
      do i0=0,nd
      a1(i0) = z1(i0) + iota*z2(i0)
      end do
                  
      end subroutine
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine back_cheb_trans (a1)                         !works beautifully
      
      integer :: i0
      complex*8, dimension(0:nd) :: a1
      real*8, dimension(0:nd)    :: z1
      real*8, dimension(0:nd)    :: z2
      
      c1(:) = 2.d0
      c1(0) = 1.d0
      c1(nd) = 1.d0
                 
      do i0=0,nd
      z1(i0) = real(a1(i0))
      z2(i0) = aimag(a1(i0))
      end do
      
      call dfftw_plan_r2r_1d (plan, nd+1, z1, z1, &
     &     fftw_redft00, fftw_estimate)
      call dfftw_execute_r2r (plan, z1, z1)
      call dfftw_destroy_plan (plan)
      
      call dfftw_plan_r2r_1d (plan, nd+1, z2, z2, &
     &     fftw_redft00, fftw_estimate)
      call dfftw_execute_r2r (plan, z2, z2)
      call dfftw_destroy_plan (plan)
      
      do i0=0,nd
      a1(i0) = (z1(i0) + iota*z2(i0))*c1(i0)/dble(2*nd)
      end do
            
      end subroutine
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine compute_nonlinear(a1,a2,a3,a4,a5,a6,i0)       
      
      integer :: i, j, i0
      complex*8, dimension(0:ld,-md-1:md) :: a1, a2, a3, a4, a5, a6
      complex*8, dimension(0:ld,-md-1:md) :: b1, b2, b3
      
      do j=-md-1,md
      do i=0,ld
      b1(i,j) = a2(i,j)*a6(i,j) - a3(i,j)*a5(i,j)
      b2(i,j) = a3(i,j)*a4(i,j) - a1(i,j)*a6(i,j)
      b3(i,j) = a1(i,j)*a5(i,j) - a4(i,j)*a2(i,j)
      end do
      end do
      
      a1(:,:) = rad(i0)*b1(:,:)
      a2(:,:) = rad(i0)*b2(:,:)
      a3(:,:) = rad(i0)*b3(:,:)
      
      end subroutine
      
      end module
