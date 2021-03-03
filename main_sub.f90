      program main
      
      use para
      use func
      use sub
      
      implicit none
      
      integer :: i, j, k, info, t, x, i0, i1
      real*8  :: tau=0.01d0
      integer, dimension(0:n) :: ipiv
      complex*8, dimension(1:l,0:nd,-1-m:m) :: f1, f2, f3
      complex*8, dimension(0:n,0:m)         :: am, ap
      complex*8, dimension(1:l,0:n,0:m)     :: fm, fp
      complex*8, dimension(0:ld,-md-1:md), target :: a1, a2, a3,  &
     &                                               a4, a5, a6
      complex*8, dimension(:), pointer :: a1p, a2p, a3p, a4p, a5p, a6p
      complex*8, dimension(0:ld,-md-1:md) :: v1, v2, v3
      complex*16, allocatable             :: qlm(:), slm(:), tlm(:)
      real*8, dimension(0:nd,0:ld,-md-1:md) :: g1, g2, g3
      real*8, dimension(0:nd,0:ld,-md-1:md) :: u1, u2, u3
      real*8, dimension(0:nd,0:ld,-md-1:md) :: res1, res2, res3
      real*8 :: ma1, ma2, ma3, loc1, loc2, loc3
      
      include "cheb_cons.f"
      include "cheb_arrays.f"
      include "precomp.f"
                          
      call shtns_calc_nlm (nlmd, ld, md, mres)
      layout = SHT_PHI_CONTIGUOUS
      call shtns_init_sh_gauss(layout, ld, md, mres, ld+1, 2*(md+1))
      allocate(qlm(1:nlmd),slm(1:nlmd),tlm(1:nlmd))
      
      open(20,file="test1.dat") 
      open(22,file="test2.dat")
      open(24,file="test3.dat")
      
      print *, 'nd,ld,md-->', nd, ld, md
                  
      vbc(:,:) = 0.d0
      vbc(0,0) = pi*.5d0*iota*sqrt(2.d0*pi/3.d0)*(1.d0-km)
      vbc(2,0) = -pi*.25d0*iota*sqrt(2.d0*pi/3.d0)*(1.d0-km)
       
      fm(:,:,:) = 0.d0
      fp(:,:,:) = 0.d0
      am(:,:) = 0.d0
      ap(:,:) = 0.d0
      tband(:,:) = 0.d0
                              
      !first time-step
      f1(:,:,:) = 0.d0      !radial in spatial, Xlm in spectral
      f2(:,:,:) = 0.d0      !polar in spatial, Vlm in spectral
      f3(:,:,:) = 0.d0      !azimuthal in spatial, Wlm in spectral
      
      !coefficients of velocity and vorticity
      
      f1(1,0,0) = - iota*sqrt(2.d0*pi/3.d0)       !velocity
      f1(1,1,0) = - iota*sqrt(2.d0*pi/3.d0)
      
      f2(1,0,-1) = sqrt(2.d0*pi)*(.5d0*km-1.d0)/3.d0  !vorticity
      f3(1,0,-1) = sqrt(pi)*2.d0*(.5d0*km+2.d0)/3.d0
      f2(1,1,-1) = - sqrt(2.d0*pi)/6.d0
      f3(1,1,-1) = - sqrt(pi)*5.d0/3.d0
           
      do k=0,min(m,j)
      do j=1,l
      call for_cheb_trans (f1(j,:,k))
      call for_cheb_trans (f2(j,:,k))
      call for_cheb_trans (f3(j,:,k))
      call for_cheb_trans (f1(j,:,-1-k))
      call for_cheb_trans (f2(j,:,-1-k))
      call for_cheb_trans (f3(j,:,-1-k))
      end do
      end do
      
      !****************************************************************      
      
      do i=0,nd
      
      i0 = i                        !computing the de-alised non-linear term from coefficients 
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      a4(:,:) = 0.d0
      a5(:,:) = 0.d0
      a6(:,:) = 0.d0
                                    !of velocity and vorticity
      do k=0,min(m,j)
      do j=1,l
      a1(j,k) = f1(j,i,k)
      a2(j,k) = f2(j,i,k)
      a3(j,k) = f3(j,i,k)
      a4(j,k) = f1(j,i,-1-k)
      a5(j,k) = f2(j,i,-1-k)
      a6(j,k) = f3(j,i,-1-k)
      end do
      end do
      
      call for_vsh_trans (a1,a2,a3)
      call for_vsh_trans (a4,a5,a6)
      
      u1(i,:,:) = real(a1(:,:))
      u2(i,:,:) = real(a2(:,:))
      u3(i,:,:) = real(a3(:,:))
      
      call compute_nonlinear(a1,a2,a3,a4,a5,a6,i0)      !non-linear product stored in a1,a2,a3
           
      call back_vsh_trans (a1,a2,a3)
               
      f1(:,i,:) = 0.d0
      f2(:,i,:) = 0.d0
      f3(:,i,:) = 0.d0
      
      do k=0,min(m,j)            
      do j=1,l
      f1(j,i,k) = a1(j,k)
      f2(j,i,k) = a2(j,k)
      f3(j,i,k) = a3(j,k)
      end do
      end do
      
      end do
     
      print *, "first time-step"
      
      !****************************************************************
      
      do j=1,l
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
           
      do k=0,min(m,j)             !copy the VSH-transformed coefficients of non-linear term 
      do i=0,nd                   !computed and saved in the previous time-step
      a1(i,k) = f1(j,i,k)
      a2(i,k) = f2(j,i,k)
      a3(i,k) = f3(j,i,k)
      end do
      end do 
      
      do k=0,min(m,j)               !backward chebyshev transform the coefficients copied above
      a1p => a1(0:nd,k)
      call back_cheb_trans (a1p)
      a2p => a2(0:nd,k)
      call back_cheb_trans (a2p)
      a3p => a3(0:nd,k)
      call back_cheb_trans (a3p)
      end do
      
      if (j.eq.(1)) am(:,:) = reinv*vbc(:,:) 
      
      do k=0,min(m,j)               !compute right hand side of the matrix equation
      do i=0,n
      do i1=0,nd
      ap(i,k) = ap(i,k) + iota*nlv(j,i,i1)*a2(i1,k) &
     &                  + iota*nlw(j,i,i1)*a3(i1,k)
      am(i,k) = am(i,k) + nlx(j,i,i1)*a1(i1,k)
      end do
      fp(j,i,k) = tau*ap(i,k)
      fm(j,i,k) = tau*am(i,k)
      end do
      end do
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      a4(:,:) = 0.d0
      a5(:,:) = 0.d0
      a6(:,:) = 0.d0
           
      !----------------------------------
      do k=0,n
      do i=max(0,k-6),min(n,k+6)
      tband(12+i-k,k) = t1m(j,i,k)
      end do
      end do
      
      call cgbsv(n+1,6,6,m+1,tband,n+1,ipiv,fm(j,:,:),n+1,info)!    !solve matrix equation to obtain coefficients of velocity
      
      tband(:,:) = 0.d0
      do k=0,n
      do i=max(0,k-10),min(n,k+10)
      tband(20+i-k,k) = t1p(j,i,k)
      end do
      end do
      
      call cgbsv(n+1,10,10,m+1,tband,n+1,ipiv,fp(j,:,:),n+1,info)! 
      !----------------------------------
      
      do k=0,min(m,j)                    !compute part of right hand side of 
      do i=0,n+4                         !the matrix equation for next time-step
      if (i.le.n) then
      am(i,k) = -(.5d0*tau)*am(i,k)
      ap(i,k) = -(.5d0*tau)*ap(i,k)
      end if
      do i1=0,n
      if (i.le.n) then
      am(i,k) = am(i,k) + t2m(j,i,i1)*fm(j,i1,k)
      ap(i,k) = ap(i,k) + t2p(j,i,i1)*fp(j,i1,k)
      end if
      a1(i,k) = a1(i,k) + cux(j,i,i1)*fm(j,i1,k)          !compute coefficients of velocity and 
      a2(i,k) = a2(i,k) + iota*cuv(j,i,i1)*fp(j,i1,k)     !vorticity in standard representaion
      a3(i,k) = a3(i,k) + iota*cuw(j,i,i1)*fp(j,i1,k)
      a4(i,k) = a4(i,k) + csx(j,i,i1)*fp(j,i1,k)          !a1,a2,a3 => velocity
      a5(i,k) = a5(i,k) + iota*csv(j,i,i1)*fm(j,i1,k)     !a4,a5,a6 => vorticity
      a6(i,k) = a6(i,k) + iota*csw(j,i,i1)*fm(j,i1,k)
      end do
      end do
      end do
      
      if(j.eq.(1)) then
      a1(0,0) = a1(0,0) - iota*sqrt(2.d0*pi/3.d0)
      a1(1,0) = a1(0,0) - iota*sqrt(2.d0*pi/3.d0)
      
      a5(0,0) = a5(0,0) + sqrt(2.d0*pi)*(.5d0*km-1.d0)/3.d0  !vorticity
      a6(0,0) = a6(0,0) + sqrt(pi)*2.d0*(.5d0*km+2.d0)/3.d0
      a5(1,0) = a5(1,0) - sqrt(2.d0*pi)/6.d0
      a6(1,0) = a6(1,0) - sqrt(pi)*5.d0/3.d0
      end if
      
      do k=0,min(m,j)                   !forward chebyshev transform the coefficients computed above
      a1p => a1(0:nd,k)
      call for_cheb_trans (a1p)
      a2p => a2(0:nd,k)
      call for_cheb_trans (a2p)
      a3p => a3(0:nd,k)
      call for_cheb_trans (a3p)
      a4p => a4(0:nd,k)
      call for_cheb_trans (a4p)
      a5p => a5(0:nd,k)
      call for_cheb_trans (a5p)
      a6p => a6(0:nd,k)
      call for_cheb_trans (a6p)
      
      f1(j,:,k) = a1p(:)                !transfer the VSH-transformed coefficients from
      f1(j,:,-1-k) = a4p(:)             !a1,a2,a3,a4,a5,a6 to f1,f2,f3
      f2(j,:,k) = a2p(:)
      f2(j,:,-1-k) = a5p(:)
      f3(j,:,k) = a3p(:)
      f3(j,:,-1-k) = a6p(:)
      end do
      
      fm(j,:,:) = am(:,:)               !save the part of right hand side for next time-step
      fp(j,:,:) = ap(:,:)
      
      am(:,:) = 0.d0
      ap(:,:) = 0.d0
      
      end do
      
      !****************************************************************
      
      do i=0,nd                        !computing the de-alised non-linear term from coefficients 
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      a4(:,:) = 0.d0
      a5(:,:) = 0.d0
      a6(:,:) = 0.d0
      
      i0 = i                           !of velocity and vorticity
      
      do k=0,min(m,j)
      do j=1,l
      a1(j,k) = f1(j,i,k)
      a2(j,k) = f2(j,i,k)
      a3(j,k) = f3(j,i,k)
      a4(j,k) = f1(j,i,-1-k)
      a5(j,k) = f2(j,i,-1-k)
      a6(j,k) = f3(j,i,-1-k)
      end do
      end do
      
      f1(:,i,:) = 0.d0
      f2(:,i,:) = 0.d0
      f3(:,i,:) = 0.d0
      
      call for_vsh_trans (a1,a2,a3)
      call for_vsh_trans (a4,a5,a6)
      
      res1(i,:,:) = sqrt(real(a1(:,:))*real(a1(:,:))  &
     &              - u1(i,:,:)*u1(i,:,:))
      res2(i,:,:) = sqrt(real(a2(:,:))*real(a2(:,:))  &
     &              - u2(i,:,:)*u2(i,:,:))
      res3(i,:,:) = sqrt(real(a3(:,:))*real(a3(:,:))  &
     &              - u3(i,:,:)*u3(i,:,:))
      
      u1(i,:,:) = real(a1(:,:))
      u2(i,:,:) = real(a2(:,:))
      u3(i,:,:) = real(a3(:,:))
      
      call compute_nonlinear(a1,a2,a3,a4,a5,a6,i0)      !non-linear product stored in a1,a2,a3
            
      call back_vsh_trans (a1,a2,a3)
                  
      do j=1,l
      do k=0,m            
      f1(j,i,k) = a1(j,k)
      f2(j,i,k) = a2(j,k)
      f3(j,i,k) = a3(j,k)
      end do
      end do
      
      end do
      
      print *, "end of first time-step"      
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      tau=0.001d0
      do t=1,100000
      !if(t.le.100) 
      !if(t.gt.100.and.t.le.1000) tau=0.0001d0
      !if(t.gt.1000) tau=0.00001d0
      
      !******************time-marching*********************************
      
      do j=1,l
      
      do k=0,min(m,j)
      do i=0,nd                   !computed and saved in the previous time-step
      a1(i,k) = f1(j,i,k)
      a2(i,k) = f2(j,i,k)
      a3(i,k) = f3(j,i,k)
      end do
      end do      
      
      do k=0,min(m,j)               !backward chebyshev transform the coefficients copied above
      a1p => a1(0:nd,k)
      call back_cheb_trans (a1p)
      a2p => a2(0:nd,k)
      call back_cheb_trans (a2p)
      a3p => a3(0:nd,k)
      call back_cheb_trans (a3p)
      end do
      
      if (j.eq.(1)) am(:,:) = reinv*vbc(:,:) 
      
      do k=0,min(m,j)               !compute right hand side of the matrix equation
      do i=0,n
      do i1=0,nd
      ap(i,k) = ap(i,k) + iota*nlv(j,i,i1)*a2(i1,k) &
     &                  + iota*nlw(j,i,i1)*a3(i1,k)
      am(i,k) = am(i,k) + nlx(j,i,i1)*a1(i1,k)
      end do
      fp(j,i,k) = fp(j,i,k) + 1.5d0*tau*ap(i,k)
      fm(j,i,k) = fm(j,i,k) + 1.5d0*tau*am(i,k)
      end do
      end do
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      a4(:,:) = 0.d0
      a5(:,:) = 0.d0
      a6(:,:) = 0.d0
      
      !----------------------------------------------------------------
      
      tband(:,:) = 0.d0
      do k=0,n
      do i=max(0,k-6),min(n,k+6)
      tband(12+i-k,k) = t1m(j,i,k)
      end do
      end do
      
      call cgbsv(n+1,6,6,m+1,tband,n+1,ipiv,fm(j,:,:),n+1,info)    !solve matrix equation to obtain coefficients of velocity
      
      tband(:,:) = 0.d0
      do k=0,n
      do i=max(0,k-10),min(n,k+10)
      tband(20+i-k,k) = t1p(j,i,k)
      end do
      end do
      
      call cgbsv(n+1,10,10,m+1,tband,n+1,ipiv,fp(j,:,:),n+1,info)
       
      !----------------------------------------------------------------
      
      do k=0,min(m,j)                    !compute part of right hand side of 
      do i=0,n+4                         !the matrix equation for next time-step
      if (i.le.n) then
      am(i,k) = - .5d0*tau*am(i,k)
      ap(i,k) = - .5d0*tau*ap(i,k)
      end if
      do i1=0,n
      if (i.le.n) then
      am(i,k) = am(i,k) + t2m(j,i,i1)*fm(j,i1,k)
      ap(i,k) = ap(i,k) + t2p(j,i,i1)*fp(j,i1,k)
      end if
      a1(i,k) = a1(i,k) + cux(j,i,i1)*fm(j,i1,k)          !compute coefficients of velocity and 
      a2(i,k) = a2(i,k) + iota*cuv(j,i,i1)*fp(j,i1,k)     !vorticity in standard representaion
      a3(i,k) = a3(i,k) + iota*cuw(j,i,i1)*fp(j,i1,k)
      a4(i,k) = a4(i,k) + csx(j,i,i1)*fp(j,i1,k)          !a1,a2,a3 => velocity
      a5(i,k) = a5(i,k) + iota*csv(j,i,i1)*fm(j,i1,k)     !a4,a5,a6 => vorticity
      a6(i,k) = a6(i,k) + iota*csw(j,i,i1)*fm(j,i1,k)
      end do
      end do
      end do
      
      if(j.eq.(1)) then
      a1(0,0) = a1(0,0) - iota*sqrt(2.d0*pi/3.d0)
      a1(1,0) = a1(0,0) - iota*sqrt(2.d0*pi/3.d0)
      
      a5(0,0) = a5(0,0) + sqrt(2.d0*pi)*(.5d0*km-1.d0)/3.d0  !vorticity
      a6(0,0) = a6(0,0) + sqrt(pi)*2.d0*(.5d0*km+2.d0)/3.d0
      a5(1,0) = a5(1,0) - sqrt(2.d0*pi)/6.d0
      a6(1,0) = a6(1,0) - sqrt(pi)*5.d0/3.d0
      end if
      
      do k=0,min(m,j)                   !forward chebyshev transform the coefficients computed above
      a1p => a1(0:nd,k)
      call for_cheb_trans (a1p)
      a2p => a2(0:nd,k)
      call for_cheb_trans (a2p)
      a3p => a3(0:nd,k)
      call for_cheb_trans (a3p)
      a4p => a4(0:nd,k)
      call for_cheb_trans (a4p)
      a5p => a5(0:nd,k)
      call for_cheb_trans (a5p)
      a6p => a6(0:nd,k)
      call for_cheb_trans (a6p)
      
      f1(j,:,k) = a1p(:)                !transfer the VSH-transformed coefficients from
      f1(j,:,-1-k) = a4p(:)             !a1,a2,a3,a4,a5,a6 to f1,f2,f3
      f2(j,:,k) = a2p(:)
      f2(j,:,-1-k) = a5p(:)
      f3(j,:,k) = a3p(:)
      f3(j,:,-1-k) = a6p(:)
      end do
      
      fm(j,:,:) = am(:,:)               !save the part of right hand side for next time-step
      fp(j,:,:) = ap(:,:)
      
      am(:,:) = 0.d0
      ap(:,:) = 0.d0
      
      end do
      
      !****************computing non-linear term************************
      
      do i=0,nd                        !computing the de-alised non-linear term from coefficients 
      
      a1(:,:) = 0.d0
      a2(:,:) = 0.d0
      a3(:,:) = 0.d0
      a4(:,:) = 0.d0
      a5(:,:) = 0.d0
      a6(:,:) = 0.d0
      
      i0 = i 
      
      do k=0,min(m,j)                                !of velocity and vorticity
      do j=1,l
      a1(j,k) = f1(j,i,k)
      a2(j,k) = f2(j,i,k)
      a3(j,k) = f3(j,i,k)
      a4(j,k) = f1(j,i,-1-k)
      a5(j,k) = f2(j,i,-1-k)
      a6(j,k) = f3(j,i,-1-k)
      end do
      end do
      
      f1(:,i,:) = 0.d0
      f2(:,i,:) = 0.d0
      f3(:,i,:) = 0.d0
      
      call for_vsh_trans (a1,a2,a3)
      call for_vsh_trans (a4,a5,a6)
      
      res1(i,:,:) = sqrt(real(a1(:,:))*real(a1(:,:))  &
     &              - u1(i,:,:)*u1(i,:,:))
      res2(i,:,:) = sqrt(real(a2(:,:))*real(a2(:,:))  &
     &              - u2(i,:,:)*u2(i,:,:))
      res3(i,:,:) = sqrt(real(a3(:,:))*real(a3(:,:))  &
     &              - u3(i,:,:)*u3(i,:,:))
     
      u1(i,:,:) = real(a1(:,:))
      u2(i,:,:) = real(a2(:,:))
      u3(i,:,:) = real(a3(:,:))
      
      if(i.eq.(20)) then
      !print *, t, real(a1(i,m)), real(a2(i,m)), real(a3(i,m))
      write(20,*) t, real(a1(i,m))
      write(22,*) t, real(a2(i,m))
      write(24,*) t, real(a3(i,m))
      !g1(i,:,:) = real(a1(:,:))
      !g2(i,:,:) = real(a2(:,:))
      !g3(i,:,:) = real(a3(:,:))
      end if
      
      call compute_nonlinear(a1,a2,a3,a4,a5,a6,i0)      !non-linear product stored in a1,a2,a3
            
      call back_vsh_trans (a1,a2,a3)
                  
      do k=0,min(m,j)
      do j=1,l
      f1(j,i,k) = a1(j,k)
      f2(j,i,k) = a2(j,k)
      f3(j,i,k) = a3(j,k)
      end do
      end do
      
      end do
      
      ma1 = maxval(res1)
      ma2 = maxval(res2)
      ma3 = maxval(res3)
      
!      print *, t, max(ma1, ma2, ma3), maxloc(res1), &
!     & maxloc(res2), maxloc(res3)

      print*,'Max u1,u2,u3=',maxval(abs(real(u1))), &
       maxval(abs(real(u2))),maxval(abs(real(u3)))
      
      !****************************************************************
      
      end do
      
      open(100,file="test1.txt")
      do j=0,ld
      do i=0,nd
      write(100,"(200F10.7)",advance="no") g1(i,j,0)
      end do
      write(100,*) ""
      end do
      
   10 end program
