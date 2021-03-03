      !velocity and vorticity conversion matrices
      do k=-5,5
      do i=0,n+5
      do j=1,l
      if((i-k).ge.(0).and.(i-k).le.(n)) then      
      cux(j,i,i-k) = - .25d0*b13(i-k,k)
      cuv(j,i,i-k) = sqrt(dble(j)/(2.d0*dble(j)+1.d0))*(.0625d0*     &
     &               dble(1-j)*b23(i-k,k) + .25d0*b33(i-k,k))
      cuw(j,i,i-k) = sqrt(dble(j+1)/(2.d0*dble(j) + 1.d0))*          &
     &               (.0625d0*dble(2+j)*b23(i-k,k) + .25d0*b33(i-k,k))
      csx(j,i,i-k) = - .0625d0*dble((2+j)*(1-j))*b23(i-k,k) -        &
     &               b33(i-k,k) - b34(i-k,k)
      csv(j,i,i-k) = sqrt(dble(j)/(2.d0*dble(j) + 1.d0))*(-          &
     &               .5d0*b24(i-k,k) + .25d0*b26(i-k,k) +            &
     &               .25d0*dble(j)*b13(i-k,k))
      csw(j,i,i-k) = sqrt(dble(j+1)/(2.d0*dble(j) + 1.d0))*(-        &
     &               .5d0*b24(i-k,k) + .25d0*b26(i-k,k) -            & 
     &               .25d0*dble(1+j)*b13(i-k,k))
      end if
      end do
      end do
      end do
      
      !matrices for non-linear terms
      do j=0,nd
      do i=0,n
      do x=1,l
      do k=j-4,j+4	
      if((k-i).ge.(-4).and.(k-i).le.(4)) then
      nlx(x,i,j) = nlx(x,i,j) - .125d0*b1(j,k-j)*b13(i,k-i)*        &
     &             .5d0*pi*c(k) 
      nlv(x,i,j) = nlv(x,i,j) - .5d0*sqrt(dble(x)/(2.d0*dble(x) +   &
     &             1.d0))*(.0625d0*(1.d0 - dble(x))*b13(j,k-j)*     &
     &             b13(i,k-i) - .0625d0*b24(j,k-j)*b13(i,k-i) +     &
     &             .25d0*b1(j,k-j)*b33(i,k-i))*.5d0*pi*c(k)
      nlw(x,i,j) = nlw(x,i,j) - .5d0*sqrt((1.d0 + dble(x))/         &
     &             (2.d0*dble(x) + 1.d0))*(.0625d0*(2.d0+dble(x))   &
     &             *b13(j,k-j)*b13(i,k-i) -.0625d0*b24(j,k-j)       &
     &             *b13(i,k-i) +.25d0*b1(j,k-j)*b33(i,k-i))*        &
     &             .5d0*pi*c(k)
      end if
      end do
      end do
      end do
      end do
      
      !inertia and diffusion matrices
      do j=1,l
      p1m(:,:) = 0.d0
      p1m =  .125d0*.125d0*.125d0*mult1(b27,b27) - (.5d0*tau*reinv  &
     &       )*(.0625d0*mult1(b25,b13) - .125d0*mult1(b26,b24)      &     
     &       - .03125d0*mult1(b16,b25) + .125d0*mult1(b24,b13) -    &
     &       .0625d0*mult1(b26,b13) - .03125d0*dble(j*(j+1))*       &
     &       mult1(b13,b13))
      t1m(j,:,:) = p1m(:,:)
      
      p1p(:,:) = 0.d0
      p1p = .125d0*dble(j*j+j+4)*mult1(b35,b35) + .25d0*(.03125d0*  & 
     &      .03125d0*mult1(b35,b36) + .03125d0*.125d0*              &
     &      mult1(b36,b37) + .03125d0*.125d0*mult1(b37,b35)) +      &
     &      .125d0*(- .125d0*.03125d0*mult1(b37,b36) + .125d0*      &
     &      .125d0*mult1(b37,b37)) - (.5d0*tau*reinv)*(.25d0*(      &
     &      .125d0*.03125d0*mult2(b20,b35) + 0.75d0*.0625d0*        &
     &      mult1(b32,b28) +9.d0*.125d0*.125d0*mult2(b21,b35) +     &
     &      3.d0*.03125d0*mult1(b36,b27))+ .125d0*(- .75d0*         &
     &      .125d0*mult3(b20,b24) + .25d0*.125d0*mult3(b20,b26)     &
     &      - 9.d0*.125d0*.25d0*mult3(b21,b28) + 9.d0*.125d0*       &
     &      mult3(b21,b33)- 3.d0*.125d0*mult2(b24,b28) + 3.d0*      &
     &      mult2(b24,b33) + 9.d0*.0625d0*mult4(b32,b24) - 3.d0*    &
     &      .0625d0*mult4(b32,b26)) + dble(j*j+j+12)*.125d0*        &
     &      .0625d0*mult1(b34,b23) + dble(j*j+j+6)*.125d0*.25d0*    &
     &      .0625d0*mult1(b33,b35) + dble(j*(1-j)*(j+1)*(j+2))*     &
     &      .0625d0*.0625d0*.125d0*mult1(b23,b23) + 6.d0*.125d0*    &
     &      (- .0625d0*mult1(b34,b28) + .25d0*mult1(b34,b33)) +     &
     &      .125d0*dble((2-j)*(3+j))*(- .25d0*.0625d0*              &  
     &      mult1(b33,b28) + .25d0*.25d0*mult1(b33,b33)) )
      t1p(j,:,:) = p1p(:,:)
      
      p2m(:,:) = 0.d0
      p2m =  .125d0*.125d0*.125d0*mult1(b27,b27) + (.5d0*tau*reinv  &
     &       )*(.0625d0*mult1(b25,b13) - .125d0*mult1(b26,b24)      &     
     &       - .03125d0*mult1(b16,b25) + .125d0*mult1(b24,b13) -    &
     &       .0625d0*mult1(b26,b13) - .03125d0*dble(j*(j+1))*       &
     &       mult1(b13,b13))
      t2m(j,:,:) = p2m(:,:)
      
      p2p(:,:) = 0.d0
      p2p = .125d0*dble(j*j+j+4)*mult1(b35,b35) + .25d0*(.03125d0*  & 
     &      .03125d0*mult1(b35,b36) + .03125d0*.125d0*              &
     &      mult1(b36,b37) + .03125d0*.125d0*mult1(b37,b35)) +      &
     &      .125d0*(- .125d0*.03125d0*mult1(b37,b36) + .125d0*      &
     &      .125d0*mult1(b37,b37)) + (.5d0*tau*reinv)*(.25d0*(      &
     &      .125d0*.03125d0*mult2(b20,b35) + 0.75d0*.0625d0*        &
     &      mult1(b32,b28) +9.d0*.125d0*.125d0*mult2(b21,b35) +     &
     &      3.d0*.03125d0*mult1(b36,b27))+ .125d0*(- .75d0*         &
     &      .125d0*mult3(b20,b24) + .25d0*.125d0*mult3(b20,b26)     &
     &      - 9.d0*.125d0*.25d0*mult3(b21,b28) + 9.d0*.125d0*       &
     &      mult3(b21,b33)- 3.d0*.125d0*mult2(b24,b28) + 3.d0*      &
     &      mult2(b24,b33) + 9.d0*.0625d0*mult4(b32,b24) - 3.d0*    &
     &      .0625d0*mult4(b32,b26)) + dble(j*j+j+12)*.125d0*        &
     &      .0625d0*mult1(b34,b23) + dble(j*j+j+6)*.125d0*.25d0*    &
     &      .0625d0*mult1(b33,b35) + dble(j*(1-j)*(j+1)*(j+2))*     &
     &      .0625d0*.0625d0*.125d0*mult1(b23,b23) + 6.d0*.125d0*    &
     &      (- .0625d0*mult1(b34,b28) + .25d0*mult1(b34,b33)) +     &
     &      .125d0*dble((2-j)*(3+j))*(- .25d0*.0625d0*              &  
     &      mult1(b33,b28) + .25d0*.25d0*mult1(b33,b33)) )
      t2p(j,:,:) = p2p(:,:)
      
      end do
