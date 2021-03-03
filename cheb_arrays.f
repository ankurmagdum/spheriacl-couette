      b1(:,:) = 0.d0
      b13(:,:) = 0.d0
      b16(:,:) = 0.d0
      b20(:,:) = 0.d0
      b21(:,:) = 0.d0
      b22(:,:) = 0.d0
      b23(:,:) = 0.d0
      b24(:,:) = 0.d0
      b25(:,:) = 0.d0
      b26(:,:) = 0.d0
      b27(:,:) = 0.d0
      b28(:,:) = 0.d0
      b32(:,:) = 0.d0
      b33(:,:) = 0.d0
      b34(:,:) = 0.d0
      b35(:,:) = 0.d0
      b36(:,:) = 0.d0
      b37(:,:) = 0.d0
      auxb14(:,:) = 0.d0
      auxb24(:,:) = 0.d0
      auxb25(:,:) = 0.d0
      
      do i=0,n
      b1(i,0) = 1.d0
      
      b13(i,-2) = d(i-2)
      b13(i,0) = e(i)
      b13(i,2) = c(i)
      
      b16(i,-2) = -d(i-2)*f(i)
      b16(i,0) = c(i)*f(i-1) + d(i-2)*f(i)
      b16(i,2) = -c(i)*f(i-1)
      
      b20(i,-3) = d(i-3)*dble(i)*f(i+1)
      b20(i,-1) = -d(i-3)*(dble(i+2)*f(i-2) + dble(2*i-2)*f(i+1))
      b20(i,1) = d(i-3)*(dble(2*i+2)*f(i-2) + dble(i-2)*f(i+1))
      b20(i,3) = -d(i-3)*dble(i)*f(i-2)
      
      b21(i,-3) = dble(i)*d(i-3)
      b21(i,-1) = dble(i)*(e(i-1) + 3.d0*d(i-1) - 4.d0/3.d0)
      b21(i,1) = dble(i)*(c(i-1) - c(i) + 1.d0/3.d0)
      b21(i,3) = -dble(i)
      
      b22(i,-3) = d(i-3)*g(i)
      b22(i,-1) = c(i-2) + d(i-1)*e(i) + .25d0*dble(i)*(1.d0-e(i-1))
      b22(i,1) = c(i) + c(i)*e(i) - .25d0*dble(i)*(2.d0 + c(i-1))
      b22(i,3) = c(i) + .25d0*dble(i)
      
      b23(i,-4) = d(i-4)
      b23(i,-2) = e(i-2) - 2.d0*d(i-2)
      b23(i,0) = c(i) + c(i-2) + e(i)*e(i)
      b23(i,2) = c(i)*(e(i)-2)
      b23(i,4) = c(i)
      
      b24(i,-2) = d(i-2)
      b24(i,-1) = 2.d0*km*d(i-1)
      b24(i,0) = c(i) + c(i-1)
      b24(i,1) = 2.d0*km*c(i)
      b24(i,2) = c(i)
      
      b25(i,-2) = d(i-2)
      b25(i,-1) = 4.d0*km*d(i-1)
      b25(i,0) = c(i) + c(i-1) + 4.d0*km*km
      b25(i,1) = 4.d0*km*c(i)
      b25(i,2) = c(i)
      
      b26(i,-2) = dble(i)*d(i-2)
      b26(i,-1) = 2.d0*dble(i)*km
      b26(i,0) = dble(i)*(c(i-1) - 1.d0)
      b26(i,1) = -2.d0*km*dble(i)
      b26(i,2) = -dble(i)
      
      b27(i,-3) = d(i-3)
      b27(i,-2) = 2.d0*km*d(i-2)
      b27(i,-1) = d(i-1)*e(i) + c(i-2)
      b27(i,0) = 2.d0*km*e(i)
      b27(i,1) = c(i)*(1.d0+ e(i))
      b27(i,2) = 2.d0*km*c(i)
      b27(i,3) = c(i)
      
      b28(i,-4) = d(i-4)
      b28(i,-3) = 2.d0*km*d(i-3)
      b28(i,-2) = d(i-2)*e(i) + c(i-2) + c(i-3)
      b28(i,-1) = 2.d0*km*(d(i-1)*e(i) + c(i-2))
      b28(i,0) = c(i)*(1.d0+e(i)) + c(i-1)*e(i) + c(i-2)
      b28(i,1) = 2.d0*km*c(i)*(1.d0+e(i))
      b28(i,2) = c(i)*(2.d0+e(i))
      b28(i,3) = 2.d0*km*c(i)
      b28(i,4) = c(i)
      
      b32(i,-4) = d(i-4)*f(i)
      b32(i,-3) = 4.d0*km*d(i-3)*f(i)
      b32(i,-2) = f(i)*(c(i-2) + c(i-3) + d(i-2)*(4.d0*km*km-1.d0)) &
     &            - f(i-1)*d(i-2)
      b32(i,-1) = 4.d0*km*(f(i)*(c(i-2)-d(i-2))-f(i-1)*d(i-1))
      b32(i,0) = f(i)*(c(i-2)-2.d0*d(i-2)*(1.d0 + 2.d0*km*km)) + &
     &           f(i-1)*c(i)*(1.d0 - 4.d0*km*km-c(i)-c(i-1))
      b32(i,1) = 4.d0*km*(f(i-1)*c(i)*(1.d0-c(i))-f(i)*d(i-2))
      b32(i,2) = f(i-1)*c(i)*(2.d0+4.d0*km*km-c(i)) - f(i)*d(i-2)
      b32(i,3) = 4.d0*km*f(i-1)*c(i)
      b32(i,4) = f(i-1)*c(i)
      
      b33(i,-4) = d(i-4)*g(i)
      b33(i,-3) = 2.d0*km*d(i-3)*g(i)
      b33(i,-2) = g(i)*(c(i-2)+c(i-3)) + d(i-2)*(i-2)
      b33(i,-1) = 2.d0*km*(c(i-2) + d(i-1)*e(i) + (1.d0-e(i-1))*   &
     &            dble(i)*.25) 
      b33(i,0) = c(i)*(1.d0 + e(i)) - (2.d0+c(i-1))*i*.25d0 +  &
     &           c(i-1)*(c(i-2) + e(i) + (1.d0-e(i-1))*.25d0*dble(i))
      b33(i,1) = 2.d0*km*(c(i)*(1.d0+e(i)) - (2.d0+c(i-1))*    &
     &           .25d0*dble(i)) 
      b33(i,2) = c(i)*(e(i)+2.d0) - (1.d0+c(i-1))*.25d0*dble(i)
      b33(i,3) = 2.d0*km*(c(i) + .25d0*dble(i))
      b33(i,4) = c(i) + .25d0*dble(i)
      
      b34(i,-4) = d(i-4)*q(i)*.25d0
      b34(i,-3) = km*d(i-3)*q(i)
      b34(i,-2) = .25d0*(r(i)*d(i-2) + q(i)*(c(i-2)+c(i-3) + &
     &            .25d0*km*km)) 
      b34(i,-1) = km*(r(i)*d(i-1) + q(i)*c(i-2))
      b34(i,0) = .25d0*(p(i) + q(i)*c(i-2) + r(i)*(c(i) + c(i-1)+ &
     &           .25d0*km*km))
      b34(i,1) = km*(p(i) + r(i)*c(i))
      b34(i,2) = p(i)*(km*km + .5d0) + r(i)*c(i)*.25d0
      b34(i,3) = km*p(i)
      b34(i,4) = .25d0*p(i)
      
      b35(i,-5) = d(i-5)
      b35(i,-4) = 2.d0*km*d(i-4)
      b35(i,-3) = d(i-3)*(c(i-3)-5.d0) + c(i-4)
      b35(i,-2) = 2.d0*km*(e(i-2) - 2.d0*d(i-2))
      b35(i,-1) = d(i-1)*(c(i) + c(i-2) + e(i)*e(i)) + &
     &            c(i-2)*(e(i-2) -2.d0)
      b35(i,0) = 2.d0*km*(c(i) + c(i-2) + e(i)*e(i))
      b35(i,1) = c(i)*(2.d0*c(i) + c(i-1) + c(i-2) - 6.d0 + e(i)*e(i))
      b35(i,2) = 2.d0*km*c(i)*(e(i)-2.d0)
      b35(i,3) = c(i)*(e(i)-1.d0)
      b35(i,4) = 2.d0*km*c(i)
      b35(i,5) = c(i)
      
      b36(i,-5) = d(i-5)
      b36(i,-4) = 4.d0*km*d(i-4)
      b36(i,-3) = d(i-3)*(1.d0+4.d0*km*km+e(i)) + c(i-3) + c(i-4)
      b36(i,-2) = 4.d0*km*(d(i-2)*e(i) + c(i-2) + c(i-3))
      b36(i,-1) = d(i-1)*(1.d0+(1.d0+4.d0*km*km)*e(i)) +       &
     &            c(i-1)*e(i) + c(i-2)*(1.d0 + 4.d0*km*km + &
     &            e(i) + c(i-2) + c(i-3))
      b36(i,0) = 4.d0*km*(c(i) + c(i-2) + (c(i)+c(i-1))*e(i))
      b36(i,1) = c(i)*(2.d0 + 4.d0*km*km + c(i) + c(i-2) + (1.d0 &
     &            + 4.d0*km*km + e(i) + c(i-2) + c(i-3))*e(i))
      b36(i,2) = 4.d0*km*c(i)*(e(i) + 2.d0)
      b36(i,3) = c(i)*(e(i) + 3.d0 + 4.d0*km*km)
      b36(i,4) = 4.d0*km*c(i)
      b36(i,5) = c(i)
      
      b37(i,-5) = d(i-5)*g(i)
      b37(i,-4) = km*d(i-4)*(4-i)
      b37(i,-3) = d(i-3)*dble(i-2) + (d(i-3)*(1.d0 + 4.d0*km*km) &
     &             + c(i-3) + c(i-4))*g(i)
      b37(i,-2) = 4.d0*km*(g(i)*(c(i-2) + c(i-3)) + d(i-2)*dble(i-2))
      b37(i,-1) = g(i)*(c(i-3) + c(i-2)*c(i-2) + c(i-1) +     &
     &            c(i-2)*dble(i-2) - d(i-1)*(2.d0 + .5d0*i) +  &
     &            (c(i-1) + 4.d0*km*km)*(c(i-2) + c(i-1) -  &
     &            3.d0*d(i-1) + (5.d0 - c(i-1) -c(i-2))*dble(i)*.25d0))
      b37(i,0) = 4.d0*km*(c(i-1)*(c(i-1)*g(i) + dble(i) - 2.d0)     &
     &           + c(i)*(c(i) - 3.d0) + c(i-2)*g(i) - .5d0*dble(i))
      b37(i,1) = c(i)*c(i)*(c(i) + 2.d0*c(i-1) - 2.d0 + 4.d0*km*km) + &
     &           c(i-1)*(c(i-1) + c(i-2) - e(i-1)*dble(i)*.25d0 - &
     &           3.d0 + 4.d0*km*km) - c(i)*(2.d0 + .5d0*dble(i) +   &
     &           12.d0*km*km) - c(i-1)*(1.d0 + 4.d0*km*km)*dble(i)* &
     &           .25d0 - (1.d0 + 8.d0*km*km)*dble(i)*.25d0
      b37(i,2) = 4.d0*km*(c(i)*(e(i) + 2.d0) - (1.d0 + c(i-1))*  &
     &           dble(i)*.25d0)
      b37(i,3) = (1.d0 + 4.d0*km*km)*(c(i) + dble(i)*.25d0) +  &
     &           c(i)*(e(i) + 2.d0) - (1.d0 + c(i-1))*.25d0*dble(i)
      b37(i,4) = km*(4.d0*c(i) + dble(i))
      b37(i,5) = c(i) + .25d0*dble(i)    
      end do 
      
      do i=-5,n+5
      auxb14(i,-1) = d(i-1)
      auxb14(i,0) = 2.d0*km
      auxb14(i,1) = c(i)
      
      auxb24(i,-2) = d(i-2)
      auxb24(i,-1) = 2.d0*km*d(i-1)
      auxb24(i,0) = c(i) + c(i-1)
      auxb24(i,1) = 2.d0*km*c(i)
      auxb24(i,2) = c(i)
      
      auxb25(i,-2) = d(i-2)
      auxb25(i,-1) = 4.d0*km*d(i-1)
      auxb25(i,0) = c(i) + c(i-1) + 4.d0*km*km
      auxb25(i,1) = 4.d0*km*c(i)
      auxb25(i,2) = c(i) 
      end do
