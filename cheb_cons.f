      c(:) = 0.d0
      d(:) = 0.d0
      c(0) = 2.d0
      d(0) = 1.d0
      do i=-5,n+5
      if (i.gt.0) then
      c(i) = 1.d0
      d(i) = 1.d0
      end if
      e(i) = c(i) + c(i-1) - 4.d0*d(i)
      f(i) = dble(i*(i+1))
      g(i) = 1.d0 - dble(i)/4.d0
      p(i) = c(i)*(3.d0 + .25d0*f(i-1)) + dble(2*i)
      q(i) = d(i-2)*(3.d0 - dble(2*i) + .25d0*f(i))
      r(i) = c(i)*(3.d0 - .25d0*f(i-1)) + c(i-1)*(3.d0 - dble(2*i)) &
     &       + dble(2*i) - 4.d0 - .25d0*d(i-2)*f(i)
      end do
      
      do i=0,nd
      rad(i) = .5d0*(radi+rado) + .5d0*(rado-radi)*                 &
     &         cos(dble(i)*pi/dble(nd))
      end do
      
      
