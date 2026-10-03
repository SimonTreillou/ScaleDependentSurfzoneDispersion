      subroutine step3d_t (tile)
      implicit none
      integer*4 tile
      return
      end
      function flux5_weno(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
!$acc routine(flux5_weno) seq
      implicit none
      REAL    :: flux5_weno
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua
      REAL    :: IS0, IS1, IS2
      REAl    :: d0, d1, d2
      REAl    :: a0, a1, a2
      REAL    :: w0, w1, w2
      REAL    :: p0, p1, p2
      REAL    :: Eps, cff1, cff2, T5
      Eps = 1.D-40
      d0=1.D0/10.D0
      d1=6.D0/10.D0
      d2=3.D0/10.D0
      cff1=13.D0/12.D0
      cff2=1.D0/6.D0
      if (ua .ge. 0.D0) then
        IS0 = cff1*(q_im3 - 2.D0*q_im2 +   q_im1)**2
     &      + 0.25D0*(q_im3 - 4.D0*q_im2 + 3*q_im1)**2
        IS1 = cff1*(q_im2 - 2.D0*q_im1 +   q_i  )**2
     &      + 0.25D0*(q_im2            -   q_i  )**2
        IS2 = cff1*(q_im1 - 2.D0*q_i   +   q_ip1)**2
     &   + 0.25D0*(3.D0*q_im1 - 4.D0*q_i   +   q_ip1)**2
        T5  = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
        p0  = cff2*(2.D0*q_im3 - 7.D0*q_im2 + 11.D0*q_im1)
        p1  = cff2*(  -q_im2 + 5.D0*q_im1 +  2.D0*q_i  )
        p2  = cff2*(2.D0*q_im1 + 5.D0*q_i   -     q_ip1)
        flux5_weno = w0*p0 + w1*p1 + w2*p2
      else
        IS0 = cff1*(q_ip2 - 2.D0*q_ip1 +   q_i  )**2
     &      + 0.25D0*(q_ip2 - 4.D0*q_ip1 + 3*q_i  )**2
        IS1 = cff1*(q_ip1 - 2.D0*q_i   +   q_im1)**2
     &      + 0.25D0*(q_ip1            -   q_im1)**2
        IS2 = cff1*(q_i   - 2.D0*q_im1 +   q_im2)**2
     &   + 0.25D0*(3.D0*q_i   - 4.D0*q_im1 +   q_im2)**2
        T5  = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
        p0  = cff2*(2.D0*q_ip2 - 7.D0*q_ip1 + 11.D0*q_i  )
        p1  = cff2*(  -q_ip1 + 5.D0*q_i   +  2.D0*q_im1)
        p2  = cff2*(2.D0*q_i   + 5.D0*q_im1 -     q_im2)
        flux5_weno = w0*p0 + w1*p1 + w2*p2
      endif
      return
      end
      function flux3_weno( q_im2, q_im1, q_i, q_ip1, ua)
!$acc routine(flux3_weno) seq
      implicit none
      REAL    :: flux3_weno
      REAL    :: q_im2, q_im1, q_i, q_ip1, ua
      REAL    :: IS0, IS1
      REAl    :: a0, a1
      REAL    :: w0, w1
      REAL    :: p0, p1
      REAL    :: Eps, d0,d1, T3
      Eps = 1.D-40
      d0=1.D0/3.D0
      d1=2.D0/3.D0
      if (ua .ge. 0.D0) then
        IS0 = (q_im1-q_im2)**2
        IS1 = (q_im1-q_i)**2
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1.D0/2.D0*(3.D0*q_im1-q_im2)
        p1  = 1.D0/2.D0*(q_im1+q_i)
        flux3_weno = w0*p0 + w1*p1
      else
        IS0 = (q_i-q_ip1)**2
        IS1 = (q_im1-q_i)**2
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1.D0/2.D0*(3.D0*q_i-q_ip1)
        p1  = 1.D0/2.D0*(q_im1+q_i)
        flux3_weno = w0*p0 + w1*p1
      endif
      return
      end
