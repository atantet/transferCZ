      SUBROUTINE BASIS_I(n,f_0,f_1,f_2,f_inv,x_i)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_1(0:n,0:n),f_2(0:n,0:n),f_inv(0:n,0:n)
      real x_i(0:n)
*     LOCAL
      integer mi,mx
      real pi,t_i
*     EXTERNAL
      real tsj
      EXTERNAL tsj
*
* definition of the grid of the variable y  & the chebs b
*

      pi = 4.*atan(1.)

      DO mx = 0,n
       t_i = pi*real(2*mx+1)/real(2*(n+1))
       x_i(mx) = - cos( t_i ) 
      END DO

      DO mx = 0,n
       DO mi = 0,n
        f_0(mi,mx) = tsj(x_i(mx),mi,0)
        f_1(mi,mx) = tsj(x_i(mx),mi,1)
        f_2(mi,mx) = tsj(x_i(mx),mi,2)
        f_inv(mi,mx) = f_0(mi,mx)
       END DO
      END DO

*     call normalize(n,f_0,f_1,f_2,f_inv)

      END 
*****************************************************************
      SUBROUTINE BASIS_E(n,f_0,f_1,f_2,f_inv,x_e)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_1(0:n,0:n),f_2(0:n,0:n),f_inv(0:n,0:n)
      real x_e(0:n)
*     LOCAL
      integer mi,mx
      real pi,t_e
*     EXTERNAL
      real tsj
      EXTERNAL tsj
*
* definition of the grid of the variable y  & the chebs b
*

      pi = 4.*atan(1.)

      DO mx = 0,n
       t_e = pi*real(mx)/real(n)
       x_e(mx) = - cos( t_e ) 
      END DO

      DO mx = 0,n
       DO mi = 0,n
        f_0(mi,mx) = tsj(x_e(mx),mi,0)
        f_1(mi,mx) = tsj(x_e(mx),mi,1)
        f_2(mi,mx) = tsj(x_e(mx),mi,2)
        f_inv(mi,mx) = f_0(mi,mx)
       END DO
      END DO

      call normalize_e(n,f_0,f_inv)

      END 
*****************************************************************
      SUBROUTINE BASIS_L(n,f_0,f_1,f_2,f_inv,y_l,L)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_1(0:n,0:n),f_2(0:n,0:n),f_inv(0:n,0:n)
      real L,y_l(0:n)
*     LOCAL
      integer mi,my
      real pi
      real t_l,x_l
*     EXTERNAL
      real tl
      EXTERNAL tl
*
* definition of the grid of the variable y  & the chebs b
*

      pi = 4.*atan(1.0)

      DO my = 0,n
          t_l = pi*real(2*my+1)/real(2*(n+1))
          x_l = - cos( t_l ) 
          y_l(my) = L*( 1. + x_l )/( 1. - x_l )
*         y_l(my) = L*cot(t_l/2.)**2
      END DO

      DO my = 0,n
       DO mi = 0,n
        f_0(mi,my) = tl(y_l(my),L,mi,0)
        f_1(mi,my) = tl(y_l(my),L,mi,1)
        f_2(mi,my) = tl(y_l(my),L,mi,2)
       END DO
      END DO

      call normalize(n,f_0,f_1,f_2,f_inv)

      END 
*****************************************************************
      SUBROUTINE BASIS_TB(n,f_0,f_1,f_2,f_inv,y_b,L)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_1(0:n,0:n),f_2(0:n,0:n),f_inv(0:n,0:n)
      real L,y_b(0:n)
*     LOCAL
      integer mi,my
      real pi
      real t_b,x_b
*     EXTERNAL
      real tb
      EXTERNAL tb
*
* definition of the grid of the variable y  & the chebs b
*

      pi = 4.*atan(1.0)

      DO my = 0,n
          t_b = pi*real(2*my+1)/real(2*(n+1))
          x_b = -cos( t_b ) 
          y_b(my) = L*x_b/sqrt(1. - x_b**2)
      END DO

      DO my = 0,n
       DO mi = 0,n
        f_0(mi,my) = tb(y_b(my),L,mi,0)
        f_1(mi,my) = tb(y_b(my),L,mi,1)
        f_2(mi,my) = tb(y_b(my),L,mi,2)
       END DO
      END DO

      call normalize(n,f_0,f_1,f_2,f_inv)

      END 
***********************************************************
      SUBROUTINE BASIS_SB(n,f_0,f_1,f_2,f_inv,y_b,L)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_1(0:n,0:n),f_2(0:n,0:n),f_inv(0:n,0:n)
      real L,y_b(0:n)
*     LOCAL
      integer mi,my
      real pi
      real t_b,x_b
*     EXTERNAL
      real sb
      EXTERNAL sb
*
* definition of the grid of the variable y  & the chebs b
*

      pi = 4.*atan(1.0)

      DO my = 0,n
          t_b = pi*real(2*my+1)/real(2*(n+1))
          x_b = -cos( t_b ) 
          y_b(my) = L*x_b/sqrt(1. - x_b**2)
      END DO

      DO my = 0,n
       DO mi = 0,n
        f_0(mi,my) = sb(y_b(my),L,mi,0)
        f_1(mi,my) = sb(y_b(my),L,mi,1)
        f_2(mi,my) = sb(y_b(my),L,mi,2)
       END DO
      END DO

      call normalize(n,f_0,f_1,f_2,f_inv)

      END 
***********************************************************
      SUBROUTINE BASIS_HE(n,f_0,f_1,f_i,w,y,L)

      implicit none
      
*     LOCAL
      integer my,mi,n,ifail
      real f_0(0:n,0:n),f_1(0:n,0:n),f_i(0:n,0:n),w(0:n),y(0:n),L
      real s_1,s_2,wd(n+1),yd(n+1)
*     EXTERNAL
      real phi
      EXTERNAL phi
*
      ifail=0
      CALL D01bcf(-4,0.,L**2,0.,0.,n+1,wd,yd,ifail)

      DO my=0,n
       w(my) = wd(my+1)
       y(my) = yd(my+1)
 	DO mi=0,n
        f_0(mi,my) = phi(L*y(my),mi)
       END DO
      END DO 

      DO my=0,n
      DO mi=0,n
       s_1 = sqrt(.5*real(mi))
       s_2 = sqrt(.5*real(mi+1))
       IF ((mi.NE.0).and.(mi.NE.n)) THEN
        f_1(mi,my) = s_1*f_0(mi-1,my)-s_2*f_0(mi+1,my)
       END IF
       IF (mi.EQ.0) THEN
        f_1(mi,my) = -s_2*f_0(mi+1,my)
       END IF
       IF (mi.EQ.n) THEN
        f_1(mi,my) = s_1*f_0(mi-1,my)
       END IF
      END DO
      END DO
      DO my=0,n
      DO mi=0,n
       f_1(mi,my) = L*f_1(mi,my)
       f_i(mi,my) = f_0(mi,my)
      END DO
      END DO
*

      END 
*****************************************************************
      FUNCTION phi(y,m)
*
      pi = 4.*atan(1.)
      s2 = sqrt(2.)
*
        y2 = y*y
        Hn = exp(-.5*y2)/pi**.25
        phi = Hn
        if (m.gt.0) then
          Hn1     = s2 * y * Hn
          if (m.gt.1) then
            do n=1,m-1
              r    = real(n+1)
              Hdum = Hn1
              Hn1  = sqrt(2./r) * y * Hdum -
     -                    sqrt(real(n)/r) * Hn
              Hn   = Hdum
            end do
          endif
          phi = Hn1
        endif
*
      end
*******************************************************
      SUBROUTINE normalize(n,f_0,f_1,f_2,f_inv)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_1(0:n,0:n),f_2(0:n,0:n),f_inv(0:n,0:n)
*     LOCAL
      integer mi,mx


      DO mx=0,n
      DO mi=0,n
       f_inv(mi,mx) = f_0(mi,mx)
       f_0(mi,mx) = f_0(mi,mx)/real(n+1)
       f_1(mi,mx) = f_1(mi,mx)/real(n+1)
       f_2(mi,mx) = f_2(mi,mx)/real(n+1)
      END DO
      END DO
      DO mx=0,n
      DO mi=1,n
       f_0(mi,mx) = 2.*f_0(mi,mx)
       f_1(mi,mx) = 2.*f_1(mi,mx)
       f_2(mi,mx) = 2.*f_2(mi,mx)
      END DO
      END DO

      END 
***********************************************************
      SUBROUTINE normalize_e(n,f_0,f_inv)

      implicit none
      
*     IMPORT?EXPORT
      integer n
      real f_0(0:n,0:n),f_inv(0:n,0:n)
*     LOCAL
      integer mi,mx

      DO mx=0,n
      DO mi=0,n
       f_inv(mi,mx) = 2.*f_0(mi,mx)/real(n)
      END DO
      END DO
      DO mx=0,n
       DO mi=0,n
        IF (mx.EQ.0) THEN
         f_inv(mi,mx) = f_inv(mi,mx)/2.
        END IF
        IF (mx.EQ.n) THEN
         f_inv(mi,mx) = f_inv(mi,mx)/2.
        END IF
        IF (mi.EQ.n) THEN
         f_inv(mi,mx) = f_inv(mi,mx)/2.
        END IF
        IF (mi.EQ.0) THEN
         f_inv(mi,mx) = f_inv(mi,mx)/2.
        END IF
       END DO
      END DO
      END 
*****************************************************************
      FUNCTION sb(y,L,p,der)

      implicit none
*     INPUT/OUTPUT
      real y,L,sb
      integer p,der
*     LOCAL
      real a,t,s,c,n

      a = y/sqrt(y*y + L*L)
      t = acos(a)
      s = sin(t)
      c = cos(t)
      n = real(p)

      IF (der.EQ.0) THEN
	sb = sin(n*t)
      END IF
*
      IF (der.EQ.1) THEN
	sb = -n*cos(n*t)*s**2/L
      END IF
*
      IF (der.EQ.2) THEN
	sb = -2.*c*cos(n*t) + n*s*sin(n*t)
       sb = -n*sb*s**3/L**2
      END IF

      END 
*****************************************************************
      FUNCTION tb(y,L,p,der)

      implicit none
*     INPUT/OUTPUT
      real y,L,tb
      integer p,der
*     LOCAL
      real a,t,s,c,n

      a = y/sqrt(y*y + L*L)
      t = acos(a)
      s = sin(t)
      c = cos(t)
      n = real(p)

      IF (der.EQ.0) THEN
	tb = cos(n*t)
      END IF
*
      IF (der.EQ.1) THEN
	tb = n*sin(n*t)*s**2/L
      END IF
*
      IF (der.EQ.2) THEN
	tb = 2.*c*sin(n*t) + n*s*cos(n*t)
       tb = -n*tb*s**3/L**2
      END IF

      END 
*****************************************************************
      FUNCTION tl(y,L,p,der)

      implicit none
*     INPUT/OUTPUT
      real y,L,tl
      integer p,der
*     LOCAL
      real t,s,c,n
*     EXTERNAL 
      real acot

      t = 2.*acot(sqrt(y/L))
      s = sin(t/2.)
      c = cos(t/2.)
      n = real(p)

      IF (der.EQ.0) THEN
	tl = cos(n*t)
      END IF
*
      IF (der.EQ.1) THEN
	tl = n*sin(n*t)*s**3/(c*L)
      END IF
*
      IF (der.EQ.2) THEN
	tl = 2.*c*s*n*cos(n*t) + (3. - 2.*s*s )*sin(n*t)
       tl = - n*tl*s**5/(2.*L**2*c**3)
      END IF

      END 
***************************************************************
      function tsj(x,l,n)

      implicit none
*     INPUT/OUTPUT
      real x,tsj
      integer l,n
*     LOCAL
      real pi,s,c,tnt,t

      pi = 4.*atan(1.0)
      t = acos(x)
      s = sin(t)
      c = cos(t)

      IF (n.EQ.0) THEN
       tsj=cos(l*t)
      END IF

      IF (n.EQ.1) THEN

       IF (x.eq.1.0) THEN
        tsj = real(l*l)
       ELSE IF (x.EQ.-1.0) THEN
        tsj = real(l*l)*cos((l+1)*pi)
       ELSE
        tnt=-l*sin(l*t)
        tsj=-tnt/s
       END IF

      END IF

      IF (n.EQ.2) THEN

       IF (x.eq.-1.0) THEN
        tsj = (l**4 - l**2)*cos(l*pi)/3.
       ELSE IF (x.eq.1.0) THEN
        tsj = (l**4 - l**2)/3.
       ELSE
        tnt = l*sin(l*t)*c - l*l*cos(l*t)*s
        tsj = tnt/s**3
       END IF

      END IF

      END
*****************************************************************
      FUNCTION acot(x)
      implicit none
*     INPUT/OUTPUT
      real    x,acot
*     LOCAL
      acot = atan(1/x)
*
      END
