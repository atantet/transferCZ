      SUBROUTINE init()
     
      call basis
      call force
      call surf_ext
      call atmosphere
      call surfmix
      call ocean

      END 
********************************************************************
      SUBROUTINE parameters
      
      include 'const.com'
      
      CHARACTER*6 DSTR

! Read experiment design
      OPEN(UNIT=70, FILE='design.par')
      READ(70, 30) DSTR, mu_0 ! Coupling
      READ(70, 30) DSTR, epsn ! Noise level
      CLOSE(70)
 30   FORMAT(A6, F18.10)
      PRINT *, 'Experiment design:'
      PRINT *, 'mu_0 = ', mu_0
      PRINT *, 'epsn = ', epsn
      
      pi=4.*atan(1.0)
*	*** Dimensional parameters ***
*	Basic parameters
      dpar(1) = 1.5e7 				!L [m]
      dpar(2) = 200.				!H [m]
      dpar(3) = 50.				!H1 [m]
      dpar(4) = dpar(2)-dpar(3)		        !H2=H-H1 [m]
      dpar(5) = 2.				!c_o [m/s]
      dpar(6) = 30.				!c_a [m/s]
      dpar(7) = 1.0922666667e-2		        !tau_0 [N/m2]
      dpar(8) = 1024.				!rho [kg/m3]		
      dpar(9) = 2.2e-11				!beta_0 [1/ms]
      dpar(10) = SQRT(dpar(5)/dpar(9))          ! Ly=sqrt(c_o/beta_0) [m]
      dpar(11) = 1.507555e7			! La=Ly/alpha, with alpha=0.2 (par(29)) [m]
      dpar(12) = 1.				!delta_T [K]
      dpar(13) = 30.				!T_0 [K]
*	Damping parameters
      dpar(14) = 1.333333333333e-8		!a_m, ocean damping [1/s]
      dpar(15) = 5e-6				!a_s, ocean surface layer damping [1/s]
      dpar(16) = 9.253333333333e-8		!a_T, ocean SST damping [1/s]
      dpar(17) = 2.5e-6				!a_M, atmosphere damping [1/s]	
*	T_sub - parametrization 
      dpar(18) = 39.92015968038723		!h1 [m]
      dpar(19) = 19.960079840319361		!h0 [m]
      dpar(20) = 23.				!T_s0 [K]
      dpar(21) = 27.5				!T_c [K]
*	Linear relation between SST and heat flux, surface wind and wind stress
      dpar(22) = 5.4e-3				! alpha_T, [m2/K s3]
      dpar(23) = 1.0				! ???? gamma_tau [1/s]

*	*** Nondimensional parameters ***  
*   parameters ocean component ( 1 - 4 )
      par(1)  = (dpar(1)*dpar(14))/dpar(5)
      par(2)  = 1.0						! delta
!      par(3)  = dpar(22)*dpar(12)*dpar(1)/(dpar(6)**3)	        !mu_0 = alpha_T delta_T L / ca^3, standard =3.0
!      par(3)  = 3*(dpar(1)/1.5e7)**2				! mu scaled with L^2
      par(3) = mu_0
      par(4)  = (dpar(1)*dpar(7))/(dpar(5)**2*dpar(8)*dpar(2))  ! F0=L tau0/(co^2 rho H)
!     par(4)  = 0.2*(dpar(1)/1.5e7)
*   parameters surface layer ( 5 - 8 )
      par(5)  = 0.3						! delta_s
      par(6)  = dpar(4)/dpar(3)				        ! gamma_s=H2/H1
      par(7)  = (dpar(1)*dpar(15))/(dpar(5))		        ! eps_s
      par(8)  = dpar(10)/dpar(1)			        ! Lambda_s
*   parameters SST component ( 9 - 17 )
      par(9)  = 1.0
      par(10) = 1.0
      par(11) = 5.0e-02						! delta_w
      par(12) = (dpar(1)*dpar(16))/dpar(5)			! eps_T
      par(13) = 1.0						! alpha_w = H1/Hu
*   T_sub - parametrization
      par(14) = dpar(19)/dpar(18)				! eta2 = h0/h1
      par(15) = dpar(2)/dpar(18)				! eta1 = H/h1
      par(16) = dpar(20)					! Ts0
      par(17) = dpar(13)					! T0
      par(18) = 1.0						! delta_sst
      par(19) = 0.0
      par(20) = 0.0
      par(21) = dpar(21)					!Tc
      par(22) = 0.76						!umin
      par(23) = 0.0
      par(24) = 1./par(22)
*   parameters atmosphere ( 26 - 30 )
      par(25) = 0.0
      par(26) = 0.0
      par(27) = 0.0
      par(28) = 1.0						!alpha_ext
      par(29) = 0.2						! alpha = Ly/La
!     par(29) = apar(3)/apar(4)
      par(30) = (dpar(1)*dpar(17))/(dpar(6))		        ! eps_a
      par(31) = epsn ! stochastic wind forcing noise level (A. Tantet)
         
      END
********************************************************************

********************************************************************
      SUBROUTINE read_tau
      
      include 'const.com'
      
      COMMON /EXT_FILE/ F0_t, F0_t_dx, F0_t_dy
      real(kind=8) F0Read(nx+1,ny+1)
      integer i_time

*     Read means
c$$$      print *, 'Reading mean fields...'
c$$$      open (20,
c$$$     &     file='../init/tau_mean_amp30.bin',
c$$$     &     status='old', form='unformatted')
c$$$      Do i_time = 1, maxstep
c$$$         read (20) F0_t(:, :, i_time)
c$$$      END DO
c$$$      close(20)
c$$$      open (20,file='../init/dtaudx_mean_amp30.bin',
c$$$     &     status='old', form='unformatted')
c$$$      Do i_time = 1, maxstep
c$$$         read (20) F0_t_dx(:, :, i_time)
c$$$      END DO
c$$$      close(20)
c$$$      open (20,file='../init/dtaudy_mean_amp30.bin',
c$$$     &     status='old', form='unformatted')
c$$$      Do i_time = 1, maxstep
c$$$         read (20) F0_t_dy(:, :, i_time)
c$$$      END DO
c$$$      close(20)

*     Read fluctuations
      print *, 'Reading fluctuations...'
      open (20,
     &     file='../init/tauWhite_1eofs_amp30.bin',
     &     status='old', form='unformatted')
      Do i_time = 1, maxstep
         read (20) F0Read
         F0_t(:, :, i_time) = par(31) * F0Read
      END DO
      close(20)
      open (20,file='../init/dtaudxWhite_1eofs_amp30.bin',
     &     status='old', form='unformatted')
      Do i_time = 1, maxstep
         read (20) F0Read
         F0_t_dx(:, :, i_time) = par(31) * F0Read
      END DO
      close(20)
      open (20,file='../init/dtaudyWhite_1eofs_amp30.bin',
     &     status='old', form='unformatted')
      Do i_time = 1, maxstep
         read (20) F0Read
         F0_t_dy(:, :, i_time) = par(31) * F0Read
      END DO
      close(20)
      print *, 'Done.'
           
      END
********************************************************************
      SUBROUTINE writepar
      include 'const.com'

      write(1,*) nx+1,ny+1,nbm,ndim,na+1
      write(1,*) 'NEW',maxstep,dstep
      write(1,*) ' ********************************************* '
      write(1,*) ' Dimensional parameters'
      write(1,*) ' ********************************************* '
      write(1,*) ' Basic parameters '
      write(1,*) ' L		   = ',01,dpar(1),' m'
      write(1,*) ' H           = ',02,dpar(2),' m'
      write(1,*) ' H1          = ',03,dpar(3),' m'
      write(1,*) ' H2          = ',04,dpar(4),' m'
      write(1,*) ' c_o         = ',05,dpar(5),' m/s'
      write(1,*) ' c_a         = ',06,dpar(6),' m/s'
      write(1,*) ' tau_0       = ',07,dpar(7),' N/m^2'
      write(1,*) ' rho         = ',08,dpar(8),' kg/m^3'
      write(1,*) ' beta_0      = ',09,dpar(09),' 1/(ms)'
      write(1,*) ' Ly          = ',10,dpar(10),' m'
      write(1,*) ' La          = ',11,dpar(11),' m'
      write(1,*) ' deltaT      = ',12,dpar(12),' K'
      write(1,*) ' T0          = ',13,dpar(13),' K'
      write(1,*) ' Damping parameters '	  
      write(1,*) ' a_m         = ',14,dpar(14),' 1/s'
      write(1,*) ' a_s         = ',15,dpar(15),' 1/s'
      write(1,*) ' a_T         = ',16,dpar(16),' 1/s'
      write(1,*) ' a_M         = ',17,dpar(17),' 1/s'
      write(1,*) ' T_sub - paramterization '	  	  
      write(1,*) ' h1          = ',18,dpar(18),' m'
      write(1,*) ' h0          = ',19,dpar(19),' m'
      write(1,*) ' T_s0        = ',20,dpar(20),' K'
      write(1,*) ' Tc          = ',21,dpar(21),' K'
      write(1,*) ' Linear relation between SST and heat flux, alpha_T'
      write(1,*) ' alpha_T     = ',22,dpar(22),' m^2/(K s^3)'
      write(1,*) ' Linear relation between surface wind and wind stress' 
      write(1,*) ' gamma_tau   = ',23,par(23),' 1/s'
      write(1,*) ' ********************************************* '
      write(1,*) ' Non-dimensional parameters'
      write(1,*) ' ********************************************* '
      write(1,*) ' parameters of ocean component '
      write(1,*) ' eps_o       = ',01,par(1)
      write(1,*) ' delta       = ',02,par(2)
      write(1,*) ' mu          = ',03,par(3)
      write(1,*) ' F_0         = ',04,par(4)
      write(1,*) ' parameters of surface layer component '
      write(1,*) ' delta_s     = ',05,par(5)
      write(1,*) ' gamma_s     = ',06,par(6)
      write(1,*) ' eps_s       = ',07,par(7)
      write(1,*) ' Lambda_s    = ',08,par(8)
      write(1,*) ' parameters of SST component '
      write(1,*) ' ad_u        = ',09,par(09)
      write(1,*) ' ad_v        = ',10,par(10)
      write(1,*) ' delta_w     = ',11,par(11)
      write(1,*) ' eps_T       = ',12,par(12)
      write(1,*) ' alpha_w     = ',13,par(13)
      write(1,*) ' H0          = ',14,par(14)
      write(1,*) ' href        = ',15,par(15)
      write(1,*) ' Ts0         = ',16,par(16)
      write(1,*) ' T0          = ',17,par(17)
      write(1,*) ' delta_sst   = ',18,par(18)
      write(1,*) ' free        = ',19,par(19)
      write(1,*) ' free        = ',20,par(20)
      write(1,*) ' Tc          = ',21,par(21)
      write(1,*) ' umin        = ',22,par(22)
      write(1,*) ' sw_E        = ',23,par(23)
      write(1,*) ' eps_E       = ',24,par(24)
      write(1,*) ' parameters of atmosphere component '
      write(1,*) ' k_Q         = ',25,par(25)
      write(1,*) ' sw_Q        = ',26,par(26)
      write(1,*) ' F_1         = ',27,par(27)
      write(1,*) ' alpha_ext   = ',28,par(28)
      write(1,*) ' alpha       = ',29,par(29)
      write(1,*) ' eps_a       = ',30,par(30)

	 	  
      END
*******************************************************
      SUBROUTINE PARMZ(w,th,mhf,tsf)

      implicit none
      include 'const.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real w(0:nx,0:ny),mhf(0:nx,0:ny,0:3)
      real th(0:nx,0:ny),tsf(0:nx,0:ny,0:3)
*     LOCAL
      integer mx,my
      real tanh1,eps,x,dtemp
      real Ts0,h0,hhat,f0,f1,f2,f3

      eps = par(11)
      mhf = 0.0

      DO mx=0,nx
      DO my=0,ny
       x = w(mx,my)
       tanh1 = tanh(x/eps)
       f0 = .5*(1.+tanh1)
       f1 = .5*(1.-tanh1**2)/eps
       f2 = -tanh1*(1.-tanh1**2)/eps**2
       f3 = -(1.-tanh1**2)*(1.-3.*tanh1**2)/eps**3
       mhf(mx,my,1) = f0 + x*f1
       mhf(mx,my,0) = x*f0
       mhf(mx,my,2) = 2.*f1 + x*f2
       mhf(mx,my,3) = 3.*f2 + x*f3
      END DO
      END DO

      Ts0  = par(16)
      hhat = 1.0/par(15) 
      h0   = hhat*par(14)  
*     tanh1=tanh((z + h0)/hhat)
      DO mx=0,nx
      DO my=0,ny
       tanh1=tanh(par(15)*th(mx,my) + par(14))
       dtemp = T0(0,mx,my)-Ts0
       tsf(mx,my,0) = Ts0 + dtemp*tanh1
       tsf(mx,my,1) = dtemp*(1.-tanh1**2)/hhat
       tsf(mx,my,2) = -2.*dtemp*(1.-tanh1**2)*tanh1/hhat**2
       tsf(mx,my,3) = -2.*dtemp*(1.-tanh1**2)*(1.-3.*tanh1**2.)/hhat**3
      END DO
      END DO

      END

*************************************************************************
      SUBROUTINE force
      
      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'tau.com'

*     LOCAL
      integer mx,my,j,i,k
      real a2
*     COMMON
      integer flag_nt, k_solve

      COMMON /EXT_FILE/ F0_t, F0_t_dx, F0_t_dy
      COMMON /EXT_FLAG/flag_nt, k_solve
      
*     Set tau to 0, otherwise, we cumulate the external forcing
      taux = 0.
      tauy = 0.
      call xstress(taux)
*      call ystress(tauy)

      if (k_solve.eq.0) then
*     Temperature at Radiative Equilibrium
         print *, "Initializing temperature at radiative equilibrium."
         a2 = par(29)**2/2.
         TRE(0, :, :) = par(17)
         TRE(1, :, :) = 0.
         TRE(2, :, :) = 0.
      END IF

      
      if (flag_nt.eq.1) then
         if (k_solve.EQ.1) then
*     Read record of residual red-noise wind stress
            print *, "Reading forcing records..."
            call read_tau
         endif

*     Define the external wind stress as the sum
*     of the artificial meridional one and the stochastic residual one
         print *, 'Assigning F0_t at ', k_solve, ' to taux'
         DO mx = 1, nx+1
            DO my = 1, ny+1
               taux(mx, my, 0) = taux(mx,my,0)+F0_t(mx, my, k_solve)
               taux(mx, my, 1) = taux(mx,my,1)+F0_t_dx(mx, my, k_solve)
               taux(mx, my, 2) = taux(mx,my,2)+F0_t_dy(mx, my, k_solve)
            END DO
         END DO

         if (ioute.eq.1) then
            write(99,*) 'wind fields, y(my),taux(5,my,0),taux(5,my,2) 
     &written to fort.76'
            DO my=0,ny
               write(76,001) y(my),taux(5,my,0),taux(5,my,2)
            END DO
         endif 

*     Define tau_x ?
         tau_x = 0.0
         DO j=0,ny
            DO mx=0,nx
               DO my=0,ny
                  tau_x(mx,j) = tau_x(mx,j)
     &                 + taux(mx,my,0)*psi(j,my)*wh(my)
               END DO
            END DO
         END DO
      END IF

001   format(3(1x,f13.6))
      END
*****************************************************
      SUBROUTINE surf_ext
*     Calculation of the surface velocities in terms of wind stress forcing,
*     see VdVaart & Dijkstra (2000), equation A10, gamma_s is applied after.

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'tau.com'
*     LOCAL
      integer mx,my
      real u_sx,v_sy,L_s,eps_s,time
*     EXTERNAL
      real di,dm

*     eps_s = L * a_s / c_0
      eps_s = par(7)
*     Ls = \Lambda_s = L_y / L = sqrt(c_0 / \beta_0) / L
      L_s = par(8)

      u_s_ext = 0.0
      v_s_ext = 0.0
      w_s_ext = 0.0
      u_sx = 0.0
      v_sy = 0.0
      
      DO my=0,ny
         di = 1./((L_s*eps_s)**2 + y(my)**2)
         dm = (L_s*eps_s)**2 - y(my)**2
         DO mx=0,nx
*     Ekman zonal velocity
            u_s_ext(mx,my) = L_s * di * (L_s * eps_s * taux(mx,my,0)
     &           + y(my) * tauy(mx,my,0))
*     Ekman meridional velocity
            v_s_ext(mx,my) = di * (L_s * eps_s * tauy(mx,my,0) 
     &           - y(my) * taux(mx,my,0))
*     x derivative of u_s
            u_sx = L_s * di * (L_s * eps_s * taux(mx,my,1) 
     &           + y(my) * tauy(mx,my,1))
*     y derivative of v_s
            v_sy = di * (-2. * y(my) * L_s * eps_s * di * tauy(mx,my,0) 
     &           - di*dm*taux(mx,my,0) + L_s * eps_s * tauy(mx,my,2)
     &           - y(my) * taux(mx,my,2))
*     Ekman vertical velocity, there should be a H1 / H somewhere?
            w_s_ext(mx,my) = u_sx + v_sy
         END DO
      END DO

      END
*****************************************************
      SUBROUTINE xstress(xf)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      real xf(0:nx,0:ny,0:2)
*     LOCAL
      integer mx,my
      real f_x,yb,a2,alpha,alpha_ext,F_1,F_0,eg

      alpha     = par(29)
      alpha_ext = 1.
      F_1       = par(27)
      F_0       = par(4)
      a2 = alpha_ext*alpha**2
      
      DO my=0,ny
         yb = y(my)
         eg = exp(-yb*yb*a2/2.)
         DO mx=0,nx
            xf(mx,my,0) = -(F_0 + F_1 * yb**2) * eg
            xf(mx,my,1) =   0.0
            xf(mx,my,2) = -eg*((2.*F_1 -a2*F_0)*yb - a2*F_1*yb**3) 
         END DO
      END DO
      
      END
*****************************************************

      SUBROUTINE ystress(yf)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      real yf(0:nx,0:ny,0:2),a2,yb,alpha,alpha_ext
      integer mx,my

      yf = 0.0
      a2 = par(29)**2
      DO my=0,ny
       yb = y(my)
       DO mx=0,nx
        yf(mx,my,0) = 0.
        yf(mx,my,1) = 0.
        yf(mx,my,2) = 0.
       END DO
      END DO

      END
*******************************************************************
      SUBROUTINE BASIS

      implicit none
      include 'const.com'
      include 'coeff.com'
      
*     LOCAL
      integer mx,my,ix,iy
      real qxx(0:nx,0:nx)
      real dx
*
      call BASIS_E(nx,q,qx,qxx,q_i,x)
      call BASIS_HE(ny,psi,psi_y,psi_i,wh,y,1.0)
*
      x = (x + 1.)/2.

      dx  = 2.
      qx  = dx*qx
      qxx = dx*dx*qxx

      DO mx=0,nx
         DO my=0,ny
            DO ix=0,nx
               DO iy=0,ny
                  f(mx,my,ix,iy)   = q(ix,mx)*psi(iy,my)
                  fx(mx,my,ix,iy)  = qx(ix,mx)*psi(iy,my)
*     fxx(mx,my,ix,iy) = qxx(ix,mx)*psi(iy,my)
                  fy(mx,my,ix,iy)  = q(ix,mx)*psi_y(iy,my)
                  f_i(mx,my,ix,iy) = q_i(ix,mx)*psi_i(iy,my)
               END DO
            END DO
         END DO
      END DO
      if (ioute.eq.1) then
         write(99,*) 'basis coefficients written to fort.21'
         write(21,*) nx+1,ny+1
         DO ix=0,nx
            DO iy=0,ny
               write(21,001) real(ix),real(iy),f(nx-3,ny/2,ix,iy)
            END DO
         END DO
      endif 

      EB = 0.
      EB(0) = 1.
      EB(1) = 0.
      DO iy=2,ny
         EB(iy) = EB(iy-2)*sqrt(real(iy-1)/real(iy))
      END DO
 002  format(4(f10.4))
 001  format(4(1x,e12.5))
      END 
***************************************************************
      SUBROUTINE atmosphere

      implicit none
      include 'const.com'
      include 'coeff.com'
      
*     LOCAL
      integer mx,ms,my,ix,iy,n,j,merws
      real c(0:ny,0:na+2)
      real cinv(0:ny,0:na+2),pck(0:ny,0:ny)
      real AI(0:nx,0:nx,0:na+2),delta(0:nx+1)
      real R(0:nx,0:nx,0:na+2,0:ny),arr(0:ny)
      real RX(0:nx,0:nx,0:na+2,0:ny)
      real S,SX,sqn,fcn,sumair,sumairx,sumairy
      real vn,vnx,sn,sumbir,sumbirx,sumbiry
      real eps_a,alpha,s_1,s_2,expx,exp0,exp1
      real psia(0:na+2,0:ny),psia_y(0:na,0:ny)
      real phi,dlt1,dlt2,qdl
*
      alpha = par(29)
      eps_a = par(30)
      DO my = 0,ny
      DO n  = 0,na+2
       psia(n,my) = phi(alpha*y(my),n)
      END DO
      END DO
*     
      DO my =0,ny
         DO n =0,na
            s_1 = sqrt(.5*real(n))
            s_2 = sqrt(.5*real(n+1))
            IF ((n.NE.0).and.(n.NE.na)) THEN
               psia_y(n,my) = s_1*psia(n-1,my)-s_2*psia(n+1,my)
            END IF
            IF ((n.EQ.na).and.(n.NE.0)) THEN
               psia_y(n,my) = s_1*psia(n-1,my)
            END IF
            IF (n.EQ.0) THEN
               psia_y(n,my) = -s_2*psia(n+1,my)
            END IF
         END DO
      END DO
      psia_y = alpha*psia_y
*     
      c = 0.0
      DO iy = 0,ny
         DO n  = 0,na+2
            DO my=0,ny
               c(iy,n) = c(iy,n) + psi(iy,my)*psia(n,my)*wh(my)
            END DO
         END DO
      END DO
      cinv = alpha*c
*
      DO mx=1,nx
         delta(mx) = x(mx) - x(mx-1)
      END DO
*     
      AI = 0.0
      DO ix = 0,nx
*     creation of I(m,0)(x)
         AI(0,ix,0) = 0.0
         DO mx = 1,nx
            exp1 = exp( eps_a*x(mx))
            exp0 = exp( eps_a*x(mx-1))
            AI(mx,ix,0) = AI(mx-1,ix,0) 
     &           + delta(mx)*(q(ix,mx)*exp1 + q(ix,mx-1)*exp0)/2.
         END DO
      END DO
      
      DO ix = 0,nx
*     creation of I(m,1)(x) = 0.0
         DO mx = 1,nx
            AI(mx,ix,1) = 0.0
         END DO
      END DO
      
*     creation of I(m,n)(x)
      DO n = 2,na+2
         DO ix = 0,nx
*     creation of I(m,n)(x)
            AI(nx,ix,n) = 0.0
            DO mx = 1,nx
               ms = nx - mx
               exp1 = exp(-real(2*n-1)*eps_a*x(ms+1))
               exp0 = exp(-real(2*n-1)*eps_a*x(ms))
               AI(ms,ix,n) = AI(ms+1,ix,n) 
     &              + delta(ms+1)*(q(ix,ms+1)*exp1 + q(ix,ms)*exp0)/2.
            END DO
         END DO
      END DO

*     creation of R_n(x) in u_a = \sum( R_n ... )\psi_n
      R = 0.0
      RX = 0.0
      DO mx = 0,nx
         expx = exp( -eps_a*x(mx))
         DO ix = 0,nx
            DO iy = 0,ny
               R(mx,ix,0,iy)  = -expx*AI(mx,ix,0)*cinv(iy,0)/2.
               R(mx,ix,1,iy)  = 0.0
               RX(mx,ix,0,iy) = -R(mx,ix,0,iy)*eps_a 
     &              -q(ix,mx)*cinv(iy,0)/2.
               RX(mx,ix,1,iy) = 0.0
            END DO
         END DO
      END DO
*
      DO n = 2,na+2
         sqn = sqrt(real(n)/real(n-1))
         DO mx = 0,nx
            expx = exp( real(2*n-1)*eps_a*x(mx))
            DO ix = 0,nx
               S  = expx*AI(mx,ix,n)
*        SX = S*eps_a*real(2*n-1) - q(ix,mx)
               DO iy = 0,ny
                  fcn = real(n-1)*(cinv(iy,n) + sqn*cinv(iy,n-2))/2.
                  R(mx,ix,n,iy)   = -S*fcn
*     RX(mx,ix,n,iy)  = -SX*fcn
               END DO
            END DO
         END DO
      END DO
*
      DO n = 2,na+2
         DO ix = 0,nx
            DO iy = 0,ny
               dlt1 = x(1)-x(0)
               dlt2 = x(2)-x(0)
               qdl = dlt1/dlt2
               RX(0,ix,n,iy) = (R(1,ix,n,iy)/qdl - qdl*R(2,ix,n,iy) 
     &              + (qdl - 1./qdl)*R(0,ix,n,iy))/(x(2)-x(1))
               DO mx = 1,nx-1
                  dlt1 = x(mx+1)-x(mx)
                  dlt2 = x(mx)-x(mx-1)
                  qdl = dlt1/dlt2
                  RX(mx,ix,n,iy) = (R(mx+1,ix,n,iy)/qdl
     &                 - qdl*R(mx-1,ix,n,iy) 
     &                 + (qdl - 1./qdl)*R(mx,ix,n,iy))/(x(mx+1)-x(mx-1))
               END DO
               dlt1 = x(nx)-x(nx-1)
               dlt2 = x(nx)-x(nx-2)
               qdl = dlt1/dlt2
               RX(nx,ix,n,iy) = -(R(nx-1,ix,n,iy)/qdl 
     &              - qdl*R(nx-2,ix,n,iy) 
     &              + (qdl - 1./qdl)*R(nx,ix,n,iy))/(x(nx-1)-x(nx-2))
            END DO
         END DO
      END DO

*     creation of air_mn(x,y) in u_a = \sum( air_mn(x,y) T_mn)
      air  = 0.0
      airx = 0.0
      DO mx = 0,nx
         DO my = 0,ny
            DO ix = 0,nx
*     DO iy = 0,ny
               DO iy = 0,5
                  sumair  = 0.0
                  sumairx  = 0.0
                  sumairy  = 0.0
                  DO n = 0,na
                     sqn = sqrt(real(n+2)/real(n+1))
                     sumair  = sumair  + (R(mx,ix,n,iy) 
     &                    - sqn*R(mx,ix,n+2,iy))*psia(n,my) 
                     sumairx = sumairx + (RX(mx,ix,n,iy) 
     &                    - sqn*RX(mx,ix,n+2,iy))*psia(n,my) 
                     sumairy = sumairy + (R(mx,ix,n,iy) 
     &                    - sqn*R(mx,ix,n+2,iy))*psia_y(n,my) 
                  END DO
                  air(mx,my,ix,iy)   = sumair
                  airx(mx,my,ix,iy)  = sumairx
                  airy(mx,my,ix,iy)  = sumairy
               END DO 
            END DO 
         END DO 
      END DO 
      if (ioute.eq.1) then 
         write(99,*) 'atm. basis coefficients written to fort.20'
         write(20,*) nx+1,ny+1
         DO ix=0,nx
            DO iy=0,ny
               write(20,001) real(ix),real(iy),air(nx-3,ny/2,ix,iy)
            END DO
         END DO
      endif 

*     creation of bir_mn(x,y) in v_a = \sum( bir_mn(x,y) T_mn)
      bir  = 0.0
      birx = 0.0
      biry = 0.0
      merws = 0
      IF (merws.EQ.1) THEN
         DO mx = 0,nx
            DO my = 0,ny
               DO ix = 0,nx
                  DO iy = 0,ny
                     sn = sqrt(1./2.)
                     vn = q(ix,mx)*cinv(iy,1)
                     vnx = qx(ix,mx)*cinv(iy,1)
                     sumbir  = sn*vn*psia(0,my)
                     sumbirx = sn*vnx*psia(0,my)
                     sumbiry = sn*vn*psia_y(0,my)
                     DO n = 1,na
                        sqn = sqrt(real(n)/real(n+1))
                        sn  = 2.*sqrt(2.*(n+1))
                        vn  = eps_a*R(mx,ix,n+1,iy) 
     &                       + q(ix,mx)*(cinv(iy,n+1)
     &                       +sqn*cinv(iy,n-1))/4.
                        vnx = eps_a*RX(mx,ix,n+1,iy) 
     &                       + qx(ix,mx)*(cinv(iy,n+1)
     &                       +sqn*cinv(iy,n-1))/4.
                        sumbir   = sumbir   + sn*vn*psia(n,my)
                        sumbirx  = sumbirx  + sn*vnx*psia(n,my)
                        sumbiry  = sumbiry  + sn*vn*psia_y(n,my)
                     END DO
                     bir(mx,my,ix,iy)   = sumbir
                     birx(mx,my,ix,iy)  = sumbirx
                     biry(mx,my,ix,iy)  = sumbiry
                  END DO 
               END DO 
            END DO 
         END DO 
      END IF

*     creation of ca_mn_j(x,y) in u_a = \sum( ca_mn_j(x)\psi_j(y) q_mn)
*     creation of cb_mn_j(x,y) in v_a = \sum( ba_mn_j(x)\psi_j(y) q_mn)
      ca = 0.0
      cb = 0.0
      DO j = 0,ny
         DO mx = 0,nx
            DO ix = 0,nx
               DO iy = 0,ny
                  DO my = 0,ny
                     ca(mx,j,ix,iy) = ca(mx,j,ix,iy) 
     &                    + air(mx,my,ix,iy)*psi(j,my)*wh(my)
*     cb(mx,j,ix,iy) = cb(mx,j,ix,iy) 
     &                    + bir(mx,my,ix,iy)*psi(j,my)*wh(my)
                  END DO
               END DO
            END DO
         END DO
      END DO
      
001   format(4(1x,e12.5))
      END
***************************************************************
      SUBROUTINE surfmix 

      implicit none
      include 'const.com'
      include 'coeff.com'
      
*     LOCAL
      integer mx,my,ix,iy
      real usx,vsy,di,dm,L_s,eps_s
*
      eps_s = par(7)
      L_s = par(8)
*     creation of us_mn(x,y) in u_smn = \sum( us_mn(x,y) T_mn)
*     creation of vs_mn(x,y) in v_smn = \sum( vs_mn(x,y) T_mn)
*     creation of ws_mn(x,y) in w_smn = \sum( ws_mn(x,y) T_mn)
      us = 0.0
      vs = 0.0
      ws = 0.0
      usx = 0.0
      vsy = 0.0
      DO mx = 0,nx
      DO my = 0,ny
       di = 1./((L_s*eps_s)**2 + y(my)**2)
       dm = (L_s*eps_s)**2 - y(my)**2
       DO ix = 0,nx
       DO iy = 0,ny
        us(mx,my,ix,iy)  = di*L_s*(eps_s*L_s*air(mx,my,ix,iy) 
     &                     + y(my)*bir(mx,my,ix,iy))
        vs(mx,my,ix,iy)  = di*(eps_s*L_s*bir(mx,my,ix,iy) 
     &                     - y(my)*air(mx,my,ix,iy))
        usx = di*L_s*(eps_s*L_s*airx(mx,my,ix,iy) 
     &                     + y(my)*birx(mx,my,ix,iy))
        vsy = di*(
     &  - 2.*y(my)*eps_s*L_s*di*bir(mx,my,ix,iy)
     &  - di*dm*air(mx,my,ix,iy)
     &  + eps_s*L_s*biry(mx,my,ix,iy)
     &  - y(my)*airy(mx,my,ix,iy) )

        ws(mx,my,ix,iy) = usx + vsy

       END DO 
       END DO 
      END DO 
      END DO 

      END
**************************************************
      SUBROUTINE ocean

      implicit none
      include 'const.com'
      include 'tau.com'
      include 'coeff.com'
      
*     LOCAL
      integer mx,my,ix,iy,j
      real sq,cr,fux,fvy

      fu = 0.0
      fh = 0.0
      fv = 0.0
      fvt = 0.0
      vtau = 0.0
      DO mx=0,nx
      DO my=0,ny
      DO ix=0,nx
*     coefficients : u(x,y) = \sum fu_mn(x,y) r_mn
*     coefficients : h(x,y) = \sum fh_mn(x,y) r_mn
*     coefficients : v(x,y) = \sum v_mn(x,y) r_mn + \sum vt_mn(x,y) q_mn
       fu(mx,my,ix,0) = f(mx,my,ix,0) 
       fu(mx,my,ix,1) = f(mx,my,ix,1)
       fh(mx,my,ix,0) = f(mx,my,ix,0) 
       fh(mx,my,ix,1) = f(mx,my,ix,1)
       fv(mx,my,ix,0) = 0.0
       fv(mx,my,ix,1) = 0.0
       DO iy=2,ny
        sq = sqrt(real(iy)/real(iy-1))
        cr = sqrt(8.*iy/real(2*iy-1))
        fu(mx,my,ix,iy) = f(mx,my,ix,iy) - sq*f(mx,my,ix,iy-2)
        fh(mx,my,ix,iy) = f(mx,my,ix,iy) + sq*f(mx,my,ix,iy-2)
        fv(mx,my,ix,iy) = fx(mx,my,ix,iy-1)*cr
       END DO
       fu(mx,my,ix,ny-1)=-sqrt(real(ny-1)/real(ny-2))*f(mx,my,ix,ny-3)
       fu(mx,my,ix,ny)  =-sqrt(real(ny)/real(ny-1))*f(mx,my,ix,ny-2)
      END DO
      END DO
      END DO

      DO mx=0,nx
      DO my=0,ny
      DO ix=0,nx
*     coefficients vt: v(x,y) = \sum v_mn(x,y) r_mn + \sum vt_mn(x,y) q_mn + vtau
      DO iy=0,ny
       fvt(mx,my,ix,iy) = ca(mx,1,ix,iy)*psi(0,my)
        DO j=1,ny
         fvt(mx,my,ix,iy) = fvt(mx,my,ix,iy) 
     &  + psi(j,my)*(sqrt(real(j+1))*ca(mx,j+1,ix,iy) 
     &  + sqrt(real(j))*ca(mx,j-1,ix,iy))/real(2*j+1)
        END DO
       END DO
      END DO
      END DO
      END DO
      fvt = fvt/sqrt(2.)
*    
      DO mx=0,nx
      DO my=0,ny
*     coefficients vt: v(x,y) = \sum v_mn(x,y) r_mn + \sum vt_mn(x,y) q_mn + vtau
       vtau(mx,my) = tau_x(mx,1)*psi(0,my)
        DO j=1,ny
         vtau(mx,my) = vtau(mx,my) 
     &  + psi(j,my)*(sqrt(real(j+1))*tau_x(mx,j+1) 
     &  + sqrt(real(j))*tau_x(mx,j-1))/real(2*j+1)
        END DO
      END DO
      END DO
      vtau = vtau/sqrt(2.)
*    
      END
