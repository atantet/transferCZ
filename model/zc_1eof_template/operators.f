*************************************************************************
* In this routine the operators defining the right hand 
* side of the equations and the Jacobian matrix are defined. 
* The state vector is indicated by z (dimension ndim) and the 
* right hand side by fvec (also dimension ndim). The state 
* vector contains the coefficients r_{ij} and q_{ij} obtained
* (see thesis VdVaart, page 102-105) 
* from collocation of the equations at collocation points of the 
* Chebychev (0, ..., nx) and Hermite functions (0, ..., ny) 
* and hence ndim = 2*(nx+1)*(ny+1); 
*  
* subroutines: 
*              rhs: right hand side of the equations 
*              matA: Jacobian matrix 
*              matB: coefficients before time derivative (stiffness matrix)
* 
*************************************************************************

      SUBROUTINE fields(z,zt,itime,time)
      
      implicit none
      include 'const.com'
      include 'outfld.com'
      include 'coeff.com'
      include 'tau.com'
*     INPUT/OUTPUT
      integer itime
      real zt(ndim),z(ndim),time
*     LOCAL
      integer j,mx,my,mj
      real delta_s,mu,gamma_s,alpha_w,F_0,eps_o,eps_E,adv,adu,eps_d
      real delta,sw_E,sw_Q,k_Q,a_h,delta_sst
      real u_s(0:nx,0:ny),v_s(0:nx,0:ny)
      real omega,gt,dgt,amp,ampfo,Mo
      real qf(ndim)

      eps_o     = par(1)
      delta     = par(2)
      mu        = par(3)
      F_0       = par(4)
      delta_s   = par(5)
      gamma_s   = par(6)
      adu       = par(9)
      adv       = par(10)
      alpha_w   = par(13)
      eps_d     = par(12)
      delta_sst = par(18)
      sw_E      = par(23)
      eps_E     = par(24)

      qf  = 0.0
      r_t = 0.0
      h_t = 0.0
      T_t = 0.0
      amp = 0.131
      omega = 6.28/4.2
      T0 = TRE
      T0_t = 0.0
      IF (itime.EQ.1) THEN
         gt  = amp*sin(omega*time)
         dgt = amp*omega*cos(omega*time)
         call meanft(zt,r_t,h_t,T_t)
c$$$         DO mx=0,nx
c$$$            DO my=0,ny
c$$$               T0(0,mx,my) = par(17) 
c$$$*     T0(0,mx,my)= par(17) + sin(omega*time)
c$$$*     T0(0,mx,my) = par(17)*cos(par(28)*y(my) - gt)
c$$$*     T0(1,mx,my) = 0.0
c$$$*     T0(2,mx,my) = -par(28)*par(17)*sin(par(28)*y(my) - gt)
c$$$*     T0_t(mx,my) = par(17)*sin(par(28)*y(my) - gt)*dgt
c$$$            END DO
c$$$         END DO
      END IF
*     
      call OBC(z,b0,b1,bcT)
      call meanf(z,r,u,vrs,vq,h,T)
      T   = T + T0
      T_t = T_t + T0_t
*
      IF (sw_A.EQ.1) THEN
         call force
         call surf_ext
         call atmosphere
         call surfmix
      END IF
*
      qf = z 
      call sky(qf,ua,va,u_A,v_A)
      call surface(qf,u_s,v_s,w_s)
      DO mx=0,nx
         DO j=0,ny
            atm(mx,j)  = tau_x(mx,j) + mu*ua(mx,j)
*     tau^x = \tau_{ext}^x + \tau_c^x
*     with \tau_c^x = \delta_\tau u_a (= mu u_a?, Heydt 2011)
            wind(mx,j) = mu*u_A(mx,j) + taux(mx,j,0)
         END DO
      END DO
      
*
      w = -(delta*h_t + eps_o*h)
      v = 0.0
      u = 0.0
      u_1 = u + gamma_s*( u_s_ext + delta_s*mu*u_s )
      v_1 = v + gamma_s*( v_s_ext + delta_s*mu*v_s )
      w_1 = w + gamma_s*( w_s_ext + delta_s*mu*w_s )
      w_1 = alpha_w*w_1

      call PARMZ(w_1,h,mheav,tsub)
*
      eps_T  = eps_d
      deps_T = 0.0

 001  format(2(1x,I3),6(1x,f10.4))
      END
***************************************************************
      SUBROUTINE rhs(fvec)

      implicit none
      include 'const.com'
      include 'outfld.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real fvec(ndim)
*     LOCAL
      integer j,mx,my,mj
      real eps_o,adv,adu
      real MW,TS,dn,sn,delta,delta_sst
 
      fvec = 0.0

      eps_o     = par(1)
      delta     = par(2)
      adu       = par(9)
      adv       = par(10)
      delta_sst = par(18)
*
*   equations for r : wave amplitude 
*
      j = 0
*
*     boundary conditions at x = 0
      mx = 0
      mj = j*(nx+1)+mx+1
      fvec(1) = b0
*
      DO mx = 1,nx
         mj = mx+1
         fvec(mj) = delta*r_t(mx,j) + eps_o*r(0,mx,j) 
     &        + r(1,mx,j) - atm(mx,j)/2.
      END DO
*
      j = 1
      DO mx = 0,nx
         mj = (nx+1)+mx+1
         fvec(mj) = r(0,mx,j)
      END DO
*     
      DO j = 2,ny
         dn = 2.*real(j) - 1.
         sn = sqrt(real(j)/real(j-1))
         DO mx = 0,nx-1
            mj = j*(nx+1)+mx+1
            fvec(mj) = dn*(delta*r_t(mx,j) + eps_o*r(0,mx,j))
     &           - r(1,mx,j) - real(j-1)*(atm(mx,j) - sn*atm(mx,j-2))/2.
         END DO
*     
*     boundary conditions at x = 1
         mx = nx
         mj = j*(nx+1)+mx+1
         fvec(mj) = b1(j)
      END DO
*
*     equations for T, coefficients q 
*     
      DO mx = 0,nx
         DO my = 0,ny
            MW = mheav(mx,my,0)
            Ts = tsub(mx,my,0)
            mj = my*(nx+1)+mx+1
            
            fvec(nbm + mj) = delta_sst*T_t(mx,my) 
     &           + eps_T(mx,my)*(T(0,mx,my)-T0(0,mx,my))
     &           + MW*(T(0,mx,my) - Ts )
     &           + adu*u_1(mx,my)*T(1,mx,my)
     &           + adv*v_1(mx,my)*T(2,mx,my)

         END DO
      END DO

      fvec = -fvec

 001  format(2(1x,I3),6(1x,f10.4))
      END
***************************************************************
      SUBROUTINE matB

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'tau.com'
      include 'outfld.com'
      include 'mat.com'
*     INPUT/OUTPUT
*     LOCAL
      integer mx,my,ix,iy,j,ij,mj
      real delta,alpha_w,eps_o,DMW,Ts,SST,delta_sst

      A = 0.0

      eps_o     = par(1)
      delta     = par(2)
      alpha_w   = par(13)
      delta_sst = par(18)

      DO mx=1,nx
       mj = mx + 1
       DO ix=0,nx
        ij = ix + 1
        A(mj,ij) = delta*q(ix,mx)
       END DO
      END DO
*
      DO j=2,ny
      DO mx=0,nx-1
       mj = j*(nx+1) + mx + 1
       DO ix=0,nx
        ij = j*(nx+1) + ix + 1
        A(mj,ij) = delta*(2.*real(j)-1.)*q(ix,mx)
       END DO
      END DO
      END DO
*
      DO mx=0,nx
      DO my=0,ny
       mj = my*(nx+1) + mx + 1
       DMW   = mheav(mx,my,1)
       Ts    = tsub(mx,my,0)
       SST   = T(0,mx,my)
       DO ix=0,nx
       DO iy=0,ny
        ij = iy*(nx+1) + ix + 1
        A(nbm + mj,ij) = -delta*fh(mx,my,ix,iy)*alpha_w*DMW*(SST-Ts)
        A(nbm + mj,nbm + ij) = delta_sst*f(mx,my,ix,iy)
       END DO
       END DO
      END DO
      END DO
      A = -A
*
      END
********************************************8
      SUBROUTINE matA

      implicit none
      include 'const.com'
      include 'tau.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'mat.com'
*     INPUT/OUTPUT
*     LOCAL
      real MW,DMW,DTs,Ts
      real SST,SSTx,SSTy,dn,sn
*     LOCAL
      integer j,mx,my,mj,ix,iy,ij,i0,tj
      real delta_s,mu,gamma_s,F_0,alpha_w,eps_o,eps_E,adv,adu
      real sw_E,sw_Q,k_Q,delta
      real duadT(0:nx,0:ny+1,0:nx,0:ny),dvadT(0:nx,0:ny+1,0:nx,0:ny)
      real du_AdT(0:nx,0:ny,0:nx,0:ny),dv_AdT(0:nx,0:ny,0:nx,0:ny)
      real dusdT(0:nx,0:ny,0:nx,0:ny),dvsdT(0:nx,0:ny,0:nx,0:ny)
      real dwsdT(0:nx,0:ny,0:nx,0:ny)
      real AdT(0:nx,0:ny+1,0:nx,0:ny),BdT(0:nx,0:ny+1,0:nx,0:ny)
      real airdT(0:nx,0:ny,0:nx,0:ny),birdT(0:nx,0:ny,0:nx,0:ny)
*
      eps_o     = par(1)
      delta     = par(2)
      mu        = par(3)
      delta_s   = par(5)
      gamma_s   = par(6)
      adu       = par(9)
      adv       = par(10)
      alpha_w   = par(13)
      k_Q       = par(25)
      sw_Q      = par(26)
*
      A = 0.0
      duadT = 0.0
      dvadT = 0.0
      du_AdT = 0.0
      dv_AdT = 0.0
      dusdT = 0.0
      dvsdT = 0.0
      dwsdT = 0.0
      IF (abs(sw_Q).GT.1.0e-04) THEN
       call dadT(T,duadT,dvadT,du_AdT,dv_AdT,dusdT,dvsdT,dwsdT)
      END IF
      AdT   = mu*((1.-sw_Q)*ca  + sw_Q*k_Q*duadT)
      BdT   = mu*((1.-sw_Q)*cb  + sw_Q*k_Q*dvadT)
      airdT = mu*((1.-sw_Q)*air + sw_Q*k_Q*du_AdT)
      birdT = mu*((1.-sw_Q)*bir + sw_Q*k_Q*dv_AdT)
*
*   equations for r : wave amplitude 
*
      j = 0
      mx = 0
      DO ix = 0,nx
       ij = ix+1
       A(1,ij) = q(ix,0)
      END DO

      DO ix = 0,nx
       DO j=2,nr
        ij = j*(nx+1)+ix+1
        A(1,ij) = - EB(j)*q(ix,0)/real(j-1)
       END DO
      END DO
*
      DO mx = 1,nx
       mj = mx+1
       DO ix = 0,nx
        i0 = ix+1
        A(mj,i0) = eps_o*q(ix,mx) + qx(ix,mx)
        DO iy=0,ny
         ij = iy*(nx+1)+ix+1
         A(mj,nbm + ij) = -AdT(mx,0,ix,iy)/2.
        END DO
       END DO
      END DO
*
      j = 1
      DO mx = 0,nx
       mj = (nx+1)+mx+1
       DO ix = 0,nx
        ij = (nx+1)+ix+1
        A(mj,ij) = q(ix,mx)
       END DO
      END DO
*
      DO j = 2,ny
       dn = 2.*real(j) - 1.
       sn = sqrt(real(j)/real(j-1))
       DO mx = 0,nx-1
        mj = j*(nx+1)+mx+1
        DO ix = 0,nx
         ij = j*(nx+1)+ix+1
         A(mj,ij) = dn*eps_o*q(ix,mx) - qx(ix,mx) 
         DO iy=0,ny
          tj = iy*(nx+1)+ix+1
          A(mj,nbm + tj) = - real(j-1)*(AdT(mx,j,ix,iy) 
     &                     - sn*AdT(mx,j-2,ix,iy))/2.
         END DO
        END DO
       END DO
      END DO
*
*     boundary conditions at x = 1
      mx = nx
      DO j=2,ny
       mj = j*(nx+1)+mx+1
       DO ix=0,nx
        i0 = ix + 1
        ij = j*(nx+1) + ix + 1
        A(mj,i0) = -EB(j)*q(ix,nx)
        A(mj,ij) = q(ix,nx)
       END DO
      END DO
*
*   equations for T
*
      DO mx = 0,nx
      DO my = 0,ny
       mj    = my*(nx+1)+mx + 1
       MW    = mheav(mx,my,0)
       DMW   = mheav(mx,my,1)
       Ts    = tsub(mx,my,0)
       DTs   = tsub(mx,my,1)
       SST   = T(0,mx,my)
       SSTx  = T(1,mx,my)
       SSTy  = T(2,mx,my)
       DO ix = 0,nx
       DO iy = 0,ny
        ij   = iy*(nx+1)+ix + 1

        A(nbm + mj,ij) = 
     &  + adu*SSTx*fu(mx,my,ix,iy) 
     &  - MW*DTs*fh(mx,my,ix,iy) 
     &  - alpha_w*eps_o*DMW*(SST-Ts)*fh(mx,my,ix,iy)

        A(nbm + mj,nbm + ij) =
     &  + (eps_T(mx,my) + MW)*f(mx,my,ix,iy)
     &  + (T(0,mx,my)-T0(0,mx,my))*deps_T(mx,my)*airdT(mx,my,ix,iy)
     &  + adu*u_1(mx,my)*fx(mx,my,ix,iy)
     &  + adv*v_1(mx,my)*fy(mx,my,ix,iy)
     &  + adu*mu*delta_s*gamma_s*SSTx*( (1.-sw_Q)*us(mx,my,ix,iy) 
     &  + sw_Q*k_Q*dusdT(mx,my,ix,iy) )
     &  + adv*mu*delta_s*gamma_s*SSTy*( (1.-sw_Q)*vs(mx,my,ix,iy) 
     &  + sw_Q*k_Q*dvsdT(mx,my,ix,iy) )
     &  + alpha_w*mu*delta_s*gamma_s*DMW*(SST - Ts)
     &    *( (1.-sw_Q)*ws(mx,my,ix,iy) 
     &  + sw_Q*k_Q*dwsdT(mx,my,ix,iy) )

       END DO
      END DO
      END DO
      END DO

001   format(2(1x,I3),6(1x,f10.4))
      END
**********************************************************
      SUBROUTINE meanf(z,r,u,vr,vq,h,T)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      real z(ndim),u(0:nx,0:ny),r(0:1,0:nx,0:ny)
      real vr(0:nx,0:ny),vq(0:nx,0:ny)
      real h(0:nx,0:ny)
      real T(0:2,0:nx,0:ny)
*     LOCAL 
      integer mx,my,ix,iy,ij,j

      r = 0.0
      DO j=0,ny
      DO mx=0,nx
       DO ix=0,nx
        ij = j*(nx+1) + ix + 1
        r(0,mx,j) = r(0,mx,j) + z(ij)*q(ix,mx)
        r(1,mx,j) = r(1,mx,j) + z(ij)*qx(ix,mx)
       END DO
      END DO
      END DO
    
      u  = 0.0
      vr = 0.0
      vq = 0.0
      h  = 0.0
      T  = 0.0
      DO mx=0,nx
      DO my=0,ny
       DO ix=0,nx
       DO iy=0,ny
        ij = iy*(nx+1) + ix + 1
        u(mx,my)   = u(mx,my)   + z(ij)*fu(mx,my,ix,iy)
        h(mx,my)   = h(mx,my)   + z(ij)*fh(mx,my,ix,iy)
        vr(mx,my)  = vr(mx,my)  + z(ij)*fv(mx,my,ix,iy)
        vq(mx,my)  = vq(mx,my)  + z(nbm+ij)*fvt(mx,my,ix,iy)
        T(0,mx,my) = T(0,mx,my) + z(nbm+ij)*f(mx,my,ix,iy)
        T(1,mx,my) = T(1,mx,my) + z(nbm+ij)*fx(mx,my,ix,iy)
        T(2,mx,my) = T(2,mx,my) + z(nbm+ij)*fy(mx,my,ix,iy)
       END DO
       END DO
      END DO
      END DO

      END
**********************************************************
      SUBROUTINE meanft(z,r,h,T)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      real z(ndim)
      real r(0:nx,0:ny)
      real h(0:nx,0:ny)
      real T(0:nx,0:ny)
*     LOCAL 
      integer mx,my,ix,iy,ij,j

      r = 0.0
      DO j=0,ny
         DO mx=0,nx
            DO ix=0,nx
               ij = j*(nx+1) + ix + 1
               r(mx,j) = r(mx,j) + z(ij)*q(ix,mx)
            END DO
         END DO
      END DO

      h = 0.0
      T = 0.0
      DO mx=0,nx
         DO my=0,ny
            DO ix=0,nx
               DO iy=0,ny
                  ij = iy*(nx+1) + ix + 1
                  h(mx,my) = h(mx,my) + z(ij)*fh(mx,my,ix,iy)
                  T(mx,my) = T(mx,my) + z(nbm+ij)*f(mx,my,ix,iy)
               END DO
            END DO
         END DO
      END DO
      
      END
**********************************************************
      SUBROUTINE sky(z,ua,va,u_A,v_A)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      real z(ndim),ua(0:nx,0:ny+1),va(0:nx,0:ny+1)
      real u_A(0:nx,0:ny),v_A(0:nx,0:ny)
*     LOCAL 
      integer mx,j,ix,iy,ij,my

      ua = 0.0
      va = 0.0
      DO mx=0,nx
      DO j=0,ny
       DO ix=0,nx
       DO iy=0,ny
        ij = iy*(nx+1) + ix + 1
        ua(mx,j) = ua(mx,j) + z(nbm+ij)*ca(mx,j,ix,iy)
        va(mx,j) = va(mx,j) + z(nbm+ij)*cb(mx,j,ix,iy)
       END DO
       END DO
      END DO
      END DO

      u_A = 0.0
      v_A = 0.0
      DO mx=0,nx
      DO my=0,ny
       DO ix=0,nx
       DO iy=0,ny
        ij = iy*(nx+1) + ix + 1
        u_A(mx,my) = u_A(mx,my) + z(nbm+ij)*air(mx,my,ix,iy)
        v_A(mx,my) = v_A(mx,my) + z(nbm+ij)*bir(mx,my,ix,iy)
       END DO
       END DO
      END DO
      END DO

      END
**********************************************************
      SUBROUTINE surface(z,u_s,v_s,w_s)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      real z(ndim),u_s(0:nx,0:ny)
      real v_s(0:nx,0:ny),w_s(0:nx,0:ny)
*     LOCAL 
      integer mx,my,ix,iy,ij

      u_s = 0.0
      v_s = 0.0
      w_s = 0.0
      DO mx=0,nx
      DO my=0,ny
       DO ix=0,nx
       DO iy=0,ny
        ij = iy*(nx+1) + ix + 1
        u_s(mx,my) = u_s(mx,my) + z(nbm+ij)*us(mx,my,ix,iy)
        v_s(mx,my) = v_s(mx,my) + z(nbm+ij)*vs(mx,my,ix,iy)
        w_s(mx,my) = w_s(mx,my) + z(nbm+ij)*ws(mx,my,ix,iy)
       END DO
       END DO
      END DO
      END DO

      END
**********************************************************
      SUBROUTINE OBC(z,b0,b1,bcT)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real z(ndim),b0,b1(0:ny),bcT(0:ny)
*     LOCAL 
      integer ix,j,ij,i0,my,iy
      real r00,r0n

*     boundary condition at x = 0
      r00 = 0.0
      DO ix=0,nx
       i0 = ix + 1
       r00 = r00 + q(ix,0)*z(i0)
      END DO
      r0n = 0.0
      DO ix=0,nx
       DO j= 2,nr
        ij = j*(nx+1) + ix + 1
        r0n = r0n + q(ix,0)*EB(j)*z(ij)/real(j-1)
       END DO
      END DO
      b0 = r00 - r0n

      bcT = 0.0
      DO my=0,ny
       DO ix=0,nx
       DO iy=0,ny
        ij = iy*(nx+1) + ix + 1
        bcT(my) = bcT(my) + f(0,my,ix,iy)*z(ij+nbm)
       END DO
       END DO
*      bcT(my) = bcT(my) - T0(0,0,my)
      END DO

*     boundary condition at x = 1
      b1 = 0.0
      DO j=2,ny
       DO ix=0,nx
        i0 = ix + 1
        ij = j*(nx+1) + ix + 1
        b1(j) = b1(j) + (z(ij)-EB(j)*z(i0))*q(ix,nx)
       END DO
      END DO

      END
**********************************************************
      SUBROUTINE dadT(SST,duadT,dvadT,du_AdT,dv_AdT,dusdT,dvsdT,dwsdT)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real SST(0:2,0:nx,0:ny),duadT(0:nx,0:ny+1,0:nx,0:ny)
      real dvadT(0:nx,0:ny+1,0:nx,0:ny)
      real du_AdT(0:nx,0:ny,0:nx,0:ny),dv_AdT(0:nx,0:ny,0:nx,0:ny)
      real dusdT(0:nx,0:ny,0:nx,0:ny),dvsdT(0:nx,0:ny,0:nx,0:ny)
      real dwsdT(0:nx,0:ny,0:nx,0:ny)
*     LOCAL 
      integer mx,my,ix,iy,jx,jy
      real Tc,T,tanh1,f0,f1,eps,dqdT(0:nx,0:ny,0:nx,0:ny)
      real dmhf(0:nx,0:ny),mtr

      Tc = par(21)
      eps = par(11)
      dqdT = 0.0
      DO mx=0,nx
      DO my=0,ny
*      Tc = par(17) - 2.*x(mx)
       T = SST(0,mx,my) - Tc
       mtr = 1. - (tanh(y(my) + 2.) - tanh(y(my)-2.))/2.
       tanh1 = tanh(T/eps)
       f0 = 1.+tanh1
       f1 = (1.-tanh1**2)/eps
       dmhf(mx,my) = mtr*(f0 + t*f1)
      END DO
      END DO
      dmhf = dmhf/2
      DO ix=0,nx
      DO iy=0,ny
      DO jx=0,nx
      DO jy=0,ny
       DO mx=0,nx
       DO my=0,ny
        dqdT(ix,iy,jx,jy) = dqdT(ix,iy,jx,jy) + f(mx,my,ix,iy)
     &                     *f_i(mx,my,jx,jy)*wh(my)*dmhf(mx,my)
       END DO
       END DO
      END DO
      END DO
      END DO
      END DO

      duadT = 0.0
      dvadT = 0.0
      du_AdT = 0.0
      dv_AdT = 0.0
      dusdT = 0.0
      dvsdT = 0.0
      dwsdT = 0.0
      DO mx=0,nx
      DO my=0,ny
      DO ix=0,nx
      DO iy=0,ny
       DO jx=0,nx
       DO jy=0,ny
        duadT(mx,my,ix,iy)  = duadT(mx,my,ix,iy)  
     &                       +  ca(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
        dvadT(mx,my,ix,iy)  = dvadT(mx,my,ix,iy)  
     &                       +  cb(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
        du_AdT(mx,my,ix,iy) = du_AdT(mx,my,ix,iy) 
     &                       + air(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
        dv_AdT(mx,my,ix,iy) = dv_AdT(mx,my,ix,iy) 
     &                       + bir(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
        dusdT(mx,my,ix,iy)  = dusdT(mx,my,ix,iy)  
     &                       +  us(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
        dvsdT(mx,my,ix,iy)  = dvsdT(mx,my,ix,iy)  
     &                       +  vs(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
        dwsdT(mx,my,ix,iy)  = dwsdT(mx,my,ix,iy)  
     &                       +  ws(mx,my,jx,jy)*dqdT(ix,iy,jx,jy)
       END DO
       END DO
      END DO
      END DO
      END DO
      END DO

001   format(4(1x,i4),2(1x,e13.5))      
      END
**********************************************************
      SUBROUTINE heatf(SST,qf)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real SST(0:2,0:nx,0:ny)
      real qnl(0:nx,0:ny),qf(ndim)
*     LOCAL 
      integer mx,my,ix,iy,ij
      real Tc,eps,tanh1,t,mtr

      Tc = par(21)
      eps = par(11)
      qnl = 0.0
      DO mx=0,nx
       DO my=0,ny
*       Tc = par(17) - 2.*x(mx)
        t = SST(0,mx,my)-Tc
        mtr = 1. - (tanh(y(my) + 2.) - tanh(y(my)-2.))/2.
        tanh1 = tanh(t/eps)
        qnl(mx,my) = mtr*t*(1.+tanh1)/2.
       END DO
      END DO
      qf = 0.0
      DO ix=0,nx
      DO iy=0,ny
       ij = iy*(nx+1) + ix + 1
       DO mx=0,nx
       DO my=0,ny
        qf(nbm + ij) = qf(nbm + ij) + wh(my)*f_i(mx,my,ix,iy)*qnl(mx,my)
       END DO
       END DO
      END DO
      END DO

001   format(3(1x,e13.5))
002   format(2(1x,i4),(1x,e13.5))
      END
**********************************************************
      SUBROUTINE wef(wind,epst,depst)

      implicit none
      include 'const.com'
      include 'tau.com'
*     INPUT/OUTPUT
      integer my,mx
      real wind(0:nx,0:ny),depst(0:nx,0:ny),epst(0:nx,0:ny)
      real eps,umin,ured1,ured2
      real tanh1,tanh2,ff1,ff2,df1,df2

      eps = par(11)
      umin = par(22)
      epst = 0.0
      depst = 0.0
      DO mx=0,nx
       DO my=0,ny
        ured1 = wind(mx,my) - umin
        tanh1 = tanh(ured1/eps)
        ff1 = 1.+tanh1
        df1 = (1.-tanh1**2)/eps
        ured2 = -wind(mx,my) - umin
        tanh2 = tanh(ured2/eps)
        ff2 = 1.+tanh2
        df2 = (1.-tanh2**2)/eps
        epst(mx,my)  = umin + (ured1*ff1/2.) + (ured2*ff2/2.)
        depst(mx,my) = ( (ff1+ured1*df1) - (ff2+ured2*df2) )/2.
       END DO
      END DO

      END
**********************************************************
      SUBROUTINE monitor(z,itp,xln,det,nfile)
      implicit none
      include 'const.com'
      include 'outfld.com'
      include 'tau.com'
*     INPUT/OUTPUT
      integer ibr,itp,nfile
      real    z(ndim),det,xln
*
      write(nfile,999) itp,xln,det
     &                    ,T(0,2,ny/2)  
     &                    ,T(0,nx/2,ny/2)  
     &                    ,T(0,nx-2,ny/2) 
     &                    ,T(0,nx-2,ny/2+2)
     &                    ,T(0,nx-2,ny/2-2)

 999  format(2x,i6,7e14.5)
 998  format(7e14.5)
      END
********************************************************
      SUBROUTINE anpnt(lab,ibr,ipnt,itp,icp,deto,sigo,xlo,xlpo)

      implicit none
      include 'const.com'
*     INPUT?OUTPUT
      integer lab,ibr,ipnt,itp,icp
      real deto,sigo,xlo,xlpo
*     COMMON
      real     uo,un,upo,w
      common /FLD/ uo(ndim),un(ndim),upo(ndim),w(ndim,nf)
*     FUNCTION
      real     g05saf,x
*     LOCAL
      integer  L,n
*
      lab=0
      ibr=1
      ipnt=1
      itp=9
      icp=0
      deto=0.0
      sigo=0.0
      xlpo=1.0

      uo=0.0
      DO n = 1,ndim
        DO L = 1,nf
*          w(n,L)=2*g05saf(x)-1.0
           w(n,L)=0.
         ENDDO
      ENDDO
      xlo=0.0

*999   format(4(1x,f13.5))
      end
