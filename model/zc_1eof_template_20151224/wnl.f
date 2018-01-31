**********************************************************************
      SUBROUTINE wnl(sigma,bs,z)
*     Not used!

      implicit none
      include 'const.com'
*     INPUT?OUTPUT
      complex sigma,z(ndim),dsdp,landau
      real bs(ndim),bst(ndim)
      integer n,itime
      
      itime = 0
      bst = 0.0
      landau = cmplx(0.,.0)
      dsdp = cmplx(0.,.0)
      call fields(bs,bst,itime,0.)
*     call singular(sigma,bs,z,dsdp)
      call svd(sigma)
      call linco(sigma,bs,z,dsdp)
      call nonl(sigma,bs,z,landau)

      write(39,*) '*************************************************************'
      write(39,*) 
      write(39,001) real(sigma),aimag(sigma),par(icp1)
     &             ,real(dsdp),aimag(dsdp),real(landau),aimag(landau)
      write(39,*) 
      write(39,*) ' dadt = ',dsdp,' a -',landau,' a|a|^2 '
      write(39,*) 
      write(39,*) ' amplitude  = ',sqrt(abs(real(dsdp)/real(landau)))
      write(39,*) 
      write(39,*) ' phase      = ',aimag(dsdp) - aimag(landau)*abs(real(dsdp)/real(landau))
      write(39,*) 
      write(39,*) '*************************************************************'
      
 001  format(7(1x,e13.5))
 
      END 
**********************************************************************
      SUBROUTINE singular(sigma,bs,z,dsdp)

      implicit none
      include 'const.com'
      include 'mat.com'
      include 'wnl.com'
*     INPUT?OUTPUT
      complex sigma,z(ndim),dsdp,zdp(ndim),zb(ndim)
      real bs(ndim)
*     SVD
      integer ifail
      logical wantp,wantq
      real S(ndim),wkr(5*ndim),dlta
      complex wkc(2*ndim+ndim**2),VH(ndim,ndim),Q,alpha(ndim,2)
*     VECTORS
      integer m
      complex norm

      dlta = 1.0e-04

      call matop1(bs,z,zb,zdp,sigma,dlta)

      DO m=1,ndim
       alpha(m,1) = zb(m)
       alpha(m,2) = zdp(m)
      END DO

      call matB
      L = A
      call matA
      L = A - sigma*L
      wantp = .false.
      wantq = .true.
      ifail = -1
      call f02xef(ndim,ndim,L,ndim,2,alpha,ndim,wantq,Q,ndim,S,
     &            wantp,VH,ndim,wkr,wkc,ifail)

      norm = alpha(ndim,1)

      DO m=1,ndim
	adjoint(m) = conjg(L(m,ndim))/norm
      END DO

      write(39,*) 'norm of svd = ',norm
      write(39,*) 'singular value = ',S(ndim),S(ndim-1)
      write(39,*) 'ifail svd = ',ifail

      dsdp = -alpha(ndim,2)/norm

      END
*****************************************************************
      SUBROUTINE linco(sigma,bs,z,dsdp)

      implicit none
      include 'const.com'
      include 'wnl.com'
*     INPUT?OUTPUT
      complex sigma,z(ndim),dsdp,zb(ndim),zdp(ndim)
      real bs(ndim),dlta
      integer m
      complex norm

      dlta = 1.0e-04

      call matop1(bs,z,zb,zdp,sigma,dlta)

      norm = cmplx(0.0,0.0)
      DO m=1,ndim
       norm    = norm + adjoint(m)*zb(m)
      END DO
      write(39,*) 'norm svd',norm
      adjoint = adjoint/norm
      DO m=1,ndim
       dsdp = dsdp + adjoint(m)*zdp(m)
      END DO
      dsdp = -dsdp

      END
*****************************************************************
      SUBROUTINE matop1(bs,z,zb,zdp,sigma,dlta)

      implicit none
      include 'const.com'
      include 'mat.com'
      include 'wnl.com'
*     INPUT?OUTPUT
      real bs(ndim),bst(ndim)
      complex sigma,zb(ndim),zdp(ndim),z(ndim)
      real dlta
*     LOCAL
      integer itime
      real p,pplus,di

      itime = 0
      bst   = 0.
      p = par(icp1)
      di = 1./dlta
      pplus  =  p + dlta

*     construction of B
      call matB
      zb = matmul(A,z)

*     construction of (A+sB)dp
      call matA
      L = A
      call matB
      L = L - sigma*A
      zdp = matmul(L,z)
      par(icp1) = pplus
      call fields(bs,bst,itime,0.)
      call matA
      L = A
      call matB
      L = L - sigma*A
      zdp = zdp - matmul(L,z)
      zdp = zdp*di

      par(icp1) = p
      
*     zb = Bu
*     dzdp = Adpu + sigma*Bdpu

      END
      
*****************************************************************
      SUBROUTINE nonl(sigma,bs,z,landau)

      implicit none
      include 'const.com'
      include 'mat.com'
      include 'wnl.com'
*     INPUT?OUTPUT
      complex sigma,z(ndim),landau
      real bs(ndim)
*     LOCAL
      integer m
      complex zstar(ndim),lambda(ndim)
      complex z02(ndim),z22(ndim),r02(ndim),r22(ndim)
      complex cnul

      cnul = cmplx(0.,.0)

      DO m=1,ndim
	zstar(m) = conjg(z(m))
      END DO

*     Nonlinear interactions at order(E**0 epsilon**2 )

      call N_2(sigma,z,conjg(sigma))
      r02 = matmul(L,zstar)
      z02 = -r02
*     call truncrhs(z02)

      call matA
      L = A
      call solvd(z02,0)

      Do m=1,ndim
        z02(m) = 2.*real(z02(m))
      END DO

*     Nonlinear interactions at order(E**2 epsilon**2 )

      call N_2(sigma,z,sigma)
      r22 = matmul(L,z)
      z22 = -r22
*     call truncrhs(z22)

      call matA
      L = A
      call matB
      L = L - 2.*sigma*A
      call solvd(z22,0)
 
*     Nonlinear interactions at order(E epsilon**3 )

*     Contribution of N2

      lambda = cnul
      call N_2(sigma,z,cnul)
      lambda = lambda + matmul(L,z02)
      call N_2(cnul,z02,sigma)
      lambda = lambda + matmul(L,z)
      call N_2(conjg(sigma),zstar,2.*sigma)
      lambda = lambda + matmul(L,z22)
      call N_2(2.*sigma,z22,conjg(sigma))
      lambda = lambda + matmul(L,zstar)

*     Contribution of N3

      call N_3(sigma,z,sigma,z,conjg(sigma))
      lambda = lambda + matmul(L,zstar)
      call N_3(sigma,z,conjg(sigma),zstar,sigma)
      lambda = lambda + matmul(L,z)
      call N_3(conjg(sigma),zstar,sigma,z,sigma)
      lambda = lambda + matmul(L,z)

*     call truncrhs(lambda)
      landau = cnul
      DO m=1,ndim
       write(60,*) m,lambda(m)
       landau = landau + lambda(m)*adjoint(m)
      END DO
       
      call monl(sigma,z,z02,r02,z22,r22,lambda)

      END
******************************************************************
      SUBROUTINE truncrhs(rhs)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT?OUTPUT
      complex rhs(ndim)
*     LOCAL
      complex cnul
      real trf,trl
      integer mx,my,mj
*     
      trl = 2.
      DO my = 0,ny
       trf = (tanh((y(my)+trl)*100.) - tanh((y(my)-trl)*100.))/2.
       DO mx = 0,nx
        mj = my*(nx+1) + mx + 1
        rhs(nbm+mj) = trf*rhs(nbm+mj)
       END DO
      END DO

      END
******************************************************************
      SUBROUTINE monl(sigma,z,z02,r02,z22,r22,lambda)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'wnl.com'
*     INPUT?OUTPUT
      complex lambda(ndim),z(ndim),sigma
      complex z02(ndim),z22(ndim),r02(ndim),r22(ndim)
*     LOCAL
      complex cnul
      integer mx,my,mj

      cnul = cmplx(0.,.0)
*     print z22
      call outper(z22,2.*sigma,52)
*     print z02
      call outper(z02,cnul,50)
*     
      write(54,*) nx+1,ny+1
      DO mx = 0,nx
      DO my = 0,ny
       mj = my*(nx+1) + mx + 1
       write(54,002) x(mx),y(my),real(r02(mj+nbm)),real(r22(mj+nbm)),aimag(r22(mj+nbm))
      END DO
      END DO
*
      write(55,*) nx+1,ny+1
      DO mx = 0,nx
      DO my = 0,ny
       mj = my*(nx+1) + mx + 1
       write(55,002) x(mx),y(my)
     & ,real(lambda(mj+nbm)),aimag(lambda(mj+nbm))
      END DO
      END DO
*
      write(56,*) nx+1,ny+1
      DO mx = 0,nx
      DO my = 0,ny
       mj = my*(nx+1) + mx + 1
       write(56,002) x(mx),y(my)
     & ,real(adjoint(mj+nbm)),aimag(adjoint(mj+nbm))
      END DO
      END DO
*
      write(58,*) nx+1,ny+1
      DO mx = 0,nx
      DO my = 0,ny
       write(58,002) x(mx),y(my),mheav(mx,my,0),mheav(mx,my,1),mheav(mx,my,2),mheav(mx,my,3)
       write(58,002) x(mx),y(my),tsub(mx,my,0),tsub(mx,my,1),tsub(mx,my,2),tsub(mx,my,3)
      END DO
      END DO
*
002   format(9(1x,e13.5))
      END
******************************************************************
      SUBROUTINE svd(sigma)

      implicit none
      include 'const.com'
      include 'mat.com'
      include 'wnl.com'
*     INPUT?OUTPUT
      complex sigma
*     LOCAL
      integer m,ni
      real norm,g05caf,x
      complex s(ndim)

      ni = 1
      DO m=1,ndim
       s(m) = g05caf(x)-.5
      END DO

      call matA
      L = A
      call matB
      L = L - sigma*A
      L = matmul(L,transpose(conjg(L)))

      call solvd(s,0)
      DO m=1,ni
       norm = dot_product(s,conjg(s))
       write(39,*) m,norm
       s = s/norm
       call solvd(s,1)
      END DO
  
      norm = dot_product(s,conjg(s))
      write(39,*) m,norm
      s = s/norm

      adjoint = conjg(s)

      call matA
      L = A
      call matB
      L = L - sigma*A
      s = matmul(adjoint,L)
      write(39,*) 'norm residu = ',dot_product(s,conjg(s))
      write(39,*) '*********************************'
*
      END
*************************************************************
      SUBROUTINE N_2(sigma1,z1,sigma2)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'tau.com'
      include 'wnl.com'
*     INPUT/OUTPUT
      complex sigma1,sigma2,z1(ndim)
*     LOCAL
      integer j,mx,my,mj,ix,iy,ij
      real MW,DMW,DDMW,DDTs,DTs,Ts,SST
      real delta_s,mu,gamma_s,alpha_w,eps_o,adv,adu,delta
      real zr(ndim),zi(ndim)
      complex u1(0:nx,0:ny),v1(0:nx,0:ny),w1(0:nx,0:ny),h1(0:nx,0:ny)
      real hr(0:nx,0:ny),hi(0:nx,0:ny),Tr(0:nx,0:ny,0:2),Ti(0:nx,0:ny,0:2)
      real ur(0:nx,0:ny),ui(0:nx,0:ny),vr(0:nx,0:ny),vi(0:nx,0:ny)
      real vqr(0:1,0:nx,0:ny),vqi(0:1,0:nx,0:ny)
      real wr(0:1,0:nx,0:ny),wi(0:1,0:nx,0:ny)
      real rr(0:1,0:nx,0:ny),ri(0:1,0:nx,0:ny)
      real u_sr(0:nx,0:ny),v_sr(0:nx,0:ny),w_sr(0:nx,0:ny)
      real u_si(0:nx,0:ny),v_si(0:nx,0:ny),w_si(0:nx,0:ny)
*
      L = cmplx(0.,.0)

      eps_o   = par(1)
      delta   = par(2)
      mu      = par(3)
      delta_s = par(5)
      gamma_s = par(6)
      adu     = par(9)
      adv     = par(10)
      alpha_w = par(13)

      DO j=1,ndim
       zr(j) = real(z1(j))
       zi(j) = aimag(z1(j))
      END DO

      call meanf(zr,rr,ur,vr,vqr,hr,Tr)
      call meanf(zi,ri,ui,vi,vqi,hi,Ti)
      call surface(zr,u_sr,v_sr,w_sr)
      call surface(zi,u_si,v_si,w_si)
*
      h1 = cmplx(hr,hi)
      u1 = gamma_s*delta_s*mu*cmplx(u_sr,u_si)
      v1 = gamma_s*delta_s*mu*cmplx(v_sr,v_si)
      w1 = -(delta*sigma1 + eps_o)*h1 + gamma_s*delta_s*mu*cmplx(w_sr,w_si)
      w1 = alpha_w*w1
*
*   equations for r : wave amplitude : no nonlinear interactions
*
*
*   equations for T
*
      DO mx = 0,nx
      DO my = 0,ny
       mj    = my*(nx+1)+mx + 1
       MW    = mheav(mx,my,0)
       DMW   = mheav(mx,my,1)*alpha_w
       DDMW  = mheav(mx,my,2)*alpha_w*alpha_w*.5
       Ts    = tsub(mx,my,0)
       DTs   = tsub(mx,my,1)
       DDTs  = tsub(mx,my,2)*.5
       SST   = T(0,mx,my)+T0(0,mx,my)
       DO ix = 0,nx
       DO iy = 0,ny
        ij   = iy*(nx+1)+ix + 1

        L(nbm + mj,ij) = 
     &  - MW*DDTs*h1(mx,my)*fh(mx,my,ix,iy)
     &  - DMW*DTs*w1(mx,my)*fh(mx,my,ix,iy)
     &  - DDMW*(SST - Ts)*w1(mx,my)*(delta*sigma2+eps_o)*fh(mx,my,ix,iy)

        L(nbm + mj,nbm + ij) = 
     &  + adu*u1(mx,my)*fx(mx,my,ix,iy) 
     &  + adu*v1(mx,my)*fy(mx,my,ix,iy) 
     &  + DMW*w1(mx,my)*f(mx,my,ix,iy)
     &  + mu*delta_s*gamma_s*DDMW*(SST - Ts)*w1(mx,my)*ws(mx,my,ix,iy) 

       END DO
      END DO
      END DO
      END DO

      END
**********************************************************
      SUBROUTINE N_3(sigma1,z1,sigma2,z2,sigma3)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'tau.com'
      include 'wnl.com'
*     INPUT/OUTPUT
      complex sigma1,sigma2,sigma3,z1(ndim),z2(ndim)
*     LOCAL
      integer j,mx,my,mj,ix,iy,ij
      real MW,DMW,DDMW,DDDMW,DDDTs,DDTs,DTs,Ts,SST
      real delta_s,mu,gamma_s,alpha_w,eps_o,delta
      real zr(ndim),zi(ndim)
      real hr(0:nx,0:ny),hi(0:nx,0:ny),Tr(0:nx,0:ny,0:2),Ti(0:nx,0:ny,0:2)
      real ur(0:nx,0:ny),ui(0:nx,0:ny),vr(0:nx,0:ny),vi(0:nx,0:ny)
      real rr(0:1,0:nx,0:ny),ri(0:1,0:nx,0:ny)
      real vqr(0:nx,0:ny),vqi(0:1,0:nx,0:ny)
      real wr(0:nx,0:ny),wi(0:nx,0:ny)
      complex w1(0:nx,0:ny),h1(0:nx,0:ny),w2(0:nx,0:ny),h2(0:nx,0:ny)
      real u_sr(0:nx,0:ny),v_sr(0:nx,0:ny),w_sr(0:nx,0:ny)
      real u_si(0:nx,0:ny),v_si(0:nx,0:ny),w_si(0:nx,0:ny)
*
      L = cmplx(0.,.0)

      eps_o   = par(1)
      delta   = par(2)
      mu      = par(3)
      delta_s = par(5)
      gamma_s = par(6)
      alpha_w = par(13)

      DO j=1,ndim
       zr(j) = real(z1(j))
       zi(j) = aimag(z1(j))
      END DO

      call meanf(zr,rr,ur,vr,vqr,hr,Tr)
      call meanf(zi,ri,ui,vi,vqi,hi,Ti)
      call surface(zr,u_sr,v_sr,w_sr)
      call surface(zi,u_si,v_si,w_si)
*
      h1 = cmplx(hr,hi)
      w1 = -(delta*sigma1 + eps_o)*h1 + gamma_s*delta_s*mu*cmplx(w_sr,w_si)
      w1 = alpha_w*w1
*
      DO j=1,ndim
       zr(j) = real(z2(j))
       zi(j) = aimag(z2(j))
      END DO

      call meanf(zr,rr,ur,vr,vqr,hr,Tr)
      call meanf(zi,ri,ui,vi,vqi,hi,Ti)
      call surface(zr,u_sr,v_sr,w_sr)
      call surface(zi,u_si,v_si,w_si)
*
      h2 = cmplx(hr,hi)
      w2 = -(delta*sigma2 + eps_o)*h2 + gamma_s*delta_s*mu*cmplx(w_sr,w_si)
      w2 = alpha_w*w2
*
*   equations for r : wave amplitude : no nonlinear interactions
*
*
*   equations for T
*
      DO mx = 0,nx
      DO my = 0,ny
       mj    = my*(nx+1)+mx + 1
       MW    = mheav(mx,my,0)
       DMW   = mheav(mx,my,1)*alpha_w
       DDMW  = mheav(mx,my,2)*alpha_w*alpha_w/2.
       DDDMW = mheav(mx,my,3)*alpha_w*alpha_w*alpha_w/6.
       Ts    = tsub(mx,my,0)
       DTs   = tsub(mx,my,1)
       DDTs  = tsub(mx,my,2)/2.
       DDDTs = tsub(mx,my,3)/6.
       SST   = T(0,mx,my)+T0(0,mx,my)
       DO ix = 0,nx
       DO iy = 0,ny
        ij   = iy*(nx+1)+ix + 1

        L(nbm + mj,ij) = 
     &  - MW*DDDTs*h1(mx,my)*h2(mx,my)*fh(mx,my,ix,iy)
     &  - DMW*DDTs*w1(mx,my)*h2(mx,my)*fh(mx,my,ix,iy)
     &  - DDMW*DTs*w1(mx,my)*w2(mx,my)*fh(mx,my,ix,iy)
     &  - DDDMW*(SST - Ts)*w1(mx,my)*w2(mx,my)*(delta*sigma3+eps_o)*fh(mx,my,ix,iy)

        L(nbm + mj,nbm + ij) = 
     &  + DDMW*w1(mx,my)*w2(mx,my)*f(mx,my,ix,iy)
     &  + mu*delta_s*gamma_s*DDDMW*(SST - Ts)*w1(mx,my)*w2(mx,my)*ws(mx,my,ix,iy) 

       END DO
      END DO
      END DO
      END DO

001   format(2(1x,I3),6(1x,f10.4))
      END
**********************************************************
