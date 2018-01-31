      SUBROUTINE outclim(z0)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real z0(ndim)
*     LOCAL
      integer ij,mx,my,nev,n
*
      if (ioute.eq.1) then 
         DO mx=0,nx
            DO my=0,ny
               write(30,001) T(0,mx,my),atm(mx,my)
     &              ,u_1(mx,my) - u(mx,my),v_1(mx,my) - v(mx,my)
     &              ,w_1(mx,my) - w(mx,my)
            END DO
         END DO

         write(32,*) nx+1,ny+1
         DO mx=0,nx
            DO my=0,ny
               ij = my*(nx+1) + mx + 1
               write(32,002) real(mx),real(my)
     &              ,abs(z0(ij)),abs(z0(nbm+ij))
            END DO
         END DO
      endif                     ! ioute 
*
! non-dimensional output
*
      if (ioutnd.eq.1) then 
!     write steady state fields
         write(33,*) nx+1,ny+1
         DO mx=0,nx
            DO my=0,ny
               write(33,002) x(mx),y(my)
     &              ,u_1(mx,my),v_1(mx,my),w_1(mx,my)
     &              ,h(mx,my),T(0,mx,my)-T0(0,mx,my),u_A(mx,my)
            END DO
         END DO
!     write steady equatorial profiles 
         my = ny/2
         DO mx=0,nx
            write(34,003) x(mx),w_1(mx,my),h(mx,my),T(0,mx,my),
     &           u_A(mx,my)
         END DO
      endif                     ! ioutnd
      
*D*
!     dimensional output
*     D*
      if (ioutd.eq.1) then 
         call dimen(u_1,v_1,w_1,h,T,u_A)
!     write steady state fields
         write(133,*) nx+1,ny+1
         DO mx=0,nx
            DO my=0,ny
               write(133,002) x(mx)*dpar(1),y(my)*dpar(10)
     &              ,u_1(mx,my),v_1(mx,my),w_1(mx,my)
     &              ,h(mx,my),T(0,mx,my)-T0(0,mx,my),u_A(mx,my)
            END DO
         END DO
!     write steady equatorial profiles 
         my = ny/2
         DO mx=0,nx
            write(134,003) x(mx)*dpar(1),w_1(mx,my),h(mx,my),
     &           T(0,mx,my),u_A(mx,my)
         END DO
      endif                     ! ioutd
      
001   format(7(1x,e13.5))
002   format(8(1x,e13.5))
003   format(8(1x,f13.5))
      END
***************************************************
      SUBROUTINE outper(z,sigma,nfile)

      implicit none
      include 'const.com'
      include 'coeff.com'
*     INPUT/OUTPUT
      integer nfile
      complex z(ndim,nf),sigma(nf)
*     LOCAL
      integer ij,mx,my,nev,n
      real delta_s,mu,gamma_s,F_0,alpha_w,eps_o,adv,adu
      real z_r(ndim),z_i(ndim),re,im
      real ur(0:nx,0:ny),vr(0:nx,0:ny),wr(0:nx,0:ny),hr(0:nx,0:ny)
      real Tr(0:2,0:nx,0:ny)
      real ui(0:nx,0:ny),vi(0:nx,0:ny),wi(0:nx,0:ny),hi(0:nx,0:ny)
      real Ti(0:2,0:nx,0:ny)
      real vqr(0:nx,0:ny),vqi(0:nx,0:ny)
      real rr(0:1,0:nx,0:ny),ri(0:1,0:nx,0:ny)
      real uar(0:nx,0:ny+1),var(0:nx,0:ny+1),u_Ar(0:nx,0:ny)
      real v_Ar(0:nx,0:ny),u_Ai(0:nx,0:ny)
      real w_sr(0:nx,0:ny),u_sr(0:nx,0:ny),v_sr(0:nx,0:ny)
      real w_si(0:nx,0:ny),u_si(0:nx,0:ny),v_si(0:nx,0:ny)
 
*     select eigenvector
      nev = 1

      re = real(sigma(nev))
      im = aimag(sigma(nev))
      write(31,*) nev,re,im
      eps_o   = par(1)
      mu      = par(3)
      delta_s = par(5)
      gamma_s = par(6)
      F_0     = par(4)
      adu     = par(9)
      adv     = par(10)
      alpha_w = par(13)
*
      DO n=1,ndim
       z_r(n) = real(z(n,nev))
       z_i(n) = aimag(z(n,nev))
      END DO
      call sky(z_r,uar,var,u_Ar,v_Ar)
      call sky(z_i,uar,var,u_Ai,v_Ar)
      call meanf(z_r,rr,ur,vr,vqr,hr,Tr)
      call meanf(z_i,ri,ui,vi,vqi,hi,Ti)
      call surface(z_r,u_sr,v_sr,w_sr)
      call surface(z_i,u_si,v_si,w_si)
     
      wr =  -(re + eps_o)*hr + im*hi
      wi =  -(re + eps_o)*hi - re*hr
      ur =  ur + adu*gamma_s*delta_s*mu*u_sr
      vr =  vr - mu*vqr + adv*gamma_s*delta_s*mu*v_sr
      wr =  wr + gamma_s*delta_s*mu*w_sr
      ui =  ui + adu*gamma_s*delta_s*mu*u_sr
      vi =  vi - mu*vqi + adv*gamma_s*delta_s*mu*v_sr
      wi =  wi + gamma_s*delta_s*mu*w_sr
      wi =  alpha_w*wi
      wr =  alpha_w*wr

      call dimen(ur,vr,wr,hr,Tr,u_Ar)
      call dimen(ui,vi,wi,hi,Ti,u_Ai)

      write(nfile,*) nx+1,ny+1
      write(nfile+1,*) nx+1,ny+1
      write(nfile+2,*) nx+1,ny+1
      DO mx=0,nx
      DO my=0,ny
       write(nfile,002) x(mx),y(my)
     &    ,ur(mx,my),vr(mx,my),wr(mx,my)
     &    ,hr(mx,my),Tr(0,mx,my),u_Ar(mx,my)
       write(nfile+1,002) x(mx),y(my)
     &    ,ui(mx,my),vi(mx,my),wi(mx,my)
     &    ,hi(mx,my),Ti(0,mx,my),u_Ai(mx,my)
       write(nfile+2,002) x(mx),real(my),rr(0,mx,my),ri(0,mx,my)
       END DO
      END DO

*
002   format(8(1x,e13.5))
003   format(8(1x,f13.5))
      END
***************************************************
      SUBROUTINE outhov(z,time)
      ! Called during time integration

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'tau.com'
*     INPUT/OUTPUT
      real time,z(ndim)
*     LOCAL
      integer mx,my,n,ij,ix,iy
	  
*
! 1) non-dimensional output
*
      if (ioutnd.eq.1) then
!     write  equatorial fields 
         my = ny/2
         DO mx=0,nx
            write(41,002) x(mx),time,
     &           w_1(mx,my),h(mx,my),T(0,mx,my)-T0(0,mx,my),
     &           u_A(mx,my),wind(mx,my)
            
            if (ioute.eq.1) then 
               write(42,002) x(mx),time,
     &              w_1(mx,my),w(mx,my),w_s(mx,my),h_t(mx,my),h(mx,my)
               write(43,002) x(mx),time,
     &              T(0,mx,my),tsub(mx,my,0),tsub(mx,my,1),
     &              T(0,mx,my)-tsub(mx,my,0),mheav(mx,my,0),
     &              mheav(mx,my,1)
               write(44,002) x(mx),time,
     &              T_t(mx,my),
     &              mheav(mx,my,0)*(T(0,mx,my)-tsub(mx,my,0)),
     &              T(0,mx,my)-T0(0,mx,my),v_1(mx,my)*T(2,mx,my)
            endif               ! ioute
         END DO
         
         if (ioute.eq.1) then
!     write eastern pacific sections 
            mx = nx-2
            DO my=0,ny
               write(45,002) y(my),time,
     &              w_1(mx,my),h(mx,my),T(0,mx,my)-T0(0,mx,my),
     &              u_A(mx,my),wind(mx,my)
               write(46,002) y(my),time,
     &              w_1(mx,my),w(mx,my),w_s(mx,my),h_t(mx,my),h(mx,my)
               write(47,002) y(my),time,
     &              T(0,mx,my),tsub(mx,my,0),tsub(mx,my,1),
     &              T(0,mx,my)-tsub(mx,my,0),
     &              mheav(mx,my,0),mheav(mx,my,1)
               write(48,002) y(my),time, T_t(mx,my),
     &              mheav(mx,my,0)*(T(0,mx,my)-tsub(mx,my,0)),
     &              T(0,mx,my)-T0(0,mx,my),
     &              v_1(mx,my)*T(2,mx,my)
            END DO
         endif                  ! ioute
*     
!     write complete fields
!     total wind stress = wind,
!     external wind stress = taux, 
!     wind stress due to coupling = u_A 
!     radiative-equilibrium temperature T0
!     Temperature of ocean T
!     Thermocline depth h
         DO mx=0,nx
            DO my=0,ny
!               write(49,002) time,x(mx),y(my),h(mx,my),T(0,mx,my),
!     &              u_A(mx,my),T0(0,mx,my),wind(mx,my),taux(mx,my,0)
               write(49,002) h(mx,my),T(0,mx,my),u_A(mx,my),
     &              taux(mx,my,0)
               ij = my*(nx+1) + mx + 1
!               write(50,002) time,real(mx),real(my),z(ij),z(nbm+ij)
            END DO
         END DO     
* 
         
!     write equatorial time series eastern and western pacific 
         my=ny/2
*     mx=3*nx/4
         mx=nx-3                ! only for check with Malika's results
         write(51,002) time,T(0,4,my),T(0,mx,my),h(4,my),
     &        h(mx,my)
      endif                     ! ioutnd
      
*     D*
! 2) dimensional output
*D*
      if (ioutd.eq.1) then
         call dimen(u,v,w,h,T,u_A)
         call dimen2(time,taux,wind)
         my = ny/2
         DO mx=0,nx
            write(141,002) x(mx)*dpar(1),time
     &           ,w_1(mx,my),h(mx,my),T(0,mx,my)-T0(0,mx,my),
     &           u_A(mx,my),wind(mx,my)
         END DO
         
!     write complete fields
!     total wind stress = wind,
!     external wind stress = taux, 
!     wind stress due to coupling = u_A 
         DO mx=0,nx
            DO my=0,ny
               write(149,002) time,x(mx)*dpar(1),y(my)*dpar(10),
     &              h(mx,my),T(0,mx,my),u_A(mx,my),
     &              T0(0,mx,my),wind(mx,my),taux(mx,my,0)
               ij = my*(nx+1) + mx + 1
               write(150,002) time,real(mx),real(my),z(ij),z(nbm+ij)
            END DO
         END DO     
*     
!     write equatorial time series eastern and western pacific 
         my=ny/2
*     mx=3*nx/4
         mx=nx-3                ! only for check with Malika's results
         write(151,002) time,T(0,4,my),T(0,mx,my),h(4,my),h(mx,my)
         
      endif                     !ioutd
*     
 002  format(10(1x,e13.5))
003   format(8(1x,f10.5))
      

      END
************************************************************
      SUBROUTINE dimen(u,v,w,h,T,u_A)

      implicit none
      include 'const.com'
*     INPUT/OUTPUT
      real Hbar,H1,Astar,c0,ca,Lxdim,Lydim,u0,v0,w0
      real u(0:nx,0:ny),v(0:nx,0:ny),w(0:nx,0:ny),h(0:nx,0:ny)
      real T(0:2,0:nx,0:ny)
      real u_A(0:nx,0:ny)

      c0 = dpar(5) 
      Lxdim  = dpar(1)
      Lydim = dpar(10)
      Hbar = dpar(2)
      H1 = dpar(3)
      Astar = 1.0
      ca  = dpar(6)
	  
      v0 = c0*Lydim/Lxdim
      u0 = c0
      w0 = c0*H1/Lxdim
*     conversion to m/day
      w0 = w0*24.*3600.
      u = u0*u
      v = v0*v
      w = w0*w
      h = Hbar*h
      u_A = Astar*u_A*ca

      END
************************************************************
      SUBROUTINE dimen2(time,taux,wind)

      implicit none
      include 'const.com'
*     INPUT/OUTPUT
      real timedim,tau0,ca,mu
      real time,taux(0:nx,0:ny,0:2),wind(0:nx,0:ny)
      
      timedim = dpar(1)/dpar(5)    ! L/c_o in seconds
	  tau0    = dpar(7)
	  ca   = dpar(6)
	  mu      = par(3)
	  
	  time = time*timedim
	  taux = taux*tau0
	  wind = wind*tau0	  
      END
***************************************************
      SUBROUTINE analysis(z,sigma,nfile)

      implicit none
      include 'const.com'
      include 'coeff.com'
      include 'outfld.com'
      include 'tau.com'
*     INPUT/OUTPUT
      integer nfile
      complex z(ndim,nf),sigma(nf)
*     LOCAL
      integer ij,mx,my,nev,n
      real delta_s,mu,gamma_s,F_0,alpha_w,eps_o,adv,adu
      real z_r(ndim),z_i(ndim),re,im
      real ur(0:nx,0:ny),vr(0:nx,0:ny),wr(0:nx,0:ny),hr(0:nx,0:ny)
      real Tr(0:2,0:nx,0:ny)
      real ui(0:nx,0:ny),vi(0:nx,0:ny),wi(0:nx,0:ny),hi(0:nx,0:ny)
      real Ti(0:2,0:nx,0:ny)
      real vqr(0:nx,0:ny),vqi(0:nx,0:ny)
      real rr(0:1,0:nx,0:ny),ri(0:1,0:nx,0:ny)
      real uar(0:nx,0:ny+1),var(0:nx,0:ny+1),u_Ar(0:nx,0:ny)
      real v_Ar(0:nx,0:ny),u_Ai(0:nx,0:ny)
      real w_sr(0:nx,0:ny),u_sr(0:nx,0:ny),v_sr(0:nx,0:ny)
      real w_si(0:nx,0:ny),u_si(0:nx,0:ny),v_si(0:nx,0:ny)
 
*     select eigenvector
*     nev = 58
      nev = 1

      re = real(sigma(nev))
      im = aimag(sigma(nev))
      eps_o   = par(1)
      mu      = par(3)
      delta_s = par(5)
      gamma_s = par(6)
      F_0     = par(4)
      adu     = par(9)
      adv     = par(10)
      alpha_w = par(13)
*
      DO n=1,ndim
       z_r(n) = real(z(n,nev))
       z_i(n) = aimag(z(n,nev))
      END DO
      call meanf(z_r,rr,ur,vr,vqr,hr,Tr)
      call meanf(z_i,ri,ui,vi,vqi,hi,Ti)
      call surface(z_r,u_sr,v_sr,w_sr)
      call surface(z_i,u_si,v_si,w_si)
     
      wr =  -(re + eps_o)*hr + im*hi
      wi =  -(re + eps_o)*hi - im*hr
      vr =  gamma_s*delta_s*mu*v_sr
      wr =  gamma_s*delta_s*mu*w_sr + wr
      vi =  gamma_s*delta_s*mu*v_si
      wi =  gamma_s*delta_s*mu*w_si + wi

      write(nfile,*) nx+1,ny+1
      write(nfile,*) re,im,par(12)
      DO mx=0,nx
      DO my=0,ny
       write(nfile,002) x(mx),y(my)
     &    ,alpha_w*mheav(mx,my,1)*wr(mx,my),vr(mx,my),Tr(0,mx,my),
     &     Tr(2,mx,my),Tr(0,mx,my) - tsub(mx,my,1)*hr(mx,my)
     &    ,alpha_w*mheav(mx,my,1)*wi(mx,my),vi(mx,my),Ti(0,mx,my),
     &     Ti(2,mx,my),Ti(0,mx,my) - tsub(mx,my,1)*hi(mx,my)
     &    ,alpha_w*mheav(mx,my,0),v_1(mx,my),T(2,mx,my),
     &     T(0,mx,my) - tsub(mx,my,0)
       END DO
      END DO
*
      write(nfile+1,*) nx+1,ny+1
      write(nfile+1,*) par(12)
      DO mx=0,nx
      DO my=0,ny
       write(nfile+1,002) x(mx),y(my)
     &    ,alpha_w*mheav(mx,my,0),v_1(mx,my),T(2,mx,my),
     &    T(0,mx,my) - tsub(mx,my,0),T0(0,mx,my),T(0,mx,my)
       END DO
      END DO
002   format(16(1x,e13.5))
      END
***************************************************
