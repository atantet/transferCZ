*      #ifndef CONST_H
*      #define CONST_H

      integer nx,ny,na,ndim,nbm,nf,nr,maxstep
      real dstep
      REAL mu_0, epsn
*     PARAMETER (nx = 16,ny = 40,na = 7,nf = 7,nr = ny)
*     PARAMETER (nx = 4,ny = 4,na = ny,nf = 7,nr=ny)
*     PARAMETER (nx = 10,ny = 20,na = ny,nf = 7,nr=ny)
*     PARAMETER (nx = 14,ny = 30,na = ny,nf = 7,nr=ny)
*     PARAMETER (nx = 28,ny = 41,na = ny,nf = 7,nr=ny)
      PARAMETER (nx = 29,ny = 30,na = ny,nf = 7,nr=ny,maxstep=75000)
*      PARAMETER (dstep=0.35) 	 
*      PARAMETER (dstep=0.175) 	 
      PARAMETER (dstep=0.060) 	 
      PARAMETER (nbm = (nx+1)*(ny+1),ndim = 2*nbm)
      integer icp1,sw_A,sw_flux
      integer ioute,ioutnd,ioutd
      real(kind=8) F0_t(nx+1,ny+1,maxstep)
      real(kind=8) F0_t_dx(nx+1,ny+1,maxstep)
      real(kind=8) F0_t_dy(nx+1,ny+1,maxstep)
      real pi,par(31)
      real dpar(23)
      COMMON /CONST/ par,dpar,pi,icp1,sw_A,sw_flux
      COMMON /OUTP/ ioute,ioutnd,ioutd

*      #end // CONST_H
