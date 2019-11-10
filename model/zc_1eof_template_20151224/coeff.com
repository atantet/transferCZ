      real x(0:nx),y(0:ny),wh(0:ny)
      real fu(0:nx,0:ny,0:nx,0:ny)
      real fv(0:nx,0:ny,0:nx,0:ny)
      real fvt(0:nx,0:ny,0:nx,0:ny)
      real fh(0:nx,0:ny,0:nx,0:ny)
      real f(0:nx,0:ny,0:nx,0:ny)
      real fx(0:nx,0:ny,0:nx,0:ny)
      real fy(0:nx,0:ny,0:nx,0:ny)
      real us(0:nx,0:ny,0:nx,0:ny)
      real vs(0:nx,0:ny,0:nx,0:ny)
      real ws(0:nx,0:ny,0:nx,0:ny)
      real ca(0:nx,0:ny+1,0:nx,0:ny)
      real cb(0:nx,0:ny+1,0:nx,0:ny)
      real psi(0:ny,0:ny),psi_y(0:ny,0:ny),EB(0:ny)
      real q(0:nx,0:nx),qx(0:nx,0:nx)
      real q_i(0:nx,0:nx),psi_i(0:ny,0:ny)
      COMMON /XY/ x,y,wh
      COMMON /BOUNDC/ psi,q,qx,EB,q_i,psi_i
      COMMON /COEFF/ f,fx,fy,fu,fv,fvt,fh
      COMMON /COEFS/ ca,cb,us,vs,ws
      real f_i(0:nx,0:ny,0:nx,0:ny)
      real air(0:nx,0:ny,0:nx,0:ny)
      real airx(0:nx,0:ny,0:nx,0:ny)
      real airy(0:nx,0:ny,0:nx,0:ny)
      real bir(0:nx,0:ny,0:nx,0:ny)
      real birx(0:nx,0:ny,0:nx,0:ny)
      real biry(0:nx,0:ny,0:nx,0:ny)
      COMMON /COEFP/ f_i,psi_y
      COMMON /COEFA/ air,airx,airy,bir,birx,biry
