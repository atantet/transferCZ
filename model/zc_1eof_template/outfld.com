      real u_A(0:nx,0:ny),T(0:2,0:nx,0:ny),h(0:nx,0:ny)
      real wind(0:nx,0:ny)
      real u_1(0:nx,0:ny),w_1(0:nx,0:ny),v_1(0:nx,0:ny)
      real w(0:nx,0:ny),w_s(0:nx,0:ny),h_t(0:nx,0:ny)
      real u(0:nx,0:ny),v(0:nx,0:ny),vrs(0:nx,0:ny),vq(0:nx,0:ny)
      COMMON /FIELD1/ u_1,v_1,w_1,h,T,u_A,v_A,w,w_s,h_t,u,v,wind
      real r(0:1,0:nx,0:ny),r_t(0:nx,0:ny)
      real T_t(0:nx,0:ny),eps_T(0:nx,0:ny),deps_T(0:nx,0:ny)
      real b0,b1(0:ny),bcT(0:ny)
      COMMON /FIELD2/ r,r_t,T_t,eps_T,deps_T,b0,b1,bcT
      real ua(0:nx,0:ny+1),va(0:nx,0:ny+1)
      real v_A(0:nx,0:ny),atm_A(0:nx,0:ny)
      real epst(0:nx,0:ny),depst(0:nx,0:ny),atm(0:nx,0:ny+1)
      COMMON /FIELD3/ ua,va,epst,depst,atm,atm_A
      real mheav(0:nx,0:ny,0:3),tsub(0:nx,0:ny,0:3)
      COMMON /PARMTZ/ mheav,tsub
