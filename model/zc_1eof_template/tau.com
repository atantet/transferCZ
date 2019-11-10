      real tau_x(0:nx,0:ny+1),taux(0:nx,0:ny,0:2),tauy(0:nx,0:ny,0:2)
      real u_s_ext(0:nx,0:ny),v_s_ext(0:nx,0:ny),w_s_ext(0:nx,0:ny)
      real vtau(0:nx,0:ny)
      real TRE(0:2,0:nx,0:ny)
      real T0(0:2,0:nx,0:ny),T0_t(0:nx,0:ny)
      real zf(ndim)
      real Tbar(0:nx,0:ny)
      real tau_obs(0:nx,0:ny),u_s_obs(0:nx,0:ny)
      real v_s_obs(0:nx,0:ny),w_s_obs(0:nx,0:ny)
      COMMON /FORCING/ tau_x,T0,T0_t,TRE,taux,tauy,vtau
      COMMON /SURFF/ u_s_ext,v_s_ext,w_s_ext
      COMMON /FLUXX/ zf,Tbar,tau_obs,u_s_obs,v_s_obs,w_s_obs
