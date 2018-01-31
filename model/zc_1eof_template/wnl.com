      complex adjoint(ndim)
      complex L(ndim,ndim)
      COMMON /MATRIC/ L,adjoint
