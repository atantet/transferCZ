      SUBROUTINE solve(x)
      implicit none
      include 'const.com'
      include 'mat.com'
*     IMPORT/EXPORT
      real      x(ndim)
*     COMMON
      integer  ipiv
      real     LU
      common /LU/ LU(ndim,ndim),ipiv(ndim)
*     LOCAL
      integer  info,iw(ndim)
      real y(ndim),w1(3*ndim),ferr,berr
      character*1 trns
*
      y=x
      LU=A
      trns = 'N'
      call f07adf(ndim,ndim,LU,ndim,ipiv,info)
      call f07aef(trns,ndim,1,LU,ndim,ipiv,x,ndim,info)
      call f07ahf(trns,ndim,1,A,ndim,LU,ndim,ipiv,y,ndim
     &           ,x,ndim,ferr,berr,w1,iw,info)
      END
********************************************************
      SUBROUTINE determinant(det)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      real      det
*     COMMON
      integer  ipiv
      real     LU
      common /LU/ LU(ndim,ndim),ipiv(ndim)
*     LOCAL
      integer   i,iscl,isclmx
*
      isclmx = 100
      det = 1.0
      DO i=1,ndim
        det = det*LU(i,i)/10.
        IF (ipiv(i).NE.i) det = -det
      ENDDO
      DO iscl = 1,isclmx
         IF (abs(det).GE.1.e+10) THEN
            det = det/1.e+10
         ELSE
            goto 100
         ENDIF
      ENDDO
 100  continue
      write(99,*) 'iscl: ', iscl
*
      END
*********************************************************
      SUBROUTINE solvd(x,id)
      implicit none
      save ipiv,LU
*     save A
      include 'const.com'
      include 'wnl.com'
*     IMPORT/EXPORT
      integer id
      complex x(ndim)
*     LOCAL
      integer ipiv(ndim),inf
      complex LU(ndim,ndim)
*
      IF (id.EQ.0) THEN
*      A = L
       call f07arf(ndim,ndim,L,ndim,ipiv,inf)
       LU = L
       IF (inf.NE.0) THEN
         write(99,*) 'info f07arf', inf
         STOP
       ENDIF
*      id = 1
      ENDIF
      call f07asf('N',ndim,1,LU,ndim,ipiv,x,ndim,inf)
      IF (inf.NE.0) THEN
       write(99,*) 'info f07asf', inf
       STOP
      ENDIF
*
      END   
*********************************************************
      SUBROUTINE backs(w,sig,ub,err,res)
      implicit none
      include 'const.com'
      include 'mat.com'
*     IMPORT/EXPORT
      complex   w(ndim,nf),sig(nf)
      real      err(nf),ub(ndim),ut(ndim),res
*     LOCAL     
      integer   i,j,L,itime
      complex   ev(ndim)
      real      e(ndim,2),B(ndim,ndim)
*     (A-SIG*B)*W 
      itime = 0
      ut = 0.0

      call matB
      B = A

      DO L = 1,nf
       call matA
       A = A - real(sig(L))*B
       e = 0.0
       DO i=1,ndim
        DO j=1,ndim
         e(i,1) = e(i,1) + A(i,j)*real(w(j,L))
         e(i,2) = e(i,2) + A(i,j)*aimag(w(j,L))
        END DO
       END DO
       DO i=1,ndim
        DO j=1,ndim
         e(i,1) = e(i,1) + aimag(sig(L))*B(i,j)*aimag(w(j,L))
         e(i,2) = e(i,2) - aimag(sig(L))*B(i,j)*real(w(j,L))
        END DO
        ev(i) = cmplx(e(i,1),e(i,2))
       END DO
       err(L) = real(dot_product(ev,ev))
      ENDDO
      write(7,999) (err(L),L=1,10)
      res=err(1)
 999  format(10(1x,e12.4))
*
      END
*********************************************************
      SUBROUTINE qz(ub,ev,ew)

      implicit none
      include 'const.com'
      include 'mat.com'
*     INPUT?OUTPUT
*     LOCAL
      integer ifail,n,j,k,jmax,itime,iter(ndim),info
      real B(ndim,ndim),vec(ndim,ndim),alfr(ndim),alfi(ndim),z(ndim)
      real eps1,err(nf),res
      real ub(ndim),ut(ndim)
      logical matv
      complex ev(ndim,nf)
      complex ew(nf),r(ndim),v0(ndim),r0,v(ndim,ndim)
      character*1 JOBVR
      Real (Kind=8) :: dummy(1,1)
      
       eps1 = 1.0e-6
      ifail = -1
      matv = .true.

      ut = 0.0
      itime = 0
      call fields(ub,ut,itime,0.)
      call matB
      B = A
      call matA
      IF (matv) THEN
       JOBVR = 'vec'
      ELSE
       JOBVR = 'ndim'
      ENDIF
      CALL DGGEV('ndim',JOBVR,ndim,A,ndim,B,ndim,alfr,alfi,z,dummy,1,
     & vec,ndim,dummy,-1,info)
      IF (INFO.EQ.0) THEN

*      call f02bjf(ndim,A,ndim,B,ndim,eps1,alfr,alfi,z,
*     & matv,vec,ndim,iter,ifail)

      write(99,*) ' in eigen  qz ',ifail
      k = 0
      DO n=1,ndim
	IF (abs(z(n)).LE.(1.0e-8)) THEN
        r(n) = -cmplx(1.0e+10,0.0)
       ELSE
        IF (alfi(n).EQ.0.0) THEN
         r(n) = cmplx(alfr(n),0.0)/z(n)
         DO j=1,ndim
          v(j,n) = cmplx(vec(j,n),0.0)
         END DO
        ELSE
         r(n) = cmplx(alfr(n),alfi(n))/z(n)
         DO j=1,ndim
          v(j,n)   = cmplx(vec(j,n-k),(-1)**(k+2)*vec(j,n+1-k))
         END DO
         k = 1 - k
        END IF 
       END IF 
      END DO
*
      DO k=1,ndim-1
       r0 = r(k)
       jmax=k
       DO j=k+1,ndim
        IF (real(r(j)).GT.real(r0)) THEN
         r0 = r(j)
         jmax=j
        ENDIF
       ENDDO
       r(jmax)=r(k)
       r(k)=r0
       DO j=1,ndim
        v0(j) = v(j,jmax)
        v(j,jmax)=v(j,k)
        v(j,k)=v0(j)
       ENDDO
      ENDDO

      write(80,003) icp1,par(icp1)
      DO j=1,ndim
       write(80,003) j,real(r(j)),aimag(r(j))
      END DO
      write(81,003) icp1,par(icp1),(real(r(j)),j=1,12)
      write(82,003) icp1,par(icp1),(aimag(r(j)),j=1,12)

      DO j=1,10
       write(7,001) j,real(r(j)),aimag(r(j))
      END DO
      DO j=1,nf
       ew(j) = r(j)
       DO n=1,ndim
        ev(n,j) = v(n,j)
       END DO
      END DO
      call backs(ev,ew,ub,err,res)
      
      ENDIF
*
001   format(i4,1x,3(1x,e13.5))
002   format(13(f10.4))
003   format(i4,1x,21(1x,f10.4))
004   format(3(1x,e13.5))

      END
*****************************************************************
      SUBROUTINE time_int(z,dt,nt)
*     z=uo, dt=dstep, nt=maxstep
      implicit none
      include 'const.com'
      include 'mat.com'
*     IMPORT/EXPORT
      integer   nt
      real      z(ndim),dt
*     LOCAL     
      integer   itime,i,k,ni,j
      real      zold(ndim),zt(ndim),zm(ndim),dz(ndim),f(ndim)
      real      eps,res,time,B(ndim,ndim)
      integer flag_nt, k_solve
      common /EXT_FLAG/ flag_nt, k_solve
*
      if (ioute.eq.1) then 
         write(41,*) nx+1,nt
         write(42,*) nx+1,nt
         write(43,*) nx+1,nt
         write(44,*) nx+1,nt
         write(45,*) ny+1,nt
         write(46,*) ny+1,nt
         write(47,*) ny+1,nt
         write(48,*) ny+1,nt
         write(49,*) nt,nx+1,ny+1
         write(50,*) nt,nx+1,ny+1
      endif 
      itime = 1
      eps = 1.0e-05
!      ni = 9
      time = 0.0
      dz = 0.0
      zold = 0.0

      write(51,*) '      time        TW        TE         hW        hE'
      print *, "Starting integration"
      DO k=1,nt
         print *, "Iteration ", k, ", time-step ", time
         time = time + dt
         res = 1.
         write(40,*)
         flag_nt = 1
!         DO i=1,ni
         i = 1
         DO WHILE (res.GT.eps)
            k_solve = k
            zt = (z - zold)/dt  
            zm = (z + zold)/2.  
            call fields(zm,zt,itime,time-dt/2.)
            flag_nt = 0
            call matB
            B = -A
            call matA
            A = .5*dt*A + B
            call rhs(f)
            dz = dt*f
            call solve(dz)
            res = sqrt(dot_product(dz,dz))
            write(40,001) time,i,res
            z = z + dz
            i = i + 1
            print *, res
         END DO
         print *, 'Exciting integration loop after ', i, ' iterations.'
         call outhov(zm,time-dt/2.)
         zold = z
      END DO
      
 999  format(11(1x,e12.4))
 001  format(1x,'time =',f12.4,'step =',i3,1x,'residue =',e13.5)
*
      END
      
