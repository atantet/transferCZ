      PROGRAM continu
*************************************************************************
* 
* This code is an implicit implementation of the Zebiak-Cane model
* The user can continue steady-states of the model versus parameters 
* and calculate the stability of these states. However, also only 
* time-stepping can be performed. 
* 
* bstep: main continuation routine
* time_int: main time stepping routine
* 
*************************************************************************
      implicit none
      include 'const.com'
*     COMMON
      real     uo,un,upo,w,time
      integer itime, k
      integer flag_nt, k_solve
      common /FLD/ uo(ndim),un(ndim),upo(ndim),w(ndim,nf)
      common /EXT_FLAG/ flag_nt, k_solve
      integer  i,ibr,ipnt,itp,icpo,irs,isw,lab
      integer  ise,istime,isf,timestep
      real     dtime 
      real     xlo,xln,xlpo,deto,detn
      real     sigo,sign,s0,s1
      real     sig(nf,2),g05caf,x
      complex  ww(ndim,nf),sigma(nf)
*
      ww = cmplx(0.,0.)
      sigma = cmplx(0.,0.)

************************************************************************
* Modify only the lines within this box                                *
************************************************************************
*
!      call random_seed()
      istime = 1   ! 0: continuation, 1: time integration
      ioute  = 0   ! 0: basic output, 1: extensive output 
      ioutnd  = 1   ! 0: dimensionless output off, 1: dimensionless output on
      ioutd  = 0   ! 0: dimensional output off, 1: dimensional output on
! control for continuation 
      if (istime.eq.0) then
         irs   = 0              ! starting label 
!     dstep = 0.0005   ! continuation step 
!     maxstep = 1      ! number of steps 
         isw   = 0              ! switch (yes = 1, no = 0) 
         icp1  = 4              ! control parameter, see parameters.f
         sw_A  = 0              ! new atmospheric features (yes = 1, no = 0) 
         ise = 0                ! stability at every point 
         isf = 0                ! stability only at end point 
      else 
!     control for time integration
         irs= 0 
!     Note: see const.com for the timestep parameters
!     maxstep = 5     ! number of time steps
!     dstep = 0.4      ! time step  dt = 1 -> dt* = 0.237 yr or 0.244?
         sw_A = 1               ! new atmospheric features (yes = 1, no = 0)
      endif 
************************************************************************
*     Indicate initialization step
      k_solve = 0
      flag_nt = 0

      IF ((icp1.EQ.4).OR.(icp1.EQ.30).OR.(icp1.EQ.29).OR.(icp1.EQ.27)) 
     &     sw_A = 1
*     
!     Initial state 
      IF (irs.EQ.0) THEN
         call parameters
         call init
         call anpnt(lab,ibr,ipnt,itp,icpo,deto,sigo,xlo,xlpo)
      ELSE
         call repnt(irs,ibr,ipnt,itp,icpo,deto,sigo,xlo,xlpo)
         call init
      END IF
      s0=0.
      lab=irs
      

* 
! Start of Continuation/Time Integration
*
      if (istime.eq.0) then 
!     Tangent Calculation 
         call prpar(dstep,maxstep,isw,irs)
         IF (icp1.NE.icpo) THEN
            call newdir(uo,upo,xlpo)
            xlo=par(icp1)
            write(99,*) xlo,xlpo
         END IF
!     Main Continuation Step
         DO k=1,maxstep
            call bstep(lab,ibr,itp,dstep,s0,s1,deto,detn,xlo,xln,uo,un,
     +           ipnt,isw,sigo,sign,upo,xlpo,ww)
            xlo=xln
            s0=s1
            deto=detn
            sigo=sign
            uo=un
            if (ise.eq.1) then
               call qz(uo,ww,sigma)
            endif 
            IF (itp.EQ.4) GOTO 100
         ENDDO
 100     CONTINUE
         if (isf.eq.1) then 
            call qz(uo,ww,sigma)
            call outper(ww,sigma,36)
            call analysis(ww,sigma,70)
         endif 
!     Write steady field to fort.33 + Eq. fields to fort.34 
         call outclim(uo)
         call writepar
      endif		! istime

* 
! Start of Time Integration
*
      if (istime.eq.1) then 
!     choose initial condition 
*
*     call wnl(sigma(1),uo,ww(1,1))
*     sigo = sig(1,1)
*     DO i=1,ndim
*      uo(i) = uo(i) + 0.2*(g05caf(x)/2.-1. )
*      uo(i) = g05caf(x)/2.-1. 
*     uo(i) = uo(i) + 0.2*real(ww(i,1))
*     END DO
         uo = 0.0
         call time_int(uo,dstep,maxstep)
         call writepar
      endif                     ! time loop 
*     
!     Restarting information
*     
      call wrpnt(lab,ibr,ipnt,itp,icp1,deto,sigo,xlo,xlpo)
 999  format(2x,i4,4x,4e14.4)
      
      END
********************************************************
      SUBROUTINE bstep(lab,ibr,itp,ds,s0,s1,deto,detn,xlo,xln,uo,un,
     +                 ipnt,isw,sigo,sign,upo,xlpo,ww)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      integer   lab,ipnt,ibr,itp,isw     
      real      s0,s1,ds
      real      deto,detn,xlo,xln,uo(ndim),un(ndim)
      real      xlpo,upo(ndim),sigo,sign
      complex   ww(ndim,nf)
*     LOCAL
      integer   i
      real      xlpa,xlpb,f0,f1,uza,uzb,sigi,epsi,sf,dsw
      complex   sigma(nf)
      logical   syes
*     FUNCTIONS
      real      uszr
* 
      itp=9
      syes=.false.
      xlpb=xlpo
      uzb=uszr(uo)
      call step(detn,ds,xlo,xln,uo,un,upo,xlpo)
*     call qz(un,ww,sigma)
      uza=uszr(un)
      xlpa=(xln-xlo)/ds
*     sign=sig(1)
      sign=0.0
      sigi=0.0
      epsi=1.0e-06
      s1=s0+ds
*
*        CHECK FOR SINGULARITY
*
      IF ((detn*deto.LT.0.0).AND.(detn.GT.0.0)) THEN
         write(99,'(A)') ' Bifurcation point '
         itp=6
         f0=deto
         f1=detn
         syes=.true.
      END IF
      IF (xlpa*xlpb.LT.0.0) THEN
         write(99,'(A)') ' Limit point '
         itp=5
         f0=xlpb
         f1=xlpa
*        syes=.true.
      END IF
      IF ((sign*sigo.LT.0.0).AND.(sigi.GT.epsi)) THEN
         write(99,'(A)') ' Hopf bifurcation point '
         itp=7
         f0=sigo
         f1=sign
*        syes=.true.
       END IF
       IF (sign*sigo.LT.0.0) THEN
         write(99,'(A)') ' EW - Bifurcation point '
         itp=8
         f0=sigo
         f1=sign
*        syes=.true.
       END IF
       IF (uza*uzb.LT.0.0) THEN
         write(99,'(A)') ' User defined point '
         itp=4
         f0=uzb
         f1=uza
         syes=.true.
       END IF
*      syes=.false.
       IF (syes) THEN
          write(99,'(A)') '*****************************************'
          call detect(s0,s1,f0,f1,sf,xlo,xln,uo,un,itp,detn,
     +                ww,sigma,upo,xlpo)
          ipnt=ipnt+1
*         call power(un,ww,sigma)
*         call qz(un,ww,sigma)
*         sign=sig(1,1)
*
*  BRANCHSWITCHING IF ISW > 0 AND ITP = 6 (BIFURCATION POINT)
*
          IF ((itp.EQ.8).AND.(abs(isw).GT.0)) THEN
               DO i=1,ndim
                 uo(i)=un(i)
               ENDDO
	        dsw=ds
               write(99,999) dsw,isw
               call switch(uo,isw,upo,xlpo)
               call step(detn,dsw,xlo,xln,uo,un,upo,xlpo)
*              call power(un,ww,sig)
*              call qz(un,ww,sig)
               sign=real(sigma(1))
               ibr=ibr+1
               ipnt=1
               s1=sf+dsw
               itp=9
*
*       SOLUTION AT THE NEW BRANCH
*
               ds=dsw
               isw=0 
           ELSE
               IF (itp.EQ.6) THEN
                   write(99,'(A)') ' Bifurcation point, no switch '
               END IF
           END IF
           write(99,'(A)') ' **************************************'
      ENDIF
      write(99,'(A)') ' Regular point '
      xlpo=(xln-xlo)/ds
      DO i=1,ndim
         upo(i)=(un(i)-uo(i))/ds
      ENDDO
      IF (itp.eq.4) RETURN
*
*       END CHECKING FOR SINGULARITY
*
*     call monitor(un,itp,xln,detn,7)
      par(icp1)=xln
      ipnt=ipnt+1
 999  format(/1x,' Bifurcation point ',/1x,
     +           ' Branch switch with DS ',  e14.5,/1x,
     +           ' ISW = ',i8)
      END
******************************************************************
      SUBROUTINE detect(s0,s1,fs0,fs1,sf,xlo,xls,uo,un,itp,detc,
     +                  ww,sigma,upo,xlpo)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      integer  itp
      real     s0,s1,fs0,fs1,sf,xlo,xls,uo(ndim),un(ndim),detc,
     +         upo(ndim),xlpo
      complex  sigma(nf),ww(ndim,nf)
*     LOCAL
      integer itmx,k
      real    sa,q0,eps,q1,q2,fq0,fq1,fq2
*     FUNCTION
      real    uszr
*    
      write(99,'(A)') ' Detection singularity  '
      eps=1.0e-04
      sa=s0
      q0=s0
      q1=s1
      fq0=fs0
      fq1=fs1
      itmx=5
      DO k=1,itmx
         IF (abs(q0-q1).GT.eps) THEN
            q2=(q0*fq1-q1*fq0)/(fq1-fq0)
            call step(detc,q2-sa,xlo,xls,uo,un,upo,xlpo)
            IF (itp.EQ.5) THEN
               fq2=(xls-xlo)/(q2-sa)
            END IF
            IF (itp.EQ.6) THEN
               fq2=detc
            END IF
            IF ((itp.EQ.7).OR.(itp.EQ.8)) THEN
*              call power(un,ww,sigma)
*              call qz(un,ww,sigma)
               fq2=real(sigma(1))
*              fq2=sig(2,1)
            END IF
            IF (itp.EQ.4) THEN
               fq2=uszr(un)
            END IF
            write(99,999) k,q2,fq2
            q0=q1
            fq0=fq1
            q1=q2
            fq1=fq2
            write(99,998) abs(q0-q1)
         END IF
      ENDDO
      sf=q2
      write(99,997) sf,xls
 999  format(' Information secant proces ',
     +       /2x,' Iteration :',i6,' s :',e12.5,2x,'F(s) : ',e12.4)
 998  format(2X,'Residu  :',e12.4,/)
 997  format(/1X,'Singularity at s :',e12.5,2x,'xl : ',e12.5)
      RETURN
      END
********************************************************************
      SUBROUTINE kern(unc,phi1)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      real    unc(ndim),phi1(ndim)
*     LOCAL
      integer i,L,itmx,itime
      real    phi0(ndim),res,eps,xnpt,ut(ndim)
*     FUNCTIONS
      real    f06ejf
 
	  
      itime=0
      ut = 0.0
      DO i=1,ndim
         phi1(i)= 1.0
      ENDDO
      WRITE(99,'(A)') ' INFORMATION INVERSE ITERATION EIGENVECTOR J'
      res= 1.0
      eps=1.0e-06
      itmx=15
      call matA
      DO L=1,itmx
         IF (res.GT.eps) THEN
            DO i=1,ndim
               phi0(i)=phi1(i)
            ENDDO
            call solve(phi1)
            xnpt=f06ejf(ndim,phi1,1)
            DO i=1,ndim
               phi1(i)=phi1(i)/xnpt
               phi0(i)=abs(phi1(i))-abs(phi0(i))
            ENDDO
            res=f06ejf(ndim,phi0,1)
            write(99,999) res,xnpt
         END IF
      ENDDO
 999  format(1x,'Residu :',e12.4,2x,' Norm :',e12.4)
      END
*********************************************************
      SUBROUTINE maxr(vect,maxv,ndim,dl)
      implicit none
*     IMPORT/EXPORT
      integer   ndim
      real      vect(ndim),maxv,dl
*     LOCAL
      integer   j,jmax
*
      maxv=0.0
      DO j=1,ndim
         IF (abs(vect(j)).GT.maxv) THEN
            maxv=abs(vect(j))
            jmax=j
         END IF
      ENDDO
      IF (abs(dl).GT.maxv) THEN
         maxv=abs(dl)
      END IF
      write(99,*) 'element no:',jmax,' max = ',maxv
      END
********************************************************************
      SUBROUTINE newdir(uo,upo,xlpo)
*     Calculates a new branch direction after change of primary bifur-
*     cation parameter.
*
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      real     uo(ndim),upo(ndim),xlpo
*     LOCAL
      integer  i,itime
      real     rb(ndim,2),xnorm,ut(ndim),rhulp(ndim),dt
*     FUNCTIONS
      real    f06ejf

      itime=0
      ut=0.0
      call fields(uo,ut,itime,0.)
      call matA
      call rlbrs(uo,ut,rb)
      DO i=1,ndim
       rhulp(i)=-rb(i,2)
      ENDDO
      call solve(rhulp)
      xnorm=sqrt(f06ejf(ndim,rhulp,1)**2+1.0)
      DO i=1,ndim
       upo(i)=rhulp(i)/xnorm
      ENDDO
      xlpo=1.0/xnorm
      write(99,*) 'Newdir done'
      END
********************************************************************************
      SUBROUTINE prpar(ds,nst,isw,irs)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      integer  nst,isw,irs
      real     ds
*     LOCAL
      integer  i
*
      write(7,996)
      write(7,999) (par(i),i=1,30)
      write(7,996)
      write(7,998) irs,nst,isw,ds
      write(7,996)
      write(7,997) icp1
 999  format('*',5e15.5)
 998  format('*',4x,'irs = ',i4,4x,' nstp = ',i4,3x,' isw = ',i4,
     +       5x,'ds = ',e12.4)
 997  format('*',3x,'typ',8x,'par(',i2,')',6x,
     +       'Testfun ',7x,'TWP  ',7x,'TCP ',7x,'TEP',7x,'TN',7x,'TS')
 996  format('*')
      END
********************************************************************************
      SUBROUTINE repnt(irs,ibr,ipnt,itp,icpo,deto,sigo,xlo,xlpo)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      integer  irs,ibr,ipnt,itp,lab,icpo
      real     deto,sigo,xlo,xlpo
      real     uo,un,upo,w        
      common /FLD/ uo(ndim),un(ndim),upo(ndim),w(ndim,nf)
*     LOCAL
      integer  i,L,nskip,nf1,ndim1
      logical  eof4
*     FUNCTION
      real     g05caf,x
*
      rewind 4
 100  CONTINUE   
         read(4,*,end=200) ibr,ipnt,itp,icpo,lab,nf1,ndim1,nskip
         IF (lab.EQ.irs) THEN
	    irs=lab
            IF (ndim.NE.ndim1) THEN
               write(99,*) 'dimension of inital point wrong'
               STOP
            ENDIF
            IF (nf.NE.nf1) THEN
               write(99,*) 'Eigenv. No of initial point wrong'
*              STOP
            ENDIF
            read(4,999,end=200) xlo,xlpo,deto,sigo
            DO i=1,ndim
               read(4,999,end=200) uo(i),upo(i),(w(i,L),L=1,nf)
*              read(4,999,end=200) uo(i)
*              upo(i)=0.
*              DO L=1,nf
*                 w(i,L)=2*g05caf(x)-1.0
*              ENDDO
            ENDDO
            read(4,999,end=200) (par(i),i=1,30)
            write(99,*) 'readpnt done'
         ELSE
            call skip4(nskip,eof4)
            IF (eof4) GOTO 200
            GOTO 100
         END IF
      RETURN
 200  STOP 'repnt failed'
*
 999  format(2X,7e18.10)
      END
********************************************************************
      SUBROUTINE rlbrs(un,ut,rb)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      real     un(ndim),ut(ndim),rb(ndim,2)
*     LOCAL
      integer  i,itime
      real     delta
*
      itime = 0
      delta=1.0e-08
      call rhs(rb(1,1))
      par(icp1)=par(icp1)+delta
      call fields(un,ut,itime,0.)
      call rhs(rb(1,2))
      par(icp1)=par(icp1)-delta
*
*     DERIVATIVE WITH RESPECT TO PRIMARY PARAMETER
*
      DO i=1,ndim
       rb(i,2)=-(rb(i,2)-rb(i,1))/delta
      ENDDO
      END
******************************************************************
      SUBROUTINE skip4(nskip,eof4)
      implicit none
*     IMPORT/EXPORT
      integer  nskip
      logical  eof4
*     LOCAL
      integer  i
*
      eof4=.false.
      DO i=1,nskip
         read(4,999,end=100)
      ENDDO 
 999  format(1x)
      return
 100  CONTINUE
      eof4=.true.
      return
      END
***************************************************************
      SUBROUTINE solvbr(rl,upo,xlpo,ndim,rlnew,dl)
      implicit none
*     IMPORT/EXPORT
      integer  ndim
      real     rl(ndim,2),upo(ndim),xlpo,dl,rlnew
*     LOCAL
      integer  i
      real     res(2)
*
      res(1)=0.0
      res(2)=0.0
      DO i=1,ndim
         res(1)=res(1)+upo(i)*rl(i,1)
         res(2)=res(2)+upo(i)*rl(i,2)
      ENDDO
      dl=(rlnew-res(1))/(xlpo-res(2))
      DO i=1,ndim
         rl(i,1)=rl(i,1)-dl*rl(i,2)
      ENDDO
      END
******************************************************************
      SUBROUTINE step(det,ds,xlo,xln,uo,un,upo,xlpo)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      real     det,ds,xlo,xln,uo(ndim),un(ndim),upo(ndim),xlpo
*     COMMON
      real     phl
      common /DD10/ phl(ndim)
*     LOCAL
      integer  i,k,itmx,itime
      real     rb(ndim,2),dl,res,reseps,rlo,rln,ut(ndim)
*     
      itmx=10
      reseps=1.0e-6
      res=1.0e+00
      write(99,*) ' Information newton process '
      write(99,999) itmx,reseps
      call updbr(upo,uo,ndim,ds,xlo,xlpo,rlo,rln,0)
      xln=xlo+xlpo*ds
      un = uo + ds*upo
      par(icp1)=xln
*
      itime=0
      ut=0.0
      DO k=1,itmx
         IF (res.GT.reseps) THEN
            call fields(un,ut,itime,0.)
            call matA
            call rlbrs(un,ut,rb)
            call updbr(upo,un,ndim,ds,xln,xlpo,rlo,rln,1)
            DO i=1,2
             call solve(rb(1,i))
            ENDDO 
            call solvbr(rb,upo,xlpo,ndim,rln,dl)
            call upd(un,rb,ndim,dl,xln)
            call maxr(rb,res,ndim,dl)
            par(icp1)=xln
            write(99,998) res,xln
         END IF
         IF ((k.EQ.itmx).AND.(res.GT.reseps)) THEN
            write(99,*) ' *** Newton process not converged'
*           call wrfld(un)
            STOP
         END IF
      ENDDO
      DO i=1,ndim
	 phl(i)=rb(i,2)
      ENDDO
*     call determinant(det)
      call monitor(un,0,xln,det,7)
 999  format(1x,' max iterations :',2x,i8,2x,' Conv. crit. :',e12.4)
 998  format(9x,'Residu :',e12.5,1x,'Par :',e12.5)
      END
********************************************************************
      SUBROUTINE switch(unc,isw,upo,xlpo)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      integer  isw
      real     unc(ndim),upo(ndim),xlpo
*     COMMON
      real     phl
      common /DD10/ phl(ndim)
*     LOCAL
      integer  i
      real     gmax,eps,phi1(ndim),rf1,rf2
*
      eps=1.0e-06
      gmax=0.0
      DO i=1,ndim
         IF (abs(phl(i)).GT.gmax) THEN
            gmax=abs(phl(i))
         END IF
      ENDDO
      write(99,999) gmax
      call kern(unc,phi1)
      IF (gmax.LT.eps) THEN
         DO i=1,ndim
            upo(i)=isw*phi1(i)
         ENDDO
         xlpo=0.0
      ELSE
         rf1=0.0
         rf2=0.0
         DO i=1,ndim
            rf1=rf1+upo(i)*phi1(i)
            rf2=rf2+upo(i)*phl(i)
         ENDDO
         xlpo=-rf1/(xlpo-rf2)
         DO i=1,ndim
            upo(i)=isw*(phi1(i)-xlpo*phl(i))
         ENDDO
      END IF
 999  format(/1x,' gmax :',e14.7)
      END
********************************************************************
      SUBROUTINE upd(un,rb,ndim,dl,xln)
      implicit none
*     IMPORT/EXPORT
      integer  ndim
      real     rb(ndim,2),un(ndim),xln,dl
*     LOCAL
      integer  i
*
      DO i=1,ndim
         un(i)=un(i)+rb(i,1)
      ENDDO
      xln=xln+dl
      END
**************************************************************
      SUBROUTINE updbr(upo,uo,ndim,ds,xl,xlp,rlold,rlnew,inew)
      implicit none
*     IMPORT/EXPORT
      integer  ndim,inew
      real     upo(ndim),uo(ndim),ds,xl,xlp,rlold,rlnew
*     LOCAL
      integer  i
      real     res
*
      res=0.0
      DO i=1,ndim
         res=res+upo(i)*uo(i)
      ENDDO
      IF (inew.EQ.0) THEN
         rlold=res+xlp*xl+ds
      ELSE
         rlnew=rlold-res-xl*xlp
      END IF
      END
********************************************************************
      real FUNCTION uszr(un)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      real     un(ndim)
*
      IF (icp1.EQ.4) THEN
       uszr=par(icp1)-0.2
      END IF
      IF (icp1.EQ.3) THEN
*      g=par(icp1)-6.5538
      END IF
      IF (icp1.EQ.5) THEN
       uszr=par(icp1)-0.6
      END IF
      IF (icp1.EQ.14) THEN
       uszr=par(icp1)-0.5
      END IF
      IF (icp1.EQ.29) THEN
       uszr=par(icp1)-0.3
      END IF
      IF (icp1.EQ.27) THEN
       uszr=par(icp1)-0.1
      END IF
      IF (icp1.EQ.15) THEN
       uszr=par(icp1)-6.67
      END IF
      IF (icp1.EQ.30) THEN
       uszr=par(icp1)-3.05
      END IF

      END
******************************************************************
      SUBROUTINE wrpnt(lab,ibr,ipnt,itp,icp,deto,sigo,xlo,xlpo)
      implicit none
      include 'const.com'
*     IMPORT/EXPORT
      integer  lab,ibr,ipnt,itp,icp
      real     deto,sigo,xlo,xlpo
*     COMMON
      real     uo,un,upo,w        
      common /FLD/ uo(ndim),un(ndim),upo(ndim),w(ndim,nf)
*     LOCAL
      integer  i,L,nskip,inpri
*
      lab=lab+1
      inpri=int((2+nf)/7)+1
      nskip=1+ndim*inpri+5
      write(3,999) ibr,ipnt,itp,icp,lab,nf,ndim,nskip
      write(3,998) xlo,xlpo,deto,sigo
      DO i=1,ndim
        write(3,998) uo(i),upo(i),(w(i,L),L=1,nf)
      ENDDO
      write(3,998) (par(i),i=1,30)
*
 999  format(2X,7i6,2i10)
 998  format(2X,7e18.10)
      END
*****************************************************************
      real FUNCTION xrm(u,ndim,iopt)
      implicit none
*     IMPORT/EXPORT
      integer  ndim,iopt
      real     u(ndim)
*     LOCAL
      integer  i
      real     x
*     FUNCTION
      real     f06ejf
*
*     IOPT=0 : MAX NORM
*     IOPT=1 : L1 NORM
*     IOPT=2 : L2 NORM
*
      x=0.0e+00
      IF (iopt.EQ.0) THEN
         DO i=1,ndim
            IF (abs(u(i)).GT.x) THEN
               x=abs(u(i))
            END IF
         ENDDO
         xrm=x
      END IF
      IF (iopt.EQ.1) THEN
         DO i=1,ndim
            x=x+abs(u(i))
         ENDDO 
         xrm=x
      END IF
      IF (IOPT.EQ.2) THEN
         xrm=f06ejf(ndim,u,1)
      END IF
      END
********************************************************************
