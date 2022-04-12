!**************** SUBROUTINE KARDIR ************************************
! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE KARDIR(NMAX,V,BASIS)!{{{
      USE constants
      !IMPLICIT REAL(dp) (A-H,O-Z)
      INTEGER(I4B),INTENT(IN)  :: NMAX
      REAL(DP)  :: V(3,NMAX),BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE!}}}

!**************** SUBROUTINE DIRKAR ************************************
! transform a set of vectors from
! ) direct lattice      (BASIS must be equal to A direct lattice)
! ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
! to cartesian coordinates
!***********************************************************************

!      SUBROUTINE DIRKAR(NMAX,V,BASIS)!{{{
!      USE constants
!      INTEGER(I4B),INTENT(IN) :: NMAX
!      REAL(DP)  :: V(3,NMAX),BASIS(3,3)
!
!      DO N=1,NMAX
!        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(1,2)+V(3,N)*BASIS(1,3)
!        V2=V(1,N)*BASIS(2,1)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(2,3)
!        V3=V(1,N)*BASIS(3,1)+V(2,N)*BASIS(3,2)+V(3,N)*BASIS(3,3)
!        V(1,N)=V1
!        V(2,N)=V2
!        V(3,N)=V3
!      ENDDO!}}}


!*********************************************************************
! subroutine IONCGR
! subroutine performes a steepest descent/ conjugate gradient step
! on the ions and the cell
! it contains considerable heuristic to make the minimization
! efficient and
! uses a varient of Brents algorithm from numerical recipies
! for the line minimization
! especially the question how long a line minimization should be
! continued is of requires a lot of fiddling
!
! IFLAG
!   on call
!           0  initial trial step (steepest descent)
!           1  initial trial step (conjugate gradient)
!           2  move ions to minimum
!   on exit
!           0  currently line optimization is done
!           1  new trial step has been performed, main can break
!           2  energy is possibly converged
! F         contains scaled forces on ions  divided by mass of ions
! FACT      scaling factor of forces
! FSIF      stress (in cartesian units) scaled by FACTSI
! FACTSI    scaling factor for stress
! POSION    coordinates of ions
! A,B       direct and reciprocal lattice
!
! FL        used to store old forces  (set by IONCGR)
! S         step  into which search is performed    (set by IONCGR)
! POSIONC   old coordinates (set by IONCGR)
! POSIONOLD coordinates that were not recalculated to be in the unit cell
! some things which are used in our heuristics:
! EBREAK    absolut accuracy of energies (should be set by caller)
!           determines whether cubic interpolation is used
! EDIFFG    required accuracy for energy (>0)
! LSTOP2    external routine can signal that a break condition was met
!
! routine modified by Robin Hirschl in Oct. 2000 to avoid
! too large step when search was started near a saddle point
!*********************************************************************
      SUBROUTINE IONCGR(IFLAG,NIONS,TOTEN,A,B,NFREE,POSION,POSIOC, &!{{{
     &      FACT,F,FACTSI,FSIF,FL,S,DISMAX,IU6,IU0, &
     &      EBREAK,EDIFFG,E1TEST,LSTOP2)
      USE constants
      USE parameters , only : Lpbc
#ifdef MPI
      use smpi_math_module
#endif

      IMPLICIT REAL(dp) (A-H,O-Z)

      DIMENSION F(3,NIONS),FL(3,NIONS),S(3,NIONS)
      DIMENSION POSION(3,NIONS),POSIOC(3,NIONS)
      DIMENSION A(3,3),B(3,3)
      DIMENSION TMP(3)
      LOGICAL   LRESET,LTRIAL,LBRENT
      LOGICAL   LSTOP2

      SAVE  AC,FSIFL,SSIF
      DIMENSION AC(3,3),FSIFL(3,3),SSIF(3,3),FSIF(3,3)

      SAVE TOTEN1,DMOVED,E1ORD1,ORTH,GAMMA,GNORM,GNORMF,GNORML
      SAVE GNORM1,GNORM2,STEP,CURVET,SNORM,SNORMOLD
      SAVE DMOVEL,LTRIAL,LBRENT
      DATA ICOUNT/0/, ICOUNT2/0/, STEP /1_dp/,SNORM/1E-10_dp/
      DATA LTRIAL/.FALSE./

!=======================================================================
!  if IFLAG =0 initialize everything
!=======================================================================
      IF (IFLAG==0) THEN
        DO NI=1,NIONS
        DO I=1,3
          S(I,NI) =0
          FL(I,NI)=0
        ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          SSIF (I,J)=0
          FSIFL(I,J)=0
        ENDDO
        ENDDO
        LTRIAL=.FALSE.
        ICOUNT=0
        SNORM =1E-1_dp
        SNORMOLD =1E10_dp
        CURVET=0
      ENDIF

! if previous step was a trial step then continue with line minimization
!print *,'in cg method',ltrial 
      !> for confined contain the real pos
      IF (LTRIAL) THEN
        GOTO 400
      ENDIF
!=======================================================================!{{{
!  calculate quantities necessary to conjugate directions
!=======================================================================
      GNORML=0
      GNORM =0
      GNORMF=0
      ORTH  =0
      IF (FACT/=0) THEN
      DO NI=1,NIONS
         GNORML = GNORML+1_dp/FACT* &
     &    (FL (1,NI)*FL (1,NI)+FL(2,NI) *FL(2,NI) +FL(3,NI) *FL(3,NI))
         GNORM  = GNORM+ 1_dp/FACT*( &
     &    (F(1,NI)-FL(1,NI))*F(1,NI) &
     &   +(F(2,NI)-FL(2,NI))*F(2,NI) &
     &   +(F(3,NI)-FL(3,NI))*F(3,NI))
         GNORMF = GNORMF+1_dp/FACT* &
     &    (F(1,NI)*F(1,NI)+F(2,NI)*F(2,NI)+F(3,NI)*F(3,NI))
         ! S  step  into which search is performed    (set by IONCGR)
         ORTH   = ORTH+  1_dp/FACT* &
     &    (F(1,NI)*S(1,NI)+F(2,NI)*S(2,NI)+F(3,NI)*S(3,NI))
      ENDDO
      ENDIF

      GNORM1=GNORMF

      IF (FACTSI/=0) THEN
      DO I=1,3
      DO J=1,3
         GNORML=GNORML+ FSIFL(I,J)*FSIFL(I,J)/FACTSI
         GNORM =GNORM +(FSIF(I,J)-FSIFL(I,J))*FSIF(I,J)/FACTSI
         GNORMF=GNORMF+ FSIF(I,J)* FSIF(I,J)/FACTSI
         ORTH  =ORTH  + FSIF(I,J)* SSIF(I,J)/FACTSI
      ENDDO
      ENDDO
      ENDIF
      GNORM2=GNORMF-GNORM1
!print *,'in cg gnorml',gnorml
!print *,'in cg gnorm',gnorm
!print *,'in cg gnormF',gnormF
!print *,'in cg orth',orth

      !CALLMPI_C( sum_chain( GNORM ))
      !CALLMPI_C( sum_chain( GNORML ))
      !CALLMPI_C( sum_chain( GNORMF))
      !CALLMPI_C( sum_chain( GNORM1))
      !CALLMPI_C( sum_chain( GNORM2))
      !CALLMPI_C( sum_chain( ORTH ))

!=======================================================================
!  calculate Gamma
!  improve line optimization if necessary
!=======================================================================
      IF (IFLAG==0) THEN
        ICOUNT=0
        GAMMA=0
      ELSE
!       this statement for Polak-Ribiere
        GAMMA=GNORM /GNORML
!       this statement for Fletcher-Reeves
!        GAMMA=GNORMF/GNORML
      ENDIF
!print *,'in cg method gamma',gamma

      GAMMIN=1_dp
      IFLAG=1
#ifdef MPI
      IF (IU0>=0 .and. parallel%isroot) &
#else
      IF (IU0>=0 ) &
#endif
      WRITE(6,30)CURVET,CURVET*GNORMF,CURVET*(ORTH/SQRT(SNORM))**2
      ! if (parallel%isroot) WRITE(IOUNITS(11),30)CURVET,CURVET*GNORMF,CURVET*(ORTH/SQRT(SNORM))**2

   30 FORMAT(' curvature: ',F12.2,' expect dE=',E10.3, &
                & ' dE for cont linesearch ',E10.3)!}}}
! required accuracy not reached in line minimization
! improve line minimization
! several conditions must be met:

! orthonormality not sufficient
!      WRITE(0,*) ORTH, MAX(GAMMA,GAMMIN),GNORMF/5, ABS(CURVET*(ORTH/SQRT(SNORM))**2), LSTOP2
      IF (ABS(ORTH)*MAX(GAMMA,GAMMIN)>ABS(GNORMF)/5 &
! expected energy change along line search must be larger then required accuracy
     &    .AND. &
     &   ( (EDIFFG>0 .AND. &
     &       ABS(CURVET*(ORTH/SQRT(SNORM))**2)>EDIFFG) &
! or force must be large enough that break condition is not met
     &    .OR. &
     &     (EDIFFG<0 .AND..NOT.LSTOP2) &
     &   ) &
! last call must have been a line minimization
     &    .AND. LBRENT &
     &  ) GOTO 400
!print *,'not goto 400'
!---- improve the trial step by adding some amount of the optimum step!{{{
      IF (ICOUNT/=0) STEP=STEP+0.2_dp*STEP*(DMOVEL-1)
!---- set GAMMA to zero if line minimization was not sufficient
      IF (5*ABS(ORTH)*GAMMA>ABS(GNORMF)) THEN
         GAMMA=0
         ICOUNT=0
      ENDIF
!---- if GNORM is very small signal calling routine to stop
      IF (CURVET/=0 .AND.ABS((GNORMF)*CURVET*2)<EDIFFG) THEN
        IFLAG=2
      ENDIF

      ICOUNT=ICOUNT+1
      ICOUNT2=ICOUNT2+1
!-----------------------------------------------------------------------
! performe trial step
!-----------------------------------------------------------------------
      E1ORD1=0
      DMOVED=0
      SNORM =1E-10_dp

!----- set GAMMA to zero for initial steepest descent steps (Robin Hirschl)
!      have to discuss this
!      IF (ICOUNT2<=NFREE) GAMMA=0
      
      DO NI=1,NIONS
!----- store last gradient
        FL(1,NI)=F(1,NI)
        FL(2,NI)=F(2,NI)
        FL(3,NI)=F(3,NI)
!----- conjugate the direction to the last direction
        S(1,NI) = F(1,NI)+ GAMMA   * S(1,NI)
        S(2,NI) = F(2,NI)+ GAMMA   * S(2,NI)
        S(3,NI) = F(3,NI)+ GAMMA   * S(3,NI)

        IF (FACT/=0) THEN
        SNORM = SNORM  +  1/FACT * &
     &  (S(1,NI)*S(1,NI)+ S(2,NI)*S(2,NI) + S(3,NI)*S(3,NI))
        ENDIF
      ENDDO
!print *,'in cg method snorm fl',snorm

      DO I=1,3
         DO J=1,3
            FSIFL(I,J)=FSIF(I,J)
            AC(I,J)   = A(I,J)
        SSIF(I,J) = FSIF(I,J)+ GAMMA* SSIF(I,J)
      ENDDO
      ENDDO

      IF (FACTSI/=0) THEN
         DO I=1,3
            DO J=1,3
               SNORM = SNORM  + 1/FACTSI *      SSIF(I,J)* SSIF(I,J)
            ENDDO
         ENDDO
      ENDIF
!print *,'in cg method snorm S'
!print *,'fsif',fsif
!print *,'in cg method snorm S'
!print *,'a',A
!print *,'factsi',factsi
!print *,'in cg method snorm',snorm,'snormld',snormold

      !CALLMPI_C( sum_chain( SNORM))

!----- if SNORM increased, rescale STEP (to avoid too large trial steps) 
!      (Robin Hirschl)
      IF (SNORM>SNORMOLD) THEN
         STEP=STEP*(SNORMOLD/SNORM)
      ENDIF
      SNORMOLD=SNORM
      
!print *,'in cg method step',step,'snorm',snorm,'snormold',snormold
!print *,'in cg method S'
!print *,'s',S
!print *,'in cg method position before'
!print *,'posion',posion
      DO NI=1,NIONS
!----- search vector from cartesian to direct lattice
         TMP(1) = S(1,NI)
         TMP(2) = S(2,NI)
         TMP(3) = S(3,NI)
         CALL KARDIR(1,TMP,B)

!----- trial step in direct grid
         POSIOC(1,NI)= POSION(1,NI)
         POSIOC(2,NI)= POSION(2,NI)
         POSIOC(3,NI)= POSION(3,NI)

         POSION(1,NI)= TMP(1)*STEP+POSIOC(1,NI)
         POSION(2,NI)= TMP(2)*STEP+POSIOC(2,NI)
         POSION(3,NI)= TMP(3)*STEP+POSIOC(3,NI)
         DMOVED= MAX( DMOVED,S(1,NI)*STEP,S(2,NI)*STEP,S(3,NI)*STEP)

!----- keep ions in unit cell (Robin Hirschl)
         POSION(1,NI)= MOD(POSION(1,NI),1.d0)
         POSION(2,NI)= MOD(POSION(2,NI),1.d0)
         POSION(3,NI)= MOD(POSION(3,NI),1.d0)

!----- force * trial step = 1. order energy change
         IF (FACT/=0) THEN
            E1ORD1= E1ORD1 - 1_dp * STEP / FACT * &
     &           (S(1,NI)*F(1,NI)+ S(2,NI)*F(2,NI) + S(3,NI)*F(3,NI))
         ENDIF
      ENDDO
!print *,'in cg method new position, e1ord1',e1ord1,'step',STEP,'fact',fact
!print *,'posion',posion

      IF (FACTSI/=0) THEN
          DO I=1,3
          DO J=1,3
             E1ORD1= E1ORD1 - STEP / FACTSI * SSIF(I,J)* FSIF(I,J)
          ENDDO
          ENDDO
      ENDIF

      !> comment for isolate
      if(Lpbc)then
        DO J=1,3
           DO I=1,3
              A(I,J)=AC(I,J)
              DO K=1,3
                 A(I,J)=A(I,J) + SSIF(I,K)*AC(K,J)*STEP
              ENDDO
           ENDDO
        ENDDO
      endif

      !CALLMPI_C( sum_chain( E1ORD1 ))

      LRESET = .TRUE.
      X=0
      Y=TOTEN
      FP=E1ORD1
      IFAIL=0

      CALL ZBRENT(IU0,LRESET,EBREAK,X,Y,FP,XNEW,XNEWH,YNEW,YD,IFAIL)
      DMOVEL=1

#ifdef MPI
      IF (IU0>=0 .and. parallel%isroot) THEN
#else
      IF (IU0>=0 ) THEN
#endif
         WRITE(6,10) GAMMA,GNORM1,GNORM2,ORTH,STEP
         ! WRITE(IOUNITS(11),10) GAMMA,GNORM1,GNORM2,ORTH,STEP
 10      FORMAT(' trial: gam=',F8.5,' g(F)= ',E10.3, &
     &        ' g(S)= ',E10.3,' ort =',E10.3,' (trialstep =',E10.3,')')
         WRITE(6,11) SNORM
         ! WRITE(IOUNITS(11),11) SNORM
 11      FORMAT(' search vector abs. value= ',E10.3)
      ENDIF
      TOTEN1= TOTEN
      E1TEST=E1ORD1
      LTRIAL=.TRUE.
      RETURN!}}}

!=======================================================================!{{{
! calculate optimal step-length and go to the minimum
!=======================================================================
!-----------------------------------------------------------------------
!  1. order energy change due to displacement at the new position
!-----------------------------------------------------------------------
  400 CONTINUE
!print *,'here is 400'
      E1ORD2=0
      IF (FACT/=0) THEN
      DO NI=1,NIONS
        E1ORD2= E1ORD2 - 1_dp * STEP / FACT * &
     &  (S(1,NI)*F(1,NI)+ S(2,NI)*F(2,NI) + S(3,NI)*F(3,NI))
      ENDDO
      ENDIF

      IF (FACTSI/=0) THEN
      DO I=1,3
      DO J=1,3
        E1ORD2= E1ORD2 - STEP / FACTSI * SSIF(I,J)* FSIF(I,J)
      ENDDO
      ENDDO
      ENDIF

      !CALLMPI_C( sum_chain( E1ORD2 ))
!-----------------------------------------------------------------------
!  calculate position of minimum
!-----------------------------------------------------------------------
      CONTINUE
      LRESET = .FALSE.
      X=DMOVEL
      Y=TOTEN
      FP=E1ORD2
      IFAIL=0

!print *,'in cg method do 400'
!print *,'parameter before zbrent,','lrest,',lrest,'ebreak,',ebreak,'x,',x,'y,',y,'fp,',fp,'xnew,',xnew,'xnewh,',xnewh,'ynew,',ynew,'yd,',yd,'ifail',ifail

!#ifdef MPI
!if (parallel%isroot) print*,'20190903','e1ord2',e1ord2,'step',step
!#else
!print*,'20190903','e1ord2',e1ord2,'step',step
!#endif
!if (parallel%isroot) print*,'s',s
!if (parallel%isroot) print*,'f',f
      CALL ZBRENT(IU0,LRESET,EBREAK,X,Y,FP,XNEW,XNEWH,YNEW,YD,IFAIL)
!print *,'parameter after zbrent,'
!print *,'lrest,',lrest
!print *,'ebreak,',ebreak
!print *,'x,',x,'y,',y
!print *,'fp,',fp
!print *,'xnew,',xnew
!print *,'xnewh,',xnewh
!print *,'ynew,',ynew
!print *,'yd,',yd
!print *,'ifail',ifail
!     estimate curvature
      CURVET=YD/(E1ORD2/STEP/SQRT(SNORM))**2

      DMOVE =XNEW
      DMOVEH=XNEWH

!print *,'in cg method do 400',ltrial
!    previous step was trial step than give long output
      IF (LTRIAL) THEN
        LBRENT=.TRUE.
        E2ORD  = TOTEN-TOTEN1
        E2ORD2 = (E1ORD1+E1ORD2)/2
#ifdef MPI
        IF (IU0>=0 .and. parallel%isroot) &
#else
        IF (IU0>=0 ) &
#endif
        WRITE(6,45) E2ORD,E2ORD2,E1ORD1,E1ORD2
        ! WRITE(IOUNITS(11),45) E2ORD,E2ORD2,E1ORD1,E1ORD2
   45   FORMAT(' trial-energy change:',F12.6,'  1 .order',3F12.6)

        E1TEST=TOTEN1-YNEW
        DISMAX=DMOVE*DMOVED
#ifdef MPI
        IF (IU0>=0 .and. parallel%isroot) &
#else
        IF (IU0>=0 ) &
#endif
        WRITE(6,20) DMOVE*STEP,DMOVEH*STEP,DISMAX,YNEW,YNEW-TOTEN1
        ! WRITE(IOUNITS(11),20) DMOVE*STEP,DMOVEH*STEP,DISMAX,YNEW,YNEW-TOTEN1
   20   FORMAT(' step: ',F8.4,'(harm=',F8.4,')', &
     &         '  dis=',F8.5,'  next Energy=',F13.6, &
     &         ' (dE=',E10.3,')')

        IF (GAMMA==0) THEN
#ifdef MPI
         IF (IU6>=0 .and. parallel%isroot) &
#else
         IF (IU6>=0 ) &
#endif
          WRITE(6,*)'Steepest descent step on ions:'
          ! WRITE(IOUNITS(11),*)'Steepest descent step on ions:'
        ELSE
#ifdef MPI
        IF (IU6>=0 .and. parallel%isroot) &
#else
        IF (IU6>=0 ) &
#endif
          WRITE(6,*)'Conjugate gradient step on ions:'
          ! WRITE(IOUNITS(11),*)'Conjugate gradient step on ions:'
        ENDIF

#ifdef MPI
        IF (IU6>=0 .and. parallel%isroot) &
#else
        IF (IU6>=0 ) &
#endif
        WRITE(6,40) E2ORD,(E1ORD1+E1ORD2)/2,E1ORD1,E1ORD2, &
     &         GNORM,GNORMF,GNORML, &
     &         GNORM1,GNORM2,ORTH,GAMMA,STEP, &
     &         DMOVE*STEP,DMOVEH*STEP,DISMAX,YNEW,YNEW-TOTEN1
     !    WRITE(IOUNITS(11),40) E2ORD,(E1ORD1+E1ORD2)/2,E1ORD1,E1ORD2, &
     ! &         GNORM,GNORMF,GNORML, &
     ! &         GNORM1,GNORM2,ORTH,GAMMA,STEP, &
     ! &         DMOVE*STEP,DMOVEH*STEP,DISMAX,YNEW,YNEW-TOTEN1

   40   FORMAT(' trial-energy change:',F12.6,'  1 .order',3F12.6/ &
     &       '  (g-gl).g =',E10.3,'      g.g   =',E10.3, &
     &       '  gl.gl    =',E10.3,/ &
     &       ' g(Force)  =',E10.3,'   g(Stress)=',E10.3, &
     &       ' ortho     =',E10.3,/ &
     &       ' gamma     =',F10.5,/ &
     &       ' trial     =',F10.5,/ &
     &       ' opt step  =',F10.5,'  (harmonic =',F10.5,')', &
     &       ' maximal distance =',F10.8/ &
     &       ' next E    =',F13.6,'   (d E  =',F10.5,')')
      ELSE
#ifdef MPI
        IF (IU0>=0 .and. parallel%isroot) &
#else
        IF (IU0>=0 ) &
#endif
        WRITE(6,25) DMOVE*STEP,YNEW,YNEW-TOTEN1
        ! WRITE(IOUNITS(11),25) DMOVE*STEP,YNEW,YNEW-TOTEN1
   25   FORMAT(' opt : ',F8.4,'  next Energy=',F13.6, &
     &         ' (dE=',E10.3,')')
!    d  o not make another call to ZBRENT if reuqired accuracy was reached
        IF (ABS(YD)<EDIFFG) LBRENT=.FALSE.
      ENDIF
!-----------------------------------------------------------------------
!    move ions to the minimum
!-----------------------------------------------------------------------
!print *,'in cg method do 400'
!print *,'posion before move','dmove',dmove,'dmovel',dmovel,'step',step
!print *,'s',s
      DO NI=1,NIONS
!----- search vector from cartesian to direct lattice
        TMP(1) = S(1,NI)
        TMP(2) = S(2,NI)
        TMP(3) = S(3,NI)
        CALL KARDIR(1,TMP,B)

        POSIOC(1,NI)= POSION(1,NI)
        POSIOC(2,NI)= POSION(2,NI)
        POSIOC(3,NI)= POSION(3,NI)

        POSION(1,NI)= TMP(1)*(DMOVE-DMOVEL)*STEP+POSIOC(1,NI)
        POSION(2,NI)= TMP(2)*(DMOVE-DMOVEL)*STEP+POSIOC(2,NI)
        POSION(3,NI)= TMP(3)*(DMOVE-DMOVEL)*STEP+POSIOC(3,NI)

!----- keep ions in unit cell (Robin Hirschl)
        POSION(1,NI)=MOD(POSION(1,NI),1.d0)
        POSION(2,NI)=MOD(POSION(2,NI),1.d0)
        POSION(3,NI)=MOD(POSION(3,NI),1.d0)
        
      ENDDO

    !> FOR confined system ,DO NOT change the cell 
    if(Lpbc)then
    !if(parallel%isroot)then
    !    print *, 'lat_mat vasp a',A*ang2bohr
    !endif
      DO J=1,3
      DO I=1,3
      A(I,J)=AC(I,J)
      DO K=1,3
        A(I,J)=A(I,J) + SSIF(I,K)*AC(K,J)*DMOVE*STEP
      ENDDO
      ENDDO
      ENDDO
    !if(parallel%isroot)then
    !    print *, 'lat_mat vasp b',A*ang2bohr
    !endif
    endif
      DMOVEL=DMOVE
      IFLAG=0
      LTRIAL=.FALSE.

      RETURN
!}}}
      END!}}}



