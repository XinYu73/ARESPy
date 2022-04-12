# 1 "brent.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "brent.f90"
!#include "symbol.inc"
!*********************************************************************
! RCS:  $Id: brent.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
! subroutine ZBRENT
!  variant of Brents root finding method used for
!  minimization of a function
!
!  derivatives and function values have to be supplied
!  the derivatives are considered to be more accurate than the
!  function values !
! ICALL determines behavior and is incremented by the routine
! until a external reset is done
! ICALL
!  0      unit step is done i.e. XNEW=1
!  1      cubic interpolation / extrapolation to minimum
!         hopefully 2. call will give sufficient precision
!  2      init Brents method
!  3      variant of Brendts method
!          bracketing, bisectioning + inverse quadratic interpolation
!          using derivatives only
!
!  on call
!  LRESET    restart algorithm
!  EBREAK    break condition for using the cubic interpolation
!  X,Y,F     position and values for function and derivatives at X
!  on exit
!  XNEW,YNEW position where function has to be evalueated
!            and expected function value for this position
!  XNEWH     harmonic step with calculated for ICALL=1
!  YD        expected change of energy in this step relative to the
!             previous step
!
!*********************************************************************

      SUBROUTINE ZBRENT(IU,LRESET,EBREAK,X,Y,F,XNEW,XNEWH,YNEW,YD,IFAIL)
      use constants

      use smpi_math_module

      ! use parameters , only : iounits

      IMPLICIT REAL(dp) (A-H,O-Z)
      REAL(dp) KOEFF1,KOEFF2,KOEFF3

      LOGICAL LRESET,LBRACK
      SAVE ICALL,LBRACK,A,B,YA,YB,FA,FB,C,FC,YC,FSTART
      SAVE D,E,XM
      PARAMETER (EPS=1E-8_dp,TOL=1E-8_dp)

!     parameter controls the maximum step width
!     DMAX is maximum step for 1. step
!     DMAXBR is maximum step for Brents algorithm, we use golden ratio
      PARAMETER (DMAX=4,DMAXBR=2)
      DATA ICALL/0/
!     minimum not found up to now
      IFAIL=1
      IDUMP=1
      IF (IU<0) IDUMP=0

      IF (LRESET) ICALL=0
!print *,'in zbrent part, icall ',icall
!-----------------------------------------------------------------------
!    case 0: trial step 1, into trial direction
!-----------------------------------------------------------------------
      IF (ICALL==0) THEN!{{{
        A   =X
        YA  =Y
        FA  =F
        YD  =F
        YNEW=Y+F
        XNEW=X+1
        ICALL=1
        IFAIL=1
      RETURN
      ENDIF!}}}
!-----------------------------------------------------------------------
!    case 1:  cubic interpolation / extrapolation
!    if precision of energy is not sufficient use secant method
!-----------------------------------------------------------------------
      IF (ICALL==1) THEN!{{{
        B   =X
        YB  =Y
        FB  =F

        FFIN= YB-YA
        FAB = (FA+FB)/2

        KOEFF3=3.d0*(FB+   FA-2.d0*FFIN)
        KOEFF2=    FB+2.d0*FA-3.d0*FFIN
        KOEFF1=    FA
! if (parallel%isroot) print *,'20190903','dmove',dmove,'koeff3/koeff1',koeff1
!----cubic extrapolation - secant method break condition
        IF ( (ABS(KOEFF3/KOEFF1)< 1E-2_dp)  .OR. &!{{{
     &     (ABS(KOEFF3)       < EBREAK*24).OR. &
     &     ( (KOEFF1*KOEFF3/KOEFF2/KOEFF2) >= 1.d0 )) THEN

!----- harmonic  case
             DMOVE  = -FA/(FAB-FA)/2.d0
             DMOVEH = -FA/(FAB-FA)/2.d0
             IF (DMOVE<0) THEN
               DMOVE =DMAX
               DMOVEH=DMAX
             ENDIF
             DY=(FA+(FAB-FA)*DMOVE)*DMOVE
             YNEW =YA+DY
         ELSE

!----- anharmonic interpolation (3rd order in energy) (jF)
             DMOVE1=KOEFF2/KOEFF3*(1.d0-SQRT(1.d0-KOEFF1*KOEFF3/KOEFF2/KOEFF2))
             DMOVE2=KOEFF2/KOEFF3*(1.d0+SQRT(1.d0-KOEFF1*KOEFF3/KOEFF2/KOEFF2))
!----- 3rd order polynomial has one minimum and one maximum ... :
             DY1=-(KOEFF1-(KOEFF2-KOEFF3*DMOVE1/3.d0)*DMOVE1)*DMOVE1
             DY2=-(KOEFF1-(KOEFF2-KOEFF3*DMOVE2/3.d0)*DMOVE2)*DMOVE2
!----- select the correct extremum:
             IF (DY1>DY2) THEN
               DMOVE=DMOVE1
               DY=DY1
             ELSE
               DMOVE=DMOVE2
               DY=DY2
             ENDIF
             YNEW=YA-DY

             DMOVEH = -FA/(FAB-FA)/2.d0
             IF (DMOVEH<0) DMOVEH=DMAX
!-----  extremely unharmonic and/or
!       minimum should be on the right  side of B but it is due to
!       3. order interpolation on left side
!     -> take the save harmonic extrapolation
             IF (DMOVE>(2*DMOVEH).OR. DMOVEH>(2*DMOVE) &
     &       .OR. (FA*FB>0 .AND.DMOVE<1.d0 ) &
     &       .OR. (FA*FB<0 .AND.DMOVE>1.d0 ))THEN
               DMOVE=DMOVEH
               DY   =(FA+(FAB-FA)*DMOVE)*DMOVE
               YNEW =YA+DY
             ENDIF
        ENDIF!}}}

        IF (DMOVE>DMAX) DMOVE=DMAX

        IF (DMOVE== DMAX .AND. IDUMP>=1 .and. parallel%isroot) &
     &       WRITE(6,*)'ZBRENT: can''t locate minimum, use default step'
        IF (DMOVE== DMAX .AND. IDUMP>=1 .and. parallel%isroot) &
     ! &       WRITE(IOUNITS(11),*) & 
     &       WRITE(6,*) & 
     & 'ZBRENT: can''t locate minimum, use default step'
# 155 "brent.f90"
        XNEW =DMOVE+A
        XNEWH=DMOVEH+A
!     estimated change relative to B
        YD   =(FB+(FAB-FB)/(A-B)*(DMOVE-B))*(DMOVE-B)
        ICALL=2
        IFAIL=1
      RETURN
      ENDIF!}}}
!-----------------------------------------------------------------------
!    case 2:  cubic interpolation / extrapolation failed
!      and was not accurate enough
!      start to use Brents algorithm
!-----------------------------------------------------------------------
      IF (ICALL==2) THEN!{{{
        LBRACK=.TRUE.
!  1.  X > B    start intervall [B,X]
        IF (X>=B) THEN!{{{
            A =B
            YA=YB
            FA=FB
            B =X
            YB=Y
            FB=F
!     check for bracketing
            IF (FA*FB>0) LBRACK=.FALSE.
            FSTART=FA
!  2.  X < B
        ELSE
            IF (FA*F<=0) THEN
!  2a. minimum between [A,X]
                B   =X
                YB  =Y
                FB  =F
!  2b. minimum between [X,B]
            ELSE IF (FB*F<=0) THEN
                A =B
                YA=YB
                FA=FB
                B   =X
                YB  =Y
                FB  =F
!  2c. here we have some serious problems no miniumum between [A,B]
!      but X (search lies between) [A,B] -> complete mess
!      happens only beacuse of cubic interpolations
!      work-around search between [X,B] LBRACK=.FALSE.
             ELSE

                IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*)'ZBRENT:  & 
                                  &  no minimum in in bracket'
                ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) & 
                IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) & 
                     & 'ZBRENT:  no minimum in in bracket'







                A =B
                YA=YB
                FA=FB
                B   =X
                YB  =Y
                FB  =F
                LBRACK=.FALSE.
            ENDIF
        ENDIF!}}}
        ICALL=3
        C =B
        FC=FB
        YC=YB
      !RETURN
      ENDIF!}}}
!-----------------------------------------------------------------------
!  fall trough to this line for ICALL > 4 set: B=X
!-----------------------------------------------------------------------
      IF (ICALL>=4) THEN!{{{
        B =X
        FB=F
        YB=Y
!   maybe a bracketing interval was found ?
        IF (.NOT. LBRACK.AND. FSTART*F < 0) THEN
!   if bracketing is started forget C
          LBRACK=.TRUE.
           C =B
           FC=FB
           YC=YB

           IF (IDUMP/=0 .and. parallel%isroot) WRITE(6,*)'ZBRENT: bracketing found'



           ! IF (IDUMP/=0 .and. parallel%isroot) WRITE(IOUNITS(11),*)'ZBRENT: bracketing found'
        ENDIF
      ENDIF!}}}
!-----------------------------------------------------------------------
! modified brent algorithm if no bracketing intervall exists
! here we have three points [C,A] and B, where B is the new guess
! and lies at the right  hand side of A
! A is the last best guess
!-----------------------------------------------------------------------
      IF (ICALL>=3) THEN!{{{
          ICALL=ICALL+1
          IF (ICALL==20) THEN

             IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*)'ZBRENT:  can not reach accuracy'



          ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*)'ZBRENT:  can not reach accuracy'
          IFAIL=2
          RETURN
        ENDIF


        IF (IDUMP>=2 .and. parallel%isroot) THEN



          WRITE(6,*)
          WRITE(6,'(A5,3E14.7)') 'A',A,YA,FA
          WRITE(6,'(A5,3E14.7)') 'B',B,YB,FB
          WRITE(6,'(A5,3E14.7)') 'C',C,YC,FC
          WRITE(6,*) LBRACK,D,E
          ! WRITE(IOUNITS(11),*)
          ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'A',A,YA,FA
          ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'B',B,YB,FB
          ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'C',C,YC,FC
          ! WRITE(IOUNITS(11),*) LBRACK,D,E
        ENDIF

        IF (.NOT.LBRACK) THEN!{{{
! ABS(FC) < ABS(FB) or AC < BA
            IF(ABS(FC)<=ABS(FB).OR.(A-C)<(B-A)) THEN!{{{
                 C=A
                 FC=FA
                 YC=YA
                 D=B-A
                 E=B-A
            ENDIF!}}}

            TOL1=2.d0*EPS*ABS(B)+0.d0*TOL
            XM=0.5d0*(C-B)


            IF (IDUMP>=2 .and. parallel%isroot) THEN!{{{



                WRITE(6,'(A5,3E14.7)') 'A',A,YA,FA
                WRITE(6,'(A5,3E14.7)') 'B',B,YB,FB
                WRITE(6,'(A5,3E14.7)') 'C',C,YC,FC
                ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'A',A,YA,FA
                ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'B',B,YB,FB
                ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'C',C,YC,FC
            ENDIF!}}}
! just for savety check for correct ordering
            IF ( .NOT.(C<=A .AND. A<=B) ) THEN!{{{

               IF (IU>0 .and. parallel%isroot) THEN



                WRITE(6,*)'ZBRENT: fatal error: bracketing interval incorrect'
                WRITE(6,*)'    please rerun with smaller EDIFF, or copy CONTCAR'
                WRITE(6,*)'    to POSCAR and continue'
                ! WRITE(IOUNITS(11),*)'ZBRENT: fatal error: bracketing interval incorrect'
                ! WRITE(IOUNITS(11),*)'    please rerun with smaller EDIFF, or copy CONTCAR'
                ! WRITE(IOUNITS(11),*)'    to POSCAR and continue'
                ENDIF
                STOP
            ENDIF!}}}
!     return if accuracy ca not be improved
            IF(ABS(XM)<=TOL1 .OR. FB==0.d0)THEN!{{{

                IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT:  accuracy reached'



                ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) 'ZBRENT:  accuracy reached'
                IFAIL=0
                RETURN
            ENDIF!}}}

            IF(ABS(E)>=TOL1 .AND. ABS(FA)>ABS(FB)) THEN  !{{{
                S=FB/FA
!  A.eq.C secant method (only one information)
                IF(A==C) THEN!{{{
                    P=2.d0*XM*S
                    QQ=1.d0-S
!  A.ne.C at        tempt inverse quadratic interpolation
                    ELSE
                    QQ=FA/FC
                    R=FB/FC
                    P=S*(2.d0*XM*QQ*(QQ-R)-(B-A)*(R-1.d0))
                    QQ=(QQ-1.d0)*(R-1.d0)*(S-1.d0)
                ENDIF!}}}
                IF(P>0.d0) QQ=-QQ!{{{
                P=ABS(P)
!  are we within the bounds; correct but tricky :-<

                IF (IDUMP>=2 .and. parallel%isroot) WRITE(6,*)'would go to ',A+P/QQ,P/QQ,E!}}}



                ! IF (IDUMP>=2 .and. parallel%isroot) WRITE(IOUNITS(11),*)'would go to ',A+P/QQ,P/QQ,E
                IF(P < MIN(DMAXBR*(B-A)*QQ-ABS(TOL1*QQ),ABS(E*QQ)/2)) THEN!{{{

                  IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT: extrapolating'



                  ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) 'ZBRENT: extrapolating'
                  E=D
                  D=P/QQ
                ELSE
! interpol      ation increase intervall

                  IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT: increasing intervall'



                  ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) 'ZBRENT: increasing intervall'
                  D=DMAXBR*(B-A)
                  E=D
                ENDIF!}}}
            ELSE
! bounds decrease to slowly increase intervall

               IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT: increasing intervall'!{{{



                ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) 'ZBRENT: increasing intervall'!{{{
                D=DMAXBR*(B-A)
                E=D
                ENDIF!}}}
! estimate new function value using B and A
                FAB = (FA+FB)/2

                YD   =(FB+(FAB-FB)/(A-B)*D)*D
                YNEW =YB+YD
! move A to C and last best guess (B) to A
                C =A
                FC=FA
                YC=YA
                A =B
                FA=FB
                YA=YB
                IF(ABS(D) > TOL1) THEN
                  B=B+D
                ELSE
                  B=B+TOL1
                ENDIF
                XNEW = B
                IFAIL=1
                RETURN
            ELSE
!-----------------------------------------------------------------------
! original brents algorithm take from numberical recipies
! (to me this is absolute mess and not
!  a genius pice of code, but I do not want to change a single line...)
! here we have three points [A,B,C] or [C,B,A] where B is the new guess
!   A is the last best guess
! if the new guess is no improvement or min between [A,B] A is set to C
!    -> secant method is applied
! the new intervall is stored in  [B,C]  B is set to best guess
!-----------------------------------------------------------------------
!  just for savety check whether intervall is correct
                IF ( .NOT.(A<=B .AND. B<=C .OR. C<=B .AND. B<=A)) THEN!{{{

                  IF (IU>0 .and. parallel%isroot) THEN



                  WRITE(6,*)'ZBRENT: fatal error in bracketing'
                  WRITE(6,*)'    please rerun with smaller EDIFF, or copy CONTCAR'
                  WRITE(6,*)'    to POSCAR and continue'
                  ! WRITE(IOUNITS(11),*)'ZBRENT: fatal error in bracketing'
                  ! WRITE(IOUNITS(11),*)'    please rerun with smaller EDIFF, or copy CONTCAR'
                  ! WRITE(IOUNITS(11),*)'    to POSCAR and continue'
                  ENDIF
                 STOP
                ENDIF!}}}

                IF((FB>0 .AND. FC>0).OR. (FB<0 .AND. FC<0))THEN!{{{
                  C=A
                  YC=YA
                  FC=FA
                  D=B-A
                  E=D
                ENDIF!}}}
! C and B are rearanged so that ABS(FC)>=ABS(FB)
                IF(ABS(FC)<ABS(FB)) THEN!{{{
                  A=B
                  B=C
                  C=A
                  FA=FB
                  FB=FC
                  FC=FA
                  YA=YB
                  YB=YC
                  YC=YA
                ENDIF!}}}
                TOL1=2.d0*EPS*ABS(B)+0.5d0*TOL
                XM=0.5d0*(C-B)

            IF (IDUMP>=2 .and. parallel%isroot) THEN!{{{



            WRITE(6,'(A5,3E14.7)') 'A',A,YA,FA
            WRITE(6,'(A5,3E14.7)') 'B',B,YB,FB
            WRITE(6,'(A5,3E14.7)') 'C',C,YC,FC
            ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'A',A,YA,FA
            ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'B',B,YB,FB
            ! WRITE(IOUNITS(11),'(A5,3E14.7)') 'C',C,YC,FC
            ENDIF!}}}
!     return if accuracy is ok
                IF(ABS(XM)<=TOL1 .OR. FB==0.d0)THEN

                IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT:  accuracy reached'



                ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) 'ZBRENT:  accuracy reached'
                IFAIL=0
                RETURN
                ENDIF
                IF(ABS(E)>=TOL1 .AND. ABS(FA)>ABS(FB)) THEN
                S=FB/FA
!  A.eq.C secant method (only one information)
            IF(A==C) THEN!{{{
                P=2.d0*XM*S
                QQ=1.d0-S
!  A.ne.C attempt inverse quadratic interpolation
            ELSE
                QQ=FA/FC
                R=FB/FC
                P=S*(2.d0*XM*QQ*(QQ-R)-(B-A)*(R-1.d0))
                QQ=(QQ-1.d0)*(R-1.d0)*(S-1.d0)
            ENDIF!}}}
                IF(P>0.d0) QQ=-QQ
                P=ABS(P)
!  this is the strangest line but it is ok (trust me)
                IF(2.d0*P < MIN(3.d0*XM*QQ-ABS(TOL1*QQ),ABS(E*QQ))) THEN

                     IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT: interpolating'



                     ! IF (IDUMP>=1 .and. parallel%isroot) WRITE(IOUNITS(11),*) 'ZBRENT: interpolating'
                     E=D
                     D=P/QQ
                ELSE
! interpolation failed bisectioning

                    IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT: bisectioning'



                    ! IF (IDUMP>=1 .and. parallel%isroot) & 
                    !             & WRITE(IOUNITS(11),*) 'ZBRENT: bisectioning'
                    D=XM
                    E=D
                ENDIF
            ELSE
! bounds decrease to slowly bisectioning

                 IF (IDUMP>=1 .and. parallel%isroot) WRITE(6,*) 'ZBRENT: bisectioning'



                 ! IF (IDUMP>=1 .and. parallel%isroot) & 
                 !                 & WRITE(IOUNITS(11),*) 'ZBRENT: bisectioning'
                 D=XM
                 E=D
            ENDIF!}}}
! estimate new function value using B and A
        FAB = (FA+FB)/2

        YD   =(FB+(FAB-FB)/(A-B)*D)*D
        YNEW =YB+YD
! move last best guess to A
        A=B
        FA=FB
        YA=YB
        IF(ABS(D) > TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
         XNEW = B
         IFAIL=1
         RETURN
        ENDIF  !}}}
      ENDIF!}}}

      END
