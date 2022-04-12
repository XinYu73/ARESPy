# 1 "Relax_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Relax_module.f90"
!#############################################################!
!---            Structual Relaxation  Module               ---!
!---Author:              Qiang Xu                          ---!
!---Date:               2018/01/17                         ---!
!#############################################################!
MODULE relax_module
  USE constants
  USE parameters , ONLY : IOPM,IGOAL,TIMES,PRESS,TOLF,TOLP
  IMPLICIT NONE
  REAL(DP) :: pstress
  LOGICAL  :: Ldone  &
           &, Lfirst
  TYPE relax_type
       !FIRE data structural
       REAL(DP)                 :: TIMSmax &
                               &,  TIMS    &
                               &,  alpha
       REAL(DP),ALLOCATABLE     :: FACTion(:)  &
                               &,  Fiond(:,:)  &
                               &,  Fceld(:,:)
       REAL(DP),ALLOCATABLE     :: VELion(:,:)
       REAL(DP),ALLOCATABLE     :: VELcel(:,:)
       REAL(DP)                 :: FACTcel
       INTEGER(I4B)             :: NEG
       LOGICAL                  :: lNEG
  ENDTYPE relax_type
  !
  TYPE(relax_type) :: relax
CONTAINS
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !------------------initialize relaxation--------------------
  SUBROUTINE initialize_relax()
     USE struct_module , ONLY : natom,naty,struct
     USE m_time_evaluate, ONLY : memory_sum
     IMPLICIT NONE
     REAL(DP) :: totMass
     INTEGER(I4B) :: Ity,Ia
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL destroy_relax()
     !
     Ldone=.FALSE.
     Lfirst=.TRUE.
     IF(IGOAL==3)THEN
        !FIRE data
        ALLOCATE(relax%FACTion(natom))
        ALLOCATE(relax%Fiond(3,natom))
        ALLOCATE(relax%Fceld(3,3))
        ALLOCATE(relax%VELion(3,natom))
        ALLOCATE(relax%VELcel(3,3))
        call memory_sum('relax_Fire',(7*real(natom,DP)+18)*DP)
        !initial data
        relax%NEG=0
        relax%lNEG=.TRUE.
        relax%TIMSmax=TIMES*10.d0
        relax%TIMS=TIMES
        relax%alpha=1.d0
        pstress=PRESS/au2gpa
        !initial vel
        relax%VELion(:,:)=0.d0
        relax%VELcel(:,:)=0.d0
        !initial FACT
        totMass=0.d0
        DO Ity=1,naty
           DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
              relax%FACTion(Ia)=CONST_MA_AU/struct%mass(Ity)
              totMASS=totMASS+struct%mass(Ity)/CONST_MA_AU
           ENDDO
        ENDDO

        IF(IGOAL==2) relax%FACTion(:)=0.d0
        IF(IGOAL==1)THEN
            relax%FACTcel=0.d0
        ELSE
            relax%FACTcel=(natom/totMASS)
        ENDIF

     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE initialize_relax
  !-------------------destroy relaxation----------------------
  SUBROUTINE destroy_relax()
     USE m_time_evaluate, ONLY: memory_free
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(relax%FACTion))THEN
        call memory_free('relax_fire%FACTion',real(size(relax%FACTion),DP)*DP)
        DEALLOCATE(relax%FACTion)
     ENDIF
     IF(ALLOCATED(relax%Fiond))THEN
        call memory_free('relax_fire%Fiond',real(size(relax%Fiond),DP)*DP)
        DEALLOCATE(relax%Fiond)
     ENDIF
     IF(ALLOCATED(relax%Fceld))THEN
        call memory_free('relax_fire%Fceld',real(size(relax%Fceld),DP)*DP)
        DEALLOCATE(relax%Fceld)
     ENDIF
     IF(ALLOCATED(relax%VELion))THEN
        call memory_free('relax_fire%VELion',real(size(relax%VELion),DP)*DP)
        DEALLOCATE(relax%VELion)
     ENDIF
     IF(ALLOCATED(relax%VELcel))THEN
        call memory_free('relax_fire%VELcel',real(size(relax%VELcel),DP)*DP)
        DEALLOCATE(relax%VELcel)
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_relax
  !-------------------------relaxer---------------------------
  SUBROUTINE relaxer(nstep)
     USE parameters ,     ONLY : Lpbc     !##ADDED FOR CONFINED SYSTEM
     USE struct_module , ONLY : struct,lat_mat
     !> add for CG
     USE cg_relax , ONLY : cg_relax_vasp_interface


     IMPLICIT NONE
     REAL(DP) :: Fcel(3,3)
     !> add for CG
     INTEGER(I4B),intent(in) :: nstep
     logical :: lfopt=.false.,lopt,lexceed
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(Lpbc)THEN
       CALL fix_direction(struct%stress)
       CALL force_on_cell(struct%stress,lat_mat,Fcel)
     ELSE
       Fcel=0.d0
     ENDIF
     !
     IF(IOPM==3)THEN
        CALL FIRE_relax(lat_mat,struct%pos,struct%forces/2.d0,Fcel)
     ELSEIF(IOPM==4)THEN
        CALL cg_relax_vasp_interface(nstep,lfopt,lopt,lexceed)
        if(lfopt)then
           Ldone=.true.
        endif


     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE relaxer
  !---------------------FIRE relaxation-----------------------
  SUBROUTINE FIRE_relax(LAT,POS,Fion,Fcel)
     USE struct_module , ONLY : natom
     USE math , ONLY : car2dir

     USE smpi_math_module, ONLY:parallel

     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(INOUT) :: LAT(:,:) &
                         &, POS(:,:)
     REAL(DP),INTENT(IN) :: Fion(:,:)  &  !force on ion
                         &, Fcel(:,:)     !force on cell
     !LOCAL
     REAL(DP) :: Norm2Fion,Norm2Fcel,Norm2Vion,Norm2Vcel &
             & , NORMF,NORMV &
             & , POWERion,POWERcel
     REAL(DP) :: VELdir(3,natom)
     INTEGER(I4B) :: Ia
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     Norm2Fion=SUM(Fion*Fion)
     Norm2Fcel=SUM(Fcel*Fcel)
     !========================================================
     IF((SQRT(Norm2Fion/(3*natom))*au2ev_force<TOLF.OR.IGOAL==2).AND. &
      & (SQRT(Norm2Fcel/9)*au2gpa<TOLP.OR.IGOAL==1) )THEN
         Ldone=.TRUE.
         RETURN
     ENDIF
     !
     IF(relax%TIMS==0.d0)THEN
          WRITE(6,*)'FIRE:stop!!!please change times step and try again'
          STOP
     ENDIF
     !========================================================
     IF(Lfirst)THEN
        relax%Fiond(:,:)=Fion(:,:)
        relax%Fceld(:,:)=Fcel(:,:)
        Lfirst=.FALSE.
     ENDIF
     !---MD integrator----------------------------------------
     !velocity
     DO Ia=1,natom
        relax%VELion(:,Ia)=relax%VELion(:,Ia) &
                 &  +0.5d0*relax%TIMS*relax%Fiond(:,Ia)*relax%FACTion(Ia)
     ENDDO
     relax%VELcel(:,:)=relax%VELcel(:,:) + &
                &   0.5d0*relax%TIMS*relax%Fceld(:,:)*relax%FACTcel
     !move the position
     !1.move ion
     CALL car2dir(relax%VELion,VELdir,LAT)
     POS(:,:)=POS(:,:)+relax%TIMS*VELdir(:,:)
     !2.move cell
     LAT(:,:)=LAT(:,:)+relax%TIMS*relax%VELcel(:,:)
     !velocity
     DO Ia=1,natom
        relax%VELion(:,Ia)=relax%VELion(:,Ia) &
                 &  +0.5d0*relax%TIMS*Fion(:,Ia)*relax%FACTion(Ia)
     ENDDO
     relax%VELcel(:,:)=relax%VELcel(:,:) + &
                &   0.5d0*relax%TIMS*Fcel(:,:)*relax%FACTcel
     !store Force
     relax%Fiond(:,:)=Fion(:,:)
     relax%Fceld(:,:)=Fcel(:,:)
     !correct the velocity------------------------------------
     !F1
     !--------------------------------------------------------
     !power=INPRODUCT(Fion,VEKion)+INPRODUCT(Fcel,VELcel)
     IF(IGOAL/=2)THEN
        POWERion=SUM(Fion*relax%VELion)
     ELSE
        POWERion=0.d0
     ENDIF
     !-----------------
     IF(IGOAL/=1)THEN
        POWERcel=SUM(Fcel*relax%VELcel)
     ELSE
        POWERcel=0.d0
     ENDIF
     !------------------
     IF(IGOAL/=2)THEN
         NORM2Vion=SUM(relax%VELion*relax%VELion)
     ELSE
         NORM2Fion=0.d0
         NORM2Vion=0.d0
     ENDIF
     !-------------------
     IF(IGOAL/=1)THEN
         NORM2Vcel=SUM(relax%VELcel*relax%VELcel)
     ELSE
         NORM2Fcel=0.d0
         NORM2Vcel=0.d0
     ENDIF
     !--------------------
     NORMF=SQRT(NORM2Fion+NORM2Fcel)
     NORMV=SQRT(NORM2Vion+NORM2Vcel)

     IF(IGOAL/=2) relax%VELion(:,:)=(1-relax%alpha)*relax%VELion(:,:) &
                      &   +relax%alpha*Fion(:,:)/NORMF*NORMV
     IF(IGOAL/=1) relax%VELcel(:,:)=(1-relax%alpha)*relax%VELcel(:,:) &
                      &   +relax%alpha*Fcel(:,:)/NORMF*NORMV
     !--F3--F4
     !-----------------
     IF((POWERion<0.d0.AND.IGOAL/=2).OR.(POWERcel<0.d0.AND.IGOAL/=1))THEN
         relax%TIMS=relax%TIMS*0.5
         relax%alpha=0.1
         !set 0
         relax%VELion(:,:)=0.d0
         relax%VELcel(:,:)=0.d0
         !
         relax%NEG=0
         relax%lNEG=.TRUE.
     ELSEIF(relax%NEG>5)THEN
         relax%TIMS=MIN(relax%TIMS*1.1,relax%TIMSmax)
         relax%alpha=relax%alpha*0.99
         relax%lNEG=.FALSE.
     ENDIF

     IF(relax%lNEG) relax%NEG=relax%NEG+1
     !--------------------------------------------------------

     if(parallel%isroot)then

     WRITE(6,*) '**************FIRE RELAX****************'
     WRITE(6,*) 'root-mean-squre:','force (eV/ang):',SQRT(NORM2Fion/(3.d0*natom))*au2ev_force
     WRITE(6,*) 'root-mean-squre:','stress (GPa):',SQRT(NORM2Fcel/9.d0)*au2gpa
     WRITE(6,*) 'present times step:',relax%TIMS
     WRITE(6,*) 'present power (ion,cel):','(',POWERion,',',POWERcel,')'
     WRITE(6,*) '**************FIRE RELAX****************'

     endif

     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE FIRE_relax
  !-------------------------FIRE XG---------------------------
  SUBROUTINE force_on_cell(stress,LAT,Fcel)
     USE math , ONLY : inv_33
     USE struct_module , ONLY : volume
     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(IN)  :: stress(:,:),LAT(:,:)
     REAL(DP),INTENT(OUT) :: Fcel(:,:)
     !LOCAL
     REAL(DP) :: as(3,3)
     INTEGER(I4B) :: i
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !force on cell
     as(:,:)=-stress(:,:)
     DO i=1,3
        as(i,i)=as(i,i)-pstress
     ENDDO
     Fcel(:,:)=volume*MATMUL(as,TRANSPOSE(inv_33(LAT)))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE force_on_cell
  !------------fix some direction in relaxing-----------------
  SUBROUTINE fix_direction(stress)
    USE parameters ,ONLY:fix_xyz
    IMPLICIT NONE
    !> in/out
    REAL(DP) :: stress(:,:)

    if(fix_xyz(1)/=0)then
       stress(1,1)=0.d0
    elseif(fix_xyz(2)/=0)then
       stress(2,2)=0.d0
    elseif(fix_xyz(3)/=0)then
       stress(3,3)=0.d0
    elseif(fix_xyz(4)/=0)then
       stress(2,1)=0.d0
       stress(1,2)=0.d0
    elseif(fix_xyz(5)/=0)then
       stress(3,1)=0.d0
       stress(1,3)=0.d0
    elseif(fix_xyz(6)/=0)then
       stress(2,3)=0.d0
       stress(3,2)=0.d0
    endif
  ENDSUBROUTINE fix_direction
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE relax_module
