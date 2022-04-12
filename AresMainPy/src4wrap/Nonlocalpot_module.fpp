# 1 "Nonlocalpot_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Nonlocalpot_module.f90"
MODULE nlpot_module
  !######################################################!
  !*For    : nonlocal pseudo-potential data and function !
  !*Author : Qiang Xu                                    !
  !*Date   : 2017-12-20                                  !
  !######################################################!
  USE constants
  IMPLICIT NONE
  !defined the type
  TYPE nol_type
     INTEGER(I4B) :: npts
     INTEGER(I4B),ALLOCATABLE :: Id(:)   !index for atom nonlpots
     REAL(DP),ALLOCATABLE :: rRvec(:,:)  !vector for r-ra
     REAL(DP),ALLOCATABLE :: proj0(:,:)  !|dv\phi>
     REAL(DP),ALLOCATABLE :: proj(:,:)   !|dv\phi_lm>
     COMPLEX(DP),ALLOCATABLE :: proj_phs(:,:,:)  !e^-ikr|dv\phi_lm>
     !Double-Grid
     REAL(DP),ALLOCATABLE :: proj0_dg(:,:)  !|dv\phi>
     REAL(DP),ALLOCATABLE :: proj_dg(:,:)   !|dv\phi_lm>
     COMPLEX(DP),ALLOCATABLE :: proj_phs_dg(:,:,:)  !e^-ikr|dv\phi_lm>
     !ISO_id
     INTEGER(I4B),ALLOCATABLE :: Id_iso(:)   !index for atom nonlpots in ISO sphere
     REAL(DP),ALLOCATABLE :: rRvec_iso(:,:)  !vector for r-ra in ISO sphere
  ENDTYPE nol_type
  !
  TYPE(nol_type),ALLOCATABLE :: nlpot(:)
  !
  INTEGER(I4B) :: max_nlnpts !max points
  !Structural factor
  !COMPLEX(DP),ALLOCATABLE :: nlpot_Sfact(:)
CONTAINS
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  SUBROUTINE initialize_nlpot ()
    USE parameters, ONLY: Lpbc
    IMPLICIT NONE

    if(Lpbc)then
       call initialize_nlpot_per()
    else
       call initialize_nlpot_iso()
    endif

  END SUBROUTINE initialize_nlpot

  !------------------- part line ---------------------
  SUBROUTINE initialize_nlpot_per()
     USE parameters , ONLY : LDG,LbvK,LGamma
     USE pspot_module , ONLY : max_nproj,max_rcut,psp
     USE struct_module , ONLY : struct,natom,lat_mat,naty
     USE grid_module ,ONLY : gap,n,nk,grid , n1,n2,n3
     USE math , ONLY : mat2thr
      USE m_time_evaluate, ONLY: memory_sum,memory_free

      USE smpi_math_module, ONLY:parallel

      IMPLICIT NONE
      !LOCAL
      REAL(DP) :: rcut2
      REAL(DP) :: ra(3),dR(3),dist2
      !!
      INTEGER(I4B) :: nrepl
      INTEGER(I4B) ::  nlnpts(natom)  !num of points
      !!index
      INTEGER(I4B) :: Ia,m,Ity  &
                &  ,  ip        &
                &  ,icx,icy,icz,ix,iy,iz
      REAL(DP)     :: rat(3)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(max_nproj==0)THEN

         if(parallel%isroot)WRITE(6,*) 'Local pseudopotential has been used'



         RETURN
      ENDIF


      if(parallel%isroot)WRITE(6,*) 'Non-Local pseudopotential has been used'



      !
      nrepl=1
      rcut2=max_rcut**2
      !inqure the number of points we need
      !tes=0.d0
      nlnpts=0
      DO Ia=1,natom
         ra(:)=struct%poscar(:,Ia)
         m=0
         DO ip=1,n
            !
            DO icz=-nrepl,nrepl
               DO icy=-nrepl,nrepl
                  DO icx=-nrepl,nrepl
                     !ret=rcell + rvec - ra
                     rat(:)= icx*lat_mat(:,1)   &
                          & +  icy*lat_mat(:,2)   &
                          & +  icz*lat_mat(:,3)   &
                          & +   grid%rVec(1:3,ip) &
                          & -  ra(:)
                     dist2=rat(1)**2+rat(2)**2+rat(3)**2
                     !
                     IF(dist2 <= rcut2) m=m+1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         nlnpts(Ia)=m
      ENDDO
      !===================store the data================
      max_nlnpts=MAXVAL(nlnpts(:))

      !
      CALL destroy_nlpot()
      !for projectors
      ALLOCATE(nlpot(natom))
      DO Ia=1,natom
         !>>>
         ALLOCATE(nlpot(Ia)%Id(nlnpts(Ia)))
         ALLOCATE(nlpot(Ia)%rRvec(4,nlnpts(Ia)))
         call memory_sum('nlpot_Id,rRvec',real(size(nlpot(Ia)%Id),DP)*I4B&
              &+size(nlpot(Ia)%rRvec)*DP)
         !store
         nlpot(Ia)%npts=nlnpts(Ia)
      ENDDO
      !projector
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            IF(.not.LDG)THEN
               ALLOCATE(nlpot(Ia)%proj0(nlnpts(Ia),psp(Ity)%nproj))
               ALLOCATE(nlpot(Ia)%proj(nlnpts(Ia),psp(Ity)%nproj))
               call memory_sum('nlpot_proj0,proj',&
                    & real(size(nlpot(Ia)%proj0),DP)*DP&
                    &+size(nlpot(Ia)%proj)*DP)
               IF(.NOT.(LGamma.and.LbvK))THEN
                  ALLOCATE(nlpot(Ia)%proj_phs(nlnpts(Ia),psp(Ity)%nproj,nk))
                  call memory_sum('nlpot_proj_phs',&
                       & real(size(nlpot(Ia)%proj_phs),DP)*DP*2)
               ENDIF
               !SF
               !ALLOCATE(nlpot(Ia)%Sfact(nk))
            ELSE
               !DG
               ALLOCATE(nlpot(Ia)%proj0_dg(nlnpts(Ia),psp(Ity)%nproj))
               ALLOCATE(nlpot(Ia)%proj_dg(nlnpts(Ia),psp(Ity)%nproj))
               call memory_sum('nlpot_proj0_dg,proj_dg',&
                    & real(size(nlpot(Ia)%proj0_dg),DP)*DP&
                    &+size(nlpot(Ia)%proj_dg)*DP)
               IF(.NOT.(LGamma.and.LbvK))THEN
                  ALLOCATE(nlpot(Ia)%proj_phs_dg(nlnpts(Ia),psp(Ity)%nproj,nk))
                  call memory_sum('nlpot_proj_phs_dg',&
                       & real(size(nlpot(Ia)%proj_phs_dg),DP)*DP*2)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !===================store data====================
      DO Ia=1,natom
         ra(:)=struct%poscar(:,Ia)
         m=0
         DO ip=1,n
            !
            DO icz=-nrepl,nrepl
               DO icy=-nrepl,nrepl
                  DO icx=-nrepl,nrepl
                     !ret=rcell + rvec - ra
                     rat(:)= icx*lat_mat(:,1)   &
                          & +  icy*lat_mat(:,2)   &
                          & +  icz*lat_mat(:,3)   &
                          & +   grid%rVec(1:3,ip) &
                          & -  ra(:)
                     dist2=rat(1)**2+rat(2)**2+rat(3)**2
                     !
                     IF(dist2 <= rcut2) THEN
                        m=m+1
                        nlpot(Ia)%Id(m)=ip
                        nlpot(Ia)%rRvec(1:3,m)=rat(1:3)
                        nlpot(Ia)%rRvec(4,m)=SQRT(dist2)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !check
         IF(m/=nlpot(Ia)%npts)THEN
            WRITE(6,*) 'Error: Initialize nolpot',m,nlpot(Ia)%npts
            STOP
         ENDIF

      ENDDO
      !CALL nlp_struct_fact(nlpot_Sfact)
      !apply Ylm phase
      IF(LDG)THEN
         !Use time-saving method
         CALL set_beta_real_dg()
      ELSE
         !Use uniform grid
         CALL set_beta_real()
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE initialize_nlpot_per
    !------------------- initialize_nlpot---------------------
    SUBROUTINE initialize_nlpot_iso()
      USE parameters , ONLY : LDG,LBvK
      USE pspot_module , ONLY : max_nproj,max_rcut,psp
      USE struct_module , ONLY : struct,natom,lat_mat,naty
      USE grid_module ,ONLY : gap,n,nk,grid , n1,n2,n3 , rho_calc
      USE math , ONLY : mat2thr

      USE smpi_math_module , ONLY : atom_split,parallel,smpi_reduce_sum_real_2d,smpi_reduce_sum_int_1d

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !LOCAL
      REAL(DP) :: rcut2
      REAL(DP) :: ra(3),dR(3),dist2
      !!
      INTEGER(I4B) :: nrepl
      INTEGER(I4B) ::  nlnpts(natom)  !num of points
      !!index
      INTEGER(I4B) :: Ia,m,Ity  &
           &  ,  ip        &
           &  ,icx,icy,icz,ix,iy,iz
      REAL(DP)     :: rat(3)

      integer(I4B) :: mysize,atom_id=0,id_core
      INTEGER(I4B),allocatable :: atom_index(:)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(max_nproj==0)THEN



         if(parallel%isroot)WRITE(6,*) 'Local pseudopotential has been used'

         RETURN
      ENDIF




      if(parallel%isroot)WRITE(6,*) 'Non-Local pseudopotential has been used'

      !
      nrepl=0
      rcut2=max_rcut**2
      !inqure the number of points we need
      !tes=0.d0
      nlnpts=0
      DO Ia=1,natom
         ra(:)=struct%poscar(:,Ia)
         m=0



            DO ip=1,parallel%mygrid_range(3)

               !
               DO icz=-nrepl,nrepl
                  DO icy=-nrepl,nrepl
                     DO icx=-nrepl,nrepl
                        !ret=rcell + rvec - ra
                        rat(:)= icx*lat_mat(:,1)   &
                             & +  icy*lat_mat(:,2)   &
                             & +  icz*lat_mat(:,3)   &



                             & +   grid%rVec(1:3,rho_calc%n(ip)) &

                             & -  ra(:)
                        dist2=rat(1)**2+rat(2)**2+rat(3)**2
                        !
                        IF(dist2 <= rcut2) m=m+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            nlnpts(Ia)=m
         ENDDO
         !===================store the data================
         max_nlnpts=MAXVAL(nlnpts(:))

         !
         CALL destroy_nlpot()
         !for projectors
         ALLOCATE(nlpot(natom))
         DO Ia=1,natom
            !>>>
            ALLOCATE(nlpot(Ia)%Id_iso(nlnpts(Ia)))
            ALLOCATE(nlpot(Ia)%rRvec_iso(4,nlnpts(Ia)))
            call memory_sum('nlpot_Id_iso,rRvec_iso',&
                 & real(size(nlpot(Ia)%Id_iso),DP)*I4B&
                 &+size(nlpot(Ia)%rRvec_iso)*DP)
            !store
            nlpot(Ia)%npts=nlnpts(Ia)
         ENDDO
         !projector
         DO Ity=1,naty
            IF(psp(Ity)%nproj==0) CYCLE
            DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
               IF(.not.LDG)THEN
                  ALLOCATE(nlpot(Ia)%proj0(nlnpts(Ia),psp(Ity)%nproj))
                  ALLOCATE(nlpot(Ia)%proj(nlnpts(Ia),psp(Ity)%nproj))
                  call memory_sum('nlpot_proj0,proj',&
                       & real(size(nlpot(Ia)%proj0),DP)*DP&
                       &+size(nlpot(Ia)%proj)*DP)
                  !SF
                  !ALLOCATE(nlpot(Ia)%Sfact(nk))
               ELSE
                  !DG
                  ALLOCATE(nlpot(Ia)%proj0_dg(nlnpts(Ia),psp(Ity)%nproj))
                  ALLOCATE(nlpot(Ia)%proj_dg(nlnpts(Ia),psp(Ity)%nproj))
                  call memory_sum('nlpot_proj0_dg,proj_dg',&
                       & real(size(nlpot(Ia)%proj0_dg),DP)*DP&
                       &+size(nlpot(Ia)%proj_dg)*DP)
               ENDIF
            ENDDO
         ENDDO
         !===================store data====================
         !##USED IN ISO FORCE CALCULATE,scf could use too {{
         DO Ia=1,natom

            nlpot(Ia)%rRvec_iso=0
            nlpot(Ia)%Id_iso=0

            ra(:)=struct%poscar(:,Ia)
            m=0



               DO ip=1,parallel%mygrid_range(3)

                  !
                  !ret=rcell + rvec - ra
                  rat(1)=rho_calc%x(ip)*gap(1)-ra(1)
                  rat(2)=rho_calc%y(ip)*gap(2)-ra(2)
                  rat(3)=rho_calc%z(ip)*gap(3)-ra(3)
                  dist2=rat(1)**2+rat(2)**2+rat(3)**2
                  !
                  IF(dist2 <= rcut2) THEN
                     m=m+1
                     nlpot(Ia)%Id_iso(m)=ip
                     nlpot(Ia)%rRvec_iso(1:3,m)=rat(1:3)
                     nlpot(Ia)%rRvec_iso(4,m)=SQRT(dist2)
                  ENDIF
               ENDDO
               !check
               IF(m/=nlpot(Ia)%npts)THEN
                  WRITE(6,*) 'Error: Initialize nolpot:may increase the radius of sphere would ok',m,nlpot(Ia)%npts
                  STOP
               ENDIF
            ENDDO
            !CALL nlp_struct_fact(nlpot_Sfact)
            !apply Ylm phase
            IF(LDG)THEN
               !Use time-saving method
               CALL set_beta_real_dg()
            ELSE
               !Use uniform grid
               CALL set_beta_real()
            ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE initialize_nlpot_iso
   !-----------------destroy array----------------------
   SUBROUTINE destroy_nlpot()
      USE struct_module , ONLY : natom
      USE m_time_evaluate, ONLY: memory_free
      IMPLICIT NONE
      INTEGER(I4B) :: Ia
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      IF(ALLOCATED(nlpot))THEN
         DO Ia=1,natom
            IF(ALLOCATED(nlpot(Ia)%Id))THEN
               call memory_free('nlpot_Id',real(size(nlpot(Ia)%Id),DP)*I4B)
               DEALLOCATE(nlpot(Ia)%Id)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%rRvec))THEN
               call memory_free('nlpot_rRvec',real(size(nlpot(Ia)%rRvec),DP)*DP)
               DEALLOCATE(nlpot(Ia)%rRvec)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%Id_iso))THEN
               call memory_free('nlpot_Id_iso',real(size(nlpot(Ia)%Id_iso),DP)*I4B)
               DEALLOCATE(nlpot(Ia)%Id_iso)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%rRvec_iso))THEN
               call memory_free('nlpot_rRvec_iso',real(size(nlpot(Ia)%rRvec_iso),DP)*DP)
               DEALLOCATE(nlpot(Ia)%rRvec_iso)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%proj0))THEN
               call memory_free('nlpot_proj0',real(size(nlpot(Ia)%proj0),DP)*DP)
               DEALLOCATE(nlpot(Ia)%proj0)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%proj))THEN
               call memory_free('nlpot_proj',real(size(nlpot(Ia)%proj),DP)*DP)
               DEALLOCATE(nlpot(Ia)%proj)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%proj_phs))THEN
               call memory_free('nlpot_proj_phs',real(size(nlpot(Ia)%proj_phs),DP)*DP)
               DEALLOCATE(nlpot(Ia)%proj_phs)
            ENDIF
               !DG
            IF(ALLOCATED(nlpot(Ia)%proj0_dg))THEN
               call memory_free('nlpot_proj0_dg',real(size(nlpot(Ia)%proj0_dg),DP)*DP)
               DEALLOCATE(nlpot(Ia)%proj0_dg)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%proj_dg))THEN
               call memory_free('nlpot_proj_dg',real(size(nlpot(Ia)%proj_dg),DP)*DP)
               DEALLOCATE(nlpot(Ia)%proj_dg)
            ENDIF
            IF(ALLOCATED(nlpot(Ia)%proj_phs_dg))THEN
               call memory_free('nlpot_proj_phs_dg',real(size(nlpot(Ia)%proj_phs_dg),DP)*DP)
               DEALLOCATE(nlpot(Ia)%proj_phs_dg)
            ENDIF
         ENDDO
      ENDIF
      !
      IF(ALLOCATED(nlpot))    DEALLOCATE(nlpot)
      !destroy structure factor
      !IF(ALLOCATED(nlpot_Sfact))    DEALLOCATE(nlpot_Sfact)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_nlpot
   !-------------------build Project to real space------------------
   SUBROUTINE set_beta_real()
      !for parpare in real space  |e^(-ikr)\beta_lm>
      USE parameters , ONLY : LBvK,Lpbc
      USE grid_module , ONLY : dvol   ,n1,n2,n3
      USE struct_module , ONLY : struct,naty,natom
      USE pspot_module , ONLY : psp,max_nproj
      USE math , ONLY : mat2thr

      USE smpi_math_module ! , ONLY : atom_split,parallel,smpi_reduce_sum_real_2d,smpi_reduce_sum_int_1d

      USE array_io
      IMPLICIT NONE
      !LOCAL
      INTEGER(I4B) :: Ity,Ia

      integer(I4B) :: mysize,atom_id=0,id_core
      INTEGER(I4B),allocatable :: atom_index(:)

      !======================
      !##PLOT
      INTEGER(I4B) :: Ipj,l,m
      !> debug
      CHARACTER(len=1) :: str1,str2
      INTEGER(I4B)     :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(max_nproj==0) RETURN
      !============apply Ylm in real space===============
! #ifdef 1
!       !> initial parallel config
!       ALLOCATE(atom_index(natom/parallel%numprocs+1))
!       CALL atom_split(mysize,atom_index)
!       id_core=1
! #endif
      !cycle the all species
      DO Ity=1,naty
         !
         IF(psp(Ity)%nproj==0) CYCLE
         !cycle all atom in this species
         DO Ia = struct%eleid(Ity),struct%eleid(Ity+1)-1
! #ifdef 1
!          If(Ia==atom_index(id_core))THEN
!             if(id_core<mysize)id_core=id_core+1
! #endif
            !interpolates data to mesh
            CALL nlp_beta_interp_r(Ity,Ia,nlpot(Ia)%proj0)
            !apply_ylm
            nlpot(Ia)%proj(:,:)=nlpot(Ia)%proj0(:,:)
            CALL nlp_beta_ylm_r(Ity,Ia,nlpot(Ia)%proj)
            IF((.NOT.LBvK).AND.Lpbc)THEN
               !apply phase
               CALL nlp_beta_phase_r(Ity,Ia,nlpot(Ia)%proj,nlpot(Ia)%proj_phs)
            ENDIF
! #ifdef 1
!          else
!             nlpot(Ia)%proj=0
!             nlpot(Ia)%proj0=0
!          endif
! #endif
         ENDDO
      ENDDO
! #ifdef 1
!       DO Ia=1,natom
!          CALL smpi_reduce_sum_real_2d(nlpot(Ia)%proj)
!          CALL smpi_reduce_sum_real_2d(nlpot(Ia)%proj0)
!       ENDDO
! #endif


# 510 "Nonlocalpot_module.f90"
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE set_beta_real
   !-----------------------interpolate data----------------------------
   SUBROUTINE nlp_beta_interp_r(Ity,Ia,beta_init)
      !USE math , ONLY : CubicSplineInterp
      USE math , ONLY : CubicHermiteInterp,interp
      USE pspot_module , ONLY : psp
      USE MathSplines  , ONLY : polynom
      USE parameters   , ONLY : Lpbc,PP_identifer
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia
      REAL(DP),INTENT(OUT) :: beta_init(:,:)
      !LOCAL
      INTEGER(I4B) :: Ipj,Ip,Id  &
                  &,  l,m,pos
      REAL(DP) :: rRnorm,betarR
      REAL(DP) :: c(3)
      INTEGER(I4B) :: max_m
      c=0.d0
      DO Ipj=1,psp(Ity)%nproj
         !(l,m)
         l=psp(Ity)%proj_l(Ipj)
         m=psp(Ity)%proj_m(Ipj)
         !> cycle all point os pot
         DO Ip=1,nlpot(Ia)%npts
            !index
            IF(Lpbc)THEN
            Id=nlpot(Ia)%Id(Ip)
            !
            rRnorm=nlpot(Ia)%rRvec(4,Ip)
            !interpole the value of nonlocal pseudopotential
            if(PP_identifer==0)then
              betarR=CubicHermiteInterp(psp(Ity)%beta_r(:,Ipj),psp(Ity)%dbeta_dr(:,Ipj) &
                   &    ,psp(Ity)%rmax, psp(Ity)%rspacing , rRnorm)
              else
               ! pos=0
               ! DO WHILE ( psp(Ity)%r_real(pos+1).lt.rRnorm .and. (pos.le.psp(Ity)%numps-2) )
               !    pos=pos+1
               ! ENDDO
               ! !##CALCULATE POTENTIAL IN "RR"
               !   betarR=polynom(0,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,rRnorm)
                 betarR=interp(psp(Ity)%numps&
                      &,psp(Ity)%beta_r(:,Ipj),psp(Ity)%r_real&
                      &,rRnorm)
                 if(isnan(betarR))print*,'betarR err',betarR,'r,',psp(Ity)%r_real(pos-1:pos+1),'f',psp(Ity)%beta_r(pos-1:pos+1,Ipj)
              ! betarR=CubicHermiteInterp(psp(Ity)%beta_r(:,Ipj),psp(Ity)%dbeta_dr(:,Ipj) &
              !      &    ,psp(Ity)%rmax, psp(Ity)%rspacing , rRnorm)
              endif
            ELSE
               !##INTERPOLATING BY THREE POINT
               !##FIND THE THREE POINT
               Id=nlpot(Ia)%Id_iso(Ip)
               rRnorm=nlpot(Ia)%rRvec_iso(4,Ip)
               pos=0
               DO WHILE ( psp(Ity)%r_real(pos+1).lt.rRnorm .and. (pos.le.psp(Ity)%numps) )
                  pos=pos+1
               ENDDO
               !##CALCULATE POTENTIAL IN "RR"
               betarR=polynom(0,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,rRnorm)
            ENDIF
            !store data
            beta_init(Ip,Ipj)=betarR
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_beta_interp_r
   !---------------------apply Ylm to beta--------------------------
   SUBROUTINE nlp_beta_ylm_r(Ity,Ia,beta_Ylm)
      USE struct_module , ONLY : struct,volume
      USE pspot_module , ONLY : psp
      USE parameters ,   ONLY : Lpbc
      USE array_io
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ia,Ity
      REAL(DP),INTENT(INOUT) :: beta_Ylm(:,:)
      !LOCAL
      INTEGER(I4B) :: Ip,Ipj &
                   &, l,m
      REAL(DP) :: rvt(4),x,y,z,fac
      !> debug
      CHARACTER(len=1) :: str1,str2
      INTEGER(I4B)     :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      fac=1.d0
      DO Ip=1,nlpot(Ia)%npts
         IF(Lpbc)THEN
            rvt(:)=nlpot(Ia)%rRvec(:,Ip)
         ELSE
            rvt(:)=nlpot(Ia)%rRvec_iso(:,Ip)
         ENDIF
         IF(rvt(4)>0.d0)THEN
            x=rvt(1)/rvt(4)
            y=rvt(2)/rvt(4)
            z=rvt(3)/rvt(4)
         ELSE
            x=0.d0
            y=0.d0
            z=0.d0
         ENDIF
         !apply Ylm
         DO Ipj=1,psp(Ity)%nproj
            !(l,m)
            l=psp(Ity)%proj_l(Ipj)
            m=psp(Ity)%proj_m(Ipj)
            !apply ylm
            CALL apply_ylm(l,m,fac,x,y,z,beta_ylm(Ip,Ipj))
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_beta_ylm_r
   !--------------------apply phase to beta-------------------------
   SUBROUTINE nlp_beta_phase_r(Ity,Ia,beta_ylm,beta_phase)
      ! <r|e^-ikr|beta>
      USE grid_module ,ONLY : nk,KPT,grid
      USE struct_module , ONLY : struct,recip_lat
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia
      REAL(DP),INTENT(IN)     :: beta_ylm(:,:)
      COMPLEX(DP),INTENT(OUT)    :: beta_phase(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Ik,Ip !,Id
      REAL(DP) :: kvec(3),rvec(3)
      REAL(DP) :: kdotr
      COMPLEX(DP) :: e_ikr
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Ik=1,nk
         kvec(:)=KPT%vcar(:,Ik)
         !
         DO Ip=1,nlpot(Ia)%npts
            rvec(1:3)=nlpot(Ia)%rRvec(1:3,Ip)+struct%poscar(1:3,Ia)
            !> kdotr =  k \cdot r
            kdotr=DOT_PRODUCT(kvec,rvec)
            !> e^-ikr
            e_ikr=EXP(-1.d0*IMAG*kdotr)
            beta_phase(Ip,:,Ik)=e_ikr*beta_ylm(Ip,:)
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_beta_phase_r
   !-------------------apply Ylm on project-------------------------
   SUBROUTINE apply_ylm(l,m,fac,x,y,z,f)
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: l
      INTEGER(I4B),INTENT(IN) :: m
      REAL(DP),INTENT(IN)     :: fac
      REAL(DP),INTENT(IN) :: x,y,z
      REAL(DP),INTENT(INOUT) :: f
      !LOCAL
      REAL(DP) :: scal
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      scal= fac /SQRT(4.0_dp*pi) !*(-IMAG)**l
      !Apply Y_lm* (spherical harmonic)
      select case(l)
      case(0)
         select case(m)
         case(0)
            f = f*scal*(1)
         case default
            WRITE(6,*) 'Apply Ylm:abs(m)>l'
            STOP
         end select
      case(1)
         select case(m)
         case(-1)
            f = f*scal*(sqrt(3.0_dp)*x)
         case(0)
            f = f*scal*(sqrt(3.0_dp)*z)
         case(1)
            f = f*scal*(-(sqrt(3.0_dp)*y))
         case default
            WRITE(6,*) 'Apply Ylm:abs(m)>l'
            STOP
         end select
      case(2)
         select case(m)
         case(-2)
            f = f*scal*(-(sqrt(15.0_dp)*x*y))
         case(-1)
            f = f*scal*(sqrt(15.0_dp)*x*z)
         case(0)
            f = f*scal*((sqrt(5.0_dp)*(-1 + 3*z**2))/2)
         case(1)
            f = f*scal*(-(sqrt(15.0_dp)*y*z))
         case(2)
            f = f*scal*((sqrt(15.0_dp)*(-x**2 + y**2))/2)  !!changed following PARSEC by xlt in 2018.7.11
         case default
            WRITE(6,*) 'Apply Ylm:abs(m)>l'
            STOP
         end select
      case(3)
         select case(m)
         case(-3)
            f = f*scal*((sqrt(17.5_dp)*x*(-1 + 4*y**2 + z**2))/2)
         case(-2)
            f = f*scal*(-(sqrt(105.0_dp)*x*y*z))
         case(-1)
            f = f*scal*((sqrt(10.5_dp)*x*(-1 + 5*z**2))/2)
         case(0)
            f = f*scal*((sqrt(7.0_dp)*z*(-3 + 5*z**2))/2)
         case(1)
            f = f*scal*((sqrt(10.5_dp)*y*(1 - 5*z**2))/2)
         case(2)
            f = f*scal*((sqrt(105.0_dp)*(-x**2 + y**2)*z)/2)
         case(3)
            f = f*scal*(-(sqrt(17.5_dp)*y*(-3*x**2 + y**2))/2)
         case default
            WRITE(6,*) 'Apply Ylm:abs(m)>l'
            STOP
         end select
      case default
         WRITE(6,*) 'Apply Ylm:l>3 not programmed'
         STOP
      end select

      !if(full_trace) call trace_exit('ion_apply_ylm',status)
      return
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE apply_ylm
   !------------------------double grid------------------------------
   SUBROUTINE nlp_wijk_dg(ndg,inpol,numdg,lft,rit,wijk,drijk)
      USE math , ONLY: lagrange_interpolation_coe,&
           lagrange_interpolation_x
      USE grid_module , ONLY : n1,n2,n3
      USE struct_module , ONLY : lat_mat
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: ndg   & !dense-grid order
                              ! &,n_near & !dense-grid cal one side near a corse grid
                              &, inpol & !the order for interpolate of dense-grid
                              &, numdg & !numdg=2ndg-1
                              &, lft , rit  !box
      !the weight of this point
      REAL(DP),INTENT(OUT) :: wijk(numdg),drijk(3,numdg)
      !LOCAL
      !the small
      REAL(DP) :: smlat(3,3),pt(3),invf,scal  &
            &  ,pijk(numdg)
      INTEGER(I4B) :: ix,iy,iz,I,izi,iyi,ixi
      REAL(DP)     :: x_set(inpol+1)
      REAL(DP)     :: coe(inpol+1)
      INTEGER(I4B) :: map(lft:rit),j
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      smlat(:,1) = lat_mat(:,1)/n1/(ndg+1)
      smlat(:,2) = lat_mat(:,2)/n2/(ndg+1)
      smlat(:,3) = lat_mat(:,3)/n3/(ndg+1)
      !
      invf=1.d0/(ndg+1)
      ! scal=1.d0/(n_near+1)**3
      ! scal=1.d0/(2*n_near+1)**3
      scal=1.d0/(ndg+1)**3
      !start
      ! print*,'lft',lft,'rit',rit
      do I=lft,rit
         if(I<0)then
            map(I)=(abs(I)-1)/(ndg+1)+1+rit/(ndg+1)
         elseif(I>0)then
            map(I)=-(I-1)/(ndg+1)+rit/(ndg+1)
         elseif(I==0)then
            map(I)=0+rit/(ndg+1)
         endif
      enddo
      ! print*,'maps',map
      !> high order interpolate
      x_set=(/(i,i=0,inpol)/)
      coe=0.d0
      CALL lagrange_interpolation_coe(inpol+1,x_set,coe)
      ! print*,'coe',coe
      I=0
      DO iz=lft,rit
      DO iy=lft,rit
      DO ix=lft,rit
         I=I+1
         !relactive coo of coarse grid
         drijk(:,I)= smlat(:,1)*ix + &
              & smlat(:,2)*iy + &
              & smlat(:,3)*iz
         !frac coo
         pt(1)=get_lag_2(inpol,lft,ndg,invf,coe,map(ix),ix)
         pt(2)=get_lag_2(inpol,lft,ndg,invf,coe,map(iy),iy)
         pt(3)=get_lag_2(inpol,lft,ndg,invf,coe,map(iz),iz)
         !weight
         pijk(I)=pt(1)*pt(2)*pt(3)
         ! print*,'pt1',pt(1)
         ! if(ix==iy.and.ix==iz.and.ix==0)then
         !    print *,'pt',pt
         !    print*,'map',map(0)
         !    print *,'coe',coe
         ! endif
      ENDDO
      ENDDO
      ENDDO
      wijk(:)=pijk(:)*scal
       ! print *,'sum_wijk',wijk
       ! stop
   CONTAINS
      FUNCTION get_lag_2(order,lftt,ndgt,hdens,coet,coeid,ig)
         USE math , ONLY: lagrange_interpolation_x
         IMPLICIT NONE
         !IN/OUT
         INTEGER(I4B),INTENT(IN) :: ig    & !id for dense-grid
                               &,   ndgt  &   ! the number of dense grid between corse grid
                               &,   order &   ! the power of interpolation funtion
                               &,   lftt  &   ! lft for shift the interpolate points
                               &,   coeid    ! id of coe
         REAL(DP),INTENT(IN) :: hdens     &  !1/ndg
                               &, coet(order+1)
         REAL(DP) :: get_lag_2
         !LOCAL
         INTEGER(I4B) :: aig,i
         REAL(DP)     :: x_set(order+1)
         REAL(DP)     :: coe_x(order+1)
         INTEGER(I4B) :: shift,ig_shift
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         x_set=(/(i,i=0,order)/)
         if(ig<0)then
            ! shift=(abs(ig)-1)/(ndg+1)+1
            ! print *,'shift_1',shift
            shift=order+lftt/(ndg+1)+(abs(ig)-1)/(ndg+1)+1
            ! print *,'shift_2',shift
            ig_shift=ig+shift*(ndg+1)
         elseif(ig>0)then
            shift=lftt/(ndg+1)+order-(ig-1)/(ndg+1)
            ig_shift=ig+shift*(ndg+1)
         else
            shift=lftt/(ndg+1)+order-ig/(ndg+1)
            ig_shift=ig+shift*(ndg+1)
         endif
         coe_x=lagrange_interpolation_x(order+1,x_set,ig_shift*hdens)
         ! if(ig==0)print*,'coe_x',coe_x
         get_lag_2=coe_x(coeid)*coet(coeid)
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       ENDFUNCTION get_lag_2
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_wijk_dg
   SUBROUTINE nlp_wijk_dg_old(ndg,n_near,numdg,lft,rit,wijk,drijk)
      USE grid_module , ONLY : n1,n2,n3
      USE struct_module , ONLY : lat_mat
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: ndg   & !dense-grid order
                              &,n_near & !dense-grid cal one side near a corse grid
                              &, numdg & !numdg=2ndg-1
                              &, lft , rit  !box
      !the weight of this point
      REAL(DP),INTENT(OUT) :: wijk(numdg),drijk(3,numdg)
      !LOCAL
      !the small
      REAL(DP) :: smlat(3,3),pt(3),invf,scal  &
            &  ,pijk(numdg)
      INTEGER(I4B) :: ix,iy,iz,I,izi,iyi,ixi
      INTEGER(I4B) :: map(4*n_near+3),len_map
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      smlat(:,1) = lat_mat(:,1)/n1/(ndg+1)
      smlat(:,2) = lat_mat(:,2)/n2/(ndg+1)
      smlat(:,3) = lat_mat(:,3)/n3/(ndg+1)
      !
      invf=1.d0/(ndg+1)
      ! scal=1.d0/(n_near+1)**3
      scal=1.d0/(2*n_near+1)**3
      !start
      do I=1,n_near+1
         map(I)=lft+I-1
      enddo
      do I=n_near+2,3*n_near+2
         map(I)=I-2*n_near-2
      enddo
      do I=3*n_near+3,4*n_near+3
         map(I)=rit-4*n_near-3+I
      enddo
      len_map=4*n_near+3
      !> just a average
      ! I=0
      ! DO iz=-n_near,n_near
      ! DO iy=-n_near,n_near
      ! DO ix=-n_near,n_near
      !    I=I+1
      !    !relactive coo of coarse grid
      !    drijk(:,I)= smlat(:,1)*ix + &
      !         & smlat(:,2)*iy + &
      !         & smlat(:,3)*iz
      !    !frac coo
      !    pt(1)=get_lag(n_near+1,invf,ix)
      !    pt(2)=get_lag(n_near+1,invf,iy)
      !    pt(3)=get_lag(n_near+1,invf,iz)
      !    !weight
      !    pijk(I)=pt(1)*pt(2)*pt(3)
      ! ENDDO
      ! ENDDO
      ! ENDDO
      I=0
      DO izi=2,len_map-1
      DO iyi=2,len_map-1
      DO ixi=2,len_map-1
         iz=map(izi)
         iy=map(iyi)
         ix=map(ixi)
         I=I+1
         !relactive coo of coarse grid
         drijk(:,I)= smlat(:,1)*ix + &
              & smlat(:,2)*iy + &
              & smlat(:,3)*iz
         !frac coo
         pt(1)=get_lag(ndg,invf,ix)
         pt(2)=get_lag(ndg,invf,iy)
         pt(3)=get_lag(ndg,invf,iz)
         !weight
         pijk(I)=pt(1)*pt(2)*pt(3)
      ENDDO
      ENDDO
      ENDDO
      !
      wijk(:)=pijk(:)*scal
   CONTAINS
      FUNCTION get_lag(ndgt,hdens,ig)
         IMPLICIT NONE
         !IN/OUT
         INTEGER(I4B),INTENT(IN) :: ig    & !id for dense-grid
                               &,   ndgt    !ndg
         REAL(DP),INTENT(IN) :: hdens       !1/ndg
         REAL(DP) :: get_lag
         !LOCAL
         INTEGER(I4B) :: aig
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         aig=ABS(ig)
         IF (aig <= n_near.or. aig>=(ndgt-n_near)) then
            get_lag = (1.d0-aig*hdens) !&
         ! get_lag=(1.d0-aig/real(ndgt,DP))
         !&   ( 1.d0 + aig*hdens)*(1.d0 - aig*hdens)*(2.d0 - aig*hdens) / 2.d0
            !get_lag=1.d0-aig*hdens
         ELSE
            print *,'error : aig > ndg+1'
            !get_lag = &
      !   &   ( 1.d0 - aig*hdens)*(2.d0 - aig*hdens)*(3.d0 - aig*hdens) / 6.d0
            !get_lag=0.d0
         ENDIF
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ENDFUNCTION get_lag
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE nlp_wijk_dg_old
   !---------------------Nlp_Interp_pot------------------------------
   SUBROUTINE nlp_Interp_betalm_dg(Ity,Ia,Ip,numdg,lft,rit,wijk,drijk,beta0_IJK,betalm_IJK)
      USE math , ONLY : CubicHermiteInterp
      USE pspot_module , ONLY : psp,max_rcut
      USE struct_module , ONLY : struct
      USE math , ONLY : Norm
      !=============================
      !##FOR ISO
      USE parameters   , ONLY : Lpbc
      USE MathSplines
      !=============================
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia    &  !atom
                           & ,   Ip        &  !the point
                           & ,   numdg     &  !total number of dens-grid
                           & ,   lft,rit      ! dens-box
      REAL(DP),INTENT(IN) :: wijk(numdg),drijk(3,numdg)
      REAL(DP),INTENT(OUT) :: beta0_IJK(:),betalm_IJK(:)
      !LOCAL
      INTEGER(I4B) :: ix,iy,iz,I  &
                  & ,Ipj,l,m
      REAL(DP) :: rt(4,numdg) !,r(3,numdg)
      !|betaYlm> only
      REAL(DP) :: x,y,z,rmod
      REAL(DP),DIMENSION(numdg,psp(Ity)%nproj) ::  &
                &   beta0_dg   & !|beta>
                & , betalm_dg    !|beta_lm>
      !================================
      !##FOR ISO
      INTEGER(I4B) :: pos
      REAL(DP)     :: c(3)
      LOGICAL      :: L_isBoundary
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      L_isBoundary=.false.
      c=0.d0
      !stort all dens-grid
      I=0
      DO iz=lft,rit
      DO iy=lft,rit
      DO ix=lft,rit
         I=I+1
         !coo center of the atom
         rt(1:3,I)=drijk(1:3,I)+nlpot(Ia)%rRvec(1:3,Ip)
         rt(4,I)=Norm(rt(1:3,I))
         if(rt(4,I)>max_rcut)L_isBoundary=.true.
      ENDDO
      ENDDO
      ENDDO
      !calculate the dens-grid potential
      beta0_dg=0.d0
      betalm_dg=0.d0
      I=0
      DO iz=lft,rit
      DO iy=lft,rit
      DO ix=lft,rit
         I=I+1
         !cosin
         IF(rt(4,I)>0.d0)THEN
            x=rt(1,I)/rt(4,I)
            y=rt(2,I)/rt(4,I)
            z=rt(3,I)/rt(4,I)
            rmod=rt(4,I)
         ELSE
            x=0.d0
            y=0.d0
            z=0.d0
            rmod=0.d0
         ENDIF
         !Interploate the value of desen-grid
         DO Ipj=1,psp(Ity)%nproj
            !cal beta only
            l=psp(Ity)%proj_l(Ipj)
            m=psp(Ity)%proj_m(Ipj)
            !============================================================================================
            if(L_isBoundary.and.(rmod>nlpot(Ia)%rRvec(4,Ip)))cycle
            IF(Lpbc)THEN
              beta0_dg(I,Ipj)=CubicHermiteInterp(psp(Ity)%beta_r(:,Ipj),psp(Ity)%dbeta_dr(:,Ipj) &
                &    ,psp(Ity)%rmax, psp(Ity)%rspacing , rmod )
            ELSE
              !##INTERPOLATING BY THREE POINT
              !##FIND THE THREE POINT
              pos=0
              DO WHILE ( psp(Ity)%r_real(pos+1).lt.rmod .and. (pos.le.psp(Ity)%numps))
                pos=pos+1
              ENDDO
              !##CALCULATE POTENTIAL IN "RR"
              beta0_dg(I,Ipj)=polynom(0,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,rmod)
            ENDIF
            !============================================================================================
            !apply Ylm
            betalm_dg(I,Ipj)=beta0_dg(I,Ipj)
            CALL apply_ylm(l,m,1.d0,x,y,z,betalm_dg(I,Ipj))
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !
      DO Ipj=1,psp(Ity)%nproj
         beta0_IJK(Ipj)=SUM(beta0_dg(:,Ipj)*wijk(:))
         betalm_IJK(Ipj)=SUM(betalm_dg(:,Ipj)*wijk(:))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_Interp_betalm_dg
   !---------------------nlp_Interp_betalm_dg_iso---------------------------
   SUBROUTINE nlp_Interp_betalm_dg_iso(Ity,Ia,Ip,numdg,lft,rit,wijk,drijk,beta0_IJK,betalm_IJK)
      USE math , ONLY : CubicHermiteInterp
      USE pspot_module , ONLY : psp,max_rcut
      USE struct_module , ONLY : struct
      USE math , ONLY : Norm
      !=============================
      !##FOR ISO
      USE parameters   , ONLY : Lpbc
      USE MathSplines
      !=============================
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia    &  !atom
                           & ,   Ip        &  !the point
                           & ,   numdg     &  !total number of dens-grid
                           & ,   lft,rit      ! dens-box
      REAL(DP),INTENT(IN) :: wijk(numdg),drijk(3,numdg)
      REAL(DP),INTENT(OUT) :: beta0_IJK(:),betalm_IJK(:)
      !LOCAL
      INTEGER(I4B) :: ix,iy,iz,I  &
                  & ,Ipj,l,m
      REAL(DP) :: rt(4,numdg) !,r(3,numdg)
      !|betaYlm> only
      REAL(DP) :: x,y,z,rmod
      REAL(DP),DIMENSION(numdg,psp(Ity)%nproj) ::  &
                &   beta0_dg   & !|beta>
                & , betalm_dg    !|beta_lm>
      !================================
      !##FOR ISO
      INTEGER(I4B) :: pos
      REAL(DP)     :: c(3)
      ! LOGICAL      :: L_isBoundary
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! L_isBoundary=.false.
      c=0.d0
      !stort all dens-grid
      I=0
      !DO iz=lft,rit
      !DO iy=lft,rit
      !DO ix=lft,rit
      DO I=1,numdg
         !I=I+1
         !coo center of the atom
         rt(1:3,I)=drijk(1:3,I)+nlpot(Ia)%rRvec_iso(1:3,Ip)
         ! rt(1:3,I)=nlpot(Ia)%rRvec_iso(1:3,Ip)
         rt(4,I)=Norm(rt(1:3,I))
         ! if(rt(4,I)>max_rcut)L_isBoundary=.true.
      !ENDDO
      !ENDDO
      !ENDDO
      ENDDO
      !calculate the dens-grid potential
      beta0_dg=0.d0
      betalm_dg=0.d0
      I=0
      !DO iz=lft,rit
      !DO iy=lft,rit
      !DO ix=lft,rit
      DO I=1,numdg
         !I=I+1
         !>cosin
         IF(rt(4,I)>0.d0)THEN
            x=rt(1,I)/rt(4,I)
            y=rt(2,I)/rt(4,I)
            z=rt(3,I)/rt(4,I)
            rmod=rt(4,I)
         ELSE
            x=0.d0
            y=0.d0
            z=0.d0
            rmod=0.d0
         ENDIF
         !Interploate the value of desen-grid
         DO Ipj=1,psp(Ity)%nproj
            !cal beta only
            l=psp(Ity)%proj_l(Ipj)
            m=psp(Ity)%proj_m(Ipj)
            !============================================================================================
            ! if(L_isBoundary.and.(rmod>nlpot(Ia)%rRvec_iso(4,Ip)))cycle
            IF(Lpbc)THEN
              beta0_dg(I,Ipj)=CubicHermiteInterp(psp(Ity)%beta_r(:,Ipj),psp(Ity)%dbeta_dr(:,Ipj) &
                &    ,psp(Ity)%rmax, psp(Ity)%rspacing , rmod )
            ELSE
              !##INTERPOLATING BY THREE POINT
              !##FIND THE THREE POINT
              pos=0
              DO WHILE ( psp(Ity)%r_real(pos+1).lt.rmod .and. (pos.le.psp(Ity)%numps))
                pos=pos+1
              ENDDO
              !##CALCULATE POTENTIAL IN "RR"
              beta0_dg(I,Ipj)=polynom(0,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,rmod)
            ENDIF
            !============================================================================================
            !apply Ylm
            betalm_dg(I,Ipj)=beta0_dg(I,Ipj)
            CALL apply_ylm(l,m,1.d0,x,y,z,betalm_dg(I,Ipj))
         ENDDO
      !ENDDO
      !ENDDO
      !ENDDO
      ENDDO
      !
      DO Ipj=1,psp(Ity)%nproj
         beta0_IJK(Ipj)=SUM(beta0_dg(:,Ipj)*wijk(:))
         betalm_IJK(Ipj)=SUM(betalm_dg(:,Ipj)*wijk(:))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE nlp_Interp_betalm_dg_iso
   !---------------------nlp_beta_ylm_r_dg---------------------------
   SUBROUTINE nlp_beta_ylm_r_dg(Ity,Ia,numdg,lft,rit,wijk,drijk,beta0,beta_Ylm)
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia,numdg &
                               &,lft,rit
      REAL(DP),INTENT(IN)    :: wijk(:),drijk(:,:)
      REAL(DP),INTENT(OUT) :: beta0(:,:) &   !|beta>
                         &,   beta_Ylm(:,:) !|beta Ylm>
      !LOCAL
      INTEGER(I4B) :: Ip
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Ip=1,nlpot(Ia)%npts
         CALL nlp_Interp_betalm_dg(Ity,Ia,Ip,numdg,lft,rit,wijk,drijk &
                 &     , beta0(Ip,:), beta_Ylm(Ip,:))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_beta_ylm_r_dg
   !---------------------nlp_beta_ylm_r_dg_iso---------------------------
   SUBROUTINE nlp_beta_ylm_r_dg_iso(Ity,Ia,numdg,lft,rit,wijk,drijk,beta0,beta_Ylm)

     USE smpi_math_module

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia,numdg &
                               &,lft,rit
      REAL(DP),INTENT(IN)    :: wijk(:),drijk(:,:)
      REAL(DP),INTENT(OUT) :: beta0(:,:) &   !|beta>
                         &,   beta_Ylm(:,:) !|beta Ylm>
      !LOCAL
      INTEGER(I4B) :: Ip
      !##PLOT
      INTEGER(I4B) :: Ipj,l,m

      !> plot
      ! if(parallel%isroot)then
      !    print*,'ipj'
      !    read(*,*),Ipj
      !    open(1173,file='f.dat')
      !    open(1174,file='r.dat')
      ! endif
      !>
      DO Ip=1,nlpot(Ia)%npts
         CALL nlp_Interp_betalm_dg_iso(Ity,Ia,Ip,numdg,lft,rit,wijk,drijk &
              &     , beta0(Ip,:), beta_Ylm(Ip,:))
         !> plot
         ! if(parallel%isroot)then
         !    write(1173,*)beta0(Ip,Ipj)
         !    write(1174,*)nlpot(Ia)%rRvec_iso(4,Ip)
         ! endif
      ENDDO

      !> plot
      ! if(parallel%isroot)then
      ! close(1173)
      ! close(1174)
      ! endif
      ! CALL MPI_BARRIER(parallel%comm,mpinfo)
      ! stop
    ENDSUBROUTINE nlp_beta_ylm_r_dg_iso
    !-------------------build Project to real space------------------
    SUBROUTINE set_beta_real_dg()
      !for parpare in real space double grid  |e^(-ikr)\beta_lm>
      USE parameters , ONLY : NDG,LBvK,Lpbc,n_near,inpol
      USE grid_module , ONLY : dvol   ,n1,n2,n3
      USE struct_module , ONLY : struct,naty
      USE pspot_module , ONLY : psp,max_nproj
      USE math , ONLY : mat2thr
      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !LOCAL
      INTEGER(I4B) :: Ity,Ia,Ilft,Irit
      INTEGER(I4B) :: numdg  !for double grid
      REAL(DP),ALLOCATABLE :: wijk(:)  & !weight
                          &, drijk(:,:)      !d vector
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(max_nproj==0) RETURN
      !get the box for dense-grid
      ! print*,'inpol',inpol
      ! print*,'(inpol+2)/2',(inpol+2)/2
      ! print*,'(inpol+1)/2',(inpol+1)/2
      Ilft=-(ndg+1)*((inpol+2)/2) !1-2*(ndg+1)
      Irit= (ndg+1)*((inpol+1)/2) !2*(ndg+1)-1
      !total points
      numdg=((inpol+1)*(ndg+1)+1)**3
      !numdg=(n_near*2+1+2*n_near)**3  !> ndg and n_near,
      ! numdg=(2*n_near+1)**3          !> just an average
      ALLOCATE(wijk(numdg),drijk(3,numdg))
      call memory_sum('dg_wijk,drijk',real(size(wijk),DP)*DP&
           &+size(drijk)*DP)
      !calculate the weight
      CALL nlp_wijk_dg(NDG,inpol,numdg,Ilft,Irit,wijk,drijk)
      ! print *,'sum(wijk)',sum(wijk)
      !cycle the all species
      DO Ity=1,naty
         !
         IF(psp(Ity)%nproj==0) CYCLE
         !cycle all atom in this species
         DO Ia = struct%eleid(Ity),struct%eleid(Ity+1)-1
            !Interpolate projector and apply_ylm
            if(lpbc)then
               CALL nlp_beta_ylm_r_dg(Ity,Ia,numdg,Ilft,Irit,wijk,drijk &
                    &  ,nlpot(Ia)%proj0_dg,nlpot(Ia)%proj_dg)
            else
               CALL nlp_beta_ylm_r_dg_iso(Ity,Ia,numdg,Ilft,Irit,wijk,drijk &
                    &  ,nlpot(Ia)%proj0_dg,nlpot(Ia)%proj_dg)
            endif
            IF((.NOT.LBvK).AND.Lpbc)THEN
               !apply phase
               CALL nlp_beta_phase_r(Ity,Ia,nlpot(Ia)%proj_dg,nlpot(Ia)%proj_phs_dg)
            ENDIF
         ENDDO
      ENDDO
      !
      call memory_free('dg_wijk,drijk',real(size(wijk),DP)*DP&
           &+size(drijk)*DP)
      DEALLOCATE(wijk,drijk)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE set_beta_real_dg

  SUBROUTINE initialize_nlpot_band()
     USE parameters , ONLY : LDG,LBvK
     USE pspot_module , ONLY : max_nproj,max_rcut,psp
     USE struct_module , ONLY : struct,natom,lat_mat,naty

     USE grid_module ,ONLY : gap,n=>global_n,nk,grid , n1,n2,n3
     USE math , ONLY : mat2thr
      USE m_time_evaluate, ONLY: memory_sum,memory_free
      USE smpi_math_module, ONLY:parallel,MPI_REAL8,mpinfo
      IMPLICIT NONE
      !LOCAL
      REAL(DP) :: rcut2
      REAL(DP) :: ra(3),dR(3),dist2
      !!
      INTEGER(I4B) :: nrepl
      INTEGER(I4B) ::  nlnpts(natom)  !num of points
      !!index
      INTEGER(I4B) :: Ia,m,Ity  &
                &  ,  ip        &
                &  ,icx,icy,icz,ix,iy,iz
      REAL(DP)     :: rat(3)
      REAL(DP)     :: rVec(4,n)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call MPI_ALLGATHERV(grid%rvec,parallel%mygrid_range(3)*4,&
          & MPI_REAL8,rVec,4*parallel%recvcounts,4*parallel%displs&
          &, MPI_REAL8, parallel%commx,mpinfo)
      IF(max_nproj==0)THEN
         WRITE(6,*) 'Local pseudopotential has been used'
         RETURN
      ENDIF

      if(parallel%isroot)WRITE(6,*) 'Non-Local pseudopotential has been used'
      !
      nrepl=1
      rcut2=max_rcut**2
      !inqure the number of points we need
      !tes=0.d0
      nlnpts=0
      DO Ia=1,natom
         ra(:)=struct%poscar(:,Ia)
         m=0
         DO ip=1,n
            !
            DO icz=-nrepl,nrepl
               DO icy=-nrepl,nrepl
                  DO icx=-nrepl,nrepl
                     !ret=rcell + rvec - ra
                     rat(:)= icx*lat_mat(:,1)   &
                          & +  icy*lat_mat(:,2)   &
                          & +  icz*lat_mat(:,3)   &
                          & +   rVec(1:3,ip) &
                          & -  ra(:)
                     dist2=rat(1)**2+rat(2)**2+rat(3)**2
                     !
                     IF(dist2 <= rcut2) m=m+1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         nlnpts(Ia)=m
      ENDDO
      !===================store the data================
      max_nlnpts=MAXVAL(nlnpts(:))

      !
      CALL destroy_nlpot()
      !for projectors
      ALLOCATE(nlpot(natom))
      DO Ia=1,natom
         !>>>
         ALLOCATE(nlpot(Ia)%Id(nlnpts(Ia)))
         ALLOCATE(nlpot(Ia)%rRvec(4,nlnpts(Ia)))
         call memory_sum('nlpot_Id,rRvec',real(size(nlpot(Ia)%Id),DP)*I4B&
              &+size(nlpot(Ia)%rRvec)*DP)
         !store
         nlpot(Ia)%npts=nlnpts(Ia)
      ENDDO
      !projector
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            IF(.not.LDG)THEN
               ALLOCATE(nlpot(Ia)%proj0(nlnpts(Ia),psp(Ity)%nproj))
               ALLOCATE(nlpot(Ia)%proj(nlnpts(Ia),psp(Ity)%nproj))
               call memory_sum('nlpot_proj0,proj',&
                    & real(size(nlpot(Ia)%proj0),DP)*DP&
                    &+size(nlpot(Ia)%proj)*DP)
               IF(.NOT.LBvK)THEN
                  ALLOCATE(nlpot(Ia)%proj_phs(nlnpts(Ia),psp(Ity)%nproj,nk))
                  call memory_sum('nlpot_proj_phs',&
                       & real(size(nlpot(Ia)%proj_phs),DP)*DP*2)
               ENDIF
               !SF
               !ALLOCATE(nlpot(Ia)%Sfact(nk))
            ELSE
               !DG
               ALLOCATE(nlpot(Ia)%proj0_dg(nlnpts(Ia),psp(Ity)%nproj))
               ALLOCATE(nlpot(Ia)%proj_dg(nlnpts(Ia),psp(Ity)%nproj))
               call memory_sum('nlpot_proj0_dg,proj_dg',&
                    & real(size(nlpot(Ia)%proj0_dg),DP)*DP&
                    &+size(nlpot(Ia)%proj_dg)*DP)
               IF(.NOT.LBvK)THEN
                  ALLOCATE(nlpot(Ia)%proj_phs_dg(nlnpts(Ia),psp(Ity)%nproj,nk))
                  call memory_sum('nlpot_proj_phs_dg',&
                       & real(size(nlpot(Ia)%proj_phs_dg),DP)*DP*2)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !===================store data====================
      DO Ia=1,natom
         ra(:)=struct%poscar(:,Ia)
         m=0
         DO ip=1,n
            !
            DO icz=-nrepl,nrepl
               DO icy=-nrepl,nrepl
                  DO icx=-nrepl,nrepl
                     !ret=rcell + rvec - ra
                     rat(:)= icx*lat_mat(:,1)   &
                          & +  icy*lat_mat(:,2)   &
                          & +  icz*lat_mat(:,3)   &
                          & +   rVec(1:3,ip) &
                          & -  ra(:)
                     dist2=rat(1)**2+rat(2)**2+rat(3)**2
                     !
                     IF(dist2 <= rcut2) THEN
                        m=m+1
                        nlpot(Ia)%Id(m)=ip
                        nlpot(Ia)%rRvec(1:3,m)=rat(1:3)
                        nlpot(Ia)%rRvec(4,m)=SQRT(dist2)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !check
         IF(m/=nlpot(Ia)%npts)THEN
            WRITE(6,*) 'Error: Initialize nolpot',m,nlpot(Ia)%npts
            STOP
         ENDIF

      ENDDO
      !CALL nlp_struct_fact(nlpot_Sfact)
      !apply Ylm phase
      IF(LDG)THEN
         !Use time-saving method
         CALL set_beta_real_dg()
      ELSE
         !Use uniform grid
         CALL set_beta_real()
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE initialize_nlpot_band

   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE nlpot_module
