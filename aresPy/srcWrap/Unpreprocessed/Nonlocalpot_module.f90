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
        COMPLEX(DCP),ALLOCATABLE :: proj_phs(:,:,:)  !e^-ikr|dv\phi_lm>
        !Double-Grid
        REAL(DP),ALLOCATABLE :: proj0_dg(:,:)  !|dv\phi>
        REAL(DP),ALLOCATABLE :: proj_dg(:,:)   !|dv\phi_lm>
        COMPLEX(DCP),ALLOCATABLE :: proj_phs_dg(:,:,:)  !e^-ikr|dv\phi_lm>
   ENDTYPE nol_type
   !
   TYPE(nol_type),ALLOCATABLE :: nlpot(:)
   !
   INTEGER(I4B) :: max_nlnpts !max points
CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   SUBROUTINE initialize_nlpot()
      USE parameters , ONLY : lcore_val
      USE pspot_module , ONLY : max_nproj,max_rcut,psp
      USE struct_module , ONLY : struct,natom,lat_mat,naty
      USE grid_module ,ONLY : gap,n,nk,grid,n1,n2,n3,dvol
      USE math , ONLY : mat2thr
#ifdef MPI
      USE smpi_math_module
#endif
      IMPLICIT NONE
      !LOCAL
      REAL(DP) :: rcut2
      REAL(DP) :: ra(3),dR(3),dist2
      !!
      INTEGER(I4B) :: nrepl=1
      INTEGER(I4B) ::  nlnpts(natom)  !num of points
      !!index
      INTEGER(I4B) :: Ia,m,Ity  &
                &  ,  ip        &
                &  ,icx,icy,icz,ix,iy,iz
      REAL(DP)     :: rat(3)
#ifdef MPI
      REAL(DP) :: ncore_loc,ncore
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(max_nproj==0)THEN
#ifdef MPI
         IF(parallel%isroot)THEN
#endif
         WRITE(*,*) 'Local pseudopotential has been used'
#ifdef MPI
         ENDIF
#endif
         RETURN
      ENDIF
      
#ifdef MPI
      IF(parallel%isroot)THEN
#endif
      WRITE(*,*) 'Non-Local pseudopotential has been used'
#ifdef MPI
      ENDIF
#endif
      !
      IF(lcore_val)THEN
         CALL nlp_init_PartialCore(grid%rhoc)
#ifdef MPI
      ncore_loc=SUM(grid%rhoc)*dvol
      CALL MPI_ALLREDUCE(ncore_loc,ncore,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo) 

      IF(parallel%isroot)THEN
#else
         ncore=SUM(grid%rhoc)*dvol
#endif
         WRITE(*,*) 'Non-Linear core-valance correlation has been used'
         PRINT*,'Num. of core charges',ncore
#ifdef MPI
      ENDIF
#endif
      ENDIF
      !
      !rcut2=max_rcut**2
      !inqure the number of points we need
!tes=0.d0
      nlnpts(:)=0
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
         !the r cut of psp
         rcut2=psp(Ity)%rcut**2
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
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
      ENDDO
      !===================store the data================
      max_nlnpts=MAXVAL(nlnpts(:))
      !
      CALL destroy_nlpot()
      !For structure factor
      !ALLOCATE(nlpot_Sfact(nk))
      !for projectors
      ALLOCATE(nlpot(natom))
      DO Ia=1,natom
         !>>>
         ALLOCATE(nlpot(Ia)%Id(nlnpts(Ia)))
         ALLOCATE(nlpot(Ia)%rRvec(4,nlnpts(Ia)))
         !store
         nlpot(Ia)%npts=nlnpts(Ia)
      ENDDO
      !projector
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            ALLOCATE(nlpot(Ia)%proj0(nlnpts(Ia),psp(Ity)%nproj))
            ALLOCATE(nlpot(Ia)%proj(nlnpts(Ia),psp(Ity)%nproj))
            ALLOCATE(nlpot(Ia)%proj_phs(nlnpts(Ia),psp(Ity)%nproj,nk))
         ENDDO
      ENDDO
      !===================store data====================
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
         !rcut
         rcut2=psp(Ity)%rcut**2
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
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
                WRITE(*,*) 'Error: Initialize nolpot',m,nlpot(Ia)%npts
                STOP
            ENDIF

         ENDDO
      ENDDO
      !==================parpare Data==================
      !Use uniform grid
      CALL set_beta_real()
!OPEN(118,FILE='pj.dat')
!   WRITE(118,*) nlpot(1)%proj0
!CLOSE(118)
!OPEN(119,FILE='pjdg.dat')
!   WRITE(119,*) nlpot(1)%proj0_dg
!CLOSE(119)
!print*,'SIZE',SIZE(nlpot(1)%proj0),SIZE(nlpot(1)%proj0_dg)
!STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE initialize_nlpot
   !-----------------destroy array----------------------
   SUBROUTINE destroy_nlpot()
      USE struct_module , ONLY : natom
      IMPLICIT NONE
      INTEGER(I4B) :: Ia
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      IF(ALLOCATED(nlpot))THEN
         DO Ia=1,natom
               IF(ALLOCATED(nlpot(Ia)%Id))     DEALLOCATE(nlpot(Ia)%Id)
               IF(ALLOCATED(nlpot(Ia)%rRvec))  DEALLOCATE(nlpot(Ia)%rRvec)
               IF(ALLOCATED(nlpot(Ia)%proj0))    DEALLOCATE(nlpot(Ia)%proj0)
               IF(ALLOCATED(nlpot(Ia)%proj))     DEALLOCATE(nlpot(Ia)%proj)
               IF(ALLOCATED(nlpot(Ia)%proj_phs)) DEALLOCATE(nlpot(Ia)%proj_phs)
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
      !USE parameters , ONLY : LBvK
      USE grid_module , ONLY : dvol   ,n1,n2,n3
      USE struct_module , ONLY : struct,naty
      USE pspot_module , ONLY : psp,max_nproj
      USE math , ONLY : mat2thr
      IMPLICIT NONE
      !LOCAL
      INTEGER(I4B) :: Ity,Ia
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(max_nproj==0) RETURN
      !============apply Ylm in real space===============
      !cycle the all species
      DO Ity=1,naty
         !
         IF(psp(Ity)%nproj==0) CYCLE
         !cycle all atom in this species
         DO Ia = struct%eleid(Ity),struct%eleid(Ity+1)-1
            !interpolates data to mesh
            CALL nlp_beta_interp_r(Ity,Ia,nlpot(Ia)%proj0)
            !apply_ylm
            nlpot(Ia)%proj(:,:)=nlpot(Ia)%proj0(:,:)
            CALL nlp_beta_ylm_r(Ity,Ia,nlpot(Ia)%proj)
            !IF(.NOT.LBvK)THEN
            !apply phase
            CALL nlp_beta_phase_r(Ity,Ia,nlpot(Ia)%proj,nlpot(Ia)%proj_phs)
            !ENDIF
         ENDDO
      ENDDO
      !===================apply phase=====================
      
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE set_beta_real
   !-----------------------interpolate data----------------------------
   SUBROUTINE nlp_beta_interp_r(Ity,Ia,beta_init)
      !USE math , ONLY : CubicSplineInterp 
      USE math , ONLY : interp
      USE pspot_module , ONLY : psp
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ity,Ia
      REAL(DP),INTENT(OUT) :: beta_init(:,:)
      !LOCAL
      INTEGER(I4B) :: Ipj,Ip,Id  &
                  &,  l,m
      REAL(DP) :: rRnorm,betarR
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Ipj=1,psp(Ity)%nproj
         !(l,m)
         l=psp(Ity)%proj_l(Ipj)
         m=psp(Ity)%proj_m(Ipj)
         !cycle all point os pot
         DO Ip=1,nlpot(Ia)%npts
            !index
            Id=nlpot(Ia)%Id(Ip)
            !
            rRnorm=nlpot(Ia)%rRvec(4,Ip)
            !interpole the value
            !betarR=CubicSplineInterp(psp(Ity)%beta_r(:,Ipj),psp(Ity)%ddbeta_dr2(:,Ipj), &
            !          & psp(Ity)%rmax, psp(Ity)%rspacing , rRnorm)

            betarR=interp(psp(Ity)%numps,psp(Ity)%beta_r(:,Ipj),psp(Ity)%r,rRnorm)
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
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ia,Ity
      REAL(DP),INTENT(INOUT) :: beta_Ylm(:,:)
      !LOCAL
      INTEGER(I4B) :: Ip,Ipj &
                   &, l,m
      REAL(DP) :: rvt(4),x,y,z,fac
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      fac=1.d0
      DO Ip=1,nlpot(Ia)%npts
         rvt(:)=nlpot(Ia)%rRvec(:,Ip)
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
!            print*,'r',rvt(4)
!            print*,'vw',beta_Ylm(Ip,Ipj)*2
            CALL apply_ylm(l,m,fac,x,y,z,beta_ylm(Ip,Ipj))
!            print*,'betalm',beta_Ylm(Ip,Ipj)*2
!pause
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
      COMPLEX(DCP),INTENT(OUT)    :: beta_phase(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Ik,Ip !,Id
      REAL(DP) :: kvec(3),rvec(3)
      REAL(DP) :: kdotr
      COMPLEX(DCP) :: e_ikr
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Ik=1,nk
         kvec(:)=KPT%vcar(:,Ik)
         !
         DO Ip=1,nlpot(Ia)%npts
            !Id=nlpot(Ia)%Id(Ip)
            rvec(1:3)=nlpot(Ia)%rRvec(1:3,Ip)+struct%poscar(1:3,Ia)
            !rvec(1:3)=grid%rVec(1:3,Id)
            !rvec(1:3)=0.d0
            !kdotr =  k \cdot r
            kdotr=DOT_PRODUCT(kvec,rvec)
            !e^-ikr
            e_ikr=EXP(-1.d0*IMAG*kdotr)
            !test
            !e_ikr=1.d0
            !e^(-ikr)|beta>
            beta_phase(Ip,:,Ik)=e_ikr*beta_ylm(Ip,:)
!print*,'betaylm',beta_ylm(Ip,:)
!print*,'beta_phase',beta_phase(Ip,:,Ik)
!pause

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
            WRITE(*,*) 'Apply Ylm:abs(m)>l'
            STOP
         end select
      case(1)
         select case(m)
         case(-1)
            f = f*scal*(sqrt(3.0_dp)*x)
         case(0)
            f = f*scal*(sqrt(3.0_dp)*z)
         case(1)
            f = f*scal*((-sqrt(3.0_dp)*y))
         case default
            WRITE(*,*) 'Apply Ylm:abs(m)>l'
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
            f = f*scal*((sqrt(15.0_dp)*(-x**2 + y**2))/2)
         case default
            WRITE(*,*) 'Apply Ylm:abs(m)>l'
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
            WRITE(*,*) 'Apply Ylm:abs(m)>l'
            STOP
         end select
      case default
         WRITE(*,*) 'Apply Ylm:l>3 not programmed'
         STOP
      end select

      !if(full_trace) call trace_exit('ion_apply_ylm',status)
      return
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE apply_ylm
   !-----------------------nlp_init_PartialCore---------------------
   SUBROUTINE nlp_init_PartialCore(rhoc)
      USE math , ONLY : interp
      USE grid_module , ONLY : grid,n !,dvol
      USE struct_module , ONLY : naty,struct,lat_mat
      USE pspot_module , ONLY : psp
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(OUT) :: rhoc(n)
      !LOCAL
      INTEGER(I4B) :: nrepl=1
      INTEGER(I4B) :: Ity,Ia,ip,icx,icy,icz
      REAL(DP) :: ra(3),rat(3),rRnorm,corec
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      rhoc(:)=0._DP
      !start
      !cycle all atom type
      DO Ity=1,naty
         !cycle all atoms
         IF(.NOT.psp(Ity)%lcore)cycle
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            ra(:)=struct%poscar(:,Ia)
            !cycle all mesh
            DO ip=1,n
               DO icz=-nrepl,nrepl
               DO icy=-nrepl,nrepl
               DO icx=-nrepl,nrepl
                  !ret=rcell + rvec - ra
                  rat(:)= icx*lat_mat(:,1)   &
                      & +  icy*lat_mat(:,2)   &
                      & +  icz*lat_mat(:,3)   &
                      & +   grid%rVec(1:3,ip) &
                      & -  ra(:)
                  !distance from core
                  rRnorm=SQRT(DOT_PRODUCT(rat,rat))
                  IF(rRnorm>=psp(Ity)%rmax)THEN
                     corec=0._DP
                  ELSE
                     !core density
                     corec=interp(psp(Ity)%numps,psp(Ity)%denc,psp(Ity)%r,rRnorm)
                  ENDIF
                  rhoc(ip)=rhoc(ip)+corec
               ENDDO
               ENDDO
               ENDDO
            ENDDO !cycle all mesh
         ENDDO !cycle all atoms
      ENDDO !cycle all type
      !
      !PRINT*,'Num. of core charges',SUM(rhoc)*dvol
      !OPEN(181,FILE='core.dat')
      !   WRITE(181,*) rhoc
      !CLOSE(181)
      !STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nlp_init_PartialCore
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE nlpot_module
