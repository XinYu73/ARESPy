!#############################################################!
!---           Calculate the force and stress              ---!
!---Author:              Qiang Xu                          ---!
!---Date:               2018/01/10                         ---!
!#############################################################!
MODULE forcestress_module
   USE constants
   IMPLICIT NONE
   REAL(DP) :: cellpress
!INTERFACE cal_force
!  module procedure cal_force_c
!  module procedure cal_force_r
!END INTERFACE cal_force
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   SUBROUTINE cal_force_c(rhoS,Uik,force)
      USE struct_module, ONLY: natom,struct,lat_mat
      USE ewald , ONLY : ewald_forces,ISO_ewald_forces
      USE pspot_module ,ONLY: psp
      USE parameters   ,ONLY: LBvK

      USE smpi_math_module, ONLY: parallel

      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      COMPLEX(DCP),INTENT(IN) :: Uik(:,:,:,:)
      REAL(DP),INTENT(OUT)   :: force(:,:)
!LOCAL
      REAL(DP) :: tmpf(3,natom)
      REAL(DP) :: fat(3),fsum(3)
      INTEGER(I4B) :: Ia
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      force(:,:)=0.d0
!##local potential part {{
      IF(.NOT.LBvK)THEN
        CALL locforce(rhoS,tmpf)
      ELSE
        CALL locforce_r(rhoS,tmpf)
      ENDIF
!##local potential part }}
! print*,'loc',tmpf
      force(:,:)=force(:,:) + tmpf(:,:)
!##ion-ion interaction {{
      tmpf=ewald_forces(lat_mat,struct%pos,struct%eleid,psp(:)%Zion)
! print*,'ewald',tmpf
      force(:,:)=force(:,:) + tmpf(:,:)
!===================================
!nlnlocal part
      CALL nonlforce(Uik,tmpf)
      force(:,:)=force(:,:) + tmpf(:,:)
! print*,'nloc',tmpf
!===================================
! PRINT*,'total force is :'
! PRINT*,force(:,:)
!net force
      fsum(1)=SUM(force(1,:))
      fsum(2)=SUM(force(2,:))
      fsum(3)=SUM(force(3,:))
!
      fat(:)=fsum(:)/REAL(natom,DP)
!
      DO Ia=1,natom
         force(:,Ia)=force(:,Ia) - fat(:)
      ENDDO
! WRITE(6,*) '[Total net force] (hartree/bohr)'

      if(parallel%isroot)THEN

      WRITE(6,*) '[TOTAL-FORCE] (hartree/bohr)'
      WRITE(6,*) '-----------------------------------------------------------------------------'
      WRITE(6,11)
      DO Ia=1,natom
         WRITE(6,'(I4,1X,3(f11.6,1X),2X,3(f11.6,1X))') &
              &    Ia,struct%pos(:,Ia),force(:,Ia)
      ENDDO
      WRITE(6,*) '-----------------------------------------------------------------------------'
      WRITE(6,*) '[Total press] (Gpa)  0.000000'
      WRITE(6,*) 'in kB  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000'
      WRITE(6,*) 'external pressure = 0.0 kB  Pullay stress = 0.00 kB'
      WRITE(6,*) '-----------------------------------------------------------------------------'

      endif

      11    FORMAT(' atom  ----x----   ----y----   ----z----' &
        &       ,  '    ----Fx----  ----Fy----  ----Fz----')
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_force_c
!-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE cal_force_r(rhoS,Uik,force)
      USE struct_module, ONLY : natom,struct,lat_mat
      USE ewald , ONLY : ewald_forces,ISO_ewald_forces
      USE pspot_module ,ONLY : psp
!xlt for judge parameters
      USE parameters, ONLY : Lpbc, LBvK,LDG

      USE smpi_math_module , ONLY:parallel

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)

      REAL(DP),INTENT(IN) :: Uik(:,:,:)

      REAL(DP),INTENT(OUT)   :: force(:,:)
!LOCAL
      REAL(DP) :: tmpf(3,natom)
      REAL(DP) :: fat(3),fsum(3)
      INTEGER(I4B) :: Ia
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('cal_force_local',real(size(tmpf),DP)*DP+6*DP)
      force(:,:)=0.d0
!> local potential part
      IF(Lpbc.AND.(.NOT.LBvK))THEN
        CALL locforce(rhoS,tmpf)
      ELSEIF(.NOT.lpbc)THEN
        CALL locforce_r(rhoS,tmpf)
      ENDIF
!>====--------=====-=====-====
! IF(parallel%isroot)PRINT*,'local part force is :'
! IF(parallel%isroot)PRINT*,tmpf(:,:)
!>====--------=====-=====-====
      force(:,:)=force(:,:) + tmpf(:,:)

!> CORRECT
! CALL locforce_correct_consistant(tmpf)
! force(:,:)=force(:,:) + tmpf(:,:)

!> ion-ion interaction
      IF(Lpbc.AND.(.NOT.LBvK))THEN
         tmpf=ewald_forces(lat_mat,struct%pos,struct%eleid,psp(:)%Zion)
      ELSEIF(.NOT.lpbc)THEN
         tmpf=ISO_ewald_forces(lat_mat,struct%poscar,struct%eleid,psp(:)%Zion)*0.5d0
      ELSE
         print*,"CHECK THE INPUT.DAT, PLEASE"
         STOP
      ENDIF
!>====--------=====-=====-====
! IF(parallel%isroot) PRINT*,'ewald part force is :'
! IF(parallel%isroot) PRINT*,tmpf(:,:)
!>====--------=====-=====-====
      force(:,:)=force(:,:) + tmpf(:,:)
!===================================
!nlnlocal part
      if(LDG)then
         CALL nonlforce_r_dg(Uik,tmpf)
      else
         CALL nonlforce_r(Uik,tmpf)
      endif
!>====--------=====-=====-====
! IF(parallel%isroot) PRINT*,'nonlocal part force is :'
! IF(parallel%isroot) PRINT*,tmpf(:,:)
!>====--------=====-=====-====
      force(:,:)=force(:,:) + tmpf(:,:)
!===================================
!net force
      fsum(1)=SUM(force(1,:))
      fsum(2)=SUM(force(2,:))
      fsum(3)=SUM(force(3,:))
!
      fat(:)=fsum(:)/REAL(natom,DP)
!
      DO Ia=1,natom
         force(:,Ia)=force(:,Ia) - fat(:)
      ENDDO

      IF(parallel%isroot)THEN
! WRITE(6,*) '[Total net force] (hartree/bohr)'
! WRITE(6,*) '[TOTAL-FORCE] (hartree/bohr)'
         WRITE(6,*) '[TOTAL-FORCE] (eV/angs)'
! WRITE(6,*) '-----------------------------------------------------------------------------'
         WRITE(6,11)
         DO Ia=1,natom
            WRITE(6,'(3(f11.6,1X),2X,3(f11.6,1X))') &
                 &    struct%pos(:,Ia),force(:,Ia)*au2ev_force
         ENDDO
         WRITE(6,*) '-----------------------------------------------------------------------------'
         WRITE(6,*) '[Total press] (Gpa)  0.000000'
         WRITE(6,*) 'in kB  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000'
         WRITE(6,*) 'external pressure = 0.0 kB  Pullay stress = 0.00 kB'
         WRITE(6,*) '-----------------------------------------------------------------------------'
      ENDIF

      11    FORMAT(' ----x----   ----y----   ----z----' &
        &       ,  '    ----Fx----  ----Fy----  ----Fz----')
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_free('cal_force_local',real(size(tmpf),DP)*DP+6*DP)
   ENDSUBROUTINE cal_force_r
!-----------------------PARTING-LINE--------------------------
!----------------------------------------------------------
   SUBROUTINE locforce(rhoS,lforce)

      USE smpi_math_module,  ONLY: parallel,mpinfo,mpi_real8&
           & ,mpi_sum
      USE grid_module, ONLY: global_n1,global_n2,global_n3

      USE parameters , ONLY : Nspin,PP_identifer
      USE grid_module , ONLY : ng1,ng2,ng3,grid,n1,n2,n3
      USE struct_module , ONLY : recip_lat , struct , naty
      USE pspot_module , ONLY : psp
      USE math , ONLY : CubicSplineInterp
      USE mathsplines , ONLY : polynom
      USE fourier,   only : fft
      IMPLICIT NONE
!IN/OUT
      REAL(DP),INTENT(IN)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: lforce(:,:)
!LOCAL
      REAL(DP) :: rho(n1,n2,n3)

      COMPLEX(DCP) :: rhog(ng1,ng2,parallel%local_z) , sfactor

      INTEGER(I4B) :: I1,I2,I3,Is,Ity,Ia,Ig
      REAL(DP) :: ra(3),qvec(3) &
             &    ,qnorm,vg,&
             & threedot(3),c(3) !interplote
      INTEGER(I4B) :: m

! REAL(DP) :: rho_global(global_n1,global_n2,global_n3)
      REAL(DP) :: rho_FFT(global_n1,global_n2,parallel%local_z)
      INTEGER(I4B) :: ix,iy
      REAL(DP) :: lforce_global(size(lforce,1),size(lforce,2))

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      rho(:,:,:)=0.d0
      DO Is=1,Nspin
         rho(:,:,:)=rho(:,:,:) + rhoS(:,:,:,Is)
      ENDDO

         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         CALL MPI_ALLTOALLV(rho(ix,iy,1),parallel%fft_scount&
              & ,parallel%fft_sdispls,MPI_REAL8&
              & ,rho_FFT,parallel%fft_rcount&
              & ,parallel%fft_rdispls,MPI_REAL8&
              & ,parallel%commx,mpinfo)
! call MPI_ALLGATHERV(rho(ix,iy,1),parallel%mygrid_range(3)&
!      &,MPI_REAL8,rho_global,parallel%recvcounts&
!      &,parallel%displs,MPI_REAL8&
!      & ,parallel%commx,mpinfo)

!
!open(11111,file="VlocqS1")
!write(11111,*)psp(1)%VlocqS
!close(11111)

         rhog = FFT(rho_FFT)

!cal force
      lforce(:,:)=0.d0
!cal
!===========================
!#PLOT
!open(23301,file="frxyz_0.2")
!write(23301,*)psp(1)%VlocqS
!close(23301)
!open(23302,file="rxyz_0.2")
!print*,"qmax/qspacing",psp(1)%qmax/psp(1)%qspacing,nint(psp(1)%qmax/psp(1)%qspacing),int(psp(1)%qmax/psp(1)%qspacing)
!DO I1=0,nint(psp(1)%qmax/psp(1)%qspacing)
!    WRITE(23302,*)I1*psp(1)%qspacing
!ENDDO
!WRITE(23302,*)psp(1)%r_real
!close(23302)
!open(23303,file="rxyz_0.3")
!open(23304,file="frxyz_0.3")
!===========================
      DO Ity=1,naty
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            ra(:)=struct%poscar(:,Ia)
            Ig=0

            DO I3=1,parallel%local_z
               Ig=(parallel%local_z_start+I3-1)*ng2*ng1

            DO I2=1,ng2
            DO I1=1,ng1
               Ig=Ig+1
               IF(grid%gMask(Ig))THEN
                  qvec(:)=grid%gVec(1:3,Ig)
                  qnorm=grid%gVec(4,Ig)
                  if(PP_identifer==0)then
                  vg = CubicSplineInterp(psp(Ity)%VlocqS,psp(Ity)%ddVl_dq2 , &
                    & psp(Ity)%qmax,psp(Ity)%qspacing,qnorm,psp(Ity)%Zion)
                  else
                  m=0
                  DO WHILE (psp(1)%qspacing*(m+1).lt.qnorm &
                       &.and. (m.le.psp(Ity)%qnumps))
                    m=m+1
                  ENDDO
                  IF(qnorm.lt.psp(Ity)%qmax)THEN
!##CALCULATE POTENTIAL IN "dist"
!===================================
!#uniform
                    IF((m+1)*psp(Ity)%qspacing.lt.psp(Ity)%qmax)THEN
                       threedot=(/(m-1)*psp(Ity)%qspacing,m*psp(Ity)%qspacing,(m+1)*psp(Ity)%qspacing/)
!threedot=(/(m-1)*0.001d0,m*0.001d0,(m+1)*0.001d0/)
                       vg=polynom(0,3,threedot,psp(Ity)%Vlocq(m-1:m+1),c,qnorm)
                    ELSE
                       threedot=(/(m-2)*psp(Ity)%qspacing,(m-1)*psp(Ity)%qspacing,(m)*psp(Ity)%qspacing/)
                       vg=0 !polynom(0,3,threedot,psp(Ity)%VlocqS(m-2:m),c,qnorm)
                    ENDIF
                  ELSE
                    vg=0.d0
                  ENDIF
                  endif !> PP_identifer==0
!======================================================================================
                  sfactor=EXP(IMAG*DOT_PRODUCT(qvec,ra))
!force
                  lforce(1,Ia) = lforce(1,Ia) + qvec(1)*vg*aimag(rhog(I1,I2,I3)*sfactor)
                  lforce(2,Ia) = lforce(2,Ia) + qvec(2)*vg*aimag(rhog(I1,I2,I3)*sfactor)
                  lforce(3,Ia) = lforce(3,Ia) + qvec(3)*vg*aimag(rhog(I1,I2,I3)*sfactor)
!============
!##PLOT
!IF(Ia==1)THEN
!  write(23303,*)qnorm
!  write(23304,*)vg
!ENDIF
!============
               ENDIF
            ENDDO
            ENDDO
            ENDDO
         ENDDO
      ENDDO
!=========================
!#PLOT
!close(23303)
!WRITE(23304,*)rho_t
!close(23304)
!=========================
!

      CALL MPI_ALLREDUCE(lforce,lforce_global,size(lforce),&
           &mpi_REAL8,mpi_sum,parallel%commx,mpinfo)
      lforce(:,:) =2.d0*lforce_global(:,:)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE locforce
!-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE locforce_r(rhoS,lforce)
      USE parameters , ONLY : Nspin
      USE grid_module , ONLY : n1,n2,n3,rho_calc, gap, dvol
      USE struct_module , ONLY : struct , naty,natom
      USE pspot_module , ONLY : psp
      USE mathsplines , ONLY : polynom

      USE smpi_math_module , ONLY : atom_split,parallel,smpi_reduce_sum_real_2d

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
!> IN/OUT
      REAL(DP),INTENT(IN)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: lforce(:,:)
!> LOCAL
      REAL(DP) :: rhoS_t(n1,n2,n3,NSPIN)
      REAL(DP) :: rho(n1,n2,n3)
      REAL(DP) :: ra(3),qvec(3) &
             &    ,qnorm,vg
      INTEGER(I4B) :: Ity, Ia, Is, k, i
      REAL(DP) :: partial_V !##{\frac{\partial V}{\partial r}}
      REAL(DP) :: c(5),r(5),f(5)   !##needed by polynom
      INTEGER(I4B) :: m   !##needed by polynom
      REAL(DP) :: temp_force, rik, rxy &
                  & , sinphi, sintheta, cosphi, costheta

!> For parallel
      integer(I4B) :: mysize,atom_id=0,id_core
      INTEGER(I4B),allocatable :: atom_index(:)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('locforce_r_local',(real(size(rhoS_t),DP)+size(rho))*DP)

! !==================================================
! !> initial parallel config
! IF(.not.allocated(atom_index))ALLOCATE(atom_index(natom/parallel%numprocs+1))
! ! CALL start_time('init_density_0')
! CALL atom_split(mysize,atom_index)
! ! print *,'atom_index',atom_index,'id',parallel%myid
! ! CALL end_time('init_density_0')
! ! CALL write_time('init_density_0')
! id_core=1
! !>=================================================

      rho(:,:,:)=0.d0
      DO Is=1,Nspin
         rho(:,:,:)=rho(:,:,:) + rhoS(:,:,:,Is)
      ENDDO
!##{cal force}
      lforce=0.0d0
! i=0  !##{natoms}
      DO Ity=1,naty
      DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
! #ifdef 1
!          !>=================================
!          If(Ia==atom_index(id_core))THEN
!             if(id_core<mysize)id_core=id_core+1
!          !>=================================
! #endif
! i=i+1
         temp_force=0.0d0
! open(unit=1001,file="frxyz_0.3")
! open(unit=1002,file="rxyz_0.3")

         DO k=1,parallel%mygrid_range(3)

!IF(i/=k)THEN
            rik=sqrt((rho_calc%x(k)*gap(1)-struct%poscar(1,Ia))**2+&
                 & (rho_calc%y(k)*gap(2)-struct%poscar(2,Ia))**2+&
                 & (rho_calc%z(k)*gap(3)-struct%poscar(3,Ia))**2)
            rxy=sqrt((rho_calc%x(k)*gap(1)-struct%poscar(1,Ia))**2+&
                 &(rho_calc%y(k)*gap(2)-struct%poscar(2,Ia))**2)
            IF(rxy.lt.0.0d0)THEN
!print*,"rxy",rxy
               sinphi=0.d0
               cosphi=0.d0
               rxy=0.d0
            ELSE
               sinphi=(rho_calc%y(k)*gap(2)-struct%poscar(2,Ia))/rxy
               cosphi=(rho_calc%x(k)*gap(1)-struct%poscar(1,Ia))/rxy
            ENDIF
            IF(rik.lt.psp(Ity)%r_real(8))cycle
            sintheta=rxy/rik
            costheta=(rho_calc%z(k)*gap(3)-struct%poscar(3,Ia))/rik
!##Cal force =========================
            IF(rik.gt.psp(Ity)%r_real(size(psp(Ity)%r_real)-8))THEN
               partial_V=struct%Zion(Ity)/rik**2
            ELSE
               m=1
               DO WHILE ( psp(Ity)%r_real(m).lt.rik .and. (m.le.psp(Ity)%numps))
                  m=m+1
               ENDDO
!##CALCULATE POTENTIAL IN "RR"
               r=psp(Ity)%r_real(m-2:m+2)
               f=psp(Ity)%V_loc(m-2:m+2)
               partial_V=polynom(1,5,r,f,c,rik) !*2.d0
            ENDIF
! write(1001,*)partial_V!rho(rho_calc%x(k),rho_calc%y(k),rho_calc%z(k))
! write(1002,*)rik
!=====================================
            temp_force=partial_V*rho(rho_calc%x(k),rho_calc%y(k),rho_calc%z(k))*dvol
            lforce(1,Ia)=lforce(1,Ia)+temp_force*sintheta*cosphi
            lforce(2,Ia)=lforce(2,Ia)+temp_force*sintheta*sinphi
            lforce(3,Ia)=lforce(3,Ia)+temp_force*costheta
!==============================================
!ENDIF
         ENDDO !> k
! close(1001)
! close(1002)
! #ifdef 1
!          endif
! #endif
      ENDDO !> Ia
      ENDDO !> Ity

      CALL smpi_reduce_sum_real_2d(lforce)

      call memory_free('locforce_r_local',(real(size(rhoS_t),DP)+size(rho))*DP)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE locforce_r
!-----------------------PARTING-LINE--------------------------
!    SUBROUTINE locforce_correct_consistant(locf_term3)
!      USE parameters , ONLY : Nspin
!      USE grid_module , ONLY : n1,n2,n3,n, gap, dvol, grid
!      USE struct_module , ONLY : struct , naty,natom
!      USE pspot_module , ONLY : psp
!      USE mathsplines , ONLY : polynom
!      USE potential_module , ONLY : V_hxc_old,V_hxc_new
! #ifdef 1
!      USE smpi_math_module , ONLY : atom_split,parallel,smpi_reduce_sum_real_2d
! #endif
!      !>variables
!      REAL(DP),INTENT(OUT)  :: locf_term3(:,:)
!      INTEGER(I4B)          :: Ity, Ia, i, k
!      INTEGER(I4B)          :: a,b,c
!      REAL(DP)              :: partial_rho
!      REAL(DP)              :: deltaV(n)
!      REAL(DP)              :: temp_force,rik,rxy,sinphi,sintheta,cosphi,&
!           & costheta,m,c3(3)
! #ifdef 1
!      integer(I4B) :: mysize,atom_id=0,id_core
!      INTEGER(I4B),allocatable :: atom_index(:)
! #endif
!      !======================================
!      !>delta V
!      i=0
!      DO c=1,n3,1
!      DO b=1,n2,1
!      DO a=1,n1,1
!        i=i+1
!        deltaV(i)=V_hxc_new(a,b,c)-V_hxc_old(a,b,c)
!      ENDDO
!      ENDDO
!      ENDDO
!      !===============================================================
! #ifdef 1
!      !> initial parallel config
!      IF(.not.allocated(atom_index))ALLOCATE(atom_index(natom/parallel%numprocs+1))
!      ! CALL start_time('init_density_0')
!      CALL atom_split(mysize,atom_index)
!      ! print *,'atom_index',atom_index,'id',parallel%myid
!      ! CALL end_time('init_density_0')
!      ! CALL write_time('init_density_0')
!      id_core=1
! #endif
!      !>calculate locf_term3 = {\delta V_hxc}*{\frac{\partial Rho}{\partial r}}
!      !##{cal force}
!      locf_term3=0.0d0
!      ! i=0  !##{natoms}
!      DO Ity=1,naty
!      DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
! #ifdef 1
!         If(Ia==atom_index(id_core))THEN
!            if(id_core<mysize)id_core=id_core+1
! #endif
!         ! i=i+1
!         temp_force=0.0d0
!         !open(unit=1001,file="frxyz_0.3")
!         !open(unit=1002,file="rxyz_0.3")
!         DO k=1,n
!            !IF(i/=k)THEN
!            rik=sqrt((grid%rVec(1,k)-struct%poscar(1,Ia))**2+&
!                 & (grid%rVec(2,k)-struct%poscar(2,Ia))**2+&
!                 & (grid%rVec(3,k)-struct%poscar(3,Ia))**2)
!            rxy=sqrt((grid%rVec(1,k)-struct%poscar(1,Ia))**2+&
!                 &(grid%rVec(2,k)-struct%poscar(2,Ia))**2)
!            IF(rxy.lt.0.000001d0)THEN
!              !print*,"rxy",rxy
!              sinphi=0.d0
!              cosphi=0.d0
!              rxy=0.d0
!            ELSE
!              sinphi=(grid%rVec(2,k)-struct%poscar(2,Ia))/rxy
!              cosphi=(grid%rVec(1,k)-struct%poscar(1,Ia))/rxy
!            ENDIF
!            IF(rik.lt.3.d0**0.5*gap(1))cycle
!            sintheta=rxy/rik
!            costheta=(grid%rVec(3,k)-struct%poscar(3,Ia))/rik
!            !##Cal force =========================
!            IF(rik.gt.psp(Ity)%r_real(size(psp(Ity)%r_real)))THEN
!              partial_rho=0.d0
!            ELSE
!              m=0
!              DO WHILE ( psp(Ity)%r_real(m+1).lt.rik )
!                m=m+1
!              ENDDO
!              !##CALCULATE POTENTIAL IN "RR"
!              partial_rho=polynom(1,3,psp(Ity)%r_real(m-1:m+1),psp(Ity)%denr(m-1:m+1),c3,rik)*2.d0
!            ENDIF
!            !write(1001,*)partial_V
!            !write(1002,*)rik
!            !=====================================
!            temp_force=partial_rho*deltaV(k)*dvol
!            locf_term3(1,Ia)=locf_term3(1,Ia)+temp_force*sintheta*cosphi
!            locf_term3(2,Ia)=locf_term3(2,Ia)+temp_force*sintheta*sinphi
!            locf_term3(3,Ia)=locf_term3(3,Ia)+temp_force*costheta
!         ENDDO !> k
! #ifdef 1
!         endif
! #endif
!      ENDDO !> Ia
!      ENDDO !> Ity
! #ifdef 1
!      CALL smpi_reduce_sum_real_2d(locf_term3)
! #endif
!    END SUBROUTINE

!------------------------nonlocal force--------------------
   SUBROUTINE nonlforce(Uik,nlforce)
      USE math , ONLY : SimpleInterp,r_dYlm
      USE MathSplines, ONLY: polynom
      USE nlpot_module , ONLY : nlpot,max_nlnpts
      USE struct_module ,ONLY : naty,struct,recip_lat
      USE pspot_module , ONLY : psp,max_nproj
      USE parameters , ONLY : Nspin,Nstates,PP_identifer
      USE grid_module , ONLY : nk,KPT,grid,dvol
      USE smearing_module , ONLY : wke

      USE smpi_math_module , ONLY : parallel ,mpinfo,mpi_real8,mpi_complex16,mpi_sum
      USE struct_module ,ONLY : natom
      USE grid_module, ONLY: n1,n2

      IMPLICIT NONE
!IN/OUT
      COMPLEX(DCP),INTENT(IN) :: Uik(:,:,:,:) !wave_function
      REAL(DP),INTENT(OUT) :: nlforce(:,:)
!LOCAL
      REAL(DP) :: xyzr(4,max_nlnpts)
      REAL(DP) :: beta,dbeta
!dYlm/dR  and  dvYlm/dR
      REAL(DP) :: dYlm(0:3) &
           & ,    dbeta_lm(0:3,max_nlnpts,max_nproj)
      COMPLEX(DCP) :: dGnlm(3),Gnlm,tmp0(3)

      COMPLEX(DCP) :: tmp0_global(3)

!psi_ik=e^ikr U_ik
      COMPLEX(DCP)  :: psi_ik(max_nlnpts,Nstates,nk,NSPIN) &
                &,    eikr
      REAL(DP)     :: kdotr , kvec(3),rvec(3)
!
      INTEGER(I4B) :: Is,Ik,Ii
      INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id,l,m
      REAL(DP) :: tmpf(3)
!> inteplote
      REAL(DP) :: c(3)
      INTEGER(I4B) :: pos

      COMPLEX(DCP) :: Gnlm_global
      REAL(DP) :: nlforce_local(3,natom)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      nlforce(:,:)=0.d0
!Cycle all species
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
!cycle all atoms
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!1.Parpare some data
!(1) Find the nlpot points' distant and cosine
!(2) calculate the \psi=eikr*Uik(r)
!cycle all points around this atom
            DO Ip=1,nlpot(Ia)%npts
!1.store points' data
               xyzr(4,Ip)=nlpot(Ia)%rRvec(4,Ip)
               IF(xyzr(4,Ip)>0.d0)THEN
                  xyzr(1,Ip)=nlpot(Ia)%rRvec(1,Ip)/xyzr(4,Ip)
                  xyzr(2,Ip)=nlpot(Ia)%rRvec(2,Ip)/xyzr(4,Ip)
                  xyzr(3,Ip)=nlpot(Ia)%rRvec(3,Ip)/xyzr(4,Ip)
               ELSE
                  xyzr(1:4,Ip)=0.d0
               ENDIF
!calculate psi_ik at this point
               Id=nlpot(Ia)%Id(Ip)
               rvec(:)=nlpot(Ia)%rRvec(1:3,Ip)+struct%poscar(:,Ia)
               DO Ik=1,nk
!e^ikr
                  kvec(:)=KPT%vcar(:,Ik)
                  kdotr=DOT_PRODUCT(kvec(:),rvec(:))
                  eikr=EXP(IMAG*kdotr)
                  psi_ik(Ip,:,Ik,:)=eikr*Uik(Id,:,Ik,:)

               ENDDO

            ENDDO
!2.Build the dbeta_lm/dr for this atom
!cycle all projectors
            DO Ipj=1,psp(Ity)%nproj
!get (l,m)
               l=psp(Ity)%proj_l(Ipj)
               m=psp(Ity)%proj_m(Ipj)
!---------------------------------
               DO Ip=1,nlpot(Ia)%npts
!> beta
                  beta=nlpot(Ia)%proj0(Ip,Ipj)
!> dbeta
                  if(PP_identifer==0)then
                  dbeta=SimpleInterp( psp(Ity)%dbeta_dr(:,Ipj),psp(Ity)%rmax, &
                       & psp(Ity)%rspacing,xyzr(4,Ip) )
                  else
                  pos=0
                  DO WHILE ( psp(Ity)%r_real(pos+1).lt.xyzr(4,Ip) .and. (pos.le.psp(Ity)%numps))
                     pos=pos+1
                  ENDDO
                  dbeta=polynom(1,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,xyzr(4,Ip))
                  endif
!dylm/dr_R
!CALL cal_dylm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip),dYlm)
                  CALL r_dYlm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip),dYlm)
!
                  dbeta_lm(0,Ip,Ipj)   = beta*dYlm(0)
                  dbeta_lm(1:3,Ip,Ipj) = beta*dYlm(1:3)+dYlm(0)*dbeta*xyzr(1:3,Ip)
               ENDDO
!---------------------------------
            ENDDO
!3.the dGlm=\SUM_ik{<dbeta_lm|\psi_ik>} and Glm=\SUM_ik{<beta_lm|\psi_ik>}
            tmpf(:)=0.d0
            DO Ipj=1,psp(Ity)%nproj
!cycle all states n: Is,Ik,Ii
               tmp0(:)=0.d0
               DO Is=1,NSPIN
               DO Ik=1,nk
               DO Ii=1,Nstates
                  Gnlm=0.d0
                  dGnlm(:)=0.d0
                  DO Ip=1,nlpot(Ia)%npts
                     Gnlm=Gnlm + dbeta_lm(0,Ip,Ipj)*psi_ik(Ip,Ii,Ik,Is)
                     dGnlm(:)=dGnlm(:) + dbeta_lm(1:3,Ip,Ipj)*CONJG(psi_ik(Ip,Ii,Ik,Is))
                  ENDDO

                  CALL MPI_ALLREDUCE(Gnlm,Gnlm_global,1,mpi_complex16,mpi_sum,parallel%commx,mpinfo)
                  tmp0(:)=tmp0(:)+wke(parallel%sub2sum(Ii,parallel%ranky+1),Ik,Is)*Gnlm_global*dGnlm(:)

               ENDDO
               ENDDO
               ENDDO
! #ifndef 1
               tmpf(:)=tmpf(:)+psp(Ity)%D0(Ipj,Ipj)*REAL(tmp0(:),DP)
! #else
! CALL MPI_ALLREDUCE(tmp0,tmp0_global,1,mpi_complex16,mpi_sum,parallel%commy,mpinfo)
! tmpf(:)=tmpf(:)+psp(Ity)%D0(Ipj,Ipj)*REAL(tmp0_global(:),DP)
! #endif
            ENDDO
!4.calculate the nonlocal force (dvol for scale)

            nlforce_local(:,Ia)=tmpf(:)*dvol*2.d0
!> communication in comm for the grid and states
            CALL MPI_ALLREDUCE(nlforce_local,nlforce,3*natom,mpi_real8,mpi_sum,parallel%comm,mpinfo)

         ENDDO
      ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nonlforce
!-----------------------DIVIDER-LINE--------------------------
!------------------------nonlocal force--------------------
   SUBROUTINE nonlforce_r_dg(Uik,nlforce)
     USE math , ONLY : SimpleInterp,r_dYlm
     USE nlpot_module , ONLY : nlpot,max_nlnpts
     USE struct_module ,ONLY : naty,struct,recip_lat,natom
     USE pspot_module , ONLY : psp,max_nproj
     USE parameters , ONLY : Nspin,Nstates
     USE grid_module , ONLY : nk,KPT,grid,dvol
     USE smearing_module , ONLY : wke
     USE mathsplines , ONLY : polynom

     USE smpi_math_module , ONLY : parallel ,mpinfo,mpi_real8,mpi_sum

     IMPLICIT NONE
!IN/OUT

     REAL(DP),INTENT(IN) :: Uik(:,:,:) !wave_function

     REAL(DP),INTENT(OUT) :: nlforce(:,:)
!LOCAL
     REAL(DP) :: xyzr(4,max_nlnpts)
     REAL(DP) :: beta,dbeta
!dYlm/dR  and  dvYlm/dR
     REAL(DP) :: dYlm(0:3) &
          & ,    dbeta_lm(0:3,max_nlnpts,max_nproj)
     REAL(DP) :: dGnlm(3),Gnlm,tmp0(3)
!psi_ik=e^ikr U_ik
! #ifndef 1
     REAL(DP)  :: psi_ik(max_nlnpts,Nstates,nk,NSPIN)
! #else
!      REAL(DP)  :: psi_ik(max_nlnpts,parallel%nstate_proc,nk,NSPIN)
! #endif
     REAL(DP)  :: eikr
     REAL(DP)     :: kdotr , kvec(3),rvec(3)
!
     INTEGER(I4B) :: Is,Ik,Ii
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id,l,m
     REAL(DP) :: tmpf(3)
!##ISO_interplote{{
     REAL(DP) :: c(3)
     INTEGER(I4B) :: pos
!##ISO_interplote}}

     REAL(DP) :: Gnlm_global
     REAL(DP) :: nlforce_local(3,natom)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     nlforce(:,:)=0.d0
!Cycle all species
     DO Ity=1,naty
     IF(psp(Ity)%nproj==0) CYCLE
!cycle all atoms
     DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!1.Parpare some data
!(1) Find the nlpot points' distant and cosine
!(2) calculate the \psi=eikr*Uik(r)
!cycle all points around this atom
! psi_ik=0.d0
        DO Ip=1,nlpot(Ia)%npts
!1.store points' data
           xyzr(4,Ip)=nlpot(Ia)%rRvec_iso(4,Ip)
           IF(xyzr(4,Ip)>0.000001d0)THEN
              xyzr(1,Ip)=nlpot(Ia)%rRvec_iso(1,Ip)/xyzr(4,Ip)
              xyzr(2,Ip)=nlpot(Ia)%rRvec_iso(2,Ip)/xyzr(4,Ip)
              xyzr(3,Ip)=nlpot(Ia)%rRvec_iso(3,Ip)/xyzr(4,Ip)
           ELSE
              print*,"atom put on a grid point, xyzr(4,Ip)",xyzr(4,Ip)
              xyzr(1:4,Ip)=0.d0
           ENDIF
!calculate psi_ik at this point
           Id=nlpot(Ia)%Id_iso(Ip)
           rvec(:)=nlpot(Ia)%rRvec_iso(1:3,Ip)+struct%poscar(:,Ia)
           Ik=1

           psi_ik(Ip,:,Ik,:)=Uik(Id,:,:)!eikr*Uik(Id,:,Ik,:)

        ENDDO !> Ip

!2.Build the dbeta_lm/dr for this atom
!cycle all projectors
        DO Ipj=1,psp(Ity)%nproj
!get (l,m)
           l=psp(Ity)%proj_l(Ipj)
           m=psp(Ity)%proj_m(Ipj)
!---------------------------------
           DO Ip=1,nlpot(Ia)%npts
!beta,dbeta
              beta=nlpot(Ia)%proj0_dg(Ip,Ipj)
!##Cal dbeta {{
!##CALCULATE POTENTIAL IN "RR"
              pos=0
              DO WHILE ( psp(Ity)%r_real(pos+1).lt.xyzr(4,Ip) .and. (pos.le.psp(Ity)%numps))
                 pos=pos+1
              ENDDO
              dbeta=polynom(1,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,xyzr(4,Ip))
!##Cal dbeta }}
!dbeta=SimpleInterp( psp(Ity)%dbeta_dr(:,Ipj),psp(Ity)%rmax, &
!             & psp(Ity)%rspacing,xyzr(4,Ip) )
!dylm/dr_R
!CALL cal_dylm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip),dYlm)
              CALL r_dYlm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip),dYlm)
!
              dbeta_lm(0,Ip,Ipj)   = beta*dYlm(0)
              dbeta_lm(1:3,Ip,Ipj) = beta*dYlm(1:3)+dYlm(0)*dbeta*xyzr(1:3,Ip)
           ENDDO
!---------------------------------
        ENDDO
!3.the dGlm=\SUM_ik{<dbeta_lm|\psi_ik>} and Glm=\SUM_ik{<beta_lm|\psi_ik>}
        tmpf(:)=0.d0
        DO Ipj=1,psp(Ity)%nproj
!cycle all states n: Is,Ik,Ii
           tmp0(:)=0.d0
           DO Is=1,NSPIN
           DO Ik=1,nk
! #ifndef 1
           DO Ii=1,Nstates
! #else
! DO Ii=1,parallel%nstate_proc
! #endif
              Gnlm=0.d0
              dGnlm(:)=0.d0
              DO Ip=1,nlpot(Ia)%npts
                 Gnlm=Gnlm + dbeta_lm(0,Ip,Ipj)*psi_ik(Ip,Ii,Ik,Is)
                 dGnlm(:)=dGnlm(:) + dbeta_lm(1:3,Ip,Ipj)*psi_ik(Ip,Ii,Ik,Is)
              ENDDO

              CALL MPI_ALLREDUCE(Gnlm,Gnlm_global,1,mpi_real8,mpi_sum,parallel%comm,mpinfo)
              tmp0(:)=tmp0(:)+wke(Ii,Ik,Is)*Gnlm_global*dGnlm(:)
! tmp0(:)=tmp0(:)+wke(parallel%sub2sum(Ii,parallel%myid+1),Ik,Is)*Gnlm*dGnlm(:)

           ENDDO
           ENDDO
           ENDDO
           tmpf(:)=tmpf(:)+psp(Ity)%D0(Ipj,Ipj)*REAL(tmp0(:),DP)
        ENDDO
!4.calculate the nonlocal force (dvol for scale)

        nlforce_local(:,Ia)=tmpf(:)*dvol*2.d0*2.d0
        CALL MPI_ALLREDUCE(nlforce_local,nlforce,3*natom,mpi_real8,mpi_sum,parallel%comm,mpinfo)

      ENDDO
      ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE nonlforce_r_dg
!------------------------nonlocal force--------------------
   SUBROUTINE nonlforce_r(Uik,nlforce)
     USE math , ONLY : SimpleInterp,r_dYlm
     USE nlpot_module , ONLY : nlpot,max_nlnpts
     USE struct_module ,ONLY : naty,struct,recip_lat,natom
     USE pspot_module , ONLY : psp,max_nproj
     USE parameters , ONLY : Nspin,Nstates
     USE grid_module , ONLY : nk,KPT,grid,dvol
     USE smearing_module , ONLY : wke
     USE mathsplines , ONLY : polynom

     USE smpi_math_module , ONLY : parallel ,mpinfo,mpi_real8,mpi_sum

     USE m_time_evaluate, ONLY: memory_sum,memory_free
     IMPLICIT NONE
!IN/OUT

     REAL(DP),INTENT(IN) :: Uik(:,:,:) !wave_function

     REAL(DP),INTENT(OUT) :: nlforce(:,:)
!LOCAL
     REAL(DP) :: xyzr(4,max_nlnpts)
     REAL(DP) :: beta,dbeta
!dYlm/dR  and  dvYlm/dR
     REAL(DP) :: dYlm(0:3) &
          & ,    dbeta_lm(0:3,max_nlnpts,max_nproj)
     REAL(DP) :: dGnlm(3),Gnlm,tmp0(3)
!psi_ik=e^ikr U_ik
! #ifndef 1
     REAL(DP)  :: psi_ik(max_nlnpts,Nstates,nk,NSPIN)
! #else
!      REAL(DP)  :: psi_ik(max_nlnpts,parallel%nstate_proc,nk,NSPIN)
! #endif
     REAL(DP)  :: eikr
     REAL(DP)     :: kdotr , kvec(3),rvec(3)
!
     INTEGER(I4B) :: Is,Ik,Ii
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id,l,m
     REAL(DP) :: tmpf(3)
!##ISO_interplote{{
     REAL(DP) :: c(3)
     INTEGER(I4B) :: pos
!##ISO_interplote}}

     REAL(DP) :: Gnlm_global
     REAL(DP) :: nlforce_local(3,natom)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     call memory_sum('nolocal_force',(real(size(xyzr),DP)+size(psi_ik)+size(nlforce_local))*DP)

     nlforce(:,:)=0.d0
!Cycle all species
     DO Ity=1,naty
     IF(psp(Ity)%nproj==0) CYCLE
!cycle all atoms
     DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!1.Parpare some data
!(1) Find the nlpot points' distant and cosine
!(2) calculate the \psi=eikr*Uik(r)
!cycle all points around this atom
! psi_ik=0.d0
        DO Ip=1,nlpot(Ia)%npts
!1.store points' data
           xyzr(4,Ip)=nlpot(Ia)%rRvec_iso(4,Ip)
           IF(xyzr(4,Ip)>0.000001d0)THEN
              xyzr(1,Ip)=nlpot(Ia)%rRvec_iso(1,Ip)/xyzr(4,Ip)
              xyzr(2,Ip)=nlpot(Ia)%rRvec_iso(2,Ip)/xyzr(4,Ip)
              xyzr(3,Ip)=nlpot(Ia)%rRvec_iso(3,Ip)/xyzr(4,Ip)
           ELSE
              print*,"atom put on a grid point, xyzr(4,Ip)",xyzr(4,Ip)
              xyzr(1:4,Ip)=0.d0
           ENDIF
!calculate psi_ik at this point
           Id=nlpot(Ia)%Id_iso(Ip)
           rvec(:)=nlpot(Ia)%rRvec_iso(1:3,Ip)+struct%poscar(:,Ia)
           Ik=1

           psi_ik(Ip,:,Ik,:)=Uik(Id,:,:)!eikr*Uik(Id,:,Ik,:)

        ENDDO !> Ip

!2.Build the dbeta_lm/dr for this atom
!cycle all projectors
        DO Ipj=1,psp(Ity)%nproj
!get (l,m)
           l=psp(Ity)%proj_l(Ipj)
           m=psp(Ity)%proj_m(Ipj)
!---------------------------------
           DO Ip=1,nlpot(Ia)%npts
!beta,dbeta
              beta=nlpot(Ia)%proj0(Ip,Ipj)
!##Cal dbeta {{
!##CALCULATE POTENTIAL IN "RR"
              pos=0
              DO WHILE ( psp(Ity)%r_real(pos+1).lt.xyzr(4,Ip) .and. (pos.le.psp(Ity)%numps))
                 pos=pos+1
              ENDDO
              dbeta=polynom(1,3,psp(Ity)%r_real(pos-1:pos+1),psp(Ity)%beta_r(pos-1:pos+1,Ipj),c,xyzr(4,Ip))
!##Cal dbeta }}
!dbeta=SimpleInterp( psp(Ity)%dbeta_dr(:,Ipj),psp(Ity)%rmax, &
!             & psp(Ity)%rspacing,xyzr(4,Ip) )
!dylm/dr_R
!CALL cal_dylm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip),dYlm)
              CALL r_dYlm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip),dYlm)
!
              dbeta_lm(0,Ip,Ipj)   = beta*dYlm(0)
              dbeta_lm(1:3,Ip,Ipj) = beta*dYlm(1:3)+dYlm(0)*dbeta*xyzr(1:3,Ip)
           ENDDO
!---------------------------------
        ENDDO
!3.the dGlm=\SUM_ik{<dbeta_lm|\psi_ik>} and Glm=\SUM_ik{<beta_lm|\psi_ik>}
        tmpf(:)=0.d0
        DO Ipj=1,psp(Ity)%nproj
!cycle all states n: Is,Ik,Ii
           tmp0(:)=0.d0
           DO Is=1,NSPIN
           DO Ik=1,nk
! #ifndef 1
           DO Ii=1,Nstates
! #else
! DO Ii=1,parallel%nstate_proc
! #endif
              Gnlm=0.d0
              dGnlm(:)=0.d0
              DO Ip=1,nlpot(Ia)%npts
                 Gnlm=Gnlm + dbeta_lm(0,Ip,Ipj)*psi_ik(Ip,Ii,Ik,Is)
                 dGnlm(:)=dGnlm(:) + dbeta_lm(1:3,Ip,Ipj)*psi_ik(Ip,Ii,Ik,Is)
              ENDDO

              CALL MPI_ALLREDUCE(Gnlm,Gnlm_global,1,mpi_real8,mpi_sum,parallel%comm,mpinfo)
              tmp0(:)=tmp0(:)+wke(Ii,Ik,Is)*Gnlm_global*dGnlm(:)
! tmp0(:)=tmp0(:)+wke(parallel%sub2sum(Ii,parallel%myid+1),Ik,Is)*Gnlm*dGnlm(:)

           ENDDO
           ENDDO
           ENDDO
           tmpf(:)=tmpf(:)+psp(Ity)%D0(Ipj,Ipj)*REAL(tmp0(:),DP)
        ENDDO
!4.calculate the nonlocal force (dvol for scale)

        nlforce_local(:,Ia)=tmpf(:)*dvol*2.d0!*2.d0
        CALL MPI_ALLREDUCE(nlforce_local,nlforce,3*natom,mpi_real8,mpi_sum,parallel%comm,mpinfo)

      ENDDO
      ENDDO

     call memory_free('nolocal_force',(real(size(xyzr),DP)+size(psi_ik)+size(nlforce_local))*DP)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE nonlforce_r
!-----------------------PARTING-LINE--------------------------
!-----------------------cal_stress-------------------------
SUBROUTINE cal_stress(rhoS,Uik,stress)
   USE struct_module, ONLY : struct,lat_mat
   USE ewald , ONLY : ewald_stress
   USE grid_module , ONLY: n1,n2,n3,ng1,ng2,ng3
   USE parameters , ONLY : Nspin
   USE FOURIER
   USE pspot_module , ONLY : psp,max_nproj
   USE potential_module , ONLY : vlda
   USE energy_module , ONLY : Exc,Ehart,Eie

   USE smpi_math_module,  ONLY: parallel,mpinfo,mpi_real8
   USE grid_module, ONLY: global_n1,global_n2,global_n3

   IMPLICIT NONE
!IN/OUT
   REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
   COMPLEX(DCP),INTENT(IN) :: Uik(:,:,:,:)
   REAL(DP),INTENT(OUT) :: stress(3,3)
!LOCAL
   REAL(DP) :: tmps(3,3)
   REAL(DP),DIMENSION(n1,n2,n3) :: rho,vxc

   COMPLEX(DCP) :: rhoq(ng1,ng2,parallel%local_z)

   INTEGER(I4B) :: Is,i,j

   REAL(DP) :: rho_FFT(global_n1,global_n2,parallel%local_z)
   INTEGER(I4B) :: ix,iy

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!total rho
   rho=0.d0
   DO Is=1,Nspin
      rho(:,:,:)=rho(:,:,:)+rhoS(:,:,:,Is)
   ENDDO

   ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
   iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
   CALL MPI_ALLTOALLV(rho(ix,iy,1),parallel%fft_scount&
        & ,parallel%fft_sdispls,MPI_REAL8&
        & ,rho_FFT,parallel%fft_rcount&
        & ,parallel%fft_rdispls,MPI_REAL8&
        & ,parallel%commx,mpinfo)
! call MPI_ALLGATHERV(rho(ix,iy,1),parallel%mygrid_range(3)&
!      &,MPI_REAL8,rho_global,parallel%recvcounts&
!      &,parallel%displs,MPI_REAL8&
!      & ,parallel%commx,mpinfo)

!rhoq

   rhoq = FFT(rho_FFT)

!cal_stress
   stress(:,:)=0.d0
!Ewald stress
   tmps(:,:)=ewald_stress(lat_mat,struct%pos,struct%eleid,psp(:)%Zion)
   stress(:,:)=stress(:,:) + tmps(:,:)
! if(parallel%isroot)print*,'ewald stress'
! if(parallel%isroot)print*,tmps
!Ion-e local part stress
   tmps(:,:)=IonEle_stress(rhoq,Eie)
   stress(:,:)=stress(:,:) + tmps(:,:)
! if(parallel%isroot)print*,'local stress'
! if(parallel%isroot)print*,tmps
!nonlocal stress
   if(max_nproj/=0)then
   CALL ion_nl_stress(Uik,tmps)
   stress(:,:)=stress(:,:) + tmps(:,:)
   endif
! if(parallel%isroot)print*,'nonlocal stress'
! if(parallel%isroot)print*,tmps
!kinetic stress
   CALL Kin_stress(Uik,tmps)
   stress(:,:)=stress(:,:) + tmps(:,:)
! if(parallel%isroot)print*,'Kinetic stress'
! if(parallel%isroot)print*,tmps
!XC
   CALL vlda(rho,vxc)
   tmps(:,:)=LDA_stress(vxc,rho,Exc)
   stress(:,:)=stress(:,:) + tmps(:,:)
! if(parallel%isroot)print*,'LDA stress'
! if(parallel%isroot)print*,tmps
!Hartree stress
   tmps(:,:)=Hart_stress(rhoq,Ehart)
   stress(:,:)=stress(:,:) + tmps(:,:)
! if(parallel%isroot)print*,'Hart stress'
! if(parallel%isroot)print*,tmps
!for total stress
   DO i=1,3
   DO j=i+1,3
      stress(j,i)=stress(i,j)
   ENDDO
   ENDDO
!total stress
   cellpress=-(stress(1,1)+stress(2,2)+stress(3,3))/3*au2gpa
!print*,'Hartree stress'
!print*,tmps

   if(parallel%isroot)Then

   print*,'[Cell Lattice parameters (ang)]'
   WRITE(6,*) '--------------------------------------------'
   WRITE(6,'(3(E15.7))') lat_mat(:,1)*bohr2ang
   WRITE(6,'(3(E15.7))') lat_mat(:,2)*bohr2ang
   WRITE(6,'(3(E15.7))') lat_mat(:,3)*bohr2ang
   WRITE(6,*) '--------------------------------------------'
   print*,'[Total stress (GPa)]'
   WRITE(6,*) '--------------------------------------------'
   WRITE(6,'(3(E15.7))') stress(:,1) * au2gpa
   WRITE(6,'(3(E15.7))') stress(:,2) * au2gpa
   WRITE(6,'(3(E15.7))') stress(:,3) * au2gpa
   WRITE(6,*) '--------------------------------------------'
   print*,'[Total press] (Gpa)',cellpress
   WRITE(6,*) '--------------------------------------------'

   endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE cal_stress
!--------------------LDA XC stress-------------------------
   FUNCTION LDA_stress(vxc,rho,Elda)
      USE struct_module , ONLY : volume

      USE smpi_math_module, ONLY: mpi_real8,mpi_sum,parallel,mpinfo
      USE grid_module ,ONLY: global_n,n,n1,n2,n3,rho_trans1D
! USE Math ,ONLY: kahan_sum

      IMPLICIT NONE
!IN/OUT
      REAL(DP),INTENT(IN) :: rho(:,:,:)
      REAL(DP),INTENT(IN) :: vxc(:,:,:)
      REAL(DP),INTENT(IN) :: Elda
      REAL(DP)   :: LDA_stress(3,3)
!LOCAL
      INTEGER(I4B)   :: i,j
      real(DP)                     :: sum_R

      real(DP) :: sum_R_local,temp(n),tempxc(n)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!

      call rho_trans1D(rho,temp)
      call rho_trans1D(vxc,tempxc)
      sum_R_local=SUM(temp*tempxc)/global_n
! temp=temp*tempxc
! sum_R_local=kahan_sum(size(temp),temp)/global_n
      CALL MPI_ALLREDUCE(sum_R_local,sum_R,1,mpi_real8,mpi_sum,parallel%commx,mpinfo)

      do i = 1, 3
         do j = i, 3
            if (i==j) then
               LDA_stress(i,j) = Elda / volume - sum_R
            else
               LDA_stress(i,j) = 0.d0
            end if
         end do
      end do
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDFUNCTION LDA_stress
!----------------------------------------------------------
   FUNCTION Hart_stress(rhoRecip,Eh)
      USE struct_module,  ONLY : volume,recip_lat
      USE grid_module , ONLY : ng1,ng2,ng3,grid

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

      IMPLICIT NONE
!IN/OUT
      REAL(DP),INTENT(IN)    :: Eh
      COMPLEX(DCP),INTENT(IN) :: rhoRecip(:,:,:)
      REAL(DP)               :: Hart_stress(3,3)
!LOCAL
      INTEGER(I4B)                 :: i,j,ix,iy,iz,Ig
      real(DP)                     :: q(3)
      real(DP)                     :: temp
      REAL(DP)                     :: J_temp(3,3)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Hart_stress = 0.d0
      J_temp = 0.d0
      Ig=0

      DO iz = 1, parallel%local_z
         Ig=(parallel%local_z_start+iz-1)*ng2*ng1

      DO iy = 1, ng2
      DO ix = 1, ng1
         Ig=Ig+1
         IF(grid%gMask(Ig))THEN
            q = -1.0_DP * grid%gVec(1:3,Ig)
! Now compute the product of the components of interest.
            temp = CONJG(rhoRecip(ix,iy,iz)) * rhoRecip(ix,iy,iz) / grid%gVec(4,Ig)**4
            DO i=1,3
               DO j=i,3
                  J_temp(i,j) = J_temp(i,j) + temp * q(i) * q(j)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!

      Hart_stress= Hart_stress + J_temp
      CALL MPI_ALLREDUCE(J_temp,Hart_stress,9,mpi_REAL8,mpi_sum,parallel%commx,mpinfo)


      Hart_stress = Hart_stress * 8._DP * pi

      DO i = 1,3
         DO j = i,3
            IF (i==j) Hart_stress(i,j) = Hart_stress(i,j) - Eh / volume
         END DO
      END DO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDFUNCTION Hart_stress
!---------------------kinetic stress-----------------------
   SUBROUTINE Kin_stress(Uik,Kinstress)
      USE parameters , ONLY : Nspin,Nstates,norder=>finite_order
      USE grid_module , ONLY : n1,n2,n3,n,nk,grid,KPT
      USE struct_module , ONLY : volume
      USE finite_module , ONLY : nabla
      USE smearing_module , ONLY : wke

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

      IMPLICIT NONE
!IN/OUT
      COMPLEX(DCP),INTENT(IN) :: Uik(:,:,:,:)
      REAL(DP),INTENT(OUT)   :: Kinstress(3,3)
!LOCAL
      COMPLEX(DCP)  :: dUik(3,n1,n2,n3),tUik(n1,n2,n3)
      INTEGER(I4B) :: Is,Ik,Ii,Ix,Iy,Iz,I

      INTEGER(I4B) :: in,j,k
      REAL(DP) :: Kinstress_local(3,3)
      INTEGER(I4B) :: d_z,y_down

      REAL(DP) :: kdotr,kvec(3),kk(6)
      COMPLEX(DCP) :: tmps(6),tmp0(6) &
                 &,sum_u2,sum_udu(3),sum_dudu(6)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Kinstress(:,:)=0.d0
      tmps(:)=0.d0
!cycle spin
      DO Is=1,Nspin
!cycle k-points
         DO Ik=1,nk
!parpare k index
            kvec(:)=KPT%vcar(:,Ik)
!kk=k_i * k_j
            kk(1)=kvec(1)*kvec(1)
            kk(2)=kvec(1)*kvec(2)
            kk(3)=kvec(1)*kvec(3)
            kk(4)=kvec(2)*kvec(2)
            kk(5)=kvec(2)*kvec(3)
            kk(6)=kvec(3)*kvec(3)
!cycle all states
            DO Ii=1,Nstates
!store psi=e^ikr*Uik(r)

               tUik=cmplx(0.d0,0.d0)
               ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
               iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
               iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
               in=0 !parallel%mygrid_range(1)-1
               do i=ix,n1
                  in=in+1
                  tUik(i,iy,iz)=Uik(in,Ii,Ik,Is)
               enddo
               d_z=ceiling((n3-iz)/(n3-iz+1.d0))
               do j=iy+1,n2*d_z
                  do i=1,n1
                     in=in+1
                     tUik(i,j,iz)= Uik(in,Ii,Ik,Is)
                  enddo
               enddo
               do k=2,n3-1
                  do j=1,n2
                     do i=1,n1
                        in=in+1
                        tUik(i,j,k)= Uik(in,Ii,Ik,Is)
                     enddo
                  enddo
               enddo
               ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
               y_down=(1-iy-1)*d_z+iy+1
               iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
               iz=n3 !(parallel%mygrid_range(2)-1)/n1/n2+1
               do j=y_down,iy-1
                  do i=1,n1
                     in=in+1
                     tUik(i,j,iz)= Uik(in,Ii,Ik,Is)
                  enddo
               enddo
               do i=1,ix
                  in=in+1
                  tUik(i,iy,iz)= Uik(in,Ii,Ik,Is)
               enddo

!calculate the gridient of psi
               CALL nabla(tUik,norder,dUik)
!store something
               sum_u2 = SUM(CONJG(tUik)*tUik)
!udui
               sum_udu(1)=SUM(CONJG(tUik(:,:,:))*dUik(1,:,:,:))
               sum_udu(2)=SUM(CONJG(tUik(:,:,:))*dUik(2,:,:,:))
               sum_udu(3)=SUM(CONJG(tUik(:,:,:))*dUik(3,:,:,:))
!dudu
               sum_dudu(1)=SUM(CONJG(dUik(1,:,:,:))*dUik(1,:,:,:))
               sum_dudu(2)=SUM(CONJG(dUik(1,:,:,:))*dUik(2,:,:,:))
               sum_dudu(3)=SUM(CONJG(dUik(1,:,:,:))*dUik(3,:,:,:))
               sum_dudu(4)=SUM(CONJG(dUik(2,:,:,:))*dUik(2,:,:,:))
               sum_dudu(5)=SUM(CONJG(dUik(2,:,:,:))*dUik(3,:,:,:))
               sum_dudu(6)=SUM(CONJG(dUik(3,:,:,:))*dUik(3,:,:,:))
!calculate the temp stress
!1:11
               tmp0(1)=  sum_u2*kk(1) + sum_dudu(1)      &
                   &   + IMAG*kvec(1)*CONJG(sum_udu(1))  &
                   &   - IMAG*kvec(1)*sum_udu(1)
!2:12
               tmp0(2)=  sum_u2*kk(2) + sum_dudu(2)      &
                   &   + IMAG*kvec(2)*CONJG(sum_udu(1))  &
                   &   - IMAG*kvec(1)*sum_udu(2)
!3:13
               tmp0(3)=  sum_u2*kk(3) + sum_dudu(3)      &
                   &   + IMAG*kvec(3)*CONJG(sum_udu(1))  &
                   &   - IMAG*kvec(1)*sum_udu(3)
!4:22
               tmp0(4)=  sum_u2*kk(4) + sum_dudu(4)      &
                   &   + IMAG*kvec(2)*CONJG(sum_udu(2))  &
                   &   - IMAG*kvec(2)*sum_udu(2)
!5:23
               tmp0(5)=  sum_u2*kk(5) + sum_dudu(5)      &
                   &   + IMAG*kvec(3)*CONJG(sum_udu(2))  &
                   &   - IMAG*kvec(2)*sum_udu(3)
!6:33
               tmp0(6)=  sum_u2*kk(6) + sum_dudu(6)      &
                   &   + IMAG*kvec(3)*CONJG(sum_udu(3))  &
                   &   - IMAG*kvec(3)*sum_udu(3)
!sum

               tmps(:) = tmps(:) + tmp0(:)*wke(parallel%sub2sum(Ii,parallel%ranky+1),Ik,Is)

            ENDDO

         ENDDO
      ENDDO
!calculate the stress
      tmps(:)=tmps(:)/volume
!out
      Kinstress(1,1)=REAL(tmps(1),DP)
      Kinstress(1,2)=REAL(tmps(2),DP)
      Kinstress(2,1)=REAL(tmps(2),DP)
      Kinstress(1,3)=REAL(tmps(3),DP)
      Kinstress(3,1)=REAL(tmps(3),DP)
      Kinstress(2,2)=REAL(tmps(4),DP)
      Kinstress(2,3)=REAL(tmps(5),DP)
      Kinstress(3,2)=REAL(tmps(5),DP)
      Kinstress(3,3)=REAL(tmps(6),DP)
!

      Kinstress_local(:,:)=-Kinstress(:,:)
      CALL MPI_ALLREDUCE(Kinstress_local,kinstress,9,mpi_real8,mpi_sum,parallel%comm,mpinfo)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Kin_stress
!-----------------------nl stress--------------------------
   SUBROUTINE ion_nl_stress(Uik,nlstress)
      USE parameters , ONLY : Nstates,Nspin,PP_identifer
      USE nlpot_module , ONLY : max_nlnpts,nlpot
      USE pspot_module , ONLY : max_nproj,psp
      USE struct_module , ONLY : naty,struct
      USE grid_module , ONLY : nk,n,KPT,grid
      USE math , ONLY : SimpleInterp,r_dYlm
      USE MathSplines , ONLY: polynom
      USE smearing_module , ONLY : wke

      USE smpi_math_module, ONLY: mpi_complex16&
           & ,mpinfo,parallel,mpi_sum,mpi_real8
      USE grid_module, ONLY: global_n

      IMPLICIT NONE
!IN/OUT
      COMPLEX(DCP),INTENT(IN) :: Uik(:,:,:,:)
      REAL(DP),INTENT(OUT)   :: nlstress(3,3)
!LOCAL
      REAL(DP),DIMENSION(0:3,max_nlnpts,max_nproj) :: &
           & dYlm      !sphereic harmonic function
!&  dbetalm & !dbeta
      REAL(DP) :: dbeta(0:1,max_nlnpts,max_nproj)
      REAL(DP) :: xyzr(4,max_nlnpts) &!cosine and distent
               &, rvec(3)
!psi_ik=e^ikr U_ik
      REAL(DP) :: kdotr
      COMPLEX(DCP)  :: psi_ik(max_nlnpts,Nstates,nk,Nspin)
!index
      INTEGER(I4B) :: Ity,Ia,Ip,Id,Is,Ik,Ii,Ipj,l,m
      COMPLEX(DCP)  :: eikr,tmps(6),tmp0(6),tG0,tGs

      COMPLEX(DCP)  :: tGs_global

      COMPLEX(DCP)  :: dGnlm(6),Gnlm,coe(6),tylm(6)   !

      COMPLEX(DCP)  :: dGnlm_local(6),Gnlm_local,tmps_local(6)
      REAL(DP)     :: nlstress_local(3,3)

!> interplote
      REAL(DP) :: c(3)
      INTEGER(I4B) :: pos
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!cycle all species
      tGs=0.d0
      tmps(:)=0.d0
      DO Ity=1,naty
         IF(psp(Ity)%nproj==0) CYCLE
!cycle all atom
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
!1.parpare the necessary data
!(1)find the nlpot points
!(2)\psi=eikr*Uik(r)
!cycle all points
            DO Ip=1,nlpot(Ia)%npts
!(1)store the distents
               xyzr(:,Ip)=nlpot(Ia)%rRvec(:,Ip)
               IF(xyzr(4,Ip)>0.d0)THEN
                  xyzr(1,Ip)=xyzr(1,Ip)/xyzr(4,Ip)
                  xyzr(2,Ip)=xyzr(2,Ip)/xyzr(4,Ip)
                  xyzr(3,Ip)=xyzr(3,Ip)/xyzr(4,Ip)
               ELSE
                  xyzr(:,Ip)=0.d0
               ENDIF
!(2)psi
               Id=nlpot(Ia)%Id(Ip)
               rvec(:)=nlpot(Ia)%rRvec(1:3,Ip)+struct%poscar(:,Ia)
               DO Ik=1,nk
!eikr
                  kdotr=DOT_PRODUCT(KPT%vcar(:,Ik),rvec(:))
                  eikr=EXP(IMAG*kdotr)
                  psi_ik(Ip,:,Ik,:)=eikr*Uik(Id,:,Ik,:)
               ENDDO

            ENDDO
!2.Build dbeta/dr ,dYlm
            DO Ipj=1,psp(Ity)%nproj
!get (l,m)
               l=psp(Ity)%proj_l(Ipj)
               m=psp(Ity)%proj_m(Ipj)
!----------------------------------
               DO Ip=1,nlpot(Ia)%npts
!dbeta
                  dbeta(0,Ip,Ipj)  = nlpot(Ia)%proj0(Ip,Ipj)
                  if(PP_identifer==0)then
                  dbeta(1,Ip,Ipj)  = SimpleInterp(         &
                              &  psp(Ity)%dbeta_dr(:,Ipj), &
                              &  psp(Ity)%rmax,            &
                              &  psp(Ity)%rspacing,        &
                              &  xyzr(4,Ip) )
                  else
                  pos=0
                  DO WHILE ( psp(Ity)%r_real(pos+1).lt.xyzr(4,Ip)&
                       & .and. (pos.le.psp(Ity)%numps))
                     pos=pos+1
                  ENDDO
                  dbeta(1,Ip,Ipj)=polynom(1,3&
                       &,psp(Ity)%r_real(pos-1:pos+1)&
                       &,psp(Ity)%beta_r(pos-1:pos+1,Ipj)&
                       &,c,xyzr(4,Ip))
                  endif
!dYlm
                  CALL r_dYlm(l,m,xyzr(1,Ip),xyzr(2,Ip),xyzr(3,Ip),xyzr(4,Ip)  &
                           &   , dYlm(:,Ip,Ipj) )
               ENDDO
!----------------------------------
            ENDDO
!3.calculate
            DO Ipj=1,psp(Ity)%nproj
               tG0=0.d0
               tmp0(:)=0.d0
!cycle all states n: Is,Ik,Ii
               DO Is=1,Nspin
               DO Ik=1,nk
               DO Ii=1,Nstates
                  Gnlm=0.d0
                  dGnlm(:)=0.d0
                  DO Ip=1,nlpot(Ia)%npts
!coe for dbeta
                     coe(1)=xyzr(1,Ip)*xyzr(1,Ip)*xyzr(4,Ip)
                     coe(2)=xyzr(1,Ip)*xyzr(2,Ip)*xyzr(4,Ip)
                     coe(3)=xyzr(1,Ip)*xyzr(3,Ip)*xyzr(4,Ip)
                     coe(4)=xyzr(2,Ip)*xyzr(2,Ip)*xyzr(4,Ip)
                     coe(5)=xyzr(2,Ip)*xyzr(3,Ip)*xyzr(4,Ip)
                     coe(6)=xyzr(3,Ip)*xyzr(3,Ip)*xyzr(4,Ip)
!for dYlm
                     tylm(1)=xyzr(4,Ip)*xyzr(1,Ip)*dYlm(1,Ip,Ipj)
                     tylm(2)=0.5d0*xyzr(4,Ip)*(xyzr(1,Ip)*dYlm(2,Ip,Ipj)+xyzr(2,Ip)*dYlm(1,Ip,Ipj))
                     tylm(3)=0.5d0*xyzr(4,Ip)*(xyzr(1,Ip)*dYlm(3,Ip,Ipj)+xyzr(3,Ip)*dYlm(1,Ip,Ipj))
                     tylm(4)=xyzr(4,Ip)*xyzr(2,Ip)*dYlm(2,Ip,Ipj)
                     tylm(5)=0.5d0*xyzr(4,Ip)*(xyzr(2,Ip)*dYlm(3,Ip,Ipj)+xyzr(3,Ip)*dYlm(2,Ip,Ipj))
                     tylm(4)=xyzr(4,Ip)*xyzr(3,Ip)*dYlm(3,Ip,Ipj)
!Gnlm,dGnlm/dE_ij
                     Gnlm=Gnlm+dbeta(0,Ip,Ipj)*dYlm(0,Ip,Ipj)*psi_ik(Ip,Ii,Ik,Is)
!d
                     dGnlm(1)=dGnlm(1)+( coe(1)*dbeta(1,Ip,Ipj)*dYlm(0,Ip,Ipj)  &
                           & + dbeta(0,Ip,Ipj)*tylm(1)  ) * psi_ik(Ip,Ii,Ik,Is)
                     dGnlm(2)=dGnlm(2)+( coe(2)*dbeta(1,Ip,Ipj)*dYlm(0,Ip,Ipj)  &
                           & + dbeta(0,Ip,Ipj)*tylm(2)  ) * psi_ik(Ip,Ii,Ik,Is)
                     dGnlm(3)=dGnlm(3)+( coe(3)*dbeta(1,Ip,Ipj)*dYlm(0,Ip,Ipj)  &
                           & + dbeta(0,Ip,Ipj)*tylm(3)  ) * psi_ik(Ip,Ii,Ik,Is)
                     dGnlm(4)=dGnlm(4)+( coe(4)*dbeta(1,Ip,Ipj)*dYlm(0,Ip,Ipj)  &
                           & + dbeta(0,Ip,Ipj)*tylm(4)  ) * psi_ik(Ip,Ii,Ik,Is)
                     dGnlm(5)=dGnlm(5)+( coe(5)*dbeta(1,Ip,Ipj)*dYlm(0,Ip,Ipj)  &
                           & + dbeta(0,Ip,Ipj)*tylm(5)  ) * psi_ik(Ip,Ii,Ik,Is)
                     dGnlm(6)=dGnlm(6)+( coe(6)*dbeta(1,Ip,Ipj)*dYlm(0,Ip,Ipj)  &
                           & + dbeta(0,Ip,Ipj)*tylm(6)  ) * psi_ik(Ip,Ii,Ik,Is)
                  ENDDO

                  Gnlm_local=Gnlm
                  CALL MPI_ALLREDUCE(Gnlm_local,Gnlm,1,mpi_complex16,mpi_sum,parallel%commx,mpinfo)
! dGnlm_local(:)=dGnlm(:)
! CALL MPI_ALLREDUCE(dGnlm_local,dGnlm,6,mpi_complex16,mpi_sum,parallel%commx,mpinfo)
                  tmp0(:)=tmp0(:) + wke(parallel%sub2sum(Ii,parallel%ranky+1),Ik,Is)*CONJG(Gnlm)*dGnlm(:)
                  tG0 = tG0 + wke(parallel%sub2sum(Ii,parallel%ranky+1),Ik,Is)*CONJG(Gnlm)*Gnlm

               ENDDO
               ENDDO
               ENDDO
               tmps(:)=tmps(:) + psp(Ity)%D0(Ipj,Ipj)*tmp0
               tGs = tGs + psp(Ity)%D0(Ipj,Ipj)*tG0
            ENDDO

         ENDDO
      ENDDO

      tmps_local=tmps
      CALL MPI_ALLREDUCE(tmps_local,tmps,6,mpi_complex16,mpi_sum,parallel%comm,mpinfo)
!> tGs need to be reduced in commy while tmps have reduced in comm
      CALL MPI_ALLREDUCE(tGs,tGs_global,6,mpi_complex16,mpi_sum,parallel%commy,mpinfo)
      tGs=tGs_global

!out the stress
      nlstress(1,1)=REAL( tGs+2*tmps(1), DP )
      nlstress(1,2)=REAL(     2*tmps(2), DP )
      nlstress(1,3)=REAL(     2*tmps(3), DP )
      nlstress(2,2)=REAL( tGs+2*tmps(4), DP )
      nlstress(2,3)=REAL(     2*tmps(5), DP )
      nlstress(3,3)=REAL( tGs+2*tmps(6), DP )
      nlstress(2,1)=nlstress(1,2)
      nlstress(3,1)=nlstress(1,3)
      nlstress(3,2)=nlstress(2,3)
!

! print *,'what',nlstress
      nlstress_local=nlstress
      nlstress_local=nlstress_local/real(global_n)
      nlstress=nlstress_local
! nlstress(:,:)=nlstress(:,:)/REAL(global_n,DP)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ion_nl_stress
!----------------------local stress------------------------
   function IonEle_stress(rhoRecip,energy)!{{{
      USE struct_module , ONLY : struct,naty,volume
      USE grid_module , ONLY : grid,ng
      USE pspot_module , ONLY : psp

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

!$ use omp_lib
      implicit none
!type (struct_type)           :: struct
      real(DP),intent(in)          :: energy
      complex(DCP),dimension(:,:,:),intent(in) :: rhoRecip
      INTEGER(I4B)                 :: sizex,sizey,sizez
      INTEGER(I4B)                 :: i,j,ix,iy,iz,Ig
      real(DP)                     :: IonEle_stress(3,3)
      real(DP)                     :: qPoint(3)
      real(DP)                     :: qNorm
      complex(DCP)                  :: pspPart
      complex(DCP)                  :: tepc
      REAL(DP)                     :: IE_stress(3,3)
!>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      sizeX = SIZE(rhoRecip, 1)
      sizeY = SIZE(rhoRecip, 2)
      sizeZ = SIZE(rhoRecip, 3)
      IonEle_stress = 0.d0
      IE_stress = 0.d0
!PRINT *,"rho111"
!PRINT *,rhoRecip(1:3,1,1)
!$OMP PARALLEL FIRSTPRIVATE(IE_stress)
!$OMP DO PRIVATE(iz,iy,ix,qPoint,qNorm,i,pspPart,j,tepc)
      Ig=0
      do iz = 1, sizeZ

         Ig=(parallel%local_z_start+Iz-1)*sizeY*sizeX

      do iy = 1, sizeY
      do ix = 1, sizeX
         Ig=Ig+1
         if(grid%gMask(Ig)) then
!qPoint = -1.0_DP * qVectors(ix, iy, iz, :)
            qPoint = grid%gVec(1:3,Ig)
            qNorm =  grid%gVec(4,Ig)
!PRINT *,"ii",ix,iy,iz
!PRINT *,"q",qPoint
!PRINT *,"qNorm",qNorm
            pspPart = (0._DP, 0._DP)
            do i=1,naty
               tepc=(0.d0,0.d0)
               do j=struct%eleid(i),struct%eleid(i+1)-1
                  tepc= tepc+ EXP(IMAG*CMPLX(DOT_PRODUCT(qPoint,struct%poscar(:,j)),0.0_DP))
               ENDDO
               pspPart = pspPart + tepc * PseudoPotDiffLookup(i,qNorm)
!PRINT *,"tepc",tepc
!PRINT *,"psu",PseudoPotDiffLookup(pp,i,qNorm)
            end do
!PRINT *,"pspPart",pspPart
            do i = 1, 3
               do j = i, 3
                  IE_stress(i,j) = IE_stress(i,j) &
                     + rhoRecip(ix, iy, iz) * qPoint(i) * qPoint(j) * pspPart / qNorm
               end do
            end do
!PRINT *,"1",IE_stress
!pause
         end if
      end do
      end do
!!$ print *, omp_get_thread_num(),iz,IE_stress(1,1)
!PRINT *,"IE iz", iz
!PRINT *,IE_stress
      end do
!$OMP END DO
!$OMP CRITICAL

      CALL MPI_ALLREDUCE(IE_stress,IonEle_stress,9,mpi_REAL8,mpi_sum,parallel%commx,mpinfo)

!$OMP END CRITICAL
!$OMP END PARALLEL
!print *,"I0"
!print *,IonEle_stress
      IonEle_stress = -2.d0 * IonEle_stress
      do i = 1, 3
         IonEle_stress(i,i) = IonEle_stress(i,i) - energy
      end do
!print *,"I2"
!print *,IonEle_stress
      IonEle_stress = IonEle_stress /volume
!print *,"I3"
!print *,IonEle_stress
!stop
   end function IonEle_stress!}}}
!----------------------------------------------------------
   FUNCTION PseudoPotDiffLookup(Ity,qNorm)  !{{{
      use constants
      use parameters, only: PP_identifer
      USE pspot_module , ONLY : tknots,psp
      USE mathsplines , ONLY : spline_cubic_val
      implicit none
!type(pseudo_potentialfile),intent(in)  :: psp

      real(kind=DP), intent(in) :: &
         qNorm            ! The norm of the q-vector (or r) to be looked up
      INTEGER(I4B)      ,intent(in)   :: Ity
      real(kind=DP) ::  PseudoPotDiffLookup


      REAL(KIND=DP) :: & ! mohan add (KIND=DP) on 02-14-2013
        pspSpacing, &    ! The spacing between values of the pseudopotential.
        rawIndex, &      ! The projected index (may be fractional) at which we
        yval, &
        ypval, &
        yppval
!>> INITIALIZATION <<!

! Calculate pspSpacing and rawIndex.
      pspSpacing =psp(Ity)%qspacing !psp%maxq / REAL(psp%numPoints(ntype) - 1, kind=DP)
      rawIndex = qNorm / pspSpacing  + 1.d0
!PRINT *,"rawIndex",rawIndex
! First take care of the possiblity that our q is out of range.
      if ( qNorm >= psp(Ity)%qmax ) then
         PseudoPotDiffLookup = 0.0_DP
! Otherwise, do a pseudo-cubic-spline interpolation.
      else
         if(PP_identifer==0)then
         call spline_cubic_val( psp(Ity)%numps, tknots &
              &, psp(Ity)%VlocqS,psp(Ity)%ddVl_dq2 &
              &, rawIndex, yval, ypval, yppval )
         else
         call spline_cubic_val( psp(Ity)%numps, tknots &
              &, psp(Ity)%VlocqS,psp(Ity)%ddVl_dq2 &
              &, rawIndex, yval, ypval, yppval )
         endif
!PRINT *,"ypval",ypval
         PseudoPotDiffLookup = ypval/pspSpacing + &
              & 8.0_DP*pi*psp(Ity)%Zion/qNorm**3
!PRINT *,"pse",PseudoPotDiffLookup
      END IF
   end function !}}}
!---------------------cal_force_stress---------------------
   SUBROUTINE cal_force_stress()

      USE smpi_math_module,  ONLY: parallel,start_time,end_time,write_time

      USE parameters , ONLY : LBvK,Idiag,NPRR,Lpbc,CalForce&
           & ,hartree_method
      USE begin_module , ONLY : Initial_grid,Initial_density,init_ISO_density
      USE grid_module , ONLY : grid,eigen,eigen_r,rho_calc,n
      USE potential_module , ONLY : vlpp
      USE nlpot_module , ONLY : initialize_nlpot

      USE scf_module , ONLY :  ISO_CheFSI, CheFSI

      USE struct_module , ONLY : struct,natom
      USE out_module , ONLY : KSout_list
      USE GetVLocalPseudoPotential , ONLY:CalVlpp
      USE poisson_isf
      USE parameters , ONLY : Nspin,Nstates

      USE grid_module , ONLY : gap,n1,n2,n3,global_n1,global_n2,global_n3
      USE array_io    , ONLY : output
      USE succeed     , ONLY : store_rho_at,store_rho,store_r,store_psi,Llastrho

      USE m_time_evaluate, ONLY: memory_sum
      IMPLICIT NONE
!LOCAL

      INTEGER(I4B),save   :: counter=0
      INTEGER(I4B)   :: i
      LOGICAL :: l1=.false.

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      counter=counter+1
! print*,counter
!initialize grid mesh
      CALL start_time('init_grid',l1)

!initialize grid mesh
      CALL Initial_grid()

      CALL end_time('init_grid',l1)
      CALL write_time('init_grid',l1)

!nonlocal part pseudo-potential

! CALL MPI_BARRIER(parallel%comm,mpinfo)
      CALL start_time('init_nlpot',l1)

      CALL initialize_nlpot()

      CALL end_time('init_nlpot',l1)
      CALL write_time('init_nlpot',l1)

!local part pseudo-potential

      CALL start_time('cal-vlpp',l1)

      IF(Lpbc)THEN
        CALL vlpp()
      ELSE
        CALL CalVlpp()
      ENDIF

      CALL start_time('cal-vlpp',l1)
      CALL end_time('cal-vlpp',l1)
      CALL write_time('cal-vlpp',l1)

!initialize charge density
      IF(Lpbc)THEN
        CALL initial_density()
      ELSE

! if(parallel%isroot)print *,"POSCAR",struct%poscar
         CALL start_time('init_density',l1)
         CALL init_ISO_density()
         CALL store_rho_at(parallel%mygrid_range(3),Nspin,rho_calc%OneDSphere)
         CALL end_time('init_density',l1)
         CALL write_time('init_density',l1)

      ENDIF
!> =======================================
!> just temporary                       !=
      IF((.not.Lpbc).and.hartree_method==6)CALL karray_set()  !=
!> =======================================
!====================
!##cal force directly
!goto 10011
!====================
!Self-consistent

!IF (parallel%isroot)THEN
!   print*,"!===========temporary============!"
!   print*,"sphere size",size(rho_calc%OneDwvf,1)
!   print*,"gap",gap
!   print*,"!===========temporary============!"
!ENDIF
      if(Lpbc)then
         CALL CheFSI(grid%rhoS,eigen%wvf,eigen%val)
      else
         CALL ISO_CheFSI(rho_calc%OneDSphere,rho_calc%OneDwvf,rho_calc%OneDval)
      endif






!output some data
      CALL KSout_list()
      if(.not.Lpbc)then

!> store the rho and atom position
         CALL store_r(natom,struct%poscar)
         CALL store_rho(parallel%mygrid_range(3),Nspin,rho_calc%OneDSphere)
         CALL store_psi(parallel%mygrid_range(3),Nstates,Nspin,rho_calc%OneDwvf)

      endif

!> cal force and stress
10011 IF(.NOT.LBvK.AND.Lpbc)THEN  !##periodic
!calculate the forces
         CALL cal_force_c(grid%rhoS,eigen%wvf,struct%forces)
!calculate the stress
! eigen%wvf=(0.001,0.001)
         CALL cal_stress(grid%rhoS,eigen%wvf,struct%stress)
      ELSEIF(.NOT.LBvK.AND.(.NOT.Lpbc).AND.CalForce)THEN  !##confined
!ALLOCATE(rho_grid(n,size(rho_calc%OneDwvf,2),size(rho_calc%OneDwvf,3),size(rho_calc%OneDwvf,4)))

         CALL start_time('cal_force',l1)
         CALL cal_force_r(grid%rhoS,rho_calc%OneDwvf,struct%forces)
         CALL end_time('cal_force',l1)
         CALL write_time('cal_force',l1)

         struct%stress=0.d0
      ENDIF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_force_stress
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE forcestress_module
