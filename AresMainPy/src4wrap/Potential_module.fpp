# 1 "Potential_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Potential_module.f90"
MODULE potential_module
   !########################################!
   !*For    : Potential calculation Module  !
   !*Author : Qiang Xu                      !
   !*Date   : 2017/07/05                    !
   !########################################!
   USE constants
   IMPLICIT NONE



   REAL(DP),ALLOCATABLE  :: V_accelerate(:)

   REAL(DP),ALLOCATABLE  :: V_hxc_old(:,:,:),V_hxc_new(:,:,:)
   LOGICAL               :: ACCELERATE=.FALSE.
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !-----------------------------------------------------------------
   SUBROUTINE cal_veff(rhoS,veffS)
      USE parameters   , ONLY : NSPIN,Lpbc
      USE libxc_module , ONLY : LDAlib_potential
      USE IsolateSet   , ONLY : IsolateVCoulomb_a,rho_calc

      USE m_time_evaluate, ONLY: memory_sum,memory_free,filename
      USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
      USE Grid_module  , ONLY : grid,n1,n2,n3,global_n1,global_n2,global_n3




      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: veffS(:,:,:,:)
      !LOCAL
      REAL(DP),DIMENSION(n1,n2,n3)       :: vh & !hartree potential
                                    &  ,    rho  !total density
      REAL(DP),DIMENSION(n1,n2,n3,NSPIN) :: vxcS    ! xc potential
      INTEGER(I4B) :: Is
      REAL(DP),DIMENSION(n1,n2,n3) :: out_grid
      ! REAL(DP),DIMENSION(rho_calc%OneDLength) :: out_sphere

      INTEGER(I4B) :: ix,iy
      REAL(DP),allocatable :: veff_global(:,:,:)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('cal_veff_local',(real(size(vh),DP)+size(rho)+size(vxcS)&
           &+size(out_grid))*DP)
      !total density
      CALL sumrhoS(rhoS,rho)
      !hartree term
      CALL vhartree(rho,vh)
      !CALL vlda(rho,grid%vxc)
      CALL LDAlib_potential(rhoS,vxcS)
      !veff
      DO Is=1,NSPIN
         veffS(:,:,:,Is)=grid%vlpp(:,:,:)+vh(:,:,:)+vxcS(:,:,:,Is)
      ENDDO
      call memory_free('cal_veff_local',(real(size(vh),DP)+size(rho)+size(vxcS)&
           &+size(out_grid))*DP)
! #ifdef 1
!      ! offset=mod(parallel%mygrid_range(1),local_n1*local_n2)
!      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!      allocate(veff_global(global_n1,global_n2,global_n3))
!      call MPI_ALLGATHERV(vh(ix,iy,1),parallel%mygrid_range(3),&
!           & MPI_REAL8,veff_global,parallel%recvcounts,parallel%displs&
!           &, MPI_REAL8, parallel%commx,mpinfo)
!      open(1111+parallel%myid,file="init_density"//filename(10:11))
!      write(1111+parallel%myid,*)shape(veff_global)
!      write(1111+parallel%myid,*)veff_global
!      close(1111+parallel%myid)
!      call smpi_exit()
! #endif
!      open(1111,file="init_density")
!      write(1111,*)shape(rhoS(:,:,:,1))
!      write(1111,*)rhoS(:,:,:,1)
!      close(1111)
!      stop
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_veff
   !-----------------------------------------------------------------
   SUBROUTINE cal_veff_iso(rhoS,veffS)
      USE parameters   , ONLY : NSPIN,Lpbc
      USE libxc_module , ONLY : LDAlib_potential_iso
      USE Grid_module  , ONLY : grid,n1,n2,n3
      USE IsolateSet   , ONLY : IsolateVCoulomb_a,rho_calc

      USE smpi_math_module     , ONLY : parallel
      USE IsolateSet   , ONLY : rho_calc
      USE Grid_module  , ONLY : parallel_s2g

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !IN/OUT




      REAL(DP),INTENT(IN) :: rhoS(parallel%mygrid_range(3),NSPIN)
      REAL(DP),INTENT(OUT) :: veffS(parallel%mygrid_range(3),NSPIN)

      !LOCAL





      REAL(DP),DIMENSION(parallel%mygrid_range(3))       :: vh & !hartree potential
                                    &  ,    rho  !total density
      REAL(DP),DIMENSION(parallel%mygrid_range(3),NSPIN) :: vxcS    ! xc potential

      INTEGER(I4B) :: Is

      INTEGER(I4B) :: i,j,x,y,z

      REAL(DP),DIMENSION(n1,n2,n3) :: out_grid
      REAL(DP),DIMENSION(rho_calc%OneDLength) :: out_sphere
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('cal_veff_local',(real(size(vh),DP)+size(rho)+size(vxcS)&
           &+size(out_grid)+size(out_sphere))*DP)
      !total density
      CALL sumrhoS_iso(rhoS,rho)
      !hartree term
      IF(ACCELERATE)THEN
         vh=V_accelerate
         ! print *,'ACCELERATE ok'
      ELSE
         CALL IsolateVCoulomb_a(rho,vh)
         ! print *,"asdfqwer"
      ENDIF

      !CALL vlda(rho,grid%vxc)
      CALL LDAlib_potential_iso(rhoS,vxcS)
      !veff
      DO Is=1,NSPIN



         ! do i=1,parallel%mygrid_range(3),1
         !    x=rho_calc%x(i)
         !    y=rho_calc%y(i)
         !    z=rho_calc%z(i)
         !    veffS(i,Is)=rho_calc%vlpp(i)+vh(x,y,z)+vxcS(x,y,z,Is)
         ! enddo
         veffS(:,Is)=rho_calc%vlpp(:)+vh(:)+vxcS(:,Is)

      ENDDO
# 187 "Potential_module.f90"
      call memory_free('cal_veff_local',(real(size(vh),DP)+size(rho)+size(vxcS)&
           &+size(out_grid)+size(out_sphere))*DP)
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_veff_iso
   !-----------------------------------------------------------------
   SUBROUTINE vhartree(rho,vhart)

     USE grid_module , ONLY : ng1,ng2,ng3,n1,n2,n3,grid,global_n1&
          & ,global_n2,global_n3,global_n,fft_sph
     USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
     USE m_time_evaluate, ONLY :filename



      USE FOURIER
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN) :: rho(:,:,:)
      REAL(DP),INTENT(OUT) :: vhart(:,:,:)
      !



      COMPLEX(DP),DIMENSION(ng1,ng2,parallel%local_z) :: rho_recip&
           &,vh_recip

      INTEGER(I4B) :: Ig,I1,I2,I3

      REAL(DP)  :: rho_FFT(global_n1,global_n2,parallel%local_z)
      REAL(DP)  :: vhart_FFT(global_n1,global_n2,parallel%local_z)
      INTEGER(I4B) :: ix,iy
      REAL(DP),allocatable :: vhart_global(:,:,:)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! ! call MPI_ALLGATHERV(rho(ix,iy,1),parallel%mygrid_range(3)&
      ! !      &,MPI_REAL8,rho_global,parallel%recvcounts&
      ! !      &,parallel%displs,MPI_REAL8&
      ! !      & ,parallel%commx,mpinfo)
      ! CALL MPI_ALLTOALLV(rho(ix,iy,1),parallel%fft_scount&
      !      & ,parallel%fft_sdispls,MPI_REAL8&
      !      & ,rho_FFT,parallel%fft_rcount&
      !      & ,parallel%fft_rdispls,MPI_REAL8&
      !      & ,parallel%commx,mpinfo)
      ! rho_recip(:,:,:)=FFT(rho_FFT)
      rho_recip=fft_sph(rho)



      !q=0
      grid%gVec(4,1)=1.d0
      !v(q)=4*pi*rho(q)/(q^2)
      Ig=0

      DO I3=1,parallel%local_z
         Ig=(parallel%local_z_start+I3-1)*ng2*ng1



      DO I2=1,ng2
      DO I1=1,ng1
         Ig=Ig+1
         vh_recip(I1,I2,I3)=(4.d0*pi)*rho_recip(I1,I2,I3)/(grid%gVec(4,Ig))**2
      ENDDO
      ENDDO
      ENDDO
      !
      grid%gVec(4,1)=0.d0
      !remove q=0 because of netural system



      if(parallel%rankx==0)vh_recip(1,1,1)=cmplx(0.d0,0.d0)

      !



      ! vhart_FFT=FFT(vh_recip)
      ! vhart=0.d0
      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! CALL MPI_ALLTOALLV(vhart_FFT,parallel%fft_rcount&
      !      & ,parallel%fft_rdispls,MPI_REAL8&
      !      & ,vhart(ix,iy,1),parallel%fft_scount&
      !      & ,parallel%fft_sdispls,MPI_REAL8&
      !      & ,parallel%commx,mpinfo)
      vhart=fft_sph(vh_recip)
      ! call MPI_SCATTERV(vhart_global,parallel%recvcounts,parallel%displs,MPI_REAL8&
      !      & ,vhart(ix,iy,1),parallel%mygrid_range(3),MPI_REAL8,parallel%rootid&
      !      & ,parallel%commx,mpinfo)

! #ifdef 1
!      ! offset=mod(parallel%mygrid_range(1),local_n1*local_n2)
!      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!      ! allocate(vhart_global(global_n1,global_n2,global_n3))
!      ! call MPI_ALLGATHERV(vhart(ix,iy,1),parallel%mygrid_range(3),&
!      !      & MPI_REAL8,vhart_global,parallel%recvcounts,parallel%displs&
!      !      &, MPI_REAL8, parallel%commx,mpinfo)
!      open(1111+parallel%myid,file="init_density"//filename(10:11))
!      write(1111+parallel%myid,*)shape(rho_recip)
!      write(1111+parallel%myid,*)rho_recip
!      close(1111+parallel%myid)
!      call smpi_exit()
! #endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE vhartree
   !-------------------------PARTING LINE---------------------------
   SUBROUTINE vlpp()
      USE math , ONLY : interp, CubicSplineInterp

      USE grid_module , ONLY : grid,ng1,ng2,ng3,n1,n2,n3,n,gap&
           & , global_n1,global_n2,global_n3,fft_sph
      USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo



      USE struct_module , ONLY : volume,struct,naty
      USE pspot_module , ONLY : psp
      USE FOURIER
      USE parameters, ONLY : PP_identifer
      IMPLICIT NONE
      !IN/OUT
      !LOCAL
      INTEGER(I4B) :: Ix,Iy,Iz,Ity,Ia,Ig

      COMPLEX(DP),DIMENSION(ng1,ng2,parallel%local_z) :: vlpp_recip



      REAL(DP) :: qPoint(3),qNorm,Vgatom
      COMPLEX(DP) :: tepc
      REAL(DP) :: gmax

      REAL(DP) :: vlpp_global(global_n1,global_n2,parallel%local_z)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! gmax=pi*pio2/maxval(gap)**2
      vlpp_recip(:,:,:)=(0.d0,0.d0)
      Ig=0

      DO Iz=1,parallel%local_z
         Ig=(parallel%local_z_start+iz-1)*ng2*ng1



      DO Iy=1,ng2
      DO Ix=1,ng1
         Ig=Ig+1
         qPoint = -1.0_DP * grid%gVec(1:3,Ig)
         !PRINT *,"q",qPoint
         qNorm = grid%gVec(4,Ig)
         DO Ity=1,naty
            IF(PP_identifer==0)THEN
               vgatom=CubicSplineInterp(psp(Ity)%VlocqS&
                    &,psp(Ity)%ddVl_dq2 , &
                    & psp(Ity)%qmax,psp(Ity)%qspacing&
                    &,qnorm,psp(Ity)%Zion)
            ELSE
               ! vgatom=interp(psp(Ity)%qnumps,psp(Ity)%Vlocq&
                    ! &,psp(Ity)%qmesh,qNorm)
               vgatom=CubicSplineInterp(psp(Ity)%VlocqS&
                    &,psp(Ity)%ddVl_dq2 , &
                    & psp(Ity)%qmax,psp(Ity)%qspacing&
                    &,qnorm,psp(Ity)%Zion)
               ! if(vgatom>gmax)vgatom=0.d0
            ENDIF
            tepc=(0.d0,0.d0)
            DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
               tepc= tepc+ EXP(IMAG*CMPLX(DOT_PRODUCT(qPoint,struct%poscar(:,Ia)),0.0_DP))
            ENDDO
            vlpp_recip(ix,iy,iz)=vlpp_recip(ix,iy,iz) + vgatom* tepc
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !
      vlpp_recip(:,:,:)=vlpp_recip/volume


      ! grid%vlpp=0
      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! vlpp_global=FFT(vlpp_recip)
      ! CALL MPI_ALLTOALLV(vlpp_global,parallel%fft_rcount&
      !      & ,parallel%fft_rdispls,MPI_REAL8&
      !      & ,grid%vlpp(ix,iy,1),parallel%fft_scount&
      !      & ,parallel%fft_sdispls,MPI_REAL8&
      !      & ,parallel%commx,mpinfo)
      grid%vlpp=fft_sph(vlpp_recip)
      ! call MPI_SCATTERV(vlpp_global,parallel%recvcounts,parallel%displs,MPI_REAL8&
      !      & ,grid%vlpp(ix,iy,1),parallel%mygrid_range(3),MPI_REAL8,parallel%rootid&
      !      & ,parallel%commx,mpinfo)



   ENDSUBROUTINE vlpp
   !----------------------------vlpp_real--------------------------------
   subroutine vlpp_real()
     use grid_module , only : grid,n1,n2,n3,n
     use struct_module , only : struct,naty
     use pspot_module , only : psp
     use MathSplines   , only : polynom
     implicit none

     INTEGER(I4B) :: max_m,m
     INTEGER(I4B) :: igrid,ix,iy,iz,Ity,Iatom
     REAL(DP) :: rr,c(3),temp

      open(1212,file='vloc')
      write(1212,*)psp(1)%V_loc
      close(1212)
     !>
     grid%vlpp=0.d0
      do Ity= 1,naty,1
      max_m=size(psp(Ity)%r_real)-1
      do Iatom= struct%eleid(Ity),struct%eleid(Ity+1)-1
      igrid=0
      do Iz= 1,n3,1
      do Iy= 1,n2,1
      do Ix= 1,n1,1
         igrid=igrid+1
         rr=sqrt((grid%rvec(1,igrid)-struct%poscar(1,Iatom))**2 + &
            & (grid%rvec(2,igrid)-struct%poscar(2,Iatom))**2 + &
            & (grid%rvec(3,igrid)-struct%poscar(3,Iatom))**2 )
         !> interpolte
         m=1
         DO WHILE ( (psp(Ity)%r_real(m).lt.rr).and.(m.le.max_m) )
            !DO WHILE ( r_fv(m+1).lt.rr )
            m=m+1
         ENDDO
         if(m==1)m=2
         if(m.le.max_m)then
            temp=polynom(0,3,psp(Ity)%r_real(m-1:m+1),psp(Ity)%V_loc(m-1:m+1),c,rr)
         else
            temp=-psp(Ity)%Zion/rr
         endif
         grid%vlpp(Ix,Iy,Iz)=grid%vlpp(Ix,Iy,Iz)+temp
      enddo !> Ix
      enddo !> Iy
      enddo !> Iz
      enddo !> Iatom
      enddo !> Ity
      open(1212,file='vlpp_upf')
      write(1212,*)grid%vlpp
      close(1212)

   end subroutine vlpp_real
   !----------------------------vlda--------------------------------
   SUBROUTINE vlda(rho,LDAPotential)
      !
      !!!!! WARNING:
      !     This subroutine is from Carter's PROFESS. Later,we Need Change
      !
      !------------------------------------------------------------------------------
      ! DESCRIPTION:
      !   This function computes the exchange-correlation potential in the Local
      !   Density Approximation (LDA) based on the real-space electron density.
      !
      ! GLOBAL/MODULE VARIABLES CHANGED:
      !
      ! CONDITIONS AND ASSUMPTIONS:
      !   We can use LDAPointPot for moudularity.
      !
      ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
      !   Complete calculation to handle spin-polarized cases.
      !
      ! REFERENCES:
      !   [1]  Perdew J.P. and Zunger A., Phys. Rev. B 23(10), 1981, 5048-79
      !------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: &
         rho                    ! Electron density in real space, spin DEPENDENT
      REAL(kind=DP), DIMENSION(:,:,:),INTENT(OUT) :: &
         LDAPotential           ! The XC potential.
      REAL(kind=DP), PARAMETER :: &
         mot = -1._DP/3._DP, &                ! minus one-third.
         ot = 1.0_DP/3.0_DP, &                ! one third
         cX = -0.98474502184269641_DP, &      ! Exchange constant:  (4/3) * -0.75_DP * (3._DP/pi)**ot
         cC = 0.62035049089940009_DP, &       ! Correlation constant
         cC1 = -5.83666666666666666E-2_DP, &  ! (b(1) - a(1) * ot)
         cC2 = 1.33333333333333333E-3_DP, &   ! 2*ot * c(1)
         cC3 = -8.4E-3_DP, &                  ! ot * (2*d(1) - c(1))
         cC4 = -0.174798948333333333_DP, &    ! g(1) * (7._DP/6._DP) * b1(1)
         cC5 = -6.325709333333333333E-2_DP    ! g(1) * 4._DP * ot * b2(1)
      REAL(kind=DP) :: &
         !  eCU, &                          ! Correlation energy at rs for 1 spin
      !  eCP, &                          ! Corr. E. at rs if fully spin polarized
      rs, &                              ! (3/(4.pi.rho))^1/3
         sqrtRS
      !  zeta, &                         ! Spin-polarization
      REAL(kind=DP), DIMENSION(2), PARAMETER :: &
         a = (3.11E-2_DP,1.555E-2_DP), &      ! Parameter A for rs < 1 in [1]
         !b = (-4.8E-2_DP,-2.69E-2_DP), &      ! Parameter B for rs < 1 in [1]
      !c = (2E-3_DP,7E-4_DP), &             ! Parameter C (Ceperley-Alder) [1]
      !d = (-1.16E-2_DP,-4.8E-3_DP), &      ! Parameter D (C-A) from [1]
      g = (-1.423E-1_DP,-8.43E-2_DP), &    ! Parameter gamma for rs > 1 [1]
         b1 = (1.0529_DP,1.3981_DP), &        ! Parameter beta1 for rs > 1 [1]
         b2 = (3.334E-1_DP,2.611E-1_DP)       ! Parameter beta2 for rs > 1 [1]

      INTEGER(I4B) :: &
         ix, iy, iz                           ! Dummy counters
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      DO iz=1, SIZE(rho,3)
         LDAPotential(:,:,iz) = rho(:,:,iz)**ot
         DO iy=1, SIZE(rho,2)
            DO ix=1, SIZE(rho,1)
               IF(rho(ix,iy,iz) > 0._DP) THEN
                  ! Calculate r subscript s [1] in atomic units.
                  rs = cC / LDAPotential(ix,iy,iz)
                  IF (rs<1) THEN
                     LDAPotential(ix,iy,iz) = cX * LDAPotential(ix,iy,iz) + &
                        LOG(rs) * (a(1) + cC2 * rs) + cC1 + cC3*rs
                  ELSE
                     sqrtRS = SQRT(rs)
                     LDAPotential(ix,iy,iz) = cX * LDAPotential(ix,iy,iz) + &
                        (g(1) + cC4 * sqrtRS + cC5 * rs) / &
                        (1._DP + b1(1) * sqrtRS + b2(1) * rs)**2
                  END IF
               ELSE
                  LDAPotential(ix,iy,iz) = 0._DP
               END IF
            END DO
         END DO
      END DO
      !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE vlda
   !-----------------------total rho--------------------------------
   SUBROUTINE sumrhoS(rhoS,rho)
      USE parameters , ONLY : NSPIN
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: rho(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      rho(:,:,:)=0.d0
      DO Is=1,NSPIN
         rho(:,:,:)=rho(:,:,:) + rhoS(:,:,:,Is)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE sumrhoS
   !-------------------partline ---------------------
   SUBROUTINE sumrhoS_iso(rhoS,rho)
      USE parameters , ONLY : NSPIN
      IMPLICIT NONE
      !IN/OUT




      REAL(DP),INTENT(IN)  :: rhoS(:,:)
      REAL(DP),INTENT(OUT) :: rho(:)

      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






      rho=0.d0
      DO Is=1,NSPIN
         rho(:)=rho(:) + rhoS(:,Is)
      ENDDO

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE sumrhoS_iso
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE potential_module
