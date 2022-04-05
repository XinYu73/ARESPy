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

!   REAL(DP),ALLOCATABLE  :: V_hxc_old(:,:,:),V_hxc_new(:,:,:)
!   LOGICAL               :: ACCELERATE=.FALSE.
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-----------------------------------------------------------------
   SUBROUTINE Calveff(nps,rhoS,rho,veffS)
      USE parameters   , ONLY : NSPIN
      USE xc_module, ONLY : XC_functional

      USE m_time_evaluate, ONLY: memory_sum,memory_free,filename
      USE smpi_math_module
      USE Grid_module  , ONLY : grid,n1,n2,n3,global_n1,global_n2,global_n3

      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN) :: rhoS(nps,nspin),rho(nps)
      REAL(DP),INTENT(OUT) :: veffS(nps,nspin)
!LOCAL
      INTEGER(I4B) :: Is
!REAL(DP) :: tmp,tmpl
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!tmpl=SUM(ABS(grid%vlpp))
!CALL MPI_ALLREDUCE(tmpl,tmp,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
!print*,'vlpp',tmp

!tmpl=SUM(rho**2)
!CALL MPI_ALLREDUCE(tmpl,tmp,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
!print*,'rho2',tmp
!hartree term
      CALL vhartree(nps,rho,grid%vh)
!tmpl=SUM(ABS(grid%vh))
!CALL MPI_ALLREDUCE(tmpl,tmp,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
!print*,'vh',tmp
!CALL vlda(rho,grid%vxc)
      CALL XC_functional(nps,rhoS,grid%vxcS)
!tmpl=SUM(ABS(grid%vxcS))
!CALL MPI_ALLREDUCE(tmpl,tmp,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
!print*,'vxc',SUM(grid%vcvS*rhoS)*dvol
!veff
      DO Is=1,NSPIN
         veffS(:,Is)=grid%vlpp(:)+grid%vh(:) +grid%vxcS(:,Is)
      ENDDO
!veffS(:,:)=0._DP
!STOP
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Calveff
!-----------------------------------------------------------------
   SUBROUTINE vhartree(nps,rho,vhart)

     USE grid_module , ONLY : ng1,ng2,ng3,grid,global_n1&
          & ,global_n2,global_n3,global_n,sphere2cubic_fft&
          & ,cubic2sphere_fft,fft_sph
     USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
     USE m_time_evaluate, ONLY :filename

      USE FOURIER
      IMPLICIT NONE
!
      INTEGER(I4B),INTENT(In) :: nps
      REAL(DP),INTENT(IN) :: rho(nps)
      REAL(DP),INTENT(OUT) :: vhart(nps)
!

      COMPLEX(DCP),DIMENSION(ng1,ng2,parallel%local_z) :: rho_recip&
           &,vh_recip

      INTEGER(I4B) :: Ig,I1,I2,I3

      REAL(DP)  :: rho_FFT(global_n1,global_n2,parallel%local_z)
      REAL(DP)  :: vhart_FFT(global_n1,global_n2,parallel%local_z)
      REAL(DP),allocatable  :: array_fft1d(:)
      INTEGER(I4B) :: fft1d_size

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! fft1d_size=parallel%fft_grid_range(2,parallel%rankx+1)&
!      & -parallel%fft_grid_range(1,parallel%rankx+1)+1
! allocate(array_fft1d(fft1d_size))
! CALL MPI_ALLTOALLV(rho,parallel%fft_scount&
!      & ,parallel%fft_sdispls,MPI_REAL8&
!      & ,array_fft1d,parallel%fft_rcount&
!      & ,parallel%fft_rdispls,MPI_REAL8&
!      & ,parallel%commx,mpinfo)
! !> sphere to cubic
! rho_FFT=0.d0
! CALL sphere2cubic_fft(fft1d_size,array_fft1d,rho_FFT&
!      & ,parallel%fft_grid_range(1,parallel%rankx+1)-1&
!      & ,parallel%local_z_start)

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

      if(parallel%rankx==0) vh_recip(1,1,1)=cmplx(0.d0,0.d0)

!

      vhart=fft_sph(vh_recip)
! vhart_FFT=FFT(vh_recip)
! CALL cubic2sphere_fft(fft1d_size,vhart_FFT,array_fft1d&
!      & ,parallel%fft_grid_range(1,parallel%rankx+1)-1 &
!      & ,parallel%local_z_start)
! CALL MPI_ALLTOALLV(array_fft1d,parallel%fft_rcount&
!      & ,parallel%fft_rdispls,MPI_REAL8&
!      & ,vhart,parallel%fft_scount&
!      & ,parallel%fft_sdispls,MPI_REAL8&
!      & ,parallel%commx,mpinfo)
! deallocate(array_fft1d)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE vhartree
!-------------------------PARTING LINE----------------------------
   SUBROUTINE vlpp()
      USE math , ONLY : interp, CubicSplineInterp

      USE grid_module , ONLY : grid,ng1,ng2,ng3,n,gap&
           & , global_n1,global_n2,global_n3 ,dvol,cubic2sphere_fft,fft_sph
      USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo

      USE ewald, ONLY : ewald_energy
      USE struct_module , ONLY : volume,struct,naty,lat_mat,Eionion
      USE pspot_module , ONLY : psp
      USE FOURIER
      IMPLICIT NONE
!IN/OUT
!LOCAL
      INTEGER(I4B) :: Ix,Iy,Iz,Ity,Ia,Ig

      COMPLEX(DCP),DIMENSION(ng1,ng2,parallel%local_z) :: vlpp_recip

      REAL(DP) :: qPoint(3),qNorm,Vgatom
      COMPLEX(DCP) :: tepc
      REAL(DP) :: gmax

      REAL(DP) :: vlpp_global(global_n1,global_n2,parallel%local_z)
      REAL(DP),allocatable :: vlpp_fft(:)

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
            vgatom=CubicSplineInterp(psp(Ity)%VlocqS&
                 &,psp(Ity)%ddVl_dq2 , &
                 & psp(Ity)%qmax,psp(Ity)%qspacing&
                 &,qnorm,psp(Ity)%Zion)
            tepc=(0.d0,0.d0)
            DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
               tepc= tepc+ EXP(IMAG*CMPLX(DOT_PRODUCT(qPoint,struct%poscar(:,Ia)),0.0_DP))
            ENDDO
            vlpp_recip(ix,iy,iz)=vlpp_recip(ix,iy,iz) + vgatom*tepc
         ENDDO
      ENDDO
      ENDDO
      ENDDO

!
      vlpp_recip(:,:,:)=vlpp_recip/volume


! vlpp_global=FFT(vlpp_recip)
! allocate(vlpp_fft(parallel%fft_grid_range(2,parallel%rankx+1)&
!      & -parallel%fft_grid_range(1,parallel%rankx+1)+1))
! CALL cubic2sphere_fft(size(vlpp_fft),vlpp_global,vlpp_fft&
!      & ,parallel%fft_grid_range(1,parallel%rankx+1)-1&
!      & ,parallel%local_z_start)
! CALL MPI_ALLTOALLV(vlpp_fft,parallel%fft_rcount&
!      & ,parallel%fft_rdispls,MPI_REAL8&
!      & ,grid%vlpp,parallel%fft_scount&
!      & ,parallel%fft_sdispls,MPI_REAL8&
!      & ,parallel%commx,mpinfo)
      grid%vlpp=fft_sph(vlpp_recip)

!STOP
!calculation the Eewald
      Eionion=ewald_energy(lat_mat,struct%pos,struct%eleid,struct%Zion)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE vlpp
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE potential_module
