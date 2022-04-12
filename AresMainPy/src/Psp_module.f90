!##########################################!
!*For    : store read pesudopotential data !
!*Author : Qiang Xu                        !
!*Date   : 2017-7-12                       !
!##########################################!
MODULE pspot_module
   USE constants
   IMPLICIT NONE
   !defined
   TYPE pspot
      REAL(DP) :: Zion
      INTEGER(I4B) :: numps
      INTEGER(I4B) :: qnumps
      !> for broadcast denr
      INTEGER(I4B) :: numps_den
      !> local part
      REAL(DP)     :: qmax
      REAL(DP)     :: qspacing
      REAL(DP),ALLOCATABLE :: qmesh(:)
      REAL(DP),ALLOCATABLE :: Vlocq(:)
      REAL(DP),ALLOCATABLE :: VlocqS(:)
      REAL(DP),ALLOCATABLE :: ddVl_dq2(:)
      !> nonlocal real space
      INTEGER(I4B) :: nproj
      INTEGER(I4B),ALLOCATABLE :: proj_l(:)
      INTEGER(I4B),ALLOCATABLE :: proj_m(:)
      REAL(DP) :: rcut
      REAL(DP) :: rmax
      REAL(DP) :: rspacing
      REAL(DP),ALLOCATABLE :: D0(:,:)
      REAL(DP),ALLOCATABLE :: beta_r(:,:)
      REAL(DP),ALLOCATABLE :: dbeta_dr(:,:)
      !> radial charge density
      REAL(DP),ALLOCATABLE :: denr(:)
      REAL(DP),ALLOCATABLE :: ddden_dr2(:)
      !> for real pseudopotential used in confined system
      REAL(DP),ALLOCATABLE :: r_real(:)
      REAL(DP),ALLOCATABLE :: V_loc(:)
   ENDTYPE pspot
   !
   TYPE(pspot),ALLOCATABLE :: psp(:)
   !
   INTEGER(I4B) :: max_nproj
   REAL(DP)     :: max_rcut
   REAL(DP),ALLOCATABLE :: tknots(:)
ENDMODULE pspot_module
