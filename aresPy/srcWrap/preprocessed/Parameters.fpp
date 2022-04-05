# 1 "Parameters.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Parameters.f90"
MODULE parameters
   USE constants
   IMPLICIT NONE
!---------------------------------------------------------
!basic parameters
   INTEGER(I4B)  :: ixc= 1   &    !type of XCDF
                 &, nspin=1       !switch for spin
   INTEGER(I4B)  :: finite_order=8  !the order of finite difference
   INTEGER(I4B)  :: ntype          !the number of elements' type
   INTEGER(I4B)  :: Naddstates=-1 , &  !number of eigenstate we need to add
                &   IGamma=0 , &
                &   ISTART=0 , &
                &   Isym=0   , &
                &   Idiag=1 ,  &
                &   CheM=-1 ,  &      !the order of chebyshev we use
                &   CheM0=-1,  &       !the free diag
                &   Nstates ,  &  !total states we calculate
                &   gridn(3)=-1 , & !grid mesh
                &   kgrid(3)=-1 , & !k-grid mesh
                &   Nssp=0           !number of simulation steps
!
   REAL(DP)      :: dcharge=0.0_DP
   REAL(DP)      :: init_gap=0.20_DP & !input gap
               &  , Ecut=-100.d0     & !Ecut
               &  , kspacing=0.2d0   & !kspacing
               &  , AtomRc=-1.d0   & !kspacing
               &  , Snlcc=0.5_DP
!
   LOGICAL       ::  & 
                &   LFIRST=.TRUE. , &
                &   LINRHO=.FALSE. ,&
                &   LAtomRho=.TRUE. , & !the density in psp file?
                &   LRROrthNorm=.FALSE. , & !chebyshev orthnormal basis set ?
                &   Lrandom=.FALSE. ,& !random chebyshev subspace
                &   Lcore_val=.FALSE.   !core-valance correlation
!
   CHARACTER(30) :: system_name='BEAST WAR' !system name,no use
   CHARACTER(30) :: cellfile_name='POSCAR'  !structre file name
   CHARACTER(30) :: ppfile_name(10)         !pseudopotential name
   CHARACTER(30) :: elements(10)            !input elements
!smearing
   INTEGER(I4B) :: Nsmear=-1         !The order for smearing
   REAL(DP)     :: Wsmear=0.01d0*ev2hart  !The width of broadening
!mixing parameter
   INTEGER(I4B)  :: IMIXER=2 &
                 &, NMITER=100 !max iteration
   INTEGER(I4B)  :: NSMIX=0 & !simple mixing number
               &,   NHMAX=5 &  !history for Anderson mixing used
               &,   NHMIN=2
   REAL(DP)      :: MALPHA=0.4d0 &  !simple mixing factor
               &,   MBETA=0.3d0  &  !Anderson mixing factor
               &, AMIX=0.4_DP    &
               &, BMIX=1*bohr2ang
   REAL(DP)      :: RESTA(3)=(/0.68,10.0,6.61/)!q0=0.68,epsilon0=10
! Rs=6.61 !>resta preconditioning
   REAL(DP)      :: W0AM=0.01 !for ill matrix
   REAL(DP)      :: RTOL=1e-4 & !for tolerance of density
               &,   ETOL=1e-3  !for tolerance of energy
!For output band
   LOGICAL       :: Lband=.TRUE.  !for bands
!>>>BCAST END<<<
!structural relax
   INTEGER(I4B)  :: IOPM=3  &
               &,   IGOAL=3 
   REAL(DP)      :: PRESS=0.d0  &
               &,   TOLF=1e-2   &
               &,   TOLP=0.5d0  &
               &,   TIMES=0.5d0 
!Outputing
   LOGICAL :: LWAVE=.FALSE. & !wave functions
           &, LCHARGE=.FALSE.
!Max Overlap Method
   LOGICAL :: LMOM=.FALSE. !Max overlap method?
   REAL(DP) :: MOMsigma=-1.0_DP
   INTEGER(I4B) :: nwf0=0
!For parallel data
   INTEGER(I4B)  :: BLOCK_MBNB=1  &      ! The NB and MB in 'pdsygvx'
               &,   Nstates_global
   INTEGER(I4B) :: NPARA(2)=-1
!---------------------------------------------------------
ENDMODULE parameters
