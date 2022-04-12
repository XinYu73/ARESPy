MODULE parameters
  USE constants
  IMPLICIT NONE
  !---------------------------------------------------------
  !basic parameters
  INTEGER(I4B)  :: iKEDF=1        !type of KEDF
  INTEGER(I4B)  :: iXCDF=1        !type of XCDF
  INTEGER(I4B)  :: finite_order=4 !the order of finite difference
  INTEGER(I4B)  :: ntype          !the number of elements' type
  INTEGER(I4B)  :: NADDSTATES=8 , &  !number of eigenstate we need to add
       &   ISTART=0 , &
       &   Idiag=1 ,  &
       &   CheM=12 ,  &      !the order of chebyshev we use
       &   CheM0=16,  &       !the free diag
       &   Nstates ,  &  !total states we calculate
       &   Nstates_global ,  &  !total states in parallel
       &   NPRR=-1 ,  &  !number of fraction occupation
       &   Nssp=0           !number of simulation steps
  !
  REAL(DP)      :: KSPACING=0.5   !input k-spacing
  REAL(DP)      :: KSHIFT(3)=0.d0   !input k-spacing
  REAL(DP)      :: init_gap=0.18 & !input gap
       &  , Ecut=-100.d0    !Ecut
  !
  LOGICAL       :: LKP=.FALSE.  !switch for kpoints file
  INTEGER(I4B)  :: PP_identifer=0
  !> identifer of Pseudopotential format:0:legacy; 1:UPF |xlt
  !For Isolate
  LOGICAL       :: Lpbc=.TRUE.        !logical of period boundary condition--xlt
  LOGICAL       :: CalForce=.True.        !logical of force calculated flag
  !>Identifier of poisson solver
  !> 0 :CG-GMG
  !> 1 :FMM
  !> 2:direct
  !> 3:MCM
  !> 4:cutoff-method
  INTEGER(I4B)  :: hartree_method=6
  REAL(DP)      :: Lcell=18, Lcellbohr             !Isolate cell length
  REAL(DP)      :: IsoRmax=9, IsoRmaxbohr            !Isolate sphere radius(ang)
  REAL(DP)      :: RadiusMax          !Isolate sphere radius(bohr)
  INTEGER(I4B)  ::NVC=6,ISOnorder=8   !grid order?,finite order of Isolate
  REAL(DP)      :: TOLCG=0.d-2        !the tolerence of poisson slover CG
  INTEGER(I4B)  :: NFCD=6             !compact finite difference order
  INTEGER(I4B)  :: ISOLmax=9          !the max l in multi-pole
  INTEGER(I4B)  :: iprec_fmm          !THE PRECSICE FOR FAST MULTIPOLE METHOD
  INTEGER(I4B)  :: ADDcharge=0        !Electric system : the number of charges
  INTEGER(I4B)  :: cell_shape=0       ! default value 0 (normal cubic cell) 1 (2d cell)
  REAL(DP)      :: cell_thick=8.d0    ! the thickness for 2d cell
  LOGICAL       :: Lpbc2iso           ! read the boundary charge from period boundary condition
  LOGICAL       :: Lradius_auto       ! Auto set the cutoff radius
  !For isolated parallel
  INTEGER(I4B)  :: BLOCK_MBNB=1       ! The NB and MB in 'pdsygvx'
  INTEGER(I4B)  :: debug_out=0       ! default value 0 (do nothing)
  !For Chebyshev
  REAL(DP)      :: Wexict=1.0d0
  !
  LOGICAL       :: LBvK=.FALSE. , &  !For BvK cell
       &   LFIRST=.FALSE. , &
       &   LINRHO=.FALSE. ,&
       &   LRadRho=.FALSE. , & !the density in psp file?
       &   LRROrthNorm=.FALSE. , & !chebyshev orthnormal basis set ?
       &   Lrandom=.FALSE.  !random chebyshev subspace
  !
  CHARACTER(30) :: system_name='ARES' !system name,no use
  CHARACTER(30) :: cellfile_name='POSCAR'  !structre file name
  CHARACTER(30) :: outfile='screen'  !structre file name
  CHARACTER(30) :: ppfile_name(10)         !pseudopotential name
  CHARACTER(30) :: elements(10)            !input elements
  !switch for spin
  INTEGER(I4B)  :: Nspin=1   !switch for spin
  !smearing
  INTEGER(I4B) :: Nsmear=0         !The order for smearing
  REAL(DP)     :: Wsmear=0.1d0     !The width of broadening
  !mixing parameter
  INTEGER(I4B)  :: IMIXER=1 &
       &, NMITER=100 !max iteration
  INTEGER(I4B)  :: NSMIX=10 & !simple mixing number
       &, NHMIX=5 &  !history for Anderson mixing used
       &, NHMIN=2
  REAL(DP)      :: MALPHA=0.3d0 & !simple mixing factor
       &,   MBETA=0.8d0  &  !Anderson mixing factor
       &,   AMIX=0.4d0   &  !for kerker
       &,   BMIX=bohr2ang       !for kerker
  REAL(DP)      :: RESTA(3)=(/0.68,10.0,6.61/)!q0=0.68,epsilon0=10
                       ! Rs=6.61 !>resta preconditioning
  REAL(DP)      :: W0AM=0.01 !for ill matrix
  REAL(DP)      :: RTOL=1e-4 & !for tolerance of density
       &,   ETOL=1e-3  !for tolerance of energy
  !for subsystem dft
  LOGICAL       :: LSUB=.FALSE.   !For Subsystem DFT
  LOGICAL       :: LFAT=.FALSE.   !For FAT-cycle
  LOGICAL       :: LONE=.FALSE.   !For one-susystem
  INTEGER(I4B)  :: Nsub=1         !number of subsystems
  INTEGER(I4B)  :: NFAT=0         !number of FAT
  REAL(DP)      :: TFVW(2)        !aTFbvW coe
  !For OFDFT
  LOGICAL       :: LOFDFT=.FALSE. !
  !For output band
  LOGICAL       :: Lband=.FALSE.  !for bands
  !structural relax
  INTEGER(I4B)  :: IOPM=3  &
       &,   IGOAL=3
  REAL(DP)      :: PRESS=0.d0  &
       &,   TOLF=1e-2   &
       &,   TOLP=0.5d0  &
       &,   TIMES=0.5d0
  INTEGER(I4B)  :: ISP
  !For double grid
  LOGICAL      :: LDG=.false.
  INTEGER(I4B) :: NDG=0 !number grid points betweed coarse grid
  INTEGER(I4B) :: n_near=0 !number grid points near coarse grid to cal
  INTEGER(I4B) :: inpol=1 !the order of interpolation for dense grid
  !For initMO
  INTEGER(I4B) :: IDinit=0
  CHARACTER(30):: MO_file
  !> optimize
  REAL(DP) :: POTIM = 0.1
  REAL(DP) :: EDIFFG = 3E-3
  INTEGER(I4B) :: step_fixrho = 0
  !> IO INFO
  TYPE IO_INFO_TYPE
     integer(i4b)  ::  IU0 = 509000 ! OUTPUT FILE OF ARES, MAYBE ANYTHING
     integer(i4b)  ::  IU1 = 509001 ! OUTPUT FILE OF 
     integer(i4b)  ::  IU2 = 509002 ! OUTPUT FILE OF ANY MD SIMULATION INFORMATION
     integer(i4b)  ::  IU3 = 509003 ! OUTPUT FILE OF 
     integer(i4b)  ::  IU4 = 509004 ! OUTPUT FILE OF 
     integer(i4b)  ::  IU5 = 509005 ! OUTPUT FILE OF 
     integer(i4b)  ::  IU6 = 509006 ! OUTPUT FILE OF GEOMETRY OPTIMIZATION (CONTCAR)
     integer(i4b)  ::  IU8 = 509008 ! OUTPUT FILE OF MD SIMULATION (XDATCAR)
  END TYPE IO_INFO_TYPE

  !> fix for stress
  INTEGER(I4B) :: fix_xyz(6)
  ! TYPE atom_fix
  !    INTEGER(I4B) :: x
  !    INTEGER(I4B) :: y
  !    INTEGER(I4B) :: z
  !    INTEGER(I4B) :: xy
  !    INTEGER(I4B) :: xz
  !    INTEGER(I4B) :: yz
  ! ENDTYPE atom_fix
  ! type(atom_fix) :: xyz_fix
  !> for CG
  REAL(DP) :: pstress
  real(dp) :: ramax
  !> for periodic parallel code
  INTEGER(I4B)  :: Grad_order=16 
  INTEGER(I4B)  :: Maxnpts = 0
  INTEGER(I4B)  :: nevshift = 0
  INTEGER(I4B)  :: gridn(3)= -1  ! grid mesh
  LOGICAL       :: LGamma=.false. 
  LOGICAL       :: Lcore_val = .FALSE.
    INTEGER(I4B)  :: IOUNITS(20)
    INTEGER(I4B)  :: Isym=0
    LOGICAL      :: Lfinite_full=.true.
    LOGICAL      :: lforce = .false.
    LOGICAL      :: lstress = .false.
    REAL(DP)     :: temper_elec = -1.d0
    CHARACTER(LEN=20)   :: lwstyle = "vasp"
    LOGICAL       :: LOPT=.FALSE.  
    INTEGER(I4B)  :: IBRON=3
    INTEGER(I4B)    :: maxsave = 20
    INTEGER(I4B)    :: IGamma=0
    INTEGER(I4B)    :: kgrid(3) = -1
    !--------------- Md Parameters --------------------------
    LOGICAL       :: LMD=.FALSE.  
    ! input
    INTEGER(I4B)           :: nhis = 200
    REAL(DP)               :: rdfr = 20.d0 
    CHARACTER(LEN=30)      :: thermostat = "nose-hoover"
    CHARACTER(LEN=30)      :: integrator = "velocity-verlet"
    ! INTEGER(I4B)           :: nssp = 0 
    INTEGER(I4B)           :: sfreq = 5    
    REAL(DP)               :: cfreq = 0.1d0
    REAL(DP)               :: temperature =297.d0
    REAL(DP)               :: mste = 0.d0
    REAL(DP)               :: delt = 0.2 
    REAL(DP)               :: relaxt = 1.d0
    CHARACTER(10)          :: ensemble = "nvt"
    INTEGER(I4B)           :: dof    
    REAL(DP)               :: pext = 0.d0
    INTEGER(I4B)           :: iresmd = 0
    ! ----------------------nose-hoover chain ----------------------
    !
    INTEGER(I4B)           :: nresn = 3
    INTEGER(I4B)           :: nyosh = 3
    INTEGER(I4B)           :: nnhc = 2
    REAL(DP),ALLOCATABLE   :: wdti2(:), wdti4(:), wdti8(:)
    REAL(DP),ALLOCATABLE   :: qmass(:)
    REAL(DP),ALLOCATABLE   :: syin_coeff(:)
    REAL(DP)               :: bmass = -1.d0
    REAL(DP)               :: fthermo = -1.d0
    REAL(DP)               :: fbaro = -1.d0
    REAL(DP)               :: fthermown = -1.d0
    REAL(DP)               :: fbarown = -1.d0
    !----------------- write ----------------------------------
    LOGICAL                  :: lbin = .false.
    !----------------- new md ----------------------------------
    INTEGER(I4B)                     :: p_flag(3)
    REAL(DP), dimension(6) :: p_start,p_stop, p_freq,p_target
    REAL(DP)                         :: t_start,t_stop,t_freq,t_target
    CHARACTER(LEN=CLen)      :: press_control='xyz'
    CHARACTER(LEN=CLen)      :: pcontrol='mtk_xyz'
    REAL(DP)                         :: pmass=-1.d0
    
    REAL(DP)                         :: pdrag = 0.0 
    REAL(DP)                         :: tdrag = 0.0 
    REAL(DP)                         :: erate = 1D10
    REAL(DP)                         :: mstrain=-1.d0
    CHARACTER(LEN=CLen)              :: sdir='z'
    LOGICAL                :: Lke_pot=.false.
  
  !---------------------------------------------------------
CONTAINS
  SUBROUTINE generate_infile(infile_name)
    implicit none
    character(len=*) :: infile_name
    !>
    open(1173,file=infile_name)
    write(1173,'(A)')"##########################BASIS DATA##############################"
    write(1173,'(A)')"#This file need more custom configure"
    write(1173,'(A)')"SYSTEM = OPTIMUS_PRIME"

    write(1173,'(A)')"#Born-von Karman super-cell or k-space representation"
    write(1173,'(A)')"LBvK = F"

    write(1173,'(A)')"#Cell file name"
    write(1173,'(A)')"CELLFILE = ARES.CONT"

    write(1173,'(A)')"#output file name"
    write(1173,'(A)')"OUTFILE = screen"

    write(1173,'(A)')"#Ecut (eV)"
    write(1173,'(A)')"Ecut = -4680"

    write(1173,'(A)')"#real-space gridsize (angs)"
    write(1173,'(A)')"GRIDSIZE = 0.20"

    write(1173,'(A)')"#k-points"
    write(1173,'(A)')"KSPACING = -0.35"

    write(1173,'(A)')"#Use k-points file"
    write(1173,'(A)')"LKPT=F"

    write(1173,'(A)')"#elements type"
    write(1173,'(A)')"ELEMENTS=C # B N"

    write(1173,'(A)')"#pseudo-potential"
    write(1173,'(A)')"#PPFILE = C_6000.rrpot"
    write(1173,'(A)')"#PPFILE = Al_lda.rrpot"
    write(1173,'(A)')"PPFILE = C_iso.realpot #B_iso.realpot N_iso.realpot"
    write(1173,'(A)')"#PPFILE = out.realpot"
    write(1173,'(A)')"#PPFILE = Al02.rrpot"
    write(1173,'(A)')"#exchange correlation type"
    write(1173,'(A)')"XCDF = 1"
    write(1173,'(A)')"#Finite difference order"
    write(1173,'(A)')"NORDER = 8"
    write(1173,'(A)')"#Switch for spin"
    write(1173,'(A)')"LSPIN= F"

    write(1173,'(A)')"#The order of smearing (2-5) , <=0 is fermi-dirac-like smearing"
    write(1173,'(A)')"Nsmear= -1"

    write(1173,'(A)')"#Width of broadening(eV.)"
    write(1173,'(A)')"Wsmear= 0.01"

    write(1173,'(A)')"#Usu double-grid technique"
    write(1173,'(A)')"LDG= T"
    write(1173,'(A)')"#Order of Double-grid"
    write(1173,'(A)')"NDG = 1520"
    write(1173,'(A)')"#dense grid calculate per side"
    write(1173,'(A)')"N_NEAR = 1"
    write(1173,'(A)')"#the order of interpolation for dense grid"
    write(1173,'(A)')"INPOL = 1"

    write(1173,'(A)')"#Number of simulate steps"
    write(1173,'(A)')"Nssp=1 #000"
    write(1173,'(A)')"#Optimiz Method"
    write(1173,'(A)')"IOPM = 4"
    write(1173,'(A)')"#Steps of SCF with fixed charge rho"
    write(1173,'(A)')"step_fixrho = 0"
    write(1173,'(A)')"#tolerance of optimize force (ipom=4)"
    write(1173,'(A)')"ediffg = 3E-3"
    write(1173,'(A)')"#ISTART"
    write(1173,'(A)')"ISTART= 0"
    write(1173,'(A)')"##############################Diag################################"
    write(1173,'(A)')"#0:ARPACK 1:Chebyshev"
    write(1173,'(A)')"Idiag=1"
    write(1173,'(A)')"#Chebyshev order [8,20]"
    write(1173,'(A)')"CheM=16"
    write(1173,'(A)')"CheM0=34"
    write(1173,'(A)')"#For add eigenstates to diag(H)"
    write(1173,'(A)')"NADST=8"
    write(1173,'(A)')"#first step use diag(H)"
    write(1173,'(A)')"Lfirst=F"
    write(1173,'(A)')"#partial RayleighRitz ()"
    write(1173,'(A)')"NPRR= -10"
    write(1173,'(A)')"#RayleighRitz need OrthNorm?"
    write(1173,'(A)')"LOrthNorm=T"
    write(1173,'(A)')"#initialize by random subspace"
    write(1173,'(A)')"Lrandom=F"
    write(1173,'(A)')"#LINRHO"
    write(1173,'(A)')"LINRHO=F"
    write(1173,'(A)')"##########################Isolate system Calculate#################"
    write(1173,'(A)')"#Lattice in periodic boundary condition"
    write(1173,'(A)')"Lpbc=F"
    write(1173,'(A)')"#Auto set the radius"
    write(1173,'(A)')"Lradius_auto = T"
    write(1173,'(A)')"#Length of Cell(cell>2*ISOrmax), Lcell, Lcellbohr"
    write(1173,'(A)')"Lcell = 14"
    write(1173,'(A)')"#The max spere radius, IsoRmax, IsoRmaxbohr"
    write(1173,'(A)')"IsoRmax = 7"
    write(1173,'(A)')"#The order of multi-grids"
    write(1173,'(A)')"NVC=6"
    write(1173,'(A)')"#the order of finite difference"
    write(1173,'(A)')"ISOnorder=8"
    write(1173,'(A)')"#For relative tolerance of CG"
    write(1173,'(A)')"TOLCG=0.1d-2"
    write(1173,'(A)')"#The order of compact finite difference:(4/6)"
    write(1173,'(A)')"NFCD=6"
    write(1173,'(A)')"#The order of  multi-pole"
    write(1173,'(A)')"ISOLmax=9"
    write(1173,'(A)')"#The precision of FastMultipoleMethod[-2,5]~{0.5-0.5E-14}"
    write(1173,'(A)')"FMMPREC=0"
    write(1173,'(A)')"#Poisson solver 0:CG_amg 1:fmm 2:direct 3: MCM 4:cutoff-method 5:MG-CG-method "
    write(1173,'(A)')"#               6:isf-method"
    write(1173,'(A)')"#xxFMM=F"
    write(1173,'(A)')"#xxHartree_direct=F"
    write(1173,'(A)')"#xxMCM_flag=T"
    write(1173,'(A)')"HARTREE_METHOD=6"
    write(1173,'(A)')"#Number of charges for electric system"
    write(1173,'(A)')"#ADDcharge = -1"
    write(1173,'(A)')"#pseudopotential formats option:{PP_OPT -> 0:rrpot,realpot;1:upf}"
    write(1173,'(A)')"PP_OPT=0"
    write(1173,'(A)')"#init method:{0:Slater orbital;1:random orbital; 1:Molecular orbital}"
    write(1173,'(A)')"IDinit=0"
    write(1173,'(A)')"#MO info file"
    write(1173,'(A)')"MO_file=C4.bind.out"
    write(1173,'(A)')"#DEBUG_OUT (INTEGER) value: 0~nothing 1~parallel-chebyshev-for-debug"
    write(1173,'(A)')"DEBUG_OUT=0"
    write(1173,'(A)')"#size of scalapack block"
    write(1173,'(A)')"BLOCK_MBNB = 1"
    write(1173,'(A)')"#Force calculate"
    write(1173,'(A)')"CALFORCE = T"
    write(1173,'(A)')"##########################mixing data#############################"
    write(1173,'(A)')"IMIXER=1"
    write(1173,'(A)')"#the max iter step"
    write(1173,'(A)')"NMITER=200"
    write(1173,'(A)')"#number of simple mixing"
    write(1173,'(A)')"NSMIX=10"
    write(1173,'(A)')"#simple mixing factor"
    write(1173,'(A)')"MALPHA=0.3"
    write(1173,'(A)')"#number of history for Anderson mixing used (2-6)"
    write(1173,'(A)')"NHMIX=8"
    write(1173,'(A)')"#Anderson mixing factor"
    write(1173,'(A)')"MBETA=0.5 #1"
    write(1173,'(A)')"#for ill matrix"
    write(1173,'(A)')"W0AM=1e-3"
    write(1173,'(A)')"#tolarence of drho in scf"
    write(1173,'(A)')"RTOL=1e-5"
    write(1173,'(A)')"#tolerance of total energy(Bohr)"
    write(1173,'(A)')"ETOL=1e-4"
    write(1173,'(A)')"############################OUT-PUT###############################"
    write(1173,'(A)')"#OUTPUT bands"
    write(1173,'(A)')"LBAND=F"
    write(1173,'(A)')"############################THE END###############################"
  ENDSUBROUTINE generate_infile
ENDMODULE parameters
