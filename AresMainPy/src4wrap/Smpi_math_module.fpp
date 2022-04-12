# 1 "Smpi_math_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Smpi_math_module.f90"
MODULE smpi_math_module

   USE CONSTANTS
   !USE mpi

   IMPLICIT NONE

   INCLUDE 'mpif.h'
   type parallel_type  !{{{
      INTEGER(I4B)              :: comm
      INTEGER(I4B)              :: myid
      INTEGER(I4B)              :: numprocs
      INTEGER(I4B)              :: rootid
      LOGICAL                   :: isroot
      !>> number of states needed by per process
      INTEGER(I4B)              :: nstate_proc
      !>> array_map('id in alignment of n_states','process id')
      INTEGER(I4B),ALLOCATABLE  :: sub2sum(:,:)
      INTEGER(I4B)              :: mygrid_range(3)
      INTEGER(I4B),ALLOCATABLE  :: recvcounts(:)
      INTEGER(I4B),ALLOCATABLE  :: displs(:)
      INTEGER(I4B),ALLOCATABLE  :: global_gridrange(:,:)
      INTEGER(I4B)              :: comm2d,commx,commy,rankx,ranky,&
           & periods(2),reorder,remainX(2),remainY(2),ndims,dims(2)
      !> fft information
      INTEGER(I4B)              :: commfft,local_z,local_z_start
      INTEGER(I4B),ALLOCATABLE  :: fft_grid_range(:,:)
      INTEGER(I4B),ALLOCATABLE  :: fft_rcount(:),fft_rdispls(:)
      INTEGER(I4B),ALLOCATABLE  :: fft_scount(:),fft_sdispls(:)
!#ifdef Pencil
      ! INTEGER(I4B)              :: comm_cart
      ! INTEGER(I4B)              :: ndims=2
      ! INTEGER(I4B)              :: dims(2)
      ! INTEGER(I4B)              :: coords(2)
      ! INTEGER(I4B)              :: nprow,npcol
      ! INTEGER(I4B)              :: color
      ! INTEGER(I4B)              :: key
      ! LOGICAL                   :: lperiods(2)
! #ifdef Pencil
      !INTEGER(I4B)              :: buffcoor(
      !INTEGER(I4B),parameter        :: ndims=2
      !INTEGER(I4B)                 :: dims(ndims)
      !INTEGER(I4B)                 :: coords(ndims)
! #endif
   end type parallel_type
   type (parallel_type)         :: parallel
   !-----------------------------------------------------------------------
   TYPE smpi_root_type
      INTEGER(I4B), allocatable :: natom_group(:, :)
   end type smpi_root_type
   type(smpi_root_type)         :: smpi_root
   !-----------------------------------------------------------------------
   TYPE smpi_comm_type
      INTEGER(I4B), allocatable :: atoms(:)
      INTEGER(I4B), allocatable :: displs(:)
   end type smpi_comm_type
   type( smpi_comm_type )       :: smpi_comm
   !-----------------------------------------------------------------------
   type time_type
      character(len=100)       :: label
      real(dp)                  :: tic
      real(dp)                  :: toc
      real(dp)                  :: total
      real(dp)                  :: sum_total
      INTEGER(I4B)              :: num
   end type time_type

   type mem_type
      character(len=100)       :: label
      real(dp)                  :: memic
      real(dp)                  :: total
      INTEGER(I4B)              :: num
   end type mem_type

   type(time_type)              :: timedat(100)
   type(mem_type)              :: memdat(100)
   INTEGER(I4B),private,save    :: ntime=0
   REAL(DP)                     :: rtic,rtoc
   !-----------------------------------------------------------------------
   TYPE grid_diff_map_type
      INTEGER(I4B),allocatable :: nz_map(:,:) !> the up id and down id for per nz
      INTEGER(I4B) :: mycomm_cores(2)         !> number of cores for communcation
      INTEGER(I4B),allocatable :: mycomm_size(:,:)
      INTEGER(I4B),allocatable :: mysend_size(:,:)
      INTEGER(I4B),allocatable :: local_map(:,:)
      INTEGER(I4B),allocatable :: local_map1d(:,:,:)
      INTEGER(I4B)             :: boundary(2,3)
   ENDTYPE grid_diff_map_type
   TYPE sphere_type
      INTEGER(I4B)             :: Length !> sphere size
      INTEGER(I4B),ALLOCATABLE :: map3d(:,:) !> (4,length)
   ENDTYPE sphere_type
   type(sphere_type) :: sphere
   type(grid_diff_map_type) :: diff_map
   !-----------------------------------------------------------------------
   INTEGER(I4B)                 :: mpinfo
   INTEGER(I4B)                 :: smpi_status(MPI_STATUS_SIZE)  !}}}
   LOGICAL  :: Lall_grid=.false.
   !------------------------------- DIVIDER LINE --------------------------------

INTERFACE SompSum  !{{{
   MODULE PROCEDURE sum_real_1d
   MODULE PROCEDURE sum_real_2d
   MODULE PROCEDURE sum_real_3d
   MODULE PROCEDURE sum_cplx_1d
   MODULE PROCEDURE sum_cplx_2d
   MODULE PROCEDURE sum_cplx_3d
END INTERFACE  

INTERFACE SmpiSum  
   MODULE PROCEDURE smpi_sum_int_1s
   MODULE PROCEDURE smpi_sum_cplx_1s
   MODULE PROCEDURE smpi_sum_real_1s
   MODULE PROCEDURE smpi_sum_real_1d
   MODULE PROCEDURE smpi_sum_real_2d
   MODULE PROCEDURE smpi_sum_real_3d
END INTERFACE  

! INTERFACE SmpiSumSub 
!   MODULE PROCEDURE smpi_sum_int_1s_sub
!   MODULE PROCEDURE smpi_sum_cplx_1s_sub
   ! MODULE PROCEDURE smpi_sum_real_1s_sub 
!    MODULE PROCEDURE smpi_sum_real_1d_sub
!    MODULE PROCEDURE smpi_sum_real_2d_sub
!    MODULE PROCEDURE smpi_sum_real_3d_sub
! END INTERFACE  

! INTERFACE SmpiSumSlab
!    !MODULE PROCEDURE smpi_sum_int_1s
!    !MODULE PROCEDURE smpi_sum_cplx_1s
!    MODULE PROCEDURE smpi_sum_real_1s_slab
!    MODULE PROCEDURE smpi_sum_real_1d
!    !MODULE PROCEDURE smpi_sum_real_2d
!    !MODULE PROCEDURE smpi_sum_real_3d
! END INTERFACE  

INTERFACE SmpiSumMem
    MODULE PROCEDURE smpi_sum_mem_1d
    MODULE PROCEDURE smpi_sum_mem_2d
    MODULE PROCEDURE smpi_sum_mem_3d
END INTERFACE
       

INTERFACE SmpiReduceSum 
   MODULE PROCEDURE smpi_reduce_sum_real_1d
   MODULE PROCEDURE smpi_reduce_sum_int_1d
   MODULE PROCEDURE smpi_reduce_sum_cplx_1d
   MODULE PROCEDURE smpi_reduce_sum_real_2d
END INTERFACE  

! INTERFACE SmpiReduceSumSub
!    MODULE PROCEDURE smpi_reduce_sum_real_1d_sub
!    MODULE PROCEDURE smpi_reduce_sum_cplx_1d_sub
!    MODULE PROCEDURE smpi_reduce_sum_real_2d_sub
! !   MODULE PROCEDURE smpi_reduce_sum_cplx_1d
! !   MODULE PROCEDURE smpi_reduce_sum_real_2d
! END INTERFACE  

! INTERFACE SmpiReduceSumSlab
!    MODULE PROCEDURE smpi_reduce_sum_real_1d_slab
!    MODULE PROCEDURE smpi_reduce_sum_cplx_1d_slab
! !   MODULE PROCEDURE smpi_reduce_sum_cplx_1d
! !   MODULE PROCEDURE smpi_reduce_sum_real_2d
! END INTERFACE  

INTERFACE SompSumPow
   module procedure sum_pow_int
   module procedure sum_pow_real
end interface

INTERFACE SmpiSumPow
   module procedure smpi_sum_pow_int
   module procedure smpi_sum_pow_real
end interface!}}}
!---------------------------- DIVIDER LINE -----------------------------
CONTAINS
  SUBROUTINE smpi_init()
    IMPLICIT NONE

    CALL MPI_INIT(mpinfo)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, parallel%myid, mpinfo)
    CALL MPI_Comm_SIZE(MPI_COMM_WORLD, parallel%numprocs, mpinfo)

    parallel%comm  = MPI_COMM_WORLD

    parallel%rootid = 0
    if(parallel%rootid == parallel%myid)then
       parallel%isroot = .true.
    else
       parallel%isroot = .false.
    endif
  END SUBROUTINE smpi_init
  !-----------------------divided line---------------------------
  !-----------------------divided line---------------------------
  SUBROUTINE smpi_init_comm(lpbc)
    IMPLICIT NONE
    LOGICAL,intent(in) :: lpbc

    if(Lpbc)then
       call smpi_init_per()
    else
       call smpi_init_iso()
    endif

  END SUBROUTINE smpi_init_comm
  !-----------------------divided line---------------------------
  SUBROUTINE smpi_init_per()!{{{
    use m_time_evaluate,only:memory_sum,memory_free,filename,&
         &file_unit

    IMPLICIT NONE
    INTEGER(I4B)      :: i

    ! CALL MPI_INIT(mpinfo)

    ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, parallel%myid, mpinfo)
    ! CALL MPI_Comm_SIZE(MPI_COMM_WORLD, parallel%numprocs, mpinfo)

    ! parallel%comm  = MPI_COMM_WORLD

    ! parallel%rootid = 0
    ! if(parallel%rootid == parallel%myid)then
    !    parallel%isroot = .true.
    ! else
    !    parallel%isroot = .false.
    ! endif

    !> 2d split
    parallel%ndims=2
    parallel%dims=(/0,0/)
    parallel%periods=(/1,1/)
    parallel%reorder=0
    parallel%remainX=(/0,1/)
    parallel%remainY=(/1,0/)
    !> get 2d layout
    CALL MPI_DIMS_CREATE(parallel%numprocs,parallel%ndims,parallel%dims,mpinfo)

    if(parallel%isroot)print*,'parallel dims',parallel%dims
    !> get a new communicator attributed with cartesian topology
    CALL MPI_CART_CREATE(parallel%comm, parallel%ndims, parallel%dims, &
         &parallel%periods, parallel%reorder, parallel%comm2d, mpinfo)

    !> set two sub communicator
    CALL MPI_CART_SUB(parallel%comm2d, parallel%remainX, parallel%commX,mpinfo)
    CALL MPI_CART_SUB(parallel%comm2d, parallel%remainY, parallel%commY,mpinfo)

    !> get rank information
    CALL MPI_COMM_RANK(parallel%commx,parallel%rankx,mpinfo)
    CALL MPI_COMM_RANK(parallel%commy,parallel%ranky,mpinfo)

    !> reset memory count file
    write(filename,*)parallel%myid
    filename='memory_sum'//trim(adjustl(filename))
    open(1367,FILE=filename)
    close(1367,STATUS='delete')
    file_unit=file_unit+parallel%myid

    !> grid split init
    parallel%mygrid_range=-2
    allocate(parallel%recvcounts(parallel%dims(2)))
    allocate(parallel%displs(parallel%dims(2)))
    allocate(parallel%global_gridrange(3,parallel%dims(2)))
    CALL memory_sum('smpi_init',(real(size(parallel%recvcounts),DP)+size(parallel%displs)+size(parallel%global_gridrange))*I4B)

  END SUBROUTINE smpi_init_per!}}}
  !-----------------------divided line---------------------------
  SUBROUTINE smpi_init_iso()!{{{
    use m_time_evaluate,only:memory_sum,memory_free,filename&
         &, file_unit

    IMPLICIT NONE
    INTEGER(I4B)      :: i

    ! CALL MPI_INIT(mpinfo)

    ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, parallel%myid, mpinfo)
    ! CALL MPI_Comm_SIZE(MPI_COMM_WORLD, parallel%numprocs, mpinfo)

    ! parallel%comm  = MPI_COMM_WORLD

    ! parallel%rootid = 0
    ! if(parallel%rootid == parallel%myid)then
    !    parallel%isroot = .true.
    ! else
    !    parallel%isroot = .false.
    ! endif
    !> suitable to periodic
    parallel%commx=parallel%comm
    parallel%mygrid_range=-2
    parallel%dims=(/1,parallel%numprocs/)
    allocate(parallel%recvcounts(parallel%numprocs))
    allocate(parallel%displs(parallel%numprocs))
    allocate(parallel%global_gridrange(3,parallel%numprocs))
    CALL memory_sum('smpi_init',(real(size(parallel%recvcounts),DP)+size(parallel%displs)+size(parallel%global_gridrange))*I4B)

    !> reset memory count file
    write(filename,*)parallel%myid
    filename='memory_sum'//trim(adjustl(filename))
    open(1367,FILE=filename)
    close(1367,STATUS='delete')
    file_unit=file_unit+parallel%myid

  END SUBROUTINE smpi_init_iso!}}}
  !-----------------------divided line---------------------------
  ! SUBROUTINE smpi_bcast_input()!{{{
  !     use parameters
  !     implicit none
  !     INTEGER(I4B)                 :: nsend
  !     nsend = 1
  !     call MPI_BCAST(init_gap,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(nstates,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(wsmear,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(malpha,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(kspacing,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(maxnpts,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(LBVK,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(LINRHO,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(Lfirst,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(nssp,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(chem,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(chem0,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(Isym,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(Istart,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(kgrid,3,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(Finite_order,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(lfinite_full,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(finite_order,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(nsmix,1,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     call MPI_BCAST(mbeta,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)

  !     !-----------------------------------------------------------------------
  !     !nsend = size(p_target,1)
  !     !call MPI_BCAST(p_target,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(p_freq,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     nsend = 1
  !     !call MPI_BCAST(t_freq,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(t_target,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(pmass,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(pdrag,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(tdrag,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(mstrain,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(erate,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !-----------------------------------------------------------------------
  !     !nsend = CLen
  !     !call MPI_BCAST(input_file,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(system_name,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(denfile_name,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(thermostat,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(integrator,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(ensemble,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(lwstyle,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(press_control,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(pcontrol,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(sdir,nsend,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
  !     !print *, 'debug lwstyle', parallel%myid, lwstyle
  !     !IF (.NOT.ALLOCATED(qmass)) ALLOCATE(qmass(nnhc))
  !     !IF (.NOT.ALLOCATED(bmass)) ALLOCATE(bmass(nnhc))
  !     !nsend=nnhc
  !     !call MPI_BCAST(qmass,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(bmass,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     nsend = 1
  !     call MPI_BCAST(LOPT,nsend,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(IOPM,nsend,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(IGOAL,nsend,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(IBRON,nsend,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(maxsave,nsend,MPI_INTEGER,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(TOLF,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(TOLP,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(Press,nsend,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !call MPI_BCAST(Ptensor,6,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
  !     !print *, 'debug qmass', parallel%myid, qmass

  !  end subroutine smpi_bcast_input!}}}

   SUBROUTINE smpi_exit()!{{{
      IMPLICIT NONE
      call MPI_finalize( mpinfo )
      STOP
      ! call exit()
      RETURN
   END SUBROUTINE!}}}

   SUBROUTINE smpi_stop(message)!{{{
      IMPLICIT NONE

      CHARACTER (LEN=*) message

      WRITE (*,*) message

      call MPI_abort(parallel%comm, 1, mpinfo )
      STOP

      RETURN
   END SUBROUTINE!}}}

   SUBROUTINE smpi_stop_info(message)!{{{
      IMPLICIT NONE

      CHARACTER (LEN=*) message

      WRITE (*,*) message, mpinfo

      call MPI_abort(MPI_comm_world , 1, mpinfo )
      STOP

      RETURN
   END SUBROUTINE!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   Subroutine  nstates_split(m,np)  !{{{
   implicit none 
   integer(i4b) , intent(inout) :: m
   integer(i4b) :: ms,np
   !-------------------------------------
   ! ms = mod(m ,np)
   ! if (ms == 0 ) then 
   !   m = m /np
   ! else 
   !   !if ( ms >= (parallel%myid + 1) ) then 
   !   !  m = m/parallel%numprocs+1
   !   !else 
   !   !  m = m / parallel%numprocs
   !   !end if 
   !   !if (parallel%myid == 0 ) then 
   !   if (parallel%coords(1) == 0 ) then 
   !     m = m/np + ms 
   !   else 
   !     m = m / np
   !   end if 
   ! end if 
   
   End Subroutine nstates_split  !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   Subroutine  nstates_split_2(m,np)  !{{{
   implicit none 
   integer(i4b) , intent(inout) :: m
   integer(i4b) :: ms,np
   !-------------------------------------
   ! ms = mod(m ,np)
   ! if (ms == 0 ) then 
   !   m = m /np
   ! else 
   !   !if ( ms >= (parallel%myid + 1) ) then 
   !   !  m = m/parallel%numprocs+1
   !   !else 
   !   !  m = m / parallel%numprocs
   !   !end if 
   !   !if (parallel%myid == 0 ) then 
   !   if (parallel%coords(2) == 0 ) then 
   !     m = m/np + ms 
   !   else 
   !     m = m / np
   !   end if 
   ! end if 
   
   End Subroutine nstates_split_2  !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_real_1d(amat) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_cplx_1d(amat) result(totals)!{{{
      complex(DP)                      :: totals
      complex(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_real_2d(amat, bmat) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_cplx_2d(amat, bmat) result(totals)!{{{
      complex(DP)                      :: totals
      complex(DP)                      :: amat(:,:,:)
      complex(DP)                      :: bmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_real_3d( amat, bmat, cmat) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:),cmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)*cmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_cplx_3d( amat, bmat, cmat) result(totals)!{{{
      complex(DP)                      :: totals
      complex(DP)                      :: amat(:,:,:)
      complex(DP)                      :: bmat(:,:,:),cmat(:,:,:)
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)*bmat(i,j,k)*cmat(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_pow_real(amat, pow) result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: pow
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)** pow
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION sum_pow_int(amat, pow)  result(totals)!{{{
      REAL(DP)                      :: totals
      REAL(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: pow
      INTEGER(I4B)                  :: n1,n2,n3
      INTEGER(I4B)                  :: i,j,k
      n1 = size(amat,1)
      n2 = size(amat,2)
      n3 = size(amat,3)
      totals = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:totals) PRIVATE(k,j,i)
      DO k = 1, n3
         DO j = 1, n2
            DO i = 1, n1
               totals = totals + amat(i,j,k)** pow
            ENDDO
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function smpi_sum_int_1s(x) result(sumx)!{{{

      IMPLICIT NONE
      INTEGER(I4B)  :: x,sumx

      CALL MPI_ALLREDUCE(x, sumx, 1, MPI_INTEGER4, MPI_SUM, parallel%comm, mpinfo)

   end function smpi_sum_int_1s!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function smpi_sum_cplx_1s(x) result(sumx)!{{{

      IMPLICIT NONE
      complex(DP)                     :: x,sumx

      CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%comm, mpinfo)

   end function smpi_sum_cplx_1s!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   function smpi_sum_real_1s(x) result(sumx)!{{{

      IMPLICIT NONE
      REAL(DP)                     :: x,sumx

      CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%comm, mpinfo)

   end function smpi_sum_real_1s!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! function smpi_sum_real_1s_sub(x) result(sumx)!{{{

   !    IMPLICIT NONE
   !    REAL(DP)                     :: x,sumx

   !    CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%subcomm, mpinfo)

   ! end function smpi_sum_real_1s_sub!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! function smpi_sum_real_1s_slab(x) result(sumx)!{{{

   !    IMPLICIT NONE
   !    REAL(DP)                     :: x,sumx

   !    CALL MPI_ALLREDUCE(x, sumx, 1, MPI_REAL8, MPI_SUM, parallel%slabcomm, mpinfo)

   ! end function smpi_sum_real_1s_slab!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_real_1d(amat) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      totals = sum_real_1d(amat)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! FUNCTION smpi_sum_real_1d_sub(amat) result(suma)!{{{
   !    REAL(DP)                      :: suma, totals
   !    REAL(DP)                      :: amat(:,:,:)
   !    totals = sum_real_1d(amat)
   !    call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   ! END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_real_2d(amat, bmat) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:)
      totals = sum_real_2d(amat,bmat)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! FUNCTION smpi_sum_real_2d_sub(amat, bmat) result(suma)!{{{
   !    REAL(DP)                      :: suma, totals
   !    REAL(DP)                      :: amat(:,:,:)
   !    REAL(DP)                      :: bmat(:,:,:)
   !    totals = sum_real_2d(amat,bmat)
   !    call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   ! END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_real_3d( amat, bmat, cmat) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: bmat(:,:,:),cmat(:,:,:)
      totals = sum_real_2d(amat,bmat)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! FUNCTION smpi_sum_real_3d_sub( amat, bmat, cmat) result(suma)!{{{
   !    REAL(DP)                      :: suma, totals
   !    REAL(DP)                      :: amat(:,:,:)
   !    REAL(DP)                      :: bmat(:,:,:),cmat(:,:,:)
   !    totals = sum_real_2d(amat,bmat)
   !    call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   ! END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_pow_real(amat, pow) result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      REAL(DP)                      :: pow
      totals = sum_pow_real(amat,pow)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_pow_int(amat, pow)  result(suma)!{{{
      REAL(DP)                      :: suma, totals
      REAL(DP)                      :: amat(:,:,:)
      INTEGER(I4B)                  :: pow
      totals = sum_pow_int(amat,pow)
      call mpi_allreduce(totals, suma, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)
   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
!   SUBROUTINE smpi_split_natom(natom, atomid_local)!{{{
!      IMPLICIT NONE
!      INTEGER(I4B)                  :: na, nm, i
!      INTEGER(I4B)                  :: natom, atomid_local(2)
!      if (parallel%isroot) then
!         allocate(smpi_root%natom_group(2, parallel%numprocs))
!      endif
!      na = natom / parallel%numprocs
!      nm = natom - na * parallel%numprocs
!      if (parallel%isroot) then
!         atomid_local(1) = 1
!         atomid_local(2) = na 
!      else if (parallel%myid < nm + 1) then
!         atomid_local(1) = na * parallel%myid + parallel%myid
!         atomid_local(2) = atomid_local(1) + na
!      else
!         atomid_local(1) = na * parallel%myid + nm + 1
!         atomid_local(2) = atomid_local(1) + na -1
!      endif
!      if (parallel%isroot) then
!         smpi_root%natom_group(1, 1) = 1
!         smpi_root%natom_group(2, 1) = na
!         if (nm==0 ) then
!            do i = 2, parallel%numprocs
!               smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
!               smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na - 1
!            enddo
!         else
!            do i = 2, nm
!               smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
!               smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na
!            enddo
!            do i = nm + 1, parallel%numprocs
!               smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
!               smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na - 1
!            enddo
!         endif
!      endif
!   end SUBROUTINE!}}}
   !---------------------------- DIVIDER LINE -----------------------------
!    SUBROUTINE smpi_split_natom(natom, atomid_local,np)!{{{
!       IMPLICIT NONE
!       INTEGER(I4B)                  :: na, nm, i,np
!       INTEGER(I4B)                  :: natom, atomid_local(2)
!       !if (parallel%isroot) then
!       if (parallel%coords(2) == 0) then
!          allocate(smpi_root%natom_group(2, np))
!       endif
!       na = natom / np
!       nm = natom - na * np
!       !if (parallel%isroot) then
!       if (parallel%coords(2) == 0) then
!          atomid_local(1) = 1
!          atomid_local(2) = na 
!       else if (parallel%mysubid < nm + 1) then
!          atomid_local(1) = na * parallel%mysubid + parallel%mysubid
!          atomid_local(2) = atomid_local(1) + na
!       else
!          atomid_local(1) = na * parallel%mysubid + nm + 1
!          atomid_local(2) = atomid_local(1) + na -1
!       endif
!       !if (parallel%isroot) then
!       if (parallel%coords(2) == 0) then
!          smpi_root%natom_group(1, 1) = 1
!          smpi_root%natom_group(2, 1) = na
!          if (nm==0 ) then
!             do i = 2, np
!                smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
!                smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na - 1
!             enddo
!          else
!             do i = 2, nm
!                smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
!                smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na
!             enddo
!             do i = nm + 1, np
!                smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
!                smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na - 1
!             enddo
!          endif
!       endif
! !print *,'atomid local',atomid_local,parallel%coords
!    end SUBROUTINE!}}}
   ! SUBROUTINE smpi_split_init()!{{{
   !   IMPLICIT NONE
   !   INTEGER(I4B)                  :: na, nm, i,np
   !   INTEGER(I4B)                  :: natom, atomid_local(2)
   !   !if (parallel%isroot) then
   !   if (parallel%coords(2) == 0) then
   !      allocate(smpi_root%natom_group(2, np))
   !   endif
   !   na = natom / np
   !   nm = natom - na * np
   !   !if (parallel%isroot) then
   !   if (parallel%coords(2) == 0) then
   !      atomid_local(1) = 1
   !      atomid_local(2) = na 
   !   else if (parallel%mysubid < nm + 1) then
   !      atomid_local(1) = na * parallel%mysubid + parallel%mysubid
   !      atomid_local(2) = atomid_local(1) + na
   !   else
   !      atomid_local(1) = na * parallel%mysubid + nm + 1
   !      atomid_local(2) = atomid_local(1) + na -1
   !   endif
   !   !if (parallel%isroot) then
   !   if (parallel%coords(2) == 0) then
   !      smpi_root%natom_group(1, 1) = 1
   !      smpi_root%natom_group(2, 1) = na
   !      if (nm==0 ) then
   !         do i = 2, np
   !            smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
   !            smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na - 1
   !         enddo
   !      else
   !         do i = 2, nm
   !            smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
   !            smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na
   !         enddo
   !         do i = nm + 1, np
   !            smpi_root%natom_group(1, i) = smpi_root%natom_group(2, i - 1)+1
   !            smpi_root%natom_group(2, i) = smpi_root%natom_group(1, i)+ na - 1
   !         enddo
   !      endif
   !   endif
   !   !print *,'atomid local',atomid_local,parallel%coords
   ! end SUBROUTINE smpi_split_init!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_int_1d(amat,ramat) !{{{
      INTEGER(I4B)                      :: amat(:)
      INTEGER(I4B),OPTIONAL             :: ramat(:)
      INTEGER(I4B), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      na = size(amat)
      call mpi_allreduce(amat, asuma, na, mpi_integer4, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_real_1d(amat,ramat) !{{{
      REAL(DP)                      :: amat(:)
      REAL(DP),OPTIONAL            :: ramat(:)
      REAL(DP), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      na = size(amat)
      call mpi_allreduce(amat, asuma, na, mpi_real8, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   SUBROUTINE smpi_reduce_sum_real(amat,na,ramat) !{{{
      REAL(DP)                      :: amat(na)
      REAL(DP),OPTIONAL            :: ramat(na)
      REAL(DP), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      call mpi_allreduce(amat, asuma, na, mpi_real8, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! SUBROUTINE smpi_reduce_sum_real_1d_sub(amat,ramat) !{{{
   !    REAL(DP)                      :: amat(:)
   !    REAL(DP),OPTIONAL            :: ramat(:)
   !    REAL(DP), dimension(size(amat)) :: asuma
   !    INTEGER(I4B)                 :: na
   !    na = size(amat)
   !    call mpi_allreduce(amat, asuma, na, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   !    if (present(ramat)) then
   !       ramat=asuma
   !    else
   !       amat=asuma
   !    endif
   ! END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! SUBROUTINE smpi_reduce_sum_cplx_1d_sub(amat,ramat) !{{{
   !    COMPLEX(DP)                      :: amat(:)
   !    COMPLEX(DP),OPTIONAL            :: ramat(:)
   !    COMPLEX(DP), dimension(size(amat)) :: asuma
   !    INTEGER(I4B)                 :: na
   !    na = size(amat)
   !    call mpi_allreduce(amat, asuma, na, mpi_complex16, mpi_sum, parallel%subcomm, mpinfo)
   !    if (present(ramat)) then
   !       ramat=asuma
   !    else
   !       amat=asuma
   !    endif
   ! END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! SUBROUTINE smpi_reduce_sum_real_1d_slab(amat,ramat) !{{{
   !    REAL(DP)                      :: amat(:)
   !    REAL(DP),OPTIONAL            :: ramat(:)
   !    REAL(DP), dimension(size(amat)) :: asuma
   !    INTEGER(I4B)                 :: na
   !    na = size(amat)
   !    call mpi_allreduce(amat, asuma, na, mpi_real8, mpi_sum, parallel%slabcomm, mpinfo)
   !    if (present(ramat)) then
   !       ramat=asuma
   !    else
   !       amat=asuma
   !    endif
   ! END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
!   SUBROUTINE smpi_reduce_sum_cplx_1d_slab(amat,ramat) !{{{
!      COMPLEX(DP)                      :: amat(:)
!      COMPLEX(DP),OPTIONAL            :: ramat(:)
!      COMPLEX(DP), dimension(size(amat)) :: asuma
!      INTEGER(I4B)                 :: na
!      na = size(amat)
!      call mpi_allreduce(amat, asuma, na, mpi_complex16, mpi_sum, parallel%slabcomm, mpinfo)
!      if (present(ramat)) then
!         ramat=asuma
!      else
!         amat=asuma
!      endif
!   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! SUBROUTINE smpi_reduce_sum_cplx_1d_slab(amat,ramat) !{{{
   !    COMPLEX(DP)                      :: amat(:)
   !    COMPLEX(DP),OPTIONAL            :: ramat(:)
   !    COMPLEX(DP), dimension(size(amat)) :: asuma
   !    INTEGER(I4B)                 :: na
   !    na = size(amat)
   !    call mpi_allreduce(amat, asuma, na, mpi_complex16, mpi_sum, parallel%slabcomm, mpinfo)
   !    if (present(ramat)) then
   !       ramat=asuma
   !    else
   !       amat=asuma
   !    endif
   ! END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_cplx_1d(amat,ramat) !{{{
      complex(DP)                      :: amat(:)
      complex(DP),OPTIONAL            :: ramat(:)
      complex(DP), dimension(size(amat)) :: asuma
      INTEGER(I4B)                 :: na
      na = size(amat)
      call mpi_allreduce(amat, asuma, na, MPI_COMPLEX16, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   SUBROUTINE smpi_reduce_sum_real_2d(amat,ramat) !{{{
      REAL(DP)                      :: amat(:,:)
      REAL(DP),OPTIONAL            :: ramat(:,:)
      REAL(DP), dimension(size(amat,1),size(amat,2)) :: asuma
      INTEGER(I4B)                 :: nx,ny,nxy
      nxy = size(amat,1) * size(amat,2)
      call mpi_allreduce(amat, asuma, nxy, mpi_real8, mpi_sum, parallel%comm, mpinfo)
      if (present(ramat)) then
         ramat=asuma
      else
         amat=asuma
      endif
   END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   ! SUBROUTINE smpi_reduce_sum_real_2d_sub(amat,ramat) !{{{
   !    REAL(DP)                      :: amat(:,:)
   !    REAL(DP),OPTIONAL            :: ramat(:,:)
   !    REAL(DP), dimension(size(amat,1),size(amat,2)) :: asuma
   !    INTEGER(I4B)                 :: nx,ny,nxy
   !    nxy = size(amat,1) * size(amat,2)
   !    call mpi_allreduce(amat, asuma, nxy, mpi_real8, mpi_sum, parallel%subcomm, mpinfo)
   !    if (present(ramat)) then
   !       ramat=asuma
   !    else
   !       amat=asuma
   !    endif
   ! END SUBROUTINE !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine start_time(inlabel,flag,tic)!{{{
      implicit none
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL  :: flag
      REAL(DP),optional :: tic
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            timedat(i)%tic= mpi_wtime()
            if (present(tic)) tic = timedat(i)%tic
            timedat(i)%num=timedat(i)%num+1
            return
         endif
      enddo
      ntime=ntime+1
      timedat(ntime)%label = label
      timedat(ntime)%tic= mpi_wtime()
      timedat(ntime)%sum_total= 0.d0
      if (present(tic)) tic = timedat(ntime)%tic
   end subroutine !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine end_time(inlabel,flag,toc)!{{{
      implicit none
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL  :: flag
      REAL(DP),optional :: toc
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            timedat(i)%toc= mpi_wtime()
            if (present(toc)) toc = timedat(i)%toc
            timedat(i)%total = timedat(i)%toc - timedat(i)%tic
            timedat(i)%sum_total = timedat(i)%sum_total + timedat(i)%total
            return
         endif
      enddo
      if (parallel%isroot) then
         WRITE(6,*) "ERROR: Haven't found ", inlabel
      endif
   end subroutine !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine write_time(inlabel,flag)!{{{
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL            :: flag
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      ! if (.not. parallel%isroot) return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            WRITE(6,'(A30,2X,ES12.4,2X,I10,I3)') inlabel, timedat(i)%total, timedat(i)%num,parallel%myid
            return
         endif
      enddo
   end subroutine write_time!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine write_sum_time(inlabel,flag)!{{{
      character(len=*) :: inlabel
      character(len=100) :: label
      LOGICAL            :: flag
      INTEGER(I4B)                 :: i

      if(.not.flag)return
      ! if (.not. parallel%isroot) return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then
            WRITE(6,*) inlabel," >> ",timedat(i)%sum_total
            return
         endif
      enddo
    end subroutine write_sum_time!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   subroutine print_time(inlabel,t)!{{{
      character(len=*) :: inlabel
      character(len=100) :: label
      real(dp) :: t
      INTEGER(I4B)                 :: i
      if (.not. parallel%isroot) return
      label = trim(inlabel)
      do i = 1, ntime
         if (label==timedat(i)%label) then

            t = timedat(i)%toc-timedat(i)%tic
            return
         endif
      enddo
   end subroutine print_time!}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_mem_1d(munit,amat) result(summem)!{{{
      REAL(DP)                      :: summem, mem1
      REAL(DP)                      :: amat(:)
      character(8) :: munit
      if (munit == 'G') then
        mem1 = size(amat,1)*8/1024/1024/1024
      else if (munit == 'M') then
        mem1 = size(amat,1)*8/1024/1024
      end if 

      call mpi_allreduce(mem1, summem, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_mem_2d(munit,amat) result(summem)!{{{
      REAL(DP)                      :: summem, mem1
      REAL(DP)                      :: amat(:,:)
      character(8) :: munit
      if (munit == 'G' .or. munit == 'g' ) then
        mem1 = size(amat,1)*size(amat,2)*8/1024/1024/1024
      else if (munit == 'M' .or. munit == 'm' ) then
        mem1 = size(amat,1)*size(amat,2)*8/1024/1024
      end if 

      call mpi_allreduce(mem1, summem, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

   END FUNCTION !}}}
   !---------------------------- DIVIDER LINE -----------------------------
   FUNCTION smpi_sum_mem_3d(munit,amat) result(summem)!{{{
      REAL(DP)                      :: summem, mem1
      REAL(DP)                      :: amat(:,:,:)
      character(8) :: munit
      if (munit == 'G' .or. munit == 'g' ) then
        mem1 = size(amat,1)*size(amat,2)*size(amat,3)*8/1024/1024/1024
      else if (munit == 'M' .or. munit == 'm' ) then
        mem1 = size(amat,1)*size(amat,2)*size(amat,3)*8/1024/1024
      end if 

      call mpi_allreduce(mem1, summem, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

   END FUNCTION smpi_sum_mem_3d !}}}
   !-------------------------------------------------------------------------
   Subroutine states_split(nev)
     USE parameters , ONLY: BLOCK_MBNB
     IMPLICIT NONE
     !>> total states to calculate
     INTEGER(I4B)             :: ncore
     INTEGER(I4B),INTENT(INOUT)  ::  nev
     !>
     INTEGER(I4B)             :: average,redund
     !>=================================
     !>> spili eigen solver
     average=nev/(parallel%dims(1)*BLOCK_MBNB)
     IF(average==0)THEN
        IF(parallel%isroot)WRITE(6,*)"[ ERROR ]: states not enough, more states or less cores"
        stop
     ENDIF
     redund=mod(nev,parallel%dims(1)*BLOCK_MBNB)
     IF(parallel%myid*BLOCK_MBNB.lt.redund)THEN
        nev=average*BLOCK_MBNB+min(BLOCK_MBNB,redund-parallel%myid*BLOCK_MBNB)
     ELSE
        nev=average*BLOCK_MBNB
     ENDIF
   END Subroutine states_split
   !---------------------------------------------------------------------
   Subroutine array_split(nev)
     USE parameters , ONLY: BLOCK_MBNB
     IMPLICIT NONE
     !>> total states to calculate
     INTEGER(I4B),INTENT(IN)  ::  nev
     INTEGER(I4B)             :: average,redund
     ! INTEGER(I4B)             :: i,j,counter,n_temp
     INTEGER(I4B) :: i,j,n_temp,seq_local,seq_global
     !>=================================
     !>> spili eigen solver
     average=nev/(parallel%dims(1)*BLOCK_MBNB)
     IF(average==0)THEN
        IF(parallel%isroot)WRITE(6,*)"[ ERROR ]: smaller block or parallel cores"
        stop
     ENDIF
     redund=mod(nev,parallel%dims(1)*BLOCK_MBNB)
     IF(parallel%ranky*BLOCK_MBNB.lt.redund)THEN
        parallel%nstate_proc=average*BLOCK_MBNB+min(BLOCK_MBNB,redund-parallel%ranky*BLOCK_MBNB)
     ELSE
        parallel%nstate_proc=average*BLOCK_MBNB
     ENDIF
     !>> set array map associated processes
     ALLOCATE(parallel%sub2sum((average+1)*BLOCK_MBNB,parallel%dims(1)))
     !>> set initial value
     parallel%sub2sum=-2
     seq_local=0
     seq_global=0
     setmap:DO j=1,nev,1
        seq_local = seq_local + 1
        DO i=1,parallel%dims(1),1
           seq_global = seq_global + 1
           parallel%sub2sum(seq_local,i) = seq_global
           if(seq_global == nev)exit setmap
        ENDDO
     ENDDO setmap
     ! counter=0
     ! !>> assignment
     ! DO j=1,parallel%numprocs,1
     !    IF(parallel%myid.eq.j-1)n_temp=parallel%nstate_proc
     !    CALL MPI_BCAST(n_temp,1,MPI_INTEGER4,j-1,parallel%comm,mpinfo)
     ! DO i=1,(average+1)*BLOCK_MBNB,1
     !    IF(i.le.n_temp)THEN
     !       counter=counter+1
     !       parallel%sub2sum(i,j)=counter
     !    ENDIF
     ! ENDDO
     ! ENDDO
     ! print *,"---->",parallel%nstate_proc,parallel%myid
     ! IF(parallel%isroot)print *,"sub2sum",parallel%sub2sum
   ENDSUBROUTINE array_split
   !---------------------------------------------------------------------
   SUBROUTINE grid_split(ngrid,ncore,comm,id,grid_range,recvcounts,displs,gridrange_sum,n1,n2,n3,n)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN)  :: ngrid,ncore,comm,id
     INTEGER(I4B),INTENT(OUT) :: grid_range(3)
     INTEGER(I4B),INTENT(OUT) :: recvcounts(ncore)
     INTEGER(I4B),INTENT(OUT) :: displs(ncore)
     INTEGER(I4B),INTENT(OUT),optional :: gridrange_sum(3,ncore)
    INTEGER(I4B),INTENT(INOUT),optional :: n1,n2,n3,n
     !> local
     INTEGER(I4B) :: i,j
     INTEGER(I4B) :: average,redund
     INTEGER(I4B) :: recvcounts_temp(ncore)
     INTEGER(I4B) :: displs_temp(ncore)

     !> split the grid into "ncore" cores
     average=ngrid/ncore
     redund=mod(ngrid,ncore)
     if(id.lt.redund)then
        !> id*(l+1) + 1, ..., (id*(l+1)+1) + (l+1) - 1
        grid_range(1)=id*(average+1)+1
        grid_range(2)=(id+1)*(average+1)
        grid_range(3)=average+1
     else
        !> redund*(l+1)+(id-redund)*l+1, ..., (redund*(l+1)+(id-redund)*l+1) + l - 1
        grid_range(1)=redund+id*average+1
        grid_range(2)=redund+(id+1)*average
        grid_range(3)=average
     endif

     !> set the recvcounts and displs
     do i=0,ncore-1,1
        j=i+1
        if(i.lt.redund)then
           !> id*(l+1) + 1, ..., (id*(l+1)+1) + (l+1) - 1
           recvcounts(j)=average+1
           displs(j)=i*(average+1)+1-1
        else
           !> redund*(l+1)+(id-redund)*l+1, ..., (redund*(l+1)+(id-redund)*l+1) + l - 1
           recvcounts(j)=average
           displs(j)=redund+i*average+1-1
        endif
     enddo

     !> set the global_gridrange
     if(present(gridrange_sum))then
        recvcounts_temp=3
        displs_temp=(/(i*3,i=0,ncore-1,1)/)
        CALL MPI_ALLGATHERV(grid_range,3,MPI_INTEGER4,gridrange_sum,recvcounts_temp,&
             displs_temp,MPI_INTEGER4,comm,mpinfo)
     endif

     !> set the local cubic grid
     if(present(n1).and.present(n2).and.present(n3))then
        n3=(grid_range(2)-1)/n1/n2+1-(grid_range(1)-1)/n1/n2-1+1
        n=grid_range(3)
     endif

     !> small test
     ! print *,"myid",id,"start",grid_range(1),"end",grid_range(2),"size",grid_range(3)
     ! if(parallel%isroot)print*,"displs",displs
     ! if(parallel%isroot)print*,"recvcounts",recvcounts
   ENDSUBROUTINE grid_split
   !-------------------------------------------------------------------
   Subroutine atom_split(mysize,natom,atom_index)
     IMPLICIT NONE
     INTEGER(I4B)             :: mysize,natom
     INTEGER(I4B)             :: atom_index(:)
     !> local
     INTEGER(I4B)             :: average,redund
     INTEGER(I4B)             :: i,j,counter,n_temp
     !>=================================
     n_temp=0
     !>> splicit atoms
     average=natom/(parallel%numprocs)
     IF(average==0)THEN
        IF(parallel%isroot)WRITE(6,*)"[ HINT ]: exist zero in the number of atoms of cores"
        ! stop
     ENDIF
     redund=mod(natom,parallel%numprocs)
     IF(parallel%myid.lt.redund)THEN
        mysize=average+1
     ELSE
        mysize=average
     ENDIF
     !>> set array map associated processes
     !ALLOCATE(atom_index(average+1))
     !>> set initial value
     atom_index=-2
     counter=0
     !>> assignment
     DO j=1,parallel%numprocs,1
        IF(parallel%myid.eq.j-1)n_temp=mysize
        CALL MPI_BCAST(n_temp,1,MPI_INTEGER4,j-1,parallel%comm,mpinfo)
        DO i=1,n_temp,1
           counter=counter+1
           IF(parallel%myid.eq.j-1)atom_index(i)=counter
        ENDDO
     ENDDO
     ! print *,"---->",parallel%nstate_proc,parallel%myid
     ! IF(parallel%isroot)print *,"sub2sum",parallel%sub2sum
   ENDSUBROUTINE atom_split

   ! SUBROUTINE grid_sphere_init(n1,n2,n3,norder,grid_pos,orign,radius,sphere,diff_map)
   SUBROUTINE grid_sphere_init(n1,n2,n3,norder)
     USE m_time_evaluate, ONLY: memory_sum
     IMPLICIT NONE
     INTEGER(I4B) :: n1,n2,n3
     ! REAL(DP)     :: grid_pos(n1*n2*n3,4)
     ! REAL(DP)     :: orign(3),radius
     INTEGER(I4B) :: norder
     !> local
     INTEGER(I4B) :: i,j,k,in
     REAL(DP)     :: vec_r(3)
     ! REAL(DP)     :: map(4,n1*n2*n3)
     INTEGER(I4B) :: local_z,up_z,down_z
     ! REAL(DP)     :: radius_g

     ! allocate(diff_map%nz_map(2,n3))
     ! in=0
     ! sphere%Length=0
     ! do k=1,n3,1
     !    diff_map%nz_map(1,k)=sphere%length+1
     !    do j=1,n2,1
     !       do i=1,n1,1
     !          in=in+1
     !          vec_r=grid_pos(in,1:3)-orign
     !          ! radius_g=ddot(3,vec_r,1,vec_r,1)
     !          radius_g=sqrt(sum(vec_r**2))
     !          if(radius_g < radius) then
     !             sphere%length=sphere%Length+1
     !             map(1,sphere%length)=i
     !             map(2,sphere%length)=j
     !             map(3,sphere%length)=k
     !             map(4,sphere%length)=in
     !          endif
     !       enddo
     !    enddo
     !    diff_map%nz_map(2,k)=sphere%length
     !    if(diff_map%nz_map(1,k) > diff_map%nz_map(2,k))diff_map%nz_map(1,k)=diff_map%nz_map(2,k)
     ! enddo

     ! !> assign the sphere
     ! allocate(sphere%map3d(4,sphere%Length))
     ! sphere%map3d=map(:,1:sphere%Length)
     ! CALL grid_split(sphere%Length,parallel%numprocs,parallel%myid,parallel%mygrid_range,&
     !      & parallel%recvcounts,parallel%displs,parallel%global_gridrange)
     ! allocate(sphere%rhoS(parallel%mygrid_range(3),1))
     ! sphere%rhoS=parallel%myid

     !> set the difference communcation list
     !> find up cores
     local_z=sphere%map3d(3,parallel%mygrid_range(2))
     do i=parallel%myid+1,parallel%dims(2),1
        up_z=sphere%map3d(3,parallel%global_gridrange(1,i))
        if( local_z+norder >= up_z )then
           diff_map%mycomm_cores(1)=i
        endif
     enddo
     !> find dowm cores
     local_z=sphere%map3d(3,parallel%mygrid_range(1))
     do i=parallel%myid+1,1,-1
        down_z=sphere%map3d(3,parallel%global_gridrange(2,i))
        if( local_z-norder <= down_z )then
           diff_map%mycomm_cores(2)=i
        endif
     enddo
     ! if(parallel%isroot)WRITE(6,'(A6,X,100(I8,2X))')'nz_map',mydiff_map%nz_map
     ! WRITE(6,'(A6,X,100(I8,2X))')'_cores',mydiff_map%mycomm_cores,parallel%myid
     ! CALL MPI_BARRIER(parallel%comm,mpinfo)
     ! stop

     !> set the size for receive
     down_z=sphere%map3d(3,parallel%mygrid_range(1))-norder
     up_z=sphere%map3d(3,parallel%mygrid_range(2))+norder
     if(down_z<sphere%map3d(3,1))down_z=sphere%map3d(3,1)
     if(up_z>sphere%map3d(3,sphere%length))up_z=n3
     allocate(diff_map%mycomm_size(2,diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     call memory_sum('diff_map_mycomm_size',real(size(diff_map%mycomm_size),DP)*I4B)
     do i=diff_map%mycomm_cores(2),diff_map%mycomm_cores(1),1
        if( parallel%global_gridrange(1,i) < diff_map%nz_map(1,down_z) )then
           diff_map%mycomm_size(1,i)=diff_map%nz_map(1,down_z)
           diff_map%mycomm_size(2,i)=parallel%global_gridrange(2,i)
        elseif( parallel%global_gridrange(2,i) > diff_map%nz_map(2,up_z) )then
           diff_map%mycomm_size(1,i)=parallel%global_gridrange(1,i)
           diff_map%mycomm_size(2,i)=diff_map%nz_map(2,up_z)
        else
           diff_map%mycomm_size(1,i)=parallel%global_gridrange(1,i)
           diff_map%mycomm_size(2,i)=parallel%global_gridrange(2,i)
        endif
     enddo

     !> set the size for send
     allocate(diff_map%mysend_size(2,diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     call memory_sum('diff_map_mysend_size',real(size(diff_map%mysend_size),DP)*I4B)
     do i=diff_map%mycomm_cores(2),diff_map%mycomm_cores(1),1
        down_z=sphere%map3d(3,parallel%global_gridrange(2,i))+norder
        up_z=sphere%map3d(3,parallel%global_gridrange(1,i))-norder
        if(down_z>sphere%map3d(3,sphere%length))down_z=sphere%map3d(3,sphere%length)
        if(up_z<sphere%map3d(3,1))up_z=sphere%map3d(3,1)
        if( diff_map%nz_map(2,down_z) < parallel%mygrid_range(2) )then
           diff_map%mysend_size(1,i)=1!parallel%mygrid_range(1)
           diff_map%mysend_size(2,i)=diff_map%nz_map(2,down_z)-parallel%mygrid_range(1)+1
        elseif( diff_map%nz_map(1,up_z) > parallel%mygrid_range(1) )then
           diff_map%mysend_size(1,i)=diff_map%nz_map(1,up_z)-parallel%mygrid_range(1)+1
           diff_map%mysend_size(2,i)=parallel%mygrid_range(2)-parallel%mygrid_range(1)+1
        else
           diff_map%mysend_size(1,i)=1!parallel%mygrid_range(1)
           diff_map%mysend_size(2,i)=parallel%mygrid_range(3)
        endif
     enddo

     !> save the 3d position for wrapped box
     down_z=sphere%map3d(3,parallel%mygrid_range(1))-norder
     up_z=sphere%map3d(3,parallel%mygrid_range(2))+norder
     if(down_z<sphere%map3d(3,1))down_z=sphere%map3d(3,1)
     if(up_z>sphere%map3d(3,sphere%length))up_z=sphere%map3d(3,sphere%length)
     ! diff_map%len_comm=diff_map%nz_map(2,up_z)-diff_map%nz_map(1,down_z)
     allocate(diff_map%local_map(4,diff_map%nz_map(1,down_z):diff_map%nz_map(2,up_z)))
     call memory_sum('diff_map_local_map',real(size(diff_map%local_map),DP)*I4B)
     diff_map%local_map=sphere%map3d(:,diff_map%nz_map(1,down_z):diff_map%nz_map(2,up_z))

     !> get the wrap-grid boundary
     diff_map%boundary(1,1)=minval(diff_map%local_map(1,:))-norder
     diff_map%boundary(2,1)=maxval(diff_map%local_map(1,:))+norder
     diff_map%boundary(1,2)=minval(diff_map%local_map(2,:))-norder
     diff_map%boundary(2,2)=maxval(diff_map%local_map(2,:))+norder
     diff_map%boundary(1,3)=sphere%map3d(3,parallel%mygrid_range(1))-norder
     diff_map%boundary(2,3)=sphere%map3d(3,parallel%mygrid_range(2))+norder
   ENDSUBROUTINE grid_sphere_init

   SUBROUTINE set_wrap_grid_iso(myrho,wrap_box)
     IMPLICIT NONE
     REAL(DP) :: wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1),&
           & diff_map%boundary(1,2):diff_map%boundary(2,2),&
           & diff_map%boundary(1,3):diff_map%boundary(2,3))
     REAL(DP) :: myrho(:)
     !> local
     INTEGER(I4B) :: i,j,len
     REAL(DP),allocatable :: rho_send(:)
     REAL(DP),allocatable :: rho_recv(:)
     integer(I4B),save                 :: tag1=100
     integer(I4B)                 :: tag2
     integer(I4B),allocatable     :: rq(:),req(:)
     integer(I4B),allocatable     :: recflag(:)
     INTEGER(I4B) :: status1(MPI_STATUS_SIZE,diff_map%mycomm_cores(2):parallel%myid+1-1)
     INTEGER(I4B) :: status2(MPI_STATUS_SIZE,parallel%myid+2:diff_map%mycomm_cores(1))

     !> data communcation
     tag1=tag1+1
     tag2=200
     if(tag1>tag2)tag1=100
     allocate(rq(diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     do i=diff_map%mycomm_cores(2),diff_map%mycomm_cores(1),1
        j=i-1
        if( parallel%myid == j )cycle
        len=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        ! allocate(rho_send(len))
        ! rho_send=myrho(diff_map%mysend_size(1,i):diff_map%mysend_size(2,i))
        call MPI_ISEND(myrho(diff_map%mysend_size(1,i)),len,MPI_REAL8,j,&
             & tag1,parallel%comm,rq(i),mpinfo)
        ! deallocate(rho_send)
     enddo
     allocate(rho_recv(diff_map%mycomm_size(1,diff_map%mycomm_cores(2)):diff_map%mycomm_size(2,diff_map%mycomm_cores(1))))
     allocate(req(diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     allocate(recflag(diff_map%mycomm_cores(2):diff_map%mycomm_cores(1)))
     do i=diff_map%mycomm_cores(1),diff_map%mycomm_cores(2),-1
        j=i-1
        if( parallel%myid == j )cycle
        len=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        call MPI_IRECV(rho_recv(diff_map%mycomm_size(1,i)),len,MPI_REAL8,j,&
             & tag1,parallel%comm,req(i),mpinfo)
     enddo

     !> deal with the local part
     wrap_box=0
     rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

     if(parallel%myid+1-diff_map%mycomm_cores(2) /= 0)then
        call MPI_WAITALL(parallel%myid+1-diff_map%mycomm_cores(2),&
             & rq(diff_map%mycomm_cores(2):parallel%myid+1-1),status1,mpinfo)
        call MPI_WAITALL(parallel%myid+1-diff_map%mycomm_cores(2),&
             & req(diff_map%mycomm_cores(2):parallel%myid+1-1),status1,mpinfo)
     endif

     if(diff_map%mycomm_cores(1)-parallel%myid-1 /= 0)then
        call MPI_WAITALL(diff_map%mycomm_cores(1)-parallel%myid-1,&
             & rq(parallel%myid+2:diff_map%mycomm_cores(1)),status2,mpinfo)
        call MPI_WAITALL(diff_map%mycomm_cores(1)-parallel%myid-1,&
             & req(parallel%myid+2:diff_map%mycomm_cores(1)),status2,mpinfo)
     endif
     ! !> deal with the local part
     ! wrap_box=0
     ! rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     ! do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
     !    wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     ! enddo

     !> sphere segment change to wrap box
     do i=diff_map%mycomm_size(1,diff_map%mycomm_cores(2)),parallel%mygrid_range(1)-1,1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

     do i=parallel%mygrid_range(2)+1,diff_map%mycomm_size(2,diff_map%mycomm_cores(1)),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

   ENDSUBROUTINE set_wrap_grid_iso

   SUBROUTINE destroy_diff_map
     USE m_time_evaluate, ONLY: memory_free
     IMPLICIT NONE
     diff_map%mycomm_cores=0
     diff_map%boundary=0
     if(allocated(diff_map%nz_map))then
        call memory_free('diff_map_nz_map',real(size(diff_map%nz_map),DP)*I4B)
        deallocate(diff_map%nz_map)
     endif
     if(allocated(diff_map%mycomm_size))then
        call memory_free('diff_map_mycomm_size',real(size(diff_map%mycomm_size),DP)*I4B)
        deallocate(diff_map%mycomm_size)
     endif
     if(allocated(diff_map%mysend_size))then
        call memory_free('diff_map_mysend_size',real(size(diff_map%mysend_size),DP)*I4B)
        deallocate(diff_map%mysend_size)
     endif
     if(allocated(diff_map%local_map1d))then
        call memory_free('diff_map_local_map1d',real(size(diff_map%local_map1d),DP)*I4B)
        deallocate(diff_map%local_map1d)
     endif
     if(allocated(diff_map%local_map))then
        call memory_free('diff_map_local_map',real(size(diff_map%local_map),DP)*I4B)
        deallocate(diff_map%local_map)
     endif
     if(allocated(sphere%map3d))then
        call memory_free('diff_map_sphere_map3d',real(size(sphere%map3d),DP)*I4B)
        deallocate(sphere%map3d)
     endif
   ENDSUBROUTINE destroy_diff_map

   ! SUBROUTINE grid_cubic_init(n1,n2,n3,n,norder)
   !   IMPLICIT NONE
   !   INTEGER(I4B) :: n1,n2,n3,n
   !   ! type(sphere_type) :: sphere
   !   INTEGER(I4B) :: norder
   !   !> local
   !   INTEGER(I4B) :: i,j,k,in,id
   !   INTEGER(I4B) :: local_z,up_z,down_z,z_translated
   !   INTEGER(I4B) :: comm_flag,sum_flag

   !   ! i=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
   !   ! j=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
   !   ! k=(parallel%mygrid_range(1)-1)/n1/n2+1
   !   ! print *,'rank',parallel%rankx,'n1,n2,n3',n1,n2,n3,'grid_range(a,b,lab)'&
   !   !      &, parallel%mygrid_range,'grid(ijk)',i,j,k
   !   ! i=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
   !   ! j=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
   !   ! k=(parallel%mygrid_range(2)-1)/n1/n2+1
   !   ! print *,'rank',parallel%rankx,'grid(ijk)',i,j,k

   !   !> set the difference communcation list
   !   !> find up cores
   !   !> local down z: (parallel%mygrid_range(1)-1)/n1/n2+1
   !   !> local up z: (parallel%mygrid_range(2)-1)/n1/n2+1
   !   local_z=(parallel%mygrid_range(2)-1)/n1/n2+1
   !   do i=parallel%rankx+1,parallel%dims(2),1
   !      up_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1
   !      if( local_z+norder >= up_z )then
   !         !> different with isolated code, there is counts not limit
   !         diff_map%mycomm_cores(1)=diff_map%mycomm_cores(1)+1
   !      endif
   !   enddo
   !   z_translated=local_z+norder
   !   do while ( z_translated >= up_z )
   !      z_translated=z_translated-n3
   !      do i=1,parallel%dims(2),1
   !         up_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1
   !         if( z_translated >= up_z )then
   !            diff_map%mycomm_cores(1)=diff_map%mycomm_cores(1)+1
   !         endif
   !      enddo
   !   enddo

   !   !> find dowm cores
   !   local_z=(parallel%mygrid_range(1)-1)/n1/n2+1
   !   do i=parallel%rankx+1,1,-1
   !      down_z=(parallel%global_gridrange(2,i)-1)/n1/n2+1
   !      if( local_z-norder <= down_z )then
   !         diff_map%mycomm_cores(2)=diff_map%mycomm_cores(2)+1
   !      endif
   !   enddo
   !   z_translated=local_z-norder
   !   do while ( z_translated <= down_z )
   !      z_translated=z_translated+n3
   !      do i=parallel%dims(2),1,-1
   !         down_z=(parallel%global_gridrange(2,i)-1)/n1/n2+1
   !         if( z_translated <= down_z )then
   !            diff_map%mycomm_cores(2)=diff_map%mycomm_cores(2)+1
   !         endif
   !      enddo
   !   enddo
   !   ! if(parallel%isroot)write(*,'(A6,X,100(I8,2X))')'nz_map',mydiff_map%nz_map
   !   ! write(*,'(A6,X,100(I8,2X))')'_cores',mydiff_map%mycomm_cores,parallel%myid
   !   ! CALL MPI_BARRIER(parallel%comm,mpinfo)
   !   ! stop
   !   ! print *,'diff_map%mycomm_cores',diff_map%mycomm_cores,'rank',parallel%rankx

   !   !> set the size for receive
   !   !> local down z: (parallel%mygrid_range(1)-1)/n1/n2+1
   !   !> local up z: (parallel%mygrid_range(2)-1)/n1/n2+1
   !   down_z=(parallel%mygrid_range(1)-1)/n1/n2+1-norder
   !   do while(down_z<=0)
   !      down_z=down_z+n3
   !   enddo
   !   up_z=(parallel%mygrid_range(2)-1)/n1/n2+1+norder
   !   do while(up_z>n3)
   !      up_z=up_z-n3
   !   enddo
   !   allocate(diff_map%mycomm_size(2,&
   !        &parallel%rankx+2-diff_map%mycomm_cores(2)&
   !        & :parallel%rankx+diff_map%mycomm_cores(1)))
   !   ! print *,'shape of mycomm_size',shape(diff_map%mycomm_size),'rank',parallel%rankx,'lower and upper boundaries',parallel%rankx+2-diff_map%mycomm_cores(2),parallel%rankx+diff_map%mycomm_cores(1)
   !   do id=parallel%rankx+1,&
   !        &parallel%rankx+diff_map%mycomm_cores(1),1
   !      i=id
   !      do while(i<=0)
   !         i=i+parallel%dims(2)
   !      enddo
   !      do while(i>parallel%dims(2))
   !         i=i-parallel%dims(2)
   !      enddo
   !      if(id == parallel%rankx+1)then
   !         diff_map%mycomm_size(1,id)=parallel%mygrid_range(1)
   !         diff_map%mycomm_size(2,id)=parallel%mygrid_range(2)
   !      elseif( id == parallel%rankx+diff_map%mycomm_cores(1) )then
   !         if(parallel%global_gridrange(2,i) > diff_map%nz_map(2,up_z))then
   !            diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id-1)+1
   !            diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id)+&
   !                 & diff_map%nz_map(2,up_z)-parallel%global_gridrange(1,i)
   !         else
   !            diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id-1)+1
   !            diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id)+parallel%global_gridrange(3,i)-1
   !         endif
   !      else
   !         diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id-1)+1
   !         diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id)+parallel%global_gridrange(3,i)-1
   !      endif
   !   enddo
   !   do id=parallel%rankx+1,&
   !        &parallel%rankx+2-diff_map%mycomm_cores(2),-1
   !      i=id
   !      do while(i<=0)
   !         i=i+parallel%dims(2)
   !      enddo
   !      do while(i>parallel%dims(2))
   !         i=i-parallel%dims(2)
   !      enddo
   !      if(id == parallel%rankx+1)then
   !         diff_map%mycomm_size(1,id)=parallel%mygrid_range(1)
   !         diff_map%mycomm_size(2,id)=parallel%mygrid_range(2)
   !      elseif( id == parallel%rankx+2-diff_map%mycomm_cores(2) )then
   !         if(parallel%global_gridrange(1,i) < diff_map%nz_map(1,down_z))then
   !            diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id+1)-1
   !            diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id)&
   !                 & -parallel%global_gridrange(2,i)+diff_map%nz_map(1,down_z)
   !         else
   !            diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id+1)-1
   !            diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id)-parallel%global_gridrange(3,i)+1
   !         endif
   !      else
   !         diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id+1)-1
   !         diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id)-parallel%global_gridrange(3,i)+1
   !      endif
   !   enddo

   !   !> set the size for send
   !   allocate(diff_map%mysend_size(2,&
   !        &parallel%rankx+2-diff_map%mycomm_cores(2)&
   !        & :parallel%rankx+diff_map%mycomm_cores(1)))
   !   do id=parallel%rankx+2-diff_map%mycomm_cores(2),&
   !        &parallel%rankx,1
   !      i=id
   !      do while(i<=0)
   !         i=i+parallel%dims(2)
   !      enddo
   !      do while(i>parallel%dims(2))
   !         i=i-parallel%dims(2)
   !      enddo
   !      !> judge if i was the border
   !      up_z=(parallel%global_gridrange(2,i)-1)/n1/n2+1+norder
   !      do while(up_z>n3)
   !         up_z=up_z-n3
   !      enddo
   !      down_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1-norder
   !      do while(down_z<=0)
   !         down_z=down_z+n3
   !      enddo
   !      if( parallel%mygrid_range(1) < diff_map%nz_map(2,up_z) .and. &
   !           & diff_map%nz_map(2,up_z) < parallel%mygrid_range(2))then
   !         !> i am upper bound
   !         diff_map%mysend_size(1,id)=1!parallel%mygrid_range(1)
   !         diff_map%mysend_size(2,id)=diff_map%nz_map(2,up_z)-parallel%mygrid_range(1)+1
   !      else
   !         diff_map%mysend_size(1,id)=1!parallel%mygrid_range(1)
   !         diff_map%mysend_size(2,id)=parallel%mygrid_range(3)
   !      endif
   !      ! print*,'rank',parallel%rankx,'send_size',diff_map%mysend_size(:,id)
   !   enddo
   !   do id=parallel%rankx+1,&
   !        &parallel%rankx+diff_map%mycomm_cores(1),1
   !      i=id
   !      do while(i<=0)
   !         i=i+parallel%dims(2)
   !      enddo
   !      do while(i>parallel%dims(2))
   !         i=i-parallel%dims(2)
   !      enddo
   !      !> judge if i was the border
   !      down_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1-norder
   !      do while(down_z<=0)
   !         down_z=down_z+n3
   !      enddo
   !      if(diff_map%nz_map(1,down_z)<parallel%mygrid_range(2)&
   !           &.and.diff_map%nz_map(1,down_z) > parallel%mygrid_range(1))then
   !         diff_map%mysend_size(1,id)=diff_map%nz_map(1,down_z)-parallel%mygrid_range(1)+1
   !         diff_map%mysend_size(2,id)=parallel%mygrid_range(2)-parallel%mygrid_range(1)+1
   !      else
   !         diff_map%mysend_size(1,id)=1!parallel%mygrid_range(1)
   !         diff_map%mysend_size(2,id)=parallel%mygrid_range(3)
   !      endif
   !      ! print*,'rank',parallel%rankx,'send_size',diff_map%mysend_size(:,id)
   !   enddo

   !   !> save the 3d position for wrapped box
   !   ! print *,'rank',parallel%rankx,'start cores',parallel%rankx-diff_map%mycomm_cores(2)+1,'start_z',(parallel%mygrid_range(1)-1)/n1/n2+1-norder
   !   ! print *,'rank',parallel%rankx,'end core',parallel%rankx+diff_map%mycomm_cores(1)-1,'end_z',(parallel%mygrid_range(2)-1)/n1/n2+1+norder
   !   down_z=(parallel%mygrid_range(1)-1)/n1/n2+1-norder
   !   up_z=(parallel%mygrid_range(2)-1)/n1/n2+1+norder
   !   if(sum(diff_map%mycomm_cores)>parallel%numprocs)then
   !      ! if(parallel%isroot)print*,'communicate all grids'
   !      comm_flag=1
   !      ! diff_map%len_comm=diff_map%nz_map(2,up_z)-diff_map%nz_map(1,down_z)
   !      allocate(diff_map%local_map(4,n))
   !      in=0
   !      do k=1,n1,1
   !         do j=1,n2,1
   !            do i=1,n1,1
   !               in=in+1
   !               diff_map%local_map(1,in)=i
   !               diff_map%local_map(2,in)=j
   !               diff_map%local_map(3,in)=k
   !               diff_map%local_map(4,in)=in
   !            enddo
   !         enddo
   !      enddo
   !      ! i=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
   !      ! j=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
   !      ! k=(parallel%mygrid_range(2)-1)/n1/n2+1
   !      !sphere%map3d(:,diff_map%nz_map(1,down_z):diff_map%nz_map(2,up_z))

   !      !> get the wrap-grid boundary
   !      diff_map%boundary(1,1)=1-norder
   !      diff_map%boundary(2,1)=n1+norder
   !      diff_map%boundary(1,2)=1-norder
   !      diff_map%boundary(2,2)=n2+norder
   !      diff_map%boundary(1,3)=down_z
   !      diff_map%boundary(2,3)=up_z
   !      ! print*,'boundary',diff_map%boundary
   !      ! if(diff_map%boundary(1,3)==1)diff_map%boundary(1,3)=-norder+1
   !      ! if(diff_map%boundary(2,3)==sphere%map3d(3,sphere%Length))diff_map%boundary(2,3)=&
   !      ! & diff_map%boundary(2,3)+norder
   !   else
   !      ! if(parallel%isroot)print*,'communicate partial grids'
   !      comm_flag=2
   !      allocate(diff_map%local_map(3,(down_z-1)*n1*n2+1:up_z*n1*n2))
   !      in=(down_z-1)*n1*n2
   !      do k=down_z,up_z,1
   !         do j=1,n2,1
   !            do i=1,n1,1
   !               in=in+1
   !               diff_map%local_map(1,in)=i
   !               diff_map%local_map(2,in)=j
   !               diff_map%local_map(3,in)=k
   !               ! diff_map%local_map(4,in)=in
   !            enddo
   !         enddo
   !      enddo
   !      allocate(diff_map%local_map1d(1-norder:n1+norder&
   !           &,1-norder:n2+norder,down_z:up_z))
   !      in=0!(down_z-1)*n1*n2
   !      do k=down_z,up_z,1
   !         do j=1,n2,1
   !            do i=1,n1,1
   !               in=in+1
   !               diff_map%local_map1d(i,j,k)=in
   !            enddo
   !         enddo
   !      enddo
   !      ! print*,'max_local_map1d',in
   !      diff_map%local_map1d(n1+1:n1+norder,1:n2,:)=&
   !           & diff_map%local_map1d(1:norder,1:n2,:)
   !      diff_map%local_map1d(1-norder:0,1:n2,:)=&
   !           & diff_map%local_map1d(n1-norder+1:n1,1:n2,:)
   !      diff_map%local_map1d(:,n2+1:n2+norder,:)=&
   !           & diff_map%local_map1d(:,1:norder,:)
   !      diff_map%local_map1d(:,1-norder:0,:)=&
   !           & diff_map%local_map1d(:,n2-norder+1:n2,:)

   !      !> get the wrap-grid boundary
   !      ! diff_map%boundary(1,1)=1
   !      ! diff_map%boundary(2,1)=n1
   !      ! diff_map%boundary(1,2)=1
   !      ! diff_map%boundary(2,2)=n2
   !      diff_map%boundary(1,1)=1-norder
   !      diff_map%boundary(2,1)=n1+norder
   !      diff_map%boundary(1,2)=1-norder
   !      diff_map%boundary(2,2)=n2+norder
   !      diff_map%boundary(1,3)=down_z
   !      diff_map%boundary(2,3)=up_z
   !   endif
   !   ! CALL MPI_BCAST(Lall_grid ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
   !   CALL MPI_ALLREDUCE(comm_flag, sum_flag, 1, MPI_INTEGER4, MPI_SUM, parallel%commx, mpinfo)
   !   if(sum_flag==parallel%numprocs*2)then
   !      Lall_grid=.false.
   !      if(parallel%isroot)print*,'communicate partial grids'
   !   else
   !      Lall_grid=.true.
   !      if(parallel%isroot)print*,'communicate all grids'
   !   endif
   !   ! Lall_grid=.true.
   ! ENDSUBROUTINE grid_cubic_init
   SUBROUTINE grid_cubic_init(n1,n2,n3,n,norder)
     IMPLICIT NONE
     INTEGER(I4B) :: n1,n2,n3,n
     ! type(sphere_type) :: sphere
     INTEGER(I4B) :: norder
     !> local
     INTEGER(I4B) :: i,j,k,in,id
     INTEGER(I4B) :: local_z,up_z,down_z,z_translated
     INTEGER(I4B) :: comm_flag,sum_flag

     ! i=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! j=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! k=(parallel%mygrid_range(1)-1)/n1/n2+1
     ! print *,'rank',parallel%rankx,'n1,n2,n3',n1,n2,n3,'grid_range(a,b,lab)'&
     !      &, parallel%mygrid_range,'grid(ijk)',i,j,k
     ! i=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
     ! j=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
     ! k=(parallel%mygrid_range(2)-1)/n1/n2+1
     ! print *,'rank',parallel%rankx,'grid(ijk)',i,j,k

     !> set the difference communcation list
     !> find up cores
     !> local down z: (parallel%mygrid_range(1)-1)/n1/n2+1
     !> local up z: (parallel%mygrid_range(2)-1)/n1/n2+1
     local_z=(parallel%mygrid_range(2)-1)/n1/n2+1
     do i=parallel%rankx+1,parallel%dims(2),1
        up_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1
        if( local_z+norder >= up_z )then
           !> different with isolated code, there is counts not limit
           diff_map%mycomm_cores(1)=diff_map%mycomm_cores(1)+1
        endif
     enddo
     z_translated=local_z+norder
     do while ( z_translated >= up_z )
        z_translated=z_translated-n3
        do i=1,parallel%dims(2),1
           up_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1
           if( z_translated >= up_z )then
              diff_map%mycomm_cores(1)=diff_map%mycomm_cores(1)+1
           endif
        enddo
     enddo

     !> find dowm cores
     local_z=(parallel%mygrid_range(1)-1)/n1/n2+1
     do i=parallel%rankx+1,1,-1
        down_z=(parallel%global_gridrange(2,i)-1)/n1/n2+1
        if( local_z-norder <= down_z )then
           diff_map%mycomm_cores(2)=diff_map%mycomm_cores(2)+1
        endif
     enddo
     z_translated=local_z-norder
     do while ( z_translated <= down_z )
        z_translated=z_translated+n3
        do i=parallel%dims(2),1,-1
           down_z=(parallel%global_gridrange(2,i)-1)/n1/n2+1
           if( z_translated <= down_z )then
              diff_map%mycomm_cores(2)=diff_map%mycomm_cores(2)+1
           endif
        enddo
     enddo
     ! if(parallel%isroot)write(*,'(A6,X,100(I8,2X))')'nz_map',mydiff_map%nz_map
     ! write(*,'(A6,X,100(I8,2X))')'_cores',mydiff_map%mycomm_cores,parallel%myid
     ! CALL MPI_BARRIER(parallel%comm,mpinfo)
     ! stop
     ! print *,'diff_map%mycomm_cores',diff_map%mycomm_cores,'rank',parallel%rankx

     !> set the size for receive
     !> local down z: (parallel%mygrid_range(1)-1)/n1/n2+1
     !> local up z: (parallel%mygrid_range(2)-1)/n1/n2+1
     down_z=(parallel%mygrid_range(1)-1)/n1/n2+1-norder
     do while(down_z<=0)
        down_z=down_z+n3
     enddo
     up_z=(parallel%mygrid_range(2)-1)/n1/n2+1+norder
     do while(up_z>n3)
        up_z=up_z-n3
     enddo
     allocate(diff_map%mycomm_size(2,&
          &parallel%rankx+2-diff_map%mycomm_cores(2)&
          & :parallel%rankx+diff_map%mycomm_cores(1)))
     ! print *,'shape of mycomm_size',shape(diff_map%mycomm_size),'rank',parallel%rankx,'lower and upper boundaries',parallel%rankx+2-diff_map%mycomm_cores(2),parallel%rankx+diff_map%mycomm_cores(1)
     do id=parallel%rankx+1,&
          &parallel%rankx+diff_map%mycomm_cores(1),1
        i=id
        do while(i<=0)
           i=i+parallel%dims(2)
        enddo
        do while(i>parallel%dims(2))
           i=i-parallel%dims(2)
        enddo
        if(id == parallel%rankx+1)then
           diff_map%mycomm_size(1,id)=parallel%mygrid_range(1)
           diff_map%mycomm_size(2,id)=parallel%mygrid_range(2)
        elseif( id == parallel%rankx+diff_map%mycomm_cores(1) )then
           if(parallel%global_gridrange(2,i) > diff_map%nz_map(2,up_z))then
              diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id-1)+1
              diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id)+&
                   & diff_map%nz_map(2,up_z)-parallel%global_gridrange(1,i)
           else
              diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id-1)+1
              diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id)+parallel%global_gridrange(3,i)-1
           endif
        else
           diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id-1)+1
           diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id)+parallel%global_gridrange(3,i)-1
        endif
     enddo
     do id=parallel%rankx+1,&
          &parallel%rankx+2-diff_map%mycomm_cores(2),-1
        i=id
        do while(i<=0)
           i=i+parallel%dims(2)
        enddo
        do while(i>parallel%dims(2))
           i=i-parallel%dims(2)
        enddo
        if(id == parallel%rankx+1)then
           diff_map%mycomm_size(1,id)=parallel%mygrid_range(1)
           diff_map%mycomm_size(2,id)=parallel%mygrid_range(2)
        elseif( id == parallel%rankx+2-diff_map%mycomm_cores(2) )then
           if(parallel%global_gridrange(1,i) < diff_map%nz_map(1,down_z))then
              diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id+1)-1
              diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id)&
                   & -parallel%global_gridrange(2,i)+diff_map%nz_map(1,down_z)
           else
              diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id+1)-1
              diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id)-parallel%global_gridrange(3,i)+1
           endif
        else
           diff_map%mycomm_size(2,id)=diff_map%mycomm_size(1,id+1)-1
           diff_map%mycomm_size(1,id)=diff_map%mycomm_size(2,id)-parallel%global_gridrange(3,i)+1
        endif
     enddo

     !> set the size for send
     allocate(diff_map%mysend_size(2,&
          &parallel%rankx+2-diff_map%mycomm_cores(2)&
          & :parallel%rankx+diff_map%mycomm_cores(1)))
     do id=parallel%rankx+2-diff_map%mycomm_cores(2),&
          &parallel%rankx,1
        i=id
        do while(i<=0)
           i=i+parallel%dims(2)
        enddo
        do while(i>parallel%dims(2))
           i=i-parallel%dims(2)
        enddo
        !> judge if i was the border
        up_z=(parallel%global_gridrange(2,i)-1)/n1/n2+1+norder
        do while(up_z>n3)
           up_z=up_z-n3
        enddo
        down_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1-norder
        do while(down_z<=0)
           down_z=down_z+n3
        enddo
        if( parallel%mygrid_range(1) < diff_map%nz_map(2,up_z) .and. &
             & diff_map%nz_map(2,up_z) < parallel%mygrid_range(2))then
           !> i am upper bound
           diff_map%mysend_size(1,id)=1!parallel%mygrid_range(1)
           diff_map%mysend_size(2,id)=diff_map%nz_map(2,up_z)-parallel%mygrid_range(1)+1
        else
           diff_map%mysend_size(1,id)=1!parallel%mygrid_range(1)
           diff_map%mysend_size(2,id)=parallel%mygrid_range(3)
        endif
        ! print*,'rank',parallel%rankx,'send_size',diff_map%mysend_size(:,id)
     enddo
     do id=parallel%rankx+1,&
          &parallel%rankx+diff_map%mycomm_cores(1),1
        i=id
        do while(i<=0)
           i=i+parallel%dims(2)
        enddo
        do while(i>parallel%dims(2))
           i=i-parallel%dims(2)
        enddo
        !> judge if i was the border
        down_z=(parallel%global_gridrange(1,i)-1)/n1/n2+1-norder
        do while(down_z<=0)
           down_z=down_z+n3
        enddo
        if(diff_map%nz_map(1,down_z)<parallel%mygrid_range(2)&
             &.and.diff_map%nz_map(1,down_z) > parallel%mygrid_range(1))then
           diff_map%mysend_size(1,id)=diff_map%nz_map(1,down_z)-parallel%mygrid_range(1)+1
           diff_map%mysend_size(2,id)=parallel%mygrid_range(2)-parallel%mygrid_range(1)+1
        else
           diff_map%mysend_size(1,id)=1!parallel%mygrid_range(1)
           diff_map%mysend_size(2,id)=parallel%mygrid_range(3)
        endif
        ! print*,'rank',parallel%rankx,'send_size',diff_map%mysend_size(:,id)
     enddo

     !> save the 3d position for wrapped box
     ! print *,'rank',parallel%rankx,'start cores',parallel%rankx-diff_map%mycomm_cores(2)+1,'start_z',(parallel%mygrid_range(1)-1)/n1/n2+1-norder
     ! print *,'rank',parallel%rankx,'end core',parallel%rankx+diff_map%mycomm_cores(1)-1,'end_z',(parallel%mygrid_range(2)-1)/n1/n2+1+norder
     down_z=(parallel%mygrid_range(1)-1)/n1/n2+1-norder
     up_z=(parallel%mygrid_range(2)-1)/n1/n2+1+norder
     if(sum(diff_map%mycomm_cores)>parallel%numprocs)then
        comm_flag=1
     else
        comm_flag=2
     endif
     CALL MPI_ALLREDUCE(comm_flag, sum_flag, 1, MPI_INTEGER4, MPI_SUM, parallel%commx, mpinfo)
     if(sum_flag==parallel%numprocs*2)then
        Lall_grid=.false.
        if(parallel%isroot)print*,'communicate partial grids'
     else
        Lall_grid=.true.
        if(parallel%isroot)print*,'communicate all grids'
     endif

     if(Lall_grid)then
        allocate(diff_map%local_map(4,n))
        in=0
        do k=1,n3,1
           do j=1,n2,1
              do i=1,n1,1
                 in=in+1
                 diff_map%local_map(1,in)=i
                 diff_map%local_map(2,in)=j
                 diff_map%local_map(3,in)=k
                 diff_map%local_map(4,in)=in
              enddo
           enddo
        enddo
        allocate(diff_map%local_map1d(1-norder:n1+norder&
             &,1-norder:n2+norder,1-norder:n3+norder))
        in=0!(down_z-1)*n1*n2
        do k=1,n3,1
           do j=1,n2,1
              do i=1,n1,1
                 in=in+1
                 diff_map%local_map1d(i,j,k)=in
              enddo
           enddo
        enddo
        !> x
        diff_map%local_map1d(n1+1:n1+norder,1:n2,:)=&
             & diff_map%local_map1d(1:norder,1:n2,:)
        diff_map%local_map1d(1-norder:0,1:n2,:)=&
             & diff_map%local_map1d(n1-norder+1:n1,1:n2,:)
        !> y
        diff_map%local_map1d(:,n2+1:n2+norder,:)=&
             & diff_map%local_map1d(:,1:norder,:)
        diff_map%local_map1d(:,1-norder:0,:)=&
             & diff_map%local_map1d(:,n2-norder+1:n2,:)
        !> z
        diff_map%local_map1d(:,:,n3+1:n3+norder)=&
             & diff_map%local_map1d(:,:,1:norder)
        diff_map%local_map1d(:,:,1-norder:0)=&
             & diff_map%local_map1d(:,:,n3-norder+1:n3)
        ! i=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
        ! j=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
        ! k=(parallel%mygrid_range(2)-1)/n1/n2+1
        !sphere%map3d(:,diff_map%nz_map(1,down_z):diff_map%nz_map(2,up_z))

        !> get the wrap-grid boundary
        diff_map%boundary(1,1)=1-norder
        diff_map%boundary(2,1)=n1+norder
        diff_map%boundary(1,2)=1-norder
        diff_map%boundary(2,2)=n2+norder
        diff_map%boundary(1,3)=down_z
        diff_map%boundary(2,3)=up_z
        ! print*,'boundary',diff_map%boundary
        ! if(diff_map%boundary(1,3)==1)diff_map%boundary(1,3)=-norder+1
        ! if(diff_map%boundary(2,3)==sphere%map3d(3,sphere%Length))diff_map%boundary(2,3)=&
        ! & diff_map%boundary(2,3)+norder
     else
        ! if(parallel%isroot)print*,'communicate partial grids'
        comm_flag=2
        allocate(diff_map%local_map(3,(down_z-1)*n1*n2+1:up_z*n1*n2))
        in=(down_z-1)*n1*n2
        do k=down_z,up_z,1
           do j=1,n2,1
              do i=1,n1,1
                 in=in+1
                 diff_map%local_map(1,in)=i
                 diff_map%local_map(2,in)=j
                 diff_map%local_map(3,in)=k
                 ! diff_map%local_map(4,in)=in
              enddo
           enddo
        enddo
        allocate(diff_map%local_map1d(1-norder:n1+norder&
             &,1-norder:n2+norder,down_z:up_z))
        in=0!(down_z-1)*n1*n2
        do k=down_z,up_z,1
           do j=1,n2,1
              do i=1,n1,1
                 in=in+1
                 diff_map%local_map1d(i,j,k)=in
              enddo
           enddo
        enddo
        ! print*,'max_local_map1d',in
        diff_map%local_map1d(n1+1:n1+norder,1:n2,:)=&
             & diff_map%local_map1d(1:norder,1:n2,:)
        diff_map%local_map1d(1-norder:0,1:n2,:)=&
             & diff_map%local_map1d(n1-norder+1:n1,1:n2,:)
        diff_map%local_map1d(:,n2+1:n2+norder,:)=&
             & diff_map%local_map1d(:,1:norder,:)
        diff_map%local_map1d(:,1-norder:0,:)=&
             & diff_map%local_map1d(:,n2-norder+1:n2,:)

        !> get the wrap-grid boundary
        ! diff_map%boundary(1,1)=1
        ! diff_map%boundary(2,1)=n1
        ! diff_map%boundary(1,2)=1
        ! diff_map%boundary(2,2)=n2
        diff_map%boundary(1,1)=1-norder
        diff_map%boundary(2,1)=n1+norder
        diff_map%boundary(1,2)=1-norder
        diff_map%boundary(2,2)=n2+norder
        diff_map%boundary(1,3)=down_z
        diff_map%boundary(2,3)=up_z
     endif
     ! CALL MPI_BCAST(Lall_grid ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
     ! Lall_grid=.true.
   ENDSUBROUTINE grid_cubic_init

   SUBROUTINE set_wrap_grid_per(myrho,wrap_box,global_n,global_n1&
        &, global_n2)
     USE parameters, ONLY: norder=>finite_order
     IMPLICIT NONE
     INTEGER(I4B) :: global_n,global_n1,global_n2
     COMPLEX(DP) :: wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1),&
          & diff_map%boundary(1,2):diff_map%boundary(2,2),&
          & diff_map%boundary(1,3):diff_map%boundary(2,3))
     COMPLEX(DP) :: myrho(:)
     !> local
     INTEGER(I4B) :: i,j,len
     COMPLEX(DP),allocatable :: rho_recv(:)
     integer(I4B),allocatable     :: rq(:),req(:)
     integer(I4B),allocatable     :: recflag(:)
     INTEGER(I4B) :: status(MPI_STATUS_SIZE,-10:10)
     INTEGER(I4B) :: core_start,core_end
     integer(I4B) :: id,icount1,icount2
     ! integer(I4B) :: tag(parallel%dims(2),parallel%dims(2)),tag1
     integer(I4B),save :: tag_i=100
     integer(I4B) :: tag1,tag2
     logical :: flag

     !> set tag
     ! id=0
     ! do j=1,parallel%dims(2),1
     !    do i=1,j,1
     !       id=id+1
     !       tag(i,j)=id
     !       tag(j,i)=id
     !    enddo
     ! enddo
     tag_i=tag_i+1
     if(tag_i>200)tag_i=100

     core_start=parallel%rankx+2-diff_map%mycomm_cores(2)
     core_end=parallel%rankx+diff_map%mycomm_cores(1)

     !> receive
     allocate(rho_recv(diff_map%mycomm_size(1,core_start):diff_map%mycomm_size(2,core_end)))
     allocate(req(core_start:core_end))
     allocate(recflag(core_start:core_end))
     ! do i=diff_map%mycomm_cores(2),diff_map%mycomm_cores(1),1
     icount2=core_start
     do i=core_end,parallel%rankx+1,-1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        len=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        ! tag1=100+tag(parallel%rankx+1,j+1)*2
        tag1=tag_i*2
        call MPI_IRECV(rho_recv(diff_map%mycomm_size(1,i)),len,MPI_COMPLEX16,j,&
             & tag1,parallel%commx,req(icount2),mpinfo)
         ! print *,'rank',parallel%rankx,'start receive',j,'tag',tag1,'size',len
        ! call MPI_WAIT(req(icount2),status,mpinfo)
        ! call MPI_TEST(req(icount2),flag,status,mpinfo)
         ! print *,'rank',parallel%rankx,'end receive',j,'tag',tag1,'size',len,'flag',flag
        ! print*,'icount2',icount2,'rank',parallel%rankx
        icount2=icount2+1
        ! call MPI_WAIT(req(icount2),status,mpinfo)
     enddo
     !> data communcation send
     allocate(rq(core_start:core_end))
     icount1=core_start
     do i=core_start,parallel%rankx+1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        len=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        ! allocate(rho_send(len))
        ! rho_send=myrho(diff_map%mysend_size(1,i):diff_map%mysend_size(2,i))
        ! tag1=100+tag(parallel%rankx+1,j+1)*2
        tag1=tag_i*2
        ! print*,'rank',parallel%rankx,'to',j,'tag',tag1,'size',len
        call MPI_ISEND(myrho(diff_map%mysend_size(1,i)),len,MPI_COMPLEX16,j,&
             & tag1,parallel%commx,rq(icount1),mpinfo)
        icount1=icount1+1
        ! deallocate(rho_send)
     enddo
     call MPI_WAITALL(diff_map%mycomm_cores(1)-1,&
          &req(core_start:icount2-1),status,mpinfo)
     !> another wait
     ! if(icount2-core_start /= 0)then
     !    call MPI_WAITALL(diff_map%mycomm_cores(2)-1,&
     !         & req(core_start:icount2),status,mpinfo)
     ! endif

     ! if(diff_map%mycomm_cores(1)-parallel%rankx-1 /= 0)then
     !    call MPI_WAITALL(diff_map%mycomm_cores(1)-1,&
     !         & req(parallel%rankx+2:core_end),status,mpinfo)
     ! endif
     ! print*,'all finish',parallel%myid
     ! print*,'req size',size(req),icount2,core_start,core_end,'len',diff_map%mycomm_cores(1:2)
          ! &req(parallel%rankx+1:core_end-1),status,mpinfo)

     !> deal with the local part
     wrap_box=0
     rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo
     do i=parallel%mygrid_range(2),diff_map%mycomm_size(2,core_end),1
        id=i
        do while(id<=0)
           id=id+global_n
        enddo
        do while(id>global_n)
           id=id-global_n
        enddo
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

     do i=parallel%rankx+1,core_start,-1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        len=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        ! tag1=100+tag(parallel%rankx+1,j+1)*2+1
        tag1=tag_i*2+1
        ! print *,'rank',parallel%rankx,'start receive',j,'tag',tag1,'size',len,'rho_recv',shape(rho_recv(diff_map%mycomm_size(1,i):))
        call MPI_IRECV(rho_recv(diff_map%mycomm_size(1,i)),len,MPI_COMPLEX16,j,&
             & tag1,parallel%commx,req(icount2),mpinfo)
        ! call MPI_WAIT(req(icount2),status,mpinfo)
        ! call MPI_TEST(req(icount2),flag,status,mpinfo)
        ! print *,'rank',parallel%rankx,'end receive',j,'tag',tag1,'size',len,'flag',flag
        ! print*,'icount2',icount2,'rank',parallel%rankx
        icount2=icount2+1
        ! call MPI_WAIT(req(icount2),status,mpinfo)
     enddo

     do i=parallel%rankx+1,core_end
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        len=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        ! allocate(rho_send(len))
        ! rho_send=myrho(diff_map%mysend_size(1,i):diff_map%mysend_size(2,i))
        ! tag1=100+tag(parallel%rankx+1,j+1)*2+1
        tag1=tag_i*2+1
        ! print*,'rank',parallel%rankx,'to',j,'tag',tag1,'size',len
        call MPI_ISEND(myrho(diff_map%mysend_size(1,i)),len,MPI_COMPLEX16,j,&
             & tag1,parallel%commx,rq(icount1),mpinfo)
        icount1=icount1+1
        ! deallocate(rho_send)
     enddo



     call MPI_WAITALL(diff_map%mycomm_cores(2)-1,&
          &req(parallel%rankx+1:core_end-1),status,mpinfo)
          ! &req(core_start:parallel%rankx+1-1),status,mpinfo)
    ! CALL MPI_BARRIER(parallel%commx,mpinfo)
     !> sphere segment change to wrap box
     do i=diff_map%mycomm_size(1,core_start),parallel%mygrid_range(1),1
        id=i
        do while(id<=0)
           id=id+global_n
        enddo
        do while(id>global_n)
           id=id-global_n
        enddo
        wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=rho_recv(i)
     enddo

     ! print*,'shape',shape(wrap_box)

     !> set periodic grid
     wrap_box(global_n1+1:global_n1+norder,1:global_n2,:)=&
          & wrap_box(1:norder,1:global_n2,:)
     wrap_box(1-norder:0,1:global_n2,:)=&
          & wrap_box(global_n1-norder+1:global_n1,1:global_n2,:)
     wrap_box(1:global_n1,global_n2+1:global_n2+norder,:)=&
          & wrap_box(1:global_n1,1:norder,:)
     wrap_box(1:global_n1,1-norder:0,:)=&
          & wrap_box(1:global_n1,global_n2-norder+1:global_n2,:)
     wrap_box(global_n1+1:global_n1+norder&
          & ,global_n2+1:global_n2+norder,:)=&
          & wrap_box(1:norder,1:norder,:)
     wrap_box(1-norder:0,global_n2+1:global_n2+norder,:)=&
          & wrap_box(global_n1-norder+1:global_n1,1:norder,:)
     wrap_box(global_n1+1:global_n1+norder,1-norder:0,:)=&
          & wrap_box(1:norder,global_n2-norder+1:global_n2,:)
     wrap_box(1-norder:0,1-norder:0,:)=&
          & wrap_box(global_n1-norder+1:global_n1&
          & ,global_n2-norder+1:global_n2,:)

   ENDSUBROUTINE set_wrap_grid_per

   SUBROUTINE set_wrap_grid_per_ata(myrho,wrap_box1d,global_n,global_n1&
        &, global_n2)
     USE parameters, ONLY: norder=>finite_order
     IMPLICIT NONE
     INTEGER(I4B) :: global_n,global_n1,global_n2
     COMPLEX(DP) :: wrap_box1d(global_n1*global_n2*&
          & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1))
     ! COMPLEX(DP) :: wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1),&
     !      & diff_map%boundary(1,2):diff_map%boundary(2,2),&
     !      & diff_map%boundary(1,3):diff_map%boundary(2,3))
     COMPLEX(DP) :: myrho(:)
     !> local
     INTEGER(I4B) :: i,j,len
     COMPLEX(DP),allocatable :: rho_recv(:,:,:)
     integer(I4B),allocatable     :: rq(:),req(:)
     integer(I4B),allocatable     :: recflag(:)
     INTEGER(I4B) :: status(MPI_STATUS_SIZE,-10:10)
     INTEGER(I4B) :: core_start,core_end
     integer(I4B) :: id,icount1,icount2
     ! integer(I4B) :: tag(parallel%dims(2),parallel%dims(2)),tag1
     integer(I4B),save :: tag_i=100
     integer(I4B) :: tag1,tag2
     logical :: flag
     !> alltoallv setting
     INTEGER(I4B),dimension(parallel%dims(2)) :: scount,sdispls&
          &, rcount,rdispls
     integer(I4B) :: recv_i(3),i3_1,i3_2
     INTEGER(I4B) :: offset1d,offset1da
     COMPLEX(DP),allocatable :: temp(:,:,:)

     call start_time('comm_ata',.true.)
     !> init setting
     core_start=parallel%rankx+2-diff_map%mycomm_cores(2)
     core_end=parallel%rankx+diff_map%mycomm_cores(1)
     ! allocate(rho_recv(diff_map%mycomm_size(1,core_start):diff_map%mycomm_size(2,core_end)))
     ! allocate(rho_recv(1:global_n1,1:global_n2,diff_map%boundary(1,3):diff_map%boundary(2,3)))
     recv_i(1)=diff_map%local_map(1,diff_map%mycomm_size(1,core_start))
     recv_i(2)=diff_map%local_map(2,diff_map%mycomm_size(1,core_start))
     recv_i(3)=diff_map%local_map(3,diff_map%mycomm_size(1,core_start))
     ! print*,'there should be equal',recv_i(1)+global_n1*(recv_i(2)-1)+global_n1*global_n2*(recv_i(3)-1),diff_map%mycomm_size(1,core_start)
     ! print*,'two array size','rho_recv',size(rho_recv),'wrap_box1d',size(wrap_box1d)
     offset1d=(recv_i(3)-1)*global_n1*global_n2
     offset1da=recv_i(1)+(recv_i(2)-1)*global_n1

     !> send left and receive right
     rcount=0
     rdispls=0
     do i=core_end,core_start,-1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        rcount(j+1)=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        rdispls(j+1)=diff_map%mycomm_size(1,i)-diff_map%mycomm_size(1,core_start)
     enddo

     scount=0
     sdispls=0
     do i=core_start,core_end !parallel%rankx+1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        scount(j+1)=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        sdispls(j+1)=diff_map%mysend_size(1,i)-1
     enddo

     CALL MPI_ALLTOALLV(myrho,scount&
          & ,sdispls,MPI_COMPLEX16&
          & ,wrap_box1d(offset1da),rcount&
          & ,rdispls,MPI_COMPLEX16&
          & ,parallel%commx,mpinfo)
     call end_time('comm_ata',.true.)
     call start_time('comm_assign1',.true.)

!     !> send right and receive left
!     rcount=0
!     rdispls=0
!     do i=parallel%rankx+1,core_start,-1
!        j=i-1
!        do while (j>=parallel%dims(2))
!           j=j-parallel%dims(2)
!        enddo
!        do while (j<0)
!           j=j+parallel%dims(2)
!        enddo
!        if( parallel%rankx == j )cycle
!        rcount(j+1)=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
!        rdispls(j+1)=diff_map%mycomm_size(1,i)-diff_map%mycomm_size(1,core_start)+1
!     enddo
!
!     scount=0
!     sdispls=0
!     do i=parallel%rankx+1,core_end
!        j=i-1
!        do while (j>=parallel%dims(2))
!           j=j-parallel%dims(2)
!        enddo
!        do while (j<0)
!           j=j+parallel%dims(2)
!        enddo
!        if( parallel%rankx == j )cycle
!        scount(j+1)=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
!        sdispls(j+1)=diff_map%mysend_size(1,i)
!     enddo
!
!     ! print*,'rank',parallel%rankx,'|',scount,'|',sdispls,'|',shape(myrho)
!     ! print*,'rank',parallel%rankx,'|',rcount,'|',rdispls,'|',shape(rho_recv)
!     CALL MPI_ALLTOALLV(myrho,scount&
!          & ,sdispls,MPI_COMPLEX16&
!          & ,rho_recv,rcount&
!          & ,rdispls,MPI_COMPLEX16&
!          & ,parallel%commx,mpinfo)
     !> deal with the local part
     ! rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     j=0
     do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
        j=j+1
        ! wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=myrho(j)
        wrap_box1d(i-offset1d)=myrho(j)
     enddo

     call end_time('comm_assign1',.true.)
     call start_time('comm_assign2',.true.)
     ! !> copy
     ! do i=diff_map%mycomm_size(1,core_start),parallel%mygrid_range(1)-1
     !    wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=&
     !         & rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     !    wrap_box1d(i-offset1d)=rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     ! enddo
     ! do i=parallel%mygrid_range(2)+1,diff_map%mycomm_size(2,core_end),1
     !    wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=&
     !         & rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     !    wrap_box1d(i-offset1d)=rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     ! enddo

     call end_time('comm_assign2',.true.)
     call start_time('comm_assign3',.true.)
     ! print *,'lbound',lbound(wrap_box),'ubound',ubound(wrap_box)
     ! !> set periodic grid
     ! wrap_box(global_n1+1:global_n1+norder,1:global_n2,:)=&
     !      & wrap_box(1:norder,1:global_n2,:)
     ! wrap_box(1-norder:0,1:global_n2,:)=&
     !      & wrap_box(global_n1-norder+1:global_n1,1:global_n2,:)
     ! wrap_box(:,global_n2+1:global_n2+norder,:)=&
     !      & wrap_box(:,1:norder,:)
     ! wrap_box(:,1-norder:0,:)=&
     !      & wrap_box(:,global_n2-norder+1:global_n2,:)

     !> deallocate
     ! deallocate(rho_recv)
     call end_time('comm_assign3',.true.)
   ENDSUBROUTINE set_wrap_grid_per_ata
   SUBROUTINE set_wrap_grid_per_ata_real(myrho,wrap_box1d,global_n,global_n1&
        &, global_n2)
     USE parameters, ONLY: norder=>finite_order
     IMPLICIT NONE
     INTEGER(I4B) :: global_n,global_n1,global_n2
     REAL(DP) :: wrap_box1d(global_n1*global_n2*&
          & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1))
     ! COMPLEX(DP) :: wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1),&
     !      & diff_map%boundary(1,2):diff_map%boundary(2,2),&
     !      & diff_map%boundary(1,3):diff_map%boundary(2,3))
     REAL(DP) :: myrho(:)
     !> local
     INTEGER(I4B) :: i,j,len
     REAL(DP),allocatable :: rho_recv(:,:,:)
     integer(I4B),allocatable     :: rq(:),req(:)
     integer(I4B),allocatable     :: recflag(:)
     INTEGER(I4B) :: status(MPI_STATUS_SIZE,-10:10)
     INTEGER(I4B) :: core_start,core_end
     integer(I4B) :: id,icount1,icount2
     integer(I4B),save :: tag_i=100
     integer(I4B) :: tag1,tag2
     logical :: flag
     !> alltoallv setting
     INTEGER(I4B),dimension(parallel%dims(2)) :: scount,sdispls&
          &, rcount,rdispls
     integer(I4B) :: recv_i(3),i3_1,i3_2
     INTEGER(I4B) :: offset1d,offset1da
     REAL(DP),allocatable :: temp(:,:,:)

     call start_time('comm_ata',.true.)
     !> init setting
     core_start=parallel%rankx+2-diff_map%mycomm_cores(2)
     core_end=parallel%rankx+diff_map%mycomm_cores(1)
     ! allocate(rho_recv(diff_map%mycomm_size(1,core_start):diff_map%mycomm_size(2,core_end)))
     ! allocate(rho_recv(1:global_n1,1:global_n2,diff_map%boundary(1,3):diff_map%boundary(2,3)))
     recv_i(1)=diff_map%local_map(1,diff_map%mycomm_size(1,core_start))
     recv_i(2)=diff_map%local_map(2,diff_map%mycomm_size(1,core_start))
     recv_i(3)=diff_map%local_map(3,diff_map%mycomm_size(1,core_start))
     ! print*,'there should be equal',recv_i(1)+global_n1*(recv_i(2)-1)+global_n1*global_n2*(recv_i(3)-1),diff_map%mycomm_size(1,core_start)
     ! print*,'two array size','rho_recv',size(rho_recv),'wrap_box1d',size(wrap_box1d)
     offset1d=(recv_i(3)-1)*global_n1*global_n2
     offset1da=recv_i(1)+(recv_i(2)-1)*global_n1

     !> send left and receive right
     rcount=0
     rdispls=0
     do i=core_end,core_start,-1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        rcount(j+1)=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
        rdispls(j+1)=diff_map%mycomm_size(1,i)-diff_map%mycomm_size(1,core_start)
     enddo

     scount=0
     sdispls=0
     do i=core_start,core_end !parallel%rankx+1
        j=i-1
        do while (j>=parallel%dims(2))
           j=j-parallel%dims(2)
        enddo
        do while (j<0)
           j=j+parallel%dims(2)
        enddo
        if( parallel%rankx == j )cycle
        scount(j+1)=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
        sdispls(j+1)=diff_map%mysend_size(1,i)-1
     enddo

     CALL MPI_ALLTOALLV(myrho,scount&
          & ,sdispls,MPI_REAL8&
          & ,wrap_box1d(offset1da),rcount&
          & ,rdispls,MPI_REAL8&
          & ,parallel%commx,mpinfo)
     call end_time('comm_ata',.true.)
     call start_time('comm_assign1',.true.)

!     !> send right and receive left
!     rcount=0
!     rdispls=0
!     do i=parallel%rankx+1,core_start,-1
!        j=i-1
!        do while (j>=parallel%dims(2))
!           j=j-parallel%dims(2)
!        enddo
!        do while (j<0)
!           j=j+parallel%dims(2)
!        enddo
!        if( parallel%rankx == j )cycle
!        rcount(j+1)=diff_map%mycomm_size(2,i)-diff_map%mycomm_size(1,i)+1
!        rdispls(j+1)=diff_map%mycomm_size(1,i)-diff_map%mycomm_size(1,core_start)+1
!     enddo
!
!     scount=0
!     sdispls=0
!     do i=parallel%rankx+1,core_end
!        j=i-1
!        do while (j>=parallel%dims(2))
!           j=j-parallel%dims(2)
!        enddo
!        do while (j<0)
!           j=j+parallel%dims(2)
!        enddo
!        if( parallel%rankx == j )cycle
!        scount(j+1)=diff_map%mysend_size(2,i)-diff_map%mysend_size(1,i)+1
!        sdispls(j+1)=diff_map%mysend_size(1,i)
!     enddo
!
!     ! print*,'rank',parallel%rankx,'|',scount,'|',sdispls,'|',shape(myrho)
!     ! print*,'rank',parallel%rankx,'|',rcount,'|',rdispls,'|',shape(rho_recv)
!     CALL MPI_ALLTOALLV(myrho,scount&
!          & ,sdispls,MPI_COMPLEX16&
!          & ,rho_recv,rcount&
!          & ,rdispls,MPI_COMPLEX16&
!          & ,parallel%commx,mpinfo)
     !> deal with the local part
     ! rho_recv(parallel%mygrid_range(1):parallel%mygrid_range(2))=myrho
     j=0
     do i=parallel%mygrid_range(1),parallel%mygrid_range(2),1
        j=j+1
        ! wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=myrho(j)
        wrap_box1d(i-offset1d)=myrho(j)
     enddo

     call end_time('comm_assign1',.true.)
     call start_time('comm_assign2',.true.)
     ! !> copy
     ! do i=diff_map%mycomm_size(1,core_start),parallel%mygrid_range(1)-1
     !    wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=&
     !         & rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     !    wrap_box1d(i-offset1d)=rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     ! enddo
     ! do i=parallel%mygrid_range(2)+1,diff_map%mycomm_size(2,core_end),1
     !    wrap_box(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))=&
     !         & rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     !    wrap_box1d(i-offset1d)=rho_recv(diff_map%local_map(1,i),diff_map%local_map(2,i),diff_map%local_map(3,i))
     ! enddo

     call end_time('comm_assign2',.true.)
     call start_time('comm_assign3',.true.)
     ! print *,'lbound',lbound(wrap_box),'ubound',ubound(wrap_box)
     ! !> set periodic grid
     ! wrap_box(global_n1+1:global_n1+norder,1:global_n2,:)=&
     !      & wrap_box(1:norder,1:global_n2,:)
     ! wrap_box(1-norder:0,1:global_n2,:)=&
     !      & wrap_box(global_n1-norder+1:global_n1,1:global_n2,:)
     ! wrap_box(:,global_n2+1:global_n2+norder,:)=&
     !      & wrap_box(:,1:norder,:)
     ! wrap_box(:,1-norder:0,:)=&
     !      & wrap_box(:,global_n2-norder+1:global_n2,:)

     !> deallocate
     ! deallocate(rho_recv)
     call end_time('comm_assign3',.true.)
   ENDSUBROUTINE set_wrap_grid_per_ata_real

   SUBROUTINE set_fft_alltoallv(fft_grid_range_temp)
     implicit NONE
     INTEGER(I4B) :: fft_grid_range_temp(:)
     !> local
     INTEGER(I4B) :: i

     allocate(parallel%fft_grid_range(2,parallel%dims(2)))
     parallel%fft_grid_range=0
     CALL MPI_ALLGATHERV(fft_grid_range_temp,2,MPI_INTEGER4,&
          & parallel%fft_grid_range,&
          & (/(2,i=1,parallel%dims(2),1)/),&
          & (/(i*2,i=0,parallel%dims(2)-1,1)/), &
          & MPI_INTEGER4,parallel%commx,mpinfo)

     !> setting all to all
     call destroy_fft_alltoallv()
     allocate(parallel%fft_sdispls(parallel%dims(2)))
     allocate(parallel%fft_scount(parallel%dims(2)))
     allocate(parallel%fft_rdispls(parallel%dims(2)))
     allocate(parallel%fft_rcount(parallel%dims(2)))
     !>> send: if overlap segment of data exist,then log it, else
     !>> set to zero
     do i=1,parallel%dims(2)
        parallel%fft_scount(i)=min(parallel%mygrid_range(2),parallel%fft_grid_range(2,i))&
             &-max(parallel%mygrid_range(1),parallel%fft_grid_range(1,i))+1
        parallel%fft_sdispls(i)=max(parallel%mygrid_range(1),parallel%fft_grid_range(1,i))&
             &-parallel%mygrid_range(1)
        if(parallel%fft_scount(i)<0)then
           parallel%fft_scount(i)=0
           parallel%fft_sdispls(i)=0
        endif
     enddo
     !>> recv: similiar to send, judge depend on the overlap
     do i=1,parallel%dims(2)
        parallel%fft_rcount(i)=min(parallel%global_gridrange(2,i),&
             & parallel%fft_grid_range(2,parallel%rankx+1))-&
             & max(parallel%global_gridrange(1,i),&
             & parallel%fft_grid_range(1,parallel%rankx+1))+1
        parallel%fft_rdispls(i)=max(parallel%global_gridrange(1,i),&
             & parallel%fft_grid_range(1,parallel%rankx+1))&
             & -parallel%fft_grid_range(1,parallel%rankx+1)
        if(parallel%fft_rcount(i)<0)then
           parallel%fft_rcount(i)=0
           parallel%fft_rdispls(i)=0
        endif
     enddo

   END SUBROUTINE set_fft_alltoallv
   SUBROUTINE destroy_fft_alltoallv()
     IMPLICIT NONE
     if(allocated(parallel%fft_scount))then
        deallocate(parallel%fft_scount)
     endif
     if(allocated(parallel%fft_sdispls))then
        deallocate(parallel%fft_sdispls)
     endif
     if(allocated(parallel%fft_rcount))then
        deallocate(parallel%fft_rcount)
     endif
     if(allocated(parallel%fft_rdispls))then
        deallocate(parallel%fft_rdispls)
     endif
   END SUBROUTINE destroy_fft_alltoallv
END MODULE smpi_math_module
