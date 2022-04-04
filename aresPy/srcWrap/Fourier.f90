MODULE Fourier
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!    MODULE Fourier
!       |_SUBROUTINE PlanFFT
!       |_SUBROUTINE PlanFST
!       |_SUBROUTINE GetFFTDims
!       |_SUBROUTINE GetFFTComplexDims
!       |_INTERFACE FFT
!         |_FUNCTION ForwardFFT_4D (Private)
!         |_FUNCTION BackFFT_4D (Private)
!         |_FUNCTION ForwardFFT_3D (Private)
!         |_FUNCTION BackFFT_3D (Private)
!       |_FUNCTION ForwardFST
!       |_FUNCTION BackFST
!       |_SUBROUTINE CleanFFT
!
! DESCRIPTION:
!   This module interfaces with the Fastest Fourier Transform in the West
!   (FFTW) public library to provide Fourier transform facilities for our
!   quantities. Each Fourier transform has to be planned for first (PlanFFT)
!   then executed as many times as necessary (FFT() ) and finally cleaned up
!   (free up the memory) using CleanFFT.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   You are encouraged to consult the FFTW3 manual online at
!   http://www.fftw.org
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!   12/15/2003  Changed INTEGER*8 to INTEGER(KIND=8) to make the compiler happy
!               Also reformatted a bunch of stuff, made blurbs (GSH)
!
!------------------------------------------------------------------------------
    ! << GLOBAL >>
    USE CONSTANTS, ONLY: DP, DCP, I4B

#ifdef MPI
    USE, intrinsic :: iso_c_binding
#endif

    IMPLICIT NONE

#ifdef MPI
    include "fftw3-mpi.f03"
#else
    INCLUDE 'fftw3.inc'
#endif

    PRIVATE :: &
        ForwardFFT_4D, & ! Use FFT
        BackFFT_4D, &    ! Use FFT
        ForwardFFT_3D, & ! Use FFT
        BackFFT_3D       ! Use FFT

#ifndef MPI
    REAL(kind=DP), DIMENSION(:, :, :), ALLOCATABLE, SAVE, PRIVATE :: &
        realRA                        ! This is a permanent but local array to
    ! transform data without losing it.

    COMPLEX(DCP), DIMENSION(:, :, :), ALLOCATABLE, SAVE, PRIVATE :: &
        cplxRA                        ! This one is its complex counterpart.
    ! RA = array. Get it?
#else
    REAL(C_DOUBLE), pointer, DIMENSION(:, :, :) :: realRA
    complex(C_DOUBLE_COMPLEX), pointer, DIMENSION(:, :, :) :: cplxRA
#endif

#ifndef MPI
    INTEGER(kind=8), SAVE, PRIVATE :: &
        planRtoC, planCtoR, &         ! I'm not sure what this is.
        planRtoR_F, planRtoR_B
#else
    type(C_PTR) :: planRtoC, planCtoR, data
    INTEGER(kind=8), SAVE, PRIVATE :: &
        planRtoR_F, planRtoR_B
#endif

    INTEGER(I4B), SAVE :: &
        offset                        ! parity of realRA's X-size.

! This interface picks the right transform to perform based on the nature of
! the incomming array: if it's real the FFT is done forward, if complex the
! back transform is done. All the calls in OFDFT should be of this type:
! FFT(f).
    INTERFACE FFT
        MODULE PROCEDURE ForwardFFT_4D
        MODULE PROCEDURE BackFFT_4D
        MODULE PROCEDURE ForwardFFT_3D
        MODULE PROCEDURE BackFFT_3D
    END INTERFACE FFT

CONTAINS

    SUBROUTINE PlanFFT(dimX, dimY, dimZ, lsp)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the initialization procedure that first gets the system name as is
!   called as an argument to OFDFT, and turns it into the various input file
!   names.  Then, it calls all the programs necessary to set variables to
!   default values, then reads the geometry file to get all the variables sets
!   to the correct values.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   realRA, cplxRA, offset
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
        USE constants
        USE m_time_evaluate, ONLY: memory_sum
#ifdef MPI
        USE smpi_math_module, ONLY: parallel, set_fft_alltoallv, diff_map
#endif
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        INTEGER(I4B) :: iert
        INTEGER(I4B), INTENT(IN) :: &
            dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd
#ifdef MPI
        INTEGER(C_INTPTR_T) :: c_dimX, c_dimY, c_dimZ
        INTEGER(C_INTPTR_T) :: alloc_local, local_z, local_z_start
        INTEGER(I4B)   :: fft_grid_range_temp(2)
#endif
        LOGICAL        :: lsp(dimX, dimY, dimZ)
        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!
        !xq
        CALL cleanFFT()

#ifndef MPI
        ALLOCATE (realRA(dimX, dimY, dimZ))
        ALLOCATE (cplxRA(dimX/2 + 1, dimY, dimZ))
        call memory_sum('fft_realRA cplxRA', real(size(realRA) + size(cplxRA)*2, DP)*DP)
#else
        !> fft comm
        c_dimX = dimX
        c_dimY = dimY
        c_dimZ = dimZ
        ! if(parallel%numprocs>dimZ)then
        parallel%commfft = parallel%commx
        ! else
        ! parallel%commfft=parallel%comm
        ! endif
        !> allocate work arrays
        call fftw_mpi_init()
        alloc_local = fftw_mpi_local_size_3d(c_dimZ, c_dimY, c_dimX/2 + 1&
             &, parallel%commfft&
             &, local_z, local_z_start)
        parallel%local_z = local_z
        parallel%local_z_start = local_z_start
        data = fftw_alloc_real(alloc_local*2)
        call c_f_pointer(data, realRA, [2*(c_dimX/2 + 1), c_dimY, local_z])
        call c_f_pointer(data, cplxRA, [c_dimX/2 + 1, c_dimY, local_z])
        !> setting grid transform information
        ! fft_grid_range_temp(1)=(local_z_start)*c_dimX*c_dimY+1
        ! fft_grid_range_temp(2)=(local_z_start+local_z)*c_dimX*c_dimY
        fft_grid_range_temp(1) = diff_map%nz_map(1, local_z_start + 1)
!xqtest
!  fft_grid_range_temp(2)=diff_map%nz_map(2,local_z_start+local_z)
        if (local_z == 0 .and. local_z_start == 0) then
            fft_grid_range_temp(2) = fft_grid_range_temp(1) - 1
        else
            fft_grid_range_temp(2) = diff_map%nz_map(2, local_z_start + local_z)
        end if
        CALL set_fft_alltoallv(fft_grid_range_temp)
        ! allocate(parallel%fft_grid_range(2,parallel%dims(1)))
        ! parallel%fft_grid_range=0
        ! CALL MPI_ALLGATHERV(fft_grid_range_temp,2,MPI_INTEGER4,&
        !      & parallel%fft_grid_range,(/(2,i=1,dims(1),1)/),&
        !      & (/(i*2,i=0,dims(1)-1,1)/), &
        !      & MPI_INTEGER4,parallel%commx,status)
#endif

#ifndef MPI
        CALL dfftw_plan_dft_r2c_3d(planRtoC, dimX, dimY, dimZ, realRA, cplxRA, &
                                   FFTW_ESTIMATE)
        CALL dfftw_plan_dft_c2r_3d(planCtoR, dimX, dimY, dimZ, cplxRA, realRA, &
                                   FFTW_ESTIMATE)
#else
        planRtoC = fftw_mpi_plan_dft_r2c_3d(c_dimZ, c_dimY, c_dimX, realRA, cplxRA&
             &, parallel%commfft, FFTW_ESTIMATE)
        planCtoR = fftw_mpi_plan_dft_c2r_3d(c_dimZ, c_dimY, c_dimX, cplxRA, realRA&
             &, parallel%commfft, FFTW_ESTIMATE)
        ! call fftw_mpi_execute_dft_r2c(planRtoC,realRA,cplxRA)
        ! call fftw_mpi_execute_dft_c2r(planCtoR,realRA,cplxRA)
#endif

        offset = MOD(dimX, 2)

    END SUBROUTINE PlanFFT

    SUBROUTINE PlanFST(dimX, dimY, dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the same as PlanFFT, except for the Fast Sine Transform.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        INTEGER(I4B), INTENT(IN) :: &
            dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

        !>> INTERNAL VARIABLES <<!
        REAL(kind=DP), DIMENSION(dimX, dimY, dimZ) :: &
            realRA_RtoR

        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

        CALL dfftw_plan_r2r_3d(planRtoR_F, dimX, dimY, dimZ, &
                               realRA_RtoR, realRA_RtoR, &
                               FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, &
                               FFTW_ESTIMATE)

        CALL dfftw_plan_r2r_3d(planRtoR_B, dimX, dimY, dimZ, &
                               realRA_RtoR, realRA_RtoR, &
                               FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, &
                               FFTW_ESTIMATE)

    END SUBROUTINE PlanFST

    SUBROUTINE GetFFTDims(dimX, dimY, dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (real-space part)
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/25/2006  Added (GSH)
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        INTEGER(I4B), INTENT(OUT) :: &
            dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

        dimX = SIZE(realRA, 1)
        dimY = SIZE(realRA, 2)
        dimZ = SIZE(realRA, 3)

    END SUBROUTINE GetFFTDims

    SUBROUTINE GetFFTComplexDims(dimX, dimY, dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (reciprocal space part)
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/26/2006 Added (GSH)
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        INTEGER(I4B), INTENT(OUT) :: &
            dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

        dimX = SIZE(cplxRA, 1)
        dimY = SIZE(cplxRA, 2)
        dimZ = SIZE(cplxRA, 3)

    END SUBROUTINE GetFFTComplexDims

    FUNCTION ForwardFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT
!   interface instead. It performs the transformation of a real 4-dimensional
!   array into its complex 4-dimensional transform. The first dimension is
!   halved.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        REAL(kind=DP), DIMENSION(:, :, :, :) :: &
            array             ! The array to transform

        COMPLEX(DCP), DIMENSION(SIZE(array, 1)/2 + 1, SIZE(array, 2), &
                                SIZE(array, 3), SIZE(array, 4)) :: &
            transform         ! The answer

        !>> INTERNAL VARIABLES <<!
        INTEGER(I4B) :: &
            is                ! dummy counter

        REAL(kind=DP) :: &
            normalize         ! normalization constant?

        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!
        DO is = 1, SIZE(array, 4)
            realRA = array(1:SIZE(array, 1), 1:SIZE(array, 2), 1:SIZE(array, 3), is)
            CALL dfftw_execute(planRtoC)
            transform(1:SIZE(transform, 1), 1:SIZE(transform, 2), 1:SIZE(transform, 3), &
                      is) = cplxRA
        END DO !is

        ! The forward transform needs to be renormalized afterwards.
        normalize = REAL(SIZE(array, 4), kind=DP)/REAL(SIZE(array), kind=DP)
        transform = transform*normalize

    END FUNCTION ForwardFFT_4D

    FUNCTION BackFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather
!   through the FFT interface. It performs the reverse Fourier transform of
!   a complex function over the half-box in reciprocal space back to real
!   space. It acts on 4-dimensional arrays, the fourth dimension being spin.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003    File created.  (VLL)
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        COMPLEX(DCP), DIMENSION(:, :, :, :) :: &
            array             ! The array to be back FFT'd

        ! The x-size of the returned array is computed from the reciprocal-space
        ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
        ! k+1 in reciprocal space. We solve the problem by storing the parity.
        REAL(kind=DP), DIMENSION(2*(SIZE(array, 1) - 1) + offset, SIZE(array, 2), &
                                 SIZE(array, 3), SIZE(array, 4)) :: &
            transform         ! The answer

        !>> INTERNAL VARIABLES <<!
        INTEGER(I4B) :: &
            is                ! dummy counter

        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!
        DO is = 1, SIZE(array, 4)
            cplxRA = array(1:SIZE(array, 1), 1:SIZE(array, 2), 1:SIZE(array, 3), is)
            CALL dfftw_execute(planCtoR)
            transform(1:SIZE(transform, 1), 1:SIZE(transform, 2), 1:SIZE(transform, 3), &
                      is) = realRA
        END DO

    END FUNCTION BackFFT_4D

    FUNCTION ForwardFFT_3D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT
!   interface instead. It performs the transformation of a real 4-dimensional
!   array into its complex 4-dimensional transform. The first dimension is
!   halved.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
#ifdef MPI
        USE smpi_math_module, ONLY: parallel, mpinfo, mpi_integer4, mpi_sum
#endif
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        REAL(kind=DP), DIMENSION(:, :, :) :: &
            array             ! The array to transform

        COMPLEX(DCP), DIMENSION(SIZE(array, 1)/2 + 1, SIZE(array, 2), &
                                SIZE(array, 3)) :: &
            transform         ! The answer
        !> local
        INTEGER(I4B)                     :: size_local, size_global
        realRA(1:size(array, 1), :, :) = array
#ifndef MPI
        CALL dfftw_execute(planRtoC)
#else
        ! print *,array(size(array),1,11)
        CALL fftw_mpi_execute_dft_r2c(planRtoC, realRA, cplxRA)
#endif
        transform = cplxRA

        ! The forward transform needs to be renormalized afterwards.
#ifdef MPI
        size_local = REAL(SIZE(array), kind=DP)
        CALL MPI_ALLREDUCE(size_local, size_global, 1, mpi_integer4, mpi_sum, parallel%commfft, mpinfo)
        transform = transform/size_global
#else
        transform = transform/REAL(SIZE(array), kind=DP)
#endif

    END FUNCTION ForwardFFT_3D

    FUNCTION BackFFT_3D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather
!   through the FFT interface. It performs the reverse Fourier transform of a
!   complex function over the half-box in reciprocal space back to real
!   space. It acts on 3-dimensional arrays.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        COMPLEX(DCP), DIMENSION(:, :, :) :: &
            array             ! The array to be back FFT'd

        ! The x-size of the returned array is computed from the reciprocal-space
        ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
        ! k+1 in reciprocal space. We solve the problem by assuming an odd real size
        REAL(kind=DP), DIMENSION(2*(SIZE(array, 1) - 1) + offset, SIZE(array, 2), &
                                 SIZE(array, 3)) :: &
            transform         ! The answer

        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

        ! print*,'shape array',shape(array)
        ! print*,'shape realRA',shape(realRA)
        ! print*,'shape transform',shape(transform)
        cplxRA = array
#ifndef MPI
        CALL dfftw_execute(planCtoR)
#else
        CALL fftw_mpi_execute_dft_c2r(planCtoR, cplxRA, realRA)
#endif
        transform = realRA(1:size(transform, 1), :, :)

    END FUNCTION BackFFT_3D

    SUBROUTINE ForwardFST(array)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!
        REAL(kind=DP), DIMENSION(:, :, :, :) :: &
            array             ! The array to transform

        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

        CALL dfftw_execute_r2r(planRtoR_F, array(:, :, :, 1), array(:, :, :, 1))

        ! The forward transform needs to be renormalized afterwards.
        array = array/REAL(SIZE(array)*8, kind=DP)

    END SUBROUTINE ForwardFST

    SUBROUTINE BackFST(array)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
        IMPLICIT NONE
        !>> EXTERNAL VARIABLES <<!

        REAL(kind=DP), DIMENSION(:, :, :, :) :: &
            array             ! The array to transform

        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

        CALL dfftw_execute_r2r(planRtoR_B, array(:, :, :, 1), array(:, :, :, 1))

    END SUBROUTINE BackFST

    SUBROUTINE CleanFFT
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine is called at the end of the run to free the memory
!   associated with the plan.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   realRA, cplxRA
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
        USE m_time_evaluate, ONLY: memory_free
#ifdef MPI
        USE smpi_math_module, ONLY: parallel
#endif
        IMPLICIT NONE
        INTEGER(I4B)  :: iert
        !>> EXTERNAL VARIABLES <<!
        !>> INTERNAL VARIABLES <<!
        !>> INITIALIZATION <<!
        !>> FUNCTION BODY <<!

#ifndef MPI
        CALL dfftw_destroy_plan(planRtoC)
        CALL dfftw_destroy_plan(planCtoR)
        IF (ALLOCATED(realRA)) THEN
            call memory_free('fft_realRA', real(size(realRA), DP)*DP)
            DEALLOCATE (realRA)
        END IF
        IF (ALLOCATED(cplxRA)) THEN
            call memory_free('fft_cplxRA', real(size(cplxRA), DP)*DP*2)
            DEALLOCATE (cplxRA)
        END IF
#else
        if (allocated(parallel%fft_grid_range)) THEN
            deallocate (parallel%fft_grid_range)
        end if
        call fftw_destroy_plan(planRtoC)
        call fftw_destroy_plan(planCtoR)
        call fftw_free(data)
        call fftw_mpi_cleanup()
#endif

    END SUBROUTINE CleanFFT

END MODULE Fourier
