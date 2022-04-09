PROGRAM ARES
    !##########################BEAST WAR#############################!
    !*BEAST   : real space Kohn-Sham DFT Software                    !
    !*Author  : Qiang Xu                                             !
    !*E-mail  : xq@calypso.cn                                        !
    !*Date    : 2017/07/01                                           !
    !*Diag H  : ARPACK                                               !
    !################################################################!
    !USE omp_lib
    USE constants
    USE read_module
    USE parameters
    USE Begin_module, ONLY: Initial_grid_pbc
    USE potential_module, ONLY: vlpp
    USE scf_module, ONLY: electronicSCF
    USE output_module
#ifdef MPI
    USE SMPI_MATH_MODULE
    USE ScaLapack_module, ONLY: Init_scala, L_useless
    USE m_time_evaluate, ONLY: memory_sum
#endif
    IMPLICIT NONE
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI
    !> Initialize the mpi environment {comm, size, rank}
    CALL smpi_init()
#endif

#ifdef MPI
    IF (PARALLEL%ISROOT) THEN
#endif
        WRITE (6, *) '=========================== ARES-PBC =============================='
        WRITE (6, *) '>>>>>>>>> Ab initio Real space Electronic Structer solver <<<<<<<<<'
#ifdef MPI
    END IF
    CALL start_time('[Total Time]', .TRUE.)
    CALL start_time('[SCF Time]', .TRUE.)
#endif

    !> Initialize the calculating information {parameters,
    !> position, pseudopotential}
    !1) Read file
    CALL read_file('ares.in')
#ifdef MPI
    CALL Init_scala()
#endif
    !2) Build real-space grids
    CALL Initial_grid_pbc()
    !3) Build externel potential
    CALL vlpp()
    !4) SCF iterations
    CALL electronicSCF()
#ifdef MPI
    CALL end_time('[SCF Time]', .TRUE.)
#endif
    !5) Output
    CALL Output()
#ifdef MPI
    CALL end_time('[Total Time]', .TRUE.)
    IF (PARALLEL%ISROOT) THEN
        CALL write_time('[Total Time]', .TRUE.)
#endif
        WRITE (*, *) '=================================Well Done==================================='
#ifdef MPI
    END IF
    CALL MPI_Barrier(parallel%comm, mpinfo)
    CALL SMPI_EXIT()
#endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END PROGRAM ARES
