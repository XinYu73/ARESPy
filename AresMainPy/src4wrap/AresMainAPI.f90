module AresMainAPI
    use constants
    use struct_module
    use smpi_math_module
    USE parameters, ONLY: elements, Lpbc, ntype
    USE math, ONLY: dir2car, car2dir, det, change_case &
           &, atom_mass, atom_effcharge
    USE struct_module
    USE m_time_evaluate, only: memory_sum
    USE Read_module, ONLY: ResetLattice
    implicit none
    type aresOut
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: forces
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: poscar
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: pos
        REAL(DP), ALLOCATABLE, DIMENSION(:, :, :, :) :: chargeRho
        REAL(DP), DIMENSION(3, 3)  :: stress
        REAL(DP)                  :: apilat_mat(3, 3)
        REAL(DP)                  :: apilat_para(6)
        INTEGER(I4B)              :: comm
        INTEGER(I4B)              :: myid
        INTEGER(I4B)              :: numprocs
        INTEGER(I4B)              :: rootid
    end type aresOut
contains
    subroutine init_alloc_arrays(dertype, nnatom)
        implicit none
        type(aresout), INTENT(inout) :: dertype
        INTEGER(I4B), intent(in) :: nnatom
        !force
        ALLOCATE (dertype%forces(3, nnatom))
        !poscar
        ALLOCATE (dertype%poscar(3, nnatom))
        !pos
        ALLOCATE (dertype%pos(3, nnatom))
        !
    end subroutine init_alloc_arrays

    !
    subroutine assignment(dertype)
        implicit none
        type(aresout), INTENT(inout) :: dertype
        dertype%forces = struct%forces
        dertype%stress = struct%stress
        dertype%poscar = struct%poscar
        dertype%pos = struct%pos
        dertype%comm = parallel%comm
        dertype%apilat_mat = lat_mat
        dertype%apilat_para = lat_para
    end subroutine assignment

    subroutine destroy_alloc_arrays(dertype)
        type(aresout), INTENT(inout) :: dertype
        if (allocated(dertype%forces)) deallocate (dertype%forces)
        if (allocated(dertype%poscar)) deallocate (dertype%poscar)
        if (allocated(dertype%pos)) deallocate (dertype%pos)
    end subroutine destroy_alloc_arrays

    ! #! input info
    subroutine updateIons(pos, lattice)
        USE parameters, only: ntype
        ! WE SIMPLILY CONSIDER NORMAL POSCAR GET IN
        REAL(DP), INTENT(IN) :: pos(:, :) !ase default car
        REAL(DP), INTENT(IN), OPTIONAL  :: lattice(3, 3)!normal
        INTEGER(I4B) :: nty
        INTEGER(I4B)  :: i, j, natom_3, nty_4, nty_3
        nty = ntype
        IF (parallel%isroot) THEN
            !turn to atomic unit
            lat_mat(:, :) = lattice*ang2bohr !neglect lat_ratio
            CALL car2dir(pos*ang2bohr, struct%pos, lat_mat)
            IF (.NOT. Lpbc) THEN
                CALL ResetLattice()
            END IF
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!
        CALL MPI_BCAST(natom, 1, MPI_INTEGER4, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 1
        CALL MPI_BCAST(struct%nati, nty, MPI_INTEGER4, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 2
        CALL MPI_BCAST(struct%eleid, nty + 1, MPI_INTEGER4, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 3
        natom_3 = 3*natom
        CALL MPI_BCAST(struct%pos, natom_3, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 4
        CALL MPI_BCAST(struct%poscar, natom_3, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 5
        CALL MPI_BCAST(struct%mass, nty, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 6
        nty_4 = 4*nty
        nty_3 = 3*nty
        CALL MPI_BCAST(struct%zeta, nty_4, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 7
        CALL MPI_BCAST(struct%prinq, nty_4, MPI_INTEGER4, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 8
        CALL MPI_BCAST(struct%Lmax, nty, MPI_INTEGER4, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 9
        CALL MPI_BCAST(struct%elements, nty_3, MPI_CHARACTER, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 10
        CALL MPI_BCAST(lat_mat, 9, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        write(*,*) 11
        call memory_sum("read_pos_local", real(nty, DP)*I4B + (nty + 1)*I4B + 9*DP)
        write(*,*) 12
    end subroutine updateIons
end module AresMainAPI
