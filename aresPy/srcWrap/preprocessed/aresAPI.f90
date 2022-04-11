module aresAPI
    USE constants
    USE Struct_module
    implicit none
    type aresOut
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: forces
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: stress
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: poscar
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: pos
        REAL(DP), ALLOCATABLE, DIMENSION(:, :, :, :) :: chargeRho
    end type aresOut
contains
    subroutine init_alloc_arrays(dertype, nnatom)
        USE parameters, ONLY: nspin
        USE grid_module, ONLY: ni1 => global_n1, ni2 => global_n2, ni3 => global_n3 &
                    & , Sphere2Cubic, sph => grid, nps => n
        implicit none
        type(aresout), INTENT(inout) :: dertype
        INTEGER(I4B), intent(in) :: nnatom
        REAL(DP) :: rhot(ni1, ni2, ni3, 2)
        !force
        ALLOCATE (dertype%forces(3, nnatom))
        dertype%forces = struct%forces
        !stress
        ALLOCATE (dertype%stress(3, 3))
        dertype%stress = struct%stress
        !poscar
        ALLOCATE (dertype%poscar(3, nnatom))
        dertype%poscar = struct%poscar
        !pos
        ALLOCATE (dertype%pos(3, nnatom))
        dertype%pos = struct%pos
        !electron density????????????????
        ALLOCATE (dertype%chargeRho(ni1, ni2, ni3, 2))
        IF (nspin == 1) THEN
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1), rhot(:, :, :, 1))
            rhot(:, :, :, 2) = 0._DP
        ELSE
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1) + sph%rhoS(:, 2), rhot(:, :, :, 1))
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1) - sph%rhoS(:, 2), rhot(:, :, :, 2))
        END IF
        dertype%chargeRho = rhot
        !
    end subroutine init_alloc_arrays

    subroutine destroy_alloc_arrays(dertype)
        type(aresout), INTENT(inout) :: dertype
        if (allocated(dertype%forces)) deallocate (dertype%forces)
        if (allocated(dertype%stress)) deallocate (dertype%stress)
        if (allocated(dertype%poscar)) deallocate (dertype%poscar)
        if (allocated(dertype%pos)) deallocate (dertype%pos)
        if (allocated(dertype%chargeRho)) deallocate (dertype%chargeRho)
    end subroutine destroy_alloc_arrays
end module aresAPI
