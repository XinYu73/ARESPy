# 1 "AresMainAPI.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "AresMainAPI.f90"
module AresMainAPI
    use constants
    use struct_module
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
        implicit none
        type(aresout), INTENT(inout) :: dertype
        INTEGER(I4B), intent(in) :: nnatom
        !force
        ALLOCATE (dertype%forces(3, nnatom))
        !stress
        ALLOCATE (dertype%stress(3, 3))
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
    end subroutine assignment

    subroutine destroy_alloc_arrays(dertype)
        type(aresout), INTENT(inout) :: dertype
        if (allocated(dertype%forces)) deallocate (dertype%forces)
        if (allocated(dertype%stress)) deallocate (dertype%stress)
        if (allocated(dertype%poscar)) deallocate (dertype%poscar)
        if (allocated(dertype%pos)) deallocate (dertype%pos)
    end subroutine destroy_alloc_arrays
end module AresMainAPI
