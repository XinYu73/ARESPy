module aresAPI
    USE constants
    USE struct_module, ONLY: natom => nat
    implicit none
    type, public :: embed_base
        ! type(input_base)                :: input
        ! type(tddft_base)                :: tddft
        ! real(kind=dp), allocatable      :: extpot(:, :)
        ! real(kind=dp)                   :: extene = 0.0
        ! integer                         :: exttype = 0
        real(kind=dp), allocatable      :: extforces(:, :)
        !     real(kind=dp)                   :: extstress(3, 3)
        !     logical                         :: initial = .true.
        !     real(kind=dp)                   :: mix_coef = -1.0
        !     logical                         :: finish = .false.
        !     real(kind=dp)                   :: etotal = 0.0
        !     real(kind=dp)                   :: dnorm = 1.0
        !     logical                         :: lewald = .true.
        !     logical                         :: nlpp = .true.
        !     real(kind=dp)                   :: diag_conv = 1.D-2
        !     logical                         :: ldescf = .false.
        ! !! add scf correction energy
        !     logical                         :: iterative = .false.
        ! !! add correction for variational energy
        !     logical                         :: lmovecell = .false.
        ! !! allow change the cell
    CONTAINS
        !--------------------------------------------------------------------------------
        ! PROCEDURE :: allocate_extpot => allocate_extpot_class
        PROCEDURE :: allocate_extforces
    end type embed_base
    TYPE(embed_base), public :: messenger
contains
    ! For Force
    SUBROUTINE allocate_extforces(embed)
        IMPLICIT NONE
        TYPE(embed_base), INTENT(INOUT) :: embed
        !
        IF (ALLOCATED(embed%extforces)) THEN
            IF (SIZE(embed%extforces, 2) /= nat) DEALLOCATE (embed%extforces)
        END IF
        IF (.NOT. ALLOCATED(embed%extforces)) THEN
            ALLOCATE (embed%extforces(3, nat))
            embed%extforces = 0.0_DP
        END IF
    END SUBROUTINE
    !
    SUBROUTINE ares_set_extforces(embed)
        !
        use struct_module, only: struct
        IMPLICIT NONE
        TYPE(embed_base), INTENT(INOUT) :: embed
        !
        call embed%allocate_extforces()
        embed%extforces(:, :) = struct%forces(:, 1:nat)
        !
    END SUBROUTINE
    ! Force Done
end module aresAPI
