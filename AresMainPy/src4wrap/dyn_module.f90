MODULE dyn_module!{{{
  USE CONSTANTS
  IMPLICIT NONE
  type  ::  dynamics
     logical  :: lstop
     real(dp),dimension(:,:),allocatable  :: posion
     real(dp),dimension(:,:),allocatable  :: posioc
     real(dp),dimension(:,:),allocatable  :: d2
     real(dp),dimension(:,:),allocatable  :: d2c
     real(dp),dimension(:,:),allocatable  :: d3
     real(dp),dimension(3,3)   :: A
     real(dp),dimension(3,3)   :: B
     real(dp),dimension(3,3)   :: AC
     real(dp)                  :: POTIM
     real(dp)                  :: EDIFFG
     real(dp)                  :: PSTRESS
     integer(i4b)              :: IBRION
     integer(i4b)              :: ISIF
     integer(i4b)              :: NFREE
  end type dynamics
  type (dynamics),save         :: dyn
contains
  subroutine create_dyn(na,dyn)
    implicit none
    integer(i4b),intent(in)  :: na
    type(dynamics)  :: dyn
    allocate(dyn%posion(3,na))
    allocate(dyn%posioc(3,na))
    allocate(dyn%d2(3,na))
    allocate(dyn%d2c(3,na))
    allocate(dyn%d3(3,na))
  end subroutine create_dyn

  subroutine destroy_dyn(dyn)
    implicit none
    type(dynamics)  :: dyn
    if (allocated(dyn%posion) )  deallocate(dyn%posion)
    if (allocated(dyn%posioc) )  deallocate(dyn%posioc)
    if (allocated(dyn%d2)     )  deallocate(dyn%d2    )
    if (allocated(dyn%d2c)    )  deallocate(dyn%d2c   )
    if (allocated(dyn%d3)     )  deallocate(dyn%d3    )
  end subroutine destroy_dyn

end module dyn_module !}}}