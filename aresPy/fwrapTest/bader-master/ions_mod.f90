! Copyright 2009
! Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!
! Bader is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! A copy of the GNU General Public License is available at
! http://www.gnu.org/licenses/

!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module with ion data structure
!-----------------------------------------------------------------------------------!

MODULE ions_mod
  USE kind_mod
  IMPLICIT NONE

  TYPE :: ions_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir,r_lat
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: ion_chg
    REAL(q2),DIMENSION(3,3) :: lattice,dir2car,car2dir
    INTEGER,ALLOCATABLE,DIMENSION(:) :: num_ion
    INTEGER,ALLOCATABLE,DIMENSION(:) :: atomic_num
    CHARACTER*330:: name_ion
    INTEGER :: niontypes,nions
  END TYPE

  PRIVATE
  PUBLIC :: ions_obj, ions_obj_initialise, ions_obj_finalise

  !-----------------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------------!
! added constructors and destructors for Python wrapper (James Kermode)
!-----------------------------------------------------------------------------------!

    SUBROUTINE ions_obj_initialise(ions, nions, niontypes)

      TYPE(ions_obj), INTENT(INOUT) :: ions
      INTEGER, INTENT(IN) :: nions, niontypes

      call ions_obj_finalise(ions)
      ions%nions = nions
      ions%niontypes = niontypes
      allocate(ions%r_dir(nions, 3), ions%r_car(nions, 3), ions%r_lat(nions, 3))
      allocate(ions%ion_chg(nions), ions%atomic_num(nions))
      allocate(ions%num_ion(niontypes))

    END SUBROUTINE ions_obj_initialise

    SUBROUTINE ions_obj_finalise(ions)
      TYPE(ions_obj), INTENT(INOUT) :: ions

      if (allocated(ions%r_car)) deallocate(ions%r_car)
      if (allocated(ions%r_dir)) deallocate(ions%r_dir)
      if (allocated(ions%r_lat)) deallocate(ions%r_lat)
      if (allocated(ions%ion_chg)) deallocate(ions%ion_chg)
      if (allocated(ions%num_ion)) deallocate(ions%num_ion)
      if (allocated(ions%atomic_num)) deallocate(ions%atomic_num)
      ions%nions = 0
      ions%niontypes = 0
      
    END SUBROUTINE ions_obj_finalise  

END MODULE ions_mod
