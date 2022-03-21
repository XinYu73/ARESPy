# 1 "test.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "test.f90"
module module_calcul
  type type_ptmes
     integer :: y
  end type type_ptmes
  
  type array_type
     type(type_ptmes) :: x(2)
  end type array_type

  type(array_type), dimension(10), target :: xarr
  
contains
  subroutine recup_point(x)
    type(array_type) :: x
    return
  end subroutine recup_point
  subroutine xytest()
        write(*,*) "hh"
  end subroutine xytest
end module module_calcul
