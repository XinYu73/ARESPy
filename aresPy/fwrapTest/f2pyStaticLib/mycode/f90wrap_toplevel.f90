subroutine f90wrap_matrixmult(a, b, c, icase, n0, n1, n2, n3, n4, n5)
    implicit none
    external matrixmult
    
    double precision, intent(in), dimension(n0,n1) :: a
    double precision, intent(in), dimension(n2,n3) :: b
    double precision, intent(inout), dimension(n4,n5) :: c
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: icase
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    integer :: n2
    !f2py intent(hide), depend(b) :: n2 = shape(b,0)
    integer :: n3
    !f2py intent(hide), depend(b) :: n3 = shape(b,1)
    integer :: n4
    !f2py intent(hide), depend(c) :: n4 = shape(c,0)
    integer :: n5
    !f2py intent(hide), depend(c) :: n5 = shape(c,1)
    call matrixmult(a, b, c, n1, n2, n3, icase)
end subroutine f90wrap_matrixmult

