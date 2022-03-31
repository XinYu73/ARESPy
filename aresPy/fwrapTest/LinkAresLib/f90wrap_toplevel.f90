subroutine f90wrap_matrixmult(n, a, d, n0, n1, n2)
    implicit none
    external matrixmult
    
    integer, intent(in) :: n
    double precision, intent(in), dimension(n0,n1) :: a
    double precision, intent(inout), dimension(n2) :: d
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    integer :: n2
    !f2py intent(hide), depend(d) :: n2 = shape(d,0)
    call matrixmult(n, a, d)
end subroutine f90wrap_matrixmult

