subroutine f90wrap_matrixmult(a, b, c, n1, n2, n3, icase)
    implicit none
    external matrixmult
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: icase
    double precision, intent(in), dimension(n1,n2) :: a
    double precision, intent(in), dimension(n2,n3) :: b
    double precision, intent(out), dimension(n1,n3) :: c
    call matrixmult(a, b, c, n1, n2, n3, icase)
end subroutine f90wrap_matrixmult

