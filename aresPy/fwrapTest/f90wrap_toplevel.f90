subroutine f90wrap_fib(a, n)
    implicit none
    external fib
    
    real :: a
    integer :: n
    call fib(a, n)
end subroutine f90wrap_fib

