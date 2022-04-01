
program test
    implicit none
    !##########Paramwter
    real, parameter::pi = 3.1415926
    ! ##########
    INTEGER :: i = 0, j = 0, k = 0
    CHARACTER(10) :: first, last
    ! ##########
    WRITE (*, *) 'Enter the numbers to multiply: '
    READ (*, *) i, j
    WRITE (*, *) k*j
    WRITE (*, *) 'Result = ', acosd(exp(real(i))/3.0) !argument of exp is real
    WRITE (*, *) 'Result = ', acos(exp(real(i))/3.0)
    ! ##########
    stop 'succeed'
end program test
