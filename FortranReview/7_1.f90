program name
    implicit none
    INTEGER :: i
    real, dimension(5) :: a = 0
    character(10)::str = 'xinyu'
    real:: aa = 10.0, bb = 121.0, result = 0
    integer:: error = 0
    call sub3(aa, bb, result, error)
    write (*, *) 'error: ', error, 'result: ', result
    call sub1(a, 5, 5)

    write (*, *) a
    call sub2(str)
end program name
subroutine sub1(a, ndim, n)
    implicit none
    integer, intent(in) :: ndim
    real, intent(out), dimension(ndim)::  a
    integer, intent(in) :: n
    integer :: i

    do i = 1, n
        a(i) = i**2
    end do
end subroutine sub1

subroutine sub2(string)
    implicit none
    character(*), intent(in) :: string
    write (*, *) "string", len(string)
end subroutine sub2

subroutine sub3(arg1, arg2, result, error)
    implicit none
    real, intent(in) :: arg1, arg2
    real, intent(out) ::  result
    integer, intent(out) :: error
!!!!!!!!!
    if (arg1 - arg2 >= 0.0) then
        result = sqrt(arg1 - arg2)
        error = 0
    else
        result = 0
        error = 1
    end if
end subroutine sub3
