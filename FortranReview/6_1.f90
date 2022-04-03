program name
    implicit none
    integer :: i, j
    REAL, DIMENSION(5) :: voltage = [1, 2, 3, 4, 5]
    CHARACTER(20), DIMENSION(50) :: last_name
    REAL, DIMENSION(10) :: array1
    !implied do loop for initial array
    integer:: alength = 100
    integer, dimension(100):: array4
    !!!!!!!!!!!!!!!
    write (*, *) voltage

    do i = 1, 10
        array1(i) = real(i)
    end do
    write (*, *) array1
    call routine(alength)
end program name

subroutine routine(arg1)
    implicit none
    integer :: i, j
    integer, intent(in) :: arg1
    integer, dimension(arg1) :: array 
    array = [((0, i=1, 9), 10*j, j=1, 10)]
    write (*, *) array
end subroutine routine
