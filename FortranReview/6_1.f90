program name
    implicit none
    integer :: i
    REAL, DIMENSION(5) :: voltage = [1, 2, 3, 4, 5]
    CHARACTER(20), DIMENSION(50) :: last_name
    REAL, DIMENSION(10) :: array1
    !!!!!!!!!!!!!!!
    write (*, *) voltage

    do i = 1, 10
        array1(i) = real(i)
    end do
    write (*, *) array1
end program name
