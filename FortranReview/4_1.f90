program name
    implicit none
    !!!!!!!!!!
    INTEGER:: n = 0
    REAL::stdDev = 0, sumX = 0, sumX2 = 0, x = 0, xBar = 0
    LOGICAL::inputFlag = .TRUE.
    !!!!!!!!!!!!!!!!
    do
        WRITE (*, *) 'Enter number: '
        READ (*, *) x
        if (x < 0) then
            WRITE (*, *) 'Exit'
            exit
        else
            WRITE (*, *) 'The number is ', x
        end if
        n = n + 1
        sumX = sumX + x
        sumX2 = sumX2 + x**2
    end do

    call calculate(sumX, sumX2, n, xBar, stdDev)

    stdDev = 0
    sumX = 0
    sumX2 = 0
    x = 0
    xBar = 0
    n = 0
    do while (x >= 0.0)
        WRITE (*, *) 'Enter number: '
        READ (*, *) x
        if (x < 0) then
            WRITE (*, *) 'Exit'
            exit
        else
            WRITE (*, *) 'The number is ', x
        end if
        n = n + 1
        sumX = sumX + x
        sumX2 = sumX2 + x**2
    end do
    call calculate(sumX, sumX2, n, xBar, stdDev)

    !!!!!!!!!

    stdDev = 0
    sumX = 0
    sumX2 = 0
    x = 0
    xBar = 0
    n = 0
    do n = 1, 30
        WRITE (*, *) 'Enter number: '
        READ (*, *) x
        if (x < 0) then
            WRITE (*, *) 'Exit'
            exit
        else
            WRITE (*, *) 'The number is ', x
        end if
        sumX = sumX + x
        sumX2 = sumX2 + x**2
    end do
    call calculate(sumX, sumX2, n, xBar, stdDev)
    stop "job done"
end program name

subroutine calculate(arg1, arg2, arg3, arg4, arg5)
    implicit none
    real, intent(in) :: arg1, arg2 !sumX,sumX2
    integer, intent(in)::arg3 ! n
    real, intent(out) :: arg4, arg5
    if (arg3 >= 2) then
        arg4 = arg1/real(arg3)
        arg5 = sqrt((real(arg3)*arg2 - arg1**2)/(real(arg3)*real(arg3 - 1)))
        WRITE (*, *) 'The mean of this data set is:', arg4
        WRITE (*, *) 'The standard deviation is: ', arg5
        WRITE (*, *) 'The number of data points is:', arg3
    else
        write (*, *) "no need to calculate"
    end if
end subroutine calculate
