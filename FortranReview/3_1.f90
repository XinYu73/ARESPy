program name
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real:: a, b, c, discriminant, imagPart, realPart
    logical :: flag = .TRUE.
    !!!!!!!!!!!!!!!!!!!!!!!!
    write (*, *) flag
    WRITE (*, *) 'Enter the coefficients A, B, and C: '
    READ (*, *) a, b, c
    WRITE (*, *) 'The coefficients A, B, and C are: ', a, b, c
    discriminant = b**2 - 4.*a*c
    if (discriminant > 0.) then
        WRITE (*, *) 'This equation has two real roots:'
        write (*, *) (-b + sqrt(discriminant))/(2.*a), (-b - sqrt(discriminant))/(2.*a)
    elseif (discriminant < 0.) then
        WRITE (*, *) 'X1 = ', (-b)/(2.*a), ' +i ', sqrt(abs(discriminant))/(2.*a)
        WRITE (*, *) 'X2 = ', (-b)/(2.*a), ' -i ', sqrt(abs(discriminant))/(2.*a)
    else
        WRITE (*, *) 'X1 = X2 = ', (-b)/(2.*a)
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case (int(a))
    case (:-1)
        WRITE (*, *) "It's below freezing today!"
    case (0)
        WRITE (*, *) "It's exactly at the freezing point."
    case (1:20)
        WRITE (*, *) "It's cool today."
    case (21:)
        WRITE (*, *) "It's hot today."
    case default
        write (*, *) "!^!"
    end select
    stop 'job done !!!'
    
end program name
