program name
    implicit none
    !!!!!!!!!!!!!
    CHARACTER(3)::file_ext = 'f'
    CHARACTER(1)::intChar
    CHARACTER(10)::a, b, c
    INTEGER:: charInt
    write (*, *) file_ext
    file_ext = '123123'
    write (*, *) file_ext
    a = '1AAAAAA'
    b = '332242r23423'
    c = a(1:3)//b(4:5)//a(6:8) ! concatenation character
    write (*, *) c
    !!!!!character intrinsic functions
    write (*, *) 'int-char'
    read (*, *) charInt
    write (*, *) charInt, ' goes to  :', achar(charInt)
    write (*, *) 'char-int'
    read (*, *) intChar
    write (*, *) intChar, ' goes to  :', iachar(intChar)
    stop 'job done'
end program name
