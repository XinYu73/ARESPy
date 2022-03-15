! https://het.as.utexas.edu/HET/Software/Numpy/f2py/python-usage.html
SUBROUTINE ZADD(a,b,c,n)
    INTEGER,INTENT(IN),DIMENSION(n) :: a,b
    INTEGER,INTENT(OUT),DIMENSION(n) :: c
    INTEGER :: n
    c = a+b
END SUBROUTINE ZADD