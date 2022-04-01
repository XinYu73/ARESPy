# f2py 链接库

## Preparation

### install LAPACK and BLAS

1. compile with flag -fPIC，f2py needs lib to be position-independent
2. compiler: ifort  and icc
3. lib path

    ```bash
    /work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/lapack-3.10.0
    /work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/BLAS-3.10.0
    ```

### Link LAPACK to Fortran Program

```fortran
    !./test.f90
    program test
    implicit none
    integer,parameter :: n= 3
    real(kind=8), dimension(n) :: x,b
    real(kind=8), dimension(n,n) :: a
    integer :: i,info,lda,ldb ,nrhs
    integer ,dimension(n) :: ipiv
    a = reshape((/1.0,45.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
    b=(/2.0,2.0,3.5/)
    write(*,*) 'j(1) = '
    x=b
    nrhs=1
    lda=n
    ldb=n
    call dgesv(n,nrhs,a,ldb,ipiv,x,ldb,info)   ! solve Ax = B,return x
    write(*,*) x
    end program test
```

compile

```make
    LAPACKLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/lapack-3.10.0 -llapack -lrefblas  -ltmglib
    BLASLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/BLAS-3.10.0 -lblas

    build:
    $(FC) -c  test.f90
    $(FC) -o test test.o $(LAPACKLIB) $(BLASLIB)
```

test:

```bash
    ./test
    -1.744186046511630E-002   1.11821705426357      -0.350775193798449
```


## Wrap

```fortran
    subroutine my(n, a, d)
    INTEGER, intent(in) :: n
    DOUBLE PRECISION, intent(in) :: a(n,n)
    DOUBLE PRECISION, intent(inout) :: d(n)
    integer :: i,info,lda,ldb ,nrhs
    integer ,dimension(n) :: ipiv
    CALL dgesv(n,nrhs,a,ldb,ipiv,d,ldb,info)
    end subroutine my
```

Wrap
    
```make
    LAPACKLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/lapack-3.10.0 -llapack -lrefblas  -ltmglib
    BLASLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/BLAS-3.10.0 -lblas
    FC=ifort
    .PHONY: buildLib Wrap
    Wrap:
    @f2py my.f90 -m my -h my.pyf --overwrite-signature
    @f2py --fcompiler=intelem --f90exec=$(FC) --f77exec=$(FC) -c my.pyf my.f90 -m my $(LAPACKLIB) $(BLASLIB)
```

```bash
    make Wrap
```

then I have

```bash
    my.cpython-38-x86_64-linux-gnu.so
```

test with python

```python
    #./test.py
    import numpy as np
    from time import *
    import sys
    import my

    n1 = int(3)
    n2 = int(3)
    n3 = int(3)

    A = np.random.rand(n1,n2)
    B = np.random.rand(n2,n3)
    C = np.random.rand(1,3)
    D = np.zeros(n1)
    #-------
    # Case 1
    #-------
    print(my.__doc__)
    beg = time()
    print(A)
    D=my.my(A,C)
    print(D)
    end = time()
```

run test.py

```bash
    python test.py
```

result

```bash
    This module 'my' is auto-generated with f2py (version:1.21.2).
    Functions:
    my(a,d,n=shape(a,0))
    .
    [[0.15121598 0.97648753 0.59436203]
    [0.89983464 0.80702176 0.22636051]
    [0.58806946 0.84986581 0.31170033]]

    Intel MKL ERROR: Parameter 2 was incorrect on entry to DGESV . 
    None
```
