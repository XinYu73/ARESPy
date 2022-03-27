# How to wrap fortran module which include static lib some

## solution 1 with f2py

1. In testlib, we have three methods to do a matrix multiplication

    ```fortran
    matMult_func.f90
    matMult_loop1.f90
    matMult_loop2.f90 
    ```

    we fist bulid a static library

    ```bash
    ifort -c -fPIC *.f90
    ar rc libmatMultiplication.a *.o
    ```

2. In mycode ,we have a module which call function in libmutiplication

    ```bash
    matrixMult.f90
    ```

    ```bash
    f2py --f90exec=ifort -L../testlib -lmatMultiplication -c matrixMult.f90 -m matrixMult
    ```

    ```bash
    -L../testlib tells wher the lib is
    -lmatMutiplication : notice that we exculde the lib in  libmatMultiplication.a
    ```

    then we got python lib called matrixMult
3. test with python f2py_matrixMult.py

    ```bash
    [[0.28123195 0.26432832 0.13521483 0.52000107 0.48099121]
    [0.87248204 0.65864055 0.43100443 1.41080808 1.39432353]
    [1.01929922 0.54439239 0.51959664 1.3659256  1.4924507 ]
    [0.11769574 0.1473179  0.05396834 0.26364367 0.22355003]]
    Loop1: time for AB(4, 5) = A(4, 2) B(2, 5) is 0.0014388561248779297 s
    Loop2: time for AB(4, 5) = A(4, 2) B(2, 5) is 1.9788742065429688e-05 s
    matmul function: time for AB(4, 5) = A(4, 2) B(2, 5) is 7.867813110351562e-06 s
    ```

## Solution 2 with f90wrap

1. build the library as solution 1.
