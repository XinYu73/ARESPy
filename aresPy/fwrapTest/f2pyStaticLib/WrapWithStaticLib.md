# How to wrap fortran module which include static lib some
## solution 1
1. In testlib, we have three methods to do a matrix multiplication
    ```fortran
    matMult_func.f90
    matMult_loop1.f90
    matMult_loop2.f90 
    ```
    we fist bulid a static library
    ```bash
    ifort -c *.f90
    ar rc libmatMultiplication.a *.o
    ```
2. In mycode ,we have a module which call function in libmutiplication
    ```bash
    matrixMult.f90
    ```
    