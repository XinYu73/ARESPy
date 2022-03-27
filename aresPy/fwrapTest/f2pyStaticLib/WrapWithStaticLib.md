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

## Solution 2 with f90wrap (it's a little bit tricky here)

1. build the library as solution 1.
2. build outter layer with f90wrap

    ```bash
    f90wrap -m matrixMult -k kind_map  matrixMult.f90
    ```

    however in this case, f90wrap_toplevel.f90 is a mess, then we delet irrelevant vars

    ```fortran
    subroutine f90wrap_matrixmult(a, b, c, n1, n2, n3, icase)
    implicit none
    external matrixmult
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: icase
    double precision, intent(in), dimension(n1,n2) :: a
    double precision, intent(in), dimension(n2,n3) :: b
    double precision, intent(out), dimension(n1,n3) :: c
    call matrixmult(a, b, c, n1, n2, n3, icase)
    end subroutine f90wrap_matrixmult
    ```

3. bulid python lib

    ```bash
    f2py-f90wrap --f90exec=ifort -m _matrixMult -L../testlib -lmatMultiplication -c f90wrap_toplevel.f90 matrixMult.f90
    ```

    yes, we just change the f2py into f2py-f90wrap and add f90wrap_toplevel.f90

    again, the python lib is a mess, wrong input vars , not return value, modify a little

    ```python
    from __future__ import print_function, absolute_import, division
    import _matrixMult
    import f90wrap.runtime
    import logging

    def matrixmult(a, b, n1, n2, n3, icase):
    """
    matrixmult(a, b, c, n1, n2, n3, icase)
    
    
    Defined at matrixMult.f90 lines 1-14
    
    Parameters
    ----------
    a : unknown array
    b : unknown array
    c : unknown array
    n1 : int
    n2 : int
    n3 : int
    icase : int
    
    """
    return _matrixMult.f90wrap_matrixmult(a=a, b=b, n1=n1, n2=n2, n3=n3, icase=icase)
    ```

    we change delete input var c and add return before _matrixMult.f90wrap_matrixmult

4. test

    to actually run the test , we add more input vars in matrixmult, anyway, the result is

    ```bash
    [[0.49783703 0.76295641 0.95495111 1.38922937 0.70565897]
    [0.16463549 0.24609853 0.30920616 0.45152313 0.22792956]
    [0.1512668  0.12061612 0.17205578 0.28074758 0.11715796]
    [0.44023201 0.64986678 0.81810698 1.19694528 0.60231151]]
    Loop1: time for AB(4, 5) = A(4, 2) B(2, 5) is 0.0014553070068359375 s
    [[0.49783703 0.76295641 0.95495111 1.38922937 0.70565897]
    [0.16463549 0.24609853 0.30920616 0.45152313 0.22792956]
    [0.1512668  0.12061612 0.17205578 0.28074758 0.11715796]
    [0.44023201 0.64986678 0.81810698 1.19694528 0.60231151]]
    Loop2: time for AB(4, 5) = A(4, 2) B(2, 5) is 0.0007596015930175781 s
    [[0.49783703 0.76295641 0.95495111 1.38922937 0.70565897]
    [0.16463549 0.24609853 0.30920616 0.45152313 0.22792956]
    [0.1512668  0.12061612 0.17205578 0.28074758 0.11715796]
    [0.44023201 0.64986678 0.81810698 1.19694528 0.60231151]]
    matmul function: time for AB(4, 5) = A(4, 2) B(2, 5) is 0.0006937980651855469 s
    ```
