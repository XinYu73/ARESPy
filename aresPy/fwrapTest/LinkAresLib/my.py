from __future__ import print_function, absolute_import, division
import _my
import f90wrap.runtime
import logging

def matrixmult(n, a, d):
    """
    matrixmult(n, a, d)
    
    
    Defined at my.f90 lines 1-7
    
    Parameters
    ----------
    n : int
    a : unknown array
    d : unknown array
    
    """
    _my.f90wrap_matrixmult(n=n, a=a, d=d)

