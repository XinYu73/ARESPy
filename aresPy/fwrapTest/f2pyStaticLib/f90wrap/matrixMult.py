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

