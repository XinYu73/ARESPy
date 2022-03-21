from __future__ import print_function, absolute_import, division
import _Fib
import f90wrap.runtime
import logging

def fib(a, n):
    """
    fib(a, n)
    
    
    Defined at fib2.f90 lines 1-11
    
    Parameters
    ----------
    a : float
    n : int
    
    """
    _Fib.f90wrap_fib(a=a, n=n)

