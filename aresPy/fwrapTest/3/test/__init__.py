from __future__ import print_function, absolute_import, division
import _test
import f90wrap.runtime
import logging

def sayhello(comm):
    """
    sayhello(comm)
    
    
    Defined at helloworld.f90 lines 2-7
    
    Parameters
    ----------
    comm : int
    
    """
    _test.f90wrap_sayhello(comm=comm)

