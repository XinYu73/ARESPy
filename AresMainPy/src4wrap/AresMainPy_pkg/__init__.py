from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging
import AresMainPy_pkg.end_module
import AresMainPy_pkg.arpack_module
import AresMainPy_pkg.poisson_isf
import AresMainPy_pkg.math
import AresMainPy_pkg.energy_module
import AresMainPy_pkg.cg_relax
import AresMainPy_pkg.nlpot_module
import AresMainPy_pkg.band_structure
import AresMainPy_pkg.write_module
import AresMainPy_pkg.pspot_module
import AresMainPy_pkg.fourier
import AresMainPy_pkg.chebyshev_module
import AresMainPy_pkg.grid_module
import AresMainPy_pkg.begin_module
import AresMainPy_pkg.array_io
import AresMainPy_pkg.mathsplines
import AresMainPy_pkg.forcestress_module
import AresMainPy_pkg.constants
import AresMainPy_pkg.relax_module
import AresMainPy_pkg.libxc_module
import AresMainPy_pkg.lapack_module
import AresMainPy_pkg.getvlocalpseudopotential
import AresMainPy_pkg.dyn_module
import AresMainPy_pkg.parameters
import AresMainPy_pkg.potential_module
import AresMainPy_pkg.smpi_math_module
import AresMainPy_pkg.finite_module
import AresMainPy_pkg.matvec_module
import AresMainPy_pkg.struct_module
import AresMainPy_pkg.smearing_module
import AresMainPy_pkg.out_module
import AresMainPy_pkg.scalapack_module
import AresMainPy_pkg.ewald
import AresMainPy_pkg.scf_module
import AresMainPy_pkg.m_time_evaluate
import AresMainPy_pkg.isolateset
import AresMainPy_pkg.aresmainapi
import AresMainPy_pkg.mixer_module
import AresMainPy_pkg.succeed
import AresMainPy_pkg.read_module
import AresMainPy_pkg.sgfft_oct_m

def zbrent(iu, lreset, ebreak, x, y, f, xnew, xnewh, ynew, yd, ifail):
    """
    zbrent(iu, lreset, ebreak, x, y, f, xnew, xnewh, ynew, yd, ifail)
    
    
    Defined at brent.fpp lines 38-454
    
    Parameters
    ----------
    iu : int
    lreset : bool
    ebreak : float
    x : float
    y : float
    f : float
    xnew : float
    xnewh : float
    ynew : float
    yd : float
    ifail : int
    
    -----------------------------------------------------------------------
        case 0: trial step 1, into trial direction
    -----------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_zbrent(iu=iu, lreset=lreset, ebreak=ebreak, x=x, y=y, \
        f=f, xnew=xnew, xnewh=xnewh, ynew=ynew, yd=yd, ifail=ifail)

def kardir(nmax, v, basis):
    """
    kardir(nmax, v, basis)
    
    
    Defined at cg_vasp.fpp lines 10-22
    
    Parameters
    ----------
    nmax : int
    v : float array
    basis : float array
    
    """
    _AresMainPy_pkg.f90wrap_kardir(nmax=nmax, v=v, basis=basis)

def ioncgr(iflag, nions, toten, a, b, nfree, posion, posioc, fact, f, factsi, \
    fsif, fl, s, dismax, iu6, iu0, ebreak, ediffg, e1test, lstop2):
    """
    ioncgr(iflag, nions, toten, a, b, nfree, posion, posioc, fact, f, factsi, fsif, \
        fl, s, dismax, iu6, iu0, ebreak, ediffg, e1test, lstop2)
    
    
    Defined at cg_vasp.fpp lines 86-503
    
    Parameters
    ----------
    iflag : int
    nions : int
    toten : float
    a : float
    b : float
    nfree : int
    posion : float
    posioc : float
    fact : float
    f : float
    factsi : float
    fsif : float
    fl : float
    s : float
    dismax : float
    iu6 : int
    iu0 : int
    ebreak : float
    ediffg : float
    e1test : float
    lstop2 : bool
    
    =======================================================================
      if IFLAG =0 initialize everything
    =======================================================================
    """
    _AresMainPy_pkg.f90wrap_ioncgr(iflag=iflag, nions=nions, toten=toten, a=a, b=b, \
        nfree=nfree, posion=posion, posioc=posioc, fact=fact, f=f, factsi=factsi, \
        fsif=fsif, fl=fl, s=s, dismax=dismax, iu6=iu6, iu0=iu0, ebreak=ebreak, \
        ediffg=ediffg, e1test=e1test, lstop2=lstop2)

