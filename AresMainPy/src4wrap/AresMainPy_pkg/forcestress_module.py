"""
Module forcestress_module


Defined at ForceStress_module.fpp lines 10-2319

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cal_force_c(rhos, uik, force):
    """
    cal_force_c(rhos, uik, force)
    
    
    Defined at ForceStress_module.fpp lines 20-84
    
    Parameters
    ----------
    rhos : float array
    uik : complex array
    force : float array
    
    ===================================
    nlnlocal part
    """
    _AresMainPy_pkg.f90wrap_cal_force_c(rhos=rhos, uik=uik, force=force)

def cal_force_gamma(rhos, uik, force):
    """
    cal_force_gamma(rhos, uik, force)
    
    
    Defined at ForceStress_module.fpp lines 87-147
    
    Parameters
    ----------
    rhos : float array
    uik : float array
    force : float array
    
    ===================================
    nlnlocal part
    """
    _AresMainPy_pkg.f90wrap_cal_force_gamma(rhos=rhos, uik=uik, force=force)

def cal_force_r(rhos, uik, force):
    """
    cal_force_r(rhos, uik, force)
    
    
    Defined at ForceStress_module.fpp lines 150-241
    
    Parameters
    ----------
    rhos : float array
    uik : float array
    force : float array
    
    ===================================
    nlnlocal part
    """
    _AresMainPy_pkg.f90wrap_cal_force_r(rhos=rhos, uik=uik, force=force)

def locforce(rhos, lforce):
    """
    locforce(rhos, lforce)
    
    
    Defined at ForceStress_module.fpp lines 245-377
    
    Parameters
    ----------
    rhos : float array
    lforce : float array
    
    ===========================
    #PLOT
    open(23301,file="frxyz_0.2")
    write(23301,*)psp(1)%VlocqS
    close(23301)
    open(23302,file="rxyz_0.2")
    print*,"qmax/qspacing",psp(1)%qmax/psp(1)%qspacing,nint(psp(1)%qmax/psp(1)%qspacing),int(psp(1)%qmax/psp(1)%qspacing)
    DO I1=0,nint(psp(1)%qmax/psp(1)%qspacing)
        WRITE(23302,*)I1*psp(1)%qspacing
    ENDDO
    WRITE(23302,*)psp(1)%r_real
    close(23302)
    open(23303,file="rxyz_0.3")
    open(23304,file="frxyz_0.3")
    ===========================
    """
    _AresMainPy_pkg.f90wrap_locforce(rhos=rhos, lforce=lforce)

def locforce_r(rhos, lforce):
    """
    locforce_r(rhos, lforce)
    
    
    Defined at ForceStress_module.fpp lines 380-488
    
    Parameters
    ----------
    rhos : float array
    lforce : float array
    
    ==================================================
    > initial parallel config
     IF(.not.allocated(atom_index))ALLOCATE(atom_index(natom/parallel%numprocs+1))
     CALL start_time('init_density_0')
     CALL atom_split(mysize,atom_index)
     print *,'atom_index',atom_index,'id',parallel%myid
     CALL end_time('init_density_0')
     CALL write_time('init_density_0')
     id_core=1
    >=================================================
    """
    _AresMainPy_pkg.f90wrap_locforce_r(rhos=rhos, lforce=lforce)

def nonlforce(uik, nlforce):
    """
    nonlforce(uik, nlforce)
    
    
    Defined at ForceStress_module.fpp lines 597-728
    
    Parameters
    ----------
    uik : complex array
    nlforce : float array
    
    ---------------------------------
    """
    _AresMainPy_pkg.f90wrap_nonlforce(uik=uik, nlforce=nlforce)

def nonlforce_gamma(uik, nlforce):
    """
    nonlforce_gamma(uik, nlforce)
    
    
    Defined at ForceStress_module.fpp lines 731-861
    
    Parameters
    ----------
    uik : float array
    nlforce : float array
    
    ---------------------------------
    """
    _AresMainPy_pkg.f90wrap_nonlforce_gamma(uik=uik, nlforce=nlforce)

def nonlforce_r_dg(uik, nlforce):
    """
    nonlforce_r_dg(uik, nlforce)
    
    
    Defined at ForceStress_module.fpp lines 865-993
    
    Parameters
    ----------
    uik : float array
    nlforce : float array
    
    ---------------------------------
    """
    _AresMainPy_pkg.f90wrap_nonlforce_r_dg(uik=uik, nlforce=nlforce)

def nonlforce_r(uik, nlforce):
    """
    nonlforce_r(uik, nlforce)
    
    
    Defined at ForceStress_module.fpp lines 996-1127
    
    Parameters
    ----------
    uik : float array
    nlforce : float array
    
    ---------------------------------
    """
    _AresMainPy_pkg.f90wrap_nonlforce_r(uik=uik, nlforce=nlforce)

def cal_stress(rhos, uik, stress):
    """
    cal_stress(rhos, uik, stress)
    
    
    Defined at ForceStress_module.fpp lines 1131-1234
    
    Parameters
    ----------
    rhos : float array
    uik : complex array
    stress : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_stress(rhos=rhos, uik=uik, stress=stress)

def cal_stress_gamma(rhos, uik, stress):
    """
    cal_stress_gamma(rhos, uik, stress)
    
    
    Defined at ForceStress_module.fpp lines 1237-1340
    
    Parameters
    ----------
    rhos : float array
    uik : float array
    stress : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_stress_gamma(rhos=rhos, uik=uik, stress=stress)

def lda_stress(vxc, rho, elda):
    """
    lda_stress = lda_stress(vxc, rho, elda)
    
    
    Defined at ForceStress_module.fpp lines 1343-1375
    
    Parameters
    ----------
    vxc : float array
    rho : float array
    elda : float
    
    Returns
    -------
    lda_stress : float array
    
    """
    lda_stress = _AresMainPy_pkg.f90wrap_lda_stress(vxc=vxc, rho=rho, elda=elda)
    return lda_stress

def hart_stress(rhorecip, eh):
    """
    hart_stress = hart_stress(rhorecip, eh)
    
    
    Defined at ForceStress_module.fpp lines 1378-1423
    
    Parameters
    ----------
    rhorecip : complex array
    eh : float
    
    Returns
    -------
    hart_stress : float array
    
    """
    hart_stress = _AresMainPy_pkg.f90wrap_hart_stress(rhorecip=rhorecip, eh=eh)
    return hart_stress

def kin_stress(uik, kinstress):
    """
    kin_stress(uik, kinstress)
    
    
    Defined at ForceStress_module.fpp lines 1426-1565
    
    Parameters
    ----------
    uik : complex array
    kinstress : float array
    
    """
    _AresMainPy_pkg.f90wrap_kin_stress(uik=uik, kinstress=kinstress)

def kin_stress_gamma(uik, kinstress):
    """
    kin_stress_gamma(uik, kinstress)
    
    
    Defined at ForceStress_module.fpp lines 1567-1706
    
    Parameters
    ----------
    uik : float array
    kinstress : float array
    
    """
    _AresMainPy_pkg.f90wrap_kin_stress_gamma(uik=uik, kinstress=kinstress)

def ion_nl_stress(uik, nlstress):
    """
    ion_nl_stress(uik, nlstress)
    
    
    Defined at ForceStress_module.fpp lines 1709-1885
    
    Parameters
    ----------
    uik : complex array
    nlstress : float array
    
    ----------------------------------
    """
    _AresMainPy_pkg.f90wrap_ion_nl_stress(uik=uik, nlstress=nlstress)

def ion_nl_stress_gamma(uik, nlstress):
    """
    ion_nl_stress_gamma(uik, nlstress)
    
    
    Defined at ForceStress_module.fpp lines 1887-2062
    
    Parameters
    ----------
    uik : float array
    nlstress : float array
    
    ----------------------------------
    """
    _AresMainPy_pkg.f90wrap_ion_nl_stress_gamma(uik=uik, nlstress=nlstress)

def ionele_stress(rhorecip, energy):
    """
    ionele_stress = ionele_stress(rhorecip, energy)
    
    
    Defined at ForceStress_module.fpp lines 2065-2145
    
    Parameters
    ----------
    rhorecip : complex array
    energy : float
    
    Returns
    -------
    ionele_stress : float array
    
    """
    ionele_stress = _AresMainPy_pkg.f90wrap_ionele_stress(rhorecip=rhorecip, \
        energy=energy)
    return ionele_stress

def pseudopotdifflookup(ity, qnorm):
    """
    pseudopotdifflookup = pseudopotdifflookup(ity, qnorm)
    
    
    Defined at ForceStress_module.fpp lines 2149-2188
    
    Parameters
    ----------
    ity : int
    qnorm : float
    
    Returns
    -------
    pseudopotdifflookup : float
    
    """
    pseudopotdifflookup = _AresMainPy_pkg.f90wrap_pseudopotdifflookup(ity=ity, \
        qnorm=qnorm)
    return pseudopotdifflookup

def cal_force_stress():
    """
    cal_force_stress()
    
    
    Defined at ForceStress_module.fpp lines 2192-2318
    
    
    ====================
    ##cal force directly
    goto 10011
    ====================
    Self-consistent
    """
    _AresMainPy_pkg.f90wrap_cal_force_stress()

def get_cellpress():
    """
    Element cellpress ftype=real(dp) pytype=float
    
    
    Defined at ForceStress_module.fpp line 13
    
    """
    return _AresMainPy_pkg.f90wrap_forcestress_module__get__cellpress()

def set_cellpress(cellpress):
    _AresMainPy_pkg.f90wrap_forcestress_module__set__cellpress(cellpress)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "forcestress_module".')

for func in _dt_array_initialisers:
    func()
