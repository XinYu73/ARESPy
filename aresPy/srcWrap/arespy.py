from __future__ import print_function, absolute_import, division
import _arespy
import f90wrap.runtime
import logging

class Arpack_Module(f90wrap.runtime.FortranModule):
    """
    Module arpack_module
    
    
    Defined at Arpack_module.f90 lines 1-189
    
    """
    @staticmethod
    def real_diagh_arpack(n, veff, nev, evec, eval, resid_restart, nec, info, \
        maxmvs, tol):
        """
        real_diagh_arpack(n, veff, nev, evec, eval, resid_restart, nec, info, maxmvs, \
            tol)
        
        
        Defined at Arpack_module.f90 lines 36-188
        
        Parameters
        ----------
        n : int
        veff : float array
        nev : int
        evec : float array
        eval : float array
        resid_restart : float array
        nec : int
        info : int
        maxmvs : int
        tol : float
        
        ------
        """
        _arespy.f90wrap_real_diagh_arpack(n=n, veff=veff, nev=nev, evec=evec, eval=eval, \
            resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, tol=tol)
    
    @property
    def maxn(self):
        """
        Element maxn ftype=integer(i4b) pytype=int
        
        
        Defined at Arpack_module.f90 line 9
        
        """
        return _arespy.f90wrap_arpack_module__get__maxn()
    
    @maxn.setter
    def maxn(self, maxn):
        _arespy.f90wrap_arpack_module__set__maxn(maxn)
    
    @property
    def maxnev(self):
        """
        Element maxnev ftype=integer(i4b) pytype=int
        
        
        Defined at Arpack_module.f90 line 9
        
        """
        return _arespy.f90wrap_arpack_module__get__maxnev()
    
    @maxnev.setter
    def maxnev(self, maxnev):
        _arespy.f90wrap_arpack_module__set__maxnev(maxnev)
    
    @property
    def maxncv(self):
        """
        Element maxncv ftype=integer(i4b) pytype=int
        
        
        Defined at Arpack_module.f90 line 9
        
        """
        return _arespy.f90wrap_arpack_module__get__maxncv()
    
    @maxncv.setter
    def maxncv(self, maxncv):
        _arespy.f90wrap_arpack_module__set__maxncv(maxncv)
    
    @property
    def logfil(self):
        """
        Element logfil ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__logfil()
    
    @logfil.setter
    def logfil(self, logfil):
        _arespy.f90wrap_arpack_module__set__logfil(logfil)
    
    @property
    def ndigit(self):
        """
        Element ndigit ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__ndigit()
    
    @ndigit.setter
    def ndigit(self, ndigit):
        _arespy.f90wrap_arpack_module__set__ndigit(ndigit)
    
    @property
    def mgetv0(self):
        """
        Element mgetv0 ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mgetv0()
    
    @mgetv0.setter
    def mgetv0(self, mgetv0):
        _arespy.f90wrap_arpack_module__set__mgetv0(mgetv0)
    
    @property
    def msaupd(self):
        """
        Element msaupd ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__msaupd()
    
    @msaupd.setter
    def msaupd(self, msaupd):
        _arespy.f90wrap_arpack_module__set__msaupd(msaupd)
    
    @property
    def msaup2(self):
        """
        Element msaup2 ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__msaup2()
    
    @msaup2.setter
    def msaup2(self, msaup2):
        _arespy.f90wrap_arpack_module__set__msaup2(msaup2)
    
    @property
    def msaitr(self):
        """
        Element msaitr ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__msaitr()
    
    @msaitr.setter
    def msaitr(self, msaitr):
        _arespy.f90wrap_arpack_module__set__msaitr(msaitr)
    
    @property
    def mseigt(self):
        """
        Element mseigt ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mseigt()
    
    @mseigt.setter
    def mseigt(self, mseigt):
        _arespy.f90wrap_arpack_module__set__mseigt(mseigt)
    
    @property
    def msapps(self):
        """
        Element msapps ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__msapps()
    
    @msapps.setter
    def msapps(self, msapps):
        _arespy.f90wrap_arpack_module__set__msapps(msapps)
    
    @property
    def msgets(self):
        """
        Element msgets ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__msgets()
    
    @msgets.setter
    def msgets(self, msgets):
        _arespy.f90wrap_arpack_module__set__msgets(msgets)
    
    @property
    def mseupd(self):
        """
        Element mseupd ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mseupd()
    
    @mseupd.setter
    def mseupd(self, mseupd):
        _arespy.f90wrap_arpack_module__set__mseupd(mseupd)
    
    @property
    def mnaupd(self):
        """
        Element mnaupd ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mnaupd()
    
    @mnaupd.setter
    def mnaupd(self, mnaupd):
        _arespy.f90wrap_arpack_module__set__mnaupd(mnaupd)
    
    @property
    def mnaup2(self):
        """
        Element mnaup2 ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mnaup2()
    
    @mnaup2.setter
    def mnaup2(self, mnaup2):
        _arespy.f90wrap_arpack_module__set__mnaup2(mnaup2)
    
    @property
    def mnaitr(self):
        """
        Element mnaitr ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mnaitr()
    
    @mnaitr.setter
    def mnaitr(self, mnaitr):
        _arespy.f90wrap_arpack_module__set__mnaitr(mnaitr)
    
    @property
    def mneigh(self):
        """
        Element mneigh ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mneigh()
    
    @mneigh.setter
    def mneigh(self, mneigh):
        _arespy.f90wrap_arpack_module__set__mneigh(mneigh)
    
    @property
    def mnapps(self):
        """
        Element mnapps ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mnapps()
    
    @mnapps.setter
    def mnapps(self, mnapps):
        _arespy.f90wrap_arpack_module__set__mnapps(mnapps)
    
    @property
    def mngets(self):
        """
        Element mngets ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mngets()
    
    @mngets.setter
    def mngets(self, mngets):
        _arespy.f90wrap_arpack_module__set__mngets(mngets)
    
    @property
    def mneupd(self):
        """
        Element mneupd ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mneupd()
    
    @mneupd.setter
    def mneupd(self, mneupd):
        _arespy.f90wrap_arpack_module__set__mneupd(mneupd)
    
    @property
    def mcaupd(self):
        """
        Element mcaupd ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mcaupd()
    
    @mcaupd.setter
    def mcaupd(self, mcaupd):
        _arespy.f90wrap_arpack_module__set__mcaupd(mcaupd)
    
    @property
    def mcaup2(self):
        """
        Element mcaup2 ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mcaup2()
    
    @mcaup2.setter
    def mcaup2(self, mcaup2):
        _arespy.f90wrap_arpack_module__set__mcaup2(mcaup2)
    
    @property
    def mcaitr(self):
        """
        Element mcaitr ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mcaitr()
    
    @mcaitr.setter
    def mcaitr(self, mcaitr):
        _arespy.f90wrap_arpack_module__set__mcaitr(mcaitr)
    
    @property
    def mceigh(self):
        """
        Element mceigh ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mceigh()
    
    @mceigh.setter
    def mceigh(self, mceigh):
        _arespy.f90wrap_arpack_module__set__mceigh(mceigh)
    
    @property
    def mcapps(self):
        """
        Element mcapps ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mcapps()
    
    @mcapps.setter
    def mcapps(self, mcapps):
        _arespy.f90wrap_arpack_module__set__mcapps(mcapps)
    
    @property
    def mcgets(self):
        """
        Element mcgets ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mcgets()
    
    @mcgets.setter
    def mcgets(self, mcgets):
        _arespy.f90wrap_arpack_module__set__mcgets(mcgets)
    
    @property
    def mceupd(self):
        """
        Element mceupd ftype=integer  pytype=int
        
        
        Defined at Arpack_module.f90 line 21
        
        """
        return _arespy.f90wrap_arpack_module__get__mceupd()
    
    @mceupd.setter
    def mceupd(self, mceupd):
        _arespy.f90wrap_arpack_module__set__mceupd(mceupd)
    
    def __str__(self):
        ret = ['<arpack_module>{\n']
        ret.append('    maxn : ')
        ret.append(repr(self.maxn))
        ret.append(',\n    maxnev : ')
        ret.append(repr(self.maxnev))
        ret.append(',\n    maxncv : ')
        ret.append(repr(self.maxncv))
        ret.append(',\n    logfil : ')
        ret.append(repr(self.logfil))
        ret.append(',\n    ndigit : ')
        ret.append(repr(self.ndigit))
        ret.append(',\n    mgetv0 : ')
        ret.append(repr(self.mgetv0))
        ret.append(',\n    msaupd : ')
        ret.append(repr(self.msaupd))
        ret.append(',\n    msaup2 : ')
        ret.append(repr(self.msaup2))
        ret.append(',\n    msaitr : ')
        ret.append(repr(self.msaitr))
        ret.append(',\n    mseigt : ')
        ret.append(repr(self.mseigt))
        ret.append(',\n    msapps : ')
        ret.append(repr(self.msapps))
        ret.append(',\n    msgets : ')
        ret.append(repr(self.msgets))
        ret.append(',\n    mseupd : ')
        ret.append(repr(self.mseupd))
        ret.append(',\n    mnaupd : ')
        ret.append(repr(self.mnaupd))
        ret.append(',\n    mnaup2 : ')
        ret.append(repr(self.mnaup2))
        ret.append(',\n    mnaitr : ')
        ret.append(repr(self.mnaitr))
        ret.append(',\n    mneigh : ')
        ret.append(repr(self.mneigh))
        ret.append(',\n    mnapps : ')
        ret.append(repr(self.mnapps))
        ret.append(',\n    mngets : ')
        ret.append(repr(self.mngets))
        ret.append(',\n    mneupd : ')
        ret.append(repr(self.mneupd))
        ret.append(',\n    mcaupd : ')
        ret.append(repr(self.mcaupd))
        ret.append(',\n    mcaup2 : ')
        ret.append(repr(self.mcaup2))
        ret.append(',\n    mcaitr : ')
        ret.append(repr(self.mcaitr))
        ret.append(',\n    mceigh : ')
        ret.append(repr(self.mceigh))
        ret.append(',\n    mcapps : ')
        ret.append(repr(self.mcapps))
        ret.append(',\n    mcgets : ')
        ret.append(repr(self.mcgets))
        ret.append(',\n    mceupd : ')
        ret.append(repr(self.mceupd))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

arpack_module = Arpack_Module()

class Band_Structure(f90wrap.runtime.FortranModule):
    """
    Module band_structure
    
    
    Defined at Bands_module.f90 lines 6-547
    
    """
    @staticmethod
    def init_bandstruct(numk, kvec):
        """
        init_bandstruct(numk, kvec)
        
        
        Defined at Bands_module.f90 lines 25-45
        
        Parameters
        ----------
        numk : int
        kvec : float array
        
        """
        _arespy.f90wrap_init_bandstruct(numk=numk, kvec=kvec)
    
    @staticmethod
    def band_begin(numk, kvec):
        """
        band_begin(numk, kvec)
        
        
        Defined at Bands_module.f90 lines 48-123
        
        Parameters
        ----------
        numk : int
        kvec : float array
        
        """
        _arespy.f90wrap_band_begin(numk=numk, kvec=kvec)
    
    @staticmethod
    def build_bands():
        """
        build_bands()
        
        
        Defined at Bands_module.f90 lines 126-229
        
        
        """
        _arespy.f90wrap_build_bands()
    
    @staticmethod
    def read_bandpath(infile):
        """
        read_bandpath(infile)
        
        
        Defined at Bands_module.f90 lines 232-352
        
        Parameters
        ----------
        infile : str
        
        -------start input.dat-----------
        """
        _arespy.f90wrap_read_bandpath(infile=infile)
    
    @staticmethod
    def band_pathread(infile):
        """
        band_pathread(infile)
        
        
        Defined at Bands_module.f90 lines 355-398
        
        Parameters
        ----------
        infile : str
        
        """
        _arespy.f90wrap_band_pathread(infile=infile)
    
    @staticmethod
    def read_density(infile, rho):
        """
        read_density(infile, rho)
        
        
        Defined at Bands_module.f90 lines 401-472
        
        Parameters
        ----------
        infile : str
        rho : float array
        
        """
        _arespy.f90wrap_read_density(infile=infile, rho=rho)
    
    @staticmethod
    def cal_band(rhos, nev, eigval):
        """
        cal_band(rhos, nev, eigval)
        
        
        Defined at Bands_module.f90 lines 475-546
        
        Parameters
        ----------
        rhos : float array
        nev : int
        eigval : float array
        
        """
        _arespy.f90wrap_cal_band(rhos=rhos, nev=nev, eigval=eigval)
    
    _dt_array_initialisers = []
    

band_structure = Band_Structure()

class Begin_Module(f90wrap.runtime.FortranModule):
    """
    Module begin_module
    
    
    Defined at Begin_module.f90 lines 1-252
    
    """
    @staticmethod
    def initial_grid_pbc():
        """
        initial_grid_pbc()
        
        
        Defined at Begin_module.f90 lines 6-89
        
        
        """
        _arespy.f90wrap_initial_grid_pbc()
    
    @staticmethod
    def initial_density():
        """
        initial_density()
        
        
        Defined at Begin_module.f90 lines 92-182
        
        
        """
        _arespy.f90wrap_initial_density()
    
    @staticmethod
    def inichrg_sp():
        """
        inichrg_sp()
        
        
        Defined at Begin_module.f90 lines 185-251
        
        
        """
        _arespy.f90wrap_inichrg_sp()
    
    _dt_array_initialisers = []
    

begin_module = Begin_Module()

class Chebyshev_Module(f90wrap.runtime.FortranModule):
    """
    Module chebyshev_module
    
    
    Defined at Chebyshev_fliter.f90 lines 6-926
    
    """
    @staticmethod
    def buildsubspace(nps, nev, veff, eig):
        """
        buildsubspace(nps, nev, veff, eig)
        
        
        Defined at Chebyshev_fliter.f90 lines 20-129
        
        Parameters
        ----------
        nps : int
        nev : int
        veff : float array
        eig : Eigen_Type
        
        """
        _arespy.f90wrap_buildsubspace(nps=nps, nev=nev, veff=veff, eig=eig._handle)
    
    @staticmethod
    def real_pseudosubspace(nps, nev, initx):
        """
        real_pseudosubspace(nps, nev, initx)
        
        
        Defined at Chebyshev_fliter.f90 lines 132-219
        
        Parameters
        ----------
        nps : int
        nev : int
        initx : float array
        
        """
        _arespy.f90wrap_real_pseudosubspace(nps=nps, nev=nev, initx=initx)
    
    @staticmethod
    def real_first_rrstep(nps, nev, veff, x, d):
        """
        real_first_rrstep(nps, nev, veff, x, d)
        
        
        Defined at Chebyshev_fliter.f90 lines 222-316
        
        Parameters
        ----------
        nps : int
        nev : int
        veff : float array
        x : float array
        d : float array
        
        -------------------
        rotation
        """
        _arespy.f90wrap_real_first_rrstep(nps=nps, nev=nev, veff=veff, x=x, d=d)
    
    @staticmethod
    def init_uplow_real(nps, k, veff, v):
        """
        a, b, al = init_uplow_real(nps, k, veff, v)
        
        
        Defined at Chebyshev_fliter.f90 lines 319-395
        
        Parameters
        ----------
        nps : int
        k : int
        veff : float array
        v : float array
        
        Returns
        -------
        a : float
        b : float
        al : float
        
        """
        a, b, al = _arespy.f90wrap_init_uplow_real(nps=nps, k=k, veff=veff, v=v)
        return a, b, al
    
    @staticmethod
    def real_first_filter(nps, nst, veff, x, eval):
        """
        real_first_filter(nps, nst, veff, x, eval)
        
        
        Defined at Chebyshev_fliter.f90 lines 398-524
        
        Parameters
        ----------
        nps : int
        nst : int
        veff : float array
        x : float array
        eval : float array
        
        -----------------
        """
        _arespy.f90wrap_real_first_filter(nps=nps, nst=nst, veff=veff, x=x, eval=eval)
    
    @staticmethod
    def rayleigh_quotient_real(nps, nst, veff, x, xhx):
        """
        rayleigh_quotient_real(nps, nst, veff, x, xhx)
        
        
        Defined at Chebyshev_fliter.f90 lines 532-561
        
        Parameters
        ----------
        nps : int
        nst : int
        veff : float array
        x : float array
        xhx : float array
        
        """
        _arespy.f90wrap_rayleigh_quotient_real(nps=nps, nst=nst, veff=veff, x=x, \
            xhx=xhx)
    
    @staticmethod
    def cal_hx_real(nps, nst, veff, v, hv):
        """
        cal_hx_real(nps, nst, veff, v, hv)
        
        
        Defined at Chebyshev_fliter.f90 lines 564-578
        
        Parameters
        ----------
        nps : int
        nst : int
        veff : float array
        v : float array
        hv : float array
        
        """
        _arespy.f90wrap_cal_hx_real(nps=nps, nst=nst, veff=veff, v=v, hv=hv)
    
    @staticmethod
    def estupb_real(nps, k, veff, vec):
        """
        b = estupb_real(nps, k, veff, vec)
        
        
        Defined at Chebyshev_fliter.f90 lines 581-656
        
        Parameters
        ----------
        nps : int
        k : int
        veff : float array
        vec : float array
        
        Returns
        -------
        b : float
        
        """
        b = _arespy.f90wrap_estupb_real(nps=nps, k=k, veff=veff, vec=vec)
        return b
    
    @staticmethod
    def chebyshev_filter_real(nps, nst, veff, x, m, a, b):
        """
        chebyshev_filter_real(nps, nst, veff, x, m, a, b)
        
        
        Defined at Chebyshev_fliter.f90 lines 659-693
        
        Parameters
        ----------
        nps : int
        nst : int
        veff : float array
        x : float array
        m : int
        a : float
        b : float
        
        """
        _arespy.f90wrap_chebyshev_filter_real(nps=nps, nst=nst, veff=veff, x=x, m=m, \
            a=a, b=b)
    
    @staticmethod
    def chebyshev_filter_scaled_real(nps, nst, veff, x, m, a, b, al):
        """
        chebyshev_filter_scaled_real(nps, nst, veff, x, m, a, b, al)
        
        
        Defined at Chebyshev_fliter.f90 lines 696-734
        
        Parameters
        ----------
        nps : int
        nst : int
        veff : float array
        x : float array
        m : int
        a : float
        b : float
        al : float
        
        """
        _arespy.f90wrap_chebyshev_filter_scaled_real(nps=nps, nst=nst, veff=veff, x=x, \
            m=m, a=a, b=b, al=al)
    
    @staticmethod
    def grayleigh_ritz_real(nps, nev, veff, x, d):
        """
        grayleigh_ritz_real(nps, nev, veff, x, d)
        
        
        Defined at Chebyshev_fliter.f90 lines 737-806
        
        Parameters
        ----------
        nps : int
        nev : int
        veff : float array
        x : float array
        d : float array
        
        """
        _arespy.f90wrap_grayleigh_ritz_real(nps=nps, nev=nev, veff=veff, x=x, d=d)
    
    @staticmethod
    def rayleigh_ritz_real(nps, sn, veff, x, d):
        """
        rayleigh_ritz_real(nps, sn, veff, x, d)
        
        
        Defined at Chebyshev_fliter.f90 lines 809-873
        
        Parameters
        ----------
        nps : int
        sn : int
        veff : float array
        x : float array
        d : float array
        
        """
        _arespy.f90wrap_rayleigh_ritz_real(nps=nps, sn=sn, veff=veff, x=x, d=d)
    
    @staticmethod
    def cheby_filtering_grrr(nps, nev, veff, x, d):
        """
        cheby_filtering_grrr(nps, nev, veff, x, d)
        
        
        Defined at Chebyshev_fliter.f90 lines 876-925
        
        Parameters
        ----------
        nps : int
        nev : int
        veff : float array
        x : float array
        d : float array
        
        """
        _arespy.f90wrap_cheby_filtering_grrr(nps=nps, nev=nev, veff=veff, x=x, d=d)
    
    @property
    def larged(self):
        """
        Element larged ftype=real(dp) pytype=float
        
        
        Defined at Chebyshev_fliter.f90 line 14
        
        """
        return _arespy.f90wrap_chebyshev_module__get__larged()
    
    @property
    def ad(self):
        """
        Element ad ftype=real(dp) pytype=float
        
        
        Defined at Chebyshev_fliter.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_chebyshev_module__array__ad(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ad = self._arrays[array_handle]
        else:
            ad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_chebyshev_module__array__ad)
            self._arrays[array_handle] = ad
        return ad
    
    @ad.setter
    def ad(self, ad):
        self.ad[...] = ad
    
    def __str__(self):
        ret = ['<chebyshev_module>{\n']
        ret.append('    larged : ')
        ret.append(repr(self.larged))
        ret.append(',\n    ad : ')
        ret.append(repr(self.ad))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

chebyshev_module = Chebyshev_Module()

class Constants(f90wrap.runtime.FortranModule):
    """
    Module constants
    
    
    Defined at Constants.f90 lines 1-93
    
    """
    @property
    def i4b(self):
        """
        Element i4b ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 21
        
        """
        return _arespy.f90wrap_constants__get__i4b()
    
    @property
    def i2b(self):
        """
        Element i2b ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 22
        
        """
        return _arespy.f90wrap_constants__get__i2b()
    
    @property
    def i1b(self):
        """
        Element i1b ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 23
        
        """
        return _arespy.f90wrap_constants__get__i1b()
    
    @property
    def sp(self):
        """
        Element sp ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 24
        
        """
        return _arespy.f90wrap_constants__get__sp()
    
    @property
    def dp(self):
        """
        Element dp ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 25
        
        """
        return _arespy.f90wrap_constants__get__dp()
    
    @property
    def scp(self):
        """
        Element scp ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 26
        
        """
        return _arespy.f90wrap_constants__get__scp()
    
    @property
    def dcp(self):
        """
        Element dcp ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 27
        
        """
        return _arespy.f90wrap_constants__get__dcp()
    
    @property
    def lgt(self):
        """
        Element lgt ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 29
        
        """
        return _arespy.f90wrap_constants__get__lgt()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 32
        
        """
        return _arespy.f90wrap_constants__get__pi()
    
    @property
    def pio2(self):
        """
        Element pio2 ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 33
        
        """
        return _arespy.f90wrap_constants__get__pio2()
    
    @property
    def twopi(self):
        """
        Element twopi ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 34
        
        """
        return _arespy.f90wrap_constants__get__twopi()
    
    @property
    def imag(self):
        """
        Element imag ftype=complex(dcp) pytype=complex
        
        
        Defined at Constants.f90 line 35
        
        """
        return _arespy.f90wrap_constants__get__imag()
    
    @property
    def inputunit(self):
        """
        Element inputunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.f90 line 36
        
        """
        return _arespy.f90wrap_constants__get__inputunit()
    
    @property
    def errunit(self):
        """
        Element errunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.f90 line 37
        
        """
        return _arespy.f90wrap_constants__get__errunit()
    
    @property
    def outputunit(self):
        """
        Element outputunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.f90 line 38
        
        """
        return _arespy.f90wrap_constants__get__outputunit()
    
    @property
    def mdposunit(self):
        """
        Element mdposunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.f90 line 39
        
        """
        return _arespy.f90wrap_constants__get__mdposunit()
    
    @property
    def rydberg(self):
        """
        Element rydberg ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 40
        
        """
        return _arespy.f90wrap_constants__get__rydberg()
    
    @property
    def bohr2ang(self):
        """
        Element bohr2ang ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 41
        
        """
        return _arespy.f90wrap_constants__get__bohr2ang()
    
    @property
    def ang2bohr(self):
        """
        Element ang2bohr ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 42
        
        """
        return _arespy.f90wrap_constants__get__ang2bohr()
    
    @property
    def hart2ev(self):
        """
        Element hart2ev ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 43
        
        """
        return _arespy.f90wrap_constants__get__hart2ev()
    
    @property
    def ev2hart(self):
        """
        Element ev2hart ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 44
        
        """
        return _arespy.f90wrap_constants__get__ev2hart()
    
    @property
    def scf_tol(self):
        """
        Element scf_tol ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 45
        
        """
        return _arespy.f90wrap_constants__get__scf_tol()
    
    @property
    def au2ev_force(self):
        """
        Element au2ev_force ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 46
        
        """
        return _arespy.f90wrap_constants__get__au2ev_force()
    
    @property
    def au2gpa(self):
        """
        Element au2gpa ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 47
        
        """
        return _arespy.f90wrap_constants__get__au2gpa()
    
    @property
    def golden(self):
        """
        Element golden ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 48
        
        """
        return _arespy.f90wrap_constants__get__golden()
    
    @property
    def vlight(self):
        """
        Element vlight ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 49
        
        """
        return _arespy.f90wrap_constants__get__vlight()
    
    @property
    def const_me_au(self):
        """
        Element const_me_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 53
        
        """
        return _arespy.f90wrap_constants__get__const_me_au()
    
    @property
    def const_e_au(self):
        """
        Element const_e_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 55
        
        """
        return _arespy.f90wrap_constants__get__const_e_au()
    
    @property
    def const_eh_au(self):
        """
        Element const_eh_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 57
        
        """
        return _arespy.f90wrap_constants__get__const_eh_au()
    
    @property
    def const_len_au(self):
        """
        Element const_len_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 59
        
        """
        return _arespy.f90wrap_constants__get__const_len_au()
    
    @property
    def const_hbar_au(self):
        """
        Element const_hbar_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 61
        
        """
        return _arespy.f90wrap_constants__get__const_hbar_au()
    
    @property
    def const_ep_au(self):
        """
        Element const_ep_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 63
        
        """
        return _arespy.f90wrap_constants__get__const_ep_au()
    
    @property
    def const_f_au(self):
        """
        Element const_f_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 65
        
        """
        return _arespy.f90wrap_constants__get__const_f_au()
    
    @property
    def const_mt_au(self):
        """
        Element const_mt_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 67
        
        """
        return _arespy.f90wrap_constants__get__const_mt_au()
    
    @property
    def const_time_au(self):
        """
        Element const_time_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 69
        
        """
        return _arespy.f90wrap_constants__get__const_time_au()
    
    @property
    def const_i_au(self):
        """
        Element const_i_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 71
        
        """
        return _arespy.f90wrap_constants__get__const_i_au()
    
    @property
    def const_temp_au(self):
        """
        Element const_temp_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 73
        
        """
        return _arespy.f90wrap_constants__get__const_temp_au()
    
    @property
    def const_p_au(self):
        """
        Element const_p_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 75
        
        """
        return _arespy.f90wrap_constants__get__const_p_au()
    
    @property
    def const_v_au(self):
        """
        Element const_v_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 77
        
        """
        return _arespy.f90wrap_constants__get__const_v_au()
    
    @property
    def const_ke_au(self):
        """
        Element const_ke_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 79
        
        """
        return _arespy.f90wrap_constants__get__const_ke_au()
    
    @property
    def const_mu_au(self):
        """
        Element const_mu_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 81
        
        """
        return _arespy.f90wrap_constants__get__const_mu_au()
    
    @property
    def const_ma_au(self):
        """
        Element const_ma_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 83
        
        """
        return _arespy.f90wrap_constants__get__const_ma_au()
    
    @property
    def const_kb_si(self):
        """
        Element const_kb_si ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 85
        
        """
        return _arespy.f90wrap_constants__get__const_kb_si()
    
    @property
    def angs(self):
        """
        Element angs ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 88
        
        """
        return _arespy.f90wrap_constants__get__angs()
    
    @property
    def force2ev(self):
        """
        Element force2ev ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 89
        
        """
        return _arespy.f90wrap_constants__get__force2ev()
    
    @property
    def clen(self):
        """
        Element clen ftype=integer pytype=int
        
        
        Defined at Constants.f90 line 91
        
        """
        return _arespy.f90wrap_constants__get__clen()
    
    @property
    def xtiny(self):
        """
        Element xtiny ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 92
        
        """
        return _arespy.f90wrap_constants__get__xtiny()
    
    @property
    def autogpa(self):
        """
        Element autogpa ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 93
        
        """
        return _arespy.f90wrap_constants__get__autogpa()
    
    @property
    def hartree2ev(self):
        """
        Element hartree2ev ftype=real(dp) pytype=float
        
        
        Defined at Constants.f90 line 94
        
        """
        return _arespy.f90wrap_constants__get__hartree2ev()
    
    def __str__(self):
        ret = ['<constants>{\n']
        ret.append('    i4b : ')
        ret.append(repr(self.i4b))
        ret.append(',\n    i2b : ')
        ret.append(repr(self.i2b))
        ret.append(',\n    i1b : ')
        ret.append(repr(self.i1b))
        ret.append(',\n    sp : ')
        ret.append(repr(self.sp))
        ret.append(',\n    dp : ')
        ret.append(repr(self.dp))
        ret.append(',\n    scp : ')
        ret.append(repr(self.scp))
        ret.append(',\n    dcp : ')
        ret.append(repr(self.dcp))
        ret.append(',\n    lgt : ')
        ret.append(repr(self.lgt))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    pio2 : ')
        ret.append(repr(self.pio2))
        ret.append(',\n    twopi : ')
        ret.append(repr(self.twopi))
        ret.append(',\n    imag : ')
        ret.append(repr(self.imag))
        ret.append(',\n    inputunit : ')
        ret.append(repr(self.inputunit))
        ret.append(',\n    errunit : ')
        ret.append(repr(self.errunit))
        ret.append(',\n    outputunit : ')
        ret.append(repr(self.outputunit))
        ret.append(',\n    mdposunit : ')
        ret.append(repr(self.mdposunit))
        ret.append(',\n    rydberg : ')
        ret.append(repr(self.rydberg))
        ret.append(',\n    bohr2ang : ')
        ret.append(repr(self.bohr2ang))
        ret.append(',\n    ang2bohr : ')
        ret.append(repr(self.ang2bohr))
        ret.append(',\n    hart2ev : ')
        ret.append(repr(self.hart2ev))
        ret.append(',\n    ev2hart : ')
        ret.append(repr(self.ev2hart))
        ret.append(',\n    scf_tol : ')
        ret.append(repr(self.scf_tol))
        ret.append(',\n    au2ev_force : ')
        ret.append(repr(self.au2ev_force))
        ret.append(',\n    au2gpa : ')
        ret.append(repr(self.au2gpa))
        ret.append(',\n    golden : ')
        ret.append(repr(self.golden))
        ret.append(',\n    vlight : ')
        ret.append(repr(self.vlight))
        ret.append(',\n    const_me_au : ')
        ret.append(repr(self.const_me_au))
        ret.append(',\n    const_e_au : ')
        ret.append(repr(self.const_e_au))
        ret.append(',\n    const_eh_au : ')
        ret.append(repr(self.const_eh_au))
        ret.append(',\n    const_len_au : ')
        ret.append(repr(self.const_len_au))
        ret.append(',\n    const_hbar_au : ')
        ret.append(repr(self.const_hbar_au))
        ret.append(',\n    const_ep_au : ')
        ret.append(repr(self.const_ep_au))
        ret.append(',\n    const_f_au : ')
        ret.append(repr(self.const_f_au))
        ret.append(',\n    const_mt_au : ')
        ret.append(repr(self.const_mt_au))
        ret.append(',\n    const_time_au : ')
        ret.append(repr(self.const_time_au))
        ret.append(',\n    const_i_au : ')
        ret.append(repr(self.const_i_au))
        ret.append(',\n    const_temp_au : ')
        ret.append(repr(self.const_temp_au))
        ret.append(',\n    const_p_au : ')
        ret.append(repr(self.const_p_au))
        ret.append(',\n    const_v_au : ')
        ret.append(repr(self.const_v_au))
        ret.append(',\n    const_ke_au : ')
        ret.append(repr(self.const_ke_au))
        ret.append(',\n    const_mu_au : ')
        ret.append(repr(self.const_mu_au))
        ret.append(',\n    const_ma_au : ')
        ret.append(repr(self.const_ma_au))
        ret.append(',\n    const_kb_si : ')
        ret.append(repr(self.const_kb_si))
        ret.append(',\n    angs : ')
        ret.append(repr(self.angs))
        ret.append(',\n    force2ev : ')
        ret.append(repr(self.force2ev))
        ret.append(',\n    clen : ')
        ret.append(repr(self.clen))
        ret.append(',\n    xtiny : ')
        ret.append(repr(self.xtiny))
        ret.append(',\n    autogpa : ')
        ret.append(repr(self.autogpa))
        ret.append(',\n    hartree2ev : ')
        ret.append(repr(self.hartree2ev))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

constants = Constants()

class End_Module(f90wrap.runtime.FortranModule):
    """
    Module end_module
    
    
    Defined at End_module.f90 lines 1-25
    
    """
    @staticmethod
    def destroy_beast():
        """
        destroy_beast()
        
        
        Defined at End_module.f90 lines 9-24
        
        
        """
        _arespy.f90wrap_destroy_beast()
    
    _dt_array_initialisers = []
    

end_module = End_Module()

class Energy_Module(f90wrap.runtime.FortranModule):
    """
    Module energy_module
    
    
    Defined at Energy_module.f90 lines 7-97
    
    """
    @staticmethod
    def totalenergy(nps, eig, rhos, rho):
        """
        totalenergy(nps, eig, rhos, rho)
        
        
        Defined at Energy_module.f90 lines 24-67
        
        Parameters
        ----------
        nps : int
        eig : Eigen_Type
        rhos : float array
        rho : float array
        
        """
        _arespy.f90wrap_totalenergy(nps=nps, eig=eig._handle, rhos=rhos, rho=rho)
    
    @staticmethod
    def ebands(eval, wke):
        """
        eband = ebands(eval, wke)
        
        
        Defined at Energy_module.f90 lines 70-95
        
        Parameters
        ----------
        eval : float array
        wke : float array
        
        Returns
        -------
        eband : float
        
        """
        eband = _arespy.f90wrap_ebands(eval=eval, wke=wke)
        return eband
    
    @property
    def etot(self):
        """
        Element etot ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__etot()
    
    @etot.setter
    def etot(self, etot):
        _arespy.f90wrap_energy_module__set__etot(etot)
    
    @property
    def eband(self):
        """
        Element eband ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__eband()
    
    @eband.setter
    def eband(self, eband):
        _arespy.f90wrap_energy_module__set__eband(eband)
    
    @property
    def eh(self):
        """
        Element eh ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__eh()
    
    @eh.setter
    def eh(self, eh):
        _arespy.f90wrap_energy_module__set__eh(eh)
    
    @property
    def eext(self):
        """
        Element eext ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__eext()
    
    @eext.setter
    def eext(self, eext):
        _arespy.f90wrap_energy_module__set__eext(eext)
    
    @property
    def exc(self):
        """
        Element exc ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__exc()
    
    @exc.setter
    def exc(self, exc):
        _arespy.f90wrap_energy_module__set__exc(exc)
    
    @property
    def eele(self):
        """
        Element eele ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__eele()
    
    @eele.setter
    def eele(self, eele):
        _arespy.f90wrap_energy_module__set__eele(eele)
    
    @property
    def fe(self):
        """
        Element fe ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__fe()
    
    @fe.setter
    def fe(self, fe):
        _arespy.f90wrap_energy_module__set__fe(fe)
    
    @property
    def fe0(self):
        """
        Element fe0 ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.f90 line 21
        
        """
        return _arespy.f90wrap_energy_module__get__fe0()
    
    @fe0.setter
    def fe0(self, fe0):
        _arespy.f90wrap_energy_module__set__fe0(fe0)
    
    def __str__(self):
        ret = ['<energy_module>{\n']
        ret.append('    etot : ')
        ret.append(repr(self.etot))
        ret.append(',\n    eband : ')
        ret.append(repr(self.eband))
        ret.append(',\n    eh : ')
        ret.append(repr(self.eh))
        ret.append(',\n    eext : ')
        ret.append(repr(self.eext))
        ret.append(',\n    exc : ')
        ret.append(repr(self.exc))
        ret.append(',\n    eele : ')
        ret.append(repr(self.eele))
        ret.append(',\n    fe : ')
        ret.append(repr(self.fe))
        ret.append(',\n    fe0 : ')
        ret.append(repr(self.fe0))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

energy_module = Energy_Module()

class Output_Module(f90wrap.runtime.FortranModule):
    """
    Module output_module
    
    
    Defined at Energy_module.f90 lines 112-232
    
    """
    @staticmethod
    def output():
        """
        output()
        
        
        Defined at Energy_module.f90 lines 124-174
        
        
        """
        _arespy.f90wrap_output()
    
    @staticmethod
    def write_density():
        """
        write_density()
        
        
        Defined at Energy_module.f90 lines 177-204
        
        
        """
        _arespy.f90wrap_write_density()
    
    @staticmethod
    def write_band():
        """
        write_band()
        
        
        Defined at Energy_module.f90 lines 208-231
        
        
        """
        _arespy.f90wrap_write_band()
    
    @property
    def time_total0(self):
        """
        Element time_total0 ftype=integer(i4b) pytype=int
        
        
        Defined at Energy_module.f90 line 121
        
        """
        return _arespy.f90wrap_output_module__get__time_total0()
    
    @time_total0.setter
    def time_total0(self, time_total0):
        _arespy.f90wrap_output_module__set__time_total0(time_total0)
    
    @property
    def time_total1(self):
        """
        Element time_total1 ftype=integer(i4b) pytype=int
        
        
        Defined at Energy_module.f90 line 121
        
        """
        return _arespy.f90wrap_output_module__get__time_total1()
    
    @time_total1.setter
    def time_total1(self, time_total1):
        _arespy.f90wrap_output_module__set__time_total1(time_total1)
    
    @property
    def time_scf0(self):
        """
        Element time_scf0 ftype=integer(i4b) pytype=int
        
        
        Defined at Energy_module.f90 line 121
        
        """
        return _arespy.f90wrap_output_module__get__time_scf0()
    
    @time_scf0.setter
    def time_scf0(self, time_scf0):
        _arespy.f90wrap_output_module__set__time_scf0(time_scf0)
    
    @property
    def time_scf1(self):
        """
        Element time_scf1 ftype=integer(i4b) pytype=int
        
        
        Defined at Energy_module.f90 line 121
        
        """
        return _arespy.f90wrap_output_module__get__time_scf1()
    
    @time_scf1.setter
    def time_scf1(self, time_scf1):
        _arespy.f90wrap_output_module__set__time_scf1(time_scf1)
    
    def __str__(self):
        ret = ['<output_module>{\n']
        ret.append('    time_total0 : ')
        ret.append(repr(self.time_total0))
        ret.append(',\n    time_total1 : ')
        ret.append(repr(self.time_total1))
        ret.append(',\n    time_scf0 : ')
        ret.append(repr(self.time_scf0))
        ret.append(',\n    time_scf1 : ')
        ret.append(repr(self.time_scf1))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

output_module = Output_Module()

class Ewald(f90wrap.runtime.FortranModule):
    """
    Module ewald
    
    
    Defined at Ewald.f90 lines 1-749
    
    """
    @f90wrap.runtime.register_class("arespy.ion")
    class ion(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=ion)
        
        
        Defined at Ewald.f90 lines 32-34
        
        """
        def __init__(self, handle=None):
            """
            self = Ion()
            
            
            Defined at Ewald.f90 lines 32-34
            
            
            Returns
            -------
            this : Ion
            	Object to be constructed
            
            
            Automatically generated constructor for ion
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_ion_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Ion
            
            
            Defined at Ewald.f90 lines 32-34
            
            Parameters
            ----------
            this : Ion
            	Object to be destructed
            
            
            Automatically generated destructor for ion
            """
            if self._alloc:
                _arespy.f90wrap_ion_finalise(this=self._handle)
        
        @property
        def charge(self):
            """
            Element charge ftype=real(dp) pytype=float
            
            
            Defined at Ewald.f90 line 33
            
            """
            return _arespy.f90wrap_ion__get__charge(self._handle)
        
        @charge.setter
        def charge(self, charge):
            _arespy.f90wrap_ion__set__charge(self._handle, charge)
        
        @property
        def fcd(self):
            """
            Element fcd ftype=real(dp) pytype=float
            
            
            Defined at Ewald.f90 line 34
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_ion__array__fcd(self._handle)
            if array_handle in self._arrays:
                fcd = self._arrays[array_handle]
            else:
                fcd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_ion__array__fcd)
                self._arrays[array_handle] = fcd
            return fcd
        
        @fcd.setter
        def fcd(self, fcd):
            self.fcd[...] = fcd
        
        def __str__(self):
            ret = ['<ion>{\n']
            ret.append('    charge : ')
            ret.append(repr(self.charge))
            ret.append(',\n    fcd : ')
            ret.append(repr(self.fcd))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def ewald_energy(latticev, ionpositions, iontpid, ioncharges):
        """
        ewald_energy = ewald_energy(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.f90 lines 54-72
        
        Parameters
        ----------
        latticev : float array
        ionpositions : float array
        iontpid : int array
        ioncharges : float array
        
        Returns
        -------
        ewald_energy : float
        
        """
        ewald_energy = _arespy.f90wrap_ewald_energy(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return ewald_energy
    
    @staticmethod
    def iso_ewald_energy(latticev, ionpositions, iontpid, ioncharges):
        """
        iso_ewald_energy = iso_ewald_energy(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.f90 lines 76-104
        
        Parameters
        ----------
        latticev : float array
        ionpositions : float array
        iontpid : int array
        ioncharges : float array
        
        Returns
        -------
        iso_ewald_energy : float
        
        """
        iso_ewald_energy = _arespy.f90wrap_iso_ewald_energy(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return iso_ewald_energy
    
    @staticmethod
    def iso_ewald_forces(latticev, ionpositions, iontpid, ioncharges):
        """
        iso_ewald_forces = iso_ewald_forces(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.f90 lines 108-203
        
        Parameters
        ----------
        latticev : float array
        ionpositions : float array
        iontpid : int array
        ioncharges : float array
        
        Returns
        -------
        iso_ewald_forces : float array
        
        ==========================================
        """
        iso_ewald_forces = _arespy.f90wrap_iso_ewald_forces(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return iso_ewald_forces
    
    @staticmethod
    def ewald_forces(latticev, ionpositions, iontpid, ioncharges):
        """
        ewald_forces = ewald_forces(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.f90 lines 206-229
        
        Parameters
        ----------
        latticev : float array
        ionpositions : float array
        iontpid : int array
        ioncharges : float array
        
        Returns
        -------
        ewald_forces : float array
        
        """
        ewald_forces = _arespy.f90wrap_ewald_forces(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return ewald_forces
    
    @staticmethod
    def ewald_stress(latticev, ionpositions, iontpid, ioncharges):
        """
        ewald_stress = ewald_stress(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.f90 lines 232-250
        
        Parameters
        ----------
        latticev : float array
        ionpositions : float array
        iontpid : int array
        ioncharges : float array
        
        Returns
        -------
        ewald_stress : float array
        
        """
        ewald_stress = _arespy.f90wrap_ewald_stress(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return ewald_stress
    
    @staticmethod
    def ewaldrpstr(eta):
        """
        ewaldrpstr = ewaldrpstr(eta)
        
        
        Defined at Ewald.f90 lines 582-626
        
        Parameters
        ----------
        eta : float
        
        Returns
        -------
        ewaldrpstr : float array
        
        """
        ewaldrpstr = _arespy.f90wrap_ewaldrpstr(eta=eta)
        return ewaldrpstr
    
    @staticmethod
    def ewaldavstr(eta):
        """
        ewaldavstr = ewaldavstr(eta)
        
        
        Defined at Ewald.f90 lines 628-635
        
        Parameters
        ----------
        eta : float
        
        Returns
        -------
        ewaldavstr : float array
        
        """
        ewaldavstr = _arespy.f90wrap_ewaldavstr(eta=eta)
        return ewaldavstr
    
    @staticmethod
    def vectorlength(vc):
        """
        vectorlength = vectorlength(vc)
        
        
        Defined at Ewald.f90 lines 638-639
        
        Parameters
        ----------
        vc : float array
        
        Returns
        -------
        vectorlength : float
        
        """
        vectorlength = _arespy.f90wrap_vectorlength(vc=vc)
        return vectorlength
    
    @staticmethod
    def recipvector(lat):
        """
        recipvector = recipvector(lat)
        
        
        Defined at Ewald.f90 lines 642-647
        
        Parameters
        ----------
        lat : float array
        
        Returns
        -------
        recipvector : float array
        
        """
        recipvector = _arespy.f90wrap_recipvector(lat=lat)
        return recipvector
    
    @staticmethod
    def volume(lat):
        """
        volume = volume(lat)
        
        
        Defined at Ewald.f90 lines 650-652
        
        Parameters
        ----------
        lat : float array
        
        Returns
        -------
        volume : float
        
        """
        volume = _arespy.f90wrap_volume(lat=lat)
        return volume
    
    @staticmethod
    def crossp(va, vb):
        """
        crossp = crossp(va, vb)
        
        
        Defined at Ewald.f90 lines 655-659
        
        Parameters
        ----------
        va : float array
        vb : float array
        
        Returns
        -------
        crossp : float array
        
        """
        crossp = _arespy.f90wrap_crossp(va=va, vb=vb)
        return crossp
    
    @staticmethod
    def erfc(x):
        """
        erfc = erfc(x)
        
        
        Defined at Ewald.f90 lines 662-748
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        erfc : float
        
        """
        erfc = _arespy.f90wrap_erfc(x=x)
        return erfc
    
    @property
    def bohr(self):
        """
        Element bohr ftype=real(dp) pytype=float
        
        
        Defined at Ewald.f90 line 36
        
        """
        return _arespy.f90wrap_ewald__get__bohr()
    
    @property
    def hartreetoev(self):
        """
        Element hartreetoev ftype=real(dp) pytype=float
        
        
        Defined at Ewald.f90 line 36
        
        """
        return _arespy.f90wrap_ewald__get__hartreetoev()
    
    def __str__(self):
        ret = ['<ewald>{\n']
        ret.append('    bohr : ')
        ret.append(repr(self.bohr))
        ret.append(',\n    hartreetoev : ')
        ret.append(repr(self.hartreetoev))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

ewald = Ewald()

class Finite_Module(f90wrap.runtime.FortranModule):
    """
    Module finite_module
    
    
    Defined at Finite_module.f90 lines 1-614
    
    """
    @staticmethod
    def destroy_finite():
        """
        destroy_finite()
        
        
        Defined at Finite_module.f90 lines 33-38
        
        
        """
        _arespy.f90wrap_destroy_finite()
    
    @staticmethod
    def init_finite(h):
        """
        init_finite(h)
        
        
        Defined at Finite_module.f90 lines 41-107
        
        Parameters
        ----------
        h : float array
        
        """
        _arespy.f90wrap_init_finite(h=h)
    
    @staticmethod
    def trans_mat_full(mat, factor, int_miu, lad_gap, err):
        """
        trans_mat_full(mat, factor, int_miu, lad_gap, err)
        
        
        Defined at Finite_module.f90 lines 110-318
        
        Parameters
        ----------
        mat : float array
        factor : float array
        int_miu : int array
        lad_gap : float array
        err : float array
        
        """
        _arespy.f90wrap_trans_mat_full(mat=mat, factor=factor, int_miu=int_miu, \
            lad_gap=lad_gap, err=err)
    
    @staticmethod
    def cmplx_keop(uk, ik, ts):
        """
        cmplx_keop(uk, ik, ts)
        
        
        Defined at Finite_module.f90 lines 322-345
        
        Parameters
        ----------
        uk : complex array
        ik : int
        ts : complex array
        
        """
        _arespy.f90wrap_cmplx_keop(uk=uk, ik=ik, ts=ts)
    
    @staticmethod
    def real_comm(ifun):
        """
        real_comm(ifun)
        
        
        Defined at Finite_module.f90 lines 352-380
        
        Parameters
        ----------
        ifun : float array
        
        """
        _arespy.f90wrap_real_comm(ifun=ifun)
    
    @staticmethod
    def real_comm_clean():
        """
        real_comm_clean()
        
        
        Defined at Finite_module.f90 lines 383-385
        
        
        """
        _arespy.f90wrap_real_comm_clean()
    
    @staticmethod
    def real_nabla1_3d(ifun, norder, derf):
        """
        real_nabla1_3d(ifun, norder, derf)
        
        
        Defined at Finite_module.f90 lines 388-436
        
        Parameters
        ----------
        ifun : float array
        norder : int
        derf : float array
        
        """
        _arespy.f90wrap_real_nabla1_3d(ifun=ifun, norder=norder, derf=derf)
    
    @staticmethod
    def real_nabla2_3d(func, norder, ofun):
        """
        real_nabla2_3d(func, norder, ofun)
        
        
        Defined at Finite_module.f90 lines 439-506
        
        Parameters
        ----------
        func : float array
        norder : int
        ofun : float array
        
        """
        _arespy.f90wrap_real_nabla2_3d(func=func, norder=norder, ofun=ofun)
    
    @staticmethod
    def real_nabla1_1d(ifun, derf, mgfun=None):
        """
        real_nabla1_1d(ifun, derf[, mgfun])
        
        
        Defined at Finite_module.f90 lines 510-560
        
        Parameters
        ----------
        ifun : float array
        derf : float array
        mgfun : float array
        
        """
        _arespy.f90wrap_real_nabla1_1d(ifun=ifun, derf=derf, mgfun=mgfun)
    
    @staticmethod
    def real_nabla2_1d(ifun, ofun):
        """
        real_nabla2_1d(ifun, ofun)
        
        
        Defined at Finite_module.f90 lines 563-612
        
        Parameters
        ----------
        ifun : float array
        ofun : float array
        
        """
        _arespy.f90wrap_real_nabla2_1d(ifun=ifun, ofun=ofun)
    
    @property
    def lapl(self):
        """
        Element lapl ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__lapl(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lapl = self._arrays[array_handle]
        else:
            lapl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__lapl)
            self._arrays[array_handle] = lapl
        return lapl
    
    @lapl.setter
    def lapl(self, lapl):
        self.lapl[...] = lapl
    
    @property
    def grad(self):
        """
        Element grad ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.f90 line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__grad(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            grad = self._arrays[array_handle]
        else:
            grad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__grad)
            self._arrays[array_handle] = grad
        return grad
    
    @grad.setter
    def grad(self, grad):
        self.grad[...] = grad
    
    @property
    def tbmat(self):
        """
        Element tbmat ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.f90 line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__tbmat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tbmat = self._arrays[array_handle]
        else:
            tbmat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__tbmat)
            self._arrays[array_handle] = tbmat
        return tbmat
    
    @tbmat.setter
    def tbmat(self, tbmat):
        self.tbmat[...] = tbmat
    
    @property
    def lap_add(self):
        """
        Element lap_add ftype=integer(i4b) pytype=int
        
        
        Defined at Finite_module.f90 line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__lap_add(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lap_add = self._arrays[array_handle]
        else:
            lap_add = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__lap_add)
            self._arrays[array_handle] = lap_add
        return lap_add
    
    @lap_add.setter
    def lap_add(self, lap_add):
        self.lap_add[...] = lap_add
    
    @property
    def cell_mu(self):
        """
        Element cell_mu ftype=integer  pytype=int
        
        
        Defined at Finite_module.f90 line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__cell_mu(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cell_mu = self._arrays[array_handle]
        else:
            cell_mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__cell_mu)
            self._arrays[array_handle] = cell_mu
        return cell_mu
    
    @cell_mu.setter
    def cell_mu(self, cell_mu):
        self.cell_mu[...] = cell_mu
    
    @property
    def cell_factor(self):
        """
        Element cell_factor ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.f90 line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__cell_factor(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cell_factor = self._arrays[array_handle]
        else:
            cell_factor = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__cell_factor)
            self._arrays[array_handle] = cell_factor
        return cell_factor
    
    @cell_factor.setter
    def cell_factor(self, cell_factor):
        self.cell_factor[...] = cell_factor
    
    @property
    def wrap_real(self):
        """
        Element wrap_real ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.f90 line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_finite_module__array__wrap_real(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wrap_real = self._arrays[array_handle]
        else:
            wrap_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_finite_module__array__wrap_real)
            self._arrays[array_handle] = wrap_real
        return wrap_real
    
    @wrap_real.setter
    def wrap_real(self, wrap_real):
        self.wrap_real[...] = wrap_real
    
    def __str__(self):
        ret = ['<finite_module>{\n']
        ret.append('    lapl : ')
        ret.append(repr(self.lapl))
        ret.append(',\n    grad : ')
        ret.append(repr(self.grad))
        ret.append(',\n    tbmat : ')
        ret.append(repr(self.tbmat))
        ret.append(',\n    lap_add : ')
        ret.append(repr(self.lap_add))
        ret.append(',\n    cell_mu : ')
        ret.append(repr(self.cell_mu))
        ret.append(',\n    cell_factor : ')
        ret.append(repr(self.cell_factor))
        ret.append(',\n    wrap_real : ')
        ret.append(repr(self.wrap_real))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

finite_module = Finite_Module()

class Forcestress_Module(f90wrap.runtime.FortranModule):
    """
    Module forcestress_module
    
    
    Defined at ForceStress_module.f90 lines 6-1970
    
    """
    @staticmethod
    def cal_force_c(rhos, uik, force):
        """
        cal_force_c(rhos, uik, force)
        
        
        Defined at ForceStress_module.f90 lines 16-86
        
        Parameters
        ----------
        rhos : float array
        uik : complex array
        force : float array
        
        ===================================
        nlnlocal part
        """
        _arespy.f90wrap_cal_force_c(rhos=rhos, uik=uik, force=force)
    
    @staticmethod
    def cal_force_r(rhos, uik, force):
        """
        cal_force_r(rhos, uik, force)
        
        
        Defined at ForceStress_module.f90 lines 89-197
        
        Parameters
        ----------
        rhos : float array
        uik : float array
        force : float array
        
        ===================================
        nlnlocal part
        """
        _arespy.f90wrap_cal_force_r(rhos=rhos, uik=uik, force=force)
    
    @staticmethod
    def locforce(rhos, lforce):
        """
        locforce(rhos, lforce)
        
        
        Defined at ForceStress_module.f90 lines 201-355
        
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
        _arespy.f90wrap_locforce(rhos=rhos, lforce=lforce)
    
    @staticmethod
    def locforce_r(rhos, lforce):
        """
        locforce_r(rhos, lforce)
        
        
        Defined at ForceStress_module.f90 lines 358-478
        
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
        _arespy.f90wrap_locforce_r(rhos=rhos, lforce=lforce)
    
    @staticmethod
    def nonlforce(uik, nlforce):
        """
        nonlforce(uik, nlforce)
        
        
        Defined at ForceStress_module.f90 lines 587-732
        
        Parameters
        ----------
        uik : complex array
        nlforce : float array
        
        ---------------------------------
        """
        _arespy.f90wrap_nonlforce(uik=uik, nlforce=nlforce)
    
    @staticmethod
    def nonlforce_r_dg(uik, nlforce):
        """
        nonlforce_r_dg(uik, nlforce)
        
        
        Defined at ForceStress_module.f90 lines 736-884
        
        Parameters
        ----------
        uik : float array
        nlforce : float array
        
        ---------------------------------
        """
        _arespy.f90wrap_nonlforce_r_dg(uik=uik, nlforce=nlforce)
    
    @staticmethod
    def nonlforce_r(uik, nlforce):
        """
        nonlforce_r(uik, nlforce)
        
        
        Defined at ForceStress_module.f90 lines 887-1046
        
        Parameters
        ----------
        uik : float array
        nlforce : float array
        
        ---------------------------------
        """
        _arespy.f90wrap_nonlforce_r(uik=uik, nlforce=nlforce)
    
    @staticmethod
    def cal_stress(rhos, uik, stress):
        """
        cal_stress(rhos, uik, stress)
        
        
        Defined at ForceStress_module.f90 lines 1050-1171
        
        Parameters
        ----------
        rhos : float array
        uik : complex array
        stress : float array
        
        """
        _arespy.f90wrap_cal_stress(rhos=rhos, uik=uik, stress=stress)
    
    @staticmethod
    def lda_stress(vxc, rho, elda):
        """
        lda_stress = lda_stress(vxc, rho, elda)
        
        
        Defined at ForceStress_module.f90 lines 1174-1216
        
        Parameters
        ----------
        vxc : float array
        rho : float array
        elda : float
        
        Returns
        -------
        lda_stress : float array
        
        """
        lda_stress = _arespy.f90wrap_lda_stress(vxc=vxc, rho=rho, elda=elda)
        return lda_stress
    
    @staticmethod
    def hart_stress(rhorecip, eh):
        """
        hart_stress = hart_stress(rhorecip, eh)
        
        
        Defined at ForceStress_module.f90 lines 1219-1274
        
        Parameters
        ----------
        rhorecip : complex array
        eh : float
        
        Returns
        -------
        hart_stress : float array
        
        """
        hart_stress = _arespy.f90wrap_hart_stress(rhorecip=rhorecip, eh=eh)
        return hart_stress
    
    @staticmethod
    def kin_stress(uik, kinstress):
        """
        kin_stress(uik, kinstress)
        
        
        Defined at ForceStress_module.f90 lines 1277-1438
        
        Parameters
        ----------
        uik : complex array
        kinstress : float array
        
        """
        _arespy.f90wrap_kin_stress(uik=uik, kinstress=kinstress)
    
    @staticmethod
    def ion_nl_stress(uik, nlstress):
        """
        ion_nl_stress(uik, nlstress)
        
        
        Defined at ForceStress_module.f90 lines 1441-1634
        
        Parameters
        ----------
        uik : complex array
        nlstress : float array
        
        ----------------------------------
        """
        _arespy.f90wrap_ion_nl_stress(uik=uik, nlstress=nlstress)
    
    @staticmethod
    def ionele_stress(rhorecip, energy):
        """
        ionele_stress = ionele_stress(rhorecip, energy)
        
        
        Defined at ForceStress_module.f90 lines 1637-1727
        
        Parameters
        ----------
        rhorecip : complex array
        energy : float
        
        Returns
        -------
        ionele_stress : float array
        
        """
        ionele_stress = _arespy.f90wrap_ionele_stress(rhorecip=rhorecip, energy=energy)
        return ionele_stress
    
    @staticmethod
    def pseudopotdifflookup(ity, qnorm):
        """
        pseudopotdifflookup = pseudopotdifflookup(ity, qnorm)
        
        
        Defined at ForceStress_module.f90 lines 1731-1770
        
        Parameters
        ----------
        ity : int
        qnorm : float
        
        Returns
        -------
        pseudopotdifflookup : float
        
        """
        pseudopotdifflookup = _arespy.f90wrap_pseudopotdifflookup(ity=ity, qnorm=qnorm)
        return pseudopotdifflookup
    
    @staticmethod
    def cal_force_stress():
        """
        cal_force_stress()
        
        
        Defined at ForceStress_module.f90 lines 1774-1969
        
        
        ====================
        ##cal force directly
        goto 10011
        ====================
        Self-consistent
        """
        _arespy.f90wrap_cal_force_stress()
    
    @property
    def cellpress(self):
        """
        Element cellpress ftype=real(dp) pytype=float
        
        
        Defined at ForceStress_module.f90 line 9
        
        """
        return _arespy.f90wrap_forcestress_module__get__cellpress()
    
    @cellpress.setter
    def cellpress(self, cellpress):
        _arespy.f90wrap_forcestress_module__set__cellpress(cellpress)
    
    def __str__(self):
        ret = ['<forcestress_module>{\n']
        ret.append('    cellpress : ')
        ret.append(repr(self.cellpress))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

forcestress_module = Forcestress_Module()

class Fourier(f90wrap.runtime.FortranModule):
    """
    Module fourier
    
    
    Defined at Fourier.f90 lines 1-567
    
    """
    @staticmethod
    def planfft(dimx, dimy, dimz, lsp):
        """
        planfft(dimx, dimy, dimz, lsp)
        
        
        Defined at Fourier.f90 lines 89-189
        
        Parameters
        ----------
        dimx : int
        dimy : int
        dimz : int
        lsp : bool array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This is the initialization procedure that first gets the system name as is
           called as an argument to OFDFT, and turns it into the various input file
           names. Then, it calls all the programs necessary to set variables to
           default values, then reads the geometry file to get all the variables sets
           to the correct values.
         GLOBAL/MODULE VARIABLES CHANGED:
           realRA, cplxRA, offset
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003  File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        _arespy.f90wrap_planfft(dimx=dimx, dimy=dimy, dimz=dimz, lsp=lsp)
    
    @staticmethod
    def planfst(dimx, dimy, dimz):
        """
        planfst(dimx, dimy, dimz)
        
        
        Defined at Fourier.f90 lines 191-223
        
        Parameters
        ----------
        dimx : int
        dimy : int
        dimz : int
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This is the same as PlanFFT, except for the Fast Sine Transform.
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003  File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        _arespy.f90wrap_planfst(dimx=dimx, dimy=dimy, dimz=dimz)
    
    @staticmethod
    def getfftdims():
        """
        dimx, dimy, dimz = getfftdims()
        
        
        Defined at Fourier.f90 lines 225-252
        
        
        Returns
        -------
        dimx : int
        dimy : int
        dimz : int
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           Gets the dimensions of the FFT(real-space part)
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           4/25/2006  Added(GSH)
        ------------------------------------------------------------------------------
        """
        dimx, dimy, dimz = _arespy.f90wrap_getfftdims()
        return dimx, dimy, dimz
    
    @staticmethod
    def getfftcomplexdims():
        """
        dimx, dimy, dimz = getfftcomplexdims()
        
        
        Defined at Fourier.f90 lines 254-281
        
        
        Returns
        -------
        dimx : int
        dimy : int
        dimz : int
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           Gets the dimensions of the FFT(reciprocal space part)
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           4/26/2006 Added(GSH)
        ------------------------------------------------------------------------------
        """
        dimx, dimy, dimz = _arespy.f90wrap_getfftcomplexdims()
        return dimx, dimy, dimz
    
    @staticmethod
    def forwardfst(array):
        """
        forwardfst(array)
        
        
        Defined at Fourier.f90 lines 466-491
        
        Parameters
        ----------
        array : float array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
        ------------------------------------------------------------------------------
        """
        _arespy.f90wrap_forwardfst(array=array)
    
    @staticmethod
    def backfst(array):
        """
        backfst(array)
        
        
        Defined at Fourier.f90 lines 493-516
        
        Parameters
        ----------
        array : float array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
        ------------------------------------------------------------------------------
        """
        _arespy.f90wrap_backfst(array=array)
    
    @staticmethod
    def cleanfft():
        """
        cleanfft()
        
        
        Defined at Fourier.f90 lines 518-567
        
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This subroutine is called at the end of the run to free the memory
           associated with the plan.
         GLOBAL/MODULE VARIABLES CHANGED:
           realRA, cplxRA
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003  File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        _arespy.f90wrap_cleanfft()
    
    @staticmethod
    def _forwardfft_4d(array):
        """
        transform = _forwardfft_4d(array)
        
        
        Defined at Fourier.f90 lines 283-326
        
        Parameters
        ----------
        array : float array
        
        Returns
        -------
        transform : complex array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This function is not called directly from the OFDFT code. Use the FFT
           interface instead. It performs the transformation of a real 4-dimensional
           array into its complex 4-dimensional transform. The first dimension is
           halved.
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003  File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        transform = _arespy.f90wrap_forwardfft_4d(array=array)
        return transform
    
    @staticmethod
    def _backfft_4d(array):
        """
        transform = _backfft_4d(array)
        
        
        Defined at Fourier.f90 lines 328-369
        
        Parameters
        ----------
        array : complex array
        
        Returns
        -------
        transform : float array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This function is not called directly from the OFDFT code, but rather
           through the FFT interface. It performs the reverse Fourier transform of
           a complex function over the half-box in reciprocal space back to real
           space. It acts on 4-dimensional arrays, the fourth dimension being spin.
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003    File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        transform = _arespy.f90wrap_backfft_4d(array=array)
        return transform
    
    @staticmethod
    def _forwardfft_3d(array):
        """
        transform = _forwardfft_3d(array)
        
        
        Defined at Fourier.f90 lines 371-419
        
        Parameters
        ----------
        array : float array
        
        Returns
        -------
        transform : complex array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This function is not called directly from the OFDFT code. Use the FFT
           interface instead. It performs the transformation of a real 4-dimensional
           array into its complex 4-dimensional transform. The first dimension is
           halved.
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003  File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        transform = _arespy.f90wrap_forwardfft_3d(array=array)
        return transform
    
    @staticmethod
    def _backfft_3d(array):
        """
        transform = _backfft_3d(array)
        
        
        Defined at Fourier.f90 lines 421-464
        
        Parameters
        ----------
        array : complex array
        
        Returns
        -------
        transform : float array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This function is not called directly from the OFDFT code, but rather
           through the FFT interface. It performs the reverse Fourier transform of a
           complex function over the half-box in reciprocal space back to real
           space. It acts on 3-dimensional arrays.
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
         REFERENCES:
        ------------------------------------------------------------------------------
         REVISION LOG:
           11/20/2003  File created.  (VLL)
        ------------------------------------------------------------------------------
        """
        transform = _arespy.f90wrap_backfft_3d(array=array)
        return transform
    
    @staticmethod
    def fft(*args, **kwargs):
        """
        fft(*args, **kwargs)
        
        
        Defined at Fourier.f90 lines 82-86
        
        Overloaded interface containing the following procedures:
          _forwardfft_4d
          _backfft_4d
          _forwardfft_3d
          _backfft_3d
        
        """
        for proc in [Fourier._forwardfft_4d, Fourier._backfft_4d, \
            Fourier._forwardfft_3d, Fourier._backfft_3d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def offset(self):
        """
        Element offset ftype=integer(i4b) pytype=int
        
        
        Defined at Fourier.f90 line 77
        
        """
        return _arespy.f90wrap_fourier__get__offset()
    
    @offset.setter
    def offset(self, offset):
        _arespy.f90wrap_fourier__set__offset(offset)
    
    def __str__(self):
        ret = ['<fourier>{\n']
        ret.append('    offset : ')
        ret.append(repr(self.offset))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

fourier = Fourier()

class Lapack_Module(f90wrap.runtime.FortranModule):
    """
    Module lapack_module
    
    
    Defined at Lapack_module.f90 lines 1-522
    
    """
    @staticmethod
    def diagm(mat, evec, eval):
        """
        diagm(mat, evec, eval)
        
        
        Defined at Lapack_module.f90 lines 16-54
        
        Parameters
        ----------
        mat : complex array
        evec : complex array
        eval : float array
        
        """
        _arespy.f90wrap_diagm(mat=mat, evec=evec, eval=eval)
    
    @staticmethod
    def generalizeeigen(dime, mata, matb, evec, eval):
        """
        generalizeeigen(dime, mata, matb, evec, eval)
        
        
        Defined at Lapack_module.f90 lines 63-99
        
        Parameters
        ----------
        dime : int
        mata : complex array
        matb : complex array
        evec : complex array
        eval : float array
        
        """
        _arespy.f90wrap_generalizeeigen(dime=dime, mata=mata, matb=matb, evec=evec, \
            eval=eval)
    
    @staticmethod
    def orthnorm(mat):
        """
        orthnorm(mat)
        
        
        Defined at Lapack_module.f90 lines 104-142
        
        Parameters
        ----------
        mat : complex array
        
        """
        _arespy.f90wrap_orthnorm(mat=mat)
    
    @staticmethod
    def norm_2(mat, k):
        """
        norm_2 = norm_2(mat, k)
        
        
        Defined at Lapack_module.f90 lines 145-159
        
        Parameters
        ----------
        mat : complex array
        k : int
        
        Returns
        -------
        norm_2 : float
        
        """
        norm_2 = _arespy.f90wrap_norm_2(mat=mat, k=k)
        return norm_2
    
    @staticmethod
    def matmat(mata, matb, opa, opb, matc):
        """
        matmat(mata, matb, opa, opb, matc)
        
        
        Defined at Lapack_module.f90 lines 167-197
        
        Parameters
        ----------
        mata : complex array
        matb : complex array
        opa : str
        opb : str
        matc : complex array
        
        """
        _arespy.f90wrap_matmat(mata=mata, matb=matb, opa=opa, opb=opb, matc=matc)
    
    @staticmethod
    def invmat(mat):
        """
        invmat(mat)
        
        
        Defined at Lapack_module.f90 lines 200-232
        
        Parameters
        ----------
        mat : complex array
        
        """
        _arespy.f90wrap_invmat(mat=mat)
    
    @staticmethod
    def orthnorm_real(mat):
        """
        orthnorm_real(mat)
        
        
        Defined at Lapack_module.f90 lines 256-296
        
        Parameters
        ----------
        mat : float array
        
        """
        _arespy.f90wrap_orthnorm_real(mat=mat)
    
    @staticmethod
    def diagm_real(mat, evec, eval):
        """
        diagm_real(mat, evec, eval)
        
        
        Defined at Lapack_module.f90 lines 301-342
        
        Parameters
        ----------
        mat : float array
        evec : float array
        eval : float array
        
        """
        _arespy.f90wrap_diagm_real(mat=mat, evec=evec, eval=eval)
    
    @staticmethod
    def generalizeeigen_real(dime, mata, matb, evec, eval):
        """
        generalizeeigen_real(dime, mata, matb, evec, eval)
        
        
        Defined at Lapack_module.f90 lines 350-387
        
        Parameters
        ----------
        dime : int
        mata : float array
        matb : float array
        evec : float array
        eval : float array
        
        """
        _arespy.f90wrap_generalizeeigen_real(dime=dime, mata=mata, matb=matb, evec=evec, \
            eval=eval)
    
    @staticmethod
    def diagmx_real(mat, dime, num, il, iu, evec, eval):
        """
        diagmx_real(mat, dime, num, il, iu, evec, eval)
        
        
        Defined at Lapack_module.f90 lines 390-430
        
        Parameters
        ----------
        mat : float array
        dime : int
        num : int
        il : int
        iu : int
        evec : float array
        eval : float array
        
        -------------------OUT PUT-----------------------------
        eigen-values
        """
        _arespy.f90wrap_diagmx_real(mat=mat, dime=dime, num=num, il=il, iu=iu, \
            evec=evec, eval=eval)
    
    @staticmethod
    def matmat_real(mata, matb, opa, opb, matc):
        """
        matmat_real(mata, matb, opa, opb, matc)
        
        
        Defined at Lapack_module.f90 lines 433-464
        
        Parameters
        ----------
        mata : float array
        matb : float array
        opa : str
        opb : str
        matc : float array
        
        """
        _arespy.f90wrap_matmat_real(mata=mata, matb=matb, opa=opa, opb=opb, matc=matc)
    
    @staticmethod
    def cholesky_factor_real(mat):
        """
        cholesky_factor_real(mat)
        
        
        Defined at Lapack_module.f90 lines 470-486
        
        Parameters
        ----------
        mat : float array
        
        """
        _arespy.f90wrap_cholesky_factor_real(mat=mat)
    
    @staticmethod
    def invmat_real(mat):
        """
        invmat_real(mat)
        
        
        Defined at Lapack_module.f90 lines 489-521
        
        Parameters
        ----------
        mat : float array
        
        """
        _arespy.f90wrap_invmat_real(mat=mat)
    
    @property
    def mmax(self):
        """
        Element mmax ftype=integer(i4b) pytype=int
        
        
        Defined at Lapack_module.f90 line 10
        
        """
        return _arespy.f90wrap_lapack_module__get__mmax()
    
    @mmax.setter
    def mmax(self, mmax):
        _arespy.f90wrap_lapack_module__set__mmax(mmax)
    
    @property
    def nmax(self):
        """
        Element nmax ftype=integer(i4b) pytype=int
        
        
        Defined at Lapack_module.f90 line 10
        
        """
        return _arespy.f90wrap_lapack_module__get__nmax()
    
    @nmax.setter
    def nmax(self, nmax):
        _arespy.f90wrap_lapack_module__set__nmax(nmax)
    
    def __str__(self):
        ret = ['<lapack_module>{\n']
        ret.append('    mmax : ')
        ret.append(repr(self.mmax))
        ret.append(',\n    nmax : ')
        ret.append(repr(self.nmax))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

lapack_module = Lapack_Module()

class Math(f90wrap.runtime.FortranModule):
    """
    Module math
    
    
    Defined at Math.f90 lines 1-2409
    
    """
    @staticmethod
    def change_case(instr, str, fun):
        """
        change_case(instr, str, fun)
        
        
        Defined at Math.f90 lines 38-78
        
        Parameters
        ----------
        instr : str
        str : str
        fun : int
        
        """
        _arespy.f90wrap_change_case(instr=instr, str=str, fun=fun)
    
    @staticmethod
    def find_keywords(str, ch_mark, id_key, id_value):
        """
        find_keywords(str, ch_mark, id_key, id_value)
        
        
        Defined at Math.f90 lines 82-111
        
        Parameters
        ----------
        str : str
        ch_mark : str
        id_key : int
        id_value : int
        
        """
        _arespy.f90wrap_find_keywords(str=str, ch_mark=ch_mark, id_key=id_key, \
            id_value=id_value)
    
    @staticmethod
    def find_nword(str, ch_comma, nword):
        """
        find_nword(str, ch_comma, nword)
        
        
        Defined at Math.f90 lines 115-141
        
        Parameters
        ----------
        str : str
        ch_comma : str
        nword : int
        
        """
        _arespy.f90wrap_find_nword(str=str, ch_comma=ch_comma, nword=nword)
    
    @staticmethod
    def det(matrix):
        """
        det = det(matrix)
        
        
        Defined at Math.f90 lines 173-177
        
        Parameters
        ----------
        matrix : float array
        
        Returns
        -------
        det : float
        
        """
        det = _arespy.f90wrap_det(matrix=matrix)
        return det
    
    @staticmethod
    def inv_33(m):
        """
        inv_33 = inv_33(m)
        
        
        Defined at Math.f90 lines 181-211
        
        Parameters
        ----------
        m : float array
        
        Returns
        -------
        inv_33 : float array
        
        """
        inv_33 = _arespy.f90wrap_inv_33(m=m)
        return inv_33
    
    @staticmethod
    def lindg(eta, lambda_, mu):
        """
        lindg = lindg(eta, lambda_, mu)
        
        
        Defined at Math.f90 lines 215-261
        
        Parameters
        ----------
        eta : float
        lambda_ : float
        mu : float
        
        Returns
        -------
        lindg : float
        
        """
        lindg = _arespy.f90wrap_lindg(eta=eta, lambda_=lambda_, mu=mu)
        return lindg
    
    @staticmethod
    def int_to_char(int_bn):
        """
        int_to_char = int_to_char(int_bn)
        
        
        Defined at Math.f90 lines 265-286
        
        Parameters
        ----------
        int_bn : int
        
        Returns
        -------
        int_to_char : str
        
        -----------------------------------------------------------------------
        """
        int_to_char = _arespy.f90wrap_int_to_char(int_bn=int_bn)
        return int_to_char
    
    @staticmethod
    def lat2matrix(lat_para, lat_mat, flag):
        """
        lat2matrix(lat_para, lat_mat, flag)
        
        
        Defined at Math.f90 lines 290-333
        
        Parameters
        ----------
        lat_para : float array
        lat_mat : float array
        flag : int
        
        """
        _arespy.f90wrap_lat2matrix(lat_para=lat_para, lat_mat=lat_mat, flag=flag)
    
    @staticmethod
    def one2three(id, n_dens, pos):
        """
        one2three(id, n_dens, pos)
        
        
        Defined at Math.f90 lines 337-359
        
        Parameters
        ----------
        id : int
        n_dens : int array
        pos : int array
        
        """
        _arespy.f90wrap_one2three(id=id, n_dens=n_dens, pos=pos)
    
    @staticmethod
    def init_random_seed():
        """
        init_random_seed()
        
        
        Defined at Math.f90 lines 453-517
        
        
        """
        _arespy.f90wrap_init_random_seed()
    
    @staticmethod
    def atom_mass(atom_name, mass):
        """
        atom_mass(atom_name, mass)
        
        
        Defined at Math.f90 lines 525-666
        
        Parameters
        ----------
        atom_name : str
        mass : float
        
        --------------------------------------------
        """
        _arespy.f90wrap_atom_mass(atom_name=atom_name, mass=mass)
    
    @staticmethod
    def newton_inter(n, x, y, m, tx, ty):
        """
        newton_inter(n, x, y, m, tx, ty)
        
        
        Defined at Math.f90 lines 670-693
        
        Parameters
        ----------
        n : int
        x : float array
        y : float array
        m : int
        tx : float array
        ty : float array
        
        ----------------------------------------------------
        """
        _arespy.f90wrap_newton_inter(n=n, x=x, y=y, m=m, tx=tx, ty=ty)
    
    @staticmethod
    def diag(n, ina, w, q):
        """
        diag(n, ina, w, q)
        
        
        Defined at Math.f90 lines 697-831
        
        Parameters
        ----------
        n : int
        ina : float array
        w : float array
        q : float array
        
        """
        _arespy.f90wrap_diag(n=n, ina=ina, w=w, q=q)
    
    @staticmethod
    def boltzmann_distribution(rnull, width):
        """
        boltzmann_distribution = boltzmann_distribution(rnull, width)
        
        
        Defined at Math.f90 lines 835-850
        
        Parameters
        ----------
        rnull : float
        width : float
        
        Returns
        -------
        boltzmann_distribution : float
        
        """
        boltzmann_distribution = _arespy.f90wrap_boltzmann_distribution(rnull=rnull, \
            width=width)
        return boltzmann_distribution
    
    @staticmethod
    def dir2car(cry_coo, ort_coo, lat):
        """
        dir2car(cry_coo, ort_coo, lat)
        
        
        Defined at Math.f90 lines 857-867
        
        Parameters
        ----------
        cry_coo : float array
        ort_coo : float array
        lat : float array
        
        """
        _arespy.f90wrap_dir2car(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)
    
    @staticmethod
    def car2dir(ort_coo, cry_coo, lat):
        """
        car2dir(ort_coo, cry_coo, lat)
        
        
        Defined at Math.f90 lines 873-884
        
        Parameters
        ----------
        ort_coo : float array
        cry_coo : float array
        lat : float array
        
        """
        _arespy.f90wrap_car2dir(ort_coo=ort_coo, cry_coo=cry_coo, lat=lat)
    
    @staticmethod
    def thr2mat(n1, n2, n3, i, j, k):
        """
        dimnu = thr2mat(n1, n2, n3, i, j, k)
        
        
        Defined at Math.f90 lines 887-893
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        i : int
        j : int
        k : int
        
        Returns
        -------
        dimnu : int
        
        """
        dimnu = _arespy.f90wrap_thr2mat(n1=n1, n2=n2, n3=n3, i=i, j=j, k=k)
        return dimnu
    
    @staticmethod
    def mat2thr(n1, n2, n3, i, offset=None):
        """
        ix, iy, iz = mat2thr(n1, n2, n3, i[, offset])
        
        
        Defined at Math.f90 lines 896-934
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        i : int
        offset : int array
        
        Returns
        -------
        ix : int
        iy : int
        iz : int
        
        ------------------------
        """
        ix, iy, iz = _arespy.f90wrap_mat2thr(n1=n1, n2=n2, n3=n3, i=i, offset=offset)
        return ix, iy, iz
    
    @staticmethod
    def sopo(a, lda, n, b, ldb, m, w):
        """
        sopo(a, lda, n, b, ldb, m, w)
        
        
        Defined at Math.f90 lines 937-992
        
        Parameters
        ----------
        a : float array
        lda : int
        n : int
        b : float array
        ldb : int
        m : int
        w : float array
        
        """
        _arespy.f90wrap_sopo(a=a, lda=lda, n=n, b=b, ldb=ldb, m=m, w=w)
    
    @staticmethod
    def csort_eigen(nev, arr, brr):
        """
        csort_eigen(nev, arr, brr)
        
        
        Defined at Math.f90 lines 995-1033
        
        Parameters
        ----------
        nev : int
        arr : float array
        brr : complex array
        
        """
        _arespy.f90wrap_csort_eigen(nev=nev, arr=arr, brr=brr)
    
    @staticmethod
    def rsort_eigen(nev, arr, brr):
        """
        rsort_eigen(nev, arr, brr)
        
        
        Defined at Math.f90 lines 1036-1070
        
        Parameters
        ----------
        nev : int
        arr : float array
        brr : float array
        
        """
        _arespy.f90wrap_rsort_eigen(nev=nev, arr=arr, brr=brr)
    
    @staticmethod
    def realint_sort(nev, arr, brr, crr=None):
        """
        realint_sort(nev, arr, brr[, crr])
        
        
        Defined at Math.f90 lines 1073-1112
        
        Parameters
        ----------
        nev : int
        arr : float array
        brr : int array
        crr : float array
        
        """
        _arespy.f90wrap_realint_sort(nev=nev, arr=arr, brr=brr, crr=crr)
    
    @staticmethod
    def sort_eigval(n, arr):
        """
        sort_eigval(n, arr)
        
        
        Defined at Math.f90 lines 1115-1137
        
        Parameters
        ----------
        n : int
        arr : float array
        
        """
        _arespy.f90wrap_sort_eigval(n=n, arr=arr)
    
    @staticmethod
    def pgfo(omat, n1, n2, n3, ix, iy, iz):
        """
        ex, ey, ez = pgfo(omat, n1, n2, n3, ix, iy, iz)
        
        
        Defined at Math.f90 lines 1140-1170
        
        Parameters
        ----------
        omat : float array
        n1 : int
        n2 : int
        n3 : int
        ix : int
        iy : int
        iz : int
        
        Returns
        -------
        ex : int
        ey : int
        ez : int
        
        """
        ex, ey, ez = _arespy.f90wrap_pgfo(omat=omat, n1=n1, n2=n2, n3=n3, ix=ix, iy=iy, \
            iz=iz)
        return ex, ey, ez
    
    @staticmethod
    def cubicsplineinterp(fun, ddfdx2, xmax, dx, x, zion=None):
        """
        cubicsplineinterp = cubicsplineinterp(fun, ddfdx2, xmax, dx, x[, zion])
        
        
        Defined at Math.f90 lines 1173-1249
        
        Parameters
        ----------
        fun : float array
        ddfdx2 : float array
        xmax : float
        dx : float
        x : float
        zion : float
        
        Returns
        -------
        cubicsplineinterp : float
        
        """
        cubicsplineinterp = _arespy.f90wrap_cubicsplineinterp(fun=fun, ddfdx2=ddfdx2, \
            xmax=xmax, dx=dx, x=x, zion=zion)
        return cubicsplineinterp
    
    @staticmethod
    def finite_factor(fnor, norder, coe):
        """
        finite_factor(fnor, norder, coe)
        
        
        Defined at Math.f90 lines 1252-1393
        
        Parameters
        ----------
        fnor : int
        norder : int
        coe : float array
        
        """
        _arespy.f90wrap_finite_factor(fnor=fnor, norder=norder, coe=coe)
    
    @staticmethod
    def finite_factor_new(fnor, norder, coe):
        """
        finite_factor_new(fnor, norder, coe)
        
        
        Defined at Math.f90 lines 1396-1453
        
        Parameters
        ----------
        fnor : int
        norder : int
        coe : float array
        
        """
        _arespy.f90wrap_finite_factor_new(fnor=fnor, norder=norder, coe=coe)
    
    @staticmethod
    def dfdr(np, h, f, df):
        """
        dfdr(np, h, f, df)
        
        
        Defined at Math.f90 lines 1456-1490
        
        Parameters
        ----------
        np : int
        h : float
        f : float array
        df : float array
        
        """
        _arespy.f90wrap_dfdr(np=np, h=h, f=f, df=df)
    
    @staticmethod
    def cubichermiteinterp(fun, dfdx, xmax, h, x):
        """
        cubichermiteinterp = cubichermiteinterp(fun, dfdx, xmax, h, x)
        
        
        Defined at Math.f90 lines 1493-1552
        
        Parameters
        ----------
        fun : float array
        dfdx : float array
        xmax : float
        h : float
        x : float
        
        Returns
        -------
        cubichermiteinterp : float
        
        """
        cubichermiteinterp = _arespy.f90wrap_cubichermiteinterp(fun=fun, dfdx=dfdx, \
            xmax=xmax, h=h, x=x)
        return cubichermiteinterp
    
    @staticmethod
    def simpleinterp(fun, xmax, h, x):
        """
        simpleinterp = simpleinterp(fun, xmax, h, x)
        
        
        Defined at Math.f90 lines 1555-1594
        
        Parameters
        ----------
        fun : float array
        xmax : float
        h : float
        x : float
        
        Returns
        -------
        simpleinterp : float
        
        """
        simpleinterp = _arespy.f90wrap_simpleinterp(fun=fun, xmax=xmax, h=h, x=x)
        return simpleinterp
    
    @staticmethod
    def r_dylm(l, m, x, y, z, rmod, f):
        """
        r_dylm(l, m, x, y, z, rmod, f)
        
        
        Defined at Math.f90 lines 1597-1737
        
        Parameters
        ----------
        l : int
        m : int
        x : float
        y : float
        z : float
        rmod : float
        f : float array
        
        """
        _arespy.f90wrap_r_dylm(l=l, m=m, x=x, y=y, z=z, rmod=rmod, f=f)
    
    @staticmethod
    def rcsntable(n1, n2, n3, dn, h, srmax, rcs, cns, sphindx):
        """
        nspt = rcsntable(n1, n2, n3, dn, h, srmax, rcs, cns, sphindx)
        
        
        Defined at Math.f90 lines 1740-1799
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        dn : int
        h : float
        srmax : float
        rcs : float array
        cns : complex array
        sphindx : bool array
        
        Returns
        -------
        nspt : int
        
        """
        nspt = _arespy.f90wrap_rcsntable(n1=n1, n2=n2, n3=n3, dn=dn, h=h, srmax=srmax, \
            rcs=rcs, cns=cns, sphindx=sphindx)
        return nspt
    
    @staticmethod
    def rcsntable_atoms(n1, n2, n3, dn, na, h, poscar, srmax, atomr, rvec, rcs, cns, \
        sphindx):
        """
        nspt = rcsntable_atoms(n1, n2, n3, dn, na, h, poscar, srmax, atomr, rvec, rcs, \
            cns, sphindx)
        
        
        Defined at Math.f90 lines 1802-1881
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        dn : int
        na : int
        h : float
        poscar : float array
        srmax : float
        atomr : float
        rvec : float array
        rcs : float array
        cns : complex array
        sphindx : bool array
        
        Returns
        -------
        nspt : int
        
        """
        nspt = _arespy.f90wrap_rcsntable_atoms(n1=n1, n2=n2, n3=n3, dn=dn, na=na, h=h, \
            poscar=poscar, srmax=srmax, atomr=atomr, rvec=rvec, rcs=rcs, cns=cns, \
            sphindx=sphindx)
        return nspt
    
    @staticmethod
    def car2spe(orig, h, x, y, z):
        """
        r, cost, sint, cosp, sinp = car2spe(orig, h, x, y, z)
        
        
        Defined at Math.f90 lines 1884-1915
        
        Parameters
        ----------
        orig : float array
        h : float
        x : int
        y : int
        z : int
        
        Returns
        -------
        r : float
        cost : float
        sint : float
        cosp : float
        sinp : float
        
        """
        r, cost, sint, cosp, sinp = _arespy.f90wrap_car2spe(orig=orig, h=h, x=x, y=y, \
            z=z)
        return r, cost, sint, cosp, sinp
    
    @staticmethod
    def calclm(lmax, clm):
        """
        calclm(lmax, clm)
        
        
        Defined at Math.f90 lines 1918-1931
        
        Parameters
        ----------
        lmax : int
        clm : float array
        
        """
        _arespy.f90wrap_calclm(lmax=lmax, clm=clm)
    
    @staticmethod
    def c(l, m):
        """
        c = c(l, m)
        
        
        Defined at Math.f90 lines 1934-1947
        
        Parameters
        ----------
        l : int
        m : int
        
        Returns
        -------
        c : float
        
        """
        c = _arespy.f90wrap_c(l=l, m=m)
        return c
    
    @staticmethod
    def cal_plm(lmax, x, sx, plm):
        """
        cal_plm(lmax, x, sx, plm)
        
        
        Defined at Math.f90 lines 1950-1977
        
        Parameters
        ----------
        lmax : int
        x : float
        sx : float
        plm : float array
        
        """
        _arespy.f90wrap_cal_plm(lmax=lmax, x=x, sx=sx, plm=plm)
    
    @staticmethod
    def interp(np, f, r, rnorm, z=None):
        """
        interp = interp(np, f, r, rnorm[, z])
        
        
        Defined at Math.f90 lines 1980-2033
        
        Parameters
        ----------
        np : int
        f : float array
        r : float array
        rnorm : float
        z : float
        
        Returns
        -------
        interp : float
        
        """
        interp = _arespy.f90wrap_interp(np=np, f=f, r=r, rnorm=rnorm, z=z)
        return interp
    
    @staticmethod
    def interp_dnf(n, np, f, r, rnorm, z=None):
        """
        interp_dnf = interp_dnf(n, np, f, r, rnorm[, z])
        
        
        Defined at Math.f90 lines 2036-2097
        
        Parameters
        ----------
        n : int
        np : int
        f : float array
        r : float array
        rnorm : float
        z : float
        
        Returns
        -------
        interp_dnf : float
        
        """
        interp_dnf = _arespy.f90wrap_interp_dnf(n=n, np=np, f=f, r=r, rnorm=rnorm, z=z)
        return interp_dnf
    
    @staticmethod
    def polynom(m, np, xa, ya, c, x):
        """
        polynom = polynom(m, np, xa, ya, c, x)
        
        
        Defined at Math.f90 lines 2100-2280
        
        Parameters
        ----------
        m : int
        np : int
        xa : float array
        ya : float array
        c : float array
        x : float
        
        Returns
        -------
        polynom : float
        
        """
        polynom = _arespy.f90wrap_polynom(m=m, np=np, xa=xa, ya=ya, c=c, x=x)
        return polynom
    
    @staticmethod
    def fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt):
        """
        fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt)
        
        
        Defined at Math.f90 lines 2283-2341
        
        Parameters
        ----------
        nr : int
        rr : float array
        rab : float array
        vr : float array
        ll : int
        nql : int
        yp : float array
        vql : float array
        vt : float
        
        """
        _arespy.f90wrap_fourier_1d(nr=nr, rr=rr, rab=rab, vr=vr, ll=ll, nql=nql, yp=yp, \
            vql=vql, vt=vt)
    
    @staticmethod
    def integ_new(rab, y):
        """
        f = integ_new(rab, y)
        
        
        Defined at Math.f90 lines 2344-2364
        
        Parameters
        ----------
        rab : float array
        y : float array
        
        Returns
        -------
        f : float
        
        """
        f = _arespy.f90wrap_integ_new(rab=rab, y=y)
        return f
    
    @staticmethod
    def sphbess(l, x):
        """
        sphbess = sphbess(l, x)
        
        
        Defined at Math.f90 lines 2367-2406
        
        Parameters
        ----------
        l : int
        x : float
        
        Returns
        -------
        sphbess : float
        
        """
        sphbess = _arespy.f90wrap_sphbess(l=l, x=x)
        return sphbess
    
    @staticmethod
    def _norm_real(a):
        """
        norm_real = _norm_real(a)
        
        
        Defined at Math.f90 lines 145-147
        
        Parameters
        ----------
        a : float array
        
        Returns
        -------
        norm_real : float
        
        """
        norm_real = _arespy.f90wrap_norm_real(a=a)
        return norm_real
    
    @staticmethod
    def _norm_complex(a):
        """
        norm_complex = _norm_complex(a)
        
        
        Defined at Math.f90 lines 151-153
        
        Parameters
        ----------
        a : complex array
        
        Returns
        -------
        norm_complex : float
        
        """
        norm_complex = _arespy.f90wrap_norm_complex(a=a)
        return norm_complex
    
    @staticmethod
    def norm(*args, **kwargs):
        """
        norm(*args, **kwargs)
        
        
        Defined at Math.f90 lines 22-23
        
        Overloaded interface containing the following procedures:
          _norm_real
          _norm_complex
        
        """
        for proc in [Math._norm_real, Math._norm_complex]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _cross_real(a, b):
        """
        cross_real = _cross_real(a, b)
        
        
        Defined at Math.f90 lines 165-169
        
        Parameters
        ----------
        a : float array
        b : float array
        
        Returns
        -------
        cross_real : float array
        
        """
        cross_real = _arespy.f90wrap_cross_real(a=a, b=b)
        return cross_real
    
    @staticmethod
    def _cross_complex(a, b):
        """
        cross_complex = _cross_complex(a, b)
        
        
        Defined at Math.f90 lines 157-161
        
        Parameters
        ----------
        a : complex array
        b : complex array
        
        Returns
        -------
        cross_complex : complex array
        
        """
        cross_complex = _arespy.f90wrap_cross_complex(a=a, b=b)
        return cross_complex
    
    @staticmethod
    def cross(*args, **kwargs):
        """
        cross(*args, **kwargs)
        
        
        Defined at Math.f90 lines 27-28
        
        Overloaded interface containing the following procedures:
          _cross_real
          _cross_complex
        
        """
        for proc in [Math._cross_real, Math._cross_complex]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _integer_index(array, n, id):
        """
        _integer_index(array, n, id)
        
        
        Defined at Math.f90 lines 408-449
        
        Parameters
        ----------
        array : int array
        n : int
        id : int array
        
        """
        _arespy.f90wrap_integer_index(array=array, n=n, id=id)
    
    @staticmethod
    def _real_index(array, n, id):
        """
        _real_index(array, n, id)
        
        
        Defined at Math.f90 lines 363-404
        
        Parameters
        ----------
        array : float array
        n : int
        id : int array
        
        """
        _arespy.f90wrap_real_index(array=array, n=n, id=id)
    
    @staticmethod
    def sort_id(*args, **kwargs):
        """
        sort_id(*args, **kwargs)
        
        
        Defined at Math.f90 lines 32-33
        
        Overloaded interface containing the following procedures:
          _integer_index
          _real_index
        
        """
        for proc in [Math._integer_index, Math._real_index]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    _dt_array_initialisers = []
    

math = Math()

class Mathsplines(f90wrap.runtime.FortranModule):
    """
    Module mathsplines
    
    
    Defined at MathSplines.f90 lines 1-832
    
    """
    @staticmethod
    def spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp):
        """
        spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp)
        
        
        Defined at MathSplines.f90 lines 174-383
        
        Parameters
        ----------
        n : int
        t : float array
        y : float array
        ibcbeg : int
        ybcbeg : float
        ibcend : int
        ybcend : float
        ypp : float array
        
        """
        _arespy.f90wrap_spline_cubic_set(n=n, t=t, y=y, ibcbeg=ibcbeg, ybcbeg=ybcbeg, \
            ibcend=ibcend, ybcend=ybcend, ypp=ypp)
    
    @staticmethod
    def spline_cubic_val(n, t, y, ypp, tval):
        """
        yval, ypval, yppval = spline_cubic_val(n, t, y, ypp, tval)
        
        
        Defined at MathSplines.f90 lines 385-479
        
        Parameters
        ----------
        n : int
        t : float array
        y : float array
        ypp : float array
        tval : float
        
        Returns
        -------
        yval : float
        ypval : float
        yppval : float
        
        """
        yval, ypval, yppval = _arespy.f90wrap_spline_cubic_val(n=n, t=t, y=y, ypp=ypp, \
            tval=tval)
        return yval, ypval, yppval
    
    @staticmethod
    def rvec_bracket(n, x, xval):
        """
        left, right = rvec_bracket(n, x, xval)
        
        
        Defined at MathSplines.f90 lines 481-577
        
        Parameters
        ----------
        n : int
        x : float array
        xval : float
        
        Returns
        -------
        left : int
        right : int
        
        """
        left, right = _arespy.f90wrap_rvec_bracket(n=n, x=x, xval=xval)
        return left, right
    
    @staticmethod
    def s3_fs(a1, a2, a3, n, b, x):
        """
        s3_fs(a1, a2, a3, n, b, x)
        
        
        Defined at MathSplines.f90 lines 579-646
        
        Parameters
        ----------
        a1 : float array
        a2 : float array
        a3 : float array
        n : int
        b : float array
        x : float array
        
        """
        _arespy.f90wrap_s3_fs(a1=a1, a2=a2, a3=a3, n=n, b=b, x=x)
    
    @staticmethod
    def polynom(m, np, xa, ya, c, x):
        """
        polynom = polynom(m, np, xa, ya, c, x)
        
        
        Defined at MathSplines.f90 lines 654-832
        
        Parameters
        ----------
        m : int
        np : int
        xa : float array
        ya : float array
        c : float array
        x : float
        
        Returns
        -------
        polynom : float
        
        """
        polynom = _arespy.f90wrap_polynom(m=m, np=np, xa=xa, ya=ya, c=c, x=x)
        return polynom
    
    _dt_array_initialisers = []
    

mathsplines = Mathsplines()

class Matvec_Module(f90wrap.runtime.FortranModule):
    """
    Module matvec_module
    
    
    Defined at Matvec_module.f90 lines 1-223
    
    """
    @staticmethod
    def cmplx_matvec(ik, veff1d, p, q, dimen):
        """
        cmplx_matvec(ik, veff1d, p, q, dimen)
        
        
        Defined at Matvec_module.f90 lines 8-26
        
        Parameters
        ----------
        ik : int
        veff1d : float array
        p : complex array
        q : complex array
        dimen : int
        
        """
        _arespy.f90wrap_cmplx_matvec(ik=ik, veff1d=veff1d, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def cmplx_nlocmatvec(ik, p, q):
        """
        cmplx_nlocmatvec(ik, p, q)
        
        
        Defined at Matvec_module.f90 lines 29-107
        
        Parameters
        ----------
        ik : int
        p : complex array
        q : complex array
        
        """
        _arespy.f90wrap_cmplx_nlocmatvec(ik=ik, p=p, q=q)
    
    @staticmethod
    def real_matvec(veff1d, p, q, dimen):
        """
        real_matvec(veff1d, p, q, dimen)
        
        
        Defined at Matvec_module.f90 lines 110-130
        
        Parameters
        ----------
        veff1d : float array
        p : float array
        q : float array
        dimen : int
        
        """
        _arespy.f90wrap_real_matvec(veff1d=veff1d, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def real_nlocmatvec(p, q):
        """
        real_nlocmatvec(p, q)
        
        
        Defined at Matvec_module.f90 lines 133-211
        
        Parameters
        ----------
        p : float array
        q : float array
        
        """
        _arespy.f90wrap_real_nlocmatvec(p=p, q=q)
    
    @staticmethod
    def real_matvec_m(mat, p, q, dimen):
        """
        real_matvec_m(mat, p, q, dimen)
        
        
        Defined at Matvec_module.f90 lines 214-222
        
        Parameters
        ----------
        mat : float array
        p : float array
        q : float array
        dimen : int
        
        """
        _arespy.f90wrap_real_matvec_m(mat=mat, p=p, q=q, dimen=dimen)
    
    _dt_array_initialisers = []
    

matvec_module = Matvec_Module()

class Mixer_Module(f90wrap.runtime.FortranModule):
    """
    Module mixer_module
    
    
    Defined at Mixer_module.f90 lines 1-806
    
    """
    @f90wrap.runtime.register_class("arespy.mixer_data")
    class mixer_data(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=mixer_data)
        
        
        Defined at Mixer_module.f90 lines 14-23
        
        """
        def __init__(self, handle=None):
            """
            self = Mixer_Data()
            
            
            Defined at Mixer_module.f90 lines 14-23
            
            
            Returns
            -------
            this : Mixer_Data
            	Object to be constructed
            
            
            Automatically generated constructor for mixer_data
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_mixer_data_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Mixer_Data
            
            
            Defined at Mixer_module.f90 lines 14-23
            
            Parameters
            ----------
            this : Mixer_Data
            	Object to be destructed
            
            
            Automatically generated destructor for mixer_data
            """
            if self._alloc:
                _arespy.f90wrap_mixer_data_finalise(this=self._handle)
        
        @property
        def dxl(self):
            """
            Element dxl ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.f90 line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_mixer_data__array__dxl(self._handle)
            if array_handle in self._arrays:
                dxl = self._arrays[array_handle]
            else:
                dxl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_mixer_data__array__dxl)
                self._arrays[array_handle] = dxl
            return dxl
        
        @dxl.setter
        def dxl(self, dxl):
            self.dxl[...] = dxl
        
        @property
        def dfl(self):
            """
            Element dfl ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.f90 line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_mixer_data__array__dfl(self._handle)
            if array_handle in self._arrays:
                dfl = self._arrays[array_handle]
            else:
                dfl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_mixer_data__array__dfl)
                self._arrays[array_handle] = dfl
            return dfl
        
        @dfl.setter
        def dfl(self, dfl):
            self.dfl[...] = dfl
        
        @property
        def voma(self):
            """
            Element voma ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.f90 line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_mixer_data__array__voma(self._handle)
            if array_handle in self._arrays:
                voma = self._arrays[array_handle]
            else:
                voma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_mixer_data__array__voma)
                self._arrays[array_handle] = voma
            return voma
        
        @voma.setter
        def voma(self, voma):
            self.voma[...] = voma
        
        @property
        def kerker(self):
            """
            Element kerker ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.f90 line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_mixer_data__array__kerker(self._handle)
            if array_handle in self._arrays:
                kerker = self._arrays[array_handle]
            else:
                kerker = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_mixer_data__array__kerker)
                self._arrays[array_handle] = kerker
            return kerker
        
        @kerker.setter
        def kerker(self, kerker):
            self.kerker[...] = kerker
        
        @property
        def dxgl(self):
            """
            Element dxgl ftype=complex(dcp) pytype=complex
            
            
            Defined at Mixer_module.f90 line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_mixer_data__array__dxgl(self._handle)
            if array_handle in self._arrays:
                dxgl = self._arrays[array_handle]
            else:
                dxgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_mixer_data__array__dxgl)
                self._arrays[array_handle] = dxgl
            return dxgl
        
        @dxgl.setter
        def dxgl(self, dxgl):
            self.dxgl[...] = dxgl
        
        @property
        def drgl(self):
            """
            Element drgl ftype=complex(dcp) pytype=complex
            
            
            Defined at Mixer_module.f90 line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_mixer_data__array__drgl(self._handle)
            if array_handle in self._arrays:
                drgl = self._arrays[array_handle]
            else:
                drgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_mixer_data__array__drgl)
                self._arrays[array_handle] = drgl
            return drgl
        
        @drgl.setter
        def drgl(self, drgl):
            self.drgl[...] = drgl
        
        def __str__(self):
            ret = ['<mixer_data>{\n']
            ret.append('    dxl : ')
            ret.append(repr(self.dxl))
            ret.append(',\n    dfl : ')
            ret.append(repr(self.dfl))
            ret.append(',\n    voma : ')
            ret.append(repr(self.voma))
            ret.append(',\n    kerker : ')
            ret.append(repr(self.kerker))
            ret.append(',\n    dxgl : ')
            ret.append(repr(self.dxgl))
            ret.append(',\n    drgl : ')
            ret.append(repr(self.drgl))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init_mixer(nps):
        """
        init_mixer(nps)
        
        
        Defined at Mixer_module.f90 lines 31-85
        
        Parameters
        ----------
        nps : int
        
        """
        _arespy.f90wrap_init_mixer(nps=nps)
    
    @staticmethod
    def destroy_mixer():
        """
        destroy_mixer()
        
        
        Defined at Mixer_module.f90 lines 88-99
        
        
        """
        _arespy.f90wrap_destroy_mixer()
    
    @staticmethod
    def mixing(iter, xout, xin):
        """
        res = mixing(iter, xout, xin)
        
        
        Defined at Mixer_module.f90 lines 102-215
        
        Parameters
        ----------
        iter : int
        xout : float array
        xin : float array
        
        Returns
        -------
        res : float
        
        """
        res = _arespy.f90wrap_mixing(iter=iter, xout=xout, xin=xin)
        return res
    
    @staticmethod
    def anderson_mixing(iiter, xin, xout):
        """
        err = anderson_mixing(iiter, xin, xout)
        
        
        Defined at Mixer_module.f90 lines 221-297
        
        Parameters
        ----------
        iiter : int
        xin : float array
        xout : float array
        
        Returns
        -------
        err : float
        
        """
        err = _arespy.f90wrap_anderson_mixing(iiter=iiter, xin=xin, xout=xout)
        return err
    
    @staticmethod
    def om1c(nam, nuh, sp, dfp, voma):
        """
        om1c(nam, nuh, sp, dfp, voma)
        
        
        Defined at Mixer_module.f90 lines 300-359
        
        Parameters
        ----------
        nam : int
        nuh : int
        sp : float
        dfp : float
        voma : float
        
        """
        _arespy.f90wrap_om1c(nam=nam, nuh=nuh, sp=sp, dfp=dfp, voma=voma)
    
    @staticmethod
    def amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn):
        """
        amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn)
        
        
        Defined at Mixer_module.f90 lines 362-477
        
        Parameters
        ----------
        beta : float
        w0 : float
        nam : int
        nuh : int
        dxp : float
        dfp : float
        sp : float
        xl : float
        fl : float
        voma : float
        xn : float
        
        """
        _arespy.f90wrap_amst(beta=beta, w0=w0, nam=nam, nuh=nuh, dxp=dxp, dfp=dfp, \
            sp=sp, xl=xl, fl=fl, voma=voma, xn=xn)
    
    @staticmethod
    def init_kerker():
        """
        init_kerker()
        
        
        Defined at Mixer_module.f90 lines 480-512
        
        
        """
        _arespy.f90wrap_init_kerker()
    
    @staticmethod
    def rpulayk_mixing(iter, rlg, xing):
        """
        rpulayk_mixing(iter, rlg, xing)
        
        
        Defined at Mixer_module.f90 lines 518-592
        
        Parameters
        ----------
        iter : int
        rlg : complex array
        xing : complex array
        
        ----------------------
        store RL
        """
        _arespy.f90wrap_rpulayk_mixing(iter=iter, rlg=rlg, xing=xing)
    
    @staticmethod
    def rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn):
        """
        rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn)
        
        
        Defined at Mixer_module.f90 lines 595-672
        
        Parameters
        ----------
        beta : float
        w0 : float
        dime : int
        nh : int
        dxl : complex array
        drl : complex array
        xl : complex array
        rl : complex array
        xn : complex array
        
        """
        _arespy.f90wrap_rpulayk_mix(beta=beta, w0=w0, dime=dime, nh=nh, dxl=dxl, \
            drl=drl, xl=xl, rl=rl, xn=xn)
    
    @staticmethod
    def resta_mixing(iter, rlg, xing):
        """
        resta_mixing(iter, rlg, xing)
        
        
        Defined at Mixer_module.f90 lines 677-749
        
        Parameters
        ----------
        iter : int
        rlg : complex array
        xing : complex array
        
        ----------------------
        store RL
        """
        _arespy.f90wrap_resta_mixing(iter=iter, rlg=rlg, xing=xing)
    
    @staticmethod
    def init_resta():
        """
        init_resta()
        
        
        Defined at Mixer_module.f90 lines 752-803
        
        
        """
        _arespy.f90wrap_init_resta()
    
    @property
    def nam(self):
        """
        Element nam ftype=integer(i4b) pytype=int
        
        
        Defined at Mixer_module.f90 line 26
        
        """
        return _arespy.f90wrap_mixer_module__get__nam()
    
    @nam.setter
    def nam(self, nam):
        _arespy.f90wrap_mixer_module__set__nam(nam)
    
    @property
    def nuh(self):
        """
        Element nuh ftype=integer(i4b) pytype=int
        
        
        Defined at Mixer_module.f90 line 26
        
        """
        return _arespy.f90wrap_mixer_module__get__nuh()
    
    @nuh.setter
    def nuh(self, nuh):
        _arespy.f90wrap_mixer_module__set__nuh(nuh)
    
    def __str__(self):
        ret = ['<mixer_module>{\n']
        ret.append('    nam : ')
        ret.append(repr(self.nam))
        ret.append(',\n    nuh : ')
        ret.append(repr(self.nuh))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

mixer_module = Mixer_Module()

class Array_Io(f90wrap.runtime.FortranModule):
    """
    Module array_io
    
    
    Defined at MPI_array_io.f90 lines 1-92
    
    """
    @f90wrap.runtime.register_class("arespy.out_label")
    class out_label(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=out_label)
        
        
        Defined at MPI_array_io.f90 lines 4-6
        
        """
        def __init__(self, handle=None):
            """
            self = Out_Label()
            
            
            Defined at MPI_array_io.f90 lines 4-6
            
            
            Returns
            -------
            this : Out_Label
            	Object to be constructed
            
            
            Automatically generated constructor for out_label
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_out_label_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Out_Label
            
            
            Defined at MPI_array_io.f90 lines 4-6
            
            Parameters
            ----------
            this : Out_Label
            	Object to be destructed
            
            
            Automatically generated destructor for out_label
            """
            if self._alloc:
                _arespy.f90wrap_out_label_finalise(this=self._handle)
        
        @property
        def label(self):
            """
            Element label ftype=character(len=100) pytype=str
            
            
            Defined at MPI_array_io.f90 line 5
            
            """
            return _arespy.f90wrap_out_label__get__label(self._handle)
        
        @label.setter
        def label(self, label):
            _arespy.f90wrap_out_label__set__label(self._handle, label)
        
        @property
        def num(self):
            """
            Element num ftype=integer(i4b) pytype=int
            
            
            Defined at MPI_array_io.f90 line 6
            
            """
            return _arespy.f90wrap_out_label__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _arespy.f90wrap_out_label__set__num(self._handle, num)
        
        def __str__(self):
            ret = ['<out_label>{\n']
            ret.append('    label : ')
            ret.append(repr(self.label))
            ret.append(',\n    num : ')
            ret.append(repr(self.num))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init_outfile():
        """
        init_outfile()
        
        
        Defined at MPI_array_io.f90 lines 17-23
        
        
        """
        _arespy.f90wrap_init_outfile()
    
    @staticmethod
    def _output_r(size_bn, array, name):
        """
        _output_r(size_bn, array, name)
        
        
        Defined at MPI_array_io.f90 lines 25-64
        
        Parameters
        ----------
        size_bn : int
        array : float array
        name : str
        
        """
        _arespy.f90wrap_output_r(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def _output_i(size_bn, array, name):
        """
        _output_i(size_bn, array, name)
        
        
        Defined at MPI_array_io.f90 lines 75-82
        
        Parameters
        ----------
        size_bn : int
        array : int array
        name : str
        
        """
        _arespy.f90wrap_output_i(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def output(*args, **kwargs):
        """
        output(*args, **kwargs)
        
        
        Defined at MPI_array_io.f90 lines 10-11
        
        Overloaded interface containing the following procedures:
          _output_r
          _output_i
        
        """
        for proc in [Array_Io._output_r, Array_Io._output_i]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _input_r(size_bn, array, name):
        """
        _input_r(size_bn, array, name)
        
        
        Defined at MPI_array_io.f90 lines 66-73
        
        Parameters
        ----------
        size_bn : int
        array : float array
        name : str
        
        """
        _arespy.f90wrap_input_r(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def _input_i(size_bn, array, name):
        """
        _input_i(size_bn, array, name)
        
        
        Defined at MPI_array_io.f90 lines 84-91
        
        Parameters
        ----------
        size_bn : int
        array : int array
        name : str
        
        """
        _arespy.f90wrap_input_i(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def input(*args, **kwargs):
        """
        input(*args, **kwargs)
        
        
        Defined at MPI_array_io.f90 lines 13-14
        
        Overloaded interface containing the following procedures:
          _input_r
          _input_i
        
        """
        for proc in [Array_Io._input_r, Array_Io._input_i]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    _dt_array_initialisers = []
    

array_io = Array_Io()

class M_Time_Evaluate(f90wrap.runtime.FortranModule):
    """
    Module m_time_evaluate
    
    
    Defined at MPI_time_evaluate.f90 lines 1-198
    
    """
    @f90wrap.runtime.register_class("arespy.time_record")
    class time_record(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=time_record)
        
        
        Defined at MPI_time_evaluate.f90 lines 11-14
        
        """
        def __init__(self, handle=None):
            """
            self = Time_Record()
            
            
            Defined at MPI_time_evaluate.f90 lines 11-14
            
            
            Returns
            -------
            this : Time_Record
            	Object to be constructed
            
            
            Automatically generated constructor for time_record
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_time_record_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Time_Record
            
            
            Defined at MPI_time_evaluate.f90 lines 11-14
            
            Parameters
            ----------
            this : Time_Record
            	Object to be destructed
            
            
            Automatically generated destructor for time_record
            """
            if self._alloc:
                _arespy.f90wrap_time_record_finalise(this=self._handle)
        
        @property
        def str(self):
            """
            Element str ftype=character(len=100) pytype=str
            
            
            Defined at MPI_time_evaluate.f90 line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_time_record__array__str(self._handle)
            if array_handle in self._arrays:
                str = self._arrays[array_handle]
            else:
                str = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_time_record__array__str)
                self._arrays[array_handle] = str
            return str
        
        @str.setter
        def str(self, str):
            self.str[...] = str
        
        @property
        def t1(self):
            """
            Element t1 ftype=integer(i4b) pytype=int
            
            
            Defined at MPI_time_evaluate.f90 line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_time_record__array__t1(self._handle)
            if array_handle in self._arrays:
                t1 = self._arrays[array_handle]
            else:
                t1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_time_record__array__t1)
                self._arrays[array_handle] = t1
            return t1
        
        @t1.setter
        def t1(self, t1):
            self.t1[...] = t1
        
        @property
        def t2(self):
            """
            Element t2 ftype=integer(i4b) pytype=int
            
            
            Defined at MPI_time_evaluate.f90 line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_time_record__array__t2(self._handle)
            if array_handle in self._arrays:
                t2 = self._arrays[array_handle]
            else:
                t2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_time_record__array__t2)
                self._arrays[array_handle] = t2
            return t2
        
        @t2.setter
        def t2(self, t2):
            self.t2[...] = t2
        
        @property
        def length(self):
            """
            Element length ftype=integer(i4b) pytype=int
            
            
            Defined at MPI_time_evaluate.f90 line 14
            
            """
            return _arespy.f90wrap_time_record__get__length(self._handle)
        
        @length.setter
        def length(self, length):
            _arespy.f90wrap_time_record__set__length(self._handle, length)
        
        def __str__(self):
            ret = ['<time_record>{\n']
            ret.append('    str : ')
            ret.append(repr(self.str))
            ret.append(',\n    t1 : ')
            ret.append(repr(self.t1))
            ret.append(',\n    t2 : ')
            ret.append(repr(self.t2))
            ret.append(',\n    length : ')
            ret.append(repr(self.length))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def record(str, t):
        """
        record(str, t)
        
        
        Defined at MPI_time_evaluate.f90 lines 21-46
        
        Parameters
        ----------
        str : str
        t : int
        
        """
        _arespy.f90wrap_record(str=str, t=t)
    
    @staticmethod
    def reset_data(myflag, str, t):
        """
        reset_data(myflag, str, t)
        
        
        Defined at MPI_time_evaluate.f90 lines 48-77
        
        Parameters
        ----------
        myflag : bool
        str : str
        t : int
        
        """
        _arespy.f90wrap_reset_data(myflag=myflag, str=str, t=t)
    
    @staticmethod
    def init_recrod():
        """
        init_recrod()
        
        
        Defined at MPI_time_evaluate.f90 lines 79-84
        
        
        """
        _arespy.f90wrap_init_recrod()
    
    @staticmethod
    def time_start(str, flag):
        """
        time_start(str, flag)
        
        
        Defined at MPI_time_evaluate.f90 lines 86-96
        
        Parameters
        ----------
        str : str
        flag : bool
        
        """
        _arespy.f90wrap_time_start(str=str, flag=flag)
    
    @staticmethod
    def time_end(str, flag):
        """
        time_end(str, flag)
        
        
        Defined at MPI_time_evaluate.f90 lines 98-108
        
        Parameters
        ----------
        str : str
        flag : bool
        
        """
        _arespy.f90wrap_time_end(str=str, flag=flag)
    
    @staticmethod
    def time_output(str, flag):
        """
        time_output(str, flag)
        
        
        Defined at MPI_time_evaluate.f90 lines 110-125
        
        Parameters
        ----------
        str : str
        flag : bool
        
        """
        _arespy.f90wrap_time_output(str=str, flag=flag)
    
    @staticmethod
    def memory_sum(label, memory_size):
        """
        memory_sum(label, memory_size)
        
        
        Defined at MPI_time_evaluate.f90 lines 127-155
        
        Parameters
        ----------
        label : str
        memory_size : float
        
        """
        _arespy.f90wrap_memory_sum(label=label, memory_size=memory_size)
    
    @staticmethod
    def memory_free(label, memory_size):
        """
        memory_free(label, memory_size)
        
        
        Defined at MPI_time_evaluate.f90 lines 157-179
        
        Parameters
        ----------
        label : str
        memory_size : float
        
        """
        _arespy.f90wrap_memory_free(label=label, memory_size=memory_size)
    
    @property
    def total_memory_consum(self):
        """
        Element total_memory_consum ftype=real(dp) pytype=float
        
        
        Defined at MPI_time_evaluate.f90 line 17
        
        """
        return _arespy.f90wrap_m_time_evaluate__get__total_memory_consum()
    
    @total_memory_consum.setter
    def total_memory_consum(self, total_memory_consum):
        _arespy.f90wrap_m_time_evaluate__set__total_memory_consum(total_memory_consum)
    
    @property
    def filename(self):
        """
        Element filename ftype=character(len=20) pytype=str
        
        
        Defined at MPI_time_evaluate.f90 line 18
        
        """
        return _arespy.f90wrap_m_time_evaluate__get__filename()
    
    @filename.setter
    def filename(self, filename):
        _arespy.f90wrap_m_time_evaluate__set__filename(filename)
    
    @property
    def file_unit(self):
        """
        Element file_unit ftype=integer(i4b) pytype=int
        
        
        Defined at MPI_time_evaluate.f90 line 19
        
        """
        return _arespy.f90wrap_m_time_evaluate__get__file_unit()
    
    @file_unit.setter
    def file_unit(self, file_unit):
        _arespy.f90wrap_m_time_evaluate__set__file_unit(file_unit)
    
    def __str__(self):
        ret = ['<m_time_evaluate>{\n']
        ret.append('    total_memory_consum : ')
        ret.append(repr(self.total_memory_consum))
        ret.append(',\n    filename : ')
        ret.append(repr(self.filename))
        ret.append(',\n    file_unit : ')
        ret.append(repr(self.file_unit))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

m_time_evaluate = M_Time_Evaluate()

class Nlpot_Module(f90wrap.runtime.FortranModule):
    """
    Module nlpot_module
    
    
    Defined at Nonlocalpot_module.f90 lines 1-497
    
    """
    @f90wrap.runtime.register_class("arespy.nol_type")
    class nol_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=nol_type)
        
        
        Defined at Nonlocalpot_module.f90 lines 10-20
        
        """
        def __init__(self, handle=None):
            """
            self = Nol_Type()
            
            
            Defined at Nonlocalpot_module.f90 lines 10-20
            
            
            Returns
            -------
            this : Nol_Type
            	Object to be constructed
            
            
            Automatically generated constructor for nol_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_nol_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Nol_Type
            
            
            Defined at Nonlocalpot_module.f90 lines 10-20
            
            Parameters
            ----------
            this : Nol_Type
            	Object to be destructed
            
            
            Automatically generated destructor for nol_type
            """
            if self._alloc:
                _arespy.f90wrap_nol_type_finalise(this=self._handle)
        
        @property
        def npts(self):
            """
            Element npts ftype=integer(i4b) pytype=int
            
            
            Defined at Nonlocalpot_module.f90 line 11
            
            """
            return _arespy.f90wrap_nol_type__get__npts(self._handle)
        
        @npts.setter
        def npts(self, npts):
            _arespy.f90wrap_nol_type__set__npts(self._handle, npts)
        
        @property
        def id(self):
            """
            Element id ftype=integer(i4b) pytype=int
            
            
            Defined at Nonlocalpot_module.f90 line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__id(self._handle)
            if array_handle in self._arrays:
                id = self._arrays[array_handle]
            else:
                id = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__id)
                self._arrays[array_handle] = id
            return id
        
        @id.setter
        def id(self, id):
            self.id[...] = id
        
        @property
        def rrvec(self):
            """
            Element rrvec ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.f90 line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__rrvec(self._handle)
            if array_handle in self._arrays:
                rrvec = self._arrays[array_handle]
            else:
                rrvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__rrvec)
                self._arrays[array_handle] = rrvec
            return rrvec
        
        @rrvec.setter
        def rrvec(self, rrvec):
            self.rrvec[...] = rrvec
        
        @property
        def proj0(self):
            """
            Element proj0 ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.f90 line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__proj0(self._handle)
            if array_handle in self._arrays:
                proj0 = self._arrays[array_handle]
            else:
                proj0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__proj0)
                self._arrays[array_handle] = proj0
            return proj0
        
        @proj0.setter
        def proj0(self, proj0):
            self.proj0[...] = proj0
        
        @property
        def proj(self):
            """
            Element proj ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.f90 line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__proj(self._handle)
            if array_handle in self._arrays:
                proj = self._arrays[array_handle]
            else:
                proj = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__proj)
                self._arrays[array_handle] = proj
            return proj
        
        @proj.setter
        def proj(self, proj):
            self.proj[...] = proj
        
        @property
        def proj_phs(self):
            """
            Element proj_phs ftype=complex(dcp) pytype=complex
            
            
            Defined at Nonlocalpot_module.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__proj_phs(self._handle)
            if array_handle in self._arrays:
                proj_phs = self._arrays[array_handle]
            else:
                proj_phs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__proj_phs)
                self._arrays[array_handle] = proj_phs
            return proj_phs
        
        @proj_phs.setter
        def proj_phs(self, proj_phs):
            self.proj_phs[...] = proj_phs
        
        @property
        def proj0_dg(self):
            """
            Element proj0_dg ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.f90 line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__proj0_dg(self._handle)
            if array_handle in self._arrays:
                proj0_dg = self._arrays[array_handle]
            else:
                proj0_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__proj0_dg)
                self._arrays[array_handle] = proj0_dg
            return proj0_dg
        
        @proj0_dg.setter
        def proj0_dg(self, proj0_dg):
            self.proj0_dg[...] = proj0_dg
        
        @property
        def proj_dg(self):
            """
            Element proj_dg ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__proj_dg(self._handle)
            if array_handle in self._arrays:
                proj_dg = self._arrays[array_handle]
            else:
                proj_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__proj_dg)
                self._arrays[array_handle] = proj_dg
            return proj_dg
        
        @proj_dg.setter
        def proj_dg(self, proj_dg):
            self.proj_dg[...] = proj_dg
        
        @property
        def proj_phs_dg(self):
            """
            Element proj_phs_dg ftype=complex(dcp) pytype=complex
            
            
            Defined at Nonlocalpot_module.f90 line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_nol_type__array__proj_phs_dg(self._handle)
            if array_handle in self._arrays:
                proj_phs_dg = self._arrays[array_handle]
            else:
                proj_phs_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_nol_type__array__proj_phs_dg)
                self._arrays[array_handle] = proj_phs_dg
            return proj_phs_dg
        
        @proj_phs_dg.setter
        def proj_phs_dg(self, proj_phs_dg):
            self.proj_phs_dg[...] = proj_phs_dg
        
        def __str__(self):
            ret = ['<nol_type>{\n']
            ret.append('    npts : ')
            ret.append(repr(self.npts))
            ret.append(',\n    id : ')
            ret.append(repr(self.id))
            ret.append(',\n    rrvec : ')
            ret.append(repr(self.rrvec))
            ret.append(',\n    proj0 : ')
            ret.append(repr(self.proj0))
            ret.append(',\n    proj : ')
            ret.append(repr(self.proj))
            ret.append(',\n    proj_phs : ')
            ret.append(repr(self.proj_phs))
            ret.append(',\n    proj0_dg : ')
            ret.append(repr(self.proj0_dg))
            ret.append(',\n    proj_dg : ')
            ret.append(repr(self.proj_dg))
            ret.append(',\n    proj_phs_dg : ')
            ret.append(repr(self.proj_phs_dg))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def initialize_nlpot():
        """
        initialize_nlpot()
        
        
        Defined at Nonlocalpot_module.f90 lines 28-192
        
        
        ===================store the data================
        """
        _arespy.f90wrap_initialize_nlpot()
    
    @staticmethod
    def destroy_nlpot():
        """
        destroy_nlpot()
        
        
        Defined at Nonlocalpot_module.f90 lines 195-214
        
        
        """
        _arespy.f90wrap_destroy_nlpot()
    
    @staticmethod
    def set_beta_real():
        """
        set_beta_real()
        
        
        Defined at Nonlocalpot_module.f90 lines 217-248
        
        
        ============apply Ylm in real space===============
        cycle the all species
        """
        _arespy.f90wrap_set_beta_real()
    
    @staticmethod
    def nlp_beta_interp_r(ity, ia, beta_init):
        """
        nlp_beta_interp_r(ity, ia, beta_init)
        
        
        Defined at Nonlocalpot_module.f90 lines 251-282
        
        Parameters
        ----------
        ity : int
        ia : int
        beta_init : float array
        
        """
        _arespy.f90wrap_nlp_beta_interp_r(ity=ity, ia=ia, beta_init=beta_init)
    
    @staticmethod
    def nlp_beta_ylm_r(ity, ia, beta_ylm):
        """
        nlp_beta_ylm_r(ity, ia, beta_ylm)
        
        
        Defined at Nonlocalpot_module.f90 lines 285-322
        
        Parameters
        ----------
        ity : int
        ia : int
        beta_ylm : float array
        
        """
        _arespy.f90wrap_nlp_beta_ylm_r(ity=ity, ia=ia, beta_ylm=beta_ylm)
    
    @staticmethod
    def nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase):
        """
        nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase)
        
        
        Defined at Nonlocalpot_module.f90 lines 325-361
        
        Parameters
        ----------
        ity : int
        ia : int
        beta_ylm : float array
        beta_phase : complex array
        
        """
        _arespy.f90wrap_nlp_beta_phase_r(ity=ity, ia=ia, beta_ylm=beta_ylm, \
            beta_phase=beta_phase)
    
    @staticmethod
    def apply_ylm(l, m, fac, x, y, z, f):
        """
        apply_ylm(l, m, fac, x, y, z, f)
        
        
        Defined at Nonlocalpot_module.f90 lines 364-440
        
        Parameters
        ----------
        l : int
        m : int
        fac : float
        x : float
        y : float
        z : float
        f : float
        
        """
        _arespy.f90wrap_apply_ylm(l=l, m=m, fac=fac, x=x, y=y, z=z, f=f)
    
    @staticmethod
    def nlp_init_partialcore(rhoc):
        """
        nlp_init_partialcore(rhoc)
        
        
        Defined at Nonlocalpot_module.f90 lines 443-496
        
        Parameters
        ----------
        rhoc : float array
        
        """
        _arespy.f90wrap_nlp_init_partialcore(rhoc=rhoc)
    
    @property
    def max_nlnpts(self):
        """
        Element max_nlnpts ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.f90 line 25
        
        """
        return _arespy.f90wrap_nlpot_module__get__max_nlnpts()
    
    @max_nlnpts.setter
    def max_nlnpts(self, max_nlnpts):
        _arespy.f90wrap_nlpot_module__set__max_nlnpts(max_nlnpts)
    
    def __str__(self):
        ret = ['<nlpot_module>{\n']
        ret.append('    max_nlnpts : ')
        ret.append(repr(self.max_nlnpts))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

nlpot_module = Nlpot_Module()

class Parameters(f90wrap.runtime.FortranModule):
    """
    Module parameters
    
    
    Defined at Parameters.f90 lines 1-79
    
    """
    @property
    def ixc(self):
        """
        Element ixc ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 7
        
        """
        return _arespy.f90wrap_parameters__get__ixc()
    
    @ixc.setter
    def ixc(self, ixc):
        _arespy.f90wrap_parameters__set__ixc(ixc)
    
    @property
    def nspin(self):
        """
        Element nspin ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 7
        
        """
        return _arespy.f90wrap_parameters__get__nspin()
    
    @nspin.setter
    def nspin(self, nspin):
        _arespy.f90wrap_parameters__set__nspin(nspin)
    
    @property
    def finite_order(self):
        """
        Element finite_order ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 8
        
        """
        return _arespy.f90wrap_parameters__get__finite_order()
    
    @finite_order.setter
    def finite_order(self, finite_order):
        _arespy.f90wrap_parameters__set__finite_order(finite_order)
    
    @property
    def ntype(self):
        """
        Element ntype ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 9
        
        """
        return _arespy.f90wrap_parameters__get__ntype()
    
    @ntype.setter
    def ntype(self, ntype):
        _arespy.f90wrap_parameters__set__ntype(ntype)
    
    @property
    def naddstates(self):
        """
        Element naddstates ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__naddstates()
    
    @naddstates.setter
    def naddstates(self, naddstates):
        _arespy.f90wrap_parameters__set__naddstates(naddstates)
    
    @property
    def igamma(self):
        """
        Element igamma ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__igamma()
    
    @igamma.setter
    def igamma(self, igamma):
        _arespy.f90wrap_parameters__set__igamma(igamma)
    
    @property
    def istart(self):
        """
        Element istart ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__istart()
    
    @istart.setter
    def istart(self, istart):
        _arespy.f90wrap_parameters__set__istart(istart)
    
    @property
    def isym(self):
        """
        Element isym ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__isym()
    
    @isym.setter
    def isym(self, isym):
        _arespy.f90wrap_parameters__set__isym(isym)
    
    @property
    def idiag(self):
        """
        Element idiag ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__idiag()
    
    @idiag.setter
    def idiag(self, idiag):
        _arespy.f90wrap_parameters__set__idiag(idiag)
    
    @property
    def chem(self):
        """
        Element chem ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__chem()
    
    @chem.setter
    def chem(self, chem):
        _arespy.f90wrap_parameters__set__chem(chem)
    
    @property
    def chem0(self):
        """
        Element chem0 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__chem0()
    
    @chem0.setter
    def chem0(self, chem0):
        _arespy.f90wrap_parameters__set__chem0(chem0)
    
    @property
    def nstates(self):
        """
        Element nstates ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__nstates()
    
    @nstates.setter
    def nstates(self, nstates):
        _arespy.f90wrap_parameters__set__nstates(nstates)
    
    @property
    def gridn(self):
        """
        Element gridn ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_parameters__array__gridn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gridn = self._arrays[array_handle]
        else:
            gridn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_parameters__array__gridn)
            self._arrays[array_handle] = gridn
        return gridn
    
    @gridn.setter
    def gridn(self, gridn):
        self.gridn[...] = gridn
    
    @property
    def kgrid(self):
        """
        Element kgrid ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_parameters__array__kgrid(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kgrid = self._arrays[array_handle]
        else:
            kgrid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_parameters__array__kgrid)
            self._arrays[array_handle] = kgrid
        return kgrid
    
    @kgrid.setter
    def kgrid(self, kgrid):
        self.kgrid[...] = kgrid
    
    @property
    def nssp(self):
        """
        Element nssp ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 20
        
        """
        return _arespy.f90wrap_parameters__get__nssp()
    
    @nssp.setter
    def nssp(self, nssp):
        _arespy.f90wrap_parameters__set__nssp(nssp)
    
    @property
    def dcharge(self):
        """
        Element dcharge ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 22
        
        """
        return _arespy.f90wrap_parameters__get__dcharge()
    
    @dcharge.setter
    def dcharge(self, dcharge):
        _arespy.f90wrap_parameters__set__dcharge(dcharge)
    
    @property
    def init_gap(self):
        """
        Element init_gap ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 27
        
        """
        return _arespy.f90wrap_parameters__get__init_gap()
    
    @init_gap.setter
    def init_gap(self, init_gap):
        _arespy.f90wrap_parameters__set__init_gap(init_gap)
    
    @property
    def ecut(self):
        """
        Element ecut ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 27
        
        """
        return _arespy.f90wrap_parameters__get__ecut()
    
    @ecut.setter
    def ecut(self, ecut):
        _arespy.f90wrap_parameters__set__ecut(ecut)
    
    @property
    def kspacing(self):
        """
        Element kspacing ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 27
        
        """
        return _arespy.f90wrap_parameters__get__kspacing()
    
    @kspacing.setter
    def kspacing(self, kspacing):
        _arespy.f90wrap_parameters__set__kspacing(kspacing)
    
    @property
    def atomrc(self):
        """
        Element atomrc ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 27
        
        """
        return _arespy.f90wrap_parameters__get__atomrc()
    
    @atomrc.setter
    def atomrc(self, atomrc):
        _arespy.f90wrap_parameters__set__atomrc(atomrc)
    
    @property
    def snlcc(self):
        """
        Element snlcc ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 27
        
        """
        return _arespy.f90wrap_parameters__get__snlcc()
    
    @snlcc.setter
    def snlcc(self, snlcc):
        _arespy.f90wrap_parameters__set__snlcc(snlcc)
    
    @property
    def lfirst(self):
        """
        Element lfirst ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 35
        
        """
        return _arespy.f90wrap_parameters__get__lfirst()
    
    @lfirst.setter
    def lfirst(self, lfirst):
        _arespy.f90wrap_parameters__set__lfirst(lfirst)
    
    @property
    def linrho(self):
        """
        Element linrho ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 35
        
        """
        return _arespy.f90wrap_parameters__get__linrho()
    
    @linrho.setter
    def linrho(self, linrho):
        _arespy.f90wrap_parameters__set__linrho(linrho)
    
    @property
    def latomrho(self):
        """
        Element latomrho ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 35
        
        """
        return _arespy.f90wrap_parameters__get__latomrho()
    
    @latomrho.setter
    def latomrho(self, latomrho):
        _arespy.f90wrap_parameters__set__latomrho(latomrho)
    
    @property
    def lrrorthnorm(self):
        """
        Element lrrorthnorm ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 35
        
        """
        return _arespy.f90wrap_parameters__get__lrrorthnorm()
    
    @lrrorthnorm.setter
    def lrrorthnorm(self, lrrorthnorm):
        _arespy.f90wrap_parameters__set__lrrorthnorm(lrrorthnorm)
    
    @property
    def lrandom(self):
        """
        Element lrandom ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 35
        
        """
        return _arespy.f90wrap_parameters__get__lrandom()
    
    @lrandom.setter
    def lrandom(self, lrandom):
        _arespy.f90wrap_parameters__set__lrandom(lrandom)
    
    @property
    def lcore_val(self):
        """
        Element lcore_val ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 35
        
        """
        return _arespy.f90wrap_parameters__get__lcore_val()
    
    @lcore_val.setter
    def lcore_val(self, lcore_val):
        _arespy.f90wrap_parameters__set__lcore_val(lcore_val)
    
    @property
    def system_name(self):
        """
        Element system_name ftype=character(30) pytype=str
        
        
        Defined at Parameters.f90 line 37
        
        """
        return _arespy.f90wrap_parameters__get__system_name()
    
    @system_name.setter
    def system_name(self, system_name):
        _arespy.f90wrap_parameters__set__system_name(system_name)
    
    @property
    def cellfile_name(self):
        """
        Element cellfile_name ftype=character(30) pytype=str
        
        
        Defined at Parameters.f90 line 38
        
        """
        return _arespy.f90wrap_parameters__get__cellfile_name()
    
    @cellfile_name.setter
    def cellfile_name(self, cellfile_name):
        _arespy.f90wrap_parameters__set__cellfile_name(cellfile_name)
    
    @property
    def ppfile_name(self):
        """
        Element ppfile_name ftype=character(30) pytype=str
        
        
        Defined at Parameters.f90 line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_parameters__array__ppfile_name(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ppfile_name = self._arrays[array_handle]
        else:
            ppfile_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_parameters__array__ppfile_name)
            self._arrays[array_handle] = ppfile_name
        return ppfile_name
    
    @ppfile_name.setter
    def ppfile_name(self, ppfile_name):
        self.ppfile_name[...] = ppfile_name
    
    @property
    def elements(self):
        """
        Element elements ftype=character(30) pytype=str
        
        
        Defined at Parameters.f90 line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_parameters__array__elements(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            elements = self._arrays[array_handle]
        else:
            elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_parameters__array__elements)
            self._arrays[array_handle] = elements
        return elements
    
    @elements.setter
    def elements(self, elements):
        self.elements[...] = elements
    
    @property
    def nsmear(self):
        """
        Element nsmear ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 42
        
        """
        return _arespy.f90wrap_parameters__get__nsmear()
    
    @nsmear.setter
    def nsmear(self, nsmear):
        _arespy.f90wrap_parameters__set__nsmear(nsmear)
    
    @property
    def wsmear(self):
        """
        Element wsmear ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 43
        
        """
        return _arespy.f90wrap_parameters__get__wsmear()
    
    @wsmear.setter
    def wsmear(self, wsmear):
        _arespy.f90wrap_parameters__set__wsmear(wsmear)
    
    @property
    def imixer(self):
        """
        Element imixer ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 46
        
        """
        return _arespy.f90wrap_parameters__get__imixer()
    
    @imixer.setter
    def imixer(self, imixer):
        _arespy.f90wrap_parameters__set__imixer(imixer)
    
    @property
    def nmiter(self):
        """
        Element nmiter ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 46
        
        """
        return _arespy.f90wrap_parameters__get__nmiter()
    
    @nmiter.setter
    def nmiter(self, nmiter):
        _arespy.f90wrap_parameters__set__nmiter(nmiter)
    
    @property
    def nsmix(self):
        """
        Element nsmix ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 49
        
        """
        return _arespy.f90wrap_parameters__get__nsmix()
    
    @nsmix.setter
    def nsmix(self, nsmix):
        _arespy.f90wrap_parameters__set__nsmix(nsmix)
    
    @property
    def nhmax(self):
        """
        Element nhmax ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 49
        
        """
        return _arespy.f90wrap_parameters__get__nhmax()
    
    @nhmax.setter
    def nhmax(self, nhmax):
        _arespy.f90wrap_parameters__set__nhmax(nhmax)
    
    @property
    def nhmin(self):
        """
        Element nhmin ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 49
        
        """
        return _arespy.f90wrap_parameters__get__nhmin()
    
    @nhmin.setter
    def nhmin(self, nhmin):
        _arespy.f90wrap_parameters__set__nhmin(nhmin)
    
    @property
    def malpha(self):
        """
        Element malpha ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 53
        
        """
        return _arespy.f90wrap_parameters__get__malpha()
    
    @malpha.setter
    def malpha(self, malpha):
        _arespy.f90wrap_parameters__set__malpha(malpha)
    
    @property
    def mbeta(self):
        """
        Element mbeta ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 53
        
        """
        return _arespy.f90wrap_parameters__get__mbeta()
    
    @mbeta.setter
    def mbeta(self, mbeta):
        _arespy.f90wrap_parameters__set__mbeta(mbeta)
    
    @property
    def amix(self):
        """
        Element amix ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 53
        
        """
        return _arespy.f90wrap_parameters__get__amix()
    
    @amix.setter
    def amix(self, amix):
        _arespy.f90wrap_parameters__set__amix(amix)
    
    @property
    def bmix(self):
        """
        Element bmix ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 53
        
        """
        return _arespy.f90wrap_parameters__get__bmix()
    
    @bmix.setter
    def bmix(self, bmix):
        _arespy.f90wrap_parameters__set__bmix(bmix)
    
    @property
    def resta(self):
        """
        Element resta ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_parameters__array__resta(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            resta = self._arrays[array_handle]
        else:
            resta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_parameters__array__resta)
            self._arrays[array_handle] = resta
        return resta
    
    @resta.setter
    def resta(self, resta):
        self.resta[...] = resta
    
    @property
    def w0am(self):
        """
        Element w0am ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 56
        
        """
        return _arespy.f90wrap_parameters__get__w0am()
    
    @w0am.setter
    def w0am(self, w0am):
        _arespy.f90wrap_parameters__set__w0am(w0am)
    
    @property
    def rtol(self):
        """
        Element rtol ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 58
        
        """
        return _arespy.f90wrap_parameters__get__rtol()
    
    @rtol.setter
    def rtol(self, rtol):
        _arespy.f90wrap_parameters__set__rtol(rtol)
    
    @property
    def etol(self):
        """
        Element etol ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 58
        
        """
        return _arespy.f90wrap_parameters__get__etol()
    
    @etol.setter
    def etol(self, etol):
        _arespy.f90wrap_parameters__set__etol(etol)
    
    @property
    def lband(self):
        """
        Element lband ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 60
        
        """
        return _arespy.f90wrap_parameters__get__lband()
    
    @lband.setter
    def lband(self, lband):
        _arespy.f90wrap_parameters__set__lband(lband)
    
    @property
    def iopm(self):
        """
        Element iopm ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 64
        
        """
        return _arespy.f90wrap_parameters__get__iopm()
    
    @iopm.setter
    def iopm(self, iopm):
        _arespy.f90wrap_parameters__set__iopm(iopm)
    
    @property
    def igoal(self):
        """
        Element igoal ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 64
        
        """
        return _arespy.f90wrap_parameters__get__igoal()
    
    @igoal.setter
    def igoal(self, igoal):
        _arespy.f90wrap_parameters__set__igoal(igoal)
    
    @property
    def press(self):
        """
        Element press ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 68
        
        """
        return _arespy.f90wrap_parameters__get__press()
    
    @press.setter
    def press(self, press):
        _arespy.f90wrap_parameters__set__press(press)
    
    @property
    def tolf(self):
        """
        Element tolf ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 68
        
        """
        return _arespy.f90wrap_parameters__get__tolf()
    
    @tolf.setter
    def tolf(self, tolf):
        _arespy.f90wrap_parameters__set__tolf(tolf)
    
    @property
    def tolp(self):
        """
        Element tolp ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 68
        
        """
        return _arespy.f90wrap_parameters__get__tolp()
    
    @tolp.setter
    def tolp(self, tolp):
        _arespy.f90wrap_parameters__set__tolp(tolp)
    
    @property
    def times(self):
        """
        Element times ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 68
        
        """
        return _arespy.f90wrap_parameters__get__times()
    
    @times.setter
    def times(self, times):
        _arespy.f90wrap_parameters__set__times(times)
    
    @property
    def lwave(self):
        """
        Element lwave ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 71
        
        """
        return _arespy.f90wrap_parameters__get__lwave()
    
    @lwave.setter
    def lwave(self, lwave):
        _arespy.f90wrap_parameters__set__lwave(lwave)
    
    @property
    def lcharge(self):
        """
        Element lcharge ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 71
        
        """
        return _arespy.f90wrap_parameters__get__lcharge()
    
    @lcharge.setter
    def lcharge(self, lcharge):
        _arespy.f90wrap_parameters__set__lcharge(lcharge)
    
    @property
    def lmom(self):
        """
        Element lmom ftype=logical pytype=bool
        
        
        Defined at Parameters.f90 line 73
        
        """
        return _arespy.f90wrap_parameters__get__lmom()
    
    @lmom.setter
    def lmom(self, lmom):
        _arespy.f90wrap_parameters__set__lmom(lmom)
    
    @property
    def momsigma(self):
        """
        Element momsigma ftype=real(dp) pytype=float
        
        
        Defined at Parameters.f90 line 74
        
        """
        return _arespy.f90wrap_parameters__get__momsigma()
    
    @momsigma.setter
    def momsigma(self, momsigma):
        _arespy.f90wrap_parameters__set__momsigma(momsigma)
    
    @property
    def nwf0(self):
        """
        Element nwf0 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 75
        
        """
        return _arespy.f90wrap_parameters__get__nwf0()
    
    @nwf0.setter
    def nwf0(self, nwf0):
        _arespy.f90wrap_parameters__set__nwf0(nwf0)
    
    @property
    def block_mbnb(self):
        """
        Element block_mbnb ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 78
        
        """
        return _arespy.f90wrap_parameters__get__block_mbnb()
    
    @block_mbnb.setter
    def block_mbnb(self, block_mbnb):
        _arespy.f90wrap_parameters__set__block_mbnb(block_mbnb)
    
    @property
    def nstates_global(self):
        """
        Element nstates_global ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 78
        
        """
        return _arespy.f90wrap_parameters__get__nstates_global()
    
    @nstates_global.setter
    def nstates_global(self, nstates_global):
        _arespy.f90wrap_parameters__set__nstates_global(nstates_global)
    
    @property
    def npara(self):
        """
        Element npara ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.f90 line 79
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_parameters__array__npara(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            npara = self._arrays[array_handle]
        else:
            npara = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_parameters__array__npara)
            self._arrays[array_handle] = npara
        return npara
    
    @npara.setter
    def npara(self, npara):
        self.npara[...] = npara
    
    def __str__(self):
        ret = ['<parameters>{\n']
        ret.append('    ixc : ')
        ret.append(repr(self.ixc))
        ret.append(',\n    nspin : ')
        ret.append(repr(self.nspin))
        ret.append(',\n    finite_order : ')
        ret.append(repr(self.finite_order))
        ret.append(',\n    ntype : ')
        ret.append(repr(self.ntype))
        ret.append(',\n    naddstates : ')
        ret.append(repr(self.naddstates))
        ret.append(',\n    igamma : ')
        ret.append(repr(self.igamma))
        ret.append(',\n    istart : ')
        ret.append(repr(self.istart))
        ret.append(',\n    isym : ')
        ret.append(repr(self.isym))
        ret.append(',\n    idiag : ')
        ret.append(repr(self.idiag))
        ret.append(',\n    chem : ')
        ret.append(repr(self.chem))
        ret.append(',\n    chem0 : ')
        ret.append(repr(self.chem0))
        ret.append(',\n    nstates : ')
        ret.append(repr(self.nstates))
        ret.append(',\n    gridn : ')
        ret.append(repr(self.gridn))
        ret.append(',\n    kgrid : ')
        ret.append(repr(self.kgrid))
        ret.append(',\n    nssp : ')
        ret.append(repr(self.nssp))
        ret.append(',\n    dcharge : ')
        ret.append(repr(self.dcharge))
        ret.append(',\n    init_gap : ')
        ret.append(repr(self.init_gap))
        ret.append(',\n    ecut : ')
        ret.append(repr(self.ecut))
        ret.append(',\n    kspacing : ')
        ret.append(repr(self.kspacing))
        ret.append(',\n    atomrc : ')
        ret.append(repr(self.atomrc))
        ret.append(',\n    snlcc : ')
        ret.append(repr(self.snlcc))
        ret.append(',\n    lfirst : ')
        ret.append(repr(self.lfirst))
        ret.append(',\n    linrho : ')
        ret.append(repr(self.linrho))
        ret.append(',\n    latomrho : ')
        ret.append(repr(self.latomrho))
        ret.append(',\n    lrrorthnorm : ')
        ret.append(repr(self.lrrorthnorm))
        ret.append(',\n    lrandom : ')
        ret.append(repr(self.lrandom))
        ret.append(',\n    lcore_val : ')
        ret.append(repr(self.lcore_val))
        ret.append(',\n    system_name : ')
        ret.append(repr(self.system_name))
        ret.append(',\n    cellfile_name : ')
        ret.append(repr(self.cellfile_name))
        ret.append(',\n    ppfile_name : ')
        ret.append(repr(self.ppfile_name))
        ret.append(',\n    elements : ')
        ret.append(repr(self.elements))
        ret.append(',\n    nsmear : ')
        ret.append(repr(self.nsmear))
        ret.append(',\n    wsmear : ')
        ret.append(repr(self.wsmear))
        ret.append(',\n    imixer : ')
        ret.append(repr(self.imixer))
        ret.append(',\n    nmiter : ')
        ret.append(repr(self.nmiter))
        ret.append(',\n    nsmix : ')
        ret.append(repr(self.nsmix))
        ret.append(',\n    nhmax : ')
        ret.append(repr(self.nhmax))
        ret.append(',\n    nhmin : ')
        ret.append(repr(self.nhmin))
        ret.append(',\n    malpha : ')
        ret.append(repr(self.malpha))
        ret.append(',\n    mbeta : ')
        ret.append(repr(self.mbeta))
        ret.append(',\n    amix : ')
        ret.append(repr(self.amix))
        ret.append(',\n    bmix : ')
        ret.append(repr(self.bmix))
        ret.append(',\n    resta : ')
        ret.append(repr(self.resta))
        ret.append(',\n    w0am : ')
        ret.append(repr(self.w0am))
        ret.append(',\n    rtol : ')
        ret.append(repr(self.rtol))
        ret.append(',\n    etol : ')
        ret.append(repr(self.etol))
        ret.append(',\n    lband : ')
        ret.append(repr(self.lband))
        ret.append(',\n    iopm : ')
        ret.append(repr(self.iopm))
        ret.append(',\n    igoal : ')
        ret.append(repr(self.igoal))
        ret.append(',\n    press : ')
        ret.append(repr(self.press))
        ret.append(',\n    tolf : ')
        ret.append(repr(self.tolf))
        ret.append(',\n    tolp : ')
        ret.append(repr(self.tolp))
        ret.append(',\n    times : ')
        ret.append(repr(self.times))
        ret.append(',\n    lwave : ')
        ret.append(repr(self.lwave))
        ret.append(',\n    lcharge : ')
        ret.append(repr(self.lcharge))
        ret.append(',\n    lmom : ')
        ret.append(repr(self.lmom))
        ret.append(',\n    momsigma : ')
        ret.append(repr(self.momsigma))
        ret.append(',\n    nwf0 : ')
        ret.append(repr(self.nwf0))
        ret.append(',\n    block_mbnb : ')
        ret.append(repr(self.block_mbnb))
        ret.append(',\n    nstates_global : ')
        ret.append(repr(self.nstates_global))
        ret.append(',\n    npara : ')
        ret.append(repr(self.npara))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

parameters = Parameters()

class Potential_Module(f90wrap.runtime.FortranModule):
    """
    Module potential_module
    
    
    Defined at Potential_module.f90 lines 1-266
    
    """
    @staticmethod
    def calveff(nps, rhos, rho, veffs):
        """
        calveff(nps, rhos, rho, veffs)
        
        
        Defined at Potential_module.f90 lines 19-61
        
        Parameters
        ----------
        nps : int
        rhos : float array
        rho : float array
        veffs : float array
        
        """
        _arespy.f90wrap_calveff(nps=nps, rhos=rhos, rho=rho, veffs=veffs)
    
    @staticmethod
    def vhartree(nps, rho, vhart):
        """
        vhartree(nps, rho, vhart)
        
        
        Defined at Potential_module.f90 lines 64-162
        
        Parameters
        ----------
        nps : int
        rho : float array
        vhart : float array
        
        """
        _arespy.f90wrap_vhartree(nps=nps, rho=rho, vhart=vhart)
    
    @staticmethod
    def vlpp():
        """
        vlpp()
        
        
        Defined at Potential_module.f90 lines 165-265
        
        
        """
        _arespy.f90wrap_vlpp()
    
    @property
    def v_accelerate(self):
        """
        Element v_accelerate ftype=real(dp) pytype=float
        
        
        Defined at Potential_module.f90 line 10
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_potential_module__array__v_accelerate(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            v_accelerate = self._arrays[array_handle]
        else:
            v_accelerate = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_potential_module__array__v_accelerate)
            self._arrays[array_handle] = v_accelerate
        return v_accelerate
    
    @v_accelerate.setter
    def v_accelerate(self, v_accelerate):
        self.v_accelerate[...] = v_accelerate
    
    @property
    def v_accelerate(self):
        """
        Element v_accelerate ftype=real(dp) pytype=float
        
        
        Defined at Potential_module.f90 line 12
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_potential_module__array__v_accelerate(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            v_accelerate = self._arrays[array_handle]
        else:
            v_accelerate = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_potential_module__array__v_accelerate)
            self._arrays[array_handle] = v_accelerate
        return v_accelerate
    
    @v_accelerate.setter
    def v_accelerate(self, v_accelerate):
        self.v_accelerate[...] = v_accelerate
    
    def __str__(self):
        ret = ['<potential_module>{\n']
        ret.append('    v_accelerate : ')
        ret.append(repr(self.v_accelerate))
        ret.append(',\n    v_accelerate : ')
        ret.append(repr(self.v_accelerate))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

potential_module = Potential_Module()

class Read_Module(f90wrap.runtime.FortranModule):
    """
    Module read_module
    
    
    Defined at Read_module.f90 lines 1-642
    
    """
    @staticmethod
    def read_file(infile):
        """
        read_file(infile)
        
        
        Defined at Read_module.f90 lines 15-453
        
        Parameters
        ----------
        infile : str
        
        -----------------------------------------------------------
        """
        _arespy.f90wrap_read_file(infile=infile)
    
    @staticmethod
    def read_poscar(nty, filename):
        """
        read_poscar(nty, filename)
        
        
        Defined at Read_module.f90 lines 456-579
        
        Parameters
        ----------
        nty : int
        filename : str
        
        """
        _arespy.f90wrap_read_poscar(nty=nty, filename=filename)
    
    @staticmethod
    def resetpos(natom, lat, pos, poscar):
        """
        resetpos(natom, lat, pos, poscar)
        
        
        Defined at Read_module.f90 lines 582-641
        
        Parameters
        ----------
        natom : int
        lat : float array
        pos : float array
        poscar : float array
        
        """
        _arespy.f90wrap_resetpos(natom=natom, lat=lat, pos=pos, poscar=poscar)
    
    _dt_array_initialisers = []
    

read_module = Read_Module()

class Scalapack_Module(f90wrap.runtime.FortranModule):
    """
    Module scalapack_module
    
    
    Defined at ScaLapack_module.f90 lines 1-1148
    
    """
    @staticmethod
    def init_scala():
        """
        init_scala()
        
        
        Defined at ScaLapack_module.f90 lines 25-45
        
        
        """
        _arespy.f90wrap_init_scala()
    
    @staticmethod
    def sl_orthnorm(amat, m, n, mb, nb):
        """
        sl_orthnorm(amat, m, n, mb, nb)
        
        
        Defined at ScaLapack_module.f90 lines 57-76
        
        Parameters
        ----------
        amat : complex array
        m : int
        n : int
        mb : int
        nb : int
        
        """
        _arespy.f90wrap_sl_orthnorm(amat=amat, m=m, n=n, mb=mb, nb=nb)
    
    @staticmethod
    def sl_orthnorm_real(amat, m, n, mb, nb):
        """
        sl_orthnorm_real(amat, m, n, mb, nb)
        
        
        Defined at ScaLapack_module.f90 lines 79-114
        
        Parameters
        ----------
        amat : float array
        m : int
        n : int
        mb : int
        nb : int
        
        """
        _arespy.f90wrap_sl_orthnorm_real(amat=amat, m=m, n=n, mb=mb, nb=nb)
    
    @staticmethod
    def sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        cmbin=None, cnbin=None):
        """
        sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
            cmbin, cnbin])
        
        
        Defined at ScaLapack_module.f90 lines 117-173
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : complex array
        bmat : complex array
        cmat : complex array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmbin : int
        cnbin : int
        
        ---------------------------------------------------------------------
        """
        _arespy.f90wrap_sl_matmat(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, \
            am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, cmbin=cmbin, \
            cnbin=cnbin)
    
    @staticmethod
    def sl_matmat_cmplx_cn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
        bmb, bnb, cmb, cnb):
        """
        sl_matmat_cmplx_cn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
            bnb, cmb, cnb)
        
        
        Defined at ScaLapack_module.f90 lines 176-285
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : complex array
        bmat : complex array
        cmat : complex array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmb : int
        cnb : int
        
        ---------------------------------------------------------------------
        """
        _arespy.f90wrap_sl_matmat_cmplx_cn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def sl_matmat_cmplx_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
        bmb, bnb, cmb, cnb):
        """
        sl_matmat_cmplx_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
            bnb, cmb, cnb)
        
        
        Defined at ScaLapack_module.f90 lines 288-347
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : complex array
        bmat : complex array
        cmat : complex array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmb : int
        cnb : int
        
        ---------------------------------------------------------------------
        > cmat_local
        """
        _arespy.f90wrap_sl_matmat_cmplx_nn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmbin=None, cnbin=None):
        """
        sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
            cmbin, cnbin])
        
        
        Defined at ScaLapack_module.f90 lines 350-435
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : complex array
        bmat : complex array
        cmat : complex array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmbin : int
        cnbin : int
        
        ---------------------------------------------------------------------
        blacs_contxt = parallel%comm
        nprow = parallel%numprocs
        npcol = 1
        CALL BLACS_GET( -1, 0, blacs_contxt )
        CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
        CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
        call blacs_pinfo(iam,nprocs)
        am = grid%tn
        bm = grid%tn
        if( .not. init_called) then
         call init_scala_sub()
        end if
        """
        _arespy.f90wrap_sl_matmat_sub(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, \
            am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, cmbin=cmbin, \
            cnbin=cnbin)
    
    @staticmethod
    def sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
        bmb, bnb, cmbin=None, cnbin=None):
        """
        sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
            bnb[, cmbin, cnbin])
        
        
        Defined at ScaLapack_module.f90 lines 438-523
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : float array
        bmat : float array
        cmat : float array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmbin : int
        cnbin : int
        
        ---------------------------------------------------------------------
        blacs_contxt = parallel%comm
        nprow = parallel%numprocs
        npcol = 1
        CALL BLACS_GET( -1, 0, blacs_contxt )
        CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
        CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
        call blacs_pinfo(iam,nprocs)
        am = grid%tn
        bm = grid%tn
        if( .not. init_called) then
         call init_scala_sub()
        end if
        """
        _arespy.f90wrap_sl_matmat_sub_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmbin=None, cnbin=None):
        """
        sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
            cmbin, cnbin])
        
        
        Defined at ScaLapack_module.f90 lines 526-623
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : float array
        bmat : float array
        cmat : float array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmbin : int
        cnbin : int
        
        ---------------------------------------------------------------------
        """
        _arespy.f90wrap_sl_matmat_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        evec, cm, cn, eval, cmbin=None, cnbin=None):
        """
        sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, \
            cm, cn, eval[, cmbin, cnbin])
        
        
        Defined at ScaLapack_module.f90 lines 626-709
        
        Parameters
        ----------
        dime : int
        amat : complex array
        bmat : complex array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        evec : complex array
        cm : int
        cn : int
        eval : float array
        cmbin : int
        cnbin : int
        
        """
        _arespy.f90wrap_sl_generalizeeigen(dime=dime, amat=amat, bmat=bmat, am=am, \
            an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, cm=cm, \
            cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, evec, cm, cn, eval, cmbin=None, cnbin=None):
        """
        sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            evec, cm, cn, eval[, cmbin, cnbin])
        
        
        Defined at ScaLapack_module.f90 lines 712-857
        
        Parameters
        ----------
        dime : int
        amat : float array
        bmat : float array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        evec : float array
        cm : int
        cn : int
        eval : float array
        cmbin : int
        cnbin : int
        
        """
        _arespy.f90wrap_sl_generalizeeigen_real(dime=dime, amat=amat, bmat=bmat, am=am, \
            an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, cm=cm, \
            cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_diagm_real(dime, amat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, \
        eval, cmb, cnb):
        """
        sl_diagm_real(dime, amat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, \
            eval, cmb, cnb)
        
        
        Defined at ScaLapack_module.f90 lines 860-943
        
        Parameters
        ----------
        dime : int
        amat : float array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        evec : float array
        cm : int
        cn : int
        eval : float array
        cmb : int
        cnb : int
        
        """
        _arespy.f90wrap_sl_diagm_real(dime=dime, amat=amat, am=am, an=an, bm=bm, bn=bn, \
            amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, cm=cm, cn=cn, eval=eval, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def twod_map_set(nstates, nrow, ncol, twod_map):
        """
        twod_map_set(nstates, nrow, ncol, twod_map)
        
        
        Defined at ScaLapack_module.f90 lines 947-965
        
        Parameters
        ----------
        nstates : int
        nrow : int
        ncol : int
        twod_map : int array
        
        """
        _arespy.f90wrap_twod_map_set(nstates=nstates, nrow=nrow, ncol=ncol, \
            twod_map=twod_map)
    
    @staticmethod
    def sl_matmat_real_tn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb):
        """
        sl_matmat_real_tn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
            bnb, cmb, cnb)
        
        
        Defined at ScaLapack_module.f90 lines 967-1085
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : float array
        bmat : float array
        cmat : float array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmb : int
        cnb : int
        
        ---------------------------------------------------------------------
        """
        _arespy.f90wrap_sl_matmat_real_tn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def sl_matmat_real_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb):
        """
        sl_matmat_real_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
            bnb, cmb, cnb)
        
        
        Defined at ScaLapack_module.f90 lines 1088-1147
        
        Parameters
        ----------
        opa : str
        opb : str
        amat : float array
        bmat : float array
        cmat : float array
        am : int
        an : int
        bm : int
        bn : int
        amb : int
        anb : int
        bmb : int
        bnb : int
        cmb : int
        cnb : int
        
        ---------------------------------------------------------------------
        > cmat_local
        """
        _arespy.f90wrap_sl_matmat_real_nn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @property
    def blacs_contxt(self):
        """
        Element blacs_contxt ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 11
        
        """
        return _arespy.f90wrap_scalapack_module__get__blacs_contxt()
    
    @blacs_contxt.setter
    def blacs_contxt(self, blacs_contxt):
        _arespy.f90wrap_scalapack_module__set__blacs_contxt(blacs_contxt)
    
    @property
    def dlen(self):
        """
        Element dlen ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 12
        
        """
        return _arespy.f90wrap_scalapack_module__get__dlen()
    
    @property
    def myrow(self):
        """
        Element myrow ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 13
        
        """
        return _arespy.f90wrap_scalapack_module__get__myrow()
    
    @myrow.setter
    def myrow(self, myrow):
        _arespy.f90wrap_scalapack_module__set__myrow(myrow)
    
    @property
    def mycol(self):
        """
        Element mycol ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 13
        
        """
        return _arespy.f90wrap_scalapack_module__get__mycol()
    
    @mycol.setter
    def mycol(self, mycol):
        _arespy.f90wrap_scalapack_module__set__mycol(mycol)
    
    @property
    def npcol(self):
        """
        Element npcol ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 14
        
        """
        return _arespy.f90wrap_scalapack_module__get__npcol()
    
    @npcol.setter
    def npcol(self, npcol):
        _arespy.f90wrap_scalapack_module__set__npcol(npcol)
    
    @property
    def nprow(self):
        """
        Element nprow ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 14
        
        """
        return _arespy.f90wrap_scalapack_module__get__nprow()
    
    @nprow.setter
    def nprow(self, nprow):
        _arespy.f90wrap_scalapack_module__set__nprow(nprow)
    
    @property
    def my_blacs_id(self):
        """
        Element my_blacs_id ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 15
        
        """
        return _arespy.f90wrap_scalapack_module__get__my_blacs_id()
    
    @my_blacs_id.setter
    def my_blacs_id(self, my_blacs_id):
        _arespy.f90wrap_scalapack_module__set__my_blacs_id(my_blacs_id)
    
    @property
    def np(self):
        """
        Element np ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 16
        
        """
        return _arespy.f90wrap_scalapack_module__get__np()
    
    @np.setter
    def np(self, np):
        _arespy.f90wrap_scalapack_module__set__np(np)
    
    @property
    def nq(self):
        """
        Element nq ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 16
        
        """
        return _arespy.f90wrap_scalapack_module__get__nq()
    
    @nq.setter
    def nq(self, nq):
        _arespy.f90wrap_scalapack_module__set__nq(nq)
    
    @property
    def iam(self):
        """
        Element iam ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 17
        
        """
        return _arespy.f90wrap_scalapack_module__get__iam()
    
    @iam.setter
    def iam(self, iam):
        _arespy.f90wrap_scalapack_module__set__iam(iam)
    
    @property
    def nprocs(self):
        """
        Element nprocs ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 17
        
        """
        return _arespy.f90wrap_scalapack_module__get__nprocs()
    
    @nprocs.setter
    def nprocs(self, nprocs):
        _arespy.f90wrap_scalapack_module__set__nprocs(nprocs)
    
    @property
    def init_called(self):
        """
        Element init_called ftype=logical pytype=bool
        
        
        Defined at ScaLapack_module.f90 line 18
        
        """
        return _arespy.f90wrap_scalapack_module__get__init_called()
    
    @init_called.setter
    def init_called(self, init_called):
        _arespy.f90wrap_scalapack_module__set__init_called(init_called)
    
    @property
    def info(self):
        """
        Element info ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 20
        
        """
        return _arespy.f90wrap_scalapack_module__get__info()
    
    @info.setter
    def info(self, info):
        _arespy.f90wrap_scalapack_module__set__info(info)
    
    @property
    def l_useless(self):
        """
        Element l_useless ftype=logical pytype=bool
        
        
        Defined at ScaLapack_module.f90 line 21
        
        """
        return _arespy.f90wrap_scalapack_module__get__l_useless()
    
    @l_useless.setter
    def l_useless(self, l_useless):
        _arespy.f90wrap_scalapack_module__set__l_useless(l_useless)
    
    @property
    def twod_map(self):
        """
        Element twod_map ftype=integer(i4b) pytype=int
        
        
        Defined at ScaLapack_module.f90 line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_scalapack_module__array__twod_map(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            twod_map = self._arrays[array_handle]
        else:
            twod_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_scalapack_module__array__twod_map)
            self._arrays[array_handle] = twod_map
        return twod_map
    
    @twod_map.setter
    def twod_map(self, twod_map):
        self.twod_map[...] = twod_map
    
    def __str__(self):
        ret = ['<scalapack_module>{\n']
        ret.append('    blacs_contxt : ')
        ret.append(repr(self.blacs_contxt))
        ret.append(',\n    dlen : ')
        ret.append(repr(self.dlen))
        ret.append(',\n    myrow : ')
        ret.append(repr(self.myrow))
        ret.append(',\n    mycol : ')
        ret.append(repr(self.mycol))
        ret.append(',\n    npcol : ')
        ret.append(repr(self.npcol))
        ret.append(',\n    nprow : ')
        ret.append(repr(self.nprow))
        ret.append(',\n    my_blacs_id : ')
        ret.append(repr(self.my_blacs_id))
        ret.append(',\n    np : ')
        ret.append(repr(self.np))
        ret.append(',\n    nq : ')
        ret.append(repr(self.nq))
        ret.append(',\n    iam : ')
        ret.append(repr(self.iam))
        ret.append(',\n    nprocs : ')
        ret.append(repr(self.nprocs))
        ret.append(',\n    init_called : ')
        ret.append(repr(self.init_called))
        ret.append(',\n    info : ')
        ret.append(repr(self.info))
        ret.append(',\n    l_useless : ')
        ret.append(repr(self.l_useless))
        ret.append(',\n    twod_map : ')
        ret.append(repr(self.twod_map))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

scalapack_module = Scalapack_Module()

class Scf_Module(f90wrap.runtime.FortranModule):
    """
    Module scf_module
    
    
    Defined at Scf_module.f90 lines 1-305
    
    """
    @staticmethod
    def electronicscf():
        """
        electronicscf()
        
        
        Defined at Scf_module.f90 lines 16-34
        
        
        """
        _arespy.f90wrap_electronicscf()
    
    @staticmethod
    def arpackscf(nps, rhos, rho, eig):
        """
        arpackscf(nps, rhos, rho, eig)
        
        
        Defined at Scf_module.f90 lines 41-114
        
        Parameters
        ----------
        nps : int
        rhos : float array
        rho : float array
        eig : Eigen_Type
        
        ===============================================================
        """
        _arespy.f90wrap_arpackscf(nps=nps, rhos=rhos, rho=rho, eig=eig._handle)
    
    @staticmethod
    def eigensolver_real(nps, nev, veff, psi, eval, diagtol):
        """
        eigensolver_real(nps, nev, veff, psi, eval, diagtol)
        
        
        Defined at Scf_module.f90 lines 117-160
        
        Parameters
        ----------
        nps : int
        nev : int
        veff : float array
        psi : float array
        eval : float array
        diagtol : float
        
        """
        _arespy.f90wrap_eigensolver_real(nps=nps, nev=nev, veff=veff, psi=psi, \
            eval=eval, diagtol=diagtol)
    
    @staticmethod
    def chefsi(nps, rhos, rho, eig):
        """
        chefsi(nps, rhos, rho, eig)
        
        
        Defined at Scf_module.f90 lines 168-286
        
        Parameters
        ----------
        nps : int
        rhos : float array
        rho : float array
        eig : Eigen_Type
        
        ===============================================================
        """
        _arespy.f90wrap_chefsi(nps=nps, rhos=rhos, rho=rho, eig=eig._handle)
    
    @staticmethod
    def filter_spin_gamma(nps, nev, veff, x, d):
        """
        filter_spin_gamma(nps, nev, veff, x, d)
        
        
        Defined at Scf_module.f90 lines 289-304
        
        Parameters
        ----------
        nps : int
        nev : int
        veff : float array
        x : float array
        d : float array
        
        """
        _arespy.f90wrap_filter_spin_gamma(nps=nps, nev=nev, veff=veff, x=x, d=d)
    
    @property
    def iwd(self):
        """
        Element iwd ftype=integer(i4b) pytype=int
        
        
        Defined at Scf_module.f90 line 12
        
        """
        return _arespy.f90wrap_scf_module__get__iwd()
    
    @iwd.setter
    def iwd(self, iwd):
        _arespy.f90wrap_scf_module__set__iwd(iwd)
    
    @property
    def lscf(self):
        """
        Element lscf ftype=logical pytype=bool
        
        
        Defined at Scf_module.f90 line 13
        
        """
        return _arespy.f90wrap_scf_module__get__lscf()
    
    @lscf.setter
    def lscf(self, lscf):
        _arespy.f90wrap_scf_module__set__lscf(lscf)
    
    def __str__(self):
        ret = ['<scf_module>{\n']
        ret.append('    iwd : ')
        ret.append(repr(self.iwd))
        ret.append(',\n    lscf : ')
        ret.append(repr(self.lscf))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

scf_module = Scf_Module()

class Smearing_Module(f90wrap.runtime.FortranModule):
    """
    Module smearing_module
    
    
    Defined at Smearing_module.f90 lines 1-401
    
    """
    @staticmethod
    def smear_init(nev):
        """
        smear_init(nev)
        
        
        Defined at Smearing_module.f90 lines 20-33
        
        Parameters
        ----------
        nev : int
        
        """
        _arespy.f90wrap_smear_init(nev=nev)
    
    @staticmethod
    def destroy_smear():
        """
        destroy_smear()
        
        
        Defined at Smearing_module.f90 lines 36-44
        
        
        """
        _arespy.f90wrap_destroy_smear()
    
    @staticmethod
    def hp(x, n):
        """
        hp = hp(x, n)
        
        
        Defined at Smearing_module.f90 lines 47-77
        
        Parameters
        ----------
        x : float
        n : int
        
        Returns
        -------
        hp : float
        
        """
        hp = _arespy.f90wrap_hp(x=x, n=n)
        return hp
    
    @staticmethod
    def smearsn(y, y0, w):
        """
        sn = smearsn(y, y0, w)
        
        
        Defined at Smearing_module.f90 lines 80-118
        
        Parameters
        ----------
        y : float
        y0 : float
        w : float
        
        Returns
        -------
        sn : float
        
        """
        sn = _arespy.f90wrap_smearsn(y=y, y0=y0, w=w)
        return sn
    
    @staticmethod
    def fermilevel(ne, nev, nk, wk, eval, sigma):
        """
        fermilevel(ne, nev, nk, wk, eval, sigma)
        
        
        Defined at Smearing_module.f90 lines 121-220
        
        Parameters
        ----------
        ne : float
        nev : int
        nk : int
        wk : float array
        eval : float array
        sigma : float
        
        """
        _arespy.f90wrap_fermilevel(ne=ne, nev=nev, nk=nk, wk=wk, eval=eval, sigma=sigma)
    
    @staticmethod
    def enpy(e_mu, fi):
        """
        enpy = enpy(e_mu, fi)
        
        
        Defined at Smearing_module.f90 lines 223-242
        
        Parameters
        ----------
        e_mu : float
        fi : float
        
        Returns
        -------
        enpy : float
        
        """
        enpy = _arespy.f90wrap_enpy(e_mu=e_mu, fi=fi)
        return enpy
    
    @staticmethod
    def whg(x, n):
        """
        whg = whg(x, n)
        
        
        Defined at Smearing_module.f90 lines 245-272
        
        Parameters
        ----------
        x : float
        n : int
        
        Returns
        -------
        whg : float
        
        """
        whg = _arespy.f90wrap_whg(x=x, n=n)
        return whg
    
    @staticmethod
    def updaterho_pbc(nps, nev, eig, wk, rhos, rho):
        """
        updaterho_pbc(nps, nev, eig, wk, rhos, rho)
        
        
        Defined at Smearing_module.f90 lines 275-348
        
        Parameters
        ----------
        nps : int
        nev : int
        eig : Eigen_Type
        wk : float array
        rhos : float array
        rho : float array
        
        """
        _arespy.f90wrap_updaterho_pbc(nps=nps, nev=nev, eig=eig._handle, wk=wk, \
            rhos=rhos, rho=rho)
    
    @staticmethod
    def smear_updaterho(nps, nev, ne, eig, rhos, rho):
        """
        smear_updaterho(nps, nev, ne, eig, rhos, rho)
        
        
        Defined at Smearing_module.f90 lines 351-400
        
        Parameters
        ----------
        nps : int
        nev : int
        ne : float
        eig : Eigen_Type
        rhos : float array
        rho : float array
        
        """
        _arespy.f90wrap_smear_updaterho(nps=nps, nev=nev, ne=ne, eig=eig._handle, \
            rhos=rhos, rho=rho)
    
    @property
    def wke(self):
        """
        Element wke ftype=real(dp) pytype=float
        
        
        Defined at Smearing_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_smearing_module__array__wke(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wke = self._arrays[array_handle]
        else:
            wke = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_smearing_module__array__wke)
            self._arrays[array_handle] = wke
        return wke
    
    @wke.setter
    def wke(self, wke):
        self.wke[...] = wke
    
    @property
    def fme(self):
        """
        Element fme ftype=real(dp) pytype=float
        
        
        Defined at Smearing_module.f90 line 17
        
        """
        return _arespy.f90wrap_smearing_module__get__fme()
    
    @fme.setter
    def fme(self, fme):
        _arespy.f90wrap_smearing_module__set__fme(fme)
    
    @property
    def ets(self):
        """
        Element ets ftype=real(dp) pytype=float
        
        
        Defined at Smearing_module.f90 line 17
        
        """
        return _arespy.f90wrap_smearing_module__get__ets()
    
    @ets.setter
    def ets(self, ets):
        _arespy.f90wrap_smearing_module__set__ets(ets)
    
    def __str__(self):
        ret = ['<smearing_module>{\n']
        ret.append('    wke : ')
        ret.append(repr(self.wke))
        ret.append(',\n    fme : ')
        ret.append(repr(self.fme))
        ret.append(',\n    ets : ')
        ret.append(repr(self.ets))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

smearing_module = Smearing_Module()

class Smpi_Math_Module(f90wrap.runtime.FortranModule):
    """
    Module smpi_math_module
    
    
    Defined at Smpi_math_module.f90 lines 1-2101
    
    """
    @f90wrap.runtime.register_class("arespy.parallel_type")
    class parallel_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=parallel_type)
        
        
        Defined at Smpi_math_module.f90 lines 6-26
        
        """
        def __init__(self, handle=None):
            """
            self = Parallel_Type()
            
            
            Defined at Smpi_math_module.f90 lines 6-26
            
            
            Returns
            -------
            this : Parallel_Type
            	Object to be constructed
            
            
            Automatically generated constructor for parallel_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_parallel_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Parallel_Type
            
            
            Defined at Smpi_math_module.f90 lines 6-26
            
            Parameters
            ----------
            this : Parallel_Type
            	Object to be destructed
            
            
            Automatically generated destructor for parallel_type
            """
            if self._alloc:
                _arespy.f90wrap_parallel_type_finalise(this=self._handle)
        
        @property
        def comm(self):
            """
            Element comm ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 7
            
            """
            return _arespy.f90wrap_parallel_type__get__comm(self._handle)
        
        @comm.setter
        def comm(self, comm):
            _arespy.f90wrap_parallel_type__set__comm(self._handle, comm)
        
        @property
        def myid(self):
            """
            Element myid ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 8
            
            """
            return _arespy.f90wrap_parallel_type__get__myid(self._handle)
        
        @myid.setter
        def myid(self, myid):
            _arespy.f90wrap_parallel_type__set__myid(self._handle, myid)
        
        @property
        def numprocs(self):
            """
            Element numprocs ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 9
            
            """
            return _arespy.f90wrap_parallel_type__get__numprocs(self._handle)
        
        @numprocs.setter
        def numprocs(self, numprocs):
            _arespy.f90wrap_parallel_type__set__numprocs(self._handle, numprocs)
        
        @property
        def rootid(self):
            """
            Element rootid ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 10
            
            """
            return _arespy.f90wrap_parallel_type__get__rootid(self._handle)
        
        @rootid.setter
        def rootid(self, rootid):
            _arespy.f90wrap_parallel_type__set__rootid(self._handle, rootid)
        
        @property
        def isroot(self):
            """
            Element isroot ftype=logical pytype=bool
            
            
            Defined at Smpi_math_module.f90 line 11
            
            """
            return _arespy.f90wrap_parallel_type__get__isroot(self._handle)
        
        @isroot.setter
        def isroot(self, isroot):
            _arespy.f90wrap_parallel_type__set__isroot(self._handle, isroot)
        
        @property
        def nstate_proc(self):
            """
            Element nstate_proc ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 13
            
            """
            return _arespy.f90wrap_parallel_type__get__nstate_proc(self._handle)
        
        @nstate_proc.setter
        def nstate_proc(self, nstate_proc):
            _arespy.f90wrap_parallel_type__set__nstate_proc(self._handle, nstate_proc)
        
        @property
        def sub2sum(self):
            """
            Element sub2sum ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__sub2sum(self._handle)
            if array_handle in self._arrays:
                sub2sum = self._arrays[array_handle]
            else:
                sub2sum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__sub2sum)
                self._arrays[array_handle] = sub2sum
            return sub2sum
        
        @sub2sum.setter
        def sub2sum(self, sub2sum):
            self.sub2sum[...] = sub2sum
        
        @property
        def mygrid_range(self):
            """
            Element mygrid_range ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__mygrid_range(self._handle)
            if array_handle in self._arrays:
                mygrid_range = self._arrays[array_handle]
            else:
                mygrid_range = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__mygrid_range)
                self._arrays[array_handle] = mygrid_range
            return mygrid_range
        
        @mygrid_range.setter
        def mygrid_range(self, mygrid_range):
            self.mygrid_range[...] = mygrid_range
        
        @property
        def recvcounts(self):
            """
            Element recvcounts ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__recvcounts(self._handle)
            if array_handle in self._arrays:
                recvcounts = self._arrays[array_handle]
            else:
                recvcounts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__recvcounts)
                self._arrays[array_handle] = recvcounts
            return recvcounts
        
        @recvcounts.setter
        def recvcounts(self, recvcounts):
            self.recvcounts[...] = recvcounts
        
        @property
        def displs(self):
            """
            Element displs ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__displs(self._handle)
            if array_handle in self._arrays:
                displs = self._arrays[array_handle]
            else:
                displs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__displs)
                self._arrays[array_handle] = displs
            return displs
        
        @displs.setter
        def displs(self, displs):
            self.displs[...] = displs
        
        @property
        def global_gridrange(self):
            """
            Element global_gridrange ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__global_gridrange(self._handle)
            if array_handle in self._arrays:
                global_gridrange = self._arrays[array_handle]
            else:
                global_gridrange = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__global_gridrange)
                self._arrays[array_handle] = global_gridrange
            return global_gridrange
        
        @global_gridrange.setter
        def global_gridrange(self, global_gridrange):
            self.global_gridrange[...] = global_gridrange
        
        @property
        def comm2d(self):
            """
            Element comm2d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__comm2d(self._handle)
        
        @comm2d.setter
        def comm2d(self, comm2d):
            _arespy.f90wrap_parallel_type__set__comm2d(self._handle, comm2d)
        
        @property
        def commx(self):
            """
            Element commx ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__commx(self._handle)
        
        @commx.setter
        def commx(self, commx):
            _arespy.f90wrap_parallel_type__set__commx(self._handle, commx)
        
        @property
        def commy(self):
            """
            Element commy ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__commy(self._handle)
        
        @commy.setter
        def commy(self, commy):
            _arespy.f90wrap_parallel_type__set__commy(self._handle, commy)
        
        @property
        def rankx(self):
            """
            Element rankx ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__rankx(self._handle)
        
        @rankx.setter
        def rankx(self, rankx):
            _arespy.f90wrap_parallel_type__set__rankx(self._handle, rankx)
        
        @property
        def ranky(self):
            """
            Element ranky ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__ranky(self._handle)
        
        @ranky.setter
        def ranky(self, ranky):
            _arespy.f90wrap_parallel_type__set__ranky(self._handle, ranky)
        
        @property
        def periods(self):
            """
            Element periods ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__periods(self._handle)
            if array_handle in self._arrays:
                periods = self._arrays[array_handle]
            else:
                periods = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__periods)
                self._arrays[array_handle] = periods
            return periods
        
        @periods.setter
        def periods(self, periods):
            self.periods[...] = periods
        
        @property
        def reorder(self):
            """
            Element reorder ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__reorder(self._handle)
        
        @reorder.setter
        def reorder(self, reorder):
            _arespy.f90wrap_parallel_type__set__reorder(self._handle, reorder)
        
        @property
        def remainx(self):
            """
            Element remainx ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__remainx(self._handle)
            if array_handle in self._arrays:
                remainx = self._arrays[array_handle]
            else:
                remainx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__remainx)
                self._arrays[array_handle] = remainx
            return remainx
        
        @remainx.setter
        def remainx(self, remainx):
            self.remainx[...] = remainx
        
        @property
        def remainy(self):
            """
            Element remainy ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__remainy(self._handle)
            if array_handle in self._arrays:
                remainy = self._arrays[array_handle]
            else:
                remainy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__remainy)
                self._arrays[array_handle] = remainy
            return remainy
        
        @remainy.setter
        def remainy(self, remainy):
            self.remainy[...] = remainy
        
        @property
        def ndims(self):
            """
            Element ndims ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            return _arespy.f90wrap_parallel_type__get__ndims(self._handle)
        
        @ndims.setter
        def ndims(self, ndims):
            _arespy.f90wrap_parallel_type__set__ndims(self._handle, ndims)
        
        @property
        def dims(self):
            """
            Element dims ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__dims(self._handle)
            if array_handle in self._arrays:
                dims = self._arrays[array_handle]
            else:
                dims = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__dims)
                self._arrays[array_handle] = dims
            return dims
        
        @dims.setter
        def dims(self, dims):
            self.dims[...] = dims
        
        @property
        def commfft(self):
            """
            Element commfft ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 23
            
            """
            return _arespy.f90wrap_parallel_type__get__commfft(self._handle)
        
        @commfft.setter
        def commfft(self, commfft):
            _arespy.f90wrap_parallel_type__set__commfft(self._handle, commfft)
        
        @property
        def local_z(self):
            """
            Element local_z ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 23
            
            """
            return _arespy.f90wrap_parallel_type__get__local_z(self._handle)
        
        @local_z.setter
        def local_z(self, local_z):
            _arespy.f90wrap_parallel_type__set__local_z(self._handle, local_z)
        
        @property
        def local_z_start(self):
            """
            Element local_z_start ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 23
            
            """
            return _arespy.f90wrap_parallel_type__get__local_z_start(self._handle)
        
        @local_z_start.setter
        def local_z_start(self, local_z_start):
            _arespy.f90wrap_parallel_type__set__local_z_start(self._handle, local_z_start)
        
        @property
        def fft_grid_range(self):
            """
            Element fft_grid_range ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__fft_grid_range(self._handle)
            if array_handle in self._arrays:
                fft_grid_range = self._arrays[array_handle]
            else:
                fft_grid_range = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__fft_grid_range)
                self._arrays[array_handle] = fft_grid_range
            return fft_grid_range
        
        @fft_grid_range.setter
        def fft_grid_range(self, fft_grid_range):
            self.fft_grid_range[...] = fft_grid_range
        
        @property
        def fft_rcount(self):
            """
            Element fft_rcount ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__fft_rcount(self._handle)
            if array_handle in self._arrays:
                fft_rcount = self._arrays[array_handle]
            else:
                fft_rcount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__fft_rcount)
                self._arrays[array_handle] = fft_rcount
            return fft_rcount
        
        @fft_rcount.setter
        def fft_rcount(self, fft_rcount):
            self.fft_rcount[...] = fft_rcount
        
        @property
        def fft_rdispls(self):
            """
            Element fft_rdispls ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__fft_rdispls(self._handle)
            if array_handle in self._arrays:
                fft_rdispls = self._arrays[array_handle]
            else:
                fft_rdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__fft_rdispls)
                self._arrays[array_handle] = fft_rdispls
            return fft_rdispls
        
        @fft_rdispls.setter
        def fft_rdispls(self, fft_rdispls):
            self.fft_rdispls[...] = fft_rdispls
        
        @property
        def fft_scount(self):
            """
            Element fft_scount ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__fft_scount(self._handle)
            if array_handle in self._arrays:
                fft_scount = self._arrays[array_handle]
            else:
                fft_scount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__fft_scount)
                self._arrays[array_handle] = fft_scount
            return fft_scount
        
        @fft_scount.setter
        def fft_scount(self, fft_scount):
            self.fft_scount[...] = fft_scount
        
        @property
        def fft_sdispls(self):
            """
            Element fft_sdispls ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_parallel_type__array__fft_sdispls(self._handle)
            if array_handle in self._arrays:
                fft_sdispls = self._arrays[array_handle]
            else:
                fft_sdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_parallel_type__array__fft_sdispls)
                self._arrays[array_handle] = fft_sdispls
            return fft_sdispls
        
        @fft_sdispls.setter
        def fft_sdispls(self, fft_sdispls):
            self.fft_sdispls[...] = fft_sdispls
        
        def __str__(self):
            ret = ['<parallel_type>{\n']
            ret.append('    comm : ')
            ret.append(repr(self.comm))
            ret.append(',\n    myid : ')
            ret.append(repr(self.myid))
            ret.append(',\n    numprocs : ')
            ret.append(repr(self.numprocs))
            ret.append(',\n    rootid : ')
            ret.append(repr(self.rootid))
            ret.append(',\n    isroot : ')
            ret.append(repr(self.isroot))
            ret.append(',\n    nstate_proc : ')
            ret.append(repr(self.nstate_proc))
            ret.append(',\n    sub2sum : ')
            ret.append(repr(self.sub2sum))
            ret.append(',\n    mygrid_range : ')
            ret.append(repr(self.mygrid_range))
            ret.append(',\n    recvcounts : ')
            ret.append(repr(self.recvcounts))
            ret.append(',\n    displs : ')
            ret.append(repr(self.displs))
            ret.append(',\n    global_gridrange : ')
            ret.append(repr(self.global_gridrange))
            ret.append(',\n    comm2d : ')
            ret.append(repr(self.comm2d))
            ret.append(',\n    commx : ')
            ret.append(repr(self.commx))
            ret.append(',\n    commy : ')
            ret.append(repr(self.commy))
            ret.append(',\n    rankx : ')
            ret.append(repr(self.rankx))
            ret.append(',\n    ranky : ')
            ret.append(repr(self.ranky))
            ret.append(',\n    periods : ')
            ret.append(repr(self.periods))
            ret.append(',\n    reorder : ')
            ret.append(repr(self.reorder))
            ret.append(',\n    remainx : ')
            ret.append(repr(self.remainx))
            ret.append(',\n    remainy : ')
            ret.append(repr(self.remainy))
            ret.append(',\n    ndims : ')
            ret.append(repr(self.ndims))
            ret.append(',\n    dims : ')
            ret.append(repr(self.dims))
            ret.append(',\n    commfft : ')
            ret.append(repr(self.commfft))
            ret.append(',\n    local_z : ')
            ret.append(repr(self.local_z))
            ret.append(',\n    local_z_start : ')
            ret.append(repr(self.local_z_start))
            ret.append(',\n    fft_grid_range : ')
            ret.append(repr(self.fft_grid_range))
            ret.append(',\n    fft_rcount : ')
            ret.append(repr(self.fft_rcount))
            ret.append(',\n    fft_rdispls : ')
            ret.append(repr(self.fft_rdispls))
            ret.append(',\n    fft_scount : ')
            ret.append(repr(self.fft_scount))
            ret.append(',\n    fft_sdispls : ')
            ret.append(repr(self.fft_sdispls))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.smpi_root_type")
    class smpi_root_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=smpi_root_type)
        
        
        Defined at Smpi_math_module.f90 lines 30-31
        
        """
        def __init__(self, handle=None):
            """
            self = Smpi_Root_Type()
            
            
            Defined at Smpi_math_module.f90 lines 30-31
            
            
            Returns
            -------
            this : Smpi_Root_Type
            	Object to be constructed
            
            
            Automatically generated constructor for smpi_root_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_smpi_root_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Smpi_Root_Type
            
            
            Defined at Smpi_math_module.f90 lines 30-31
            
            Parameters
            ----------
            this : Smpi_Root_Type
            	Object to be destructed
            
            
            Automatically generated destructor for smpi_root_type
            """
            if self._alloc:
                _arespy.f90wrap_smpi_root_type_finalise(this=self._handle)
        
        @property
        def natom_group(self):
            """
            Element natom_group ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 31
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_smpi_root_type__array__natom_group(self._handle)
            if array_handle in self._arrays:
                natom_group = self._arrays[array_handle]
            else:
                natom_group = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_smpi_root_type__array__natom_group)
                self._arrays[array_handle] = natom_group
            return natom_group
        
        @natom_group.setter
        def natom_group(self, natom_group):
            self.natom_group[...] = natom_group
        
        def __str__(self):
            ret = ['<smpi_root_type>{\n']
            ret.append('    natom_group : ')
            ret.append(repr(self.natom_group))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.smpi_comm_type")
    class smpi_comm_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=smpi_comm_type)
        
        
        Defined at Smpi_math_module.f90 lines 35-37
        
        """
        def __init__(self, handle=None):
            """
            self = Smpi_Comm_Type()
            
            
            Defined at Smpi_math_module.f90 lines 35-37
            
            
            Returns
            -------
            this : Smpi_Comm_Type
            	Object to be constructed
            
            
            Automatically generated constructor for smpi_comm_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_smpi_comm_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Smpi_Comm_Type
            
            
            Defined at Smpi_math_module.f90 lines 35-37
            
            Parameters
            ----------
            this : Smpi_Comm_Type
            	Object to be destructed
            
            
            Automatically generated destructor for smpi_comm_type
            """
            if self._alloc:
                _arespy.f90wrap_smpi_comm_type_finalise(this=self._handle)
        
        @property
        def atoms(self):
            """
            Element atoms ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 36
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_smpi_comm_type__array__atoms(self._handle)
            if array_handle in self._arrays:
                atoms = self._arrays[array_handle]
            else:
                atoms = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_smpi_comm_type__array__atoms)
                self._arrays[array_handle] = atoms
            return atoms
        
        @atoms.setter
        def atoms(self, atoms):
            self.atoms[...] = atoms
        
        @property
        def displs(self):
            """
            Element displs ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 37
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_smpi_comm_type__array__displs(self._handle)
            if array_handle in self._arrays:
                displs = self._arrays[array_handle]
            else:
                displs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_smpi_comm_type__array__displs)
                self._arrays[array_handle] = displs
            return displs
        
        @displs.setter
        def displs(self, displs):
            self.displs[...] = displs
        
        def __str__(self):
            ret = ['<smpi_comm_type>{\n']
            ret.append('    atoms : ')
            ret.append(repr(self.atoms))
            ret.append(',\n    displs : ')
            ret.append(repr(self.displs))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.time_type")
    class time_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=time_type)
        
        
        Defined at Smpi_math_module.f90 lines 41-47
        
        """
        def __init__(self, handle=None):
            """
            self = Time_Type()
            
            
            Defined at Smpi_math_module.f90 lines 41-47
            
            
            Returns
            -------
            this : Time_Type
            	Object to be constructed
            
            
            Automatically generated constructor for time_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_time_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Time_Type
            
            
            Defined at Smpi_math_module.f90 lines 41-47
            
            Parameters
            ----------
            this : Time_Type
            	Object to be destructed
            
            
            Automatically generated destructor for time_type
            """
            if self._alloc:
                _arespy.f90wrap_time_type_finalise(this=self._handle)
        
        @property
        def label(self):
            """
            Element label ftype=character(len=100) pytype=str
            
            
            Defined at Smpi_math_module.f90 line 42
            
            """
            return _arespy.f90wrap_time_type__get__label(self._handle)
        
        @label.setter
        def label(self, label):
            _arespy.f90wrap_time_type__set__label(self._handle, label)
        
        @property
        def tic(self):
            """
            Element tic ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.f90 line 43
            
            """
            return _arespy.f90wrap_time_type__get__tic(self._handle)
        
        @tic.setter
        def tic(self, tic):
            _arespy.f90wrap_time_type__set__tic(self._handle, tic)
        
        @property
        def toc(self):
            """
            Element toc ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.f90 line 44
            
            """
            return _arespy.f90wrap_time_type__get__toc(self._handle)
        
        @toc.setter
        def toc(self, toc):
            _arespy.f90wrap_time_type__set__toc(self._handle, toc)
        
        @property
        def total(self):
            """
            Element total ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.f90 line 45
            
            """
            return _arespy.f90wrap_time_type__get__total(self._handle)
        
        @total.setter
        def total(self, total):
            _arespy.f90wrap_time_type__set__total(self._handle, total)
        
        @property
        def sum_total(self):
            """
            Element sum_total ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.f90 line 46
            
            """
            return _arespy.f90wrap_time_type__get__sum_total(self._handle)
        
        @sum_total.setter
        def sum_total(self, sum_total):
            _arespy.f90wrap_time_type__set__sum_total(self._handle, sum_total)
        
        @property
        def num(self):
            """
            Element num ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 47
            
            """
            return _arespy.f90wrap_time_type__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _arespy.f90wrap_time_type__set__num(self._handle, num)
        
        def __str__(self):
            ret = ['<time_type>{\n']
            ret.append('    label : ')
            ret.append(repr(self.label))
            ret.append(',\n    tic : ')
            ret.append(repr(self.tic))
            ret.append(',\n    toc : ')
            ret.append(repr(self.toc))
            ret.append(',\n    total : ')
            ret.append(repr(self.total))
            ret.append(',\n    sum_total : ')
            ret.append(repr(self.sum_total))
            ret.append(',\n    num : ')
            ret.append(repr(self.num))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.mem_type")
    class mem_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=mem_type)
        
        
        Defined at Smpi_math_module.f90 lines 49-53
        
        """
        def __init__(self, handle=None):
            """
            self = Mem_Type()
            
            
            Defined at Smpi_math_module.f90 lines 49-53
            
            
            Returns
            -------
            this : Mem_Type
            	Object to be constructed
            
            
            Automatically generated constructor for mem_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_mem_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Mem_Type
            
            
            Defined at Smpi_math_module.f90 lines 49-53
            
            Parameters
            ----------
            this : Mem_Type
            	Object to be destructed
            
            
            Automatically generated destructor for mem_type
            """
            if self._alloc:
                _arespy.f90wrap_mem_type_finalise(this=self._handle)
        
        @property
        def label(self):
            """
            Element label ftype=character(len=100) pytype=str
            
            
            Defined at Smpi_math_module.f90 line 50
            
            """
            return _arespy.f90wrap_mem_type__get__label(self._handle)
        
        @label.setter
        def label(self, label):
            _arespy.f90wrap_mem_type__set__label(self._handle, label)
        
        @property
        def memic(self):
            """
            Element memic ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.f90 line 51
            
            """
            return _arespy.f90wrap_mem_type__get__memic(self._handle)
        
        @memic.setter
        def memic(self, memic):
            _arespy.f90wrap_mem_type__set__memic(self._handle, memic)
        
        @property
        def total(self):
            """
            Element total ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.f90 line 52
            
            """
            return _arespy.f90wrap_mem_type__get__total(self._handle)
        
        @total.setter
        def total(self, total):
            _arespy.f90wrap_mem_type__set__total(self._handle, total)
        
        @property
        def num(self):
            """
            Element num ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 53
            
            """
            return _arespy.f90wrap_mem_type__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _arespy.f90wrap_mem_type__set__num(self._handle, num)
        
        def __str__(self):
            ret = ['<mem_type>{\n']
            ret.append('    label : ')
            ret.append(repr(self.label))
            ret.append(',\n    memic : ')
            ret.append(repr(self.memic))
            ret.append(',\n    total : ')
            ret.append(repr(self.total))
            ret.append(',\n    num : ')
            ret.append(repr(self.num))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.grid_diff_map_type")
    class grid_diff_map_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=grid_diff_map_type)
        
        
        Defined at Smpi_math_module.f90 lines 60-70
        
        """
        def __init__(self, handle=None):
            """
            self = Grid_Diff_Map_Type()
            
            
            Defined at Smpi_math_module.f90 lines 60-70
            
            
            Returns
            -------
            this : Grid_Diff_Map_Type
            	Object to be constructed
            
            
            Automatically generated constructor for grid_diff_map_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_grid_diff_map_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Grid_Diff_Map_Type
            
            
            Defined at Smpi_math_module.f90 lines 60-70
            
            Parameters
            ----------
            this : Grid_Diff_Map_Type
            	Object to be destructed
            
            
            Automatically generated destructor for grid_diff_map_type
            """
            if self._alloc:
                _arespy.f90wrap_grid_diff_map_type_finalise(this=self._handle)
        
        @property
        def nz_map(self):
            """
            Element nz_map ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 61
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__nz_map(self._handle)
            if array_handle in self._arrays:
                nz_map = self._arrays[array_handle]
            else:
                nz_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__nz_map)
                self._arrays[array_handle] = nz_map
            return nz_map
        
        @nz_map.setter
        def nz_map(self, nz_map):
            self.nz_map[...] = nz_map
        
        @property
        def mycomm_cores(self):
            """
            Element mycomm_cores ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 62
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__mycomm_cores(self._handle)
            if array_handle in self._arrays:
                mycomm_cores = self._arrays[array_handle]
            else:
                mycomm_cores = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__mycomm_cores)
                self._arrays[array_handle] = mycomm_cores
            return mycomm_cores
        
        @mycomm_cores.setter
        def mycomm_cores(self, mycomm_cores):
            self.mycomm_cores[...] = mycomm_cores
        
        @property
        def mycomm_size(self):
            """
            Element mycomm_size ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 63
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__mycomm_size(self._handle)
            if array_handle in self._arrays:
                mycomm_size = self._arrays[array_handle]
            else:
                mycomm_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__mycomm_size)
                self._arrays[array_handle] = mycomm_size
            return mycomm_size
        
        @mycomm_size.setter
        def mycomm_size(self, mycomm_size):
            self.mycomm_size[...] = mycomm_size
        
        @property
        def mysend_size(self):
            """
            Element mysend_size ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 64
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__mysend_size(self._handle)
            if array_handle in self._arrays:
                mysend_size = self._arrays[array_handle]
            else:
                mysend_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__mysend_size)
                self._arrays[array_handle] = mysend_size
            return mysend_size
        
        @mysend_size.setter
        def mysend_size(self, mysend_size):
            self.mysend_size[...] = mysend_size
        
        @property
        def local_map(self):
            """
            Element local_map ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 65
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__local_map(self._handle)
            if array_handle in self._arrays:
                local_map = self._arrays[array_handle]
            else:
                local_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__local_map)
                self._arrays[array_handle] = local_map
            return local_map
        
        @local_map.setter
        def local_map(self, local_map):
            self.local_map[...] = local_map
        
        @property
        def local_map1d(self):
            """
            Element local_map1d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 66
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__local_map1d(self._handle)
            if array_handle in self._arrays:
                local_map1d = self._arrays[array_handle]
            else:
                local_map1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__local_map1d)
                self._arrays[array_handle] = local_map1d
            return local_map1d
        
        @local_map1d.setter
        def local_map1d(self, local_map1d):
            self.local_map1d[...] = local_map1d
        
        @property
        def boundary(self):
            """
            Element boundary ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 67
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__boundary(self._handle)
            if array_handle in self._arrays:
                boundary = self._arrays[array_handle]
            else:
                boundary = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__boundary)
                self._arrays[array_handle] = boundary
            return boundary
        
        @boundary.setter
        def boundary(self, boundary):
            self.boundary[...] = boundary
        
        @property
        def boundary1d(self):
            """
            Element boundary1d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 68
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__boundary1d(self._handle)
            if array_handle in self._arrays:
                boundary1d = self._arrays[array_handle]
            else:
                boundary1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__boundary1d)
                self._arrays[array_handle] = boundary1d
            return boundary1d
        
        @boundary1d.setter
        def boundary1d(self, boundary1d):
            self.boundary1d[...] = boundary1d
        
        @property
        def rcount(self):
            """
            Element rcount ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 69
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__rcount(self._handle)
            if array_handle in self._arrays:
                rcount = self._arrays[array_handle]
            else:
                rcount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__rcount)
                self._arrays[array_handle] = rcount
            return rcount
        
        @rcount.setter
        def rcount(self, rcount):
            self.rcount[...] = rcount
        
        @property
        def rdispls(self):
            """
            Element rdispls ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 69
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__rdispls(self._handle)
            if array_handle in self._arrays:
                rdispls = self._arrays[array_handle]
            else:
                rdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__rdispls)
                self._arrays[array_handle] = rdispls
            return rdispls
        
        @rdispls.setter
        def rdispls(self, rdispls):
            self.rdispls[...] = rdispls
        
        @property
        def scount(self):
            """
            Element scount ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 70
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__scount(self._handle)
            if array_handle in self._arrays:
                scount = self._arrays[array_handle]
            else:
                scount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__scount)
                self._arrays[array_handle] = scount
            return scount
        
        @scount.setter
        def scount(self, scount):
            self.scount[...] = scount
        
        @property
        def sdispls(self):
            """
            Element sdispls ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 70
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_diff_map_type__array__sdispls(self._handle)
            if array_handle in self._arrays:
                sdispls = self._arrays[array_handle]
            else:
                sdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_diff_map_type__array__sdispls)
                self._arrays[array_handle] = sdispls
            return sdispls
        
        @sdispls.setter
        def sdispls(self, sdispls):
            self.sdispls[...] = sdispls
        
        def __str__(self):
            ret = ['<grid_diff_map_type>{\n']
            ret.append('    nz_map : ')
            ret.append(repr(self.nz_map))
            ret.append(',\n    mycomm_cores : ')
            ret.append(repr(self.mycomm_cores))
            ret.append(',\n    mycomm_size : ')
            ret.append(repr(self.mycomm_size))
            ret.append(',\n    mysend_size : ')
            ret.append(repr(self.mysend_size))
            ret.append(',\n    local_map : ')
            ret.append(repr(self.local_map))
            ret.append(',\n    local_map1d : ')
            ret.append(repr(self.local_map1d))
            ret.append(',\n    boundary : ')
            ret.append(repr(self.boundary))
            ret.append(',\n    boundary1d : ')
            ret.append(repr(self.boundary1d))
            ret.append(',\n    rcount : ')
            ret.append(repr(self.rcount))
            ret.append(',\n    rdispls : ')
            ret.append(repr(self.rdispls))
            ret.append(',\n    scount : ')
            ret.append(repr(self.scount))
            ret.append(',\n    sdispls : ')
            ret.append(repr(self.sdispls))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.sphere_type")
    class sphere_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=sphere_type)
        
        
        Defined at Smpi_math_module.f90 lines 72-74
        
        """
        def __init__(self, handle=None):
            """
            self = Sphere_Type()
            
            
            Defined at Smpi_math_module.f90 lines 72-74
            
            
            Returns
            -------
            this : Sphere_Type
            	Object to be constructed
            
            
            Automatically generated constructor for sphere_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_sphere_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Sphere_Type
            
            
            Defined at Smpi_math_module.f90 lines 72-74
            
            Parameters
            ----------
            this : Sphere_Type
            	Object to be destructed
            
            
            Automatically generated destructor for sphere_type
            """
            if self._alloc:
                _arespy.f90wrap_sphere_type_finalise(this=self._handle)
        
        @property
        def length(self):
            """
            Element length ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 73
            
            """
            return _arespy.f90wrap_sphere_type__get__length(self._handle)
        
        @length.setter
        def length(self, length):
            _arespy.f90wrap_sphere_type__set__length(self._handle, length)
        
        @property
        def map3d(self):
            """
            Element map3d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.f90 line 74
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_sphere_type__array__map3d(self._handle)
            if array_handle in self._arrays:
                map3d = self._arrays[array_handle]
            else:
                map3d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_sphere_type__array__map3d)
                self._arrays[array_handle] = map3d
            return map3d
        
        @map3d.setter
        def map3d(self, map3d):
            self.map3d[...] = map3d
        
        def __str__(self):
            ret = ['<sphere_type>{\n']
            ret.append('    length : ')
            ret.append(repr(self.length))
            ret.append(',\n    map3d : ')
            ret.append(repr(self.map3d))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def smpi_init():
        """
        smpi_init()
        
        
        Defined at Smpi_math_module.f90 lines 124-135
        
        
        """
        _arespy.f90wrap_smpi_init()
    
    @staticmethod
    def smpi_init_pbc():
        """
        smpi_init_pbc()
        
        
        Defined at Smpi_math_module.f90 lines 138-184
        
        
        """
        _arespy.f90wrap_smpi_init_pbc()
    
    @staticmethod
    def smpi_exit():
        """
        smpi_exit()
        
        
        Defined at Smpi_math_module.f90 lines 188-192
        
        
        """
        _arespy.f90wrap_smpi_exit()
    
    @staticmethod
    def smpi_stop(message):
        """
        smpi_stop(message)
        
        
        Defined at Smpi_math_module.f90 lines 195-200
        
        Parameters
        ----------
        message : str
        
        """
        _arespy.f90wrap_smpi_stop(message=message)
    
    @staticmethod
    def smpi_stop_info(message):
        """
        smpi_stop_info(message)
        
        
        Defined at Smpi_math_module.f90 lines 203-208
        
        Parameters
        ----------
        message : str
        
        """
        _arespy.f90wrap_smpi_stop_info(message=message)
    
    @staticmethod
    def nstates_split(m, np):
        """
        nstates_split(m, np)
        
        
        Defined at Smpi_math_module.f90 lines 212-231
        
        Parameters
        ----------
        m : int
        np : int
        
        -------------------------------------
         ms = mod(m ,np)
         if(ms == 0 ) then
           m = m /np
         else
        if( ms >= (parallel%myid + 1) ) then
          m = m/parallel%numprocs+1
        else
          m = m / parallel%numprocs
        end if
        if(parallel%myid == 0 ) then
           if(parallel%coords(1) == 0 ) then
             m = m/np + ms
           else
             m = m / np
           end if
         end if
        """
        _arespy.f90wrap_nstates_split(m=m, np=np)
    
    @staticmethod
    def nstates_split_2(m, np):
        """
        nstates_split_2(m, np)
        
        
        Defined at Smpi_math_module.f90 lines 235-254
        
        Parameters
        ----------
        m : int
        np : int
        
        -------------------------------------
         ms = mod(m ,np)
         if(ms == 0 ) then
           m = m /np
         else
        if( ms >= (parallel%myid + 1) ) then
          m = m/parallel%numprocs+1
        else
          m = m / parallel%numprocs
        end if
        if(parallel%myid == 0 ) then
           if(parallel%coords(2) == 0 ) then
             m = m/np + ms
           else
             m = m / np
           end if
         end if
        """
        _arespy.f90wrap_nstates_split_2(m=m, np=np)
    
    @staticmethod
    def smpi_reduce_sum_real(amat, na, ramat=None):
        """
        smpi_reduce_sum_real(amat, na[, ramat])
        
        
        Defined at Smpi_math_module.f90 lines 509-518
        
        Parameters
        ----------
        amat : float array
        na : int
        ramat : float array
        
        """
        _arespy.f90wrap_smpi_reduce_sum_real(amat=amat, na=na, ramat=ramat)
    
    @staticmethod
    def start_time(inlabel, flag, tic=None):
        """
        start_time(inlabel, flag[, tic])
        
        
        Defined at Smpi_math_module.f90 lines 550-570
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        tic : float
        
        """
        _arespy.f90wrap_start_time(inlabel=inlabel, flag=flag, tic=tic)
    
    @staticmethod
    def end_time(inlabel, flag, toc=None):
        """
        end_time(inlabel, flag[, toc])
        
        
        Defined at Smpi_math_module.f90 lines 574-593
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        toc : float
        
        """
        _arespy.f90wrap_end_time(inlabel=inlabel, flag=flag, toc=toc)
    
    @staticmethod
    def write_time(inlabel, flag):
        """
        write_time(inlabel, flag)
        
        
        Defined at Smpi_math_module.f90 lines 597-609
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        
        """
        _arespy.f90wrap_write_time(inlabel=inlabel, flag=flag)
    
    @staticmethod
    def write_sum_time(inlabel, flag):
        """
        write_sum_time(inlabel, flag)
        
        
        Defined at Smpi_math_module.f90 lines 613-625
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        
        """
        _arespy.f90wrap_write_sum_time(inlabel=inlabel, flag=flag)
    
    @staticmethod
    def print_time(inlabel, t):
        """
        print_time(inlabel, t)
        
        
        Defined at Smpi_math_module.f90 lines 629-640
        
        Parameters
        ----------
        inlabel : str
        t : float
        
        """
        _arespy.f90wrap_print_time(inlabel=inlabel, t=t)
    
    @staticmethod
    def states_split(nev):
        """
        states_split(nev)
        
        
        Defined at Smpi_math_module.f90 lines 680-700
        
        Parameters
        ----------
        nev : int
        
        """
        _arespy.f90wrap_states_split(nev=nev)
    
    @staticmethod
    def array_split(nev):
        """
        array_split(nev)
        
        
        Defined at Smpi_math_module.f90 lines 703-751
        
        Parameters
        ----------
        nev : int
        
        """
        _arespy.f90wrap_array_split(nev=nev)
    
    @staticmethod
    def grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
        gridrange_sum=None, n1=None, n2=None, n3=None, n=None):
        """
        grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs[, \
            gridrange_sum, n1, n2, n3, n])
        
        
        Defined at Smpi_math_module.f90 lines 754-809
        
        Parameters
        ----------
        ngrid : int
        ncore : int
        comm : int
        id : int
        grid_range : int array
        recvcounts : int array
        displs : int array
        gridrange_sum : int array
        n1 : int
        n2 : int
        n3 : int
        n : int
        
        """
        _arespy.f90wrap_grid_split(ngrid=ngrid, ncore=ncore, comm=comm, id=id, \
            grid_range=grid_range, recvcounts=recvcounts, displs=displs, \
            gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)
    
    @staticmethod
    def atom_split(mysize, natom, atom_index):
        """
        atom_split(mysize, natom, atom_index)
        
        
        Defined at Smpi_math_module.f90 lines 812-848
        
        Parameters
        ----------
        mysize : int
        natom : int
        atom_index : int array
        
        """
        _arespy.f90wrap_atom_split(mysize=mysize, natom=natom, atom_index=atom_index)
    
    @staticmethod
    def grid_sphere_init(n1, n2, n3, norder):
        """
        grid_sphere_init(n1, n2, n3, norder)
        
        
        Defined at Smpi_math_module.f90 lines 850-967
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        norder : int
        
        """
        _arespy.f90wrap_grid_sphere_init(n1=n1, n2=n2, n3=n3, norder=norder)
    
    @staticmethod
    def set_wrap_grid_iso(myrho, wrap_box):
        """
        set_wrap_grid_iso(myrho, wrap_box)
        
        
        Defined at Smpi_math_module.f90 lines 969-1040
        
        Parameters
        ----------
        myrho : float array
        wrap_box : float array
        
        """
        _arespy.f90wrap_set_wrap_grid_iso(myrho=myrho, wrap_box=wrap_box)
    
    @staticmethod
    def destroy_diff_map():
        """
        destroy_diff_map()
        
        
        Defined at Smpi_math_module.f90 lines 1042-1086
        
        
        """
        _arespy.f90wrap_destroy_diff_map()
    
    @staticmethod
    def smpi_diff_init_sph(n1, n2, n3, n, norder, cell_mu, lsp):
        """
        smpi_diff_init_sph(n1, n2, n3, n, norder, cell_mu, lsp)
        
        
        Defined at Smpi_math_module.f90 lines 1088-1460
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        n : int
        norder : int
        cell_mu : int array
        lsp : bool array
        
        """
        _arespy.f90wrap_smpi_diff_init_sph(n1=n1, n2=n2, n3=n3, n=n, norder=norder, \
            cell_mu=cell_mu, lsp=lsp)
    
    @staticmethod
    def smpi_diff_init(n1, n2, n3, n, norder, cell_mu):
        """
        smpi_diff_init(n1, n2, n3, n, norder, cell_mu)
        
        
        Defined at Smpi_math_module.f90 lines 1462-1794
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        n : int
        norder : int
        cell_mu : int array
        
        """
        _arespy.f90wrap_smpi_diff_init(n1=n1, n2=n2, n3=n3, n=n, norder=norder, \
            cell_mu=cell_mu)
    
    @staticmethod
    def set_wrap_sph_pbc_ata_real(myrho, wrap_box1d):
        """
        set_wrap_sph_pbc_ata_real(myrho, wrap_box1d)
        
        
        Defined at Smpi_math_module.f90 lines 1796-1815
        
        Parameters
        ----------
        myrho : float array
        wrap_box1d : float array
        
        """
        _arespy.f90wrap_set_wrap_sph_pbc_ata_real(myrho=myrho, wrap_box1d=wrap_box1d)
    
    @staticmethod
    def set_fft_alltoallv(fft_grid_range_temp):
        """
        set_fft_alltoallv(fft_grid_range_temp)
        
        
        Defined at Smpi_math_module.f90 lines 2042-2085
        
        Parameters
        ----------
        fft_grid_range_temp : int array
        
        """
        _arespy.f90wrap_set_fft_alltoallv(fft_grid_range_temp=fft_grid_range_temp)
    
    @staticmethod
    def destroy_fft_alltoallv():
        """
        destroy_fft_alltoallv()
        
        
        Defined at Smpi_math_module.f90 lines 2087-2100
        
        
        """
        _arespy.f90wrap_destroy_fft_alltoallv()
    
    @staticmethod
    def _sum_real_1d(amat):
        """
        totals = _sum_real_1d(amat)
        
        
        Defined at Smpi_math_module.f90 lines 258-274
        
        Parameters
        ----------
        amat : float array
        
        Returns
        -------
        totals : float
        
        """
        totals = _arespy.f90wrap_sum_real_1d(amat=amat)
        return totals
    
    @staticmethod
    def _sum_real_2d(amat, bmat):
        """
        totals = _sum_real_2d(amat, bmat)
        
        
        Defined at Smpi_math_module.f90 lines 298-315
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        
        Returns
        -------
        totals : float
        
        """
        totals = _arespy.f90wrap_sum_real_2d(amat=amat, bmat=bmat)
        return totals
    
    @staticmethod
    def _sum_real_3d(amat, bmat, cmat):
        """
        totals = _sum_real_3d(amat, bmat, cmat)
        
        
        Defined at Smpi_math_module.f90 lines 340-357
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        cmat : float array
        
        Returns
        -------
        totals : float
        
        """
        totals = _arespy.f90wrap_sum_real_3d(amat=amat, bmat=bmat, cmat=cmat)
        return totals
    
    @staticmethod
    def _sum_cplx_1d(amat):
        """
        totals = _sum_cplx_1d(amat)
        
        
        Defined at Smpi_math_module.f90 lines 278-294
        
        Parameters
        ----------
        amat : complex array
        
        Returns
        -------
        totals : complex
        
        """
        totals = _arespy.f90wrap_sum_cplx_1d(amat=amat)
        return totals
    
    @staticmethod
    def _sum_cplx_2d(amat, bmat):
        """
        totals = _sum_cplx_2d(amat, bmat)
        
        
        Defined at Smpi_math_module.f90 lines 319-336
        
        Parameters
        ----------
        amat : complex array
        bmat : complex array
        
        Returns
        -------
        totals : complex
        
        """
        totals = _arespy.f90wrap_sum_cplx_2d(amat=amat, bmat=bmat)
        return totals
    
    @staticmethod
    def _sum_cplx_3d(amat, bmat, cmat):
        """
        totals = _sum_cplx_3d(amat, bmat, cmat)
        
        
        Defined at Smpi_math_module.f90 lines 361-378
        
        Parameters
        ----------
        amat : complex array
        bmat : complex array
        cmat : complex array
        
        Returns
        -------
        totals : complex
        
        """
        totals = _arespy.f90wrap_sum_cplx_3d(amat=amat, bmat=bmat, cmat=cmat)
        return totals
    
    @staticmethod
    def sompsum(*args, **kwargs):
        """
        sompsum(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 83-89
        
        Overloaded interface containing the following procedures:
          _sum_real_1d
          _sum_real_2d
          _sum_real_3d
          _sum_cplx_1d
          _sum_cplx_2d
          _sum_cplx_3d
        
        """
        for proc in [Smpi_Math_Module._sum_real_1d, Smpi_Math_Module._sum_real_2d, \
            Smpi_Math_Module._sum_real_3d, Smpi_Math_Module._sum_cplx_1d, \
            Smpi_Math_Module._sum_cplx_2d, Smpi_Math_Module._sum_cplx_3d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _smpi_sum_int_1s(x):
        """
        sumx = _smpi_sum_int_1s(x)
        
        
        Defined at Smpi_math_module.f90 lines 424-426
        
        Parameters
        ----------
        x : int
        
        Returns
        -------
        sumx : int
        
        """
        sumx = _arespy.f90wrap_smpi_sum_int_1s(x=x)
        return sumx
    
    @staticmethod
    def _smpi_sum_cplx_1s(x):
        """
        sumx = _smpi_sum_cplx_1s(x)
        
        
        Defined at Smpi_math_module.f90 lines 430-432
        
        Parameters
        ----------
        x : complex
        
        Returns
        -------
        sumx : complex
        
        """
        sumx = _arespy.f90wrap_smpi_sum_cplx_1s(x=x)
        return sumx
    
    @staticmethod
    def _smpi_sum_real_1s(x):
        """
        sumx = _smpi_sum_real_1s(x)
        
        
        Defined at Smpi_math_module.f90 lines 436-438
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        sumx : float
        
        """
        sumx = _arespy.f90wrap_smpi_sum_real_1s(x=x)
        return sumx
    
    @staticmethod
    def _smpi_sum_real_1d(amat):
        """
        suma = _smpi_sum_real_1d(amat)
        
        
        Defined at Smpi_math_module.f90 lines 442-445
        
        Parameters
        ----------
        amat : float array
        
        Returns
        -------
        suma : float
        
        """
        suma = _arespy.f90wrap_smpi_sum_real_1d(amat=amat)
        return suma
    
    @staticmethod
    def _smpi_sum_real_2d(amat, bmat):
        """
        suma = _smpi_sum_real_2d(amat, bmat)
        
        
        Defined at Smpi_math_module.f90 lines 449-453
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        
        Returns
        -------
        suma : float
        
        """
        suma = _arespy.f90wrap_smpi_sum_real_2d(amat=amat, bmat=bmat)
        return suma
    
    @staticmethod
    def _smpi_sum_real_3d(amat, bmat, cmat):
        """
        suma = _smpi_sum_real_3d(amat, bmat, cmat)
        
        
        Defined at Smpi_math_module.f90 lines 457-461
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        cmat : float array
        
        Returns
        -------
        suma : float
        
        """
        suma = _arespy.f90wrap_smpi_sum_real_3d(amat=amat, bmat=bmat, cmat=cmat)
        return suma
    
    @staticmethod
    def smpisum(*args, **kwargs):
        """
        smpisum(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 91-97
        
        Overloaded interface containing the following procedures:
          _smpi_sum_int_1s
          _smpi_sum_cplx_1s
          _smpi_sum_real_1s
          _smpi_sum_real_1d
          _smpi_sum_real_2d
          _smpi_sum_real_3d
        
        """
        for proc in [Smpi_Math_Module._smpi_sum_int_1s, \
            Smpi_Math_Module._smpi_sum_cplx_1s, Smpi_Math_Module._smpi_sum_real_1s, \
            Smpi_Math_Module._smpi_sum_real_1d, Smpi_Math_Module._smpi_sum_real_2d, \
            Smpi_Math_Module._smpi_sum_real_3d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _smpi_sum_mem_1d(munit, amat):
        """
        summem = _smpi_sum_mem_1d(munit, amat)
        
        
        Defined at Smpi_math_module.f90 lines 644-652
        
        Parameters
        ----------
        munit : str
        amat : float array
        
        Returns
        -------
        summem : float
        
        """
        summem = _arespy.f90wrap_smpi_sum_mem_1d(munit=munit, amat=amat)
        return summem
    
    @staticmethod
    def _smpi_sum_mem_2d(munit, amat):
        """
        summem = _smpi_sum_mem_2d(munit, amat)
        
        
        Defined at Smpi_math_module.f90 lines 656-664
        
        Parameters
        ----------
        munit : str
        amat : float array
        
        Returns
        -------
        summem : float
        
        """
        summem = _arespy.f90wrap_smpi_sum_mem_2d(munit=munit, amat=amat)
        return summem
    
    @staticmethod
    def _smpi_sum_mem_3d(munit, amat):
        """
        summem = _smpi_sum_mem_3d(munit, amat)
        
        
        Defined at Smpi_math_module.f90 lines 668-676
        
        Parameters
        ----------
        munit : str
        amat : float array
        
        Returns
        -------
        summem : float
        
        """
        summem = _arespy.f90wrap_smpi_sum_mem_3d(munit=munit, amat=amat)
        return summem
    
    @staticmethod
    def smpisummem(*args, **kwargs):
        """
        smpisummem(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 99-102
        
        Overloaded interface containing the following procedures:
          _smpi_sum_mem_1d
          _smpi_sum_mem_2d
          _smpi_sum_mem_3d
        
        """
        for proc in [Smpi_Math_Module._smpi_sum_mem_1d, \
            Smpi_Math_Module._smpi_sum_mem_2d, Smpi_Math_Module._smpi_sum_mem_3d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _smpi_reduce_sum_real_1d(amat, ramat=None):
        """
        _smpi_reduce_sum_real_1d(amat[, ramat])
        
        
        Defined at Smpi_math_module.f90 lines 495-505
        
        Parameters
        ----------
        amat : float array
        ramat : float array
        
        """
        _arespy.f90wrap_smpi_reduce_sum_real_1d(amat=amat, ramat=ramat)
    
    @staticmethod
    def _smpi_reduce_sum_int_1d(amat, ramat=None):
        """
        _smpi_reduce_sum_int_1d(amat[, ramat])
        
        
        Defined at Smpi_math_module.f90 lines 481-491
        
        Parameters
        ----------
        amat : int array
        ramat : int array
        
        """
        _arespy.f90wrap_smpi_reduce_sum_int_1d(amat=amat, ramat=ramat)
    
    @staticmethod
    def _smpi_reduce_sum_cplx_1d(amat, ramat=None):
        """
        _smpi_reduce_sum_cplx_1d(amat[, ramat])
        
        
        Defined at Smpi_math_module.f90 lines 522-532
        
        Parameters
        ----------
        amat : complex array
        ramat : complex array
        
        """
        _arespy.f90wrap_smpi_reduce_sum_cplx_1d(amat=amat, ramat=ramat)
    
    @staticmethod
    def _smpi_reduce_sum_real_2d(amat, ramat=None):
        """
        _smpi_reduce_sum_real_2d(amat[, ramat])
        
        
        Defined at Smpi_math_module.f90 lines 536-546
        
        Parameters
        ----------
        amat : float array
        ramat : float array
        
        """
        _arespy.f90wrap_smpi_reduce_sum_real_2d(amat=amat, ramat=ramat)
    
    @staticmethod
    def smpireducesum(*args, **kwargs):
        """
        smpireducesum(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 104-108
        
        Overloaded interface containing the following procedures:
          _smpi_reduce_sum_real_1d
          _smpi_reduce_sum_int_1d
          _smpi_reduce_sum_cplx_1d
          _smpi_reduce_sum_real_2d
        
        """
        for proc in [Smpi_Math_Module._smpi_reduce_sum_real_1d, \
            Smpi_Math_Module._smpi_reduce_sum_int_1d, \
            Smpi_Math_Module._smpi_reduce_sum_cplx_1d, \
            Smpi_Math_Module._smpi_reduce_sum_real_2d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _sum_pow_int(amat, pow):
        """
        totals = _sum_pow_int(amat, pow)
        
        
        Defined at Smpi_math_module.f90 lines 403-420
        
        Parameters
        ----------
        amat : float array
        pow : int
        
        Returns
        -------
        totals : float
        
        """
        totals = _arespy.f90wrap_sum_pow_int(amat=amat, pow=pow)
        return totals
    
    @staticmethod
    def _sum_pow_real(amat, pow):
        """
        totals = _sum_pow_real(amat, pow)
        
        
        Defined at Smpi_math_module.f90 lines 382-399
        
        Parameters
        ----------
        amat : float array
        pow : float
        
        Returns
        -------
        totals : float
        
        """
        totals = _arespy.f90wrap_sum_pow_real(amat=amat, pow=pow)
        return totals
    
    @staticmethod
    def sompsumpow(*args, **kwargs):
        """
        sompsumpow(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 110-112
        
        Overloaded interface containing the following procedures:
          _sum_pow_int
          _sum_pow_real
        
        """
        for proc in [Smpi_Math_Module._sum_pow_int, Smpi_Math_Module._sum_pow_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _smpi_sum_pow_int(amat, pow):
        """
        suma = _smpi_sum_pow_int(amat, pow)
        
        
        Defined at Smpi_math_module.f90 lines 473-477
        
        Parameters
        ----------
        amat : float array
        pow : int
        
        Returns
        -------
        suma : float
        
        """
        suma = _arespy.f90wrap_smpi_sum_pow_int(amat=amat, pow=pow)
        return suma
    
    @staticmethod
    def _smpi_sum_pow_real(amat, pow):
        """
        suma = _smpi_sum_pow_real(amat, pow)
        
        
        Defined at Smpi_math_module.f90 lines 465-469
        
        Parameters
        ----------
        amat : float array
        pow : float
        
        Returns
        -------
        suma : float
        
        """
        suma = _arespy.f90wrap_smpi_sum_pow_real(amat=amat, pow=pow)
        return suma
    
    @staticmethod
    def smpisumpow(*args, **kwargs):
        """
        smpisumpow(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 114-115
        
        Overloaded interface containing the following procedures:
          _smpi_sum_pow_int
          _smpi_sum_pow_real
        
        """
        for proc in [Smpi_Math_Module._smpi_sum_pow_int, \
            Smpi_Math_Module._smpi_sum_pow_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _set_wrap_grid_pbc_ata_cmplx(myrho, wrap_box1d, global_n, global_n1, \
        global_n2):
        """
        _set_wrap_grid_pbc_ata_cmplx(myrho, wrap_box1d, global_n, global_n1, global_n2)
        
        
        Defined at Smpi_math_module.f90 lines 1817-1964
        
        Parameters
        ----------
        myrho : complex array
        wrap_box1d : complex array
        global_n : int
        global_n1 : int
        global_n2 : int
        
        """
        _arespy.f90wrap_set_wrap_grid_pbc_ata_cmplx(myrho=myrho, wrap_box1d=wrap_box1d, \
            global_n=global_n, global_n1=global_n1, global_n2=global_n2)
    
    @staticmethod
    def _set_wrap_grid_pbc_ata_real(myrho, wrap_box1d, global_n, global_n1, \
        global_n2):
        """
        _set_wrap_grid_pbc_ata_real(myrho, wrap_box1d, global_n, global_n1, global_n2)
        
        
        Defined at Smpi_math_module.f90 lines 1968-2040
        
        Parameters
        ----------
        myrho : float array
        wrap_box1d : float array
        global_n : int
        global_n1 : int
        global_n2 : int
        
        """
        _arespy.f90wrap_set_wrap_grid_pbc_ata_real(myrho=myrho, wrap_box1d=wrap_box1d, \
            global_n=global_n, global_n1=global_n1, global_n2=global_n2)
    
    @staticmethod
    def set_wrap_grid_pbc_ata(*args, **kwargs):
        """
        set_wrap_grid_pbc_ata(*args, **kwargs)
        
        
        Defined at Smpi_math_module.f90 lines 118-120
        
        Overloaded interface containing the following procedures:
          _set_wrap_grid_pbc_ata_cmplx
          _set_wrap_grid_pbc_ata_real
        
        """
        for proc in [Smpi_Math_Module._set_wrap_grid_pbc_ata_cmplx, \
            Smpi_Math_Module._set_wrap_grid_pbc_ata_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def rtic(self):
        """
        Element rtic ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.f90 line 58
        
        """
        return _arespy.f90wrap_smpi_math_module__get__rtic()
    
    @rtic.setter
    def rtic(self, rtic):
        _arespy.f90wrap_smpi_math_module__set__rtic(rtic)
    
    @property
    def rtoc(self):
        """
        Element rtoc ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.f90 line 58
        
        """
        return _arespy.f90wrap_smpi_math_module__get__rtoc()
    
    @rtoc.setter
    def rtoc(self, rtoc):
        _arespy.f90wrap_smpi_math_module__set__rtoc(rtoc)
    
    @property
    def mpinfo(self):
        """
        Element mpinfo ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.f90 line 79
        
        """
        return _arespy.f90wrap_smpi_math_module__get__mpinfo()
    
    @mpinfo.setter
    def mpinfo(self, mpinfo):
        _arespy.f90wrap_smpi_math_module__set__mpinfo(mpinfo)
    
    @property
    def smpi_status(self):
        """
        Element smpi_status ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.f90 line 80
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_smpi_math_module__array__smpi_status(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            smpi_status = self._arrays[array_handle]
        else:
            smpi_status = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_smpi_math_module__array__smpi_status)
            self._arrays[array_handle] = smpi_status
        return smpi_status
    
    @smpi_status.setter
    def smpi_status(self, smpi_status):
        self.smpi_status[...] = smpi_status
    
    @property
    def lall_grid(self):
        """
        Element lall_grid ftype=logical pytype=bool
        
        
        Defined at Smpi_math_module.f90 line 81
        
        """
        return _arespy.f90wrap_smpi_math_module__get__lall_grid()
    
    @lall_grid.setter
    def lall_grid(self, lall_grid):
        _arespy.f90wrap_smpi_math_module__set__lall_grid(lall_grid)
    
    def __str__(self):
        ret = ['<smpi_math_module>{\n']
        ret.append('    rtic : ')
        ret.append(repr(self.rtic))
        ret.append(',\n    rtoc : ')
        ret.append(repr(self.rtoc))
        ret.append(',\n    mpinfo : ')
        ret.append(repr(self.mpinfo))
        ret.append(',\n    smpi_status : ')
        ret.append(repr(self.smpi_status))
        ret.append(',\n    lall_grid : ')
        ret.append(repr(self.lall_grid))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

smpi_math_module = Smpi_Math_Module()

class Struct_Module(f90wrap.runtime.FortranModule):
    """
    Module struct_module
    
    
    Defined at Struct_module.f90 lines 7-98
    
    """
    @f90wrap.runtime.register_class("arespy.struct_type")
    class struct_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=struct_type)
        
        
        Defined at Struct_module.f90 lines 15-28
        
        """
        def __init__(self, handle=None):
            """
            self = Struct_Type()
            
            
            Defined at Struct_module.f90 lines 15-28
            
            
            Returns
            -------
            this : Struct_Type
            	Object to be constructed
            
            
            Automatically generated constructor for struct_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_struct_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Struct_Type
            
            
            Defined at Struct_module.f90 lines 15-28
            
            Parameters
            ----------
            this : Struct_Type
            	Object to be destructed
            
            
            Automatically generated destructor for struct_type
            """
            if self._alloc:
                _arespy.f90wrap_struct_type_finalise(this=self._handle)
        
        @property
        def zion(self):
            """
            Element zion ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__zion(self._handle)
            if array_handle in self._arrays:
                zion = self._arrays[array_handle]
            else:
                zion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__zion)
                self._arrays[array_handle] = zion
            return zion
        
        @zion.setter
        def zion(self, zion):
            self.zion[...] = zion
        
        @property
        def nati(self):
            """
            Element nati ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__nati(self._handle)
            if array_handle in self._arrays:
                nati = self._arrays[array_handle]
            else:
                nati = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__nati)
                self._arrays[array_handle] = nati
            return nati
        
        @nati.setter
        def nati(self, nati):
            self.nati[...] = nati
        
        @property
        def eleid(self):
            """
            Element eleid ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__eleid(self._handle)
            if array_handle in self._arrays:
                eleid = self._arrays[array_handle]
            else:
                eleid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__eleid)
                self._arrays[array_handle] = eleid
            return eleid
        
        @eleid.setter
        def eleid(self, eleid):
            self.eleid[...] = eleid
        
        @property
        def pos(self):
            """
            Element pos ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__pos(self._handle)
            if array_handle in self._arrays:
                pos = self._arrays[array_handle]
            else:
                pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__pos)
                self._arrays[array_handle] = pos
            return pos
        
        @pos.setter
        def pos(self, pos):
            self.pos[...] = pos
        
        @property
        def poscar(self):
            """
            Element poscar ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__poscar(self._handle)
            if array_handle in self._arrays:
                poscar = self._arrays[array_handle]
            else:
                poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__poscar)
                self._arrays[array_handle] = poscar
            return poscar
        
        @poscar.setter
        def poscar(self, poscar):
            self.poscar[...] = poscar
        
        @property
        def stress(self):
            """
            Element stress ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__stress(self._handle)
            if array_handle in self._arrays:
                stress = self._arrays[array_handle]
            else:
                stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__stress)
                self._arrays[array_handle] = stress
            return stress
        
        @stress.setter
        def stress(self, stress):
            self.stress[...] = stress
        
        @property
        def forces(self):
            """
            Element forces ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 22
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__forces(self._handle)
            if array_handle in self._arrays:
                forces = self._arrays[array_handle]
            else:
                forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__forces)
                self._arrays[array_handle] = forces
            return forces
        
        @forces.setter
        def forces(self, forces):
            self.forces[...] = forces
        
        @property
        def mass(self):
            """
            Element mass ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__mass(self._handle)
            if array_handle in self._arrays:
                mass = self._arrays[array_handle]
            else:
                mass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__mass)
                self._arrays[array_handle] = mass
            return mass
        
        @mass.setter
        def mass(self, mass):
            self.mass[...] = mass
        
        @property
        def zeta(self):
            """
            Element zeta ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__zeta(self._handle)
            if array_handle in self._arrays:
                zeta = self._arrays[array_handle]
            else:
                zeta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__zeta)
                self._arrays[array_handle] = zeta
            return zeta
        
        @zeta.setter
        def zeta(self, zeta):
            self.zeta[...] = zeta
        
        @property
        def prinq(self):
            """
            Element prinq ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__prinq(self._handle)
            if array_handle in self._arrays:
                prinq = self._arrays[array_handle]
            else:
                prinq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__prinq)
                self._arrays[array_handle] = prinq
            return prinq
        
        @prinq.setter
        def prinq(self, prinq):
            self.prinq[...] = prinq
        
        @property
        def lmax(self):
            """
            Element lmax ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 27
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__lmax(self._handle)
            if array_handle in self._arrays:
                lmax = self._arrays[array_handle]
            else:
                lmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__lmax)
                self._arrays[array_handle] = lmax
            return lmax
        
        @lmax.setter
        def lmax(self, lmax):
            self.lmax[...] = lmax
        
        @property
        def elements(self):
            """
            Element elements ftype=character(len=3) pytype=str
            
            
            Defined at Struct_module.f90 line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__elements(self._handle)
            if array_handle in self._arrays:
                elements = self._arrays[array_handle]
            else:
                elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_struct_type__array__elements)
                self._arrays[array_handle] = elements
            return elements
        
        @elements.setter
        def elements(self, elements):
            self.elements[...] = elements
        
        def __str__(self):
            ret = ['<struct_type>{\n']
            ret.append('    zion : ')
            ret.append(repr(self.zion))
            ret.append(',\n    nati : ')
            ret.append(repr(self.nati))
            ret.append(',\n    eleid : ')
            ret.append(repr(self.eleid))
            ret.append(',\n    pos : ')
            ret.append(repr(self.pos))
            ret.append(',\n    poscar : ')
            ret.append(repr(self.poscar))
            ret.append(',\n    stress : ')
            ret.append(repr(self.stress))
            ret.append(',\n    forces : ')
            ret.append(repr(self.forces))
            ret.append(',\n    mass : ')
            ret.append(repr(self.mass))
            ret.append(',\n    zeta : ')
            ret.append(repr(self.zeta))
            ret.append(',\n    prinq : ')
            ret.append(repr(self.prinq))
            ret.append(',\n    lmax : ')
            ret.append(repr(self.lmax))
            ret.append(',\n    elements : ')
            ret.append(repr(self.elements))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def creat_struct(numtyp, numatom):
        """
        creat_struct(numtyp, numatom)
        
        
        Defined at Struct_module.f90 lines 53-77
        
        Parameters
        ----------
        numtyp : int
        numatom : int
        
        """
        _arespy.f90wrap_creat_struct(numtyp=numtyp, numatom=numatom)
    
    @staticmethod
    def destroy_struct():
        """
        destroy_struct()
        
        
        Defined at Struct_module.f90 lines 80-96
        
        
        """
        _arespy.f90wrap_destroy_struct()
    
    @property
    def natom(self):
        """
        Element natom ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 31
        
        """
        return _arespy.f90wrap_struct_module__get__natom()
    
    @natom.setter
    def natom(self, natom):
        _arespy.f90wrap_struct_module__set__natom(natom)
    
    @property
    def nzion(self):
        """
        Element nzion ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 31
        
        """
        return _arespy.f90wrap_struct_module__get__nzion()
    
    @nzion.setter
    def nzion(self, nzion):
        _arespy.f90wrap_struct_module__set__nzion(nzion)
    
    @property
    def naty(self):
        """
        Element naty ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 32
        
        """
        return _arespy.f90wrap_struct_module__get__naty()
    
    @naty.setter
    def naty(self, naty):
        _arespy.f90wrap_struct_module__set__naty(naty)
    
    @property
    def ncharge(self):
        """
        Element ncharge ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 33
        
        """
        return _arespy.f90wrap_struct_module__get__ncharge()
    
    @ncharge.setter
    def ncharge(self, ncharge):
        _arespy.f90wrap_struct_module__set__ncharge(ncharge)
    
    @property
    def charge_ave(self):
        """
        Element charge_ave ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 34
        
        """
        return _arespy.f90wrap_struct_module__get__charge_ave()
    
    @charge_ave.setter
    def charge_ave(self, charge_ave):
        _arespy.f90wrap_struct_module__set__charge_ave(charge_ave)
    
    @property
    def volume(self):
        """
        Element volume ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 35
        
        """
        return _arespy.f90wrap_struct_module__get__volume()
    
    @volume.setter
    def volume(self, volume):
        _arespy.f90wrap_struct_module__set__volume(volume)
    
    @property
    def volsp(self):
        """
        Element volsp ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 35
        
        """
        return _arespy.f90wrap_struct_module__get__volsp()
    
    @volsp.setter
    def volsp(self, volsp):
        _arespy.f90wrap_struct_module__set__volsp(volsp)
    
    @property
    def lat_mat(self):
        """
        Element lat_mat ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__lat_mat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lat_mat = self._arrays[array_handle]
        else:
            lat_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__lat_mat)
            self._arrays[array_handle] = lat_mat
        return lat_mat
    
    @lat_mat.setter
    def lat_mat(self, lat_mat):
        self.lat_mat[...] = lat_mat
    
    @property
    def lat_para(self):
        """
        Element lat_para ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__lat_para(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lat_para = self._arrays[array_handle]
        else:
            lat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__lat_para)
            self._arrays[array_handle] = lat_para
        return lat_para
    
    @lat_para.setter
    def lat_para(self, lat_para):
        self.lat_para[...] = lat_para
    
    @property
    def recip_lat(self):
        """
        Element recip_lat ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__recip_lat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            recip_lat = self._arrays[array_handle]
        else:
            recip_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__recip_lat)
            self._arrays[array_handle] = recip_lat
        return recip_lat
    
    @recip_lat.setter
    def recip_lat(self, recip_lat):
        self.recip_lat[...] = recip_lat
    
    @property
    def reclat_para(self):
        """
        Element reclat_para ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__reclat_para(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            reclat_para = self._arrays[array_handle]
        else:
            reclat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__reclat_para)
            self._arrays[array_handle] = reclat_para
        return reclat_para
    
    @reclat_para.setter
    def reclat_para(self, reclat_para):
        self.reclat_para[...] = reclat_para
    
    @property
    def eionion(self):
        """
        Element eionion ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 40
        
        """
        return _arespy.f90wrap_struct_module__get__eionion()
    
    @eionion.setter
    def eionion(self, eionion):
        _arespy.f90wrap_struct_module__set__eionion(eionion)
    
    @property
    def eshift_ps(self):
        """
        Element eshift_ps ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 41
        
        """
        return _arespy.f90wrap_struct_module__get__eshift_ps()
    
    @eshift_ps.setter
    def eshift_ps(self, eshift_ps):
        _arespy.f90wrap_struct_module__set__eshift_ps(eshift_ps)
    
    @property
    def eshift_tot(self):
        """
        Element eshift_tot ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 41
        
        """
        return _arespy.f90wrap_struct_module__get__eshift_tot()
    
    @eshift_tot.setter
    def eshift_tot(self, eshift_tot):
        _arespy.f90wrap_struct_module__set__eshift_tot(eshift_tot)
    
    @property
    def opsym(self):
        """
        Element opsym ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__opsym(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            opsym = self._arrays[array_handle]
        else:
            opsym = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__opsym)
            self._arrays[array_handle] = opsym
        return opsym
    
    @opsym.setter
    def opsym(self, opsym):
        self.opsym[...] = opsym
    
    @property
    def otrans(self):
        """
        Element otrans ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__otrans(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            otrans = self._arrays[array_handle]
        else:
            otrans = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__otrans)
            self._arrays[array_handle] = otrans
        return otrans
    
    @otrans.setter
    def otrans(self, otrans):
        self.otrans[...] = otrans
    
    @property
    def nsym(self):
        """
        Element nsym ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 45
        
        """
        return _arespy.f90wrap_struct_module__get__nsym()
    
    @nsym.setter
    def nsym(self, nsym):
        _arespy.f90wrap_struct_module__set__nsym(nsym)
    
    @property
    def num_t(self):
        """
        Element num_t ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 46
        
        """
        return _arespy.f90wrap_struct_module__get__num_t()
    
    @num_t.setter
    def num_t(self, num_t):
        _arespy.f90wrap_struct_module__set__num_t(num_t)
    
    @property
    def c_i(self):
        """
        Element c_i ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 47
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__c_i(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            c_i = self._arrays[array_handle]
        else:
            c_i = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__c_i)
            self._arrays[array_handle] = c_i
        return c_i
    
    @c_i.setter
    def c_i(self, c_i):
        self.c_i[...] = c_i
    
    @property
    def odet(self):
        """
        Element odet ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 48
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_struct_module__array__odet(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            odet = self._arrays[array_handle]
        else:
            odet = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_struct_module__array__odet)
            self._arrays[array_handle] = odet
        return odet
    
    @odet.setter
    def odet(self, odet):
        self.odet[...] = odet
    
    def __str__(self):
        ret = ['<struct_module>{\n']
        ret.append('    natom : ')
        ret.append(repr(self.natom))
        ret.append(',\n    nzion : ')
        ret.append(repr(self.nzion))
        ret.append(',\n    naty : ')
        ret.append(repr(self.naty))
        ret.append(',\n    ncharge : ')
        ret.append(repr(self.ncharge))
        ret.append(',\n    charge_ave : ')
        ret.append(repr(self.charge_ave))
        ret.append(',\n    volume : ')
        ret.append(repr(self.volume))
        ret.append(',\n    volsp : ')
        ret.append(repr(self.volsp))
        ret.append(',\n    lat_mat : ')
        ret.append(repr(self.lat_mat))
        ret.append(',\n    lat_para : ')
        ret.append(repr(self.lat_para))
        ret.append(',\n    recip_lat : ')
        ret.append(repr(self.recip_lat))
        ret.append(',\n    reclat_para : ')
        ret.append(repr(self.reclat_para))
        ret.append(',\n    eionion : ')
        ret.append(repr(self.eionion))
        ret.append(',\n    eshift_ps : ')
        ret.append(repr(self.eshift_ps))
        ret.append(',\n    eshift_tot : ')
        ret.append(repr(self.eshift_tot))
        ret.append(',\n    opsym : ')
        ret.append(repr(self.opsym))
        ret.append(',\n    otrans : ')
        ret.append(repr(self.otrans))
        ret.append(',\n    nsym : ')
        ret.append(repr(self.nsym))
        ret.append(',\n    num_t : ')
        ret.append(repr(self.num_t))
        ret.append(',\n    c_i : ')
        ret.append(repr(self.c_i))
        ret.append(',\n    odet : ')
        ret.append(repr(self.odet))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

struct_module = Struct_Module()

class Pspot_Module(f90wrap.runtime.FortranModule):
    """
    Module pspot_module
    
    
    Defined at Struct_module.f90 lines 101-1057
    
    """
    @f90wrap.runtime.register_class("arespy.pspot")
    class pspot(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=pspot)
        
        
        Defined at Struct_module.f90 lines 105-145
        
        """
        def __init__(self, handle=None):
            """
            self = Pspot()
            
            
            Defined at Struct_module.f90 lines 105-145
            
            
            Returns
            -------
            this : Pspot
            	Object to be constructed
            
            
            Automatically generated constructor for pspot
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_pspot_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Pspot
            
            
            Defined at Struct_module.f90 lines 105-145
            
            Parameters
            ----------
            this : Pspot
            	Object to be destructed
            
            
            Automatically generated destructor for pspot
            """
            if self._alloc:
                _arespy.f90wrap_pspot_finalise(this=self._handle)
        
        @property
        def zion(self):
            """
            Element zion ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 106
            
            """
            return _arespy.f90wrap_pspot__get__zion(self._handle)
        
        @zion.setter
        def zion(self, zion):
            _arespy.f90wrap_pspot__set__zion(self._handle, zion)
        
        @property
        def numps(self):
            """
            Element numps ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 107
            
            """
            return _arespy.f90wrap_pspot__get__numps(self._handle)
        
        @numps.setter
        def numps(self, numps):
            _arespy.f90wrap_pspot__set__numps(self._handle, numps)
        
        @property
        def qnumps(self):
            """
            Element qnumps ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 108
            
            """
            return _arespy.f90wrap_pspot__get__qnumps(self._handle)
        
        @qnumps.setter
        def qnumps(self, qnumps):
            _arespy.f90wrap_pspot__set__qnumps(self._handle, qnumps)
        
        @property
        def elename(self):
            """
            Element elename ftype=character(len=3) pytype=str
            
            
            Defined at Struct_module.f90 line 109
            
            """
            return _arespy.f90wrap_pspot__get__elename(self._handle)
        
        @elename.setter
        def elename(self, elename):
            _arespy.f90wrap_pspot__set__elename(self._handle, elename)
        
        @property
        def qmax(self):
            """
            Element qmax ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 111
            
            """
            return _arespy.f90wrap_pspot__get__qmax(self._handle)
        
        @qmax.setter
        def qmax(self, qmax):
            _arespy.f90wrap_pspot__set__qmax(self._handle, qmax)
        
        @property
        def qspacing(self):
            """
            Element qspacing ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 112
            
            """
            return _arespy.f90wrap_pspot__get__qspacing(self._handle)
        
        @qspacing.setter
        def qspacing(self, qspacing):
            _arespy.f90wrap_pspot__set__qspacing(self._handle, qspacing)
        
        @property
        def qmesh(self):
            """
            Element qmesh ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 113
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__qmesh(self._handle)
            if array_handle in self._arrays:
                qmesh = self._arrays[array_handle]
            else:
                qmesh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__qmesh)
                self._arrays[array_handle] = qmesh
            return qmesh
        
        @qmesh.setter
        def qmesh(self, qmesh):
            self.qmesh[...] = qmesh
        
        @property
        def vlocq(self):
            """
            Element vlocq ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 114
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__vlocq(self._handle)
            if array_handle in self._arrays:
                vlocq = self._arrays[array_handle]
            else:
                vlocq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__vlocq)
                self._arrays[array_handle] = vlocq
            return vlocq
        
        @vlocq.setter
        def vlocq(self, vlocq):
            self.vlocq[...] = vlocq
        
        @property
        def vlocqs(self):
            """
            Element vlocqs ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 115
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__vlocqs(self._handle)
            if array_handle in self._arrays:
                vlocqs = self._arrays[array_handle]
            else:
                vlocqs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__vlocqs)
                self._arrays[array_handle] = vlocqs
            return vlocqs
        
        @vlocqs.setter
        def vlocqs(self, vlocqs):
            self.vlocqs[...] = vlocqs
        
        @property
        def ddvl_dq2(self):
            """
            Element ddvl_dq2 ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 116
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__ddvl_dq2(self._handle)
            if array_handle in self._arrays:
                ddvl_dq2 = self._arrays[array_handle]
            else:
                ddvl_dq2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__ddvl_dq2)
                self._arrays[array_handle] = ddvl_dq2
            return ddvl_dq2
        
        @ddvl_dq2.setter
        def ddvl_dq2(self, ddvl_dq2):
            self.ddvl_dq2[...] = ddvl_dq2
        
        @property
        def r(self):
            """
            Element r ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 118
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__r(self._handle)
            if array_handle in self._arrays:
                r = self._arrays[array_handle]
            else:
                r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__r)
                self._arrays[array_handle] = r
            return r
        
        @r.setter
        def r(self, r):
            self.r[...] = r
        
        @property
        def rab(self):
            """
            Element rab ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 118
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__rab(self._handle)
            if array_handle in self._arrays:
                rab = self._arrays[array_handle]
            else:
                rab = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__rab)
                self._arrays[array_handle] = rab
            return rab
        
        @rab.setter
        def rab(self, rab):
            self.rab[...] = rab
        
        @property
        def vlocr(self):
            """
            Element vlocr ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 119
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__vlocr(self._handle)
            if array_handle in self._arrays:
                vlocr = self._arrays[array_handle]
            else:
                vlocr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__vlocr)
                self._arrays[array_handle] = vlocr
            return vlocr
        
        @vlocr.setter
        def vlocr(self, vlocr):
            self.vlocr[...] = vlocr
        
        @property
        def nproj(self):
            """
            Element nproj ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 121
            
            """
            return _arespy.f90wrap_pspot__get__nproj(self._handle)
        
        @nproj.setter
        def nproj(self, nproj):
            _arespy.f90wrap_pspot__set__nproj(self._handle, nproj)
        
        @property
        def proj_l(self):
            """
            Element proj_l ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 122
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__proj_l(self._handle)
            if array_handle in self._arrays:
                proj_l = self._arrays[array_handle]
            else:
                proj_l = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__proj_l)
                self._arrays[array_handle] = proj_l
            return proj_l
        
        @proj_l.setter
        def proj_l(self, proj_l):
            self.proj_l[...] = proj_l
        
        @property
        def proj_m(self):
            """
            Element proj_m ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 123
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__proj_m(self._handle)
            if array_handle in self._arrays:
                proj_m = self._arrays[array_handle]
            else:
                proj_m = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__proj_m)
                self._arrays[array_handle] = proj_m
            return proj_m
        
        @proj_m.setter
        def proj_m(self, proj_m):
            self.proj_m[...] = proj_m
        
        @property
        def indx(self):
            """
            Element indx ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 123
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__indx(self._handle)
            if array_handle in self._arrays:
                indx = self._arrays[array_handle]
            else:
                indx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__indx)
                self._arrays[array_handle] = indx
            return indx
        
        @indx.setter
        def indx(self, indx):
            self.indx[...] = indx
        
        @property
        def rcut(self):
            """
            Element rcut ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 124
            
            """
            return _arespy.f90wrap_pspot__get__rcut(self._handle)
        
        @rcut.setter
        def rcut(self, rcut):
            _arespy.f90wrap_pspot__set__rcut(self._handle, rcut)
        
        @property
        def rmax(self):
            """
            Element rmax ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 125
            
            """
            return _arespy.f90wrap_pspot__get__rmax(self._handle)
        
        @rmax.setter
        def rmax(self, rmax):
            _arespy.f90wrap_pspot__set__rmax(self._handle, rmax)
        
        @property
        def rspacing(self):
            """
            Element rspacing ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 126
            
            """
            return _arespy.f90wrap_pspot__get__rspacing(self._handle)
        
        @rspacing.setter
        def rspacing(self, rspacing):
            _arespy.f90wrap_pspot__set__rspacing(self._handle, rspacing)
        
        @property
        def dij(self):
            """
            Element dij ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 127
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__dij(self._handle)
            if array_handle in self._arrays:
                dij = self._arrays[array_handle]
            else:
                dij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__dij)
                self._arrays[array_handle] = dij
            return dij
        
        @dij.setter
        def dij(self, dij):
            self.dij[...] = dij
        
        @property
        def beta_r(self):
            """
            Element beta_r ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 128
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__beta_r(self._handle)
            if array_handle in self._arrays:
                beta_r = self._arrays[array_handle]
            else:
                beta_r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__beta_r)
                self._arrays[array_handle] = beta_r
            return beta_r
        
        @beta_r.setter
        def beta_r(self, beta_r):
            self.beta_r[...] = beta_r
        
        @property
        def dbeta_dr(self):
            """
            Element dbeta_dr ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 129
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__dbeta_dr(self._handle)
            if array_handle in self._arrays:
                dbeta_dr = self._arrays[array_handle]
            else:
                dbeta_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__dbeta_dr)
                self._arrays[array_handle] = dbeta_dr
            return dbeta_dr
        
        @dbeta_dr.setter
        def dbeta_dr(self, dbeta_dr):
            self.dbeta_dr[...] = dbeta_dr
        
        @property
        def lden(self):
            """
            Element lden ftype=logical pytype=bool
            
            
            Defined at Struct_module.f90 line 131
            
            """
            return _arespy.f90wrap_pspot__get__lden(self._handle)
        
        @lden.setter
        def lden(self, lden):
            _arespy.f90wrap_pspot__set__lden(self._handle, lden)
        
        @property
        def denr(self):
            """
            Element denr ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 132
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__denr(self._handle)
            if array_handle in self._arrays:
                denr = self._arrays[array_handle]
            else:
                denr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__denr)
                self._arrays[array_handle] = denr
            return denr
        
        @denr.setter
        def denr(self, denr):
            self.denr[...] = denr
        
        @property
        def dden_dr(self):
            """
            Element dden_dr ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 133
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__dden_dr(self._handle)
            if array_handle in self._arrays:
                dden_dr = self._arrays[array_handle]
            else:
                dden_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__dden_dr)
                self._arrays[array_handle] = dden_dr
            return dden_dr
        
        @dden_dr.setter
        def dden_dr(self, dden_dr):
            self.dden_dr[...] = dden_dr
        
        @property
        def lcore(self):
            """
            Element lcore ftype=logical pytype=bool
            
            
            Defined at Struct_module.f90 line 135
            
            """
            return _arespy.f90wrap_pspot__get__lcore(self._handle)
        
        @lcore.setter
        def lcore(self, lcore):
            _arespy.f90wrap_pspot__set__lcore(self._handle, lcore)
        
        @property
        def denc(self):
            """
            Element denc ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 137
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__denc(self._handle)
            if array_handle in self._arrays:
                denc = self._arrays[array_handle]
            else:
                denc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__denc)
                self._arrays[array_handle] = denc
            return denc
        
        @denc.setter
        def denc(self, denc):
            self.denc[...] = denc
        
        @property
        def ddenc_dr(self):
            """
            Element ddenc_dr ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 137
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__ddenc_dr(self._handle)
            if array_handle in self._arrays:
                ddenc_dr = self._arrays[array_handle]
            else:
                ddenc_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__ddenc_dr)
                self._arrays[array_handle] = ddenc_dr
            return ddenc_dr
        
        @ddenc_dr.setter
        def ddenc_dr(self, ddenc_dr):
            self.ddenc_dr[...] = ddenc_dr
        
        @property
        def rnoverlap(self):
            """
            Element rnoverlap ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 139
            
            """
            return _arespy.f90wrap_pspot__get__rnoverlap(self._handle)
        
        @rnoverlap.setter
        def rnoverlap(self, rnoverlap):
            _arespy.f90wrap_pspot__set__rnoverlap(self._handle, rnoverlap)
        
        @property
        def eps(self):
            """
            Element eps ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 141
            
            """
            return _arespy.f90wrap_pspot__get__eps(self._handle)
        
        @eps.setter
        def eps(self, eps):
            _arespy.f90wrap_pspot__set__eps(self._handle, eps)
        
        @property
        def eae(self):
            """
            Element eae ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 141
            
            """
            return _arespy.f90wrap_pspot__get__eae(self._handle)
        
        @eae.setter
        def eae(self, eae):
            _arespy.f90wrap_pspot__set__eae(self._handle, eae)
        
        @property
        def nwfa(self):
            """
            Element nwfa ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 143
            
            """
            return _arespy.f90wrap_pspot__get__nwfa(self._handle)
        
        @nwfa.setter
        def nwfa(self, nwfa):
            _arespy.f90wrap_pspot__set__nwfa(self._handle, nwfa)
        
        @property
        def wfal(self):
            """
            Element wfal ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 144
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__wfal(self._handle)
            if array_handle in self._arrays:
                wfal = self._arrays[array_handle]
            else:
                wfal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__wfal)
                self._arrays[array_handle] = wfal
            return wfal
        
        @wfal.setter
        def wfal(self, wfal):
            self.wfal[...] = wfal
        
        @property
        def wfar(self):
            """
            Element wfar ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 145
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_pspot__array__wfar(self._handle)
            if array_handle in self._arrays:
                wfar = self._arrays[array_handle]
            else:
                wfar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_pspot__array__wfar)
                self._arrays[array_handle] = wfar
            return wfar
        
        @wfar.setter
        def wfar(self, wfar):
            self.wfar[...] = wfar
        
        def __str__(self):
            ret = ['<pspot>{\n']
            ret.append('    zion : ')
            ret.append(repr(self.zion))
            ret.append(',\n    numps : ')
            ret.append(repr(self.numps))
            ret.append(',\n    qnumps : ')
            ret.append(repr(self.qnumps))
            ret.append(',\n    elename : ')
            ret.append(repr(self.elename))
            ret.append(',\n    qmax : ')
            ret.append(repr(self.qmax))
            ret.append(',\n    qspacing : ')
            ret.append(repr(self.qspacing))
            ret.append(',\n    qmesh : ')
            ret.append(repr(self.qmesh))
            ret.append(',\n    vlocq : ')
            ret.append(repr(self.vlocq))
            ret.append(',\n    vlocqs : ')
            ret.append(repr(self.vlocqs))
            ret.append(',\n    ddvl_dq2 : ')
            ret.append(repr(self.ddvl_dq2))
            ret.append(',\n    r : ')
            ret.append(repr(self.r))
            ret.append(',\n    rab : ')
            ret.append(repr(self.rab))
            ret.append(',\n    vlocr : ')
            ret.append(repr(self.vlocr))
            ret.append(',\n    nproj : ')
            ret.append(repr(self.nproj))
            ret.append(',\n    proj_l : ')
            ret.append(repr(self.proj_l))
            ret.append(',\n    proj_m : ')
            ret.append(repr(self.proj_m))
            ret.append(',\n    indx : ')
            ret.append(repr(self.indx))
            ret.append(',\n    rcut : ')
            ret.append(repr(self.rcut))
            ret.append(',\n    rmax : ')
            ret.append(repr(self.rmax))
            ret.append(',\n    rspacing : ')
            ret.append(repr(self.rspacing))
            ret.append(',\n    dij : ')
            ret.append(repr(self.dij))
            ret.append(',\n    beta_r : ')
            ret.append(repr(self.beta_r))
            ret.append(',\n    dbeta_dr : ')
            ret.append(repr(self.dbeta_dr))
            ret.append(',\n    lden : ')
            ret.append(repr(self.lden))
            ret.append(',\n    denr : ')
            ret.append(repr(self.denr))
            ret.append(',\n    dden_dr : ')
            ret.append(repr(self.dden_dr))
            ret.append(',\n    lcore : ')
            ret.append(repr(self.lcore))
            ret.append(',\n    denc : ')
            ret.append(repr(self.denc))
            ret.append(',\n    ddenc_dr : ')
            ret.append(repr(self.ddenc_dr))
            ret.append(',\n    rnoverlap : ')
            ret.append(repr(self.rnoverlap))
            ret.append(',\n    eps : ')
            ret.append(repr(self.eps))
            ret.append(',\n    eae : ')
            ret.append(repr(self.eae))
            ret.append(',\n    nwfa : ')
            ret.append(repr(self.nwfa))
            ret.append(',\n    wfal : ')
            ret.append(repr(self.wfal))
            ret.append(',\n    wfar : ')
            ret.append(repr(self.wfar))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.attribute")
    class attribute(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=attribute)
        
        
        Defined at Struct_module.f90 lines 158-159
        
        """
        def __init__(self, handle=None):
            """
            self = Attribute()
            
            
            Defined at Struct_module.f90 lines 158-159
            
            
            Returns
            -------
            this : Attribute
            	Object to be constructed
            
            
            Automatically generated constructor for attribute
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_attribute_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Attribute
            
            
            Defined at Struct_module.f90 lines 158-159
            
            Parameters
            ----------
            this : Attribute
            	Object to be destructed
            
            
            Automatically generated destructor for attribute
            """
            if self._alloc:
                _arespy.f90wrap_attribute_finalise(this=self._handle)
        
        @property
        def value(self):
            """
            Element value ftype=character(len=150) pytype=str
            
            
            Defined at Struct_module.f90 line 159
            
            """
            return _arespy.f90wrap_attribute__get__value(self._handle)
        
        @value.setter
        def value(self, value):
            _arespy.f90wrap_attribute__set__value(self._handle, value)
        
        def __str__(self):
            ret = ['<attribute>{\n']
            ret.append('    value : ')
            ret.append(repr(self.value))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def read_pspot(nty, filenames):
        """
        read_pspot(nty, filenames)
        
        
        Defined at Struct_module.f90 lines 165-364
        
        Parameters
        ----------
        nty : int
        filenames : str array
        
        """
        _arespy.f90wrap_read_pspot(nty=nty, filenames=filenames)
    
    @staticmethod
    def read_psupf_atom(ity, filename, ps):
        """
        read_psupf_atom(ity, filename, ps)
        
        
        Defined at Struct_module.f90 lines 367-584
        
        Parameters
        ----------
        ity : int
        filename : str
        ps : Pspot
        
        -------->Search for Local potential
        """
        _arespy.f90wrap_read_psupf_atom(ity=ity, filename=filename, ps=ps._handle)
    
    @staticmethod
    def scan_head(file_unit, title, start_old):
        """
        scan_head(file_unit, title, start_old)
        
        
        Defined at Struct_module.f90 lines 587-645
        
        Parameters
        ----------
        file_unit : int
        title : str
        start_old : bool
        
        """
        _arespy.f90wrap_scan_head(file_unit=file_unit, title=title, start_old=start_old)
    
    @staticmethod
    def scan_tail(file_unit, title):
        """
        scan_tail(file_unit, title)
        
        
        Defined at Struct_module.f90 lines 648-671
        
        Parameters
        ----------
        file_unit : int
        title : str
        
        """
        _arespy.f90wrap_scan_tail(file_unit=file_unit, title=title)
    
    @staticmethod
    def read_pseudo_header():
        """
        ps, nproj = read_pseudo_header()
        
        
        Defined at Struct_module.f90 lines 674-789
        
        
        Returns
        -------
        ps : Pspot
        nproj : int
        
        """
        ps, nproj = _arespy.f90wrap_read_pseudo_header()
        ps = f90wrap.runtime.lookup_class("arespy.pspot").from_handle(ps, alloc=True)
        return ps, nproj
    
    @staticmethod
    def exist_in(string1, string2):
        """
        exist_in = exist_in(string1, string2)
        
        
        Defined at Struct_module.f90 lines 912-935
        
        Parameters
        ----------
        string1 : str
        string2 : str
        
        Returns
        -------
        exist_in : bool
        
        ====== ====== ======
        """
        exist_in = _arespy.f90wrap_exist_in(string1=string1, string2=string2)
        return exist_in
    
    @staticmethod
    def exist_ibegin(string1, string2):
        """
        exist_ibegin = exist_ibegin(string1, string2)
        
        
        Defined at Struct_module.f90 lines 938-961
        
        Parameters
        ----------
        string1 : str
        string2 : str
        
        Returns
        -------
        exist_ibegin : int
        
        ====== ====== ======
        """
        exist_ibegin = _arespy.f90wrap_exist_ibegin(string1=string1, string2=string2)
        return exist_ibegin
    
    @staticmethod
    def read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l):
        """
        read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l)
        
        
        Defined at Struct_module.f90 lines 964-996
        
        Parameters
        ----------
        unit_upf : int
        nl : int
        beta_r : float array
        d0 : float array
        rcut : float
        proj_l : int array
        
        """
        _arespy.f90wrap_read_pseudo_nonlocal(unit_upf=unit_upf, nl=nl, beta_r=beta_r, \
            d0=d0, rcut=rcut, proj_l=proj_l)
    
    @staticmethod
    def read_pseudo_pswfc(unit_upf, nwfc, nps, wfcl, wfcr):
        """
        read_pseudo_pswfc(unit_upf, nwfc, nps, wfcl, wfcr)
        
        
        Defined at Struct_module.f90 lines 999-1024
        
        Parameters
        ----------
        unit_upf : int
        nwfc : int
        nps : int
        wfcl : int array
        wfcr : float array
        
        ---------------------------------------------------------------------
        """
        _arespy.f90wrap_read_pseudo_pswfc(unit_upf=unit_upf, nwfc=nwfc, nps=nps, \
            wfcl=wfcl, wfcr=wfcr)
    
    @staticmethod
    def aep_generator(nz):
        """
        ps = aep_generator(nz)
        
        
        Defined at Struct_module.f90 lines 1027-1054
        
        Parameters
        ----------
        nz : int
        
        Returns
        -------
        ps : Pspot
        
        """
        ps = _arespy.f90wrap_aep_generator(nz=nz)
        ps = f90wrap.runtime.lookup_class("arespy.pspot").from_handle(ps, alloc=True)
        return ps
    
    @staticmethod
    def _get_value_int(char_in, char_find, variable, find_flag):
        """
        _get_value_int(char_in, char_find, variable, find_flag)
        
        
        Defined at Struct_module.f90 lines 792-819
        
        Parameters
        ----------
        char_in : str
        char_find : str
        variable : int
        find_flag : bool
        
        """
        _arespy.f90wrap_get_value_int(char_in=char_in, char_find=char_find, \
            variable=variable, find_flag=find_flag)
    
    @staticmethod
    def _get_value_real(char_in, char_find, variable, find_flag):
        """
        _get_value_real(char_in, char_find, variable, find_flag)
        
        
        Defined at Struct_module.f90 lines 822-849
        
        Parameters
        ----------
        char_in : str
        char_find : str
        variable : float
        find_flag : bool
        
        """
        _arespy.f90wrap_get_value_real(char_in=char_in, char_find=char_find, \
            variable=variable, find_flag=find_flag)
    
    @staticmethod
    def _get_value_char(char_in, char_find, variable, find_flag):
        """
        _get_value_char(char_in, char_find, variable, find_flag)
        
        
        Defined at Struct_module.f90 lines 852-879
        
        Parameters
        ----------
        char_in : str
        char_find : str
        variable : str
        find_flag : bool
        
        """
        _arespy.f90wrap_get_value_char(char_in=char_in, char_find=char_find, \
            variable=variable, find_flag=find_flag)
    
    @staticmethod
    def _get_value_logic(char_in, char_find, variable, find_flag):
        """
        _get_value_logic(char_in, char_find, variable, find_flag)
        
        
        Defined at Struct_module.f90 lines 882-909
        
        Parameters
        ----------
        char_in : str
        char_find : str
        variable : bool
        find_flag : bool
        
        """
        _arespy.f90wrap_get_value_logic(char_in=char_in, char_find=char_find, \
            variable=variable, find_flag=find_flag)
    
    @staticmethod
    def get_value(*args, **kwargs):
        """
        get_value(*args, **kwargs)
        
        
        Defined at Struct_module.f90 lines 154-156
        
        Overloaded interface containing the following procedures:
          _get_value_int
          _get_value_real
          _get_value_char
          _get_value_logic
        
        """
        for proc in [Pspot_Module._get_value_int, Pspot_Module._get_value_real, \
            Pspot_Module._get_value_char, Pspot_Module._get_value_logic]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def max_nproj(self):
        """
        Element max_nproj ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 150
        
        """
        return _arespy.f90wrap_pspot_module__get__max_nproj()
    
    @max_nproj.setter
    def max_nproj(self, max_nproj):
        _arespy.f90wrap_pspot_module__set__max_nproj(max_nproj)
    
    @property
    def max_nwfa(self):
        """
        Element max_nwfa ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 150
        
        """
        return _arespy.f90wrap_pspot_module__get__max_nwfa()
    
    @max_nwfa.setter
    def max_nwfa(self, max_nwfa):
        _arespy.f90wrap_pspot_module__set__max_nwfa(max_nwfa)
    
    @property
    def max_rcut(self):
        """
        Element max_rcut ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 151
        
        """
        return _arespy.f90wrap_pspot_module__get__max_rcut()
    
    @max_rcut.setter
    def max_rcut(self, max_rcut):
        _arespy.f90wrap_pspot_module__set__max_rcut(max_rcut)
    
    def __str__(self):
        ret = ['<pspot_module>{\n']
        ret.append('    max_nproj : ')
        ret.append(repr(self.max_nproj))
        ret.append(',\n    max_nwfa : ')
        ret.append(repr(self.max_nwfa))
        ret.append(',\n    max_rcut : ')
        ret.append(repr(self.max_rcut))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

pspot_module = Pspot_Module()

class Grid_Module(f90wrap.runtime.FortranModule):
    """
    Module grid_module
    
    
    Defined at Struct_module.f90 lines 1060-2168
    
    """
    @f90wrap.runtime.register_class("arespy.grid_type")
    class grid_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=grid_type)
        
        
        Defined at Struct_module.f90 lines 1069-1088
        
        """
        def __init__(self, handle=None):
            """
            self = Grid_Type()
            
            
            Defined at Struct_module.f90 lines 1069-1088
            
            
            Returns
            -------
            this : Grid_Type
            	Object to be constructed
            
            
            Automatically generated constructor for grid_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_grid_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Grid_Type
            
            
            Defined at Struct_module.f90 lines 1069-1088
            
            Parameters
            ----------
            this : Grid_Type
            	Object to be destructed
            
            
            Automatically generated destructor for grid_type
            """
            if self._alloc:
                _arespy.f90wrap_grid_type_finalise(this=self._handle)
        
        @property
        def rhos(self):
            """
            Element rhos ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1071
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__rhos(self._handle)
            if array_handle in self._arrays:
                rhos = self._arrays[array_handle]
            else:
                rhos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__rhos)
                self._arrays[array_handle] = rhos
            return rhos
        
        @rhos.setter
        def rhos(self, rhos):
            self.rhos[...] = rhos
        
        @property
        def rho(self):
            """
            Element rho ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1072
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__rho(self._handle)
            if array_handle in self._arrays:
                rho = self._arrays[array_handle]
            else:
                rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__rho)
                self._arrays[array_handle] = rho
            return rho
        
        @rho.setter
        def rho(self, rho):
            self.rho[...] = rho
        
        @property
        def vxcs(self):
            """
            Element vxcs ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1073
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__vxcs(self._handle)
            if array_handle in self._arrays:
                vxcs = self._arrays[array_handle]
            else:
                vxcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__vxcs)
                self._arrays[array_handle] = vxcs
            return vxcs
        
        @vxcs.setter
        def vxcs(self, vxcs):
            self.vxcs[...] = vxcs
        
        @property
        def vhxcd(self):
            """
            Element vhxcd ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1074
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__vhxcd(self._handle)
            if array_handle in self._arrays:
                vhxcd = self._arrays[array_handle]
            else:
                vhxcd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__vhxcd)
                self._arrays[array_handle] = vhxcd
            return vhxcd
        
        @vhxcd.setter
        def vhxcd(self, vhxcd):
            self.vhxcd[...] = vhxcd
        
        @property
        def vlpp(self):
            """
            Element vlpp ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1075
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__vlpp(self._handle)
            if array_handle in self._arrays:
                vlpp = self._arrays[array_handle]
            else:
                vlpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__vlpp)
                self._arrays[array_handle] = vlpp
            return vlpp
        
        @vlpp.setter
        def vlpp(self, vlpp):
            self.vlpp[...] = vlpp
        
        @property
        def vh(self):
            """
            Element vh ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1076
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__vh(self._handle)
            if array_handle in self._arrays:
                vh = self._arrays[array_handle]
            else:
                vh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__vh)
                self._arrays[array_handle] = vh
            return vh
        
        @vh.setter
        def vh(self, vh):
            self.vh[...] = vh
        
        @property
        def eval(self):
            """
            Element eval ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1077
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__eval(self._handle)
            if array_handle in self._arrays:
                eval = self._arrays[array_handle]
            else:
                eval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__eval)
                self._arrays[array_handle] = eval
            return eval
        
        @eval.setter
        def eval(self, eval):
            self.eval[...] = eval
        
        @property
        def rhoc(self):
            """
            Element rhoc ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1079
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__rhoc(self._handle)
            if array_handle in self._arrays:
                rhoc = self._arrays[array_handle]
            else:
                rhoc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__rhoc)
                self._arrays[array_handle] = rhoc
            return rhoc
        
        @rhoc.setter
        def rhoc(self, rhoc):
            self.rhoc[...] = rhoc
        
        @property
        def gvec(self):
            """
            Element gvec ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1081
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__gvec(self._handle)
            if array_handle in self._arrays:
                gvec = self._arrays[array_handle]
            else:
                gvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__gvec)
                self._arrays[array_handle] = gvec
            return gvec
        
        @gvec.setter
        def gvec(self, gvec):
            self.gvec[...] = gvec
        
        @property
        def gmask(self):
            """
            Element gmask ftype=logical pytype=bool
            
            
            Defined at Struct_module.f90 line 1082
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__gmask(self._handle)
            if array_handle in self._arrays:
                gmask = self._arrays[array_handle]
            else:
                gmask = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__gmask)
                self._arrays[array_handle] = gmask
            return gmask
        
        @gmask.setter
        def gmask(self, gmask):
            self.gmask[...] = gmask
        
        @property
        def rvec(self):
            """
            Element rvec ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1084
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__rvec(self._handle)
            if array_handle in self._arrays:
                rvec = self._arrays[array_handle]
            else:
                rvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__rvec)
                self._arrays[array_handle] = rvec
            return rvec
        
        @rvec.setter
        def rvec(self, rvec):
            self.rvec[...] = rvec
        
        @property
        def lsp(self):
            """
            Element lsp ftype=logical pytype=bool
            
            
            Defined at Struct_module.f90 line 1086
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_grid_type__array__lsp(self._handle)
            if array_handle in self._arrays:
                lsp = self._arrays[array_handle]
            else:
                lsp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_grid_type__array__lsp)
                self._arrays[array_handle] = lsp
            return lsp
        
        @lsp.setter
        def lsp(self, lsp):
            self.lsp[...] = lsp
        
        @property
        def onedlength(self):
            """
            Element onedlength ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.f90 line 1088
            
            """
            return _arespy.f90wrap_grid_type__get__onedlength(self._handle)
        
        @onedlength.setter
        def onedlength(self, onedlength):
            _arespy.f90wrap_grid_type__set__onedlength(self._handle, onedlength)
        
        def __str__(self):
            ret = ['<grid_type>{\n']
            ret.append('    rhos : ')
            ret.append(repr(self.rhos))
            ret.append(',\n    rho : ')
            ret.append(repr(self.rho))
            ret.append(',\n    vxcs : ')
            ret.append(repr(self.vxcs))
            ret.append(',\n    vhxcd : ')
            ret.append(repr(self.vhxcd))
            ret.append(',\n    vlpp : ')
            ret.append(repr(self.vlpp))
            ret.append(',\n    vh : ')
            ret.append(repr(self.vh))
            ret.append(',\n    eval : ')
            ret.append(repr(self.eval))
            ret.append(',\n    rhoc : ')
            ret.append(repr(self.rhoc))
            ret.append(',\n    gvec : ')
            ret.append(repr(self.gvec))
            ret.append(',\n    gmask : ')
            ret.append(repr(self.gmask))
            ret.append(',\n    rvec : ')
            ret.append(repr(self.rvec))
            ret.append(',\n    lsp : ')
            ret.append(repr(self.lsp))
            ret.append(',\n    onedlength : ')
            ret.append(repr(self.onedlength))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.kgrid_type")
    class kgrid_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=kgrid_type)
        
        
        Defined at Struct_module.f90 lines 1091-1095
        
        """
        def __init__(self, handle=None):
            """
            self = Kgrid_Type()
            
            
            Defined at Struct_module.f90 lines 1091-1095
            
            
            Returns
            -------
            this : Kgrid_Type
            	Object to be constructed
            
            
            Automatically generated constructor for kgrid_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_kgrid_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Kgrid_Type
            
            
            Defined at Struct_module.f90 lines 1091-1095
            
            Parameters
            ----------
            this : Kgrid_Type
            	Object to be destructed
            
            
            Automatically generated destructor for kgrid_type
            """
            if self._alloc:
                _arespy.f90wrap_kgrid_type_finalise(this=self._handle)
        
        @property
        def vec(self):
            """
            Element vec ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1093
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_kgrid_type__array__vec(self._handle)
            if array_handle in self._arrays:
                vec = self._arrays[array_handle]
            else:
                vec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_kgrid_type__array__vec)
                self._arrays[array_handle] = vec
            return vec
        
        @vec.setter
        def vec(self, vec):
            self.vec[...] = vec
        
        @property
        def vcar(self):
            """
            Element vcar ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1094
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_kgrid_type__array__vcar(self._handle)
            if array_handle in self._arrays:
                vcar = self._arrays[array_handle]
            else:
                vcar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_kgrid_type__array__vcar)
                self._arrays[array_handle] = vcar
            return vcar
        
        @vcar.setter
        def vcar(self, vcar):
            self.vcar[...] = vcar
        
        @property
        def wk(self):
            """
            Element wk ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1095
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_kgrid_type__array__wk(self._handle)
            if array_handle in self._arrays:
                wk = self._arrays[array_handle]
            else:
                wk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_kgrid_type__array__wk)
                self._arrays[array_handle] = wk
            return wk
        
        @wk.setter
        def wk(self, wk):
            self.wk[...] = wk
        
        def __str__(self):
            ret = ['<kgrid_type>{\n']
            ret.append('    vec : ')
            ret.append(repr(self.vec))
            ret.append(',\n    vcar : ')
            ret.append(repr(self.vcar))
            ret.append(',\n    wk : ')
            ret.append(repr(self.wk))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("arespy.eigen_type")
    class eigen_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=eigen_type)
        
        
        Defined at Struct_module.f90 lines 1098-1103
        
        """
        def __init__(self, handle=None):
            """
            self = Eigen_Type()
            
            
            Defined at Struct_module.f90 lines 1098-1103
            
            
            Returns
            -------
            this : Eigen_Type
            	Object to be constructed
            
            
            Automatically generated constructor for eigen_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _arespy.f90wrap_eigen_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Eigen_Type
            
            
            Defined at Struct_module.f90 lines 1098-1103
            
            Parameters
            ----------
            this : Eigen_Type
            	Object to be destructed
            
            
            Automatically generated destructor for eigen_type
            """
            if self._alloc:
                _arespy.f90wrap_eigen_type_finalise(this=self._handle)
        
        @property
        def wvf(self):
            """
            Element wvf ftype=complex(dcp) pytype=complex
            
            
            Defined at Struct_module.f90 line 1100
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_eigen_type__array__wvf(self._handle)
            if array_handle in self._arrays:
                wvf = self._arrays[array_handle]
            else:
                wvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_eigen_type__array__wvf)
                self._arrays[array_handle] = wvf
            return wvf
        
        @wvf.setter
        def wvf(self, wvf):
            self.wvf[...] = wvf
        
        @property
        def val(self):
            """
            Element val ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1102
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_eigen_type__array__val(self._handle)
            if array_handle in self._arrays:
                val = self._arrays[array_handle]
            else:
                val = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_eigen_type__array__val)
                self._arrays[array_handle] = val
            return val
        
        @val.setter
        def val(self, val):
            self.val[...] = val
        
        @property
        def wvfg(self):
            """
            Element wvfg ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.f90 line 1103
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_eigen_type__array__wvfg(self._handle)
            if array_handle in self._arrays:
                wvfg = self._arrays[array_handle]
            else:
                wvfg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _arespy.f90wrap_eigen_type__array__wvfg)
                self._arrays[array_handle] = wvfg
            return wvfg
        
        @wvfg.setter
        def wvfg(self, wvfg):
            self.wvfg[...] = wvfg
        
        def __str__(self):
            ret = ['<eigen_type>{\n']
            ret.append('    wvf : ')
            ret.append(repr(self.wvf))
            ret.append(',\n    val : ')
            ret.append(repr(self.val))
            ret.append(',\n    wvfg : ')
            ret.append(repr(self.wvfg))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def build_rgrid():
        """
        build_rgrid()
        
        
        Defined at Struct_module.f90 lines 1129-1192
        
        
        """
        _arespy.f90wrap_build_rgrid()
    
    @staticmethod
    def destroy_rgrid():
        """
        destroy_rgrid()
        
        
        Defined at Struct_module.f90 lines 1195-1212
        
        
        """
        _arespy.f90wrap_destroy_rgrid()
    
    @staticmethod
    def build_kgrid():
        """
        build_kgrid()
        
        
        Defined at Struct_module.f90 lines 1215-1247
        
        
        """
        _arespy.f90wrap_build_kgrid()
    
    @staticmethod
    def destroy_kpt():
        """
        destroy_kpt()
        
        
        Defined at Struct_module.f90 lines 1250-1257
        
        
        """
        _arespy.f90wrap_destroy_kpt()
    
    @staticmethod
    def build_eigen():
        """
        build_eigen()
        
        
        Defined at Struct_module.f90 lines 1260-1283
        
        
        """
        _arespy.f90wrap_build_eigen()
    
    @staticmethod
    def destroy_eigen():
        """
        destroy_eigen()
        
        
        Defined at Struct_module.f90 lines 1286-1292
        
        
        """
        _arespy.f90wrap_destroy_eigen()
    
    @staticmethod
    def fillqtable():
        """
        fillqtable()
        
        
        Defined at Struct_module.f90 lines 1295-1353
        
        
        """
        _arespy.f90wrap_fillqtable()
    
    @staticmethod
    def fillrtable():
        """
        fillrtable()
        
        
        Defined at Struct_module.f90 lines 1356-1411
        
        
        """
        _arespy.f90wrap_fillrtable()
    
    @staticmethod
    def symm_kgrid():
        """
        symm_kgrid()
        
        
        Defined at Struct_module.f90 lines 1414-1588
        
        
        -----------------------
        """
        _arespy.f90wrap_symm_kgrid()
    
    @staticmethod
    def symm_density(rho):
        """
        symm_density(rho)
        
        
        Defined at Struct_module.f90 lines 1591-1668
        
        Parameters
        ----------
        rho : float array
        
        """
        _arespy.f90wrap_symm_density(rho=rho)
    
    @staticmethod
    def isymm_apply(isy, nsize, vin, vout):
        """
        lsymm = isymm_apply(isy, nsize, vin, vout)
        
        
        Defined at Struct_module.f90 lines 1671-1730
        
        Parameters
        ----------
        isy : int
        nsize : int array
        vin : int array
        vout : int array
        
        Returns
        -------
        lsymm : bool
        
        """
        lsymm = _arespy.f90wrap_isymm_apply(isy=isy, nsize=nsize, vin=vin, vout=vout)
        return lsymm
    
    @staticmethod
    def build_parallel_3d_grid():
        """
        build_parallel_3d_grid()
        
        
        Defined at Struct_module.f90 lines 1733-1753
        
        
        """
        _arespy.f90wrap_build_parallel_3d_grid()
    
    @staticmethod
    def build_parallel_sph_grid():
        """
        build_parallel_sph_grid()
        
        
        Defined at Struct_module.f90 lines 1757-1770
        
        
        """
        _arespy.f90wrap_build_parallel_sph_grid()
    
    @staticmethod
    def grid_split_sph(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
        gridrange_sum, n1, n2, n3, n):
        """
        grid_split_sph(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
            gridrange_sum, n1, n2, n3, n)
        
        
        Defined at Struct_module.f90 lines 1773-1872
        
        Parameters
        ----------
        ngrid : int
        ncore : int
        comm : int
        id : int
        grid_range : int array
        recvcounts : int array
        displs : int array
        gridrange_sum : int array
        n1 : int
        n2 : int
        n3 : int
        n : int
        
        """
        _arespy.f90wrap_grid_split_sph(ngrid=ngrid, ncore=ncore, comm=comm, id=id, \
            grid_range=grid_range, recvcounts=recvcounts, displs=displs, \
            gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)
    
    @staticmethod
    def sphere_region(n1, n2, n3, lsphere, nz_map, nsphere):
        """
        sphere_region(n1, n2, n3, lsphere, nz_map, nsphere)
        
        
        Defined at Struct_module.f90 lines 1877-1966
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        lsphere : bool array
        nz_map : int array
        nsphere : int
        
        """
        _arespy.f90wrap_sphere_region(n1=n1, n2=n2, n3=n3, lsphere=lsphere, \
            nz_map=nz_map, nsphere=nsphere)
    
    @staticmethod
    def sphere2cubic(nps, f1d, f3d, rfill=None):
        """
        sphere2cubic(nps, f1d, f3d[, rfill])
        
        
        Defined at Struct_module.f90 lines 1970-2018
        
        Parameters
        ----------
        nps : int
        f1d : float array
        f3d : float array
        rfill : float
        
        """
        _arespy.f90wrap_sphere2cubic(nps=nps, f1d=f1d, f3d=f3d, rfill=rfill)
    
    @staticmethod
    def cubic2sphere(nps, f3d, f1d):
        """
        cubic2sphere(nps, f3d, f1d)
        
        
        Defined at Struct_module.f90 lines 2021-2046
        
        Parameters
        ----------
        nps : int
        f3d : float array
        f1d : float array
        
        """
        _arespy.f90wrap_cubic2sphere(nps=nps, f3d=f3d, f1d=f1d)
    
    @staticmethod
    def cubic2sphere_fft(nps, f3d, f1d, shiftn, shiftz):
        """
        cubic2sphere_fft(nps, f3d, f1d, shiftn, shiftz)
        
        
        Defined at Struct_module.f90 lines 2050-2068
        
        Parameters
        ----------
        nps : int
        f3d : float array
        f1d : float array
        shiftn : int
        shiftz : int
        
        """
        _arespy.f90wrap_cubic2sphere_fft(nps=nps, f3d=f3d, f1d=f1d, shiftn=shiftn, \
            shiftz=shiftz)
    
    @staticmethod
    def sphere2cubic_fft(nps, f1d, f3d, shiftn, shiftz):
        """
        sphere2cubic_fft(nps, f1d, f3d, shiftn, shiftz)
        
        
        Defined at Struct_module.f90 lines 2071-2089
        
        Parameters
        ----------
        nps : int
        f1d : float array
        f3d : float array
        shiftn : int
        shiftz : int
        
        """
        _arespy.f90wrap_sphere2cubic_fft(nps=nps, f1d=f1d, f3d=f3d, shiftn=shiftn, \
            shiftz=shiftz)
    
    @staticmethod
    def sumrhos(nps, rhos, rho):
        """
        sumrhos(nps, rhos, rho)
        
        
        Defined at Struct_module.f90 lines 2150-2167
        
        Parameters
        ----------
        nps : int
        rhos : float array
        rho : float array
        
        """
        _arespy.f90wrap_sumrhos(nps=nps, rhos=rhos, rho=rho)
    
    @staticmethod
    def _fft_sph_r2c(array_r):
        """
        array_c = _fft_sph_r2c(array_r)
        
        
        Defined at Struct_module.f90 lines 2092-2118
        
        Parameters
        ----------
        array_r : float array
        
        Returns
        -------
        array_c : complex array
        
        """
        array_c = _arespy.f90wrap_fft_sph_r2c(array_r=array_r)
        return array_c
    
    @staticmethod
    def _fft_sph_c2r(array_c):
        """
        array_r = _fft_sph_c2r(array_c)
        
        
        Defined at Struct_module.f90 lines 2121-2145
        
        Parameters
        ----------
        array_c : complex array
        
        Returns
        -------
        array_r : float array
        
        """
        array_r = _arespy.f90wrap_fft_sph_c2r(array_c=array_c)
        return array_r
    
    @staticmethod
    def fft_sph(*args, **kwargs):
        """
        fft_sph(*args, **kwargs)
        
        
        Defined at Struct_module.f90 lines 1122-1124
        
        Overloaded interface containing the following procedures:
          _fft_sph_r2c
          _fft_sph_c2r
        
        """
        for proc in [Grid_Module._fft_sph_r2c, Grid_Module._fft_sph_c2r]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def n1(self):
        """
        Element n1 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1106
        
        """
        return _arespy.f90wrap_grid_module__get__n1()
    
    @n1.setter
    def n1(self, n1):
        _arespy.f90wrap_grid_module__set__n1(n1)
    
    @property
    def n2(self):
        """
        Element n2 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1106
        
        """
        return _arespy.f90wrap_grid_module__get__n2()
    
    @n2.setter
    def n2(self, n2):
        _arespy.f90wrap_grid_module__set__n2(n2)
    
    @property
    def n3(self):
        """
        Element n3 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1106
        
        """
        return _arespy.f90wrap_grid_module__get__n3()
    
    @n3.setter
    def n3(self, n3):
        _arespy.f90wrap_grid_module__set__n3(n3)
    
    @property
    def n(self):
        """
        Element n ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1106
        
        """
        return _arespy.f90wrap_grid_module__get__n()
    
    @n.setter
    def n(self, n):
        _arespy.f90wrap_grid_module__set__n(n)
    
    @property
    def ng1(self):
        """
        Element ng1 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1107
        
        """
        return _arespy.f90wrap_grid_module__get__ng1()
    
    @ng1.setter
    def ng1(self, ng1):
        _arespy.f90wrap_grid_module__set__ng1(ng1)
    
    @property
    def ng2(self):
        """
        Element ng2 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1107
        
        """
        return _arespy.f90wrap_grid_module__get__ng2()
    
    @ng2.setter
    def ng2(self, ng2):
        _arespy.f90wrap_grid_module__set__ng2(ng2)
    
    @property
    def ng3(self):
        """
        Element ng3 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1107
        
        """
        return _arespy.f90wrap_grid_module__get__ng3()
    
    @ng3.setter
    def ng3(self, ng3):
        _arespy.f90wrap_grid_module__set__ng3(ng3)
    
    @property
    def ng(self):
        """
        Element ng ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1107
        
        """
        return _arespy.f90wrap_grid_module__get__ng()
    
    @ng.setter
    def ng(self, ng):
        _arespy.f90wrap_grid_module__set__ng(ng)
    
    @property
    def gap(self):
        """
        Element gap ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 1108
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_grid_module__array__gap(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gap = self._arrays[array_handle]
        else:
            gap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_grid_module__array__gap)
            self._arrays[array_handle] = gap
        return gap
    
    @gap.setter
    def gap(self, gap):
        self.gap[...] = gap
    
    @property
    def dvol(self):
        """
        Element dvol ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 1109
        
        """
        return _arespy.f90wrap_grid_module__get__dvol()
    
    @dvol.setter
    def dvol(self, dvol):
        _arespy.f90wrap_grid_module__set__dvol(dvol)
    
    @property
    def nk1(self):
        """
        Element nk1 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1111
        
        """
        return _arespy.f90wrap_grid_module__get__nk1()
    
    @nk1.setter
    def nk1(self, nk1):
        _arespy.f90wrap_grid_module__set__nk1(nk1)
    
    @property
    def nk2(self):
        """
        Element nk2 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1111
        
        """
        return _arespy.f90wrap_grid_module__get__nk2()
    
    @nk2.setter
    def nk2(self, nk2):
        _arespy.f90wrap_grid_module__set__nk2(nk2)
    
    @property
    def nk3(self):
        """
        Element nk3 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1111
        
        """
        return _arespy.f90wrap_grid_module__get__nk3()
    
    @nk3.setter
    def nk3(self, nk3):
        _arespy.f90wrap_grid_module__set__nk3(nk3)
    
    @property
    def nk(self):
        """
        Element nk ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1111
        
        """
        return _arespy.f90wrap_grid_module__get__nk()
    
    @nk.setter
    def nk(self, nk):
        _arespy.f90wrap_grid_module__set__nk(nk)
    
    @property
    def kdispl(self):
        """
        Element kdispl ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.f90 line 1112
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_grid_module__array__kdispl(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kdispl = self._arrays[array_handle]
        else:
            kdispl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_grid_module__array__kdispl)
            self._arrays[array_handle] = kdispl
        return kdispl
    
    @kdispl.setter
    def kdispl(self, kdispl):
        self.kdispl[...] = kdispl
    
    @property
    def global_n1(self):
        """
        Element global_n1 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1121
        
        """
        return _arespy.f90wrap_grid_module__get__global_n1()
    
    @global_n1.setter
    def global_n1(self, global_n1):
        _arespy.f90wrap_grid_module__set__global_n1(global_n1)
    
    @property
    def global_n2(self):
        """
        Element global_n2 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1121
        
        """
        return _arespy.f90wrap_grid_module__get__global_n2()
    
    @global_n2.setter
    def global_n2(self, global_n2):
        _arespy.f90wrap_grid_module__set__global_n2(global_n2)
    
    @property
    def global_n3(self):
        """
        Element global_n3 ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1121
        
        """
        return _arespy.f90wrap_grid_module__get__global_n3()
    
    @global_n3.setter
    def global_n3(self, global_n3):
        _arespy.f90wrap_grid_module__set__global_n3(global_n3)
    
    @property
    def global_n(self):
        """
        Element global_n ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.f90 line 1121
        
        """
        return _arespy.f90wrap_grid_module__get__global_n()
    
    @global_n.setter
    def global_n(self, global_n):
        _arespy.f90wrap_grid_module__set__global_n(global_n)
    
    def __str__(self):
        ret = ['<grid_module>{\n']
        ret.append('    n1 : ')
        ret.append(repr(self.n1))
        ret.append(',\n    n2 : ')
        ret.append(repr(self.n2))
        ret.append(',\n    n3 : ')
        ret.append(repr(self.n3))
        ret.append(',\n    n : ')
        ret.append(repr(self.n))
        ret.append(',\n    ng1 : ')
        ret.append(repr(self.ng1))
        ret.append(',\n    ng2 : ')
        ret.append(repr(self.ng2))
        ret.append(',\n    ng3 : ')
        ret.append(repr(self.ng3))
        ret.append(',\n    ng : ')
        ret.append(repr(self.ng))
        ret.append(',\n    gap : ')
        ret.append(repr(self.gap))
        ret.append(',\n    dvol : ')
        ret.append(repr(self.dvol))
        ret.append(',\n    nk1 : ')
        ret.append(repr(self.nk1))
        ret.append(',\n    nk2 : ')
        ret.append(repr(self.nk2))
        ret.append(',\n    nk3 : ')
        ret.append(repr(self.nk3))
        ret.append(',\n    nk : ')
        ret.append(repr(self.nk))
        ret.append(',\n    kdispl : ')
        ret.append(repr(self.kdispl))
        ret.append(',\n    global_n1 : ')
        ret.append(repr(self.global_n1))
        ret.append(',\n    global_n2 : ')
        ret.append(repr(self.global_n2))
        ret.append(',\n    global_n3 : ')
        ret.append(repr(self.global_n3))
        ret.append(',\n    global_n : ')
        ret.append(repr(self.global_n))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

grid_module = Grid_Module()

class Succeed(f90wrap.runtime.FortranModule):
    """
    Module succeed
    
    
    Defined at Succeed_module.f90 lines 1-404
    
    """
    @staticmethod
    def init_succeed_rho(n_rho, n_r, n_s, nspin, dv):
        """
        init_succeed_rho(n_rho, n_r, n_s, nspin, dv)
        
        
        Defined at Succeed_module.f90 lines 21-60
        
        Parameters
        ----------
        n_rho : int
        n_r : int
        n_s : int
        nspin : int
        dv : float
        
        """
        _arespy.f90wrap_init_succeed_rho(n_rho=n_rho, n_r=n_r, n_s=n_s, nspin=nspin, \
            dv=dv)
    
    @staticmethod
    def destory_succeed():
        """
        destory_succeed()
        
        
        Defined at Succeed_module.f90 lines 62-69
        
        
        """
        _arespy.f90wrap_destory_succeed()
    
    @staticmethod
    def store_rho(n_rho, nspin, rho):
        """
        store_rho(n_rho, nspin, rho)
        
        
        Defined at Succeed_module.f90 lines 71-84
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        rho : float array
        
        """
        _arespy.f90wrap_store_rho(n_rho=n_rho, nspin=nspin, rho=rho)
    
    @staticmethod
    def store_rho_at(n_rho, nspin, rho_in):
        """
        store_rho_at(n_rho, nspin, rho_in)
        
        
        Defined at Succeed_module.f90 lines 86-91
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        rho_in : float array
        
        """
        _arespy.f90wrap_store_rho_at(n_rho=n_rho, nspin=nspin, rho_in=rho_in)
    
    @staticmethod
    def store_r(nr, r):
        """
        store_r(nr, r)
        
        
        Defined at Succeed_module.f90 lines 93-102
        
        Parameters
        ----------
        nr : int
        r : float array
        
        """
        _arespy.f90wrap_store_r(nr=nr, r=r)
    
    @staticmethod
    def store_psi(n_rho, n_s, nspin, psi):
        """
        store_psi(n_rho, n_s, nspin, psi)
        
        
        Defined at Succeed_module.f90 lines 104-115
        
        Parameters
        ----------
        n_rho : int
        n_s : int
        nspin : int
        psi : float array
        
        """
        _arespy.f90wrap_store_psi(n_rho=n_rho, n_s=n_s, nspin=nspin, psi=psi)
    
    @staticmethod
    def get_rho(nr, r_new, nrho, nspin, rho_new):
        """
        get_rho(nr, r_new, nrho, nspin, rho_new)
        
        
        Defined at Succeed_module.f90 lines 117-207
        
        Parameters
        ----------
        nr : int
        r_new : float array
        nrho : int
        nspin : int
        rho_new : float array
        
        """
        _arespy.f90wrap_get_rho(nr=nr, r_new=r_new, nrho=nrho, nspin=nspin, \
            rho_new=rho_new)
    
    @staticmethod
    def get_psi(nrho, n_s, nspin, psi_new):
        """
        get_psi(nrho, n_s, nspin, psi_new)
        
        
        Defined at Succeed_module.f90 lines 209-230
        
        Parameters
        ----------
        nrho : int
        n_s : int
        nspin : int
        psi_new : float array
        
        """
        _arespy.f90wrap_get_psi(nrho=nrho, n_s=n_s, nspin=nspin, psi_new=psi_new)
    
    @staticmethod
    def cal_trans_phase(nr, nspin, r_new, n1, n2, n3, ng1, ng2, ng3, gvec, \
        trans_phase):
        """
        cal_trans_phase(nr, nspin, r_new, n1, n2, n3, ng1, ng2, ng3, gvec, trans_phase)
        
        
        Defined at Succeed_module.f90 lines 240-292
        
        Parameters
        ----------
        nr : int
        nspin : int
        r_new : float array
        n1 : int
        n2 : int
        n3 : int
        ng1 : int
        ng2 : int
        ng3 : int
        gvec : float array
        trans_phase : complex array
        
        """
        _arespy.f90wrap_cal_trans_phase(nr=nr, nspin=nspin, r_new=r_new, n1=n1, n2=n2, \
            n3=n3, ng1=ng1, ng2=ng2, ng3=ng3, gvec=gvec, trans_phase=trans_phase)
    
    @staticmethod
    def get_new_rho_psi(nr, r_new, nrho, n1, n2, n3, nspin, rho_new, n_s, psi_new, \
        gvec):
        """
        get_new_rho_psi(nr, r_new, nrho, n1, n2, n3, nspin, rho_new, n_s, psi_new, gvec)
        
        
        Defined at Succeed_module.f90 lines 294-353
        
        Parameters
        ----------
        nr : int
        r_new : float array
        nrho : int
        n1 : int
        n2 : int
        n3 : int
        nspin : int
        rho_new : float array
        n_s : int
        psi_new : float array
        gvec : float array
        
        """
        _arespy.f90wrap_get_new_rho_psi(nr=nr, r_new=r_new, nrho=nrho, n1=n1, n2=n2, \
            n3=n3, nspin=nspin, rho_new=rho_new, n_s=n_s, psi_new=psi_new, gvec=gvec)
    
    @staticmethod
    def store_rho_fft_trans(n_rho, nspin, rho):
        """
        store_rho_fft_trans(n_rho, nspin, rho)
        
        
        Defined at Succeed_module.f90 lines 355-363
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        rho : float array
        
        """
        _arespy.f90wrap_store_rho_fft_trans(n_rho=n_rho, nspin=nspin, rho=rho)
    
    @staticmethod
    def store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2):
        """
        store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2)
        
        
        Defined at Succeed_module.f90 lines 365-380
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        na : int
        rho_in : float array
        rho_in2 : float array
        
        """
        _arespy.f90wrap_store_rho_at_fft_trans(n_rho=n_rho, nspin=nspin, na=na, \
            rho_in=rho_in, rho_in2=rho_in2)
    
    @staticmethod
    def store_r_fft_trans(nr, r):
        """
        store_r_fft_trans(nr, r)
        
        
        Defined at Succeed_module.f90 lines 382-390
        
        Parameters
        ----------
        nr : int
        r : float array
        
        """
        _arespy.f90wrap_store_r_fft_trans(nr=nr, r=r)
    
    @staticmethod
    def store_psi_fft_trans(n_rho, n_s, nspin, psi):
        """
        store_psi_fft_trans(n_rho, n_s, nspin, psi)
        
        
        Defined at Succeed_module.f90 lines 392-404
        
        Parameters
        ----------
        n_rho : int
        n_s : int
        nspin : int
        psi : float array
        
        """
        _arespy.f90wrap_store_psi_fft_trans(n_rho=n_rho, n_s=n_s, nspin=nspin, psi=psi)
    
    @property
    def llastrho(self):
        """
        Element llastrho ftype=logical pytype=bool
        
        
        Defined at Succeed_module.f90 line 10
        
        """
        return _arespy.f90wrap_succeed__get__llastrho()
    
    @llastrho.setter
    def llastrho(self, llastrho):
        _arespy.f90wrap_succeed__set__llastrho(llastrho)
    
    @property
    def lsr(self):
        """
        Element lsr ftype=logical pytype=bool
        
        
        Defined at Succeed_module.f90 line 10
        
        """
        return _arespy.f90wrap_succeed__get__lsr()
    
    @lsr.setter
    def lsr(self, lsr):
        _arespy.f90wrap_succeed__set__lsr(lsr)
    
    @property
    def lsrho(self):
        """
        Element lsrho ftype=logical pytype=bool
        
        
        Defined at Succeed_module.f90 line 10
        
        """
        return _arespy.f90wrap_succeed__get__lsrho()
    
    @lsrho.setter
    def lsrho(self, lsrho):
        _arespy.f90wrap_succeed__set__lsrho(lsrho)
    
    @property
    def lspsi(self):
        """
        Element lspsi ftype=logical pytype=bool
        
        
        Defined at Succeed_module.f90 line 10
        
        """
        return _arespy.f90wrap_succeed__get__lspsi()
    
    @lspsi.setter
    def lspsi(self, lspsi):
        _arespy.f90wrap_succeed__set__lspsi(lspsi)
    
    @property
    def lsrho_at(self):
        """
        Element lsrho_at ftype=logical pytype=bool
        
        
        Defined at Succeed_module.f90 line 10
        
        """
        return _arespy.f90wrap_succeed__get__lsrho_at()
    
    @lsrho_at.setter
    def lsrho_at(self, lsrho_at):
        _arespy.f90wrap_succeed__set__lsrho_at(lsrho_at)
    
    @property
    def rho1(self):
        """
        Element rho1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rho1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho1 = self._arrays[array_handle]
        else:
            rho1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rho1)
            self._arrays[array_handle] = rho1
        return rho1
    
    @rho1.setter
    def rho1(self, rho1):
        self.rho1[...] = rho1
    
    @property
    def rho2(self):
        """
        Element rho2 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rho2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho2 = self._arrays[array_handle]
        else:
            rho2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rho2)
            self._arrays[array_handle] = rho2
        return rho2
    
    @rho2.setter
    def rho2(self, rho2):
        self.rho2[...] = rho2
    
    @property
    def rho3(self):
        """
        Element rho3 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rho3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho3 = self._arrays[array_handle]
        else:
            rho3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rho3)
            self._arrays[array_handle] = rho3
        return rho3
    
    @rho3.setter
    def rho3(self, rho3):
        self.rho3[...] = rho3
    
    @property
    def r1(self):
        """
        Element r1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__r1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r1 = self._arrays[array_handle]
        else:
            r1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__r1)
            self._arrays[array_handle] = r1
        return r1
    
    @r1.setter
    def r1(self, r1):
        self.r1[...] = r1
    
    @property
    def r2(self):
        """
        Element r2 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__r2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r2 = self._arrays[array_handle]
        else:
            r2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__r2)
            self._arrays[array_handle] = r2
        return r2
    
    @r2.setter
    def r2(self, r2):
        self.r2[...] = r2
    
    @property
    def r3(self):
        """
        Element r3 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__r3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r3 = self._arrays[array_handle]
        else:
            r3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__r3)
            self._arrays[array_handle] = r3
        return r3
    
    @r3.setter
    def r3(self, r3):
        self.r3[...] = r3
    
    @property
    def psi1(self):
        """
        Element psi1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__psi1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi1 = self._arrays[array_handle]
        else:
            psi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__psi1)
            self._arrays[array_handle] = psi1
        return psi1
    
    @psi1.setter
    def psi1(self, psi1):
        self.psi1[...] = psi1
    
    @property
    def psi2(self):
        """
        Element psi2 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__psi2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi2 = self._arrays[array_handle]
        else:
            psi2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__psi2)
            self._arrays[array_handle] = psi2
        return psi2
    
    @psi2.setter
    def psi2(self, psi2):
        self.psi2[...] = psi2
    
    @property
    def psi3(self):
        """
        Element psi3 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__psi3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi3 = self._arrays[array_handle]
        else:
            psi3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__psi3)
            self._arrays[array_handle] = psi3
        return psi3
    
    @psi3.setter
    def psi3(self, psi3):
        self.psi3[...] = psi3
    
    @property
    def rho_at(self):
        """
        Element rho_at ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rho_at(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho_at = self._arrays[array_handle]
        else:
            rho_at = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rho_at)
            self._arrays[array_handle] = rho_at
        return rho_at
    
    @rho_at.setter
    def rho_at(self, rho_at):
        self.rho_at[...] = rho_at
    
    @property
    def rhoi(self):
        """
        Element rhoi ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rhoi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rhoi = self._arrays[array_handle]
        else:
            rhoi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rhoi)
            self._arrays[array_handle] = rhoi
        return rhoi
    
    @rhoi.setter
    def rhoi(self, rhoi):
        self.rhoi[...] = rhoi
    
    @property
    def rho_at1(self):
        """
        Element rho_at1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rho_at1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho_at1 = self._arrays[array_handle]
        else:
            rho_at1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rho_at1)
            self._arrays[array_handle] = rho_at1
        return rho_at1
    
    @rho_at1.setter
    def rho_at1(self, rho_at1):
        self.rho_at1[...] = rho_at1
    
    @property
    def rhoi1(self):
        """
        Element rhoi1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy.f90wrap_succeed__array__rhoi1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rhoi1 = self._arrays[array_handle]
        else:
            rhoi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _arespy.f90wrap_succeed__array__rhoi1)
            self._arrays[array_handle] = rhoi1
        return rhoi1
    
    @rhoi1.setter
    def rhoi1(self, rhoi1):
        self.rhoi1[...] = rhoi1
    
    @property
    def dvol(self):
        """
        Element dvol ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 17
        
        """
        return _arespy.f90wrap_succeed__get__dvol()
    
    @dvol.setter
    def dvol(self, dvol):
        _arespy.f90wrap_succeed__set__dvol(dvol)
    
    @property
    def counter1(self):
        """
        Element counter1 ftype=integer(i4b) pytype=int
        
        
        Defined at Succeed_module.f90 line 18
        
        """
        return _arespy.f90wrap_succeed__get__counter1()
    
    @counter1.setter
    def counter1(self, counter1):
        _arespy.f90wrap_succeed__set__counter1(counter1)
    
    @property
    def counter2(self):
        """
        Element counter2 ftype=integer(i4b) pytype=int
        
        
        Defined at Succeed_module.f90 line 18
        
        """
        return _arespy.f90wrap_succeed__get__counter2()
    
    @counter2.setter
    def counter2(self, counter2):
        _arespy.f90wrap_succeed__set__counter2(counter2)
    
    @property
    def counter3(self):
        """
        Element counter3 ftype=integer(i4b) pytype=int
        
        
        Defined at Succeed_module.f90 line 18
        
        """
        return _arespy.f90wrap_succeed__get__counter3()
    
    @counter3.setter
    def counter3(self, counter3):
        _arespy.f90wrap_succeed__set__counter3(counter3)
    
    @property
    def alpha(self):
        """
        Element alpha ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 19
        
        """
        return _arespy.f90wrap_succeed__get__alpha()
    
    @alpha.setter
    def alpha(self, alpha):
        _arespy.f90wrap_succeed__set__alpha(alpha)
    
    @property
    def beta(self):
        """
        Element beta ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.f90 line 19
        
        """
        return _arespy.f90wrap_succeed__get__beta()
    
    @beta.setter
    def beta(self, beta):
        _arespy.f90wrap_succeed__set__beta(beta)
    
    def __str__(self):
        ret = ['<succeed>{\n']
        ret.append('    llastrho : ')
        ret.append(repr(self.llastrho))
        ret.append(',\n    lsr : ')
        ret.append(repr(self.lsr))
        ret.append(',\n    lsrho : ')
        ret.append(repr(self.lsrho))
        ret.append(',\n    lspsi : ')
        ret.append(repr(self.lspsi))
        ret.append(',\n    lsrho_at : ')
        ret.append(repr(self.lsrho_at))
        ret.append(',\n    rho1 : ')
        ret.append(repr(self.rho1))
        ret.append(',\n    rho2 : ')
        ret.append(repr(self.rho2))
        ret.append(',\n    rho3 : ')
        ret.append(repr(self.rho3))
        ret.append(',\n    r1 : ')
        ret.append(repr(self.r1))
        ret.append(',\n    r2 : ')
        ret.append(repr(self.r2))
        ret.append(',\n    r3 : ')
        ret.append(repr(self.r3))
        ret.append(',\n    psi1 : ')
        ret.append(repr(self.psi1))
        ret.append(',\n    psi2 : ')
        ret.append(repr(self.psi2))
        ret.append(',\n    psi3 : ')
        ret.append(repr(self.psi3))
        ret.append(',\n    rho_at : ')
        ret.append(repr(self.rho_at))
        ret.append(',\n    rhoi : ')
        ret.append(repr(self.rhoi))
        ret.append(',\n    rho_at1 : ')
        ret.append(repr(self.rho_at1))
        ret.append(',\n    rhoi1 : ')
        ret.append(repr(self.rhoi1))
        ret.append(',\n    dvol : ')
        ret.append(repr(self.dvol))
        ret.append(',\n    counter1 : ')
        ret.append(repr(self.counter1))
        ret.append(',\n    counter2 : ')
        ret.append(repr(self.counter2))
        ret.append(',\n    counter3 : ')
        ret.append(repr(self.counter3))
        ret.append(',\n    alpha : ')
        ret.append(repr(self.alpha))
        ret.append(',\n    beta : ')
        ret.append(repr(self.beta))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

succeed = Succeed()

class Write_Module(f90wrap.runtime.FortranModule):
    """
    Module write_module
    
    
    Defined at Write_module.f90 lines 1-248
    
    """
    @staticmethod
    def write_poscar(filename, lat_mat, eleid, pos, atom_symbol=None, fixpos=None):
        """
        write_poscar(filename, lat_mat, eleid, pos[, atom_symbol, fixpos])
        
        
        Defined at Write_module.f90 lines 20-64
        
        Parameters
        ----------
        filename : str
        lat_mat : float array
        eleid : int array
        pos : float array
        atom_symbol : str array
        fixpos : bool array
        
        """
        _arespy.f90wrap_write_poscar(filename=filename, lat_mat=lat_mat, eleid=eleid, \
            pos=pos, atom_symbol=atom_symbol, fixpos=fixpos)
    
    @staticmethod
    def write_cif(filename, lat_para, nati, pos, atom_symbol=None):
        """
        write_cif(filename, lat_para, nati, pos[, atom_symbol])
        
        
        Defined at Write_module.f90 lines 68-133
        
        Parameters
        ----------
        filename : str
        lat_para : float array
        nati : int array
        pos : float array
        atom_symbol : str array
        
        """
        _arespy.f90wrap_write_cif(filename=filename, lat_para=lat_para, nati=nati, \
            pos=pos, atom_symbol=atom_symbol)
    
    @staticmethod
    def write3dat(filename, rho):
        """
        write3dat(filename, rho)
        
        
        Defined at Write_module.f90 lines 138-167
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        -----------------------------------------------------------------------
        """
        _arespy.f90wrap_write3dat(filename=filename, rho=rho)
    
    @staticmethod
    def write3d(filename, rho):
        """
        write3d(filename, rho)
        
        
        Defined at Write_module.f90 lines 171-195
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        -----------------------------------------------------------------------
        """
        _arespy.f90wrap_write3d(filename=filename, rho=rho)
    
    @staticmethod
    def write3dpencil(filename, rho):
        """
        write3dpencil(filename, rho)
        
        
        Defined at Write_module.f90 lines 200-233
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        -----------------------------------------------------------------------
        """
        _arespy.f90wrap_write3dpencil(filename=filename, rho=rho)
    
    @staticmethod
    def writedensity(filename, rho):
        """
        writedensity(filename, rho)
        
        
        Defined at Write_module.f90 lines 238-245
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        """
        _arespy.f90wrap_writedensity(filename=filename, rho=rho)
    
    _dt_array_initialisers = []
    

write_module = Write_Module()

class Xc_Module(f90wrap.runtime.FortranModule):
    """
    Module xc_module
    
    
    Defined at XC_functional.f90 lines 5-486
    
    """
    @staticmethod
    def xc_functional(nps, rhos, vxc, exc=None):
        """
        xc_functional(nps, rhos, vxc[, exc])
        
        
        Defined at XC_functional.f90 lines 19-44
        
        Parameters
        ----------
        nps : int
        rhos : float array
        vxc : float array
        exc : float
        
        """
        _arespy.f90wrap_xc_functional(nps=nps, rhos=rhos, vxc=vxc, exc=exc)
    
    @staticmethod
    def libxc_lda_set(nps, rhos, vxcs, exc=None):
        """
        libxc_lda_set(nps, rhos, vxcs[, exc])
        
        
        Defined at XC_functional.f90 lines 49-120
        
        Parameters
        ----------
        nps : int
        rhos : float array
        vxcs : float array
        exc : float
        
        """
        _arespy.f90wrap_libxc_lda_set(nps=nps, rhos=rhos, vxcs=vxcs, exc=exc)
    
    @staticmethod
    def libxc_lda_x(nps, rhos, vxs=None):
        """
        ex = libxc_lda_x(nps, rhos[, vxs])
        
        
        Defined at XC_functional.f90 lines 125-178
        
        Parameters
        ----------
        nps : int
        rhos : float array
        vxs : float array
        
        Returns
        -------
        ex : float
        
        """
        ex = _arespy.f90wrap_libxc_lda_x(nps=nps, rhos=rhos, vxs=vxs)
        return ex
    
    @staticmethod
    def libxc_vwn1rpa_c(nps, rhos, vcs=None):
        """
        ec = libxc_vwn1rpa_c(nps, rhos[, vcs])
        
        
        Defined at XC_functional.f90 lines 183-236
        
        Parameters
        ----------
        nps : int
        rhos : float array
        vcs : float array
        
        Returns
        -------
        ec : float
        
        """
        ec = _arespy.f90wrap_libxc_vwn1rpa_c(nps=nps, rhos=rhos, vcs=vcs)
        return ec
    
    @staticmethod
    def libxc_gga_set(nps, rhos, vxcs, exc=None):
        """
        libxc_gga_set(nps, rhos, vxcs[, exc])
        
        
        Defined at XC_functional.f90 lines 241-357
        
        Parameters
        ----------
        nps : int
        rhos : float array
        vxcs : float array
        exc : float
        
        """
        _arespy.f90wrap_libxc_gga_set(nps=nps, rhos=rhos, vxcs=vxcs, exc=exc)
    
    @staticmethod
    def libxc_b88_x(nps, rhos, sigma, vxs=None):
        """
        ex = libxc_b88_x(nps, rhos, sigma[, vxs])
        
        
        Defined at XC_functional.f90 lines 362-419
        
        Parameters
        ----------
        nps : int
        rhos : float array
        sigma : float array
        vxs : float array
        
        Returns
        -------
        ex : float
        
        """
        ex = _arespy.f90wrap_libxc_b88_x(nps=nps, rhos=rhos, sigma=sigma, vxs=vxs)
        return ex
    
    @staticmethod
    def libxc_lyp_c(nps, rhos, sigma, vcs=None):
        """
        ec = libxc_lyp_c(nps, rhos, sigma[, vcs])
        
        
        Defined at XC_functional.f90 lines 424-482
        
        Parameters
        ----------
        nps : int
        rhos : float array
        sigma : float array
        vcs : float array
        
        Returns
        -------
        ec : float
        
        """
        ec = _arespy.f90wrap_libxc_lyp_c(nps=nps, rhos=rhos, sigma=sigma, vcs=vcs)
        return ec
    
    _dt_array_initialisers = []
    

xc_module = Xc_Module()

