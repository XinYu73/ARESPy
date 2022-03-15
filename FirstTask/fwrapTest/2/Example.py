from __future__ import print_function, absolute_import, division
import _Example
import f90wrap.runtime
import logging

class Mcyldnad(f90wrap.runtime.FortranModule):
    """
    Module mcyldnad
    
    
    Defined at cyldnad.fpp lines 5-14
    
    """
    @staticmethod
    def cyldnad(self, height):
        """
        vol = cyldnad(self, height)
        
        
        Defined at cyldnad.fpp lines 8-13
        
        Parameters
        ----------
        radius : Dual_Num
        height : Dual_Num
        
        Returns
        -------
        vol : Dual_Num
        
        """
        vol = _Example.f90wrap_cyldnad(radius=self._handle, height=height._handle)
        vol = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(vol, \
            alloc=True)
        return vol
    
    _dt_array_initialisers = []
    

mcyldnad = Mcyldnad()

class Dual_Num_Auto_Diff(f90wrap.runtime.FortranModule):
    """
    Module dual_num_auto_diff
    
    
    Defined at DNAD.fpp lines 112-1609
    
    """
    @f90wrap.runtime.register_class("Example.DUAL_NUM")
    class DUAL_NUM(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=dual_num)
        
        
        Defined at DNAD.fpp lines 119-124
        
        """
        def __init__(self, handle=None):
            """
            self = Dual_Num()
            
            
            Defined at DNAD.fpp lines 119-124
            
            
            Returns
            -------
            this : Dual_Num
            	Object to be constructed
            
            
            Automatically generated constructor for dual_num
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _Example.f90wrap_dual_num_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Dual_Num
            
            
            Defined at DNAD.fpp lines 119-124
            
            Parameters
            ----------
            this : Dual_Num
            	Object to be destructed
            
            
            Automatically generated destructor for dual_num
            """
            if self._alloc:
                _Example.f90wrap_dual_num_finalise(this=self._handle)
        
        @property
        def x_ad_(self):
            """
            Element x_ad_ ftype=real(dbl_ad) pytype=float
            
            
            Defined at DNAD.fpp line 123
            
            """
            return _Example.f90wrap_dual_num__get__x_ad_(self._handle)
        
        @x_ad_.setter
        def x_ad_(self, x_ad_):
            _Example.f90wrap_dual_num__set__x_ad_(self._handle, x_ad_)
        
        @property
        def xp_ad_(self):
            """
            Element xp_ad_ ftype=real(dbl_ad) pytype=float
            
            
            Defined at DNAD.fpp line 124
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _Example.f90wrap_dual_num__array__xp_ad_(self._handle)
            if array_handle in self._arrays:
                xp_ad_ = self._arrays[array_handle]
            else:
                xp_ad_ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _Example.f90wrap_dual_num__array__xp_ad_)
                self._arrays[array_handle] = xp_ad_
            return xp_ad_
        
        @xp_ad_.setter
        def xp_ad_(self, xp_ad_):
            self.xp_ad_[...] = xp_ad_
        
        def __str__(self):
            ret = ['<dual_num>{\n']
            ret.append('    x_ad_ : ')
            ret.append(repr(self.x_ad_))
            ret.append(',\n    xp_ad_ : ')
            ret.append(repr(self.xp_ad_))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def abs_d(self):
        """
        res = abs_d(self)
        
        
        Defined at DNAD.fpp lines 1223-1235
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_abs_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def acos_d(self):
        """
        res = acos_d(self)
        
        
        Defined at DNAD.fpp lines 1241-1251
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_acos_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def asin_d(self):
        """
        res = asin_d(self)
        
        
        Defined at DNAD.fpp lines 1257-1267
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_asin_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def cos_d(self):
        """
        res = cos_d(self)
        
        
        Defined at DNAD.fpp lines 1273-1279
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_cos_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def exp_d(self):
        """
        res = exp_d(self)
        
        
        Defined at DNAD.fpp lines 1298-1304
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_exp_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def int_d(self):
        """
        res = int_d(self)
        
        
        Defined at DNAD.fpp lines 1310-1315
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : int
        
        """
        res = _Example.f90wrap_int_d(u=self._handle)
        return res
    
    @staticmethod
    def log_d(self):
        """
        res = log_d(self)
        
        
        Defined at DNAD.fpp lines 1323-1329
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_log_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def log10_d(self):
        """
        res = log10_d(self)
        
        
        Defined at DNAD.fpp lines 1337-1343
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_log10_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def nint_d(self):
        """
        res = nint_d(self)
        
        
        Defined at DNAD.fpp lines 1533-1536
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : int
        
        """
        res = _Example.f90wrap_nint_d(u=self._handle)
        return res
    
    @staticmethod
    def sin_d(self):
        """
        res = sin_d(self)
        
        
        Defined at DNAD.fpp lines 1569-1575
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_sin_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def sqrt_d(self):
        """
        res = sqrt_d(self)
        
        
        Defined at DNAD.fpp lines 1581-1592
        
        Parameters
        ----------
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_sqrt_d(u=self._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _max_dd(self, val2, val3=None, val4=None, val5=None):
        """
        res = _max_dd(self, val2[, val3, val4, val5])
        
        
        Defined at DNAD.fpp lines 1391-1408
        
        Parameters
        ----------
        val1 : Dual_Num
        val2 : Dual_Num
        val3 : Dual_Num
        val4 : Dual_Num
        val5 : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_max_dd(val1=self._handle, val2=val2._handle, val3=None if \
            val3 is None else val3._handle, val4=None if val4 is None else val4._handle, \
            val5=None if val5 is None else val5._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _max_di(self, n):
        """
        res = _max_di(self, n)
        
        
        Defined at DNAD.fpp lines 1413-1421
        
        Parameters
        ----------
        u : Dual_Num
        n : int
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_max_di(u=self._handle, n=n)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _max_dr(self, n):
        """
        res = _max_dr(self, n)
        
        
        Defined at DNAD.fpp lines 1426-1434
        
        Parameters
        ----------
        u : Dual_Num
        n : float
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_max_dr(u=self._handle, n=n)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _max_ds(self, n):
        """
        res = _max_ds(self, n)
        
        
        Defined at DNAD.fpp lines 1439-1447
        
        Parameters
        ----------
        u : Dual_Num
        n : float
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_max_ds(u=self._handle, n=n)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _max_rd(r, u):
        """
        res = _max_rd(r, u)
        
        
        Defined at DNAD.fpp lines 1455-1463
        
        Parameters
        ----------
        r : float
        u : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_max_rd(r=r, u=u._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def max(*args, **kwargs):
        """
        max(*args, **kwargs)
        
        
        Defined at DNAD.fpp lines 334-339
        
        Overloaded interface containing the following procedures:
          _max_dd
          _max_di
          _max_dr
          _max_ds
          _max_rd
        
        """
        for proc in [Dual_Num_Auto_Diff._max_dd, Dual_Num_Auto_Diff._max_di, \
            Dual_Num_Auto_Diff._max_dr, Dual_Num_Auto_Diff._max_ds, \
            Dual_Num_Auto_Diff._max_rd]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _min_dd(self, val2, val3=None, val4=None):
        """
        res = _min_dd(self, val2[, val3, val4])
        
        
        Defined at DNAD.fpp lines 1478-1492
        
        Parameters
        ----------
        val1 : Dual_Num
        val2 : Dual_Num
        val3 : Dual_Num
        val4 : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_min_dd(val1=self._handle, val2=val2._handle, val3=None if \
            val3 is None else val3._handle, val4=None if val4 is None else val4._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _min_dr(self, n):
        """
        res = _min_dr(self, n)
        
        
        Defined at DNAD.fpp lines 1497-1505
        
        Parameters
        ----------
        u : Dual_Num
        n : float
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_min_dr(u=self._handle, n=n)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _min_ds(self, n):
        """
        res = _min_ds(self, n)
        
        
        Defined at DNAD.fpp lines 1510-1518
        
        Parameters
        ----------
        u : Dual_Num
        n : float
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_min_ds(u=self._handle, n=n)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def min(*args, **kwargs):
        """
        min(*args, **kwargs)
        
        
        Defined at DNAD.fpp lines 350-353
        
        Overloaded interface containing the following procedures:
          _min_dd
          _min_dr
          _min_ds
        
        """
        for proc in [Dual_Num_Auto_Diff._min_dd, Dual_Num_Auto_Diff._min_dr, \
            Dual_Num_Auto_Diff._min_ds]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _sign_dd(self, val2):
        """
        res = _sign_dd(self, val2)
        
        
        Defined at DNAD.fpp lines 1542-1549
        
        Parameters
        ----------
        val1 : Dual_Num
        val2 : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_sign_dd(val1=self._handle, val2=val2._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def _sign_rd(val1, val2):
        """
        res = _sign_rd(val1, val2)
        
        
        Defined at DNAD.fpp lines 1555-1563
        
        Parameters
        ----------
        val1 : float
        val2 : Dual_Num
        
        Returns
        -------
        res : Dual_Num
        
        """
        res = _Example.f90wrap_sign_rd(val1=val1, val2=val2._handle)
        res = f90wrap.runtime.lookup_class("Example.DUAL_NUM").from_handle(res, \
            alloc=True)
        return res
    
    @staticmethod
    def sign(*args, **kwargs):
        """
        sign(*args, **kwargs)
        
        
        Defined at DNAD.fpp lines 368-370
        
        Overloaded interface containing the following procedures:
          _sign_dd
          _sign_rd
        
        """
        for proc in [Dual_Num_Auto_Diff._sign_dd, Dual_Num_Auto_Diff._sign_rd]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def ndv_ad(self):
        """
        Element ndv_ad ftype=integer(2) pytype=int
        
        
        Defined at DNAD.fpp line 114
        
        """
        return _Example.f90wrap_dual_num_auto_diff__get__ndv_ad()
    
    def __str__(self):
        ret = ['<dual_num_auto_diff>{\n']
        ret.append('    ndv_ad : ')
        ret.append(repr(self.ndv_ad))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

dual_num_auto_diff = Dual_Num_Auto_Diff()

