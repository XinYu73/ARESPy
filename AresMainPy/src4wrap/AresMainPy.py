from __future__ import print_function, absolute_import, division
import _AresMainPy
import f90wrap.runtime
import logging

class Constants(f90wrap.runtime.FortranModule):
    """
    Module constants
    
    
    Defined at Constants.fpp lines 5-97
    
    """
    @property
    def i4b(self):
        """
        Element i4b ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 25
        
        """
        return _AresMainPy.f90wrap_constants__get__i4b()
    
    @property
    def i2b(self):
        """
        Element i2b ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 26
        
        """
        return _AresMainPy.f90wrap_constants__get__i2b()
    
    @property
    def i1b(self):
        """
        Element i1b ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 27
        
        """
        return _AresMainPy.f90wrap_constants__get__i1b()
    
    @property
    def sp(self):
        """
        Element sp ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 28
        
        """
        return _AresMainPy.f90wrap_constants__get__sp()
    
    @property
    def dp(self):
        """
        Element dp ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 29
        
        """
        return _AresMainPy.f90wrap_constants__get__dp()
    
    @property
    def scp(self):
        """
        Element scp ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 30
        
        """
        return _AresMainPy.f90wrap_constants__get__scp()
    
    @property
    def dcp(self):
        """
        Element dcp ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 31
        
        """
        return _AresMainPy.f90wrap_constants__get__dcp()
    
    @property
    def lgt(self):
        """
        Element lgt ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 33
        
        """
        return _AresMainPy.f90wrap_constants__get__lgt()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 36
        
        """
        return _AresMainPy.f90wrap_constants__get__pi()
    
    @property
    def pio2(self):
        """
        Element pio2 ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 37
        
        """
        return _AresMainPy.f90wrap_constants__get__pio2()
    
    @property
    def twopi(self):
        """
        Element twopi ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 38
        
        """
        return _AresMainPy.f90wrap_constants__get__twopi()
    
    @property
    def imag(self):
        """
        Element imag ftype=complex(dp) pytype=complex
        
        
        Defined at Constants.fpp line 39
        
        """
        return _AresMainPy.f90wrap_constants__get__imag()
    
    @property
    def inputunit(self):
        """
        Element inputunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.fpp line 40
        
        """
        return _AresMainPy.f90wrap_constants__get__inputunit()
    
    @property
    def errunit(self):
        """
        Element errunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.fpp line 41
        
        """
        return _AresMainPy.f90wrap_constants__get__errunit()
    
    @property
    def outputunit(self):
        """
        Element outputunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.fpp line 42
        
        """
        return _AresMainPy.f90wrap_constants__get__outputunit()
    
    @property
    def mdposunit(self):
        """
        Element mdposunit ftype=integer(i4b) pytype=int
        
        
        Defined at Constants.fpp line 43
        
        """
        return _AresMainPy.f90wrap_constants__get__mdposunit()
    
    @property
    def rydberg(self):
        """
        Element rydberg ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 44
        
        """
        return _AresMainPy.f90wrap_constants__get__rydberg()
    
    @property
    def bohr2ang(self):
        """
        Element bohr2ang ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 45
        
        """
        return _AresMainPy.f90wrap_constants__get__bohr2ang()
    
    @property
    def ang2bohr(self):
        """
        Element ang2bohr ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 46
        
        """
        return _AresMainPy.f90wrap_constants__get__ang2bohr()
    
    @property
    def hart2ev(self):
        """
        Element hart2ev ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 47
        
        """
        return _AresMainPy.f90wrap_constants__get__hart2ev()
    
    @property
    def ev2hart(self):
        """
        Element ev2hart ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 48
        
        """
        return _AresMainPy.f90wrap_constants__get__ev2hart()
    
    @property
    def scf_tol(self):
        """
        Element scf_tol ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 49
        
        """
        return _AresMainPy.f90wrap_constants__get__scf_tol()
    
    @property
    def au2ev_force(self):
        """
        Element au2ev_force ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 50
        
        """
        return _AresMainPy.f90wrap_constants__get__au2ev_force()
    
    @property
    def au2gpa(self):
        """
        Element au2gpa ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 51
        
        """
        return _AresMainPy.f90wrap_constants__get__au2gpa()
    
    @property
    def golden(self):
        """
        Element golden ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 52
        
        """
        return _AresMainPy.f90wrap_constants__get__golden()
    
    @property
    def vlight(self):
        """
        Element vlight ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 53
        
        """
        return _AresMainPy.f90wrap_constants__get__vlight()
    
    @property
    def const_me_au(self):
        """
        Element const_me_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 57
        
        """
        return _AresMainPy.f90wrap_constants__get__const_me_au()
    
    @property
    def const_e_au(self):
        """
        Element const_e_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 59
        
        """
        return _AresMainPy.f90wrap_constants__get__const_e_au()
    
    @property
    def const_eh_au(self):
        """
        Element const_eh_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 61
        
        """
        return _AresMainPy.f90wrap_constants__get__const_eh_au()
    
    @property
    def const_len_au(self):
        """
        Element const_len_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 63
        
        """
        return _AresMainPy.f90wrap_constants__get__const_len_au()
    
    @property
    def const_hbar_au(self):
        """
        Element const_hbar_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 65
        
        """
        return _AresMainPy.f90wrap_constants__get__const_hbar_au()
    
    @property
    def const_ep_au(self):
        """
        Element const_ep_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 67
        
        """
        return _AresMainPy.f90wrap_constants__get__const_ep_au()
    
    @property
    def const_f_au(self):
        """
        Element const_f_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 69
        
        """
        return _AresMainPy.f90wrap_constants__get__const_f_au()
    
    @property
    def const_mt_au(self):
        """
        Element const_mt_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 71
        
        """
        return _AresMainPy.f90wrap_constants__get__const_mt_au()
    
    @property
    def const_time_au(self):
        """
        Element const_time_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 73
        
        """
        return _AresMainPy.f90wrap_constants__get__const_time_au()
    
    @property
    def const_i_au(self):
        """
        Element const_i_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 75
        
        """
        return _AresMainPy.f90wrap_constants__get__const_i_au()
    
    @property
    def const_temp_au(self):
        """
        Element const_temp_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 77
        
        """
        return _AresMainPy.f90wrap_constants__get__const_temp_au()
    
    @property
    def const_p_au(self):
        """
        Element const_p_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 79
        
        """
        return _AresMainPy.f90wrap_constants__get__const_p_au()
    
    @property
    def const_v_au(self):
        """
        Element const_v_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 81
        
        """
        return _AresMainPy.f90wrap_constants__get__const_v_au()
    
    @property
    def const_ke_au(self):
        """
        Element const_ke_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 83
        
        """
        return _AresMainPy.f90wrap_constants__get__const_ke_au()
    
    @property
    def const_mu_au(self):
        """
        Element const_mu_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 85
        
        """
        return _AresMainPy.f90wrap_constants__get__const_mu_au()
    
    @property
    def const_ma_au(self):
        """
        Element const_ma_au ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 87
        
        """
        return _AresMainPy.f90wrap_constants__get__const_ma_au()
    
    @property
    def const_kb_si(self):
        """
        Element const_kb_si ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 89
        
        """
        return _AresMainPy.f90wrap_constants__get__const_kb_si()
    
    @property
    def angs(self):
        """
        Element angs ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 92
        
        """
        return _AresMainPy.f90wrap_constants__get__angs()
    
    @property
    def force2ev(self):
        """
        Element force2ev ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 93
        
        """
        return _AresMainPy.f90wrap_constants__get__force2ev()
    
    @property
    def clen(self):
        """
        Element clen ftype=integer pytype=int
        
        
        Defined at Constants.fpp line 95
        
        """
        return _AresMainPy.f90wrap_constants__get__clen()
    
    @property
    def xtiny(self):
        """
        Element xtiny ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 96
        
        """
        return _AresMainPy.f90wrap_constants__get__xtiny()
    
    @property
    def autogpa(self):
        """
        Element autogpa ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 97
        
        """
        return _AresMainPy.f90wrap_constants__get__autogpa()
    
    @property
    def hartree2ev(self):
        """
        Element hartree2ev ftype=real(dp) pytype=float
        
        
        Defined at Constants.fpp line 98
        
        """
        return _AresMainPy.f90wrap_constants__get__hartree2ev()
    
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

class M_Time_Evaluate(f90wrap.runtime.FortranModule):
    """
    Module m_time_evaluate
    
    
    Defined at m_time_evaluate.fpp lines 5-202
    
    """
    @f90wrap.runtime.register_class("AresMainPy.time_record")
    class time_record(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=time_record)
        
        
        Defined at m_time_evaluate.fpp lines 15-18
        
        """
        def __init__(self, handle=None):
            """
            self = Time_Record()
            
            
            Defined at m_time_evaluate.fpp lines 15-18
            
            
            Returns
            -------
            this : Time_Record
            	Object to be constructed
            
            
            Automatically generated constructor for time_record
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_time_record_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Time_Record
            
            
            Defined at m_time_evaluate.fpp lines 15-18
            
            Parameters
            ----------
            this : Time_Record
            	Object to be destructed
            
            
            Automatically generated destructor for time_record
            """
            if self._alloc:
                _AresMainPy.f90wrap_time_record_finalise(this=self._handle)
        
        @property
        def str(self):
            """
            Element str ftype=character(len=100) pytype=str
            
            
            Defined at m_time_evaluate.fpp line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_time_record__array__str(self._handle)
            if array_handle in self._arrays:
                str = self._arrays[array_handle]
            else:
                str = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_time_record__array__str)
                self._arrays[array_handle] = str
            return str
        
        @str.setter
        def str(self, str):
            self.str[...] = str
        
        @property
        def t1(self):
            """
            Element t1 ftype=integer(i4b) pytype=int
            
            
            Defined at m_time_evaluate.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_time_record__array__t1(self._handle)
            if array_handle in self._arrays:
                t1 = self._arrays[array_handle]
            else:
                t1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_time_record__array__t1)
                self._arrays[array_handle] = t1
            return t1
        
        @t1.setter
        def t1(self, t1):
            self.t1[...] = t1
        
        @property
        def t2(self):
            """
            Element t2 ftype=integer(i4b) pytype=int
            
            
            Defined at m_time_evaluate.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_time_record__array__t2(self._handle)
            if array_handle in self._arrays:
                t2 = self._arrays[array_handle]
            else:
                t2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_time_record__array__t2)
                self._arrays[array_handle] = t2
            return t2
        
        @t2.setter
        def t2(self, t2):
            self.t2[...] = t2
        
        @property
        def length(self):
            """
            Element length ftype=integer(i4b) pytype=int
            
            
            Defined at m_time_evaluate.fpp line 18
            
            """
            return _AresMainPy.f90wrap_time_record__get__length(self._handle)
        
        @length.setter
        def length(self, length):
            _AresMainPy.f90wrap_time_record__set__length(self._handle, length)
        
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
        
        
        Defined at m_time_evaluate.fpp lines 25-50
        
        Parameters
        ----------
        str : str
        t : int
        
        """
        _AresMainPy.f90wrap_record(str=str, t=t)
    
    @staticmethod
    def reset_data(myflag, str, t):
        """
        reset_data(myflag, str, t)
        
        
        Defined at m_time_evaluate.fpp lines 52-81
        
        Parameters
        ----------
        myflag : bool
        str : str
        t : int
        
        """
        _AresMainPy.f90wrap_reset_data(myflag=myflag, str=str, t=t)
    
    @staticmethod
    def init_recrod():
        """
        init_recrod()
        
        
        Defined at m_time_evaluate.fpp lines 83-88
        
        
        """
        _AresMainPy.f90wrap_init_recrod()
    
    @staticmethod
    def time_start(str, flag):
        """
        time_start(str, flag)
        
        
        Defined at m_time_evaluate.fpp lines 90-100
        
        Parameters
        ----------
        str : str
        flag : bool
        
        """
        _AresMainPy.f90wrap_time_start(str=str, flag=flag)
    
    @staticmethod
    def time_end(str, flag):
        """
        time_end(str, flag)
        
        
        Defined at m_time_evaluate.fpp lines 102-112
        
        Parameters
        ----------
        str : str
        flag : bool
        
        """
        _AresMainPy.f90wrap_time_end(str=str, flag=flag)
    
    @staticmethod
    def time_output(str, flag):
        """
        time_output(str, flag)
        
        
        Defined at m_time_evaluate.fpp lines 114-129
        
        Parameters
        ----------
        str : str
        flag : bool
        
        """
        _AresMainPy.f90wrap_time_output(str=str, flag=flag)
    
    @staticmethod
    def memory_sum(label, memory_size):
        """
        memory_sum(label, memory_size)
        
        
        Defined at m_time_evaluate.fpp lines 131-159
        
        Parameters
        ----------
        label : str
        memory_size : float
        
        """
        _AresMainPy.f90wrap_memory_sum(label=label, memory_size=memory_size)
    
    @staticmethod
    def memory_free(label, memory_size):
        """
        memory_free(label, memory_size)
        
        
        Defined at m_time_evaluate.fpp lines 161-183
        
        Parameters
        ----------
        label : str
        memory_size : float
        
        """
        _AresMainPy.f90wrap_memory_free(label=label, memory_size=memory_size)
    
    @property
    def total_memory_consum(self):
        """
        Element total_memory_consum ftype=real(dp) pytype=float
        
        
        Defined at m_time_evaluate.fpp line 21
        
        """
        return _AresMainPy.f90wrap_m_time_evaluate__get__total_memory_consum()
    
    @total_memory_consum.setter
    def total_memory_consum(self, total_memory_consum):
        _AresMainPy.f90wrap_m_time_evaluate__set__total_memory_consum(total_memory_consum)
    
    @property
    def filename(self):
        """
        Element filename ftype=character(len=20) pytype=str
        
        
        Defined at m_time_evaluate.fpp line 22
        
        """
        return _AresMainPy.f90wrap_m_time_evaluate__get__filename()
    
    @filename.setter
    def filename(self, filename):
        _AresMainPy.f90wrap_m_time_evaluate__set__filename(filename)
    
    @property
    def file_unit(self):
        """
        Element file_unit ftype=integer(i4b) pytype=int
        
        
        Defined at m_time_evaluate.fpp line 23
        
        """
        return _AresMainPy.f90wrap_m_time_evaluate__get__file_unit()
    
    @file_unit.setter
    def file_unit(self, file_unit):
        _AresMainPy.f90wrap_m_time_evaluate__set__file_unit(file_unit)
    
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

class Parameters(f90wrap.runtime.FortranModule):
    """
    Module parameters
    
    
    Defined at Parameters.fpp lines 5-354
    
    """
    @f90wrap.runtime.register_class("AresMainPy.IO_INFO_TYPE")
    class IO_INFO_TYPE(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=io_info_type)
        
        
        Defined at Parameters.fpp lines 125-133
        
        """
        def __init__(self, handle=None):
            """
            self = Io_Info_Type()
            
            
            Defined at Parameters.fpp lines 125-133
            
            
            Returns
            -------
            this : Io_Info_Type
            	Object to be constructed
            
            
            Automatically generated constructor for io_info_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_io_info_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Io_Info_Type
            
            
            Defined at Parameters.fpp lines 125-133
            
            Parameters
            ----------
            this : Io_Info_Type
            	Object to be destructed
            
            
            Automatically generated destructor for io_info_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_io_info_type_finalise(this=self._handle)
        
        @property
        def iu0(self):
            """
            Element iu0 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 126
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu0(self._handle)
        
        @iu0.setter
        def iu0(self, iu0):
            _AresMainPy.f90wrap_io_info_type__set__iu0(self._handle, iu0)
        
        @property
        def iu1(self):
            """
            Element iu1 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 127
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu1(self._handle)
        
        @iu1.setter
        def iu1(self, iu1):
            _AresMainPy.f90wrap_io_info_type__set__iu1(self._handle, iu1)
        
        @property
        def iu2(self):
            """
            Element iu2 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 128
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu2(self._handle)
        
        @iu2.setter
        def iu2(self, iu2):
            _AresMainPy.f90wrap_io_info_type__set__iu2(self._handle, iu2)
        
        @property
        def iu3(self):
            """
            Element iu3 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 129
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu3(self._handle)
        
        @iu3.setter
        def iu3(self, iu3):
            _AresMainPy.f90wrap_io_info_type__set__iu3(self._handle, iu3)
        
        @property
        def iu4(self):
            """
            Element iu4 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 130
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu4(self._handle)
        
        @iu4.setter
        def iu4(self, iu4):
            _AresMainPy.f90wrap_io_info_type__set__iu4(self._handle, iu4)
        
        @property
        def iu5(self):
            """
            Element iu5 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 131
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu5(self._handle)
        
        @iu5.setter
        def iu5(self, iu5):
            _AresMainPy.f90wrap_io_info_type__set__iu5(self._handle, iu5)
        
        @property
        def iu6(self):
            """
            Element iu6 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 132
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu6(self._handle)
        
        @iu6.setter
        def iu6(self, iu6):
            _AresMainPy.f90wrap_io_info_type__set__iu6(self._handle, iu6)
        
        @property
        def iu8(self):
            """
            Element iu8 ftype=integer(i4b) pytype=int
            
            
            Defined at Parameters.fpp line 133
            
            """
            return _AresMainPy.f90wrap_io_info_type__get__iu8(self._handle)
        
        @iu8.setter
        def iu8(self, iu8):
            _AresMainPy.f90wrap_io_info_type__set__iu8(self._handle, iu8)
        
        def __str__(self):
            ret = ['<io_info_type>{\n']
            ret.append('    iu0 : ')
            ret.append(repr(self.iu0))
            ret.append(',\n    iu1 : ')
            ret.append(repr(self.iu1))
            ret.append(',\n    iu2 : ')
            ret.append(repr(self.iu2))
            ret.append(',\n    iu3 : ')
            ret.append(repr(self.iu3))
            ret.append(',\n    iu4 : ')
            ret.append(repr(self.iu4))
            ret.append(',\n    iu5 : ')
            ret.append(repr(self.iu5))
            ret.append(',\n    iu6 : ')
            ret.append(repr(self.iu6))
            ret.append(',\n    iu8 : ')
            ret.append(repr(self.iu8))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def generate_infile(infile_name):
        """
        generate_infile(infile_name)
        
        
        Defined at Parameters.fpp lines 216-354
        
        Parameters
        ----------
        infile_name : str
        
        """
        _AresMainPy.f90wrap_generate_infile(infile_name=infile_name)
    
    @property
    def ikedf(self):
        """
        Element ikedf ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 10
        
        """
        return _AresMainPy.f90wrap_parameters__get__ikedf()
    
    @ikedf.setter
    def ikedf(self, ikedf):
        _AresMainPy.f90wrap_parameters__set__ikedf(ikedf)
    
    @property
    def ixcdf(self):
        """
        Element ixcdf ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 11
        
        """
        return _AresMainPy.f90wrap_parameters__get__ixcdf()
    
    @ixcdf.setter
    def ixcdf(self, ixcdf):
        _AresMainPy.f90wrap_parameters__set__ixcdf(ixcdf)
    
    @property
    def finite_order(self):
        """
        Element finite_order ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 12
        
        """
        return _AresMainPy.f90wrap_parameters__get__finite_order()
    
    @finite_order.setter
    def finite_order(self, finite_order):
        _AresMainPy.f90wrap_parameters__set__finite_order(finite_order)
    
    @property
    def ntype(self):
        """
        Element ntype ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 13
        
        """
        return _AresMainPy.f90wrap_parameters__get__ntype()
    
    @ntype.setter
    def ntype(self, ntype):
        _AresMainPy.f90wrap_parameters__set__ntype(ntype)
    
    @property
    def naddstates(self):
        """
        Element naddstates ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__naddstates()
    
    @naddstates.setter
    def naddstates(self, naddstates):
        _AresMainPy.f90wrap_parameters__set__naddstates(naddstates)
    
    @property
    def istart(self):
        """
        Element istart ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__istart()
    
    @istart.setter
    def istart(self, istart):
        _AresMainPy.f90wrap_parameters__set__istart(istart)
    
    @property
    def idiag(self):
        """
        Element idiag ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__idiag()
    
    @idiag.setter
    def idiag(self, idiag):
        _AresMainPy.f90wrap_parameters__set__idiag(idiag)
    
    @property
    def chem(self):
        """
        Element chem ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__chem()
    
    @chem.setter
    def chem(self, chem):
        _AresMainPy.f90wrap_parameters__set__chem(chem)
    
    @property
    def chem0(self):
        """
        Element chem0 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__chem0()
    
    @chem0.setter
    def chem0(self, chem0):
        _AresMainPy.f90wrap_parameters__set__chem0(chem0)
    
    @property
    def nstates(self):
        """
        Element nstates ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__nstates()
    
    @nstates.setter
    def nstates(self, nstates):
        _AresMainPy.f90wrap_parameters__set__nstates(nstates)
    
    @property
    def nstates_global(self):
        """
        Element nstates_global ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__nstates_global()
    
    @nstates_global.setter
    def nstates_global(self, nstates_global):
        _AresMainPy.f90wrap_parameters__set__nstates_global(nstates_global)
    
    @property
    def nprr(self):
        """
        Element nprr ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__nprr()
    
    @nprr.setter
    def nprr(self, nprr):
        _AresMainPy.f90wrap_parameters__set__nprr(nprr)
    
    @property
    def nssp(self):
        """
        Element nssp ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 22
        
        """
        return _AresMainPy.f90wrap_parameters__get__nssp()
    
    @nssp.setter
    def nssp(self, nssp):
        _AresMainPy.f90wrap_parameters__set__nssp(nssp)
    
    @property
    def kspacing(self):
        """
        Element kspacing ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 24
        
        """
        return _AresMainPy.f90wrap_parameters__get__kspacing()
    
    @kspacing.setter
    def kspacing(self, kspacing):
        _AresMainPy.f90wrap_parameters__set__kspacing(kspacing)
    
    @property
    def kshift(self):
        """
        Element kshift ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__kshift(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kshift = self._arrays[array_handle]
        else:
            kshift = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__kshift)
            self._arrays[array_handle] = kshift
        return kshift
    
    @kshift.setter
    def kshift(self, kshift):
        self.kshift[...] = kshift
    
    @property
    def init_gap(self):
        """
        Element init_gap ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 27
        
        """
        return _AresMainPy.f90wrap_parameters__get__init_gap()
    
    @init_gap.setter
    def init_gap(self, init_gap):
        _AresMainPy.f90wrap_parameters__set__init_gap(init_gap)
    
    @property
    def ecut(self):
        """
        Element ecut ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 27
        
        """
        return _AresMainPy.f90wrap_parameters__get__ecut()
    
    @ecut.setter
    def ecut(self, ecut):
        _AresMainPy.f90wrap_parameters__set__ecut(ecut)
    
    @property
    def lkp(self):
        """
        Element lkp ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 29
        
        """
        return _AresMainPy.f90wrap_parameters__get__lkp()
    
    @lkp.setter
    def lkp(self, lkp):
        _AresMainPy.f90wrap_parameters__set__lkp(lkp)
    
    @property
    def pp_identifer(self):
        """
        Element pp_identifer ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 30
        
        """
        return _AresMainPy.f90wrap_parameters__get__pp_identifer()
    
    @pp_identifer.setter
    def pp_identifer(self, pp_identifer):
        _AresMainPy.f90wrap_parameters__set__pp_identifer(pp_identifer)
    
    @property
    def lpbc(self):
        """
        Element lpbc ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 33
        
        """
        return _AresMainPy.f90wrap_parameters__get__lpbc()
    
    @lpbc.setter
    def lpbc(self, lpbc):
        _AresMainPy.f90wrap_parameters__set__lpbc(lpbc)
    
    @property
    def calforce(self):
        """
        Element calforce ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 34
        
        """
        return _AresMainPy.f90wrap_parameters__get__calforce()
    
    @calforce.setter
    def calforce(self, calforce):
        _AresMainPy.f90wrap_parameters__set__calforce(calforce)
    
    @property
    def hartree_method(self):
        """
        Element hartree_method ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 41
        
        """
        return _AresMainPy.f90wrap_parameters__get__hartree_method()
    
    @hartree_method.setter
    def hartree_method(self, hartree_method):
        _AresMainPy.f90wrap_parameters__set__hartree_method(hartree_method)
    
    @property
    def lcell(self):
        """
        Element lcell ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 42
        
        """
        return _AresMainPy.f90wrap_parameters__get__lcell()
    
    @lcell.setter
    def lcell(self, lcell):
        _AresMainPy.f90wrap_parameters__set__lcell(lcell)
    
    @property
    def lcellbohr(self):
        """
        Element lcellbohr ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 42
        
        """
        return _AresMainPy.f90wrap_parameters__get__lcellbohr()
    
    @lcellbohr.setter
    def lcellbohr(self, lcellbohr):
        _AresMainPy.f90wrap_parameters__set__lcellbohr(lcellbohr)
    
    @property
    def isormax(self):
        """
        Element isormax ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 43
        
        """
        return _AresMainPy.f90wrap_parameters__get__isormax()
    
    @isormax.setter
    def isormax(self, isormax):
        _AresMainPy.f90wrap_parameters__set__isormax(isormax)
    
    @property
    def isormaxbohr(self):
        """
        Element isormaxbohr ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 43
        
        """
        return _AresMainPy.f90wrap_parameters__get__isormaxbohr()
    
    @isormaxbohr.setter
    def isormaxbohr(self, isormaxbohr):
        _AresMainPy.f90wrap_parameters__set__isormaxbohr(isormaxbohr)
    
    @property
    def radiusmax(self):
        """
        Element radiusmax ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 44
        
        """
        return _AresMainPy.f90wrap_parameters__get__radiusmax()
    
    @radiusmax.setter
    def radiusmax(self, radiusmax):
        _AresMainPy.f90wrap_parameters__set__radiusmax(radiusmax)
    
    @property
    def nvc(self):
        """
        Element nvc ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 45
        
        """
        return _AresMainPy.f90wrap_parameters__get__nvc()
    
    @nvc.setter
    def nvc(self, nvc):
        _AresMainPy.f90wrap_parameters__set__nvc(nvc)
    
    @property
    def isonorder(self):
        """
        Element isonorder ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 45
        
        """
        return _AresMainPy.f90wrap_parameters__get__isonorder()
    
    @isonorder.setter
    def isonorder(self, isonorder):
        _AresMainPy.f90wrap_parameters__set__isonorder(isonorder)
    
    @property
    def tolcg(self):
        """
        Element tolcg ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 46
        
        """
        return _AresMainPy.f90wrap_parameters__get__tolcg()
    
    @tolcg.setter
    def tolcg(self, tolcg):
        _AresMainPy.f90wrap_parameters__set__tolcg(tolcg)
    
    @property
    def nfcd(self):
        """
        Element nfcd ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 47
        
        """
        return _AresMainPy.f90wrap_parameters__get__nfcd()
    
    @nfcd.setter
    def nfcd(self, nfcd):
        _AresMainPy.f90wrap_parameters__set__nfcd(nfcd)
    
    @property
    def isolmax(self):
        """
        Element isolmax ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 48
        
        """
        return _AresMainPy.f90wrap_parameters__get__isolmax()
    
    @isolmax.setter
    def isolmax(self, isolmax):
        _AresMainPy.f90wrap_parameters__set__isolmax(isolmax)
    
    @property
    def iprec_fmm(self):
        """
        Element iprec_fmm ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 49
        
        """
        return _AresMainPy.f90wrap_parameters__get__iprec_fmm()
    
    @iprec_fmm.setter
    def iprec_fmm(self, iprec_fmm):
        _AresMainPy.f90wrap_parameters__set__iprec_fmm(iprec_fmm)
    
    @property
    def addcharge(self):
        """
        Element addcharge ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 50
        
        """
        return _AresMainPy.f90wrap_parameters__get__addcharge()
    
    @addcharge.setter
    def addcharge(self, addcharge):
        _AresMainPy.f90wrap_parameters__set__addcharge(addcharge)
    
    @property
    def cell_shape(self):
        """
        Element cell_shape ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 51
        
        """
        return _AresMainPy.f90wrap_parameters__get__cell_shape()
    
    @cell_shape.setter
    def cell_shape(self, cell_shape):
        _AresMainPy.f90wrap_parameters__set__cell_shape(cell_shape)
    
    @property
    def cell_thick(self):
        """
        Element cell_thick ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 52
        
        """
        return _AresMainPy.f90wrap_parameters__get__cell_thick()
    
    @cell_thick.setter
    def cell_thick(self, cell_thick):
        _AresMainPy.f90wrap_parameters__set__cell_thick(cell_thick)
    
    @property
    def lpbc2iso(self):
        """
        Element lpbc2iso ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 53
        
        """
        return _AresMainPy.f90wrap_parameters__get__lpbc2iso()
    
    @lpbc2iso.setter
    def lpbc2iso(self, lpbc2iso):
        _AresMainPy.f90wrap_parameters__set__lpbc2iso(lpbc2iso)
    
    @property
    def lradius_auto(self):
        """
        Element lradius_auto ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 54
        
        """
        return _AresMainPy.f90wrap_parameters__get__lradius_auto()
    
    @lradius_auto.setter
    def lradius_auto(self, lradius_auto):
        _AresMainPy.f90wrap_parameters__set__lradius_auto(lradius_auto)
    
    @property
    def block_mbnb(self):
        """
        Element block_mbnb ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 56
        
        """
        return _AresMainPy.f90wrap_parameters__get__block_mbnb()
    
    @block_mbnb.setter
    def block_mbnb(self, block_mbnb):
        _AresMainPy.f90wrap_parameters__set__block_mbnb(block_mbnb)
    
    @property
    def debug_out(self):
        """
        Element debug_out ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 57
        
        """
        return _AresMainPy.f90wrap_parameters__get__debug_out()
    
    @debug_out.setter
    def debug_out(self, debug_out):
        _AresMainPy.f90wrap_parameters__set__debug_out(debug_out)
    
    @property
    def wexict(self):
        """
        Element wexict ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 59
        
        """
        return _AresMainPy.f90wrap_parameters__get__wexict()
    
    @wexict.setter
    def wexict(self, wexict):
        _AresMainPy.f90wrap_parameters__set__wexict(wexict)
    
    @property
    def lbvk(self):
        """
        Element lbvk ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 66
        
        """
        return _AresMainPy.f90wrap_parameters__get__lbvk()
    
    @lbvk.setter
    def lbvk(self, lbvk):
        _AresMainPy.f90wrap_parameters__set__lbvk(lbvk)
    
    @property
    def lfirst(self):
        """
        Element lfirst ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 66
        
        """
        return _AresMainPy.f90wrap_parameters__get__lfirst()
    
    @lfirst.setter
    def lfirst(self, lfirst):
        _AresMainPy.f90wrap_parameters__set__lfirst(lfirst)
    
    @property
    def linrho(self):
        """
        Element linrho ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 66
        
        """
        return _AresMainPy.f90wrap_parameters__get__linrho()
    
    @linrho.setter
    def linrho(self, linrho):
        _AresMainPy.f90wrap_parameters__set__linrho(linrho)
    
    @property
    def lradrho(self):
        """
        Element lradrho ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 66
        
        """
        return _AresMainPy.f90wrap_parameters__get__lradrho()
    
    @lradrho.setter
    def lradrho(self, lradrho):
        _AresMainPy.f90wrap_parameters__set__lradrho(lradrho)
    
    @property
    def lrrorthnorm(self):
        """
        Element lrrorthnorm ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 66
        
        """
        return _AresMainPy.f90wrap_parameters__get__lrrorthnorm()
    
    @lrrorthnorm.setter
    def lrrorthnorm(self, lrrorthnorm):
        _AresMainPy.f90wrap_parameters__set__lrrorthnorm(lrrorthnorm)
    
    @property
    def lrandom(self):
        """
        Element lrandom ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 66
        
        """
        return _AresMainPy.f90wrap_parameters__get__lrandom()
    
    @lrandom.setter
    def lrandom(self, lrandom):
        _AresMainPy.f90wrap_parameters__set__lrandom(lrandom)
    
    @property
    def system_name(self):
        """
        Element system_name ftype=character(30) pytype=str
        
        
        Defined at Parameters.fpp line 68
        
        """
        return _AresMainPy.f90wrap_parameters__get__system_name()
    
    @system_name.setter
    def system_name(self, system_name):
        _AresMainPy.f90wrap_parameters__set__system_name(system_name)
    
    @property
    def cellfile_name(self):
        """
        Element cellfile_name ftype=character(30) pytype=str
        
        
        Defined at Parameters.fpp line 69
        
        """
        return _AresMainPy.f90wrap_parameters__get__cellfile_name()
    
    @cellfile_name.setter
    def cellfile_name(self, cellfile_name):
        _AresMainPy.f90wrap_parameters__set__cellfile_name(cellfile_name)
    
    @property
    def outfile(self):
        """
        Element outfile ftype=character(30) pytype=str
        
        
        Defined at Parameters.fpp line 70
        
        """
        return _AresMainPy.f90wrap_parameters__get__outfile()
    
    @outfile.setter
    def outfile(self, outfile):
        _AresMainPy.f90wrap_parameters__set__outfile(outfile)
    
    @property
    def ppfile_name(self):
        """
        Element ppfile_name ftype=character(30) pytype=str
        
        
        Defined at Parameters.fpp line 71
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__ppfile_name(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ppfile_name = self._arrays[array_handle]
        else:
            ppfile_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__ppfile_name)
            self._arrays[array_handle] = ppfile_name
        return ppfile_name
    
    @ppfile_name.setter
    def ppfile_name(self, ppfile_name):
        self.ppfile_name[...] = ppfile_name
    
    @property
    def elements(self):
        """
        Element elements ftype=character(30) pytype=str
        
        
        Defined at Parameters.fpp line 72
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__elements(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            elements = self._arrays[array_handle]
        else:
            elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__elements)
            self._arrays[array_handle] = elements
        return elements
    
    @elements.setter
    def elements(self, elements):
        self.elements[...] = elements
    
    @property
    def nspin(self):
        """
        Element nspin ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 74
        
        """
        return _AresMainPy.f90wrap_parameters__get__nspin()
    
    @nspin.setter
    def nspin(self, nspin):
        _AresMainPy.f90wrap_parameters__set__nspin(nspin)
    
    @property
    def nsmear(self):
        """
        Element nsmear ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 76
        
        """
        return _AresMainPy.f90wrap_parameters__get__nsmear()
    
    @nsmear.setter
    def nsmear(self, nsmear):
        _AresMainPy.f90wrap_parameters__set__nsmear(nsmear)
    
    @property
    def wsmear(self):
        """
        Element wsmear ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 77
        
        """
        return _AresMainPy.f90wrap_parameters__get__wsmear()
    
    @wsmear.setter
    def wsmear(self, wsmear):
        _AresMainPy.f90wrap_parameters__set__wsmear(wsmear)
    
    @property
    def imixer(self):
        """
        Element imixer ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 80
        
        """
        return _AresMainPy.f90wrap_parameters__get__imixer()
    
    @imixer.setter
    def imixer(self, imixer):
        _AresMainPy.f90wrap_parameters__set__imixer(imixer)
    
    @property
    def nmiter(self):
        """
        Element nmiter ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 80
        
        """
        return _AresMainPy.f90wrap_parameters__get__nmiter()
    
    @nmiter.setter
    def nmiter(self, nmiter):
        _AresMainPy.f90wrap_parameters__set__nmiter(nmiter)
    
    @property
    def nsmix(self):
        """
        Element nsmix ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 83
        
        """
        return _AresMainPy.f90wrap_parameters__get__nsmix()
    
    @nsmix.setter
    def nsmix(self, nsmix):
        _AresMainPy.f90wrap_parameters__set__nsmix(nsmix)
    
    @property
    def nhmix(self):
        """
        Element nhmix ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 83
        
        """
        return _AresMainPy.f90wrap_parameters__get__nhmix()
    
    @nhmix.setter
    def nhmix(self, nhmix):
        _AresMainPy.f90wrap_parameters__set__nhmix(nhmix)
    
    @property
    def nhmin(self):
        """
        Element nhmin ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 83
        
        """
        return _AresMainPy.f90wrap_parameters__get__nhmin()
    
    @nhmin.setter
    def nhmin(self, nhmin):
        _AresMainPy.f90wrap_parameters__set__nhmin(nhmin)
    
    @property
    def malpha(self):
        """
        Element malpha ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 87
        
        """
        return _AresMainPy.f90wrap_parameters__get__malpha()
    
    @malpha.setter
    def malpha(self, malpha):
        _AresMainPy.f90wrap_parameters__set__malpha(malpha)
    
    @property
    def mbeta(self):
        """
        Element mbeta ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 87
        
        """
        return _AresMainPy.f90wrap_parameters__get__mbeta()
    
    @mbeta.setter
    def mbeta(self, mbeta):
        _AresMainPy.f90wrap_parameters__set__mbeta(mbeta)
    
    @property
    def amix(self):
        """
        Element amix ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 87
        
        """
        return _AresMainPy.f90wrap_parameters__get__amix()
    
    @amix.setter
    def amix(self, amix):
        _AresMainPy.f90wrap_parameters__set__amix(amix)
    
    @property
    def bmix(self):
        """
        Element bmix ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 87
        
        """
        return _AresMainPy.f90wrap_parameters__get__bmix()
    
    @bmix.setter
    def bmix(self, bmix):
        _AresMainPy.f90wrap_parameters__set__bmix(bmix)
    
    @property
    def resta(self):
        """
        Element resta ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 88
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__resta(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            resta = self._arrays[array_handle]
        else:
            resta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__resta)
            self._arrays[array_handle] = resta
        return resta
    
    @resta.setter
    def resta(self, resta):
        self.resta[...] = resta
    
    @property
    def w0am(self):
        """
        Element w0am ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 90
        
        """
        return _AresMainPy.f90wrap_parameters__get__w0am()
    
    @w0am.setter
    def w0am(self, w0am):
        _AresMainPy.f90wrap_parameters__set__w0am(w0am)
    
    @property
    def rtol(self):
        """
        Element rtol ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 92
        
        """
        return _AresMainPy.f90wrap_parameters__get__rtol()
    
    @rtol.setter
    def rtol(self, rtol):
        _AresMainPy.f90wrap_parameters__set__rtol(rtol)
    
    @property
    def etol(self):
        """
        Element etol ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 92
        
        """
        return _AresMainPy.f90wrap_parameters__get__etol()
    
    @etol.setter
    def etol(self, etol):
        _AresMainPy.f90wrap_parameters__set__etol(etol)
    
    @property
    def lsub(self):
        """
        Element lsub ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 94
        
        """
        return _AresMainPy.f90wrap_parameters__get__lsub()
    
    @lsub.setter
    def lsub(self, lsub):
        _AresMainPy.f90wrap_parameters__set__lsub(lsub)
    
    @property
    def lfat(self):
        """
        Element lfat ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 95
        
        """
        return _AresMainPy.f90wrap_parameters__get__lfat()
    
    @lfat.setter
    def lfat(self, lfat):
        _AresMainPy.f90wrap_parameters__set__lfat(lfat)
    
    @property
    def lone(self):
        """
        Element lone ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 96
        
        """
        return _AresMainPy.f90wrap_parameters__get__lone()
    
    @lone.setter
    def lone(self, lone):
        _AresMainPy.f90wrap_parameters__set__lone(lone)
    
    @property
    def nsub(self):
        """
        Element nsub ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 97
        
        """
        return _AresMainPy.f90wrap_parameters__get__nsub()
    
    @nsub.setter
    def nsub(self, nsub):
        _AresMainPy.f90wrap_parameters__set__nsub(nsub)
    
    @property
    def nfat(self):
        """
        Element nfat ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 98
        
        """
        return _AresMainPy.f90wrap_parameters__get__nfat()
    
    @nfat.setter
    def nfat(self, nfat):
        _AresMainPy.f90wrap_parameters__set__nfat(nfat)
    
    @property
    def tfvw(self):
        """
        Element tfvw ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 99
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__tfvw(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tfvw = self._arrays[array_handle]
        else:
            tfvw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__tfvw)
            self._arrays[array_handle] = tfvw
        return tfvw
    
    @tfvw.setter
    def tfvw(self, tfvw):
        self.tfvw[...] = tfvw
    
    @property
    def lofdft(self):
        """
        Element lofdft ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 101
        
        """
        return _AresMainPy.f90wrap_parameters__get__lofdft()
    
    @lofdft.setter
    def lofdft(self, lofdft):
        _AresMainPy.f90wrap_parameters__set__lofdft(lofdft)
    
    @property
    def lband(self):
        """
        Element lband ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 103
        
        """
        return _AresMainPy.f90wrap_parameters__get__lband()
    
    @lband.setter
    def lband(self, lband):
        _AresMainPy.f90wrap_parameters__set__lband(lband)
    
    @property
    def iopm(self):
        """
        Element iopm ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 106
        
        """
        return _AresMainPy.f90wrap_parameters__get__iopm()
    
    @iopm.setter
    def iopm(self, iopm):
        _AresMainPy.f90wrap_parameters__set__iopm(iopm)
    
    @property
    def igoal(self):
        """
        Element igoal ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 106
        
        """
        return _AresMainPy.f90wrap_parameters__get__igoal()
    
    @igoal.setter
    def igoal(self, igoal):
        _AresMainPy.f90wrap_parameters__set__igoal(igoal)
    
    @property
    def press(self):
        """
        Element press ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 110
        
        """
        return _AresMainPy.f90wrap_parameters__get__press()
    
    @press.setter
    def press(self, press):
        _AresMainPy.f90wrap_parameters__set__press(press)
    
    @property
    def tolf(self):
        """
        Element tolf ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 110
        
        """
        return _AresMainPy.f90wrap_parameters__get__tolf()
    
    @tolf.setter
    def tolf(self, tolf):
        _AresMainPy.f90wrap_parameters__set__tolf(tolf)
    
    @property
    def tolp(self):
        """
        Element tolp ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 110
        
        """
        return _AresMainPy.f90wrap_parameters__get__tolp()
    
    @tolp.setter
    def tolp(self, tolp):
        _AresMainPy.f90wrap_parameters__set__tolp(tolp)
    
    @property
    def times(self):
        """
        Element times ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 110
        
        """
        return _AresMainPy.f90wrap_parameters__get__times()
    
    @times.setter
    def times(self, times):
        _AresMainPy.f90wrap_parameters__set__times(times)
    
    @property
    def isp(self):
        """
        Element isp ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 111
        
        """
        return _AresMainPy.f90wrap_parameters__get__isp()
    
    @isp.setter
    def isp(self, isp):
        _AresMainPy.f90wrap_parameters__set__isp(isp)
    
    @property
    def ldg(self):
        """
        Element ldg ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 113
        
        """
        return _AresMainPy.f90wrap_parameters__get__ldg()
    
    @ldg.setter
    def ldg(self, ldg):
        _AresMainPy.f90wrap_parameters__set__ldg(ldg)
    
    @property
    def ndg(self):
        """
        Element ndg ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 114
        
        """
        return _AresMainPy.f90wrap_parameters__get__ndg()
    
    @ndg.setter
    def ndg(self, ndg):
        _AresMainPy.f90wrap_parameters__set__ndg(ndg)
    
    @property
    def n_near(self):
        """
        Element n_near ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 115
        
        """
        return _AresMainPy.f90wrap_parameters__get__n_near()
    
    @n_near.setter
    def n_near(self, n_near):
        _AresMainPy.f90wrap_parameters__set__n_near(n_near)
    
    @property
    def inpol(self):
        """
        Element inpol ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 116
        
        """
        return _AresMainPy.f90wrap_parameters__get__inpol()
    
    @inpol.setter
    def inpol(self, inpol):
        _AresMainPy.f90wrap_parameters__set__inpol(inpol)
    
    @property
    def idinit(self):
        """
        Element idinit ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 118
        
        """
        return _AresMainPy.f90wrap_parameters__get__idinit()
    
    @idinit.setter
    def idinit(self, idinit):
        _AresMainPy.f90wrap_parameters__set__idinit(idinit)
    
    @property
    def mo_file(self):
        """
        Element mo_file ftype=character(30) pytype=str
        
        
        Defined at Parameters.fpp line 119
        
        """
        return _AresMainPy.f90wrap_parameters__get__mo_file()
    
    @mo_file.setter
    def mo_file(self, mo_file):
        _AresMainPy.f90wrap_parameters__set__mo_file(mo_file)
    
    @property
    def potim(self):
        """
        Element potim ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 121
        
        """
        return _AresMainPy.f90wrap_parameters__get__potim()
    
    @potim.setter
    def potim(self, potim):
        _AresMainPy.f90wrap_parameters__set__potim(potim)
    
    @property
    def ediffg(self):
        """
        Element ediffg ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 122
        
        """
        return _AresMainPy.f90wrap_parameters__get__ediffg()
    
    @ediffg.setter
    def ediffg(self, ediffg):
        _AresMainPy.f90wrap_parameters__set__ediffg(ediffg)
    
    @property
    def step_fixrho(self):
        """
        Element step_fixrho ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 123
        
        """
        return _AresMainPy.f90wrap_parameters__get__step_fixrho()
    
    @step_fixrho.setter
    def step_fixrho(self, step_fixrho):
        _AresMainPy.f90wrap_parameters__set__step_fixrho(step_fixrho)
    
    @property
    def fix_xyz(self):
        """
        Element fix_xyz ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 136
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__fix_xyz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fix_xyz = self._arrays[array_handle]
        else:
            fix_xyz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__fix_xyz)
            self._arrays[array_handle] = fix_xyz
        return fix_xyz
    
    @fix_xyz.setter
    def fix_xyz(self, fix_xyz):
        self.fix_xyz[...] = fix_xyz
    
    @property
    def pstress(self):
        """
        Element pstress ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 147
        
        """
        return _AresMainPy.f90wrap_parameters__get__pstress()
    
    @pstress.setter
    def pstress(self, pstress):
        _AresMainPy.f90wrap_parameters__set__pstress(pstress)
    
    @property
    def ramax(self):
        """
        Element ramax ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 148
        
        """
        return _AresMainPy.f90wrap_parameters__get__ramax()
    
    @ramax.setter
    def ramax(self, ramax):
        _AresMainPy.f90wrap_parameters__set__ramax(ramax)
    
    @property
    def grad_order(self):
        """
        Element grad_order ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 150
        
        """
        return _AresMainPy.f90wrap_parameters__get__grad_order()
    
    @grad_order.setter
    def grad_order(self, grad_order):
        _AresMainPy.f90wrap_parameters__set__grad_order(grad_order)
    
    @property
    def maxnpts(self):
        """
        Element maxnpts ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 151
        
        """
        return _AresMainPy.f90wrap_parameters__get__maxnpts()
    
    @maxnpts.setter
    def maxnpts(self, maxnpts):
        _AresMainPy.f90wrap_parameters__set__maxnpts(maxnpts)
    
    @property
    def nevshift(self):
        """
        Element nevshift ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 152
        
        """
        return _AresMainPy.f90wrap_parameters__get__nevshift()
    
    @nevshift.setter
    def nevshift(self, nevshift):
        _AresMainPy.f90wrap_parameters__set__nevshift(nevshift)
    
    @property
    def gridn(self):
        """
        Element gridn ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 153
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__gridn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gridn = self._arrays[array_handle]
        else:
            gridn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__gridn)
            self._arrays[array_handle] = gridn
        return gridn
    
    @gridn.setter
    def gridn(self, gridn):
        self.gridn[...] = gridn
    
    @property
    def lgamma(self):
        """
        Element lgamma ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 154
        
        """
        return _AresMainPy.f90wrap_parameters__get__lgamma()
    
    @lgamma.setter
    def lgamma(self, lgamma):
        _AresMainPy.f90wrap_parameters__set__lgamma(lgamma)
    
    @property
    def lcore_val(self):
        """
        Element lcore_val ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 155
        
        """
        return _AresMainPy.f90wrap_parameters__get__lcore_val()
    
    @lcore_val.setter
    def lcore_val(self, lcore_val):
        _AresMainPy.f90wrap_parameters__set__lcore_val(lcore_val)
    
    @property
    def iounits(self):
        """
        Element iounits ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 156
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__iounits(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iounits = self._arrays[array_handle]
        else:
            iounits = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__iounits)
            self._arrays[array_handle] = iounits
        return iounits
    
    @iounits.setter
    def iounits(self, iounits):
        self.iounits[...] = iounits
    
    @property
    def isym(self):
        """
        Element isym ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 157
        
        """
        return _AresMainPy.f90wrap_parameters__get__isym()
    
    @isym.setter
    def isym(self, isym):
        _AresMainPy.f90wrap_parameters__set__isym(isym)
    
    @property
    def lfinite_full(self):
        """
        Element lfinite_full ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 158
        
        """
        return _AresMainPy.f90wrap_parameters__get__lfinite_full()
    
    @lfinite_full.setter
    def lfinite_full(self, lfinite_full):
        _AresMainPy.f90wrap_parameters__set__lfinite_full(lfinite_full)
    
    @property
    def lforce(self):
        """
        Element lforce ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 159
        
        """
        return _AresMainPy.f90wrap_parameters__get__lforce()
    
    @lforce.setter
    def lforce(self, lforce):
        _AresMainPy.f90wrap_parameters__set__lforce(lforce)
    
    @property
    def lstress(self):
        """
        Element lstress ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 160
        
        """
        return _AresMainPy.f90wrap_parameters__get__lstress()
    
    @lstress.setter
    def lstress(self, lstress):
        _AresMainPy.f90wrap_parameters__set__lstress(lstress)
    
    @property
    def temper_elec(self):
        """
        Element temper_elec ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 161
        
        """
        return _AresMainPy.f90wrap_parameters__get__temper_elec()
    
    @temper_elec.setter
    def temper_elec(self, temper_elec):
        _AresMainPy.f90wrap_parameters__set__temper_elec(temper_elec)
    
    @property
    def lwstyle(self):
        """
        Element lwstyle ftype=character(len=20) pytype=str
        
        
        Defined at Parameters.fpp line 162
        
        """
        return _AresMainPy.f90wrap_parameters__get__lwstyle()
    
    @lwstyle.setter
    def lwstyle(self, lwstyle):
        _AresMainPy.f90wrap_parameters__set__lwstyle(lwstyle)
    
    @property
    def lopt(self):
        """
        Element lopt ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 163
        
        """
        return _AresMainPy.f90wrap_parameters__get__lopt()
    
    @lopt.setter
    def lopt(self, lopt):
        _AresMainPy.f90wrap_parameters__set__lopt(lopt)
    
    @property
    def ibron(self):
        """
        Element ibron ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 164
        
        """
        return _AresMainPy.f90wrap_parameters__get__ibron()
    
    @ibron.setter
    def ibron(self, ibron):
        _AresMainPy.f90wrap_parameters__set__ibron(ibron)
    
    @property
    def maxsave(self):
        """
        Element maxsave ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 165
        
        """
        return _AresMainPy.f90wrap_parameters__get__maxsave()
    
    @maxsave.setter
    def maxsave(self, maxsave):
        _AresMainPy.f90wrap_parameters__set__maxsave(maxsave)
    
    @property
    def igamma(self):
        """
        Element igamma ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 166
        
        """
        return _AresMainPy.f90wrap_parameters__get__igamma()
    
    @igamma.setter
    def igamma(self, igamma):
        _AresMainPy.f90wrap_parameters__set__igamma(igamma)
    
    @property
    def kgrid(self):
        """
        Element kgrid ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 167
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__kgrid(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kgrid = self._arrays[array_handle]
        else:
            kgrid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__kgrid)
            self._arrays[array_handle] = kgrid
        return kgrid
    
    @kgrid.setter
    def kgrid(self, kgrid):
        self.kgrid[...] = kgrid
    
    @property
    def lmd(self):
        """
        Element lmd ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 169
        
        """
        return _AresMainPy.f90wrap_parameters__get__lmd()
    
    @lmd.setter
    def lmd(self, lmd):
        _AresMainPy.f90wrap_parameters__set__lmd(lmd)
    
    @property
    def nhis(self):
        """
        Element nhis ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 171
        
        """
        return _AresMainPy.f90wrap_parameters__get__nhis()
    
    @nhis.setter
    def nhis(self, nhis):
        _AresMainPy.f90wrap_parameters__set__nhis(nhis)
    
    @property
    def rdfr(self):
        """
        Element rdfr ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 172
        
        """
        return _AresMainPy.f90wrap_parameters__get__rdfr()
    
    @rdfr.setter
    def rdfr(self, rdfr):
        _AresMainPy.f90wrap_parameters__set__rdfr(rdfr)
    
    @property
    def thermostat(self):
        """
        Element thermostat ftype=character(len=30) pytype=str
        
        
        Defined at Parameters.fpp line 173
        
        """
        return _AresMainPy.f90wrap_parameters__get__thermostat()
    
    @thermostat.setter
    def thermostat(self, thermostat):
        _AresMainPy.f90wrap_parameters__set__thermostat(thermostat)
    
    @property
    def integrator(self):
        """
        Element integrator ftype=character(len=30) pytype=str
        
        
        Defined at Parameters.fpp line 174
        
        """
        return _AresMainPy.f90wrap_parameters__get__integrator()
    
    @integrator.setter
    def integrator(self, integrator):
        _AresMainPy.f90wrap_parameters__set__integrator(integrator)
    
    @property
    def sfreq(self):
        """
        Element sfreq ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 176
        
        """
        return _AresMainPy.f90wrap_parameters__get__sfreq()
    
    @sfreq.setter
    def sfreq(self, sfreq):
        _AresMainPy.f90wrap_parameters__set__sfreq(sfreq)
    
    @property
    def cfreq(self):
        """
        Element cfreq ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 177
        
        """
        return _AresMainPy.f90wrap_parameters__get__cfreq()
    
    @cfreq.setter
    def cfreq(self, cfreq):
        _AresMainPy.f90wrap_parameters__set__cfreq(cfreq)
    
    @property
    def temperature(self):
        """
        Element temperature ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 178
        
        """
        return _AresMainPy.f90wrap_parameters__get__temperature()
    
    @temperature.setter
    def temperature(self, temperature):
        _AresMainPy.f90wrap_parameters__set__temperature(temperature)
    
    @property
    def mste(self):
        """
        Element mste ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 179
        
        """
        return _AresMainPy.f90wrap_parameters__get__mste()
    
    @mste.setter
    def mste(self, mste):
        _AresMainPy.f90wrap_parameters__set__mste(mste)
    
    @property
    def delt(self):
        """
        Element delt ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 180
        
        """
        return _AresMainPy.f90wrap_parameters__get__delt()
    
    @delt.setter
    def delt(self, delt):
        _AresMainPy.f90wrap_parameters__set__delt(delt)
    
    @property
    def relaxt(self):
        """
        Element relaxt ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 181
        
        """
        return _AresMainPy.f90wrap_parameters__get__relaxt()
    
    @relaxt.setter
    def relaxt(self, relaxt):
        _AresMainPy.f90wrap_parameters__set__relaxt(relaxt)
    
    @property
    def ensemble(self):
        """
        Element ensemble ftype=character(10) pytype=str
        
        
        Defined at Parameters.fpp line 182
        
        """
        return _AresMainPy.f90wrap_parameters__get__ensemble()
    
    @ensemble.setter
    def ensemble(self, ensemble):
        _AresMainPy.f90wrap_parameters__set__ensemble(ensemble)
    
    @property
    def dof(self):
        """
        Element dof ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 183
        
        """
        return _AresMainPy.f90wrap_parameters__get__dof()
    
    @dof.setter
    def dof(self, dof):
        _AresMainPy.f90wrap_parameters__set__dof(dof)
    
    @property
    def pext(self):
        """
        Element pext ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 184
        
        """
        return _AresMainPy.f90wrap_parameters__get__pext()
    
    @pext.setter
    def pext(self, pext):
        _AresMainPy.f90wrap_parameters__set__pext(pext)
    
    @property
    def iresmd(self):
        """
        Element iresmd ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 185
        
        """
        return _AresMainPy.f90wrap_parameters__get__iresmd()
    
    @iresmd.setter
    def iresmd(self, iresmd):
        _AresMainPy.f90wrap_parameters__set__iresmd(iresmd)
    
    @property
    def nresn(self):
        """
        Element nresn ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 188
        
        """
        return _AresMainPy.f90wrap_parameters__get__nresn()
    
    @nresn.setter
    def nresn(self, nresn):
        _AresMainPy.f90wrap_parameters__set__nresn(nresn)
    
    @property
    def nyosh(self):
        """
        Element nyosh ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 189
        
        """
        return _AresMainPy.f90wrap_parameters__get__nyosh()
    
    @nyosh.setter
    def nyosh(self, nyosh):
        _AresMainPy.f90wrap_parameters__set__nyosh(nyosh)
    
    @property
    def nnhc(self):
        """
        Element nnhc ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 190
        
        """
        return _AresMainPy.f90wrap_parameters__get__nnhc()
    
    @nnhc.setter
    def nnhc(self, nnhc):
        _AresMainPy.f90wrap_parameters__set__nnhc(nnhc)
    
    @property
    def wdti2(self):
        """
        Element wdti2 ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 191
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__wdti2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wdti2 = self._arrays[array_handle]
        else:
            wdti2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__wdti2)
            self._arrays[array_handle] = wdti2
        return wdti2
    
    @wdti2.setter
    def wdti2(self, wdti2):
        self.wdti2[...] = wdti2
    
    @property
    def wdti4(self):
        """
        Element wdti4 ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 191
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__wdti4(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wdti4 = self._arrays[array_handle]
        else:
            wdti4 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__wdti4)
            self._arrays[array_handle] = wdti4
        return wdti4
    
    @wdti4.setter
    def wdti4(self, wdti4):
        self.wdti4[...] = wdti4
    
    @property
    def wdti8(self):
        """
        Element wdti8 ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 191
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__wdti8(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wdti8 = self._arrays[array_handle]
        else:
            wdti8 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__wdti8)
            self._arrays[array_handle] = wdti8
        return wdti8
    
    @wdti8.setter
    def wdti8(self, wdti8):
        self.wdti8[...] = wdti8
    
    @property
    def qmass(self):
        """
        Element qmass ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 192
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__qmass(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            qmass = self._arrays[array_handle]
        else:
            qmass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__qmass)
            self._arrays[array_handle] = qmass
        return qmass
    
    @qmass.setter
    def qmass(self, qmass):
        self.qmass[...] = qmass
    
    @property
    def syin_coeff(self):
        """
        Element syin_coeff ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 193
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__syin_coeff(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            syin_coeff = self._arrays[array_handle]
        else:
            syin_coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__syin_coeff)
            self._arrays[array_handle] = syin_coeff
        return syin_coeff
    
    @syin_coeff.setter
    def syin_coeff(self, syin_coeff):
        self.syin_coeff[...] = syin_coeff
    
    @property
    def bmass(self):
        """
        Element bmass ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 194
        
        """
        return _AresMainPy.f90wrap_parameters__get__bmass()
    
    @bmass.setter
    def bmass(self, bmass):
        _AresMainPy.f90wrap_parameters__set__bmass(bmass)
    
    @property
    def fthermo(self):
        """
        Element fthermo ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 195
        
        """
        return _AresMainPy.f90wrap_parameters__get__fthermo()
    
    @fthermo.setter
    def fthermo(self, fthermo):
        _AresMainPy.f90wrap_parameters__set__fthermo(fthermo)
    
    @property
    def fbaro(self):
        """
        Element fbaro ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 196
        
        """
        return _AresMainPy.f90wrap_parameters__get__fbaro()
    
    @fbaro.setter
    def fbaro(self, fbaro):
        _AresMainPy.f90wrap_parameters__set__fbaro(fbaro)
    
    @property
    def fthermown(self):
        """
        Element fthermown ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 197
        
        """
        return _AresMainPy.f90wrap_parameters__get__fthermown()
    
    @fthermown.setter
    def fthermown(self, fthermown):
        _AresMainPy.f90wrap_parameters__set__fthermown(fthermown)
    
    @property
    def fbarown(self):
        """
        Element fbarown ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 198
        
        """
        return _AresMainPy.f90wrap_parameters__get__fbarown()
    
    @fbarown.setter
    def fbarown(self, fbarown):
        _AresMainPy.f90wrap_parameters__set__fbarown(fbarown)
    
    @property
    def lbin(self):
        """
        Element lbin ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 200
        
        """
        return _AresMainPy.f90wrap_parameters__get__lbin()
    
    @lbin.setter
    def lbin(self, lbin):
        _AresMainPy.f90wrap_parameters__set__lbin(lbin)
    
    @property
    def p_flag(self):
        """
        Element p_flag ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 202
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__p_flag(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            p_flag = self._arrays[array_handle]
        else:
            p_flag = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__p_flag)
            self._arrays[array_handle] = p_flag
        return p_flag
    
    @p_flag.setter
    def p_flag(self, p_flag):
        self.p_flag[...] = p_flag
    
    @property
    def p_start(self):
        """
        Element p_start ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 203
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__p_start(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            p_start = self._arrays[array_handle]
        else:
            p_start = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__p_start)
            self._arrays[array_handle] = p_start
        return p_start
    
    @p_start.setter
    def p_start(self, p_start):
        self.p_start[...] = p_start
    
    @property
    def p_stop(self):
        """
        Element p_stop ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 203
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__p_stop(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            p_stop = self._arrays[array_handle]
        else:
            p_stop = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__p_stop)
            self._arrays[array_handle] = p_stop
        return p_stop
    
    @p_stop.setter
    def p_stop(self, p_stop):
        self.p_stop[...] = p_stop
    
    @property
    def p_freq(self):
        """
        Element p_freq ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 203
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__p_freq(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            p_freq = self._arrays[array_handle]
        else:
            p_freq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__p_freq)
            self._arrays[array_handle] = p_freq
        return p_freq
    
    @p_freq.setter
    def p_freq(self, p_freq):
        self.p_freq[...] = p_freq
    
    @property
    def p_target(self):
        """
        Element p_target ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 203
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_parameters__array__p_target(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            p_target = self._arrays[array_handle]
        else:
            p_target = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_parameters__array__p_target)
            self._arrays[array_handle] = p_target
        return p_target
    
    @p_target.setter
    def p_target(self, p_target):
        self.p_target[...] = p_target
    
    @property
    def t_start(self):
        """
        Element t_start ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 204
        
        """
        return _AresMainPy.f90wrap_parameters__get__t_start()
    
    @t_start.setter
    def t_start(self, t_start):
        _AresMainPy.f90wrap_parameters__set__t_start(t_start)
    
    @property
    def t_stop(self):
        """
        Element t_stop ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 204
        
        """
        return _AresMainPy.f90wrap_parameters__get__t_stop()
    
    @t_stop.setter
    def t_stop(self, t_stop):
        _AresMainPy.f90wrap_parameters__set__t_stop(t_stop)
    
    @property
    def t_freq(self):
        """
        Element t_freq ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 204
        
        """
        return _AresMainPy.f90wrap_parameters__get__t_freq()
    
    @t_freq.setter
    def t_freq(self, t_freq):
        _AresMainPy.f90wrap_parameters__set__t_freq(t_freq)
    
    @property
    def t_target(self):
        """
        Element t_target ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 204
        
        """
        return _AresMainPy.f90wrap_parameters__get__t_target()
    
    @t_target.setter
    def t_target(self, t_target):
        _AresMainPy.f90wrap_parameters__set__t_target(t_target)
    
    @property
    def press_control(self):
        """
        Element press_control ftype=character(len=clen) pytype=str
        
        
        Defined at Parameters.fpp line 205
        
        """
        return _AresMainPy.f90wrap_parameters__get__press_control()
    
    @press_control.setter
    def press_control(self, press_control):
        _AresMainPy.f90wrap_parameters__set__press_control(press_control)
    
    @property
    def pcontrol(self):
        """
        Element pcontrol ftype=character(len=clen) pytype=str
        
        
        Defined at Parameters.fpp line 206
        
        """
        return _AresMainPy.f90wrap_parameters__get__pcontrol()
    
    @pcontrol.setter
    def pcontrol(self, pcontrol):
        _AresMainPy.f90wrap_parameters__set__pcontrol(pcontrol)
    
    @property
    def pmass(self):
        """
        Element pmass ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 207
        
        """
        return _AresMainPy.f90wrap_parameters__get__pmass()
    
    @pmass.setter
    def pmass(self, pmass):
        _AresMainPy.f90wrap_parameters__set__pmass(pmass)
    
    @property
    def pdrag(self):
        """
        Element pdrag ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 208
        
        """
        return _AresMainPy.f90wrap_parameters__get__pdrag()
    
    @pdrag.setter
    def pdrag(self, pdrag):
        _AresMainPy.f90wrap_parameters__set__pdrag(pdrag)
    
    @property
    def tdrag(self):
        """
        Element tdrag ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 209
        
        """
        return _AresMainPy.f90wrap_parameters__get__tdrag()
    
    @tdrag.setter
    def tdrag(self, tdrag):
        _AresMainPy.f90wrap_parameters__set__tdrag(tdrag)
    
    @property
    def erate(self):
        """
        Element erate ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 210
        
        """
        return _AresMainPy.f90wrap_parameters__get__erate()
    
    @erate.setter
    def erate(self, erate):
        _AresMainPy.f90wrap_parameters__set__erate(erate)
    
    @property
    def mstrain(self):
        """
        Element mstrain ftype=real(dp) pytype=float
        
        
        Defined at Parameters.fpp line 211
        
        """
        return _AresMainPy.f90wrap_parameters__get__mstrain()
    
    @mstrain.setter
    def mstrain(self, mstrain):
        _AresMainPy.f90wrap_parameters__set__mstrain(mstrain)
    
    @property
    def sdir(self):
        """
        Element sdir ftype=character(len=clen) pytype=str
        
        
        Defined at Parameters.fpp line 212
        
        """
        return _AresMainPy.f90wrap_parameters__get__sdir()
    
    @sdir.setter
    def sdir(self, sdir):
        _AresMainPy.f90wrap_parameters__set__sdir(sdir)
    
    @property
    def lke_pot(self):
        """
        Element lke_pot ftype=logical pytype=bool
        
        
        Defined at Parameters.fpp line 213
        
        """
        return _AresMainPy.f90wrap_parameters__get__lke_pot()
    
    @lke_pot.setter
    def lke_pot(self, lke_pot):
        _AresMainPy.f90wrap_parameters__set__lke_pot(lke_pot)
    
    def __str__(self):
        ret = ['<parameters>{\n']
        ret.append('    ikedf : ')
        ret.append(repr(self.ikedf))
        ret.append(',\n    ixcdf : ')
        ret.append(repr(self.ixcdf))
        ret.append(',\n    finite_order : ')
        ret.append(repr(self.finite_order))
        ret.append(',\n    ntype : ')
        ret.append(repr(self.ntype))
        ret.append(',\n    naddstates : ')
        ret.append(repr(self.naddstates))
        ret.append(',\n    istart : ')
        ret.append(repr(self.istart))
        ret.append(',\n    idiag : ')
        ret.append(repr(self.idiag))
        ret.append(',\n    chem : ')
        ret.append(repr(self.chem))
        ret.append(',\n    chem0 : ')
        ret.append(repr(self.chem0))
        ret.append(',\n    nstates : ')
        ret.append(repr(self.nstates))
        ret.append(',\n    nstates_global : ')
        ret.append(repr(self.nstates_global))
        ret.append(',\n    nprr : ')
        ret.append(repr(self.nprr))
        ret.append(',\n    nssp : ')
        ret.append(repr(self.nssp))
        ret.append(',\n    kspacing : ')
        ret.append(repr(self.kspacing))
        ret.append(',\n    kshift : ')
        ret.append(repr(self.kshift))
        ret.append(',\n    init_gap : ')
        ret.append(repr(self.init_gap))
        ret.append(',\n    ecut : ')
        ret.append(repr(self.ecut))
        ret.append(',\n    lkp : ')
        ret.append(repr(self.lkp))
        ret.append(',\n    pp_identifer : ')
        ret.append(repr(self.pp_identifer))
        ret.append(',\n    lpbc : ')
        ret.append(repr(self.lpbc))
        ret.append(',\n    calforce : ')
        ret.append(repr(self.calforce))
        ret.append(',\n    hartree_method : ')
        ret.append(repr(self.hartree_method))
        ret.append(',\n    lcell : ')
        ret.append(repr(self.lcell))
        ret.append(',\n    lcellbohr : ')
        ret.append(repr(self.lcellbohr))
        ret.append(',\n    isormax : ')
        ret.append(repr(self.isormax))
        ret.append(',\n    isormaxbohr : ')
        ret.append(repr(self.isormaxbohr))
        ret.append(',\n    radiusmax : ')
        ret.append(repr(self.radiusmax))
        ret.append(',\n    nvc : ')
        ret.append(repr(self.nvc))
        ret.append(',\n    isonorder : ')
        ret.append(repr(self.isonorder))
        ret.append(',\n    tolcg : ')
        ret.append(repr(self.tolcg))
        ret.append(',\n    nfcd : ')
        ret.append(repr(self.nfcd))
        ret.append(',\n    isolmax : ')
        ret.append(repr(self.isolmax))
        ret.append(',\n    iprec_fmm : ')
        ret.append(repr(self.iprec_fmm))
        ret.append(',\n    addcharge : ')
        ret.append(repr(self.addcharge))
        ret.append(',\n    cell_shape : ')
        ret.append(repr(self.cell_shape))
        ret.append(',\n    cell_thick : ')
        ret.append(repr(self.cell_thick))
        ret.append(',\n    lpbc2iso : ')
        ret.append(repr(self.lpbc2iso))
        ret.append(',\n    lradius_auto : ')
        ret.append(repr(self.lradius_auto))
        ret.append(',\n    block_mbnb : ')
        ret.append(repr(self.block_mbnb))
        ret.append(',\n    debug_out : ')
        ret.append(repr(self.debug_out))
        ret.append(',\n    wexict : ')
        ret.append(repr(self.wexict))
        ret.append(',\n    lbvk : ')
        ret.append(repr(self.lbvk))
        ret.append(',\n    lfirst : ')
        ret.append(repr(self.lfirst))
        ret.append(',\n    linrho : ')
        ret.append(repr(self.linrho))
        ret.append(',\n    lradrho : ')
        ret.append(repr(self.lradrho))
        ret.append(',\n    lrrorthnorm : ')
        ret.append(repr(self.lrrorthnorm))
        ret.append(',\n    lrandom : ')
        ret.append(repr(self.lrandom))
        ret.append(',\n    system_name : ')
        ret.append(repr(self.system_name))
        ret.append(',\n    cellfile_name : ')
        ret.append(repr(self.cellfile_name))
        ret.append(',\n    outfile : ')
        ret.append(repr(self.outfile))
        ret.append(',\n    ppfile_name : ')
        ret.append(repr(self.ppfile_name))
        ret.append(',\n    elements : ')
        ret.append(repr(self.elements))
        ret.append(',\n    nspin : ')
        ret.append(repr(self.nspin))
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
        ret.append(',\n    nhmix : ')
        ret.append(repr(self.nhmix))
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
        ret.append(',\n    lsub : ')
        ret.append(repr(self.lsub))
        ret.append(',\n    lfat : ')
        ret.append(repr(self.lfat))
        ret.append(',\n    lone : ')
        ret.append(repr(self.lone))
        ret.append(',\n    nsub : ')
        ret.append(repr(self.nsub))
        ret.append(',\n    nfat : ')
        ret.append(repr(self.nfat))
        ret.append(',\n    tfvw : ')
        ret.append(repr(self.tfvw))
        ret.append(',\n    lofdft : ')
        ret.append(repr(self.lofdft))
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
        ret.append(',\n    isp : ')
        ret.append(repr(self.isp))
        ret.append(',\n    ldg : ')
        ret.append(repr(self.ldg))
        ret.append(',\n    ndg : ')
        ret.append(repr(self.ndg))
        ret.append(',\n    n_near : ')
        ret.append(repr(self.n_near))
        ret.append(',\n    inpol : ')
        ret.append(repr(self.inpol))
        ret.append(',\n    idinit : ')
        ret.append(repr(self.idinit))
        ret.append(',\n    mo_file : ')
        ret.append(repr(self.mo_file))
        ret.append(',\n    potim : ')
        ret.append(repr(self.potim))
        ret.append(',\n    ediffg : ')
        ret.append(repr(self.ediffg))
        ret.append(',\n    step_fixrho : ')
        ret.append(repr(self.step_fixrho))
        ret.append(',\n    fix_xyz : ')
        ret.append(repr(self.fix_xyz))
        ret.append(',\n    pstress : ')
        ret.append(repr(self.pstress))
        ret.append(',\n    ramax : ')
        ret.append(repr(self.ramax))
        ret.append(',\n    grad_order : ')
        ret.append(repr(self.grad_order))
        ret.append(',\n    maxnpts : ')
        ret.append(repr(self.maxnpts))
        ret.append(',\n    nevshift : ')
        ret.append(repr(self.nevshift))
        ret.append(',\n    gridn : ')
        ret.append(repr(self.gridn))
        ret.append(',\n    lgamma : ')
        ret.append(repr(self.lgamma))
        ret.append(',\n    lcore_val : ')
        ret.append(repr(self.lcore_val))
        ret.append(',\n    iounits : ')
        ret.append(repr(self.iounits))
        ret.append(',\n    isym : ')
        ret.append(repr(self.isym))
        ret.append(',\n    lfinite_full : ')
        ret.append(repr(self.lfinite_full))
        ret.append(',\n    lforce : ')
        ret.append(repr(self.lforce))
        ret.append(',\n    lstress : ')
        ret.append(repr(self.lstress))
        ret.append(',\n    temper_elec : ')
        ret.append(repr(self.temper_elec))
        ret.append(',\n    lwstyle : ')
        ret.append(repr(self.lwstyle))
        ret.append(',\n    lopt : ')
        ret.append(repr(self.lopt))
        ret.append(',\n    ibron : ')
        ret.append(repr(self.ibron))
        ret.append(',\n    maxsave : ')
        ret.append(repr(self.maxsave))
        ret.append(',\n    igamma : ')
        ret.append(repr(self.igamma))
        ret.append(',\n    kgrid : ')
        ret.append(repr(self.kgrid))
        ret.append(',\n    lmd : ')
        ret.append(repr(self.lmd))
        ret.append(',\n    nhis : ')
        ret.append(repr(self.nhis))
        ret.append(',\n    rdfr : ')
        ret.append(repr(self.rdfr))
        ret.append(',\n    thermostat : ')
        ret.append(repr(self.thermostat))
        ret.append(',\n    integrator : ')
        ret.append(repr(self.integrator))
        ret.append(',\n    sfreq : ')
        ret.append(repr(self.sfreq))
        ret.append(',\n    cfreq : ')
        ret.append(repr(self.cfreq))
        ret.append(',\n    temperature : ')
        ret.append(repr(self.temperature))
        ret.append(',\n    mste : ')
        ret.append(repr(self.mste))
        ret.append(',\n    delt : ')
        ret.append(repr(self.delt))
        ret.append(',\n    relaxt : ')
        ret.append(repr(self.relaxt))
        ret.append(',\n    ensemble : ')
        ret.append(repr(self.ensemble))
        ret.append(',\n    dof : ')
        ret.append(repr(self.dof))
        ret.append(',\n    pext : ')
        ret.append(repr(self.pext))
        ret.append(',\n    iresmd : ')
        ret.append(repr(self.iresmd))
        ret.append(',\n    nresn : ')
        ret.append(repr(self.nresn))
        ret.append(',\n    nyosh : ')
        ret.append(repr(self.nyosh))
        ret.append(',\n    nnhc : ')
        ret.append(repr(self.nnhc))
        ret.append(',\n    wdti2 : ')
        ret.append(repr(self.wdti2))
        ret.append(',\n    wdti4 : ')
        ret.append(repr(self.wdti4))
        ret.append(',\n    wdti8 : ')
        ret.append(repr(self.wdti8))
        ret.append(',\n    qmass : ')
        ret.append(repr(self.qmass))
        ret.append(',\n    syin_coeff : ')
        ret.append(repr(self.syin_coeff))
        ret.append(',\n    bmass : ')
        ret.append(repr(self.bmass))
        ret.append(',\n    fthermo : ')
        ret.append(repr(self.fthermo))
        ret.append(',\n    fbaro : ')
        ret.append(repr(self.fbaro))
        ret.append(',\n    fthermown : ')
        ret.append(repr(self.fthermown))
        ret.append(',\n    fbarown : ')
        ret.append(repr(self.fbarown))
        ret.append(',\n    lbin : ')
        ret.append(repr(self.lbin))
        ret.append(',\n    p_flag : ')
        ret.append(repr(self.p_flag))
        ret.append(',\n    p_start : ')
        ret.append(repr(self.p_start))
        ret.append(',\n    p_stop : ')
        ret.append(repr(self.p_stop))
        ret.append(',\n    p_freq : ')
        ret.append(repr(self.p_freq))
        ret.append(',\n    p_target : ')
        ret.append(repr(self.p_target))
        ret.append(',\n    t_start : ')
        ret.append(repr(self.t_start))
        ret.append(',\n    t_stop : ')
        ret.append(repr(self.t_stop))
        ret.append(',\n    t_freq : ')
        ret.append(repr(self.t_freq))
        ret.append(',\n    t_target : ')
        ret.append(repr(self.t_target))
        ret.append(',\n    press_control : ')
        ret.append(repr(self.press_control))
        ret.append(',\n    pcontrol : ')
        ret.append(repr(self.pcontrol))
        ret.append(',\n    pmass : ')
        ret.append(repr(self.pmass))
        ret.append(',\n    pdrag : ')
        ret.append(repr(self.pdrag))
        ret.append(',\n    tdrag : ')
        ret.append(repr(self.tdrag))
        ret.append(',\n    erate : ')
        ret.append(repr(self.erate))
        ret.append(',\n    mstrain : ')
        ret.append(repr(self.mstrain))
        ret.append(',\n    sdir : ')
        ret.append(repr(self.sdir))
        ret.append(',\n    lke_pot : ')
        ret.append(repr(self.lke_pot))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

parameters = Parameters()

class Smpi_Math_Module(f90wrap.runtime.FortranModule):
    """
    Module smpi_math_module
    
    
    Defined at Smpi_math_module.fpp lines 5-2616
    
    """
    @f90wrap.runtime.register_class("AresMainPy.parallel_type")
    class parallel_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=parallel_type)
        
        
        Defined at Smpi_math_module.fpp lines 10-45
        
        """
        def __init__(self, handle=None):
            """
            self = Parallel_Type()
            
            
            Defined at Smpi_math_module.fpp lines 10-45
            
            
            Returns
            -------
            this : Parallel_Type
            	Object to be constructed
            
            
            Automatically generated constructor for parallel_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_parallel_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Parallel_Type
            
            
            Defined at Smpi_math_module.fpp lines 10-45
            
            Parameters
            ----------
            this : Parallel_Type
            	Object to be destructed
            
            
            Automatically generated destructor for parallel_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_parallel_type_finalise(this=self._handle)
        
        @property
        def comm(self):
            """
            Element comm ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 11
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__comm(self._handle)
        
        @comm.setter
        def comm(self, comm):
            _AresMainPy.f90wrap_parallel_type__set__comm(self._handle, comm)
        
        @property
        def myid(self):
            """
            Element myid ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 12
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__myid(self._handle)
        
        @myid.setter
        def myid(self, myid):
            _AresMainPy.f90wrap_parallel_type__set__myid(self._handle, myid)
        
        @property
        def numprocs(self):
            """
            Element numprocs ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 13
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__numprocs(self._handle)
        
        @numprocs.setter
        def numprocs(self, numprocs):
            _AresMainPy.f90wrap_parallel_type__set__numprocs(self._handle, numprocs)
        
        @property
        def rootid(self):
            """
            Element rootid ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 14
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__rootid(self._handle)
        
        @rootid.setter
        def rootid(self, rootid):
            _AresMainPy.f90wrap_parallel_type__set__rootid(self._handle, rootid)
        
        @property
        def isroot(self):
            """
            Element isroot ftype=logical pytype=bool
            
            
            Defined at Smpi_math_module.fpp line 15
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__isroot(self._handle)
        
        @isroot.setter
        def isroot(self, isroot):
            _AresMainPy.f90wrap_parallel_type__set__isroot(self._handle, isroot)
        
        @property
        def nstate_proc(self):
            """
            Element nstate_proc ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 17
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__nstate_proc(self._handle)
        
        @nstate_proc.setter
        def nstate_proc(self, nstate_proc):
            _AresMainPy.f90wrap_parallel_type__set__nstate_proc(self._handle, nstate_proc)
        
        @property
        def sub2sum(self):
            """
            Element sub2sum ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__sub2sum(self._handle)
            if array_handle in self._arrays:
                sub2sum = self._arrays[array_handle]
            else:
                sub2sum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__sub2sum)
                self._arrays[array_handle] = sub2sum
            return sub2sum
        
        @sub2sum.setter
        def sub2sum(self, sub2sum):
            self.sub2sum[...] = sub2sum
        
        @property
        def mygrid_range(self):
            """
            Element mygrid_range ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__mygrid_range(self._handle)
            if array_handle in self._arrays:
                mygrid_range = self._arrays[array_handle]
            else:
                mygrid_range = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__mygrid_range)
                self._arrays[array_handle] = mygrid_range
            return mygrid_range
        
        @mygrid_range.setter
        def mygrid_range(self, mygrid_range):
            self.mygrid_range[...] = mygrid_range
        
        @property
        def recvcounts(self):
            """
            Element recvcounts ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__recvcounts(self._handle)
            if array_handle in self._arrays:
                recvcounts = self._arrays[array_handle]
            else:
                recvcounts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__recvcounts)
                self._arrays[array_handle] = recvcounts
            return recvcounts
        
        @recvcounts.setter
        def recvcounts(self, recvcounts):
            self.recvcounts[...] = recvcounts
        
        @property
        def displs(self):
            """
            Element displs ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 22
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__displs(self._handle)
            if array_handle in self._arrays:
                displs = self._arrays[array_handle]
            else:
                displs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__displs)
                self._arrays[array_handle] = displs
            return displs
        
        @displs.setter
        def displs(self, displs):
            self.displs[...] = displs
        
        @property
        def global_gridrange(self):
            """
            Element global_gridrange ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__global_gridrange(self._handle)
            if array_handle in self._arrays:
                global_gridrange = self._arrays[array_handle]
            else:
                global_gridrange = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__global_gridrange)
                self._arrays[array_handle] = global_gridrange
            return global_gridrange
        
        @global_gridrange.setter
        def global_gridrange(self, global_gridrange):
            self.global_gridrange[...] = global_gridrange
        
        @property
        def comm2d(self):
            """
            Element comm2d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__comm2d(self._handle)
        
        @comm2d.setter
        def comm2d(self, comm2d):
            _AresMainPy.f90wrap_parallel_type__set__comm2d(self._handle, comm2d)
        
        @property
        def commx(self):
            """
            Element commx ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__commx(self._handle)
        
        @commx.setter
        def commx(self, commx):
            _AresMainPy.f90wrap_parallel_type__set__commx(self._handle, commx)
        
        @property
        def commy(self):
            """
            Element commy ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__commy(self._handle)
        
        @commy.setter
        def commy(self, commy):
            _AresMainPy.f90wrap_parallel_type__set__commy(self._handle, commy)
        
        @property
        def rankx(self):
            """
            Element rankx ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__rankx(self._handle)
        
        @rankx.setter
        def rankx(self, rankx):
            _AresMainPy.f90wrap_parallel_type__set__rankx(self._handle, rankx)
        
        @property
        def ranky(self):
            """
            Element ranky ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__ranky(self._handle)
        
        @ranky.setter
        def ranky(self, ranky):
            _AresMainPy.f90wrap_parallel_type__set__ranky(self._handle, ranky)
        
        @property
        def periods(self):
            """
            Element periods ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__periods(self._handle)
            if array_handle in self._arrays:
                periods = self._arrays[array_handle]
            else:
                periods = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__periods)
                self._arrays[array_handle] = periods
            return periods
        
        @periods.setter
        def periods(self, periods):
            self.periods[...] = periods
        
        @property
        def reorder(self):
            """
            Element reorder ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__reorder(self._handle)
        
        @reorder.setter
        def reorder(self, reorder):
            _AresMainPy.f90wrap_parallel_type__set__reorder(self._handle, reorder)
        
        @property
        def remainx(self):
            """
            Element remainx ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__remainx(self._handle)
            if array_handle in self._arrays:
                remainx = self._arrays[array_handle]
            else:
                remainx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__remainx)
                self._arrays[array_handle] = remainx
            return remainx
        
        @remainx.setter
        def remainx(self, remainx):
            self.remainx[...] = remainx
        
        @property
        def remainy(self):
            """
            Element remainy ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__remainy(self._handle)
            if array_handle in self._arrays:
                remainy = self._arrays[array_handle]
            else:
                remainy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__remainy)
                self._arrays[array_handle] = remainy
            return remainy
        
        @remainy.setter
        def remainy(self, remainy):
            self.remainy[...] = remainy
        
        @property
        def ndims(self):
            """
            Element ndims ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__ndims(self._handle)
        
        @ndims.setter
        def ndims(self, ndims):
            _AresMainPy.f90wrap_parallel_type__set__ndims(self._handle, ndims)
        
        @property
        def dims(self):
            """
            Element dims ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__dims(self._handle)
            if array_handle in self._arrays:
                dims = self._arrays[array_handle]
            else:
                dims = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__dims)
                self._arrays[array_handle] = dims
            return dims
        
        @dims.setter
        def dims(self, dims):
            self.dims[...] = dims
        
        @property
        def commfft(self):
            """
            Element commfft ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 27
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__commfft(self._handle)
        
        @commfft.setter
        def commfft(self, commfft):
            _AresMainPy.f90wrap_parallel_type__set__commfft(self._handle, commfft)
        
        @property
        def local_z(self):
            """
            Element local_z ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 27
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__local_z(self._handle)
        
        @local_z.setter
        def local_z(self, local_z):
            _AresMainPy.f90wrap_parallel_type__set__local_z(self._handle, local_z)
        
        @property
        def local_z_start(self):
            """
            Element local_z_start ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 27
            
            """
            return _AresMainPy.f90wrap_parallel_type__get__local_z_start(self._handle)
        
        @local_z_start.setter
        def local_z_start(self, local_z_start):
            _AresMainPy.f90wrap_parallel_type__set__local_z_start(self._handle, \
                local_z_start)
        
        @property
        def fft_grid_range(self):
            """
            Element fft_grid_range ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__fft_grid_range(self._handle)
            if array_handle in self._arrays:
                fft_grid_range = self._arrays[array_handle]
            else:
                fft_grid_range = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__fft_grid_range)
                self._arrays[array_handle] = fft_grid_range
            return fft_grid_range
        
        @fft_grid_range.setter
        def fft_grid_range(self, fft_grid_range):
            self.fft_grid_range[...] = fft_grid_range
        
        @property
        def fft_rcount(self):
            """
            Element fft_rcount ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__fft_rcount(self._handle)
            if array_handle in self._arrays:
                fft_rcount = self._arrays[array_handle]
            else:
                fft_rcount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__fft_rcount)
                self._arrays[array_handle] = fft_rcount
            return fft_rcount
        
        @fft_rcount.setter
        def fft_rcount(self, fft_rcount):
            self.fft_rcount[...] = fft_rcount
        
        @property
        def fft_rdispls(self):
            """
            Element fft_rdispls ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__fft_rdispls(self._handle)
            if array_handle in self._arrays:
                fft_rdispls = self._arrays[array_handle]
            else:
                fft_rdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__fft_rdispls)
                self._arrays[array_handle] = fft_rdispls
            return fft_rdispls
        
        @fft_rdispls.setter
        def fft_rdispls(self, fft_rdispls):
            self.fft_rdispls[...] = fft_rdispls
        
        @property
        def fft_scount(self):
            """
            Element fft_scount ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__fft_scount(self._handle)
            if array_handle in self._arrays:
                fft_scount = self._arrays[array_handle]
            else:
                fft_scount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__fft_scount)
                self._arrays[array_handle] = fft_scount
            return fft_scount
        
        @fft_scount.setter
        def fft_scount(self, fft_scount):
            self.fft_scount[...] = fft_scount
        
        @property
        def fft_sdispls(self):
            """
            Element fft_sdispls ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_parallel_type__array__fft_sdispls(self._handle)
            if array_handle in self._arrays:
                fft_sdispls = self._arrays[array_handle]
            else:
                fft_sdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_parallel_type__array__fft_sdispls)
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
        
    
    @f90wrap.runtime.register_class("AresMainPy.smpi_root_type")
    class smpi_root_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=smpi_root_type)
        
        
        Defined at Smpi_math_module.fpp lines 49-50
        
        """
        def __init__(self, handle=None):
            """
            self = Smpi_Root_Type()
            
            
            Defined at Smpi_math_module.fpp lines 49-50
            
            
            Returns
            -------
            this : Smpi_Root_Type
            	Object to be constructed
            
            
            Automatically generated constructor for smpi_root_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_smpi_root_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Smpi_Root_Type
            
            
            Defined at Smpi_math_module.fpp lines 49-50
            
            Parameters
            ----------
            this : Smpi_Root_Type
            	Object to be destructed
            
            
            Automatically generated destructor for smpi_root_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_smpi_root_type_finalise(this=self._handle)
        
        @property
        def natom_group(self):
            """
            Element natom_group ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 50
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_smpi_root_type__array__natom_group(self._handle)
            if array_handle in self._arrays:
                natom_group = self._arrays[array_handle]
            else:
                natom_group = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_smpi_root_type__array__natom_group)
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
        
    
    @f90wrap.runtime.register_class("AresMainPy.smpi_comm_type")
    class smpi_comm_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=smpi_comm_type)
        
        
        Defined at Smpi_math_module.fpp lines 54-56
        
        """
        def __init__(self, handle=None):
            """
            self = Smpi_Comm_Type()
            
            
            Defined at Smpi_math_module.fpp lines 54-56
            
            
            Returns
            -------
            this : Smpi_Comm_Type
            	Object to be constructed
            
            
            Automatically generated constructor for smpi_comm_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_smpi_comm_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Smpi_Comm_Type
            
            
            Defined at Smpi_math_module.fpp lines 54-56
            
            Parameters
            ----------
            this : Smpi_Comm_Type
            	Object to be destructed
            
            
            Automatically generated destructor for smpi_comm_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_smpi_comm_type_finalise(this=self._handle)
        
        @property
        def atoms(self):
            """
            Element atoms ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 55
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_smpi_comm_type__array__atoms(self._handle)
            if array_handle in self._arrays:
                atoms = self._arrays[array_handle]
            else:
                atoms = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_smpi_comm_type__array__atoms)
                self._arrays[array_handle] = atoms
            return atoms
        
        @atoms.setter
        def atoms(self, atoms):
            self.atoms[...] = atoms
        
        @property
        def displs(self):
            """
            Element displs ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 56
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_smpi_comm_type__array__displs(self._handle)
            if array_handle in self._arrays:
                displs = self._arrays[array_handle]
            else:
                displs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_smpi_comm_type__array__displs)
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
        
    
    @f90wrap.runtime.register_class("AresMainPy.time_type")
    class time_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=time_type)
        
        
        Defined at Smpi_math_module.fpp lines 60-66
        
        """
        def __init__(self, handle=None):
            """
            self = Time_Type()
            
            
            Defined at Smpi_math_module.fpp lines 60-66
            
            
            Returns
            -------
            this : Time_Type
            	Object to be constructed
            
            
            Automatically generated constructor for time_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_time_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Time_Type
            
            
            Defined at Smpi_math_module.fpp lines 60-66
            
            Parameters
            ----------
            this : Time_Type
            	Object to be destructed
            
            
            Automatically generated destructor for time_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_time_type_finalise(this=self._handle)
        
        @property
        def label(self):
            """
            Element label ftype=character(len=100) pytype=str
            
            
            Defined at Smpi_math_module.fpp line 61
            
            """
            return _AresMainPy.f90wrap_time_type__get__label(self._handle)
        
        @label.setter
        def label(self, label):
            _AresMainPy.f90wrap_time_type__set__label(self._handle, label)
        
        @property
        def tic(self):
            """
            Element tic ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.fpp line 62
            
            """
            return _AresMainPy.f90wrap_time_type__get__tic(self._handle)
        
        @tic.setter
        def tic(self, tic):
            _AresMainPy.f90wrap_time_type__set__tic(self._handle, tic)
        
        @property
        def toc(self):
            """
            Element toc ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.fpp line 63
            
            """
            return _AresMainPy.f90wrap_time_type__get__toc(self._handle)
        
        @toc.setter
        def toc(self, toc):
            _AresMainPy.f90wrap_time_type__set__toc(self._handle, toc)
        
        @property
        def total(self):
            """
            Element total ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.fpp line 64
            
            """
            return _AresMainPy.f90wrap_time_type__get__total(self._handle)
        
        @total.setter
        def total(self, total):
            _AresMainPy.f90wrap_time_type__set__total(self._handle, total)
        
        @property
        def sum_total(self):
            """
            Element sum_total ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.fpp line 65
            
            """
            return _AresMainPy.f90wrap_time_type__get__sum_total(self._handle)
        
        @sum_total.setter
        def sum_total(self, sum_total):
            _AresMainPy.f90wrap_time_type__set__sum_total(self._handle, sum_total)
        
        @property
        def num(self):
            """
            Element num ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 66
            
            """
            return _AresMainPy.f90wrap_time_type__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _AresMainPy.f90wrap_time_type__set__num(self._handle, num)
        
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
        
    
    @f90wrap.runtime.register_class("AresMainPy.mem_type")
    class mem_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=mem_type)
        
        
        Defined at Smpi_math_module.fpp lines 68-72
        
        """
        def __init__(self, handle=None):
            """
            self = Mem_Type()
            
            
            Defined at Smpi_math_module.fpp lines 68-72
            
            
            Returns
            -------
            this : Mem_Type
            	Object to be constructed
            
            
            Automatically generated constructor for mem_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_mem_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Mem_Type
            
            
            Defined at Smpi_math_module.fpp lines 68-72
            
            Parameters
            ----------
            this : Mem_Type
            	Object to be destructed
            
            
            Automatically generated destructor for mem_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_mem_type_finalise(this=self._handle)
        
        @property
        def label(self):
            """
            Element label ftype=character(len=100) pytype=str
            
            
            Defined at Smpi_math_module.fpp line 69
            
            """
            return _AresMainPy.f90wrap_mem_type__get__label(self._handle)
        
        @label.setter
        def label(self, label):
            _AresMainPy.f90wrap_mem_type__set__label(self._handle, label)
        
        @property
        def memic(self):
            """
            Element memic ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.fpp line 70
            
            """
            return _AresMainPy.f90wrap_mem_type__get__memic(self._handle)
        
        @memic.setter
        def memic(self, memic):
            _AresMainPy.f90wrap_mem_type__set__memic(self._handle, memic)
        
        @property
        def total(self):
            """
            Element total ftype=real(dp) pytype=float
            
            
            Defined at Smpi_math_module.fpp line 71
            
            """
            return _AresMainPy.f90wrap_mem_type__get__total(self._handle)
        
        @total.setter
        def total(self, total):
            _AresMainPy.f90wrap_mem_type__set__total(self._handle, total)
        
        @property
        def num(self):
            """
            Element num ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 72
            
            """
            return _AresMainPy.f90wrap_mem_type__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _AresMainPy.f90wrap_mem_type__set__num(self._handle, num)
        
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
        
    
    @f90wrap.runtime.register_class("AresMainPy.grid_diff_map_type")
    class grid_diff_map_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=grid_diff_map_type)
        
        
        Defined at Smpi_math_module.fpp lines 79-86
        
        """
        def __init__(self, handle=None):
            """
            self = Grid_Diff_Map_Type()
            
            
            Defined at Smpi_math_module.fpp lines 79-86
            
            
            Returns
            -------
            this : Grid_Diff_Map_Type
            	Object to be constructed
            
            
            Automatically generated constructor for grid_diff_map_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_grid_diff_map_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Grid_Diff_Map_Type
            
            
            Defined at Smpi_math_module.fpp lines 79-86
            
            Parameters
            ----------
            this : Grid_Diff_Map_Type
            	Object to be destructed
            
            
            Automatically generated destructor for grid_diff_map_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_grid_diff_map_type_finalise(this=self._handle)
        
        @property
        def nz_map(self):
            """
            Element nz_map ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 80
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__nz_map(self._handle)
            if array_handle in self._arrays:
                nz_map = self._arrays[array_handle]
            else:
                nz_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__nz_map)
                self._arrays[array_handle] = nz_map
            return nz_map
        
        @nz_map.setter
        def nz_map(self, nz_map):
            self.nz_map[...] = nz_map
        
        @property
        def mycomm_cores(self):
            """
            Element mycomm_cores ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 81
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__mycomm_cores(self._handle)
            if array_handle in self._arrays:
                mycomm_cores = self._arrays[array_handle]
            else:
                mycomm_cores = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__mycomm_cores)
                self._arrays[array_handle] = mycomm_cores
            return mycomm_cores
        
        @mycomm_cores.setter
        def mycomm_cores(self, mycomm_cores):
            self.mycomm_cores[...] = mycomm_cores
        
        @property
        def mycomm_size(self):
            """
            Element mycomm_size ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 82
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__mycomm_size(self._handle)
            if array_handle in self._arrays:
                mycomm_size = self._arrays[array_handle]
            else:
                mycomm_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__mycomm_size)
                self._arrays[array_handle] = mycomm_size
            return mycomm_size
        
        @mycomm_size.setter
        def mycomm_size(self, mycomm_size):
            self.mycomm_size[...] = mycomm_size
        
        @property
        def mysend_size(self):
            """
            Element mysend_size ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 83
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__mysend_size(self._handle)
            if array_handle in self._arrays:
                mysend_size = self._arrays[array_handle]
            else:
                mysend_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__mysend_size)
                self._arrays[array_handle] = mysend_size
            return mysend_size
        
        @mysend_size.setter
        def mysend_size(self, mysend_size):
            self.mysend_size[...] = mysend_size
        
        @property
        def local_map(self):
            """
            Element local_map ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 84
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__local_map(self._handle)
            if array_handle in self._arrays:
                local_map = self._arrays[array_handle]
            else:
                local_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__local_map)
                self._arrays[array_handle] = local_map
            return local_map
        
        @local_map.setter
        def local_map(self, local_map):
            self.local_map[...] = local_map
        
        @property
        def local_map1d(self):
            """
            Element local_map1d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 85
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__local_map1d(self._handle)
            if array_handle in self._arrays:
                local_map1d = self._arrays[array_handle]
            else:
                local_map1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__local_map1d)
                self._arrays[array_handle] = local_map1d
            return local_map1d
        
        @local_map1d.setter
        def local_map1d(self, local_map1d):
            self.local_map1d[...] = local_map1d
        
        @property
        def boundary(self):
            """
            Element boundary ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 86
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_diff_map_type__array__boundary(self._handle)
            if array_handle in self._arrays:
                boundary = self._arrays[array_handle]
            else:
                boundary = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_diff_map_type__array__boundary)
                self._arrays[array_handle] = boundary
            return boundary
        
        @boundary.setter
        def boundary(self, boundary):
            self.boundary[...] = boundary
        
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
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("AresMainPy.sphere_type")
    class sphere_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=sphere_type)
        
        
        Defined at Smpi_math_module.fpp lines 88-90
        
        """
        def __init__(self, handle=None):
            """
            self = Sphere_Type()
            
            
            Defined at Smpi_math_module.fpp lines 88-90
            
            
            Returns
            -------
            this : Sphere_Type
            	Object to be constructed
            
            
            Automatically generated constructor for sphere_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_sphere_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Sphere_Type
            
            
            Defined at Smpi_math_module.fpp lines 88-90
            
            Parameters
            ----------
            this : Sphere_Type
            	Object to be destructed
            
            
            Automatically generated destructor for sphere_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_sphere_type_finalise(this=self._handle)
        
        @property
        def length(self):
            """
            Element length ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 89
            
            """
            return _AresMainPy.f90wrap_sphere_type__get__length(self._handle)
        
        @length.setter
        def length(self, length):
            _AresMainPy.f90wrap_sphere_type__set__length(self._handle, length)
        
        @property
        def map3d(self):
            """
            Element map3d ftype=integer(i4b) pytype=int
            
            
            Defined at Smpi_math_module.fpp line 90
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_sphere_type__array__map3d(self._handle)
            if array_handle in self._arrays:
                map3d = self._arrays[array_handle]
            else:
                map3d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_sphere_type__array__map3d)
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
        
        
        Defined at Smpi_math_module.fpp lines 165-176
        
        
        """
        _AresMainPy.f90wrap_smpi_init()
    
    @staticmethod
    def smpi_init_comm(lpbc):
        """
        smpi_init_comm(lpbc)
        
        
        Defined at Smpi_math_module.fpp lines 180-187
        
        Parameters
        ----------
        lpbc : bool
        
        """
        _AresMainPy.f90wrap_smpi_init_comm(lpbc=lpbc)
    
    @staticmethod
    def smpi_init_per():
        """
        smpi_init_per()
        
        
        Defined at Smpi_math_module.fpp lines 190-234
        
        
        """
        _AresMainPy.f90wrap_smpi_init_per()
    
    @staticmethod
    def smpi_init_iso():
        """
        smpi_init_iso()
        
        
        Defined at Smpi_math_module.fpp lines 238-265
        
        
        """
        _AresMainPy.f90wrap_smpi_init_iso()
    
    @staticmethod
    def smpi_exit():
        """
        smpi_exit()
        
        
        Defined at Smpi_math_module.fpp lines 336-340
        
        
        """
        _AresMainPy.f90wrap_smpi_exit()
    
    @staticmethod
    def smpi_stop(message):
        """
        smpi_stop(message)
        
        
        Defined at Smpi_math_module.fpp lines 343-348
        
        Parameters
        ----------
        message : str
        
        """
        _AresMainPy.f90wrap_smpi_stop(message=message)
    
    @staticmethod
    def smpi_stop_info(message):
        """
        smpi_stop_info(message)
        
        
        Defined at Smpi_math_module.fpp lines 351-356
        
        Parameters
        ----------
        message : str
        
        """
        _AresMainPy.f90wrap_smpi_stop_info(message=message)
    
    @staticmethod
    def nstates_split(m, np):
        """
        nstates_split(m, np)
        
        
        Defined at Smpi_math_module.fpp lines 360-379
        
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
        _AresMainPy.f90wrap_nstates_split(m=m, np=np)
    
    @staticmethod
    def nstates_split_2(m, np):
        """
        nstates_split_2(m, np)
        
        
        Defined at Smpi_math_module.fpp lines 383-402
        
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
        _AresMainPy.f90wrap_nstates_split_2(m=m, np=np)
    
    @staticmethod
    def smpi_reduce_sum_real(amat, na, ramat=None):
        """
        smpi_reduce_sum_real(amat, na[, ramat])
        
        
        Defined at Smpi_math_module.fpp lines 818-827
        
        Parameters
        ----------
        amat : float array
        na : int
        ramat : float array
        
        """
        _AresMainPy.f90wrap_smpi_reduce_sum_real(amat=amat, na=na, ramat=ramat)
    
    @staticmethod
    def start_time(inlabel, flag, tic=None):
        """
        start_time(inlabel, flag[, tic])
        
        
        Defined at Smpi_math_module.fpp lines 943-963
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        tic : float
        
        """
        _AresMainPy.f90wrap_start_time(inlabel=inlabel, flag=flag, tic=tic)
    
    @staticmethod
    def end_time(inlabel, flag, toc=None):
        """
        end_time(inlabel, flag[, toc])
        
        
        Defined at Smpi_math_module.fpp lines 967-986
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        toc : float
        
        """
        _AresMainPy.f90wrap_end_time(inlabel=inlabel, flag=flag, toc=toc)
    
    @staticmethod
    def write_time(inlabel, flag):
        """
        write_time(inlabel, flag)
        
        
        Defined at Smpi_math_module.fpp lines 990-1002
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        
        """
        _AresMainPy.f90wrap_write_time(inlabel=inlabel, flag=flag)
    
    @staticmethod
    def write_sum_time(inlabel, flag):
        """
        write_sum_time(inlabel, flag)
        
        
        Defined at Smpi_math_module.fpp lines 1006-1018
        
        Parameters
        ----------
        inlabel : str
        flag : bool
        
        """
        _AresMainPy.f90wrap_write_sum_time(inlabel=inlabel, flag=flag)
    
    @staticmethod
    def print_time(inlabel, t):
        """
        print_time(inlabel, t)
        
        
        Defined at Smpi_math_module.fpp lines 1022-1033
        
        Parameters
        ----------
        inlabel : str
        t : float
        
        """
        _AresMainPy.f90wrap_print_time(inlabel=inlabel, t=t)
    
    @staticmethod
    def states_split(nev):
        """
        states_split(nev)
        
        
        Defined at Smpi_math_module.fpp lines 1073-1093
        
        Parameters
        ----------
        nev : int
        
        """
        _AresMainPy.f90wrap_states_split(nev=nev)
    
    @staticmethod
    def array_split(nev):
        """
        array_split(nev)
        
        
        Defined at Smpi_math_module.fpp lines 1096-1144
        
        Parameters
        ----------
        nev : int
        
        """
        _AresMainPy.f90wrap_array_split(nev=nev)
    
    @staticmethod
    def grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
        gridrange_sum=None, n1xy=None, n2xy=None, n3xy=None, n=None):
        """
        grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs[, \
            gridrange_sum, n1xy, n2xy, n3xy, n])
        
        
        Defined at Smpi_math_module.fpp lines 1147-1202
        
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
        n1xy : int
        n2xy : int
        n3xy : int
        n : int
        
        """
        _AresMainPy.f90wrap_grid_split(ngrid=ngrid, ncore=ncore, comm=comm, id=id, \
            grid_range=grid_range, recvcounts=recvcounts, displs=displs, \
            gridrange_sum=gridrange_sum, n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, n=n)
    
    @staticmethod
    def atom_split(mysize, natom, atom_index):
        """
        atom_split(mysize, natom, atom_index)
        
        
        Defined at Smpi_math_module.fpp lines 1205-1241
        
        Parameters
        ----------
        mysize : int
        natom : int
        atom_index : int array
        
        """
        _AresMainPy.f90wrap_atom_split(mysize=mysize, natom=natom, \
            atom_index=atom_index)
    
    @staticmethod
    def grid_sphere_init(n1, n2, n3, norder):
        """
        grid_sphere_init(n1, n2, n3, norder)
        
        
        Defined at Smpi_math_module.fpp lines 1244-1361
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        norder : int
        
        """
        _AresMainPy.f90wrap_grid_sphere_init(n1=n1, n2=n2, n3=n3, norder=norder)
    
    @staticmethod
    def set_wrap_grid_iso(myrho, wrap_box):
        """
        set_wrap_grid_iso(myrho, wrap_box)
        
        
        Defined at Smpi_math_module.fpp lines 1363-1434
        
        Parameters
        ----------
        myrho : float array
        wrap_box : float array
        
        """
        _AresMainPy.f90wrap_set_wrap_grid_iso(myrho=myrho, wrap_box=wrap_box)
    
    @staticmethod
    def destroy_diff_map():
        """
        destroy_diff_map()
        
        
        Defined at Smpi_math_module.fpp lines 1436-1464
        
        
        """
        _AresMainPy.f90wrap_destroy_diff_map()
    
    @staticmethod
    def grid_cubic_init(n1, n2, n3, n, norder):
        """
        grid_cubic_init(n1, n2, n3, n, norder)
        
        
        Defined at Smpi_math_module.fpp lines 1749-2059
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        n : int
        norder : int
        
        """
        _AresMainPy.f90wrap_grid_cubic_init(n1=n1, n2=n2, n3=n3, n=n, norder=norder)
    
    @staticmethod
    def set_wrap_grid_per(myrho, wrap_box, global_n, global_n1, global_n2):
        """
        set_wrap_grid_per(myrho, wrap_box, global_n, global_n1, global_n2)
        
        
        Defined at Smpi_math_module.fpp lines 2061-2252
        
        Parameters
        ----------
        myrho : complex array
        wrap_box : complex array
        global_n : int
        global_n1 : int
        global_n2 : int
        
        """
        _AresMainPy.f90wrap_set_wrap_grid_per(myrho=myrho, wrap_box=wrap_box, \
            global_n=global_n, global_n1=global_n1, global_n2=global_n2)
    
    @staticmethod
    def set_wrap_grid_per_ata(myrho, wrap_box1d, global_n, global_n1, global_n2):
        """
        set_wrap_grid_per_ata(myrho, wrap_box1d, global_n, global_n1, global_n2)
        
        
        Defined at Smpi_math_module.fpp lines 2254-2405
        
        Parameters
        ----------
        myrho : complex array
        wrap_box1d : complex array
        global_n : int
        global_n1 : int
        global_n2 : int
        
        """
        _AresMainPy.f90wrap_set_wrap_grid_per_ata(myrho=myrho, wrap_box1d=wrap_box1d, \
            global_n=global_n, global_n1=global_n1, global_n2=global_n2)
    
    @staticmethod
    def set_wrap_grid_per_ata_real(myrho, wrap_box1d, global_n, global_n1, \
        global_n2):
        """
        set_wrap_grid_per_ata_real(myrho, wrap_box1d, global_n, global_n1, global_n2)
        
        
        Defined at Smpi_math_module.fpp lines 2407-2556
        
        Parameters
        ----------
        myrho : float array
        wrap_box1d : float array
        global_n : int
        global_n1 : int
        global_n2 : int
        
        """
        _AresMainPy.f90wrap_set_wrap_grid_per_ata_real(myrho=myrho, \
            wrap_box1d=wrap_box1d, global_n=global_n, global_n1=global_n1, \
            global_n2=global_n2)
    
    @staticmethod
    def set_fft_alltoallv(fft_grid_range_temp):
        """
        set_fft_alltoallv(fft_grid_range_temp)
        
        
        Defined at Smpi_math_module.fpp lines 2558-2601
        
        Parameters
        ----------
        fft_grid_range_temp : int array
        
        """
        _AresMainPy.f90wrap_set_fft_alltoallv(fft_grid_range_temp=fft_grid_range_temp)
    
    @staticmethod
    def destroy_fft_alltoallv():
        """
        destroy_fft_alltoallv()
        
        
        Defined at Smpi_math_module.fpp lines 2603-2616
        
        
        """
        _AresMainPy.f90wrap_destroy_fft_alltoallv()
    
    @staticmethod
    def _sum_real_1d(amat):
        """
        totals = _sum_real_1d(amat)
        
        
        Defined at Smpi_math_module.fpp lines 406-422
        
        Parameters
        ----------
        amat : float array
        
        Returns
        -------
        totals : float
        
        """
        totals = _AresMainPy.f90wrap_sum_real_1d(amat=amat)
        return totals
    
    @staticmethod
    def _sum_real_2d(amat, bmat):
        """
        totals = _sum_real_2d(amat, bmat)
        
        
        Defined at Smpi_math_module.fpp lines 446-463
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        
        Returns
        -------
        totals : float
        
        """
        totals = _AresMainPy.f90wrap_sum_real_2d(amat=amat, bmat=bmat)
        return totals
    
    @staticmethod
    def _sum_real_3d(amat, bmat, cmat):
        """
        totals = _sum_real_3d(amat, bmat, cmat)
        
        
        Defined at Smpi_math_module.fpp lines 488-505
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        cmat : float array
        
        Returns
        -------
        totals : float
        
        """
        totals = _AresMainPy.f90wrap_sum_real_3d(amat=amat, bmat=bmat, cmat=cmat)
        return totals
    
    @staticmethod
    def _sum_cplx_1d(amat):
        """
        totals = _sum_cplx_1d(amat)
        
        
        Defined at Smpi_math_module.fpp lines 426-442
        
        Parameters
        ----------
        amat : complex array
        
        Returns
        -------
        totals : complex
        
        """
        totals = _AresMainPy.f90wrap_sum_cplx_1d(amat=amat)
        return totals
    
    @staticmethod
    def _sum_cplx_2d(amat, bmat):
        """
        totals = _sum_cplx_2d(amat, bmat)
        
        
        Defined at Smpi_math_module.fpp lines 467-484
        
        Parameters
        ----------
        amat : complex array
        bmat : complex array
        
        Returns
        -------
        totals : complex
        
        """
        totals = _AresMainPy.f90wrap_sum_cplx_2d(amat=amat, bmat=bmat)
        return totals
    
    @staticmethod
    def _sum_cplx_3d(amat, bmat, cmat):
        """
        totals = _sum_cplx_3d(amat, bmat, cmat)
        
        
        Defined at Smpi_math_module.fpp lines 509-526
        
        Parameters
        ----------
        amat : complex array
        bmat : complex array
        cmat : complex array
        
        Returns
        -------
        totals : complex
        
        """
        totals = _AresMainPy.f90wrap_sum_cplx_3d(amat=amat, bmat=bmat, cmat=cmat)
        return totals
    
    @staticmethod
    def sompsum(*args, **kwargs):
        """
        sompsum(*args, **kwargs)
        
        
        Defined at Smpi_math_module.fpp lines 99-105
        
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
        
        
        Defined at Smpi_math_module.fpp lines 572-574
        
        Parameters
        ----------
        x : int
        
        Returns
        -------
        sumx : int
        
        """
        sumx = _AresMainPy.f90wrap_smpi_sum_int_1s(x=x)
        return sumx
    
    @staticmethod
    def _smpi_sum_cplx_1s(x):
        """
        sumx = _smpi_sum_cplx_1s(x)
        
        
        Defined at Smpi_math_module.fpp lines 578-580
        
        Parameters
        ----------
        x : complex
        
        Returns
        -------
        sumx : complex
        
        """
        sumx = _AresMainPy.f90wrap_smpi_sum_cplx_1s(x=x)
        return sumx
    
    @staticmethod
    def _smpi_sum_real_1s(x):
        """
        sumx = _smpi_sum_real_1s(x)
        
        
        Defined at Smpi_math_module.fpp lines 584-586
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        sumx : float
        
        """
        sumx = _AresMainPy.f90wrap_smpi_sum_real_1s(x=x)
        return sumx
    
    @staticmethod
    def _smpi_sum_real_1d(amat):
        """
        suma = _smpi_sum_real_1d(amat)
        
        
        Defined at Smpi_math_module.fpp lines 602-605
        
        Parameters
        ----------
        amat : float array
        
        Returns
        -------
        suma : float
        
        """
        suma = _AresMainPy.f90wrap_smpi_sum_real_1d(amat=amat)
        return suma
    
    @staticmethod
    def _smpi_sum_real_2d(amat, bmat):
        """
        suma = _smpi_sum_real_2d(amat, bmat)
        
        
        Defined at Smpi_math_module.fpp lines 616-620
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        
        Returns
        -------
        suma : float
        
        """
        suma = _AresMainPy.f90wrap_smpi_sum_real_2d(amat=amat, bmat=bmat)
        return suma
    
    @staticmethod
    def _smpi_sum_real_3d(amat, bmat, cmat):
        """
        suma = _smpi_sum_real_3d(amat, bmat, cmat)
        
        
        Defined at Smpi_math_module.fpp lines 632-636
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        cmat : float array
        
        Returns
        -------
        suma : float
        
        """
        suma = _AresMainPy.f90wrap_smpi_sum_real_3d(amat=amat, bmat=bmat, cmat=cmat)
        return suma
    
    @staticmethod
    def smpisum(*args, **kwargs):
        """
        smpisum(*args, **kwargs)
        
        
        Defined at Smpi_math_module.fpp lines 107-113
        
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
        
        
        Defined at Smpi_math_module.fpp lines 1037-1045
        
        Parameters
        ----------
        munit : str
        amat : float array
        
        Returns
        -------
        summem : float
        
        """
        summem = _AresMainPy.f90wrap_smpi_sum_mem_1d(munit=munit, amat=amat)
        return summem
    
    @staticmethod
    def _smpi_sum_mem_2d(munit, amat):
        """
        summem = _smpi_sum_mem_2d(munit, amat)
        
        
        Defined at Smpi_math_module.fpp lines 1049-1057
        
        Parameters
        ----------
        munit : str
        amat : float array
        
        Returns
        -------
        summem : float
        
        """
        summem = _AresMainPy.f90wrap_smpi_sum_mem_2d(munit=munit, amat=amat)
        return summem
    
    @staticmethod
    def _smpi_sum_mem_3d(munit, amat):
        """
        summem = _smpi_sum_mem_3d(munit, amat)
        
        
        Defined at Smpi_math_module.fpp lines 1061-1069
        
        Parameters
        ----------
        munit : str
        amat : float array
        
        Returns
        -------
        summem : float
        
        """
        summem = _AresMainPy.f90wrap_smpi_sum_mem_3d(munit=munit, amat=amat)
        return summem
    
    @staticmethod
    def smpisummem(*args, **kwargs):
        """
        smpisummem(*args, **kwargs)
        
        
        Defined at Smpi_math_module.fpp lines 131-134
        
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
        
        
        Defined at Smpi_math_module.fpp lines 805-815
        
        Parameters
        ----------
        amat : float array
        ramat : float array
        
        """
        _AresMainPy.f90wrap_smpi_reduce_sum_real_1d(amat=amat, ramat=ramat)
    
    @staticmethod
    def _smpi_reduce_sum_int_1d(amat, ramat=None):
        """
        _smpi_reduce_sum_int_1d(amat[, ramat])
        
        
        Defined at Smpi_math_module.fpp lines 791-801
        
        Parameters
        ----------
        amat : int array
        ramat : int array
        
        """
        _AresMainPy.f90wrap_smpi_reduce_sum_int_1d(amat=amat, ramat=ramat)
    
    @staticmethod
    def _smpi_reduce_sum_cplx_1d(amat, ramat=None):
        """
        _smpi_reduce_sum_cplx_1d(amat[, ramat])
        
        
        Defined at Smpi_math_module.fpp lines 901-911
        
        Parameters
        ----------
        amat : complex array
        ramat : complex array
        
        """
        _AresMainPy.f90wrap_smpi_reduce_sum_cplx_1d(amat=amat, ramat=ramat)
    
    @staticmethod
    def _smpi_reduce_sum_real_2d(amat, ramat=None):
        """
        _smpi_reduce_sum_real_2d(amat[, ramat])
        
        
        Defined at Smpi_math_module.fpp lines 915-925
        
        Parameters
        ----------
        amat : float array
        ramat : float array
        
        """
        _AresMainPy.f90wrap_smpi_reduce_sum_real_2d(amat=amat, ramat=ramat)
    
    @staticmethod
    def smpireducesum(*args, **kwargs):
        """
        smpireducesum(*args, **kwargs)
        
        
        Defined at Smpi_math_module.fpp lines 136-140
        
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
        
        
        Defined at Smpi_math_module.fpp lines 551-568
        
        Parameters
        ----------
        amat : float array
        pow : int
        
        Returns
        -------
        totals : float
        
        """
        totals = _AresMainPy.f90wrap_sum_pow_int(amat=amat, pow=pow)
        return totals
    
    @staticmethod
    def _sum_pow_real(amat, pow):
        """
        totals = _sum_pow_real(amat, pow)
        
        
        Defined at Smpi_math_module.fpp lines 530-547
        
        Parameters
        ----------
        amat : float array
        pow : float
        
        Returns
        -------
        totals : float
        
        """
        totals = _AresMainPy.f90wrap_sum_pow_real(amat=amat, pow=pow)
        return totals
    
    @staticmethod
    def sompsumpow(*args, **kwargs):
        """
        sompsumpow(*args, **kwargs)
        
        
        Defined at Smpi_math_module.fpp lines 155-157
        
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
        
        
        Defined at Smpi_math_module.fpp lines 656-660
        
        Parameters
        ----------
        amat : float array
        pow : int
        
        Returns
        -------
        suma : float
        
        """
        suma = _AresMainPy.f90wrap_smpi_sum_pow_int(amat=amat, pow=pow)
        return suma
    
    @staticmethod
    def _smpi_sum_pow_real(amat, pow):
        """
        suma = _smpi_sum_pow_real(amat, pow)
        
        
        Defined at Smpi_math_module.fpp lines 648-652
        
        Parameters
        ----------
        amat : float array
        pow : float
        
        Returns
        -------
        suma : float
        
        """
        suma = _AresMainPy.f90wrap_smpi_sum_pow_real(amat=amat, pow=pow)
        return suma
    
    @staticmethod
    def smpisumpow(*args, **kwargs):
        """
        smpisumpow(*args, **kwargs)
        
        
        Defined at Smpi_math_module.fpp lines 159-160
        
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
    
    @property
    def rtic(self):
        """
        Element rtic ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 77
        
        """
        return _AresMainPy.f90wrap_smpi_math_module__get__rtic()
    
    @rtic.setter
    def rtic(self, rtic):
        _AresMainPy.f90wrap_smpi_math_module__set__rtic(rtic)
    
    @property
    def rtoc(self):
        """
        Element rtoc ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 77
        
        """
        return _AresMainPy.f90wrap_smpi_math_module__get__rtoc()
    
    @rtoc.setter
    def rtoc(self, rtoc):
        _AresMainPy.f90wrap_smpi_math_module__set__rtoc(rtoc)
    
    @property
    def mpinfo(self):
        """
        Element mpinfo ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 95
        
        """
        return _AresMainPy.f90wrap_smpi_math_module__get__mpinfo()
    
    @mpinfo.setter
    def mpinfo(self, mpinfo):
        _AresMainPy.f90wrap_smpi_math_module__set__mpinfo(mpinfo)
    
    @property
    def smpi_status(self):
        """
        Element smpi_status ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 96
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_smpi_math_module__array__smpi_status(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            smpi_status = self._arrays[array_handle]
        else:
            smpi_status = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_smpi_math_module__array__smpi_status)
            self._arrays[array_handle] = smpi_status
        return smpi_status
    
    @smpi_status.setter
    def smpi_status(self, smpi_status):
        self.smpi_status[...] = smpi_status
    
    @property
    def lall_grid(self):
        """
        Element lall_grid ftype=logical pytype=bool
        
        
        Defined at Smpi_math_module.fpp line 97
        
        """
        return _AresMainPy.f90wrap_smpi_math_module__get__lall_grid()
    
    @lall_grid.setter
    def lall_grid(self, lall_grid):
        _AresMainPy.f90wrap_smpi_math_module__set__lall_grid(lall_grid)
    
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

class Array_Io(f90wrap.runtime.FortranModule):
    """
    Module array_io
    
    
    Defined at array_io.fpp lines 5-96
    
    """
    @f90wrap.runtime.register_class("AresMainPy.out_label")
    class out_label(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=out_label)
        
        
        Defined at array_io.fpp lines 8-10
        
        """
        def __init__(self, handle=None):
            """
            self = Out_Label()
            
            
            Defined at array_io.fpp lines 8-10
            
            
            Returns
            -------
            this : Out_Label
            	Object to be constructed
            
            
            Automatically generated constructor for out_label
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_out_label_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Out_Label
            
            
            Defined at array_io.fpp lines 8-10
            
            Parameters
            ----------
            this : Out_Label
            	Object to be destructed
            
            
            Automatically generated destructor for out_label
            """
            if self._alloc:
                _AresMainPy.f90wrap_out_label_finalise(this=self._handle)
        
        @property
        def label(self):
            """
            Element label ftype=character(len=100) pytype=str
            
            
            Defined at array_io.fpp line 9
            
            """
            return _AresMainPy.f90wrap_out_label__get__label(self._handle)
        
        @label.setter
        def label(self, label):
            _AresMainPy.f90wrap_out_label__set__label(self._handle, label)
        
        @property
        def num(self):
            """
            Element num ftype=integer(i4b) pytype=int
            
            
            Defined at array_io.fpp line 10
            
            """
            return _AresMainPy.f90wrap_out_label__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _AresMainPy.f90wrap_out_label__set__num(self._handle, num)
        
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
        
        
        Defined at array_io.fpp lines 21-27
        
        
        """
        _AresMainPy.f90wrap_init_outfile()
    
    @staticmethod
    def _output_r(size_bn, array, name):
        """
        _output_r(size_bn, array, name)
        
        
        Defined at array_io.fpp lines 29-68
        
        Parameters
        ----------
        size_bn : int
        array : float array
        name : str
        
        """
        _AresMainPy.f90wrap_output_r(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def _output_i(size_bn, array, name):
        """
        _output_i(size_bn, array, name)
        
        
        Defined at array_io.fpp lines 79-86
        
        Parameters
        ----------
        size_bn : int
        array : int array
        name : str
        
        """
        _AresMainPy.f90wrap_output_i(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def output(*args, **kwargs):
        """
        output(*args, **kwargs)
        
        
        Defined at array_io.fpp lines 14-15
        
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
        
        
        Defined at array_io.fpp lines 70-77
        
        Parameters
        ----------
        size_bn : int
        array : float array
        name : str
        
        """
        _AresMainPy.f90wrap_input_r(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def _input_i(size_bn, array, name):
        """
        _input_i(size_bn, array, name)
        
        
        Defined at array_io.fpp lines 88-95
        
        Parameters
        ----------
        size_bn : int
        array : int array
        name : str
        
        """
        _AresMainPy.f90wrap_input_i(size_bn=size_bn, array=array, name=name)
    
    @staticmethod
    def input(*args, **kwargs):
        """
        input(*args, **kwargs)
        
        
        Defined at array_io.fpp lines 17-18
        
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

class Mathsplines(f90wrap.runtime.FortranModule):
    """
    Module mathsplines
    
    
    Defined at MathSplines.fpp lines 5-836
    
    """
    @staticmethod
    def spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp):
        """
        spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp)
        
        
        Defined at MathSplines.fpp lines 178-387
        
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
        _AresMainPy.f90wrap_spline_cubic_set(n=n, t=t, y=y, ibcbeg=ibcbeg, \
            ybcbeg=ybcbeg, ibcend=ibcend, ybcend=ybcend, ypp=ypp)
    
    @staticmethod
    def spline_cubic_val(n, t, y, ypp, tval):
        """
        yval, ypval, yppval = spline_cubic_val(n, t, y, ypp, tval)
        
        
        Defined at MathSplines.fpp lines 389-483
        
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
        yval, ypval, yppval = _AresMainPy.f90wrap_spline_cubic_val(n=n, t=t, y=y, \
            ypp=ypp, tval=tval)
        return yval, ypval, yppval
    
    @staticmethod
    def rvec_bracket(n, x, xval):
        """
        left, right = rvec_bracket(n, x, xval)
        
        
        Defined at MathSplines.fpp lines 485-581
        
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
        left, right = _AresMainPy.f90wrap_rvec_bracket(n=n, x=x, xval=xval)
        return left, right
    
    @staticmethod
    def s3_fs(a1, a2, a3, n, b, x):
        """
        s3_fs(a1, a2, a3, n, b, x)
        
        
        Defined at MathSplines.fpp lines 583-650
        
        Parameters
        ----------
        a1 : float array
        a2 : float array
        a3 : float array
        n : int
        b : float array
        x : float array
        
        """
        _AresMainPy.f90wrap_s3_fs(a1=a1, a2=a2, a3=a3, n=n, b=b, x=x)
    
    @staticmethod
    def polynom(m, np, xa, ya, c, x):
        """
        polynom = polynom(m, np, xa, ya, c, x)
        
        
        Defined at MathSplines.fpp lines 658-836
        
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
        polynom = _AresMainPy.f90wrap_polynom(m=m, np=np, xa=xa, ya=ya, c=c, x=x)
        return polynom
    
    _dt_array_initialisers = []
    

mathsplines = Mathsplines()

class Math(f90wrap.runtime.FortranModule):
    """
    Module math
    
    
    Defined at Math.fpp lines 5-3439
    
    """
    @staticmethod
    def change_case(instr, str, fun):
        """
        change_case(instr, str, fun)
        
        
        Defined at Math.fpp lines 55-95
        
        Parameters
        ----------
        instr : str
        str : str
        fun : int
        
        """
        _AresMainPy.f90wrap_change_case(instr=instr, str=str, fun=fun)
    
    @staticmethod
    def find_keywords(str, ch_mark, id_key, id_value):
        """
        find_keywords(str, ch_mark, id_key, id_value)
        
        
        Defined at Math.fpp lines 99-128
        
        Parameters
        ----------
        str : str
        ch_mark : str
        id_key : int
        id_value : int
        
        """
        _AresMainPy.f90wrap_find_keywords(str=str, ch_mark=ch_mark, id_key=id_key, \
            id_value=id_value)
    
    @staticmethod
    def find_nword(str, ch_comma, nword):
        """
        find_nword(str, ch_comma, nword)
        
        
        Defined at Math.fpp lines 132-158
        
        Parameters
        ----------
        str : str
        ch_comma : str
        nword : int
        
        """
        _AresMainPy.f90wrap_find_nword(str=str, ch_comma=ch_comma, nword=nword)
    
    @staticmethod
    def det(matrix):
        """
        det = det(matrix)
        
        
        Defined at Math.fpp lines 190-194
        
        Parameters
        ----------
        matrix : float array
        
        Returns
        -------
        det : float
        
        """
        det = _AresMainPy.f90wrap_det(matrix=matrix)
        return det
    
    @staticmethod
    def inv_33(m):
        """
        inv_33 = inv_33(m)
        
        
        Defined at Math.fpp lines 198-228
        
        Parameters
        ----------
        m : float array
        
        Returns
        -------
        inv_33 : float array
        
        """
        inv_33 = _AresMainPy.f90wrap_inv_33(m=m)
        return inv_33
    
    @staticmethod
    def lindg(eta, lambda_, mu):
        """
        lindg = lindg(eta, lambda_, mu)
        
        
        Defined at Math.fpp lines 232-278
        
        Parameters
        ----------
        eta : float
        lambda_ : float
        mu : float
        
        Returns
        -------
        lindg : float
        
        """
        lindg = _AresMainPy.f90wrap_lindg(eta=eta, lambda_=lambda_, mu=mu)
        return lindg
    
    @staticmethod
    def int_to_char(int_bn):
        """
        int_to_char = int_to_char(int_bn)
        
        
        Defined at Math.fpp lines 282-303
        
        Parameters
        ----------
        int_bn : int
        
        Returns
        -------
        int_to_char : str
        
        -----------------------------------------------------------------------
        """
        int_to_char = _AresMainPy.f90wrap_int_to_char(int_bn=int_bn)
        return int_to_char
    
    @staticmethod
    def lat2matrix(lat_para, lat_mat, flag):
        """
        lat2matrix(lat_para, lat_mat, flag)
        
        
        Defined at Math.fpp lines 307-350
        
        Parameters
        ----------
        lat_para : float array
        lat_mat : float array
        flag : int
        
        """
        _AresMainPy.f90wrap_lat2matrix(lat_para=lat_para, lat_mat=lat_mat, flag=flag)
    
    @staticmethod
    def one2three(id, n_dens, pos):
        """
        one2three(id, n_dens, pos)
        
        
        Defined at Math.fpp lines 354-376
        
        Parameters
        ----------
        id : int
        n_dens : int array
        pos : int array
        
        """
        _AresMainPy.f90wrap_one2three(id=id, n_dens=n_dens, pos=pos)
    
    @staticmethod
    def init_random_seed():
        """
        init_random_seed()
        
        
        Defined at Math.fpp lines 470-534
        
        
        """
        _AresMainPy.f90wrap_init_random_seed()
    
    @staticmethod
    def atom_mass(atom_name, mass):
        """
        atom_mass(atom_name, mass)
        
        
        Defined at Math.fpp lines 694-835
        
        Parameters
        ----------
        atom_name : str
        mass : float
        
        --------------------------------------------
        """
        _AresMainPy.f90wrap_atom_mass(atom_name=atom_name, mass=mass)
    
    @staticmethod
    def newton_inter(n, x, y, m, tx, ty):
        """
        newton_inter(n, x, y, m, tx, ty)
        
        
        Defined at Math.fpp lines 839-862
        
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
        _AresMainPy.f90wrap_newton_inter(n=n, x=x, y=y, m=m, tx=tx, ty=ty)
    
    @staticmethod
    def diag(n, ina, w, q):
        """
        diag(n, ina, w, q)
        
        
        Defined at Math.fpp lines 866-1000
        
        Parameters
        ----------
        n : int
        ina : float array
        w : float array
        q : float array
        
        """
        _AresMainPy.f90wrap_diag(n=n, ina=ina, w=w, q=q)
    
    @staticmethod
    def boltzmann_distribution(rnull, width):
        """
        boltzmann_distribution = boltzmann_distribution(rnull, width)
        
        
        Defined at Math.fpp lines 1004-1019
        
        Parameters
        ----------
        rnull : float
        width : float
        
        Returns
        -------
        boltzmann_distribution : float
        
        """
        boltzmann_distribution = _AresMainPy.f90wrap_boltzmann_distribution(rnull=rnull, \
            width=width)
        return boltzmann_distribution
    
    @staticmethod
    def dir2car(cry_coo, ort_coo, lat):
        """
        dir2car(cry_coo, ort_coo, lat)
        
        
        Defined at Math.fpp lines 1026-1036
        
        Parameters
        ----------
        cry_coo : float array
        ort_coo : float array
        lat : float array
        
        """
        _AresMainPy.f90wrap_dir2car(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)
    
    @staticmethod
    def car2dir(ort_coo, cry_coo, lat):
        """
        car2dir(ort_coo, cry_coo, lat)
        
        
        Defined at Math.fpp lines 1042-1053
        
        Parameters
        ----------
        ort_coo : float array
        cry_coo : float array
        lat : float array
        
        """
        _AresMainPy.f90wrap_car2dir(ort_coo=ort_coo, cry_coo=cry_coo, lat=lat)
    
    @staticmethod
    def thr2mat(n1, n2, n3, i, j, k):
        """
        dimnu = thr2mat(n1, n2, n3, i, j, k)
        
        
        Defined at Math.fpp lines 1056-1062
        
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
        dimnu = _AresMainPy.f90wrap_thr2mat(n1=n1, n2=n2, n3=n3, i=i, j=j, k=k)
        return dimnu
    
    @staticmethod
    def mat2thr(n1, n2, n3, i):
        """
        ix, iy, iz = mat2thr(n1, n2, n3, i)
        
        
        Defined at Math.fpp lines 1065-1096
        
        Parameters
        ----------
        n1 : int
        n2 : int
        n3 : int
        i : int
        
        Returns
        -------
        ix : int
        iy : int
        iz : int
        
        ------------------------
        """
        ix, iy, iz = _AresMainPy.f90wrap_mat2thr(n1=n1, n2=n2, n3=n3, i=i)
        return ix, iy, iz
    
    @staticmethod
    def sopo(a, lda, n, b, ldb, m, w):
        """
        sopo(a, lda, n, b, ldb, m, w)
        
        
        Defined at Math.fpp lines 1099-1155
        
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
        _AresMainPy.f90wrap_sopo(a=a, lda=lda, n=n, b=b, ldb=ldb, m=m, w=w)
    
    @staticmethod
    def csort_eigen(nev, arr, brr):
        """
        csort_eigen(nev, arr, brr)
        
        
        Defined at Math.fpp lines 1158-1191
        
        Parameters
        ----------
        nev : int
        arr : float array
        brr : complex array
        
        """
        _AresMainPy.f90wrap_csort_eigen(nev=nev, arr=arr, brr=brr)
    
    @staticmethod
    def rsort_eigen(nev, arr, brr):
        """
        rsort_eigen(nev, arr, brr)
        
        
        Defined at Math.fpp lines 1194-1227
        
        Parameters
        ----------
        nev : int
        arr : float array
        brr : float array
        
        """
        _AresMainPy.f90wrap_rsort_eigen(nev=nev, arr=arr, brr=brr)
    
    @staticmethod
    def sort_eigval(n, arr):
        """
        sort_eigval(n, arr)
        
        
        Defined at Math.fpp lines 1230-1252
        
        Parameters
        ----------
        n : int
        arr : float array
        
        """
        _AresMainPy.f90wrap_sort_eigval(n=n, arr=arr)
    
    @staticmethod
    def realint_sort(nev, arr, brr, crr=None):
        """
        realint_sort(nev, arr, brr[, crr])
        
        
        Defined at Math.fpp lines 1255-1293
        
        Parameters
        ----------
        nev : int
        arr : float array
        brr : int array
        crr : float array
        
        """
        _AresMainPy.f90wrap_realint_sort(nev=nev, arr=arr, brr=brr, crr=crr)
    
    @staticmethod
    def pgfo(omat, n1, n2, n3, ix, iy, iz):
        """
        ex, ey, ez = pgfo(omat, n1, n2, n3, ix, iy, iz)
        
        
        Defined at Math.fpp lines 1297-1327
        
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
        ex, ey, ez = _AresMainPy.f90wrap_pgfo(omat=omat, n1=n1, n2=n2, n3=n3, ix=ix, \
            iy=iy, iz=iz)
        return ex, ey, ez
    
    @staticmethod
    def cubicsplineinterp(fun, ddfdx2, xmax, dx, x, zion=None):
        """
        cubicsplineinterp = cubicsplineinterp(fun, ddfdx2, xmax, dx, x[, zion])
        
        
        Defined at Math.fpp lines 1330-1406
        
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
        cubicsplineinterp = _AresMainPy.f90wrap_cubicsplineinterp(fun=fun, \
            ddfdx2=ddfdx2, xmax=xmax, dx=dx, x=x, zion=zion)
        return cubicsplineinterp
    
    @staticmethod
    def finite_factor(fnor, norder, coe):
        """
        finite_factor(fnor, norder, coe)
        
        
        Defined at Math.fpp lines 1409-1550
        
        Parameters
        ----------
        fnor : int
        norder : int
        coe : float array
        
        """
        _AresMainPy.f90wrap_finite_factor(fnor=fnor, norder=norder, coe=coe)
    
    @staticmethod
    def finite_factor_new(fnor, norder, coe):
        """
        finite_factor_new(fnor, norder, coe)
        
        
        Defined at Math.fpp lines 1552-1608
        
        Parameters
        ----------
        fnor : int
        norder : int
        coe : float array
        
        """
        _AresMainPy.f90wrap_finite_factor_new(fnor=fnor, norder=norder, coe=coe)
    
    @staticmethod
    def dfdr(np, h, f, df):
        """
        dfdr(np, h, f, df)
        
        
        Defined at Math.fpp lines 1612-1646
        
        Parameters
        ----------
        np : int
        h : float
        f : float array
        df : float array
        
        """
        _AresMainPy.f90wrap_dfdr(np=np, h=h, f=f, df=df)
    
    @staticmethod
    def cubichermiteinterp(fun, dfdx, xmax, h, x):
        """
        cubichermiteinterp = cubichermiteinterp(fun, dfdx, xmax, h, x)
        
        
        Defined at Math.fpp lines 1649-1708
        
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
        cubichermiteinterp = _AresMainPy.f90wrap_cubichermiteinterp(fun=fun, dfdx=dfdx, \
            xmax=xmax, h=h, x=x)
        return cubichermiteinterp
    
    @staticmethod
    def simpleinterp(fun, xmax, h, x):
        """
        simpleinterp = simpleinterp(fun, xmax, h, x)
        
        
        Defined at Math.fpp lines 1711-1750
        
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
        simpleinterp = _AresMainPy.f90wrap_simpleinterp(fun=fun, xmax=xmax, h=h, x=x)
        return simpleinterp
    
    @staticmethod
    def r_dylm(l, m, x, y, z, rmod, f):
        """
        r_dylm(l, m, x, y, z, rmod, f)
        
        
        Defined at Math.fpp lines 1753-1893
        
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
        _AresMainPy.f90wrap_r_dylm(l=l, m=m, x=x, y=y, z=z, rmod=rmod, f=f)
    
    @staticmethod
    def atom_effcharge(atom_name, nquan, zeta):
        """
        lmax = atom_effcharge(atom_name, nquan, zeta)
        
        
        Defined at Math.fpp lines 2067-2161
        
        Parameters
        ----------
        atom_name : str
        nquan : int array
        zeta : float array
        
        Returns
        -------
        lmax : int
        
        --------------------------------------------
        """
        lmax = _AresMainPy.f90wrap_atom_effcharge(atom_name=atom_name, nquan=nquan, \
            zeta=zeta)
        return lmax
    
    @staticmethod
    def atom_sto(p, l, m, zeta, r):
        """
        f = atom_sto(p, l, m, zeta, r)
        
        
        Defined at Math.fpp lines 2446-2647
        
        Parameters
        ----------
        p : int
        l : int
        m : int
        zeta : float
        r : float array
        
        Returns
        -------
        f : float
        
        """
        f = _AresMainPy.f90wrap_atom_sto(p=p, l=l, m=m, zeta=zeta, r=r)
        return f
    
    @staticmethod
    def gammp(a, x):
        """
        gammp = gammp(a, x)
        
        
        Defined at Math.fpp lines 2805-2822
        
        Parameters
        ----------
        a : float
        x : float
        
        Returns
        -------
        gammp : float
        
        """
        gammp = _AresMainPy.f90wrap_gammp(a=a, x=x)
        return gammp
    
    @staticmethod
    def gser(gamser, a, x, gln):
        """
        gser(gamser, a, x, gln)
        
        
        Defined at Math.fpp lines 2825-2852
        
        Parameters
        ----------
        gamser : float
        a : float
        x : float
        gln : float
        
        """
        _AresMainPy.f90wrap_gser(gamser=gamser, a=a, x=x, gln=gln)
    
    @staticmethod
    def gcf(gammcf, a, x, gln):
        """
        gcf(gammcf, a, x, gln)
        
        
        Defined at Math.fpp lines 2855-2891
        
        Parameters
        ----------
        gammcf : float
        a : float
        x : float
        gln : float
        
        """
        _AresMainPy.f90wrap_gcf(gammcf=gammcf, a=a, x=x, gln=gln)
    
    @staticmethod
    def gammln(xx):
        """
        gammln = gammln(xx)
        
        
        Defined at Math.fpp lines 2894-2914
        
        Parameters
        ----------
        xx : float
        
        Returns
        -------
        gammln : float
        
        """
        gammln = _AresMainPy.f90wrap_gammln(xx=xx)
        return gammln
    
    @staticmethod
    def integral(l, dx, x):
        """
        integral = integral(l, dx, x)
        
        
        Defined at Math.fpp lines 2917-2930
        
        Parameters
        ----------
        l : int
        dx : float
        x : float
        
        Returns
        -------
        integral : float
        
        """
        integral = _AresMainPy.f90wrap_integral(l=l, dx=dx, x=x)
        return integral
    
    @staticmethod
    def myfun(l, t):
        """
        myfun = myfun(l, t)
        
        
        Defined at Math.fpp lines 2933-2938
        
        Parameters
        ----------
        l : int
        t : float
        
        Returns
        -------
        myfun : float
        
        """
        myfun = _AresMainPy.f90wrap_myfun(l=l, t=t)
        return myfun
    
    @staticmethod
    def plgndr(l, m, x):
        """
        plgndr = plgndr(l, m, x)
        
        
        Defined at Math.fpp lines 2941-2974
        
        Parameters
        ----------
        l : int
        m : int
        x : float
        
        Returns
        -------
        plgndr : float
        
        """
        plgndr = _AresMainPy.f90wrap_plgndr(l=l, m=m, x=x)
        return plgndr
    
    @staticmethod
    def lagrange_interpolation_coe(npoint, scatter_x, coe):
        """
        lagrange_interpolation_coe(npoint, scatter_x, coe)
        
        
        Defined at Math.fpp lines 2977-2994
        
        Parameters
        ----------
        npoint : int
        scatter_x : float array
        coe : float array
        
        """
        _AresMainPy.f90wrap_lagrange_interpolation_coe(npoint=npoint, \
            scatter_x=scatter_x, coe=coe)
    
    @staticmethod
    def lagrange_interpolation_x(npoint, x_sample, x_in):
        """
        lagrange_interpolation_x = lagrange_interpolation_x(npoint, x_sample, x_in)
        
        
        Defined at Math.fpp lines 2996-3011
        
        Parameters
        ----------
        npoint : int
        x_sample : float array
        x_in : float
        
        Returns
        -------
        lagrange_interpolation_x : float array
        
        """
        lagrange_interpolation_x = \
            _AresMainPy.f90wrap_lagrange_interpolation_x(npoint=npoint, \
            x_sample=x_sample, x_in=x_in)
        return lagrange_interpolation_x
    
    @staticmethod
    def interpolation_test():
        """
        interpolation_test()
        
        
        Defined at Math.fpp lines 3013-3037
        
        
        """
        _AresMainPy.f90wrap_interpolation_test()
    
    @staticmethod
    def direct_productlm(nll, nml, index_ll, index_ml, mat_in, mat_out):
        """
        direct_productlm(nll, nml, index_ll, index_ml, mat_in, mat_out)
        
        
        Defined at Math.fpp lines 3039-3116
        
        Parameters
        ----------
        nll : int
        nml : int
        index_ll : int array
        index_ml : int array
        mat_in : float array
        mat_out : float array
        
        """
        _AresMainPy.f90wrap_direct_productlm(nll=nll, nml=nml, index_ll=index_ll, \
            index_ml=index_ml, mat_in=mat_in, mat_out=mat_out)
    
    @staticmethod
    def fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt):
        """
        fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt)
        
        
        Defined at Math.fpp lines 3161-3219
        
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
        _AresMainPy.f90wrap_fourier_1d(nr=nr, rr=rr, rab=rab, vr=vr, ll=ll, nql=nql, \
            yp=yp, vql=vql, vt=vt)
    
    @staticmethod
    def invfourier_1d(g, fg, ll, r, fr):
        """
        invfourier_1d(g, fg, ll, r, fr)
        
        
        Defined at Math.fpp lines 3222-3253
        
        Parameters
        ----------
        g : float array
        fg : float array
        ll : int
        r : float array
        fr : float array
        
        """
        _AresMainPy.f90wrap_invfourier_1d(g=g, fg=fg, ll=ll, r=r, fr=fr)
    
    @staticmethod
    def integ_new(rab, y):
        """
        f = integ_new(rab, y)
        
        
        Defined at Math.fpp lines 3256-3276
        
        Parameters
        ----------
        rab : float array
        y : float array
        
        Returns
        -------
        f : float
        
        """
        f = _AresMainPy.f90wrap_integ_new(rab=rab, y=y)
        return f
    
    @staticmethod
    def interp(np, f, r, rnorm, z=None):
        """
        interp = interp(np, f, r, rnorm[, z])
        
        
        Defined at Math.fpp lines 3279-3333
        
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
        interp = _AresMainPy.f90wrap_interp(np=np, f=f, r=r, rnorm=rnorm, z=z)
        return interp
    
    @staticmethod
    def getfileunit():
        """
        getfileunit = getfileunit()
        
        
        Defined at Math.fpp lines 3408-3421
        
        
        Returns
        -------
        getfileunit : int
        
        """
        getfileunit = _AresMainPy.f90wrap_getfileunit()
        return getfileunit
    
    @staticmethod
    def kahan_sum(n, array):
        """
        kahan_sum = kahan_sum(n, array)
        
        
        Defined at Math.fpp lines 3424-3439
        
        Parameters
        ----------
        n : int
        array : float array
        
        Returns
        -------
        kahan_sum : float
        
        """
        kahan_sum = _AresMainPy.f90wrap_kahan_sum(n=n, array=array)
        return kahan_sum
    
    @staticmethod
    def _gasdev_s_sp(inmu=None, insigma=None):
        """
        harvest = _gasdev_s_sp([inmu, insigma])
        
        
        Defined at Math.fpp lines 538-567
        
        Parameters
        ----------
        inmu : float
        insigma : float
        
        Returns
        -------
        harvest : float
        
        """
        harvest = _AresMainPy.f90wrap_gasdev_s_sp(inmu=inmu, insigma=insigma)
        return harvest
    
    @staticmethod
    def _gasdev_s_dp(inmu=None, insigma=None):
        """
        harvest = _gasdev_s_dp([inmu, insigma])
        
        
        Defined at Math.fpp lines 571-600
        
        Parameters
        ----------
        inmu : float
        insigma : float
        
        Returns
        -------
        harvest : float
        
        """
        harvest = _AresMainPy.f90wrap_gasdev_s_dp(inmu=inmu, insigma=insigma)
        return harvest
    
    @staticmethod
    def _gasdev_v_sp(harvest, inmu=None, insigma=None):
        """
        _gasdev_v_sp(harvest[, inmu, insigma])
        
        
        Defined at Math.fpp lines 604-645
        
        Parameters
        ----------
        harvest : float array
        inmu : float
        insigma : float
        
        """
        _AresMainPy.f90wrap_gasdev_v_sp(harvest=harvest, inmu=inmu, insigma=insigma)
    
    @staticmethod
    def _gasdev_v_dp(harvest, inmu=None, insigma=None):
        """
        _gasdev_v_dp(harvest[, inmu, insigma])
        
        
        Defined at Math.fpp lines 649-690
        
        Parameters
        ----------
        harvest : float array
        inmu : float
        insigma : float
        
        """
        _AresMainPy.f90wrap_gasdev_v_dp(harvest=harvest, inmu=inmu, insigma=insigma)
    
    @staticmethod
    def gasdev(*args, **kwargs):
        """
        gasdev(*args, **kwargs)
        
        
        Defined at Math.fpp lines 21-22
        
        Overloaded interface containing the following procedures:
          _gasdev_s_sp
          _gasdev_s_dp
          _gasdev_v_sp
          _gasdev_v_dp
        
        """
        for proc in [Math._gasdev_s_sp, Math._gasdev_s_dp, Math._gasdev_v_sp, \
            Math._gasdev_v_dp]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _norm_real(a):
        """
        norm_real = _norm_real(a)
        
        
        Defined at Math.fpp lines 162-164
        
        Parameters
        ----------
        a : float array
        
        Returns
        -------
        norm_real : float
        
        """
        norm_real = _AresMainPy.f90wrap_norm_real(a=a)
        return norm_real
    
    @staticmethod
    def _norm_complex(a):
        """
        norm_complex = _norm_complex(a)
        
        
        Defined at Math.fpp lines 168-170
        
        Parameters
        ----------
        a : complex array
        
        Returns
        -------
        norm_complex : float
        
        """
        norm_complex = _AresMainPy.f90wrap_norm_complex(a=a)
        return norm_complex
    
    @staticmethod
    def norm(*args, **kwargs):
        """
        norm(*args, **kwargs)
        
        
        Defined at Math.fpp lines 26-27
        
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
        
        
        Defined at Math.fpp lines 182-186
        
        Parameters
        ----------
        a : float array
        b : float array
        
        Returns
        -------
        cross_real : float array
        
        """
        cross_real = _AresMainPy.f90wrap_cross_real(a=a, b=b)
        return cross_real
    
    @staticmethod
    def _cross_complex(a, b):
        """
        cross_complex = _cross_complex(a, b)
        
        
        Defined at Math.fpp lines 174-178
        
        Parameters
        ----------
        a : complex array
        b : complex array
        
        Returns
        -------
        cross_complex : complex array
        
        """
        cross_complex = _AresMainPy.f90wrap_cross_complex(a=a, b=b)
        return cross_complex
    
    @staticmethod
    def cross(*args, **kwargs):
        """
        cross(*args, **kwargs)
        
        
        Defined at Math.fpp lines 31-32
        
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
        
        
        Defined at Math.fpp lines 425-466
        
        Parameters
        ----------
        array : int array
        n : int
        id : int array
        
        """
        _AresMainPy.f90wrap_integer_index(array=array, n=n, id=id)
    
    @staticmethod
    def _real_index(array, n, id):
        """
        _real_index(array, n, id)
        
        
        Defined at Math.fpp lines 380-421
        
        Parameters
        ----------
        array : float array
        n : int
        id : int array
        
        """
        _AresMainPy.f90wrap_real_index(array=array, n=n, id=id)
    
    @staticmethod
    def sort_id(*args, **kwargs):
        """
        sort_id(*args, **kwargs)
        
        
        Defined at Math.fpp lines 36-37
        
        Overloaded interface containing the following procedures:
          _integer_index
          _real_index
        
        """
        for proc in [Math._integer_index, Math._real_index]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _three2one_d_cplx(amat, bmat):
        """
        _three2one_d_cplx(amat, bmat)
        
        
        Defined at Math.fpp lines 3337-3351
        
        Parameters
        ----------
        amat : complex array
        bmat : complex array
        
        -------------------------------------------------------------
        allocate(bmat(grid%n))
        """
        _AresMainPy.f90wrap_three2one_d_cplx(amat=amat, bmat=bmat)
    
    @staticmethod
    def _three2one_d_real(amat, bmat):
        """
        _three2one_d_real(amat, bmat)
        
        
        Defined at Math.fpp lines 3355-3369
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        
        -------------------------------------------------------------
        allocate(bmat(grid%n))
        """
        _AresMainPy.f90wrap_three2one_d_real(amat=amat, bmat=bmat)
    
    @staticmethod
    def three2one_dim(*args, **kwargs):
        """
        three2one_dim(*args, **kwargs)
        
        
        Defined at Math.fpp lines 42-44
        
        Overloaded interface containing the following procedures:
          _three2one_d_cplx
          _three2one_d_real
        
        """
        for proc in [Math._three2one_d_cplx, Math._three2one_d_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _one2three_d_cplx(amat, bmat):
        """
        _one2three_d_cplx(amat, bmat)
        
        
        Defined at Math.fpp lines 3373-3387
        
        Parameters
        ----------
        amat : complex array
        bmat : complex array
        
        -------------------------------------------------------------
        allocate(bmat(grid%n1,grid%n2,grid%n3))
        """
        _AresMainPy.f90wrap_one2three_d_cplx(amat=amat, bmat=bmat)
    
    @staticmethod
    def _one2three_d_real(amat, bmat):
        """
        _one2three_d_real(amat, bmat)
        
        
        Defined at Math.fpp lines 3391-3405
        
        Parameters
        ----------
        amat : float array
        bmat : float array
        
        -------------------------------------------------------------
        allocate(bmat(grid%n1,grid%n2,grid%n3))
        """
        _AresMainPy.f90wrap_one2three_d_real(amat=amat, bmat=bmat)
    
    @staticmethod
    def one2three_dim(*args, **kwargs):
        """
        one2three_dim(*args, **kwargs)
        
        
        Defined at Math.fpp lines 46-47
        
        Overloaded interface containing the following procedures:
          _one2three_d_cplx
          _one2three_d_real
        
        """
        for proc in [Math._one2three_d_cplx, Math._one2three_d_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    _dt_array_initialisers = []
    

math = Math()

class Fourier(f90wrap.runtime.FortranModule):
    """
    Module fourier
    
    
    Defined at Fourier.fpp lines 5-497
    
    """
    @staticmethod
    def planfft(dimx, dimy, dimz):
        """
        planfft(dimx, dimy, dimz)
        
        
        Defined at Fourier.fpp lines 73-147
        
        Parameters
        ----------
        dimx : int
        dimy : int
        dimz : int
        
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
        _AresMainPy.f90wrap_planfft(dimx=dimx, dimy=dimy, dimz=dimz)
    
    @staticmethod
    def planfst(dimx, dimy, dimz):
        """
        planfst(dimx, dimy, dimz)
        
        
        Defined at Fourier.fpp lines 149-181
        
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
        _AresMainPy.f90wrap_planfst(dimx=dimx, dimy=dimy, dimz=dimz)
    
    @staticmethod
    def getfftdims():
        """
        dimx, dimy, dimz = getfftdims()
        
        
        Defined at Fourier.fpp lines 183-210
        
        
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
        dimx, dimy, dimz = _AresMainPy.f90wrap_getfftdims()
        return dimx, dimy, dimz
    
    @staticmethod
    def getfftcomplexdims():
        """
        dimx, dimy, dimz = getfftcomplexdims()
        
        
        Defined at Fourier.fpp lines 212-239
        
        
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
        dimx, dimy, dimz = _AresMainPy.f90wrap_getfftcomplexdims()
        return dimx, dimy, dimz
    
    @staticmethod
    def forwardfst(array):
        """
        forwardfst(array)
        
        
        Defined at Fourier.fpp lines 410-435
        
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
        _AresMainPy.f90wrap_forwardfst(array=array)
    
    @staticmethod
    def backfst(array):
        """
        backfst(array)
        
        
        Defined at Fourier.fpp lines 437-460
        
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
        _AresMainPy.f90wrap_backfst(array=array)
    
    @staticmethod
    def cleanfft():
        """
        cleanfft()
        
        
        Defined at Fourier.fpp lines 462-497
        
        
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
        _AresMainPy.f90wrap_cleanfft()
    
    @staticmethod
    def _forwardfft_4d(array):
        """
        transform = _forwardfft_4d(array)
        
        
        Defined at Fourier.fpp lines 241-284
        
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
        transform = _AresMainPy.f90wrap_forwardfft_4d(array=array)
        return transform
    
    @staticmethod
    def _backfft_4d(array):
        """
        transform = _backfft_4d(array)
        
        
        Defined at Fourier.fpp lines 286-327
        
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
        transform = _AresMainPy.f90wrap_backfft_4d(array=array)
        return transform
    
    @staticmethod
    def _forwardfft_3d(array):
        """
        transform = _forwardfft_3d(array)
        
        
        Defined at Fourier.fpp lines 329-367
        
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
        transform = _AresMainPy.f90wrap_forwardfft_3d(array=array)
        return transform
    
    @staticmethod
    def _backfft_3d(array):
        """
        transform = _backfft_3d(array)
        
        
        Defined at Fourier.fpp lines 369-408
        
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
        transform = _AresMainPy.f90wrap_backfft_3d(array=array)
        return transform
    
    @staticmethod
    def fft(*args, **kwargs):
        """
        fft(*args, **kwargs)
        
        
        Defined at Fourier.fpp lines 66-70
        
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
        
        
        Defined at Fourier.fpp line 61
        
        """
        return _AresMainPy.f90wrap_fourier__get__offset()
    
    @offset.setter
    def offset(self, offset):
        _AresMainPy.f90wrap_fourier__set__offset(offset)
    
    def __str__(self):
        ret = ['<fourier>{\n']
        ret.append('    offset : ')
        ret.append(repr(self.offset))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

fourier = Fourier()

class Struct_Module(f90wrap.runtime.FortranModule):
    """
    Module struct_module
    
    
    Defined at Struct_module.fpp lines 5-106
    
    """
    @f90wrap.runtime.register_class("AresMainPy.struct_type")
    class struct_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=struct_type)
        
        
        Defined at Struct_module.fpp lines 13-30
        
        """
        def __init__(self, handle=None):
            """
            self = Struct_Type()
            
            
            Defined at Struct_module.fpp lines 13-30
            
            
            Returns
            -------
            this : Struct_Type
            	Object to be constructed
            
            
            Automatically generated constructor for struct_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_struct_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Struct_Type
            
            
            Defined at Struct_module.fpp lines 13-30
            
            Parameters
            ----------
            this : Struct_Type
            	Object to be destructed
            
            
            Automatically generated destructor for struct_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_struct_type_finalise(this=self._handle)
        
        @property
        def zion(self):
            """
            Element zion ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__zion(self._handle)
            if array_handle in self._arrays:
                zion = self._arrays[array_handle]
            else:
                zion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__zion)
                self._arrays[array_handle] = zion
            return zion
        
        @zion.setter
        def zion(self, zion):
            self.zion[...] = zion
        
        @property
        def nati(self):
            """
            Element nati ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.fpp line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__nati(self._handle)
            if array_handle in self._arrays:
                nati = self._arrays[array_handle]
            else:
                nati = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__nati)
                self._arrays[array_handle] = nati
            return nati
        
        @nati.setter
        def nati(self, nati):
            self.nati[...] = nati
        
        @property
        def eleid(self):
            """
            Element eleid ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.fpp line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__eleid(self._handle)
            if array_handle in self._arrays:
                eleid = self._arrays[array_handle]
            else:
                eleid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__eleid)
                self._arrays[array_handle] = eleid
            return eleid
        
        @eleid.setter
        def eleid(self, eleid):
            self.eleid[...] = eleid
        
        @property
        def pos(self):
            """
            Element pos ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__pos(self._handle)
            if array_handle in self._arrays:
                pos = self._arrays[array_handle]
            else:
                pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__pos)
                self._arrays[array_handle] = pos
            return pos
        
        @pos.setter
        def pos(self, pos):
            self.pos[...] = pos
        
        @property
        def poscar(self):
            """
            Element poscar ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__poscar(self._handle)
            if array_handle in self._arrays:
                poscar = self._arrays[array_handle]
            else:
                poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__poscar)
                self._arrays[array_handle] = poscar
            return poscar
        
        @poscar.setter
        def poscar(self, poscar):
            self.poscar[...] = poscar
        
        @property
        def stress(self):
            """
            Element stress ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__stress(self._handle)
            if array_handle in self._arrays:
                stress = self._arrays[array_handle]
            else:
                stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__stress)
                self._arrays[array_handle] = stress
            return stress
        
        @stress.setter
        def stress(self, stress):
            self.stress[...] = stress
        
        @property
        def forces(self):
            """
            Element forces ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__forces(self._handle)
            if array_handle in self._arrays:
                forces = self._arrays[array_handle]
            else:
                forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__forces)
                self._arrays[array_handle] = forces
            return forces
        
        @forces.setter
        def forces(self, forces):
            self.forces[...] = forces
        
        @property
        def mass(self):
            """
            Element mass ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__mass(self._handle)
            if array_handle in self._arrays:
                mass = self._arrays[array_handle]
            else:
                mass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__mass)
                self._arrays[array_handle] = mass
            return mass
        
        @mass.setter
        def mass(self, mass):
            self.mass[...] = mass
        
        @property
        def zeta(self):
            """
            Element zeta ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__zeta(self._handle)
            if array_handle in self._arrays:
                zeta = self._arrays[array_handle]
            else:
                zeta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__zeta)
                self._arrays[array_handle] = zeta
            return zeta
        
        @zeta.setter
        def zeta(self, zeta):
            self.zeta[...] = zeta
        
        @property
        def prinq(self):
            """
            Element prinq ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.fpp line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__prinq(self._handle)
            if array_handle in self._arrays:
                prinq = self._arrays[array_handle]
            else:
                prinq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__prinq)
                self._arrays[array_handle] = prinq
            return prinq
        
        @prinq.setter
        def prinq(self, prinq):
            self.prinq[...] = prinq
        
        @property
        def lmax(self):
            """
            Element lmax ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__lmax(self._handle)
            if array_handle in self._arrays:
                lmax = self._arrays[array_handle]
            else:
                lmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__lmax)
                self._arrays[array_handle] = lmax
            return lmax
        
        @lmax.setter
        def lmax(self, lmax):
            self.lmax[...] = lmax
        
        @property
        def elements(self):
            """
            Element elements ftype=character(len=3) pytype=str
            
            
            Defined at Struct_module.fpp line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__elements(self._handle)
            if array_handle in self._arrays:
                elements = self._arrays[array_handle]
            else:
                elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__elements)
                self._arrays[array_handle] = elements
            return elements
        
        @elements.setter
        def elements(self, elements):
            self.elements[...] = elements
        
        @property
        def coeff(self):
            """
            Element coeff ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__coeff(self._handle)
            if array_handle in self._arrays:
                coeff = self._arrays[array_handle]
            else:
                coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__coeff)
                self._arrays[array_handle] = coeff
            return coeff
        
        @coeff.setter
        def coeff(self, coeff):
            self.coeff[...] = coeff
        
        @property
        def occupy(self):
            """
            Element occupy ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_struct_type__array__occupy(self._handle)
            if array_handle in self._arrays:
                occupy = self._arrays[array_handle]
            else:
                occupy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_struct_type__array__occupy)
                self._arrays[array_handle] = occupy
            return occupy
        
        @occupy.setter
        def occupy(self, occupy):
            self.occupy[...] = occupy
        
        @property
        def noccupy(self):
            """
            Element noccupy ftype=integer(i4b) pytype=int
            
            
            Defined at Struct_module.fpp line 30
            
            """
            return _AresMainPy.f90wrap_struct_type__get__noccupy(self._handle)
        
        @noccupy.setter
        def noccupy(self, noccupy):
            _AresMainPy.f90wrap_struct_type__set__noccupy(self._handle, noccupy)
        
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
            ret.append(',\n    coeff : ')
            ret.append(repr(self.coeff))
            ret.append(',\n    occupy : ')
            ret.append(repr(self.occupy))
            ret.append(',\n    noccupy : ')
            ret.append(repr(self.noccupy))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def creat_struct(numtyp, numatom):
        """
        creat_struct(numtyp, numatom)
        
        
        Defined at Struct_module.fpp lines 47-76
        
        Parameters
        ----------
        numtyp : int
        numatom : int
        
        """
        _AresMainPy.f90wrap_creat_struct(numtyp=numtyp, numatom=numatom)
    
    @staticmethod
    def destroy_struct():
        """
        destroy_struct()
        
        
        Defined at Struct_module.fpp lines 79-104
        
        
        """
        _AresMainPy.f90wrap_destroy_struct()
    
    @property
    def natom(self):
        """
        Element natom ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 33
        
        """
        return _AresMainPy.f90wrap_struct_module__get__natom()
    
    @natom.setter
    def natom(self, natom):
        _AresMainPy.f90wrap_struct_module__set__natom(natom)
    
    @property
    def naty(self):
        """
        Element naty ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 34
        
        """
        return _AresMainPy.f90wrap_struct_module__get__naty()
    
    @naty.setter
    def naty(self, naty):
        _AresMainPy.f90wrap_struct_module__set__naty(naty)
    
    @property
    def ncharge(self):
        """
        Element ncharge ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 35
        
        """
        return _AresMainPy.f90wrap_struct_module__get__ncharge()
    
    @ncharge.setter
    def ncharge(self, ncharge):
        _AresMainPy.f90wrap_struct_module__set__ncharge(ncharge)
    
    @property
    def charge_ave(self):
        """
        Element charge_ave ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 36
        
        """
        return _AresMainPy.f90wrap_struct_module__get__charge_ave()
    
    @charge_ave.setter
    def charge_ave(self, charge_ave):
        _AresMainPy.f90wrap_struct_module__set__charge_ave(charge_ave)
    
    @property
    def volume(self):
        """
        Element volume ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 37
        
        """
        return _AresMainPy.f90wrap_struct_module__get__volume()
    
    @volume.setter
    def volume(self, volume):
        _AresMainPy.f90wrap_struct_module__set__volume(volume)
    
    @property
    def lat_mat(self):
        """
        Element lat_mat ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_struct_module__array__lat_mat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lat_mat = self._arrays[array_handle]
        else:
            lat_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_struct_module__array__lat_mat)
            self._arrays[array_handle] = lat_mat
        return lat_mat
    
    @lat_mat.setter
    def lat_mat(self, lat_mat):
        self.lat_mat[...] = lat_mat
    
    @property
    def lat_para(self):
        """
        Element lat_para ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_struct_module__array__lat_para(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lat_para = self._arrays[array_handle]
        else:
            lat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_struct_module__array__lat_para)
            self._arrays[array_handle] = lat_para
        return lat_para
    
    @lat_para.setter
    def lat_para(self, lat_para):
        self.lat_para[...] = lat_para
    
    @property
    def recip_lat(self):
        """
        Element recip_lat ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_struct_module__array__recip_lat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            recip_lat = self._arrays[array_handle]
        else:
            recip_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_struct_module__array__recip_lat)
            self._arrays[array_handle] = recip_lat
        return recip_lat
    
    @recip_lat.setter
    def recip_lat(self, recip_lat):
        self.recip_lat[...] = recip_lat
    
    @property
    def reclat_para(self):
        """
        Element reclat_para ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_struct_module__array__reclat_para(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            reclat_para = self._arrays[array_handle]
        else:
            reclat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_struct_module__array__reclat_para)
            self._arrays[array_handle] = reclat_para
        return reclat_para
    
    @reclat_para.setter
    def reclat_para(self, reclat_para):
        self.reclat_para[...] = reclat_para
    
    @property
    def energy(self):
        """
        Element energy ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_struct_module__array__energy(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            energy = self._arrays[array_handle]
        else:
            energy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_struct_module__array__energy)
            self._arrays[array_handle] = energy
        return energy
    
    @energy.setter
    def energy(self, energy):
        self.energy[...] = energy
    
    def __str__(self):
        ret = ['<struct_module>{\n']
        ret.append('    natom : ')
        ret.append(repr(self.natom))
        ret.append(',\n    naty : ')
        ret.append(repr(self.naty))
        ret.append(',\n    ncharge : ')
        ret.append(repr(self.ncharge))
        ret.append(',\n    charge_ave : ')
        ret.append(repr(self.charge_ave))
        ret.append(',\n    volume : ')
        ret.append(repr(self.volume))
        ret.append(',\n    lat_mat : ')
        ret.append(repr(self.lat_mat))
        ret.append(',\n    lat_para : ')
        ret.append(repr(self.lat_para))
        ret.append(',\n    recip_lat : ')
        ret.append(repr(self.recip_lat))
        ret.append(',\n    reclat_para : ')
        ret.append(repr(self.reclat_para))
        ret.append(',\n    energy : ')
        ret.append(repr(self.energy))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

struct_module = Struct_Module()

class Pspot_Module(f90wrap.runtime.FortranModule):
    """
    Module pspot_module
    
    
    Defined at Psp_module.fpp lines 10-48
    
    """
    @f90wrap.runtime.register_class("AresMainPy.pspot")
    class pspot(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=pspot)
        
        
        Defined at Psp_module.fpp lines 14-42
        
        """
        def __init__(self, handle=None):
            """
            self = Pspot()
            
            
            Defined at Psp_module.fpp lines 14-42
            
            
            Returns
            -------
            this : Pspot
            	Object to be constructed
            
            
            Automatically generated constructor for pspot
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_pspot_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Pspot
            
            
            Defined at Psp_module.fpp lines 14-42
            
            Parameters
            ----------
            this : Pspot
            	Object to be destructed
            
            
            Automatically generated destructor for pspot
            """
            if self._alloc:
                _AresMainPy.f90wrap_pspot_finalise(this=self._handle)
        
        @property
        def zion(self):
            """
            Element zion ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 15
            
            """
            return _AresMainPy.f90wrap_pspot__get__zion(self._handle)
        
        @zion.setter
        def zion(self, zion):
            _AresMainPy.f90wrap_pspot__set__zion(self._handle, zion)
        
        @property
        def numps(self):
            """
            Element numps ftype=integer(i4b) pytype=int
            
            
            Defined at Psp_module.fpp line 16
            
            """
            return _AresMainPy.f90wrap_pspot__get__numps(self._handle)
        
        @numps.setter
        def numps(self, numps):
            _AresMainPy.f90wrap_pspot__set__numps(self._handle, numps)
        
        @property
        def qnumps(self):
            """
            Element qnumps ftype=integer(i4b) pytype=int
            
            
            Defined at Psp_module.fpp line 17
            
            """
            return _AresMainPy.f90wrap_pspot__get__qnumps(self._handle)
        
        @qnumps.setter
        def qnumps(self, qnumps):
            _AresMainPy.f90wrap_pspot__set__qnumps(self._handle, qnumps)
        
        @property
        def numps_den(self):
            """
            Element numps_den ftype=integer(i4b) pytype=int
            
            
            Defined at Psp_module.fpp line 19
            
            """
            return _AresMainPy.f90wrap_pspot__get__numps_den(self._handle)
        
        @numps_den.setter
        def numps_den(self, numps_den):
            _AresMainPy.f90wrap_pspot__set__numps_den(self._handle, numps_den)
        
        @property
        def qmax(self):
            """
            Element qmax ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_pspot__get__qmax(self._handle)
        
        @qmax.setter
        def qmax(self, qmax):
            _AresMainPy.f90wrap_pspot__set__qmax(self._handle, qmax)
        
        @property
        def qspacing(self):
            """
            Element qspacing ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 22
            
            """
            return _AresMainPy.f90wrap_pspot__get__qspacing(self._handle)
        
        @qspacing.setter
        def qspacing(self, qspacing):
            _AresMainPy.f90wrap_pspot__set__qspacing(self._handle, qspacing)
        
        @property
        def qmesh(self):
            """
            Element qmesh ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__qmesh(self._handle)
            if array_handle in self._arrays:
                qmesh = self._arrays[array_handle]
            else:
                qmesh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__qmesh)
                self._arrays[array_handle] = qmesh
            return qmesh
        
        @qmesh.setter
        def qmesh(self, qmesh):
            self.qmesh[...] = qmesh
        
        @property
        def vlocq(self):
            """
            Element vlocq ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__vlocq(self._handle)
            if array_handle in self._arrays:
                vlocq = self._arrays[array_handle]
            else:
                vlocq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__vlocq)
                self._arrays[array_handle] = vlocq
            return vlocq
        
        @vlocq.setter
        def vlocq(self, vlocq):
            self.vlocq[...] = vlocq
        
        @property
        def vlocqs(self):
            """
            Element vlocqs ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__vlocqs(self._handle)
            if array_handle in self._arrays:
                vlocqs = self._arrays[array_handle]
            else:
                vlocqs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__vlocqs)
                self._arrays[array_handle] = vlocqs
            return vlocqs
        
        @vlocqs.setter
        def vlocqs(self, vlocqs):
            self.vlocqs[...] = vlocqs
        
        @property
        def ddvl_dq2(self):
            """
            Element ddvl_dq2 ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__ddvl_dq2(self._handle)
            if array_handle in self._arrays:
                ddvl_dq2 = self._arrays[array_handle]
            else:
                ddvl_dq2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__ddvl_dq2)
                self._arrays[array_handle] = ddvl_dq2
            return ddvl_dq2
        
        @ddvl_dq2.setter
        def ddvl_dq2(self, ddvl_dq2):
            self.ddvl_dq2[...] = ddvl_dq2
        
        @property
        def nproj(self):
            """
            Element nproj ftype=integer(i4b) pytype=int
            
            
            Defined at Psp_module.fpp line 28
            
            """
            return _AresMainPy.f90wrap_pspot__get__nproj(self._handle)
        
        @nproj.setter
        def nproj(self, nproj):
            _AresMainPy.f90wrap_pspot__set__nproj(self._handle, nproj)
        
        @property
        def proj_l(self):
            """
            Element proj_l ftype=integer(i4b) pytype=int
            
            
            Defined at Psp_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__proj_l(self._handle)
            if array_handle in self._arrays:
                proj_l = self._arrays[array_handle]
            else:
                proj_l = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__proj_l)
                self._arrays[array_handle] = proj_l
            return proj_l
        
        @proj_l.setter
        def proj_l(self, proj_l):
            self.proj_l[...] = proj_l
        
        @property
        def proj_m(self):
            """
            Element proj_m ftype=integer(i4b) pytype=int
            
            
            Defined at Psp_module.fpp line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__proj_m(self._handle)
            if array_handle in self._arrays:
                proj_m = self._arrays[array_handle]
            else:
                proj_m = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__proj_m)
                self._arrays[array_handle] = proj_m
            return proj_m
        
        @proj_m.setter
        def proj_m(self, proj_m):
            self.proj_m[...] = proj_m
        
        @property
        def rcut(self):
            """
            Element rcut ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 31
            
            """
            return _AresMainPy.f90wrap_pspot__get__rcut(self._handle)
        
        @rcut.setter
        def rcut(self, rcut):
            _AresMainPy.f90wrap_pspot__set__rcut(self._handle, rcut)
        
        @property
        def rmax(self):
            """
            Element rmax ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 32
            
            """
            return _AresMainPy.f90wrap_pspot__get__rmax(self._handle)
        
        @rmax.setter
        def rmax(self, rmax):
            _AresMainPy.f90wrap_pspot__set__rmax(self._handle, rmax)
        
        @property
        def rspacing(self):
            """
            Element rspacing ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 33
            
            """
            return _AresMainPy.f90wrap_pspot__get__rspacing(self._handle)
        
        @rspacing.setter
        def rspacing(self, rspacing):
            _AresMainPy.f90wrap_pspot__set__rspacing(self._handle, rspacing)
        
        @property
        def d0(self):
            """
            Element d0 ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 34
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__d0(self._handle)
            if array_handle in self._arrays:
                d0 = self._arrays[array_handle]
            else:
                d0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__d0)
                self._arrays[array_handle] = d0
            return d0
        
        @d0.setter
        def d0(self, d0):
            self.d0[...] = d0
        
        @property
        def beta_r(self):
            """
            Element beta_r ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 35
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__beta_r(self._handle)
            if array_handle in self._arrays:
                beta_r = self._arrays[array_handle]
            else:
                beta_r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__beta_r)
                self._arrays[array_handle] = beta_r
            return beta_r
        
        @beta_r.setter
        def beta_r(self, beta_r):
            self.beta_r[...] = beta_r
        
        @property
        def dbeta_dr(self):
            """
            Element dbeta_dr ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 36
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__dbeta_dr(self._handle)
            if array_handle in self._arrays:
                dbeta_dr = self._arrays[array_handle]
            else:
                dbeta_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__dbeta_dr)
                self._arrays[array_handle] = dbeta_dr
            return dbeta_dr
        
        @dbeta_dr.setter
        def dbeta_dr(self, dbeta_dr):
            self.dbeta_dr[...] = dbeta_dr
        
        @property
        def denr(self):
            """
            Element denr ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 38
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__denr(self._handle)
            if array_handle in self._arrays:
                denr = self._arrays[array_handle]
            else:
                denr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__denr)
                self._arrays[array_handle] = denr
            return denr
        
        @denr.setter
        def denr(self, denr):
            self.denr[...] = denr
        
        @property
        def ddden_dr2(self):
            """
            Element ddden_dr2 ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__ddden_dr2(self._handle)
            if array_handle in self._arrays:
                ddden_dr2 = self._arrays[array_handle]
            else:
                ddden_dr2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__ddden_dr2)
                self._arrays[array_handle] = ddden_dr2
            return ddden_dr2
        
        @ddden_dr2.setter
        def ddden_dr2(self, ddden_dr2):
            self.ddden_dr2[...] = ddden_dr2
        
        @property
        def r_real(self):
            """
            Element r_real ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 41
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__r_real(self._handle)
            if array_handle in self._arrays:
                r_real = self._arrays[array_handle]
            else:
                r_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__r_real)
                self._arrays[array_handle] = r_real
            return r_real
        
        @r_real.setter
        def r_real(self, r_real):
            self.r_real[...] = r_real
        
        @property
        def v_loc(self):
            """
            Element v_loc ftype=real(dp) pytype=float
            
            
            Defined at Psp_module.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_pspot__array__v_loc(self._handle)
            if array_handle in self._arrays:
                v_loc = self._arrays[array_handle]
            else:
                v_loc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_pspot__array__v_loc)
                self._arrays[array_handle] = v_loc
            return v_loc
        
        @v_loc.setter
        def v_loc(self, v_loc):
            self.v_loc[...] = v_loc
        
        def __str__(self):
            ret = ['<pspot>{\n']
            ret.append('    zion : ')
            ret.append(repr(self.zion))
            ret.append(',\n    numps : ')
            ret.append(repr(self.numps))
            ret.append(',\n    qnumps : ')
            ret.append(repr(self.qnumps))
            ret.append(',\n    numps_den : ')
            ret.append(repr(self.numps_den))
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
            ret.append(',\n    nproj : ')
            ret.append(repr(self.nproj))
            ret.append(',\n    proj_l : ')
            ret.append(repr(self.proj_l))
            ret.append(',\n    proj_m : ')
            ret.append(repr(self.proj_m))
            ret.append(',\n    rcut : ')
            ret.append(repr(self.rcut))
            ret.append(',\n    rmax : ')
            ret.append(repr(self.rmax))
            ret.append(',\n    rspacing : ')
            ret.append(repr(self.rspacing))
            ret.append(',\n    d0 : ')
            ret.append(repr(self.d0))
            ret.append(',\n    beta_r : ')
            ret.append(repr(self.beta_r))
            ret.append(',\n    dbeta_dr : ')
            ret.append(repr(self.dbeta_dr))
            ret.append(',\n    denr : ')
            ret.append(repr(self.denr))
            ret.append(',\n    ddden_dr2 : ')
            ret.append(repr(self.ddden_dr2))
            ret.append(',\n    r_real : ')
            ret.append(repr(self.r_real))
            ret.append(',\n    v_loc : ')
            ret.append(repr(self.v_loc))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @property
    def max_nproj(self):
        """
        Element max_nproj ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 47
        
        """
        return _AresMainPy.f90wrap_pspot_module__get__max_nproj()
    
    @max_nproj.setter
    def max_nproj(self, max_nproj):
        _AresMainPy.f90wrap_pspot_module__set__max_nproj(max_nproj)
    
    @property
    def max_rcut(self):
        """
        Element max_rcut ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 48
        
        """
        return _AresMainPy.f90wrap_pspot_module__get__max_rcut()
    
    @max_rcut.setter
    def max_rcut(self, max_rcut):
        _AresMainPy.f90wrap_pspot_module__set__max_rcut(max_rcut)
    
    @property
    def tknots(self):
        """
        Element tknots ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_pspot_module__array__tknots(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tknots = self._arrays[array_handle]
        else:
            tknots = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_pspot_module__array__tknots)
            self._arrays[array_handle] = tknots
        return tknots
    
    @tknots.setter
    def tknots(self, tknots):
        self.tknots[...] = tknots
    
    def __str__(self):
        ret = ['<pspot_module>{\n']
        ret.append('    max_nproj : ')
        ret.append(repr(self.max_nproj))
        ret.append(',\n    max_rcut : ')
        ret.append(repr(self.max_rcut))
        ret.append(',\n    tknots : ')
        ret.append(repr(self.tknots))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

pspot_module = Pspot_Module()

class Read_Module(f90wrap.runtime.FortranModule):
    """
    Module read_module
    
    
    Defined at Read_module.fpp lines 5-2613
    
    """
    @f90wrap.runtime.register_class("AresMainPy.attribute")
    class attribute(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=attribute)
        
        
        Defined at Read_module.fpp lines 19-20
        
        """
        def __init__(self, handle=None):
            """
            self = Attribute()
            
            
            Defined at Read_module.fpp lines 19-20
            
            
            Returns
            -------
            this : Attribute
            	Object to be constructed
            
            
            Automatically generated constructor for attribute
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_attribute_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Attribute
            
            
            Defined at Read_module.fpp lines 19-20
            
            Parameters
            ----------
            this : Attribute
            	Object to be destructed
            
            
            Automatically generated destructor for attribute
            """
            if self._alloc:
                _AresMainPy.f90wrap_attribute_finalise(this=self._handle)
        
        @property
        def value(self):
            """
            Element value ftype=character(len=120) pytype=str
            
            
            Defined at Read_module.fpp line 20
            
            """
            return _AresMainPy.f90wrap_attribute__get__value(self._handle)
        
        @value.setter
        def value(self, value):
            _AresMainPy.f90wrap_attribute__set__value(self._handle, value)
        
        def __str__(self):
            ret = ['<attribute>{\n']
            ret.append('    value : ')
            ret.append(repr(self.value))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def read_file(infile):
        """
        read_file(infile)
        
        
        Defined at Read_module.fpp lines 26-703
        
        Parameters
        ----------
        infile : str
        
        -----------------------------------------------------------
        """
        _AresMainPy.f90wrap_read_file(infile=infile)
    
    @staticmethod
    def read_pos(nty, filename):
        """
        read_pos(nty, filename)
        
        
        Defined at Read_module.fpp lines 706-828
        
        Parameters
        ----------
        nty : int
        filename : str
        
        """
        _AresMainPy.f90wrap_read_pos(nty=nty, filename=filename)
    
    @staticmethod
    def resetlattice():
        """
        resetlattice()
        
        
        Defined at Read_module.fpp lines 831-851
        
        
        """
        _AresMainPy.f90wrap_resetlattice()
    
    @staticmethod
    def read_pspot_atom(ity, filename):
        """
        ps = read_pspot_atom(ity, filename)
        
        
        Defined at Read_module.fpp lines 854-1172
        
        Parameters
        ----------
        ity : int
        filename : str
        
        Returns
        -------
        ps : Pspot
        
        ======================read comment=========================
        """
        ps = _AresMainPy.f90wrap_read_pspot_atom(ity=ity, filename=filename)
        ps = f90wrap.runtime.lookup_class("AresMainPy.pspot").from_handle(ps, \
            alloc=True)
        return ps
    
    @staticmethod
    def read_pspot(nty, filenames):
        """
        read_pspot(nty, filenames)
        
        
        Defined at Read_module.fpp lines 1175-1456
        
        Parameters
        ----------
        nty : int
        filenames : str array
        
        """
        _AresMainPy.f90wrap_read_pspot(nty=nty, filenames=filenames)
    
    @staticmethod
    def read_realpot_atom(ity, filename):
        """
        ps = read_realpot_atom(ity, filename)
        
        
        Defined at Read_module.fpp lines 1459-1804
        
        Parameters
        ----------
        ity : int
        filename : str
        
        Returns
        -------
        ps : Pspot
        
        ======================read comment=========================
        """
        ps = _AresMainPy.f90wrap_read_realpot_atom(ity=ity, filename=filename)
        ps = f90wrap.runtime.lookup_class("AresMainPy.pspot").from_handle(ps, \
            alloc=True)
        return ps
    
    @staticmethod
    def exist_in(string1, string2):
        """
        exist_in = exist_in(string1, string2)
        
        
        Defined at Read_module.fpp lines 1808-1831
        
        Parameters
        ----------
        string1 : str
        string2 : str
        
        Returns
        -------
        exist_in : bool
        
        ====== ====== ======
        """
        exist_in = _AresMainPy.f90wrap_exist_in(string1=string1, string2=string2)
        return exist_in
    
    @staticmethod
    def exist_ibegin(string1, string2):
        """
        exist_ibegin = exist_ibegin(string1, string2)
        
        
        Defined at Read_module.fpp lines 1833-1856
        
        Parameters
        ----------
        string1 : str
        string2 : str
        
        Returns
        -------
        exist_ibegin : int
        
        ====== ====== ======
        """
        exist_ibegin = _AresMainPy.f90wrap_exist_ibegin(string1=string1, \
            string2=string2)
        return exist_ibegin
    
    @staticmethod
    def scan_head(file_unit, title, start_old):
        """
        scan_head(file_unit, title, start_old)
        
        
        Defined at Read_module.fpp lines 1859-1916
        
        Parameters
        ----------
        file_unit : int
        title : str
        start_old : bool
        
        ====== ====== ======
         print*,"read"//title
        """
        _AresMainPy.f90wrap_scan_head(file_unit=file_unit, title=title, \
            start_old=start_old)
    
    @staticmethod
    def scan_tail(file_unit, title):
        """
        scan_tail(file_unit, title)
        
        
        Defined at Read_module.fpp lines 1919-1942
        
        Parameters
        ----------
        file_unit : int
        title : str
        
        ====== ====== ======
        """
        _AresMainPy.f90wrap_scan_tail(file_unit=file_unit, title=title)
    
    @staticmethod
    def read_upf(ity, filename, ps):
        """
        read_upf(ity, filename, ps)
        
        
        Defined at Read_module.fpp lines 1945-2224
        
        Parameters
        ----------
        ity : int
        filename : str
        ps : Pspot
        
        -------->Search for Nonlocal potential
        """
        _AresMainPy.f90wrap_read_upf(ity=ity, filename=filename, ps=ps._handle)
    
    @staticmethod
    def read_pseudo_header(zion, mesh_size, nproj):
        """
        read_pseudo_header(zion, mesh_size, nproj)
        
        
        Defined at Read_module.fpp lines 2226-2289
        
        Parameters
        ----------
        zion : float
        mesh_size : int
        nproj : int
        
        """
        _AresMainPy.f90wrap_read_pseudo_header(zion=zion, mesh_size=mesh_size, \
            nproj=nproj)
    
    @staticmethod
    def read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l):
        """
        read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l)
        
        
        Defined at Read_module.fpp lines 2291-2323
        
        Parameters
        ----------
        unit_upf : int
        nl : int
        beta_r : float array
        d0 : float array
        rcut : float
        proj_l : int array
        
        """
        _AresMainPy.f90wrap_read_pseudo_nonlocal(unit_upf=unit_upf, nl=nl, \
            beta_r=beta_r, d0=d0, rcut=rcut, proj_l=proj_l)
    
    @staticmethod
    def get_mo_coefficient():
        """
        get_mo_coefficient()
        
        
        Defined at Read_module.fpp lines 2386-2493
        
        
        """
        _AresMainPy.f90wrap_get_mo_coefficient()
    
    @staticmethod
    def parse_headline(str, atom_sign, atom_id, atom_orbital):
        """
        parse_headline(str, atom_sign, atom_id, atom_orbital)
        
        
        Defined at Read_module.fpp lines 2496-2547
        
        Parameters
        ----------
        str : str
        atom_sign : str array
        atom_id : int array
        atom_orbital : str array
        
        """
        _AresMainPy.f90wrap_parse_headline(str=str, atom_sign=atom_sign, \
            atom_id=atom_id, atom_orbital=atom_orbital)
    
    @staticmethod
    def sign2lm(atom_orbital, l, m):
        """
        sign2lm(atom_orbital, l, m)
        
        
        Defined at Read_module.fpp lines 2550-2572
        
        Parameters
        ----------
        atom_orbital : str array
        l : int array
        m : int array
        
        """
        _AresMainPy.f90wrap_sign2lm(atom_orbital=atom_orbital, l=l, m=m)
    
    @staticmethod
    def destroy_moinit():
        """
        destroy_moinit()
        
        
        Defined at Read_module.fpp lines 2575-2580
        
        
        """
        _AresMainPy.f90wrap_destroy_moinit()
    
    @staticmethod
    def out_concar(filename):
        """
        out_concar(filename)
        
        
        Defined at Read_module.fpp lines 2583-2601
        
        Parameters
        ----------
        filename : str
        
        """
        _AresMainPy.f90wrap_out_concar(filename=filename)
    
    @staticmethod
    def read_chgcar():
        """
        read_chgcar()
        
        
        Defined at Read_module.fpp lines 2603-2612
        
        
        """
        _AresMainPy.f90wrap_read_chgcar()
    
    @staticmethod
    def _get_value_int(char_in, char_find, variable, find_flag):
        """
        _get_value_int(char_in, char_find, variable, find_flag)
        
        
        Defined at Read_module.fpp lines 2326-2353
        
        Parameters
        ----------
        char_in : str
        char_find : str
        variable : int
        find_flag : bool
        
        """
        _AresMainPy.f90wrap_get_value_int(char_in=char_in, char_find=char_find, \
            variable=variable, find_flag=find_flag)
    
    @staticmethod
    def _get_value_real(char_in, char_find, variable, find_flag):
        """
        _get_value_real(char_in, char_find, variable, find_flag)
        
        
        Defined at Read_module.fpp lines 2355-2382
        
        Parameters
        ----------
        char_in : str
        char_find : str
        variable : float
        find_flag : bool
        
        """
        _AresMainPy.f90wrap_get_value_real(char_in=char_in, char_find=char_find, \
            variable=variable, find_flag=find_flag)
    
    @staticmethod
    def get_value(*args, **kwargs):
        """
        get_value(*args, **kwargs)
        
        
        Defined at Read_module.fpp lines 16-17
        
        Overloaded interface containing the following procedures:
          _get_value_int
          _get_value_real
        
        """
        for proc in [Read_Module._get_value_int, Read_Module._get_value_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    _dt_array_initialisers = []
    

read_module = Read_Module()

class Succeed(f90wrap.runtime.FortranModule):
    """
    Module succeed
    
    
    Defined at Succeed_module.fpp lines 5-487
    
    """
    @staticmethod
    def init_succeed_rho_real(n_rho, n_r, n_s, nspin, dv):
        """
        init_succeed_rho_real(n_rho, n_r, n_s, nspin, dv)
        
        
        Defined at Succeed_module.fpp lines 31-70
        
        Parameters
        ----------
        n_rho : int
        n_r : int
        n_s : int
        nspin : int
        dv : float
        
        """
        _AresMainPy.f90wrap_init_succeed_rho_real(n_rho=n_rho, n_r=n_r, n_s=n_s, \
            nspin=nspin, dv=dv)
    
    @staticmethod
    def init_succeed_rho_cmplx(n_rho, n_r, n_s, n_k, nspin, dv):
        """
        init_succeed_rho_cmplx(n_rho, n_r, n_s, n_k, nspin, dv)
        
        
        Defined at Succeed_module.fpp lines 72-113
        
        Parameters
        ----------
        n_rho : int
        n_r : int
        n_s : int
        n_k : int
        nspin : int
        dv : float
        
        """
        _AresMainPy.f90wrap_init_succeed_rho_cmplx(n_rho=n_rho, n_r=n_r, n_s=n_s, \
            n_k=n_k, nspin=nspin, dv=dv)
    
    @staticmethod
    def destroy_succeed():
        """
        destroy_succeed()
        
        
        Defined at Succeed_module.fpp lines 115-124
        
        
        """
        _AresMainPy.f90wrap_destroy_succeed()
    
    @staticmethod
    def store_rho(n_rho, nspin, rho, dvol):
        """
        store_rho(n_rho, nspin, rho, dvol)
        
        
        Defined at Succeed_module.fpp lines 126-140
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        rho : float array
        dvol : float
        
        """
        _AresMainPy.f90wrap_store_rho(n_rho=n_rho, nspin=nspin, rho=rho, dvol=dvol)
    
    @staticmethod
    def store_rho_at(n_rho, nspin, rho_in):
        """
        store_rho_at(n_rho, nspin, rho_in)
        
        
        Defined at Succeed_module.fpp lines 142-147
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        rho_in : float array
        
        """
        _AresMainPy.f90wrap_store_rho_at(n_rho=n_rho, nspin=nspin, rho_in=rho_in)
    
    @staticmethod
    def store_r(nr, r):
        """
        store_r(nr, r)
        
        
        Defined at Succeed_module.fpp lines 149-158
        
        Parameters
        ----------
        nr : int
        r : float array
        
        """
        _AresMainPy.f90wrap_store_r(nr=nr, r=r)
    
    @staticmethod
    def get_rho(nr, r_new, nrho, nspin, rho_new, dvol):
        """
        get_rho(nr, r_new, nrho, nspin, rho_new, dvol)
        
        
        Defined at Succeed_module.fpp lines 186-267
        
        Parameters
        ----------
        nr : int
        r_new : float array
        nrho : int
        nspin : int
        rho_new : float array
        dvol : float
        
        """
        _AresMainPy.f90wrap_get_rho(nr=nr, r_new=r_new, nrho=nrho, nspin=nspin, \
            rho_new=rho_new, dvol=dvol)
    
    @staticmethod
    def cal_trans_phase(nr, nspin, r_new, n1xy, n2xy, n3xy, ng1, ng2, ng3, gvec, \
        trans_phase):
        """
        cal_trans_phase(nr, nspin, r_new, n1xy, n2xy, n3xy, ng1, ng2, ng3, gvec, \
            trans_phase)
        
        
        Defined at Succeed_module.fpp lines 323-375
        
        Parameters
        ----------
        nr : int
        nspin : int
        r_new : float array
        n1xy : int
        n2xy : int
        n3xy : int
        ng1 : int
        ng2 : int
        ng3 : int
        gvec : float array
        trans_phase : complex array
        
        """
        _AresMainPy.f90wrap_cal_trans_phase(nr=nr, nspin=nspin, r_new=r_new, n1xy=n1xy, \
            n2xy=n2xy, n3xy=n3xy, ng1=ng1, ng2=ng2, ng3=ng3, gvec=gvec, \
            trans_phase=trans_phase)
    
    @staticmethod
    def get_new_rho_psi(nr, r_new, nrho, n1xy, n2xy, n3xy, nspin, rho_new, n_s, \
        psi_new, gvec):
        """
        get_new_rho_psi(nr, r_new, nrho, n1xy, n2xy, n3xy, nspin, rho_new, n_s, psi_new, \
            gvec)
        
        
        Defined at Succeed_module.fpp lines 377-436
        
        Parameters
        ----------
        nr : int
        r_new : float array
        nrho : int
        n1xy : int
        n2xy : int
        n3xy : int
        nspin : int
        rho_new : float array
        n_s : int
        psi_new : float array
        gvec : float array
        
        """
        _AresMainPy.f90wrap_get_new_rho_psi(nr=nr, r_new=r_new, nrho=nrho, n1xy=n1xy, \
            n2xy=n2xy, n3xy=n3xy, nspin=nspin, rho_new=rho_new, n_s=n_s, \
            psi_new=psi_new, gvec=gvec)
    
    @staticmethod
    def store_rho_fft_trans(n_rho, nspin, rho):
        """
        store_rho_fft_trans(n_rho, nspin, rho)
        
        
        Defined at Succeed_module.fpp lines 438-446
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        rho : float array
        
        """
        _AresMainPy.f90wrap_store_rho_fft_trans(n_rho=n_rho, nspin=nspin, rho=rho)
    
    @staticmethod
    def store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2):
        """
        store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2)
        
        
        Defined at Succeed_module.fpp lines 448-463
        
        Parameters
        ----------
        n_rho : int
        nspin : int
        na : int
        rho_in : float array
        rho_in2 : float array
        
        """
        _AresMainPy.f90wrap_store_rho_at_fft_trans(n_rho=n_rho, nspin=nspin, na=na, \
            rho_in=rho_in, rho_in2=rho_in2)
    
    @staticmethod
    def store_r_fft_trans(nr, r):
        """
        store_r_fft_trans(nr, r)
        
        
        Defined at Succeed_module.fpp lines 465-473
        
        Parameters
        ----------
        nr : int
        r : float array
        
        """
        _AresMainPy.f90wrap_store_r_fft_trans(nr=nr, r=r)
    
    @staticmethod
    def store_psi_fft_trans(n_rho, n_s, nspin, psi):
        """
        store_psi_fft_trans(n_rho, n_s, nspin, psi)
        
        
        Defined at Succeed_module.fpp lines 475-487
        
        Parameters
        ----------
        n_rho : int
        n_s : int
        nspin : int
        psi : float array
        
        """
        _AresMainPy.f90wrap_store_psi_fft_trans(n_rho=n_rho, n_s=n_s, nspin=nspin, \
            psi=psi)
    
    @staticmethod
    def _store_psi_cmplx(n_rho, n_s, n_k, nspin, psi):
        """
        _store_psi_cmplx(n_rho, n_s, n_k, nspin, psi)
        
        
        Defined at Succeed_module.fpp lines 173-184
        
        Parameters
        ----------
        n_rho : int
        n_s : int
        n_k : int
        nspin : int
        psi : complex array
        
        """
        _AresMainPy.f90wrap_store_psi_cmplx(n_rho=n_rho, n_s=n_s, n_k=n_k, nspin=nspin, \
            psi=psi)
    
    @staticmethod
    def _store_psi_real(n_rho, n_s, nspin, psi):
        """
        _store_psi_real(n_rho, n_s, nspin, psi)
        
        
        Defined at Succeed_module.fpp lines 160-171
        
        Parameters
        ----------
        n_rho : int
        n_s : int
        nspin : int
        psi : float array
        
        """
        _AresMainPy.f90wrap_store_psi_real(n_rho=n_rho, n_s=n_s, nspin=nspin, psi=psi)
    
    @staticmethod
    def store_psi(*args, **kwargs):
        """
        store_psi(*args, **kwargs)
        
        
        Defined at Succeed_module.fpp lines 22-24
        
        Overloaded interface containing the following procedures:
          _store_psi_cmplx
          _store_psi_real
        
        """
        for proc in [Succeed._store_psi_cmplx, Succeed._store_psi_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _get_psi_cmplx(nrho, n_s, n_k, nspin, psi_new):
        """
        _get_psi_cmplx(nrho, n_s, n_k, nspin, psi_new)
        
        
        Defined at Succeed_module.fpp lines 292-313
        
        Parameters
        ----------
        nrho : int
        n_s : int
        n_k : int
        nspin : int
        psi_new : complex array
        
        """
        _AresMainPy.f90wrap_get_psi_cmplx(nrho=nrho, n_s=n_s, n_k=n_k, nspin=nspin, \
            psi_new=psi_new)
    
    @staticmethod
    def _get_psi_real(nrho, n_s, nspin, psi_new):
        """
        _get_psi_real(nrho, n_s, nspin, psi_new)
        
        
        Defined at Succeed_module.fpp lines 269-290
        
        Parameters
        ----------
        nrho : int
        n_s : int
        nspin : int
        psi_new : float array
        
        """
        _AresMainPy.f90wrap_get_psi_real(nrho=nrho, n_s=n_s, nspin=nspin, \
            psi_new=psi_new)
    
    @staticmethod
    def get_psi(*args, **kwargs):
        """
        get_psi(*args, **kwargs)
        
        
        Defined at Succeed_module.fpp lines 26-28
        
        Overloaded interface containing the following procedures:
          _get_psi_cmplx
          _get_psi_real
        
        """
        for proc in [Succeed._get_psi_cmplx, Succeed._get_psi_real]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def llastrho(self):
        """
        Element llastrho ftype=logical pytype=bool
        
        
        Defined at Succeed_module.fpp line 12
        
        """
        return _AresMainPy.f90wrap_succeed__get__llastrho()
    
    @llastrho.setter
    def llastrho(self, llastrho):
        _AresMainPy.f90wrap_succeed__set__llastrho(llastrho)
    
    @property
    def lsr(self):
        """
        Element lsr ftype=logical pytype=bool
        
        
        Defined at Succeed_module.fpp line 12
        
        """
        return _AresMainPy.f90wrap_succeed__get__lsr()
    
    @lsr.setter
    def lsr(self, lsr):
        _AresMainPy.f90wrap_succeed__set__lsr(lsr)
    
    @property
    def lsrho(self):
        """
        Element lsrho ftype=logical pytype=bool
        
        
        Defined at Succeed_module.fpp line 12
        
        """
        return _AresMainPy.f90wrap_succeed__get__lsrho()
    
    @lsrho.setter
    def lsrho(self, lsrho):
        _AresMainPy.f90wrap_succeed__set__lsrho(lsrho)
    
    @property
    def lspsi(self):
        """
        Element lspsi ftype=logical pytype=bool
        
        
        Defined at Succeed_module.fpp line 12
        
        """
        return _AresMainPy.f90wrap_succeed__get__lspsi()
    
    @lspsi.setter
    def lspsi(self, lspsi):
        _AresMainPy.f90wrap_succeed__set__lspsi(lspsi)
    
    @property
    def lsrho_at(self):
        """
        Element lsrho_at ftype=logical pytype=bool
        
        
        Defined at Succeed_module.fpp line 12
        
        """
        return _AresMainPy.f90wrap_succeed__get__lsrho_at()
    
    @lsrho_at.setter
    def lsrho_at(self, lsrho_at):
        _AresMainPy.f90wrap_succeed__set__lsrho_at(lsrho_at)
    
    @property
    def rho1(self):
        """
        Element rho1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rho1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho1 = self._arrays[array_handle]
        else:
            rho1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rho1)
            self._arrays[array_handle] = rho1
        return rho1
    
    @rho1.setter
    def rho1(self, rho1):
        self.rho1[...] = rho1
    
    @property
    def rho2(self):
        """
        Element rho2 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rho2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho2 = self._arrays[array_handle]
        else:
            rho2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rho2)
            self._arrays[array_handle] = rho2
        return rho2
    
    @rho2.setter
    def rho2(self, rho2):
        self.rho2[...] = rho2
    
    @property
    def rho3(self):
        """
        Element rho3 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rho3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho3 = self._arrays[array_handle]
        else:
            rho3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rho3)
            self._arrays[array_handle] = rho3
        return rho3
    
    @rho3.setter
    def rho3(self, rho3):
        self.rho3[...] = rho3
    
    @property
    def r1(self):
        """
        Element r1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__r1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r1 = self._arrays[array_handle]
        else:
            r1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__r1)
            self._arrays[array_handle] = r1
        return r1
    
    @r1.setter
    def r1(self, r1):
        self.r1[...] = r1
    
    @property
    def r2(self):
        """
        Element r2 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__r2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r2 = self._arrays[array_handle]
        else:
            r2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__r2)
            self._arrays[array_handle] = r2
        return r2
    
    @r2.setter
    def r2(self, r2):
        self.r2[...] = r2
    
    @property
    def r3(self):
        """
        Element r3 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__r3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            r3 = self._arrays[array_handle]
        else:
            r3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__r3)
            self._arrays[array_handle] = r3
        return r3
    
    @r3.setter
    def r3(self, r3):
        self.r3[...] = r3
    
    @property
    def psi1(self):
        """
        Element psi1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__psi1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi1 = self._arrays[array_handle]
        else:
            psi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__psi1)
            self._arrays[array_handle] = psi1
        return psi1
    
    @psi1.setter
    def psi1(self, psi1):
        self.psi1[...] = psi1
    
    @property
    def psi2(self):
        """
        Element psi2 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__psi2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi2 = self._arrays[array_handle]
        else:
            psi2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__psi2)
            self._arrays[array_handle] = psi2
        return psi2
    
    @psi2.setter
    def psi2(self, psi2):
        self.psi2[...] = psi2
    
    @property
    def psi3(self):
        """
        Element psi3 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__psi3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi3 = self._arrays[array_handle]
        else:
            psi3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__psi3)
            self._arrays[array_handle] = psi3
        return psi3
    
    @psi3.setter
    def psi3(self, psi3):
        self.psi3[...] = psi3
    
    @property
    def rho_at(self):
        """
        Element rho_at ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rho_at(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho_at = self._arrays[array_handle]
        else:
            rho_at = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rho_at)
            self._arrays[array_handle] = rho_at
        return rho_at
    
    @rho_at.setter
    def rho_at(self, rho_at):
        self.rho_at[...] = rho_at
    
    @property
    def rhoi(self):
        """
        Element rhoi ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rhoi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rhoi = self._arrays[array_handle]
        else:
            rhoi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rhoi)
            self._arrays[array_handle] = rhoi
        return rhoi
    
    @rhoi.setter
    def rhoi(self, rhoi):
        self.rhoi[...] = rhoi
    
    @property
    def rho_at1(self):
        """
        Element rho_at1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rho_at1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rho_at1 = self._arrays[array_handle]
        else:
            rho_at1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rho_at1)
            self._arrays[array_handle] = rho_at1
        return rho_at1
    
    @rho_at1.setter
    def rho_at1(self, rho_at1):
        self.rho_at1[...] = rho_at1
    
    @property
    def rhoi1(self):
        """
        Element rhoi1 ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__rhoi1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rhoi1 = self._arrays[array_handle]
        else:
            rhoi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__rhoi1)
            self._arrays[array_handle] = rhoi1
        return rhoi1
    
    @rhoi1.setter
    def rhoi1(self, rhoi1):
        self.rhoi1[...] = rhoi1
    
    @property
    def psi1_c(self):
        """
        Element psi1_c ftype=complex(dcp) pytype=complex
        
        
        Defined at Succeed_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__psi1_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi1_c = self._arrays[array_handle]
        else:
            psi1_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__psi1_c)
            self._arrays[array_handle] = psi1_c
        return psi1_c
    
    @psi1_c.setter
    def psi1_c(self, psi1_c):
        self.psi1_c[...] = psi1_c
    
    @property
    def psi2_c(self):
        """
        Element psi2_c ftype=complex(dcp) pytype=complex
        
        
        Defined at Succeed_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__psi2_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi2_c = self._arrays[array_handle]
        else:
            psi2_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__psi2_c)
            self._arrays[array_handle] = psi2_c
        return psi2_c
    
    @psi2_c.setter
    def psi2_c(self, psi2_c):
        self.psi2_c[...] = psi2_c
    
    @property
    def psi3_c(self):
        """
        Element psi3_c ftype=complex(dcp) pytype=complex
        
        
        Defined at Succeed_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_succeed__array__psi3_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            psi3_c = self._arrays[array_handle]
        else:
            psi3_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_succeed__array__psi3_c)
            self._arrays[array_handle] = psi3_c
        return psi3_c
    
    @psi3_c.setter
    def psi3_c(self, psi3_c):
        self.psi3_c[...] = psi3_c
    
    @property
    def counter1(self):
        """
        Element counter1 ftype=integer(i4b) pytype=int
        
        
        Defined at Succeed_module.fpp line 20
        
        """
        return _AresMainPy.f90wrap_succeed__get__counter1()
    
    @counter1.setter
    def counter1(self, counter1):
        _AresMainPy.f90wrap_succeed__set__counter1(counter1)
    
    @property
    def counter2(self):
        """
        Element counter2 ftype=integer(i4b) pytype=int
        
        
        Defined at Succeed_module.fpp line 20
        
        """
        return _AresMainPy.f90wrap_succeed__get__counter2()
    
    @counter2.setter
    def counter2(self, counter2):
        _AresMainPy.f90wrap_succeed__set__counter2(counter2)
    
    @property
    def counter3(self):
        """
        Element counter3 ftype=integer(i4b) pytype=int
        
        
        Defined at Succeed_module.fpp line 20
        
        """
        return _AresMainPy.f90wrap_succeed__get__counter3()
    
    @counter3.setter
    def counter3(self, counter3):
        _AresMainPy.f90wrap_succeed__set__counter3(counter3)
    
    @property
    def alpha(self):
        """
        Element alpha ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 21
        
        """
        return _AresMainPy.f90wrap_succeed__get__alpha()
    
    @alpha.setter
    def alpha(self, alpha):
        _AresMainPy.f90wrap_succeed__set__alpha(alpha)
    
    @property
    def beta(self):
        """
        Element beta ftype=real(dp) pytype=float
        
        
        Defined at Succeed_module.fpp line 21
        
        """
        return _AresMainPy.f90wrap_succeed__get__beta()
    
    @beta.setter
    def beta(self, beta):
        _AresMainPy.f90wrap_succeed__set__beta(beta)
    
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
        ret.append(',\n    psi1_c : ')
        ret.append(repr(self.psi1_c))
        ret.append(',\n    psi2_c : ')
        ret.append(repr(self.psi2_c))
        ret.append(',\n    psi3_c : ')
        ret.append(repr(self.psi3_c))
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

class Scalapack_Module(f90wrap.runtime.FortranModule):
    """
    Module scalapack_module
    
    
    Defined at Scala.fpp lines 5-1150
    
    """
    @staticmethod
    def init_scala():
        """
        init_scala()
        
        
        Defined at Scala.fpp lines 29-47
        
        
        """
        _AresMainPy.f90wrap_init_scala()
    
    @staticmethod
    def sl_orthnorm(amat, m, n, mb, nb):
        """
        sl_orthnorm(amat, m, n, mb, nb)
        
        
        Defined at Scala.fpp lines 59-78
        
        Parameters
        ----------
        amat : complex array
        m : int
        n : int
        mb : int
        nb : int
        
        """
        _AresMainPy.f90wrap_sl_orthnorm(amat=amat, m=m, n=n, mb=mb, nb=nb)
    
    @staticmethod
    def sl_orthnorm_real(amat, m, n, mb, nb):
        """
        sl_orthnorm_real(amat, m, n, mb, nb)
        
        
        Defined at Scala.fpp lines 81-117
        
        Parameters
        ----------
        amat : float array
        m : int
        n : int
        mb : int
        nb : int
        
        """
        _AresMainPy.f90wrap_sl_orthnorm_real(amat=amat, m=m, n=n, mb=mb, nb=nb)
    
    @staticmethod
    def sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        cmbin=None, cnbin=None):
        """
        sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
            cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 120-176
        
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
        _AresMainPy.f90wrap_sl_matmat(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, \
            am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, cmbin=cmbin, \
            cnbin=cnbin)
    
    @staticmethod
    def sl_matmat_gridcn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb):
        """
        sl_matmat_gridcn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            cmb, cnb)
        
        
        Defined at Scala.fpp lines 179-288
        
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
        _AresMainPy.f90wrap_sl_matmat_gridcn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def sl_matmat_gridnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb):
        """
        sl_matmat_gridnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            cmb, cnb)
        
        
        Defined at Scala.fpp lines 291-350
        
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
        _AresMainPy.f90wrap_sl_matmat_gridnn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmbin=None, cnbin=None):
        """
        sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
            cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 353-438
        
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
        _AresMainPy.f90wrap_sl_matmat_sub(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
        bmb, bnb, cmbin=None, cnbin=None):
        """
        sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
            bnb[, cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 441-526
        
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
        _AresMainPy.f90wrap_sl_matmat_sub_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmbin=None, cnbin=None):
        """
        sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
            cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 529-626
        
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
        _AresMainPy.f90wrap_sl_matmat_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        evec, cm, cn, eval, cmbin=None, cnbin=None):
        """
        sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, \
            cm, cn, eval[, cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 629-712
        
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
        _AresMainPy.f90wrap_sl_generalizeeigen(dime=dime, amat=amat, bmat=bmat, am=am, \
            an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, cm=cm, \
            cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_generalizeeigen_real2(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, evec, cm, cn, eval, cmbin=None, cnbin=None):
        """
        sl_generalizeeigen_real2(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            evec, cm, cn, eval[, cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 715-802
        
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
        _AresMainPy.f90wrap_sl_generalizeeigen_real2(dime=dime, amat=amat, bmat=bmat, \
            am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, \
            cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, evec, cm, cn, eval, cmbin=None, cnbin=None):
        """
        sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            evec, cm, cn, eval[, cmbin, cnbin])
        
        
        Defined at Scala.fpp lines 805-954
        
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
        _AresMainPy.f90wrap_sl_generalizeeigen_real(dime=dime, amat=amat, bmat=bmat, \
            am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, \
            cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
    
    @staticmethod
    def twod_map_set(nstates, nrow, ncol, twod_map):
        """
        twod_map_set(nstates, nrow, ncol, twod_map)
        
        
        Defined at Scala.fpp lines 958-976
        
        Parameters
        ----------
        nstates : int
        nrow : int
        ncol : int
        twod_map : int array
        
        """
        _AresMainPy.f90wrap_twod_map_set(nstates=nstates, nrow=nrow, ncol=ncol, \
            twod_map=twod_map)
    
    @staticmethod
    def sl_matmat_realtn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb):
        """
        sl_matmat_realtn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            cmb, cnb)
        
        
        Defined at Scala.fpp lines 978-1087
        
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
        _AresMainPy.f90wrap_sl_matmat_realtn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @staticmethod
    def sl_matmat_realnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb):
        """
        sl_matmat_realnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
            cmb, cnb)
        
        
        Defined at Scala.fpp lines 1090-1149
        
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
        _AresMainPy.f90wrap_sl_matmat_realnn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
            cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
            cmb=cmb, cnb=cnb)
    
    @property
    def blacs_contxt(self):
        """
        Element blacs_contxt ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 15
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__blacs_contxt()
    
    @blacs_contxt.setter
    def blacs_contxt(self, blacs_contxt):
        _AresMainPy.f90wrap_scalapack_module__set__blacs_contxt(blacs_contxt)
    
    @property
    def dlen(self):
        """
        Element dlen ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 16
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__dlen()
    
    @property
    def myrow(self):
        """
        Element myrow ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 17
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__myrow()
    
    @myrow.setter
    def myrow(self, myrow):
        _AresMainPy.f90wrap_scalapack_module__set__myrow(myrow)
    
    @property
    def mycol(self):
        """
        Element mycol ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 17
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__mycol()
    
    @mycol.setter
    def mycol(self, mycol):
        _AresMainPy.f90wrap_scalapack_module__set__mycol(mycol)
    
    @property
    def npcol(self):
        """
        Element npcol ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 18
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__npcol()
    
    @npcol.setter
    def npcol(self, npcol):
        _AresMainPy.f90wrap_scalapack_module__set__npcol(npcol)
    
    @property
    def nprow(self):
        """
        Element nprow ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 18
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__nprow()
    
    @nprow.setter
    def nprow(self, nprow):
        _AresMainPy.f90wrap_scalapack_module__set__nprow(nprow)
    
    @property
    def my_blacs_id(self):
        """
        Element my_blacs_id ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 19
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__my_blacs_id()
    
    @my_blacs_id.setter
    def my_blacs_id(self, my_blacs_id):
        _AresMainPy.f90wrap_scalapack_module__set__my_blacs_id(my_blacs_id)
    
    @property
    def np(self):
        """
        Element np ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 20
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__np()
    
    @np.setter
    def np(self, np):
        _AresMainPy.f90wrap_scalapack_module__set__np(np)
    
    @property
    def nq(self):
        """
        Element nq ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 20
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__nq()
    
    @nq.setter
    def nq(self, nq):
        _AresMainPy.f90wrap_scalapack_module__set__nq(nq)
    
    @property
    def iam(self):
        """
        Element iam ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 21
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__iam()
    
    @iam.setter
    def iam(self, iam):
        _AresMainPy.f90wrap_scalapack_module__set__iam(iam)
    
    @property
    def nprocs(self):
        """
        Element nprocs ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 21
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__nprocs()
    
    @nprocs.setter
    def nprocs(self, nprocs):
        _AresMainPy.f90wrap_scalapack_module__set__nprocs(nprocs)
    
    @property
    def init_called(self):
        """
        Element init_called ftype=logical pytype=bool
        
        
        Defined at Scala.fpp line 22
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__init_called()
    
    @init_called.setter
    def init_called(self, init_called):
        _AresMainPy.f90wrap_scalapack_module__set__init_called(init_called)
    
    @property
    def info(self):
        """
        Element info ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 24
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__info()
    
    @info.setter
    def info(self, info):
        _AresMainPy.f90wrap_scalapack_module__set__info(info)
    
    @property
    def l_useless(self):
        """
        Element l_useless ftype=logical pytype=bool
        
        
        Defined at Scala.fpp line 25
        
        """
        return _AresMainPy.f90wrap_scalapack_module__get__l_useless()
    
    @l_useless.setter
    def l_useless(self, l_useless):
        _AresMainPy.f90wrap_scalapack_module__set__l_useless(l_useless)
    
    @property
    def twod_map(self):
        """
        Element twod_map ftype=integer(i4b) pytype=int
        
        
        Defined at Scala.fpp line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_scalapack_module__array__twod_map(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            twod_map = self._arrays[array_handle]
        else:
            twod_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_scalapack_module__array__twod_map)
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

class Grid_Module(f90wrap.runtime.FortranModule):
    """
    Module grid_module
    
    
    Defined at Grid_module.fpp lines 5-1248
    
    """
    @f90wrap.runtime.register_class("AresMainPy.grid_type")
    class grid_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=grid_type)
        
        
        Defined at Grid_module.fpp lines 14-23
        
        """
        def __init__(self, handle=None):
            """
            self = Grid_Type()
            
            
            Defined at Grid_module.fpp lines 14-23
            
            
            Returns
            -------
            this : Grid_Type
            	Object to be constructed
            
            
            Automatically generated constructor for grid_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_grid_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Grid_Type
            
            
            Defined at Grid_module.fpp lines 14-23
            
            Parameters
            ----------
            this : Grid_Type
            	Object to be destructed
            
            
            Automatically generated destructor for grid_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_grid_type_finalise(this=self._handle)
        
        @property
        def rhos(self):
            """
            Element rhos ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_type__array__rhos(self._handle)
            if array_handle in self._arrays:
                rhos = self._arrays[array_handle]
            else:
                rhos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_type__array__rhos)
                self._arrays[array_handle] = rhos
            return rhos
        
        @rhos.setter
        def rhos(self, rhos):
            self.rhos[...] = rhos
        
        @property
        def vlpp(self):
            """
            Element vlpp ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_type__array__vlpp(self._handle)
            if array_handle in self._arrays:
                vlpp = self._arrays[array_handle]
            else:
                vlpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_type__array__vlpp)
                self._arrays[array_handle] = vlpp
            return vlpp
        
        @vlpp.setter
        def vlpp(self, vlpp):
            self.vlpp[...] = vlpp
        
        @property
        def eval(self):
            """
            Element eval ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_type__array__eval(self._handle)
            if array_handle in self._arrays:
                eval = self._arrays[array_handle]
            else:
                eval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_type__array__eval)
                self._arrays[array_handle] = eval
            return eval
        
        @eval.setter
        def eval(self, eval):
            self.eval[...] = eval
        
        @property
        def gvec(self):
            """
            Element gvec ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_type__array__gvec(self._handle)
            if array_handle in self._arrays:
                gvec = self._arrays[array_handle]
            else:
                gvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_type__array__gvec)
                self._arrays[array_handle] = gvec
            return gvec
        
        @gvec.setter
        def gvec(self, gvec):
            self.gvec[...] = gvec
        
        @property
        def gmask(self):
            """
            Element gmask ftype=logical pytype=bool
            
            
            Defined at Grid_module.fpp line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_type__array__gmask(self._handle)
            if array_handle in self._arrays:
                gmask = self._arrays[array_handle]
            else:
                gmask = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_type__array__gmask)
                self._arrays[array_handle] = gmask
            return gmask
        
        @gmask.setter
        def gmask(self, gmask):
            self.gmask[...] = gmask
        
        @property
        def rvec(self):
            """
            Element rvec ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_grid_type__array__rvec(self._handle)
            if array_handle in self._arrays:
                rvec = self._arrays[array_handle]
            else:
                rvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_grid_type__array__rvec)
                self._arrays[array_handle] = rvec
            return rvec
        
        @rvec.setter
        def rvec(self, rvec):
            self.rvec[...] = rvec
        
        def __str__(self):
            ret = ['<grid_type>{\n']
            ret.append('    rhos : ')
            ret.append(repr(self.rhos))
            ret.append(',\n    vlpp : ')
            ret.append(repr(self.vlpp))
            ret.append(',\n    eval : ')
            ret.append(repr(self.eval))
            ret.append(',\n    gvec : ')
            ret.append(repr(self.gvec))
            ret.append(',\n    gmask : ')
            ret.append(repr(self.gmask))
            ret.append(',\n    rvec : ')
            ret.append(repr(self.rvec))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("AresMainPy.kgrid_type")
    class kgrid_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=kgrid_type)
        
        
        Defined at Grid_module.fpp lines 26-30
        
        """
        def __init__(self, handle=None):
            """
            self = Kgrid_Type()
            
            
            Defined at Grid_module.fpp lines 26-30
            
            
            Returns
            -------
            this : Kgrid_Type
            	Object to be constructed
            
            
            Automatically generated constructor for kgrid_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_kgrid_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Kgrid_Type
            
            
            Defined at Grid_module.fpp lines 26-30
            
            Parameters
            ----------
            this : Kgrid_Type
            	Object to be destructed
            
            
            Automatically generated destructor for kgrid_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_kgrid_type_finalise(this=self._handle)
        
        @property
        def vec(self):
            """
            Element vec ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_kgrid_type__array__vec(self._handle)
            if array_handle in self._arrays:
                vec = self._arrays[array_handle]
            else:
                vec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_kgrid_type__array__vec)
                self._arrays[array_handle] = vec
            return vec
        
        @vec.setter
        def vec(self, vec):
            self.vec[...] = vec
        
        @property
        def vcar(self):
            """
            Element vcar ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_kgrid_type__array__vcar(self._handle)
            if array_handle in self._arrays:
                vcar = self._arrays[array_handle]
            else:
                vcar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_kgrid_type__array__vcar)
                self._arrays[array_handle] = vcar
            return vcar
        
        @vcar.setter
        def vcar(self, vcar):
            self.vcar[...] = vcar
        
        @property
        def wk(self):
            """
            Element wk ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_kgrid_type__array__wk(self._handle)
            if array_handle in self._arrays:
                wk = self._arrays[array_handle]
            else:
                wk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_kgrid_type__array__wk)
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
        
    
    @f90wrap.runtime.register_class("AresMainPy.eigen_type")
    class eigen_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=eigen_type)
        
        
        Defined at Grid_module.fpp lines 33-37
        
        """
        def __init__(self, handle=None):
            """
            self = Eigen_Type()
            
            
            Defined at Grid_module.fpp lines 33-37
            
            
            Returns
            -------
            this : Eigen_Type
            	Object to be constructed
            
            
            Automatically generated constructor for eigen_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_eigen_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Eigen_Type
            
            
            Defined at Grid_module.fpp lines 33-37
            
            Parameters
            ----------
            this : Eigen_Type
            	Object to be destructed
            
            
            Automatically generated destructor for eigen_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_eigen_type_finalise(this=self._handle)
        
        @property
        def wvf(self):
            """
            Element wvf ftype=complex(dp) pytype=complex
            
            
            Defined at Grid_module.fpp line 35
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_eigen_type__array__wvf(self._handle)
            if array_handle in self._arrays:
                wvf = self._arrays[array_handle]
            else:
                wvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_eigen_type__array__wvf)
                self._arrays[array_handle] = wvf
            return wvf
        
        @wvf.setter
        def wvf(self, wvf):
            self.wvf[...] = wvf
        
        @property
        def val(self):
            """
            Element val ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 37
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_eigen_type__array__val(self._handle)
            if array_handle in self._arrays:
                val = self._arrays[array_handle]
            else:
                val = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_eigen_type__array__val)
                self._arrays[array_handle] = val
            return val
        
        @val.setter
        def val(self, val):
            self.val[...] = val
        
        def __str__(self):
            ret = ['<eigen_type>{\n']
            ret.append('    wvf : ')
            ret.append(repr(self.wvf))
            ret.append(',\n    val : ')
            ret.append(repr(self.val))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("AresMainPy.eigen_type_r")
    class eigen_type_r(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=eigen_type_r)
        
        
        Defined at Grid_module.fpp lines 40-44
        
        """
        def __init__(self, handle=None):
            """
            self = Eigen_Type_R()
            
            
            Defined at Grid_module.fpp lines 40-44
            
            
            Returns
            -------
            this : Eigen_Type_R
            	Object to be constructed
            
            
            Automatically generated constructor for eigen_type_r
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_eigen_type_r_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Eigen_Type_R
            
            
            Defined at Grid_module.fpp lines 40-44
            
            Parameters
            ----------
            this : Eigen_Type_R
            	Object to be destructed
            
            
            Automatically generated destructor for eigen_type_r
            """
            if self._alloc:
                _AresMainPy.f90wrap_eigen_type_r_finalise(this=self._handle)
        
        @property
        def wvf(self):
            """
            Element wvf ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_eigen_type_r__array__wvf(self._handle)
            if array_handle in self._arrays:
                wvf = self._arrays[array_handle]
            else:
                wvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_eigen_type_r__array__wvf)
                self._arrays[array_handle] = wvf
            return wvf
        
        @wvf.setter
        def wvf(self, wvf):
            self.wvf[...] = wvf
        
        @property
        def val(self):
            """
            Element val ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 44
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_eigen_type_r__array__val(self._handle)
            if array_handle in self._arrays:
                val = self._arrays[array_handle]
            else:
                val = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_eigen_type_r__array__val)
                self._arrays[array_handle] = val
            return val
        
        @val.setter
        def val(self, val):
            self.val[...] = val
        
        def __str__(self):
            ret = ['<eigen_type_r>{\n']
            ret.append('    wvf : ')
            ret.append(repr(self.wvf))
            ret.append(',\n    val : ')
            ret.append(repr(self.val))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("AresMainPy.charge_sphere")
    class charge_sphere(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=charge_sphere)
        
        
        Defined at Grid_module.fpp lines 48-57
        
        """
        def __init__(self, handle=None):
            """
            self = Charge_Sphere()
            
            
            Defined at Grid_module.fpp lines 48-57
            
            
            Returns
            -------
            this : Charge_Sphere
            	Object to be constructed
            
            
            Automatically generated constructor for charge_sphere
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_charge_sphere_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Charge_Sphere
            
            
            Defined at Grid_module.fpp lines 48-57
            
            Parameters
            ----------
            this : Charge_Sphere
            	Object to be destructed
            
            
            Automatically generated destructor for charge_sphere
            """
            if self._alloc:
                _AresMainPy.f90wrap_charge_sphere_finalise(this=self._handle)
        
        @property
        def onedlength(self):
            """
            Element onedlength ftype=integer(i4b) pytype=int
            
            
            Defined at Grid_module.fpp line 49
            
            """
            return _AresMainPy.f90wrap_charge_sphere__get__onedlength(self._handle)
        
        @onedlength.setter
        def onedlength(self, onedlength):
            _AresMainPy.f90wrap_charge_sphere__set__onedlength(self._handle, onedlength)
        
        @property
        def volume(self):
            """
            Element volume ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 50
            
            """
            return _AresMainPy.f90wrap_charge_sphere__get__volume(self._handle)
        
        @volume.setter
        def volume(self, volume):
            _AresMainPy.f90wrap_charge_sphere__set__volume(self._handle, volume)
        
        @property
        def x(self):
            """
            Element x ftype=integer(i4b) pytype=int
            
            
            Defined at Grid_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__x(self._handle)
            if array_handle in self._arrays:
                x = self._arrays[array_handle]
            else:
                x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__x)
                self._arrays[array_handle] = x
            return x
        
        @x.setter
        def x(self, x):
            self.x[...] = x
        
        @property
        def y(self):
            """
            Element y ftype=integer(i4b) pytype=int
            
            
            Defined at Grid_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__y(self._handle)
            if array_handle in self._arrays:
                y = self._arrays[array_handle]
            else:
                y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__y)
                self._arrays[array_handle] = y
            return y
        
        @y.setter
        def y(self, y):
            self.y[...] = y
        
        @property
        def z(self):
            """
            Element z ftype=integer(i4b) pytype=int
            
            
            Defined at Grid_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__z(self._handle)
            if array_handle in self._arrays:
                z = self._arrays[array_handle]
            else:
                z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__z)
                self._arrays[array_handle] = z
            return z
        
        @z.setter
        def z(self, z):
            self.z[...] = z
        
        @property
        def n(self):
            """
            Element n ftype=integer(i4b) pytype=int
            
            
            Defined at Grid_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__n(self._handle)
            if array_handle in self._arrays:
                n = self._arrays[array_handle]
            else:
                n = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__n)
                self._arrays[array_handle] = n
            return n
        
        @n.setter
        def n(self, n):
            self.n[...] = n
        
        @property
        def onedsphere(self):
            """
            Element onedsphere ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 52
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__onedsphere(self._handle)
            if array_handle in self._arrays:
                onedsphere = self._arrays[array_handle]
            else:
                onedsphere = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__onedsphere)
                self._arrays[array_handle] = onedsphere
            return onedsphere
        
        @onedsphere.setter
        def onedsphere(self, onedsphere):
            self.onedsphere[...] = onedsphere
        
        @property
        def onedveff(self):
            """
            Element onedveff ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 53
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__onedveff(self._handle)
            if array_handle in self._arrays:
                onedveff = self._arrays[array_handle]
            else:
                onedveff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__onedveff)
                self._arrays[array_handle] = onedveff
            return onedveff
        
        @onedveff.setter
        def onedveff(self, onedveff):
            self.onedveff[...] = onedveff
        
        @property
        def onedwvf(self):
            """
            Element onedwvf ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 54
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__onedwvf(self._handle)
            if array_handle in self._arrays:
                onedwvf = self._arrays[array_handle]
            else:
                onedwvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__onedwvf)
                self._arrays[array_handle] = onedwvf
            return onedwvf
        
        @onedwvf.setter
        def onedwvf(self, onedwvf):
            self.onedwvf[...] = onedwvf
        
        @property
        def onedval(self):
            """
            Element onedval ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 55
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__onedval(self._handle)
            if array_handle in self._arrays:
                onedval = self._arrays[array_handle]
            else:
                onedval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__onedval)
                self._arrays[array_handle] = onedval
            return onedval
        
        @onedval.setter
        def onedval(self, onedval):
            self.onedval[...] = onedval
        
        @property
        def initmo(self):
            """
            Element initmo ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 56
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__initmo(self._handle)
            if array_handle in self._arrays:
                initmo = self._arrays[array_handle]
            else:
                initmo = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__initmo)
                self._arrays[array_handle] = initmo
            return initmo
        
        @initmo.setter
        def initmo(self, initmo):
            self.initmo[...] = initmo
        
        @property
        def vlpp(self):
            """
            Element vlpp ftype=real(dp) pytype=float
            
            
            Defined at Grid_module.fpp line 57
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_charge_sphere__array__vlpp(self._handle)
            if array_handle in self._arrays:
                vlpp = self._arrays[array_handle]
            else:
                vlpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_charge_sphere__array__vlpp)
                self._arrays[array_handle] = vlpp
            return vlpp
        
        @vlpp.setter
        def vlpp(self, vlpp):
            self.vlpp[...] = vlpp
        
        def __str__(self):
            ret = ['<charge_sphere>{\n']
            ret.append('    onedlength : ')
            ret.append(repr(self.onedlength))
            ret.append(',\n    volume : ')
            ret.append(repr(self.volume))
            ret.append(',\n    x : ')
            ret.append(repr(self.x))
            ret.append(',\n    y : ')
            ret.append(repr(self.y))
            ret.append(',\n    z : ')
            ret.append(repr(self.z))
            ret.append(',\n    n : ')
            ret.append(repr(self.n))
            ret.append(',\n    onedsphere : ')
            ret.append(repr(self.onedsphere))
            ret.append(',\n    onedveff : ')
            ret.append(repr(self.onedveff))
            ret.append(',\n    onedwvf : ')
            ret.append(repr(self.onedwvf))
            ret.append(',\n    onedval : ')
            ret.append(repr(self.onedval))
            ret.append(',\n    initmo : ')
            ret.append(repr(self.initmo))
            ret.append(',\n    vlpp : ')
            ret.append(repr(self.vlpp))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def build_rgrid():
        """
        build_rgrid()
        
        
        Defined at Grid_module.fpp lines 102-165
        
        
        =============================================================
        isolate system and odd grid
         IF(.NOT.Lpbc.AND.MOD(n1,2)==0 )n1=n1+1
         IF(.NOT.Lpbc.AND.MOD(n2,2)==0 )n2=n2+1
         IF(.NOT.Lpbc.AND.MOD(n3,2)==0 )n3=n3+1
        =============================================================
        grid size
        """
        _AresMainPy.f90wrap_build_rgrid()
    
    @staticmethod
    def build_rgrid_iso():
        """
        build_rgrid_iso()
        
        
        Defined at Grid_module.fpp lines 168-219
        
        
        =============================================================
        isolate system and odd grid
        """
        _AresMainPy.f90wrap_build_rgrid_iso()
    
    @staticmethod
    def destroy_rgrid():
        """
        destroy_rgrid()
        
        
        Defined at Grid_module.fpp lines 222-251
        
        
        """
        _AresMainPy.f90wrap_destroy_rgrid()
    
    @staticmethod
    def build_kgrid():
        """
        build_kgrid()
        
        
        Defined at Grid_module.fpp lines 254-353
        
        
        """
        _AresMainPy.f90wrap_build_kgrid()
    
    @staticmethod
    def destroy_kpt():
        """
        destroy_kpt()
        
        
        Defined at Grid_module.fpp lines 356-372
        
        
        """
        _AresMainPy.f90wrap_destroy_kpt()
    
    @staticmethod
    def build_eigen():
        """
        build_eigen()
        
        
        Defined at Grid_module.fpp lines 375-419
        
        
        """
        _AresMainPy.f90wrap_build_eigen()
    
    @staticmethod
    def destroy_eigen():
        """
        destroy_eigen()
        
        
        Defined at Grid_module.fpp lines 422-444
        
        
        """
        _AresMainPy.f90wrap_destroy_eigen()
    
    @staticmethod
    def fillqtable():
        """
        fillqtable()
        
        
        Defined at Grid_module.fpp lines 447-505
        
        
        """
        _AresMainPy.f90wrap_fillqtable()
    
    @staticmethod
    def fillrtable():
        """
        fillrtable()
        
        
        Defined at Grid_module.fpp lines 508-589
        
        
        """
        _AresMainPy.f90wrap_fillrtable()
    
    @staticmethod
    def fillrtable_iso():
        """
        fillrtable_iso()
        
        
        Defined at Grid_module.fpp lines 592-620
        
        
        """
        _AresMainPy.f90wrap_fillrtable_iso()
    
    @staticmethod
    def build_iso_sphere_grid():
        """
        build_iso_sphere_grid()
        
        
        Defined at Grid_module.fpp lines 623-749
        
        
        ======================================
        ##XLT OBTAIN GRID INFORMATION
        print*,"orig",orig
        print*,"struct%poscar",struct%poscar
        ======================================
        """
        _AresMainPy.f90wrap_build_iso_sphere_grid()
    
    @staticmethod
    def matchposcar_iso():
        """
        matchposcar_iso()
        
        
        Defined at Grid_module.fpp lines 752-766
        
        
        =========================================================
        ##move atoms in line with the grid coordinate system
        """
        _AresMainPy.f90wrap_matchposcar_iso()
    
    @staticmethod
    def set_car2spe(orig, x, y, z):
        """
        r, cost, sint, cosp, sinp = set_car2spe(orig, x, y, z)
        
        
        Defined at Grid_module.fpp lines 771-802
        
        Parameters
        ----------
        orig : float array
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
        r, cost, sint, cosp, sinp = _AresMainPy.f90wrap_set_car2spe(orig=orig, x=x, y=y, \
            z=z)
        return r, cost, sint, cosp, sinp
    
    @staticmethod
    def destroy_iso_sphere_grid():
        """
        destroy_iso_sphere_grid()
        
        
        Defined at Grid_module.fpp lines 806-857
        
        
        """
        _AresMainPy.f90wrap_destroy_iso_sphere_grid()
    
    @staticmethod
    def iso_vsphere(vgrid, vsphere):
        """
        iso_vsphere(vgrid, vsphere)
        
        
        Defined at Grid_module.fpp lines 861-872
        
        Parameters
        ----------
        vgrid : float array
        vsphere : float array
        
        """
        _AresMainPy.f90wrap_iso_vsphere(vgrid=vgrid, vsphere=vsphere)
    
    @staticmethod
    def iso_rho2grid(rhosphere, rhogrid):
        """
        iso_rho2grid(rhosphere, rhogrid)
        
        
        Defined at Grid_module.fpp lines 876-888
        
        Parameters
        ----------
        rhosphere : float array
        rhogrid : float array
        
        """
        _AresMainPy.f90wrap_iso_rho2grid(rhosphere=rhosphere, rhogrid=rhogrid)
    
    @staticmethod
    def parallel_s2g(p, thrq):
        """
        parallel_s2g(p, thrq)
        
        
        Defined at Grid_module.fpp lines 891-914
        
        Parameters
        ----------
        p : float array
        thrq : float array
        
        """
        _AresMainPy.f90wrap_parallel_s2g(p=p, thrq=thrq)
    
    @staticmethod
    def parallel_g2s(p, thrq):
        """
        parallel_g2s(p, thrq)
        
        
        Defined at Grid_module.fpp lines 917-932
        
        Parameters
        ----------
        p : float array
        thrq : float array
        
        """
        _AresMainPy.f90wrap_parallel_g2s(p=p, thrq=thrq)
    
    @staticmethod
    def get_z_range(my_z_range, delta_z, ngrid_z, cell_thick, cell_shape):
        """
        get_z_range(my_z_range, delta_z, ngrid_z, cell_thick, cell_shape)
        
        
        Defined at Grid_module.fpp lines 935-953
        
        Parameters
        ----------
        my_z_range : int array
        delta_z : float
        ngrid_z : int
        cell_thick : float
        cell_shape : int
        
        """
        _AresMainPy.f90wrap_get_z_range(my_z_range=my_z_range, delta_z=delta_z, \
            ngrid_z=ngrid_z, cell_thick=cell_thick, cell_shape=cell_shape)
    
    @staticmethod
    def confirm_iso_radius():
        """
        confirm_iso_radius()
        
        
        Defined at Grid_module.fpp lines 956-1006
        
        
        ========================================
        > find MAX(r)
        """
        _AresMainPy.f90wrap_confirm_iso_radius()
    
    @staticmethod
    def reshape_center():
        """
        reshape_center()
        
        
        Defined at Grid_module.fpp lines 1009-1034
        
        
        -----------------------------------------------
        > reset the direct pos
        """
        _AresMainPy.f90wrap_reshape_center()
    
    @staticmethod
    def reset_poscar(new_gap, new_ni):
        """
        reset_poscar(new_gap, new_ni)
        
        
        Defined at Grid_module.fpp lines 1038-1071
        
        Parameters
        ----------
        new_gap : float
        new_ni : int
        
        -----------------------------------------------
        > reset the direct pos
        """
        _AresMainPy.f90wrap_reset_poscar(new_gap=new_gap, new_ni=new_ni)
    
    @staticmethod
    def build_parallel_cubic_grid():
        """
        build_parallel_cubic_grid()
        
        
        Defined at Grid_module.fpp lines 1074-1111
        
        
        """
        _AresMainPy.f90wrap_build_parallel_cubic_grid()
    
    @staticmethod
    def rho_trans1d(rho3d, rho1d):
        """
        rho_trans1d(rho3d, rho1d)
        
        
        Defined at Grid_module.fpp lines 1113-1157
        
        Parameters
        ----------
        rho3d : float array
        rho1d : float array
        
        """
        _AresMainPy.f90wrap_rho_trans1d(rho3d=rho3d, rho1d=rho1d)
    
    @staticmethod
    def rho_trans3d(rho1d, rho3d):
        """
        rho_trans3d(rho1d, rho3d)
        
        
        Defined at Grid_module.fpp lines 1159-1203
        
        Parameters
        ----------
        rho1d : float array
        rho3d : float array
        
        """
        _AresMainPy.f90wrap_rho_trans3d(rho1d=rho1d, rho3d=rho3d)
    
    @staticmethod
    def _fft_sph_r2c(array_r):
        """
        array_c = _fft_sph_r2c(array_r)
        
        
        Defined at Grid_module.fpp lines 1206-1224
        
        Parameters
        ----------
        array_r : float array
        
        Returns
        -------
        array_c : complex array
        
        """
        array_c = _AresMainPy.f90wrap_fft_sph_r2c(array_r=array_r)
        return array_c
    
    @staticmethod
    def _fft_sph_c2r(array_c):
        """
        array_r = _fft_sph_c2r(array_c)
        
        
        Defined at Grid_module.fpp lines 1227-1246
        
        Parameters
        ----------
        array_c : complex array
        
        Returns
        -------
        array_r : float array
        
        """
        array_r = _AresMainPy.f90wrap_fft_sph_c2r(array_c=array_c)
        return array_r
    
    @staticmethod
    def fft_sph(*args, **kwargs):
        """
        fft_sph(*args, **kwargs)
        
        
        Defined at Grid_module.fpp lines 96-98
        
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
        
        
        Defined at Grid_module.fpp line 71
        
        """
        return _AresMainPy.f90wrap_grid_module__get__n1()
    
    @n1.setter
    def n1(self, n1):
        _AresMainPy.f90wrap_grid_module__set__n1(n1)
    
    @property
    def n2(self):
        """
        Element n2 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 71
        
        """
        return _AresMainPy.f90wrap_grid_module__get__n2()
    
    @n2.setter
    def n2(self, n2):
        _AresMainPy.f90wrap_grid_module__set__n2(n2)
    
    @property
    def n3(self):
        """
        Element n3 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 71
        
        """
        return _AresMainPy.f90wrap_grid_module__get__n3()
    
    @n3.setter
    def n3(self, n3):
        _AresMainPy.f90wrap_grid_module__set__n3(n3)
    
    @property
    def n(self):
        """
        Element n ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 71
        
        """
        return _AresMainPy.f90wrap_grid_module__get__n()
    
    @n.setter
    def n(self, n):
        _AresMainPy.f90wrap_grid_module__set__n(n)
    
    @property
    def nsn(self):
        """
        Element nsn ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 71
        
        """
        return _AresMainPy.f90wrap_grid_module__get__nsn()
    
    @nsn.setter
    def nsn(self, nsn):
        _AresMainPy.f90wrap_grid_module__set__nsn(nsn)
    
    @property
    def global_n1(self):
        """
        Element global_n1 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 72
        
        """
        return _AresMainPy.f90wrap_grid_module__get__global_n1()
    
    @global_n1.setter
    def global_n1(self, global_n1):
        _AresMainPy.f90wrap_grid_module__set__global_n1(global_n1)
    
    @property
    def global_n2(self):
        """
        Element global_n2 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 72
        
        """
        return _AresMainPy.f90wrap_grid_module__get__global_n2()
    
    @global_n2.setter
    def global_n2(self, global_n2):
        _AresMainPy.f90wrap_grid_module__set__global_n2(global_n2)
    
    @property
    def global_n3(self):
        """
        Element global_n3 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 72
        
        """
        return _AresMainPy.f90wrap_grid_module__get__global_n3()
    
    @global_n3.setter
    def global_n3(self, global_n3):
        _AresMainPy.f90wrap_grid_module__set__global_n3(global_n3)
    
    @property
    def global_n(self):
        """
        Element global_n ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 72
        
        """
        return _AresMainPy.f90wrap_grid_module__get__global_n()
    
    @global_n.setter
    def global_n(self, global_n):
        _AresMainPy.f90wrap_grid_module__set__global_n(global_n)
    
    @property
    def ng1(self):
        """
        Element ng1 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 73
        
        """
        return _AresMainPy.f90wrap_grid_module__get__ng1()
    
    @ng1.setter
    def ng1(self, ng1):
        _AresMainPy.f90wrap_grid_module__set__ng1(ng1)
    
    @property
    def ng2(self):
        """
        Element ng2 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 73
        
        """
        return _AresMainPy.f90wrap_grid_module__get__ng2()
    
    @ng2.setter
    def ng2(self, ng2):
        _AresMainPy.f90wrap_grid_module__set__ng2(ng2)
    
    @property
    def ng3(self):
        """
        Element ng3 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 73
        
        """
        return _AresMainPy.f90wrap_grid_module__get__ng3()
    
    @ng3.setter
    def ng3(self, ng3):
        _AresMainPy.f90wrap_grid_module__set__ng3(ng3)
    
    @property
    def ng(self):
        """
        Element ng ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 73
        
        """
        return _AresMainPy.f90wrap_grid_module__get__ng()
    
    @ng.setter
    def ng(self, ng):
        _AresMainPy.f90wrap_grid_module__set__ng(ng)
    
    @property
    def gap(self):
        """
        Element gap ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__gap(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            gap = self._arrays[array_handle]
        else:
            gap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__gap)
            self._arrays[array_handle] = gap
        return gap
    
    @gap.setter
    def gap(self, gap):
        self.gap[...] = gap
    
    @property
    def dvol(self):
        """
        Element dvol ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 75
        
        """
        return _AresMainPy.f90wrap_grid_module__get__dvol()
    
    @dvol.setter
    def dvol(self, dvol):
        _AresMainPy.f90wrap_grid_module__set__dvol(dvol)
    
    @property
    def nk1(self):
        """
        Element nk1 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 77
        
        """
        return _AresMainPy.f90wrap_grid_module__get__nk1()
    
    @nk1.setter
    def nk1(self, nk1):
        _AresMainPy.f90wrap_grid_module__set__nk1(nk1)
    
    @property
    def nk2(self):
        """
        Element nk2 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 77
        
        """
        return _AresMainPy.f90wrap_grid_module__get__nk2()
    
    @nk2.setter
    def nk2(self, nk2):
        _AresMainPy.f90wrap_grid_module__set__nk2(nk2)
    
    @property
    def nk3(self):
        """
        Element nk3 ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 77
        
        """
        return _AresMainPy.f90wrap_grid_module__get__nk3()
    
    @nk3.setter
    def nk3(self, nk3):
        _AresMainPy.f90wrap_grid_module__set__nk3(nk3)
    
    @property
    def nk(self):
        """
        Element nk ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 77
        
        """
        return _AresMainPy.f90wrap_grid_module__get__nk()
    
    @nk.setter
    def nk(self, nk):
        _AresMainPy.f90wrap_grid_module__set__nk(nk)
    
    @property
    def kdispl(self):
        """
        Element kdispl ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__kdispl(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            kdispl = self._arrays[array_handle]
        else:
            kdispl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__kdispl)
            self._arrays[array_handle] = kdispl
        return kdispl
    
    @kdispl.setter
    def kdispl(self, kdispl):
        self.kdispl[...] = kdispl
    
    @property
    def dr(self):
        """
        Element dr ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 88
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__dr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dr = self._arrays[array_handle]
        else:
            dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__dr)
            self._arrays[array_handle] = dr
        return dr
    
    @dr.setter
    def dr(self, dr):
        self.dr[...] = dr
    
    @property
    def cost(self):
        """
        Element cost ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__cost(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cost = self._arrays[array_handle]
        else:
            cost = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__cost)
            self._arrays[array_handle] = cost
        return cost
    
    @cost.setter
    def cost(self, cost):
        self.cost[...] = cost
    
    @property
    def sint(self):
        """
        Element sint ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 90
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__sint(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sint = self._arrays[array_handle]
        else:
            sint = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__sint)
            self._arrays[array_handle] = sint
        return sint
    
    @sint.setter
    def sint(self, sint):
        self.sint[...] = sint
    
    @property
    def cns(self):
        """
        Element cns ftype=complex(dcp) pytype=complex
        
        
        Defined at Grid_module.fpp line 91
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__cns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cns = self._arrays[array_handle]
        else:
            cns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__cns)
            self._arrays[array_handle] = cns
        return cns
    
    @cns.setter
    def cns(self, cns):
        self.cns[...] = cns
    
    @property
    def indx(self):
        """
        Element indx ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 92
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_grid_module__array__indx(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            indx = self._arrays[array_handle]
        else:
            indx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_grid_module__array__indx)
            self._arrays[array_handle] = indx
        return indx
    
    @indx.setter
    def indx(self, indx):
        self.indx[...] = indx
    
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
        ret.append(',\n    nsn : ')
        ret.append(repr(self.nsn))
        ret.append(',\n    global_n1 : ')
        ret.append(repr(self.global_n1))
        ret.append(',\n    global_n2 : ')
        ret.append(repr(self.global_n2))
        ret.append(',\n    global_n3 : ')
        ret.append(repr(self.global_n3))
        ret.append(',\n    global_n : ')
        ret.append(repr(self.global_n))
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
        ret.append(',\n    dr : ')
        ret.append(repr(self.dr))
        ret.append(',\n    cost : ')
        ret.append(repr(self.cost))
        ret.append(',\n    sint : ')
        ret.append(repr(self.sint))
        ret.append(',\n    cns : ')
        ret.append(repr(self.cns))
        ret.append(',\n    indx : ')
        ret.append(repr(self.indx))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

grid_module = Grid_Module()

class Ewald(f90wrap.runtime.FortranModule):
    """
    Module ewald
    
    
    Defined at Ewald.fpp lines 5-826
    
    """
    @f90wrap.runtime.register_class("AresMainPy.ion")
    class ion(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=ion)
        
        
        Defined at Ewald.fpp lines 36-38
        
        """
        def __init__(self, handle=None):
            """
            self = Ion()
            
            
            Defined at Ewald.fpp lines 36-38
            
            
            Returns
            -------
            this : Ion
            	Object to be constructed
            
            
            Automatically generated constructor for ion
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_ion_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Ion
            
            
            Defined at Ewald.fpp lines 36-38
            
            Parameters
            ----------
            this : Ion
            	Object to be destructed
            
            
            Automatically generated destructor for ion
            """
            if self._alloc:
                _AresMainPy.f90wrap_ion_finalise(this=self._handle)
        
        @property
        def charge(self):
            """
            Element charge ftype=real(dp) pytype=float
            
            
            Defined at Ewald.fpp line 37
            
            """
            return _AresMainPy.f90wrap_ion__get__charge(self._handle)
        
        @charge.setter
        def charge(self, charge):
            _AresMainPy.f90wrap_ion__set__charge(self._handle, charge)
        
        @property
        def fcd(self):
            """
            Element fcd ftype=real(dp) pytype=float
            
            
            Defined at Ewald.fpp line 38
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_ion__array__fcd(self._handle)
            if array_handle in self._arrays:
                fcd = self._arrays[array_handle]
            else:
                fcd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_ion__array__fcd)
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
        
        
        Defined at Ewald.fpp lines 58-76
        
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
        ewald_energy = _AresMainPy.f90wrap_ewald_energy(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return ewald_energy
    
    @staticmethod
    def iso_ewald_energy(latticev, ionpositions, iontpid, ioncharges):
        """
        iso_ewald_energy = iso_ewald_energy(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.fpp lines 80-108
        
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
        iso_ewald_energy = _AresMainPy.f90wrap_iso_ewald_energy(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return iso_ewald_energy
    
    @staticmethod
    def iso_ewald_forces(latticev, ionpositions, iontpid, ioncharges):
        """
        iso_ewald_forces = iso_ewald_forces(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.fpp lines 112-193
        
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
        iso_ewald_forces = _AresMainPy.f90wrap_iso_ewald_forces(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return iso_ewald_forces
    
    @staticmethod
    def ewald_forces(latticev, ionpositions, iontpid, ioncharges):
        """
        ewald_forces = ewald_forces(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.fpp lines 196-219
        
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
        ewald_forces = _AresMainPy.f90wrap_ewald_forces(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return ewald_forces
    
    @staticmethod
    def ewald_stress(latticev, ionpositions, iontpid, ioncharges):
        """
        ewald_stress = ewald_stress(latticev, ionpositions, iontpid, ioncharges)
        
        
        Defined at Ewald.fpp lines 222-240
        
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
        ewald_stress = _AresMainPy.f90wrap_ewald_stress(latticev=latticev, \
            ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
        return ewald_stress
    
    @staticmethod
    def vectorlength(vc):
        """
        vectorlength = vectorlength(vc)
        
        
        Defined at Ewald.fpp lines 715-716
        
        Parameters
        ----------
        vc : float array
        
        Returns
        -------
        vectorlength : float
        
        """
        vectorlength = _AresMainPy.f90wrap_vectorlength(vc=vc)
        return vectorlength
    
    @staticmethod
    def recipvector(lat):
        """
        recipvector = recipvector(lat)
        
        
        Defined at Ewald.fpp lines 719-724
        
        Parameters
        ----------
        lat : float array
        
        Returns
        -------
        recipvector : float array
        
        """
        recipvector = _AresMainPy.f90wrap_recipvector(lat=lat)
        return recipvector
    
    @staticmethod
    def volume(lat):
        """
        volume = volume(lat)
        
        
        Defined at Ewald.fpp lines 727-729
        
        Parameters
        ----------
        lat : float array
        
        Returns
        -------
        volume : float
        
        """
        volume = _AresMainPy.f90wrap_volume(lat=lat)
        return volume
    
    @staticmethod
    def crossp(va, vb):
        """
        crossp = crossp(va, vb)
        
        
        Defined at Ewald.fpp lines 732-736
        
        Parameters
        ----------
        va : float array
        vb : float array
        
        Returns
        -------
        crossp : float array
        
        """
        crossp = _AresMainPy.f90wrap_crossp(va=va, vb=vb)
        return crossp
    
    @staticmethod
    def erfc(x):
        """
        erfc = erfc(x)
        
        
        Defined at Ewald.fpp lines 739-825
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        erfc : float
        
        """
        erfc = _AresMainPy.f90wrap_erfc(x=x)
        return erfc
    
    @property
    def bohr(self):
        """
        Element bohr ftype=real(dp) pytype=float
        
        
        Defined at Ewald.fpp line 40
        
        """
        return _AresMainPy.f90wrap_ewald__get__bohr()
    
    @property
    def hartreetoev(self):
        """
        Element hartreetoev ftype=real(dp) pytype=float
        
        
        Defined at Ewald.fpp line 40
        
        """
        return _AresMainPy.f90wrap_ewald__get__hartreetoev()
    
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
    
    
    Defined at Finite_module.fpp lines 5-2216
    
    """
    @staticmethod
    def init_finite(norder, h):
        """
        init_finite(norder, h)
        
        
        Defined at Finite_module.fpp lines 39-88
        
        Parameters
        ----------
        norder : int
        h : float array
        
        """
        _AresMainPy.f90wrap_init_finite(norder=norder, h=h)
    
    @staticmethod
    def destroy_finite():
        """
        destroy_finite()
        
        
        Defined at Finite_module.fpp lines 91-103
        
        
        """
        _AresMainPy.f90wrap_destroy_finite()
    
    @staticmethod
    def trans_mat_full(mat, factor, int_miu, lad_gap, err):
        """
        trans_mat_full(mat, factor, int_miu, lad_gap, err)
        
        
        Defined at Finite_module.fpp lines 106-309
        
        Parameters
        ----------
        mat : float array
        factor : float array
        int_miu : int array
        lad_gap : float array
        err : float array
        
        """
        _AresMainPy.f90wrap_trans_mat_full(mat=mat, factor=factor, int_miu=int_miu, \
            lad_gap=lad_gap, err=err)
    
    @staticmethod
    def cmplx_nabla2(ifun, norder, ofun):
        """
        cmplx_nabla2(ifun, norder, ofun)
        
        
        Defined at Finite_module.fpp lines 312-439
        
        Parameters
        ----------
        ifun : complex array
        norder : int
        ofun : complex array
        
        """
        _AresMainPy.f90wrap_cmplx_nabla2(ifun=ifun, norder=norder, ofun=ofun)
    
    @staticmethod
    def cmplx_nabla1(ifun, norder, derf):
        """
        cmplx_nabla1(ifun, norder, derf)
        
        
        Defined at Finite_module.fpp lines 442-634
        
        Parameters
        ----------
        ifun : complex array
        norder : int
        derf : complex array
        
        """
        _AresMainPy.f90wrap_cmplx_nabla1(ifun=ifun, norder=norder, derf=derf)
    
    @staticmethod
    def nabla(ifun, norder, derf):
        """
        nabla(ifun, norder, derf)
        
        
        Defined at Finite_module.fpp lines 637-741
        
        Parameters
        ----------
        ifun : complex array
        norder : int
        derf : complex array
        
        """
        _AresMainPy.f90wrap_nabla(ifun=ifun, norder=norder, derf=derf)
    
    @staticmethod
    def kskeperiod_op(psi, ik, kepsi):
        """
        kskeperiod_op(psi, ik, kepsi)
        
        
        Defined at Finite_module.fpp lines 743-866
        
        Parameters
        ----------
        psi : complex array
        ik : int
        kepsi : complex array
        
        """
        _AresMainPy.f90wrap_kskeperiod_op(psi=psi, ik=ik, kepsi=kepsi)
    
    @staticmethod
    def ke1(psi, ts1):
        """
        ke1(psi, ts1)
        
        
        Defined at Finite_module.fpp lines 869-880
        
        Parameters
        ----------
        psi : complex array
        ts1 : complex array
        
        """
        _AresMainPy.f90wrap_ke1(psi=psi, ts1=ts1)
    
    @staticmethod
    def ke2(psi, ik, ts2):
        """
        ke2(psi, ik, ts2)
        
        
        Defined at Finite_module.fpp lines 883-909
        
        Parameters
        ----------
        psi : complex array
        ik : int
        ts2 : complex array
        
        """
        _AresMainPy.f90wrap_ke2(psi=psi, ik=ik, ts2=ts2)
    
    @staticmethod
    def realnabla2(ifun, norder, ofun):
        """
        realnabla2(ifun, norder, ofun)
        
        
        Defined at Finite_module.fpp lines 918-1008
        
        Parameters
        ----------
        ifun : float array
        norder : int
        ofun : float array
        
        """
        _AresMainPy.f90wrap_realnabla2(ifun=ifun, norder=norder, ofun=ofun)
    
    @staticmethod
    def realke(psi, kepsi):
        """
        realke(psi, kepsi)
        
        
        Defined at Finite_module.fpp lines 1011-1023
        
        Parameters
        ----------
        psi : float array
        kepsi : float array
        
        """
        _AresMainPy.f90wrap_realke(psi=psi, kepsi=kepsi)
    
    @staticmethod
    def iso_realnabla2(ifun, norder, ofun):
        """
        iso_realnabla2(ifun, norder, ofun)
        
        
        Defined at Finite_module.fpp lines 1033-1188
        
        Parameters
        ----------
        ifun : float array
        norder : int
        ofun : float array
        
        ==================================================
        """
        _AresMainPy.f90wrap_iso_realnabla2(ifun=ifun, norder=norder, ofun=ofun)
    
    @staticmethod
    def iso_realke(psi, kepsi):
        """
        iso_realke(psi, kepsi)
        
        
        Defined at Finite_module.fpp lines 1243-1255
        
        Parameters
        ----------
        psi : float array
        kepsi : float array
        
        """
        _AresMainPy.f90wrap_iso_realke(psi=psi, kepsi=kepsi)
    
    @staticmethod
    def kskeperiod_op_band(psi, ik, kepsi):
        """
        kskeperiod_op_band(psi, ik, kepsi)
        
        
        Defined at Finite_module.fpp lines 1258-1308
        
        Parameters
        ----------
        psi : complex array
        ik : int
        kepsi : complex array
        
        """
        _AresMainPy.f90wrap_kskeperiod_op_band(psi=psi, ik=ik, kepsi=kepsi)
    
    @staticmethod
    def nabla2_np(ifun, norder, ofun):
        """
        nabla2_np(ifun, norder, ofun)
        
        
        Defined at Finite_module.fpp lines 1310-1431
        
        Parameters
        ----------
        ifun : complex array
        norder : int
        ofun : complex array
        
        """
        _AresMainPy.f90wrap_nabla2_np(ifun=ifun, norder=norder, ofun=ofun)
    
    @staticmethod
    def nabla1_np(ifun, norder, derf):
        """
        nabla1_np(ifun, norder, derf)
        
        
        Defined at Finite_module.fpp lines 1434-1512
        
        Parameters
        ----------
        ifun : complex array
        norder : int
        derf : complex array
        
        """
        _AresMainPy.f90wrap_nabla1_np(ifun=ifun, norder=norder, derf=derf)
    
    @staticmethod
    def ke1_np(psi, ts1):
        """
        ke1_np(psi, ts1)
        
        
        Defined at Finite_module.fpp lines 1515-1526
        
        Parameters
        ----------
        psi : complex array
        ts1 : complex array
        
        """
        _AresMainPy.f90wrap_ke1_np(psi=psi, ts1=ts1)
    
    @staticmethod
    def ke2_np(psi, ik, ts2):
        """
        ke2_np(psi, ik, ts2)
        
        
        Defined at Finite_module.fpp lines 1529-1555
        
        Parameters
        ----------
        psi : complex array
        ik : int
        ts2 : complex array
        
        """
        _AresMainPy.f90wrap_ke2_np(psi=psi, ik=ik, ts2=ts2)
    
    @staticmethod
    def kskeperiod_op_gamma(psi, ik, kepsi):
        """
        kskeperiod_op_gamma(psi, ik, kepsi)
        
        
        Defined at Finite_module.fpp lines 1557-1654
        
        Parameters
        ----------
        psi : float array
        ik : int
        kepsi : float array
        
        """
        _AresMainPy.f90wrap_kskeperiod_op_gamma(psi=psi, ik=ik, kepsi=kepsi)
    
    @staticmethod
    def ke_gamma(psi, ts1):
        """
        ke_gamma(psi, ts1)
        
        
        Defined at Finite_module.fpp lines 1656-1668
        
        Parameters
        ----------
        psi : float array
        ts1 : float array
        
        """
        _AresMainPy.f90wrap_ke_gamma(psi=psi, ts1=ts1)
    
    @staticmethod
    def real_nabla2(ifun, norder, ofun):
        """
        real_nabla2(ifun, norder, ofun)
        
        
        Defined at Finite_module.fpp lines 1670-1795
        
        Parameters
        ----------
        ifun : float array
        norder : int
        ofun : float array
        
        """
        _AresMainPy.f90wrap_real_nabla2(ifun=ifun, norder=norder, ofun=ofun)
    
    @staticmethod
    def real_nabla2_1d(ifun, ofun):
        """
        real_nabla2_1d(ifun, ofun)
        
        
        Defined at Finite_module.fpp lines 1797-1928
        
        Parameters
        ----------
        ifun : float array
        ofun : float array
        
        """
        _AresMainPy.f90wrap_real_nabla2_1d(ifun=ifun, ofun=ofun)
    
    @staticmethod
    def nabla_gamma(ifun, norder, derf):
        """
        nabla_gamma(ifun, norder, derf)
        
        
        Defined at Finite_module.fpp lines 1930-2022
        
        Parameters
        ----------
        ifun : float array
        norder : int
        derf : float array
        
        """
        _AresMainPy.f90wrap_nabla_gamma(ifun=ifun, norder=norder, derf=derf)
    
    @staticmethod
    def real_nabla1(ifun, norder, derf):
        """
        real_nabla1(ifun, norder, derf)
        
        
        Defined at Finite_module.fpp lines 2024-2216
        
        Parameters
        ----------
        ifun : float array
        norder : int
        derf : float array
        
        """
        _AresMainPy.f90wrap_real_nabla1(ifun=ifun, norder=norder, derf=derf)
    
    @property
    def lapl(self):
        """
        Element lapl ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__lapl(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lapl = self._arrays[array_handle]
        else:
            lapl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__lapl)
            self._arrays[array_handle] = lapl
        return lapl
    
    @lapl.setter
    def lapl(self, lapl):
        self.lapl[...] = lapl
    
    @property
    def grad(self):
        """
        Element grad ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__grad(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            grad = self._arrays[array_handle]
        else:
            grad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__grad)
            self._arrays[array_handle] = grad
        return grad
    
    @grad.setter
    def grad(self, grad):
        self.grad[...] = grad
    
    @property
    def tbmat(self):
        """
        Element tbmat ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__tbmat(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tbmat = self._arrays[array_handle]
        else:
            tbmat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__tbmat)
            self._arrays[array_handle] = tbmat
        return tbmat
    
    @tbmat.setter
    def tbmat(self, tbmat):
        self.tbmat[...] = tbmat
    
    @property
    def lap_add(self):
        """
        Element lap_add ftype=integer(i4b) pytype=int
        
        
        Defined at Finite_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__lap_add(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            lap_add = self._arrays[array_handle]
        else:
            lap_add = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__lap_add)
            self._arrays[array_handle] = lap_add
        return lap_add
    
    @lap_add.setter
    def lap_add(self, lap_add):
        self.lap_add[...] = lap_add
    
    @property
    def cell_mu(self):
        """
        Element cell_mu ftype=integer  pytype=int
        
        
        Defined at Finite_module.fpp line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__cell_mu(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cell_mu = self._arrays[array_handle]
        else:
            cell_mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__cell_mu)
            self._arrays[array_handle] = cell_mu
        return cell_mu
    
    @cell_mu.setter
    def cell_mu(self, cell_mu):
        self.cell_mu[...] = cell_mu
    
    @property
    def cell_factor(self):
        """
        Element cell_factor ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 28
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__cell_factor(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cell_factor = self._arrays[array_handle]
        else:
            cell_factor = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__cell_factor)
            self._arrays[array_handle] = cell_factor
        return cell_factor
    
    @cell_factor.setter
    def cell_factor(self, cell_factor):
        self.cell_factor[...] = cell_factor
    
    @property
    def wrap_box(self):
        """
        Element wrap_box ftype=complex(dcp) pytype=complex
        
        
        Defined at Finite_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__wrap_box(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wrap_box = self._arrays[array_handle]
        else:
            wrap_box = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__wrap_box)
            self._arrays[array_handle] = wrap_box
        return wrap_box
    
    @wrap_box.setter
    def wrap_box(self, wrap_box):
        self.wrap_box[...] = wrap_box
    
    @property
    def fun_global(self):
        """
        Element fun_global ftype=complex(dcp) pytype=complex
        
        
        Defined at Finite_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__fun_global(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fun_global = self._arrays[array_handle]
        else:
            fun_global = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__fun_global)
            self._arrays[array_handle] = fun_global
        return fun_global
    
    @fun_global.setter
    def fun_global(self, fun_global):
        self.fun_global[...] = fun_global
    
    @property
    def fun_1d(self):
        """
        Element fun_1d ftype=complex(dcp) pytype=complex
        
        
        Defined at Finite_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__fun_1d(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fun_1d = self._arrays[array_handle]
        else:
            fun_1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__fun_1d)
            self._arrays[array_handle] = fun_1d
        return fun_1d
    
    @fun_1d.setter
    def fun_1d(self, fun_1d):
        self.fun_1d[...] = fun_1d
    
    @property
    def wrap_box1d(self):
        """
        Element wrap_box1d ftype=complex(dcp) pytype=complex
        
        
        Defined at Finite_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__wrap_box1d(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wrap_box1d = self._arrays[array_handle]
        else:
            wrap_box1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__wrap_box1d)
            self._arrays[array_handle] = wrap_box1d
        return wrap_box1d
    
    @wrap_box1d.setter
    def wrap_box1d(self, wrap_box1d):
        self.wrap_box1d[...] = wrap_box1d
    
    @property
    def wrap_box_real(self):
        """
        Element wrap_box_real ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__wrap_box_real(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wrap_box_real = self._arrays[array_handle]
        else:
            wrap_box_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__wrap_box_real)
            self._arrays[array_handle] = wrap_box_real
        return wrap_box_real
    
    @wrap_box_real.setter
    def wrap_box_real(self, wrap_box_real):
        self.wrap_box_real[...] = wrap_box_real
    
    @property
    def fun_global_real(self):
        """
        Element fun_global_real ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__fun_global_real(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fun_global_real = self._arrays[array_handle]
        else:
            fun_global_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__fun_global_real)
            self._arrays[array_handle] = fun_global_real
        return fun_global_real
    
    @fun_global_real.setter
    def fun_global_real(self, fun_global_real):
        self.fun_global_real[...] = fun_global_real
    
    @property
    def fun_1d_real(self):
        """
        Element fun_1d_real ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__fun_1d_real(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fun_1d_real = self._arrays[array_handle]
        else:
            fun_1d_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__fun_1d_real)
            self._arrays[array_handle] = fun_1d_real
        return fun_1d_real
    
    @fun_1d_real.setter
    def fun_1d_real(self, fun_1d_real):
        self.fun_1d_real[...] = fun_1d_real
    
    @property
    def wrap_box1d_real(self):
        """
        Element wrap_box1d_real ftype=real(dp) pytype=float
        
        
        Defined at Finite_module.fpp line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_finite_module__array__wrap_box1d_real(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wrap_box1d_real = self._arrays[array_handle]
        else:
            wrap_box1d_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_finite_module__array__wrap_box1d_real)
            self._arrays[array_handle] = wrap_box1d_real
        return wrap_box1d_real
    
    @wrap_box1d_real.setter
    def wrap_box1d_real(self, wrap_box1d_real):
        self.wrap_box1d_real[...] = wrap_box1d_real
    
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
        ret.append(',\n    wrap_box : ')
        ret.append(repr(self.wrap_box))
        ret.append(',\n    fun_global : ')
        ret.append(repr(self.fun_global))
        ret.append(',\n    fun_1d : ')
        ret.append(repr(self.fun_1d))
        ret.append(',\n    wrap_box1d : ')
        ret.append(repr(self.wrap_box1d))
        ret.append(',\n    wrap_box_real : ')
        ret.append(repr(self.wrap_box_real))
        ret.append(',\n    fun_global_real : ')
        ret.append(repr(self.fun_global_real))
        ret.append(',\n    fun_1d_real : ')
        ret.append(repr(self.fun_1d_real))
        ret.append(',\n    wrap_box1d_real : ')
        ret.append(repr(self.wrap_box1d_real))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

finite_module = Finite_Module()

class Libxc_Module(f90wrap.runtime.FortranModule):
    """
    Module libxc_module
    
    
    Defined at XC_functional.fpp lines 5-344
    
    """
    @staticmethod
    def ldalib_energy(rhos):
        """
        exc = ldalib_energy(rhos)
        
        
        Defined at XC_functional.fpp lines 19-138
        
        Parameters
        ----------
        rhos : float array
        
        Returns
        -------
        exc : float
        
        """
        exc = _AresMainPy.f90wrap_ldalib_energy(rhos=rhos)
        return exc
    
    @staticmethod
    def ldalib_potential(rhos, vxcs):
        """
        ldalib_potential(rhos, vxcs)
        
        
        Defined at XC_functional.fpp lines 141-249
        
        Parameters
        ----------
        rhos : float array
        vxcs : float array
        
        """
        _AresMainPy.f90wrap_ldalib_potential(rhos=rhos, vxcs=vxcs)
    
    @staticmethod
    def ldalib_energy_iso(rhos):
        """
        exc = ldalib_energy_iso(rhos)
        
        
        Defined at XC_functional.fpp lines 252-298
        
        Parameters
        ----------
        rhos : float array
        
        Returns
        -------
        exc : float
        
        """
        exc = _AresMainPy.f90wrap_ldalib_energy_iso(rhos=rhos)
        return exc
    
    @staticmethod
    def ldalib_potential_iso(rhos, vxcs):
        """
        ldalib_potential_iso(rhos, vxcs)
        
        
        Defined at XC_functional.fpp lines 301-343
        
        Parameters
        ----------
        rhos : float array
        vxcs : float array
        
        """
        _AresMainPy.f90wrap_ldalib_potential_iso(rhos=rhos, vxcs=vxcs)
    
    _dt_array_initialisers = []
    

libxc_module = Libxc_Module()

class Sgfft_Oct_M(f90wrap.runtime.FortranModule):
    """
    Module sgfft_oct_m
    
    
    Defined at sgfft.fpp lines 26-4425
    
    """
    @staticmethod
    def fourier_dim(n):
        """
        n_next = fourier_dim(n)
        
        
        Defined at sgfft.fpp lines 51-81
        
        Parameters
        ----------
        n : int
        
        Returns
        -------
        n_next : int
        
        """
        n_next = _AresMainPy.f90wrap_fourier_dim(n=n)
        return n_next
    
    @staticmethod
    def fft(n1xy, n2xy, n3xy, nd1, nd2, nd3, z, isign, inzee):
        """
        fft(n1xy, n2xy, n3xy, nd1, nd2, nd3, z, isign, inzee)
        
        
        Defined at sgfft.fpp lines 161-385
        
        Parameters
        ----------
        n1xy : int
        n2xy : int
        n3xy : int
        nd1 : int
        nd2 : int
        nd3 : int
        z : float array
        isign : int
        inzee : int
        
        """
        _AresMainPy.f90wrap_fft(n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, nd1=nd1, nd2=nd2, \
            nd3=nd3, z=z, isign=isign, inzee=inzee)
    
    @staticmethod
    def convolxc_off(n1xy, n2xy, n3xy, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, \
        pot, zf, scal, comm):
        """
        convolxc_off(n1xy, n2xy, n3xy, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, pot, \
            zf, scal, comm)
        
        
        Defined at sgfft.fpp lines 3504-3796
        
        Parameters
        ----------
        n1xy : int
        n2xy : int
        n3xy : int
        nd1 : int
        nd2 : int
        nd3 : int
        md1 : int
        md2 : int
        md3 : int
        nproc : int
        iproc : int
        pot : float array
        zf : float array
        scal : float
        comm : int
        
        """
        _AresMainPy.f90wrap_convolxc_off(n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, nd1=nd1, \
            nd2=nd2, nd3=nd3, md1=md1, md2=md2, md3=md3, nproc=nproc, iproc=iproc, \
            pot=pot, zf=zf, scal=scal, comm=comm)
    
    _dt_array_initialisers = []
    

sgfft_oct_m = Sgfft_Oct_M()

class Poisson_Isf(f90wrap.runtime.FortranModule):
    """
    Module poisson_isf
    
    
    Defined at poisson_isf.fpp lines 5-1430
    
    """
    @staticmethod
    def karray_set():
        """
        karray_set()
        
        
        Defined at poisson_isf.fpp lines 72-94
        
        
        """
        _AresMainPy.f90wrap_karray_set()
    
    @staticmethod
    def poisson_isf_method(rho, vh):
        """
        poisson_isf_method(rho, vh)
        
        
        Defined at poisson_isf.fpp lines 97-110
        
        Parameters
        ----------
        rho : float array
        vh : float array
        
        """
        _AresMainPy.f90wrap_poisson_isf_method(rho=rho, vh=vh)
    
    @staticmethod
    def build_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, itype_scf, \
        karrayout):
        """
        build_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, itype_scf, karrayout)
        
        
        Defined at poisson_isf.fpp lines 113-249
        
        Parameters
        ----------
        n01 : int
        n02 : int
        n03 : int
        nfft1 : int
        nfft2 : int
        nfft3 : int
        hgrid : float
        itype_scf : int
        karrayout : float array
        
        """
        _AresMainPy.f90wrap_build_kernel(n01=n01, n02=n02, n03=n03, nfft1=nfft1, \
            nfft2=nfft2, nfft3=nfft3, hgrid=hgrid, itype_scf=itype_scf, \
            karrayout=karrayout)
    
    @staticmethod
    def scaling_function(itype, nd, a, x):
        """
        nrange = scaling_function(itype, nd, a, x)
        
        
        Defined at poisson_isf.fpp lines 252-290
        
        Parameters
        ----------
        itype : int
        nd : int
        a : float array
        x : float array
        
        Returns
        -------
        nrange : int
        
        """
        nrange = _AresMainPy.f90wrap_scaling_function(itype=itype, nd=nd, a=a, x=x)
        return nrange
    
    @staticmethod
    def back_trans_8(nd, nt, x, y):
        """
        back_trans_8(nd, nt, x, y)
        
        
        Defined at poisson_isf.fpp lines 293-383
        
        Parameters
        ----------
        nd : int
        nt : int
        x : float array
        y : float array
        
        """
        _AresMainPy.f90wrap_back_trans_8(nd=nd, nt=nt, x=x, y=y)
    
    @staticmethod
    def gequad(n_gauss, p_gauss, w_gauss):
        """
        ur_gauss, dr_gauss, acc_gauss = gequad(n_gauss, p_gauss, w_gauss)
        
        
        Defined at poisson_isf.fpp lines 386-582
        
        Parameters
        ----------
        n_gauss : int
        p_gauss : float array
        w_gauss : float array
        
        Returns
        -------
        ur_gauss : float
        dr_gauss : float
        acc_gauss : float
        
        """
        ur_gauss, dr_gauss, acc_gauss = _AresMainPy.f90wrap_gequad(n_gauss=n_gauss, \
            p_gauss=p_gauss, w_gauss=w_gauss)
        return ur_gauss, dr_gauss, acc_gauss
    
    @staticmethod
    def calculate_dimensions(n01, n02, n03):
        """
        nfft1, nfft2, nfft3 = calculate_dimensions(n01, n02, n03)
        
        
        Defined at poisson_isf.fpp lines 585-619
        
        Parameters
        ----------
        n01 : int
        n02 : int
        n03 : int
        
        Returns
        -------
        nfft1 : int
        nfft2 : int
        nfft3 : int
        
        """
        nfft1, nfft2, nfft3 = _AresMainPy.f90wrap_calculate_dimensions(n01=n01, n02=n02, \
            n03=n03)
        return nfft1, nfft2, nfft3
    
    @staticmethod
    def psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot):
        """
        psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot)
        
        
        Defined at poisson_isf.fpp lines 622-663
        
        Parameters
        ----------
        n01 : int
        n02 : int
        n03 : int
        nfft1 : int
        nfft2 : int
        nfft3 : int
        hgrid : float
        karray : float array
        rhopot : float array
        
        """
        _AresMainPy.f90wrap_psolver_kernel(n01=n01, n02=n02, n03=n03, nfft1=nfft1, \
            nfft2=nfft2, nfft3=nfft3, hgrid=hgrid, karray=karray, rhopot=rhopot)
    
    @staticmethod
    def zarray_in(n01, n02, n03, nd1, nd2, nd3, density, zarray):
        """
        zarray_in(n01, n02, n03, nd1, nd2, nd3, density, zarray)
        
        
        Defined at poisson_isf.fpp lines 666-700
        
        Parameters
        ----------
        n01 : int
        n02 : int
        n03 : int
        nd1 : int
        nd2 : int
        nd3 : int
        density : float array
        zarray : float array
        
        """
        _AresMainPy.f90wrap_zarray_in(n01=n01, n02=n02, n03=n03, nd1=nd1, nd2=nd2, \
            nd3=nd3, density=density, zarray=zarray)
    
    @staticmethod
    def zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor):
        """
        zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor)
        
        
        Defined at poisson_isf.fpp lines 703-716
        
        Parameters
        ----------
        n01 : int
        n02 : int
        n03 : int
        nd1 : int
        nd2 : int
        nd3 : int
        rhopot : float array
        zarray : float array
        factor : float
        
        """
        _AresMainPy.f90wrap_zarray_out(n01=n01, n02=n02, n03=n03, nd1=nd1, nd2=nd2, \
            nd3=nd3, rhopot=rhopot, zarray=zarray, factor=factor)
    
    @staticmethod
    def scf_recursion(itype, n_iter, n_range, kernel_scf, kern_1_scf):
        """
        scf_recursion(itype, n_iter, n_range, kernel_scf, kern_1_scf)
        
        
        Defined at poisson_isf.fpp lines 729-745
        
        Parameters
        ----------
        itype : int
        n_iter : int
        n_range : int
        kernel_scf : float array
        kern_1_scf : float array
        
        """
        _AresMainPy.f90wrap_scf_recursion(itype=itype, n_iter=n_iter, n_range=n_range, \
            kernel_scf=kernel_scf, kern_1_scf=kern_1_scf)
    
    @staticmethod
    def scf_recursion_8(n_iter, n_range, kernel_scf, kern_1_scf):
        """
        scf_recursion_8(n_iter, n_range, kernel_scf, kern_1_scf)
        
        
        Defined at poisson_isf.fpp lines 758-852
        
        Parameters
        ----------
        n_iter : int
        n_range : int
        kernel_scf : float array
        kern_1_scf : float array
        
        """
        _AresMainPy.f90wrap_scf_recursion_8(n_iter=n_iter, n_range=n_range, \
            kernel_scf=kernel_scf, kern_1_scf=kern_1_scf)
    
    @staticmethod
    def karrayhalf_in(n01, n02, n03, n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, \
        nd3, kernel, karrayhalf):
        """
        karrayhalf_in(n01, n02, n03, n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, \
            kernel, karrayhalf)
        
        
        Defined at poisson_isf.fpp lines 864-905
        
        Parameters
        ----------
        n01 : int
        n02 : int
        n03 : int
        n1k : int
        n2k : int
        n3k : int
        nfft1 : int
        nfft2 : int
        nfft3 : int
        nd1 : int
        nd2 : int
        nd3 : int
        kernel : float array
        karrayhalf : float array
        
        """
        _AresMainPy.f90wrap_karrayhalf_in(n01=n01, n02=n02, n03=n03, n1k=n1k, n2k=n2k, \
            n3k=n3k, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, nd1=nd1, nd2=nd2, nd3=nd3, \
            kernel=kernel, karrayhalf=karrayhalf)
    
    @staticmethod
    def kernel_recon(n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, zarray, \
        karray):
        """
        kernel_recon(n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, zarray, karray)
        
        
        Defined at poisson_isf.fpp lines 919-969
        
        Parameters
        ----------
        n1k : int
        n2k : int
        n3k : int
        nfft1 : int
        nfft2 : int
        nfft3 : int
        nd1 : int
        nd2 : int
        nd3 : int
        zarray : float array
        karray : float array
        
        """
        _AresMainPy.f90wrap_kernel_recon(n1k=n1k, n2k=n2k, n3k=n3k, nfft1=nfft1, \
            nfft2=nfft2, nfft3=nfft3, nd1=nd1, nd2=nd2, nd3=nd3, zarray=zarray, \
            karray=karray)
    
    @staticmethod
    def norm_ind(nd1, nd2, nd3, i1, i2, i3, ind):
        """
        norm_ind(nd1, nd2, nd3, i1, i2, i3, ind)
        
        
        Defined at poisson_isf.fpp lines 981-1001
        
        Parameters
        ----------
        nd1 : int
        nd2 : int
        nd3 : int
        i1 : int
        i2 : int
        i3 : int
        ind : int
        
        """
        _AresMainPy.f90wrap_norm_ind(nd1=nd1, nd2=nd2, nd3=nd3, i1=i1, i2=i2, i3=i3, \
            ind=ind)
    
    @staticmethod
    def symm_ind(nd1, nd2, nd3, i1, i2, i3, ind):
        """
        symm_ind(nd1, nd2, nd3, i1, i2, i3, ind)
        
        
        Defined at poisson_isf.fpp lines 1014-1033
        
        Parameters
        ----------
        nd1 : int
        nd2 : int
        nd3 : int
        i1 : int
        i2 : int
        i3 : int
        ind : int
        
        """
        _AresMainPy.f90wrap_symm_ind(nd1=nd1, nd2=nd2, nd3=nd3, i1=i1, i2=i2, i3=i3, \
            ind=ind)
    
    @staticmethod
    def kernel_application(n1xy, n2xy, n3xy, nd1h, nd2, nd3, nfft1, nfft2, nfft3, \
        zarray, karray, inzee):
        """
        kernel_application(n1xy, n2xy, n3xy, nd1h, nd2, nd3, nfft1, nfft2, nfft3, \
            zarray, karray, inzee)
        
        
        Defined at poisson_isf.fpp lines 1037-1429
        
        Parameters
        ----------
        n1xy : int
        n2xy : int
        n3xy : int
        nd1h : int
        nd2 : int
        nd3 : int
        nfft1 : int
        nfft2 : int
        nfft3 : int
        zarray : float array
        karray : float array
        inzee : int
        
        --------------------------------------------
        --- Starting reconstruction half -> full ---
        --------------------------------------------
        -------------Case i3 = 1
        """
        _AresMainPy.f90wrap_kernel_application(n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, \
            nd1h=nd1h, nd2=nd2, nd3=nd3, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, \
            zarray=zarray, karray=karray, inzee=inzee)
    
    @property
    def karray(self):
        """
        Element karray ftype=real(dp) pytype=float
        
        
        Defined at poisson_isf.fpp line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_poisson_isf__array__karray(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            karray = self._arrays[array_handle]
        else:
            karray = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_poisson_isf__array__karray)
            self._arrays[array_handle] = karray
        return karray
    
    @karray.setter
    def karray(self, karray):
        self.karray[...] = karray
    
    def __str__(self):
        ret = ['<poisson_isf>{\n']
        ret.append('    karray : ')
        ret.append(repr(self.karray))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

poisson_isf = Poisson_Isf()

class Begin_Module(f90wrap.runtime.FortranModule):
    """
    Module begin_module
    
    
    Defined at Begin_module.fpp lines 5-650
    
    """
    @staticmethod
    def initial_grid():
        """
        initial_grid()
        
        
        Defined at Begin_module.fpp lines 11-18
        
        
        """
        _AresMainPy.f90wrap_initial_grid()
    
    @staticmethod
    def initial_grid_per():
        """
        initial_grid_per()
        
        
        Defined at Begin_module.fpp lines 21-80
        
        
        """
        _AresMainPy.f90wrap_initial_grid_per()
    
    @staticmethod
    def initial_grid_iso():
        """
        initial_grid_iso()
        
        
        Defined at Begin_module.fpp lines 83-143
        
        
        """
        _AresMainPy.f90wrap_initial_grid_iso()
    
    @staticmethod
    def initial_density():
        """
        initial_density()
        
        
        Defined at Begin_module.fpp lines 146-185
        
        
        =================================================
        SHIELD
        """
        _AresMainPy.f90wrap_initial_density()
    
    @staticmethod
    def inichrg():
        """
        inichrg()
        
        
        Defined at Begin_module.fpp lines 188-375
        
        
        """
        _AresMainPy.f90wrap_inichrg()
    
    @staticmethod
    def init_iso_density():
        """
        init_iso_density()
        
        
        Defined at Begin_module.fpp lines 378-405
        
        
        =====================================================================
        ##JUDGE THE INITIALIZE METHOD
        """
        _AresMainPy.f90wrap_init_iso_density()
    
    @staticmethod
    def inichrg_iso():
        """
        inichrg_iso()
        
        
        Defined at Begin_module.fpp lines 408-503
        
        
        """
        _AresMainPy.f90wrap_inichrg_iso()
    
    @staticmethod
    def init_mo_density_iso():
        """
        init_mo_density_iso()
        
        
        Defined at Begin_module.fpp lines 506-569
        
        
        ===============================================================
        ##assignment to grid
        """
        _AresMainPy.f90wrap_init_mo_density_iso()
    
    @staticmethod
    def iso_init_mo(initx_sto):
        """
        iso_init_mo(initx_sto)
        
        
        Defined at Begin_module.fpp lines 572-648
        
        Parameters
        ----------
        initx_sto : float array
        
        """
        _AresMainPy.f90wrap_iso_init_mo(initx_sto=initx_sto)
    
    _dt_array_initialisers = []
    

begin_module = Begin_Module()

class Smearing_Module(f90wrap.runtime.FortranModule):
    """
    Module smearing_module
    
    
    Defined at Smearing_module.fpp lines 5-380
    
    """
    @staticmethod
    def smear_init(nev):
        """
        smear_init(nev)
        
        
        Defined at Smearing_module.fpp lines 21-48
        
        Parameters
        ----------
        nev : int
        
        """
        _AresMainPy.f90wrap_smear_init(nev=nev)
    
    @staticmethod
    def destroy_smear():
        """
        destroy_smear()
        
        
        Defined at Smearing_module.fpp lines 51-63
        
        
        """
        _AresMainPy.f90wrap_destroy_smear()
    
    @staticmethod
    def hp(x, n):
        """
        hp = hp(x, n)
        
        
        Defined at Smearing_module.fpp lines 66-96
        
        Parameters
        ----------
        x : float
        n : int
        
        Returns
        -------
        hp : float
        
        """
        hp = _AresMainPy.f90wrap_hp(x=x, n=n)
        return hp
    
    @staticmethod
    def smearsn(y, y0, w):
        """
        sn = smearsn(y, y0, w)
        
        
        Defined at Smearing_module.fpp lines 99-137
        
        Parameters
        ----------
        y : float
        y0 : float
        w : float
        
        Returns
        -------
        sn : float
        
        """
        sn = _AresMainPy.f90wrap_smearsn(y=y, y0=y0, w=w)
        return sn
    
    @staticmethod
    def fermilevel(ne, nev, nk, wk, eval, wke):
        """
        fme, ets = fermilevel(ne, nev, nk, wk, eval, wke)
        
        
        Defined at Smearing_module.fpp lines 140-229
        
        Parameters
        ----------
        ne : float
        nev : int
        nk : int
        wk : float array
        eval : float array
        wke : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_fermilevel(ne=ne, nev=nev, nk=nk, wk=wk, \
            eval=eval, wke=wke)
        return fme, ets
    
    @staticmethod
    def fermilevel_iso(ne, nev, nk, wk, eval, wke):
        """
        fme, ets = fermilevel_iso(ne, nev, nk, wk, eval, wke)
        
        
        Defined at Smearing_module.fpp lines 232-327
        
        Parameters
        ----------
        ne : float
        nev : int
        nk : int
        wk : float array
        eval : float array
        wke : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_fermilevel_iso(ne=ne, nev=nev, nk=nk, wk=wk, \
            eval=eval, wke=wke)
        return fme, ets
    
    @staticmethod
    def enpy(e_mu, fi):
        """
        enpy = enpy(e_mu, fi)
        
        
        Defined at Smearing_module.fpp lines 330-349
        
        Parameters
        ----------
        e_mu : float
        fi : float
        
        Returns
        -------
        enpy : float
        
        """
        enpy = _AresMainPy.f90wrap_enpy(e_mu=e_mu, fi=fi)
        return enpy
    
    @staticmethod
    def whg(x, n):
        """
        whg = whg(x, n)
        
        
        Defined at Smearing_module.fpp lines 352-379
        
        Parameters
        ----------
        x : float
        n : int
        
        Returns
        -------
        whg : float
        
        """
        whg = _AresMainPy.f90wrap_whg(x=x, n=n)
        return whg
    
    @property
    def wke(self):
        """
        Element wke ftype=real(dp) pytype=float
        
        
        Defined at Smearing_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_smearing_module__array__wke(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wke = self._arrays[array_handle]
        else:
            wke = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_smearing_module__array__wke)
            self._arrays[array_handle] = wke
        return wke
    
    @wke.setter
    def wke(self, wke):
        self.wke[...] = wke
    
    @property
    def subwke(self):
        """
        Element subwke ftype=real(dp) pytype=float
        
        
        Defined at Smearing_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_smearing_module__array__subwke(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            subwke = self._arrays[array_handle]
        else:
            subwke = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_smearing_module__array__subwke)
            self._arrays[array_handle] = subwke
        return subwke
    
    @subwke.setter
    def subwke(self, subwke):
        self.subwke[...] = subwke
    
    def __str__(self):
        ret = ['<smearing_module>{\n']
        ret.append('    wke : ')
        ret.append(repr(self.wke))
        ret.append(',\n    subwke : ')
        ret.append(repr(self.subwke))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

smearing_module = Smearing_Module()

class Nlpot_Module(f90wrap.runtime.FortranModule):
    """
    Module nlpot_module
    
    
    Defined at Nonlocalpot_module.fpp lines 5-1371
    
    """
    @f90wrap.runtime.register_class("AresMainPy.nol_type")
    class nol_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=nol_type)
        
        
        Defined at Nonlocalpot_module.fpp lines 14-27
        
        """
        def __init__(self, handle=None):
            """
            self = Nol_Type()
            
            
            Defined at Nonlocalpot_module.fpp lines 14-27
            
            
            Returns
            -------
            this : Nol_Type
            	Object to be constructed
            
            
            Automatically generated constructor for nol_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_nol_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Nol_Type
            
            
            Defined at Nonlocalpot_module.fpp lines 14-27
            
            Parameters
            ----------
            this : Nol_Type
            	Object to be destructed
            
            
            Automatically generated destructor for nol_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_nol_type_finalise(this=self._handle)
        
        @property
        def npts(self):
            """
            Element npts ftype=integer(i4b) pytype=int
            
            
            Defined at Nonlocalpot_module.fpp line 15
            
            """
            return _AresMainPy.f90wrap_nol_type__get__npts(self._handle)
        
        @npts.setter
        def npts(self, npts):
            _AresMainPy.f90wrap_nol_type__set__npts(self._handle, npts)
        
        @property
        def id(self):
            """
            Element id ftype=integer(i4b) pytype=int
            
            
            Defined at Nonlocalpot_module.fpp line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__id(self._handle)
            if array_handle in self._arrays:
                id = self._arrays[array_handle]
            else:
                id = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__id)
                self._arrays[array_handle] = id
            return id
        
        @id.setter
        def id(self, id):
            self.id[...] = id
        
        @property
        def rrvec(self):
            """
            Element rrvec ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__rrvec(self._handle)
            if array_handle in self._arrays:
                rrvec = self._arrays[array_handle]
            else:
                rrvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__rrvec)
                self._arrays[array_handle] = rrvec
            return rrvec
        
        @rrvec.setter
        def rrvec(self, rrvec):
            self.rrvec[...] = rrvec
        
        @property
        def proj0(self):
            """
            Element proj0 ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.fpp line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__proj0(self._handle)
            if array_handle in self._arrays:
                proj0 = self._arrays[array_handle]
            else:
                proj0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__proj0)
                self._arrays[array_handle] = proj0
            return proj0
        
        @proj0.setter
        def proj0(self, proj0):
            self.proj0[...] = proj0
        
        @property
        def proj(self):
            """
            Element proj ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.fpp line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__proj(self._handle)
            if array_handle in self._arrays:
                proj = self._arrays[array_handle]
            else:
                proj = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__proj)
                self._arrays[array_handle] = proj
            return proj
        
        @proj.setter
        def proj(self, proj):
            self.proj[...] = proj
        
        @property
        def proj_phs(self):
            """
            Element proj_phs ftype=complex(dp) pytype=complex
            
            
            Defined at Nonlocalpot_module.fpp line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__proj_phs(self._handle)
            if array_handle in self._arrays:
                proj_phs = self._arrays[array_handle]
            else:
                proj_phs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__proj_phs)
                self._arrays[array_handle] = proj_phs
            return proj_phs
        
        @proj_phs.setter
        def proj_phs(self, proj_phs):
            self.proj_phs[...] = proj_phs
        
        @property
        def proj0_dg(self):
            """
            Element proj0_dg ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.fpp line 22
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__proj0_dg(self._handle)
            if array_handle in self._arrays:
                proj0_dg = self._arrays[array_handle]
            else:
                proj0_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__proj0_dg)
                self._arrays[array_handle] = proj0_dg
            return proj0_dg
        
        @proj0_dg.setter
        def proj0_dg(self, proj0_dg):
            self.proj0_dg[...] = proj0_dg
        
        @property
        def proj_dg(self):
            """
            Element proj_dg ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.fpp line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__proj_dg(self._handle)
            if array_handle in self._arrays:
                proj_dg = self._arrays[array_handle]
            else:
                proj_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__proj_dg)
                self._arrays[array_handle] = proj_dg
            return proj_dg
        
        @proj_dg.setter
        def proj_dg(self, proj_dg):
            self.proj_dg[...] = proj_dg
        
        @property
        def proj_phs_dg(self):
            """
            Element proj_phs_dg ftype=complex(dp) pytype=complex
            
            
            Defined at Nonlocalpot_module.fpp line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__proj_phs_dg(self._handle)
            if array_handle in self._arrays:
                proj_phs_dg = self._arrays[array_handle]
            else:
                proj_phs_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__proj_phs_dg)
                self._arrays[array_handle] = proj_phs_dg
            return proj_phs_dg
        
        @proj_phs_dg.setter
        def proj_phs_dg(self, proj_phs_dg):
            self.proj_phs_dg[...] = proj_phs_dg
        
        @property
        def id_iso(self):
            """
            Element id_iso ftype=integer(i4b) pytype=int
            
            
            Defined at Nonlocalpot_module.fpp line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__id_iso(self._handle)
            if array_handle in self._arrays:
                id_iso = self._arrays[array_handle]
            else:
                id_iso = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__id_iso)
                self._arrays[array_handle] = id_iso
            return id_iso
        
        @id_iso.setter
        def id_iso(self, id_iso):
            self.id_iso[...] = id_iso
        
        @property
        def rrvec_iso(self):
            """
            Element rrvec_iso ftype=real(dp) pytype=float
            
            
            Defined at Nonlocalpot_module.fpp line 27
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_nol_type__array__rrvec_iso(self._handle)
            if array_handle in self._arrays:
                rrvec_iso = self._arrays[array_handle]
            else:
                rrvec_iso = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_nol_type__array__rrvec_iso)
                self._arrays[array_handle] = rrvec_iso
            return rrvec_iso
        
        @rrvec_iso.setter
        def rrvec_iso(self, rrvec_iso):
            self.rrvec_iso[...] = rrvec_iso
        
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
            ret.append(',\n    id_iso : ')
            ret.append(repr(self.id_iso))
            ret.append(',\n    rrvec_iso : ')
            ret.append(repr(self.rrvec_iso))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def initialize_nlpot():
        """
        initialize_nlpot()
        
        
        Defined at Nonlocalpot_module.fpp lines 37-44
        
        
        """
        _AresMainPy.f90wrap_initialize_nlpot()
    
    @staticmethod
    def initialize_nlpot_per():
        """
        initialize_nlpot_per()
        
        
        Defined at Nonlocalpot_module.fpp lines 47-191
        
        
        ===================store the data================
        """
        _AresMainPy.f90wrap_initialize_nlpot_per()
    
    @staticmethod
    def initialize_nlpot_iso():
        """
        initialize_nlpot_iso()
        
        
        Defined at Nonlocalpot_module.fpp lines 194-326
        
        
        ===================store the data================
        """
        _AresMainPy.f90wrap_initialize_nlpot_iso()
    
    @staticmethod
    def destroy_nlpot():
        """
        destroy_nlpot()
        
        
        Defined at Nonlocalpot_module.fpp lines 329-385
        
        
        """
        _AresMainPy.f90wrap_destroy_nlpot()
    
    @staticmethod
    def set_beta_real():
        """
        set_beta_real()
        
        
        Defined at Nonlocalpot_module.fpp lines 388-451
        
        
        ======================
        ##PLOT
        """
        _AresMainPy.f90wrap_set_beta_real()
    
    @staticmethod
    def nlp_beta_interp_r(ity, ia, beta_init):
        """
        nlp_beta_interp_r(ity, ia, beta_init)
        
        
        Defined at Nonlocalpot_module.fpp lines 454-516
        
        Parameters
        ----------
        ity : int
        ia : int
        beta_init : float array
        
        """
        _AresMainPy.f90wrap_nlp_beta_interp_r(ity=ity, ia=ia, beta_init=beta_init)
    
    @staticmethod
    def nlp_beta_ylm_r(ity, ia, beta_ylm):
        """
        nlp_beta_ylm_r(ity, ia, beta_ylm)
        
        
        Defined at Nonlocalpot_module.fpp lines 519-561
        
        Parameters
        ----------
        ity : int
        ia : int
        beta_ylm : float array
        
        """
        _AresMainPy.f90wrap_nlp_beta_ylm_r(ity=ity, ia=ia, beta_ylm=beta_ylm)
    
    @staticmethod
    def nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase):
        """
        nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase)
        
        
        Defined at Nonlocalpot_module.fpp lines 564-591
        
        Parameters
        ----------
        ity : int
        ia : int
        beta_ylm : float array
        beta_phase : complex array
        
        """
        _AresMainPy.f90wrap_nlp_beta_phase_r(ity=ity, ia=ia, beta_ylm=beta_ylm, \
            beta_phase=beta_phase)
    
    @staticmethod
    def apply_ylm(l, m, fac, x, y, z, f):
        """
        apply_ylm(l, m, fac, x, y, z, f)
        
        
        Defined at Nonlocalpot_module.fpp lines 594-670
        
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
        _AresMainPy.f90wrap_apply_ylm(l=l, m=m, fac=fac, x=x, y=y, z=z, f=f)
    
    @staticmethod
    def nlp_wijk_dg(ndg, inpol, numdg, lft, rit, wijk, drijk):
        """
        nlp_wijk_dg(ndg, inpol, numdg, lft, rit, wijk, drijk)
        
        
        Defined at Nonlocalpot_module.fpp lines 673-786
        
        Parameters
        ----------
        ndg : int
        inpol : int
        numdg : int
        lft : int
        rit : int
        wijk : float array
        drijk : float array
        
        """
        _AresMainPy.f90wrap_nlp_wijk_dg(ndg=ndg, inpol=inpol, numdg=numdg, lft=lft, \
            rit=rit, wijk=wijk, drijk=drijk)
    
    @staticmethod
    def nlp_wijk_dg_old(ndg, n_near, numdg, lft, rit, wijk, drijk):
        """
        nlp_wijk_dg_old(ndg, n_near, numdg, lft, rit, wijk, drijk)
        
        
        Defined at Nonlocalpot_module.fpp lines 788-891
        
        Parameters
        ----------
        ndg : int
        n_near : int
        numdg : int
        lft : int
        rit : int
        wijk : float array
        drijk : float array
        
        """
        _AresMainPy.f90wrap_nlp_wijk_dg_old(ndg=ndg, n_near=n_near, numdg=numdg, \
            lft=lft, rit=rit, wijk=wijk, drijk=drijk)
    
    @staticmethod
    def nlp_interp_betalm_dg(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, \
        betalm_ijk):
        """
        nlp_interp_betalm_dg(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, \
            betalm_ijk)
        
        
        Defined at Nonlocalpot_module.fpp lines 894-995
        
        Parameters
        ----------
        ity : int
        ia : int
        ip : int
        numdg : int
        lft : int
        rit : int
        wijk : float array
        drijk : float array
        beta0_ijk : float array
        betalm_ijk : float array
        
        =============================
        ##FOR ISO
        """
        _AresMainPy.f90wrap_nlp_interp_betalm_dg(ity=ity, ia=ia, ip=ip, numdg=numdg, \
            lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0_ijk=beta0_ijk, \
            betalm_ijk=betalm_ijk)
    
    @staticmethod
    def nlp_interp_betalm_dg_iso(ity, ia, ip, numdg, lft, rit, wijk, drijk, \
        beta0_ijk, betalm_ijk):
        """
        nlp_interp_betalm_dg_iso(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, \
            betalm_ijk)
        
        
        Defined at Nonlocalpot_module.fpp lines 998-1104
        
        Parameters
        ----------
        ity : int
        ia : int
        ip : int
        numdg : int
        lft : int
        rit : int
        wijk : float array
        drijk : float array
        beta0_ijk : float array
        betalm_ijk : float array
        
        =============================
        ##FOR ISO
        """
        _AresMainPy.f90wrap_nlp_interp_betalm_dg_iso(ity=ity, ia=ia, ip=ip, numdg=numdg, \
            lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0_ijk=beta0_ijk, \
            betalm_ijk=betalm_ijk)
    
    @staticmethod
    def nlp_beta_ylm_r_dg(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm):
        """
        nlp_beta_ylm_r_dg(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm)
        
        
        Defined at Nonlocalpot_module.fpp lines 1107-1122
        
        Parameters
        ----------
        ity : int
        ia : int
        numdg : int
        lft : int
        rit : int
        wijk : float array
        drijk : float array
        beta0 : float array
        beta_ylm : float array
        
        """
        _AresMainPy.f90wrap_nlp_beta_ylm_r_dg(ity=ity, ia=ia, numdg=numdg, lft=lft, \
            rit=rit, wijk=wijk, drijk=drijk, beta0=beta0, beta_ylm=beta_ylm)
    
    @staticmethod
    def nlp_beta_ylm_r_dg_iso(ity, ia, numdg, lft, rit, wijk, drijk, beta0, \
        beta_ylm):
        """
        nlp_beta_ylm_r_dg_iso(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm)
        
        
        Defined at Nonlocalpot_module.fpp lines 1125-1161
        
        Parameters
        ----------
        ity : int
        ia : int
        numdg : int
        lft : int
        rit : int
        wijk : float array
        drijk : float array
        beta0 : float array
        beta_ylm : float array
        
        """
        _AresMainPy.f90wrap_nlp_beta_ylm_r_dg_iso(ity=ity, ia=ia, numdg=numdg, lft=lft, \
            rit=rit, wijk=wijk, drijk=drijk, beta0=beta0, beta_ylm=beta_ylm)
    
    @staticmethod
    def set_beta_real_dg():
        """
        set_beta_real_dg()
        
        
        Defined at Nonlocalpot_module.fpp lines 1164-1220
        
        
        """
        _AresMainPy.f90wrap_set_beta_real_dg()
    
    @staticmethod
    def initialize_nlpot_band():
        """
        initialize_nlpot_band()
        
        
        Defined at Nonlocalpot_module.fpp lines 1222-1370
        
        
        ===================store the data================
        """
        _AresMainPy.f90wrap_initialize_nlpot_band()
    
    @property
    def max_nlnpts(self):
        """
        Element max_nlnpts ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.fpp line 32
        
        """
        return _AresMainPy.f90wrap_nlpot_module__get__max_nlnpts()
    
    @max_nlnpts.setter
    def max_nlnpts(self, max_nlnpts):
        _AresMainPy.f90wrap_nlpot_module__set__max_nlnpts(max_nlnpts)
    
    def __str__(self):
        ret = ['<nlpot_module>{\n']
        ret.append('    max_nlnpts : ')
        ret.append(repr(self.max_nlnpts))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

nlpot_module = Nlpot_Module()

class Getvlocalpseudopotential(f90wrap.runtime.FortranModule):
    """
    Module getvlocalpseudopotential
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 5-940
    
    """
    @staticmethod
    def calvlpp():
        """
        calvlpp()
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 26-107
        
        
        ==========================================================
        """
        _AresMainPy.f90wrap_calvlpp()
    
    @staticmethod
    def ionpotentialassignment(ity, zion, poscar, temp):
        """
        ionpotentialassignment(ity, zion, poscar, temp)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 110-186
        
        Parameters
        ----------
        ity : int
        zion : float
        poscar : float array
        temp : float array
        
        ========================================================================
        """
        _AresMainPy.f90wrap_ionpotentialassignment(ity=ity, zion=zion, poscar=poscar, \
            temp=temp)
    
    @staticmethod
    def sphbess(l, x):
        """
        sphbess = sphbess(l, x)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 191-230
        
        Parameters
        ----------
        l : int
        x : float
        
        Returns
        -------
        sphbess : float
        
        """
        sphbess = _AresMainPy.f90wrap_sphbess(l=l, x=x)
        return sphbess
    
    @staticmethod
    def fourbess_gr(g, fg, r, fr):
        """
        fourbess_gr(g, fg, r, fr)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 234-257
        
        Parameters
        ----------
        g : float array
        fg : float array
        r : float array
        fr : float array
        
        """
        _AresMainPy.f90wrap_fourbess_gr(g=g, fg=fg, r=r, fr=fr)
    
    @staticmethod
    def dir2car_single(cry_coo, ort_coo, lat):
        """
        dir2car_single(cry_coo, ort_coo, lat)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 598-605
        
        Parameters
        ----------
        cry_coo : float array
        ort_coo : float array
        lat : float array
        
        """
        _AresMainPy.f90wrap_dir2car_single(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)
    
    @staticmethod
    def ionpotentialassignment_dg(ity, zion, poscar, temp):
        """
        ionpotentialassignment_dg(ity, zion, poscar, temp)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 800-866
        
        Parameters
        ----------
        ity : int
        zion : float
        poscar : float array
        temp : float array
        
        ========================================================================
        ==========================================================xlt test
        """
        _AresMainPy.f90wrap_ionpotentialassignment_dg(ity=ity, zion=zion, poscar=poscar, \
            temp=temp)
    
    @staticmethod
    def set_vloc_dg(ity, xyz):
        """
        vloc = set_vloc_dg(ity, xyz)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 869-915
        
        Parameters
        ----------
        ity : int
        xyz : float array
        
        Returns
        -------
        vloc : float
        
        """
        vloc = _AresMainPy.f90wrap_set_vloc_dg(ity=ity, xyz=xyz)
        return vloc
    
    @staticmethod
    def set_vcomp(n_inp, zion, rgauss, r_inp, v_inp, vcomp):
        """
        set_vcomp(n_inp, zion, rgauss, r_inp, v_inp, vcomp)
        
        
        Defined at IonLocalPotentialAssignment.fpp lines 918-938
        
        Parameters
        ----------
        n_inp : int
        zion : int
        rgauss : float
        r_inp : float array
        v_inp : float array
        vcomp : float array
        
        """
        _AresMainPy.f90wrap_set_vcomp(n_inp=n_inp, zion=zion, rgauss=rgauss, \
            r_inp=r_inp, v_inp=v_inp, vcomp=vcomp)
    
    @property
    def radiusn(self):
        """
        Element radiusn ftype=integer(i4b) pytype=int
        
        
        Defined at IonLocalPotentialAssignment.fpp line 15
        
        """
        return _AresMainPy.f90wrap_getvlocalpseudopotential__get__radiusn()
    
    @radiusn.setter
    def radiusn(self, radiusn):
        _AresMainPy.f90wrap_getvlocalpseudopotential__set__radiusn(radiusn)
    
    @property
    def fr(self):
        """
        Element fr ftype=real(dp) pytype=float
        
        
        Defined at IonLocalPotentialAssignment.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_getvlocalpseudopotential__array__fr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            fr = self._arrays[array_handle]
        else:
            fr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_getvlocalpseudopotential__array__fr)
            self._arrays[array_handle] = fr
        return fr
    
    @fr.setter
    def fr(self, fr):
        self.fr[...] = fr
    
    @property
    def wijk(self):
        """
        Element wijk ftype=real(dp) pytype=float
        
        
        Defined at IonLocalPotentialAssignment.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_getvlocalpseudopotential__array__wijk(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            wijk = self._arrays[array_handle]
        else:
            wijk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_getvlocalpseudopotential__array__wijk)
            self._arrays[array_handle] = wijk
        return wijk
    
    @wijk.setter
    def wijk(self, wijk):
        self.wijk[...] = wijk
    
    @property
    def drijk(self):
        """
        Element drijk ftype=real(dp) pytype=float
        
        
        Defined at IonLocalPotentialAssignment.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_getvlocalpseudopotential__array__drijk(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            drijk = self._arrays[array_handle]
        else:
            drijk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_getvlocalpseudopotential__array__drijk)
            self._arrays[array_handle] = drijk
        return drijk
    
    @drijk.setter
    def drijk(self, drijk):
        self.drijk[...] = drijk
    
    @property
    def ilft(self):
        """
        Element ilft ftype=integer(i4b) pytype=int
        
        
        Defined at IonLocalPotentialAssignment.fpp line 20
        
        """
        return _AresMainPy.f90wrap_getvlocalpseudopotential__get__ilft()
    
    @ilft.setter
    def ilft(self, ilft):
        _AresMainPy.f90wrap_getvlocalpseudopotential__set__ilft(ilft)
    
    @property
    def irit(self):
        """
        Element irit ftype=integer(i4b) pytype=int
        
        
        Defined at IonLocalPotentialAssignment.fpp line 20
        
        """
        return _AresMainPy.f90wrap_getvlocalpseudopotential__get__irit()
    
    @irit.setter
    def irit(self, irit):
        _AresMainPy.f90wrap_getvlocalpseudopotential__set__irit(irit)
    
    @property
    def numdg(self):
        """
        Element numdg ftype=integer(i4b) pytype=int
        
        
        Defined at IonLocalPotentialAssignment.fpp line 20
        
        """
        return _AresMainPy.f90wrap_getvlocalpseudopotential__get__numdg()
    
    @numdg.setter
    def numdg(self, numdg):
        _AresMainPy.f90wrap_getvlocalpseudopotential__set__numdg(numdg)
    
    def __str__(self):
        ret = ['<getvlocalpseudopotential>{\n']
        ret.append('    radiusn : ')
        ret.append(repr(self.radiusn))
        ret.append(',\n    fr : ')
        ret.append(repr(self.fr))
        ret.append(',\n    wijk : ')
        ret.append(repr(self.wijk))
        ret.append(',\n    drijk : ')
        ret.append(repr(self.drijk))
        ret.append(',\n    ilft : ')
        ret.append(repr(self.ilft))
        ret.append(',\n    irit : ')
        ret.append(repr(self.irit))
        ret.append(',\n    numdg : ')
        ret.append(repr(self.numdg))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

getvlocalpseudopotential = Getvlocalpseudopotential()

class Isolateset(f90wrap.runtime.FortranModule):
    """
    Module isolateset
    
    
    Defined at Isolate_module.fpp lines 5-2205
    
    """
    @staticmethod
    def isolatevcoulomb_a(rho_s, v_s):
        """
        isolatevcoulomb_a(rho_s, v_s)
        
        
        Defined at Isolate_module.fpp lines 25-93
        
        Parameters
        ----------
        rho_s : float array
        v_s : float array
        
        ===========================================================
        """
        _AresMainPy.f90wrap_isolatevcoulomb_a(rho_s=rho_s, v_s=v_s)
    
    @staticmethod
    def evaluateboundrygrid_b(rhos):
        """
        evaluateboundrygrid_b(rhos)
        
        
        Defined at Isolate_module.fpp lines 97-122
        
        Parameters
        ----------
        rhos : float array
        
        """
        _AresMainPy.f90wrap_evaluateboundrygrid_b(rhos=rhos)
    
    @staticmethod
    def endevaluate_b():
        """
        endevaluate_b()
        
        
        Defined at Isolate_module.fpp lines 126-137
        
        
        ------------
        """
        _AresMainPy.f90wrap_endevaluate_b()
    
    @staticmethod
    def calculateboundrypotential_b(rho, lmax):
        """
        calculateboundrypotential_b(rho, lmax)
        
        
        Defined at Isolate_module.fpp lines 141-177
        
        Parameters
        ----------
        rho : float array
        lmax : int
        
        """
        _AresMainPy.f90wrap_calculateboundrypotential_b(rho=rho, lmax=lmax)
    
    @staticmethod
    def laplaceequationslover_b(iter):
        """
        laplaceequationslover_b(iter)
        
        
        Defined at Isolate_module.fpp lines 181-249
        
        Parameters
        ----------
        iter : int
        
        ---------------------------------------------------------
        """
        _AresMainPy.f90wrap_laplaceequationslover_b(iter=iter)
    
    @staticmethod
    def vcoulombassignment_b(iter, niter_poisson, phi):
        """
        vcoulombassignment_b(iter, niter_poisson, phi)
        
        
        Defined at Isolate_module.fpp lines 253-263
        
        Parameters
        ----------
        iter : int
        niter_poisson : int
        phi : float array
        
        -------------------------------------------------
        """
        _AresMainPy.f90wrap_vcoulombassignment_b(iter=iter, niter_poisson=niter_poisson, \
            phi=phi)
    
    @staticmethod
    def cal_qlm(rho, lmax, q_l0, q_lm):
        """
        cal_qlm(rho, lmax, q_l0, q_lm)
        
        
        Defined at Isolate_module.fpp lines 267-310
        
        Parameters
        ----------
        rho : float array
        lmax : int
        q_l0 : float array
        q_lm : complex array
        
        """
        _AresMainPy.f90wrap_cal_qlm(rho=rho, lmax=lmax, q_l0=q_l0, q_lm=q_lm)
    
    @staticmethod
    def car2spe(orig, x, y, z):
        """
        r, cost, sint, cosp, sinp = car2spe(orig, x, y, z)
        
        
        Defined at Isolate_module.fpp lines 314-345
        
        Parameters
        ----------
        orig : float array
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
        r, cost, sint, cosp, sinp = _AresMainPy.f90wrap_car2spe(orig=orig, x=x, y=y, \
            z=z)
        return r, cost, sint, cosp, sinp
    
    @staticmethod
    def cal_plm(lmax, x, z, plm):
        """
        cal_plm(lmax, x, z, plm)
        
        
        Defined at Isolate_module.fpp lines 349-377
        
        Parameters
        ----------
        lmax : int
        x : float
        z : float
        plm : float array
        
        """
        _AresMainPy.f90wrap_cal_plm(lmax=lmax, x=x, z=z, plm=plm)
    
    @staticmethod
    def cal1pot(iex, iey, iez, orig, lmax, q_l0, q_lm):
        """
        phi = cal1pot(iex, iey, iez, orig, lmax, q_l0, q_lm)
        
        
        Defined at Isolate_module.fpp lines 381-411
        
        Parameters
        ----------
        iex : int
        iey : int
        iez : int
        orig : float array
        lmax : int
        q_l0 : float array
        q_lm : complex array
        
        Returns
        -------
        phi : float
        
        """
        phi = _AresMainPy.f90wrap_cal1pot(iex=iex, iey=iey, iez=iez, orig=orig, \
            lmax=lmax, q_l0=q_l0, q_lm=q_lm)
        return phi
    
    @staticmethod
    def apply_boundary(lmax, orig, q_l0, q_lm, bphi):
        """
        apply_boundary(lmax, orig, q_l0, q_lm, bphi)
        
        
        Defined at Isolate_module.fpp lines 415-446
        
        Parameters
        ----------
        lmax : int
        orig : float array
        q_l0 : float array
        q_lm : complex array
        bphi : float array
        
        """
        _AresMainPy.f90wrap_apply_boundary(lmax=lmax, orig=orig, q_l0=q_l0, q_lm=q_lm, \
            bphi=bphi)
    
    @staticmethod
    def apply_boundary00(rho, bphi):
        """
        apply_boundary00(rho, bphi)
        
        
        Defined at Isolate_module.fpp lines 449-483
        
        Parameters
        ----------
        rho : float array
        bphi : float array
        
        """
        _AresMainPy.f90wrap_apply_boundary00(rho=rho, bphi=bphi)
    
    @staticmethod
    def laplb(ford, bf, f):
        """
        laplb(ford, bf, f)
        
        
        Defined at Isolate_module.fpp lines 486-499
        
        Parameters
        ----------
        ford : int
        bf : float array
        f : float array
        
        """
        _AresMainPy.f90wrap_laplb(ford=ford, bf=bf, f=f)
    
    @staticmethod
    def lapla(ford, dxyz, f, af):
        """
        lapla(ford, dxyz, f, af)
        
        
        Defined at Isolate_module.fpp lines 503-520
        
        Parameters
        ----------
        ford : int
        dxyz : float array
        f : float array
        af : float array
        
        """
        _AresMainPy.f90wrap_lapla(ford=ford, dxyz=dxyz, f=f, af=af)
    
    @staticmethod
    def lapla_coe(ford, dxyz, coe):
        """
        lapla_coe(ford, dxyz, coe)
        
        
        Defined at Isolate_module.fpp lines 524-537
        
        Parameters
        ----------
        ford : int
        dxyz : float array
        coe : float array
        
        """
        _AresMainPy.f90wrap_lapla_coe(ford=ford, dxyz=dxyz, coe=coe)
    
    @staticmethod
    def damped_jacobi_iterate(ford, dxyz, func, rhs):
        """
        damped_jacobi_iterate(ford, dxyz, func, rhs)
        
        
        Defined at Isolate_module.fpp lines 541-595
        
        Parameters
        ----------
        ford : int
        dxyz : float array
        func : float array
        rhs : float array
        
        """
        _AresMainPy.f90wrap_damped_jacobi_iterate(ford=ford, dxyz=dxyz, func=func, \
            rhs=rhs)
    
    @staticmethod
    def gauss_seidel_iterate(ford, dxyz, func, rhs):
        """
        gauss_seidel_iterate(ford, dxyz, func, rhs)
        
        
        Defined at Isolate_module.fpp lines 600-656
        
        Parameters
        ----------
        ford : int
        dxyz : float array
        func : float array
        rhs : float array
        
        """
        _AresMainPy.f90wrap_gauss_seidel_iterate(ford=ford, dxyz=dxyz, func=func, \
            rhs=rhs)
    
    @staticmethod
    def refine(var, nxyz, nxyz_fine, var_fine):
        """
        refine(var, nxyz, nxyz_fine, var_fine)
        
        
        Defined at Isolate_module.fpp lines 661-674
        
        Parameters
        ----------
        var : float array
        nxyz : int array
        nxyz_fine : int array
        var_fine : float array
        
        """
        _AresMainPy.f90wrap_refine(var=var, nxyz=nxyz, nxyz_fine=nxyz_fine, \
            var_fine=var_fine)
    
    @staticmethod
    def interpolate(var1, nxyz1, nxyz2, var2):
        """
        interpolate(var1, nxyz1, nxyz2, var2)
        
        
        Defined at Isolate_module.fpp lines 679-741
        
        Parameters
        ----------
        var1 : float array
        nxyz1 : int array
        nxyz2 : int array
        var2 : float array
        
        -----------------------------------------------------------------------
        """
        _AresMainPy.f90wrap_interpolate(var1=var1, nxyz1=nxyz1, nxyz2=nxyz2, var2=var2)
    
    @staticmethod
    def smooth(f):
        """
        smooth(f)
        
        
        Defined at Isolate_module.fpp lines 746-797
        
        Parameters
        ----------
        f : float array
        
        """
        _AresMainPy.f90wrap_smooth(f=f)
    
    @staticmethod
    def restrict(var, nxyz, nxyz_coarse, var_coarse):
        """
        restrict(var, nxyz, nxyz_coarse, var_coarse)
        
        
        Defined at Isolate_module.fpp lines 802-821
        
        Parameters
        ----------
        var : float array
        nxyz : int array
        nxyz_coarse : int array
        var_coarse : float array
        
        """
        _AresMainPy.f90wrap_restrict(var=var, nxyz=nxyz, nxyz_coarse=nxyz_coarse, \
            var_coarse=var_coarse)
    
    @staticmethod
    def residual(ford, dxyz, f, rhs, res):
        """
        residual(ford, dxyz, f, rhs, res)
        
        
        Defined at Isolate_module.fpp lines 826-844
        
        Parameters
        ----------
        ford : int
        dxyz : float array
        f : float array
        rhs : float array
        res : float array
        
        """
        _AresMainPy.f90wrap_residual(ford=ford, dxyz=dxyz, f=f, rhs=rhs, res=res)
    
    @staticmethod
    def v_cycle(ford, f, rhs, dxyz, nvc):
        """
        v_cycle(ford, f, rhs, dxyz, nvc)
        
        
        Defined at Isolate_module.fpp lines 849-965
        
        Parameters
        ----------
        ford : int
        f : float array
        rhs : float array
        dxyz : float array
        nvc : int
        
        """
        _AresMainPy.f90wrap_v_cycle(ford=ford, f=f, rhs=rhs, dxyz=dxyz, nvc=nvc)
    
    @staticmethod
    def csix_coe(dxyz, coe):
        """
        csix_coe(dxyz, coe)
        
        
        Defined at Isolate_module.fpp lines 970-983
        
        Parameters
        ----------
        dxyz : float array
        coe : float array
        
        """
        _AresMainPy.f90wrap_csix_coe(dxyz=dxyz, coe=coe)
    
    @staticmethod
    def lapla_csix(f, af, coe):
        """
        lapla_csix(f, af, coe)
        
        
        Defined at Isolate_module.fpp lines 987-1021
        
        Parameters
        ----------
        f : float array
        af : float array
        coe : float array
        
        """
        _AresMainPy.f90wrap_lapla_csix(f=f, af=af, coe=coe)
    
    @staticmethod
    def laplb_csix(bff, f):
        """
        laplb_csix(bff, f)
        
        
        Defined at Isolate_module.fpp lines 1025-1066
        
        Parameters
        ----------
        bff : float array
        f : float array
        
        """
        _AresMainPy.f90wrap_laplb_csix(bff=bff, f=f)
    
    @staticmethod
    def cfour_coe(gaps, coe):
        """
        cfour_coe(gaps, coe)
        
        
        Defined at Isolate_module.fpp lines 1070-1082
        
        Parameters
        ----------
        gaps : float array
        coe : float array
        
        """
        _AresMainPy.f90wrap_cfour_coe(gaps=gaps, coe=coe)
    
    @staticmethod
    def lapla_cfour(u, au, coe):
        """
        lapla_cfour(u, au, coe)
        
        
        Defined at Isolate_module.fpp lines 1086-1106
        
        Parameters
        ----------
        u : float array
        au : float array
        coe : float array
        
        """
        _AresMainPy.f90wrap_lapla_cfour(u=u, au=au, coe=coe)
    
    @staticmethod
    def laplb_cfour(rhs, nrhs):
        """
        laplb_cfour(rhs, nrhs)
        
        
        Defined at Isolate_module.fpp lines 1110-1132
        
        Parameters
        ----------
        rhs : float array
        nrhs : float array
        
        """
        _AresMainPy.f90wrap_laplb_cfour(rhs=rhs, nrhs=nrhs)
    
    @staticmethod
    def calclm(lmax, clm):
        """
        calclm(lmax, clm)
        
        
        Defined at Isolate_module.fpp lines 1136-1145
        
        Parameters
        ----------
        lmax : int
        clm : float array
        
        """
        _AresMainPy.f90wrap_calclm(lmax=lmax, clm=clm)
    
    @staticmethod
    def c(l, m):
        """
        c = c(l, m)
        
        
        Defined at Isolate_module.fpp lines 1149-1162
        
        Parameters
        ----------
        l : int
        m : int
        
        Returns
        -------
        c : float
        
        """
        c = _AresMainPy.f90wrap_c(l=l, m=m)
        return c
    
    @staticmethod
    def rcs(n1xy, n2xy, n3xy, dr, cost, sint, cns, indx):
        """
        rcs(n1xy, n2xy, n3xy, dr, cost, sint, cns, indx)
        
        
        Defined at Isolate_module.fpp lines 1166-1195
        
        Parameters
        ----------
        n1xy : int
        n2xy : int
        n3xy : int
        dr : float array
        cost : float array
        sint : float array
        cns : complex array
        indx : int array
        
        """
        _AresMainPy.f90wrap_rcs(n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, dr=dr, cost=cost, \
            sint=sint, cns=cns, indx=indx)
    
    @staticmethod
    def vhartree_fmm(rhos, vcoulomb):
        """
        vhartree_fmm(rhos, vcoulomb)
        
        
        Defined at Isolate_module.fpp lines 1199-1371
        
        Parameters
        ----------
        rhos : float array
        vcoulomb : float array
        
        ===============================================
        ##CALCULATE VHARTREE VIA FMM CALLED FROM LFMM3D
        ===============================================
        """
        _AresMainPy.f90wrap_vhartree_fmm(rhos=rhos, vcoulomb=vcoulomb)
    
    @staticmethod
    def vhartree_direct(rhos, vcoulomb):
        """
        vhartree_direct(rhos, vcoulomb)
        
        
        Defined at Isolate_module.fpp lines 1374-1459
        
        Parameters
        ----------
        rhos : float array
        vcoulomb : float array
        
        ===============================================
        ##CALCULATE VHARTREE VIA FMM CALLED FROM LFMM3D
        ===============================================
        """
        _AresMainPy.f90wrap_vhartree_direct(rhos=rhos, vcoulomb=vcoulomb)
    
    @staticmethod
    def cal_vsrcpot(size_src, src, pos_src, size_tar, pot_tar, pos_tar):
        """
        cal_vsrcpot(size_src, src, pos_src, size_tar, pot_tar, pos_tar)
        
        
        Defined at Isolate_module.fpp lines 1462-1492
        
        Parameters
        ----------
        size_src : int
        src : float array
        pos_src : float array
        size_tar : int
        pot_tar : float array
        pos_tar : float array
        
        """
        _AresMainPy.f90wrap_cal_vsrcpot(size_src=size_src, src=src, pos_src=pos_src, \
            size_tar=size_tar, pot_tar=pot_tar, pos_tar=pos_tar)
    
    @staticmethod
    def vhart_iso(rho, vhart):
        """
        vhart_iso(rho, vhart)
        
        
        Defined at Isolate_module.fpp lines 1494-1524
        
        Parameters
        ----------
        rho : float array
        vhart : float array
        
        """
        _AresMainPy.f90wrap_vhart_iso(rho=rho, vhart=vhart)
    
    @staticmethod
    def poisson_mcm(rho, pot):
        """
        poisson_mcm(rho, pot)
        
        
        Defined at Isolate_module.fpp lines 1526-1591
        
        Parameters
        ----------
        rho : float array
        pot : float array
        
        """
        _AresMainPy.f90wrap_poisson_mcm(rho=rho, pot=pot)
    
    @staticmethod
    def mcm_calvcorr(q_l0, q_lm, v_corr):
        """
        mcm_calvcorr(q_l0, q_lm, v_corr)
        
        
        Defined at Isolate_module.fpp lines 1593-1740
        
        Parameters
        ----------
        q_l0 : float array
        q_lm : complex array
        v_corr : float array
        
        """
        _AresMainPy.f90wrap_mcm_calvcorr(q_l0=q_l0, q_lm=q_lm, v_corr=v_corr)
    
    @staticmethod
    def mcm_cal_naux(q_l0, q_lm, rhoaux):
        """
        mcm_cal_naux(q_l0, q_lm, rhoaux)
        
        
        Defined at Isolate_module.fpp lines 1742-1801
        
        Parameters
        ----------
        q_l0 : float array
        q_lm : complex array
        rhoaux : float array
        
        """
        _AresMainPy.f90wrap_mcm_cal_naux(q_l0=q_l0, q_lm=q_lm, rhoaux=rhoaux)
    
    @staticmethod
    def cal_ylm(mycost, mysint, e_phi, ylm):
        """
        cal_ylm(mycost, mysint, e_phi, ylm)
        
        
        Defined at Isolate_module.fpp lines 1804-1823
        
        Parameters
        ----------
        mycost : float
        mysint : float
        e_phi : complex
        ylm : complex array
        
        """
        _AresMainPy.f90wrap_cal_ylm(mycost=mycost, mysint=mysint, e_phi=e_phi, ylm=ylm)
    
    @staticmethod
    def rcs_rec(ng1, ng2, ng3, cost, sint, cns):
        """
        rcs_rec(ng1, ng2, ng3, cost, sint, cns)
        
        
        Defined at Isolate_module.fpp lines 1825-1858
        
        Parameters
        ----------
        ng1 : int
        ng2 : int
        ng3 : int
        cost : float array
        sint : float array
        cns : complex array
        
        """
        _AresMainPy.f90wrap_rcs_rec(ng1=ng1, ng2=ng2, ng3=ng3, cost=cost, sint=sint, \
            cns=cns)
    
    @staticmethod
    def car2spe_rec(orig, x, y, z):
        """
        cost, sint, cosp, sinp = car2spe_rec(orig, x, y, z)
        
        
        Defined at Isolate_module.fpp lines 1860-1898
        
        Parameters
        ----------
        orig : float array
        x : float
        y : float
        z : float
        
        Returns
        -------
        cost : float
        sint : float
        cosp : float
        sinp : float
        
        """
        cost, sint, cosp, sinp = _AresMainPy.f90wrap_car2spe_rec(orig=orig, x=x, y=y, \
            z=z)
        return cost, sint, cosp, sinp
    
    @staticmethod
    def integrade_i(l, x):
        """
        integrade_i = integrade_i(l, x)
        
        
        Defined at Isolate_module.fpp lines 1900-2130
        
        Parameters
        ----------
        l : int
        x : float
        
        Returns
        -------
        integrade_i : float
        
        """
        integrade_i = _AresMainPy.f90wrap_integrade_i(l=l, x=x)
        return integrade_i
    
    @staticmethod
    def gmg_cg_hart(rhos, vcoulomb):
        """
        gmg_cg_hart(rhos, vcoulomb)
        
        
        Defined at Isolate_module.fpp lines 2132-2171
        
        Parameters
        ----------
        rhos : float array
        vcoulomb : float array
        
        ===========================================================
        """
        _AresMainPy.f90wrap_gmg_cg_hart(rhos=rhos, vcoulomb=vcoulomb)
    
    @staticmethod
    def cutoff_method(rho, vh):
        """
        cutoff_method(rho, vh)
        
        
        Defined at Isolate_module.fpp lines 2174-2204
        
        Parameters
        ----------
        rho : float array
        vh : float array
        
        """
        _AresMainPy.f90wrap_cutoff_method(rho=rho, vh=vh)
    
    @property
    def center(self):
        """
        Element center ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__center(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            center = self._arrays[array_handle]
        else:
            center = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__center)
            self._arrays[array_handle] = center
        return center
    
    @center.setter
    def center(self, center):
        self.center[...] = center
    
    @property
    def celllengthn(self):
        """
        Element celllengthn ftype=integer(i4b) pytype=int
        
        
        Defined at Isolate_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__celllengthn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            celllengthn = self._arrays[array_handle]
        else:
            celllengthn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__celllengthn)
            self._arrays[array_handle] = celllengthn
        return celllengthn
    
    @celllengthn.setter
    def celllengthn(self, celllengthn):
        self.celllengthn[...] = celllengthn
    
    @property
    def thickness(self):
        """
        Element thickness ftype=integer(i4b) pytype=int
        
        
        Defined at Isolate_module.fpp line 17
        
        """
        return _AresMainPy.f90wrap_isolateset__get__thickness()
    
    @thickness.setter
    def thickness(self, thickness):
        _AresMainPy.f90wrap_isolateset__set__thickness(thickness)
    
    @property
    def cellleft(self):
        """
        Element cellleft ftype=integer(i4b) pytype=int
        
        
        Defined at Isolate_module.fpp line 17
        
        """
        return _AresMainPy.f90wrap_isolateset__get__cellleft()
    
    @cellleft.setter
    def cellleft(self, cellleft):
        _AresMainPy.f90wrap_isolateset__set__cellleft(cellleft)
    
    @property
    def cellright(self):
        """
        Element cellright ftype=integer(i4b) pytype=int
        
        
        Defined at Isolate_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__cellright(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            cellright = self._arrays[array_handle]
        else:
            cellright = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__cellright)
            self._arrays[array_handle] = cellright
        return cellright
    
    @cellright.setter
    def cellright(self, cellright):
        self.cellright[...] = cellright
    
    @property
    def bandright(self):
        """
        Element bandright ftype=integer(i4b) pytype=int
        
        
        Defined at Isolate_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__bandright(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bandright = self._arrays[array_handle]
        else:
            bandright = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__bandright)
            self._arrays[array_handle] = bandright
        return bandright
    
    @bandright.setter
    def bandright(self, bandright):
        self.bandright[...] = bandright
    
    @property
    def boundaryrhos(self):
        """
        Element boundaryrhos ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__boundaryrhos(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            boundaryrhos = self._arrays[array_handle]
        else:
            boundaryrhos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__boundaryrhos)
            self._arrays[array_handle] = boundaryrhos
        return boundaryrhos
    
    @boundaryrhos.setter
    def boundaryrhos(self, boundaryrhos):
        self.boundaryrhos[...] = boundaryrhos
    
    @property
    def boundaryvcoulomb(self):
        """
        Element boundaryvcoulomb ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__boundaryvcoulomb(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            boundaryvcoulomb = self._arrays[array_handle]
        else:
            boundaryvcoulomb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__boundaryvcoulomb)
            self._arrays[array_handle] = boundaryvcoulomb
        return boundaryvcoulomb
    
    @boundaryvcoulomb.setter
    def boundaryvcoulomb(self, boundaryvcoulomb):
        self.boundaryvcoulomb[...] = boundaryvcoulomb
    
    @property
    def ap(self):
        """
        Element ap ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__ap(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            ap = self._arrays[array_handle]
        else:
            ap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__ap)
            self._arrays[array_handle] = ap
        return ap
    
    @ap.setter
    def ap(self, ap):
        self.ap[...] = ap
    
    @property
    def aphi(self):
        """
        Element aphi ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__aphi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            aphi = self._arrays[array_handle]
        else:
            aphi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__aphi)
            self._arrays[array_handle] = aphi
        return aphi
    
    @aphi.setter
    def aphi(self, aphi):
        self.aphi[...] = aphi
    
    @property
    def res(self):
        """
        Element res ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__res(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            res = self._arrays[array_handle]
        else:
            res = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__res)
            self._arrays[array_handle] = res
        return res
    
    @res.setter
    def res(self, res):
        self.res[...] = res
    
    @property
    def res1(self):
        """
        Element res1 ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__res1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            res1 = self._arrays[array_handle]
        else:
            res1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__res1)
            self._arrays[array_handle] = res1
        return res1
    
    @res1.setter
    def res1(self, res1):
        self.res1[...] = res1
    
    @property
    def vcoulomb_old(self):
        """
        Element vcoulomb_old ftype=real(dp) pytype=float
        
        
        Defined at Isolate_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_isolateset__array__vcoulomb_old(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            vcoulomb_old = self._arrays[array_handle]
        else:
            vcoulomb_old = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_isolateset__array__vcoulomb_old)
            self._arrays[array_handle] = vcoulomb_old
        return vcoulomb_old
    
    @vcoulomb_old.setter
    def vcoulomb_old(self, vcoulomb_old):
        self.vcoulomb_old[...] = vcoulomb_old
    
    def __str__(self):
        ret = ['<isolateset>{\n']
        ret.append('    center : ')
        ret.append(repr(self.center))
        ret.append(',\n    celllengthn : ')
        ret.append(repr(self.celllengthn))
        ret.append(',\n    thickness : ')
        ret.append(repr(self.thickness))
        ret.append(',\n    cellleft : ')
        ret.append(repr(self.cellleft))
        ret.append(',\n    cellright : ')
        ret.append(repr(self.cellright))
        ret.append(',\n    bandright : ')
        ret.append(repr(self.bandright))
        ret.append(',\n    boundaryrhos : ')
        ret.append(repr(self.boundaryrhos))
        ret.append(',\n    boundaryvcoulomb : ')
        ret.append(repr(self.boundaryvcoulomb))
        ret.append(',\n    ap : ')
        ret.append(repr(self.ap))
        ret.append(',\n    aphi : ')
        ret.append(repr(self.aphi))
        ret.append(',\n    res : ')
        ret.append(repr(self.res))
        ret.append(',\n    res1 : ')
        ret.append(repr(self.res1))
        ret.append(',\n    vcoulomb_old : ')
        ret.append(repr(self.vcoulomb_old))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

isolateset = Isolateset()

class Potential_Module(f90wrap.runtime.FortranModule):
    """
    Module potential_module
    
    
    Defined at Potential_module.fpp lines 5-436
    
    """
    @staticmethod
    def cal_veff(rhos, veffs):
        """
        cal_veff(rhos, veffs)
        
        
        Defined at Potential_module.fpp lines 19-73
        
        Parameters
        ----------
        rhos : float array
        veffs : float array
        
        """
        _AresMainPy.f90wrap_cal_veff(rhos=rhos, veffs=veffs)
    
    @staticmethod
    def cal_veff_iso(rhos, veffs):
        """
        cal_veff_iso(rhos, veffs)
        
        
        Defined at Potential_module.fpp lines 76-125
        
        Parameters
        ----------
        rhos : float array
        veffs : float array
        
        """
        _AresMainPy.f90wrap_cal_veff_iso(rhos=rhos, veffs=veffs)
    
    @staticmethod
    def vhartree(rho, vhart):
        """
        vhartree(rho, vhart)
        
        
        Defined at Potential_module.fpp lines 128-205
        
        Parameters
        ----------
        rho : float array
        vhart : float array
        
        """
        _AresMainPy.f90wrap_vhartree(rho=rho, vhart=vhart)
    
    @staticmethod
    def vlpp():
        """
        vlpp()
        
        
        Defined at Potential_module.fpp lines 208-276
        
        
        """
        _AresMainPy.f90wrap_vlpp()
    
    @staticmethod
    def vlpp_real():
        """
        vlpp_real()
        
        
        Defined at Potential_module.fpp lines 279-324
        
        
        """
        _AresMainPy.f90wrap_vlpp_real()
    
    @staticmethod
    def vlda(rho, ldapotential):
        """
        vlda(rho, ldapotential)
        
        
        Defined at Potential_module.fpp lines 327-403
        
        Parameters
        ----------
        rho : float array
        ldapotential : float array
        
        ------------------------------------------------------------------------------
         DESCRIPTION:
           This function computes the exchange-correlation potential in the Local
           Density Approximation(LDA) based on the real-space electron density.
         GLOBAL/MODULE VARIABLES CHANGED:
         CONDITIONS AND ASSUMPTIONS:
           We can use LDAPointPot for moudularity.
         FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
           Complete calculation to handle spin-polarized cases.
         REFERENCES:
           [1]  Perdew J.P. and Zunger A., Phys. Rev. B 23(10), 1981, 5048-79
        ------------------------------------------------------------------------------
        """
        _AresMainPy.f90wrap_vlda(rho=rho, ldapotential=ldapotential)
    
    @staticmethod
    def sumrhos(rhos, rho):
        """
        sumrhos(rhos, rho)
        
        
        Defined at Potential_module.fpp lines 406-419
        
        Parameters
        ----------
        rhos : float array
        rho : float array
        
        """
        _AresMainPy.f90wrap_sumrhos(rhos=rhos, rho=rho)
    
    @staticmethod
    def sumrhos_iso(rhos, rho):
        """
        sumrhos_iso(rhos, rho)
        
        
        Defined at Potential_module.fpp lines 422-435
        
        Parameters
        ----------
        rhos : float array
        rho : float array
        
        """
        _AresMainPy.f90wrap_sumrhos_iso(rhos=rhos, rho=rho)
    
    @property
    def v_accelerate(self):
        """
        Element v_accelerate ftype=real(dp) pytype=float
        
        
        Defined at Potential_module.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_potential_module__array__v_accelerate(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            v_accelerate = self._arrays[array_handle]
        else:
            v_accelerate = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_potential_module__array__v_accelerate)
            self._arrays[array_handle] = v_accelerate
        return v_accelerate
    
    @v_accelerate.setter
    def v_accelerate(self, v_accelerate):
        self.v_accelerate[...] = v_accelerate
    
    @property
    def v_hxc_old(self):
        """
        Element v_hxc_old ftype=real(dp) pytype=float
        
        
        Defined at Potential_module.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_potential_module__array__v_hxc_old(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            v_hxc_old = self._arrays[array_handle]
        else:
            v_hxc_old = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_potential_module__array__v_hxc_old)
            self._arrays[array_handle] = v_hxc_old
        return v_hxc_old
    
    @v_hxc_old.setter
    def v_hxc_old(self, v_hxc_old):
        self.v_hxc_old[...] = v_hxc_old
    
    @property
    def v_hxc_new(self):
        """
        Element v_hxc_new ftype=real(dp) pytype=float
        
        
        Defined at Potential_module.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_potential_module__array__v_hxc_new(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            v_hxc_new = self._arrays[array_handle]
        else:
            v_hxc_new = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_potential_module__array__v_hxc_new)
            self._arrays[array_handle] = v_hxc_new
        return v_hxc_new
    
    @v_hxc_new.setter
    def v_hxc_new(self, v_hxc_new):
        self.v_hxc_new[...] = v_hxc_new
    
    @property
    def accelerate(self):
        """
        Element accelerate ftype=logical pytype=bool
        
        
        Defined at Potential_module.fpp line 15
        
        """
        return _AresMainPy.f90wrap_potential_module__get__accelerate()
    
    @accelerate.setter
    def accelerate(self, accelerate):
        _AresMainPy.f90wrap_potential_module__set__accelerate(accelerate)
    
    def __str__(self):
        ret = ['<potential_module>{\n']
        ret.append('    v_accelerate : ')
        ret.append(repr(self.v_accelerate))
        ret.append(',\n    v_hxc_old : ')
        ret.append(repr(self.v_hxc_old))
        ret.append(',\n    v_hxc_new : ')
        ret.append(repr(self.v_hxc_new))
        ret.append(',\n    accelerate : ')
        ret.append(repr(self.accelerate))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

potential_module = Potential_Module()

class Lapack_Module(f90wrap.runtime.FortranModule):
    """
    Module lapack_module
    
    
    Defined at Lapack_module.fpp lines 5-526
    
    """
    @staticmethod
    def diagm(mat, evec, eval):
        """
        diagm(mat, evec, eval)
        
        
        Defined at Lapack_module.fpp lines 20-58
        
        Parameters
        ----------
        mat : complex array
        evec : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_diagm(mat=mat, evec=evec, eval=eval)
    
    @staticmethod
    def generalizeeigen(dime, mata, matb, evec, eval):
        """
        generalizeeigen(dime, mata, matb, evec, eval)
        
        
        Defined at Lapack_module.fpp lines 67-103
        
        Parameters
        ----------
        dime : int
        mata : complex array
        matb : complex array
        evec : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_generalizeeigen(dime=dime, mata=mata, matb=matb, evec=evec, \
            eval=eval)
    
    @staticmethod
    def orthnorm(mat):
        """
        orthnorm(mat)
        
        
        Defined at Lapack_module.fpp lines 108-146
        
        Parameters
        ----------
        mat : complex array
        
        """
        _AresMainPy.f90wrap_orthnorm(mat=mat)
    
    @staticmethod
    def norm_2(mat, k):
        """
        norm_2 = norm_2(mat, k)
        
        
        Defined at Lapack_module.fpp lines 149-163
        
        Parameters
        ----------
        mat : complex array
        k : int
        
        Returns
        -------
        norm_2 : float
        
        """
        norm_2 = _AresMainPy.f90wrap_norm_2(mat=mat, k=k)
        return norm_2
    
    @staticmethod
    def matmat(mata, matb, opa, opb, matc):
        """
        matmat(mata, matb, opa, opb, matc)
        
        
        Defined at Lapack_module.fpp lines 171-201
        
        Parameters
        ----------
        mata : complex array
        matb : complex array
        opa : str
        opb : str
        matc : complex array
        
        """
        _AresMainPy.f90wrap_matmat(mata=mata, matb=matb, opa=opa, opb=opb, matc=matc)
    
    @staticmethod
    def invmat(mat):
        """
        invmat(mat)
        
        
        Defined at Lapack_module.fpp lines 204-236
        
        Parameters
        ----------
        mat : complex array
        
        """
        _AresMainPy.f90wrap_invmat(mat=mat)
    
    @staticmethod
    def orthnorm_real(mat):
        """
        orthnorm_real(mat)
        
        
        Defined at Lapack_module.fpp lines 260-300
        
        Parameters
        ----------
        mat : float array
        
        """
        _AresMainPy.f90wrap_orthnorm_real(mat=mat)
    
    @staticmethod
    def diagm_real(mat, evec, eval):
        """
        diagm_real(mat, evec, eval)
        
        
        Defined at Lapack_module.fpp lines 305-346
        
        Parameters
        ----------
        mat : float array
        evec : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_diagm_real(mat=mat, evec=evec, eval=eval)
    
    @staticmethod
    def generalizeeigen_real(dime, mata, matb, evec, eval):
        """
        generalizeeigen_real(dime, mata, matb, evec, eval)
        
        
        Defined at Lapack_module.fpp lines 354-391
        
        Parameters
        ----------
        dime : int
        mata : float array
        matb : float array
        evec : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_generalizeeigen_real(dime=dime, mata=mata, matb=matb, \
            evec=evec, eval=eval)
    
    @staticmethod
    def diagmx_real(mat, dime, num, il, iu, evec, eval):
        """
        diagmx_real(mat, dime, num, il, iu, evec, eval)
        
        
        Defined at Lapack_module.fpp lines 394-434
        
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
        _AresMainPy.f90wrap_diagmx_real(mat=mat, dime=dime, num=num, il=il, iu=iu, \
            evec=evec, eval=eval)
    
    @staticmethod
    def matmat_real(mata, matb, opa, opb, matc):
        """
        matmat_real(mata, matb, opa, opb, matc)
        
        
        Defined at Lapack_module.fpp lines 437-468
        
        Parameters
        ----------
        mata : float array
        matb : float array
        opa : str
        opb : str
        matc : float array
        
        """
        _AresMainPy.f90wrap_matmat_real(mata=mata, matb=matb, opa=opa, opb=opb, \
            matc=matc)
    
    @staticmethod
    def cholesky_factor_real(mat):
        """
        cholesky_factor_real(mat)
        
        
        Defined at Lapack_module.fpp lines 474-490
        
        Parameters
        ----------
        mat : float array
        
        """
        _AresMainPy.f90wrap_cholesky_factor_real(mat=mat)
    
    @staticmethod
    def invmat_real(mat):
        """
        invmat_real(mat)
        
        
        Defined at Lapack_module.fpp lines 493-525
        
        Parameters
        ----------
        mat : float array
        
        """
        _AresMainPy.f90wrap_invmat_real(mat=mat)
    
    @property
    def mmax(self):
        """
        Element mmax ftype=integer(i4b) pytype=int
        
        
        Defined at Lapack_module.fpp line 14
        
        """
        return _AresMainPy.f90wrap_lapack_module__get__mmax()
    
    @mmax.setter
    def mmax(self, mmax):
        _AresMainPy.f90wrap_lapack_module__set__mmax(mmax)
    
    @property
    def nmax(self):
        """
        Element nmax ftype=integer(i4b) pytype=int
        
        
        Defined at Lapack_module.fpp line 14
        
        """
        return _AresMainPy.f90wrap_lapack_module__get__nmax()
    
    @nmax.setter
    def nmax(self, nmax):
        _AresMainPy.f90wrap_lapack_module__set__nmax(nmax)
    
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

class Matvec_Module(f90wrap.runtime.FortranModule):
    """
    Module matvec_module
    
    
    Defined at Matvec_module.fpp lines 5-945
    
    """
    @staticmethod
    def cmatvec(veff, ik, p, q, dimen):
        """
        cmatvec(veff, ik, p, q, dimen)
        
        
        Defined at Matvec_module.fpp lines 9-163
        
        Parameters
        ----------
        veff : float array
        ik : int
        p : complex array
        q : complex array
        dimen : int
        
        """
        _AresMainPy.f90wrap_cmatvec(veff=veff, ik=ik, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def nlocmatvec(ik, p, q):
        """
        nlocmatvec(ik, p, q)
        
        
        Defined at Matvec_module.fpp lines 166-227
        
        Parameters
        ----------
        ik : int
        p : complex array
        q : complex array
        
        """
        _AresMainPy.f90wrap_nlocmatvec(ik=ik, p=p, q=q)
    
    @staticmethod
    def nlocmatvec_dg(ik, p, q):
        """
        nlocmatvec_dg(ik, p, q)
        
        
        Defined at Matvec_module.fpp lines 230-284
        
        Parameters
        ----------
        ik : int
        p : complex array
        q : complex array
        
        """
        _AresMainPy.f90wrap_nlocmatvec_dg(ik=ik, p=p, q=q)
    
    @staticmethod
    def rmatvec(veff, p, q, dimen):
        """
        rmatvec(veff, p, q, dimen)
        
        
        Defined at Matvec_module.fpp lines 287-347
        
        Parameters
        ----------
        veff : float array
        p : float array
        q : float array
        dimen : int
        
        =====================================================
        ##OUTPUT NONLOCAL TERM TO TEST ITS CORRECTNESS
        open(unit=521,file="nl_Vrho_period")
        write(521,*)n1,n2,n3
        write(521,*)q
        close(521)
        =====================================================
        """
        _AresMainPy.f90wrap_rmatvec(veff=veff, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def nlocmatvec_r(p, q):
        """
        nlocmatvec_r(p, q)
        
        
        Defined at Matvec_module.fpp lines 350-405
        
        Parameters
        ----------
        p : float array
        q : float array
        
        """
        _AresMainPy.f90wrap_nlocmatvec_r(p=p, q=q)
    
    @staticmethod
    def rmatvec_new(mat, p, q, dimen):
        """
        rmatvec_new(mat, p, q, dimen)
        
        
        Defined at Matvec_module.fpp lines 408-416
        
        Parameters
        ----------
        mat : float array
        p : float array
        q : float array
        dimen : int
        
        """
        _AresMainPy.f90wrap_rmatvec_new(mat=mat, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def iso_rmatvec(veff_3d, p, q, dimen):
        """
        iso_rmatvec(veff_3d, p, q, dimen)
        
        
        Defined at Matvec_module.fpp lines 420-457
        
        Parameters
        ----------
        veff_3d : float array
        p : float array
        q : float array
        dimen : int
        
        """
        _AresMainPy.f90wrap_iso_rmatvec(veff_3d=veff_3d, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def iso_grid2sphere(grid_v, grid_v1, sphere_v):
        """
        iso_grid2sphere(grid_v, grid_v1, sphere_v)
        
        
        Defined at Matvec_module.fpp lines 460-473
        
        Parameters
        ----------
        grid_v : float array
        grid_v1 : float array
        sphere_v : float array
        
        """
        _AresMainPy.f90wrap_iso_grid2sphere(grid_v=grid_v, grid_v1=grid_v1, \
            sphere_v=sphere_v)
    
    @staticmethod
    def iso_sphere2grid(sphere_v, grid_v):
        """
        iso_sphere2grid(sphere_v, grid_v)
        
        
        Defined at Matvec_module.fpp lines 476-488
        
        Parameters
        ----------
        sphere_v : float array
        grid_v : float array
        
        """
        _AresMainPy.f90wrap_iso_sphere2grid(sphere_v=sphere_v, grid_v=grid_v)
    
    @staticmethod
    def nlocmatvec_iso(p, q):
        """
        nlocmatvec_iso(p, q)
        
        
        Defined at Matvec_module.fpp lines 491-549
        
        Parameters
        ----------
        p : float array
        q : float array
        
        """
        _AresMainPy.f90wrap_nlocmatvec_iso(p=p, q=q)
    
    @staticmethod
    def nlocmatvec_iso_dg(p, q):
        """
        nlocmatvec_iso_dg(p, q)
        
        
        Defined at Matvec_module.fpp lines 552-609
        
        Parameters
        ----------
        p : float array
        q : float array
        
        """
        _AresMainPy.f90wrap_nlocmatvec_iso_dg(p=p, q=q)
    
    @staticmethod
    def cmatvec_band(veff, ik, p, q, dimen):
        """
        cmatvec_band(veff, ik, p, q, dimen)
        
        
        Defined at Matvec_module.fpp lines 613-699
        
        Parameters
        ----------
        veff : float array
        ik : int
        p : complex array
        q : complex array
        dimen : int
        
        """
        _AresMainPy.f90wrap_cmatvec_band(veff=veff, ik=ik, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def nlocmatvec_band(ik, p, q):
        """
        nlocmatvec_band(ik, p, q)
        
        
        Defined at Matvec_module.fpp lines 702-760
        
        Parameters
        ----------
        ik : int
        p : complex array
        q : complex array
        
        """
        _AresMainPy.f90wrap_nlocmatvec_band(ik=ik, p=p, q=q)
    
    @staticmethod
    def rmatvec_gamma(veff, ik, p, q, dimen):
        """
        rmatvec_gamma(veff, ik, p, q, dimen)
        
        
        Defined at Matvec_module.fpp lines 763-881
        
        Parameters
        ----------
        veff : float array
        ik : int
        p : float array
        q : float array
        dimen : int
        
        """
        _AresMainPy.f90wrap_rmatvec_gamma(veff=veff, ik=ik, p=p, q=q, dimen=dimen)
    
    @staticmethod
    def nlocmatvec_gamma(ik, p, q):
        """
        nlocmatvec_gamma(ik, p, q)
        
        
        Defined at Matvec_module.fpp lines 884-945
        
        Parameters
        ----------
        ik : int
        p : float array
        q : float array
        
        """
        _AresMainPy.f90wrap_nlocmatvec_gamma(ik=ik, p=p, q=q)
    
    _dt_array_initialisers = []
    

matvec_module = Matvec_Module()

class Arpack_Module(f90wrap.runtime.FortranModule):
    """
    Module arpack_module
    
    
    Defined at Arpack_module.fpp lines 5-919
    
    """
    @staticmethod
    def diagh_arpack(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, \
        tol):
        """
        diagh_arpack(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
        
        
        Defined at Arpack_module.fpp lines 23-192
        
        Parameters
        ----------
        veff : float array
        ik : int
        nev : int
        evec : complex array
        eval : float array
        resid_restart : complex array
        nec : int
        info : int
        maxmvs : int
        tol : float
        
        ------
        """
        _AresMainPy.f90wrap_diagh_arpack(veff=veff, ik=ik, nev=nev, evec=evec, \
            eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
            tol=tol)
    
    @staticmethod
    def real_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, \
        tol):
        """
        real_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
        
        
        Defined at Arpack_module.fpp lines 202-363
        
        Parameters
        ----------
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
        _AresMainPy.f90wrap_real_diagh_arpack(veff=veff, nev=nev, evec=evec, eval=eval, \
            resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, tol=tol)
    
    @staticmethod
    def real_diagm_arpk(mat, nev, evec, eval, resid_restart, nec, info, maxmvs, \
        tol):
        """
        real_diagm_arpk(mat, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
        
        
        Defined at Arpack_module.fpp lines 373-506
        
        Parameters
        ----------
        mat : float array
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
        _AresMainPy.f90wrap_real_diagm_arpk(mat=mat, nev=nev, evec=evec, eval=eval, \
            resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, tol=tol)
    
    @staticmethod
    def rdiagm_arpk(mat, dimen, nev, evec, eval):
        """
        rdiagm_arpk(mat, dimen, nev, evec, eval)
        
        
        Defined at Arpack_module.fpp lines 509-544
        
        Parameters
        ----------
        mat : float array
        dimen : int
        nev : int
        evec : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_rdiagm_arpk(mat=mat, dimen=dimen, nev=nev, evec=evec, \
            eval=eval)
    
    @staticmethod
    def iso_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, \
        tol):
        """
        iso_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
        
        
        Defined at Arpack_module.fpp lines 553-716
        
        Parameters
        ----------
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
        _AresMainPy.f90wrap_iso_diagh_arpack(veff=veff, nev=nev, evec=evec, eval=eval, \
            resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, tol=tol)
    
    @staticmethod
    def diagh_arpack_band(veff, ik, nev, evec, eval, resid_restart, nec, info, \
        maxmvs, tol):
        """
        diagh_arpack_band(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, \
            tol)
        
        
        Defined at Arpack_module.fpp lines 719-918
        
        Parameters
        ----------
        veff : float array
        ik : int
        nev : int
        evec : complex array
        eval : float array
        resid_restart : complex array
        nec : int
        info : int
        maxmvs : int
        tol : float
        
        ------
        """
        _AresMainPy.f90wrap_diagh_arpack_band(veff=veff, ik=ik, nev=nev, evec=evec, \
            eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
            tol=tol)
    
    @property
    def maxn(self):
        """
        Element maxn ftype=integer(i4b) pytype=int
        
        
        Defined at Arpack_module.fpp line 14
        
        """
        return _AresMainPy.f90wrap_arpack_module__get__maxn()
    
    @property
    def maxnev(self):
        """
        Element maxnev ftype=integer(i4b) pytype=int
        
        
        Defined at Arpack_module.fpp line 14
        
        """
        return _AresMainPy.f90wrap_arpack_module__get__maxnev()
    
    @property
    def maxncv(self):
        """
        Element maxncv ftype=integer(i4b) pytype=int
        
        
        Defined at Arpack_module.fpp line 14
        
        """
        return _AresMainPy.f90wrap_arpack_module__get__maxncv()
    
    def __str__(self):
        ret = ['<arpack_module>{\n']
        ret.append('    maxn : ')
        ret.append(repr(self.maxn))
        ret.append(',\n    maxnev : ')
        ret.append(repr(self.maxnev))
        ret.append(',\n    maxncv : ')
        ret.append(repr(self.maxncv))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

arpack_module = Arpack_Module()

class Chebyshev_Module(f90wrap.runtime.FortranModule):
    """
    Module chebyshev_module
    
    
    Defined at Chebyshev_fliter.fpp lines 10-2437
    
    """
    @staticmethod
    def cheby_filter_rr(ik, veff, x, d):
        """
        cheby_filter_rr(ik, veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 24-66
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : complex array
        d : float array
        
        """
        _AresMainPy.f90wrap_cheby_filter_rr(ik=ik, veff=veff, x=x, d=d)
    
    @staticmethod
    def grayleigh_ritz(ik, veff, x, d):
        """
        grayleigh_ritz(ik, veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 69-127
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : complex array
        d : float array
        
        """
        _AresMainPy.f90wrap_grayleigh_ritz(ik=ik, veff=veff, x=x, d=d)
    
    @staticmethod
    def chebyshev_filter(ik, veff, x, m, a, b):
        """
        chebyshev_filter(ik, veff, x, m, a, b)
        
        
        Defined at Chebyshev_fliter.fpp lines 130-164
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : complex array
        m : int
        a : float
        b : float
        
        """
        _AresMainPy.f90wrap_chebyshev_filter(ik=ik, veff=veff, x=x, m=m, a=a, b=b)
    
    @staticmethod
    def chebyshev_filter_scaled(ik, veff, x, m, a, b, al):
        """
        chebyshev_filter_scaled(ik, veff, x, m, a, b, al)
        
        
        Defined at Chebyshev_fliter.fpp lines 167-211
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : complex array
        m : int
        a : float
        b : float
        al : float
        
        """
        _AresMainPy.f90wrap_chebyshev_filter_scaled(ik=ik, veff=veff, x=x, m=m, a=a, \
            b=b, al=al)
    
    @staticmethod
    def cal_hx(ik, veff, nst, v, hv):
        """
        cal_hx(ik, veff, nst, v, hv)
        
        
        Defined at Chebyshev_fliter.fpp lines 214-231
        
        Parameters
        ----------
        ik : int
        veff : float array
        nst : int
        v : complex array
        hv : complex array
        
        """
        _AresMainPy.f90wrap_cal_hx(ik=ik, veff=veff, nst=nst, v=v, hv=hv)
    
    @staticmethod
    def rayleigh_quotient(ik, veff, nst, x, xhx):
        """
        rayleigh_quotient(ik, veff, nst, x, xhx)
        
        
        Defined at Chebyshev_fliter.fpp lines 234-281
        
        Parameters
        ----------
        ik : int
        veff : float array
        nst : int
        x : complex array
        xhx : complex array
        
        """
        _AresMainPy.f90wrap_rayleigh_quotient(ik=ik, veff=veff, nst=nst, x=x, xhx=xhx)
    
    @staticmethod
    def rayleigh_ritz(ik, veff, x, d):
        """
        rayleigh_ritz(ik, veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 284-315
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : complex array
        d : float array
        
        """
        _AresMainPy.f90wrap_rayleigh_ritz(ik=ik, veff=veff, x=x, d=d)
    
    @staticmethod
    def estupb(k, ik, veff, vec):
        """
        b = estupb(k, ik, veff, vec)
        
        
        Defined at Chebyshev_fliter.fpp lines 318-377
        
        Parameters
        ----------
        k : int
        ik : int
        veff : float array
        vec : complex array
        
        Returns
        -------
        b : float
        
        """
        b = _AresMainPy.f90wrap_estupb(k=k, ik=ik, veff=veff, vec=vec)
        return b
    
    @staticmethod
    def inituplow(k, ik, veff, v, eval):
        """
        a, b, al = inituplow(k, ik, veff, v, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 384-427
        
        Parameters
        ----------
        k : int
        ik : int
        veff : float array
        v : complex array
        eval : float array
        
        Returns
        -------
        a : float
        b : float
        al : float
        
        """
        a, b, al = _AresMainPy.f90wrap_inituplow(k=k, ik=ik, veff=veff, v=v, eval=eval)
        return a, b, al
    
    @staticmethod
    def first_scfstep_filter(ik, veff, x, eval):
        """
        first_scfstep_filter(ik, veff, x, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 430-484
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : complex array
        eval : float array
        
        -----------------
        """
        _AresMainPy.f90wrap_first_scfstep_filter(ik=ik, veff=veff, x=x, eval=eval)
    
    @staticmethod
    def randomfast(psi, nev, radius):
        """
        randomfast(psi, nev, radius)
        
        
        Defined at Chebyshev_fliter.fpp lines 487-511
        
        Parameters
        ----------
        psi : complex array
        nev : int
        radius : float
        
        """
        _AresMainPy.f90wrap_randomfast(psi=psi, nev=nev, radius=radius)
    
    @staticmethod
    def first_chescf_random(rhos, nev, psi, eval):
        """
        first_chescf_random(rhos, nev, psi, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 514-537
        
        Parameters
        ----------
        rhos : float array
        nev : int
        psi : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_first_chescf_random(rhos=rhos, nev=nev, psi=psi, eval=eval)
    
    @staticmethod
    def first_chescf(veff, psi_ran, nev, psi, eval):
        """
        first_chescf(veff, psi_ran, nev, psi, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 540-559
        
        Parameters
        ----------
        veff : float array
        psi_ran : complex array
        nev : int
        psi : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_first_chescf(veff=veff, psi_ran=psi_ran, nev=nev, psi=psi, \
            eval=eval)
    
    @staticmethod
    def cheby_init_sto(nmax, npw, initx_sto):
        """
        cheby_init_sto(nmax, npw, initx_sto)
        
        
        Defined at Chebyshev_fliter.fpp lines 565-655
        
        Parameters
        ----------
        nmax : int
        npw : int
        initx_sto : complex array
        
        """
        _AresMainPy.f90wrap_cheby_init_sto(nmax=nmax, npw=npw, initx_sto=initx_sto)
    
    @staticmethod
    def first_chescf_sto(veff, nev, psi, eval):
        """
        first_chescf_sto(veff, nev, psi, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 658-720
        
        Parameters
        ----------
        veff : float array
        nev : int
        psi : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_first_chescf_sto(veff=veff, nev=nev, psi=psi, eval=eval)
    
    @staticmethod
    def first_subspace_stopw(ik, veff, nmax, nev, initx, x, d):
        """
        first_subspace_stopw(ik, veff, nmax, nev, initx, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 723-823
        
        Parameters
        ----------
        ik : int
        veff : float array
        nmax : int
        nev : int
        initx : complex array
        x : complex array
        d : float array
        
        """
        _AresMainPy.f90wrap_first_subspace_stopw(ik=ik, veff=veff, nmax=nmax, nev=nev, \
            initx=initx, x=x, d=d)
    
    @staticmethod
    def bvk_first_chescf_sto_rand(rhos, nev, psi, eval):
        """
        bvk_first_chescf_sto_rand(rhos, nev, psi, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 832-875
        
        Parameters
        ----------
        rhos : float array
        nev : int
        psi : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_bvk_first_chescf_sto_rand(rhos=rhos, nev=nev, psi=psi, \
            eval=eval)
    
    @staticmethod
    def bvk_cheby_init_sto_rand(nmax, nrand, initx_sto):
        """
        bvk_cheby_init_sto_rand(nmax, nrand, initx_sto)
        
        
        Defined at Chebyshev_fliter.fpp lines 878-958
        
        Parameters
        ----------
        nmax : int
        nrand : int
        initx_sto : float array
        
        """
        _AresMainPy.f90wrap_bvk_cheby_init_sto_rand(nmax=nmax, nrand=nrand, \
            initx_sto=initx_sto)
    
    @staticmethod
    def first_subspace_sto_rand(veff, nmax, nev, initx, x, d):
        """
        first_subspace_sto_rand(veff, nmax, nev, initx, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 961-988
        
        Parameters
        ----------
        veff : float array
        nmax : int
        nev : int
        initx : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_first_subspace_sto_rand(veff=veff, nmax=nmax, nev=nev, \
            initx=initx, x=x, d=d)
    
    @staticmethod
    def rayleigh_quotient_real(veff, nst, x, xhx):
        """
        rayleigh_quotient_real(veff, nst, x, xhx)
        
        
        Defined at Chebyshev_fliter.fpp lines 996-1021
        
        Parameters
        ----------
        veff : float array
        nst : int
        x : float array
        xhx : float array
        
        """
        _AresMainPy.f90wrap_rayleigh_quotient_real(veff=veff, nst=nst, x=x, xhx=xhx)
    
    @staticmethod
    def cal_hx_real(veff, nst, v, hv):
        """
        cal_hx_real(veff, nst, v, hv)
        
        
        Defined at Chebyshev_fliter.fpp lines 1024-1043
        
        Parameters
        ----------
        veff : float array
        nst : int
        v : float array
        hv : float array
        
        =========================================
        print *,"nst",nst
        print *,'cal HX time -->',(ts-tt)/10000.d0
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        """
        _AresMainPy.f90wrap_cal_hx_real(veff=veff, nst=nst, v=v, hv=hv)
    
    @staticmethod
    def estupb_real(k, veff, vec):
        """
        b = estupb_real(k, veff, vec)
        
        
        Defined at Chebyshev_fliter.fpp lines 1046-1089
        
        Parameters
        ----------
        k : int
        veff : float array
        vec : float array
        
        Returns
        -------
        b : float
        
        """
        b = _AresMainPy.f90wrap_estupb_real(k=k, veff=veff, vec=vec)
        return b
    
    @staticmethod
    def chebyshev_filter_scaled_real(veff, x, m, a, b, al):
        """
        chebyshev_filter_scaled_real(veff, x, m, a, b, al)
        
        
        Defined at Chebyshev_fliter.fpp lines 1092-1135
        
        Parameters
        ----------
        veff : float array
        x : float array
        m : int
        a : float
        b : float
        al : float
        
        ===========================
        print*,"{-------------------"
        ===========================
        """
        _AresMainPy.f90wrap_chebyshev_filter_scaled_real(veff=veff, x=x, m=m, a=a, b=b, \
            al=al)
    
    @staticmethod
    def grayleigh_ritz_real(veff, x, d):
        """
        grayleigh_ritz_real(veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 1138-1181
        
        Parameters
        ----------
        veff : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_grayleigh_ritz_real(veff=veff, x=x, d=d)
    
    @staticmethod
    def cheby_filtering_grrr(veff, x, d):
        """
        cheby_filtering_grrr(veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 1184-1202
        
        Parameters
        ----------
        veff : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_cheby_filtering_grrr(veff=veff, x=x, d=d)
    
    @staticmethod
    def cheby_filtering_prrr(veff, x, nfs, cfr, efr):
        """
        cheby_filtering_prrr(veff, x, nfs, cfr, efr)
        
        
        Defined at Chebyshev_fliter.fpp lines 1210-1230
        
        Parameters
        ----------
        veff : float array
        x : float array
        nfs : int
        cfr : float array
        efr : float array
        
        """
        _AresMainPy.f90wrap_cheby_filtering_prrr(veff=veff, x=x, nfs=nfs, cfr=cfr, \
            efr=efr)
    
    @staticmethod
    def partialrayleighritz(veff, x, nfs, cfr, efr):
        """
        partialrayleighritz(veff, x, nfs, cfr, efr)
        
        
        Defined at Chebyshev_fliter.fpp lines 1233-1250
        
        Parameters
        ----------
        veff : float array
        x : float array
        nfs : int
        cfr : float array
        efr : float array
        
        """
        _AresMainPy.f90wrap_partialrayleighritz(veff=veff, x=x, nfs=nfs, cfr=cfr, \
            efr=efr)
    
    @staticmethod
    def prr_orthnorm(veff, x, nfs, cfr, efr):
        """
        prr_orthnorm(veff, x, nfs, cfr, efr)
        
        
        Defined at Chebyshev_fliter.fpp lines 1253-1271
        
        Parameters
        ----------
        veff : float array
        x : float array
        nfs : int
        cfr : float array
        efr : float array
        
        """
        _AresMainPy.f90wrap_prr_orthnorm(veff=veff, x=x, nfs=nfs, cfr=cfr, efr=efr)
    
    @staticmethod
    def iso_first_chescf_sto_rand(rhos, nev, psi, eval):
        """
        iso_first_chescf_sto_rand(rhos, nev, psi, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 1279-1375
        
        Parameters
        ----------
        rhos : float array
        nev : int
        psi : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_iso_first_chescf_sto_rand(rhos=rhos, nev=nev, psi=psi, \
            eval=eval)
    
    @staticmethod
    def iso_cheby_init_sto_rand(nmax, nrand, initx_sto):
        """
        iso_cheby_init_sto_rand(nmax, nrand, initx_sto)
        
        
        Defined at Chebyshev_fliter.fpp lines 1379-1457
        
        Parameters
        ----------
        nmax : int
        nrand : int
        initx_sto : float array
        
        """
        _AresMainPy.f90wrap_iso_cheby_init_sto_rand(nmax=nmax, nrand=nrand, \
            initx_sto=initx_sto)
    
    @staticmethod
    def iso_cheby_init_rand(nmax, nrand, initx_sto):
        """
        iso_cheby_init_rand(nmax, nrand, initx_sto)
        
        
        Defined at Chebyshev_fliter.fpp lines 1460-1498
        
        Parameters
        ----------
        nmax : int
        nrand : int
        initx_sto : float array
        
        """
        _AresMainPy.f90wrap_iso_cheby_init_rand(nmax=nmax, nrand=nrand, \
            initx_sto=initx_sto)
    
    @staticmethod
    def iso_first_subspace_sto_rand(veff_3d, nmax, nev, initx, x, d):
        """
        iso_first_subspace_sto_rand(veff_3d, nmax, nev, initx, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 1501-1628
        
        Parameters
        ----------
        veff_3d : float array
        nmax : int
        nev : int
        initx : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_iso_first_subspace_sto_rand(veff_3d=veff_3d, nmax=nmax, \
            nev=nev, initx=initx, x=x, d=d)
    
    @staticmethod
    def rayleigh_quotient_iso(veff_3d, nst, x, xhx):
        """
        rayleigh_quotient_iso(veff_3d, nst, x, xhx)
        
        
        Defined at Chebyshev_fliter.fpp lines 1632-1657
        
        Parameters
        ----------
        veff_3d : float array
        nst : int
        x : float array
        xhx : float array
        
        """
        _AresMainPy.f90wrap_rayleigh_quotient_iso(veff_3d=veff_3d, nst=nst, x=x, \
            xhx=xhx)
    
    @staticmethod
    def cal_hx_iso(veff, nst, v, hv):
        """
        cal_hx_iso(veff, nst, v, hv)
        
        
        Defined at Chebyshev_fliter.fpp lines 1661-1689
        
        Parameters
        ----------
        veff : float array
        nst : int
        v : float array
        hv : float array
        
        =========================================
        """
        _AresMainPy.f90wrap_cal_hx_iso(veff=veff, nst=nst, v=v, hv=hv)
    
    @staticmethod
    def cheby_filtering_grriso(veff_3d, x_sphe, d):
        """
        cheby_filtering_grriso(veff_3d, x_sphe, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 1694-1715
        
        Parameters
        ----------
        veff_3d : float array
        x_sphe : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_cheby_filtering_grriso(veff_3d=veff_3d, x_sphe=x_sphe, d=d)
    
    @staticmethod
    def estupb_iso(k, veff_3d, vec):
        """
        b = estupb_iso(k, veff_3d, vec)
        
        
        Defined at Chebyshev_fliter.fpp lines 1720-1783
        
        Parameters
        ----------
        k : int
        veff_3d : float array
        vec : float array
        
        Returns
        -------
        b : float
        
        """
        b = _AresMainPy.f90wrap_estupb_iso(k=k, veff_3d=veff_3d, vec=vec)
        return b
    
    @staticmethod
    def chebyshev_filter_scaled_iso(veff_3d, x, m, a, b, al):
        """
        chebyshev_filter_scaled_iso(veff_3d, x, m, a, b, al)
        
        
        Defined at Chebyshev_fliter.fpp lines 1788-1841
        
        Parameters
        ----------
        veff_3d : float array
        x : float array
        m : int
        a : float
        b : float
        al : float
        
        """
        _AresMainPy.f90wrap_chebyshev_filter_scaled_iso(veff_3d=veff_3d, x=x, m=m, a=a, \
            b=b, al=al)
    
    @staticmethod
    def grayleigh_ritz_iso(veff_3d, x, d):
        """
        grayleigh_ritz_iso(veff_3d, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 1845-1925
        
        Parameters
        ----------
        veff_3d : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_grayleigh_ritz_iso(veff_3d=veff_3d, x=x, d=d)
    
    @staticmethod
    def first_chescf_sto_gamma(veff, nev, psi, eval):
        """
        first_chescf_sto_gamma(veff, nev, psi, eval)
        
        
        Defined at Chebyshev_fliter.fpp lines 1929-1990
        
        Parameters
        ----------
        veff : float array
        nev : int
        psi : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_first_chescf_sto_gamma(veff=veff, nev=nev, psi=psi, \
            eval=eval)
    
    @staticmethod
    def cheby_init_sto_gamma(nmax, npw, initx_sto):
        """
        cheby_init_sto_gamma(nmax, npw, initx_sto)
        
        
        Defined at Chebyshev_fliter.fpp lines 1993-2083
        
        Parameters
        ----------
        nmax : int
        npw : int
        initx_sto : float array
        
        """
        _AresMainPy.f90wrap_cheby_init_sto_gamma(nmax=nmax, npw=npw, \
            initx_sto=initx_sto)
    
    @staticmethod
    def first_subspace_stopw_gamma(ik, veff, nmax, nev, initx, x, d):
        """
        first_subspace_stopw_gamma(ik, veff, nmax, nev, initx, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 2086-2185
        
        Parameters
        ----------
        ik : int
        veff : float array
        nmax : int
        nev : int
        initx : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_first_subspace_stopw_gamma(ik=ik, veff=veff, nmax=nmax, \
            nev=nev, initx=initx, x=x, d=d)
    
    @staticmethod
    def rayleigh_quotient_gamma(ik, veff, nst, x, xhx):
        """
        rayleigh_quotient_gamma(ik, veff, nst, x, xhx)
        
        
        Defined at Chebyshev_fliter.fpp lines 2188-2213
        
        Parameters
        ----------
        ik : int
        veff : float array
        nst : int
        x : float array
        xhx : float array
        
        """
        _AresMainPy.f90wrap_rayleigh_quotient_gamma(ik=ik, veff=veff, nst=nst, x=x, \
            xhx=xhx)
    
    @staticmethod
    def cal_hx_gamma(ik, veff, nst, v, hv):
        """
        cal_hx_gamma(ik, veff, nst, v, hv)
        
        
        Defined at Chebyshev_fliter.fpp lines 2216-2233
        
        Parameters
        ----------
        ik : int
        veff : float array
        nst : int
        v : float array
        hv : float array
        
        """
        _AresMainPy.f90wrap_cal_hx_gamma(ik=ik, veff=veff, nst=nst, v=v, hv=hv)
    
    @staticmethod
    def cheby_filter_rr_gamma(ik, veff, x, d):
        """
        cheby_filter_rr_gamma(ik, veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 2235-2269
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_cheby_filter_rr_gamma(ik=ik, veff=veff, x=x, d=d)
    
    @staticmethod
    def estupb_gamma(k, ik, veff, vec):
        """
        b = estupb_gamma(k, ik, veff, vec)
        
        
        Defined at Chebyshev_fliter.fpp lines 2271-2330
        
        Parameters
        ----------
        k : int
        ik : int
        veff : float array
        vec : float array
        
        Returns
        -------
        b : float
        
        """
        b = _AresMainPy.f90wrap_estupb_gamma(k=k, ik=ik, veff=veff, vec=vec)
        return b
    
    @staticmethod
    def chebyshev_filter_scaled_gamma(ik, veff, x, m, a, b, al):
        """
        chebyshev_filter_scaled_gamma(ik, veff, x, m, a, b, al)
        
        
        Defined at Chebyshev_fliter.fpp lines 2332-2376
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : float array
        m : int
        a : float
        b : float
        al : float
        
        """
        _AresMainPy.f90wrap_chebyshev_filter_scaled_gamma(ik=ik, veff=veff, x=x, m=m, \
            a=a, b=b, al=al)
    
    @staticmethod
    def grayleigh_ritz_gamma(ik, veff, x, d):
        """
        grayleigh_ritz_gamma(ik, veff, x, d)
        
        
        Defined at Chebyshev_fliter.fpp lines 2378-2436
        
        Parameters
        ----------
        ik : int
        veff : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_grayleigh_ritz_gamma(ik=ik, veff=veff, x=x, d=d)
    
    @property
    def iiii(self):
        """
        Element iiii ftype=integer(i4b) pytype=int
        
        
        Defined at Chebyshev_fliter.fpp line 16
        
        """
        return _AresMainPy.f90wrap_chebyshev_module__get__iiii()
    
    @iiii.setter
    def iiii(self, iiii):
        _AresMainPy.f90wrap_chebyshev_module__set__iiii(iiii)
    
    def __str__(self):
        ret = ['<chebyshev_module>{\n']
        ret.append('    iiii : ')
        ret.append(repr(self.iiii))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

chebyshev_module = Chebyshev_Module()

class Energy_Module(f90wrap.runtime.FortranModule):
    """
    Module energy_module
    
    
    Defined at Energy_module.fpp lines 5-763
    
    """
    @staticmethod
    def totalenergy(psi, rhos, eval, fmenergy, ets, llast):
        """
        totalenergy(psi, rhos, eval, fmenergy, ets, llast)
        
        
        Defined at Energy_module.fpp lines 32-127
        
        Parameters
        ----------
        psi : complex array
        rhos : float array
        eval : float array
        fmenergy : float
        ets : float
        llast : bool
        
        """
        _AresMainPy.f90wrap_totalenergy(psi=psi, rhos=rhos, eval=eval, \
            fmenergy=fmenergy, ets=ets, llast=llast)
    
    @staticmethod
    def out_ke_potential(efm, vie, vh, vxc):
        """
        out_ke_potential(efm, vie, vh, vxc)
        
        
        Defined at Energy_module.fpp lines 129-156
        
        Parameters
        ----------
        efm : float
        vie : float array
        vh : float array
        vxc : float array
        
        """
        _AresMainPy.f90wrap_out_ke_potential(efm=efm, vie=vie, vh=vh, vxc=vxc)
    
    @staticmethod
    def totalenergy_gamma(psi, rhos, eval, fmenergy, ets, llast):
        """
        totalenergy_gamma(psi, rhos, eval, fmenergy, ets, llast)
        
        
        Defined at Energy_module.fpp lines 159-250
        
        Parameters
        ----------
        psi : float array
        rhos : float array
        eval : float array
        fmenergy : float
        ets : float
        llast : bool
        
        """
        _AresMainPy.f90wrap_totalenergy_gamma(psi=psi, rhos=rhos, eval=eval, \
            fmenergy=fmenergy, ets=ets, llast=llast)
    
    @staticmethod
    def ehartree(rho, vhart):
        """
        ehart = ehartree(rho, vhart)
        
        
        Defined at Energy_module.fpp lines 253-270
        
        Parameters
        ----------
        rho : float array
        vhart : float array
        
        Returns
        -------
        ehart : float
        
        """
        ehart = _AresMainPy.f90wrap_ehartree(rho=rho, vhart=vhart)
        return ehart
    
    @staticmethod
    def ehartree_iso(rho, vhart):
        """
        ehart = ehartree_iso(rho, vhart)
        
        
        Defined at Energy_module.fpp lines 273-284
        
        Parameters
        ----------
        rho : float array
        vhart : float array
        
        Returns
        -------
        ehart : float
        
        """
        ehart = _AresMainPy.f90wrap_ehartree_iso(rho=rho, vhart=vhart)
        return ehart
    
    @staticmethod
    def e_ie(rho, vion):
        """
        eie = e_ie(rho, vion)
        
        
        Defined at Energy_module.fpp lines 287-304
        
        Parameters
        ----------
        rho : float array
        vion : float array
        
        Returns
        -------
        eie : float
        
        """
        eie = _AresMainPy.f90wrap_e_ie(rho=rho, vion=vion)
        return eie
    
    @staticmethod
    def ebands(eval, wke):
        """
        eband = ebands(eval, wke)
        
        
        Defined at Energy_module.fpp lines 307-338
        
        Parameters
        ----------
        eval : float array
        wke : float array
        
        Returns
        -------
        eband : float
        
        """
        eband = _AresMainPy.f90wrap_ebands(eval=eval, wke=wke)
        return eband
    
    @staticmethod
    def bvk_totalenergy(rhos, eval, fmenergy, ets, llast):
        """
        bvk_totalenergy(rhos, eval, fmenergy, ets, llast)
        
        
        Defined at Energy_module.fpp lines 344-394
        
        Parameters
        ----------
        rhos : float array
        eval : float array
        fmenergy : float
        ets : float
        llast : bool
        
        """
        _AresMainPy.f90wrap_bvk_totalenergy(rhos=rhos, eval=eval, fmenergy=fmenergy, \
            ets=ets, llast=llast)
    
    @staticmethod
    def bvk_ebands(eval, focc):
        """
        eband = bvk_ebands(eval, focc)
        
        
        Defined at Energy_module.fpp lines 397-420
        
        Parameters
        ----------
        eval : float array
        focc : float array
        
        Returns
        -------
        eband : float
        
        """
        eband = _AresMainPy.f90wrap_bvk_ebands(eval=eval, focc=focc)
        return eband
    
    @staticmethod
    def prr_ebands(veff, phi, pbar):
        """
        eband = prr_ebands(veff, phi, pbar)
        
        
        Defined at Energy_module.fpp lines 426-456
        
        Parameters
        ----------
        veff : float array
        phi : float array
        pbar : float array
        
        Returns
        -------
        eband : float
        
        """
        eband = _AresMainPy.f90wrap_prr_ebands(veff=veff, phi=phi, pbar=pbar)
        return eband
    
    @staticmethod
    def prr_totalenergy(rhos, phi, pbar, fmenergy, ets, llast):
        """
        prr_totalenergy(rhos, phi, pbar, fmenergy, ets, llast)
        
        
        Defined at Energy_module.fpp lines 459-516
        
        Parameters
        ----------
        rhos : float array
        phi : float array
        pbar : float array
        fmenergy : float
        ets : float
        llast : bool
        
        """
        _AresMainPy.f90wrap_prr_totalenergy(rhos=rhos, phi=phi, pbar=pbar, \
            fmenergy=fmenergy, ets=ets, llast=llast)
    
    @staticmethod
    def iso_totalenergy(rhos, eval, fmenergy, ets, llast):
        """
        iso_totalenergy(rhos, eval, fmenergy, ets, llast)
        
        
        Defined at Energy_module.fpp lines 522-646
        
        Parameters
        ----------
        rhos : float array
        eval : float array
        fmenergy : float
        ets : float
        llast : bool
        
        ------------------total energy calculation------------------
        """
        _AresMainPy.f90wrap_iso_totalenergy(rhos=rhos, eval=eval, fmenergy=fmenergy, \
            ets=ets, llast=llast)
    
    @staticmethod
    def iso_ebands(eval, focc):
        """
        eband = iso_ebands(eval, focc)
        
        
        Defined at Energy_module.fpp lines 650-673
        
        Parameters
        ----------
        eval : float array
        focc : float array
        
        Returns
        -------
        eband : float
        
        """
        eband = _AresMainPy.f90wrap_iso_ebands(eval=eval, focc=focc)
        return eband
    
    @staticmethod
    def iso_e_ie(rho, vionlpp):
        """
        eie = iso_e_ie(rho, vionlpp)
        
        
        Defined at Energy_module.fpp lines 676-708
        
        Parameters
        ----------
        rho : float array
        vionlpp : float array
        
        Returns
        -------
        eie : float
        
        """
        eie = _AresMainPy.f90wrap_iso_e_ie(rho=rho, vionlpp=vionlpp)
        return eie
    
    @staticmethod
    def lda_energy(rhoreal):
        """
        lda_energy = lda_energy(rhoreal)
        
        
        Defined at Energy_module.fpp lines 712-762
        
        Parameters
        ----------
        rhoreal : float array
        
        Returns
        -------
        lda_energy : float
        
        """
        lda_energy = _AresMainPy.f90wrap_lda_energy(rhoreal=rhoreal)
        return lda_energy
    
    @property
    def etot(self):
        """
        Element etot ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__etot()
    
    @etot.setter
    def etot(self, etot):
        _AresMainPy.f90wrap_energy_module__set__etot(etot)
    
    @property
    def ekine(self):
        """
        Element ekine ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__ekine()
    
    @ekine.setter
    def ekine(self, ekine):
        _AresMainPy.f90wrap_energy_module__set__ekine(ekine)
    
    @property
    def ehart(self):
        """
        Element ehart ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__ehart()
    
    @ehart.setter
    def ehart(self, ehart):
        _AresMainPy.f90wrap_energy_module__set__ehart(ehart)
    
    @property
    def exc(self):
        """
        Element exc ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__exc()
    
    @exc.setter
    def exc(self, exc):
        _AresMainPy.f90wrap_energy_module__set__exc(exc)
    
    @property
    def eband(self):
        """
        Element eband ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__eband()
    
    @eband.setter
    def eband(self, eband):
        _AresMainPy.f90wrap_energy_module__set__eband(eband)
    
    @property
    def eie(self):
        """
        Element eie ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__eie()
    
    @eie.setter
    def eie(self, eie):
        _AresMainPy.f90wrap_energy_module__set__eie(eie)
    
    @property
    def eienl(self):
        """
        Element eienl ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__eienl()
    
    @eienl.setter
    def eienl(self, eienl):
        _AresMainPy.f90wrap_energy_module__set__eienl(eienl)
    
    @property
    def eewald(self):
        """
        Element eewald ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__eewald()
    
    @eewald.setter
    def eewald(self, eewald):
        _AresMainPy.f90wrap_energy_module__set__eewald(eewald)
    
    @property
    def efm(self):
        """
        Element efm ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__efm()
    
    @efm.setter
    def efm(self, efm):
        _AresMainPy.f90wrap_energy_module__set__efm(efm)
    
    @property
    def tnad(self):
        """
        Element tnad ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__tnad()
    
    @tnad.setter
    def tnad(self, tnad):
        _AresMainPy.f90wrap_energy_module__set__tnad(tnad)
    
    @property
    def fe(self):
        """
        Element fe ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__fe()
    
    @fe.setter
    def fe(self, fe):
        _AresMainPy.f90wrap_energy_module__set__fe(fe)
    
    @property
    def fe0(self):
        """
        Element fe0 ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 23
        
        """
        return _AresMainPy.f90wrap_energy_module__get__fe0()
    
    @fe0.setter
    def fe0(self, fe0):
        _AresMainPy.f90wrap_energy_module__set__fe0(fe0)
    
    @property
    def wmaxl(self):
        """
        Element wmaxl ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 25
        
        """
        return _AresMainPy.f90wrap_energy_module__get__wmaxl()
    
    @wmaxl.setter
    def wmaxl(self, wmaxl):
        _AresMainPy.f90wrap_energy_module__set__wmaxl(wmaxl)
    
    @property
    def wminl(self):
        """
        Element wminl ftype=real(dp) pytype=float
        
        
        Defined at Energy_module.fpp line 25
        
        """
        return _AresMainPy.f90wrap_energy_module__get__wminl()
    
    @wminl.setter
    def wminl(self, wminl):
        _AresMainPy.f90wrap_energy_module__set__wminl(wminl)
    
    def __str__(self):
        ret = ['<energy_module>{\n']
        ret.append('    etot : ')
        ret.append(repr(self.etot))
        ret.append(',\n    ekine : ')
        ret.append(repr(self.ekine))
        ret.append(',\n    ehart : ')
        ret.append(repr(self.ehart))
        ret.append(',\n    exc : ')
        ret.append(repr(self.exc))
        ret.append(',\n    eband : ')
        ret.append(repr(self.eband))
        ret.append(',\n    eie : ')
        ret.append(repr(self.eie))
        ret.append(',\n    eienl : ')
        ret.append(repr(self.eienl))
        ret.append(',\n    eewald : ')
        ret.append(repr(self.eewald))
        ret.append(',\n    efm : ')
        ret.append(repr(self.efm))
        ret.append(',\n    tnad : ')
        ret.append(repr(self.tnad))
        ret.append(',\n    fe : ')
        ret.append(repr(self.fe))
        ret.append(',\n    fe0 : ')
        ret.append(repr(self.fe0))
        ret.append(',\n    wmaxl : ')
        ret.append(repr(self.wmaxl))
        ret.append(',\n    wminl : ')
        ret.append(repr(self.wminl))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

energy_module = Energy_Module()

class Mixer_Module(f90wrap.runtime.FortranModule):
    """
    Module mixer_module
    
    
    Defined at Mixer_module.fpp lines 5-1148
    
    """
    @f90wrap.runtime.register_class("AresMainPy.mixer_data")
    class mixer_data(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=mixer_data)
        
        
        Defined at Mixer_module.fpp lines 16-37
        
        """
        def __init__(self, handle=None):
            """
            self = Mixer_Data()
            
            
            Defined at Mixer_module.fpp lines 16-37
            
            
            Returns
            -------
            this : Mixer_Data
            	Object to be constructed
            
            
            Automatically generated constructor for mixer_data
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_mixer_data_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Mixer_Data
            
            
            Defined at Mixer_module.fpp lines 16-37
            
            Parameters
            ----------
            this : Mixer_Data
            	Object to be destructed
            
            
            Automatically generated destructor for mixer_data
            """
            if self._alloc:
                _AresMainPy.f90wrap_mixer_data_finalise(this=self._handle)
        
        @property
        def nhmix(self):
            """
            Element nhmix ftype=integer  pytype=int
            
            
            Defined at Mixer_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__nhmix(self._handle)
        
        @nhmix.setter
        def nhmix(self, nhmix):
            _AresMainPy.f90wrap_mixer_data__set__nhmix(self._handle, nhmix)
        
        @property
        def nhmin(self):
            """
            Element nhmin ftype=integer  pytype=int
            
            
            Defined at Mixer_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__nhmin(self._handle)
        
        @nhmin.setter
        def nhmin(self, nhmin):
            _AresMainPy.f90wrap_mixer_data__set__nhmin(self._handle, nhmin)
        
        @property
        def nam(self):
            """
            Element nam ftype=integer  pytype=int
            
            
            Defined at Mixer_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__nam(self._handle)
        
        @nam.setter
        def nam(self, nam):
            _AresMainPy.f90wrap_mixer_data__set__nam(self._handle, nam)
        
        @property
        def nitra(self):
            """
            Element nitra ftype=integer  pytype=int
            
            
            Defined at Mixer_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__nitra(self._handle)
        
        @nitra.setter
        def nitra(self, nitra):
            _AresMainPy.f90wrap_mixer_data__set__nitra(self._handle, nitra)
        
        @property
        def alpha(self):
            """
            Element alpha ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__alpha(self._handle)
        
        @alpha.setter
        def alpha(self, alpha):
            _AresMainPy.f90wrap_mixer_data__set__alpha(self._handle, alpha)
        
        @property
        def beta(self):
            """
            Element beta ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__beta(self._handle)
        
        @beta.setter
        def beta(self, beta):
            _AresMainPy.f90wrap_mixer_data__set__beta(self._handle, beta)
        
        @property
        def w0(self):
            """
            Element w0 ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 25
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__w0(self._handle)
        
        @w0.setter
        def w0(self, w0):
            _AresMainPy.f90wrap_mixer_data__set__w0(self._handle, w0)
        
        @property
        def dxl(self):
            """
            Element dxl ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_mixer_data__array__dxl(self._handle)
            if array_handle in self._arrays:
                dxl = self._arrays[array_handle]
            else:
                dxl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_mixer_data__array__dxl)
                self._arrays[array_handle] = dxl
            return dxl
        
        @dxl.setter
        def dxl(self, dxl):
            self.dxl[...] = dxl
        
        @property
        def dfl(self):
            """
            Element dfl ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_mixer_data__array__dfl(self._handle)
            if array_handle in self._arrays:
                dfl = self._arrays[array_handle]
            else:
                dfl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_mixer_data__array__dfl)
                self._arrays[array_handle] = dfl
            return dfl
        
        @dfl.setter
        def dfl(self, dfl):
            self.dfl[...] = dfl
        
        @property
        def voma(self):
            """
            Element voma ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_mixer_data__array__voma(self._handle)
            if array_handle in self._arrays:
                voma = self._arrays[array_handle]
            else:
                voma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_mixer_data__array__voma)
                self._arrays[array_handle] = voma
            return voma
        
        @voma.setter
        def voma(self, voma):
            self.voma[...] = voma
        
        @property
        def sp(self):
            """
            Element sp ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 30
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__sp(self._handle)
        
        @sp.setter
        def sp(self, sp):
            _AresMainPy.f90wrap_mixer_data__set__sp(self._handle, sp)
        
        @property
        def kerker(self):
            """
            Element kerker ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 32
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_mixer_data__array__kerker(self._handle)
            if array_handle in self._arrays:
                kerker = self._arrays[array_handle]
            else:
                kerker = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_mixer_data__array__kerker)
                self._arrays[array_handle] = kerker
            return kerker
        
        @kerker.setter
        def kerker(self, kerker):
            self.kerker[...] = kerker
        
        @property
        def dxgl(self):
            """
            Element dxgl ftype=complex(dp) pytype=complex
            
            
            Defined at Mixer_module.fpp line 35
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_mixer_data__array__dxgl(self._handle)
            if array_handle in self._arrays:
                dxgl = self._arrays[array_handle]
            else:
                dxgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_mixer_data__array__dxgl)
                self._arrays[array_handle] = dxgl
            return dxgl
        
        @dxgl.setter
        def dxgl(self, dxgl):
            self.dxgl[...] = dxgl
        
        @property
        def dfgl(self):
            """
            Element dfgl ftype=complex(dp) pytype=complex
            
            
            Defined at Mixer_module.fpp line 35
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_mixer_data__array__dfgl(self._handle)
            if array_handle in self._arrays:
                dfgl = self._arrays[array_handle]
            else:
                dfgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_mixer_data__array__dfgl)
                self._arrays[array_handle] = dfgl
            return dfgl
        
        @dfgl.setter
        def dfgl(self, dfgl):
            self.dfgl[...] = dfgl
        
        @property
        def bmix(self):
            """
            Element bmix ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 37
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__bmix(self._handle)
        
        @bmix.setter
        def bmix(self, bmix):
            _AresMainPy.f90wrap_mixer_data__set__bmix(self._handle, bmix)
        
        @property
        def amix(self):
            """
            Element amix ftype=real(dp) pytype=float
            
            
            Defined at Mixer_module.fpp line 37
            
            """
            return _AresMainPy.f90wrap_mixer_data__get__amix(self._handle)
        
        @amix.setter
        def amix(self, amix):
            _AresMainPy.f90wrap_mixer_data__set__amix(self._handle, amix)
        
        def __str__(self):
            ret = ['<mixer_data>{\n']
            ret.append('    nhmix : ')
            ret.append(repr(self.nhmix))
            ret.append(',\n    nhmin : ')
            ret.append(repr(self.nhmin))
            ret.append(',\n    nam : ')
            ret.append(repr(self.nam))
            ret.append(',\n    nitra : ')
            ret.append(repr(self.nitra))
            ret.append(',\n    alpha : ')
            ret.append(repr(self.alpha))
            ret.append(',\n    beta : ')
            ret.append(repr(self.beta))
            ret.append(',\n    w0 : ')
            ret.append(repr(self.w0))
            ret.append(',\n    dxl : ')
            ret.append(repr(self.dxl))
            ret.append(',\n    dfl : ')
            ret.append(repr(self.dfl))
            ret.append(',\n    voma : ')
            ret.append(repr(self.voma))
            ret.append(',\n    sp : ')
            ret.append(repr(self.sp))
            ret.append(',\n    kerker : ')
            ret.append(repr(self.kerker))
            ret.append(',\n    dxgl : ')
            ret.append(repr(self.dxgl))
            ret.append(',\n    dfgl : ')
            ret.append(repr(self.dfgl))
            ret.append(',\n    bmix : ')
            ret.append(repr(self.bmix))
            ret.append(',\n    amix : ')
            ret.append(repr(self.amix))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init_mixer_data():
        """
        init_mixer_data()
        
        
        Defined at Mixer_module.fpp lines 44-51
        
        
        """
        _AresMainPy.f90wrap_init_mixer_data()
    
    @staticmethod
    def init_mixer_data_per():
        """
        init_mixer_data_per()
        
        
        Defined at Mixer_module.fpp lines 54-129
        
        
        """
        _AresMainPy.f90wrap_init_mixer_data_per()
    
    @staticmethod
    def init_mixer_data_iso():
        """
        init_mixer_data_iso()
        
        
        Defined at Mixer_module.fpp lines 133-200
        
        
        """
        _AresMainPy.f90wrap_init_mixer_data_iso()
    
    @staticmethod
    def destroy_mixer():
        """
        destroy_mixer()
        
        
        Defined at Mixer_module.fpp lines 204-234
        
        
        """
        _AresMainPy.f90wrap_destroy_mixer()
    
    @staticmethod
    def mixing(iter, xout, xin):
        """
        res = mixing(iter, xout, xin)
        
        
        Defined at Mixer_module.fpp lines 238-321
        
        Parameters
        ----------
        iter : int
        xout : float array
        xin : float array
        
        Returns
        -------
        res : float
        
        """
        res = _AresMainPy.f90wrap_mixing(iter=iter, xout=xout, xin=xin)
        return res
    
    @staticmethod
    def mixing_iso(iter, xout, xin):
        """
        res = mixing_iso(iter, xout, xin)
        
        
        Defined at Mixer_module.fpp lines 325-382
        
        Parameters
        ----------
        iter : int
        xout : float array
        xin : float array
        
        Returns
        -------
        res : float
        
        """
        res = _AresMainPy.f90wrap_mixing_iso(iter=iter, xout=xout, xin=xin)
        return res
    
    @staticmethod
    def anderson_mixing(iiter, xout, xin):
        """
        err = anderson_mixing(iiter, xout, xin)
        
        
        Defined at Mixer_module.fpp lines 389-512
        
        Parameters
        ----------
        iiter : int
        xout : float array
        xin : float array
        
        Returns
        -------
        err : float
        
        """
        err = _AresMainPy.f90wrap_anderson_mixing(iiter=iiter, xout=xout, xin=xin)
        return err
    
    @staticmethod
    def om1c(nam, nuh, sp, dfp, voma):
        """
        om1c(nam, nuh, sp, dfp, voma)
        
        
        Defined at Mixer_module.fpp lines 516-571
        
        Parameters
        ----------
        nam : int
        nuh : int
        sp : float
        dfp : float array
        voma : float array
        
        """
        _AresMainPy.f90wrap_om1c(nam=nam, nuh=nuh, sp=sp, dfp=dfp, voma=voma)
    
    @staticmethod
    def amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn):
        """
        amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn)
        
        
        Defined at Mixer_module.fpp lines 575-685
        
        Parameters
        ----------
        beta : float
        w0 : float
        nam : int
        nuh : int
        dxp : float array
        dfp : float array
        sp : float
        xl : float array
        fl : float array
        voma : float array
        xn : float array
        
        """
        _AresMainPy.f90wrap_amst(beta=beta, w0=w0, nam=nam, nuh=nuh, dxp=dxp, dfp=dfp, \
            sp=sp, xl=xl, fl=fl, voma=voma, xn=xn)
    
    @staticmethod
    def rpulay_mixing(iter, xout, xin):
        """
        err = rpulay_mixing(iter, xout, xin)
        
        
        Defined at Mixer_module.fpp lines 693-757
        
        Parameters
        ----------
        iter : int
        xout : float array
        xin : float array
        
        Returns
        -------
        err : float
        
        ----------------------
        store RL
        """
        err = _AresMainPy.f90wrap_rpulay_mixing(iter=iter, xout=xout, xin=xin)
        return err
    
    @staticmethod
    def rpulay_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn):
        """
        rpulay_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn)
        
        
        Defined at Mixer_module.fpp lines 761-799
        
        Parameters
        ----------
        beta : float
        w0 : float
        dime : int
        nh : int
        dxl : float array
        drl : float array
        xl : float array
        rl : float array
        xn : float array
        
        """
        _AresMainPy.f90wrap_rpulay_mix(beta=beta, w0=w0, dime=dime, nh=nh, dxl=dxl, \
            drl=drl, xl=xl, rl=rl, xn=xn)
    
    @staticmethod
    def init_kerker():
        """
        init_kerker()
        
        
        Defined at Mixer_module.fpp lines 892-909
        
        
        """
        _AresMainPy.f90wrap_init_kerker()
    
    @staticmethod
    def rpulayk_mixing(iter, rlg, xing):
        """
        rpulayk_mixing(iter, rlg, xing)
        
        
        Defined at Mixer_module.fpp lines 917-978
        
        Parameters
        ----------
        iter : int
        rlg : complex array
        xing : complex array
        
        ----------------------
        store RL
        """
        _AresMainPy.f90wrap_rpulayk_mixing(iter=iter, rlg=rlg, xing=xing)
    
    @staticmethod
    def rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn):
        """
        rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn)
        
        
        Defined at Mixer_module.fpp lines 982-1037
        
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
        _AresMainPy.f90wrap_rpulayk_mix(beta=beta, w0=w0, dime=dime, nh=nh, dxl=dxl, \
            drl=drl, xl=xl, rl=rl, xn=xn)
    
    @staticmethod
    def resta_mixing(iter, rlg, xing):
        """
        resta_mixing(iter, rlg, xing)
        
        
        Defined at Mixer_module.fpp lines 1043-1107
        
        Parameters
        ----------
        iter : int
        rlg : complex array
        xing : complex array
        
        ----------------------
        store RL
        """
        _AresMainPy.f90wrap_resta_mixing(iter=iter, rlg=rlg, xing=xing)
    
    @staticmethod
    def init_resta():
        """
        init_resta()
        
        
        Defined at Mixer_module.fpp lines 1110-1145
        
        
        """
        _AresMainPy.f90wrap_init_resta()
    
    _dt_array_initialisers = []
    

mixer_module = Mixer_Module()

class Scf_Module(f90wrap.runtime.FortranModule):
    """
    Module scf_module
    
    
    Defined at Scf_module.fpp lines 5-1462
    
    """
    @staticmethod
    def arpackscf(rhos, psi, eval):
        """
        arpackscf(rhos, psi, eval)
        
        
        Defined at Scf_module.fpp lines 21-43
        
        Parameters
        ----------
        rhos : float array
        psi : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_arpackscf(rhos=rhos, psi=psi, eval=eval)
    
    @staticmethod
    def solver_spin(rhos, psi, eval, nev, diagtol):
        """
        solver_spin(rhos, psi, eval, nev, diagtol)
        
        
        Defined at Scf_module.fpp lines 46-66
        
        Parameters
        ----------
        rhos : float array
        psi : complex array
        eval : float array
        nev : int
        diagtol : float
        
        """
        _AresMainPy.f90wrap_solver_spin(rhos=rhos, psi=psi, eval=eval, nev=nev, \
            diagtol=diagtol)
    
    @staticmethod
    def ksolver(veff, kpsi, keval, nev, tol):
        """
        ksolver(veff, kpsi, keval, nev, tol)
        
        
        Defined at Scf_module.fpp lines 69-113
        
        Parameters
        ----------
        veff : float array
        kpsi : complex array
        keval : float array
        nev : int
        tol : float
        
        """
        _AresMainPy.f90wrap_ksolver(veff=veff, kpsi=kpsi, keval=keval, nev=nev, tol=tol)
    
    @staticmethod
    def smear_updaterho(nev, ne, psi, eval, wke_l, rhos):
        """
        fme, ets = smear_updaterho(nev, ne, psi, eval, wke_l, rhos)
        
        
        Defined at Scf_module.fpp lines 117-174
        
        Parameters
        ----------
        nev : int
        ne : float
        psi : complex array
        eval : float array
        wke_l : float array
        rhos : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_smear_updaterho(nev=nev, ne=ne, psi=psi, \
            eval=eval, wke_l=wke_l, rhos=rhos)
        return fme, ets
    
    @staticmethod
    def chefsi(rhos, psi, eval):
        """
        chefsi(rhos, psi, eval)
        
        
        Defined at Scf_module.fpp lines 182-386
        
        Parameters
        ----------
        rhos : float array
        psi : complex array
        eval : float array
        
        ===============================================================
         WRITE(6,"(1X,A8,I4,1X,A6,E15.7,1X,A4,E15.7)") &
              &      '>CheFSI:',iter,'dTOTEN',dtoten,      &
              &  'dRHO',drho
         IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
        ===============================================================
        """
        _AresMainPy.f90wrap_chefsi(rhos=rhos, psi=psi, eval=eval)
    
    @staticmethod
    def filter_spin(veff, psi, eval):
        """
        filter_spin(veff, psi, eval)
        
        
        Defined at Scf_module.fpp lines 389-411
        
        Parameters
        ----------
        veff : float array
        psi : complex array
        eval : float array
        
        """
        _AresMainPy.f90wrap_filter_spin(veff=veff, psi=psi, eval=eval)
    
    @staticmethod
    def arpackscf_r(rhos, psi, eval):
        """
        arpackscf_r(rhos, psi, eval)
        
        
        Defined at Scf_module.fpp lines 435-502
        
        Parameters
        ----------
        rhos : float array
        psi : float array
        eval : float array
        
        ===============================================================
        """
        _AresMainPy.f90wrap_arpackscf_r(rhos=rhos, psi=psi, eval=eval)
    
    @staticmethod
    def solver_spin_r(rhos, psi, eval, nev, diagtol):
        """
        solver_spin_r(rhos, psi, eval, nev, diagtol)
        
        
        Defined at Scf_module.fpp lines 505-558
        
        Parameters
        ----------
        rhos : float array
        psi : float array
        eval : float array
        nev : int
        diagtol : float
        
        """
        _AresMainPy.f90wrap_solver_spin_r(rhos=rhos, psi=psi, eval=eval, nev=nev, \
            diagtol=diagtol)
    
    @staticmethod
    def real_smear_updaterho(nev, ne, psi, eval, wke_l, rhos):
        """
        fme, ets = real_smear_updaterho(nev, ne, psi, eval, wke_l, rhos)
        
        
        Defined at Scf_module.fpp lines 561-606
        
        Parameters
        ----------
        nev : int
        ne : float
        psi : float array
        eval : float array
        wke_l : float array
        rhos : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_real_smear_updaterho(nev=nev, ne=ne, psi=psi, \
            eval=eval, wke_l=wke_l, rhos=rhos)
        return fme, ets
    
    @staticmethod
    def bvk_chefsi(rhos, psi, eval):
        """
        bvk_chefsi(rhos, psi, eval)
        
        
        Defined at Scf_module.fpp lines 613-710
        
        Parameters
        ----------
        rhos : float array
        psi : float array
        eval : float array
        
        ================================================
        ##XLT test I/O FOR DETECT ISO
        open(unit=1118,file="RHO_ISO")
        READ(1118,*)
        READ(1118,*)rhoS
        close(1118)
        goto 1188
        ================================================
        initial mixer
        """
        _AresMainPy.f90wrap_bvk_chefsi(rhos=rhos, psi=psi, eval=eval)
    
    @staticmethod
    def bvk_filter_spin(rhos, x, d):
        """
        bvk_filter_spin(rhos, x, d)
        
        
        Defined at Scf_module.fpp lines 713-732
        
        Parameters
        ----------
        rhos : float array
        x : float array
        d : float array
        
        """
        _AresMainPy.f90wrap_bvk_filter_spin(rhos=rhos, x=x, d=d)
    
    @staticmethod
    def prr_filter_spin(nfs, nfe, rhos, x, efr, pbar):
        """
        fme, ets = prr_filter_spin(nfs, nfe, rhos, x, efr, pbar)
        
        
        Defined at Scf_module.fpp lines 739-827
        
        Parameters
        ----------
        nfs : int
        nfe : float
        rhos : float array
        x : float array
        efr : float array
        pbar : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_prr_filter_spin(nfs=nfs, nfe=nfe, rhos=rhos, x=x, \
            efr=efr, pbar=pbar)
        return fme, ets
    
    @staticmethod
    def prr_chefsi(rhos, phi):
        """
        prr_chefsi(rhos, phi)
        
        
        Defined at Scf_module.fpp lines 830-908
        
        Parameters
        ----------
        rhos : float array
        phi : float array
        
        -------------------------------------------------------
        """
        _AresMainPy.f90wrap_prr_chefsi(rhos=rhos, phi=phi)
    
    @staticmethod
    def iso_solver_spin_r(rhos, psi, eval, nev, diagtol):
        """
        iso_solver_spin_r(rhos, psi, eval, nev, diagtol)
        
        
        Defined at Scf_module.fpp lines 918-971
        
        Parameters
        ----------
        rhos : float array
        psi : float array
        eval : float array
        nev : int
        diagtol : float
        
        """
        _AresMainPy.f90wrap_iso_solver_spin_r(rhos=rhos, psi=psi, eval=eval, nev=nev, \
            diagtol=diagtol)
    
    @staticmethod
    def iso_smear_updaterho(nev, ne, psi, eval, wke_l, rhos_out):
        """
        fme, ets = iso_smear_updaterho(nev, ne, psi, eval, wke_l, rhos_out)
        
        
        Defined at Scf_module.fpp lines 974-1023
        
        Parameters
        ----------
        nev : int
        ne : float
        psi : float array
        eval : float array
        wke_l : float array
        rhos_out : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_iso_smear_updaterho(nev=nev, ne=ne, psi=psi, \
            eval=eval, wke_l=wke_l, rhos_out=rhos_out)
        return fme, ets
    
    @staticmethod
    def iso_chefsi(rhos, psi, eval):
        """
        iso_chefsi(rhos, psi, eval)
        
        
        Defined at Scf_module.fpp lines 1030-1204
        
        Parameters
        ----------
        rhos : float array
        psi : float array
        eval : float array
        
        ================================================
        xlt 2018-5-18 USED FOR TEST
         USE grid_module , ONLY : ISO_Rho2grid,rho_calc
        ================================================
        """
        _AresMainPy.f90wrap_iso_chefsi(rhos=rhos, psi=psi, eval=eval)
    
    @staticmethod
    def iso_filter_spin(rhos, x, d):
        """
        iso_filter_spin(rhos, x, d)
        
        
        Defined at Scf_module.fpp lines 1207-1238
        
        Parameters
        ----------
        rhos : float array
        x : float array
        d : float array
        
        =============================================
        print *,'calculate effectial time -->',(tf-te)/10000.d0
        print *,'filter GeneralRealyRitz  time -->',(tg-tf)/10000.d0
        print *,'total chebyshev filter time -->',(tg-te)/10000.d0
        =============================================
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        """
        _AresMainPy.f90wrap_iso_filter_spin(rhos=rhos, x=x, d=d)
    
    @staticmethod
    def chefsi_gamma(rhos, psi, eval):
        """
        chefsi_gamma(rhos, psi, eval)
        
        
        Defined at Scf_module.fpp lines 1241-1378
        
        Parameters
        ----------
        rhos : float array
        psi : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_chefsi_gamma(rhos=rhos, psi=psi, eval=eval)
    
    @staticmethod
    def smear_updaterho_gamma(nev, ne, psi, eval, wke_l, rhos):
        """
        fme, ets = smear_updaterho_gamma(nev, ne, psi, eval, wke_l, rhos)
        
        
        Defined at Scf_module.fpp lines 1381-1438
        
        Parameters
        ----------
        nev : int
        ne : float
        psi : float array
        eval : float array
        wke_l : float array
        rhos : float array
        
        Returns
        -------
        fme : float
        ets : float
        
        """
        fme, ets = _AresMainPy.f90wrap_smear_updaterho_gamma(nev=nev, ne=ne, psi=psi, \
            eval=eval, wke_l=wke_l, rhos=rhos)
        return fme, ets
    
    @staticmethod
    def filter_spin_gamma(veff, psi, eval):
        """
        filter_spin_gamma(veff, psi, eval)
        
        
        Defined at Scf_module.fpp lines 1441-1461
        
        Parameters
        ----------
        veff : float array
        psi : float array
        eval : float array
        
        """
        _AresMainPy.f90wrap_filter_spin_gamma(veff=veff, psi=psi, eval=eval)
    
    _dt_array_initialisers = []
    

scf_module = Scf_Module()

class Out_Module(f90wrap.runtime.FortranModule):
    """
    Module out_module
    
    
    Defined at Out_module.fpp lines 5-130
    
    """
    @staticmethod
    def ksout_list():
        """
        ksout_list()
        
        
        Defined at Out_module.fpp lines 9-16
        
        
        """
        _AresMainPy.f90wrap_ksout_list()
    
    @staticmethod
    def ksout_list_per():
        """
        ksout_list_per()
        
        
        Defined at Out_module.fpp lines 19-69
        
        
        """
        _AresMainPy.f90wrap_ksout_list_per()
    
    @staticmethod
    def ksout_list_iso():
        """
        ksout_list_iso()
        
        
        Defined at Out_module.fpp lines 72-130
        
        
        """
        _AresMainPy.f90wrap_ksout_list_iso()
    
    _dt_array_initialisers = []
    

out_module = Out_Module()

class Forcestress_Module(f90wrap.runtime.FortranModule):
    """
    Module forcestress_module
    
    
    Defined at ForceStress_module.fpp lines 10-2319
    
    """
    @staticmethod
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
        _AresMainPy.f90wrap_cal_force_c(rhos=rhos, uik=uik, force=force)
    
    @staticmethod
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
        _AresMainPy.f90wrap_cal_force_gamma(rhos=rhos, uik=uik, force=force)
    
    @staticmethod
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
        _AresMainPy.f90wrap_cal_force_r(rhos=rhos, uik=uik, force=force)
    
    @staticmethod
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
        _AresMainPy.f90wrap_locforce(rhos=rhos, lforce=lforce)
    
    @staticmethod
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
        _AresMainPy.f90wrap_locforce_r(rhos=rhos, lforce=lforce)
    
    @staticmethod
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
        _AresMainPy.f90wrap_nonlforce(uik=uik, nlforce=nlforce)
    
    @staticmethod
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
        _AresMainPy.f90wrap_nonlforce_gamma(uik=uik, nlforce=nlforce)
    
    @staticmethod
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
        _AresMainPy.f90wrap_nonlforce_r_dg(uik=uik, nlforce=nlforce)
    
    @staticmethod
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
        _AresMainPy.f90wrap_nonlforce_r(uik=uik, nlforce=nlforce)
    
    @staticmethod
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
        _AresMainPy.f90wrap_cal_stress(rhos=rhos, uik=uik, stress=stress)
    
    @staticmethod
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
        _AresMainPy.f90wrap_cal_stress_gamma(rhos=rhos, uik=uik, stress=stress)
    
    @staticmethod
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
        lda_stress = _AresMainPy.f90wrap_lda_stress(vxc=vxc, rho=rho, elda=elda)
        return lda_stress
    
    @staticmethod
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
        hart_stress = _AresMainPy.f90wrap_hart_stress(rhorecip=rhorecip, eh=eh)
        return hart_stress
    
    @staticmethod
    def kin_stress(uik, kinstress):
        """
        kin_stress(uik, kinstress)
        
        
        Defined at ForceStress_module.fpp lines 1426-1565
        
        Parameters
        ----------
        uik : complex array
        kinstress : float array
        
        """
        _AresMainPy.f90wrap_kin_stress(uik=uik, kinstress=kinstress)
    
    @staticmethod
    def kin_stress_gamma(uik, kinstress):
        """
        kin_stress_gamma(uik, kinstress)
        
        
        Defined at ForceStress_module.fpp lines 1567-1706
        
        Parameters
        ----------
        uik : float array
        kinstress : float array
        
        """
        _AresMainPy.f90wrap_kin_stress_gamma(uik=uik, kinstress=kinstress)
    
    @staticmethod
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
        _AresMainPy.f90wrap_ion_nl_stress(uik=uik, nlstress=nlstress)
    
    @staticmethod
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
        _AresMainPy.f90wrap_ion_nl_stress_gamma(uik=uik, nlstress=nlstress)
    
    @staticmethod
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
        ionele_stress = _AresMainPy.f90wrap_ionele_stress(rhorecip=rhorecip, \
            energy=energy)
        return ionele_stress
    
    @staticmethod
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
        pseudopotdifflookup = _AresMainPy.f90wrap_pseudopotdifflookup(ity=ity, \
            qnorm=qnorm)
        return pseudopotdifflookup
    
    @staticmethod
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
        _AresMainPy.f90wrap_cal_force_stress()
    
    @property
    def cellpress(self):
        """
        Element cellpress ftype=real(dp) pytype=float
        
        
        Defined at ForceStress_module.fpp line 13
        
        """
        return _AresMainPy.f90wrap_forcestress_module__get__cellpress()
    
    @cellpress.setter
    def cellpress(self, cellpress):
        _AresMainPy.f90wrap_forcestress_module__set__cellpress(cellpress)
    
    def __str__(self):
        ret = ['<forcestress_module>{\n']
        ret.append('    cellpress : ')
        ret.append(repr(self.cellpress))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

forcestress_module = Forcestress_Module()

class Write_Module(f90wrap.runtime.FortranModule):
    """
    Module write_module
    
    
    Defined at Write_module.fpp lines 5-203
    
    """
    @staticmethod
    def write_poscar(filename, lat_mat, eleid, pos, atom_symbol=None, fixpos=None):
        """
        write_poscar(filename, lat_mat, eleid, pos[, atom_symbol, fixpos])
        
        
        Defined at Write_module.fpp lines 22-62
        
        Parameters
        ----------
        filename : str
        lat_mat : float array
        eleid : int array
        pos : float array
        atom_symbol : str array
        fixpos : bool array
        
        """
        _AresMainPy.f90wrap_write_poscar(filename=filename, lat_mat=lat_mat, \
            eleid=eleid, pos=pos, atom_symbol=atom_symbol, fixpos=fixpos)
    
    @staticmethod
    def write_cif(filename, lat_para, nati, pos, atom_symbol=None):
        """
        write_cif(filename, lat_para, nati, pos[, atom_symbol])
        
        
        Defined at Write_module.fpp lines 66-131
        
        Parameters
        ----------
        filename : str
        lat_para : float array
        nati : int array
        pos : float array
        atom_symbol : str array
        
        """
        _AresMainPy.f90wrap_write_cif(filename=filename, lat_para=lat_para, nati=nati, \
            pos=pos, atom_symbol=atom_symbol)
    
    @staticmethod
    def write3dat(filename, rho):
        """
        write3dat(filename, rho)
        
        
        Defined at Write_module.fpp lines 135-164
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        -----------------------------------------------------------------------
        """
        _AresMainPy.f90wrap_write3dat(filename=filename, rho=rho)
    
    @staticmethod
    def write3d(filename, rho):
        """
        write3d(filename, rho)
        
        
        Defined at Write_module.fpp lines 168-192
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        -----------------------------------------------------------------------
        """
        _AresMainPy.f90wrap_write3d(filename=filename, rho=rho)
    
    @staticmethod
    def writedensity(filename, rho):
        """
        writedensity(filename, rho)
        
        
        Defined at Write_module.fpp lines 198-201
        
        Parameters
        ----------
        filename : str
        rho : float array
        
        """
        _AresMainPy.f90wrap_writedensity(filename=filename, rho=rho)
    
    _dt_array_initialisers = []
    

write_module = Write_Module()

class Dyn_Module(f90wrap.runtime.FortranModule):
    """
    Module dyn_module
    
    
    Defined at dyn_module.fpp lines 5-44
    
    """
    @f90wrap.runtime.register_class("AresMainPy.dynamics")
    class dynamics(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=dynamics)
        
        
        Defined at dyn_module.fpp lines 8-23
        
        """
        def __init__(self, handle=None):
            """
            self = Dynamics()
            
            
            Defined at dyn_module.fpp lines 8-23
            
            
            Returns
            -------
            this : Dynamics
            	Object to be constructed
            
            
            Automatically generated constructor for dynamics
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_dynamics_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Dynamics
            
            
            Defined at dyn_module.fpp lines 8-23
            
            Parameters
            ----------
            this : Dynamics
            	Object to be destructed
            
            
            Automatically generated destructor for dynamics
            """
            if self._alloc:
                _AresMainPy.f90wrap_dynamics_finalise(this=self._handle)
        
        @property
        def lstop(self):
            """
            Element lstop ftype=logical pytype=bool
            
            
            Defined at dyn_module.fpp line 9
            
            """
            return _AresMainPy.f90wrap_dynamics__get__lstop(self._handle)
        
        @lstop.setter
        def lstop(self, lstop):
            _AresMainPy.f90wrap_dynamics__set__lstop(self._handle, lstop)
        
        @property
        def posion(self):
            """
            Element posion ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 10
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__posion(self._handle)
            if array_handle in self._arrays:
                posion = self._arrays[array_handle]
            else:
                posion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__posion)
                self._arrays[array_handle] = posion
            return posion
        
        @posion.setter
        def posion(self, posion):
            self.posion[...] = posion
        
        @property
        def posioc(self):
            """
            Element posioc ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__posioc(self._handle)
            if array_handle in self._arrays:
                posioc = self._arrays[array_handle]
            else:
                posioc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__posioc)
                self._arrays[array_handle] = posioc
            return posioc
        
        @posioc.setter
        def posioc(self, posioc):
            self.posioc[...] = posioc
        
        @property
        def d2(self):
            """
            Element d2 ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__d2(self._handle)
            if array_handle in self._arrays:
                d2 = self._arrays[array_handle]
            else:
                d2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__d2)
                self._arrays[array_handle] = d2
            return d2
        
        @d2.setter
        def d2(self, d2):
            self.d2[...] = d2
        
        @property
        def d2c(self):
            """
            Element d2c ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__d2c(self._handle)
            if array_handle in self._arrays:
                d2c = self._arrays[array_handle]
            else:
                d2c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__d2c)
                self._arrays[array_handle] = d2c
            return d2c
        
        @d2c.setter
        def d2c(self, d2c):
            self.d2c[...] = d2c
        
        @property
        def d3(self):
            """
            Element d3 ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__d3(self._handle)
            if array_handle in self._arrays:
                d3 = self._arrays[array_handle]
            else:
                d3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__d3)
                self._arrays[array_handle] = d3
            return d3
        
        @d3.setter
        def d3(self, d3):
            self.d3[...] = d3
        
        @property
        def a(self):
            """
            Element a ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__a(self._handle)
            if array_handle in self._arrays:
                a = self._arrays[array_handle]
            else:
                a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__a)
                self._arrays[array_handle] = a
            return a
        
        @a.setter
        def a(self, a):
            self.a[...] = a
        
        @property
        def b(self):
            """
            Element b ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__b(self._handle)
            if array_handle in self._arrays:
                b = self._arrays[array_handle]
            else:
                b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__b)
                self._arrays[array_handle] = b
            return b
        
        @b.setter
        def b(self, b):
            self.b[...] = b
        
        @property
        def ac(self):
            """
            Element ac ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_dynamics__array__ac(self._handle)
            if array_handle in self._arrays:
                ac = self._arrays[array_handle]
            else:
                ac = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_dynamics__array__ac)
                self._arrays[array_handle] = ac
            return ac
        
        @ac.setter
        def ac(self, ac):
            self.ac[...] = ac
        
        @property
        def potim(self):
            """
            Element potim ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 18
            
            """
            return _AresMainPy.f90wrap_dynamics__get__potim(self._handle)
        
        @potim.setter
        def potim(self, potim):
            _AresMainPy.f90wrap_dynamics__set__potim(self._handle, potim)
        
        @property
        def ediffg(self):
            """
            Element ediffg ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 19
            
            """
            return _AresMainPy.f90wrap_dynamics__get__ediffg(self._handle)
        
        @ediffg.setter
        def ediffg(self, ediffg):
            _AresMainPy.f90wrap_dynamics__set__ediffg(self._handle, ediffg)
        
        @property
        def pstress(self):
            """
            Element pstress ftype=real(dp) pytype=float
            
            
            Defined at dyn_module.fpp line 20
            
            """
            return _AresMainPy.f90wrap_dynamics__get__pstress(self._handle)
        
        @pstress.setter
        def pstress(self, pstress):
            _AresMainPy.f90wrap_dynamics__set__pstress(self._handle, pstress)
        
        @property
        def ibrion(self):
            """
            Element ibrion ftype=integer(i4b) pytype=int
            
            
            Defined at dyn_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_dynamics__get__ibrion(self._handle)
        
        @ibrion.setter
        def ibrion(self, ibrion):
            _AresMainPy.f90wrap_dynamics__set__ibrion(self._handle, ibrion)
        
        @property
        def isif(self):
            """
            Element isif ftype=integer(i4b) pytype=int
            
            
            Defined at dyn_module.fpp line 22
            
            """
            return _AresMainPy.f90wrap_dynamics__get__isif(self._handle)
        
        @isif.setter
        def isif(self, isif):
            _AresMainPy.f90wrap_dynamics__set__isif(self._handle, isif)
        
        @property
        def nfree(self):
            """
            Element nfree ftype=integer(i4b) pytype=int
            
            
            Defined at dyn_module.fpp line 23
            
            """
            return _AresMainPy.f90wrap_dynamics__get__nfree(self._handle)
        
        @nfree.setter
        def nfree(self, nfree):
            _AresMainPy.f90wrap_dynamics__set__nfree(self._handle, nfree)
        
        def __str__(self):
            ret = ['<dynamics>{\n']
            ret.append('    lstop : ')
            ret.append(repr(self.lstop))
            ret.append(',\n    posion : ')
            ret.append(repr(self.posion))
            ret.append(',\n    posioc : ')
            ret.append(repr(self.posioc))
            ret.append(',\n    d2 : ')
            ret.append(repr(self.d2))
            ret.append(',\n    d2c : ')
            ret.append(repr(self.d2c))
            ret.append(',\n    d3 : ')
            ret.append(repr(self.d3))
            ret.append(',\n    a : ')
            ret.append(repr(self.a))
            ret.append(',\n    b : ')
            ret.append(repr(self.b))
            ret.append(',\n    ac : ')
            ret.append(repr(self.ac))
            ret.append(',\n    potim : ')
            ret.append(repr(self.potim))
            ret.append(',\n    ediffg : ')
            ret.append(repr(self.ediffg))
            ret.append(',\n    pstress : ')
            ret.append(repr(self.pstress))
            ret.append(',\n    ibrion : ')
            ret.append(repr(self.ibrion))
            ret.append(',\n    isif : ')
            ret.append(repr(self.isif))
            ret.append(',\n    nfree : ')
            ret.append(repr(self.nfree))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def create_dyn(na, dyn):
        """
        create_dyn(na, dyn)
        
        
        Defined at dyn_module.fpp lines 27-35
        
        Parameters
        ----------
        na : int
        dyn : Dynamics
        
        """
        _AresMainPy.f90wrap_create_dyn(na=na, dyn=dyn._handle)
    
    @staticmethod
    def destroy_dyn(self):
        """
        destroy_dyn(self)
        
        
        Defined at dyn_module.fpp lines 37-44
        
        Parameters
        ----------
        dyn : Dynamics
        
        """
        _AresMainPy.f90wrap_destroy_dyn(dyn=self._handle)
    
    _dt_array_initialisers = []
    

dyn_module = Dyn_Module()

class Cg_Relax(f90wrap.runtime.FortranModule):
    """
    Module cg_relax
    
    
    Defined at cg_module.fpp lines 5-289
    
    """
    @f90wrap.runtime.register_class("AresMainPy.lattice")
    class lattice(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=lattice)
        
        
        Defined at cg_module.fpp lines 16-20
        
        """
        def __init__(self, handle=None):
            """
            self = Lattice()
            
            
            Defined at cg_module.fpp lines 16-20
            
            
            Returns
            -------
            this : Lattice
            	Object to be constructed
            
            
            Automatically generated constructor for lattice
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_lattice_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Lattice
            
            
            Defined at cg_module.fpp lines 16-20
            
            Parameters
            ----------
            this : Lattice
            	Object to be destructed
            
            
            Automatically generated destructor for lattice
            """
            if self._alloc:
                _AresMainPy.f90wrap_lattice_finalise(this=self._handle)
        
        @property
        def a(self):
            """
            Element a ftype=real(dp) pytype=float
            
            
            Defined at cg_module.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_lattice__array__a(self._handle)
            if array_handle in self._arrays:
                a = self._arrays[array_handle]
            else:
                a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_lattice__array__a)
                self._arrays[array_handle] = a
            return a
        
        @a.setter
        def a(self, a):
            self.a[...] = a
        
        @property
        def b(self):
            """
            Element b ftype=real(dp) pytype=float
            
            
            Defined at cg_module.fpp line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_lattice__array__b(self._handle)
            if array_handle in self._arrays:
                b = self._arrays[array_handle]
            else:
                b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_lattice__array__b)
                self._arrays[array_handle] = b
            return b
        
        @b.setter
        def b(self, b):
            self.b[...] = b
        
        @property
        def anorm(self):
            """
            Element anorm ftype=real(dp) pytype=float
            
            
            Defined at cg_module.fpp line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_lattice__array__anorm(self._handle)
            if array_handle in self._arrays:
                anorm = self._arrays[array_handle]
            else:
                anorm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_lattice__array__anorm)
                self._arrays[array_handle] = anorm
            return anorm
        
        @anorm.setter
        def anorm(self, anorm):
            self.anorm[...] = anorm
        
        @property
        def bnorm(self):
            """
            Element bnorm ftype=real(dp) pytype=float
            
            
            Defined at cg_module.fpp line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_lattice__array__bnorm(self._handle)
            if array_handle in self._arrays:
                bnorm = self._arrays[array_handle]
            else:
                bnorm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_lattice__array__bnorm)
                self._arrays[array_handle] = bnorm
            return bnorm
        
        @bnorm.setter
        def bnorm(self, bnorm):
            self.bnorm[...] = bnorm
        
        @property
        def omega(self):
            """
            Element omega ftype=real(dp) pytype=float
            
            
            Defined at cg_module.fpp line 20
            
            """
            return _AresMainPy.f90wrap_lattice__get__omega(self._handle)
        
        @omega.setter
        def omega(self, omega):
            _AresMainPy.f90wrap_lattice__set__omega(self._handle, omega)
        
        def __str__(self):
            ret = ['<lattice>{\n']
            ret.append('    a : ')
            ret.append(repr(self.a))
            ret.append(',\n    b : ')
            ret.append(repr(self.b))
            ret.append(',\n    anorm : ')
            ret.append(repr(self.anorm))
            ret.append(',\n    bnorm : ')
            ret.append(repr(self.bnorm))
            ret.append(',\n    omega : ')
            ret.append(repr(self.omega))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def cg_relax_vasp_interface(nstep, lfopt, lopt, lexceed):
        """
        cg_relax_vasp_interface(nstep, lfopt, lopt, lexceed)
        
        
        Defined at cg_module.fpp lines 23-196
        
        Parameters
        ----------
        nstep : int
        lfopt : bool
        lopt : bool
        lexceed : bool
        
        ----------------------------------------------------------------
        """
        _AresMainPy.f90wrap_cg_relax_vasp_interface(nstep=nstep, lfopt=lfopt, lopt=lopt, \
            lexceed=lexceed)
    
    @staticmethod
    def check_opt(na, flag):
        """
        check_opt(na, flag)
        
        
        Defined at cg_module.fpp lines 199-236
        
        Parameters
        ----------
        na : int
        flag : bool
        
        """
        _AresMainPy.f90wrap_check_opt(na=na, flag=flag)
    
    @staticmethod
    def lattic(self):
        """
        lattic(self)
        
        
        Defined at cg_module.fpp lines 239-259
        
        Parameters
        ----------
        mylatt : Lattice
        
        """
        _AresMainPy.f90wrap_lattic(mylatt=self._handle)
    
    @staticmethod
    def expro(h, u1, u2):
        """
        expro(h, u1, u2)
        
        
        Defined at cg_module.fpp lines 262-267
        
        Parameters
        ----------
        h : float array
        u1 : float array
        u2 : float array
        
        """
        _AresMainPy.f90wrap_expro(h=h, u1=u1, u2=u2)
    
    @staticmethod
    def check_distance(l, pos, lwrong):
        """
        check_distance(l, pos, lwrong)
        
        
        Defined at cg_module.fpp lines 270-287
        
        Parameters
        ----------
        l : float
        pos : float array
        lwrong : bool
        
        -----------------------------------------------
        """
        _AresMainPy.f90wrap_check_distance(l=l, pos=pos, lwrong=lwrong)
    
    @property
    def maxforce(self):
        """
        Element maxforce ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 13
        
        """
        return _AresMainPy.f90wrap_cg_relax__get__maxforce()
    
    @maxforce.setter
    def maxforce(self, maxforce):
        _AresMainPy.f90wrap_cg_relax__set__maxforce(maxforce)
    
    @property
    def maxstress(self):
        """
        Element maxstress ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 13
        
        """
        return _AresMainPy.f90wrap_cg_relax__get__maxstress()
    
    @maxstress.setter
    def maxstress(self, maxstress):
        _AresMainPy.f90wrap_cg_relax__set__maxstress(maxstress)
    
    @property
    def tstress(self):
        """
        Element tstress ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_cg_relax__array__tstress(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            tstress = self._arrays[array_handle]
        else:
            tstress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_cg_relax__array__tstress)
            self._arrays[array_handle] = tstress
        return tstress
    
    @tstress.setter
    def tstress(self, tstress):
        self.tstress[...] = tstress
    
    @property
    def latmark(self):
        """
        Element latmark ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy.f90wrap_cg_relax__array__latmark(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            latmark = self._arrays[array_handle]
        else:
            latmark = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _AresMainPy.f90wrap_cg_relax__array__latmark)
            self._arrays[array_handle] = latmark
        return latmark
    
    @latmark.setter
    def latmark(self, latmark):
        self.latmark[...] = latmark
    
    def __str__(self):
        ret = ['<cg_relax>{\n']
        ret.append('    maxforce : ')
        ret.append(repr(self.maxforce))
        ret.append(',\n    maxstress : ')
        ret.append(repr(self.maxstress))
        ret.append(',\n    tstress : ')
        ret.append(repr(self.tstress))
        ret.append(',\n    latmark : ')
        ret.append(repr(self.latmark))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

cg_relax = Cg_Relax()

class Relax_Module(f90wrap.runtime.FortranModule):
    """
    Module relax_module
    
    
    Defined at Relax_module.fpp lines 10-305
    
    """
    @f90wrap.runtime.register_class("AresMainPy.relax_type")
    class relax_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=relax_type)
        
        
        Defined at Relax_module.fpp lines 17-29
        
        """
        def __init__(self, handle=None):
            """
            self = Relax_Type()
            
            
            Defined at Relax_module.fpp lines 17-29
            
            
            Returns
            -------
            this : Relax_Type
            	Object to be constructed
            
            
            Automatically generated constructor for relax_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_relax_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Relax_Type
            
            
            Defined at Relax_module.fpp lines 17-29
            
            Parameters
            ----------
            this : Relax_Type
            	Object to be destructed
            
            
            Automatically generated destructor for relax_type
            """
            if self._alloc:
                _AresMainPy.f90wrap_relax_type_finalise(this=self._handle)
        
        @property
        def timsmax(self):
            """
            Element timsmax ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_relax_type__get__timsmax(self._handle)
        
        @timsmax.setter
        def timsmax(self, timsmax):
            _AresMainPy.f90wrap_relax_type__set__timsmax(self._handle, timsmax)
        
        @property
        def tims(self):
            """
            Element tims ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_relax_type__get__tims(self._handle)
        
        @tims.setter
        def tims(self, tims):
            _AresMainPy.f90wrap_relax_type__set__tims(self._handle, tims)
        
        @property
        def alpha(self):
            """
            Element alpha ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 21
            
            """
            return _AresMainPy.f90wrap_relax_type__get__alpha(self._handle)
        
        @alpha.setter
        def alpha(self, alpha):
            _AresMainPy.f90wrap_relax_type__set__alpha(self._handle, alpha)
        
        @property
        def faction(self):
            """
            Element faction ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_relax_type__array__faction(self._handle)
            if array_handle in self._arrays:
                faction = self._arrays[array_handle]
            else:
                faction = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_relax_type__array__faction)
                self._arrays[array_handle] = faction
            return faction
        
        @faction.setter
        def faction(self, faction):
            self.faction[...] = faction
        
        @property
        def fiond(self):
            """
            Element fiond ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_relax_type__array__fiond(self._handle)
            if array_handle in self._arrays:
                fiond = self._arrays[array_handle]
            else:
                fiond = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_relax_type__array__fiond)
                self._arrays[array_handle] = fiond
            return fiond
        
        @fiond.setter
        def fiond(self, fiond):
            self.fiond[...] = fiond
        
        @property
        def fceld(self):
            """
            Element fceld ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_relax_type__array__fceld(self._handle)
            if array_handle in self._arrays:
                fceld = self._arrays[array_handle]
            else:
                fceld = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_relax_type__array__fceld)
                self._arrays[array_handle] = fceld
            return fceld
        
        @fceld.setter
        def fceld(self, fceld):
            self.fceld[...] = fceld
        
        @property
        def velion(self):
            """
            Element velion ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_relax_type__array__velion(self._handle)
            if array_handle in self._arrays:
                velion = self._arrays[array_handle]
            else:
                velion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_relax_type__array__velion)
                self._arrays[array_handle] = velion
            return velion
        
        @velion.setter
        def velion(self, velion):
            self.velion[...] = velion
        
        @property
        def velcel(self):
            """
            Element velcel ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_relax_type__array__velcel(self._handle)
            if array_handle in self._arrays:
                velcel = self._arrays[array_handle]
            else:
                velcel = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_relax_type__array__velcel)
                self._arrays[array_handle] = velcel
            return velcel
        
        @velcel.setter
        def velcel(self, velcel):
            self.velcel[...] = velcel
        
        @property
        def factcel(self):
            """
            Element factcel ftype=real(dp) pytype=float
            
            
            Defined at Relax_module.fpp line 27
            
            """
            return _AresMainPy.f90wrap_relax_type__get__factcel(self._handle)
        
        @factcel.setter
        def factcel(self, factcel):
            _AresMainPy.f90wrap_relax_type__set__factcel(self._handle, factcel)
        
        @property
        def neg(self):
            """
            Element neg ftype=integer(i4b) pytype=int
            
            
            Defined at Relax_module.fpp line 28
            
            """
            return _AresMainPy.f90wrap_relax_type__get__neg(self._handle)
        
        @neg.setter
        def neg(self, neg):
            _AresMainPy.f90wrap_relax_type__set__neg(self._handle, neg)
        
        @property
        def lneg(self):
            """
            Element lneg ftype=logical pytype=bool
            
            
            Defined at Relax_module.fpp line 29
            
            """
            return _AresMainPy.f90wrap_relax_type__get__lneg(self._handle)
        
        @lneg.setter
        def lneg(self, lneg):
            _AresMainPy.f90wrap_relax_type__set__lneg(self._handle, lneg)
        
        def __str__(self):
            ret = ['<relax_type>{\n']
            ret.append('    timsmax : ')
            ret.append(repr(self.timsmax))
            ret.append(',\n    tims : ')
            ret.append(repr(self.tims))
            ret.append(',\n    alpha : ')
            ret.append(repr(self.alpha))
            ret.append(',\n    faction : ')
            ret.append(repr(self.faction))
            ret.append(',\n    fiond : ')
            ret.append(repr(self.fiond))
            ret.append(',\n    fceld : ')
            ret.append(repr(self.fceld))
            ret.append(',\n    velion : ')
            ret.append(repr(self.velion))
            ret.append(',\n    velcel : ')
            ret.append(repr(self.velcel))
            ret.append(',\n    factcel : ')
            ret.append(repr(self.factcel))
            ret.append(',\n    neg : ')
            ret.append(repr(self.neg))
            ret.append(',\n    lneg : ')
            ret.append(repr(self.lneg))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def initialize_relax():
        """
        initialize_relax()
        
        
        Defined at Relax_module.fpp lines 36-80
        
        
        """
        _AresMainPy.f90wrap_initialize_relax()
    
    @staticmethod
    def destroy_relax():
        """
        destroy_relax()
        
        
        Defined at Relax_module.fpp lines 83-107
        
        
        """
        _AresMainPy.f90wrap_destroy_relax()
    
    @staticmethod
    def relaxer(nstep):
        """
        relaxer(nstep)
        
        
        Defined at Relax_module.fpp lines 110-136
        
        Parameters
        ----------
        nstep : int
        
        """
        _AresMainPy.f90wrap_relaxer(nstep=nstep)
    
    @staticmethod
    def fire_relax(lat, pos, fion, fcel):
        """
        fire_relax(lat, pos, fion, fcel)
        
        
        Defined at Relax_module.fpp lines 139-261
        
        Parameters
        ----------
        lat : float array
        pos : float array
        fion : float array
        fcel : float array
        
        ========================================================
        """
        _AresMainPy.f90wrap_fire_relax(lat=lat, pos=pos, fion=fion, fcel=fcel)
    
    @staticmethod
    def force_on_cell(stress, lat, fcel):
        """
        force_on_cell(stress, lat, fcel)
        
        
        Defined at Relax_module.fpp lines 264-281
        
        Parameters
        ----------
        stress : float array
        lat : float array
        fcel : float array
        
        """
        _AresMainPy.f90wrap_force_on_cell(stress=stress, lat=lat, fcel=fcel)
    
    @staticmethod
    def fix_direction(stress):
        """
        fix_direction(stress)
        
        
        Defined at Relax_module.fpp lines 284-304
        
        Parameters
        ----------
        stress : float array
        
        """
        _AresMainPy.f90wrap_fix_direction(stress=stress)
    
    @property
    def pstress(self):
        """
        Element pstress ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 14
        
        """
        return _AresMainPy.f90wrap_relax_module__get__pstress()
    
    @pstress.setter
    def pstress(self, pstress):
        _AresMainPy.f90wrap_relax_module__set__pstress(pstress)
    
    @property
    def ldone(self):
        """
        Element ldone ftype=logical pytype=bool
        
        
        Defined at Relax_module.fpp line 16
        
        """
        return _AresMainPy.f90wrap_relax_module__get__ldone()
    
    @ldone.setter
    def ldone(self, ldone):
        _AresMainPy.f90wrap_relax_module__set__ldone(ldone)
    
    @property
    def lfirst(self):
        """
        Element lfirst ftype=logical pytype=bool
        
        
        Defined at Relax_module.fpp line 16
        
        """
        return _AresMainPy.f90wrap_relax_module__get__lfirst()
    
    @lfirst.setter
    def lfirst(self, lfirst):
        _AresMainPy.f90wrap_relax_module__set__lfirst(lfirst)
    
    def __str__(self):
        ret = ['<relax_module>{\n']
        ret.append('    pstress : ')
        ret.append(repr(self.pstress))
        ret.append(',\n    ldone : ')
        ret.append(repr(self.ldone))
        ret.append(',\n    lfirst : ')
        ret.append(repr(self.lfirst))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

relax_module = Relax_Module()

class End_Module(f90wrap.runtime.FortranModule):
    """
    Module end_module
    
    
    Defined at End_module.fpp lines 5-29
    
    """
    @staticmethod
    def destroy_beast():
        """
        destroy_beast()
        
        
        Defined at End_module.fpp lines 13-28
        
        
        """
        _AresMainPy.f90wrap_destroy_beast()
    
    _dt_array_initialisers = []
    

end_module = End_Module()

class Band_Structure(f90wrap.runtime.FortranModule):
    """
    Module band_structure
    
    
    Defined at Bands_module.fpp lines 10-405
    
    """
    @staticmethod
    def init_bandstruct(numk, kvec):
        """
        init_bandstruct(numk, kvec)
        
        
        Defined at Bands_module.fpp lines 29-49
        
        Parameters
        ----------
        numk : int
        kvec : float array
        
        """
        _AresMainPy.f90wrap_init_bandstruct(numk=numk, kvec=kvec)
    
    @staticmethod
    def band_begin(numk, kvec):
        """
        band_begin(numk, kvec)
        
        
        Defined at Bands_module.fpp lines 52-107
        
        Parameters
        ----------
        numk : int
        kvec : float array
        
        """
        _AresMainPy.f90wrap_band_begin(numk=numk, kvec=kvec)
    
    @staticmethod
    def build_bands():
        """
        build_bands()
        
        
        Defined at Bands_module.fpp lines 110-171
        
        
        """
        _AresMainPy.f90wrap_build_bands()
    
    @staticmethod
    def read_bandpath(infile):
        """
        read_bandpath(infile)
        
        
        Defined at Bands_module.fpp lines 174-268
        
        Parameters
        ----------
        infile : str
        
        -------start input.dat-----------
        """
        _AresMainPy.f90wrap_read_bandpath(infile=infile)
    
    @staticmethod
    def band_pathread(infile):
        """
        band_pathread(infile)
        
        
        Defined at Bands_module.fpp lines 271-304
        
        Parameters
        ----------
        infile : str
        
        """
        _AresMainPy.f90wrap_band_pathread(infile=infile)
    
    @staticmethod
    def read_density(infile, rho):
        """
        read_density(infile, rho)
        
        
        Defined at Bands_module.fpp lines 307-351
        
        Parameters
        ----------
        infile : str
        rho : float array
        
        """
        _AresMainPy.f90wrap_read_density(infile=infile, rho=rho)
    
    @staticmethod
    def cal_band(rhos, nev, eigval):
        """
        cal_band(rhos, nev, eigval)
        
        
        Defined at Bands_module.fpp lines 354-404
        
        Parameters
        ----------
        rhos : float array
        nev : int
        eigval : float array
        
        """
        _AresMainPy.f90wrap_cal_band(rhos=rhos, nev=nev, eigval=eigval)
    
    _dt_array_initialisers = []
    

band_structure = Band_Structure()

class Aresmainapi(f90wrap.runtime.FortranModule):
    """
    Module aresmainapi
    
    
    Defined at AresMainAPI.fpp lines 5-104
    
    """
    @f90wrap.runtime.register_class("AresMainPy.aresOut")
    class aresOut(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=aresout)
        
        
        Defined at AresMainAPI.fpp lines 16-27
        
        """
        def __init__(self, handle=None):
            """
            self = Aresout()
            
            
            Defined at AresMainAPI.fpp lines 16-27
            
            
            Returns
            -------
            this : Aresout
            	Object to be constructed
            
            
            Automatically generated constructor for aresout
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _AresMainPy.f90wrap_aresout_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Aresout
            
            
            Defined at AresMainAPI.fpp lines 16-27
            
            Parameters
            ----------
            this : Aresout
            	Object to be destructed
            
            
            Automatically generated destructor for aresout
            """
            if self._alloc:
                _AresMainPy.f90wrap_aresout_finalise(this=self._handle)
        
        @property
        def forces(self):
            """
            Element forces ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__forces(self._handle)
            if array_handle in self._arrays:
                forces = self._arrays[array_handle]
            else:
                forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__forces)
                self._arrays[array_handle] = forces
            return forces
        
        @forces.setter
        def forces(self, forces):
            self.forces[...] = forces
        
        @property
        def poscar(self):
            """
            Element poscar ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__poscar(self._handle)
            if array_handle in self._arrays:
                poscar = self._arrays[array_handle]
            else:
                poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__poscar)
                self._arrays[array_handle] = poscar
            return poscar
        
        @poscar.setter
        def poscar(self, poscar):
            self.poscar[...] = poscar
        
        @property
        def pos(self):
            """
            Element pos ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__pos(self._handle)
            if array_handle in self._arrays:
                pos = self._arrays[array_handle]
            else:
                pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__pos)
                self._arrays[array_handle] = pos
            return pos
        
        @pos.setter
        def pos(self, pos):
            self.pos[...] = pos
        
        @property
        def chargerho(self):
            """
            Element chargerho ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__chargerho(self._handle)
            if array_handle in self._arrays:
                chargerho = self._arrays[array_handle]
            else:
                chargerho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__chargerho)
                self._arrays[array_handle] = chargerho
            return chargerho
        
        @chargerho.setter
        def chargerho(self, chargerho):
            self.chargerho[...] = chargerho
        
        @property
        def stress(self):
            """
            Element stress ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__stress(self._handle)
            if array_handle in self._arrays:
                stress = self._arrays[array_handle]
            else:
                stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__stress)
                self._arrays[array_handle] = stress
            return stress
        
        @stress.setter
        def stress(self, stress):
            self.stress[...] = stress
        
        @property
        def apilat_mat(self):
            """
            Element apilat_mat ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 22
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__apilat_mat(self._handle)
            if array_handle in self._arrays:
                apilat_mat = self._arrays[array_handle]
            else:
                apilat_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__apilat_mat)
                self._arrays[array_handle] = apilat_mat
            return apilat_mat
        
        @apilat_mat.setter
        def apilat_mat(self, apilat_mat):
            self.apilat_mat[...] = apilat_mat
        
        @property
        def apilat_para(self):
            """
            Element apilat_para ftype=real(dp) pytype=float
            
            
            Defined at AresMainAPI.fpp line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _AresMainPy.f90wrap_aresout__array__apilat_para(self._handle)
            if array_handle in self._arrays:
                apilat_para = self._arrays[array_handle]
            else:
                apilat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _AresMainPy.f90wrap_aresout__array__apilat_para)
                self._arrays[array_handle] = apilat_para
            return apilat_para
        
        @apilat_para.setter
        def apilat_para(self, apilat_para):
            self.apilat_para[...] = apilat_para
        
        @property
        def comm(self):
            """
            Element comm ftype=integer(i4b) pytype=int
            
            
            Defined at AresMainAPI.fpp line 24
            
            """
            return _AresMainPy.f90wrap_aresout__get__comm(self._handle)
        
        @comm.setter
        def comm(self, comm):
            _AresMainPy.f90wrap_aresout__set__comm(self._handle, comm)
        
        @property
        def myid(self):
            """
            Element myid ftype=integer(i4b) pytype=int
            
            
            Defined at AresMainAPI.fpp line 25
            
            """
            return _AresMainPy.f90wrap_aresout__get__myid(self._handle)
        
        @myid.setter
        def myid(self, myid):
            _AresMainPy.f90wrap_aresout__set__myid(self._handle, myid)
        
        @property
        def numprocs(self):
            """
            Element numprocs ftype=integer(i4b) pytype=int
            
            
            Defined at AresMainAPI.fpp line 26
            
            """
            return _AresMainPy.f90wrap_aresout__get__numprocs(self._handle)
        
        @numprocs.setter
        def numprocs(self, numprocs):
            _AresMainPy.f90wrap_aresout__set__numprocs(self._handle, numprocs)
        
        @property
        def rootid(self):
            """
            Element rootid ftype=integer(i4b) pytype=int
            
            
            Defined at AresMainAPI.fpp line 27
            
            """
            return _AresMainPy.f90wrap_aresout__get__rootid(self._handle)
        
        @rootid.setter
        def rootid(self, rootid):
            _AresMainPy.f90wrap_aresout__set__rootid(self._handle, rootid)
        
        def __str__(self):
            ret = ['<aresout>{\n']
            ret.append('    forces : ')
            ret.append(repr(self.forces))
            ret.append(',\n    poscar : ')
            ret.append(repr(self.poscar))
            ret.append(',\n    pos : ')
            ret.append(repr(self.pos))
            ret.append(',\n    chargerho : ')
            ret.append(repr(self.chargerho))
            ret.append(',\n    stress : ')
            ret.append(repr(self.stress))
            ret.append(',\n    apilat_mat : ')
            ret.append(repr(self.apilat_mat))
            ret.append(',\n    apilat_para : ')
            ret.append(repr(self.apilat_para))
            ret.append(',\n    comm : ')
            ret.append(repr(self.comm))
            ret.append(',\n    myid : ')
            ret.append(repr(self.myid))
            ret.append(',\n    numprocs : ')
            ret.append(repr(self.numprocs))
            ret.append(',\n    rootid : ')
            ret.append(repr(self.rootid))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init_alloc_arrays(self, nnatom):
        """
        init_alloc_arrays(self, nnatom)
        
        
        Defined at AresMainAPI.fpp lines 30-40
        
        Parameters
        ----------
        dertype : Aresout
        nnatom : int
        
        """
        _AresMainPy.f90wrap_init_alloc_arrays(dertype=self._handle, nnatom=nnatom)
    
    @staticmethod
    def assignment(self):
        """
        assignment(self)
        
        
        Defined at AresMainAPI.fpp lines 43-52
        
        Parameters
        ----------
        dertype : Aresout
        
        """
        _AresMainPy.f90wrap_assignment(dertype=self._handle)
    
    @staticmethod
    def destroy_alloc_arrays(self):
        """
        destroy_alloc_arrays(self)
        
        
        Defined at AresMainAPI.fpp lines 54-58
        
        Parameters
        ----------
        dertype : Aresout
        
        """
        _AresMainPy.f90wrap_destroy_alloc_arrays(dertype=self._handle)
    
    @staticmethod
    def updateions(pos, lattice=None):
        """
        updateions(pos[, lattice])
        
        
        Defined at AresMainAPI.fpp lines 61-104
        
        Parameters
        ----------
        pos : float array
        lattice : float array
        
        """
        _AresMainPy.f90wrap_updateions(pos=pos, lattice=lattice)
    
    _dt_array_initialisers = []
    

aresmainapi = Aresmainapi()

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
    _AresMainPy.f90wrap_zbrent(iu=iu, lreset=lreset, ebreak=ebreak, x=x, y=y, f=f, \
        xnew=xnew, xnewh=xnewh, ynew=ynew, yd=yd, ifail=ifail)

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
    _AresMainPy.f90wrap_kardir(nmax=nmax, v=v, basis=basis)

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
    _AresMainPy.f90wrap_ioncgr(iflag=iflag, nions=nions, toten=toten, a=a, b=b, \
        nfree=nfree, posion=posion, posioc=posioc, fact=fact, f=f, factsi=factsi, \
        fsif=fsif, fl=fl, s=s, dismax=dismax, iu6=iu6, iu0=iu0, ebreak=ebreak, \
        ediffg=ediffg, e1test=e1test, lstop2=lstop2)

