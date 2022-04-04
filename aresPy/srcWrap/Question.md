# Question

ares程序中对mpi有下面的处理,举例:(potential_module)

```fortran
#ifndef MPI
   REAL(DP),ALLOCATABLE  :: V_accelerate(:,:,:)
#else
   REAL(DP),ALLOCATABLE  :: V_accelerate(:)
#endif
```

我用下面的命令生成 arespy.py

```bash
f90WrapFile:
    f90wrap -m arespy $(ARESSRC) -k kindMap.json
```

结果arespy.py中就有两个V_accelerate

```python
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
```

现在的问题就是我该怎么处理源文件中检测是不是有mpi的代码块

```fortran
#ifndef mpi

#endif
```
