# Question

## 问题关于derived type中Allocatable Array

为了获取ares的结果,按照师兄的指点,重点关注energy，force，stress

1. energy 是一个全局变量可以轻易获取
2. force和stress在一个结构体中

    ```fortran
    TYPE struct_type
        REAL(DP),ALLOCATABLE      :: Zion (:)
        INTEGER(I4B),ALLOCATABLE  :: nati(:)  !number of atoms for each atom type
        INTEGER(I4B),ALLOCATABLE  :: eleid(:)   !atoms for each atom type id
        REAL(DP),ALLOCATABLE      :: pos(:,:)  !direct
        REAL(DP),ALLOCATABLE      :: poscar(:,:) ! car
        REAL(DP)                  :: stress(3,3)
        REAL(DP),ALLOCATABLE      :: forces(:,:)       ! total forces
        REAL(DP),ALLOCATABLE      :: mass(:)       ! atomic mass
        REAL(DP),ALLOCATABLE      :: zeta(:,:)    ! zeta for STO
        INTEGER(I4B),ALLOCATABLE  :: prinq(:,:)   !atom's principle quantum
        INTEGER(I4B),ALLOCATABLE  :: Lmax(:)   !max L need to consider
        CHARACTER(len=3),ALLOCATABLE :: elements(:)
    END TYPE struct_type
    ```

    stress 是一个简单的大小明确的数组，但是force和其他量都是Allocatable的数组,我尝试按照stress的方法获取得到

    forece 在arespy中

    ```python
        @property
        def forces(self):
            """
            Element forces ftype=real(dp) pytype=float
            
            
            Defined at Struct_module.fpp line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _arespy.f90wrap_struct_type__array__forces(self._handle)
            if array_handle in self._arrays:
                forces = self._arrays[array_handle]
            else:
                print("IIIIIIIIII")
                print(f90wrap.runtime.sizeof_fortran_t,self._handle,_arespy.f90wrap_struct_type__array__forces)
                print("IIIIIIIIII")
                forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,self._handle,_arespy.f90wrap_struct_type__array__forces)
                self._arrays[array_handle] = forces
            return forces
        
        @forces.setter
        def forces(self, forces):
            self.forces[...] = forces
    ```

    在python程序中

    ```python
    print("Force",aresstruct_module.struct_type().forces)
    ```

    得到下面错误

    ```bash
    IIIIIIIIII
    Traceback (most recent call last):
    File "WrappedAresTest.py", line 45, in <module>
        print("Force",aresstruct_module.struct_type().forces)
    File "/work/home/xinyu/workplace/PhdProgram/aresPy/srcWrap/preprocessed/arespy.py", line 6962, in forces
        forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,self._handle,_arespy.f90wrap_struct_type__array__forces)
        ValueError: array is NULL
    ```

    然后我猜测是不是struct中的数组内存被收回了，我就在ares程序中的自洽计算部分后面补充了一struct的打印输出

    ```fortran
        write (*, *) 'struct%forces', struct%forces
        write (*, *) 'struct%Zion', struct%Zion
        write (*, *) 'struct%nati', struct%nati
        write (*, *) 'struct%eleid', struct%eleid
        write (*, *) 'struct%mass', struct%mass
    ```

    结果发现

    ```bash
    struct%forces  0.000000000000000E+000  0.000000000000000E+000
    0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    0.000000000000000E+000
    struct%Zion   4.00000000000000        1.00000000000000     
    struct%nati           1           4
    struct%eleid           1           2           6
    struct%mass  0.000000000000000E+000  0.000000000000000E+000
    ```

    这说明我对于内存被收回的想法是错的

    说明要么是我获取值的方式有问题，要么是封装的问题

    现在我打算再写一个fortran module, 在自洽计算后吧struct的值用普通数组保留一下，但这显然有点笨

    不知道师兄有没有什么解决方案推荐
