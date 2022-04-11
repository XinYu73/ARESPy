# Question 关于derived type中Allocatable array的封装

按照师兄的指点，Fortran的derived type的allocatable array是不能直接访问的，需要写一个接口类似于师兄的[qepy的api](https://gitlab.com/shaoxc/qepy/-/tree/master/src/api)

但是我没太看懂哈哈,然后我在这个[链接](https://github.com/jameskermode/f90wrap/blob/master/examples/derivedtypes/datatypes.f90)发现

```fortran
type alloc_arrays
    REAL(idp), DIMENSION(:,:), ALLOCATABLE :: chi
    REAL(idp), DIMENSION(:,:), ALLOCATABLE :: psi
    INTEGER(4) chi_shape(2)
    INTEGER(4) psi_shape(2)
end type alloc_arrays
```

这个和ares中struct类似的derived type,这个例子中对alloc_arrays配置了下面的函数

```fortran
subroutine init_alloc_arrays(dertype, m, n)
    type(alloc_arrays), INTENT(inout) :: dertype
    INTEGER(4), INTENT(in) :: m, n
    allocate(dertype%chi(m,n))
    allocate(dertype%psi(m,n))
end subroutine init_alloc_arrays

subroutine destroy_alloc_arrays(dertype)
    type(alloc_arrays), INTENT(inout) :: dertype
    if (allocated(dertype%chi)) deallocate(dertype%chi)
    if (allocated(dertype%psi)) deallocate(dertype%psi)
end subroutine destroy_alloc_arrays
```

和ares中给struct分配内存的方式有所不同，然后我就尝试了一下这个方法,写了一个类似的aresAPI.f90

```fortran
module aresAPI
    USE constants
    USE Struct_module
    implicit none
    type aresOut
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: forces
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: stress
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: poscar
        REAL(DP), ALLOCATABLE, DIMENSION(:, :)  :: pos
        REAL(DP), ALLOCATABLE, DIMENSION(:, :, :, :) :: chargeRho
    end type aresOut
contains
    subroutine init_alloc_arrays(dertype, nnatom)
        USE parameters, ONLY: nspin
        USE grid_module, ONLY: ni1 => global_n1, ni2 => global_n2, ni3 => global_n3 &
                    & , Sphere2Cubic, sph => grid, nps => n
        implicit none
        type(aresout), INTENT(inout) :: dertype
        INTEGER(I4B), intent(in) :: nnatom
        REAL(DP) :: rhot(ni1, ni2, ni3, 2)
        !force
        ALLOCATE (dertype%forces(3, nnatom))
        dertype%forces = struct%forces
        !stress
        ALLOCATE (dertype%stress(3, 3))
        dertype%stress = struct%stress
        !poscar
        ALLOCATE (dertype%poscar(3, nnatom))
        dertype%poscar = struct%poscar
        !pos
        ALLOCATE (dertype%pos(3, nnatom))
        dertype%pos = struct%pos
        !electron density????????????????
        ALLOCATE (dertype%chargeRho(ni1, ni2, ni3, 2))
        IF (nspin == 1) THEN
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1), rhot(:, :, :, 1))
            rhot(:, :, :, 2) = 0._DP
        ELSE
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1) + sph%rhoS(:, 2), rhot(:, :, :, 1))
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1) - sph%rhoS(:, 2), rhot(:, :, :, 2))
        END IF
        dertype%chargeRho = rhot
        !
    end subroutine init_alloc_arrays

    subroutine destroy_alloc_arrays(dertype)
        type(aresout), INTENT(inout) :: dertype
        if (allocated(dertype%forces)) deallocate (dertype%forces)
        if (allocated(dertype%stress)) deallocate (dertype%stress)
        if (allocated(dertype%poscar)) deallocate (dertype%poscar)
        if (allocated(dertype%pos)) deallocate (dertype%pos)
        if (allocated(dertype%chargeRho)) deallocate (dertype%chargeRho)
    end subroutine destroy_alloc_arrays
end module aresAPI
```

在调用ares的python程序中我是这样调用的

```python
Out = ares.Aresapi.aresOut()
#calculate
ares.Scalapack_Module.init_scala()
ares.Begin_Module.initial_grid_pbc()
ares.Potential_Module.vlpp()
ares.Scf_Module.electronicscf()
#calculate done
ares.Aresapi.init_alloc_arrays(Out,5)

print("Out.forces\n",Out.forces)
print("Out.Stress\n",Out.stress)
print("Out.poscar\n",Out.poscar)
print("Out.pos\n",Out.pos)
print("Out.chargeRho\n",Out.chargerho[0,0,:,0])
```

这样的方法目前倒是也能调用,但是因为看不太懂师兄的api，不知这两种方法有何区别 ,这是一个问题

另外的话今天发现不管怎么调节input文件，ares的给的force和stress都是0，后来问了问发现给我的ares根本不能算force和stress,明天又要重新包另外一个版本的，估计到时候又会有很多问题，还得麻烦师兄哈
