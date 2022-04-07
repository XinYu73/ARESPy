# Question

现在可以动态库已经可以生成了,arespy也可以导入了

```python
>>> import arespy
>>> 
```

我模仿ares的主程序main.f90写了如下python程序

```python
import sys

sys.path.insert(0, '/work/home/xinyu/workplace/PhdProgram/aresPy/srcWrap/preprocessed')

import arespy
from mpi4py import MPI
comm = MPI.COMM_WORLD
arespy.Smpi_Math_Module.smpi_init()
rank = comm.Get_rank()
print("comm:    ",comm)
print("rank:    ",rank)
arespy.Smpi_Math_Module.start_time('[Total Time]',True)
arespy.Smpi_Math_Module.start_time('[SCF Time]',True)
arespy.Read_Module.read_file('ares.in')
arespy.Scalapack_Module.init_scala()
arespy.Begin_Module.initial_grid_pbc()
arespy.Potential_Module.vlpp()
arespy.Scf_Module.electronicscf()
# arespy.Smpi_Math_Module.end_time('[SCF Time]',True)
# arespy.Smpi_Math_Module.end_time('[Total Time]',True)
# arespy.Smpi_Math_Module.write_time('[Total Time]',True)
a = arespy.Constants()

print(a.dp)
```

结果如下，出现了Segmentation fault

```bash
comm:     <mpi4py.MPI.Intracomm object at 0x2b1ac07e4ef0>
rank:     0
Task name >>> OPTIMUS_PRIME                 
[The Atomic Spherical Cutting Radius is (Angs)]  -4.00000000000000     
Ecut=   940.075301025635      eV
The order of finite difference is:           8
Spin unpolarized
Max simulate steps is:           0
Max simulate steps is:           0
Chebyshev filter is used,order:          12
first Chebyshev filter order:         -24
Chebyshev:first step diag free
RayleighRitz needn''t OrthNorm step
Initialize the subspace by Pseudo orbitals
initial density by program
[The Adding Charge]  0.000000000000000E+000
Use Simple Mixing + (r)Pulay Mixing + kerker
parallel dims           1           1
[Reading PseudoPP FILE]C.pbe-mt-cpi.UPF              
[Reading PseudoPP FILE]H.pbe-mt-cpi.UPF              
[Ion and Charge Num.]           8   8.00000000000000     
[Eigen States Num.]          19
XC: PBE-GGA (Pewdew-Burke-Ernzerhof 1996)
communicate all grids
Total Charge# in Sphere is   7.99999999999999     
Non-Local pseudopotential has been used
R-space-GRIDS:          50          50          50
K-space-GRIDS:           1           1           1
Num of K-used:           1
>----SCF Iterations for Solving KS Equations----<
CheFSI: Dim. of Pseudo Subspace:           8
CheFSI: # of random samping states:          11
Segmentation fault
```

问题可能产生的位置:

1. 源代码preprocess出错
    1. 用preprocessed代码重新编译ares fortran程序
    2. 用mpirun -n 2测试, 结果显示没问题
    3. 说明preprocessed的代码没问题
2. f90wrap的代码有问题 ，由于f90wrap直接生成的f90wrap*.f90文件编译出错，我直接改了一部分f90wrap*.f90的代码，针对的问题主要是duplicate define,和数组定义

    数组定义问题：
    f90wrap 对下面的dfp没有给dimension，源程序给的dimension属性是dimension(:,:), 我的处理方式是手工给f90wrap文件加了dimension

    ```fortran
    subroutine f90wrap_om1c(nam, nuh, sp, dfp, voma, n1, n2, n3, n4)
        use mixer_module, only: om1c
        implicit none

        integer, intent(in) :: nam
        integer, intent(in) :: nuh
        real(8) :: sp
        real(8), dimension(n1, n2) :: dfp
        real(8), dimension(n3, n4) :: voma
        integer :: n1
        !f2py intent(hide), depend(dfp) :: n1 = shape(dfp,0)
        integer :: n2
        !f2py intent(hide), depend(dfp) :: n2 = shape(dfp,1)
        integer :: n3
        !f2py intent(hide), depend(voma) :: n3 = shape(voma,0)
        integer :: n4
        !f2py intent(hide), depend(voma) :: n4 = shape(voma,1)
        call om1c(NAM=nam, NUH=nuh, sp=sp, dfp=dfp, voma=voma)
    end subroutine f90wrap_om1c
    ```

    duplicate define 问题:
    由于源程序中部分的数组的index不是从1开始，f90wrap会把index单独赋予变量，下面是一个例子，参数列表中n11, n22, n33,是我手工改的，原来是n1,n2,n3,造成了duplicate define的问题

    ```fortran
    subroutine f90wrap_rcsntable(n1, n2, n3, dn, h, srmax, rcs, cns, sphindx, nspt, n0, n11, n22, n33, n4, n5, n6, n7, n8)
        use math, only: rcsntable
        implicit none

        integer(4), intent(in) :: n1
        integer(4), intent(in) :: n2
        integer(4), intent(in) :: n3
        integer(4), intent(in) :: dn
        real(8), intent(in) :: h
        real(8), intent(in) :: srmax
        real(8), intent(inout), dimension(3, n0, n1, n2) :: rcs
        complex(8), intent(inout), dimension(n3, n4, n5) :: cns
        logical, intent(inout), dimension(n6, n7, n8) :: sphindx
        integer(4), intent(out) :: nspt
        integer :: n0
        !f2py intent(hide), depend(rcs) :: n0 = shape(rcs,1)
        integer :: n11
        !f2py intent(hide), depend(rcs) :: n11 = shape(rcs,2)
        integer :: n22
        !f2py intent(hide), depend(rcs) :: n22 = shape(rcs,3)
        integer :: n33
        !f2py intent(hide), depend(cns) :: n33 = shape(cns,0)
        integer :: n4
        !f2py intent(hide), depend(cns) :: n4 = shape(cns,1)
        integer :: n5
        !f2py intent(hide), depend(cns) :: n5 = shape(cns,2)
        integer :: n6
        !f2py intent(hide), depend(sphindx) :: n6 = shape(sphindx,0)
        integer :: n7
        !f2py intent(hide), depend(sphindx) :: n7 = shape(sphindx,1)
        integer :: n8
        !f2py intent(hide), depend(sphindx) :: n8 = shape(sphindx,2)
        call rcsntable(n1=n1, n2=n2, n3=n3, dn=dn, h=h, srmax=srmax, rcs=rcs, cns=cns, Sphindx=sphindx, nspt=nspt)
    end subroutine f90wrap_rcsntable
    ```

现在问题就是我不知道改怎么定位这个segmentation fault的问题,源代码没问题,可能还是f90wrap*.f90的问题，不清楚该怎么处理

## 下面的和问题无关

```fortran
subroutine f90wrap_amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    use mixer_module, only: amst
    implicit none

    real(8) :: beta
    real(8) :: w0
    integer, intent(in) :: nam
    integer, intent(in) :: nuh
    real(8), dimension(n1, n2) :: dxp
    real(8), dimension(n3, n4) :: dfp
    real(8) :: sp
    real(8), dimension(n5) :: xl
    real(8), dimension(n6) :: fl
    real(8), dimension(n7, n8) :: voma
    real(8), dimension(n9) :: xn
    integer :: n1
    !f2py intent(hide), depend(dxp) :: n1 = shape(dxp,0)
    integer :: n2
    !f2py intent(hide), depend(dxp) :: n2 = shape(dxp,1)
    integer :: n3
    !f2py intent(hide), depend(dfp) :: n3 = shape(dfp,0)
    integer :: n4
    !f2py intent(hide), depend(dfp) :: n4 = shape(dfp,1)
    integer :: n5
    !f2py intent(hide), depend(xl) :: n5 = shape(xl,0)
    integer :: n6
    !f2py intent(hide), depend(fl) :: n6 = shape(fl,0)
    integer :: n7
    !f2py intent(hide), depend(voma) :: n7 = shape(voma,0)
    integer :: n8
    !f2py intent(hide), depend(voma) :: n8 = shape(voma,1)
    integer :: n9
    !f2py intent(hide), depend(xn) :: n9 = shape(xn,0)
    call amst(beta=beta, w0=w0, NAM=nam, NUH=nuh, dxp=dxp, dfp=dfp, sp=sp, xl=xl, fl=fl, voma=voma, xn=xn)
end subroutine f90wrap_amst
```

```fortran
subroutine f90wrap_rcsntable_atoms(n1, n2, n3, dn, na, h, poscar, srmax, atomr, rvec, rcs, cns, sphindx, nspt, n0, n111, &
                                   n22, n33, n4, n5, n6, n7, n8, n9, n10, n11, n12)
    use math, only: rcsntable_atoms
    implicit none

    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: dn
    integer(4), intent(in) :: na
    real(8), intent(in) :: h
    real(8), intent(in), dimension(3, n0) :: poscar
    real(8), intent(in) :: srmax
    real(8), intent(in) :: atomr
    real(8), intent(in), dimension(4, n1, n2, n3) :: rvec
    real(8), intent(inout), dimension(3, n4, n5, n6) :: rcs
    complex(8), intent(inout), dimension(n7, n8, n9) :: cns
    logical, intent(inout), dimension(n10, n11, n12) :: sphindx
    integer(4), intent(out) :: nspt
    integer :: n0
    !f2py intent(hide), depend(poscar) :: n0 = shape(poscar,1)
    integer :: n111
    !f2py intent(hide), depend(rvec) :: n111 = shape(rvec,1)
    integer :: n22
    !f2py intent(hide), depend(rvec) :: n22 = shape(rvec,2)
    integer :: n33
    !f2py intent(hide), depend(rvec) :: n33 = shape(rvec,3)
    integer :: n4
    !f2py intent(hide), depend(rcs) :: n4 = shape(rcs,1)
    integer :: n5
    !f2py intent(hide), depend(rcs) :: n5 = shape(rcs,2)
    integer :: n6
    !f2py intent(hide), depend(rcs) :: n6 = shape(rcs,3)
    integer :: n7
    !f2py intent(hide), depend(cns) :: n7 = shape(cns,0)
    integer :: n8
    !f2py intent(hide), depend(cns) :: n8 = shape(cns,1)
    integer :: n9
    !f2py intent(hide), depend(cns) :: n9 = shape(cns,2)
    integer :: n10
    !f2py intent(hide), depend(sphindx) :: n10 = shape(sphindx,0)
    integer :: n11
    !f2py intent(hide), depend(sphindx) :: n11 = shape(sphindx,1)
    integer :: n12
    !f2py intent(hide), depend(sphindx) :: n12 = shape(sphindx,2)
    call rcsntable_atoms(n1=n1, n2=n2, n3=n3, dn=dn, na=na, h=h, poscar=poscar, srmax=srmax, atomR=atomr, rVec=rvec, &
                         rcs=rcs, cns=cns, Sphindx=sphindx, nspt=nspt)
end subroutine f90wrap_rcsntable_atoms
```

```fortran
subroutine f90wrap_grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs, gridrange_sum, n1, n2, n3, n, n0, &
                              n11, n22)
    use smpi_math_module, only: grid_split
    implicit none

    integer(4), intent(in) :: ngrid
    integer(4), intent(in) :: ncore
    integer(4), intent(in) :: comm
    integer(4), intent(in) :: id
    integer(4), dimension(3), intent(inout) :: grid_range
    integer(4), intent(inout), dimension(n0) :: recvcounts
    integer(4), intent(inout), dimension(n1) :: displs
    integer(4), optional, intent(inout), dimension(3, n2) :: gridrange_sum
    integer(4), intent(inout), optional :: n1
    integer(4), intent(inout), optional :: n2
    integer(4), intent(inout), optional :: n3
    integer(4), intent(inout), optional :: n
    integer :: n0
    !f2py intent(hide), depend(recvcounts) :: n0 = shape(recvcounts,0)
    integer :: n11
    !f2py intent(hide), depend(displs) :: n11 = shape(displs,0)
    integer :: n22
    !f2py intent(hide), depend(gridrange_sum) :: n2 = shape(gridrange_sum,1)
    call grid_split(ngrid=ngrid, ncore=ncore, comm=comm, id=id, grid_range=grid_range, recvcounts=recvcounts, displs=displs, &
                    gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)
end subroutine f90wrap_grid_split
```

```fortran
subroutine f90wrap_grid_split_sph(ngrid, ncore, comm, id, grid_range, recvcounts, displs, gridrange_sum, n1, n2, n3, n, &
                                  n0, n11, n22)
    use grid_module, only: grid_split_sph
    implicit none

    integer(4), intent(in) :: ngrid
    integer(4), intent(in) :: ncore
    integer(4), intent(in) :: comm
    integer(4), intent(in) :: id
    integer(4), dimension(3), intent(inout) :: grid_range
    integer(4), intent(inout), dimension(n0) :: recvcounts
    integer(4), intent(inout), dimension(n1) :: displs
    integer(4), intent(inout), dimension(3, n2) :: gridrange_sum
    integer(4), intent(inout) :: n1
    integer(4), intent(inout) :: n2
    integer(4), intent(inout) :: n3
    integer(4), intent(inout) :: n
    integer :: n0
    !f2py intent(hide), depend(recvcounts) :: n0 = shape(recvcounts,0)
    integer :: n11
    !f2py intent(hide), depend(displs) :: n11 = shape(displs,0)
    integer :: n22
    !f2py intent(hide), depend(gridrange_sum) :: n22 = shape(gridrange_sum,1)
    call grid_split_sph(ngrid=ngrid, ncore=ncore, comm=comm, id=id, grid_range=grid_range, recvcounts=recvcounts, &
                        displs=displs, gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)
end subroutine f90wrap_grid_split_sph

subroutine f90wrap_sphere_region(n1, n2, n3, lsphere, nz_map, nsphere, n0, n11, n22, n33)
    use grid_module, only: sphere_region
    implicit none

    integer(4) :: n1
    integer(4) :: n2
    integer(4) :: n3
    logical, dimension(n0, n1, n2) :: lsphere
    integer(4), dimension(2, n3) :: nz_map
    integer(4) :: nsphere
    integer :: n0
    !f2py intent(hide), depend(lsphere) :: n0 = shape(lsphere,0)
    integer :: n11
    !f2py intent(hide), depend(lsphere) :: n11 = shape(lsphere,1)
    integer :: n22
    !f2py intent(hide), depend(lsphere) :: n22 = shape(lsphere,2)
    integer :: n33
    !f2py intent(hide), depend(nz_map) :: n33 = shape(nz_map,1)
    call sphere_region(n1=n1, n2=n2, n3=n3, Lsphere=lsphere, nz_map=nz_map, nsphere=nsphere)
end subroutine f90wrap_sphere_region
```

```fortran
subroutine f90wrap_cal_trans_phase(nr, nspin, r_new, n1, n2, n3, ng1, ng2, ng3, gvec, trans_phase, n0, n11, n22, n33, n4, &
                                   n5)
    use succeed, only: cal_trans_phase
    implicit none

    integer(4), intent(in) :: nr
    integer(4), intent(in) :: nspin
    real(8), intent(in), dimension(3, n0) :: r_new
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: ng1
    integer(4), intent(in) :: ng2
    integer(4), intent(in) :: ng3
    real(8), intent(in), dimension(4, n1) :: gvec
    complex(8), intent(inout), dimension(n2, n3, n4, n5) :: trans_phase
    integer :: n0
    !f2py intent(hide), depend(r_new) :: n0 = shape(r_new,1)
    integer :: n11
    !f2py intent(hide), depend(gvec) :: n11 = shape(gvec,1)
    integer :: n22
    !f2py intent(hide), depend(trans_phase) :: n22 = shape(trans_phase,0)
    integer :: n33
    !f2py intent(hide), depend(trans_phase) :: n33 = shape(trans_phase,1)
    integer :: n4
    !f2py intent(hide), depend(trans_phase) :: n4 = shape(trans_phase,2)
    integer :: n5
    !f2py intent(hide), depend(trans_phase) :: n5 = shape(trans_phase,3)
    call cal_trans_phase(nr=nr, nspin=nspin, r_new=r_new, n1=n1, n2=n2, n3=n3, ng1=ng1, ng2=ng2, ng3=ng3, gvec=gvec, &
                         trans_phase=trans_phase)
end subroutine f90wrap_cal_trans_phase
```

```fortran
subroutine f90wrap_get_new_rho_psi(nr, r_new, nrho, n1, n2, n3, nspin, rho_new, n_s, psi_new, gvec, n0, n11, n22, n33, n44, &
                                   n5, n6)
    use succeed, only: get_new_rho_psi
    implicit none

    integer(4), intent(in) :: nr
    real(8), intent(in), dimension(3, n0) :: r_new
    integer(4), intent(in) :: nrho
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: nspin
    real(8), intent(inout), dimension(n1, n2) :: rho_new
    integer(4), intent(in) :: n_s
    real(8), intent(inout), dimension(n3, n44, n5) :: psi_new
    real(8), intent(in), dimension(4, n6) :: gvec
    integer :: n0
    !f2py intent(hide), depend(r_new) :: n0 = shape(r_new,1)
    integer :: n11
    !f2py intent(hide), depend(rho_new) :: n11 = shape(rho_new,0)
    integer :: n22
    !f2py intent(hide), depend(rho_new) :: n22 = shape(rho_new,1)
    integer :: n33
    !f2py intent(hide), depend(psi_new) :: n33 = shape(psi_new,0)
    integer :: n44
    !f2py intent(hide), depend(psi_new) :: n44 = shape(psi_new,1)
    integer :: n5
    !f2py intent(hide), depend(psi_new) :: n5 = shape(psi_new,2)
    integer :: n6
    !f2py intent(hide), depend(gvec) :: n6 = shape(gvec,1)
    call get_new_rho_psi(nr=nr, r_new=r_new, nrho=nrho, n1=n1, n2=n2, n3=n3, Nspin=nspin, rho_new=rho_new, n_s=n_s, &
                         psi_new=psi_new, gvec=gvec)
end subroutine f90wrap_get_new_rho_psi
```
