# Link libs that Ares is dependent on

## f2py

### 测试链接自己写的静态库

```fortran
subroutine matrixMult(A, B, C, n1, n2, n3, iCase)

    INTEGER, intent(in) :: n1, n2, n3
    INTEGER, intent(in) :: iCase
    DOUBLE PRECISION, intent(in) :: A(n1,n2)
    DOUBLE PRECISION, intent(in) :: B(n2,n3)
    DOUBLE PRECISION, intent(out) :: C(n1,n3)
    SELECT CASE (iCase)
        CASE (1)
            CALL matMult_loop1(A, B, C, n1, n2, n3)
        CASE (2)
            CALL matMult_loop2(A, B, C, n1, n2, n3)
        CASE (3)
            CALL matMult_func(A, B, C, n1, n2, n3)
    END SELECT
    RETURN

end subroutine matrixMult
```

matMult_loop1,matMult_loop1,matMult_loop1 三个做矩阵乘法的方法，已经用下面的方法编译为libmatMultiplication.a

```bash
ifort -c -fPIC *.f90
ar rc libmatMultiplication.a *.o
```

用下面的方式包

```bash
f2py --f90exec=ifort -L../testlib -lmatMultiplication -c matrixMult.f90 -m matrixMult
```

运行正常

### 测试链接APRACK写的静态库

```fortran
subroutine matrixMult(n1, A, D)

    INTEGER, intent(in) :: n1
    DOUBLE PRECISION, intent(in) :: A(n1,n1)
    DOUBLE PRECISION, intent(out) :: D(n1)
    CALL EVLRG(A,D)
    RETURN

end subroutine matrixMult
```

现在我添加了 EVLRG LAPACK 的函数, 想测试一下能不能链接上，运行了下面的命令

```bash
f2py --f90exec=ifort -L/work/home/xinyu/soft/ARES/ares_lib/lib -larpack_ifort -lxcf90 -lxc -c my.f90 -m my
```

f2py 正常运行结束

测试一下

```python
import my
ImportError: /work/home/xinyu/workplace/PhdProgram/aresPy/fwrapTest/LinkAresLib/my.cpython-38-x86_64-linux-gnu.so: undefined symbol: evlrg_
```

发现还是没找到这个函数
