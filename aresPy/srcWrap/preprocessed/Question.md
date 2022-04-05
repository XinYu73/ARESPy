# 今天的问题

## 麻烦师兄啦（！V !）

今天用下面的命令对ares代码中#ifdef MPI的部分做了处理

```bash
prepro:
    for item in $(ARESSRC);do $(MPIFC) -EP -P -DMPI $$item -cpp  -o ./preprocessed/$$item; done;
```

处理的结果是arespy.py中不在有重复定义的成分

接下来用继续用我以前写的makefile编译发现关于derived type的错误,然后我就去github找答案，发现example [derived type](https://github.com/jameskermode/f90wrap/tree/master/examples/derivedtypes),之后我仿照example中的makefile又写了一个makefile 命名为MakefileExample

执行后我遇到了下面的问题，我理解的是Ewald.f90的function未被识别,但是example中也是有function的，我就不太理解了

MakefileExample的思路应该是先把源文件建成一个静态库，然后f90wrap*.f90这些文件去访问

```bash
mpiifort:f90: f90wrap_Arpack_module.f90
mpiifort:f90: f90wrap_Begin_module.f90
mpiifort:f90: f90wrap_Chebyshev_fliter.f90
mpiifort:f90: f90wrap_Constants.f90
mpiifort:f90: f90wrap_Energy_module.f90
mpiifort:f90: f90wrap_Ewald.f90
f90wrap_Ewald.f90(195): error #6580: Name in only-list does not exist or is not accessible.   [EWALDRPSTR]
    use ewald, only: ewaldrpstr
---------------------^
f90wrap_Ewald.f90(201): error #6404: This name does not have a type, and must have an explicit type.   [EWALDRPSTR]
    ret_ewaldrpstr = ewaldrpstr(eta=eta)
---------------------^
f90wrap_Ewald.f90(205): error #6580: Name in only-list does not exist or is not accessible.   [EWALDAVSTR]
    use ewald, only: ewaldavstr
---------------------^
f90wrap_Ewald.f90(211): error #6404: This name does not have a type, and must have an explicit type.   [EWALDAVSTR]
    ret_ewaldavstr = ewaldavstr(eta=eta)
---------------------^
compilation aborted for f90wrap_Ewald.f90 (code 1)
error: Command "mpiifort -FR -fPIC -O3 -traceback -heap-arrays 1000 -align array32byte -align rec32byte -fpp -DMPI -I/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib/fftw/include -I/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib/include -module . -fPIC -fp-model strict -O1 -assume minus0 -qopenmp -I./src.linux-x86_64-3.8/./src.linux-x86_64-3.8 -I/work/home/xinyu/soft/XinYuEnv/lib/python3.8/site-packages/numpy/core/include -I/work/home/xinyu/soft/XinYuEnv/include/python3.8 -c -c f90wrap_Ewald.f90 -o ./f90wrap_Ewald.o" failed with exit status 1
make: *** [_arespy.so] Error 1
```
