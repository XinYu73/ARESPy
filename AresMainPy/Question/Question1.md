# Question

## 关于 Arespy 的并行

ares的并行初始化在函数内部

```fortran
  SUBROUTINE smpi_init()
    IMPLICIT NONE

    CALL MPI_INIT(mpinfo)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, parallel%myid, mpinfo)
    CALL MPI_Comm_SIZE(MPI_COMM_WORLD, parallel%numprocs, mpinfo)

    parallel%comm  = MPI_COMM_WORLD

    parallel%rootid = 0
    if(parallel%rootid == parallel%myid)then
       parallel%isroot = .true.
    else
       parallel%isroot = .false.
    endif
  END SUBROUTINE smpi_init
```

该函数并不需要外部参数，我们在python中调用smpi_init()应该就可以给ares初始化并行环境