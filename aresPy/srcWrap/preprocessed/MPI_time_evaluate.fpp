# 1 "MPI_time_evaluate.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "MPI_time_evaluate.f90"
MODULE m_time_evaluate
  USE constants
! USE parameters, ONLY: time_flag
!>> ------------------------------------------------
!> useage:
!  CALL time_start('info of the time to calculate')
!  CALL time_end('info of the time to calculate')
!  CALL time_output('info of the time to calculate')
!>> ------------------------------------------------

  IMPLICIT NONE
  TYPE time_record
     CHARACTER(len=100),allocatable :: str(:)
     INTEGER(I4B),allocatable       :: t1(:),t2(:)
     INTEGER(I4B)                   :: length
  ENDTYPE time_record
  TYPE(time_record) :: data
  REAL(DP) :: total_memory_consum=0.d0
  CHARACTER(len=20) :: filename='memory_sum'
  integer(I4B) :: file_unit=1315
CONTAINS

  SUBROUTINE record(str,t)
    IMPLICIT NONE
    CHARACTER(len=100),INTENT(IN) :: str
    INTEGER(I4B),INTENT(IN) :: t
    LOGICAL  :: add_new
    INTEGER(I4B) :: i

    add_new=.TRUE.

    IF(data%length/=0)THEN
       do i=1,data%length,1
          IF(str==data%str(i))THEN
             data%t1(i)=data%t2(i)
             data%t2(i)=t
             add_new=.FALSE.
          ENDIF
       enddo
    ELSE
       data%length=1
       allocate(data%str(100))
       allocate(data%t1(100))
       allocate(data%t2(100))
       data%str(1)=str
       data%t1(1)=data%t2(1)
       data%t2(1)=t
       add_new=.FALSE.
    ENDIF

    CALL reset_data(add_new,str,t)

  ENDSUBROUTINE record


  SUBROUTINE reset_data(myflag,str,t)
    IMPLICIT NONE
    LOGICAL,INTENT(IN) :: myflag
    CHARACTER(len=100) :: str
    INTEGER(I4B)      :: t
    type(time_record) :: data_local
    INTEGER(I4B)      :: size_new, size_old

    IF(myflag)THEN
       data%length = data%length + 1
       IF(data%length > size(data%str))THEN

          size_old=size(data%str)
          size_new=size_old+50

          allocate(data_local%str(size_old))
          allocate(data_local%t1(size_old))
          allocate(data_local%t2(size_old))

          data_local%str(:)=data%str(:)
          data_local%t1(:)=data%t1(:)
          data_local%t2(:)=data%t2(:)

          deallocate(data%str,data%t1,data%t2)
          allocate(data%str(size_new))
          allocate(data%t1(size_new))
          allocate(data%t2(size_new))

          data%str(1:size_old)=data_local%str
          data%t1(1:size_old)=data_local%t1
          data%t2(1:size_old)=data_local%t2
       ENDIF

       data%str(data%length)=str
       data%t1(data%length)=t
       data%t2(data%length)=t
    ENDIF
  ENDSUBROUTINE reset_data

  SUBROUTINE init_recrod()
    IMPLICIT NONE

! allocate(data_local%str(1))
! allocate(data_local%t1(si))
! allocate(data_local%t2(size_old))
    data%length=0
  ENDSUBROUTINE init_recrod

  SUBROUTINE time_start(str,flag)
    IMPLICIT NONE
    CHARACTER(len=*) :: str
    LOGICAL          :: flag
!> local
    CHARACTER(len=100) :: str1
    INTEGER(I4B)     :: t

    if(.not.flag)return
    str1=str
    CALL system_clock(t)
    CALL record(str1,t)
  ENDSUBROUTINE time_start

  SUBROUTINE time_end(str,flag)
    IMPLICIT NONE
    CHARACTER(len=*) :: str
    LOGICAL          :: flag
!> local
    CHARACTER(len=100) :: str1
    INTEGER(I4B)     :: t

    if(.not.flag)return
    str1=str
    CALL system_clock(t)
    CALL record(str1,t)
  ENDSUBROUTINE time_end

  SUBROUTINE time_output(str,flag)
    IMPLICIT NONE
    CHARACTER(len=*) :: str
    LOGICAL          :: flag
!> local
    CHARACTER(len=100) :: str1
    INTEGER(I4B)     :: i

    if(.not.flag)return
    str1=str
    do i=1,data%length,1
       IF(str1==data%str(i))THEN
          print *,str,' : ',(data%t2(i)-data%t1(i))/10000.d0
          return
       ENDIF
    enddo
    print*, 'time evalute output error\n',"can't find out the profile string"
  ENDSUBROUTINE time_output


  Subroutine memory_sum(label,memory_size)
    implicit none
    CHARACTER(len=*) :: label
    real(DP)      :: memory_size
!> local
    logical :: lexist
    real(DP)     :: size_GB
!>
return
    size_GB = real(memory_size,DP)/1000.d0/1000.d0/1000.d0
    total_memory_consum =total_memory_consum + size_GB

!> out
    INQUIRE(FILE=trim(adjustl(filename)),EXIST=lexist)
    if(.not. lexist)then
! print*,"mem_summary not exist, a new file is generate"
       open(file_unit,FILE=trim(adjustl(filename)))
       write(file_unit,'(A20,1X,F25.9,A2)',advance='no')&
            &trim(adjustl(label)),size_GB,'GB'
       write(file_unit,'(A20,1X,F25.9,A2)')&
            &"total_memory_consum",total_memory_consum,'GB'
       close(file_unit)
    else
       open(file_unit,FILE=trim(adjustl(filename)),POSITION='append')
       write(file_unit,'(A20,1X,F25.9,A2)',advance='no')&
            &trim(adjustl(label)),size_GB,'GB'
       write(file_unit,'(A20,1X,F25.9,A2)')&
            &"total_memory_consum",total_memory_consum,'GB'
       close(file_unit)
    endif

  end Subroutine memory_sum


  Subroutine memory_free(label,memory_size)
    implicit none
    CHARACTER(len=*) :: label
    real(DP)      :: memory_size
!> local
    logical :: lexist
    real(DP)     :: size_GB

!>
return
    size_GB = real(memory_size,DP)/1000.d0/1000.d0/1000.d0
    total_memory_consum =total_memory_consum - size_GB

!> out
    INQUIRE(FILE=trim(adjustl(filename)),EXIST=lexist)
    if(.not. lexist)then
       print *,"ERROR: can not find outfile when memory free stat"
    else
       open(file_unit,FILE=trim(adjustl(filename)),POSITION='append')
       write(file_unit,'(A20,1X,F25.9,A2)',advance='no')&
            & trim(adjustl(label)),-size_GB,'GB'
       write(file_unit,'(A20,1X,F25.9,A2)')&
            & "total_memory_consum",total_memory_consum,'GB'
       close(file_unit)
    endif

  end Subroutine memory_free


! Subroutine out_total_memory()
!   implicit none
!   !> local
!   logical :: lexist
!   integer(I4B) :: file_unit=1315

!   INQUIRE(FILE=trim(adjustl(filename)),EXIST=lexist)
!   if(.not. lexist)then
!      print*,"mem_summary not exist, a new file is generate"
!      open(file_unit,FILE=trim(adjustl(filename)))
!      write(file_unit,'(A20,1X,F25.9,A2)')&
!           &trim(adjustl(label)),total_memory_consum,'GB'
!      close(file_unit)
!   else
!      open(file_unit,FILE=trim(adjustl(filename)),POSITION='append')
!      write(file_unit,'(A20,1X,F25.9,A2)')&
!           &trim(adjustl(label)),total_memory_consum,'GB'
!      close(file_unit)
!   endif
! end Subroutine out_total_memory


ENDMODULE m_time_evaluate
