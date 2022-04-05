# 1 "MPI_array_io.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "MPI_array_io.f90"
MODULE array_io
  USE constants
  IMPLICIT NONE
  type out_label
     character(len=100)  :: label
     integer(I4B)      :: num
  endtype out_label
  INTEGER(I4B),save,private :: ntime=0

  type(out_label)      :: outfile(100)
  INTERFACE output
     MODULE PROCEDURE output_r,output_i
  END INTERFACE output
  INTERFACE input
     MODULE PROCEDURE input_r,input_i
  END INTERFACE input
CONTAINS

  SUBROUTINE init_outfile()
    IMPLICIT NONE
    INTEGER(I4B) :: i

    do i=1,100
       outfile(i)%label=""
       outfile(i)%num=0
    enddo
  ENDSUBROUTINE init_outfile

  SUBROUTINE output_r(size,array,name)
    IMPLICIT NONE
    INTEGER(I4B)        :: size
    REAL(DP),intent(in) :: array(size)
    CHARACTER(len=*)    :: name

    CHARACTER(len=100)  :: label
    CHARACTER(len=3)    :: label_suffix
    INTEGER(I4B)        :: i
    LOGICAL             :: add_more_label

    label=trim(name)
    add_more_label=.true.

!> # label classify and count
!> ## the first and when the label is old
    if(ntime==0)then
       CALL init_outfile()
       ntime=1
       outfile(1)%label=label
       outfile(1)%num=outfile(1)%num+1
       write(label_suffix,'(I3)')outfile(1)%num
       add_more_label=.false.
    else
       do i=1,ntime
          if(label==outfile(i)%label)then
             outfile(i)%num=outfile(i)%num+1
             write(label_suffix,'(I3)')outfile(i)%num
             add_more_label=.false.
          endif
       enddo
    endif

!> ## when the label is new
    if(add_more_label)then
       ntime=ntime+1
       outfile(ntime)%label=label
       outfile(ntime)%num=outfile(ntime)%num+1
       write(label_suffix,'(I3)')outfile(ntime)%num
    endif

!> # write the array
    open(222,file=trim(adjustl(label))//'_'//trim(adjustl(label_suffix)))
    write(222,*)array
    close(222)
  ENDSUBROUTINE output_r

  SUBROUTINE input_r(size,array,name)
    IMPLICIT NONE
    INTEGER(I4B)        :: size
    REAL(DP),intent(out) :: array(size)
    CHARACTER(len=*)    :: name

    open(222,file=name)
    read(222,*)array
    close(222)
  ENDSUBROUTINE input_r

  SUBROUTINE output_i(size,array,name)
    IMPLICIT NONE
    INTEGER(I4B)        :: size
    INTEGER(I4B),intent(in) :: array(size)
    CHARACTER(len=*)    :: name

    open(222,file=name)
    write(222,*)array
    close(222)
  ENDSUBROUTINE output_i

  SUBROUTINE input_i(size,array,name)
    IMPLICIT NONE
    INTEGER(I4B)        :: size
    INTEGER(I4B),intent(out) :: array(size)
    CHARACTER(len=*)    :: name

    open(222,file=name)
    read(222,*)array
    close(222)
  ENDSUBROUTINE input_i

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ENDMODULE array_io
