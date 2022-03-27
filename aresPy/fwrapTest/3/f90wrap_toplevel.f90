subroutine f90wrap_sayhello(comm)
    implicit none
    external sayhello
    
    integer :: comm
    call sayhello(comm)
end subroutine f90wrap_sayhello

