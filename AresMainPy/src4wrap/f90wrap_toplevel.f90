subroutine f90wrap_zbrent(iu, lreset, ebreak, x, y, f, xnew, xnewh, ynew, yd, ifail)
    implicit none
    external zbrent
    
    integer :: iu
    logical :: lreset
    real :: ebreak
    real :: x
    real :: y
    real :: f
    real :: xnew
    real :: xnewh
    real :: ynew
    real :: yd
    integer :: ifail
    call zbrent(iu, lreset, ebreak, x, y, f, xnew, xnewh, ynew, yd, ifail)
end subroutine f90wrap_zbrent

subroutine f90wrap_kardir(nmax, v, basis, n0)
    implicit none
    external kardir
    
    integer(4), intent(in) :: nmax
    real(8), dimension(3,n0) :: v
    real(8), dimension(3,3) :: basis
    integer :: n0
    !f2py intent(hide), depend(v) :: n0 = shape(v,1)
    call kardir(nmax, v, basis)
end subroutine f90wrap_kardir

subroutine f90wrap_ioncgr(iflag, nions, toten, a, b, nfree, posion, posioc, fact, f, factsi, fsif, fl, s, dismax, iu6, &
    iu0, ebreak, ediffg, e1test, lstop2)
    implicit none
    external ioncgr
    
    integer :: iflag
    integer :: nions
    real :: toten
    real :: a
    real :: b
    integer :: nfree
    real :: posion
    real :: posioc
    real :: fact
    real :: f
    real :: factsi
    real :: fsif
    real :: fl
    real :: s
    real :: dismax
    integer :: iu6
    integer :: iu0
    real :: ebreak
    real :: ediffg
    real :: e1test
    logical :: lstop2
    call ioncgr(iflag, nions, toten, a, b, nfree, posion, posioc, fact, f, factsi, fsif, fl, s, dismax, iu6, iu0, ebreak, &
        ediffg, e1test, lstop2)
end subroutine f90wrap_ioncgr

