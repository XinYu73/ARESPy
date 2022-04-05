! Module mathsplines defined in file MathSplines.fpp

subroutine f90wrap_spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp, n0, n1, n2)
    use mathsplines, only: spline_cubic_set
    implicit none
    
    integer(4) :: n
    real(8), dimension(n0) :: t
    real(8), dimension(n1) :: y
    integer(4) :: ibcbeg
    real(8) :: ybcbeg
    integer(4) :: ibcend
    real(8) :: ybcend
    real(8), dimension(n2) :: ypp
    integer :: n0
    !f2py intent(hide), depend(t) :: n0 = shape(t,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    integer :: n2
    !f2py intent(hide), depend(ypp) :: n2 = shape(ypp,0)
    call spline_cubic_set(n=n, t=t, y=y, ibcbeg=ibcbeg, ybcbeg=ybcbeg, ibcend=ibcend, ybcend=ybcend, ypp=ypp)
end subroutine f90wrap_spline_cubic_set

subroutine f90wrap_spline_cubic_val(n, t, y, ypp, tval, yval, ypval, yppval, n0, n1, n2)
    use mathsplines, only: spline_cubic_val
    implicit none
    
    integer(4), intent(in) :: n
    real(8), intent(in), dimension(n0) :: t
    real(8), intent(in), dimension(n1) :: y
    real(8), intent(in), dimension(n2) :: ypp
    real(8), intent(in) :: tval
    real(8), intent(out) :: yval
    real(8), intent(out) :: ypval
    real(8), intent(out) :: yppval
    integer :: n0
    !f2py intent(hide), depend(t) :: n0 = shape(t,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    integer :: n2
    !f2py intent(hide), depend(ypp) :: n2 = shape(ypp,0)
    call spline_cubic_val(n=n, t=t, y=y, ypp=ypp, tval=tval, yval=yval, ypval=ypval, yppval=yppval)
end subroutine f90wrap_spline_cubic_val

subroutine f90wrap_rvec_bracket(n, x, xval, left, right, n0)
    use mathsplines, only: rvec_bracket
    implicit none
    
    integer(4), intent(in) :: n
    real(8), intent(in), dimension(n0) :: x
    real(8), intent(in) :: xval
    integer(4), intent(out) :: left
    integer(4), intent(out) :: right
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    call rvec_bracket(n=n, x=x, xval=xval, left=left, right=right)
end subroutine f90wrap_rvec_bracket

subroutine f90wrap_s3_fs(a1, a2, a3, n, b, x, n0, n1, n2, n3, n4)
    use mathsplines, only: s3_fs
    implicit none
    
    real(8), dimension(n0) :: a1
    real(8), dimension(n1) :: a2
    real(8), dimension(n2) :: a3
    integer(4) :: n
    real(8), dimension(n3) :: b
    real(8), dimension(n4) :: x
    integer :: n0
    !f2py intent(hide), depend(a1) :: n0 = shape(a1,0)
    integer :: n1
    !f2py intent(hide), depend(a2) :: n1 = shape(a2,0)
    integer :: n2
    !f2py intent(hide), depend(a3) :: n2 = shape(a3,0)
    integer :: n3
    !f2py intent(hide), depend(b) :: n3 = shape(b,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,0)
    call s3_fs(a1=a1, a2=a2, a3=a3, n=n, b=b, x=x)
end subroutine f90wrap_s3_fs

! End of module mathsplines defined in file MathSplines.fpp

