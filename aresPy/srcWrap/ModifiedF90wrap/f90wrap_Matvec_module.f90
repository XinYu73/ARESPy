! Module matvec_module defined in file Matvec_module.fpp

subroutine f90wrap_cmplx_matvec(ik, veff1d, p, q, dimen, n0, n1, n2)
    use matvec_module, only: cmplx_matvec
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0) :: veff1d
    complex(8), intent(in), dimension(n1) :: p
    complex(8), intent(inout), dimension(n2) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff1d) :: n0 = shape(veff1d,0)
    integer :: n1
    !f2py intent(hide), depend(p) :: n1 = shape(p,0)
    integer :: n2
    !f2py intent(hide), depend(q) :: n2 = shape(q,0)
    call cmplx_matvec(Ik=ik, veff1d=veff1d, p=p, q=q, dimen=dimen)
end subroutine f90wrap_cmplx_matvec

subroutine f90wrap_cmplx_nlocmatvec(ik, p, q, n0, n1)
    use matvec_module, only: cmplx_nlocmatvec
    implicit none
    
    integer(4), intent(in) :: ik
    complex(8), intent(in), dimension(n0) :: p
    complex(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call cmplx_nlocmatvec(Ik=ik, p=p, q=q)
end subroutine f90wrap_cmplx_nlocmatvec

subroutine f90wrap_real_matvec(veff1d, p, q, dimen, n0, n1, n2)
    use matvec_module, only: real_matvec
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff1d
    real(8), intent(in), dimension(n1) :: p
    real(8), intent(inout), dimension(n2) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff1d) :: n0 = shape(veff1d,0)
    integer :: n1
    !f2py intent(hide), depend(p) :: n1 = shape(p,0)
    integer :: n2
    !f2py intent(hide), depend(q) :: n2 = shape(q,0)
    call real_matvec(veff1d=veff1d, p=p, q=q, dimen=dimen)
end subroutine f90wrap_real_matvec

subroutine f90wrap_real_nlocmatvec(p, q, n0, n1)
    use matvec_module, only: real_nlocmatvec
    implicit none
    
    real(8), intent(in), dimension(n0) :: p
    real(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call real_nlocmatvec(p=p, q=q)
end subroutine f90wrap_real_nlocmatvec

subroutine f90wrap_real_matvec_m(mat, p, q, dimen, n0, n1, n2, n3)
    use matvec_module, only: real_matvec_m
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: mat
    real(8), intent(in), dimension(n2) :: p
    real(8), intent(inout), dimension(n3) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    integer :: n2
    !f2py intent(hide), depend(p) :: n2 = shape(p,0)
    integer :: n3
    !f2py intent(hide), depend(q) :: n3 = shape(q,0)
    call real_matvec_m(mat=mat, p=p, q=q, dimen=dimen)
end subroutine f90wrap_real_matvec_m

! End of module matvec_module defined in file Matvec_module.fpp

