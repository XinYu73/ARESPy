! Module lapack_module defined in file Lapack_module.f90

subroutine f90wrap_diagm(mat, evec, eval, n0, n1, n2, n3, n4)
    use lapack_module, only: diagm
    implicit none
    
    complex(8), intent(in), dimension(n0,n1) :: mat
    complex(8), intent(inout), dimension(n2,n3) :: evec
    real(8), intent(inout), dimension(n4) :: eval
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,0)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,1)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    call diagm(mat=mat, evec=evec, eval=eval)
end subroutine f90wrap_diagm

subroutine f90wrap_generalizeeigen(dime, mata, matb, evec, eval, n0, n1, n2, n3, n4, n5, n6)
    use lapack_module, only: generalizeeigen
    implicit none
    
    integer(8), intent(in) :: dime
    complex(8), intent(in), dimension(n0,n1) :: mata
    complex(8), intent(in), dimension(n2,n3) :: matb
    complex(8), intent(inout), dimension(n4,n5) :: evec
    real(8), intent(inout), dimension(n6) :: eval
    integer :: n0
    !f2py intent(hide), depend(mata) :: n0 = shape(mata,0)
    integer :: n1
    !f2py intent(hide), depend(mata) :: n1 = shape(mata,1)
    integer :: n2
    !f2py intent(hide), depend(matb) :: n2 = shape(matb,0)
    integer :: n3
    !f2py intent(hide), depend(matb) :: n3 = shape(matb,1)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,0)
    integer :: n5
    !f2py intent(hide), depend(evec) :: n5 = shape(evec,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,0)
    call generalizeeigen(dime=dime, matA=mata, matB=matb, evec=evec, eval=eval)
end subroutine f90wrap_generalizeeigen

subroutine f90wrap_orthnorm(mat, n0, n1)
    use lapack_module, only: orthnorm
    implicit none
    
    complex(8), dimension(n0,n1) :: mat
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    call orthnorm(mat=mat)
end subroutine f90wrap_orthnorm

subroutine f90wrap_norm_2(mat, ret_norm_2, k, n0, n1)
    use lapack_module, only: norm_2
    implicit none
    
    complex(8), intent(in), dimension(n0,n1) :: mat
    real(8), intent(out) :: ret_norm_2
    integer(8), intent(in) :: k
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    ret_norm_2 = norm_2(mat=mat, k=k)
end subroutine f90wrap_norm_2

subroutine f90wrap_matmat(mata, matb, opa, opb, matc, n0, n1, n2, n3, n4, n5)
    use lapack_module, only: matmat
    implicit none
    
    complex(8), intent(in), dimension(n0,n1) :: mata
    complex(8), intent(in), dimension(n2,n3) :: matb
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    complex(8), intent(inout), dimension(n4,n5) :: matc
    integer :: n0
    !f2py intent(hide), depend(mata) :: n0 = shape(mata,0)
    integer :: n1
    !f2py intent(hide), depend(mata) :: n1 = shape(mata,1)
    integer :: n2
    !f2py intent(hide), depend(matb) :: n2 = shape(matb,0)
    integer :: n3
    !f2py intent(hide), depend(matb) :: n3 = shape(matb,1)
    integer :: n4
    !f2py intent(hide), depend(matc) :: n4 = shape(matc,0)
    integer :: n5
    !f2py intent(hide), depend(matc) :: n5 = shape(matc,1)
    call matmat(matA=mata, matB=matb, opA=opa, opB=opb, matC=matc)
end subroutine f90wrap_matmat

subroutine f90wrap_invmat(mat, n0, n1)
    use lapack_module, only: invmat
    implicit none
    
    complex(8), intent(inout), dimension(n0,n1) :: mat
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    call invmat(mat=mat)
end subroutine f90wrap_invmat

subroutine f90wrap_orthnorm_real(mat, n0, n1)
    use lapack_module, only: orthnorm_real
    implicit none
    
    real(8), dimension(n0,n1) :: mat
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    call orthnorm_real(mat=mat)
end subroutine f90wrap_orthnorm_real

subroutine f90wrap_diagm_real(mat, evec, eval, n0, n1, n2, n3, n4)
    use lapack_module, only: diagm_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: mat
    real(8), intent(inout), dimension(n2,n3) :: evec
    real(8), intent(inout), dimension(n4) :: eval
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,0)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,1)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    call diagm_real(mat=mat, evec=evec, eval=eval)
end subroutine f90wrap_diagm_real

subroutine f90wrap_generalizeeigen_real(dime, mata, matb, evec, eval, n0, n1, n2, n3, n4, n5, n6)
    use lapack_module, only: generalizeeigen_real
    implicit none
    
    integer(8), intent(in) :: dime
    real(8), intent(in), dimension(n0,n1) :: mata
    real(8), intent(in), dimension(n2,n3) :: matb
    real(8), intent(inout), dimension(n4,n5) :: evec
    real(8), intent(inout), dimension(n6) :: eval
    integer :: n0
    !f2py intent(hide), depend(mata) :: n0 = shape(mata,0)
    integer :: n1
    !f2py intent(hide), depend(mata) :: n1 = shape(mata,1)
    integer :: n2
    !f2py intent(hide), depend(matb) :: n2 = shape(matb,0)
    integer :: n3
    !f2py intent(hide), depend(matb) :: n3 = shape(matb,1)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,0)
    integer :: n5
    !f2py intent(hide), depend(evec) :: n5 = shape(evec,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,0)
    call generalizeeigen_real(dime=dime, matA=mata, matB=matb, evec=evec, eval=eval)
end subroutine f90wrap_generalizeeigen_real

subroutine f90wrap_diagmx_real(mat, dime, num, il, iu, evec, eval, n0, n1, n2, n3, n4)
    use lapack_module, only: diagmx_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: mat
    integer(8), intent(in) :: dime
    integer(8), intent(in) :: num
    integer(8), intent(in) :: il
    integer(8), intent(in) :: iu
    real(8), intent(inout), dimension(n2,n3) :: evec
    real(8), intent(inout), dimension(n4) :: eval
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,0)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,1)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    call diagmx_real(mat=mat, dime=dime, num=num, Il=il, Iu=iu, evec=evec, eval=eval)
end subroutine f90wrap_diagmx_real

subroutine f90wrap_matmat_real(mata, matb, opa, opb, matc, n0, n1, n2, n3, n4, n5)
    use lapack_module, only: matmat_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: mata
    real(8), intent(in), dimension(n2,n3) :: matb
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    real(8), intent(inout), dimension(n4,n5) :: matc
    integer :: n0
    !f2py intent(hide), depend(mata) :: n0 = shape(mata,0)
    integer :: n1
    !f2py intent(hide), depend(mata) :: n1 = shape(mata,1)
    integer :: n2
    !f2py intent(hide), depend(matb) :: n2 = shape(matb,0)
    integer :: n3
    !f2py intent(hide), depend(matb) :: n3 = shape(matb,1)
    integer :: n4
    !f2py intent(hide), depend(matc) :: n4 = shape(matc,0)
    integer :: n5
    !f2py intent(hide), depend(matc) :: n5 = shape(matc,1)
    call matmat_real(matA=mata, matB=matb, opA=opa, opB=opb, matC=matc)
end subroutine f90wrap_matmat_real

subroutine f90wrap_cholesky_factor_real(mat, n0, n1)
    use lapack_module, only: cholesky_factor_real
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: mat
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    call cholesky_factor_real(mat=mat)
end subroutine f90wrap_cholesky_factor_real

subroutine f90wrap_invmat_real(mat, n0, n1)
    use lapack_module, only: invmat_real
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: mat
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    call invmat_real(mat=mat)
end subroutine f90wrap_invmat_real

subroutine f90wrap_lapack_module__get__MMAX(f90wrap_MMAX)
    use lapack_module, only: lapack_module_MMAX => MMAX
    implicit none
    integer(8), intent(out) :: f90wrap_MMAX
    
    f90wrap_MMAX = lapack_module_MMAX
end subroutine f90wrap_lapack_module__get__MMAX

subroutine f90wrap_lapack_module__set__MMAX(f90wrap_MMAX)
    use lapack_module, only: lapack_module_MMAX => MMAX
    implicit none
    integer(8), intent(in) :: f90wrap_MMAX
    
    lapack_module_MMAX = f90wrap_MMAX
end subroutine f90wrap_lapack_module__set__MMAX

subroutine f90wrap_lapack_module__get__NMAX(f90wrap_NMAX)
    use lapack_module, only: lapack_module_NMAX => NMAX
    implicit none
    integer(8), intent(out) :: f90wrap_NMAX
    
    f90wrap_NMAX = lapack_module_NMAX
end subroutine f90wrap_lapack_module__get__NMAX

subroutine f90wrap_lapack_module__set__NMAX(f90wrap_NMAX)
    use lapack_module, only: lapack_module_NMAX => NMAX
    implicit none
    integer(8), intent(in) :: f90wrap_NMAX
    
    lapack_module_NMAX = f90wrap_NMAX
end subroutine f90wrap_lapack_module__set__NMAX

! End of module lapack_module defined in file Lapack_module.f90

