! Module xc_module defined in file XC_functional.fpp

subroutine f90wrap_xc_functional(nps, rhos, vxc, exc, n0, n1, n2, n3)
    use xc_module, only: xc_functional
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3) :: vxc
    real(8), optional, intent(inout) :: exc
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(vxc) :: n2 = shape(vxc,0)
    integer :: n3
    !f2py intent(hide), depend(vxc) :: n3 = shape(vxc,1)
    call xc_functional(nps=nps, rhoS=rhos, vxc=vxc, exc=exc)
end subroutine f90wrap_xc_functional

subroutine f90wrap_libxc_lda_set(nps, rhos, vxcs, exc, n0, n1, n2, n3)
    use xc_module, only: libxc_lda_set
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3) :: vxcs
    real(8), optional, intent(inout) :: exc
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(vxcs) :: n2 = shape(vxcs,0)
    integer :: n3
    !f2py intent(hide), depend(vxcs) :: n3 = shape(vxcs,1)
    call libxc_lda_set(nps=nps, rhoS=rhos, vxcS=vxcs, exc=exc)
end subroutine f90wrap_libxc_lda_set

subroutine f90wrap_libxc_lda_x(nps, rhos, ex, vxs, n0, n1, n2, n3)
    use xc_module, only: libxc_lda_x
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(out) :: ex
    real(8), optional, intent(inout), dimension(n2,n3) :: vxs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(vxs) :: n2 = shape(vxs,0)
    integer :: n3
    !f2py intent(hide), depend(vxs) :: n3 = shape(vxs,1)
    call libxc_lda_x(nps=nps, rhoS=rhos, ex=ex, vxS=vxs)
end subroutine f90wrap_libxc_lda_x

subroutine f90wrap_libxc_vwn1rpa_c(nps, rhos, ec, vcs, n0, n1, n2, n3)
    use xc_module, only: libxc_vwn1rpa_c
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(out) :: ec
    real(8), optional, intent(inout), dimension(n2,n3) :: vcs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(vcs) :: n2 = shape(vcs,0)
    integer :: n3
    !f2py intent(hide), depend(vcs) :: n3 = shape(vcs,1)
    call libxc_vwn1rpa_c(nps=nps, rhoS=rhos, ec=ec, vcS=vcs)
end subroutine f90wrap_libxc_vwn1rpa_c

subroutine f90wrap_libxc_gga_set(nps, rhos, vxcs, exc, n0, n1, n2, n3)
    use xc_module, only: libxc_gga_set
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3) :: vxcs
    real(8), optional, intent(inout) :: exc
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(vxcs) :: n2 = shape(vxcs,0)
    integer :: n3
    !f2py intent(hide), depend(vxcs) :: n3 = shape(vxcs,1)
    call libxc_gga_set(nps=nps, rhoS=rhos, vxcS=vxcs, exc=exc)
end subroutine f90wrap_libxc_gga_set

subroutine f90wrap_libxc_b88_x(nps, rhos, sigma, ex, vxs, n0, n1, n2, n3, n4, n5)
    use xc_module, only: libxc_b88_x
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(in), dimension(n2,n3) :: sigma
    real(8), intent(out) :: ex
    real(8), optional, intent(inout), dimension(n4,n5) :: vxs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(sigma) :: n2 = shape(sigma,0)
    integer :: n3
    !f2py intent(hide), depend(sigma) :: n3 = shape(sigma,1)
    integer :: n4
    !f2py intent(hide), depend(vxs) :: n4 = shape(vxs,0)
    integer :: n5
    !f2py intent(hide), depend(vxs) :: n5 = shape(vxs,1)
    call libxc_b88_x(nps=nps, rhoS=rhos, sigma=sigma, ex=ex, vxS=vxs)
end subroutine f90wrap_libxc_b88_x

subroutine f90wrap_libxc_lyp_c(nps, rhos, sigma, ec, vcs, n0, n1, n2, n3, n4, n5)
    use xc_module, only: libxc_lyp_c
    implicit none
    
    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(in), dimension(n2,n3) :: sigma
    real(8), intent(out) :: ec
    real(8), optional, intent(inout), dimension(n4,n5) :: vcs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(sigma) :: n2 = shape(sigma,0)
    integer :: n3
    !f2py intent(hide), depend(sigma) :: n3 = shape(sigma,1)
    integer :: n4
    !f2py intent(hide), depend(vcs) :: n4 = shape(vcs,0)
    integer :: n5
    !f2py intent(hide), depend(vcs) :: n5 = shape(vcs,1)
    call libxc_lyp_c(nps=nps, rhoS=rhos, sigma=sigma, ec=ec, vcS=vcs)
end subroutine f90wrap_libxc_lyp_c

! End of module xc_module defined in file XC_functional.fpp

