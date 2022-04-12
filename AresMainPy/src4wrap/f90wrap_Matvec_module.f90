! Module matvec_module defined in file Matvec_module.fpp

subroutine f90wrap_cmatvec(veff, ik, p, q, dimen, n0, n1, n2, n3, n4)
    use matvec_module, only: cmatvec
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: ik
    complex(8), intent(in), dimension(n3) :: p
    complex(8), intent(inout), dimension(n4) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(p) :: n3 = shape(p,0)
    integer :: n4
    !f2py intent(hide), depend(q) :: n4 = shape(q,0)
    call cmatvec(veff=veff, Ik=ik, p=p, q=q, dimen=dimen)
end subroutine f90wrap_cmatvec

subroutine f90wrap_nlocmatvec(ik, p, q, n0, n1)
    use matvec_module, only: nlocmatvec
    implicit none
    
    integer(4), intent(in) :: ik
    complex(8), intent(in), dimension(n0) :: p
    complex(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec(Ik=ik, p=p, q=q)
end subroutine f90wrap_nlocmatvec

subroutine f90wrap_nlocmatvec_dg(ik, p, q, n0, n1)
    use matvec_module, only: nlocmatvec_dg
    implicit none
    
    integer(4), intent(in) :: ik
    complex(8), intent(in), dimension(n0) :: p
    complex(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec_dg(Ik=ik, p=p, q=q)
end subroutine f90wrap_nlocmatvec_dg

subroutine f90wrap_rmatvec(veff, p, q, dimen, n0, n1, n2, n3, n4)
    use matvec_module, only: rmatvec
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(in), dimension(n3) :: p
    real(8), intent(inout), dimension(n4) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(p) :: n3 = shape(p,0)
    integer :: n4
    !f2py intent(hide), depend(q) :: n4 = shape(q,0)
    call rmatvec(veff=veff, p=p, q=q, dimen=dimen)
end subroutine f90wrap_rmatvec

subroutine f90wrap_nlocmatvec_r(p, q, n0, n1)
    use matvec_module, only: nlocmatvec_r
    implicit none
    
    real(8), intent(in), dimension(n0) :: p
    real(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec_r(p=p, q=q)
end subroutine f90wrap_nlocmatvec_r

subroutine f90wrap_rmatvec_new(mat, p, q, dimen, n0, n1, n2, n3)
    use matvec_module, only: rmatvec_new
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
    call rmatvec_new(mat=mat, p=p, q=q, dimen=dimen)
end subroutine f90wrap_rmatvec_new

subroutine f90wrap_iso_rmatvec(veff_3d, p, q, dimen, n0, n1, n2)
    use matvec_module, only: iso_rmatvec
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff_3d
    real(8), intent(in), dimension(n1) :: p
    real(8), intent(inout), dimension(n2) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(p) :: n1 = shape(p,0)
    integer :: n2
    !f2py intent(hide), depend(q) :: n2 = shape(q,0)
    call iso_rmatvec(veff_3D=veff_3d, p=p, q=q, dimen=dimen)
end subroutine f90wrap_iso_rmatvec

subroutine f90wrap_iso_grid2sphere(grid_v, grid_v1, sphere_v, n0, n1, n2, n3, n4)
    use matvec_module, only: iso_grid2sphere
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: grid_v
    real(8), intent(in), dimension(n3) :: grid_v1
    real(8), intent(inout), dimension(n4) :: sphere_v
    integer :: n0
    !f2py intent(hide), depend(grid_v) :: n0 = shape(grid_v,0)
    integer :: n1
    !f2py intent(hide), depend(grid_v) :: n1 = shape(grid_v,1)
    integer :: n2
    !f2py intent(hide), depend(grid_v) :: n2 = shape(grid_v,2)
    integer :: n3
    !f2py intent(hide), depend(grid_v1) :: n3 = shape(grid_v1,0)
    integer :: n4
    !f2py intent(hide), depend(sphere_v) :: n4 = shape(sphere_v,0)
    call iso_grid2sphere(grid_V=grid_v, grid_V1=grid_v1, sphere_V=sphere_v)
end subroutine f90wrap_iso_grid2sphere

subroutine f90wrap_iso_sphere2grid(sphere_v, grid_v, n0, n1, n2, n3)
    use matvec_module, only: iso_sphere2grid
    implicit none
    
    real(8), intent(in), dimension(n0) :: sphere_v
    real(8), intent(inout), dimension(n1,n2,n3) :: grid_v
    integer :: n0
    !f2py intent(hide), depend(sphere_v) :: n0 = shape(sphere_v,0)
    integer :: n1
    !f2py intent(hide), depend(grid_v) :: n1 = shape(grid_v,0)
    integer :: n2
    !f2py intent(hide), depend(grid_v) :: n2 = shape(grid_v,1)
    integer :: n3
    !f2py intent(hide), depend(grid_v) :: n3 = shape(grid_v,2)
    call iso_sphere2grid(sphere_V=sphere_v, grid_V=grid_v)
end subroutine f90wrap_iso_sphere2grid

subroutine f90wrap_nlocmatvec_iso(p, q, n0, n1)
    use matvec_module, only: nlocmatvec_iso
    implicit none
    
    real(8), intent(in), dimension(n0) :: p
    real(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec_iso(p=p, q=q)
end subroutine f90wrap_nlocmatvec_iso

subroutine f90wrap_nlocmatvec_iso_dg(p, q, n0, n1)
    use matvec_module, only: nlocmatvec_iso_dg
    implicit none
    
    real(8), intent(in), dimension(n0) :: p
    real(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec_iso_dg(p=p, q=q)
end subroutine f90wrap_nlocmatvec_iso_dg

subroutine f90wrap_cmatvec_band(veff, ik, p, q, dimen, n0, n1, n2, n3, n4)
    use matvec_module, only: cmatvec_band
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: ik
    complex(8), intent(in), dimension(n3) :: p
    complex(8), intent(inout), dimension(n4) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(p) :: n3 = shape(p,0)
    integer :: n4
    !f2py intent(hide), depend(q) :: n4 = shape(q,0)
    call cmatvec_band(veff=veff, Ik=ik, p=p, q=q, dimen=dimen)
end subroutine f90wrap_cmatvec_band

subroutine f90wrap_nlocmatvec_band(ik, p, q, n0, n1)
    use matvec_module, only: nlocmatvec_band
    implicit none
    
    integer(4), intent(in) :: ik
    complex(8), intent(in), dimension(n0) :: p
    complex(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec_band(Ik=ik, p=p, q=q)
end subroutine f90wrap_nlocmatvec_band

subroutine f90wrap_rmatvec_gamma(veff, ik, p, q, dimen, n0, n1, n2, n3, n4)
    use matvec_module, only: rmatvec_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n3) :: p
    real(8), intent(inout), dimension(n4) :: q
    integer(4), intent(in) :: dimen
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(p) :: n3 = shape(p,0)
    integer :: n4
    !f2py intent(hide), depend(q) :: n4 = shape(q,0)
    call rmatvec_gamma(veff=veff, Ik=ik, p=p, q=q, dimen=dimen)
end subroutine f90wrap_rmatvec_gamma

subroutine f90wrap_nlocmatvec_gamma(ik, p, q, n0, n1)
    use matvec_module, only: nlocmatvec_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0) :: p
    real(8), intent(inout), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    call nlocmatvec_gamma(Ik=ik, p=p, q=q)
end subroutine f90wrap_nlocmatvec_gamma

! End of module matvec_module defined in file Matvec_module.fpp

