! Module isolateset defined in file Isolate_module.fpp

subroutine f90wrap_isolatevcoulomb_a(rho_s, v_s, n0, n1)
    use isolateset, only: isolatevcoulomb_a
    implicit none
    
    real(8), intent(in), dimension(n0) :: rho_s
    real(8), intent(inout), dimension(n1) :: v_s
    integer :: n0
    !f2py intent(hide), depend(rho_s) :: n0 = shape(rho_s,0)
    integer :: n1
    !f2py intent(hide), depend(v_s) :: n1 = shape(v_s,0)
    call isolatevcoulomb_a(rho_s=rho_s, V_s=v_s)
end subroutine f90wrap_isolatevcoulomb_a

subroutine f90wrap_evaluateboundrygrid_b(rhos, n0, n1, n2)
    use isolateset, only: evaluateboundrygrid_b
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rhos
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    call evaluateboundrygrid_b(rhoS=rhos)
end subroutine f90wrap_evaluateboundrygrid_b

subroutine f90wrap_endevaluate_b
    use isolateset, only: endevaluate_b
    implicit none
    
    call endevaluate_b()
end subroutine f90wrap_endevaluate_b

subroutine f90wrap_calculateboundrypotential_b(rho, lmax, n0, n1, n2)
    use isolateset, only: calculateboundrypotential_b
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    integer(4) :: lmax
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    call calculateboundrypotential_b(rho=rho, Lmax=lmax)
end subroutine f90wrap_calculateboundrypotential_b

subroutine f90wrap_laplaceequationslover_b(iter)
    use isolateset, only: laplaceequationslover_b
    implicit none
    
    integer(4) :: iter
    call laplaceequationslover_b(iter=iter)
end subroutine f90wrap_laplaceequationslover_b

subroutine f90wrap_vcoulombassignment_b(iter, niter_poisson, phi, n0, n1, n2)
    use isolateset, only: vcoulombassignment_b
    implicit none
    
    integer(4), intent(in) :: iter
    integer(4), intent(in) :: niter_poisson
    real(8), intent(inout), dimension(n0,n1,n2) :: phi
    integer :: n0
    !f2py intent(hide), depend(phi) :: n0 = shape(phi,0)
    integer :: n1
    !f2py intent(hide), depend(phi) :: n1 = shape(phi,1)
    integer :: n2
    !f2py intent(hide), depend(phi) :: n2 = shape(phi,2)
    call vcoulombassignment_b(iter=iter, niter_poisson=niter_poisson, phi=phi)
end subroutine f90wrap_vcoulombassignment_b

subroutine f90wrap_cal_qlm(rho, lmax, q_l0, q_lm, n0, n1, n2, n3, n4, n5)
    use isolateset, only: cal_qlm
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    integer(4), intent(in) :: lmax
    real(8), intent(inout), dimension(n3) :: q_l0
    complex(8), intent(inout), dimension(n4,n5) :: q_lm
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(q_l0) :: n3 = shape(q_l0,0)
    integer :: n4
    !f2py intent(hide), depend(q_lm) :: n4 = shape(q_lm,0)
    integer :: n5
    !f2py intent(hide), depend(q_lm) :: n5 = shape(q_lm,1)
    call cal_qlm(rho=rho, Lmax=lmax, Q_l0=q_l0, Q_lm=q_lm)
end subroutine f90wrap_cal_qlm

subroutine f90wrap_car2spe(orig, x, y, z, r, cost, sint, cosp, sinp)
    use isolateset, only: car2spe
    implicit none
    
    real(8), intent(in), dimension(3) :: orig
    integer, intent(in) :: x
    integer, intent(in) :: y
    integer, intent(in) :: z
    real(8), intent(out) :: r
    real(8), intent(out) :: cost
    real(8), intent(out) :: sint
    real(8), intent(out) :: cosp
    real(8), intent(out) :: sinp
    call car2spe(ORIG=orig, x=x, y=y, z=z, r=r, cost=cost, sint=sint, cosp=cosp, sinp=sinp)
end subroutine f90wrap_car2spe

subroutine f90wrap_cal_plm(lmax, x, z, plm, n0, n1)
    use isolateset, only: cal_plm
    implicit none
    
    integer(4), intent(in) :: lmax
    real(8), intent(in) :: x
    real(8), intent(in) :: z
    real(8), intent(inout), dimension(n0,n1) :: plm
    integer :: n0
    !f2py intent(hide), depend(plm) :: n0 = shape(plm,0)
    integer :: n1
    !f2py intent(hide), depend(plm) :: n1 = shape(plm,1)
    call cal_plm(Lmax=lmax, x=x, z=z, plm=plm)
end subroutine f90wrap_cal_plm

subroutine f90wrap_cal1pot(iex, iey, iez, orig, lmax, q_l0, q_lm, phi, n0, n1, n2, n3)
    use isolateset, only: cal1pot
    implicit none
    
    integer(4), intent(in) :: iex
    integer(4), intent(in) :: iey
    integer(4), intent(in) :: iez
    real(8), intent(in), dimension(n0) :: orig
    integer(4), intent(in) :: lmax
    real(8), intent(in), dimension(n1) :: q_l0
    complex(8), intent(in), dimension(n2,n3) :: q_lm
    real(8), intent(out) :: phi
    integer :: n0
    !f2py intent(hide), depend(orig) :: n0 = shape(orig,0)
    integer :: n1
    !f2py intent(hide), depend(q_l0) :: n1 = shape(q_l0,0)
    integer :: n2
    !f2py intent(hide), depend(q_lm) :: n2 = shape(q_lm,0)
    integer :: n3
    !f2py intent(hide), depend(q_lm) :: n3 = shape(q_lm,1)
    call cal1pot(Iex=iex, Iey=iey, Iez=iez, orig=orig, Lmax=lmax, Q_l0=q_l0, Q_lm=q_lm, phi=phi)
end subroutine f90wrap_cal1pot

subroutine f90wrap_apply_boundary(lmax, orig, q_l0, q_lm, bphi, n0, n1, n2, n3, n4, n5)
    use isolateset, only: apply_boundary
    implicit none
    
    integer(4), intent(in) :: lmax
    real(8), intent(in), dimension(3) :: orig
    real(8), intent(in), dimension(n0) :: q_l0
    complex(8), intent(in), dimension(n1,n2) :: q_lm
    real(8), intent(inout), dimension(n3,n4,n5) :: bphi
    integer :: n0
    !f2py intent(hide), depend(q_l0) :: n0 = shape(q_l0,0)
    integer :: n1
    !f2py intent(hide), depend(q_lm) :: n1 = shape(q_lm,0)
    integer :: n2
    !f2py intent(hide), depend(q_lm) :: n2 = shape(q_lm,1)
    integer :: n3
    !f2py intent(hide), depend(bphi) :: n3 = shape(bphi,0)
    integer :: n4
    !f2py intent(hide), depend(bphi) :: n4 = shape(bphi,1)
    integer :: n5
    !f2py intent(hide), depend(bphi) :: n5 = shape(bphi,2)
    call apply_boundary(Lmax=lmax, orig=orig, Q_l0=q_l0, Q_lm=q_lm, bphi=bphi)
end subroutine f90wrap_apply_boundary

subroutine f90wrap_apply_boundary00(rho, bphi, n0, n1, n2, n3, n4, n5)
    use isolateset, only: apply_boundary00
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(inout), dimension(n3,n4,n5) :: bphi
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(bphi) :: n3 = shape(bphi,0)
    integer :: n4
    !f2py intent(hide), depend(bphi) :: n4 = shape(bphi,1)
    integer :: n5
    !f2py intent(hide), depend(bphi) :: n5 = shape(bphi,2)
    call apply_boundary00(rho=rho, bphi=bphi)
end subroutine f90wrap_apply_boundary00

subroutine f90wrap_laplb(ford, bf, f, n0, n1, n2, n3, n4, n5)
    use isolateset, only: laplb
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), intent(in), dimension(n0,n1,n2) :: bf
    real(8), intent(inout), dimension(n3,n4,n5) :: f
    integer :: n0
    !f2py intent(hide), depend(bf) :: n0 = shape(bf,0)
    integer :: n1
    !f2py intent(hide), depend(bf) :: n1 = shape(bf,1)
    integer :: n2
    !f2py intent(hide), depend(bf) :: n2 = shape(bf,2)
    integer :: n3
    !f2py intent(hide), depend(f) :: n3 = shape(f,0)
    integer :: n4
    !f2py intent(hide), depend(f) :: n4 = shape(f,1)
    integer :: n5
    !f2py intent(hide), depend(f) :: n5 = shape(f,2)
    call laplb(Ford=ford, bf=bf, f=f)
end subroutine f90wrap_laplb

subroutine f90wrap_lapla(ford, dxyz, f, af, n0, n1, n2, n3, n4, n5)
    use isolateset, only: lapla
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), intent(in), dimension(3) :: dxyz
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: af
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(af) :: n3 = shape(af,0)
    integer :: n4
    !f2py intent(hide), depend(af) :: n4 = shape(af,1)
    integer :: n5
    !f2py intent(hide), depend(af) :: n5 = shape(af,2)
    call lapla(Ford=ford, dxyz=dxyz, f=f, af=af)
end subroutine f90wrap_lapla

subroutine f90wrap_lapla_coe(ford, dxyz, coe, n0)
    use isolateset, only: lapla_coe
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), intent(in), dimension(3) :: dxyz
    real(8), dimension(n0) :: coe
    integer :: n0
    !f2py intent(hide), depend(coe) :: n0 = shape(coe,0)
    call lapla_coe(Ford=ford, dxyz=dxyz, coe=coe)
end subroutine f90wrap_lapla_coe

subroutine f90wrap_damped_jacobi_iterate(ford, dxyz, func, rhs, n0, n1, n2, n3, n4, n5)
    use isolateset, only: damped_jacobi_iterate
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), dimension(3) :: dxyz
    real(8), intent(inout), dimension(n0,n1,n2) :: func
    real(8), intent(in), dimension(n3,n4,n5) :: rhs
    integer :: n0
    !f2py intent(hide), depend(func) :: n0 = shape(func,0)
    integer :: n1
    !f2py intent(hide), depend(func) :: n1 = shape(func,1)
    integer :: n2
    !f2py intent(hide), depend(func) :: n2 = shape(func,2)
    integer :: n3
    !f2py intent(hide), depend(rhs) :: n3 = shape(rhs,0)
    integer :: n4
    !f2py intent(hide), depend(rhs) :: n4 = shape(rhs,1)
    integer :: n5
    !f2py intent(hide), depend(rhs) :: n5 = shape(rhs,2)
    call damped_jacobi_iterate(Ford=ford, dxyz=dxyz, func=func, rhs=rhs)
end subroutine f90wrap_damped_jacobi_iterate

subroutine f90wrap_gauss_seidel_iterate(ford, dxyz, func, rhs, n0, n1, n2, n3, n4, n5)
    use isolateset, only: gauss_seidel_iterate
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), dimension(3) :: dxyz
    real(8), intent(inout), dimension(n0,n1,n2) :: func
    real(8), intent(in), dimension(n3,n4,n5) :: rhs
    integer :: n0
    !f2py intent(hide), depend(func) :: n0 = shape(func,0)
    integer :: n1
    !f2py intent(hide), depend(func) :: n1 = shape(func,1)
    integer :: n2
    !f2py intent(hide), depend(func) :: n2 = shape(func,2)
    integer :: n3
    !f2py intent(hide), depend(rhs) :: n3 = shape(rhs,0)
    integer :: n4
    !f2py intent(hide), depend(rhs) :: n4 = shape(rhs,1)
    integer :: n5
    !f2py intent(hide), depend(rhs) :: n5 = shape(rhs,2)
    call gauss_seidel_iterate(Ford=ford, dxyz=dxyz, func=func, rhs=rhs)
end subroutine f90wrap_gauss_seidel_iterate

subroutine f90wrap_refine(var, nxyz, nxyz_fine, var_fine, n0, n1, n2, n3, n4, n5)
    use isolateset, only: refine
    implicit none
    
    real(8), dimension(n0,n1,n2) :: var
    integer(4), dimension(3) :: nxyz
    integer(4), dimension(3) :: nxyz_fine
    real(8), dimension(n3,n4,n5) :: var_fine
    integer :: n0
    !f2py intent(hide), depend(var) :: n0 = shape(var,0)
    integer :: n1
    !f2py intent(hide), depend(var) :: n1 = shape(var,1)
    integer :: n2
    !f2py intent(hide), depend(var) :: n2 = shape(var,2)
    integer :: n3
    !f2py intent(hide), depend(var_fine) :: n3 = shape(var_fine,0)
    integer :: n4
    !f2py intent(hide), depend(var_fine) :: n4 = shape(var_fine,1)
    integer :: n5
    !f2py intent(hide), depend(var_fine) :: n5 = shape(var_fine,2)
    call refine(var=var, nxyz=nxyz, nxyz_fine=nxyz_fine, var_fine=var_fine)
end subroutine f90wrap_refine

subroutine f90wrap_interpolate(var1, nxyz1, nxyz2, var2, n0, n1, n2, n3, n4, n5)
    use isolateset, only: interpolate
    implicit none
    
    real(8), dimension(n0,n1,n2) :: var1
    integer(4), dimension(3) :: nxyz1
    integer(4), dimension(3) :: nxyz2
    real(8), dimension(n3,n4,n5) :: var2
    integer :: n0
    !f2py intent(hide), depend(var1) :: n0 = shape(var1,0)
    integer :: n1
    !f2py intent(hide), depend(var1) :: n1 = shape(var1,1)
    integer :: n2
    !f2py intent(hide), depend(var1) :: n2 = shape(var1,2)
    integer :: n3
    !f2py intent(hide), depend(var2) :: n3 = shape(var2,0)
    integer :: n4
    !f2py intent(hide), depend(var2) :: n4 = shape(var2,1)
    integer :: n5
    !f2py intent(hide), depend(var2) :: n5 = shape(var2,2)
    call interpolate(var1=var1, nxyz1=nxyz1, nxyz2=nxyz2, var2=var2)
end subroutine f90wrap_interpolate

subroutine f90wrap_smooth(f, n0, n1, n2)
    use isolateset, only: smooth
    implicit none
    
    real(8), dimension(n0,n1,n2) :: f
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    call smooth(f=f)
end subroutine f90wrap_smooth

subroutine f90wrap_restrict(var, nxyz, nxyz_coarse, var_coarse, n0, n1, n2, n3, n4, n5)
    use isolateset, only: restrict
    implicit none
    
    real(8), dimension(n0,n1,n2) :: var
    integer(4), dimension(3) :: nxyz
    integer(4), dimension(3) :: nxyz_coarse
    real(8), dimension(n3,n4,n5) :: var_coarse
    integer :: n0
    !f2py intent(hide), depend(var) :: n0 = shape(var,0)
    integer :: n1
    !f2py intent(hide), depend(var) :: n1 = shape(var,1)
    integer :: n2
    !f2py intent(hide), depend(var) :: n2 = shape(var,2)
    integer :: n3
    !f2py intent(hide), depend(var_coarse) :: n3 = shape(var_coarse,0)
    integer :: n4
    !f2py intent(hide), depend(var_coarse) :: n4 = shape(var_coarse,1)
    integer :: n5
    !f2py intent(hide), depend(var_coarse) :: n5 = shape(var_coarse,2)
    call restrict(var=var, nxyz=nxyz, nxyz_coarse=nxyz_coarse, var_coarse=var_coarse)
end subroutine f90wrap_restrict

subroutine f90wrap_residual(ford, dxyz, f, rhs, res, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use isolateset, only: residual
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), dimension(3) :: dxyz
    real(8), dimension(n0,n1,n2) :: f
    real(8), dimension(n3,n4,n5) :: rhs
    real(8), dimension(n6,n7,n8) :: res
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(rhs) :: n3 = shape(rhs,0)
    integer :: n4
    !f2py intent(hide), depend(rhs) :: n4 = shape(rhs,1)
    integer :: n5
    !f2py intent(hide), depend(rhs) :: n5 = shape(rhs,2)
    integer :: n6
    !f2py intent(hide), depend(res) :: n6 = shape(res,0)
    integer :: n7
    !f2py intent(hide), depend(res) :: n7 = shape(res,1)
    integer :: n8
    !f2py intent(hide), depend(res) :: n8 = shape(res,2)
    call residual(Ford=ford, dxyz=dxyz, f=f, rhs=rhs, res=res)
end subroutine f90wrap_residual

subroutine f90wrap_v_cycle(ford, f, rhs, dxyz, nvc, n0, n1, n2, n3, n4, n5)
    use isolateset, only: v_cycle
    implicit none
    
    integer(4), intent(in) :: ford
    real(8), dimension(n0,n1,n2) :: f
    real(8), dimension(n3,n4,n5) :: rhs
    real(8), dimension(3) :: dxyz
    integer(4) :: nvc
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(rhs) :: n3 = shape(rhs,0)
    integer :: n4
    !f2py intent(hide), depend(rhs) :: n4 = shape(rhs,1)
    integer :: n5
    !f2py intent(hide), depend(rhs) :: n5 = shape(rhs,2)
    call v_cycle(Ford=ford, f=f, rhs=rhs, dxyz=dxyz, nvc=nvc)
end subroutine f90wrap_v_cycle

subroutine f90wrap_csix_coe(dxyz, coe, n0, n1)
    use isolateset, only: csix_coe
    implicit none
    
    real(8), intent(in), dimension(n0) :: dxyz
    real(8), intent(inout), dimension(n1) :: coe
    integer :: n0
    !f2py intent(hide), depend(dxyz) :: n0 = shape(dxyz,0)
    integer :: n1
    !f2py intent(hide), depend(coe) :: n1 = shape(coe,0)
    call csix_coe(dxyz=dxyz, coe=coe)
end subroutine f90wrap_csix_coe

subroutine f90wrap_lapla_csix(f, af, coe, n0, n1, n2, n3, n4, n5, n6)
    use isolateset, only: lapla_csix
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: af
    real(8), intent(in), dimension(n6) :: coe
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(af) :: n3 = shape(af,0)
    integer :: n4
    !f2py intent(hide), depend(af) :: n4 = shape(af,1)
    integer :: n5
    !f2py intent(hide), depend(af) :: n5 = shape(af,2)
    integer :: n6
    !f2py intent(hide), depend(coe) :: n6 = shape(coe,0)
    call lapla_csix(f=f, af=af, coe=coe)
end subroutine f90wrap_lapla_csix

subroutine f90wrap_laplb_csix(bff, f, n0, n1, n2, n3, n4, n5)
    use isolateset, only: laplb_csix
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: bff
    real(8), intent(inout), dimension(n3,n4,n5) :: f
    integer :: n0
    !f2py intent(hide), depend(bff) :: n0 = shape(bff,0)
    integer :: n1
    !f2py intent(hide), depend(bff) :: n1 = shape(bff,1)
    integer :: n2
    !f2py intent(hide), depend(bff) :: n2 = shape(bff,2)
    integer :: n3
    !f2py intent(hide), depend(f) :: n3 = shape(f,0)
    integer :: n4
    !f2py intent(hide), depend(f) :: n4 = shape(f,1)
    integer :: n5
    !f2py intent(hide), depend(f) :: n5 = shape(f,2)
    call laplb_csix(bff=bff, f=f)
end subroutine f90wrap_laplb_csix

subroutine f90wrap_cfour_coe(gaps, coe, n0)
    use isolateset, only: cfour_coe
    implicit none
    
    real(8), intent(in), dimension(3) :: gaps
    real(8), dimension(n0) :: coe
    integer :: n0
    !f2py intent(hide), depend(coe) :: n0 = shape(coe,0)
    call cfour_coe(gaps=gaps, coe=coe)
end subroutine f90wrap_cfour_coe

subroutine f90wrap_lapla_cfour(u, au, coe, n0, n1, n2, n3, n4, n5, n6)
    use isolateset, only: lapla_cfour
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: u
    real(8), intent(inout), dimension(n3,n4,n5) :: au
    real(8), intent(in), dimension(n6) :: coe
    integer :: n0
    !f2py intent(hide), depend(u) :: n0 = shape(u,0)
    integer :: n1
    !f2py intent(hide), depend(u) :: n1 = shape(u,1)
    integer :: n2
    !f2py intent(hide), depend(u) :: n2 = shape(u,2)
    integer :: n3
    !f2py intent(hide), depend(au) :: n3 = shape(au,0)
    integer :: n4
    !f2py intent(hide), depend(au) :: n4 = shape(au,1)
    integer :: n5
    !f2py intent(hide), depend(au) :: n5 = shape(au,2)
    integer :: n6
    !f2py intent(hide), depend(coe) :: n6 = shape(coe,0)
    call lapla_cfour(u=u, au=au, coe=coe)
end subroutine f90wrap_lapla_cfour

subroutine f90wrap_laplb_cfour(rhs, nrhs, n0, n1, n2, n3, n4, n5)
    use isolateset, only: laplb_cfour
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rhs
    real(8), intent(inout), dimension(n3,n4,n5) :: nrhs
    integer :: n0
    !f2py intent(hide), depend(rhs) :: n0 = shape(rhs,0)
    integer :: n1
    !f2py intent(hide), depend(rhs) :: n1 = shape(rhs,1)
    integer :: n2
    !f2py intent(hide), depend(rhs) :: n2 = shape(rhs,2)
    integer :: n3
    !f2py intent(hide), depend(nrhs) :: n3 = shape(nrhs,0)
    integer :: n4
    !f2py intent(hide), depend(nrhs) :: n4 = shape(nrhs,1)
    integer :: n5
    !f2py intent(hide), depend(nrhs) :: n5 = shape(nrhs,2)
    call laplb_cfour(rhs=rhs, nrhs=nrhs)
end subroutine f90wrap_laplb_cfour

subroutine f90wrap_calclm(lmax, clm, n0, n1)
    use isolateset, only: calclm
    implicit none
    
    integer(4) :: lmax
    real(8), dimension(n0,n1) :: clm
    integer :: n0
    !f2py intent(hide), depend(clm) :: n0 = shape(clm,0)
    integer :: n1
    !f2py intent(hide), depend(clm) :: n1 = shape(clm,1)
    call calclm(Lmax=lmax, clm=clm)
end subroutine f90wrap_calclm

subroutine f90wrap_c(l, ret_c, m)
    use isolateset, only: c
    implicit none
    
    integer(4) :: l
    real(8), intent(out) :: ret_c
    integer(4) :: m
    ret_c = c(l=l, m=m)
end subroutine f90wrap_c

subroutine f90wrap_rcs(n1, n2, n3, dr, cost, sint, cns, indx, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, &
    n13, n14)
    use isolateset, only: rcs
    implicit none
    
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    real(8), intent(inout), dimension(n0,n1,n2) :: dr
    real(8), intent(inout), dimension(n3,n4,n5) :: cost
    real(8), intent(inout), dimension(n6,n7,n8) :: sint
    complex(8), intent(inout), dimension(n9,n10,n11) :: cns
    integer(4), intent(inout), dimension(n12,n13,n14) :: indx
    integer :: n0
    !f2py intent(hide), depend(dr) :: n0 = shape(dr,0)
    integer :: n1
    !f2py intent(hide), depend(dr) :: n1 = shape(dr,1)
    integer :: n2
    !f2py intent(hide), depend(dr) :: n2 = shape(dr,2)
    integer :: n3
    !f2py intent(hide), depend(cost) :: n3 = shape(cost,0)
    integer :: n4
    !f2py intent(hide), depend(cost) :: n4 = shape(cost,1)
    integer :: n5
    !f2py intent(hide), depend(cost) :: n5 = shape(cost,2)
    integer :: n6
    !f2py intent(hide), depend(sint) :: n6 = shape(sint,0)
    integer :: n7
    !f2py intent(hide), depend(sint) :: n7 = shape(sint,1)
    integer :: n8
    !f2py intent(hide), depend(sint) :: n8 = shape(sint,2)
    integer :: n9
    !f2py intent(hide), depend(cns) :: n9 = shape(cns,0)
    integer :: n10
    !f2py intent(hide), depend(cns) :: n10 = shape(cns,1)
    integer :: n11
    !f2py intent(hide), depend(cns) :: n11 = shape(cns,2)
    integer :: n12
    !f2py intent(hide), depend(indx) :: n12 = shape(indx,0)
    integer :: n13
    !f2py intent(hide), depend(indx) :: n13 = shape(indx,1)
    integer :: n14
    !f2py intent(hide), depend(indx) :: n14 = shape(indx,2)
    call rcs(n1=n1, n2=n2, n3=n3, dr=dr, cost=cost, sint=sint, cns=cns, indx=indx)
end subroutine f90wrap_rcs

subroutine f90wrap_vhartree_fmm(rhos, vcoulomb, n0, n1, n2, n3, n4, n5)
    use isolateset, only: vhartree_fmm
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rhos
    real(8), intent(inout), dimension(n3,n4,n5) :: vcoulomb
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(vcoulomb) :: n3 = shape(vcoulomb,0)
    integer :: n4
    !f2py intent(hide), depend(vcoulomb) :: n4 = shape(vcoulomb,1)
    integer :: n5
    !f2py intent(hide), depend(vcoulomb) :: n5 = shape(vcoulomb,2)
    call vhartree_fmm(rhoS=rhos, VCoulomb=vcoulomb)
end subroutine f90wrap_vhartree_fmm

subroutine f90wrap_vhartree_direct(rhos, vcoulomb, n0, n1, n2, n3, n4, n5)
    use isolateset, only: vhartree_direct
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rhos
    real(8), intent(inout), dimension(n3,n4,n5) :: vcoulomb
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(vcoulomb) :: n3 = shape(vcoulomb,0)
    integer :: n4
    !f2py intent(hide), depend(vcoulomb) :: n4 = shape(vcoulomb,1)
    integer :: n5
    !f2py intent(hide), depend(vcoulomb) :: n5 = shape(vcoulomb,2)
    call vhartree_direct(rhoS=rhos, VCoulomb=vcoulomb)
end subroutine f90wrap_vhartree_direct

subroutine f90wrap_cal_vsrcpot(size_src, src, pos_src, size_tar, pot_tar, pos_tar, n0, n1, n2, n3)
    use isolateset, only: cal_vsrcpot
    implicit none
    
    integer(4) :: size_src
    real(8), dimension(n0) :: src
    real(8), dimension(3,n1) :: pos_src
    integer(4) :: size_tar
    real(8), dimension(n2) :: pot_tar
    real(8), dimension(3,n3) :: pos_tar
    integer :: n0
    !f2py intent(hide), depend(src) :: n0 = shape(src,0)
    integer :: n1
    !f2py intent(hide), depend(pos_src) :: n1 = shape(pos_src,1)
    integer :: n2
    !f2py intent(hide), depend(pot_tar) :: n2 = shape(pot_tar,0)
    integer :: n3
    !f2py intent(hide), depend(pos_tar) :: n3 = shape(pos_tar,1)
    call cal_vsrcpot(size_src=size_src, src=src, pos_src=pos_src, size_tar=size_tar, pot_tar=pot_tar, pos_tar=pos_tar)
end subroutine f90wrap_cal_vsrcpot

subroutine f90wrap_vhart_iso(rho, vhart, n0, n1, n2, n3, n4, n5)
    use isolateset, only: vhart_iso
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(inout), dimension(n3,n4,n5) :: vhart
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(vhart) :: n3 = shape(vhart,0)
    integer :: n4
    !f2py intent(hide), depend(vhart) :: n4 = shape(vhart,1)
    integer :: n5
    !f2py intent(hide), depend(vhart) :: n5 = shape(vhart,2)
    call vhart_iso(rho=rho, vhart=vhart)
end subroutine f90wrap_vhart_iso

subroutine f90wrap_poisson_mcm(rho, pot, n0, n1, n2, n3, n4, n5)
    use isolateset, only: poisson_mcm
    implicit none
    
    real(8), dimension(n0,n1,n2) :: rho
    real(8), dimension(n3,n4,n5) :: pot
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(pot) :: n3 = shape(pot,0)
    integer :: n4
    !f2py intent(hide), depend(pot) :: n4 = shape(pot,1)
    integer :: n5
    !f2py intent(hide), depend(pot) :: n5 = shape(pot,2)
    call poisson_mcm(rho=rho, pot=pot)
end subroutine f90wrap_poisson_mcm

subroutine f90wrap_mcm_calvcorr(q_l0, q_lm, v_corr, n0, n1, n2, n3, n4, n5)
    use isolateset, only: mcm_calvcorr
    implicit none
    
    real(8), dimension(n0) :: q_l0
    complex(8), dimension(n1,n2) :: q_lm
    real(8), dimension(n3,n4,n5) :: v_corr
    integer :: n0
    !f2py intent(hide), depend(q_l0) :: n0 = shape(q_l0,0)
    integer :: n1
    !f2py intent(hide), depend(q_lm) :: n1 = shape(q_lm,0)
    integer :: n2
    !f2py intent(hide), depend(q_lm) :: n2 = shape(q_lm,1)
    integer :: n3
    !f2py intent(hide), depend(v_corr) :: n3 = shape(v_corr,0)
    integer :: n4
    !f2py intent(hide), depend(v_corr) :: n4 = shape(v_corr,1)
    integer :: n5
    !f2py intent(hide), depend(v_corr) :: n5 = shape(v_corr,2)
    call mcm_calvcorr(Q_l0=q_l0, Q_lm=q_lm, V_corr=v_corr)
end subroutine f90wrap_mcm_calvcorr

subroutine f90wrap_mcm_cal_naux(q_l0, q_lm, rhoaux, n0, n1, n2, n3)
    use isolateset, only: mcm_cal_naux
    implicit none
    
    real(8), dimension(n0) :: q_l0
    complex(8), dimension(n1,n2) :: q_lm
    real(8), dimension(n3) :: rhoaux
    integer :: n0
    !f2py intent(hide), depend(q_l0) :: n0 = shape(q_l0,0)
    integer :: n1
    !f2py intent(hide), depend(q_lm) :: n1 = shape(q_lm,0)
    integer :: n2
    !f2py intent(hide), depend(q_lm) :: n2 = shape(q_lm,1)
    integer :: n3
    !f2py intent(hide), depend(rhoaux) :: n3 = shape(rhoaux,0)
    call mcm_cal_naux(Q_l0=q_l0, Q_lm=q_lm, rhoaux=rhoaux)
end subroutine f90wrap_mcm_cal_naux

subroutine f90wrap_cal_ylm(mycost, mysint, e_phi, ylm, n0, n1)
    use isolateset, only: cal_ylm
    implicit none
    
    real(8) :: mycost
    real(8) :: mysint
    complex(8) :: e_phi
    complex(8), dimension(n0,n1) :: ylm
    integer :: n0
    !f2py intent(hide), depend(ylm) :: n0 = shape(ylm,0)
    integer :: n1
    !f2py intent(hide), depend(ylm) :: n1 = shape(ylm,1)
    call cal_ylm(mycost=mycost, mysint=mysint, e_phi=e_phi, ylm=ylm)
end subroutine f90wrap_cal_ylm

subroutine f90wrap_rcs_rec(ng1, ng2, ng3, cost, sint, cns, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use isolateset, only: rcs_rec
    implicit none
    
    integer(4), intent(in) :: ng1
    integer(4), intent(in) :: ng2
    integer(4), intent(in) :: ng3
    real(8), intent(inout), dimension(n0,n1,n2) :: cost
    real(8), intent(inout), dimension(n3,n4,n5) :: sint
    complex(8), intent(inout), dimension(n6,n7,n8) :: cns
    integer :: n0
    !f2py intent(hide), depend(cost) :: n0 = shape(cost,0)
    integer :: n1
    !f2py intent(hide), depend(cost) :: n1 = shape(cost,1)
    integer :: n2
    !f2py intent(hide), depend(cost) :: n2 = shape(cost,2)
    integer :: n3
    !f2py intent(hide), depend(sint) :: n3 = shape(sint,0)
    integer :: n4
    !f2py intent(hide), depend(sint) :: n4 = shape(sint,1)
    integer :: n5
    !f2py intent(hide), depend(sint) :: n5 = shape(sint,2)
    integer :: n6
    !f2py intent(hide), depend(cns) :: n6 = shape(cns,0)
    integer :: n7
    !f2py intent(hide), depend(cns) :: n7 = shape(cns,1)
    integer :: n8
    !f2py intent(hide), depend(cns) :: n8 = shape(cns,2)
    call rcs_rec(ng1=ng1, ng2=ng2, ng3=ng3, cost=cost, sint=sint, cns=cns)
end subroutine f90wrap_rcs_rec

subroutine f90wrap_car2spe_rec(orig, x, y, z, cost, sint, cosp, sinp)
    use isolateset, only: car2spe_rec
    implicit none
    
    real(8), intent(in), dimension(3) :: orig
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: z
    real(8), intent(out) :: cost
    real(8), intent(out) :: sint
    real(8), intent(out) :: cosp
    real(8), intent(out) :: sinp
    call car2spe_rec(ORIG=orig, x=x, y=y, z=z, cost=cost, sint=sint, cosp=cosp, sinp=sinp)
end subroutine f90wrap_car2spe_rec

subroutine f90wrap_integrade_i(l, ret_integrade_i, x)
    use isolateset, only: integrade_i
    implicit none
    
    integer(4) :: l
    real(8), intent(out) :: ret_integrade_i
    real(8), intent(in) :: x
    ret_integrade_i = integrade_i(l=l, x=x)
end subroutine f90wrap_integrade_i

subroutine f90wrap_gmg_cg_hart(rhos, vcoulomb, n0, n1, n2, n3, n4, n5)
    use isolateset, only: gmg_cg_hart
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rhos
    real(8), intent(inout), dimension(n3,n4,n5) :: vcoulomb
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(vcoulomb) :: n3 = shape(vcoulomb,0)
    integer :: n4
    !f2py intent(hide), depend(vcoulomb) :: n4 = shape(vcoulomb,1)
    integer :: n5
    !f2py intent(hide), depend(vcoulomb) :: n5 = shape(vcoulomb,2)
    call gmg_cg_hart(rhoS=rhos, VCoulomb=vcoulomb)
end subroutine f90wrap_gmg_cg_hart

subroutine f90wrap_cutoff_method(rho, vh, n0, n1, n2, n3, n4, n5)
    use isolateset, only: cutoff_method
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(inout), dimension(n3,n4,n5) :: vh
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(vh) :: n3 = shape(vh,0)
    integer :: n4
    !f2py intent(hide), depend(vh) :: n4 = shape(vh,1)
    integer :: n5
    !f2py intent(hide), depend(vh) :: n5 = shape(vh,2)
    call cutoff_method(rho=rho, Vh=vh)
end subroutine f90wrap_cutoff_method

subroutine f90wrap_isolateset__array__Center(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_center => center
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(isolateset_Center)
    dloc = loc(isolateset_Center)
end subroutine f90wrap_isolateset__array__Center

subroutine f90wrap_isolateset__array__CellLengthN(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_celllengthn => celllengthn
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(isolateset_CellLengthN)
    dloc = loc(isolateset_CellLengthN)
end subroutine f90wrap_isolateset__array__CellLengthN

subroutine f90wrap_isolateset__get__Thickness(f90wrap_Thickness)
    use isolateset, only: isolateset_Thickness => Thickness
    implicit none
    integer(4), intent(out) :: f90wrap_Thickness
    
    f90wrap_Thickness = isolateset_Thickness
end subroutine f90wrap_isolateset__get__Thickness

subroutine f90wrap_isolateset__set__Thickness(f90wrap_Thickness)
    use isolateset, only: isolateset_Thickness => Thickness
    implicit none
    integer(4), intent(in) :: f90wrap_Thickness
    
    isolateset_Thickness = f90wrap_Thickness
end subroutine f90wrap_isolateset__set__Thickness

subroutine f90wrap_isolateset__get__CellLeft(f90wrap_CellLeft)
    use isolateset, only: isolateset_CellLeft => CellLeft
    implicit none
    integer(4), intent(out) :: f90wrap_CellLeft
    
    f90wrap_CellLeft = isolateset_CellLeft
end subroutine f90wrap_isolateset__get__CellLeft

subroutine f90wrap_isolateset__set__CellLeft(f90wrap_CellLeft)
    use isolateset, only: isolateset_CellLeft => CellLeft
    implicit none
    integer(4), intent(in) :: f90wrap_CellLeft
    
    isolateset_CellLeft = f90wrap_CellLeft
end subroutine f90wrap_isolateset__set__CellLeft

subroutine f90wrap_isolateset__array__CellRight(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_cellright => cellright
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(isolateset_CellRight)
    dloc = loc(isolateset_CellRight)
end subroutine f90wrap_isolateset__array__CellRight

subroutine f90wrap_isolateset__array__BandRight(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_bandright => bandright
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(isolateset_BandRight)
    dloc = loc(isolateset_BandRight)
end subroutine f90wrap_isolateset__array__BandRight

subroutine f90wrap_isolateset__array__BoundaryrhoS(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_boundaryrhos => boundaryrhos
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_BoundaryrhoS)) then
        dshape(1:3) = shape(isolateset_BoundaryrhoS)
        dloc = loc(isolateset_BoundaryrhoS)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__BoundaryrhoS

subroutine f90wrap_isolateset__array__BoundaryVCoulomb(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_boundaryvcoulomb => boundaryvcoulomb
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_BoundaryVCoulomb)) then
        dshape(1:3) = shape(isolateset_BoundaryVCoulomb)
        dloc = loc(isolateset_BoundaryVCoulomb)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__BoundaryVCoulomb

subroutine f90wrap_isolateset__array__ap(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_ap => ap
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_ap)) then
        dshape(1:3) = shape(isolateset_ap)
        dloc = loc(isolateset_ap)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__ap

subroutine f90wrap_isolateset__array__aphi(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_aphi => aphi
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_aphi)) then
        dshape(1:3) = shape(isolateset_aphi)
        dloc = loc(isolateset_aphi)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__aphi

subroutine f90wrap_isolateset__array__res(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_res => res
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_res)) then
        dshape(1:3) = shape(isolateset_res)
        dloc = loc(isolateset_res)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__res

subroutine f90wrap_isolateset__array__res1(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_res1 => res1
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_res1)) then
        dshape(1:3) = shape(isolateset_res1)
        dloc = loc(isolateset_res1)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__res1

subroutine f90wrap_isolateset__array__VCoulomb_old(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use struct_module
    use m_time_evaluate
    use poisson_isf
    use isolateset, only: isolateset_vcoulomb_old => vcoulomb_old
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(isolateset_VCoulomb_old)) then
        dshape(1:3) = shape(isolateset_VCoulomb_old)
        dloc = loc(isolateset_VCoulomb_old)
    else
        dloc = 0
    end if
end subroutine f90wrap_isolateset__array__VCoulomb_old

! End of module isolateset defined in file Isolate_module.fpp

