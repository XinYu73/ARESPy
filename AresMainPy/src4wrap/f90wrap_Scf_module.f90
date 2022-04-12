! Module scf_module defined in file Scf_module.fpp

subroutine f90wrap_arpackscf(rhos, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: arpackscf
    implicit none
    
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    complex(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call arpackscf(rhoS=rhos, psi=psi, eval=eval)
end subroutine f90wrap_arpackscf

subroutine f90wrap_solver_spin(rhos, psi, eval, nev, diagtol, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: solver_spin
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    complex(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer, intent(in) :: nev
    real(8), intent(in) :: diagtol
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call solver_spin(rhoS=rhos, psi=psi, eval=eval, nev=nev, diagTOL=diagtol)
end subroutine f90wrap_solver_spin

subroutine f90wrap_ksolver(veff, kpsi, keval, nev, tol, n0, n1, n2, n3, n4, n5, n6, n7)
    use scf_module, only: ksolver
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4,n5) :: kpsi
    real(8), intent(inout), dimension(n6,n7) :: keval
    integer(4), intent(in) :: nev
    real(8), intent(in) :: tol
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(kpsi) :: n3 = shape(kpsi,0)
    integer :: n4
    !f2py intent(hide), depend(kpsi) :: n4 = shape(kpsi,1)
    integer :: n5
    !f2py intent(hide), depend(kpsi) :: n5 = shape(kpsi,2)
    integer :: n6
    !f2py intent(hide), depend(keval) :: n6 = shape(keval,0)
    integer :: n7
    !f2py intent(hide), depend(keval) :: n7 = shape(keval,1)
    call ksolver(veff=veff, kpsi=kpsi, keval=keval, nev=nev, TOL=tol)
end subroutine f90wrap_ksolver

subroutine f90wrap_smear_updaterho(nev, ne, psi, eval, wke_l, fme, ets, rhos, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, &
    n10)
    use scf_module, only: smear_updaterho
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(in) :: ne
    complex(8), intent(in), dimension(n0,n1,n2,n3) :: psi
    real(8), intent(in), dimension(n4,n5,n6) :: eval
    real(8), intent(inout), dimension(n7,n8,n9) :: wke_l
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    real(8), intent(inout), dimension(n10) :: rhos
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,3)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,2)
    integer :: n7
    !f2py intent(hide), depend(wke_l) :: n7 = shape(wke_l,0)
    integer :: n8
    !f2py intent(hide), depend(wke_l) :: n8 = shape(wke_l,1)
    integer :: n9
    !f2py intent(hide), depend(wke_l) :: n9 = shape(wke_l,2)
    integer :: n10
    !f2py intent(hide), depend(rhos) :: n10 = shape(rhos,0)
    call smear_updaterho(nev=nev, ne=ne, psi=psi, eval=eval, wke_l=wke_l, fme=fme, ets=ets, rhoS=rhos)
end subroutine f90wrap_smear_updaterho

subroutine f90wrap_chefsi(rhos, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: chefsi
    implicit none
    
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    complex(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call chefsi(rhoS=rhos, psi=psi, eval=eval)
end subroutine f90wrap_chefsi

subroutine f90wrap_filter_spin(veff, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: filter_spin
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: veff
    complex(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(veff) :: n3 = shape(veff,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call filter_spin(veff=veff, psi=psi, eval=eval)
end subroutine f90wrap_filter_spin

subroutine f90wrap_arpackscf_r(rhos, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: arpackscf_r
    implicit none
    
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call arpackscf_r(rhoS=rhos, psi=psi, eval=eval)
end subroutine f90wrap_arpackscf_r

subroutine f90wrap_solver_spin_r(rhos, psi, eval, nev, diagtol, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use scf_module, only: solver_spin_r
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6) :: psi
    real(8), intent(inout), dimension(n7,n8) :: eval
    integer, intent(in) :: nev
    real(8), intent(in) :: diagtol
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(eval) :: n7 = shape(eval,0)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,1)
    call solver_spin_r(rhoS=rhos, psi=psi, eval=eval, nev=nev, diagTOL=diagtol)
end subroutine f90wrap_solver_spin_r

subroutine f90wrap_real_smear_updaterho(nev, ne, psi, eval, wke_l, fme, ets, rhos, n0, n1, n2, n3, n4, n5, n6, n7, n8, &
    n9, n10)
    use scf_module, only: real_smear_updaterho
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(in) :: ne
    real(8), intent(in), dimension(n0,n1,n2,n3) :: psi
    real(8), intent(in), dimension(n4,n5,n6) :: eval
    real(8), intent(inout), dimension(n7,n8,n9) :: wke_l
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    real(8), intent(inout), dimension(n10) :: rhos
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,3)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,2)
    integer :: n7
    !f2py intent(hide), depend(wke_l) :: n7 = shape(wke_l,0)
    integer :: n8
    !f2py intent(hide), depend(wke_l) :: n8 = shape(wke_l,1)
    integer :: n9
    !f2py intent(hide), depend(wke_l) :: n9 = shape(wke_l,2)
    integer :: n10
    !f2py intent(hide), depend(rhos) :: n10 = shape(rhos,0)
    call real_smear_updaterho(nev=nev, ne=ne, psi=psi, eval=eval, wke_l=wke_l, fme=fme, ets=ets, rhoS=rhos)
end subroutine f90wrap_real_smear_updaterho

subroutine f90wrap_bvk_chefsi(rhos, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: bvk_chefsi
    implicit none
    
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call bvk_chefsi(rhoS=rhos, psi=psi, eval=eval)
end subroutine f90wrap_bvk_chefsi

subroutine f90wrap_bvk_filter_spin(rhos, x, d, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use scf_module, only: bvk_filter_spin
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6) :: x
    real(8), intent(inout), dimension(n7,n8) :: d
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,0)
    integer :: n5
    !f2py intent(hide), depend(x) :: n5 = shape(x,1)
    integer :: n6
    !f2py intent(hide), depend(x) :: n6 = shape(x,2)
    integer :: n7
    !f2py intent(hide), depend(d) :: n7 = shape(d,0)
    integer :: n8
    !f2py intent(hide), depend(d) :: n8 = shape(d,1)
    call bvk_filter_spin(rhoS=rhos, X=x, D=d)
end subroutine f90wrap_bvk_filter_spin

subroutine f90wrap_prr_filter_spin(nfs, nfe, rhos, x, efr, fme, ets, pbar, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, &
    n11, n12)
    use scf_module, only: prr_filter_spin
    implicit none
    
    integer(4), intent(in) :: nfs
    real(8), intent(inout) :: nfe
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6) :: x
    real(8), intent(inout), dimension(n7,n8,n9) :: efr
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    real(8), intent(inout), dimension(n10,n11,n12) :: pbar
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,0)
    integer :: n5
    !f2py intent(hide), depend(x) :: n5 = shape(x,1)
    integer :: n6
    !f2py intent(hide), depend(x) :: n6 = shape(x,2)
    integer :: n7
    !f2py intent(hide), depend(efr) :: n7 = shape(efr,0)
    integer :: n8
    !f2py intent(hide), depend(efr) :: n8 = shape(efr,1)
    integer :: n9
    !f2py intent(hide), depend(efr) :: n9 = shape(efr,2)
    integer :: n10
    !f2py intent(hide), depend(pbar) :: n10 = shape(pbar,0)
    integer :: n11
    !f2py intent(hide), depend(pbar) :: n11 = shape(pbar,1)
    integer :: n12
    !f2py intent(hide), depend(pbar) :: n12 = shape(pbar,2)
    call prr_filter_spin(Nfs=nfs, Nfe=nfe, rhoS=rhos, X=x, Efr=efr, fme=fme, ets=ets, Pbar=pbar)
end subroutine f90wrap_prr_filter_spin

subroutine f90wrap_prr_chefsi(rhos, phi, n0, n1, n2, n3, n4, n5, n6, n7)
    use scf_module, only: prr_chefsi
    implicit none
    
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: phi
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(phi) :: n4 = shape(phi,0)
    integer :: n5
    !f2py intent(hide), depend(phi) :: n5 = shape(phi,1)
    integer :: n6
    !f2py intent(hide), depend(phi) :: n6 = shape(phi,2)
    integer :: n7
    !f2py intent(hide), depend(phi) :: n7 = shape(phi,3)
    call prr_chefsi(rhoS=rhos, Phi=phi)
end subroutine f90wrap_prr_chefsi

subroutine f90wrap_iso_solver_spin_r(rhos, psi, eval, nev, diagtol, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use scf_module, only: iso_solver_spin_r
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6) :: psi
    real(8), intent(inout), dimension(n7,n8) :: eval
    integer, intent(in) :: nev
    real(8), intent(in) :: diagtol
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(eval) :: n7 = shape(eval,0)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,1)
    call iso_solver_spin_r(rhoS=rhos, psi=psi, eval=eval, nev=nev, diagTOL=diagtol)
end subroutine f90wrap_iso_solver_spin_r

subroutine f90wrap_iso_smear_updaterho(nev, ne, psi, eval, wke_l, fme, ets, rhos_out, n0, n1, n2, n3, n4, n5, n6, n7, &
    n8, n9)
    use scf_module, only: iso_smear_updaterho
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(in) :: ne
    real(8), intent(in), dimension(n0,n1,n2) :: psi
    real(8), intent(in), dimension(n3,n4) :: eval
    real(8), intent(inout), dimension(n5,n6,n7) :: wke_l
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    real(8), intent(inout), dimension(n8,n9) :: rhos_out
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(eval) :: n3 = shape(eval,0)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,1)
    integer :: n5
    !f2py intent(hide), depend(wke_l) :: n5 = shape(wke_l,0)
    integer :: n6
    !f2py intent(hide), depend(wke_l) :: n6 = shape(wke_l,1)
    integer :: n7
    !f2py intent(hide), depend(wke_l) :: n7 = shape(wke_l,2)
    integer :: n8
    !f2py intent(hide), depend(rhos_out) :: n8 = shape(rhos_out,0)
    integer :: n9
    !f2py intent(hide), depend(rhos_out) :: n9 = shape(rhos_out,1)
    call iso_smear_updaterho(nev=nev, ne=ne, psi=psi, eval=eval, wke_l=wke_l, fme=fme, ets=ets, rhoS_out=rhos_out)
end subroutine f90wrap_iso_smear_updaterho

subroutine f90wrap_iso_chefsi(rhos, psi, eval, n0, n1, n2, n3, n4, n5, n6)
    use scf_module, only: iso_chefsi
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3,n4) :: psi
    real(8), intent(inout), dimension(n5,n6) :: eval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,0)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,1)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,2)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,0)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,1)
    call iso_chefsi(rhoS=rhos, psi=psi, eval=eval)
end subroutine f90wrap_iso_chefsi

subroutine f90wrap_iso_filter_spin(rhos, x, d, n0, n1, n2, n3, n4, n5, n6)
    use scf_module, only: iso_filter_spin
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3,n4) :: x
    real(8), intent(inout), dimension(n5,n6) :: d
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,0)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,1)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,2)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    integer :: n6
    !f2py intent(hide), depend(d) :: n6 = shape(d,1)
    call iso_filter_spin(rhoS=rhos, X=x, D=d)
end subroutine f90wrap_iso_filter_spin

subroutine f90wrap_chefsi_gamma(rhos, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: chefsi_gamma
    implicit none
    
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call chefsi_gamma(rhoS=rhos, psi=psi, eval=eval)
end subroutine f90wrap_chefsi_gamma

subroutine f90wrap_smear_updaterho_gamma(nev, ne, psi, eval, wke_l, fme, ets, rhos, n0, n1, n2, n3, n4, n5, n6, n7, n8, &
    n9, n10)
    use scf_module, only: smear_updaterho_gamma
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(in) :: ne
    real(8), intent(in), dimension(n0,n1,n2,n3) :: psi
    real(8), intent(in), dimension(n4,n5,n6) :: eval
    real(8), intent(inout), dimension(n7,n8,n9) :: wke_l
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    real(8), intent(inout), dimension(n10) :: rhos
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,3)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,2)
    integer :: n7
    !f2py intent(hide), depend(wke_l) :: n7 = shape(wke_l,0)
    integer :: n8
    !f2py intent(hide), depend(wke_l) :: n8 = shape(wke_l,1)
    integer :: n9
    !f2py intent(hide), depend(wke_l) :: n9 = shape(wke_l,2)
    integer :: n10
    !f2py intent(hide), depend(rhos) :: n10 = shape(rhos,0)
    call smear_updaterho_gamma(nev=nev, ne=ne, psi=psi, eval=eval, wke_l=wke_l, fme=fme, ets=ets, rhoS=rhos)
end subroutine f90wrap_smear_updaterho_gamma

subroutine f90wrap_filter_spin_gamma(veff, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use scf_module, only: filter_spin_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: veff
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9,n10) :: eval
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(veff) :: n3 = shape(veff,3)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,0)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,1)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,2)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call filter_spin_gamma(veff=veff, psi=psi, eval=eval)
end subroutine f90wrap_filter_spin_gamma

! End of module scf_module defined in file Scf_module.fpp

