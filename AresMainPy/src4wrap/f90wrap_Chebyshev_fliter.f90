! Module chebyshev_module defined in file Chebyshev_fliter.fpp

subroutine f90wrap_cheby_filter_rr(ik, veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: cheby_filter_rr
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call cheby_filter_rr(Ik=ik, veff=veff, X=x, D=d)
end subroutine f90wrap_cheby_filter_rr

subroutine f90wrap_grayleigh_ritz(ik, veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: grayleigh_ritz
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call grayleigh_ritz(Ik=ik, veff=veff, X=x, D=d)
end subroutine f90wrap_grayleigh_ritz

subroutine f90wrap_chebyshev_filter(ik, veff, x, m, a, b, n0, n1, n2, n3, n4)
    use chebyshev_module, only: chebyshev_filter
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    call chebyshev_filter(Ik=ik, veff=veff, X=x, m=m, a=a, b=b)
end subroutine f90wrap_chebyshev_filter

subroutine f90wrap_chebyshev_filter_scaled(ik, veff, x, m, a, b, al, n0, n1, n2, n3, n4)
    use chebyshev_module, only: chebyshev_filter_scaled
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: al
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    call chebyshev_filter_scaled(Ik=ik, veff=veff, X=x, m=m, a=a, b=b, al=al)
end subroutine f90wrap_chebyshev_filter_scaled

subroutine f90wrap_cal_hx(ik, veff, nst, v, hv, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: cal_hx
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nst
    complex(8), intent(in), dimension(n3,n4) :: v
    complex(8), intent(inout), dimension(n5,n6) :: hv
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(v) :: n3 = shape(v,0)
    integer :: n4
    !f2py intent(hide), depend(v) :: n4 = shape(v,1)
    integer :: n5
    !f2py intent(hide), depend(hv) :: n5 = shape(hv,0)
    integer :: n6
    !f2py intent(hide), depend(hv) :: n6 = shape(hv,1)
    call cal_hx(Ik=ik, veff=veff, nst=nst, V=v, HV=hv)
end subroutine f90wrap_cal_hx

subroutine f90wrap_rayleigh_quotient(ik, veff, nst, x, xhx, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: rayleigh_quotient
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nst
    complex(8), dimension(n3,n4) :: x
    complex(8), dimension(n5,n6) :: xhx
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(xhx) :: n5 = shape(xhx,0)
    integer :: n6
    !f2py intent(hide), depend(xhx) :: n6 = shape(xhx,1)
    call rayleigh_quotient(Ik=ik, veff=veff, nst=nst, x=x, xhx=xhx)
end subroutine f90wrap_rayleigh_quotient

subroutine f90wrap_rayleigh_ritz(ik, veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: rayleigh_ritz
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call rayleigh_ritz(Ik=ik, veff=veff, X=x, D=d)
end subroutine f90wrap_rayleigh_ritz

subroutine f90wrap_estupb(k, ik, veff, vec, b, n0, n1, n2, n3)
    use chebyshev_module, only: estupb
    implicit none
    
    integer(4), intent(in) :: k
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(in), dimension(n3) :: vec
    real(8), intent(out) :: b
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(vec) :: n3 = shape(vec,0)
    call estupb(k=k, Ik=ik, veff=veff, vec=vec, b=b)
end subroutine f90wrap_estupb

subroutine f90wrap_inituplow(k, ik, veff, v, eval, a, b, al, n0, n1, n2, n3, n4)
    use chebyshev_module, only: inituplow
    implicit none
    
    integer(4), intent(in) :: k
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3) :: v
    real(8), intent(inout), dimension(n4) :: eval
    real(8), intent(out) :: a
    real(8), intent(out) :: b
    real(8), intent(out) :: al
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(v) :: n3 = shape(v,0)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    call inituplow(k=k, Ik=ik, veff=veff, v=v, eval=eval, a=a, b=b, al=al)
end subroutine f90wrap_inituplow

subroutine f90wrap_first_scfstep_filter(ik, veff, x, eval, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: first_scfstep_filter
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: eval
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,0)
    call first_scfstep_filter(Ik=ik, veff=veff, X=x, eval=eval)
end subroutine f90wrap_first_scfstep_filter

subroutine f90wrap_randomfast(psi, nev, radius, n0, n1)
    use chebyshev_module, only: randomfast
    implicit none
    
    complex(8), intent(inout), dimension(n0,n1) :: psi
    integer(4), intent(in) :: nev
    real(8), intent(in) :: radius
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    call randomfast(psi=psi, nev=nev, radius=radius)
end subroutine f90wrap_randomfast

subroutine f90wrap_first_chescf_random(rhos, nev, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use chebyshev_module, only: first_chescf_random
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    integer(4), intent(in) :: nev
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
    call first_chescf_random(rhoS=rhos, nev=nev, psi=psi, eval=eval)
end subroutine f90wrap_first_chescf_random

subroutine f90wrap_first_chescf(veff, psi_ran, nev, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    use chebyshev_module, only: first_chescf
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    complex(8), intent(in), dimension(n3,n4) :: psi_ran
    integer(4), intent(in) :: nev
    complex(8), intent(inout), dimension(n5,n6,n7) :: psi
    real(8), intent(inout), dimension(n8,n9) :: eval
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(psi_ran) :: n3 = shape(psi_ran,0)
    integer :: n4
    !f2py intent(hide), depend(psi_ran) :: n4 = shape(psi_ran,1)
    integer :: n5
    !f2py intent(hide), depend(psi) :: n5 = shape(psi,0)
    integer :: n6
    !f2py intent(hide), depend(psi) :: n6 = shape(psi,1)
    integer :: n7
    !f2py intent(hide), depend(psi) :: n7 = shape(psi,2)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    call first_chescf(veff=veff, psi_ran=psi_ran, nev=nev, psi=psi, eval=eval)
end subroutine f90wrap_first_chescf

subroutine f90wrap_cheby_init_sto(nmax, npw, initx_sto, n0, n1)
    use chebyshev_module, only: cheby_init_sto
    implicit none
    
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: npw
    complex(8), intent(inout), dimension(n0,n1) :: initx_sto
    integer :: n0
    !f2py intent(hide), depend(initx_sto) :: n0 = shape(initx_sto,0)
    integer :: n1
    !f2py intent(hide), depend(initx_sto) :: n1 = shape(initx_sto,1)
    call cheby_init_sto(Nmax=nmax, Npw=npw, initX_sto=initx_sto)
end subroutine f90wrap_cheby_init_sto

subroutine f90wrap_first_chescf_sto(veff, nev, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use chebyshev_module, only: first_chescf_sto
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: veff
    integer(4), intent(in) :: nev
    complex(8), dimension(n4,n5,n6,n7) :: psi
    real(8), dimension(n8,n9,n10) :: eval
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
    call first_chescf_sto(veff=veff, nev=nev, psi=psi, eval=eval)
end subroutine f90wrap_first_chescf_sto

subroutine f90wrap_first_subspace_stopw(ik, veff, nmax, nev, initx, x, d, n0, n1, n2, n3, n4, n5, n6, n7)
    use chebyshev_module, only: first_subspace_stopw
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nev
    complex(8), intent(inout), dimension(n3,n4) :: initx
    complex(8), intent(inout), dimension(n5,n6) :: x
    real(8), intent(inout), dimension(n7) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(initx) :: n3 = shape(initx,0)
    integer :: n4
    !f2py intent(hide), depend(initx) :: n4 = shape(initx,1)
    integer :: n5
    !f2py intent(hide), depend(x) :: n5 = shape(x,0)
    integer :: n6
    !f2py intent(hide), depend(x) :: n6 = shape(x,1)
    integer :: n7
    !f2py intent(hide), depend(d) :: n7 = shape(d,0)
    call first_subspace_stopw(Ik=ik, veff=veff, Nmax=nmax, nev=nev, initX=initx, X=x, D=d)
end subroutine f90wrap_first_subspace_stopw

subroutine f90wrap_bvk_first_chescf_sto_rand(rhos, nev, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use chebyshev_module, only: bvk_first_chescf_sto_rand
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    integer(4), intent(in) :: nev
    real(8), dimension(n4,n5,n6) :: psi
    real(8), dimension(n7,n8) :: eval
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
    call bvk_first_chescf_sto_rand(rhoS=rhos, nev=nev, psi=psi, eval=eval)
end subroutine f90wrap_bvk_first_chescf_sto_rand

subroutine f90wrap_bvk_cheby_init_sto_rand(nmax, nrand, initx_sto, n0, n1)
    use chebyshev_module, only: bvk_cheby_init_sto_rand
    implicit none
    
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nrand
    real(8), intent(inout), dimension(n0,n1) :: initx_sto
    integer :: n0
    !f2py intent(hide), depend(initx_sto) :: n0 = shape(initx_sto,0)
    integer :: n1
    !f2py intent(hide), depend(initx_sto) :: n1 = shape(initx_sto,1)
    call bvk_cheby_init_sto_rand(Nmax=nmax, Nrand=nrand, initX_sto=initx_sto)
end subroutine f90wrap_bvk_cheby_init_sto_rand

subroutine f90wrap_first_subspace_sto_rand(veff, nmax, nev, initx, x, d, n0, n1, n2, n3, n4, n5, n6, n7)
    use chebyshev_module, only: first_subspace_sto_rand
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n3,n4) :: initx
    real(8), intent(inout), dimension(n5,n6) :: x
    real(8), intent(inout), dimension(n7) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(initx) :: n3 = shape(initx,0)
    integer :: n4
    !f2py intent(hide), depend(initx) :: n4 = shape(initx,1)
    integer :: n5
    !f2py intent(hide), depend(x) :: n5 = shape(x,0)
    integer :: n6
    !f2py intent(hide), depend(x) :: n6 = shape(x,1)
    integer :: n7
    !f2py intent(hide), depend(d) :: n7 = shape(d,0)
    call first_subspace_sto_rand(veff=veff, Nmax=nmax, nev=nev, initX=initx, X=x, D=d)
end subroutine f90wrap_first_subspace_sto_rand

subroutine f90wrap_rayleigh_quotient_real(veff, nst, x, xhx, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: rayleigh_quotient_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nst
    real(8), dimension(n3,n4) :: x
    real(8), dimension(n5,n6) :: xhx
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(xhx) :: n5 = shape(xhx,0)
    integer :: n6
    !f2py intent(hide), depend(xhx) :: n6 = shape(xhx,1)
    call rayleigh_quotient_real(veff=veff, nst=nst, x=x, xhx=xhx)
end subroutine f90wrap_rayleigh_quotient_real

subroutine f90wrap_cal_hx_real(veff, nst, v, hv, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: cal_hx_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n3,n4) :: v
    real(8), intent(inout), dimension(n5,n6) :: hv
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(v) :: n3 = shape(v,0)
    integer :: n4
    !f2py intent(hide), depend(v) :: n4 = shape(v,1)
    integer :: n5
    !f2py intent(hide), depend(hv) :: n5 = shape(hv,0)
    integer :: n6
    !f2py intent(hide), depend(hv) :: n6 = shape(hv,1)
    call cal_hx_real(veff=veff, nst=nst, V=v, HV=hv)
end subroutine f90wrap_cal_hx_real

subroutine f90wrap_estupb_real(k, veff, vec, b, n0, n1, n2, n3)
    use chebyshev_module, only: estupb_real
    implicit none
    
    integer(4), intent(in) :: k
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(in), dimension(n3) :: vec
    real(8), intent(out) :: b
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(vec) :: n3 = shape(vec,0)
    call estupb_real(k=k, veff=veff, vec=vec, b=b)
end subroutine f90wrap_estupb_real

subroutine f90wrap_chebyshev_filter_scaled_real(veff, x, m, a, b, al, n0, n1, n2, n3, n4)
    use chebyshev_module, only: chebyshev_filter_scaled_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: al
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    call chebyshev_filter_scaled_real(veff=veff, X=x, m=m, a=a, b=b, al=al)
end subroutine f90wrap_chebyshev_filter_scaled_real

subroutine f90wrap_grayleigh_ritz_real(veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: grayleigh_ritz_real
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call grayleigh_ritz_real(veff=veff, X=x, D=d)
end subroutine f90wrap_grayleigh_ritz_real

subroutine f90wrap_cheby_filtering_grrr(veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: cheby_filtering_grrr
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call cheby_filtering_grrr(veff=veff, X=x, D=d)
end subroutine f90wrap_cheby_filtering_grrr

subroutine f90wrap_cheby_filtering_prrr(veff, x, nfs, cfr, efr, n0, n1, n2, n3, n4, n5, n6, n7)
    use chebyshev_module, only: cheby_filtering_prrr
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    integer(4), intent(in) :: nfs
    real(8), intent(inout), dimension(n5,n6) :: cfr
    real(8), intent(inout), dimension(n7) :: efr
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(cfr) :: n5 = shape(cfr,0)
    integer :: n6
    !f2py intent(hide), depend(cfr) :: n6 = shape(cfr,1)
    integer :: n7
    !f2py intent(hide), depend(efr) :: n7 = shape(efr,0)
    call cheby_filtering_prrr(veff=veff, X=x, Nfs=nfs, Cfr=cfr, Efr=efr)
end subroutine f90wrap_cheby_filtering_prrr

subroutine f90wrap_partialrayleighritz(veff, x, nfs, cfr, efr, n0, n1, n2, n3, n4, n5, n6, n7)
    use chebyshev_module, only: partialrayleighritz
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    integer(4) :: nfs
    real(8), intent(inout), dimension(n5,n6) :: cfr
    real(8), intent(inout), dimension(n7) :: efr
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(cfr) :: n5 = shape(cfr,0)
    integer :: n6
    !f2py intent(hide), depend(cfr) :: n6 = shape(cfr,1)
    integer :: n7
    !f2py intent(hide), depend(efr) :: n7 = shape(efr,0)
    call partialrayleighritz(veff=veff, X=x, Nfs=nfs, Cfr=cfr, Efr=efr)
end subroutine f90wrap_partialrayleighritz

subroutine f90wrap_prr_orthnorm(veff, x, nfs, cfr, efr, n0, n1, n2, n3, n4, n5, n6, n7)
    use chebyshev_module, only: prr_orthnorm
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    integer(4) :: nfs
    real(8), intent(inout), dimension(n5,n6) :: cfr
    real(8), intent(inout), dimension(n7) :: efr
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(cfr) :: n5 = shape(cfr,0)
    integer :: n6
    !f2py intent(hide), depend(cfr) :: n6 = shape(cfr,1)
    integer :: n7
    !f2py intent(hide), depend(efr) :: n7 = shape(efr,0)
    call prr_orthnorm(veff=veff, X=x, Nfs=nfs, Cfr=cfr, Efr=efr)
end subroutine f90wrap_prr_orthnorm

subroutine f90wrap_iso_first_chescf_sto_rand(rhos, nev, psi, eval, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: iso_first_chescf_sto_rand
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    integer(4), intent(in) :: nev
    real(8), dimension(n2,n3,n4) :: psi
    real(8), dimension(n5,n6) :: eval
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
    call iso_first_chescf_sto_rand(rhoS=rhos, nev=nev, psi=psi, eval=eval)
end subroutine f90wrap_iso_first_chescf_sto_rand

subroutine f90wrap_iso_cheby_init_sto_rand(nmax, nrand, initx_sto, n0, n1)
    use chebyshev_module, only: iso_cheby_init_sto_rand
    implicit none
    
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nrand
    real(8), intent(inout), dimension(n0,n1) :: initx_sto
    integer :: n0
    !f2py intent(hide), depend(initx_sto) :: n0 = shape(initx_sto,0)
    integer :: n1
    !f2py intent(hide), depend(initx_sto) :: n1 = shape(initx_sto,1)
    call iso_cheby_init_sto_rand(Nmax=nmax, Nrand=nrand, initX_sto=initx_sto)
end subroutine f90wrap_iso_cheby_init_sto_rand

subroutine f90wrap_iso_cheby_init_rand(nmax, nrand, initx_sto, n0, n1)
    use chebyshev_module, only: iso_cheby_init_rand
    implicit none
    
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nrand
    real(8), intent(inout), dimension(n0,n1) :: initx_sto
    integer :: n0
    !f2py intent(hide), depend(initx_sto) :: n0 = shape(initx_sto,0)
    integer :: n1
    !f2py intent(hide), depend(initx_sto) :: n1 = shape(initx_sto,1)
    call iso_cheby_init_rand(Nmax=nmax, Nrand=nrand, initX_sto=initx_sto)
end subroutine f90wrap_iso_cheby_init_rand

subroutine f90wrap_iso_first_subspace_sto_rand(veff_3d, nmax, nev, initx, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: iso_first_subspace_sto_rand
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff_3d
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n1,n2) :: initx
    real(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(initx) :: n1 = shape(initx,0)
    integer :: n2
    !f2py intent(hide), depend(initx) :: n2 = shape(initx,1)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call iso_first_subspace_sto_rand(veff_3D=veff_3d, Nmax=nmax, nev=nev, initX=initx, X=x, D=d)
end subroutine f90wrap_iso_first_subspace_sto_rand

subroutine f90wrap_rayleigh_quotient_iso(veff_3d, nst, x, xhx, n0, n1, n2, n3, n4)
    use chebyshev_module, only: rayleigh_quotient_iso
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff_3d
    integer(4), intent(in) :: nst
    real(8), dimension(n1,n2) :: x
    real(8), dimension(n3,n4) :: xhx
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(xhx) :: n3 = shape(xhx,0)
    integer :: n4
    !f2py intent(hide), depend(xhx) :: n4 = shape(xhx,1)
    call rayleigh_quotient_iso(veff_3D=veff_3d, nst=nst, x=x, xhx=xhx)
end subroutine f90wrap_rayleigh_quotient_iso

subroutine f90wrap_cal_hx_iso(veff, nst, v, hv, n0, n1, n2, n3, n4)
    use chebyshev_module, only: cal_hx_iso
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n1,n2) :: v
    real(8), intent(inout), dimension(n3,n4) :: hv
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(v) :: n1 = shape(v,0)
    integer :: n2
    !f2py intent(hide), depend(v) :: n2 = shape(v,1)
    integer :: n3
    !f2py intent(hide), depend(hv) :: n3 = shape(hv,0)
    integer :: n4
    !f2py intent(hide), depend(hv) :: n4 = shape(hv,1)
    call cal_hx_iso(veff=veff, nst=nst, V=v, HV=hv)
end subroutine f90wrap_cal_hx_iso

subroutine f90wrap_cheby_filtering_grriso(veff_3d, x_sphe, d, n0, n1, n2, n3)
    use chebyshev_module, only: cheby_filtering_grriso
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff_3d
    real(8), intent(inout), dimension(n1,n2) :: x_sphe
    real(8), intent(inout), dimension(n3) :: d
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(x_sphe) :: n1 = shape(x_sphe,0)
    integer :: n2
    !f2py intent(hide), depend(x_sphe) :: n2 = shape(x_sphe,1)
    integer :: n3
    !f2py intent(hide), depend(d) :: n3 = shape(d,0)
    call cheby_filtering_grriso(veff_3D=veff_3d, X_sphe=x_sphe, D=d)
end subroutine f90wrap_cheby_filtering_grriso

subroutine f90wrap_estupb_iso(k, veff_3d, vec, b, n0, n1)
    use chebyshev_module, only: estupb_iso
    implicit none
    
    integer(4), intent(in) :: k
    real(8), intent(in), dimension(n0) :: veff_3d
    real(8), intent(in), dimension(n1) :: vec
    real(8), intent(out) :: b
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(vec) :: n1 = shape(vec,0)
    call estupb_iso(k=k, veff_3D=veff_3d, vec=vec, b=b)
end subroutine f90wrap_estupb_iso

subroutine f90wrap_chebyshev_filter_scaled_iso(veff_3d, x, m, a, b, al, n0, n1, n2)
    use chebyshev_module, only: chebyshev_filter_scaled_iso
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff_3d
    real(8), intent(inout), dimension(n1,n2) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: al
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    call chebyshev_filter_scaled_iso(veff_3D=veff_3d, X=x, m=m, a=a, b=b, al=al)
end subroutine f90wrap_chebyshev_filter_scaled_iso

subroutine f90wrap_grayleigh_ritz_iso(veff_3d, x, d, n0, n1, n2, n3)
    use chebyshev_module, only: grayleigh_ritz_iso
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff_3d
    real(8), intent(inout), dimension(n1,n2) :: x
    real(8), intent(inout), dimension(n3) :: d
    integer :: n0
    !f2py intent(hide), depend(veff_3d) :: n0 = shape(veff_3d,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(d) :: n3 = shape(d,0)
    call grayleigh_ritz_iso(veff_3D=veff_3d, X=x, D=d)
end subroutine f90wrap_grayleigh_ritz_iso

subroutine f90wrap_first_chescf_sto_gamma(veff, nev, psi, eval, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use chebyshev_module, only: first_chescf_sto_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: veff
    integer(4), intent(in) :: nev
    real(8), dimension(n4,n5,n6,n7) :: psi
    real(8), dimension(n8,n9,n10) :: eval
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
    call first_chescf_sto_gamma(veff=veff, nev=nev, psi=psi, eval=eval)
end subroutine f90wrap_first_chescf_sto_gamma

subroutine f90wrap_cheby_init_sto_gamma(nmax, npw, initx_sto, n0, n1)
    use chebyshev_module, only: cheby_init_sto_gamma
    implicit none
    
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: npw
    real(8), intent(inout), dimension(n0,n1) :: initx_sto
    integer :: n0
    !f2py intent(hide), depend(initx_sto) :: n0 = shape(initx_sto,0)
    integer :: n1
    !f2py intent(hide), depend(initx_sto) :: n1 = shape(initx_sto,1)
    call cheby_init_sto_gamma(Nmax=nmax, Npw=npw, initX_sto=initx_sto)
end subroutine f90wrap_cheby_init_sto_gamma

subroutine f90wrap_first_subspace_stopw_gamma(ik, veff, nmax, nev, initx, x, d, n0, n1, n2, n3, n4, n5, n6, n7)
    use chebyshev_module, only: first_subspace_stopw_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nmax
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n3,n4) :: initx
    real(8), intent(inout), dimension(n5,n6) :: x
    real(8), intent(inout), dimension(n7) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(initx) :: n3 = shape(initx,0)
    integer :: n4
    !f2py intent(hide), depend(initx) :: n4 = shape(initx,1)
    integer :: n5
    !f2py intent(hide), depend(x) :: n5 = shape(x,0)
    integer :: n6
    !f2py intent(hide), depend(x) :: n6 = shape(x,1)
    integer :: n7
    !f2py intent(hide), depend(d) :: n7 = shape(d,0)
    call first_subspace_stopw_gamma(Ik=ik, veff=veff, Nmax=nmax, nev=nev, initX=initx, X=x, D=d)
end subroutine f90wrap_first_subspace_stopw_gamma

subroutine f90wrap_rayleigh_quotient_gamma(ik, veff, nst, x, xhx, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: rayleigh_quotient_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nst
    real(8), dimension(n3,n4) :: x
    real(8), dimension(n5,n6) :: xhx
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(xhx) :: n5 = shape(xhx,0)
    integer :: n6
    !f2py intent(hide), depend(xhx) :: n6 = shape(xhx,1)
    call rayleigh_quotient_gamma(Ik=ik, veff=veff, nst=nst, x=x, xhx=xhx)
end subroutine f90wrap_rayleigh_quotient_gamma

subroutine f90wrap_cal_hx_gamma(ik, veff, nst, v, hv, n0, n1, n2, n3, n4, n5, n6)
    use chebyshev_module, only: cal_hx_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n3,n4) :: v
    real(8), intent(inout), dimension(n5,n6) :: hv
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(v) :: n3 = shape(v,0)
    integer :: n4
    !f2py intent(hide), depend(v) :: n4 = shape(v,1)
    integer :: n5
    !f2py intent(hide), depend(hv) :: n5 = shape(hv,0)
    integer :: n6
    !f2py intent(hide), depend(hv) :: n6 = shape(hv,1)
    call cal_hx_gamma(Ik=ik, veff=veff, nst=nst, V=v, HV=hv)
end subroutine f90wrap_cal_hx_gamma

subroutine f90wrap_cheby_filter_rr_gamma(ik, veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: cheby_filter_rr_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call cheby_filter_rr_gamma(Ik=ik, veff=veff, X=x, D=d)
end subroutine f90wrap_cheby_filter_rr_gamma

subroutine f90wrap_estupb_gamma(k, ik, veff, vec, b, n0, n1, n2, n3)
    use chebyshev_module, only: estupb_gamma
    implicit none
    
    integer(4), intent(in) :: k
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(in), dimension(n3) :: vec
    real(8), intent(out) :: b
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(vec) :: n3 = shape(vec,0)
    call estupb_gamma(k=k, Ik=ik, veff=veff, vec=vec, b=b)
end subroutine f90wrap_estupb_gamma

subroutine f90wrap_chebyshev_filter_scaled_gamma(ik, veff, x, m, a, b, al, n0, n1, n2, n3, n4)
    use chebyshev_module, only: chebyshev_filter_scaled_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: al
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    call chebyshev_filter_scaled_gamma(Ik=ik, veff=veff, X=x, m=m, a=a, b=b, al=al)
end subroutine f90wrap_chebyshev_filter_scaled_gamma

subroutine f90wrap_grayleigh_ritz_gamma(ik, veff, x, d, n0, n1, n2, n3, n4, n5)
    use chebyshev_module, only: grayleigh_ritz_gamma
    implicit none
    
    integer(4), intent(in) :: ik
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    real(8), intent(inout), dimension(n3,n4) :: x
    real(8), intent(inout), dimension(n5) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,0)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,1)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    call grayleigh_ritz_gamma(Ik=ik, veff=veff, X=x, D=d)
end subroutine f90wrap_grayleigh_ritz_gamma

subroutine f90wrap_chebyshev_module__get__iiii(f90wrap_iiii)
    use chebyshev_module, only: chebyshev_module_iiii => iiii
    implicit none
    integer(4), intent(out) :: f90wrap_iiii
    
    f90wrap_iiii = chebyshev_module_iiii
end subroutine f90wrap_chebyshev_module__get__iiii

subroutine f90wrap_chebyshev_module__set__iiii(f90wrap_iiii)
    use chebyshev_module, only: chebyshev_module_iiii => iiii
    implicit none
    integer(4), intent(in) :: f90wrap_iiii
    
    chebyshev_module_iiii = f90wrap_iiii
end subroutine f90wrap_chebyshev_module__set__iiii

! End of module chebyshev_module defined in file Chebyshev_fliter.fpp

