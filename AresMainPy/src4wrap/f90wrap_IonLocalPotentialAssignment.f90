! Module getvlocalpseudopotential defined in file IonLocalPotentialAssignment.fpp

subroutine f90wrap_calvlpp
    use getvlocalpseudopotential, only: calvlpp
    implicit none
    
    call calvlpp()
end subroutine f90wrap_calvlpp

subroutine f90wrap_ionpotentialassignment(ity, zion, poscar, temp, n0)
    use getvlocalpseudopotential, only: ionpotentialassignment
    implicit none
    
    integer(4), intent(in) :: ity
    real(8), intent(in) :: zion
    real(8), dimension(3) :: poscar
    real(8), intent(inout), dimension(n0) :: temp
    integer :: n0
    !f2py intent(hide), depend(temp) :: n0 = shape(temp,0)
    call ionpotentialassignment(Ity=ity, Zion=zion, poscar=poscar, temp=temp)
end subroutine f90wrap_ionpotentialassignment

subroutine f90wrap_sphbess(l, ret_sphbess, x)
    use getvlocalpseudopotential, only: sphbess
    implicit none
    
    integer, intent(in) :: l
    real(8), intent(out) :: ret_sphbess
    real(8), intent(in) :: x
    ret_sphbess = sphbess(l=l, x=x)
end subroutine f90wrap_sphbess

subroutine f90wrap_fourbess_gr(g, fg, r, fr, n0, n1, n2, n3)
    use getvlocalpseudopotential, only: fourbess_gr
    implicit none
    
    real(8), intent(in), dimension(n0) :: g
    real(8), intent(in), dimension(n1) :: fg
    real(8), intent(in), dimension(n2) :: r
    real(8), intent(inout), dimension(n3) :: fr
    integer :: n0
    !f2py intent(hide), depend(g) :: n0 = shape(g,0)
    integer :: n1
    !f2py intent(hide), depend(fg) :: n1 = shape(fg,0)
    integer :: n2
    !f2py intent(hide), depend(r) :: n2 = shape(r,0)
    integer :: n3
    !f2py intent(hide), depend(fr) :: n3 = shape(fr,0)
    call fourbess_gr(g=g, fg=fg, r=r, fr=fr)
end subroutine f90wrap_fourbess_gr

subroutine f90wrap_cubicsplineinterp(fun, ddfdx2, xmax, dx, x, ret_cubicsplineinterp, zion, n0, n1)
    use getvlocalpseudopotential, only: cubicsplineinterp
    implicit none
    
    real(8), intent(in), dimension(n0) :: fun
    real(8), intent(in), dimension(n1) :: ddfdx2
    real(8), intent(in) :: xmax
    real(8), intent(in) :: dx
    real(8), intent(in) :: x
    real(8), intent(out) :: ret_cubicsplineinterp
    integer(4), optional, intent(in) :: zion
    integer :: n0
    !f2py intent(hide), depend(fun) :: n0 = shape(fun,0)
    integer :: n1
    !f2py intent(hide), depend(ddfdx2) :: n1 = shape(ddfdx2,0)
    ret_cubicsplineinterp = cubicsplineinterp(fun=fun, ddfdx2=ddfdx2, xmax=xmax, dx=dx, x=x, Zion=zion)
end subroutine f90wrap_cubicsplineinterp

subroutine f90wrap_cubichermiteinterp(fun, dfdx, xmax, h, x, ret_cubichermiteinterp, zion, n0, n1)
    use getvlocalpseudopotential, only: cubichermiteinterp
    implicit none
    
    real(8), intent(in), dimension(n0) :: fun
    real(8), intent(in), dimension(n1) :: dfdx
    real(8), intent(in) :: xmax
    real(8), intent(in) :: h
    real(8), intent(in) :: x
    real(8), intent(out) :: ret_cubichermiteinterp
    integer(4) :: zion
    integer :: n0
    !f2py intent(hide), depend(fun) :: n0 = shape(fun,0)
    integer :: n1
    !f2py intent(hide), depend(dfdx) :: n1 = shape(dfdx,0)
    ret_cubichermiteinterp = cubichermiteinterp(fun=fun, dfdx=dfdx, xmax=xmax, h=h, x=x, zion=zion)
end subroutine f90wrap_cubichermiteinterp

subroutine f90wrap_dfdr(np, h, f, zion, df, n0, n1)
    use getvlocalpseudopotential, only: dfdr
    implicit none
    
    integer(4), intent(in) :: np
    real(8), intent(in) :: h
    real(8), intent(in), dimension(n0) :: f
    integer(4), intent(in) :: zion
    real(8), intent(inout), dimension(n1) :: df
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(df) :: n1 = shape(df,0)
    call dfdr(np=np, h=h, f=f, zion=zion, df=df)
end subroutine f90wrap_dfdr

subroutine f90wrap_finite_factor(fnor, norder, coe, n0)
    use getvlocalpseudopotential, only: finite_factor
    implicit none
    
    integer(4), intent(in) :: fnor
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n0) :: coe
    integer :: n0
    !f2py intent(hide), depend(coe) :: n0 = shape(coe,0)
    call finite_factor(fnor=fnor, norder=norder, coe=coe)
end subroutine f90wrap_finite_factor

subroutine f90wrap_dir2car_single(cry_coo, ort_coo, lat, n0, n1)
    use getvlocalpseudopotential, only: dir2car_single
    implicit none
    
    real(8), dimension(n0) :: cry_coo
    real(8), dimension(n1) :: ort_coo
    real(8), dimension(3,3) :: lat
    integer :: n0
    !f2py intent(hide), depend(cry_coo) :: n0 = shape(cry_coo,0)
    integer :: n1
    !f2py intent(hide), depend(ort_coo) :: n1 = shape(ort_coo,0)
    call dir2car_single(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)
end subroutine f90wrap_dir2car_single

subroutine f90wrap_polynom(m, np, xa, ya, c, ret_polynom, x, n0, n1, n2)
    use getvlocalpseudopotential, only: polynom
    implicit none
    
    integer(4), intent(in) :: m
    integer(4), intent(in) :: np
    real(8), intent(in), dimension(n0) :: xa
    real(8), intent(in), dimension(n1) :: ya
    real(8), intent(inout), dimension(n2) :: c
    real(8), intent(out) :: ret_polynom
    real(8), intent(in) :: x
    integer :: n0
    !f2py intent(hide), depend(xa) :: n0 = shape(xa,0)
    integer :: n1
    !f2py intent(hide), depend(ya) :: n1 = shape(ya,0)
    integer :: n2
    !f2py intent(hide), depend(c) :: n2 = shape(c,0)
    ret_polynom = polynom(m=m, np=np, xa=xa, ya=ya, c=c, x=x)
end subroutine f90wrap_polynom

subroutine f90wrap_ionpotentialassignment_dg(ity, zion, poscar, temp, n0)
    use getvlocalpseudopotential, only: ionpotentialassignment_dg
    implicit none
    
    integer(4), intent(in) :: ity
    real(8), intent(in) :: zion
    real(8), dimension(3) :: poscar
    real(8), intent(inout), dimension(n0) :: temp
    integer :: n0
    !f2py intent(hide), depend(temp) :: n0 = shape(temp,0)
    call ionpotentialassignment_dg(Ity=ity, Zion=zion, poscar=poscar, temp=temp)
end subroutine f90wrap_ionpotentialassignment_dg

subroutine f90wrap_set_vloc_dg(ity, xyz, vloc)
    use getvlocalpseudopotential, only: set_vloc_dg
    implicit none
    
    integer(4), intent(in) :: ity
    real(8), intent(in), dimension(3) :: xyz
    real(8), intent(out) :: vloc
    call set_vloc_dg(Ity=ity, xyz=xyz, Vloc=vloc)
end subroutine f90wrap_set_vloc_dg

subroutine f90wrap_set_vcomp(n_inp, zion, rgauss, r_inp, v_inp, vcomp, n0, n1, n2)
    use getvlocalpseudopotential, only: set_vcomp
    implicit none
    
    integer(4), intent(in) :: n_inp
    integer(4), intent(in) :: zion
    real(8), intent(in) :: rgauss
    real(8), intent(in), dimension(n0) :: r_inp
    real(8), intent(in), dimension(n1) :: v_inp
    real(8), intent(inout), dimension(n2) :: vcomp
    integer :: n0
    !f2py intent(hide), depend(r_inp) :: n0 = shape(r_inp,0)
    integer :: n1
    !f2py intent(hide), depend(v_inp) :: n1 = shape(v_inp,0)
    integer :: n2
    !f2py intent(hide), depend(vcomp) :: n2 = shape(vcomp,0)
    call set_vcomp(n_inp=n_inp, zion=zion, rgauss=rgauss, r_inp=r_inp, V_inp=v_inp, Vcomp=vcomp)
end subroutine f90wrap_set_vcomp

subroutine f90wrap_getvlocalpseudopotential__get__RadiusN(f90wrap_RadiusN)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_RadiusN => RadiusN
    implicit none
    integer(4), intent(out) :: f90wrap_RadiusN
    
    f90wrap_RadiusN = getvlocalpseudopotential_RadiusN
end subroutine f90wrap_getvlocalpseudopotential__get__RadiusN

subroutine f90wrap_getvlocalpseudopotential__set__RadiusN(f90wrap_RadiusN)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_RadiusN => RadiusN
    implicit none
    integer(4), intent(in) :: f90wrap_RadiusN
    
    getvlocalpseudopotential_RadiusN = f90wrap_RadiusN
end subroutine f90wrap_getvlocalpseudopotential__set__RadiusN

subroutine f90wrap_getvlocalpseudopotential__array__fr(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use pspot_module
    use nlpot_module
    use smpi_math_module
    use getvlocalpseudopotential, only: getvlocalpseudopotential_fr => fr
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(getvlocalpseudopotential_fr)) then
        dshape(1:1) = shape(getvlocalpseudopotential_fr)
        dloc = loc(getvlocalpseudopotential_fr)
    else
        dloc = 0
    end if
end subroutine f90wrap_getvlocalpseudopotential__array__fr

subroutine f90wrap_getvlocalpseudopotential__array__wijk(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use pspot_module
    use nlpot_module
    use smpi_math_module
    use getvlocalpseudopotential, only: getvlocalpseudopotential_wijk => wijk
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(getvlocalpseudopotential_wijk)) then
        dshape(1:1) = shape(getvlocalpseudopotential_wijk)
        dloc = loc(getvlocalpseudopotential_wijk)
    else
        dloc = 0
    end if
end subroutine f90wrap_getvlocalpseudopotential__array__wijk

subroutine f90wrap_getvlocalpseudopotential__array__drijk(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use grid_module
    use pspot_module
    use nlpot_module
    use smpi_math_module
    use getvlocalpseudopotential, only: getvlocalpseudopotential_drijk => drijk
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(getvlocalpseudopotential_drijk)) then
        dshape(1:2) = shape(getvlocalpseudopotential_drijk)
        dloc = loc(getvlocalpseudopotential_drijk)
    else
        dloc = 0
    end if
end subroutine f90wrap_getvlocalpseudopotential__array__drijk

subroutine f90wrap_getvlocalpseudopotential__get__Ilft(f90wrap_Ilft)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_Ilft => Ilft
    implicit none
    integer(4), intent(out) :: f90wrap_Ilft
    
    f90wrap_Ilft = getvlocalpseudopotential_Ilft
end subroutine f90wrap_getvlocalpseudopotential__get__Ilft

subroutine f90wrap_getvlocalpseudopotential__set__Ilft(f90wrap_Ilft)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_Ilft => Ilft
    implicit none
    integer(4), intent(in) :: f90wrap_Ilft
    
    getvlocalpseudopotential_Ilft = f90wrap_Ilft
end subroutine f90wrap_getvlocalpseudopotential__set__Ilft

subroutine f90wrap_getvlocalpseudopotential__get__Irit(f90wrap_Irit)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_Irit => Irit
    implicit none
    integer(4), intent(out) :: f90wrap_Irit
    
    f90wrap_Irit = getvlocalpseudopotential_Irit
end subroutine f90wrap_getvlocalpseudopotential__get__Irit

subroutine f90wrap_getvlocalpseudopotential__set__Irit(f90wrap_Irit)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_Irit => Irit
    implicit none
    integer(4), intent(in) :: f90wrap_Irit
    
    getvlocalpseudopotential_Irit = f90wrap_Irit
end subroutine f90wrap_getvlocalpseudopotential__set__Irit

subroutine f90wrap_getvlocalpseudopotential__get__numdg(f90wrap_numdg)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_numdg => numdg
    implicit none
    integer(4), intent(out) :: f90wrap_numdg
    
    f90wrap_numdg = getvlocalpseudopotential_numdg
end subroutine f90wrap_getvlocalpseudopotential__get__numdg

subroutine f90wrap_getvlocalpseudopotential__set__numdg(f90wrap_numdg)
    use getvlocalpseudopotential, only: getvlocalpseudopotential_numdg => numdg
    implicit none
    integer(4), intent(in) :: f90wrap_numdg
    
    getvlocalpseudopotential_numdg = f90wrap_numdg
end subroutine f90wrap_getvlocalpseudopotential__set__numdg

! End of module getvlocalpseudopotential defined in file IonLocalPotentialAssignment.fpp

