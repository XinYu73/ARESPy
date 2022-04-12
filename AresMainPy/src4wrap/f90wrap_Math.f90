! Module math defined in file Math.fpp

subroutine f90wrap_change_case(instr, str, fun)
    use math, only: change_case
    implicit none
    
    character*(*) :: instr
    character*(*) :: str
    integer(4) :: fun
    call change_case(instr=instr, str=str, fun=fun)
end subroutine f90wrap_change_case

subroutine f90wrap_find_keywords(str, ch_mark, id_key, id_value)
    use math, only: find_keywords
    implicit none
    
    character*(*) :: str
    character(1) :: ch_mark
    integer(4) :: id_key
    integer(4) :: id_value
    call find_keywords(str=str, ch_mark=ch_mark, id_key=id_key, id_value=id_value)
end subroutine f90wrap_find_keywords

subroutine f90wrap_find_nword(str, ch_comma, nword)
    use math, only: find_nword
    implicit none
    
    character*(*) :: str
    character(1) :: ch_comma
    integer(4) :: nword
    call find_nword(str=str, ch_comma=ch_comma, nword=nword)
end subroutine f90wrap_find_nword

subroutine f90wrap_det(ret_det, matrix)
    use math, only: det
    implicit none
    
    real(8), intent(out) :: ret_det
    real(8), intent(in), dimension(3,3) :: matrix
    ret_det = det(Matrix=matrix)
end subroutine f90wrap_det

subroutine f90wrap_inv_33(ret_inv_33, m)
    use math, only: inv_33
    implicit none
    
    real(8), dimension(3,3), intent(out) :: ret_inv_33
    real(8), intent(in), dimension(3,3) :: m
    ret_inv_33 = inv_33(M=m)
end subroutine f90wrap_inv_33

subroutine f90wrap_lindg(eta, lambda_, ret_lindg, mu)
    use math, only: lindg
    implicit none
    
    real(8), intent(in) :: eta
    real(8), intent(in) :: lambda_
    real(8), intent(out) :: ret_lindg
    real(8), intent(in) :: mu
    ret_lindg = lindg(eta=eta, lambda=lambda_, mu=mu)
end subroutine f90wrap_lindg

subroutine f90wrap_int_to_char(ret_int_to_char, int_bn)
    use math, only: int_to_char
    implicit none
    
    character(6), intent(out) :: ret_int_to_char
    integer(4), intent(in) :: int_bn
    ret_int_to_char = int_to_char(int=int_bn)
end subroutine f90wrap_int_to_char

subroutine f90wrap_lat2matrix(lat_para, lat_mat, flag)
    use math, only: lat2matrix
    implicit none
    
    real(8), dimension(6) :: lat_para
    real(8), dimension(3,3) :: lat_mat
    integer(4) :: flag
    call lat2matrix(lat_para=lat_para, lat_mat=lat_mat, flag=flag)
end subroutine f90wrap_lat2matrix

subroutine f90wrap_one2three(id, n_dens, pos)
    use math, only: one2three
    implicit none
    
    integer(4) :: id
    integer(4), dimension(3) :: n_dens
    integer(4), dimension(3) :: pos
    call one2three(id=id, n_dens=n_dens, pos=pos)
end subroutine f90wrap_one2three

subroutine f90wrap_init_random_seed
    use math, only: init_random_seed
    implicit none
    
    call init_random_seed()
end subroutine f90wrap_init_random_seed

subroutine f90wrap_atom_mass(atom_name, mass)
    use math, only: atom_mass
    implicit none
    
    character*(*) :: atom_name
    real(8) :: mass
    call atom_mass(atom_name=atom_name, mass=mass)
end subroutine f90wrap_atom_mass

subroutine f90wrap_newton_inter(n, x, y, m, tx, ty, n0, n1, n2, n3)
    use math, only: newton_inter
    implicit none
    
    integer(4) :: n
    real(8), dimension(n0) :: x
    real(8), dimension(n1) :: y
    integer :: m
    real(8), dimension(n2) :: tx
    real(8), dimension(n3) :: ty
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    integer :: n2
    !f2py intent(hide), depend(tx) :: n2 = shape(tx,0)
    integer :: n3
    !f2py intent(hide), depend(ty) :: n3 = shape(ty,0)
    call newton_inter(n=n, x=x, y=y, m=m, tx=tx, ty=ty)
end subroutine f90wrap_newton_inter

subroutine f90wrap_diag(n, ina, w, q, n0, n1, n2, n3, n4)
    use math, only: diag
    implicit none
    
    integer(4) :: n
    real(8), dimension(n0,n1) :: ina
    real(8), dimension(n2) :: w
    real(8), dimension(n3,n4) :: q
    integer :: n0
    !f2py intent(hide), depend(ina) :: n0 = shape(ina,0)
    integer :: n1
    !f2py intent(hide), depend(ina) :: n1 = shape(ina,1)
    integer :: n2
    !f2py intent(hide), depend(w) :: n2 = shape(w,0)
    integer :: n3
    !f2py intent(hide), depend(q) :: n3 = shape(q,0)
    integer :: n4
    !f2py intent(hide), depend(q) :: n4 = shape(q,1)
    call diag(N=n, inA=ina, W=w, Q=q)
end subroutine f90wrap_diag

subroutine f90wrap_boltzmann_distribution(rnull, ret_boltzmann_distribution, width)
    use math, only: boltzmann_distribution
    implicit none
    
    real(8) :: rnull
    real(8), intent(out) :: ret_boltzmann_distribution
    real(8) :: width
    ret_boltzmann_distribution = boltzmann_distribution(rnull=rnull, width=width)
end subroutine f90wrap_boltzmann_distribution

subroutine f90wrap_dir2car(cry_coo, ort_coo, lat, n0, n1, n2, n3)
    use math, only: dir2car
    implicit none
    
    real(8), dimension(n0,n1) :: cry_coo
    real(8), dimension(n2,n3) :: ort_coo
    real(8), dimension(3,3) :: lat
    integer :: n0
    !f2py intent(hide), depend(cry_coo) :: n0 = shape(cry_coo,0)
    integer :: n1
    !f2py intent(hide), depend(cry_coo) :: n1 = shape(cry_coo,1)
    integer :: n2
    !f2py intent(hide), depend(ort_coo) :: n2 = shape(ort_coo,0)
    integer :: n3
    !f2py intent(hide), depend(ort_coo) :: n3 = shape(ort_coo,1)
    call dir2car(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)
end subroutine f90wrap_dir2car

subroutine f90wrap_car2dir(ort_coo, cry_coo, lat, n0, n1, n2, n3)
    use math, only: car2dir
    implicit none
    
    real(8), dimension(n0,n1) :: ort_coo
    real(8), dimension(n2,n3) :: cry_coo
    real(8), dimension(3,3) :: lat
    integer :: n0
    !f2py intent(hide), depend(ort_coo) :: n0 = shape(ort_coo,0)
    integer :: n1
    !f2py intent(hide), depend(ort_coo) :: n1 = shape(ort_coo,1)
    integer :: n2
    !f2py intent(hide), depend(cry_coo) :: n2 = shape(cry_coo,0)
    integer :: n3
    !f2py intent(hide), depend(cry_coo) :: n3 = shape(cry_coo,1)
    call car2dir(ort_coo=ort_coo, cry_coo=cry_coo, lat=lat)
end subroutine f90wrap_car2dir

subroutine f90wrap_thr2mat(n1, n2, n3, i, j, k, dimnu)
    use math, only: thr2mat
    implicit none
    
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: i
    integer(4), intent(in) :: j
    integer(4), intent(in) :: k
    integer(4), intent(out) :: dimnu
    call thr2mat(n1=n1, n2=n2, n3=n3, I=i, J=j, K=k, dimnu=dimnu)
end subroutine f90wrap_thr2mat

subroutine f90wrap_mat2thr(n1, n2, n3, i, ix, iy, iz)
    use math, only: mat2thr
    implicit none
    
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: i
    integer(4), intent(out) :: ix
    integer(4), intent(out) :: iy
    integer(4), intent(out) :: iz
    call mat2thr(n1=n1, n2=n2, n3=n3, I=i, Ix=ix, Iy=iy, Iz=iz)
end subroutine f90wrap_mat2thr

subroutine f90wrap_sopo(a, lda, n, b, ldb, m, w, n0, n1, n2, n3, n4)
    use math, only: sopo
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: a
    integer, intent(in) :: lda
    integer, intent(in) :: n
    real(8), intent(inout), dimension(n2,n3) :: b
    integer, intent(in) :: ldb
    integer, intent(in) :: m
    real(8), intent(inout), dimension(n4) :: w
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    integer :: n2
    !f2py intent(hide), depend(b) :: n2 = shape(b,0)
    integer :: n3
    !f2py intent(hide), depend(b) :: n3 = shape(b,1)
    integer :: n4
    !f2py intent(hide), depend(w) :: n4 = shape(w,0)
    call sopo(A=a, LDA=lda, N=n, B=b, LDB=ldb, M=m, W=w)
end subroutine f90wrap_sopo

subroutine f90wrap_csort_eigen(nev, arr, brr, n0, n1, n2)
    use math, only: csort_eigen
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n0) :: arr
    complex(8), intent(inout), dimension(n1,n2) :: brr
    integer :: n0
    !f2py intent(hide), depend(arr) :: n0 = shape(arr,0)
    integer :: n1
    !f2py intent(hide), depend(brr) :: n1 = shape(brr,0)
    integer :: n2
    !f2py intent(hide), depend(brr) :: n2 = shape(brr,1)
    call csort_eigen(nev=nev, arr=arr, brr=brr)
end subroutine f90wrap_csort_eigen

subroutine f90wrap_rsort_eigen(nev, arr, brr, n0, n1, n2)
    use math, only: rsort_eigen
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n0) :: arr
    real(8), intent(inout), dimension(n1,n2) :: brr
    integer :: n0
    !f2py intent(hide), depend(arr) :: n0 = shape(arr,0)
    integer :: n1
    !f2py intent(hide), depend(brr) :: n1 = shape(brr,0)
    integer :: n2
    !f2py intent(hide), depend(brr) :: n2 = shape(brr,1)
    call rsort_eigen(nev=nev, arr=arr, brr=brr)
end subroutine f90wrap_rsort_eigen

subroutine f90wrap_sort_eigval(n, arr, n0)
    use math, only: sort_eigval
    implicit none
    
    integer(4), intent(in) :: n
    real(8), intent(inout), dimension(n0) :: arr
    integer :: n0
    !f2py intent(hide), depend(arr) :: n0 = shape(arr,0)
    call sort_eigval(n=n, arr=arr)
end subroutine f90wrap_sort_eigval

subroutine f90wrap_realint_sort(nev, arr, brr, crr, n0, n1, n2, n3)
    use math, only: realint_sort
    implicit none
    
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n0) :: arr
    integer(4), intent(inout), dimension(n1,n2) :: brr
    real(8), optional, intent(inout), dimension(3,n3) :: crr
    integer :: n0
    !f2py intent(hide), depend(arr) :: n0 = shape(arr,0)
    integer :: n1
    !f2py intent(hide), depend(brr) :: n1 = shape(brr,0)
    integer :: n2
    !f2py intent(hide), depend(brr) :: n2 = shape(brr,1)
    integer :: n3
    !f2py intent(hide), depend(crr) :: n3 = shape(crr,1)
    call realint_sort(nev=nev, arr=arr, brr=brr, crr=crr)
end subroutine f90wrap_realint_sort

subroutine f90wrap_pgfo(omat, n1, n2, n3, ix, iy, iz, ex, ey, ez)
    use math, only: pgfo
    implicit none
    
    real(8), intent(in), dimension(3,3) :: omat
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: ix
    integer(4), intent(in) :: iy
    integer(4), intent(in) :: iz
    integer(4), intent(out) :: ex
    integer(4), intent(out) :: ey
    integer(4), intent(out) :: ez
    call pgfo(Omat=omat, n1=n1, n2=n2, n3=n3, Ix=ix, Iy=iy, Iz=iz, Ex=ex, Ey=ey, Ez=ez)
end subroutine f90wrap_pgfo

subroutine f90wrap_cubicsplineinterp(fun, ddfdx2, xmax, dx, x, ret_cubicsplineinterp, zion, n0, n1)
    use math, only: cubicsplineinterp
    implicit none
    
    real(8), intent(in), dimension(n0) :: fun
    real(8), intent(in), dimension(n1) :: ddfdx2
    real(8), intent(in) :: xmax
    real(8), intent(in) :: dx
    real(8), intent(in) :: x
    real(8), intent(out) :: ret_cubicsplineinterp
    real(8), optional, intent(in) :: zion
    integer :: n0
    !f2py intent(hide), depend(fun) :: n0 = shape(fun,0)
    integer :: n1
    !f2py intent(hide), depend(ddfdx2) :: n1 = shape(ddfdx2,0)
    ret_cubicsplineinterp = cubicsplineinterp(fun=fun, ddfdx2=ddfdx2, xmax=xmax, dx=dx, x=x, Zion=zion)
end subroutine f90wrap_cubicsplineinterp

subroutine f90wrap_finite_factor(fnor, norder, coe, n0)
    use math, only: finite_factor
    implicit none
    
    integer(4), intent(in) :: fnor
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n0) :: coe
    integer :: n0
    !f2py intent(hide), depend(coe) :: n0 = shape(coe,0)
    call finite_factor(fnor=fnor, norder=norder, coe=coe)
end subroutine f90wrap_finite_factor

subroutine f90wrap_finite_factor_new(fnor, norder, coe, n0)
    use math, only: finite_factor_new
    implicit none
    
    integer(4), intent(in) :: fnor
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n0) :: coe
    integer :: n0
    !f2py intent(hide), depend(coe) :: n0 = shape(coe,0)
    call finite_factor_new(fnor=fnor, norder=norder, coe=coe)
end subroutine f90wrap_finite_factor_new

subroutine f90wrap_dfdr(np, h, f, df, n0, n1)
    use math, only: dfdr
    implicit none
    
    integer(4), intent(in) :: np
    real(8), intent(in) :: h
    real(8), intent(in), dimension(n0) :: f
    real(8), intent(inout), dimension(n1) :: df
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(df) :: n1 = shape(df,0)
    call dfdr(np=np, h=h, f=f, df=df)
end subroutine f90wrap_dfdr

subroutine f90wrap_cubichermiteinterp(fun, dfdx, xmax, h, ret_cubichermiteinterp, x, n0, n1)
    use math, only: cubichermiteinterp
    implicit none
    
    real(8), intent(in), dimension(n0) :: fun
    real(8), intent(in), dimension(n1) :: dfdx
    real(8), intent(in) :: xmax
    real(8), intent(in) :: h
    real(8), intent(out) :: ret_cubichermiteinterp
    real(8), intent(in) :: x
    integer :: n0
    !f2py intent(hide), depend(fun) :: n0 = shape(fun,0)
    integer :: n1
    !f2py intent(hide), depend(dfdx) :: n1 = shape(dfdx,0)
    ret_cubichermiteinterp = cubichermiteinterp(fun=fun, dfdx=dfdx, xmax=xmax, h=h, x=x)
end subroutine f90wrap_cubichermiteinterp

subroutine f90wrap_simpleinterp(fun, xmax, h, ret_simpleinterp, x, n0)
    use math, only: simpleinterp
    implicit none
    
    real(8), intent(in), dimension(n0) :: fun
    real(8), intent(in) :: xmax
    real(8), intent(in) :: h
    real(8), intent(out) :: ret_simpleinterp
    real(8), intent(in) :: x
    integer :: n0
    !f2py intent(hide), depend(fun) :: n0 = shape(fun,0)
    ret_simpleinterp = simpleinterp(fun=fun, xmax=xmax, h=h, x=x)
end subroutine f90wrap_simpleinterp

subroutine f90wrap_r_dylm(l, m, x, y, z, rmod, f, n0)
    use math, only: r_dylm
    implicit none
    
    integer(4), intent(in) :: l
    integer(4), intent(in) :: m
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: z
    real(8), intent(in) :: rmod
    real(8), intent(inout), dimension(n0) :: f
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    call r_dylm(l=l, m=m, x=x, y=y, z=z, rmod=rmod, f=f)
end subroutine f90wrap_r_dylm

subroutine f90wrap_atom_effcharge(atom_name, lmax, nquan, zeta)
    use math, only: atom_effcharge
    implicit none
    
    character*(*) :: atom_name
    integer, intent(out) :: lmax
    integer, dimension(4), intent(inout) :: nquan
    real(8), dimension(4), intent(inout) :: zeta
    call atom_effcharge(atom_name=atom_name, Lmax=lmax, nquan=nquan, zeta=zeta)
end subroutine f90wrap_atom_effcharge

subroutine f90wrap_atom_sto(p, l, m, zeta, r, f, n0)
    use math, only: atom_sto
    implicit none
    
    integer(4), intent(in) :: p
    integer(4), intent(in) :: l
    integer(4), intent(in) :: m
    real(8), intent(in) :: zeta
    real(8), intent(in), dimension(n0) :: r
    real(8), intent(out) :: f
    integer :: n0
    !f2py intent(hide), depend(r) :: n0 = shape(r,0)
    call atom_sto(p=p, l=l, m=m, zeta=zeta, r=r, f=f)
end subroutine f90wrap_atom_sto

subroutine f90wrap_gammp(a, ret_gammp, x)
    use math, only: gammp
    implicit none
    
    real(8) :: a
    real(8), intent(out) :: ret_gammp
    real(8) :: x
    ret_gammp = gammp(a=a, x=x)
end subroutine f90wrap_gammp

subroutine f90wrap_gser(gamser, a, x, gln)
    use math, only: gser
    implicit none
    
    real(8) :: gamser
    real(8) :: a
    real(8) :: x
    real(8) :: gln
    call gser(gamser=gamser, a=a, x=x, gln=gln)
end subroutine f90wrap_gser

subroutine f90wrap_gcf(gammcf, a, x, gln)
    use math, only: gcf
    implicit none
    
    real(8) :: gammcf
    real(8) :: a
    real(8) :: x
    real(8) :: gln
    call gcf(gammcf=gammcf, a=a, x=x, gln=gln)
end subroutine f90wrap_gcf

subroutine f90wrap_gammln(ret_gammln, xx)
    use math, only: gammln
    implicit none
    
    real(8), intent(out) :: ret_gammln
    real(8) :: xx
    ret_gammln = gammln(xx=xx)
end subroutine f90wrap_gammln

subroutine f90wrap_integral(l, dx, ret_integral, x)
    use math, only: integral
    implicit none
    
    integer(4) :: l
    real(8) :: dx
    real(8), intent(out) :: ret_integral
    real(8) :: x
    ret_integral = integral(l=l, dx=dx, x=x)
end subroutine f90wrap_integral

subroutine f90wrap_myfun(l, ret_myfun, t)
    use math, only: myfun
    implicit none
    
    integer(4) :: l
    real(8), intent(out) :: ret_myfun
    real(8) :: t
    ret_myfun = myfun(l=l, t=t)
end subroutine f90wrap_myfun

subroutine f90wrap_plgndr(l, m, ret_plgndr, x)
    use math, only: plgndr
    implicit none
    
    integer(4) :: l
    integer(4) :: m
    real(8), intent(out) :: ret_plgndr
    real(8) :: x
    ret_plgndr = plgndr(l=l, m=m, x=x)
end subroutine f90wrap_plgndr

subroutine f90wrap_lagrange_interpolation_coe(npoint, scatter_x, coe, n0, n1)
    use math, only: lagrange_interpolation_coe
    implicit none
    
    integer(4) :: npoint
    real(8), dimension(n0) :: scatter_x
    real(8), dimension(n1) :: coe
    integer :: n0
    !f2py intent(hide), depend(scatter_x) :: n0 = shape(scatter_x,0)
    integer :: n1
    !f2py intent(hide), depend(coe) :: n1 = shape(coe,0)
    call lagrange_interpolation_coe(npoint=npoint, scatter_x=scatter_x, coe=coe)
end subroutine f90wrap_lagrange_interpolation_coe

subroutine f90wrap_lagrange_interpolation_x(npoint, x_sample, ret_lagrange_interpolation_x, x_in, n0, n1)
    use math, only: lagrange_interpolation_x
    implicit none
    
    integer(4), intent(in) :: npoint
    real(8), intent(in), dimension(n0) :: x_sample
    real(8), intent(out), dimension(n1) :: ret_lagrange_interpolation_x
    real(8), intent(in) :: x_in
    integer :: n0
    !f2py intent(hide), depend(x_sample) :: n0 = shape(x_sample,0)
    integer :: n1
    ret_lagrange_interpolation_x = lagrange_interpolation_x(npoint=npoint, x_sample=x_sample, x_in=x_in)
end subroutine f90wrap_lagrange_interpolation_x

subroutine f90wrap_interpolation_test
    use math, only: interpolation_test
    implicit none
    
    call interpolation_test()
end subroutine f90wrap_interpolation_test

subroutine f90wrap_direct_productlm(nll, nml, index_ll, index_ml, mat_in, mat_out, n0, n1, n2, n3, n4, n5)
    use math, only: direct_productlm
    implicit none
    
    integer(4), intent(in) :: nll
    integer(4), intent(in) :: nml
    integer(4), intent(in), dimension(n0) :: index_ll
    integer(4), intent(in), dimension(n1) :: index_ml
    real(8), intent(in), dimension(n2,n3) :: mat_in
    real(8), intent(inout), dimension(n4,n5) :: mat_out
    integer :: n0
    !f2py intent(hide), depend(index_ll) :: n0 = shape(index_ll,0)
    integer :: n1
    !f2py intent(hide), depend(index_ml) :: n1 = shape(index_ml,0)
    integer :: n2
    !f2py intent(hide), depend(mat_in) :: n2 = shape(mat_in,0)
    integer :: n3
    !f2py intent(hide), depend(mat_in) :: n3 = shape(mat_in,1)
    integer :: n4
    !f2py intent(hide), depend(mat_out) :: n4 = shape(mat_out,0)
    integer :: n5
    !f2py intent(hide), depend(mat_out) :: n5 = shape(mat_out,1)
    call direct_productlm(nll=nll, nml=nml, index_ll=index_ll, index_ml=index_ml, mat_in=mat_in, mat_out=mat_out)
end subroutine f90wrap_direct_productlm

subroutine f90wrap_fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt, n0, n1, n2, n3, n4)
    use math, only: fourier_1d
    implicit none
    
    integer, intent(in) :: nr
    real(8), intent(in), dimension(n0) :: rr
    real(8), intent(in), dimension(n1) :: rab
    real(8), intent(in), dimension(n2) :: vr
    integer, intent(in) :: ll
    integer, intent(in) :: nql
    real(8), intent(in), dimension(n3) :: yp
    real(8), intent(inout), dimension(n4) :: vql
    real(8), intent(in) :: vt
    integer :: n0
    !f2py intent(hide), depend(rr) :: n0 = shape(rr,0)
    integer :: n1
    !f2py intent(hide), depend(rab) :: n1 = shape(rab,0)
    integer :: n2
    !f2py intent(hide), depend(vr) :: n2 = shape(vr,0)
    integer :: n3
    !f2py intent(hide), depend(yp) :: n3 = shape(yp,0)
    integer :: n4
    !f2py intent(hide), depend(vql) :: n4 = shape(vql,0)
    call fourier_1d(nr=nr, rr=rr, rab=rab, vr=vr, ll=ll, nql=nql, yp=yp, vql=vql, vt=vt)
end subroutine f90wrap_fourier_1d

subroutine f90wrap_invfourier_1d(g, fg, ll, r, fr, n0, n1, n2, n3)
    use math, only: invfourier_1d
    implicit none
    
    real(8), intent(in), dimension(n0) :: g
    real(8), intent(in), dimension(n1) :: fg
    integer(4), intent(in) :: ll
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
    call invfourier_1d(g=g, fg=fg, ll=ll, r=r, fr=fr)
end subroutine f90wrap_invfourier_1d

subroutine f90wrap_integ_new(rab, y, f, n0, n1)
    use math, only: integ_new
    implicit none
    
    real(8), intent(in), dimension(n0) :: rab
    real(8), intent(in), dimension(n1) :: y
    real(8), intent(out) :: f
    integer :: n0
    !f2py intent(hide), depend(rab) :: n0 = shape(rab,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    call integ_new(rab=rab, y=y, f=f)
end subroutine f90wrap_integ_new

subroutine f90wrap_interp(np, f, r, rnorm, ret_interp, z, n0, n1)
    use math, only: interp
    implicit none
    
    integer(4), intent(in) :: np
    real(8), intent(in), dimension(n0) :: f
    real(8), intent(in), dimension(n1) :: r
    real(8), intent(in) :: rnorm
    real(8), intent(out) :: ret_interp
    real(8), optional :: z
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(r) :: n1 = shape(r,0)
    ret_interp = interp(np=np, f=f, r=r, rnorm=rnorm, Z=z)
end subroutine f90wrap_interp

subroutine f90wrap_getfileunit(ret_getfileunit)
    use math, only: getfileunit
    implicit none
    
    integer(4), intent(out) :: ret_getfileunit
    ret_getfileunit = getfileunit()
end subroutine f90wrap_getfileunit

subroutine f90wrap_kahan_sum(n, ret_kahan_sum, array, n0)
    use math, only: kahan_sum
    implicit none
    
    integer(4) :: n
    real(8), intent(out) :: ret_kahan_sum
    real(8), dimension(n0) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    ret_kahan_sum = kahan_sum(n=n, array=array)
end subroutine f90wrap_kahan_sum

subroutine f90wrap_gasdev_s_sp(harvest, inmu, insigma)
    use math, only: gasdev
    implicit none
    
    real(8), intent(out) :: harvest
    real(8), optional, intent(in) :: inmu
    real(8), optional, intent(in) :: insigma
    call gasdev(harvest=harvest, inmu=inmu, insigma=insigma)
end subroutine f90wrap_gasdev_s_sp

subroutine f90wrap_gasdev_s_dp(harvest, inmu, insigma)
    use math, only: gasdev
    implicit none
    
    real(8), intent(out) :: harvest
    real(8), optional, intent(in) :: inmu
    real(8), optional, intent(in) :: insigma
    call gasdev(harvest=harvest, inmu=inmu, insigma=insigma)
end subroutine f90wrap_gasdev_s_dp

subroutine f90wrap_gasdev_v_sp(harvest, inmu, insigma, n0)
    use math, only: gasdev
    implicit none
    
    real(8), intent(inout), dimension(n0) :: harvest
    real(8), optional, intent(in) :: inmu
    real(8), optional, intent(in) :: insigma
    integer :: n0
    !f2py intent(hide), depend(harvest) :: n0 = shape(harvest,0)
    call gasdev(harvest=harvest, inmu=inmu, insigma=insigma)
end subroutine f90wrap_gasdev_v_sp

subroutine f90wrap_gasdev_v_dp(harvest, inmu, insigma, n0)
    use math, only: gasdev
    implicit none
    
    real(8), intent(inout), dimension(n0) :: harvest
    real(8), optional, intent(in) :: inmu
    real(8), optional, intent(in) :: insigma
    integer :: n0
    !f2py intent(hide), depend(harvest) :: n0 = shape(harvest,0)
    call gasdev(harvest=harvest, inmu=inmu, insigma=insigma)
end subroutine f90wrap_gasdev_v_dp

subroutine f90wrap_norm_real(ret_norm_real, a)
    use math, only: norm
    implicit none
    
    real(8), intent(out) :: ret_norm_real
    real(8), intent(in), dimension(3) :: a
    ret_norm_real = norm(a=a)
end subroutine f90wrap_norm_real

subroutine f90wrap_norm_complex(ret_norm_complex, a)
    use math, only: norm
    implicit none
    
    real(8), intent(out) :: ret_norm_complex
    complex(8), intent(in), dimension(3) :: a
    ret_norm_complex = norm(a=a)
end subroutine f90wrap_norm_complex

subroutine f90wrap_cross_real(a, ret_cross_real, b)
    use math, only: cross
    implicit none
    
    real(8), intent(in), dimension(3) :: a
    real(8), dimension(3), intent(out) :: ret_cross_real
    real(8), intent(in), dimension(3) :: b
    ret_cross_real = cross(a=a, b=b)
end subroutine f90wrap_cross_real

subroutine f90wrap_cross_complex(a, ret_cross_complex, b)
    use math, only: cross
    implicit none
    
    complex(8), intent(in), dimension(3) :: a
    complex(8), dimension(3), intent(out) :: ret_cross_complex
    complex(8), intent(in), dimension(3) :: b
    ret_cross_complex = cross(a=a, b=b)
end subroutine f90wrap_cross_complex

subroutine f90wrap_integer_index(array, n, id, n0, n1)
    use math, only: sort_id
    implicit none
    
    integer(4), dimension(n0) :: array
    integer(4) :: n
    integer(4), dimension(n1) :: id
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(id) :: n1 = shape(id,0)
    call sort_id(array=array, n=n, id=id)
end subroutine f90wrap_integer_index

subroutine f90wrap_real_index(array, n, id, n0, n1)
    use math, only: sort_id
    implicit none
    
    real(8), dimension(n0) :: array
    integer(4) :: n
    integer(4), dimension(n1) :: id
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(id) :: n1 = shape(id,0)
    call sort_id(array=array, n=n, id=id)
end subroutine f90wrap_real_index

subroutine f90wrap_three2one_d_cplx(amat, bmat, n0, n1, n2, n3)
    use math, only: three2one_dim
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: amat
    complex(8), intent(inout), dimension(n3) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    call three2one_dim(amat=amat, bmat=bmat)
end subroutine f90wrap_three2one_d_cplx

subroutine f90wrap_three2one_d_real(amat, bmat, n0, n1, n2, n3)
    use math, only: three2one_dim
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: amat
    real(8), intent(inout), dimension(n3) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    call three2one_dim(amat=amat, bmat=bmat)
end subroutine f90wrap_three2one_d_real

subroutine f90wrap_one2three_d_cplx(amat, bmat, n0, n1, n2, n3)
    use math, only: one2three_dim
    implicit none
    
    complex(8), intent(in), dimension(n0) :: amat
    complex(8), intent(inout), dimension(n1,n2,n3) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(bmat) :: n1 = shape(bmat,0)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,1)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,2)
    call one2three_dim(amat=amat, bmat=bmat)
end subroutine f90wrap_one2three_d_cplx

subroutine f90wrap_one2three_d_real(amat, bmat, n0, n1, n2, n3)
    use math, only: one2three_dim
    implicit none
    
    real(8), intent(in), dimension(n0) :: amat
    real(8), intent(inout), dimension(n1,n2,n3) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(bmat) :: n1 = shape(bmat,0)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,1)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,2)
    call one2three_dim(amat=amat, bmat=bmat)
end subroutine f90wrap_one2three_d_real

! End of module math defined in file Math.fpp

