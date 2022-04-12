! Module poisson_isf defined in file poisson_isf.fpp

subroutine f90wrap_karray_set
    use poisson_isf, only: karray_set
    implicit none
    
    call karray_set()
end subroutine f90wrap_karray_set

subroutine f90wrap_poisson_isf_method(rho, vh, n0, n1, n2, n3, n4, n5)
    use poisson_isf, only: poisson_isf_method
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
    call poisson_isf_method(rho=rho, vh=vh)
end subroutine f90wrap_poisson_isf_method

subroutine f90wrap_build_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, itype_scf, karrayout, n0, n1, n2)
    use poisson_isf, only: build_kernel
    implicit none
    
    integer, intent(in) :: n01
    integer, intent(in) :: n02
    integer, intent(in) :: n03
    integer, intent(in) :: nfft1
    integer, intent(in) :: nfft2
    integer, intent(in) :: nfft3
    real(8), intent(in) :: hgrid
    integer, intent(in) :: itype_scf
    real(8), intent(inout), dimension(n0,n1,n2) :: karrayout
    integer :: n0
    !f2py intent(hide), depend(karrayout) :: n0 = shape(karrayout,0)
    integer :: n1
    !f2py intent(hide), depend(karrayout) :: n1 = shape(karrayout,1)
    integer :: n2
    !f2py intent(hide), depend(karrayout) :: n2 = shape(karrayout,2)
    call build_kernel(n01=n01, n02=n02, n03=n03, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, hgrid=hgrid, itype_scf=itype_scf, &
        karrayout=karrayout)
end subroutine f90wrap_build_kernel

subroutine f90wrap_scaling_function(itype, nd, nrange, a, x, n0, n1)
    use poisson_isf, only: scaling_function
    implicit none
    
    integer, intent(in) :: itype
    integer, intent(in) :: nd
    integer, intent(out) :: nrange
    real(8), intent(inout), dimension(n0) :: a
    real(8), intent(inout), dimension(n1) :: x
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    call scaling_function(itype=itype, nd=nd, nrange=nrange, a=a, x=x)
end subroutine f90wrap_scaling_function

subroutine f90wrap_back_trans_8(nd, nt, x, y, n0, n1)
    use poisson_isf, only: back_trans_8
    implicit none
    
    integer, intent(in) :: nd
    integer, intent(in) :: nt
    real(8), intent(in), dimension(n0) :: x
    real(8), intent(inout), dimension(n1) :: y
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    call back_trans_8(nd=nd, nt=nt, x=x, y=y)
end subroutine f90wrap_back_trans_8

subroutine f90wrap_gequad(n_gauss, p_gauss, w_gauss, ur_gauss, dr_gauss, acc_gauss, n0, n1)
    use poisson_isf, only: gequad
    implicit none
    
    integer, intent(in) :: n_gauss
    real(8), intent(inout), dimension(n0) :: p_gauss
    real(8), intent(inout), dimension(n1) :: w_gauss
    real(8), intent(out) :: ur_gauss
    real(8), intent(out) :: dr_gauss
    real(8), intent(out) :: acc_gauss
    integer :: n0
    !f2py intent(hide), depend(p_gauss) :: n0 = shape(p_gauss,0)
    integer :: n1
    !f2py intent(hide), depend(w_gauss) :: n1 = shape(w_gauss,0)
    call gequad(n_gauss=n_gauss, p_gauss=p_gauss, w_gauss=w_gauss, ur_gauss=ur_gauss, dr_gauss=dr_gauss, &
        acc_gauss=acc_gauss)
end subroutine f90wrap_gequad

subroutine f90wrap_calculate_dimensions(n01, n02, n03, nfft1, nfft2, nfft3)
    use poisson_isf, only: calculate_dimensions
    implicit none
    
    integer, intent(in) :: n01
    integer, intent(in) :: n02
    integer, intent(in) :: n03
    integer, intent(out) :: nfft1
    integer, intent(out) :: nfft2
    integer, intent(out) :: nfft3
    call calculate_dimensions(n01=n01, n02=n02, n03=n03, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3)
end subroutine f90wrap_calculate_dimensions

subroutine f90wrap_psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot, n0, n1, n2, n3, n4, n5)
    use poisson_isf, only: psolver_kernel
    implicit none
    
    integer, intent(in) :: n01
    integer, intent(in) :: n02
    integer, intent(in) :: n03
    integer, intent(in) :: nfft1
    integer, intent(in) :: nfft2
    integer, intent(in) :: nfft3
    real(8), intent(in) :: hgrid
    real(8), intent(in), dimension(n0,n1,n2) :: karray
    real(8), intent(inout), dimension(n3,n4,n5) :: rhopot
    integer :: n0
    !f2py intent(hide), depend(karray) :: n0 = shape(karray,0)
    integer :: n1
    !f2py intent(hide), depend(karray) :: n1 = shape(karray,1)
    integer :: n2
    !f2py intent(hide), depend(karray) :: n2 = shape(karray,2)
    integer :: n3
    !f2py intent(hide), depend(rhopot) :: n3 = shape(rhopot,0)
    integer :: n4
    !f2py intent(hide), depend(rhopot) :: n4 = shape(rhopot,1)
    integer :: n5
    !f2py intent(hide), depend(rhopot) :: n5 = shape(rhopot,2)
    call psolver_kernel(n01=n01, n02=n02, n03=n03, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, hgrid=hgrid, karray=karray, &
        rhopot=rhopot)
end subroutine f90wrap_psolver_kernel

subroutine f90wrap_zarray_in(n01, n02, n03, nd1, nd2, nd3, density, zarray, n0, n1, n2, n3, n4, n5)
    use poisson_isf, only: zarray_in
    implicit none
    
    integer :: n01
    integer :: n02
    integer :: n03
    integer :: nd1
    integer :: nd2
    integer :: nd3
    real(8), dimension(n0,n1,n2) :: density
    real(8), dimension(2,n3,n4,n5) :: zarray
    integer :: n0
    !f2py intent(hide), depend(density) :: n0 = shape(density,0)
    integer :: n1
    !f2py intent(hide), depend(density) :: n1 = shape(density,1)
    integer :: n2
    !f2py intent(hide), depend(density) :: n2 = shape(density,2)
    integer :: n3
    !f2py intent(hide), depend(zarray) :: n3 = shape(zarray,1)
    integer :: n4
    !f2py intent(hide), depend(zarray) :: n4 = shape(zarray,2)
    integer :: n5
    !f2py intent(hide), depend(zarray) :: n5 = shape(zarray,3)
    call zarray_in(n01=n01, n02=n02, n03=n03, nd1=nd1, nd2=nd2, nd3=nd3, density=density, zarray=zarray)
end subroutine f90wrap_zarray_in

subroutine f90wrap_zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor, n0, n1, n2, n3, n4, n5)
    use poisson_isf, only: zarray_out
    implicit none
    
    integer, intent(in) :: n01
    integer, intent(in) :: n02
    integer, intent(in) :: n03
    integer, intent(in) :: nd1
    integer, intent(in) :: nd2
    integer, intent(in) :: nd3
    real(8), intent(inout), dimension(n0,n1,n2) :: rhopot
    real(8), intent(in), dimension(n3,n4,n5) :: zarray
    real(8), intent(in) :: factor
    integer :: n0
    !f2py intent(hide), depend(rhopot) :: n0 = shape(rhopot,0)
    integer :: n1
    !f2py intent(hide), depend(rhopot) :: n1 = shape(rhopot,1)
    integer :: n2
    !f2py intent(hide), depend(rhopot) :: n2 = shape(rhopot,2)
    integer :: n3
    !f2py intent(hide), depend(zarray) :: n3 = shape(zarray,0)
    integer :: n4
    !f2py intent(hide), depend(zarray) :: n4 = shape(zarray,1)
    integer :: n5
    !f2py intent(hide), depend(zarray) :: n5 = shape(zarray,2)
    call zarray_out(n01=n01, n02=n02, n03=n03, nd1=nd1, nd2=nd2, nd3=nd3, rhopot=rhopot, zarray=zarray, factor=factor)
end subroutine f90wrap_zarray_out

subroutine f90wrap_scf_recursion(itype, n_iter, n_range, kernel_scf, kern_1_scf, n0, n1)
    use poisson_isf, only: scf_recursion
    implicit none
    
    integer, intent(in) :: itype
    integer, intent(in) :: n_iter
    integer, intent(in) :: n_range
    real(8), intent(inout), dimension(n0) :: kernel_scf
    real(8), intent(inout), dimension(n1) :: kern_1_scf
    integer :: n0
    !f2py intent(hide), depend(kernel_scf) :: n0 = shape(kernel_scf,0)
    integer :: n1
    !f2py intent(hide), depend(kern_1_scf) :: n1 = shape(kern_1_scf,0)
    call scf_recursion(itype=itype, n_iter=n_iter, n_range=n_range, kernel_scf=kernel_scf, kern_1_scf=kern_1_scf)
end subroutine f90wrap_scf_recursion

subroutine f90wrap_scf_recursion_8(n_iter, n_range, kernel_scf, kern_1_scf, n0, n1)
    use poisson_isf, only: scf_recursion_8
    implicit none
    
    integer, intent(in) :: n_iter
    integer, intent(in) :: n_range
    real(8), intent(inout), dimension(n0) :: kernel_scf
    real(8), intent(inout), dimension(n1) :: kern_1_scf
    integer :: n0
    !f2py intent(hide), depend(kernel_scf) :: n0 = shape(kernel_scf,0)
    integer :: n1
    !f2py intent(hide), depend(kern_1_scf) :: n1 = shape(kern_1_scf,0)
    call scf_recursion_8(n_iter=n_iter, n_range=n_range, kernel_scf=kernel_scf, kern_1_scf=kern_1_scf)
end subroutine f90wrap_scf_recursion_8

subroutine f90wrap_karrayhalf_in(n01, n02, n03, n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, kernel, karrayhalf, &
    n0, n1, n2, n3, n4, n5)
    use poisson_isf, only: karrayhalf_in
    implicit none
    
    integer, intent(in) :: n01
    integer, intent(in) :: n02
    integer, intent(in) :: n03
    integer, intent(in) :: n1k
    integer, intent(in) :: n2k
    integer, intent(in) :: n3k
    integer, intent(in) :: nfft1
    integer, intent(in) :: nfft2
    integer, intent(in) :: nfft3
    integer, intent(in) :: nd1
    integer, intent(in) :: nd2
    integer, intent(in) :: nd3
    real(8), intent(in), dimension(n0,n1,n2) :: kernel
    real(8), intent(inout), dimension(2,n3,n4,n5) :: karrayhalf
    integer :: n0
    !f2py intent(hide), depend(kernel) :: n0 = shape(kernel,0)
    integer :: n1
    !f2py intent(hide), depend(kernel) :: n1 = shape(kernel,1)
    integer :: n2
    !f2py intent(hide), depend(kernel) :: n2 = shape(kernel,2)
    integer :: n3
    !f2py intent(hide), depend(karrayhalf) :: n3 = shape(karrayhalf,1)
    integer :: n4
    !f2py intent(hide), depend(karrayhalf) :: n4 = shape(karrayhalf,2)
    integer :: n5
    !f2py intent(hide), depend(karrayhalf) :: n5 = shape(karrayhalf,3)
    call karrayhalf_in(n01=n01, n02=n02, n03=n03, n1k=n1k, n2k=n2k, n3k=n3k, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, nd1=nd1, &
        nd2=nd2, nd3=nd3, kernel=kernel, karrayhalf=karrayhalf)
end subroutine f90wrap_karrayhalf_in

subroutine f90wrap_kernel_recon(n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, zarray, karray, n0, n1, n2, n3)
    use poisson_isf, only: kernel_recon
    implicit none
    
    integer, intent(in) :: n1k
    integer, intent(in) :: n2k
    integer, intent(in) :: n3k
    integer, intent(in) :: nfft1
    integer, intent(in) :: nfft2
    integer, intent(in) :: nfft3
    integer, intent(in) :: nd1
    integer, intent(in) :: nd2
    integer, intent(in) :: nd3
    real(8), intent(in), dimension(2,n0) :: zarray
    real(8), intent(inout), dimension(n1,n2,n3) :: karray
    integer :: n0
    !f2py intent(hide), depend(zarray) :: n0 = shape(zarray,1)
    integer :: n1
    !f2py intent(hide), depend(karray) :: n1 = shape(karray,0)
    integer :: n2
    !f2py intent(hide), depend(karray) :: n2 = shape(karray,1)
    integer :: n3
    !f2py intent(hide), depend(karray) :: n3 = shape(karray,2)
    call kernel_recon(n1k=n1k, n2k=n2k, n3k=n3k, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, nd1=nd1, nd2=nd2, nd3=nd3, &
        zarray=zarray, karray=karray)
end subroutine f90wrap_kernel_recon

subroutine f90wrap_norm_ind(nd1, nd2, nd3, i1, i2, i3, ind)
    use poisson_isf, only: norm_ind
    implicit none
    
    integer :: nd1
    integer :: nd2
    integer :: nd3
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: ind
    call norm_ind(nd1=nd1, nd2=nd2, nd3=nd3, i1=i1, i2=i2, i3=i3, ind=ind)
end subroutine f90wrap_norm_ind

subroutine f90wrap_symm_ind(nd1, nd2, nd3, i1, i2, i3, ind)
    use poisson_isf, only: symm_ind
    implicit none
    
    integer :: nd1
    integer :: nd2
    integer :: nd3
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: ind
    call symm_ind(nd1=nd1, nd2=nd2, nd3=nd3, i1=i1, i2=i2, i3=i3, ind=ind)
end subroutine f90wrap_symm_ind

subroutine f90wrap_kernel_application(n1xy, n2xy, n3xy, nd1h, nd2, nd3, nfft1, nfft2, nfft3, zarray, karray, inzee, n0, &
    n1, n2, n3, n4, n5, n6, n7)
    use poisson_isf, only: kernel_application
    implicit none
    
    integer, intent(in) :: n1xy
    integer, intent(in) :: n2xy
    integer, intent(in) :: n3xy
    integer, intent(in) :: nd1h
    integer, intent(in) :: nd2
    integer, intent(in) :: nd3
    integer, intent(in) :: nfft1
    integer, intent(in) :: nfft2
    integer, intent(in) :: nfft3
    real(8), intent(inout), dimension(n0,n1,n2,n3,n4) :: zarray
    real(8), intent(in), dimension(n5,n6,n7) :: karray
    integer, intent(in) :: inzee
    integer :: n0
    !f2py intent(hide), depend(zarray) :: n0 = shape(zarray,0)
    integer :: n1
    !f2py intent(hide), depend(zarray) :: n1 = shape(zarray,1)
    integer :: n2
    !f2py intent(hide), depend(zarray) :: n2 = shape(zarray,2)
    integer :: n3
    !f2py intent(hide), depend(zarray) :: n3 = shape(zarray,3)
    integer :: n4
    !f2py intent(hide), depend(zarray) :: n4 = shape(zarray,4)
    integer :: n5
    !f2py intent(hide), depend(karray) :: n5 = shape(karray,0)
    integer :: n6
    !f2py intent(hide), depend(karray) :: n6 = shape(karray,1)
    integer :: n7
    !f2py intent(hide), depend(karray) :: n7 = shape(karray,2)
    call kernel_application(n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, nd1h=nd1h, nd2=nd2, nd3=nd3, nfft1=nfft1, nfft2=nfft2, &
        nfft3=nfft3, zarray=zarray, karray=karray, inzee=inzee)
end subroutine f90wrap_kernel_application

subroutine f90wrap_poisson_isf__array__karray(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use sgfft_oct_m
    use grid_module
    use poisson_isf, only: poisson_isf_karray => karray
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(poisson_isf_karray)) then
        dshape(1:3) = shape(poisson_isf_karray)
        dloc = loc(poisson_isf_karray)
    else
        dloc = 0
    end if
end subroutine f90wrap_poisson_isf__array__karray

! End of module poisson_isf defined in file poisson_isf.fpp

