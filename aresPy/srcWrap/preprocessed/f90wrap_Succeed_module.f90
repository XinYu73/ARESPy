! Module succeed defined in file Succeed_module.fpp

subroutine f90wrap_init_succeed_rho(n_rho, n_r, n_s, nspin, dv)
    use succeed, only: init_succeed_rho
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: n_r
    integer(4) :: n_s
    integer(4) :: nspin
    real(8) :: dv
    call init_succeed_rho(n_rho=n_rho, n_r=n_r, n_s=n_s, nspin=nspin, dv=dv)
end subroutine f90wrap_init_succeed_rho

subroutine f90wrap_destory_succeed
    use succeed, only: destory_succeed
    implicit none
    
    call destory_succeed()
end subroutine f90wrap_destory_succeed

subroutine f90wrap_store_rho(n_rho, nspin, rho, n0, n1)
    use succeed, only: store_rho
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: nspin
    real(8), dimension(n0,n1) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    call store_rho(n_rho=n_rho, Nspin=nspin, rho=rho)
end subroutine f90wrap_store_rho

subroutine f90wrap_store_rho_at(n_rho, nspin, rho_in, n0, n1)
    use succeed, only: store_rho_at
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: nspin
    real(8), dimension(n0,n1) :: rho_in
    integer :: n0
    !f2py intent(hide), depend(rho_in) :: n0 = shape(rho_in,0)
    integer :: n1
    !f2py intent(hide), depend(rho_in) :: n1 = shape(rho_in,1)
    call store_rho_at(n_rho=n_rho, Nspin=nspin, rho_in=rho_in)
end subroutine f90wrap_store_rho_at

subroutine f90wrap_store_r(nr, r, n0)
    use succeed, only: store_r
    implicit none
    
    integer(4) :: nr
    real(8), dimension(3,n0) :: r
    integer :: n0
    !f2py intent(hide), depend(r) :: n0 = shape(r,1)
    call store_r(nr=nr, r=r)
end subroutine f90wrap_store_r

subroutine f90wrap_store_psi(n_rho, n_s, nspin, psi, n0, n1, n2)
    use succeed, only: store_psi
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: n_s
    integer(4) :: nspin
    real(8), dimension(n0,n1,n2) :: psi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    call store_psi(n_rho=n_rho, n_s=n_s, nspin=nspin, psi=psi)
end subroutine f90wrap_store_psi

subroutine f90wrap_get_rho(nr, r_new, nrho, nspin, rho_new, n0, n1, n2)
    use succeed, only: get_rho
    implicit none
    
    integer(4) :: nr
    real(8), intent(in), dimension(3,n0) :: r_new
    integer(4) :: nrho
    integer(4) :: nspin
    real(8), intent(inout), dimension(n1,n2) :: rho_new
    integer :: n0
    !f2py intent(hide), depend(r_new) :: n0 = shape(r_new,1)
    integer :: n1
    !f2py intent(hide), depend(rho_new) :: n1 = shape(rho_new,0)
    integer :: n2
    !f2py intent(hide), depend(rho_new) :: n2 = shape(rho_new,1)
    call get_rho(nr=nr, r_new=r_new, nrho=nrho, Nspin=nspin, rho_new=rho_new)
end subroutine f90wrap_get_rho

subroutine f90wrap_get_psi(nrho, n_s, nspin, psi_new, n0, n1, n2)
    use succeed, only: get_psi
    implicit none
    
    integer(4), intent(in) :: nrho
    integer(4), intent(in) :: n_s
    integer(4), intent(in) :: nspin
    real(8), intent(inout), dimension(n0,n1,n2) :: psi_new
    integer :: n0
    !f2py intent(hide), depend(psi_new) :: n0 = shape(psi_new,0)
    integer :: n1
    !f2py intent(hide), depend(psi_new) :: n1 = shape(psi_new,1)
    integer :: n2
    !f2py intent(hide), depend(psi_new) :: n2 = shape(psi_new,2)
    call get_psi(nrho=nrho, n_s=n_s, nspin=nspin, psi_new=psi_new)
end subroutine f90wrap_get_psi

subroutine f90wrap_cal_trans_phase(nr, nspin, r_new, n1, n2, n3, ng1, ng2, ng3, gvec, trans_phase, n0, n4, &
    n5)
    use succeed, only: cal_trans_phase
    implicit none
    
    integer(4), intent(in) :: nr
    integer(4), intent(in) :: nspin
    real(8), intent(in), dimension(3,n0) :: r_new
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: ng1
    integer(4), intent(in) :: ng2
    integer(4), intent(in) :: ng3
    real(8), intent(in), dimension(4,n1) :: gvec
    complex(8), intent(inout), dimension(n2,n3,n4,n5) :: trans_phase
    integer :: n0
    !f2py intent(hide), depend(r_new) :: n0 = shape(r_new,1)
    integer :: n4
    !f2py intent(hide), depend(trans_phase) :: n4 = shape(trans_phase,2)
    integer :: n5
    !f2py intent(hide), depend(trans_phase) :: n5 = shape(trans_phase,3)
    call cal_trans_phase(nr=nr, nspin=nspin, r_new=r_new, n1=n1, n2=n2, n3=n3, ng1=ng1, ng2=ng2, ng3=ng3, gvec=gvec, &
        trans_phase=trans_phase)
end subroutine f90wrap_cal_trans_phase

subroutine f90wrap_get_new_rho_psi(nr, r_new, nrho, n1, n2, n3, nspin, rho_new, n_s, psi_new, gvec, n0, n4, &
    n5, n6)
    use succeed, only: get_new_rho_psi
    implicit none
    
    integer(4), intent(in) :: nr
    real(8), intent(in), dimension(3,n0) :: r_new
    integer(4), intent(in) :: nrho
    integer(4), intent(in) :: n1
    integer(4), intent(in) :: n2
    integer(4), intent(in) :: n3
    integer(4), intent(in) :: nspin
    real(8), intent(inout), dimension(n1,n2) :: rho_new
    integer(4), intent(in) :: n_s
    real(8), intent(inout), dimension(n3,n4,n5) :: psi_new
    real(8), intent(in), dimension(4,n6) :: gvec
    integer :: n0
    !f2py intent(hide), depend(r_new) :: n0 = shape(r_new,1)
    integer :: n4
    !f2py intent(hide), depend(psi_new) :: n4 = shape(psi_new,1)
    integer :: n5
    !f2py intent(hide), depend(psi_new) :: n5 = shape(psi_new,2)
    integer :: n6
    !f2py intent(hide), depend(gvec) :: n6 = shape(gvec,1)
    call get_new_rho_psi(nr=nr, r_new=r_new, nrho=nrho, n1=n1, n2=n2, n3=n3, Nspin=nspin, rho_new=rho_new, n_s=n_s, &
        psi_new=psi_new, gvec=gvec)
end subroutine f90wrap_get_new_rho_psi

subroutine f90wrap_store_rho_fft_trans(n_rho, nspin, rho, n0, n1)
    use succeed, only: store_rho_fft_trans
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: nspin
    real(8), dimension(n0,n1) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    call store_rho_fft_trans(n_rho=n_rho, Nspin=nspin, rho=rho)
end subroutine f90wrap_store_rho_fft_trans

subroutine f90wrap_store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2, n0, n1, n2, n3, n4)
    use succeed, only: store_rho_at_fft_trans
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: nspin
    integer(4) :: na
    real(8), dimension(n0,n1) :: rho_in
    real(8), dimension(n2,n3,n4) :: rho_in2
    integer :: n0
    !f2py intent(hide), depend(rho_in) :: n0 = shape(rho_in,0)
    integer :: n1
    !f2py intent(hide), depend(rho_in) :: n1 = shape(rho_in,1)
    integer :: n2
    !f2py intent(hide), depend(rho_in2) :: n2 = shape(rho_in2,0)
    integer :: n3
    !f2py intent(hide), depend(rho_in2) :: n3 = shape(rho_in2,1)
    integer :: n4
    !f2py intent(hide), depend(rho_in2) :: n4 = shape(rho_in2,2)
    call store_rho_at_fft_trans(n_rho=n_rho, Nspin=nspin, na=na, rho_in=rho_in, rho_in2=rho_in2)
end subroutine f90wrap_store_rho_at_fft_trans

subroutine f90wrap_store_r_fft_trans(nr, r, n0)
    use succeed, only: store_r_fft_trans
    implicit none
    
    integer(4) :: nr
    real(8), dimension(3,n0) :: r
    integer :: n0
    !f2py intent(hide), depend(r) :: n0 = shape(r,1)
    call store_r_fft_trans(nr=nr, r=r)
end subroutine f90wrap_store_r_fft_trans

subroutine f90wrap_store_psi_fft_trans(n_rho, n_s, nspin, psi, n0, n1, n2)
    use succeed, only: store_psi_fft_trans
    implicit none
    
    integer(4) :: n_rho
    integer(4) :: n_s
    integer(4) :: nspin
    real(8), dimension(n0,n1,n2) :: psi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    call store_psi_fft_trans(n_rho=n_rho, n_s=n_s, nspin=nspin, psi=psi)
end subroutine f90wrap_store_psi_fft_trans

subroutine f90wrap_succeed__get__Llastrho(f90wrap_Llastrho)
    use succeed, only: succeed_Llastrho => Llastrho
    implicit none
    logical, intent(out) :: f90wrap_Llastrho
    
    f90wrap_Llastrho = succeed_Llastrho
end subroutine f90wrap_succeed__get__Llastrho

subroutine f90wrap_succeed__set__Llastrho(f90wrap_Llastrho)
    use succeed, only: succeed_Llastrho => Llastrho
    implicit none
    logical, intent(in) :: f90wrap_Llastrho
    
    succeed_Llastrho = f90wrap_Llastrho
end subroutine f90wrap_succeed__set__Llastrho

subroutine f90wrap_succeed__get__Lsr(f90wrap_Lsr)
    use succeed, only: succeed_Lsr => Lsr
    implicit none
    logical, intent(out) :: f90wrap_Lsr
    
    f90wrap_Lsr = succeed_Lsr
end subroutine f90wrap_succeed__get__Lsr

subroutine f90wrap_succeed__set__Lsr(f90wrap_Lsr)
    use succeed, only: succeed_Lsr => Lsr
    implicit none
    logical, intent(in) :: f90wrap_Lsr
    
    succeed_Lsr = f90wrap_Lsr
end subroutine f90wrap_succeed__set__Lsr

subroutine f90wrap_succeed__get__Lsrho(f90wrap_Lsrho)
    use succeed, only: succeed_Lsrho => Lsrho
    implicit none
    logical, intent(out) :: f90wrap_Lsrho
    
    f90wrap_Lsrho = succeed_Lsrho
end subroutine f90wrap_succeed__get__Lsrho

subroutine f90wrap_succeed__set__Lsrho(f90wrap_Lsrho)
    use succeed, only: succeed_Lsrho => Lsrho
    implicit none
    logical, intent(in) :: f90wrap_Lsrho
    
    succeed_Lsrho = f90wrap_Lsrho
end subroutine f90wrap_succeed__set__Lsrho

subroutine f90wrap_succeed__get__Lspsi(f90wrap_Lspsi)
    use succeed, only: succeed_Lspsi => Lspsi
    implicit none
    logical, intent(out) :: f90wrap_Lspsi
    
    f90wrap_Lspsi = succeed_Lspsi
end subroutine f90wrap_succeed__get__Lspsi

subroutine f90wrap_succeed__set__Lspsi(f90wrap_Lspsi)
    use succeed, only: succeed_Lspsi => Lspsi
    implicit none
    logical, intent(in) :: f90wrap_Lspsi
    
    succeed_Lspsi = f90wrap_Lspsi
end subroutine f90wrap_succeed__set__Lspsi

subroutine f90wrap_succeed__get__Lsrho_at(f90wrap_Lsrho_at)
    use succeed, only: succeed_Lsrho_at => Lsrho_at
    implicit none
    logical, intent(out) :: f90wrap_Lsrho_at
    
    f90wrap_Lsrho_at = succeed_Lsrho_at
end subroutine f90wrap_succeed__get__Lsrho_at

subroutine f90wrap_succeed__set__Lsrho_at(f90wrap_Lsrho_at)
    use succeed, only: succeed_Lsrho_at => Lsrho_at
    implicit none
    logical, intent(in) :: f90wrap_Lsrho_at
    
    succeed_Lsrho_at = f90wrap_Lsrho_at
end subroutine f90wrap_succeed__set__Lsrho_at

subroutine f90wrap_succeed__array__rho1(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rho1 => rho1
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_rho1)) then
        dshape(1:2) = shape(succeed_rho1)
        dloc = loc(succeed_rho1)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rho1

subroutine f90wrap_succeed__array__rho2(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rho2 => rho2
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_rho2)) then
        dshape(1:2) = shape(succeed_rho2)
        dloc = loc(succeed_rho2)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rho2

subroutine f90wrap_succeed__array__rho3(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rho3 => rho3
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_rho3)) then
        dshape(1:2) = shape(succeed_rho3)
        dloc = loc(succeed_rho3)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rho3

subroutine f90wrap_succeed__array__r1(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_r1 => r1
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_r1)) then
        dshape(1:2) = shape(succeed_r1)
        dloc = loc(succeed_r1)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__r1

subroutine f90wrap_succeed__array__r2(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_r2 => r2
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_r2)) then
        dshape(1:2) = shape(succeed_r2)
        dloc = loc(succeed_r2)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__r2

subroutine f90wrap_succeed__array__r3(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_r3 => r3
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_r3)) then
        dshape(1:2) = shape(succeed_r3)
        dloc = loc(succeed_r3)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__r3

subroutine f90wrap_succeed__array__psi1(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_psi1 => psi1
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(succeed_psi1)) then
        dshape(1:3) = shape(succeed_psi1)
        dloc = loc(succeed_psi1)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__psi1

subroutine f90wrap_succeed__array__psi2(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_psi2 => psi2
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(succeed_psi2)) then
        dshape(1:3) = shape(succeed_psi2)
        dloc = loc(succeed_psi2)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__psi2

subroutine f90wrap_succeed__array__psi3(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_psi3 => psi3
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(succeed_psi3)) then
        dshape(1:3) = shape(succeed_psi3)
        dloc = loc(succeed_psi3)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__psi3

subroutine f90wrap_succeed__array__rho_at(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rho_at => rho_at
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_rho_at)) then
        dshape(1:2) = shape(succeed_rho_at)
        dloc = loc(succeed_rho_at)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rho_at

subroutine f90wrap_succeed__array__rhoi(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rhoi => rhoi
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(succeed_rhoi)) then
        dshape(1:3) = shape(succeed_rhoi)
        dloc = loc(succeed_rhoi)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rhoi

subroutine f90wrap_succeed__array__rho_at1(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rho_at1 => rho_at1
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(succeed_rho_at1)) then
        dshape(1:2) = shape(succeed_rho_at1)
        dloc = loc(succeed_rho_at1)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rho_at1

subroutine f90wrap_succeed__array__rhoi1(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module
    use smpi_math_module
    use succeed, only: succeed_rhoi1 => rhoi1
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(succeed_rhoi1)) then
        dshape(1:3) = shape(succeed_rhoi1)
        dloc = loc(succeed_rhoi1)
    else
        dloc = 0
    end if
end subroutine f90wrap_succeed__array__rhoi1

subroutine f90wrap_succeed__get__dvol(f90wrap_dvol)
    use succeed, only: succeed_dvol => dvol
    implicit none
    real(8), intent(out) :: f90wrap_dvol
    
    f90wrap_dvol = succeed_dvol
end subroutine f90wrap_succeed__get__dvol

subroutine f90wrap_succeed__set__dvol(f90wrap_dvol)
    use succeed, only: succeed_dvol => dvol
    implicit none
    real(8), intent(in) :: f90wrap_dvol
    
    succeed_dvol = f90wrap_dvol
end subroutine f90wrap_succeed__set__dvol

subroutine f90wrap_succeed__get__counter1(f90wrap_counter1)
    use succeed, only: succeed_counter1 => counter1
    implicit none
    integer(4), intent(out) :: f90wrap_counter1
    
    f90wrap_counter1 = succeed_counter1
end subroutine f90wrap_succeed__get__counter1

subroutine f90wrap_succeed__set__counter1(f90wrap_counter1)
    use succeed, only: succeed_counter1 => counter1
    implicit none
    integer(4), intent(in) :: f90wrap_counter1
    
    succeed_counter1 = f90wrap_counter1
end subroutine f90wrap_succeed__set__counter1

subroutine f90wrap_succeed__get__counter2(f90wrap_counter2)
    use succeed, only: succeed_counter2 => counter2
    implicit none
    integer(4), intent(out) :: f90wrap_counter2
    
    f90wrap_counter2 = succeed_counter2
end subroutine f90wrap_succeed__get__counter2

subroutine f90wrap_succeed__set__counter2(f90wrap_counter2)
    use succeed, only: succeed_counter2 => counter2
    implicit none
    integer(4), intent(in) :: f90wrap_counter2
    
    succeed_counter2 = f90wrap_counter2
end subroutine f90wrap_succeed__set__counter2

subroutine f90wrap_succeed__get__counter3(f90wrap_counter3)
    use succeed, only: succeed_counter3 => counter3
    implicit none
    integer(4), intent(out) :: f90wrap_counter3
    
    f90wrap_counter3 = succeed_counter3
end subroutine f90wrap_succeed__get__counter3

subroutine f90wrap_succeed__set__counter3(f90wrap_counter3)
    use succeed, only: succeed_counter3 => counter3
    implicit none
    integer(4), intent(in) :: f90wrap_counter3
    
    succeed_counter3 = f90wrap_counter3
end subroutine f90wrap_succeed__set__counter3

subroutine f90wrap_succeed__get__alpha(f90wrap_alpha)
    use succeed, only: succeed_alpha => alpha
    implicit none
    real(8), intent(out) :: f90wrap_alpha
    
    f90wrap_alpha = succeed_alpha
end subroutine f90wrap_succeed__get__alpha

subroutine f90wrap_succeed__set__alpha(f90wrap_alpha)
    use succeed, only: succeed_alpha => alpha
    implicit none
    real(8), intent(in) :: f90wrap_alpha
    
    succeed_alpha = f90wrap_alpha
end subroutine f90wrap_succeed__set__alpha

subroutine f90wrap_succeed__get__beta(f90wrap_beta)
    use succeed, only: succeed_beta => beta
    implicit none
    real(8), intent(out) :: f90wrap_beta
    
    f90wrap_beta = succeed_beta
end subroutine f90wrap_succeed__get__beta

subroutine f90wrap_succeed__set__beta(f90wrap_beta)
    use succeed, only: succeed_beta => beta
    implicit none
    real(8), intent(in) :: f90wrap_beta
    
    succeed_beta = f90wrap_beta
end subroutine f90wrap_succeed__set__beta

! End of module succeed defined in file Succeed_module.fpp

