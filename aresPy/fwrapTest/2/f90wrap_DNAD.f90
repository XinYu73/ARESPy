! Module dual_num_auto_diff defined in file DNAD.fpp

subroutine f90wrap_dual_num__get__x_ad_(this, f90wrap_x_ad_)
    use dual_num_auto_diff, only: dual_num
    implicit none
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(in)   :: this(2)
    type(dual_num_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x_ad_
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x_ad_ = this_ptr%p%x_ad_
end subroutine f90wrap_dual_num__get__x_ad_

subroutine f90wrap_dual_num__set__x_ad_(this, f90wrap_x_ad_)
    use dual_num_auto_diff, only: dual_num
    implicit none
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(in)   :: this(2)
    type(dual_num_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x_ad_
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x_ad_ = f90wrap_x_ad_
end subroutine f90wrap_dual_num__set__x_ad_

subroutine f90wrap_dual_num__array__xp_ad_(this, nd, dtype, dshape, dloc)
    use dual_num_auto_diff, only: dual_num
    implicit none
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(in) :: this(2)
    type(dual_num_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%xp_ad_)
    dloc = loc(this_ptr%p%xp_ad_)
end subroutine f90wrap_dual_num__array__xp_ad_

subroutine f90wrap_dual_num_initialise(this)
    use dual_num_auto_diff, only: dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_dual_num_initialise

subroutine f90wrap_dual_num_finalise(this)
    use dual_num_auto_diff, only: dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_dual_num_finalise

subroutine f90wrap_abs_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, abs
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = abs(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_abs_d

subroutine f90wrap_acos_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, acos
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = acos(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_acos_d

subroutine f90wrap_asin_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, asin
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = asin(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_asin_d

subroutine f90wrap_cos_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, cos
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = cos(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_cos_d

subroutine f90wrap_exp_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, exp
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = exp(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_exp_d

subroutine f90wrap_int_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, int
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(out) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    ret_res = int(u=u_ptr%p)
end subroutine f90wrap_int_d

subroutine f90wrap_log_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, log
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = log(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_log_d

subroutine f90wrap_log10_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, log10
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = log10(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_log10_d

subroutine f90wrap_nint_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, nint
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(out) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    ret_res = nint(u=u_ptr%p)
end subroutine f90wrap_nint_d

subroutine f90wrap_sin_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, sin
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = sin(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_sin_d

subroutine f90wrap_sqrt_d(ret_res, u)
    use dual_num_auto_diff, only: dual_num, sqrt
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = sqrt(u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_sqrt_d

subroutine f90wrap_max_dd(val1, val2, ret_res, val3, val4, val5)
    use dual_num_auto_diff, only: max, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: val1_ptr
    integer, intent(in), dimension(2) :: val1
    type(dual_num_ptr_type) :: val2_ptr
    integer, intent(in), dimension(2) :: val2
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: val3_ptr
    integer, optional, intent(in), dimension(2) :: val3
    type(dual_num_ptr_type) :: val4_ptr
    integer, optional, intent(in), dimension(2) :: val4
    type(dual_num_ptr_type) :: val5_ptr
    integer, optional, intent(in), dimension(2) :: val5
    val1_ptr = transfer(val1, val1_ptr)
    val2_ptr = transfer(val2, val2_ptr)
    if (present(val3)) then
        val3_ptr = transfer(val3, val3_ptr)
    else
        val3_ptr%p => null()
    end if
    if (present(val4)) then
        val4_ptr = transfer(val4, val4_ptr)
    else
        val4_ptr%p => null()
    end if
    if (present(val5)) then
        val5_ptr = transfer(val5, val5_ptr)
    else
        val5_ptr%p => null()
    end if
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = max(val1=val1_ptr%p, val2=val2_ptr%p, val3=val3_ptr%p, val4=val4_ptr%p, val5=val5_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_max_dd

subroutine f90wrap_max_di(u, ret_res, n)
    use dual_num_auto_diff, only: max, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    integer, intent(in) :: n
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = max(u=u_ptr%p, n=n)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_max_di

subroutine f90wrap_max_dr(u, ret_res, n)
    use dual_num_auto_diff, only: max, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    real(8), intent(in) :: n
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = max(u=u_ptr%p, n=n)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_max_dr

subroutine f90wrap_max_ds(u, ret_res, n)
    use dual_num_auto_diff, only: max, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    real(4), intent(in) :: n
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = max(u=u_ptr%p, n=n)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_max_ds

subroutine f90wrap_max_rd(r, ret_res, u)
    use dual_num_auto_diff, only: max, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    real(8), intent(in) :: r
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = max(r=r, u=u_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_max_rd

subroutine f90wrap_min_dd(val1, val2, ret_res, val3, val4)
    use dual_num_auto_diff, only: min, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: val1_ptr
    integer, intent(in), dimension(2) :: val1
    type(dual_num_ptr_type) :: val2_ptr
    integer, intent(in), dimension(2) :: val2
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: val3_ptr
    integer, optional, intent(in), dimension(2) :: val3
    type(dual_num_ptr_type) :: val4_ptr
    integer, optional, intent(in), dimension(2) :: val4
    val1_ptr = transfer(val1, val1_ptr)
    val2_ptr = transfer(val2, val2_ptr)
    if (present(val3)) then
        val3_ptr = transfer(val3, val3_ptr)
    else
        val3_ptr%p => null()
    end if
    if (present(val4)) then
        val4_ptr = transfer(val4, val4_ptr)
    else
        val4_ptr%p => null()
    end if
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = min(val1=val1_ptr%p, val2=val2_ptr%p, val3=val3_ptr%p, val4=val4_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_min_dd

subroutine f90wrap_min_dr(u, ret_res, n)
    use dual_num_auto_diff, only: min, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    real(8), intent(in) :: n
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = min(u=u_ptr%p, n=n)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_min_dr

subroutine f90wrap_min_ds(u, ret_res, n)
    use dual_num_auto_diff, only: min, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: u_ptr
    integer, intent(in), dimension(2) :: u
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    real(4), intent(in) :: n
    u_ptr = transfer(u, u_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = min(u=u_ptr%p, n=n)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_min_ds

subroutine f90wrap_sign_dd(val1, ret_res, val2)
    use dual_num_auto_diff, only: sign, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: val1_ptr
    integer, intent(in), dimension(2) :: val1
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: val2_ptr
    integer, intent(in), dimension(2) :: val2
    val1_ptr = transfer(val1, val1_ptr)
    val2_ptr = transfer(val2, val2_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = sign(val1=val1_ptr%p, val2=val2_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_sign_dd

subroutine f90wrap_sign_rd(val1, ret_res, val2)
    use dual_num_auto_diff, only: sign, dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    real(8), intent(in) :: val1
    type(dual_num_ptr_type) :: ret_res_ptr
    integer, intent(out), dimension(2) :: ret_res
    type(dual_num_ptr_type) :: val2_ptr
    integer, intent(in), dimension(2) :: val2
    val2_ptr = transfer(val2, val2_ptr)
    allocate(ret_res_ptr%p)
    ret_res_ptr%p = sign(val1=val1, val2=val2_ptr%p)
    ret_res = transfer(ret_res_ptr, ret_res)
end subroutine f90wrap_sign_rd

subroutine f90wrap_dual_num_auto_diff__get__NDV_AD(f90wrap_NDV_AD)
    use dual_num_auto_diff, only: dual_num_auto_diff_NDV_AD => NDV_AD
    implicit none
    integer(4), intent(out) :: f90wrap_NDV_AD
    
    f90wrap_NDV_AD = dual_num_auto_diff_NDV_AD
end subroutine f90wrap_dual_num_auto_diff__get__NDV_AD

! End of module dual_num_auto_diff defined in file DNAD.fpp

