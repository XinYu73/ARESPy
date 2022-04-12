! Module sgfft_oct_m defined in file sgfft.fpp

subroutine f90wrap_fourier_dim(n, n_next)
    use sgfft_oct_m, only: fourier_dim
    implicit none
    
    integer, intent(in) :: n
    integer, intent(out) :: n_next
    call fourier_dim(n=n, n_next=n_next)
end subroutine f90wrap_fourier_dim

subroutine f90wrap_fft(n1, n2, n3, nd1, nd2, nd3, z, isign, inzee, n0, n1, n2)
    use sgfft_oct_m, only: fft
    implicit none
    
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: nd1
    integer, intent(in) :: nd2
    integer, intent(in) :: nd3
    real(8), intent(inout), dimension(n0,n1,n2) :: z
    integer, intent(in) :: isign
    integer, intent(inout) :: inzee
    integer :: n0
    !f2py intent(hide), depend(z) :: n0 = shape(z,0)
    integer :: n1
    !f2py intent(hide), depend(z) :: n1 = shape(z,1)
    integer :: n2
    !f2py intent(hide), depend(z) :: n2 = shape(z,2)
    call fft(n1=n1, n2=n2, n3=n3, nd1=nd1, nd2=nd2, nd3=nd3, z=z, isign=isign, inzee=inzee)
end subroutine f90wrap_fft

subroutine f90wrap_convolxc_off(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, pot, zf, scal, comm, n0, n1, n2, &
    n3, n4, n5)
    use sgfft_oct_m, only: convolxc_off
    implicit none
    
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: nd1
    integer, intent(in) :: nd2
    integer, intent(in) :: nd3
    integer, intent(in) :: md1
    integer, intent(in) :: md2
    integer, intent(in) :: md3
    integer, intent(in) :: nproc
    integer, intent(in) :: iproc
    real(8), intent(in), dimension(n0,n1,n2) :: pot
    real(8), intent(inout), dimension(n3,n4,n5) :: zf
    real(8), intent(in) :: scal
    integer, intent(in) :: comm
    integer :: n0
    !f2py intent(hide), depend(pot) :: n0 = shape(pot,0)
    integer :: n1
    !f2py intent(hide), depend(pot) :: n1 = shape(pot,1)
    integer :: n2
    !f2py intent(hide), depend(pot) :: n2 = shape(pot,2)
    integer :: n3
    !f2py intent(hide), depend(zf) :: n3 = shape(zf,0)
    integer :: n4
    !f2py intent(hide), depend(zf) :: n4 = shape(zf,1)
    integer :: n5
    !f2py intent(hide), depend(zf) :: n5 = shape(zf,2)
    call convolxc_off(n1=n1, n2=n2, n3=n3, nd1=nd1, nd2=nd2, nd3=nd3, md1=md1, md2=md2, md3=md3, nproc=nproc, iproc=iproc, &
        pot=pot, zf=zf, scal=scal, comm=comm)
end subroutine f90wrap_convolxc_off

! End of module sgfft_oct_m defined in file sgfft.fpp

