! Module fourier defined in file Fourier.fpp

subroutine f90wrap_planfft(dimx, dimy, dimz)
    use fourier, only: planfft
    implicit none
    
    integer(4), intent(in) :: dimx
    integer(4), intent(in) :: dimy
    integer(4), intent(in) :: dimz
    call planfft(dimX=dimx, dimY=dimy, dimZ=dimz)
end subroutine f90wrap_planfft

subroutine f90wrap_planfst(dimx, dimy, dimz)
    use fourier, only: planfst
    implicit none
    
    integer(4), intent(in) :: dimx
    integer(4), intent(in) :: dimy
    integer(4), intent(in) :: dimz
    call planfst(dimX=dimx, dimY=dimy, dimZ=dimz)
end subroutine f90wrap_planfst

subroutine f90wrap_getfftdims(dimx, dimy, dimz)
    use fourier, only: getfftdims
    implicit none
    
    integer(4), intent(out) :: dimx
    integer(4), intent(out) :: dimy
    integer(4), intent(out) :: dimz
    call getfftdims(dimX=dimx, dimY=dimy, dimZ=dimz)
end subroutine f90wrap_getfftdims

subroutine f90wrap_getfftcomplexdims(dimx, dimy, dimz)
    use fourier, only: getfftcomplexdims
    implicit none
    
    integer(4), intent(out) :: dimx
    integer(4), intent(out) :: dimy
    integer(4), intent(out) :: dimz
    call getfftcomplexdims(dimX=dimx, dimY=dimy, dimZ=dimz)
end subroutine f90wrap_getfftcomplexdims

subroutine f90wrap_forwardfst(array, n0, n1, n2, n3)
    use fourier, only: forwardfst
    implicit none
    
    real(8), dimension(n0,n1,n2,n3) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(array) :: n1 = shape(array,1)
    integer :: n2
    !f2py intent(hide), depend(array) :: n2 = shape(array,2)
    integer :: n3
    !f2py intent(hide), depend(array) :: n3 = shape(array,3)
    call forwardfst(array=array)
end subroutine f90wrap_forwardfst

subroutine f90wrap_backfst(array, n0, n1, n2, n3)
    use fourier, only: backfst
    implicit none
    
    real(8), dimension(n0,n1,n2,n3) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(array) :: n1 = shape(array,1)
    integer :: n2
    !f2py intent(hide), depend(array) :: n2 = shape(array,2)
    integer :: n3
    !f2py intent(hide), depend(array) :: n3 = shape(array,3)
    call backfst(array=array)
end subroutine f90wrap_backfst

subroutine f90wrap_cleanfft
    use fourier, only: cleanfft
    implicit none
    
    call cleanfft()
end subroutine f90wrap_cleanfft

subroutine f90wrap_forwardfft_4d(ret_transform, array, n0, n1, n2, n3, n4, n5, n6, n7)
    use fourier, only: fft
    implicit none
    
    complex(8), intent(out), dimension(n0,n1,n2,n3) :: ret_transform
    real(8), dimension(n4,n5,n6,n7) :: array
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    integer :: n4
    !f2py intent(hide), depend(array) :: n4 = shape(array,0)
    integer :: n5
    !f2py intent(hide), depend(array) :: n5 = shape(array,1)
    integer :: n6
    !f2py intent(hide), depend(array) :: n6 = shape(array,2)
    integer :: n7
    !f2py intent(hide), depend(array) :: n7 = shape(array,3)
    ret_transform = fft(array=array)
end subroutine f90wrap_forwardfft_4d

subroutine f90wrap_backfft_4d(ret_transform, array, n0, n1, n2, n3, n4, n5, n6, n7)
    use fourier, only: fft
    implicit none
    
    real(8), intent(out), dimension(n0,n1,n2,n3) :: ret_transform
    complex(8), dimension(n4,n5,n6,n7) :: array
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    integer :: n4
    !f2py intent(hide), depend(array) :: n4 = shape(array,0)
    integer :: n5
    !f2py intent(hide), depend(array) :: n5 = shape(array,1)
    integer :: n6
    !f2py intent(hide), depend(array) :: n6 = shape(array,2)
    integer :: n7
    !f2py intent(hide), depend(array) :: n7 = shape(array,3)
    ret_transform = fft(array=array)
end subroutine f90wrap_backfft_4d

subroutine f90wrap_forwardfft_3d(ret_transform, array, n0, n1, n2, n3, n4, n5)
    use fourier, only: fft
    implicit none
    
    complex(8), intent(out), dimension(n0,n1,n2) :: ret_transform
    real(8), dimension(n3,n4,n5) :: array
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    !f2py intent(hide), depend(array) :: n3 = shape(array,0)
    integer :: n4
    !f2py intent(hide), depend(array) :: n4 = shape(array,1)
    integer :: n5
    !f2py intent(hide), depend(array) :: n5 = shape(array,2)
    ret_transform = fft(array=array)
end subroutine f90wrap_forwardfft_3d

subroutine f90wrap_backfft_3d(ret_transform, array, n0, n1, n2, n3, n4, n5)
    use fourier, only: fft
    implicit none
    
    real(8), intent(out), dimension(n0,n1,n2) :: ret_transform
    complex(8), dimension(n3,n4,n5) :: array
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    !f2py intent(hide), depend(array) :: n3 = shape(array,0)
    integer :: n4
    !f2py intent(hide), depend(array) :: n4 = shape(array,1)
    integer :: n5
    !f2py intent(hide), depend(array) :: n5 = shape(array,2)
    ret_transform = fft(array=array)
end subroutine f90wrap_backfft_3d

subroutine f90wrap_fourier__get__offset(f90wrap_offset)
    use fourier, only: fourier_offset => offset
    implicit none
    integer(4), intent(out) :: f90wrap_offset
    
    f90wrap_offset = fourier_offset
end subroutine f90wrap_fourier__get__offset

subroutine f90wrap_fourier__set__offset(f90wrap_offset)
    use fourier, only: fourier_offset => offset
    implicit none
    integer(4), intent(in) :: f90wrap_offset
    
    fourier_offset = f90wrap_offset
end subroutine f90wrap_fourier__set__offset

! End of module fourier defined in file Fourier.fpp

