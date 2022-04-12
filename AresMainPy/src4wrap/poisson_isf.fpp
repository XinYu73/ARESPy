# 1 "poisson_isf.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "poisson_isf.f90"
MODULE poisson_isf
  USE constants
  USE sgfft_oct_m
  USE grid_module ,ONLY : n1,n2,n3,n,gap
  IMPLICIT NONE
  ! include "constant.h"
  !> ----------------- copy ------------------
  !! Copies a vector, x, to a vector, y.
  interface blas_copy
    subroutine scopy(n_lib, dx, incx, dy, incy)
      implicit none
      integer,    intent(in)  :: n_lib, incx, incy
      real(KIND(1.0)),    intent(in)  :: dx !< dx(n)
      real(KIND(1.0)),    intent(out) :: dy !< dy(n)
    end subroutine scopy

    subroutine dcopy(n_lib, dx, incx, dy, incy)
      implicit none
      integer,    intent(in)  :: n_lib, incx, incy
      real(KIND(1.d0)),    intent(in)  :: dx !< dx(n)
      real(KIND(1.d0)),    intent(out) :: dy !< dy(n)
    end subroutine dcopy

    subroutine ccopy(n_lib, dx, incx, dy, incy)
      implicit none
      integer,    intent(in)  :: n_lib, incx, incy
      complex(KIND((1.0,1.0))), intent(in)  :: dx !< dx(n)
      complex(KIND((1.0,1.0))), intent(out) :: dy !< dy(n)
    end subroutine ccopy

    subroutine zcopy(n_lib, dx, incx, dy, incy)
      implicit none
      integer,    intent(in)  :: n_lib, incx, incy
      complex(KIND((1.d0,1.d0))), intent(in)  :: dx !< dx(n)
      complex(KIND((1.d0,1.d0))), intent(out) :: dy !< dy(n)
    end subroutine zcopy
  end interface blas_copy
  REAL(DP),allocatable :: karray(:,:,:)

  ! REAL(8) :: x_scf(0:1024)
  ! REAL(8) :: y_scf(0:1024)
  ! INTEGER(4) :: n_range
  ! ! Indices for the cnf array
  ! integer, parameter :: SERIAL = 1
  ! integer, parameter :: WORLD = 2
  ! integer, parameter :: DOMAIN = 3
  ! integer, parameter :: N_CNF = 3
  ! Datatype to store kernel values to solve Poisson equation
  ! on different communicators (configurations).

  ! type isf_cnf_t
  !   real(8), pointer  :: kernel(:, :, :)
  !   integer           :: nfft1, nfft2, nfft3
  !   ! type(mpi_grp_t)   :: mpi_grp
  !   logical           :: all_nodes
  ! end type isf_cnf_t

  ! type poisson_isf_t
  !   integer         :: all_nodes_comm
  !   type(isf_cnf_t) :: cnf(1:N_CNF)
  ! end type poisson_isf_t

  ! integer, parameter :: order_scaling_function = 8

  !   call scaling_function(8,1024,n_range,x_scf,y_scf)
  !   !> print scaling function
  !   print*,"n_range",n_range
  !   open(666,file="rxyz")
  !   write(666,*)x_scf
  !   close(666)
  !   open(666,file="frxyz")
  !   write(666,*)y_scf
  !   close(666)

CONTAINS
  !-----------------------------divided line--------------------------------
  SUBROUTINE karray_set()
    USE m_time_evaluate, ONLY: memory_sum,memory_free
    IMPLICIT NONE
    !> local
    INTEGER(I4B) :: nfft1,nfft2,nfft3
    INTEGER(I4B) :: karray_n1,karray_n2,karray_n3
    INTEGER(I4B) :: isf_order

    !> init recip dimensions
    CALL calculate_dimensions(n1,n2,n3,nfft1,nfft2,nfft3)
    isf_order=8

    !> set karray
    IF(allocated(karray))THEN
       call memory_free('hartree_karray',real(size(karray),DP)*DP)
       deallocate(karray)
    ENDIF
    karray_n1=nfft1/2+1
    karray_n2=nfft2/2+1
    karray_n3=nfft3/2+1
    allocate(karray(1:karray_n1,1:karray_n2,1:karray_n3))
    call memory_sum('hartree_karray',real(size(karray),DP)*DP)

    !> subroutine from octopus
    karray=0.d0
    CALL build_kernel(n1,n2,n3,nfft1,nfft2,nfft3,gap(1),isf_order,karray)

  ENDSUBROUTINE karray_set
  !-----------------------------divided line--------------------------------
  subroutine poisson_isf_method(rho,vh)
    implicit none
    REAL(DP),INTENT(IN)   :: rho(n1,n2,n3)
    REAL(DP),INTENT(OUT)  :: vh(n1,n2,n3)

    !> local
    REAL(DP)     :: rhopot(n1,n2,n3)
    INTEGER(I4B)     :: nfft1,nfft2,nfft3

    !> init setting
    rhopot=rho
    CALL calculate_dimensions(n1,n2,n3,nfft1,nfft2,nfft3)

    !> subroutine form octopus
    CALL psolver_kernel(n1,n2,n3,nfft1,nfft2,nfft3,gap(1),karray,rhopot)

    !> end
    vh=rhopot

  ENDSUBROUTINE poisson_isf_method
  !-----------------------------divided line--------------------------------
  subroutine build_kernel(n01,n02,n03,nfft1,nfft2,nfft3,hgrid,itype_scf,karrayout)
    USE m_time_evaluate, ONLY: memory_sum,memory_free
    integer, intent(in) :: n01,n02,n03,nfft1,nfft2,nfft3,itype_scf
    real(DP), intent(in) :: hgrid
    real(DP), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1), intent(out) :: karrayout

    !Do not touch !!!!
    integer, parameter :: N_GAUSS = 89
    !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
    integer, parameter :: n_points = 2**6

    !Better p_gauss for calculation
    !(the support of the exponential should be inside [-n_range/2,n_range/2])
    real(DP), parameter :: p0_ref = 1.d0
    real(DP), dimension(N_GAUSS) :: p_gauss,w_gauss

    real(DP), allocatable :: kernel_scf(:), kern_1_scf(:)
    real(DP), allocatable :: x_scf(:), y_scf(:)
    real(DP), allocatable :: karrayhalf(:, :, :)

    real(DP) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range
    real(DP) :: factor,factor2,dx,absci,p0gauss,p0_cell
    real(DP) :: a1,a2,a3
    integer :: nd1,nd2,nd3,n1k,n2k,n3k,n_scf
    integer :: i_gauss,n_range,n_cell
    integer :: i,n_iter,i1,i2,i3,i_kern
    integer :: i01,i02,i03,inkee,n1h,n2h,n3h,nd1h

    !Number of integration points : 2*itype_scf*n_points
    n_scf=2*itype_scf*n_points
    !Dimensions of Kernel
    n1k=nfft1/2+1 ; n2k=nfft2/2+1 ; n3k=nfft3/2+1
    n1h=nfft1/2  ; n2h=nfft2/2 ; n3h=nfft3/2
    nd1 = nfft1 + modulo(nfft1+1,2)
    nd2 = nfft2 + modulo(nfft2+1,2)
    nd3 = nfft3 + modulo(nfft3+1,2)

    !Half size for the half FFT
    nd1h=(nd1+1)/2

    !Allocations
    ALLOCATE(x_scf(0:n_scf))
    ALLOCATE(y_scf(0:n_scf))
    call memory_sum('hartree_karray_xy',real(size(x_scf),DP)*DP*2)

    !Build the scaling function
    call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
    !Step grid for the integration
    dx = real(n_range,8)/real(n_scf,8)
    !Extend the range (no more calculations because fill in by 0.0_8)
    n_cell = max(n01,n02,n03)
    n_range = max(n_cell,n_range)

    !Allocations
    ALLOCATE(kernel_scf(-n_range:n_range))
    ALLOCATE(kern_1_scf(-n_range:n_range))
    call memory_sum('hartree_karray_kernel',real(size(kernel_scf),DP)*DP*2)

    !Lengthes of the box (use FFT dimension)
    a1 = hgrid * real(n01,8)
    a2 = hgrid * real(n02,8)
    a3 = hgrid * real(n03,8)

    x_scf(:) = hgrid * x_scf(:)
    y_scf(:) = 1.d0/hgrid * y_scf(:)
    dx = hgrid * dx
    !To have a correct integration
    p0_cell = p0_ref/(hgrid*hgrid)

    !Initialisation of the gaussian (Beylkin)
    call gequad(N_GAUSS,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)

    !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3) 
    !(biggest length in the cube)
    !We divide the p_gauss by a_range**2 and a_gauss by a_range
    a_range = sqrt(a1*a1+a2*a2+a3*a3)
    factor = 1.d0/a_range
    !factor2 = factor*factor
    factor2 = 1.d0/(a1*a1+a2*a2+a3*a3)
    do i_gauss=1,N_GAUSS
      p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
    end do
    do i_gauss=1,N_GAUSS
      w_gauss(i_gauss) = factor*w_gauss(i_gauss)
    end do

    karrayout(:,:,:) = 0.d0
    ! karrayout(:,:,:) = 0.0_8

    !Use in this order (better for accuracy).
    loop_gauss: do i_gauss=N_GAUSS,1,-1
      !Gaussian
      pgauss = p_gauss(i_gauss)

      !We calculate the number of iterations to go from pgauss to p0_ref
      n_iter = nint((log(pgauss) - log(p0_cell))/log(4.d0))
      if (n_iter <= 0)then
        n_iter = 0
        p0gauss = pgauss
      else
        p0gauss = pgauss/4.d0**n_iter
      end if

      !Stupid integration
      !Do the integration with the exponential centered in i_kern
      kernel_scf(:) = 0.0_8
      do i_kern=0,n_range
        kern = 0.0_8
        do i=0,n_scf
          absci = x_scf(i) - real(i_kern,8)*hgrid
          absci = absci*absci
          kern = kern + y_scf(i)*exp(-p0gauss*absci)*dx
        end do
        kernel_scf(i_kern) = kern
        kernel_scf(-i_kern) = kern
        if (abs(kern) < 1.d-18) then
          !Too small not useful to calculate
          exit
        end if
      end do

      !Start the iteration to go from p0gauss to pgauss
      call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

      !Add to the kernel.
      do i3=1,n03
        i03 = i3-1
        do i2=1,n02
          i02 = i2-1
          do i1=1,n01
            i01 = i1-1
            karrayout(i1,i2,i3) = karrayout(i1,i2,i3) + w_gauss(i_gauss)* &
              kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          end do
        end do
      end do

    end do loop_gauss

    call memory_free('hartree_karray',real(size(kernel_scf),DP)*DP*2+&
         &size(x_scf)*DP*2)
    DEALLOCATE(kernel_scf)
    DEALLOCATE(kern_1_scf)
    DEALLOCATE(x_scf)
    DEALLOCATE(y_scf)

    !Set karray
    ALLOCATE(karrayhalf(1:2, 1:nd1h*nd2*nd3, 1:2))
    call memory_sum('hartree_karray',real(size(karrayhalf),DP)*DP)

    !Set karray : use mirror symmetries
    inkee=1
    call karrayhalf_in(n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
      karrayout,karrayhalf)
    call fft(n1h,nfft2,nfft3,nd1h,nd2,nd3,karrayhalf,1,inkee)
    !Reconstruct the real kernel
    call kernel_recon(n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
      karrayhalf(1,1,inkee),karrayout)

    call memory_free('hartree_karray',real(size(karrayhalf),DP)*DP)
    DEALLOCATE(karrayhalf)
  end subroutine Build_Kernel
  !-----------------------------divided line--------------------------------
  subroutine scaling_function(itype,nd,nrange,a,x)
    integer, intent(in) :: itype
    integer, intent(in) :: nd
    integer, intent(out) :: nrange
    real(DP), dimension(0:nd), intent(out) :: a,x

    real(DP), dimension(:), allocatable :: y
    integer :: i,nt,ni

    !Only itype=8,14,16,20,24,30,40,50,60,100
!!$  write(unit=*,fmt="(1x,a,i0,a)") &
!!$       "Use interpolating scaling functions of ",itype," order"

    !Give the range of the scaling function
    !from -itype to itype
    ni=2*itype
    nrange = ni
    allocate(y(0:nd))

    ! plot scaling function
    x = 0.d0
    y = 0.d0

    nt=ni
    x(nt/2-1)=1.d0
    loop1: do
      nt=2*nt
      !	write(6,*) 'nd,nt',nd,nt
      select case(itype)
      case(8)
        call back_trans_8(nd,nt,x,y)
      end select
      call blas_copy(nt, y(0), 1, x(0) ,1)
      if (nt == nd) then
        exit loop1
      end if
    end do loop1

    !open (unit=1,file='scfunction',status='unknown')
    do i=0,nd
      a(i) = 1.d0*i*ni/nd-(.5d0*ni-1.d0)
      !write(1,*) 1.d0*i*ni/nd-(.5d0*ni-1.d0),x(i)
    end do
    !close(1)

    deallocate(y)
  end subroutine scaling_function
  !-----------------------------divided line--------------------------------
  subroutine back_trans_8(nd,nt,x,y)
    integer, intent(in) :: nd,nt
    real(kind=8), intent(in) :: x(0:nd-1)
    real(kind=8), intent(out) :: y(0:nd-1)

    integer :: i,j,ind

    ! include "lazy_8_inc.F90"
!! Copyright (C) 2011 X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

!!****h* BigDFT/lazy_8
!! NAME
!!   lazy_8
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 8)
!!
!! SOURCE
!!
integer, parameter :: m=10
real(kind=8), dimension(-m:m) ::  ch,cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
   ch(i)=0.d0
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
end do

! The normalization is chosen such that a constant function remains the same constant 
! on each level of the transform

ch(-7)=-5.d0/2048.d0
ch(-6)=0.d0
ch(-5)=49.d0/2048.d0
ch(-4)=0.d0
ch(-3)=-245.d0/2048.d0
ch(-2)=0.d0
ch(-1)=1225.d0/2048.d0
ch( 0)=1.d0
ch( 1)=1225.d0/2048.d0
ch( 2)=0.d0
ch( 3)=-245.d0/2048.d0
ch( 4)=0.d0
ch( 5)=49.d0/2048.d0
ch( 6)=0.d0
ch( 7)=-5.d0/2048.d0
! 
cht( 0)=1.d0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1)**(i+1)
   cgt(i+1)=ch(-i)*(-1)**(i+1)
end do
!!***

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

    do i=0,nt/2-1
      y(2*i+0)=0.d0
      y(2*i+1)=0.d0

      do j=-m/2,m/2-1

        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
          if (ind < 0) then 
            ind=ind+nt/2
            cycle loop99
          end if
          if (ind >= nt/2) then 
            ind=ind-nt/2
            cycle loop99
          end if
          exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
      end do

    end do

  end subroutine back_trans_8
  !-----------------------------divided line--------------------------------
  subroutine gequad(n_gauss, p_gauss, w_gauss, ur_gauss, dr_gauss, acc_gauss)
    integer, intent(in)    :: n_gauss
    real(DP), intent(out)   :: p_gauss(:)
    real(DP), intent(out)   :: w_gauss(:)
    real(DP), intent(out)   :: ur_gauss
    real(DP), intent(out)   :: dr_gauss
    real(DP), intent(out)   :: acc_gauss

    integer :: iunit, i, idx

    ur_gauss = 1.0_8
    dr_gauss = 1.0e-08_8
    acc_gauss = 1.0e-08_8

    ! iunit=777
    ! open(iunit,file="gequad.data")
    ! do i = 1, n_gauss
    !   read(iunit, *) idx, p_gauss(i), w_gauss(i)
    ! end do
    ! close(iunit)
    p_gauss(1:n_gauss)=(/4.96142640560223518720E+19 ,&
         1.37454269147978055680E+19 ,&
         7.58610013441204633600E+18 ,&
         4.42040691347806976000E+18 ,&
         2.61986077948367872000E+18 ,&
         1.56320138155496678400E+18 ,&
         9.35645215863028352000E+17 ,&
         5.60962910452691712000E+17 ,&
         3.36662251196867584000E+17 ,&
         2.02182531979478656000E+17 ,&
         1.21477756091902016000E+17 ,&
         7.30129825136084960000E+16 ,&
         4.38951893556421120000E+16 ,&
         2.63949482512262320000E+16 ,&
         1.58742054072786180000E+16 ,&
         9.54806587737665600000E+15 ,&
         5.74353712364571700000E+15 ,&
         3.45521487738944500000E+15 ,&
         2.07871658520326800000E+15 ,&
         1.25064667315629925000E+15 ,&
         7.52469429541933750000E+14 ,&
         4.52746033372531750000E+14 ,&
         2.72414006900059562500E+14 ,&
         1.63912168349216750000E+14 ,&
         9.86275802590865781250E+13 ,&
         5.93457701624975000000E+13 ,&
         3.57095543222962968750E+13 ,&
         2.14872890367310468750E+13 ,&
         1.29294719957726894531E+13 ,&
         7.78003375426361035156E+12 ,&
         4.68148199759876660156E+12 ,&
         2.81699550248298681641E+12 ,&
         1.69507790481958471680E+12 ,&
         1.01998486064607580566E+12 ,&
         6.13759486539856445312E+11 ,&
         3.69320183828682556152E+11 ,&
         2.22232783898905090332E+11 ,&
         1.33725247623668685913E+11 ,&
         8.04671927390362854004E+10 ,&
         4.84199582415144119263E+10 ,&
         2.91360091170559577942E+10 ,&
         1.75321747475309219360E+10 ,&
         1.05497355522109947205E+10 ,&
         6.34815321079006576538E+09 ,&
         3.81991113733594226837E+09 ,&
         2.29857747533101129532E+09 ,&
         1.38313653595483684540E+09 ,&
         8.32282908580025315285E+08 ,&
         5.00814519374587476254E+08 ,&
         3.01358090773319005966E+08 ,&
         1.81337994217503547668E+08 ,&
         1.09117589961086824536E+08 ,&
         6.56599771718640327454E+07 ,&
         3.95099693638497143984E+07 ,&
         2.37745694710665978491E+07 ,&
         1.43060135285912808031E+07 ,&
         8.60844290313506685197E+06 ,&
         5.18000974075383413583E+06 ,&
         3.11699819305746583268E+06 ,&
         1.87560993870024033822E+06 ,&
         1.12862197183979558758E+06 ,&
         6.79132441326077212580E+05 ,&
         4.08658421279877948109E+05 ,&
         2.45904473450669785962E+05 ,&
         1.47969568088321015239E+05 ,&
         8.90386123573111399310E+04 ,&
         5.35777362552358899848E+04 ,&
         3.22396513926914667536E+04 ,&
         1.93997580852362807491E+04 ,&
         1.16735323603058641311E+04 ,&
         7.02438438577707802324E+03 ,&
         4.22682479307686026004E+03 ,&
         2.54343254175354286417E+03 ,&
         1.53047486269122669000E+03 ,&
         9.20941785160749532224E+02 ,&
         5.54163803906291605017E+02 ,&
         3.33460297407856955942E+02 ,&
         2.00655057533504106004E+02 ,&
         1.20741366914147278067E+02 ,&
         7.26544243200329873389E+01 ,&
         4.37187810415471034275E+01 ,&
         2.63071631447061058395E+01 ,&
         1.58299486353816334372E+01 ,&
         9.52493152341243920489E+00 ,&
         5.72200417067776001545E+00 ,&
         3.36242234070940915203E+00 ,&
         1.75371394604499464265E+00 ,&
         6.47059326506589704842E-01 ,&
         7.27659059437082422761E-02 &
         /)
    w_gauss(1:n_gauss)=(/4.76744548452830444336E+11 ,&
         1.13748577475044212341E+10 ,&
         7.86434097688018989563E+09 ,&
         4.62733578875959014893E+09 ,&
         2.47380464827152967453E+09 ,&
         1.36290411643898773193E+09 ,&
         9.27956002904588317871E+09 ,&
         5.21593121625466060638E+09 ,&
         3.16701801106166601181E+09 ,&
         1.29291036801493048668E+08 ,&
         1.00139319988015860319E+08 ,&
         7.75892350510188341141E+07 ,&
         6.01333567950731292367E+07 ,&
         4.66141178654796853662E+07 ,&
         3.61398903394911438227E+07 ,&
         2.80225846672956384718E+07 ,&
         2.17305091809302456677E+07 ,&
         1.68524482625876963139E+07 ,&
         1.30701489345870334655E+07 ,&
         1.01371784832269288599E+07 ,&
         7.86264116300379298627E+06 ,&
         6.09861667912273760885E+06 ,&
         4.73045784039455652237E+06 ,&
         3.66928949951594183221E+06 ,&
         2.84620508362302603200E+06 ,&
         2.20777394798527006060E+06 ,&
         1.71256191589205525815E+06 ,&
         1.32843556197737087496E+06 ,&
         1.03047312759559892584E+06 ,&
         7.99345206572271417826E+05 ,&
         6.20059354143595322967E+05 ,&
         4.80986704107449331786E+05 ,&
         3.73107167700228514150E+05 ,&
         2.89424083374121342786E+05 ,&
         2.24510248231581790606E+05 ,&
         1.74155825690028956160E+05 ,&
         1.35095256919654057128E+05 ,&
         1.04795442776800307911E+05 ,&
         8.12914458222430403112E+04 ,&
         6.30590493649328709580E+04 ,&
         4.89159040455329668475E+04 ,&
         3.79448484018048766302E+04 ,&
         2.94344290473253968230E+04 ,&
         2.28327622054490057053E+04 ,&
         1.77117439501512344577E+04 ,&
         1.37392878671041762573E+04 ,&
         1.06577895710752582090E+04 ,&
         8.26742141053961859143E+03 ,&
         6.41317397520136455569E+03 ,&
         4.97480402838654299558E+03 ,&
         3.85903698188553062209E+03 ,&
         2.99351824493299136520E+03 ,&
         2.32212119668117520632E+03 ,&
         1.80130750964719641161E+03 ,&
         1.39730379659817049287E+03 ,&
         1.08391149143250686393E+03 ,&
         8.40807939169209134889E+02 ,&
         6.52228524366749411456E+02 ,&
         5.05944376983506117540E+02 ,&
         3.92469362317941090623E+02 ,&
         3.04444930257324301692E+02 ,&
         2.36162932842453614057E+02 ,&
         1.83195466078603516280E+02 ,&
         1.42107732186551459108E+02 ,&
         1.10235302157239914322E+02 ,&
         8.55113346705382326718E+01 ,&
         6.63325469806696617070E+01 ,&
         5.14552463353841389448E+01 ,&
         3.99146798429449276568E+01 ,&
         3.09624728409162095488E+01 ,&
         2.40180988122150118613E+01 ,&
         1.86312338024296586525E+01 ,&
         1.44525541233150498499E+01 ,&
         1.12110836519105934173E+01 ,&
         8.69662175848497120967E+00 ,&
         6.74611236165731931180E+00 ,&
         5.23307018057530015653E+00 ,&
         4.05937850501539543302E+00 ,&
         3.14892659076635705873E+00 ,&
         2.44267408211071623825E+00 ,&
         1.89482240522855271969E+00 ,&
         1.46984505907050078122E+00 ,&
         1.14019261330526999743E+00 ,&
         8.84791217422925324598E-01 ,&
         6.92686387080616472467E-01 ,&
         5.85244576897023249806E-01 ,&
         5.76182522545327535646E-01 ,&
         5.96688817388997150282E-01 ,&
         6.07879901151108792412E-01 &
         /)
    
  end subroutine gequad
  !-----------------------------divided line--------------------------------
  subroutine calculate_dimensions(n01,n02,n03,nfft1,nfft2,nfft3)
    integer, intent(in) :: n01,n02,n03
    integer, intent(out) :: nfft1,nfft2,nfft3
    
    integer :: i1,i2,i3,l1

    !Test 2*n01, 2*n02, 2*n03
    !WRITE(6,*) 'in dimensions_fft',n01,n02,n03
    i1=2*n01
    i2=2*n02
    i3=2*n03
    do
      call fourier_dim(i1,nfft1)
      call fourier_dim(nfft1/2,l1)
      if (modulo(nfft1,2) == 0 .and. modulo(l1,2) == 0 .and. 2*l1 == nfft1) then
        exit
      end if
      i1=i1+1
    end do
    do
      call fourier_dim(i2,nfft2)
      if (modulo(nfft2,2) == 0) then
        exit
      end if
      i2=i2+1
    end do
    do
      call fourier_dim(i3,nfft3)
      if (modulo(nfft3,2) == 0) then
        exit
      end if
      i3=i3+1
    end do
    !nd1 = nfft1 + modulo(nfft1+1,2)
    !nd2 = nfft2 + modulo(nfft2+1,2)
    !nd3 = nfft3 + modulo(nfft3+1,2)
    !WRITE(6,*) 'out dimensions_fft',nfft1,nfft2,nfft3

  end subroutine calculate_dimensions
  !-----------------------------divided line--------------------------------
  subroutine psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot)
    integer, intent(in)    :: n01
    integer, intent(in)    :: n02
    integer, intent(in)    :: n03
    integer, intent(in)    :: nfft1
    integer, intent(in)    :: nfft2
    integer, intent(in)    :: nfft3
    real(DP), intent(in)    :: hgrid
    real(DP), intent(in)    :: karray(nfft1/2 + 1,nfft2/2 + 1, nfft3/2 + 1)
    real(DP), intent(inout) :: rhopot(n01, n02, n03)

    real(DP), allocatable :: zarray(:,:,:)
    real(DP) :: factor
    integer :: n1, n2, n3, nd1, nd2, nd3, n1h, nd1h
    integer :: inzee, i_sign

    !Dimension of the FFT
    call calculate_dimensions(n01, n02, n03, n1, n2, n3)

    !Half size of nd1
    n1h=n1/2
    nd1 = n1 + modulo(n1+1,2)
    nd2 = n2 + modulo(n2+1,2)
    nd3 = n3 + modulo(n3+1,2)
    nd1h=(nd1+1)/2

    ALLOCATE(zarray(1:2, 1:nd1h*nd2*nd3, 1:2))

    !Set zarray
    call zarray_in(n01,n02,n03,nd1h,nd2,nd3,rhopot,zarray)

    !FFT
    !print *,"Do a 3D HalFFT for the density"
    i_sign=1
    inzee=1
    call fft(n1h,n2,n3,nd1h,nd2,nd3,zarray,i_sign,inzee)

    !print *, "Apply the kernel"
    call kernel_application(n1,n2,n3,nd1h,nd2,nd3,nfft1,nfft2,nfft3,zarray,karray,inzee)

    !Inverse FFT
    i_sign=-1
    !print *,"Do a 3D inverse HalFFT"
    call fft(n1h,n2,n3,nd1h,nd2,nd3,zarray,i_sign,inzee)

    !Recollect the result
    !We have to multiply by a factor
    factor = hgrid**3/(n1*n2*n3)

    ! Calling this routine gives only the Hartree potential
    call zarray_out(n01, n02, n03, nd1h, nd2, nd3, rhopot, zarray(1, 1, inzee), factor)

    DEALLOCATE(zarray)
  end subroutine psolver_kernel
  !-----------------------------divided line--------------------------------
  subroutine zarray_in(n01,n02,n03,nd1,nd2,nd3,density,zarray)
    integer :: n01,n02,n03,nd1,nd2,nd3
    real(DP), dimension(n01,n02,n03) :: density
    real(DP), dimension(2,nd1,nd2,nd3) :: zarray

    integer :: i1,i2,i3,n01h,nd1hm,nd3hm,nd2hm

    !Half the size of n01
    n01h=n01/2
    nd1hm=(nd1-1)/2
    nd2hm=(nd2-1)/2
    nd3hm=(nd3-1)/2
    !Set to zero
    do i3=1,nd3
      do i2=1,nd2
        do i1=1,nd1
          zarray(1,i1,i2,i3) = 0.0_8
          zarray(2,i1,i2,i3) = 0.0_8
        end do
      end do
    end do
    !Set zarray
    do i3=1,n03
      do i2=1,n02
        do i1=1,n01h
          zarray(1,i1+nd1hm,i2+nd2hm,i3+nd3hm) = density(2*i1-1,i2,i3)
          zarray(2,i1+nd1hm,i2+nd2hm,i3+nd3hm) = density(2*i1,i2,i3)
        end do
      end do
    end do
    if(modulo(n01,2) == 1) then
      do i3=1,n03
        do i2=1,n02
          zarray(1,n01h+1+nd1hm,i2+nd2hm,i3+nd3hm) = density(n01,i2,i3)
        end do
      end do
    end if

  end subroutine zarray_in
  !-----------------------------divided line--------------------------------
  subroutine zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor)
    integer, intent(in)    :: n01,n02,n03,nd1,nd2,nd3
    real(DP), intent(out)   :: rhopot(n01,n02,n03)
    real(DP), intent(in)    :: zarray(2*nd1,nd2,nd3)  ! Convert zarray(2,nd1,nd2,nd3) -> zarray(2*nd1,nd2,nd3) to use i1=1,n01
                                                     ! instead of i1=1,n1h + special case for modulo(n01,2)
    real(DP), intent(in)    :: factor

    integer :: i1,i2,i3

    do i3=1,n03
      do i2=1,n02
        do i1=1,n01
          rhopot(i1, i2, i3) = factor*zarray(i1,i2,i3)
        end do
      end do
    end do

  end subroutine zarray_out
  !-----------------------------divided line--------------------------------
  !!****h* BigDFT/scf_recursion
  !! NAME
  !!   scf_recursion
  !!
  !! FUNCTION
  !!   Do iterations to go from p0gauss to pgauss
  !!   order interpolating scaling function
  !!
  !! SOURCE
  !!
  subroutine scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
    integer, intent(in) :: itype,n_iter,n_range
    real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
    real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)

    !Only itype=8,14,16,20,24,30,40,50,60,100
    select case(itype)
    case(8)
      !O.K.
    case default
      ! message(1) = "Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100."
      print*,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100."
      stop
    end select

    select case(itype)
    case(8)
      call scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
    end select

  end subroutine scf_recursion
  !-----------------------------divided line--------------------------------
  !!****h* BigDFT/scf_recursion_8
  !! NAME
  !!   scf_recursion_8
  !!
  !! FUNCTION
  !!   Do iterations to go from p0gauss to pgauss
  !!   8th-order interpolating scaling function
  !!
  !! SOURCE
  !!
  subroutine scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
    integer, intent(in) :: n_iter,n_range
    real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
    real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)

    real(kind=8) :: kern,kern_tot
    integer :: i_iter,i,j,ind

! #include "lazy_8_inc.F90"
!! Copyright (C) 2011 X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

!!****h* BigDFT/lazy_8
!! NAME
!!   lazy_8
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 8)
!!
!! SOURCE
!!
integer, parameter :: m=10
real(kind=8), dimension(-m:m) ::  ch,cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
   ch(i)=0.d0
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
end do

! The normalization is chosen such that a constant function remains the same constant 
! on each level of the transform

ch(-7)=-5.d0/2048.d0
ch(-6)=0.d0
ch(-5)=49.d0/2048.d0
ch(-4)=0.d0
ch(-3)=-245.d0/2048.d0
ch(-2)=0.d0
ch(-1)=1225.d0/2048.d0
ch( 0)=1.d0
ch( 1)=1225.d0/2048.d0
ch( 2)=0.d0
ch( 3)=-245.d0/2048.d0
ch( 4)=0.d0
ch( 5)=49.d0/2048.d0
ch( 6)=0.d0
ch( 7)=-5.d0/2048.d0
! 
cht( 0)=1.d0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1)**(i+1)
   cgt(i+1)=ch(-i)*(-1)**(i+1)
end do
!!***

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

    !Start the iteration to go from p0gauss to pgauss
    loop_iter_scf: do i_iter=1,n_iter
      kern_1_scf(:) = kernel_scf(:)
      kernel_scf(:) = 0.d0
      loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
          ind = 2*i-j
          if (abs(ind) > n_range) then
            kern = 0.d0
          else
            kern = kern_1_scf(ind)
          end if
          kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
          !zero after (be sure because strictly == 0.d0)
          exit loop_iter_i
        else
          kernel_scf( i) = 0.5d0*kern_tot
          kernel_scf(-i) = kernel_scf(i)
        end if
      end do loop_iter_i
    end do loop_iter_scf
  end subroutine scf_recursion_8
  !-----------------------------divided line--------------------------------
  !!****h* BigDFT/karrayhalf_in
  !! NAME
  !!   karrayhalf_in
  !!
  !! FUNCTION
  !!    Put in the array for4446666444 FFT
  !!
  !! SOURCE
  !!
  subroutine karrayhalf_in(n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3, kernel,karrayhalf)
    integer, intent(in) :: n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3
    real(DP), dimension(n1k,n2k,n3k), intent(in) :: kernel
    real(DP), dimension(2,(nd1+1)/2,nd2,nd3), intent(out) :: karrayhalf

    real(DP), dimension(:), allocatable :: karray
    integer :: i1,i2,i3,nd1h,n1h,n2h,n3h

    !Body
    n1h=nfft1/2
    n2h=nfft2/2
    n3h=nfft3/2

    ALLOCATE(karray(1:nfft1))

    nd1h=(nd1+1)/2
    karrayhalf(:,:,:,:) = 0.0_8
    do i3=1,n03
      do i2=1,n02
        karray(:) = 0.0_8
        do i1=1,n01
          karray(i1+n1h) = kernel(i1,i2,i3)
        end do
        do i1=2,n01
          karray(n1h-i1+1+nd1-nfft1) = kernel(i1,i2,i3)
        end do
        do i1=1,n1h
          karrayhalf(1,i1,i2+n2h,i3+n3h) = karray(2*i1-1)
          karrayhalf(2,i1,i2+n2h,i3+n3h) = karray(2*i1)
        end do
      end do
      do i2=2,n02
        do i1=1,nd1h
          karrayhalf(:,i1,n2h-i2+1+nd2-nfft2,i3+n3h) = &
            karrayhalf(:,i1,i2+n2h,i3+n3h)
        end do
      end do
    end do
    do i3=2,n03
      do i2=1,nd2
        do i1=1,nd1h
          karrayhalf(:,i1,i2,n3h-i3+1+nd3-nfft3) = karrayhalf(:,i1,i2,i3+n3h)
        end do
      end do
    end do

    DEALLOCATE(karray)
  end subroutine karrayhalf_in
  !!***
  !-----------------------------divided line--------------------------------
  !!****h* BigDFT/kernel_recon
  !! NAME
  !!   kernel_recon
  !!
  !! FUNCTION
  !!    Reconstruction of the kernel from the FFT array zarray
  !!    We keep only the half kernel in each direction (x,y,z).
  !!
  !! SOURCE
  !!
  subroutine kernel_recon(n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,zarray,karray)
    integer, intent(in) :: n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3
    real(DP), dimension(2,(nd1+1)/2*nd2*nd3), intent(in) :: zarray
    real(DP), dimension(n1k,n2k,n3k), intent(out) :: karray

    real(DP), dimension(:), allocatable :: cos_array,sin_array
    integer :: i1,i2,i3,ind1,ind2,nd1h,n1h,n2h,n3h
    real(DP) :: rfe,ife,rfo,ifo,cp,sp,rk,ik,a,b,c,d,pi2

    !Body
    n1h=nfft1/2
    n2h=nfft2/2
    n3h=nfft3/2
    nd1h=(nd1+1)/2
    pi2=8.d0*datan(1.d0)
    pi2=pi2/real(nfft1,8)

    ALLOCATE(cos_array(1:nd1h))
    ALLOCATE(sin_array(1:nd1h))

    do i1=1,nd1h
      cos_array(i1)= dcos(pi2*(i1-1))
      sin_array(i1)=-dsin(pi2*(i1-1))
    end do
    do i3=1,n3h+1
      do i2=1,n2h+1
        do i1=1,nd1h
          call norm_ind(nd1h,nd2,nd3,i1,i2,i3,ind1)
          call symm_ind(nd1h,nd2,nd3,i1,i2,i3,ind2)
          a=zarray(1,ind1)
          b=zarray(2,ind1)
          c=zarray(1,ind2)
          d=zarray(2,ind2)
          rfe=0.5d0*(a+c)
          ife=0.5d0*(b-d)
          rfo=0.5d0*(a-c)
          ifo=0.5d0*(b+d) 
          cp=cos_array(i1)
          sp=sin_array(i1)
          rk=rfe+cp*ifo-sp*rfo
          ik=ife-cp*rfo-sp*ifo
          !For big dimension 1.d-9 otherwise 1.d-10
          !Remove the test
          !if(abs(ik) >= 1.d-10) then
          !   print *,"non real kernel FFT",i1,i2,i3,ik  
          !   stop
          !end if
          !Build the intermediate FFT convolution (full)
          !call norm_ind(nd1,nd2,nd3,i1,i2,i3,indA)
          karray(i1,i2,i3)=rk 
        end do
      end do
    end do

    DEALLOCATE(cos_array)
    DEALLOCATE(sin_array)

  end subroutine kernel_recon
  !-----------------------------divided line--------------------------------
  !!****h* BigDFT/norm_ind
  !! NAME
  !!   norm_ind
  !!
  !! FUNCTION
  !!   Index in zarray
  !!
  !! SOURCE
  !!
  subroutine norm_ind(nd1,nd2,nd3,i1,i2,i3,ind)
    integer :: nd1,nd2,nd3,i1,i2,i3
    integer :: ind

    !Local variables
    integer :: a1,a2,a3
    if ( i1 == nd1 ) then
      a1=1
    else
      a1=i1
    end if
    if ( i2 == nd2 ) then
      a2=1
    else
      a2=i2
    end if
    if ( i3 == nd3 ) then
      a3=1
    else
      a3=i3
    end if
    ind=a1+nd1*(a2-1)+nd1*nd2*(a3-1)
  end subroutine norm_ind
  !!***
  !-----------------------------divided line--------------------------------
  !!****h* BigDFT/symm_ind
  !! NAME
  !!   symm_ind
  !!
  !! FUNCTION
  !!   Index in zarray for -g vector
  !!
  !! SOURCE
  !!
  subroutine symm_ind(nd1,nd2,nd3,i1,i2,i3,ind)
    integer :: nd1,nd2,nd3,i1,i2,i3
    integer :: ind

    integer ::  a1,a2,a3
    if (i1 /= 1) then 
      a1=nd1+1-i1
    else
      a1=i1
    end if
    if (i2 /= 1) then 
      a2=nd2+1-i2
    else
      a2=i2
    end if
    if (i3 /= 1) then 
      a3=nd3+1-i3
    else
      a3=i3
    end if
    ind=a1+nd1*(a2-1)+nd1*nd2*(a3-1)
  end subroutine symm_ind
  !!***
  !-----------------------------divided line--------------------------------
  subroutine kernel_application( n1xy, n2xy, n3xy,nd1h,nd2,nd3,nfft1,nfft2,nfft3,zarray,karray,inzee)
    integer, intent(in)    ::  n1xy, n2xy, n3xy,nd1h,nd2,nd3,nfft1,nfft2,nfft3,inzee
    real(DP), intent(in)    :: karray(1:nfft1/2 + 1, 1:nfft2/2 + 1, 1:nfft3/2 + 1)
    real(DP), intent(inout) :: zarray(1:2, 1:nd1h, 1:nd2, 1:nd3, 1:2)

    real(DP), dimension(:), allocatable :: cos_array,sin_array
    real(DP) :: a,b,c,d,pi2,g1,cp,sp
    real(DP) :: rfe,ife,rfo,ifo,rk,ik,rk2,ik2,re,ro,ie,io,rhk,ihk
    integer :: i1,i2,i3,j1,j2,j3, ouzee,n1h,n2h,n3h
    integer :: si1,si2,si3

    n1h =  n1xy/2
    n2h =  n2xy/2
    n3h =  n3xy/2

    ALLOCATE(cos_array(1:n1h + 1))
    ALLOCATE(sin_array(1:n1h + 1))

    pi2=8.d0*datan(1.d0)
    pi2=pi2/real( n1xy,8)
    do i1=1,n1h+1
      cos_array(i1)=dcos(pi2*(i1-1))
      sin_array(i1)=-dsin(pi2*(i1-1))
    end do

    ouzee=3-inzee

    !--------------------------------------------!
    !--- Starting reconstruction half -> full ---!
    !--------------------------------------------!   

    !-------------Case i3 = 1
    i3=1
    j3=1
    si3=1

    !-------------Case i2 = 1, i3 = 1
    i2=1
    j2=1
    si2=1

    !Case i1 == 1
    i1=1
    si1=1
    a=zarray(1,i1,i2,i3,inzee)
    b=zarray(2,i1,i2,i3,inzee)
    c=zarray(1,si1,si2,si3,inzee)
    d=zarray(2,si1,si2,si3,inzee)
    rfe=.5d0*(a+c)
    ife=.5d0*(b-d)
    rfo=.5d0*(a-c)
    ifo=.5d0*(b+d) 
    cp=cos_array(i1)
    sp=sin_array(i1)
    rk=rfe+cp*ifo-sp*rfo
    ik=ife-cp*rfo-sp*ifo
    g1=karray(i1,j2,j3)
    rk2=rk*g1
    ik2=ik*g1

    zarray(1,1,i2,i3,ouzee) = rk2
    zarray(2,1,i2,i3,ouzee) = ik2

    !Case i1=2,n1h
    do i1=2,n1h
      si1=n1h+2-i1

      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
    end do

    !Case i1=n1h+1
    i1=n1h+1
    si1=n1h+2-i1

    a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
    b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
    c=zarray(1,si1,si2,si3,inzee)
    d=zarray(2,si1,si2,si3,inzee)
    rfe=.5d0*(a+c)
    ife=.5d0*(b-d)
    rfo=.5d0*(a-c)
    ifo=.5d0*(b+d) 
    cp=cos_array(i1)
    sp=sin_array(i1)
    rk=rfe+cp*ifo-sp*rfo
    ik=ife-cp*rfo-sp*ifo
    g1=karray(i1,j2,j3)
    rk2=rk*g1
    ik2=ik*g1

    zarray(1,i1,i2,i3,ouzee) = rk2
    zarray(2,i1,i2,i3,ouzee) = ik2
    !-------------END case i2 = 1 , i3=1

    !case i2 >=2
    do i2=2, n2xy
      j2=n2h+1-abs(n2h+1-i2)
      si2= n2xy+2-i2 !if i2 /=1, otherwise si2=1

      !Case i1 == 1
      i1=1
      si1=1
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,1,i2,i3,ouzee) = rk2
      zarray(2,1,i2,i3,ouzee) = ik2

      !Case i1=2,n1h
      do i1=2,n1h
        si1=n1h+2-i1

        a=zarray(1,i1,i2,i3,inzee)
        b=zarray(2,i1,i2,i3,inzee)
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,i1,i2,i3,ouzee) = rk2
        zarray(2,i1,i2,i3,ouzee) = ik2
      end do

      !Case i1=n1h+1
      i1=n1h+1
      si1=n1h+2-i1

      a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
      b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
    end do
    !-------------END Case i3 = 1

    !case i3 >=2
    do i3=2, n3xy
      j3=n3h+1-abs(n3h+1-i3)
      si3= n3xy+2-i3 !if i3 /=1, otherwise si3=1

      !-------------Case i2 = 1
      i2=1
      j2=1
      si2=1

      !Case i1 == 1
      i1=1
      si1=1
      a=zarray(1,i1,i2,i3,inzee)
      b=zarray(2,i1,i2,i3,inzee)
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,1,i2,i3,ouzee) = rk2
      zarray(2,1,i2,i3,ouzee) = ik2

      !Case i1=2,n1h
      do i1=2,n1h
        si1=n1h+2-i1

        a=zarray(1,i1,i2,i3,inzee)
        b=zarray(2,i1,i2,i3,inzee)
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,i1,i2,i3,ouzee) = rk2
        zarray(2,i1,i2,i3,ouzee) = ik2
      end do

      !Case i1=n1h+1
      i1=n1h+1
      si1=n1h+2-i1

      a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
      b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
      c=zarray(1,si1,si2,si3,inzee)
      d=zarray(2,si1,si2,si3,inzee)
      rfe=.5d0*(a+c)
      ife=.5d0*(b-d)
      rfo=.5d0*(a-c)
      ifo=.5d0*(b+d) 
      cp=cos_array(i1)
      sp=sin_array(i1)
      rk=rfe+cp*ifo-sp*rfo
      ik=ife-cp*rfo-sp*ifo
      g1=karray(i1,j2,j3)
      rk2=rk*g1
      ik2=ik*g1

      zarray(1,i1,i2,i3,ouzee) = rk2
      zarray(2,i1,i2,i3,ouzee) = ik2
      !-------------END case i2 = 1      

      !case i2 >=2
      do i2=2, n2xy
        j2=n2h+1-abs(n2h+1-i2)
        si2= n2xy+2-i2 !if i2 /=1, otherwise si2=1

        !Case i1 == 1
        i1=1
        si1=1
        a=zarray(1,i1,i2,i3,inzee)
        b=zarray(2,i1,i2,i3,inzee)
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,1,i2,i3,ouzee) = rk2
        zarray(2,1,i2,i3,ouzee) = ik2

        !Case i1=2,n1h
        do i1=2,n1h
          si1=n1h+2-i1

          a=zarray(1,i1,i2,i3,inzee)
          b=zarray(2,i1,i2,i3,inzee)
          c=zarray(1,si1,si2,si3,inzee)
          d=zarray(2,si1,si2,si3,inzee)
          rfe=.5d0*(a+c)
          ife=.5d0*(b-d)
          rfo=.5d0*(a-c)
          ifo=.5d0*(b+d) 
          cp=cos_array(i1)
          sp=sin_array(i1)
          rk=rfe+cp*ifo-sp*rfo
          ik=ife-cp*rfo-sp*ifo
          g1=karray(i1,j2,j3)
          rk2=rk*g1
          ik2=ik*g1

          zarray(1,i1,i2,i3,ouzee) = rk2
          zarray(2,i1,i2,i3,ouzee) = ik2
        end do

        !Case i1=n1h+1
        i1=n1h+1
        si1=n1h+2-i1

        a=zarray(1,1,i2,i3,inzee) !beware here i1 -> 1
        b=zarray(2,1,i2,i3,inzee) !beware here i1 -> 1
        c=zarray(1,si1,si2,si3,inzee)
        d=zarray(2,si1,si2,si3,inzee)
        rfe=.5d0*(a+c)
        ife=.5d0*(b-d)
        rfo=.5d0*(a-c)
        ifo=.5d0*(b+d) 
        cp=cos_array(i1)
        sp=sin_array(i1)
        rk=rfe+cp*ifo-sp*rfo
        ik=ife-cp*rfo-sp*ifo
        g1=karray(i1,j2,j3)
        rk2=rk*g1
        ik2=ik*g1

        zarray(1,i1,i2,i3,ouzee) = rk2
        zarray(2,i1,i2,i3,ouzee) = ik2
      end do

    end do


    !--------------------------------------------!
    !--- Starting reconstruction full -> half ---!
    !--------------------------------------------!   

    !case i3=1
    i3=1
    j3=1
    !case i2=1
    i2=1
    j2=1
    do i1=1,n1h
      j1=n1h+2-i1

      a=zarray(1,i1,i2,i3,ouzee)
      b=zarray(2,i1,i2,i3,ouzee)
      c=zarray(1,j1,j2,j3,ouzee)
      d=-zarray(2,j1,j2,j3,ouzee)
      cp=cos_array(i1)
      sp=sin_array(i1)
      re=(a+c)
      ie=(b+d)
      ro=(a-c)*cp-(b-d)*sp
      io=(a-c)*sp+(b-d)*cp
      rhk=re-io 
      ihk=ie+ro

      zarray(1,i1,i2,i3,inzee)=rhk
      zarray(2,i1,i2,i3,inzee)=ihk
    end do
    !case i2 >= 2
    do i2=2, n2xy
      j2=nd2+1-i2
      do i1=1,n1h
        j1=n1h+2-i1

        a=zarray(1,i1,i2,i3,ouzee)
        b=zarray(2,i1,i2,i3,ouzee)
        c=zarray(1,j1,j2,j3,ouzee)
        d=-zarray(2,j1,j2,j3,ouzee)
        cp=cos_array(i1)
        sp=sin_array(i1)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rhk=re-io 
        ihk=ie+ro

        zarray(1,i1,i2,i3,inzee)=rhk
        zarray(2,i1,i2,i3,inzee)=ihk
      end do
    end do


    !case i3 >=2
    do i3=2, n3xy
      j3=nd3+1-i3
      !case i2=1
      i2=1
      j2=1
      do i1=1,n1h
        j1=n1h+2-i1

        a=zarray(1,i1,i2,i3,ouzee)
        b=zarray(2,i1,i2,i3,ouzee)
        c=zarray(1,j1,j2,j3,ouzee)
        d=-zarray(2,j1,j2,j3,ouzee)
        cp=cos_array(i1)
        sp=sin_array(i1)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rhk=re-io 
        ihk=ie+ro

        zarray(1,i1,i2,i3,inzee)=rhk
        zarray(2,i1,i2,i3,inzee)=ihk
      end do
      !case i2 >= 2
      do i2=2, n2xy
        j2=nd2+1-i2
        do i1=1,n1h
          j1=n1h+2-i1

          a=zarray(1,i1,i2,i3,ouzee)
          b=zarray(2,i1,i2,i3,ouzee)
          c=zarray(1,j1,j2,j3,ouzee)
          d=-zarray(2,j1,j2,j3,ouzee)
          cp=cos_array(i1)
          sp=sin_array(i1)
          re=(a+c)
          ie=(b+d)
          ro=(a-c)*cp-(b-d)*sp
          io=(a-c)*sp+(b-d)*cp
          rhk=re-io 
          ihk=ie+ro

          zarray(1,i1,i2,i3,inzee)=rhk
          zarray(2,i1,i2,i3,inzee)=ihk
        end do
      end do

    end do

    !De-allocations
    DEALLOCATE(cos_array)
    DEALLOCATE(sin_array)

  end subroutine kernel_application
  !-----------------------------divided line--------------------------------

END MODULE poisson_isf
