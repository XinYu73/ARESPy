MODULE succeed
  USE constants
  USE struct_module , ONLY : ncharge
#ifdef MPI
  USE smpi_math_module
#endif
  IMPLICIT NONE
  LOGICAL :: Llastrho = .false. ,Lsr =.false. &
       ,Lsrho = .false. ,Lspsi = .false. &
       , Lsrho_at=.false.
  REAL(DP),allocatable :: rho1(:,:),rho2(:,:),rho3(:,:)&
       ,r1(:,:),r2(:,:),r3(:,:)&
       ,psi1(:,:,:),psi2(:,:,:),psi3(:,:,:) &
       ,rho_at(:,:),rhoi(:,:,:)&
       ,rho_at1(:,:),rhoi1(:,:,:)
  !> grid info
  REAL(DP),save :: dvol=0.d0
  INTEGER(I4B) :: counter1=0,counter2=0,counter3=0
  REAL(DP) :: alpha,beta

CONTAINS
  SUBROUTINE init_succeed_rho(n_rho,n_r,n_s,nspin,dv)
    USE m_time_evaluate, ONLY: memory_sum
    implicit none
    !> number of local rho grid
    INTEGER(I4B) :: n_rho
    !> number of atom
    INTEGER(I4B) :: n_r
    !> number of state
    INTEGER(I4B) :: n_s
    !> number of spin
    INTEGER(I4B) :: nspin
    !> volume of per grid point
    REAL(DP)     :: dv
    !> local

    if(counter1>0)return
    Llastrho=.FALSE.
    Lsr=.FALSE.
    Lsrho=.FALSE.
    allocate(rho1(n_rho,nspin))
    allocate(rho2(n_rho,nspin))
    allocate(rho3(n_rho,nspin))
    allocate(rho_at(n_rho,nspin))
    rho1=0.d0
    rho2=0.d0
    rho3=0.d0
    rho_at=0.d0
    allocate(r1(3,n_r))
    allocate(r2(3,n_r))
    allocate(r3(3,n_r))
    r1=0.d0
    r2=0.d0
    r3=0.d0
    allocate(psi1(n_rho,n_s,nspin))
    allocate(psi2(n_rho,n_s,nspin))
    allocate(psi3(n_rho,n_s,nspin))
    psi1=0.d0
    psi2=0.d0
    psi3=0.d0
    dvol=dv
    call memory_sum('succeed_rho',real(size(rho1),DP)*4*DP+size(r1)*3*DP+size(psi1)*3*DP)
  ENDSUBROUTINE init_succeed_rho

  SUBROUTINE destory_succeed()
    USE m_time_evaluate, ONLY: memory_free
    implicit none

    call memory_free('succeed_rho',real(size(rho1),DP)*4*DP+size(r1)*3*DP+size(psi1)*3*DP)
    deallocate(rho1,rho2,rho3)
    deallocate(r1,r2,r3)
    deallocate(psi1,psi2,psi3)
    deallocate(rho_at)

  ENDSUBROUTINE destory_succeed

  SUBROUTINE store_rho(n_rho,Nspin,rho)
    implicit none
    INTEGER(I4B) :: n_rho,Nspin
    REAL(DP) :: rho(n_rho,Nspin)

    rho3=rho2
    rho2=rho1
    if(.not.Lsrho_at)then
       rho1=rho
    else
       rho1=rho-rho_at
    endif

    counter1=counter1+1
    if(counter1==1)Lsrho = .true.
    if(Lsrho.and.Lsr.and.Lspsi)Llastrho=.true.

  ENDSUBROUTINE store_rho

  SUBROUTINE store_rho_at(n_rho,Nspin,rho_in)
    implicit none
    INTEGER(I4B) :: n_rho,Nspin
    REAL(DP) :: rho_in(n_rho,Nspin)

    rho_at=rho_in
    Lsrho_at=.true.

  ENDSUBROUTINE store_rho_at

  SUBROUTINE store_r(nr,r)
    implicit none
    INTEGER(I4B) :: nr
    REAL(DP)     :: r(3,nr)

    r3=r2
    r2=r1
    r1=r
    counter2=counter2+1
    if(counter2==1)Lsr = .true.
    if(Lsrho.and.Lsr.and.Lspsi)Llastrho=.true.

  ENDSUBROUTINE store_r

  SUBROUTINE store_psi(n_rho,n_s,nspin,psi)
    implicit none
    !> array dimension
    INTEGER(I4B) :: n_rho,n_s,nspin
    !> psi
    REAL(DP) :: psi(n_rho,n_s,nspin)

    psi3=psi2
    psi2=psi1
    psi1=psi
    counter3=counter3+1
    if(counter3==1)Lspsi = .true.
    if(Lsrho.and.Lsr.and.Lspsi)Llastrho=.true.

  ENDSUBROUTINE store_psi

  SUBROUTINE get_rho(nr,r_new,nrho,nspin,rho_new)
    USE math , ONLY: Norm
    implicit none
    INTEGER(I4B) :: nr, nrho, Nspin
    REAL(DP),intent(IN)     :: r_new(3,nr)
    REAL(DP),intent(OUT)    :: rho_new(nrho,Nspin)
    !> local
    ! REAL(DP) :: alpha,beta
    REAL(DP) :: a11,a22,a12,b1,b2
    REAL(DP) :: detA
    REAL(DP) :: new_max,old_max,delta_r(nr),delta_r1(nr)
    INTEGER(I4B) :: dimension_rho
    INTEGER(I4B) :: i
    !>> total charge
    REAL(DP) :: tele,tele_local
    REAL(DP) :: alpha2,beta2

    a11=0.d0
    a22=0.d0
    a12=0.d0
    b1=0.d0
    b2=0.d0
    delta_r=0.d0
    delta_r1=0.d0

    do i =1,nr
       a11=a11+sum((r1(:,i)-r2(:,i))**2)
       a22=a22+sum((r2(:,i)-r3(:,i))**2)
       a12=a12+sum( (r1(:,i)-r2(:,i)) * (r2(:,i)-r3(:,i)) )
       b1 =b1 +sum( (r_new(:,i)-r1(:,i)) * (r1(:,i)-r2(:,i)) )
       b2 =b2 +sum( (r_new(:,i)-r1(:,i)) * (r2(:,i)-r3(:,i)) )
    enddo
    detA=a11*a22-a12**2
    if(detA<1.d0.or.counter1<4)then
       do i =1,nr
          delta_r(i)= Norm(r_new(:,i)-r1(:,i))
          delta_r1(i)= Norm(r1(:,i)-r2(:,i))
       enddo
       old_max=sqrt(maxval(delta_r1))
       new_max=sqrt(maxval(delta_r))
       alpha=new_max/old_max
       beta=0.d0
       ! print*,'old_max,new_max',old_max,new_max
       if(old_max<1.d-5)alpha=0
       if(alpha>10)alpha=0
    else
       alpha2=(b1*a22-b2*a12)**2/detA**2
       beta2=(b2*a11-b1*a12)**2/detA**2
       alpha=sqrt(alpha2/(alpha2+beta2))
       if(alpha>10)alpha=0
       beta=1.d0-alpha**2
    endif
    alpha=0
    beta=0

    ! print *,'r1-r2',r1-r2
    ! print *,'r2-r3',r2-r3
#ifdef MPI
    ! if(parallel%isroot)print *,'alpha,beta',alpha,beta
    ! if(parallel%isroot)print *,'detA',detA
#else
    ! print *,'alpha,beta',alpha,beta
#endif
    ! print*,'detA',detA

    if(counter1>3)then
       dimension_rho = nrho*nspin
       CALL dcopy(dimension_rho, rho1,1,rho_new,1)
       CALL daxpy(dimension_rho, alpha, rho1,1,rho_new,1)
       CALL daxpy(dimension_rho, beta-alpha, rho2,1,rho_new,1)
       CALL daxpy(dimension_rho, -beta, rho3,1,rho_new,1)
       if(Lsrho_at)rho_new=rho_new+rho_at
    else
       dimension_rho = nrho*nspin
       CALL dcopy(dimension_rho, rho1,1,rho_new,1)
       CALL daxpy(dimension_rho, alpha, rho1,1,rho_new,1)
       !rho_new=rho1
       if(Lsrho_at)rho_new = rho_new + rho_at
    endif

    !> unified the charge_density
    where(rho_new<0)
       rho_new=0
    endwhere

#ifdef MPI
    tele_local=SUM(rho_new)*dvol
    CALL MPI_ALLREDUCE(tele_local,tele,1,MPI_REAL8,MPI_SUM,parallel%comm,mpinfo)
    if(parallel%isroot)print *,'ncharge,tele >',ncharge,tele
#else
    print *,"marvalus, ",dvol,sum(rho_new)
    tele=SUM(rho_new)*dvol
    print *,'ncharge,tele >',ncharge,tele
#endif
    rho_new=rho_new*ncharge/tele

    ! rho_new=rho1
    ! if(parallel%isroot)print *,'get_psi',rho_new(1:10,1)

  ENDSUBROUTINE get_rho

  SUBROUTINE get_psi(nrho,n_s,nspin,psi_new)
    implicit none
    INTEGER(I4B),intent(IN)     :: nrho,n_s,nspin
    REAL(DP),intent(OUT)    :: psi_new(nrho,n_s,nspin)
    !> local
    ! REAL(DP) :: alpha,beta
    ! REAL(DP) :: a11,a22,a12,b1,b2
    ! REAL(DP) :: detA
    INTEGER(I4B) :: dimension_psi
    INTEGER(I4B) :: i

    !print *,'alpha,beta',alpha,beta

    if(counter3>3)then
       dimension_psi = nrho*n_s*nspin
       CALL dcopy(dimension_psi, psi1,1,psi_new,1)
       CALL daxpy(dimension_psi, alpha,psi1,1,psi_new,1)
       CALL daxpy(dimension_psi, beta-alpha,psi2,1,psi_new,1)
       CALL daxpy(dimension_psi, -beta,psi3,1,psi_new,1)
    else
       psi_new=psi1
    endif

    ! psi_new=psi1
    ! print *,'get_psi',psi_new(1:10,1,1)

  ENDSUBROUTINE get_psi

  !>===============================================================================
  !> inherit by fft translation
  !> SUBROUTINE cal_trans_phase(nr,nspin,r_new,n1,n2,n3,ng1,ng2,ng3,gvec,trans_phase)
  !> SUBROUTINE get_new_rho_psi(nr,r_new,nrho,n1,n2,n3,nspin,rho_new,n_s,psi_new,gvec)
  !> SUBROUTINE store_rho_fft_trans(n_rho,Nspin,rho)
  !> SUBROUTINE store_rho_at_fft_trans(n_rho,Nspin,na,rho_in,rho_in2)
  !> SUBROUTINE store_r_fft_trans(nr,r)
  !> SUBROUTINE store_psi_fft_trans(n_rho,n_s,nspin,psi)

  SUBROUTINE cal_trans_phase(nr,nspin,r_new,n1,n2,n3,ng1,ng2,ng3,gvec,trans_phase)
    USE FOURIER , ONLY: FFT
    IMPLICIT NONE
    !> reciprocal grid
    INTEGER(I4B),intent(in) :: nr,nspin
    REAL(DP),intent(in)     :: r_new(3,nr)
    INTEGER(I4B),intent(in) :: ng1,ng2,ng3,n1,n2,n3
    REAL(DP),intent(in)     :: gvec(4,ng1*ng2*ng3)
    COMPLEX(DCP),intent(out)  :: trans_phase(ng1,ng2,ng3,nspin)
    !> local
    INTEGER(I4B) :: i,ig,j,i1,i2,i3
    REAL(DP)     :: len_trans(3) !> the distance of translation
    COMPLEX(DCP)  :: eigr
    COMPLEX(DCP)  :: rhoi_reci(ng1,ng2,ng3,nspin,nr),rho_at_reci(ng1,ng2,ng3,nspin)
    COMPLEX(DCP)  :: rhoi_over_rhoat(ng1,ng2,ng3)

    !> initialize the reciprocal array of rho_i, rho_atom and their quotient
    ! allocate(rhoi_reci(ng1,ng2,ng3,nspin,nr))
    ! allocate(rho_at_reci(ng1,ng2,ng3,nspin))
    trans_phase=(0,0)

    !> rho_atom for per atom in reciprocal
    do i=1,nr,1   !> per atom position
       do j=1,Nspin,1 !> per spin
          rhoi_reci(:,:,:,j,i)=FFT(reshape(rhoi1(:,j,i),(/n1,n2,n3/)))
       enddo
    enddo

    !> rho_atom for all atom in reciprocal
    do i=1,nspin,1 !> per spin
       rho_at_reci(:,:,:,i)=FFT( reshape(rho_at1(:,i),(/n1,n2,n3/)) )
    enddo
    ! do i=1,nspin,1 !> per spin
    !    rhoi_reci(:,:,:,1,i)=FFT( reshape(rho_at(:,i),(/n1,n2,n3/)) )
    ! enddo

   !> cal
    do i=1,nr,1   !> per atom position
    len_trans=r_new(:,i)-r1(:,i)
    ! print *,'len_trans',len_trans
    ! len_trans=1.d0
    do j=1,nspin,1
       rhoi_over_rhoat = rhoi_reci(:,:,:,j,i)/rho_at_reci(:,:,:,j)
       ig=0
       do i3=1,ng3,1
       do i2=1,ng2,1
       do i1=1,ng1,1  !> per reciprocal grid
          ig=ig+1
          eigr=exp(-imag*sum(gvec(1:3,ig)*len_trans))
          ! eigr=exp(-imag*2*pi*(real(i1-1)/real(n1)*len_trans(1) + real(i2-1)/real(n2)*len_trans(2) + real(i3-1)/real(n3)*len_trans(3)))
          trans_phase(i1,i2,i3,j)=trans_phase(i1,i2,i3,j)+rhoi_over_rhoat(i1,i2,i3)*eigr
       enddo
       enddo
       enddo  !> per reciprocal grid
    enddo
    enddo     !> per atom position
    ! trans_phase=rhoi_reci(:,:,:,1,:)/rho_at_reci

  ENDSUBROUTINE cal_trans_phase

  SUBROUTINE get_new_rho_psi(nr,r_new,nrho,n1,n2,n3,nspin,rho_new,n_s,psi_new,gvec)
    !> get the new rho and new psi by fourier translation
    USE Fourier, ONLY: FFT
    implicit none
    INTEGER(I4B),intent(IN) :: nr, nrho,n1,n2,n3, Nspin, n_s
    REAL(DP),intent(IN)     :: r_new(3,nr),gvec(4,nrho)
    REAL(DP),intent(OUT)    :: rho_new(nrho,Nspin)
    REAL(DP),intent(OUT)    :: psi_new(nrho,n_s,nspin)
    !> local
    INTEGER(I4B) :: ng1,ng2,ng3
    COMPLEX(DCP),allocatable :: trans_phase(:,:,:,:),rho1_reci(:,:,:)
    REAL(DP),allocatable    ::  nrho_l(:,:,:),nrho_new_l(:,:,:)
    INTEGER(I4B) :: i,j
    CHARACTER(len=10) :: nu

    !> initialize
    ng1=n1/2+1
    ng2=n2
    ng3=n3
    allocate(trans_phase(ng1,ng2,ng3,nspin))
    nu=''
    write(nu,'(I4)')counter1

    !> get translation phase
    CALL cal_trans_phase(nr,nspin,r_new,n1,n2,n3,ng1,ng2,ng3,gvec,trans_phase)
    open(1256,file='tran_phase')
    write(1256,*)trans_phase
    close(1256)

    !> cal new rho !>> nrho_l,nrho_new_l
    allocate(rho1_reci(ng1,ng2,ng3))
    allocate(nrho_l(n1,n2,n3))
    allocate(nrho_new_l(n1,n2,n3))
    do i=1,Nspin,1
       ! CALL dcopy(nrho,rho_at1(:,i),1,nrho_l,1)
       CALL dcopy(nrho,rho1(:,i),1,nrho_l,1)
    open(1256,file='nrho_l'//trim(adjustl(nu)))
    write(1256,*)n1,n2,n3
    write(1256,*)nrho_l
    close(1256)
       !> cal array2 which "FFT(array2)=FFT(array1)*trans"
       rho1_reci=FFT(nrho_l)
       rho1_reci=rho1_reci*trans_phase(:,:,:,i)
       nrho_new_l=FFT(rho1_reci)
    open(1256,file='nrho_new_l'//trim(adjustl(nu)))
    write(1256,*)n1,n2,n3
    write(1256,*)nrho_new_l
    close(1256)
       CALL dcopy(nrho,nrho_new_l,1,rho_new(:,i),1)
    enddo
    ! rho_new=rho1

    !> cal new psi
    ! do i=1,Nspin,1
    ! do j=1,n_s,1
    !    CALL dcopy(nrho,psi1(:,j,i),1,nrho_l,1)
    !    rho1_reci=FFT(nrho_l)
    !    rho1_reci=rho1_reci*trans_phase(:,:,:,i)
    !    nrho_new_l=FFT(rho1_reci)
    !    CALL dcopy(nrho,nrho_new_l,1,psi_new(:,j,i),1)
    ! enddo
    ! enddo
    psi_new=psi1
    deallocate(rho1_reci,nrho_l,nrho_new_l)

  ENDSUBROUTINE get_new_rho_psi

  SUBROUTINE store_rho_fft_trans(n_rho,Nspin,rho)
    implicit none
    INTEGER(I4B) :: n_rho,Nspin
    REAL(DP) :: rho(n_rho,Nspin)

    if(counter1==0)allocate(rho1(n_rho,Nspin))
    rho1=rho
    counter1=counter1+1
    if(counter1==1)Lsrho = .true.
    if(Lsrho.and.Lsr.and.Lspsi)Llastrho=.true.

  ENDSUBROUTINE store_rho_fft_trans

  SUBROUTINE store_rho_at_fft_trans(n_rho,Nspin,na,rho_in,rho_in2)
    implicit none
    INTEGER(I4B) :: n_rho,Nspin,na
    REAL(DP) :: rho_in(n_rho,Nspin)
    REAL(DP) :: rho_in2(n_rho,Nspin,na)


    if(.not.allocated(rho_at))then
       allocate(rho_at(n_rho,Nspin))
       allocate(rhoi(n_rho,Nspin,na))
       allocate(rho_at1(n_rho,Nspin))
       allocate(rhoi1(n_rho,Nspin,na))
    endif
    rho_at1=rho_at
    rhoi1=rhoi
    rho_at=rho_in
    rhoi=rho_in2
    Lsrho_at=.true.

  ENDSUBROUTINE store_rho_at_fft_trans

  SUBROUTINE store_r_fft_trans(nr,r)
    implicit none
    INTEGER(I4B) :: nr
    REAL(DP)     :: r(3,nr)

    if(counter2==0)allocate(r1(3,nr))
    r1=r
    counter2=counter2+1
    if(counter2==1)Lsr = .true.
    if(Lsrho.and.Lsr.and.Lspsi)Llastrho=.true.

  ENDSUBROUTINE store_r_fft_trans

  SUBROUTINE store_psi_fft_trans(n_rho,n_s,nspin,psi)
    implicit none
    !> array dimension
    INTEGER(I4B) :: n_rho,n_s,nspin
    !> psi
    REAL(DP) :: psi(n_rho,n_s,nspin)

    if(counter3==0)allocate(psi1(n_rho,n_s,nspin))
    ! print*,'shape psi',shape(psi)
    ! print*,'shape psi1',shape(psi1)
    psi1=psi
    counter3=counter3+1
    if(counter3==1)Lspsi = .true.
    if(Lsrho.and.Lsr.and.Lspsi)Llastrho=.true.

  ENDSUBROUTINE store_psi_fft_trans

ENDMODULE succeed
