MODULE CG_RELAX
  USE CONSTANTS
  ! USE PARAMETERS , only : iounits
  ! use system_info

  ! use struct_module
  ! use dyn_module
  ! USE SMPI_MATH_MODULE

  IMPLICIT NONE
  REAL(DP)  :: MAXFORCE , MAXSTRESS
  REAL(DP)  :: TSTRESS(3,3) = 0.D0
  REAL(DP)  :: Latmark(3,3) = 1.d0
  type  :: lattice
     real(dp)  :: A(3,3)
     real(dp)  :: B(3,3)
     real(dP)  :: anorm(3),bnorm(3)
     real(dp)  :: omega
  end type lattice
CONTAINS
  subroutine cg_relax_vasp_interface(nstep,lfopt,lopt,lexceed )!{{{
    use struct_module , ONLY : struct,lat_mat,recip_lat,natom,naty,energy
    use parameters  , only : potim,ediffg,pstress,lcell,Lpbc
#ifdef MPI
    USE SMPI_MATH_MODULE , ONLY: parallel,mpi_real8,mpinfo
#endif
    use dyn_module
    use math ,        only: inv_33,lat2matrix
    use write_module , only : write_poscar
    logical,intent(inout)  :: lopt , lfopt
    integer(i4b),intent(in) :: nstep
    integer(i4b)    :: na
    integer(i4b)    :: iflag
    real(dp),save        :: toten=0.d0 , toteng=0.d0
    real(dp)        :: fact,factsi,dismax,ebreak!,ediffg
    real(dp),save   :: e1test
    real(dp)        :: d2sif(3,3)
    real(dp)        :: tsif(3,3)
    real(dp),allocatable :: posion(:,:) ,posioc(:,:), force(:,:)
    real(dp)        :: lat_para(6)
    logical,save         :: lstop2
    real(dp)        :: econv = 0.d0
    integer(i4b)    :: iu0, iu6
    integer(i4b)    :: ia,ityp,istep,i,j
    integer(i4b)    :: ii,jj,kk
    logical :: lexist
    character(5) :: txt
    type(lattice)  :: latt
    logical         :: lexceed 
    !----------------------------------------------------------------
#ifdef MPI
    if(parallel%isroot)then
#else
    if(.true.)then
#endif
      ! if ( nstep == 1 ) then
      !  open(8888,file="file_8888")
      !  open(11111,file="file_11111")
      ! else
      !  open(8888,file="file_8888",access='append')
      !  open(11111,file="file_11111",access='append')
      ! end if
    endif
    if(.not.Lpbc)struct%stress = 0.d0 
    na = natom
    iflag = 1   
    if (nstep == 1) iflag = 0
    !ediffg = -1E-3_dp 
    dyn%ediffg = ediffg
    dyn%pstress = pstress
    econv = 0.d0
    ebreak = max(abs(ediffg),abs(econv))
    iu0 = 6
    iu6 = 8
    
    dyn%potim = 0.1
    dyn%potim = potim
    fact = 0.d0
    fact = 10.d0*dyn%potim*CONST_E_AU/CONST_MU_AU * 1E-10
    dyn%nfree = 1
    ebreak = 0.25d0 * ediffg
    print *,"nstep",NSTEP
    if (nstep == 1) then 
      call destroy_dyn(dyn)
      call create_dyn(na,dyn)
      dyn%d2c = 0.d0
      dyn%d2 = 0.d0
      dyn%d3 = 0.d0
      e1test = 0.d0
      !lstop2 = .false.
    end if 
    if (allocated(force)) deallocate(force)
    allocate(force(3,na))
    dyn%posion = struct%pos
    if ( nstep == 1) dyn%posioc = struct%pos

    force (:,:) = struct%forces(:,:) * force2ev
    lstop2 = .true.
    do ityp = 1,naty
      do ia = struct%eleid(ityp),struct%eleid(ityp+1)-1
        dyn%d2c(1,ia) = force(1,ia) * fact
        dyn%d2c(2,ia) = force(2,ia) * fact
        dyn%d2c(3,ia) = force(3,ia) * fact
        if (sqrt(force(1,ia)**2+force(2,ia)**2+force(3,ia)**2) &
            & > abs(dyn%ediffg)) lstop2 = .false.
      end do
    end do

    dismax = 0.d0
    latt%A = lat_mat / ang2bohr
    call lattic(latt)
    factsi = 0.d0 
    factsi = 10*dyn%potim*CONST_E_AU/CONST_MU_AU/NA*1E-10_dp
    toten = energy(1)
    if(Lpbc)then
    do ii = 1,3
    do jj = 1,3
      d2sif(ii,jj) = -1.d0 * struct%stress(ii,jj)*factsi*autogpa
    end do
    !d2sif(ii,ii)=d2sif(ii,ii)-dyn%pstress/(CONST_E_AU*1e22_dp)*latt%omega*factsi
    d2sif(ii,ii)=d2sif(ii,ii)-dyn%pstress/ 10.d0 * factsi
    end do
    else
    do ii = 1,3
    do jj = 1,3
      d2sif(ii,jj) = 0.d0
    end do
    !d2sif(ii,ii)=d2sif(ii,ii)-dyn%pstress/(CONST_E_AU*1e22_dp)*latt%omega*factsi
    d2sif(ii,ii)=0.d0
    end do
    endif

    if ( factsi /= 0 ) then 
       do ii = 1 , 3
       do jj = 1 , 3 
          if ( factsi /= 0 ) then 
             if ( ( abs( d2sif (jj,ii) / factsi / na ) ) > abs(dyn%ediffg) ) lstop2 = .false.
          end if
       enddo
       enddo
    endif

    call ioncgr(iflag,na,toten,latt%A,latt%B,dyn%nfree, &
        &  dyn%posion,dyn%posioc,fact,dyn%d2c,factsi,d2sif,dyn%d2,dyn%d3,dismax,iu6,iu0,ebreak,dyn%ediffg, &
        &  e1test,lstop2)

#ifdef MPI
    if (parallel%isroot) then
#else
    if (.true.) then
#endif
       !inquire(file = 'xdatcar',exist=lexist)
       ! write(8888,'(A18)') 'local optimisation'
       ! write(8888,*) '1.d0'
       ! write(8888,'(1x,3f12.6)') (latt%A(1,i)/1, latt%A(2,i)/1,latt%A(3,i)/1,i = 1,3)
       ! write(8888,'(6x,A10)') struct%elements(:)
       ! write(8888,'(I8)')  natom
       ! write(8888,'(A,I6)') 'Direct configuration=',nstep
       ! write(8888,1221) ((dyn%posioc(I,j),i=1,3),j=1,na)
       ! write(8888,'(A,I6)')
    endif

    lopt = .false.

    call check_opt(na,lfopt)

    if ( iflag == 1 )then
       lopt = ( abs( toten - toteng) < dyn%ediffg )
       toteng = toten
    endif

    if ( iflag == 2 ) lopt = .true.
    if ( iflag == 1 ) lopt = ( lopt .or. ( abs(e1test) < 0.1d0 * dyn%ediffg ) )

1221 FORMAT(1X,3F12.8)

    latt%B = transpose(inv_33(latt%A))
    call lat2matrix(lat_para,latt%A,2) 

    struct%pos = dyn%posion
    lat_mat = latt%A * ang2bohr
    recip_lat = 2.d0*pi*transpose(inv_33(lat_mat))

    call check_distance(lcell , struct%pos,lexceed )

#ifdef MPI
    if (parallel%isroot) then
#else
    if (.true.) then
#endif
       ! write(11111,*)
       ! write(11111,'(10x,2A30)') 'real lattice','reciprocal lattice'
       ! WRITE(11111,'(3(F15.7),10x,3(F15.7))') latt%A(:,1) , 2.d0 * pi * latt%B(:,1)
       ! WRITE(11111,'(3(F15.7),10x,3(F15.7))') latt%A(:,2) , 2.d0 * pi * latt%B(:,2)
       ! WRITE(11111,'(3(F15.7),10x,3(F15.7))') latt%A(:,3) , 2.d0 * pi * latt%B(:,3)
       ! !WRITE(IOUNITS(11),'(3(F15.7),10x,3(F15.7))') struct%lat_mat(:,1) , struct%recip_lat(:,1)
       ! !WRITE(IOUNITS(11),'(3(F15.7),10x,3(F15.7))') struct%lat_mat(:,2) , struct%recip_lat(:,2)
       ! !WRITE(IOUNITS(11),'(3(F15.7),10x,3(F15.7))') struct%lat_mat(:,3) , struct%recip_lat(:,3)
       ! write(11111,*)
       ! write(11111,'(10x,2A30)') 'Lattice Parameters','Cell Angles'
       ! write(11111,'(20x,A3,3x,F13.6,16x,A8,F13.6)') 'a =',lat_para(1) , 'alpha =',lat_para(4) * 57.29578
       ! write(11111,'(20x,A3,3x,F13.6,16x,A8,F13.6)') 'b =',lat_para(2) , 'beta =', lat_para(5) * 57.29578
       ! write(11111,'(20x,A3,3x,F13.6,16x,A8,F13.6)') 'c =',lat_para(3) , 'gamma =',lat_para(6) * 57.29578
       ! write(11111,*)
       ! write(11111,*)
    endif


#ifdef MPI
    call MPI_BCAST(struct%pos,na*3,MPI_REAL8,PARALLEL%ROOTID,PARALLEL%COMM,MPINFO)
    call MPI_BCAST(LAT_MAT,9,MPI_REAL8,PARALLEL%ROOTID,PARALLEL%COMM,MPINFO)
    call MPI_BCAST(RECIP_LAT,9,MPI_REAL8,PARALLEL%ROOTID,PARALLEL%COMM,MPINFO)
#endif

#ifdef MPI
    if (parallel%isroot) then
#else
    if (.true.) then
#endif
       do ityp = 1,naty
          do ia = struct%eleid(ityp),struct%eleid(ityp+1)-1
             ! write(11111,'(20x,A6,I4,3x,3F13.6)') 'posion',ia,struct%pos(:,ia)
          enddo
       enddo
    endif

    ! struct_local%pos = struct%pos(:,struct_local%atomid_local(1) : struct_local%atomid_local(2) ) 

#ifdef MPI
    if (parallel%isroot) then
#else
    if (.true.) then
#endif
       CALL write_poscar('CONTCAR', lat_mat,struct%eleid,struct%pos)
       ! close(11111)
       ! close(8888)
    endif
#ifdef MPI
    ! call MPI_BARRIER(PARALLEL%COMM,MPINFO)
#endif

    RETURN
  End Subroutine cg_relax_vasp_interface!}}}

  SUBROUTINE check_opt(na,flag)!{{{
    use struct_module , only : struct
    use parameters , only : igoal , tolf , tolp , pstress,Lpbc
#ifdef MPI
    use smpi_math_module , only: parallel
#endif
    !
    IMPLICIT NONE
    LOGICAL                           :: flag
    integer(i4b), intent(in) :: na
    REAL(DP)       :: rtmp
    REAL(DP)       :: press
    REAL(DP),save      :: savetolf(3)=1E30
    REAL(DP),save      :: savetolE(3)=1E30
    !LOGICAL,optional                  :: reflag
    LOGICAL :: flagE=.false.
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    if (igoal == 1 ) then
       maxforce = maxval(abs(struct%forces))
       maxstress = 0.d0
    elseif (igoal == 2 ) then
       maxforce = 0.d0
       maxstress= maxval(abs(struct%stress+tstress)*latmark)
    elseif (igoal == 3 ) then
       maxforce = maxval(abs(struct%forces(:,:)))
       maxstress= maxval(abs(struct%stress+tstress)*latmark)
    endif
    maxforce = maxforce * force2ev
    if(Lpbc)then
       maxstress = maxstress * autogpa
    else
       maxstress = 0.d0!maxstress * autogpa
    endif
    maxstress = abs( maxstress - pstress / 10.d0 )
#ifdef MPI
    if (parallel%isroot) print *,'maxforce',maxforce
    if (parallel%isroot) print *,'maxstress',maxstress
#else
    print *,'maxforce',maxforce
    print *,'maxstress',maxstress
#endif

    if ( maxforce < tolf .and. maxstress < tolp ) then
      flag = .true.
    else
      flag = .false.
    endif

  END SUBROUTINE check_opt !}}}

  subroutine lattic(mylatt)!{{{
    IMPLICIT NONE
    type(lattice)  :: mylatt
    real(dp)  :: omega
    integer(i4b)  ::  i , j
    integer(i4b)  ::  sum

    CALL EXPRO(Mylatt%B(1:3,1),Mylatt%A(1:3,2),Mylatt%A(1:3,3))
    CALL EXPRO(Mylatt%B(1:3,2),Mylatt%A(1:3,3),Mylatt%A(1:3,1))
    CALL EXPRO(Mylatt%B(1:3,3),Mylatt%A(1:3,1),Mylatt%A(1:3,2))

    Omega =Mylatt%B(1,1)*Mylatt%A(1,1)+Mylatt%B(2,1)*Mylatt%A(2,1) &
         &      +Mylatt%B(3,1)*Mylatt%A(3,1)

    DO I=1,3
    DO J=1,3
       Mylatt%B(I,J)=Mylatt%B(I,J)/Omega
    ENDDO
    ENDDO

    DO I=1,3
       Mylatt%ANORM(I)=SQRT(SUM(Mylatt%A(:,I)*Mylatt%A(:,I)))
       Mylatt%BNORM(I)=SQRT(SUM(Mylatt%B(:,I)*Mylatt%B(:,I)))
    ENDDO
    Mylatt%Omega=Omega
    RETURN
  END SUBROUTINE lattic!}}}

  SUBROUTINE EXPRO(H,U1,U2)!{{{
    implicit none
    real(dp)  ::  H(3),U1(3),U2(3)

    H(1)=U1(2)*U2(3)-U1(3)*U2(2)
    H(2)=U1(3)*U2(1)-U1(1)*U2(3)
    H(3)=U1(1)*U2(2)-U1(2)*U2(1)

    RETURN
  END SUBROUTINE EXPRO!}}}


  subroutine check_distance(l , pos , lwrong)!{{{
    use math
    implicit none
    real(dp),intent(in)  ::  l , pos(:,:)
    real(dp)             ::  r(size(pos,2)), rmax, pos1(size(pos,1),size(pos,2))
    integer(i4b)         ::  i
    logical              ::  lwrong
    !-----------------------------------------------
    lwrong = .false.
    do i = 1, 3
      pos1(i,:) = abs(pos(i,:) - 0.5d0)
    end do

    do i = 1, size(pos,2)
      r(i) =  norm(pos1(:,i) )
    end do

    rmax = maxval(r)

    if ( rmax >= 0.4999999999 ) then
      lwrong = .true.
    end if

  end subroutine check_distance !}}}

END MODULE CG_RELAX

