# 1 "Matvec_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Matvec_module.f90"
MODULE matvec_module
  USE constants
CONTAINS
  !-------------------------------------------------------------
  SUBROUTINE cmatvec(veff,Ik,p,q,dimen)
     USE parameters , ONLY : LDG
     USE finite_module , ONLY : ksKEperiod_op
     USE grid_module , ONLY : n1,n2,n3
     use m_time_evaluate, only: memory_sum,memory_free,filename

     ! USE smpi_math_module, ONLY: parallel
     USE smpi_math_module , ONLY : parallel,mpi_complex16,mpi_sum,mpinfo,mpi_real8

     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     INTEGER(I4B),INTENT(IN) :: Ik
     REAL(DP),INTENT(IN) :: veff(:,:,:)
     COMPLEX(DP),INTENT(IN) :: p(dimen)
     COMPLEX(DP),INTENT(OUT) :: q(dimen)
     !LOCAL
     COMPLEX(DP),DIMENSION(n1,n2,n3) :: thrq  &
                                    &   ,Tsp
     INTEGER(I4B) :: I,Ix,Iy,Iz

     INTEGER(I4B) :: In,j,k
     INTEGER(I4B) :: d_z,y_down
     ! COMPLEX(DCP)  :: temp_global(global_n)
     ! REAL(DP)  :: temp_globalr(global_n)

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("cmatvec",real(size(thrq),DP)*DP+size(Tsp)*DP)
     !transfrom vec to 3D

     thrq=cmplx(0.d0,0.d0)
     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
     in=0 !parallel%mygrid_range(1)-1
     do i=ix,n1
        in=in+1
        thrq(i,iy,iz)=p(in)
     enddo
     d_z=ceiling((n3-iz)/(n3-iz+1.d0))
     do j=iy+1,n2*d_z
        do i=1,n1
           in=in+1
           thrq(i,j,iz)= p(in)
        enddo
     enddo
     do k=2,n3-1
        do j=1,n2
           do i=1,n1
              in=in+1
              thrq(i,j,k)= p(in)
           enddo
        enddo
     enddo
     ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
     y_down=(1-iy-1)*d_z+iy+1
     iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
     iz=n3 !(parallel%mygrid_range(2)-1)/n1/n2+1
     do j=y_down,iy-1
        do i=1,n1
           in=in+1
           thrq(i,j,iz)= p(in)
        enddo
     enddo
     do i=1,ix
        in=in+1
        thrq(i,iy,iz)= p(in)
     enddo
# 84 "Matvec_module.f90"

     !nonlocal part Sum_lm(D0_lm |\beta_lm><\beta_lm| )
     IF(LDG)THEN
         CALL nlocmatvec_dg(Ik,p,q)
     ELSE
         CALL nlocmatvec(Ik,p,q)
     ENDIF
     !KE+local part
     CALL ksKEperiod_op(thrq,Ik,Tsp)
     ! CALL MPI_BARRIER(parallel%commx,mpinfo)
     thrq(:,:,:)=Tsp(:,:,:)+thrq(:,:,:)*veff(:,:,:)

     ! print*,'dimen',dimen,parallel%mygrid_range(3),global_n
     ! print*,'recvcounts',parallel%recvcounts
     ! print*,'displs',parallel%displs
     ! print*,'size',size(veff),'rank',parallel%myid
     ! open(1111+parallel%myid,file="q"//filename(10:11))
     ! call MPI_ALLGATHERV(q,parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! open(1111+parallel%myid,file="veff"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(veff(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_REAL8,temp_globalr,parallel%recvcounts&
     !      &,parallel%displs,MPI_REAL8&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_globalr)
     ! close(1111+parallel%myid)
     ! open(1111+parallel%myid,file="Tsp"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(Tsp(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! stop

     !out q : thrq  -> vec q

     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
     in=0 !parallel%mygrid_range(1)-1
     do i=ix,n1
        in=in+1
        q(in)=q(in)+thrq(i,iy,iz)
     enddo
     do j=iy+1,n2*d_z
        do i=1,n1
           in=in+1
           q(in)=q(in)+thrq(i,j,iz)
        enddo
     enddo
     do k=2,n3-1
        do j=1,n2
           do i=1,n1
              in=in+1
              q(in)=q(in)+thrq(i,j,k)
           enddo
        enddo
     enddo
     ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
     y_down=(1-iy-1)*d_z+iy+1
     iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
     iz=n3 !(parallel%mygrid_range(2)-1)/n1/n2+1
     do j=y_down,iy-1
        do i=1,n1
           in=in+1
           q(in)=q(in)+thrq(i,j,iz)
        enddo
     enddo
     do i=1,ix
        in=in+1
        q(in)=q(in)+thrq(i,iy,iz)
     enddo
# 178 "Matvec_module.f90"
     ! open(1111+parallel%myid,file="q"//filename(10:11))
     ! call MPI_ALLGATHERV(q,parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%comm,mpinfo)
     ! stop 'parallel'
     ! call memory_free("cmatvec",real(size(thrq),DP)*DP+size(Tsp)*DP)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE cmatvec
  !--------------------------------------------------------------
  SUBROUTINE nlocmatvec(Ik,p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot  !,nlpot_Sfact
     USE grid_module ,ONLY : dvol
     use m_time_evaluate ,only:memory_sum,memory_free

     USE smpi_math_module , ONLY : parallel,mpi_complex16,mpi_sum,mpinfo

     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN) :: Ik
     COMPLEX(DP),INTENT(IN)     :: p(:)
     COMPLEX(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     COMPLEX(DP) :: dots(max_nproj,natom),tmp0

     COMPLEX(DP) :: dots_local(max_nproj,natom)

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("nlocmatvec",real(size(dots),DP)*DP)
     !
     q(:)=0.d0
     IF(max_nproj==0) RETURN
     !



     dots_local(:,:)=0.d0

     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0) CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0

              if( nlpot(Ia)%npts==0 )cycle

              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+CONJG(nlpot(Ia)%proj_phs(Ip,Ipj,Ik))*p(Id)
              ENDDO
              !save dots




              dots_local(:,Ia)=dots_local(:,Ia)+tmp0*psp(Ity)%D0(:,Ipj)

           ENDDO
        ENDDO
     ENDDO

     CALL MPI_ALLREDUCE(dots_local,dots,size(dots),mpi_complex16,mpi_sum,parallel%commx,mpinfo)

     !scale
     dots(:,:)=dots(:,:)*dvol !* nlpot_Sfact(Ik)
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj_phs(Ip,Ipj,Ik)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     call memory_free("nlocmatvec",real(size(dots),DP)*DP)
  ENDSUBROUTINE nlocmatvec
  !------------------------double grid----------------------------
  SUBROUTINE nlocmatvec_dg(Ik,p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot
     USE grid_module ,ONLY : dvol
     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN) :: Ik
     COMPLEX(DP),INTENT(IN)     :: p(:)
     COMPLEX(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     COMPLEX(DP) :: dots(max_nproj,natom),tmp0
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !
     q(:)=0.d0
     IF(max_nproj==0) RETURN
     !
     dots(:,:)=0.d0
     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0) CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+CONJG(nlpot(Ia)%proj_phs_dg(Ip,Ipj,Ik))*p(Id)
              ENDDO
              !save dots
              dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)
           ENDDO
        ENDDO
     ENDDO
     !scale
     dots(:,:)=dots(:,:)*dvol
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj_phs(Ip,Ipj,Ik)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE nlocmatvec_dg
  !------------------------for real matver------------------------
  SUBROUTINE rmatvec(veff,p,q,dimen)
     USE grid_module , ONLY : n1,n2,n3
     USE finite_module , ONLY : realKE
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     REAL(DP),INTENT(IN) :: veff(:,:,:)
     REAL(DP),INTENT(IN) :: p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
     !
     REAL(DP),DIMENSION(n1,n2,n3)   :: thrq  &
                                    &   ,Tsp
     INTEGER(I4B) :: I,Ix,Iy,Iz,t1,t2,t3,t4,t5,t6
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL system_clock(t1)
     !transfrom vec to 3D
     I=0
     DO Iz=1,n3
     DO Iy=1,n2
     DO Ix=1,n1
        I=I+1
        thrq(Ix,Iy,Iz)=p(I)
     ENDDO
     ENDDO
     ENDDO
     !nonlocal part
     CALL system_clock(t2)
     CALL nlocmatvec_r(p,q)
     CALL system_clock(t3)
     !=====================================================
     !##OUTPUT NONLOCAL TERM TO TEST ITS CORRECTNESS
     !open(unit=521,file="nl_Vrho_period")
     !write(521,*)n1,n2,n3
     !write(521,*)q
     !close(521)
     !=====================================================
     !
     CALL realKE(thrq,Tsp)
     CALL system_clock(t4)
     thrq(:,:,:)=Tsp(:,:,:)+thrq(:,:,:)*veff(:,:,:)
     !out q : thrq  -> vec q
     CALL system_clock(t5)
     I=0
     DO Iz=1,n3
     DO Iy=1,n2
     DO Ix=1,n1
        I=I+1
        q(I)=q(I)+thrq(Ix,Iy,Iz)
     ENDDO
     ENDDO
     ENDDO
     CALL system_clock(t6)
     !============================================
     !print *,'transform time -->',(t2-t1)/10000.d0
     !print *,'nonloc time -->',(t3-t2)/10000.d0
     !print *,'kinetic time -->',(t4-t3)/10000.d0
     !print *,'VlocalPsi time -->',(t5-t4)/10000.d0
     !print *,'add time -->',(t6-t5)/10000.d0
     !print *,'TOTAL time -->',(t6-t1)/10000.d0
     !============================================
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE rmatvec
  !--------------------------------------------------------------
  SUBROUTINE nlocmatvec_r(p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot  !,nlpot_Sfact
     USE grid_module ,ONLY : dvol
     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(IN)     :: p(:)
     REAL(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     REAL(DP) :: dots(max_nproj,natom),tmp0
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !
     q(:)=0.d0
     IF(max_nproj==0)RETURN
     !
     dots(:,:)=0.d0
     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0)CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+ nlpot(Ia)%proj(Ip,Ipj) *p(Id)
              ENDDO
              !save dots
              ! dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)
              !>xlt: for off-diagnal element of D0
              dots(:,Ia)=dots(:,Ia)+tmp0*psp(Ity)%D0(:,Ipj)
           ENDDO
        ENDDO
     ENDDO
     !scale
     dots(:,:)=dots(:,:)*dvol
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj(Ip,Ipj)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE nlocmatvec_r
  !--------------------------------------------------------------
  SUBROUTINE rmatvec_new(mat,p,q,dimen)
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     REAL(DP),INTENT(IN) :: mat(dimen,dimen),p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     q(:)=MATMUL(mat,p)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE rmatvec_new

!-------------------------------------------------------------------
!------------------------for real matver------------------------
  SUBROUTINE ISO_rmatvec(veff_3D,p,q,dimen)
  !SUBROUTINE ISO_rmatvec(veff,p,q,dimen)
     USE parameters , ONLY : LDG
     USE grid_module , ONLY : rho_calc,n1,n2,n3,n
     USE finite_module , ONLY : ISO_realKE

     USE smpi_math_module

     USE m_time_evaluate, ONLY: memory_sum,memory_free
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen



     REAL(DP),INTENT(IN) :: veff_3D(dimen)

     REAL(DP),INTENT(IN) :: p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
     !> local






     REAL(DP),DIMENSION(n1,n2,n3) :: thrq,thrtsp
     REAL(DP)                     :: Tsp(dimen)

     INTEGER(I4B) :: I,Ix,Iy,Iz
     REAL(DP)     :: gridp(n)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



     call memory_sum('ISO_rmatvec',real(size(thrq),DP)*2*DP+size(Tsp)*DP+n*DP)







     IF(LDG)THEN
        CALL nlocmatvec_iso_dg(p,q)
     ELSE
        CALL nlocmatvec_iso(p,q)
     ENDIF
# 522 "Matvec_module.f90"

     !> calculate the kinetic part



     CALL ISO_realKE(p,Tsp)



# 544 "Matvec_module.f90"
     do I=1,parallel%mygrid_range(3)
        Ix=rho_calc%x(I)
        Iy=rho_calc%y(I)
        Iz=rho_calc%z(I)
        q(I)=q(I)+Tsp(I)+p(I)*Veff_3D(I)
     enddo



! #endif



     call memory_free('ISO_rmatvec',real(size(thrq),DP)*2*DP+size(Tsp)*DP+n*DP)

     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE ISO_rmatvec
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE ISO_grid2sphere(grid_V,grid_V1,sphere_V)
    USE parameters ,  ONLY : NSPIN
    USE grid_module , ONLY : rho_calc,n1,n2,n3
    IMPLICIT NONE
    INTEGER(I4B)         :: i, j
    REAL(DP),INTENT(IN)  :: grid_V(1:,1:,1:),grid_V1(1:)
    REAL(DP),INTENT(OUT) :: sphere_V(1:)
    sphere_V=0.d0
    !DO j=1,NSPIN,1
    DO i=1,rho_calc%OneDLength
      sphere_V(i)=grid_V(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))+ &
      &grid_V1(rho_calc%x(i)+n1*(rho_calc%y(i)-1)+n1*n2*(rho_calc%z(i)-1))
    ENDDO
    !ENDDO
  END SUBROUTINE ISO_grid2sphere
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE ISO_sphere2grid(sphere_V,grid_V)
    USE parameters ,  ONLY : NSPIN
    USE grid_module , ONLY : rho_calc
    IMPLICIT NONE
    INTEGER(I4B)         :: i, j
    REAL(DP),INTENT(IN)  :: sphere_V(:)
    REAL(DP),INTENT(OUT) :: grid_V(:,:,:)
    grid_V=0.d0
    !DO j=1,NSPIN,1
    DO i=1,rho_calc%OneDLength
      grid_V(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))=sphere_V(i)
    ENDDO
    !ENDDO
  END SUBROUTINE ISO_sphere2grid
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE nlocmatvec_iso(p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot  !,nlpot_Sfact
     USE grid_module ,ONLY : dvol

     USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo

     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(IN)     :: p(:)
     REAL(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     REAL(DP) :: dots(max_nproj,natom),tmp0

     REAL(DP) :: dots_local(max_nproj,natom)

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !
     q(:)=0.d0
     IF(max_nproj==0)RETURN
     !



     dots_local(:,:)=0.d0

     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0)CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0

              if( nlpot(Ia)%npts==0 )cycle

              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id_iso(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+nlpot(Ia)%proj(Ip,Ipj)*p(Id)
              ENDDO
              !save dots





              dots_local(:,Ia)=dots_local(:,Ia)+tmp0*psp(Ity)%D0(:,Ipj)

           ENDDO
        ENDDO
     ENDDO

     CALL MPI_ALLREDUCE(dots_local,dots,size(dots),mpi_real8,mpi_sum,parallel%comm,mpinfo)

     !scale
     dots(:,:)=dots(:,:)*dvol
     ! print *,dots
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id_iso(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj(Ip,Ipj)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE nlocmatvec_iso
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE nlocmatvec_iso_dg(p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot  !,nlpot_Sfact
     USE grid_module ,ONLY : dvol

     USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo

     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(IN)     :: p(:)
     REAL(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     REAL(DP) :: dots(max_nproj,natom),tmp0

     REAL(DP) :: dots_local(max_nproj,natom)

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !
     q(:)=0.d0
     IF(max_nproj==0)RETURN
     !



     dots_local(:,:)=0.d0

     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0)CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0

              if( nlpot(Ia)%npts==0 )cycle

              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id_iso(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+nlpot(Ia)%proj_dg(Ip,Ipj)*p(Id)
              ENDDO
              !save dots



              dots_local(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)

           ENDDO
        ENDDO
     ENDDO

     CALL MPI_ALLREDUCE(dots_local,dots,size(dots),mpi_real8,mpi_sum,parallel%comm,mpinfo)

     !scale
     dots(:,:)=dots(:,:)*dvol
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id_iso(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj_dg(Ip,Ipj)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE nlocmatvec_iso_dg
  !-----------------------PARTING-LINE--------------------------
  !--------------------------------------------------------------

  SUBROUTINE cmatvec_band(veff,Ik,p,q,dimen)
     USE parameters , ONLY : LDG
     USE finite_module , ONLY : ksKEperiod_op_band
     USE grid_module , ONLY :n1=> global_n1,n2=>global_n2,n3=>global_n3
     use m_time_evaluate, only: memory_sum,memory_free,filename
     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     INTEGER(I4B),INTENT(IN) :: Ik
     REAL(DP),INTENT(IN) :: veff(:,:,:)
     COMPLEX(DP),INTENT(IN) :: p(dimen)
     COMPLEX(DP),INTENT(OUT) :: q(dimen)
     !LOCAL
     COMPLEX(DP),DIMENSION(n1,n2,n3) :: thrq  &
                                    &   ,Tsp
     INTEGER(I4B) :: I,Ix,Iy,Iz
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("cmatvec",real(size(thrq),DP)*DP+size(Tsp)*DP)
     !transfrom vec to 3D
     I=0
     DO Iz=1,n3
     DO Iy=1,n2
     DO Ix=1,n1
        I=I+1
        thrq(Ix,Iy,Iz)=p(I)
     ENDDO
     ENDDO
     ENDDO

     !nonlocal part Sum_lm(D0_lm |\beta_lm><\beta_lm| )
     CALL nlocmatvec_band(Ik,p,q)
     !KE+local part
     CALL ksKEperiod_op_band(thrq,Ik,Tsp)
     ! CALL MPI_BARRIER(parallel%commx,mpinfo)
     thrq(:,:,:)=Tsp(:,:,:)+thrq(:,:,:)*veff(:,:,:)

     ! print*,'dimen',dimen,parallel%mygrid_range(3),global_n
     ! print*,'recvcounts',parallel%recvcounts
     ! print*,'displs',parallel%displs
     ! print*,'size',size(veff),'rank',parallel%myid
     ! open(1111+parallel%myid,file="q"//filename(10:11))
     ! call MPI_ALLGATHERV(q,parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! open(1111+parallel%myid,file="veff"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(veff(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_REAL8,temp_globalr,parallel%recvcounts&
     !      &,parallel%displs,MPI_REAL8&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_globalr)
     ! close(1111+parallel%myid)
     ! open(1111+parallel%myid,file="Tsp"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(Tsp(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! stop

     !out q : thrq  -> vec q
     I=0
     DO Iz=1,n3
     DO Iy=1,n2
     DO Ix=1,n1
        I=I+1
        q(I)=q(I)+thrq(Ix,Iy,Iz)
     ENDDO
     ENDDO
     ENDDO
     ! open(1111+parallel%myid,file="q"//filename(10:11))
     ! call MPI_ALLGATHERV(q,parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%comm,mpinfo)
     ! stop 'parallel'
     ! call memory_free("cmatvec",real(size(thrq),DP)*DP+size(Tsp)*DP)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmatvec_band
  !--------------------------------------------------------------
  SUBROUTINE nlocmatvec_band(Ik,p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot  !,nlpot_Sfact
     USE grid_module ,ONLY : dvol
     use m_time_evaluate ,only:memory_sum,memory_free
     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN) :: Ik
     COMPLEX(DP),INTENT(IN)     :: p(:)
     COMPLEX(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     COMPLEX(DP) :: dots(max_nproj,natom),tmp0
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("nlocmatvec",real(size(dots),DP)*DP)
     !
     q(:)=0.d0
     IF(max_nproj==0) RETURN
     !
     dots(:,:)=0.d0
     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0) CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+CONJG(nlpot(Ia)%proj_phs(Ip,Ipj,Ik))*p(Id)
              ENDDO
              !save dots
              ! dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)
              dots(:,Ia)=dots(:,Ia)+tmp0*psp(Ity)%D0(:,Ipj)
           ENDDO
        ENDDO
     ENDDO
     !scale
     dots(:,:)=dots(:,:)*dvol !* nlpot_Sfact(Ik)
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj_phs(Ip,Ipj,Ik)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     call memory_free("nlocmatvec",real(size(dots),DP)*DP)
   ENDSUBROUTINE nlocmatvec_band

  !-------------------------------------------------------------
  SUBROUTINE rmatvec_gamma(veff,Ik,p,q,dimen)
     USE parameters , ONLY : LDG
     USE finite_module , ONLY : ksKEperiod_op_gamma
     USE grid_module , ONLY : n1,n2,n3
     use m_time_evaluate, only: memory_sum,memory_free,filename

     ! USE smpi_math_module, ONLY: parallel
     USE smpi_math_module , ONLY : parallel,mpi_complex16,mpi_sum,mpinfo,mpi_real8

     IMPLICIT NONE
     !
     INTEGER(I4B),INTENT(IN) :: dimen
     INTEGER(I4B),INTENT(IN) :: Ik
     REAL(DP),INTENT(IN) :: veff(:,:,:)
     REAL(DP),INTENT(IN) :: p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
     !LOCAL
     REAL(DP),DIMENSION(n1,n2,n3) :: thrq  &
                                    &   ,Tsp
     INTEGER(I4B) :: I,Ix,Iy,Iz

     INTEGER(I4B) :: In,j,k
     INTEGER(I4B) :: d_z,y_down
     ! COMPLEX(DCP)  :: temp_global(global_n)
     ! REAL(DP)  :: temp_globalr(global_n)

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("cmatvec",real(size(thrq),DP)*DP+size(Tsp)*DP)
     !transfrom vec to 3D

     thrq=cmplx(0.d0,0.d0)
     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
     in=0 !parallel%mygrid_range(1)-1
     do i=ix,n1
        in=in+1
        thrq(i,iy,iz)=p(in)
     enddo
     d_z=ceiling((n3-iz)/(n3-iz+1.d0))
     do j=iy+1,n2*d_z
        do i=1,n1
           in=in+1
           thrq(i,j,iz)= p(in)
        enddo
     enddo
     do k=2,n3-1
        do j=1,n2
           do i=1,n1
              in=in+1
              thrq(i,j,k)= p(in)
           enddo
        enddo
     enddo
     ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
     y_down=(1-iy-1)*d_z+iy+1
     iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
     iz=n3 !(parallel%mygrid_range(2)-1)/n1/n2+1
     do j=y_down,iy-1
        do i=1,n1
           in=in+1
           thrq(i,j,iz)= p(in)
        enddo
     enddo
     do i=1,ix
        in=in+1
        thrq(i,iy,iz)= p(in)
     enddo
# 983 "Matvec_module.f90"

     !nonlocal part Sum_lm(D0_lm |\beta_lm><\beta_lm| )
     CALL nlocmatvec_gamma(Ik,p,q)
     !KE+local part
     CALL ksKEperiod_op_gamma(thrq,Ik,Tsp)
     ! CALL MPI_BARRIER(parallel%commx,mpinfo)
     thrq(:,:,:)=Tsp(:,:,:)+thrq(:,:,:)*veff(:,:,:)

     !out q : thrq  -> vec q

     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
     in=0 !parallel%mygrid_range(1)-1
     do i=ix,n1
        in=in+1
        q(in)=q(in)+thrq(i,iy,iz)
     enddo
     do j=iy+1,n2*d_z
        do i=1,n1
           in=in+1
           q(in)=q(in)+thrq(i,j,iz)
        enddo
     enddo
     do k=2,n3-1
        do j=1,n2
           do i=1,n1
              in=in+1
              q(in)=q(in)+thrq(i,j,k)
           enddo
        enddo
     enddo
     ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
     y_down=(1-iy-1)*d_z+iy+1
     iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
     iz=n3 !(parallel%mygrid_range(2)-1)/n1/n2+1
     do j=y_down,iy-1
        do i=1,n1
           in=in+1
           q(in)=q(in)+thrq(i,j,iz)
        enddo
     enddo
     do i=1,ix
        in=in+1
        q(in)=q(in)+thrq(i,iy,iz)
     enddo
# 1040 "Matvec_module.f90"
     ! open(1111+parallel%myid,file="q"//filename(10:11))
     ! call MPI_ALLGATHERV(q,parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%comm,mpinfo)
     ! stop 'parallel'
     ! call memory_free("cmatvec",real(size(thrq),DP)*DP+size(Tsp)*DP)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rmatvec_gamma
  !--------------------------------------------------------------
  SUBROUTINE nlocmatvec_gamma(Ik,p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot  !,nlpot_Sfact
     USE grid_module ,ONLY : dvol
     use m_time_evaluate ,only:memory_sum,memory_free

     USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo

     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN) :: Ik
     REAL(DP),INTENT(IN)     :: p(:)
     REAL(DP),INTENT(OUT)    :: q(:)
     !LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     REAL(DP) :: dots(max_nproj,natom),tmp0

     REAL(DP) :: dots_local(max_nproj,natom)

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("nlocmatvec",real(size(dots),DP)*DP)
     !
     q(:)=0.d0
     IF(max_nproj==0) RETURN
     !



     dots_local(:,:)=0.d0

     !all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0) CYCLE
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0

              if( nlpot(Ia)%npts==0 )cycle

              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !<proj_lm|p>
                 tmp0=tmp0+nlpot(Ia)%proj(Ip,Ipj)*p(Id)
              ENDDO
              !save dots




              dots_local(:,Ia)=dots_local(:,Ia)+tmp0*psp(Ity)%D0(:,Ipj)

           ENDDO
        ENDDO
     ENDDO

     CALL MPI_ALLREDUCE(dots_local,dots,size(dots),mpi_real8,mpi_sum,parallel%commx,mpinfo)

     !scale
     dots(:,:)=dots(:,:)*dvol !* nlpot_Sfact(Ik)
     !Multiply the nonlocal vectors
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
                 !q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj(Ip,Ipj)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     call memory_free("nlocmatvec",real(size(dots),DP)*DP)
   ENDSUBROUTINE nlocmatvec_gamma
ENDMODULE matvec_module
