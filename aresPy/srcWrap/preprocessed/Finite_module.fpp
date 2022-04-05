# 1 "Finite_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Finite_module.f90"
MODULE finite_module
!##########################################################!{{{
!* CREATED_TIME  : 2017-07-07
!* AUTHOR        : Qiang Xu
!* DESCRIPTION   :
!     ------
!* REFERENCES    :
!     ------
!* LOG           :
!##########################################################!}}}
   USE constants
   USE grid_module , ONLY : n1,n2,n3
   IMPLICIT NONE
!REAL(DP),ALLOCATABLE :: Laplcoe(:) !laplace coe
   REAL(DP),ALLOCATABLE :: Lapl(:,:)  !laplace total coe
!REAL(DP),ALLOCATABLE :: Gradcoe(:) !gradient coe
   REAL(DP),ALLOCATABLE :: Grad(:,:)  !total gradient coe
   REAL(DP) :: tBmat(3,3) ! for non-orthrog grids
!REAL(DP) :: lap_gap(3)
   INTEGER(I4B) :: lap_add(3)
!transform full
   INTEGER :: cell_mu(3,3)
   REAL(DP) :: cell_factor(3)

   REAL(DP),ALLOCATABLE :: wrap_real(:)

CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!###################################################################!
!       For complex wave function finite difference operator        !
!###################################################################!
!------------------------destroy  finite------------------------------
   SUBROUTINE destroy_finite()
      IMPLICIT NONE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ALLOCATED(Lapl)) DEALLOCATE(Lapl)
      IF(ALLOCATED(Grad)) DEALLOCATE(Grad)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_finite
!---------------------------initial finite----------------------------
   SUBROUTINE init_finite(h)
      USE parameters , ONLY : finite_order
      USE math , ONLY : finite_factor,finite_factor_new,inv_33,Det
      USE struct_module , ONLY : lat_mat,lat_para
      IMPLICIT NONE
!IN/OUT
      REAL(DP) :: h(3)  !gridsize
!LOCAL
      REAL(DP),ALLOCATABLE,DIMENSION(:) ::  &
            &   Laplcoe  &
            & , Gradcoe
      REAL(DP) :: Amat(3,3),Bmat(3,3)
      REAL(DP) :: lap_gap(3),factor(6)
      INTEGER(I4B) :: i
      INTEGER(I4B) :: norder
      REAL(DP) :: err(6)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL destroy_finite()
!exit
      IF(finite_order/=0)THEN
         norder=finite_order
      ELSE
         norder=8
         finite_order=8
      ENDIF
!
      ALLOCATE(Laplcoe(-norder:norder))
      ALLOCATE(Gradcoe(-norder:norder))
      ALLOCATE(Lapl(-norder:norder,6))
      ALLOCATE(Grad(-norder:norder,3))
!set data
!      CALL finite_factor(1,norder,Gradcoe)
!      CALL finite_factor(2,norder,Laplcoe)
!print*,'Grad coe'
!print*,Gradcoe
!print*,'Laplcoe'
!print*,Laplcoe
      CALL finite_factor_new(1,norder,Gradcoe)
      CALL finite_factor_new(2,norder,Laplcoe)
!print*,'our'
!print*,'Grad coe'
!print*,Gradcoe
!print*,'Laplcoe'
!print*,Laplcoe


!transform to cartensien coor
      Amat(:,1)=lat_mat(:,1)/lat_para(1) !SQRT(SUM(lat_mat(:,1)**2))
      Amat(:,2)=lat_mat(:,2)/lat_para(2)!SQRT(SUM(lat_mat(:,2)**2))
      Amat(:,3)=lat_mat(:,3)/lat_para(3)!SQRT(SUM(lat_mat(:,3)**2))
      Bmat=inv_33(Amat)
      tBmat=TRANSPOSE(Bmat)
!tBmat=Bmat
!grad
      Grad(:,1)=Gradcoe(:)/h(1)
      Grad(:,2)=Gradcoe(:)/h(2)
      Grad(:,3)=Gradcoe(:)/h(3)
!Lapl
!total coefficient in grid
!full
      CALL trans_mat_full(lat_mat,factor,cell_mu,lap_gap,err)
      DO i=1,3
         Lapl(:,i)=factor(i)*Laplcoe(:)/h(i)**2
         Lapl(:,i+3)=factor(i+3)*Laplcoe(:)/lap_gap(i)**2
      ENDDO
!deallocate
      DEALLOCATE(Laplcoe,Gradcoe)
!print*,'gap'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_finite
!-----------------------trans_mat_full_new----------------------------
   subroutine trans_mat_full(mat,factor,int_miu,lad_gap,err)!{{{
!##########################################################
!* CREATED_TIME  : 2018-08-31
!* AUTHOR        : Xuecheng Shao & Qiang Xu
!* CHANGE        : Xuecheng Shao
!* ADD           : Xuecheng Shao
!* DESCRIPTION   :
!     ------
!* REFERENCES    :
!     1.PHYSICAL REVIEW B 78, 075109 (2008)
!       "Real-space pseudopotential method for first principles calculations
!       of general periodic and partially periodic systems"
!* LOG           :
!     2015-05-10
!* LAST MODIFIED : 2015-05-10 08:39:18 PM
!##########################################################
      USE constants,   ONLY : DP,I4B
      use math
!USE Lapack_module , ONLY : invmat_real
      USE grid_module ,ONLY : gap
      USE struct_module , ONLY : lat_para
      implicit none
!IN/OUT
      REAL(DP),INTENT(IN) :: mat(3,3) !lattice matrix
      REAL(DP),INTENT(OUT) :: factor(6),lad_gap(3)
      INTEGER(I4B),INTENT(OUT) :: int_miu(3,3)
      REAL(DP),INTENT(OUT) :: err(6)
!LOCAL
      INTEGER(I4B)                     :: i,j,k,i1,i2,i3,ucount
!INTEGER(I4B)                     :: offset(3),bound(3),n_one,idx(27)
      INTEGER(I4B) :: intmiu(3,27),n_one
      INTEGER(I4B)                     :: sum_miu(3) !,int_tmp
      real(DP),dimension(3,3)     :: A,F_mat,inva,M_mat,invM,new_r,miu
      real(DP)                    :: f_vec(3),vec_tmp(3,27),scal_one(27)
      real(DP)                    :: b_vec(3),r(3),lgap(3)
      real(DP)                    :: vol,rmod
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do i=1,3
         A(:,i)=mat(:,i)/lat_para(i)
      ENDDO
      invA=inv_33(A)
      F_mat=matmul(invA,transpose(invA))
      f_vec(1)=F_mat(1,2)+F_mat(2,1)
      f_vec(2)=F_mat(1,3)+F_mat(3,1)
      f_vec(3)=F_mat(2,3)+F_mat(3,2)
      k=0 !count total number
      ucount=0 !count useful number
      do i3=-1,1
      do i2=-1,1
      do i1=-1,1
         k=k+1
         intmiu(1,k)=i1
         intmiu(2,k)=i2
         intmiu(3,k)=i3
!delete the cell lattice vector
         if ((abs(i1)+abs(i2)+abs(i3))<=1) then
            vec_tmp(:,k)=0.d0
            scal_one(k)=10000.d0
         else
!r(1:3)=i3*a(:,3)+i2*a(:,2)+i1*a(:,1)
            ucount=ucount+1
            r(1:3)=gap(1)*i1*A(:,1)+gap(2)*i2*A(:,2)+gap(3)*i3*A(:,3)
            vec_tmp(:,k)=r
            scal_one(k)=SQRT(SUM(r**2))
         endif
      enddo
      enddo
      enddo
      n_one=27    ! n_one=27
!
!print*,'lgap'
!print*,scal_one
!print*,vec_tmp(1,:)
!sort by length
      CALL realInt_sort(n_one,scal_one,intmiu,vec_tmp)
!      print*,'sort'
!      print*,scal_one
!print*,vec_tmp(1,:)
!STOP
!sort the directions by length
!call sort_id(scal_one,n_one,idx)
!select the direction
!changed by Qiang Xu
!offset(:)=2
!bound=3
      findit:DO i=1,n_one
         new_r(:,1)=vec_tmp(:,i)
         int_miu(1,:)=intmiu(:,i)
         lad_gap(1)=scal_one(i)
      DO j=1,n_one
         new_r(:,2)=vec_tmp(:,j)
         int_miu(2,:)=intmiu(:,j)
         lad_gap(2)=scal_one(j)
      DO k=1,n_one
         new_r(:,3)=vec_tmp(:,k)
         int_miu(3,:)=intmiu(:,k)
         lad_gap(3)=scal_one(k)
!get real miu
         miu(:,:)=REAL(int_miu(:,:),DP)
         r = miu(1,1)*A(:,1) + miu(1,2) * A(:,2) + miu(1,3)*A(:,3)
         rmod=SQRT(SUM(r*r))
         IF(rmod<xtiny) CYCLE
         miu(1,:)=miu(1,:)/rmod
         r = miu(2,1)*A(:,1) + miu(2,2) * A(:,2) + miu(2,3)*A(:,3)
         rmod=SQRT(SUM(r*r))
         IF(rmod<xtiny) CYCLE
         miu(2,:)=miu(2,:)/rmod
         r = miu(3,1)*A(:,1) + miu(3,2) * A(:,2) + miu(3,3)*A(:,3)
         rmod=SQRT(SUM(r*r))
         IF(rmod<xtiny) CYCLE
         miu(3,:)=miu(3,:)/rmod
!get M_mat=2*M_mu
         do i1=1,3
            M_mat(1,i1)=2*miu(i1,1)*miu(i1,2)
            M_mat(2,i1)=2*miu(i1,1)*miu(i1,3)
            M_mat(3,i1)=2*miu(i1,2)*miu(i1,3)
         enddo
         
         vol=ABS(Det(M_mat))
!print*,int_miu(:,:)
!print*,miu(:,:)
!print*,'vol',vol
!pause
         IF(vol<1e-7) CYCLE
!satisfy
         EXIT findit

      ENDDO
      ENDDO
      ENDDO findit
!print*,miu(1,:)
!print*,miu(2,:)
!print*,miu(3,:)
!print*,vol
!STOP

!get factors ang gaps

!print*,i,j,k
!!print*,
!vol=ABS(Det(M_mat))
      IF(vol>xtiny)THEN
!print *,"M",M_mat
!call gauss_u(M_mat,f_vec,b_vec)
         invM=inv_33(M_mat)
         b_vec(:)=MATMUL(invM,f_vec)
         DO i=1,3
            IF(ABS(b_vec(i))<5e-7) b_vec(i)=0._DP
         ENDDO
!print *,"b",b_vec
         factor(1)=F_mat(1,1) - sum(b_vec(:)*(miu(:,1)**2))
         factor(2)=F_mat(2,2) - sum(b_vec(:)*(miu(:,2)**2))
         factor(3)=F_mat(3,3) - sum(b_vec(:)*(miu(:,3)**2))
         factor(4)=b_vec(1)
         factor(5)=b_vec(2)
         factor(6)=b_vec(3)
!gap
!DO i=1,3
!   !r(:)= gap(1)*a(:,1)*int_miu(i,1) + &
!   !    & gap(2)*a(:,2)*int_miu(i,2) + &
!   !    & gap(3)*a(:,3)*int_miu(i,3)
!   r(:)=new_r(:,i)
!   lad_gap(i)=SQRT(SUM(r*r))
!ENDDO

      ELSE
!volum is zero
         print*,'transmat_full : invert is not exit!!!'
         STOP
!CALL transmat(factor,lad_gap,lad)
      ENDIF
!err analysis
      err(1)=factor(1)*gap(1)**2
      err(2)=factor(2)*gap(2)**2
      err(3)=factor(3)*gap(3)**2
      err(4)=factor(4)*lad_gap(1)**2
      err(5)=factor(5)*lad_gap(2)**2
      err(6)=factor(6)*lad_gap(3)**2

!where(factor<0.1) factor=0.d0
!print *,"fac",factor
!pause
!print *,"mu",(int_miu(i,:),i=1,3)
!print *,"mu",int_miu(:,:)
!print*,gap(1:3)
!DO i=1,3
!   r(:)= gap(1)*a(:,1)*int_miu(i,1) + &
!       & gap(2)*a(:,2)*int_miu(i,2) + &
!       & gap(3)*a(:,3)*int_miu(i,3)
!   lad_gap(i)=SQRT(SUM(r*r))
!ENDDO
!lad_gap(3)=lad_gap(2)
!print*,'mu'
!print*,det(REAL(int_miu,8))
!print*,int_miu(1,:)
!print*,int_miu(2,:)
!print*,int_miu(3,:)
!!!int_miu(3,1)=0
!!!int_miu(3,2)=1
!!!int_miu(3,3)=1
!!!!print*,int_miu(:,:)
!factor(1:4)=sum(factor(1:4))/4
!gap(1:3)=(sum(gap(1:3))+lad_gap(1))/4
!lad_gap(4)=gap(1)
!print*,'factor'
!!factor(4)=0._DP
!print *, factor(1:6)
!print*,'gap'
!!lad_gap(4)=gap(1)
!print*,gap(:)
!print*,lad_gap(:)
!stop
!!
!pause
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   END subroutine trans_mat_full
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!----------------------------PARTING LINE-----------------------------
   SUBROUTINE cmplx_keop(Uk,Ik,Ts)
!#################################################################!
!ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
!#################################################################!
      USE struct_module, ONLY : recip_lat
      USE m_time_evaluate, only: memory_sum,memory_free,filename

      USE grid_module , ONLY :KPT,global_n1,global_n2,global_n3,&
           & global_n,n1,n2,n3
      USE smpi_math_module

      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)  :: Ik
      COMPLEX(DCP),INTENT(IN)  :: Uk(n1,n2,n3)
      COMPLEX(DCP),INTENT(OUT) :: Ts(n1,n2,n3)
!
      REAL(DP) :: kmod2_2
      COMPLEX(DCP) :: Ts1(n1,n2,n3)
      COMPLEX(DCP) :: Ts2(n1,n2,n3)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Ts=0._DP
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_keop
!####################################################################!
!       For real wave function finite difference operator            !
!*Author : Qiang Xu                                                  !
!*Date   : 2017/08/28                                                !
!####################################################################!
   SUBROUTINE real_comm(ifun)
      USE parameters, ONLY : norder=>finite_order
      USE grid_module , ONLY : global_n,global_n1,global_n2,grid
      USE smpi_math_module , ONLY: parallel,diff_map,MPI_REAL8&
           & ,mpinfo,set_wrap_sph_pbc_ata_real,Lall_grid
      IMPLICIT NONE
!IN/OUT
      REAL(DP),INTENT(IN)  :: ifun(:)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL real_comm_clean()
!> comm setting init
      IF(Lall_grid)THEN
         ALLOCATE(wrap_real(grid%OneDlength+1))
         wrap_real(grid%OneDlength+1)=0.d0
         CALL MPI_ALLGATHERV(ifun,parallel%mygrid_range(3),&
              & MPI_REAL8,wrap_real,parallel%recvcounts&
              & ,parallel%displs, MPI_REAL8, parallel%commx&
              & ,mpinfo)

      ELSE
! ALLOCATE(wrap_real(global_n1*global_n2* &
!  & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1)))
! !> set wrap_grid
! CALL set_wrap_grid_pbc_ata(ifun,wrap_real,global_n ,global_n1, global_n2)
         ALLOCATE(wrap_real(diff_map%boundary1d(2)-diff_map%boundary1d(1)+2))
         wrap_real(diff_map%boundary1d(2)-diff_map%boundary1d(1)+2)=0.d0
!> set wrap_grid
         CALL set_wrap_sph_pbc_ata_real(ifun,wrap_real)
      ENDIF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_comm
!--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_comm_clean()
      IMPLICIT NONE
      IF(allocated(wrap_real)) DEALLOCATE(wrap_real)
   ENDSUBROUTINE real_comm_clean
!--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_nabla1_3d(ifun,norder,derf)
!
      IMPLICIT NONE
!
      REAL(DP),INTENT(IN)  :: ifun(:,:,:)
      INTEGER(I4B),INTENT(IN)  :: norder 
      REAL(DP),INTENT(OUT) :: derf(:,:)
!
      COMPLEX(DCP) :: vec(-norder+1:norder+n1,-norder+1:norder+n2,-norder+1:norder+n3)
      INTEGER(I4B) :: i,ish,ix,iy,iz,ip
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO iz = 1, n3
!! xm ym zm
         vec(1:n1,1:n2,iz)           = ifun(1:n1,1:n2,iz)
!! xl ym zm
         vec(-norder+1:0,1:n2,iz)    = ifun(n1-norder+1:n1,1:n2,iz)
!! xr ym zm
         vec(n1+1:n1+norder,1:n2,iz) = ifun(1:norder,1:n2,iz)
      ENDDO
!! xm ym zm
      vec(1:n1,1:n2,1:n3)=ifun(1:n1,1:n2,1:n3)
!! xl ym zm
      vec(-norder+1:0,1:n2,1:n3)=ifun(n1-norder+1:n1,1:n2,1:n3)
!! xr ym zm
      vec(n1+1:n1+norder,1:n2,1:n3)=ifun(1:norder,1:n2,1:n3)
!! xa yl zm
      vec(:,-norder+1:0,1:n3)=vec(:,n2-norder+1:n2,1:n3)
!! xa yr zm
      vec(:,n2+1:n2+norder,1:n3)=vec(:,1:norder,1:n3)
!! xa ya zl
      vec(:,:,-norder+1:0)=vec(:,:,n3-norder+1:n3)
!! xa ya zr
      vec(:,:,n3+1:n3+norder)=vec(:,:,1:norder)
      derf(:,:) = 0.d0
!df/dx
      ip=0
      DO iz=1,n3
      DO iy=1,n2
      DO ix=1,n1
         ip=ip+1
         DO ish=-norder,norder
            derf(ip,1)=derf(ip,1)+Grad(ish,1)*vec(ix+ish,iy,iz)
            derf(ip,2)=derf(ip,2)+Grad(ish,2)*vec(ix,iy+ish,iz)
            derf(ip,3)=derf(ip,3)+Grad(ish,3)*vec(ix,iy,iz+ish)
         ENDDO
      ENDDO
      ENDDO
      ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_nabla1_3d
!--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_nabla2_3d(func,norder,ofun)!{{{
!##########################################################
!* CREATED_TIME  : 2013-03-12
!* AUTHOR        : Yanchao Wang
!* CHANGE        : Xuecheng Shao
!* ADD           : Xuecheng Shao
!* DESCRIPTION   :
!     ------
!* REFERENCES    :
!     ------
!* LOG           :
!     2015-05-08 :
!* LAST MODIFIED : 2015-05-10 07:47:12 PM
!##########################################################
      use constants
!USE parameters,        only : norder=>finite_order
      USE grid_module , ONLY : gap
!
      implicit none
!IN/OUT
      REAL(DP),INTENT(IN) :: func(n1,n2,n3)
      INTEGER(I4B),INTENT(IN) :: norder
      REAL(DP),INTENT(OUT) :: ofun(n1,n2,n3)
!LOCAL
      REAL(DP) :: vec(-norder+1:norder+n1,-norder+1:norder+n2,-norder+1:norder+n3)
!real(dp)             :: coeke(-norder:norder,6)
!
      INTEGER(I4B)              :: i,ish,ix,iy,iz
      real(dp)             :: gaps(6)
!>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      do iz = 1, n3
!! xm ym zm
         vec(1:n1,1:n2,iz)           = func(1:n1,1:n2,iz)
!! xl ym zm
         vec(-norder+1:0,1:n2,iz)    = func(n1-norder+1:n1,1:n2,iz)
!! xr ym zm
         vec(n1+1:n1+norder,1:n2,iz) = func(1:norder,1:n2,iz)
      enddo

      do iz = 1, n3
! xa yl zm
         vec(:,-norder+1:0,iz)=vec(:,n2-norder+1:n2,iz)
! xa yr zm
         vec(:,n2+1:n2+norder,iz)=vec(:,1:norder,iz)
      enddo

!! xa ya zl
      do iz = -norder+1, 0
         vec(:,:, iz)=vec(:,:, iz+n3)
      enddo

!! xa ya zr
      do iz = 1, norder
         vec(:,:, n3+iz)=vec(:,:, iz)
      enddo

      do iz=1,n3
         do iy=1,n2
            do ix=1,n1
               ofun(ix,iy,iz) = 0.d0
               do ish = -norder,norder
                  ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                   &  Lapl(ish,1)*(vec(ix+ish,iy,iz))+ &
                   &  Lapl(ish,2)*(vec(ix,iy+ish,iz))+ &
                   &  Lapl(ish,3)*(vec(ix,iy,iz+ish))+ &
                   &  Lapl(ish,4)*(vec(ix+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish))+&
                   &  Lapl(ish,5)*(vec(ix+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish))+&
                   &  Lapl(ish,6)*(vec(ix+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish))
               ENDDO
            enddo
         enddo
      ENDDO
!>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
   ENDSUBROUTINE real_nabla2_3d!}}}
!--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_nabla1_1d(ifun,derf,mgfun)
      USE parameters , ONLY : norder=>finite_order
      USE grid_module, ONLY : n,global_n1,global_n2,global_n3

      USE smpi_math_module

      IMPLICIT NONE
!IN/OUT
      REAL(DP),INTENT(IN)  :: ifun(n)
      REAL(DP),INTENT(OUT) :: derf(n,3)
      REAL(DP),OPTIONAL :: mgfun(n)
!LOCAL
      INTEGER(I4B) :: ix,iy,iz,ip,ish,id1,id2,id3
      LOGICAL :: lmgfun

      INTEGER(I4B) :: shiftn

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lmgfun=PRESENT(mgfun)
!by finite difference
!zero boundary

      CALL real_comm(ifun)

!transform to cartensian
      DO ip=1,n

!set 0
         derf(ip,:)=0._DP
!
         DO ish=-norder,norder
!def id
            id1=diff_map%local_map1d(1,ish,ip)
            id2=diff_map%local_map1d(2,ish,ip)
            id3=diff_map%local_map1d(3,ish,ip)

            derf(ip,1)= derf(ip,1) + Grad(ish,1)*wrap_real(id1)
            derf(ip,2)= derf(ip,2) + Grad(ish,2)*wrap_real(id2)
            derf(ip,3)= derf(ip,3) + Grad(ish,3)*wrap_real(id3)

         ENDDO

 
         derf(ip,:)=MATMUL(tBmat,derf(ip,:))
      ENDDO
!mod of grad
      IF(lmgfun)THEN
         mgfun(:)=SQRT(derf(:,1)**2+derf(:,2)**2+derf(:,3)**2)
      ENDIF
!clean
      CALL real_comm_clean()
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_nabla1_1d
!--------------------------PARTING LINE-------------------------------
   SUBROUTINE real_nabla2_1d(ifun,ofun)
!##########################################################
!*For    : \nabla^2(f)
!*Author : Xuecheng Shao & Qiang Xu
!*Date   : 2017/08/08
!##########################################################
      USE parameters , ONLY : norder=>finite_order

      USE grid_module, ONLY : n,global_n1,global_n2,global_n3
      USE smpi_math_module

      IMPLICIT NONE
!
      REAL(DP),INTENT(IN)  :: ifun(n) 
      REAL(DP),INTENT(OUT) :: ofun(n) 
!

      INTEGER(I4B) :: ip,ish,ix,iy,iz &
         &, id1,id2,id3,id4,id5,id6   &
         &, shiftn
      INTEGER(I4B) :: i,j,k,in
      INTEGER(I4B) :: d_z,y_down

!>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
!zero boundary

      CALL real_comm(ifun)
!
! call start_time('diff',.true.)
      DO ip=1,n

         ofun(ip) = 0._DP

         DO ish = -norder,norder


            ofun(ip)=ofun(ip)+&
             &  Lapl(ish,1)*wrap_real(diff_map%local_map1d(1,ish,ip))+ &
             &  Lapl(ish,2)*wrap_real(diff_map%local_map1d(2,ish,ip))+ &
             &  Lapl(ish,3)*wrap_real(diff_map%local_map1d(3,ish,ip))+ &
             &  Lapl(ish,4)*wrap_real(diff_map%local_map1d(4,ish,ip))+ &
             &  Lapl(ish,5)*wrap_real(diff_map%local_map1d(5,ish,ip))+ &
             &  Lapl(ish,6)*wrap_real(diff_map%local_map1d(6,ish,ip))
         ENDDO 

      ENDDO
      CALL real_comm_clean()

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! call end_time('diff',.true.)
! call write_time('diff',.true.)
   ENDSUBROUTINE real_nabla2_1d
!----------------------------PARTING LINE-----------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE finite_module
