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
   !COMPLEX(DP) :: KEgradcore(3,3)
   !REAL(DP) :: lap_gap(3)
   INTEGER(I4B) :: lap_add(3)
   !transform full
   INTEGER :: cell_mu(3,3)
   REAL(DP) :: cell_factor(3)

   COMPLEX(DCP),allocatable :: wrap_box(:,:,:),fun_global(:,:,:)
   COMPLEX(DCP),allocatable :: fun_1d(:),wrap_box1d(:)
   REAL(DP),allocatable :: wrap_box_real(:,:,:),fun_global_real(:,:,:)
   REAL(DP),allocatable :: fun_1d_real(:),wrap_box1d_real(:)

CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !###################################################################!
   !       For complex wave function finite difference operator        !
   !###################################################################!
   !---------------------------initial finite----------------------------
   SUBROUTINE init_finite(norder,h)
      USE parameters , ONLY : finite_order,Grad_order
      USE math , ONLY : finite_factor,finite_factor_new,inv_33,Det
      USE struct_module , ONLY : lat_mat,lat_para
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B) :: norder
      REAL(DP) :: h(3)  !gridsize
      !LOCAL
      REAL(DP),ALLOCATABLE,DIMENSION(:) ::  &
            &   Laplcoe  &
            & , Gradcoe
      REAL(DP) :: Amat(3,3),Bmat(3,3)
      REAL(DP) :: lap_gap(3),factor(6),err(6)
      INTEGER(I4B) :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL destroy_finite()
      !exit
      !
      ALLOCATE(Laplcoe(-norder:norder))
      ALLOCATE(Gradcoe(-norder:norder))
      ALLOCATE(Lapl(-norder:norder,6))
      ALLOCATE(Grad(-norder:norder,3))
      !set data
      CALL finite_factor_new(1,norder,Gradcoe)
      CALL finite_factor_new(2,norder,Laplcoe)


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
      !full finite difference
      !CALL trans_mat_full(lat_mat,cell_factor,cell_mu,lap_gap)
      CALL trans_mat_full(lat_mat,factor,cell_mu,lap_gap,err)
      DO i=1,3
         Lapl(:,i)=factor(i)*Laplcoe(:)/h(i)**2
         Lapl(:,i+3)=factor(i+3)*Laplcoe(:)/lap_gap(i)**2
      ENDDO
      !deallocate
      DEALLOCATE(Laplcoe,Gradcoe)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_finite
   !------------------------destroy  finite------------------------------
   SUBROUTINE destroy_finite()
     USE m_time_evaluate, ONLY: memory_free
      IMPLICIT NONE
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ALLOCATED(Lapl))THEN
         call memory_free('finite_Lapl',real(size(Lapl),DP)*DP)
         DEALLOCATE(Lapl)
      ENDIF
      IF(ALLOCATED(Grad))THEN
         call memory_free('finite_Drad',real(size(Grad),DP)*DP)
         DEALLOCATE(Grad)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_finite
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
            IF(ABS(b_vec(i))<xtiny) b_vec(i)=0._DP
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
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE cmplx_nabla2(ifun,norder,ofun)
      !##########################################################
      !*For    : -1/2*\nabla^2(f)
      !*Author : Xuecheng Shao & Qiang Xu
      !*Date   : 2017/08/08
      !##########################################################

    USE smpi_math_module, ONLY: parallel,start_time,end_time,diff_map
      ! use m_time_evaluate, only: filename
      ! USE grid_module , ONLY :global_n
      ! USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per

      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: norder
      COMPLEX(DP),INTENT(IN)     :: ifun(:,:,:)
      COMPLEX(DP),INTENT(OUT)    :: ofun(:,:,:)
      !



      !
      INTEGER(I4B) :: i,ish,ix,iy,iz
      REAL(DP)     :: kcoe

      INTEGER(I4B) :: j,k
      ! COMPLEX(DCP) :: temp_global(global_n)
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>

      call start_time('nabla2',.true.)

      kcoe = SUM( Lapl(0,1:3) )
      !$OMP PARALLEL PRIVATE(iz,iy,ix,ish)
      !$OMP DO
# 379 "Finite_module.f90"
      !OMP END DO
      !
      !$OMP DO
# 400 "Finite_module.f90"
      ix=(parallel%mygrid_range(2)-1)/n1/n2
      iy=(parallel%mygrid_range(1)-1)/n1/n2
      d_z=ceiling((ix-iy)/(ix-iy+1.d0))
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      ! iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
      do i=ix,n1
         ofun(i,iy,1) = 0.d0
         DO ish = -norder,norder
            ofun(i,iy,1)=ofun(i,iy,1)+&
                 Lapl(ish,1)*(wrap_box1d(diff_map%local_map1d(i+ish,iy,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d(diff_map%local_map1d(i,iy+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d(diff_map%local_map1d(i,iy,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         ofun(i,j,1) = 0
         DO ish = -norder,norder
            ofun(i,j,1)=ofun(i,j,1)+&
                 Lapl(ish,1)*(wrap_box1d(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         ofun(i,j,k) = 0
         DO ish = -norder,norder
            ofun(i,j,k)=ofun(i,j,k)+&
                 Lapl(ish,1)*(wrap_box1d(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      ! iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         ofun(i,j,n3) = 0
         DO ish = -norder,norder
            ofun(i,j,n3)=ofun(i,j,n3)+&
                 Lapl(ish,1)*(wrap_box1d(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      do i=1,ix
         ofun(i,iy,n3) = 0
         DO ish = -norder,norder
            ofun(i,iy,n3)=ofun(i,iy,n3)+&
                 Lapl(ish,1)*(wrap_box1d(diff_map%local_map1d(i+ish,iy,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d(diff_map%local_map1d(i,iy+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d(diff_map%local_map1d(i,iy,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d(diff_map%local_map1d(i+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo

     ! print*,'kcoe',kcoe
     ! open(1111+parallel%myid,file="Ts1"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(ofun(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! stop
      !OMP END DO

      call end_time('nabla2',.true.)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_nabla2
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE cmplx_nabla1(ifun,norder,derf)

      USE smpi_math_module, ONLY: parallel,start_time,end_time,diff_map
      USE Grid_module, ONLY:global_n1,global_n2

      !
      IMPLICIT NONE
      !
      COMPLEX(DP),INTENT(IN)  :: ifun(:,:,:)
      INTEGER(I4B),INTENT(IN)  :: norder
      COMPLEX(DP),INTENT(OUT) :: derf(:,:,:,:)
      !



      INTEGER(I4B) :: i,ish,ix,iy,iz

      INTEGER(I4B) :: j,k,in,i_temp
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !$OMP PARALLEL PRIVATE(iz,iy,ix,ish)
      !$OMP DO
# 550 "Finite_module.f90"

      !df/dx
      !$OMP DO
# 565 "Finite_module.f90"
      call start_time('nabla1',.true.)
      derf=(0.d0,0.d0)
      ix=(parallel%mygrid_range(2)-1)/n1/n2
      iy=(parallel%mygrid_range(1)-1)/n1/n2
      d_z=ceiling((ix-iy)/(ix-iy+1.d0))
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0 !parallel%mygrid_range(1)-1
      do i=ix,n1
         DO ish=-norder,norder
            derf(1,i,iy,1)=derf(1,i,iy,1)+Grad(ish,1)*wrap_box1d(diff_map%local_map1d(i+ish,iy,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,iy,1)=derf(1,i,iy,1)+Grad(ish,1)*wrap_box(i_temp,iy,iz)
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         DO ish=-norder,norder
            derf(1,i,j,1)=derf(1,i,j,1)+Grad(ish,1)*wrap_box1d(diff_map%local_map1d(i+ish,j,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,j,1)=derf(1,i,j,1)+Grad(ish,1)*wrap_box(i_temp,j,iz)
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         DO ish=-norder,norder
            derf(1,i,j,k)=derf(1,i,j,k)+Grad(ish,1)*wrap_box1d(diff_map%local_map1d(i+ish,j,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,j,k)=derf(1,i,j,k)+Grad(ish,1)*wrap_box(i_temp,j,iz)
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         DO ish=-norder,norder
            derf(1,i,j,n3)=derf(1,i,j,n3)+Grad(ish,1)*wrap_box1d(diff_map%local_map1d(i+ish,j,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,j,n3)=derf(1,i,j,n3)+Grad(ish,1)*wrap_box(i_temp,j,iz)
         ENDDO
      enddo
      enddo
      do i=1,ix
         DO ish=-norder,norder
            derf(1,i,iy,n3)=derf(1,i,iy,n3)+Grad(ish,1)*wrap_box1d(diff_map%local_map1d(i+ish,iy,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,iy,n3)=derf(1,i,iy,n3)+Grad(ish,1)*wrap_box(i_temp,iy,iz)
         ENDDO
      enddo

      !$OMP END DO

      !df/dy
      !$OMP DO
# 639 "Finite_module.f90"
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0 !parallel%mygrid_range(1)-1
      do i=ix,n1
         DO ish=-norder,norder
            derf(2,i,iy,1)=derf(2,i,iy,1)+Grad(ish,2)*wrap_box1d(diff_map%local_map1d(i,iy+ish,iz))
            ! i_temp=mod(iy+ish+global_n2,global_n2)
            ! derf(2,i,iy,1)=derf(2,i,iy,1)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         DO ish=-norder,norder
            derf(2,i,j,1)=derf(2,i,j,1)+Grad(ish,2)*wrap_box1d(diff_map%local_map1d(i,j+ish,iz))
            ! i_temp=mod(j+ish+global_n2,global_n2)
            ! derf(2,i,j,1)=derf(2,i,j,1)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         DO ish=-norder,norder
            derf(2,i,j,k)=derf(2,i,j,k)+Grad(ish,2)*wrap_box1d(diff_map%local_map1d(i,j+ish,iz))
            ! i_temp=mod(j+ish+global_n2,global_n2)
            ! derf(2,i,j,k)=derf(2,i,j,k)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         DO ish=-norder,norder
            derf(2,i,j,n3)=derf(2,i,j,n3)+Grad(ish,2)*wrap_box1d(diff_map%local_map1d(i,j+ish,iz))
            ! i_temp=mod(j+ish+global_n2,global_n2)
            ! derf(2,i,j,n3)=derf(2,i,j,n3)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      enddo
      do i=1,ix
         DO ish=-norder,norder
            derf(2,i,iy,n3)=derf(2,i,iy,n3)+Grad(ish,2)*wrap_box1d(diff_map%local_map1d(i,iy+ish,iz))
            ! i_temp=mod(iy+ish+global_n2,global_n2)
            ! derf(2,i,iy,n3)=derf(2,i,iy,n3)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo

      !$OMP END DO

      !df/dz
      !$OMP DO
# 708 "Finite_module.f90"
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0 !parallel%mygrid_range(1)-1
      do i=ix,n1
         DO ish=-norder,norder
            derf(3,i,iy,1)=derf(3,i,iy,1)+Grad(ish,3)*wrap_box1d(diff_map%local_map1d(i,iy,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,iy,1)=derf(3,i,iy,1)+Grad(ish,3)*wrap_box(i,iy,i_temp)
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         DO ish=-norder,norder
            derf(3,i,j,1)=derf(3,i,j,1)+Grad(ish,3)*wrap_box1d(diff_map%local_map1d(i,j,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,j,1)=derf(3,i,j,1)+Grad(ish,3)*wrap_box(i,j,i_temp)
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         DO ish=-norder,norder
            derf(3,i,j,k)=derf(3,i,j,k)+Grad(ish,3)*wrap_box1d(diff_map%local_map1d(i,j,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,j,k)=derf(3,i,j,k)+Grad(ish,3)*wrap_box(i,j,i_temp)
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         DO ish=-norder,norder
            derf(3,i,j,n3)=derf(3,i,j,n3)+Grad(ish,3)*wrap_box1d(diff_map%local_map1d(i,j,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,j,n3)=derf(3,i,j,n3)+Grad(ish,3)*wrap_box(i,j,i_temp)
         ENDDO
      enddo
      enddo
      do i=1,ix
         DO ish=-norder,norder
            derf(3,i,iy,n3)=derf(3,i,iy,n3)+Grad(ish,3)*wrap_box1d(diff_map%local_map1d(i,iy,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,iy,n3)=derf(3,i,iy,n3)+Grad(ish,3)*wrap_box(i,iy,i_temp)
         ENDDO
      enddo
      call end_time('nabla1',.true.)

      !$OMP END DO
      !$OMP END PARALLEL
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cmplx_nabla1
   !----------------------------parting line-----------------------------
   SUBROUTINE nabla(ifun,norder,derf)

      USE grid_module , ONLY :KPT,global_n1,global_n2,global_n3,&
           & global_n
      ! USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per
      USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per_ata

      IMPLICIT NONE
      !IN/OUT
      COMPLEX(DP),INTENT(IN)  :: ifun(:,:,:)
      INTEGER(I4B),INTENT(IN)  :: norder
      COMPLEX(DP),INTENT(OUT) :: derf(:,:,:,:)
      !LOCAL
      INTEGER(I4B) :: ix,iy,iz

      INTEGER(I4B) :: i,j,k,in,id
      INTEGER(I4B) :: y_down,d_z

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> comm setting init
      if(Lall_grid)then
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         allocate(wrap_box1d(global_n1*global_n2*&
          & global_n3))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         call MPI_ALLGATHERV(ifun(ix,iy,1),parallel%mygrid_range(3),&
              & MPI_COMPLEX16,wrap_box1d,parallel%recvcounts&
              & ,parallel%displs, MPI_COMPLEX16, parallel%commx&
              & ,mpinfo)
         ! !> set wrap_grid
         ! j=0
         ! do i=diff_map%boundary(1,3),diff_map%boundary(2,3),1
         !    id=i
         !    j=j+1
         !    do while(id<1)
         !       id=id+global_n3
         !    enddo
         !    do while(id>global_n3)
         !       id=id-global_n3
         !    enddo
         !    wrap_box1d((j-1)*global_n1*global_n2+1:j*global_n1*global_n2)=reshape(fun_global(:,:,id),(/global_n1*global_n2/))
         ! enddo
      else
         !> set wrap_grid
         allocate(wrap_box1d(global_n1*global_n2*&
          & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1) ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         allocate(fun_1d(parallel%mygrid_range(3)))!> oneD rho smaller then 3d

         ix=(parallel%mygrid_range(2)-1)/n1/n2
         iy=(parallel%mygrid_range(1)-1)/n1/n2
         d_z=ceiling((ix-iy)/(ix-iy+1.d0))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
         in=0 !parallel%mygrid_range(1)-1
         do i=ix,n1
            in=in+1
            fun_1d(in)=ifun(i,iy,iz)
         enddo
         do j=iy+1,n2*d_z
            do i=1,n1
               in=in+1
               fun_1d(in)=ifun(i,j,iz)
            enddo
         enddo
         do k=2,n3-1
            do j=1,n2
               do i=1,n1
                  in=in+1
                  fun_1d(in)=ifun(i,j,k)
               enddo
            enddo
         enddo
         ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
         y_down=(1-iy-1)*d_z+iy+1
         iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
         iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
         do j=y_down,iy-1
            do i=1,n1
               in=in+1
               fun_1d(in)=ifun(i,j,iz)
            enddo
         enddo
         do i=1,ix
            in=in+1
            fun_1d(in)=ifun(i,iy,iz)
         enddo
         call set_wrap_grid_per_ata(fun_1d,wrap_box1d,global_n&
              & ,global_n1, global_n2)
      endif

      CALL cmplx_nabla1(ifun,norder,derf)
      !transform to cartensian
      DO iz=1,n3
      DO iy=1,n2
      DO ix=1,n1
         derf(:,ix,iy,iz)=MATMUL(tBmat,derf(:,ix,iy,iz))
      ENDDO
      ENDDO
      ENDDO

      if(allocated(wrap_box))deallocate(wrap_box)
      if(allocated(wrap_box1d))deallocate(wrap_box1d)
      if(allocated(fun_global))deallocate(fun_global)
      if(allocated(fun_1d))deallocate(fun_1d)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nabla
   SUBROUTINE ksKEperiod_op(psi,Ik,KEpsi)
      !#################################################################!
      !ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
      !#################################################################!
      USE struct_module, ONLY : recip_lat
      use m_time_evaluate, only: memory_sum,memory_free,filename

      USE grid_module , ONLY :KPT,global_n1,global_n2,global_n3,&
           & global_n
      ! USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per
      USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per_ata



      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)  :: Ik
      COMPLEX(DP),INTENT(IN)  :: psi(:,:,:)
      COMPLEX(DP),INTENT(OUT) :: KEpsi(:,:,:)
      !
      REAL(DP) :: kmod2_2
      COMPLEX(DP) :: Ts1(n1,n2,n3)
      COMPLEX(DP) :: Ts2(n1,n2,n3)

      INTEGER(I4B) :: ix,iy,iz,id,i,j,k,in
      ! COMPLEX(DCP) :: temp_global(global_n)
      INTEGER(I4B) :: y_down,d_z


      call memory_sum("ksKEperiod_op",real(size(Ts1),DP)*DP+size(Ts2)*DP)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !kvec(:)=MATMUL(recip_lat,kp)
      !|\vec{k}|/2
      kmod2_2=SUM(KPT%vcar(:,Ik)*KPT%vcar(:,Ik))/2


      if(Lall_grid)then
         allocate(wrap_box1d(global_n1*global_n2*&
          & global_n3 ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         call MPI_ALLGATHERV(psi(ix,iy,1),parallel%mygrid_range(3),&
              & MPI_COMPLEX16,wrap_box1d,parallel%recvcounts&
              & ,parallel%displs, MPI_COMPLEX16, parallel%commx&
              & ,mpinfo)
      else
         !> set wrap_grid
         allocate(wrap_box1d(global_n1*global_n2*&
          & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1) ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         allocate(fun_1d(parallel%mygrid_range(3)))!> oneD rho smaller then 3d
         ix=(parallel%mygrid_range(2)-1)/n1/n2
         iy=(parallel%mygrid_range(1)-1)/n1/n2
         d_z=ceiling((ix-iy)/(ix-iy+1.d0))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
         in=0 !parallel%mygrid_range(1)-1
         do i=ix,n1
            in=in+1
            fun_1d(in)=psi(i,iy,iz)
         enddo
         do j=iy+1,n2*d_z
            do i=1,n1
               in=in+1
               fun_1d(in)=psi(i,j,iz)
            enddo
         enddo
         do k=2,n3-1
            do j=1,n2
               do i=1,n1
                  in=in+1
                  fun_1d(in)=psi(i,j,k)
               enddo
            enddo
         enddo
         ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
         y_down=(1-iy-1)*d_z+iy+1
         iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
         iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
         do j=y_down,iy-1
            do i=1,n1
               in=in+1
               fun_1d(in)=psi(i,j,iz)
            enddo
         enddo
         do i=1,ix
            in=in+1
            fun_1d(in)=psi(i,iy,iz)
         enddo
         call set_wrap_grid_per_ata(fun_1d,wrap_box1d,global_n&
              & ,global_n1, global_n2)
         ! print*,'communite over rank->',parallel%myid
      endif

      !Ts1=-0.5*\nabla^2 \psi
      CALL KE1(psi,Ts1)
         ! print*,'KE1 over rank->',parallel%myid
      ! open(1111+parallel%myid,file="Ts1"//filename(10:11))
      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! call MPI_ALLGATHERV(Ts1(ix,iy,1),parallel%mygrid_range(3)&
      !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
      !      &,parallel%displs,MPI_COMPLEX16&
      !      & ,parallel%commx,mpinfo)
      ! write(1111+parallel%myid,*)real(temp_global)
      ! close(1111+parallel%myid)
      !Ts2=-i\vec{k}\cdot\nabla\psi
      CALL KE2(psi,Ik,Ts2)
         ! print*,'KE2 over rank->',parallel%myid
      ! open(1111+parallel%myid,file="Ts2"//filename(10:11))
      ! call MPI_ALLGATHERV(Ts2(ix,iy,1),parallel%mygrid_range(3)&
      !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
      !      &,parallel%displs,MPI_COMPLEX16&
      !      & ,parallel%commx,mpinfo)
      ! write(1111+parallel%myid,*)real(temp_global)
      ! close(1111+parallel%myid)

      !>
      !Ts3=|\vec{k}|^2/
      !Ts=Ts1+Ts2+Ts3
      KEpsi=Ts1(:,:,:)+Ts2(:,:,:)+kmod2_2*psi(:,:,:)
      !  !KEpsi=Ts2+kmod2_2*psi(:,:,:)
      ! call MPI_BARRIER(parallel%commx,mpinfo)
      ! stop

      if(allocated(wrap_box))deallocate(wrap_box)
      if(allocated(wrap_box1d))deallocate(wrap_box1d)
      if(allocated(fun_global))deallocate(fun_global)
      if(allocated(fun_1d))deallocate(fun_1d)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! call memory_free("ksKEperiod_op",real(size(Ts1),DP)*DP+size(Ts2)*DP)
   ENDSUBROUTINE ksKEperiod_op
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE KE1(psi,Ts1)
      !#################################################################!
      !Ts1=-1/2*(\nabla)^2\psi
      !#################################################################!
      USE parameters , ONLY : finite_order
      IMPLICIT NONE
      COMPLEX(DP),INTENT(IN)  :: psi(:,:,:)
      COMPLEX(DP),INTENT(OUT) :: Ts1(:,:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      CALL cmplx_nabla2(psi,finite_order,Ts1)

      Ts1(:,:,:)=-0.5d0*Ts1(:,:,:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE KE1
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE KE2(psi,Ik,Ts2)
      !#################################!
      !  Ts2=-i\kvec\cdot\nabla         !
      !#################################!
      USE parameters , ONLY :finite_order
      USE grid_module , ONLY :KPT
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Ik
      COMPLEX(DP),INTENT(IN)  :: psi(:,:,:)
      COMPLEX(DP),INTENT(OUT) :: Ts2(:,:,:)
      !
      INTEGER(I4B) :: ix,iy,iz
      COMPLEX(DP) :: Rcoe(3)

      COMPLEX(DP) :: napsi(3,n1,n2,n3)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Rcoe(:)=-IMAG*MATMUL(KPT%vcar(:,Ik),tBmat(:,:))
      !
      CALL cmplx_nabla1(psi,finite_order,napsi)
      !
      DO iz=1,n3
      DO iy=1,n2
      DO ix=1,n1
         Ts2(ix,iy,iz)=SUM(napsi(:,ix,iy,iz)*Rcoe(:))
      ENDDO
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE KE2
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !####################################################################!
   !       For real wave function finite difference operator            !
   !*Author : Qiang Xu                                                  !
   !*Date   : 2017/08/28                                                !
   !####################################################################!
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE realnabla2(ifun,norder,ofun)
      !##########################################################
      !*For    : -1/2*\nabla^2(f)
      !*Author : Xuecheng Shao & Qiang Xu
      !*Date   : 2017/08/08
      !##########################################################
      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: norder
      REAL(DP),INTENT(IN)     :: ifun(:,:,:)
      REAL(DP),INTENT(OUT)    :: ofun(:,:,:)
      !
      REAL(DP) :: vec(-norder+1:norder+n1,-norder+1:norder+n2,-norder+1:norder+n3)
      !
      INTEGER(I4B) :: i,ish,ix,iy,iz
      REAL(DP)     :: kcoe
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      kcoe = SUM( Lapl(0,1:3) )
      DO iz = 1,n3
         vec(1:n1,1:n2,iz)           = ifun(1:n1,1:n2,iz)
         vec(-norder+1:0,1:n2,iz)    = ifun(n1-norder+1:n1,1:n2,iz)
         vec(n1+1:n1+norder,1:n2,iz) = ifun(1:norder,1:n2,iz)
      ENDDO

      DO iz = 1, n3
         vec(:,-norder+1:0,iz)=vec(:,n2-norder+1:n2,iz)
         vec(:,n2+1:n2+norder,iz)=vec(:,1:norder,iz)
      ENDDO

      DO iz = -norder+1, 0
         vec(:,:, iz)=vec(:,:, iz+n3)
      ENDDO

      DO iz = 1, norder
         vec(:,:, n3+iz)=vec(:,:, iz)
      ENDDO
      !
      DO iz=1,n3
         DO iy=1,n2
            DO ix=1,n1
               ofun(ix,iy,iz) = kcoe * vec(ix,iy,iz)
               DO ish = 1,norder
                  ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                     Lapl(ish,1)*(vec(ix-ish,iy,iz) + vec(ix+ish,iy,iz) ) + &
                     Lapl(ish,2)*(vec(ix,iy-ish,iz) + vec(ix,iy+ish,iz) ) + &
                     Lapl(ish,3)*(vec(ix,iy,iz-ish) + vec(ix,iy,iz+ish) )
               ENDDO
            ENDDO
         enddo
      ENDDO

      IF ( lap_add(1) /= 0 ) then
         kcoe = Lapl(0,4)
         DO iz=1,n3
            DO iy=1,n2
               DO ix=1,n1
                  ofun(ix,iy,iz)=ofun(ix,iy,iz)+kcoe * vec(ix,iy,iz)
                  DO ish = 1,norder
                     ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                        Lapl(ish,4)*( vec(ix-ish,iy-lap_add(1)*ish,iz) + vec(ix+ish,iy+lap_add(1)*ish,iz) )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF ( lap_add(2) /= 0 ) then
         kcoe = Lapl(0,5)
         DO iz=1,n3
            DO iy=1,n2
               DO ix=1,n1
               ofun(ix,iy,iz)=ofun(ix,iy,iz)+kcoe * vec(ix,iy,iz)
                  DO ish = 1,norder
                     ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                        Lapl(ish,5)*( vec(ix,iy-ish,iz-lap_add(2)*ish) + vec(ix,iy+ish,iz+lap_add(2)*ish))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF ( lap_add(3) /= 0 ) then
         kcoe = Lapl(0,6)
         DO iz=1,n3
            DO iy=1,n2
               DO ix=1,n1
               ofun(ix,iy,iz)=ofun(ix,iy,iz)+kcoe * vec(ix,iy,iz)
                  DO ish = 1,norder
                     ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                        Lapl(ish,6)*( vec(ix-ish,iy,iz-lap_add(3)*ish) + vec(ix+ish,iy,iz+lap_add(3)*ish))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE realnabla2
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE realKE(psi,KEpsi)
      !#################################################################!
      !ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
      !#################################################################!
      USE parameters , ONLY : finite_order
      REAL(DP),INTENT(IN)  :: psi(:,:,:)
      REAL(DP),INTENT(OUT) :: KEpsi(:,:,:)
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !Ts1=-0.5*\nabla^2 \psi
      CALL realnabla2(psi,finite_order,KEpsi)
      KEpsi(:,:,:)=-0.5d0*KEpsi(:,:,:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE realKE
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !####################################################################!
   !       For real wave function finite difference operator            !
   !*Author : Qiang Xu                                                  !
   !*Date   : 2017/08/28                                                !
   !####################################################################!
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE ISO_realnabla2(ifun,norder,ofun)
      !##########################################################
      !*For    : -1/2*\nabla^2(f)
      !*Author : Xuecheng Shao & Qiang Xu
      !*Date   : 2017/08/08
      !##########################################################
      USE grid_module , ONLY :rho_calc

      USE smpi_math_module, ONLY :parallel,start_time,end_time,mpi_real8,&
           & mpi_sum,mpinfo,set_wrap_grid_iso,diff_map

      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: norder





      REAL(DP),INTENT(IN)     :: ifun(:)
      REAL(DP),INTENT(OUT)    :: ofun(:)
      ! REAL(DP) :: vec_local(-norder+1:norder+n1,-norder+1:norder+n2,-norder+1:norder+n3)
      REAL(DP) :: vec(diff_map%boundary(1,1):diff_map%boundary(2,1),&
           & diff_map%boundary(1,2):diff_map%boundary(2,2),&
           & diff_map%boundary(1,3):diff_map%boundary(2,3))

      !
      !
      INTEGER(I4B) :: i,ish,ix,iy,iz,ir
      REAL(DP)     :: kcoe

      REAL(DP),DIMENSION(parallel%mygrid_range(3)) :: out_sphere

      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      kcoe = SUM( Lapl(0,1:3) )
      !DO iz = 1,n3
      !   vec(1:n1,1:n2,iz)           = ifun(1:n1,1:n2,iz)
      !   vec(-norder+1:0,1:n2,iz)    = 0
      !   vec(n1+1:n1+norder,1:n2,iz) = 0
      !ENDDO

      !DO iz = 1, n3
      !   vec(:,-norder+1:0,iz)=0
      !   vec(:,n2+1:n2+norder,iz)=0
      !ENDDO

      !DO iz = -norder+1, 0
      !   vec(:,:, iz)=0
      !ENDDO

      !DO iz = 1, norder
      !   vec(:,:, n3+iz)=0
      !ENDDO
      !==================================================
# 1258 "Finite_module.f90"
      CALL set_wrap_grid_iso(ifun,vec)

! #ifndef 1
!       DO I=1,rho_calc%OneDLength
!          out_sphere(I)=ifun(rho_calc%x(I),rho_calc%y(I),rho_calc%z(I))
!       ENDDO
!       open(111,file="ifun_s")
!       write(111,*)vec
!       close(111)
! #else
!       if(parallel%myid==1)then
!          print *,'parallel%mygrid_range',parallel%mygrid_range(3)
!       DO I=1,parallel%mygrid_range(3)
!          out_sphere(I)=vec(rho_calc%x(I),rho_calc%y(I),rho_calc%z(I))
!       ENDDO
!       ! CALL parallel_s2g(vh,out_grid)
!          open(111,file="ifun_p")
!          write(111,*)ifun
!          close(111)
!          open(111,file="vec_p")
!          write(111,*)out_sphere
!          close(111)

!       endif
! #endif
!       CALL MPI_BARRIER(parallel%comm,mpinfo)
!       stop
      !==================================================
      !
      !==================================================
      !##cal only sphere



      DO ir=1,parallel%mygrid_range(3)

        ix=rho_calc%x(ir)
        iy=rho_calc%y(ir)
        iz=rho_calc%z(ir)
# 1305 "Finite_module.f90"
        ofun(ir) = kcoe * vec(ix,iy,iz)
        DO ish = 1,norder
          ofun(ir)=ofun(ir)+&
          Lapl(ish,1)*(vec(ix-ish,iy,iz) + vec(ix+ish,iy,iz) ) + &
          Lapl(ish,2)*(vec(ix,iy-ish,iz) + vec(ix,iy+ish,iz) ) + &
          Lapl(ish,3)*(vec(ix,iy,iz-ish) + vec(ix,iy,iz+ish) )

        ENDDO
      ENDDO
      !==================================================
      !DO iz=1,n3
      !   DO iy=1,n2
      !      DO ix=1,n1
      !         ofun(ix,iy,iz) = kcoe * vec(ix,iy,iz)
      !         DO ish = 1,norder
      !            ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
      !               Lapl(ish,1)*(vec(ix-ish,iy,iz) + vec(ix+ish,iy,iz) ) + &
      !               Lapl(ish,2)*(vec(ix,iy-ish,iz) + vec(ix,iy+ish,iz) ) + &
      !               Lapl(ish,3)*(vec(ix,iy,iz-ish) + vec(ix,iy,iz+ish) )
      !         ENDDO
      !      ENDDO
      !   enddo
      !ENDDO

      IF ( lap_add(1) /= 0 ) then
         kcoe = Lapl(0,4)



         DO ir=1,parallel%mygrid_range(3)

           ix=rho_calc%x(ir)
           iy=rho_calc%y(ir)
           iz=rho_calc%z(ir)






           ofun(ir)=ofun(ir)+kcoe * vec(ix,iy,iz)
           DO ish = 1,norder
             ofun(ir)=ofun(ir)+&
                Lapl(ish,4)*( vec(ix-ish,iy-lap_add(1)*ish,iz) + vec(ix+ish,iy+lap_add(1)*ish,iz) )

           ENDDO
         ENDDO
      ENDIF

      IF ( lap_add(2) /= 0 ) then
         kcoe = Lapl(0,5)



         DO ir=1,parallel%mygrid_range(3)

           ix=rho_calc%x(ir)
           iy=rho_calc%y(ir)
           iz=rho_calc%z(ir)






           ofun(ir)=ofun(ir)+kcoe * vec(ix,iy,iz)
           DO ish = 1,norder
             ofun(ir)=ofun(ir)+&
                Lapl(ish,5)*( vec(ix,iy-ish,iz-lap_add(2)*ish) + vec(ix,iy+ish,iz+lap_add(2)*ish))

           ENDDO
         ENDDO
      ENDIF

      IF ( lap_add(3) /= 0 ) then
         kcoe = Lapl(0,6)



         DO ir=1,parallel%mygrid_range(3)

           ix=rho_calc%x(ir)
           iy=rho_calc%y(ir)
           iz=rho_calc%z(ir)






           ofun(ir)=ofun(ir)+kcoe * vec(ix,iy,iz)
           DO ish = 1,norder
             ofun(ir)=ofun(ir)+&
                Lapl(ish,6)*( vec(ix-ish,iy,iz-lap_add(3)*ish) + vec(ix+ish,iy,iz+lap_add(3)*ish))

           ENDDO
         ENDDO
      ENDIF
! #ifdef 1
!       if(parallel%isroot)then
!          ! print *,'parallel%mygrid_range',parallel%mygrid_range(3)
!       ! DO I=1,parallel%mygrid_range(3)
!          ! out_sphere(I)=ofun(rho_calc%x(I),rho_calc%y(I),rho_calc%z(I))
!       ! ENDDO
!       ! CALL parallel_s2g(vh,out_grid)
!          open(111,file="ofun_p")
!          write(111,*)ofun
!          close(111)
!          ! open(111,file="vec_p")
!          ! write(111,*)out_sphere
!          ! close(111)

!       endif
! #endif
!       CALL MPI_BARRIER(parallel%comm,mpinfo)
!       stop
      !
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_realnabla2
   !----------------------------PARTING LINE-----------------------------
   ! SUBROUTINE real_nabla2_wb(vec,norder,ofun)!{{{
   !    !##########################################################
   !    !* CREATED_TIME  : 2019-05-06
   !    !* AUTHOR        : Qiang Xu
   !    !* DESCRIPTION   :              
   !    !     ------
   !    !* REFERENCES    :              
   !    !     ------
   !    !* LOG           :              
   !    !     2015-05-08 :              
   !    !* LAST MODIFIED : 2015-05-10 07:47:12 PM
   !    !##########################################################
   !    use constants
   !    !USE parameters,        only : norder=>finite_order
   !    USE grid_module , ONLY : gap
   !    !
   !    implicit none
   !    !IN/OUT
   !    REAL(DP),INTENT(IN) :: vec(-norder+1:,-norder+1:,-norder+1:)
   !    INTEGER(I4B),INTENT(IN) :: norder
   !    REAL(DP),INTENT(OUT) :: ofun(:,:,:)
   !    !real(dp)             :: coeke(-norder:norder,6)
   !    !
   !    INTEGER(I4B)              :: i,ish,ix,iy,iz
   !    INTEGER(I4B) :: nn1,nn2,nn3
   !    !real(dp)             :: gaps(6)
   !    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
   !    !zero boundary
   !    !vec=0._DP
   !    !vec(1:n1,1:n2,1:n3)=func(1:n1,1:n2,1:n3)
   !    nn1=SIZE(ofun,1)
   !    nn2=SIZE(ofun,2)
   !    nn3=SIZE(ofun,3)

   !    do iz=1,nn3
   !       do iy=1,nn2
   !          do ix=1,nn1
   !             ofun(ix,iy,iz) = 0.d0
   !             do ish = -norder,norder
   !                ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
   !                 &  Lapl(ish,1)*(vec(ix+ish,iy,iz))+ &
   !                 &  Lapl(ish,2)*(vec(ix,iy+ish,iz))+ &
   !                 &  Lapl(ish,3)*(vec(ix,iy,iz+ish))+ &
   !                 &  Lapl(ish,4)*(vec(ix+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish))+&
   !                 &  Lapl(ish,5)*(vec(ix+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish))+&
   !                 &  Lapl(ish,6)*(vec(ix+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish))
   !             ENDDO
   !          enddo
   !       enddo
   !    ENDDO
   !    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
   ! ENDSUBROUTINE real_nabla2_wb!}}}
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE ISO_realKE(psi,KEpsi)
      !#################################################################!
      !ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
      !#################################################################!
      USE parameters , ONLY : finite_order




      REAL(DP),INTENT(IN)  :: psi(:)
      REAL(DP),INTENT(OUT) :: KEpsi(:)

      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !Ts1=-0.5*\nabla^2 \psi
      CALL ISO_realnabla2(psi,finite_order,KEpsi)



      KEpsi(:)=-0.5d0*KEpsi(:)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_realKE
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   SUBROUTINE ksKEperiod_op_band(psi,Ik,KEpsi)
      !#################################################################!
      !ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
      !#################################################################!
      USE struct_module, ONLY : recip_lat
      use m_time_evaluate, only: memory_sum,memory_free,filename
      USE grid_module , ONLY :KPT,n1=>global_n1,n2=>global_n2,n3=>global_n3
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)  :: Ik
      COMPLEX(DP),INTENT(IN)  :: psi(:,:,:)
      COMPLEX(DP),INTENT(OUT) :: KEpsi(:,:,:)
      !
      REAL(DP) :: kmod2_2
      COMPLEX(DP) :: Ts1(n1,n2,n3)
      COMPLEX(DP) :: Ts2(n1,n2,n3)
      call memory_sum("ksKEperiod_op",real(size(Ts1),DP)*DP+size(Ts2)*DP)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !kvec(:)=MATMUL(recip_lat,kp)
      !|\vec{k}|/2
      kmod2_2=SUM(KPT%vcar(:,Ik)*KPT%vcar(:,Ik))/2

      !Ts1=-0.5*\nabla^2 \psi
      CALL KE1_np(psi,Ts1)
         ! print*,'KE1 over rank->',parallel%myid
      ! open(1111+parallel%myid,file="Ts1"//filename(10:11))
      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! call MPI_ALLGATHERV(Ts1(ix,iy,1),parallel%mygrid_range(3)&
      !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
      !      &,parallel%displs,MPI_COMPLEX16&
      !      & ,parallel%commx,mpinfo)
      ! write(1111+parallel%myid,*)real(temp_global)
      ! close(1111+parallel%myid)
      !Ts2=-i\vec{k}\cdot\nabla\psi
      CALL KE2_np(psi,Ik,Ts2)
         ! print*,'KE2 over rank->',parallel%myid
      ! open(1111+parallel%myid,file="Ts2"//filename(10:11))
      ! call MPI_ALLGATHERV(Ts2(ix,iy,1),parallel%mygrid_range(3)&
      !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
      !      &,parallel%displs,MPI_COMPLEX16&
      !      & ,parallel%commx,mpinfo)
      ! write(1111+parallel%myid,*)real(temp_global)
      ! close(1111+parallel%myid)

      !>
      !Ts3=|\vec{k}|^2/
      !Ts=Ts1+Ts2+Ts3
      KEpsi=Ts1(:,:,:)+Ts2(:,:,:)+kmod2_2*psi(:,:,:)
      !  !KEpsi=Ts2+kmod2_2*psi(:,:,:)
      ! call MPI_BARRIER(parallel%commx,mpinfo)
      ! stop
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_free("ksKEperiod_op",real(size(Ts1),DP)*DP+size(Ts2)*DP)
    ENDSUBROUTINE ksKEperiod_op_band
   SUBROUTINE nabla2_np(ifun,norder,ofun)
      !##########################################################
      !*For    : -1/2*\nabla^2(f)
      !*Author : Xuecheng Shao & Qiang Xu
      !*Date   : 2017/08/08
      !##########################################################
     USE grid_module, only:n1=>global_n1,n2=>global_n2,n3=>global_n3
      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: norder
      COMPLEX(DP),INTENT(IN)     :: ifun(:,:,:)
      COMPLEX(DP),INTENT(OUT)    :: ofun(:,:,:)
      !
      COMPLEX(DP) :: vec(-norder+1:norder+n1,-norder+1:norder+n2,-norder+1:norder+n3)
      !
      INTEGER(I4B) :: i,ish,ix,iy,iz
      REAL(DP)     :: kcoe
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      kcoe = SUM( Lapl(0,1:3) )
      !$OMP PARALLEL PRIVATE(iz,iy,ix,ish)
      !$OMP DO
      DO iz = 1,n3
         vec(1:n1,1:n2,iz)           = ifun(1:n1,1:n2,iz)
         vec(-norder+1:0,1:n2,iz)    = ifun(n1-norder+1:n1,1:n2,iz)
         vec(n1+1:n1+norder,1:n2,iz) = ifun(1:norder,1:n2,iz)
      ENDDO
      !OMP END DO

      !$OMP DO
      DO iz = 1, n3
         vec(:,-norder+1:0,iz)=vec(:,n2-norder+1:n2,iz)
         vec(:,n2+1:n2+norder,iz)=vec(:,1:norder,iz)
      ENDDO
      !OMP END DO

      !$OMP DO
      DO iz = -norder+1, 0
         vec(:,:, iz)=vec(:,:, iz+n3)
      ENDDO
      !OMP END DO

      !$OMP DO
      DO iz = 1, norder
         vec(:,:, n3+iz)=vec(:,:, iz)
      ENDDO
      !OMP END DO
      !
      !$OMP DO
      DO iz=1,n3
         DO iy=1,n2
            DO ix=1,n1
               ofun(ix,iy,iz) = kcoe * vec(ix,iy,iz)
               DO ish = 1,norder
                  ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                     Lapl(ish,1)*(vec(ix-ish,iy,iz) + vec(ix+ish,iy,iz) ) + &
                     Lapl(ish,2)*(vec(ix,iy-ish,iz) + vec(ix,iy+ish,iz) ) + &
                     Lapl(ish,3)*(vec(ix,iy,iz-ish) + vec(ix,iy,iz+ish) )
               ENDDO
            ENDDO
         enddo
      ENDDO
     ! print*,'kcoe',kcoe
     ! open(1111+parallel%myid,file="Ts1"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(ofun(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! stop
      !OMP END DO

      IF ( lap_add(1) /= 0 ) then
         kcoe = Lapl(0,4)
      !$OMP DO
         DO iz=1,n3
            DO iy=1,n2
               DO ix=1,n1
                  ofun(ix,iy,iz)=ofun(ix,iy,iz)+kcoe * vec(ix,iy,iz)
                  DO ish = 1,norder
                     ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                        Lapl(ish,4)*( vec(ix-ish,iy-lap_add(1)*ish,iz) + vec(ix+ish,iy+lap_add(1)*ish,iz) )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      !OMP END DO
      ENDIF

      IF ( lap_add(2) /= 0 ) then
         kcoe = Lapl(0,5)
      !$OMP DO
         DO iz=1,n3
            DO iy=1,n2
               DO ix=1,n1
               ofun(ix,iy,iz)=ofun(ix,iy,iz)+kcoe * vec(ix,iy,iz)
                  DO ish = 1,norder
                     ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                        Lapl(ish,5)*( vec(ix,iy-ish,iz-lap_add(2)*ish) + vec(ix,iy+ish,iz+lap_add(2)*ish))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      !OMP END DO
      ENDIF

      IF ( lap_add(3) /= 0 ) then
         kcoe = Lapl(0,6)
      !$OMP DO
         DO iz=1,n3
            DO iy=1,n2
               DO ix=1,n1
               ofun(ix,iy,iz)=ofun(ix,iy,iz)+kcoe * vec(ix,iy,iz)
                  DO ish = 1,norder
                     ofun(ix,iy,iz)=ofun(ix,iy,iz)+&
                        Lapl(ish,6)*( vec(ix-ish,iy,iz-lap_add(3)*ish) + vec(ix+ish,iy,iz+lap_add(3)*ish))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      !OMP END DO
      ENDIF
      !$OMP END PARALLEL
      !
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE nabla2_np
   !--------------------------PARTING LINE-------------------------------
   SUBROUTINE nabla1_np(ifun,norder,derf)
      !
     USE grid_module, only:n1=>global_n1,n2=>global_n2,n3=>global_n3
      IMPLICIT NONE
      !
      COMPLEX(DP),INTENT(IN)  :: ifun(:,:,:)
      INTEGER(I4B),INTENT(IN)  :: norder
      COMPLEX(DP),INTENT(OUT) :: derf(:,:,:,:)
      !
      COMPLEX(DP) :: vec(-norder+1:norder+n1,-norder+1:norder+n2,-norder+1:norder+n3)
      INTEGER(I4B) :: i,ish,ix,iy,iz
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !$OMP PARALLEL PRIVATE(iz,iy,ix,ish)
      !$OMP DO
      DO iz = 1, n3
         !! xm ym zm
         vec(1:n1,1:n2,iz)           = ifun(1:n1,1:n2,iz)
         !! xl ym zm
         vec(-norder+1:0,1:n2,iz)    = ifun(n1-norder+1:n1,1:n2,iz)
         !! xr ym zm
         vec(n1+1:n1+norder,1:n2,iz) = ifun(1:norder,1:n2,iz)
      ENDDO
      !$OMP END DO
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
      !!$OMP END DO

      !df/dx
      !$OMP DO
      DO iz=1,n3
         derf(1,:,:,iz) = (0.d0,0.d0)
      DO iy=1,n2
      DO ix=1,n1
         DO ish=-norder,norder
            derf(1,ix,iy,iz)=derf(1,ix,iy,iz)+Grad(ish,1)*vec(ix+ish,iy,iz)
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END DO

      !df/dy
      !$OMP DO
      DO iz=1,n3
         derf(2,:,:,iz) = (0.d0,0.d0)
      DO iy=1,n2
      DO ish = -norder,norder
      DO ix=1,n1
         derf(2,ix,iy,iz)=derf(2,ix,iy,iz)+ Grad(ish,2)*(vec(ix,iy+ish,iz))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END DO

      !df/dz
      !$OMP DO
      DO iz=1,n3
         derf(3,:,:,iz) = (0.d0,0.d0)
      DO ish = -norder,norder
      DO iy=1,n2
      DO ix=1,n1
         derf(3,ix,iy,iz)=derf(3,ix,iy,iz)+ Grad(ish,3)*(vec(ix,iy,iz+ish))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nabla1_np
   !----------------------------parting line-----------------------------
   SUBROUTINE KE1_np(psi,Ts1)
      !#################################################################!
      !Ts1=-1/2*(\nabla)^2\psi
      !#################################################################!
      USE parameters , ONLY : finite_order
      IMPLICIT NONE
      COMPLEX(DP),INTENT(IN)  :: psi(:,:,:)
      COMPLEX(DP),INTENT(OUT) :: Ts1(:,:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      CALL nabla2_np(psi,finite_order,Ts1)

      Ts1(:,:,:)=-0.5d0*Ts1(:,:,:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE KE1_np
   !----------------------------PARTING LINE-----------------------------
   SUBROUTINE KE2_np(psi,Ik,Ts2)
      !#################################!
      !  Ts2=-i\kvec\cdot\nabla         !
      !#################################!
      USE parameters , ONLY :finite_order
      USE grid_module , ONLY :KPT,n1=>global_n1,n2=>global_n2,n3=>global_n3
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Ik
      COMPLEX(DP),INTENT(IN)  :: psi(:,:,:)
      COMPLEX(DP),INTENT(OUT) :: Ts2(:,:,:)
      !
      INTEGER(I4B) :: ix,iy,iz
      COMPLEX(DP) :: Rcoe(3)

      COMPLEX(DP) :: napsi(3,n1,n2,n3)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Rcoe(:)=-IMAG*MATMUL(KPT%vcar(:,Ik),tBmat(:,:))
      !
      CALL nabla1_np(psi,finite_order,napsi)
      !
      DO iz=1,n3
      DO iy=1,n2
      DO ix=1,n1
         Ts2(ix,iy,iz)=SUM(napsi(:,ix,iy,iz)*Rcoe(:))
      ENDDO
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE KE2_np

   SUBROUTINE ksKEperiod_op_gamma(psi,Ik,KEpsi)
      !#################################################################!
      !ksKE_op\psi=Ts1+Ts2+k^2/2\psi                                    !
      !#################################################################!
      USE struct_module, ONLY : recip_lat
      use m_time_evaluate, only: memory_sum,memory_free,filename

      USE grid_module , ONLY :KPT,global_n1,global_n2,global_n3,&
           & global_n
      ! USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per
      USE smpi_math_module, ONLY: parallel,diff_map,mpi_real8,mpinfo,Lall_grid,set_wrap_grid_per_ata_real



      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)  :: Ik
      REAL(DP),INTENT(IN)  :: psi(:,:,:)
      REAL(DP),INTENT(OUT) :: KEpsi(:,:,:)
      !
      REAL(DP) :: kmod2_2
      REAL(DP) :: Ts1(n1,n2,n3)
      REAL(DP) :: Ts2(n1,n2,n3)

      INTEGER(I4B) :: ix,iy,iz,id,i,j,k,in
      ! COMPLEX(DCP) :: temp_global(global_n)
      INTEGER(I4B) :: y_down,d_z


      call memory_sum("ksKEperiod_op",real(size(Ts1),DP)*DP+size(Ts2)*DP)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !kvec(:)=MATMUL(recip_lat,kp)
      !|\vec{k}|/2
      kmod2_2=SUM(KPT%vcar(:,Ik)*KPT%vcar(:,Ik))/2


      if(Lall_grid)then
         allocate(wrap_box1d_real(global_n1*global_n2*&
          & global_n3 ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         call MPI_ALLGATHERV(psi(ix,iy,1),parallel%mygrid_range(3),&
              & MPI_real8,wrap_box1d_real,parallel%recvcounts&
              & ,parallel%displs, MPI_real8, parallel%commx&
              & ,mpinfo)
      else
         !> set wrap_grid
         allocate(wrap_box1d_real(global_n1*global_n2*&
          & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1) ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         allocate(fun_1d_real(parallel%mygrid_range(3)))!> oneD rho smaller then 3d
         ix=(parallel%mygrid_range(2)-1)/n1/n2
         iy=(parallel%mygrid_range(1)-1)/n1/n2
         d_z=ceiling((ix-iy)/(ix-iy+1.d0))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
         in=0 !parallel%mygrid_range(1)-1
         do i=ix,n1
            in=in+1
            fun_1d_real(in)=psi(i,iy,iz)
         enddo
         do j=iy+1,n2*d_z
            do i=1,n1
               in=in+1
               fun_1d_real(in)=psi(i,j,iz)
            enddo
         enddo
         do k=2,n3-1
            do j=1,n2
               do i=1,n1
                  in=in+1
                  fun_1d_real(in)=psi(i,j,k)
               enddo
            enddo
         enddo
         ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
         y_down=(1-iy-1)*d_z+iy+1
         iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
         iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
         do j=y_down,iy-1
            do i=1,n1
               in=in+1
               fun_1d_real(in)=psi(i,j,iz)
            enddo
         enddo
         do i=1,ix
            in=in+1
            fun_1d_real(in)=psi(i,iy,iz)
         enddo
         call set_wrap_grid_per_ata_real(fun_1d_real,wrap_box1d_real,global_n&
              & ,global_n1, global_n2)
         ! print*,'communite over rank->',parallel%myid
      endif

      !Ts1=-0.5*\nabla^2 \psi
      CALL KE_gamma(psi,KEpsi)
      ! stop

      if(allocated(wrap_box_real))deallocate(wrap_box_real)
      if(allocated(wrap_box1d_real))deallocate(wrap_box1d_real)
      if(allocated(fun_global_real))deallocate(fun_global_real)
      if(allocated(fun_1d_real))deallocate(fun_1d_real)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! call memory_free("ksKEperiod_op",real(size(Ts1),DP)*DP+size(Ts2)*DP)
    ENDSUBROUTINE ksKEperiod_op_gamma
   SUBROUTINE KE_gamma(psi,Ts1)
      !#################################################################!
      !Ts1=-1/2*(\nabla)^2\psi
      !#################################################################!
      USE parameters , ONLY : finite_order
      IMPLICIT NONE
      REAL(DP),INTENT(IN)  :: psi(:,:,:)
      REAL(DP),INTENT(OUT) :: Ts1(:,:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      ! CALL real_nabla2_1d(psi(1:,1:,1:),Ts1)
      CALL real_nabla2(psi(1:,1:,1:),finite_order,Ts1)

      Ts1(:,:,:)=-0.5d0*Ts1(:,:,:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE KE_GAMMA
   SUBROUTINE real_nabla2(ifun,norder,ofun)
      !##########################################################
      !*For    : -1/2*\nabla^2(f)
      !*Author : Xuecheng Shao & Qiang Xu
      !*Date   : 2017/08/08
      !##########################################################

    USE smpi_math_module, ONLY: parallel,start_time,end_time,diff_map,write_time

      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: norder
      REAL(DP),INTENT(IN)     :: ifun(:,:,:)
      REAL(DP),INTENT(OUT)    :: ofun(:,:,:)
      !



      !
      INTEGER(I4B) :: i,ish,ix,iy,iz
      REAL(DP)     :: kcoe

      INTEGER(I4B) :: j,k
      ! COMPLEX(DCP) :: temp_global(global_n)
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>

      call start_time('nabla2',.true.)

      kcoe = SUM( Lapl(0,1:3) )
      !$OMP PARALLEL PRIVATE(iz,iy,ix,ish)
      !$OMP DO
# 2004 "Finite_module.f90"
      !OMP END DO
      !
      !$OMP DO
# 2025 "Finite_module.f90"
      ix=(parallel%mygrid_range(2)-1)/n1/n2
      iy=(parallel%mygrid_range(1)-1)/n1/n2
      d_z=ceiling((ix-iy)/(ix-iy+1.d0))
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      ! iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
      do i=ix,n1
         ofun(i,iy,1) = 0.d0
         DO ish = -norder,norder
            ofun(i,iy,1)=ofun(i,iy,1)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,iy,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,iy+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,iy,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         ofun(i,j,1) = 0
         DO ish = -norder,norder
            ofun(i,j,1)=ofun(i,j,1)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         ofun(i,j,k) = 0
         DO ish = -norder,norder
            ofun(i,j,k)=ofun(i,j,k)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      ! iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         ofun(i,j,n3) = 0
         DO ish = -norder,norder
            ofun(i,j,n3)=ofun(i,j,n3)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      do i=1,ix
         ofun(i,iy,n3) = 0
         DO ish = -norder,norder
            ofun(i,iy,n3)=ofun(i,iy,n3)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,iy,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,iy+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,iy,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo

     ! print*,'kcoe',kcoe
     ! open(1111+parallel%myid,file="Ts1"//filename(10:11))
     ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     ! call MPI_ALLGATHERV(ofun(ix,iy,1),parallel%mygrid_range(3)&
     !      &,MPI_COMPLEX16,temp_global,parallel%recvcounts&
     !      &,parallel%displs,MPI_COMPLEX16&
     !      & ,parallel%commx,mpinfo)
     ! write(1111+parallel%myid,*)real(temp_global)
     ! close(1111+parallel%myid)
     ! call MPI_BARRIER(parallel%commx,mpinfo)
     ! stop
      !OMP END DO

      call end_time('nabla2',.true.)
      ! call write_time('nabla2',.true.)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE real_nabla2
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
      REAL(DP),INTENT(IN)  :: ifun(*) 
      REAL(DP),INTENT(OUT) :: ofun(n) 
      !

      INTEGER(I4B) :: ip,ish,ix,iy,iz &
         &, id1,id2,id3,id4,id5,id6   &
         &, shiftn
      INTEGER(I4B) :: i,j,k,in
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      !zero boundary
      ! shiftn=parallel%mygrid_range(1)-1
      ! !
      ! DO ip=1,n

      !    iz=diff_map%local_map(3,shiftn+ip)
      !    iy=diff_map%local_map(2,shiftn+ip)
      !    ix=diff_map%local_map(1,shiftn+ip)

      !    ofun(ip) = 0._DP

      !    DO ish = -norder,norder

      !       id1 = diff_map%local_map1d(ix+ish,iy,iz)
      !       id2 = diff_map%local_map1d(ix,iy+ish,iz)
      !       id3 = diff_map%local_map1d(ix,iy,iz+ish)
      !       id4 = diff_map%local_map1d(ix+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)
      !       id5 = diff_map%local_map1d(ix+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)
      !       id6 = diff_map%local_map1d(ix+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)

      !       ofun(ip)=ofun(ip)+&
      !        &  Lapl(ish,1)*wrap_box1d_real(id1)+ &
      !        &  Lapl(ish,2)*wrap_box1d_real(id2)+ &
      !        &  Lapl(ish,3)*wrap_box1d_real(id3)+ &
      !        &  Lapl(ish,4)*wrap_box1d_real(id4)+ &
      !        &  Lapl(ish,5)*wrap_box1d_real(id5)+ &
      !        &  Lapl(ish,6)*wrap_box1d_real(id6)
      !    ENDDO 

      ! ENDDO
      ix=(parallel%mygrid_range(2)-1)/n1/n2
      iy=(parallel%mygrid_range(1)-1)/n1/n2
      d_z=ceiling((ix-iy)/(ix-iy+1.d0))
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0
      ! iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
      do i=ix,n1
         in=in+1
         ofun(in) = 0.d0
         DO ish = -norder,norder
            ofun(in)=ofun(in)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,iy,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,iy+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,iy,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         in=in+1
         ofun(in) = 0
         DO ish = -norder,norder
            ofun(in)=ofun(in)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         in=in+1
         ofun(in) = 0
         DO ish = -norder,norder
            ofun(in)=ofun(in)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      ! iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         in=in+1
         ofun(in) = 0
         DO ish = -norder,norder
            ofun(in)=ofun(in)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,j+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,j+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,j+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      enddo
      do i=1,ix
         in=in+1
         ofun(in) = 0
         DO ish = -norder,norder
            ofun(in)=ofun(in)+&
                 Lapl(ish,1)*(wrap_box1d_real(diff_map%local_map1d(i+ish,iy,iz)) ) + &
                 Lapl(ish,2)*(wrap_box1d_real(diff_map%local_map1d(i,iy+ish,iz)) ) + &
                 Lapl(ish,3)*(wrap_box1d_real(diff_map%local_map1d(i,iy,iz+ish)) ) + &
                 Lapl(ish,4)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(1,1)*ish,iy+cell_mu(1,2)*ish,iz+cell_mu(1,3)*ish)))+&
                 Lapl(ish,5)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(2,1)*ish,iy+cell_mu(2,2)*ish,iz+cell_mu(2,3)*ish)))+&
                 Lapl(ish,6)*(wrap_box1d_real(diff_map%local_map1d(i+cell_mu(3,1)*ish,iy+cell_mu(3,2)*ish,iz+cell_mu(3,3)*ish)))
         ENDDO
      enddo
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_nabla2_1d
   SUBROUTINE nabla_gamma(ifun,norder,derf)

      USE grid_module , ONLY :KPT,global_n1,global_n2,global_n3,&
           & global_n
      ! USE smpi_math_module, ONLY: parallel,diff_map,mpi_complex16,mpinfo,Lall_grid,set_wrap_grid_per
      USE smpi_math_module, ONLY: parallel,diff_map,mpi_real8,mpinfo,Lall_grid,set_wrap_grid_per_ata_real

      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN)  :: ifun(:,:,:)
      INTEGER(I4B),INTENT(IN)  :: norder
      REAL(DP),INTENT(OUT) :: derf(:,:,:,:)
      !LOCAL
      INTEGER(I4B) :: ix,iy,iz

      INTEGER(I4B) :: i,j,k,in,id
      INTEGER(I4B) :: y_down,d_z

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> comm setting init
      if(Lall_grid)then
         allocate(wrap_box1d_real(global_n1*global_n2*&
          & global_n3 ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         call MPI_ALLGATHERV(ifun(ix,iy,1),parallel%mygrid_range(3),&
              & MPI_real8,wrap_box1d_real,parallel%recvcounts&
              & ,parallel%displs, MPI_real8, parallel%commx&
              & ,mpinfo)
      else
         !> set wrap_grid
         allocate(wrap_box1d_real(global_n1*global_n2*&
          & (diff_map%boundary(2,3)-diff_map%boundary(1,3)+1) ))
         ! allocate(wrap_box(diff_map%boundary(1,1):diff_map%boundary(2,1)&
         !      & ,diff_map%boundary(1,2):diff_map%boundary(2,2)&
         !      & ,diff_map%boundary(1,3):diff_map%boundary(2,3) ))
         allocate(fun_1d_real(parallel%mygrid_range(3)))!> oneD rho smaller then 3d
         ix=(parallel%mygrid_range(2)-1)/n1/n2
         iy=(parallel%mygrid_range(1)-1)/n1/n2
         d_z=ceiling((ix-iy)/(ix-iy+1.d0))
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
         in=0 !parallel%mygrid_range(1)-1
         do i=ix,n1
            in=in+1
            fun_1d_real(in)=ifun(i,iy,iz)
         enddo
         do j=iy+1,n2*d_z
            do i=1,n1
               in=in+1
               fun_1d_real(in)=ifun(i,j,iz)
            enddo
         enddo
         do k=2,n3-1
            do j=1,n2
               do i=1,n1
                  in=in+1
                  fun_1d_real(in)=ifun(i,j,k)
               enddo
            enddo
         enddo
         ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
         y_down=(1-iy-1)*d_z+iy+1
         iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
         iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
         do j=y_down,iy-1
            do i=1,n1
               in=in+1
               fun_1d_real(in)=ifun(i,j,iz)
            enddo
         enddo
         do i=1,ix
            in=in+1
            fun_1d_real(in)=ifun(i,iy,iz)
         enddo
         call set_wrap_grid_per_ata_real(fun_1d_real,wrap_box1d_real,global_n&
              & ,global_n1, global_n2)
         ! print*,'communite over rank->',parallel%myid
      endif

      CALL real_nabla1(ifun,norder,derf)
      !transform to cartensian
      DO iz=1,n3
      DO iy=1,n2
      DO ix=1,n1
         derf(:,ix,iy,iz)=MATMUL(tBmat,derf(:,ix,iy,iz))
      ENDDO
      ENDDO
      ENDDO

      if(allocated(wrap_box_real))deallocate(wrap_box_real)
      if(allocated(wrap_box1d_real))deallocate(wrap_box1d_real)
      if(allocated(fun_global_real))deallocate(fun_global_real)
      if(allocated(fun_1d_real))deallocate(fun_1d_real)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE nabla_gamma
   SUBROUTINE real_nabla1(ifun,norder,derf)

      USE smpi_math_module, ONLY: parallel,start_time,end_time,diff_map
      USE Grid_module, ONLY:global_n1,global_n2

      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN)  :: ifun(:,:,:)
      INTEGER(I4B),INTENT(IN)  :: norder
      REAL(DP),INTENT(OUT) :: derf(:,:,:,:)
      !



      INTEGER(I4B) :: i,ish,ix,iy,iz

      INTEGER(I4B) :: j,k,in,i_temp
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !$OMP PARALLEL PRIVATE(iz,iy,ix,ish)
      !$OMP DO
# 2422 "Finite_module.f90"

      !df/dx
      !$OMP DO
# 2437 "Finite_module.f90"
      call start_time('nabla1',.true.)
      derf=(0.d0,0.d0)
      ix=(parallel%mygrid_range(2)-1)/n1/n2
      iy=(parallel%mygrid_range(1)-1)/n1/n2
      d_z=ceiling((ix-iy)/(ix-iy+1.d0))
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0 !parallel%mygrid_range(1)-1
      do i=ix,n1
         DO ish=-norder,norder
            derf(1,i,iy,1)=derf(1,i,iy,1)+Grad(ish,1)*wrap_box1d_real(diff_map%local_map1d(i+ish,iy,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,iy,1)=derf(1,i,iy,1)+Grad(ish,1)*wrap_box(i_temp,iy,iz)
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         DO ish=-norder,norder
            derf(1,i,j,1)=derf(1,i,j,1)+Grad(ish,1)*wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,j,1)=derf(1,i,j,1)+Grad(ish,1)*wrap_box(i_temp,j,iz)
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         DO ish=-norder,norder
            derf(1,i,j,k)=derf(1,i,j,k)+Grad(ish,1)*wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,j,k)=derf(1,i,j,k)+Grad(ish,1)*wrap_box(i_temp,j,iz)
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         DO ish=-norder,norder
            derf(1,i,j,n3)=derf(1,i,j,n3)+Grad(ish,1)*wrap_box1d_real(diff_map%local_map1d(i+ish,j,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,j,n3)=derf(1,i,j,n3)+Grad(ish,1)*wrap_box(i_temp,j,iz)
         ENDDO
      enddo
      enddo
      do i=1,ix
         DO ish=-norder,norder
            derf(1,i,iy,n3)=derf(1,i,iy,n3)+Grad(ish,1)*wrap_box1d_real(diff_map%local_map1d(i+ish,iy,iz))
            ! i_temp=mod(i+ish+global_n1,global_n1)
            ! derf(1,i,iy,n3)=derf(1,i,iy,n3)+Grad(ish,1)*wrap_box(i_temp,iy,iz)
         ENDDO
      enddo

      !$OMP END DO

      !df/dy
      !$OMP DO
# 2511 "Finite_module.f90"
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0 !parallel%mygrid_range(1)-1
      do i=ix,n1
         DO ish=-norder,norder
            derf(2,i,iy,1)=derf(2,i,iy,1)+Grad(ish,2)*wrap_box1d_real(diff_map%local_map1d(i,iy+ish,iz))
            ! i_temp=mod(iy+ish+global_n2,global_n2)
            ! derf(2,i,iy,1)=derf(2,i,iy,1)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         DO ish=-norder,norder
            derf(2,i,j,1)=derf(2,i,j,1)+Grad(ish,2)*wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz))
            ! i_temp=mod(j+ish+global_n2,global_n2)
            ! derf(2,i,j,1)=derf(2,i,j,1)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         DO ish=-norder,norder
            derf(2,i,j,k)=derf(2,i,j,k)+Grad(ish,2)*wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz))
            ! i_temp=mod(j+ish+global_n2,global_n2)
            ! derf(2,i,j,k)=derf(2,i,j,k)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         DO ish=-norder,norder
            derf(2,i,j,n3)=derf(2,i,j,n3)+Grad(ish,2)*wrap_box1d_real(diff_map%local_map1d(i,j+ish,iz))
            ! i_temp=mod(j+ish+global_n2,global_n2)
            ! derf(2,i,j,n3)=derf(2,i,j,n3)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo
      enddo
      do i=1,ix
         DO ish=-norder,norder
            derf(2,i,iy,n3)=derf(2,i,iy,n3)+Grad(ish,2)*wrap_box1d_real(diff_map%local_map1d(i,iy+ish,iz))
            ! i_temp=mod(iy+ish+global_n2,global_n2)
            ! derf(2,i,iy,n3)=derf(2,i,iy,n3)+Grad(ish,2)*wrap_box(i,i_temp,iz)
         ENDDO
      enddo

      !$OMP END DO

      !df/dz
      !$OMP DO
# 2580 "Finite_module.f90"
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(1)-1)/n1/n2+1
      in=0 !parallel%mygrid_range(1)-1
      do i=ix,n1
         DO ish=-norder,norder
            derf(3,i,iy,1)=derf(3,i,iy,1)+Grad(ish,3)*wrap_box1d_real(diff_map%local_map1d(i,iy,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,iy,1)=derf(3,i,iy,1)+Grad(ish,3)*wrap_box(i,iy,i_temp)
         ENDDO
      enddo
      do j=iy+1,n2*d_z
      do i=1,n1
         DO ish=-norder,norder
            derf(3,i,j,1)=derf(3,i,j,1)+Grad(ish,3)*wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,j,1)=derf(3,i,j,1)+Grad(ish,3)*wrap_box(i,j,i_temp)
         ENDDO
      enddo
      enddo
      do k=2,n3-1
         iz=iz+1
      do j=1,n2
      do i=1,n1
         DO ish=-norder,norder
            derf(3,i,j,k)=derf(3,i,j,k)+Grad(ish,3)*wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,j,k)=derf(3,i,j,k)+Grad(ish,3)*wrap_box(i,j,i_temp)
         ENDDO
      enddo
      enddo
      enddo
      ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
      y_down=(1-iy-1)*d_z+iy+1
      iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
      iz=(parallel%mygrid_range(2)-1)/n1/n2+1
      do j=y_down,iy-1
      do i=1,n1
         DO ish=-norder,norder
            derf(3,i,j,n3)=derf(3,i,j,n3)+Grad(ish,3)*wrap_box1d_real(diff_map%local_map1d(i,j,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,j,n3)=derf(3,i,j,n3)+Grad(ish,3)*wrap_box(i,j,i_temp)
         ENDDO
      enddo
      enddo
      do i=1,ix
         DO ish=-norder,norder
            derf(3,i,iy,n3)=derf(3,i,iy,n3)+Grad(ish,3)*wrap_box1d_real(diff_map%local_map1d(i,iy,iz+ish))
            ! i_temp=mod(iz+ish+global_n3,global_n3)
            ! derf(3,i,iy,n3)=derf(3,i,iy,n3)+Grad(ish,3)*wrap_box(i,iy,i_temp)
         ENDDO
      enddo
      call end_time('nabla1',.true.)

      !$OMP END DO
      !$OMP END PARALLEL
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_nabla1
ENDMODULE finite_module
