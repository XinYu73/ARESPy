# 1 "XC_functional.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "XC_functional.f90"
MODULE libxc_module
   USE constants
   USE parameters , ONLY : NSPIN,iXCDF
   USE grid_module , ONLY : nsn,dvol,n,n1,n2,n3
   !libxc
   USE xc_f90_types_m
   USE xc_f90_lib_m
   USE libxc_funcs_m

   USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

   USE m_time_evaluate, ONLY: memory_sum,memory_free
   IMPLICIT NONE
   !
CONTAINS
   !---------------------------energy--------------------------
   SUBROUTINE LDAlib_energy(rhoS,Exc)



      use grid_module, only:nsn,dvol,n,n1,n2,n3,global_n1,global_n2,global_n3

      use m_time_evaluate, only:filename
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) ::Exc
      !LOCAL
      TYPE(xc_f90_pointer_t) :: x_func,c_func
      TYPE(xc_f90_pointer_t) :: x_info,c_info
      REAL(DP) :: ext,ect,exct(n1,n2,n3)
      REAL(DP) :: rhot(NSPIN)
      INTEGER(I4B) :: Ix,Iy,Iz,Is

      REAL(DP),external :: ddot
      REAL(DP) :: Exc_local
      INTEGER(I4B) :: i,j,k,in
      REAL(DP) :: rho_global(global_n1,global_n2,global_n3)
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('totalenergy_xc',real(size(exct),DP)*DP)
      SELECT CASE(iXCDF)
      CASE(XC_FAMILY_LDA)
         !initial exchange and correlation
         CALL xc_f90_func_init(x_func,x_info,XC_LDA_X,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_LDA_C_PZ,NSPIN)

         !> TEST TEST TEST
! #ifdef 1
!          ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!          iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!          call MPI_ALLGATHERV(rhoS(ix,iy,1,1),parallel%mygrid_range(3)&
!               &,MPI_REAL8,rho_global,parallel%recvcounts&
!               &,parallel%displs,MPI_REAL8&
!               & ,parallel%commx,mpinfo)
!          open(1111+parallel%myid,file="rho"//filename(10:11))
!          ! do i=1,size(psi,2),1
!          !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
!          ! enddo
!          write(1111+parallel%myid,*)rho_global
!          close(1111+parallel%myid)
!          ! stop
!          print*,'sum rho_global',sum(rho_global)
! #endif
         !calculate the exc
# 79 "XC_functional.f90"
        ! print*,'z2-z1',(parallel%mygrid_range(2)-1)/n1/n2+1-(parallel%mygrid_range(1)-1)/n1/n2-1+1,'local_n3',local_n3
        ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
        iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
        iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
        in=0 !parallel%mygrid_range(1)-1
        do i=ix,n1
           in=in+1
           rhot(:)=rhoS(i,iy,iz,:)
           CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)
           CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)
           exct(I,Iy,Iz)=ext+ect
        enddo
        d_z=ceiling((n3-iz)/(n3-iz+1.d0))
        do j=iy+1,n2*d_z
           do i=1,n1
              in=in+1
              rhot(:)=rhoS(i,j,iz,:)
              CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)
              CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)
              exct(I,J,Iz)=ext+ect
           enddo
        enddo
        do k=2,n3-1
           do j=1,n2
              do i=1,n1
                 in=in+1
                 rhot(:)=rhoS(i,j,k,:)
                 CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)
                 CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)
                 exct(I,J,K)=ext+ect
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
              rhot(:)=rhoS(i,j,iz,:)
              CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)
              CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)
              exct(I,J,Iz)=ext+ect
           enddo
        enddo
        do i=1,ix
           in=in+1
           rhot(:)=rhoS(i,iy,iz,:)
           CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)
           CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)
           exct(I,Iy,Iz)=ext+ect
        enddo

      CASE(XC_FAMILY_GGA)
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      CASE default
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      ENDSELECT

      !total xc energy
      Exc=0.d0





      Exc_local=0.d0
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=1
      DO Is=1,NSPIN
         Exc_local=Exc_local+ddot(parallel%mygrid_range(3), rhoS(ix,iy,iz,Is),1,exct(ix,iy,iz),1 )
      ENDDO
      call mpi_allreduce(Exc_local, Exc, 1, mpi_real8, mpi_sum, parallel%commx, mpinfo)

      Exc=Exc*dvol
      !exit xc call
      CALL xc_f90_func_end(x_func)
      CALL xc_f90_func_end(c_func)
      call memory_free('totalenergy_xc',real(size(exct),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE  LDAlib_energy
   !-----------------------LDA potential------------------------
   SUBROUTINE LDAlib_potential(rhoS,vxcS)
     ! use grid_module, only: global_n1,global_n2,global_n3
     ! use m_time_evaluate, only: filename
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: vxcS(:,:,:,:)
      !LOCAL
      TYPE(xc_f90_pointer_t) :: x_func,c_func
      TYPE(xc_f90_pointer_t) :: x_info,c_info
      REAL(DP) :: vxt(NSPIN),vct(NSPIN),rhot(NSPIN) !temp array
      INTEGER(I4B) :: Ix,Iy,Iz

      INTEGER(I4B) :: i,j,k,in
      ! REAL(DP)     :: vxc_global(global_n1,global_n2,global_n3)
      INTEGER(I4B) :: d_z,y_down

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !select a functional
      SELECT CASE(iXCDF)
      CASE(XC_FAMILY_LDA)
         !initial exchange and correlation
         CALL xc_f90_func_init(x_func,x_info,XC_LDA_X,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_LDA_C_PZ,NSPIN)

         !calculate the exc
# 204 "XC_functional.f90"
         vxcs=0.d0
         ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
         iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
         iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
         in=0 !parallel%mygrid_range(1)-1
         do i=ix,n1
            in=in+1
            rhot(:)=rhoS(i,iy,iz,:)
            CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
            CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
            vxcS(I,Iy,Iz,:)=vxt+vct
         enddo
         d_z=ceiling((n3-iz)/(n3-iz+1.d0))
         do j=iy+1,n2*d_z
            do i=1,n1
               in=in+1
               rhot(:)=rhoS(i,j,iz,:)
               CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
               CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
               vxcS(I,J,Iz,:)=vxt+vct
            enddo
         enddo
         do k=2,n3-1
            do j=1,n2
               do i=1,n1
                  in=in+1
                  rhot(:)=rhoS(i,j,k,:)
                  CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
                  CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
                  vxcS(I,J,K,:)=vxt+vct
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
               rhot(:)=rhoS(i,j,iz,:)
               CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
               CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
               vxcS(I,J,Iz,:)=vxt+vct
            enddo
         enddo
         do i=1,ix
            in=in+1
            rhot(:)=rhoS(i,iy,iz,:)
            CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
            CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
            vxcS(I,Iy,Iz,:)=vxt+vct
         enddo

      CASE(XC_FAMILY_GGA)
         !!initial exchange and correlation
         !CALL xc_f90_func_init(x_func,x_info,XC_GGA_X_PBE,NSPIN)
         !CALL xc_f90_func_init(c_func,c_info,XC_GGA_C_PBE,NSPIN)
         !!calculate the exc
         !CALL xc_f90_lda_vxc(x_func,nsn,rhoS(1),sigma(1),vxt(1))
         !CALL xc_f90_lda_vxc(c_func,nsn,rhoS(1),sigma(1),vct(1))
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      CASE default
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      ENDSELECT
      !exit xc call
! #ifdef 1
!       ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!       iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!       call MPI_ALLGATHERV(vxcs(ix,iy,1,1),parallel%mygrid_range(3)&
!            &,MPI_REAL8,vxc_global,parallel%recvcounts&
!            &,parallel%displs,MPI_REAL8&
!            & ,parallel%commx,mpinfo)
!       open(1111+parallel%myid,file="rho"//filename(10:11))
!       ! do i=1,size(psi,2),1
!       !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
!       ! enddo
!       write(1111+parallel%myid,*)vxc_global
!       close(1111+parallel%myid)
!       stop
! #endif
      CALL xc_f90_func_end(x_func)
      CALL xc_f90_func_end(c_func)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE LDAlib_potential
   !---------------------------energy--------------------------
   SUBROUTINE LDAlib_energy_iso(rhoS,Exc)
      IMPLICIT NONE
      !IN/OUT



      REAL(DP),INTENT(IN)  :: rhoS(:,:)

      REAL(DP),INTENT(OUT) ::Exc
      !LOCAL
      TYPE(xc_f90_pointer_t) :: x_func,c_func
      TYPE(xc_f90_pointer_t) :: x_info,c_info



      REAL(DP) :: ext,ect,exct(parallel%mygrid_range(3))

      REAL(DP) :: rhot(NSPIN)
      INTEGER(I4B) :: Ix,Iy,Iz,Is

      REAL(DP),external :: ddot
      REAL(DP) :: Exc_local

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('iso_totalenergy_xc',real(size(exct),DP)*DP)
      SELECT CASE(iXCDF)
      CASE(XC_FAMILY_LDA)
         !initial exchange and correlation
         CALL xc_f90_func_init(x_func,x_info,XC_LDA_X,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_LDA_C_PZ,NSPIN)

         !calculate the exc
# 336 "XC_functional.f90"
         DO Ix=1,parallel%mygrid_range(3)
            CALL xc_f90_lda_exc(x_func,1,rhoS(Ix,1),ext)
            CALL xc_f90_lda_exc(c_func,1,rhoS(Ix,1),ect)
            exct(Ix)=ext+ect
         ENDDO

      CASE(XC_FAMILY_GGA)
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      CASE default
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      ENDSELECT

      !total xc energy
      Exc=0.d0





      Exc_local=0.d0
      DO Is=1,NSPIN
         Exc_local=Exc_local+ddot(parallel%mygrid_range(3), rhoS(:,Is),1,exct,1 )
      ENDDO
      call mpi_allreduce(Exc_local, Exc, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

      Exc=Exc*dvol
      !exit xc call
      CALL xc_f90_func_end(x_func)
      CALL xc_f90_func_end(c_func)
      call memory_free('iso_totalenergy_xc',real(size(exct),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE  LDAlib_energy_iso
   !-----------------------LDA potential------------------------
   SUBROUTINE LDAlib_potential_iso(rhoS,vxcS)
      IMPLICIT NONE
      !IN/OUT




      REAL(DP),INTENT(IN)  :: rhoS(:,:)
      REAL(DP),INTENT(OUT) :: vxcS(:,:)

      !LOCAL
      TYPE(xc_f90_pointer_t) :: x_func,c_func
      TYPE(xc_f90_pointer_t) :: x_info,c_info
      REAL(DP) :: vxt(NSPIN),vct(NSPIN),rhot(NSPIN) !temp array
      INTEGER(I4B) :: Ix,Iy,Iz
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !select a functional
      SELECT CASE(iXCDF)
      CASE(XC_FAMILY_LDA)
         !initial exchange and correlation
         CALL xc_f90_func_init(x_func,x_info,XC_LDA_X,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_LDA_C_PZ,NSPIN)

         !calculate the exc
# 408 "XC_functional.f90"
         DO Ix=1,parallel%mygrid_range(3)
            rhot(:)=rhoS(Ix,:)
            CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
            CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
            !store
            vxcS(Ix,:)=vxt(:)+vct(:)
         ENDDO

      CASE(XC_FAMILY_GGA)
         !!initial exchange and correlation
         !CALL xc_f90_func_init(x_func,x_info,XC_GGA_X_PBE,NSPIN)
         !CALL xc_f90_func_init(c_func,c_info,XC_GGA_C_PBE,NSPIN)
         !!calculate the exc
         !CALL xc_f90_lda_vxc(x_func,nsn,rhoS(1),sigma(1),vxt(1))
         !CALL xc_f90_lda_vxc(c_func,nsn,rhoS(1),sigma(1),vct(1))
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      CASE default
         WRITE(6,*) "XC Functional: now haven't consider this case"
         STOP
      ENDSELECT
      !exit xc call
      CALL xc_f90_func_end(x_func)
      CALL xc_f90_func_end(c_func)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE LDAlib_potential_iso
   !------------------------------------------------------------
ENDMODULE libxc_module
