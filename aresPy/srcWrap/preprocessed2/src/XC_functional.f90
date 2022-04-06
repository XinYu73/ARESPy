!#################################!
!Qiang XU Sep.1st 2021            !
!#################################!
!---------------------------------------------------------------
MODULE xc_module
   USE constants
   USE parameters , ONLY : nspin,iXCDF=>ixc,Lcore_val,Snlcc
   USE grid_module , ONLY : sph=>grid,dvol
!libxc
   USE xc_f90_types_m
   USE xc_f90_lib_m
   USE libxc_funcs_m

   USE smpi_math_module

   IMPLICIT NONE
CONTAINS
!------------------------XC potential------------------------
   SUBROUTINE XC_functional(nps,rhoS,vxc,exc)
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps !num. of grids
      REAL(DP),INTENT(IN) :: rhoS(nps,nspin)
      REAL(DP),INTENT(OUT) :: vxc(nps,nspin)
      REAL(DP),OPTIONAL,INTENT(OUT) :: exc
!LOCAL
      LOGICAL :: lexc
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lexc=PRESENT(exc)
!spin-unpolarization
      IF(iXCDF>0)THEN !LDA family
         IF(lexc)THEN
            CALL LibXC_lda_set(nps,rhoS,vxc,exc)
         ELSE
            CALL LibXC_lda_set(nps,rhoS,vxc)
         ENDIF
      ELSE !GGA family
         IF(lexc)THEN
            CALL Libxc_gga_set(nps,rhoS,vxc,exc)
         ELSE
            CALL Libxc_gga_set(nps,rhoS,vxc)
         ENDIF
      ENDIF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE XC_functional
!------------------------------------------------------------
!##########################LDA-level########################!
!------------------------ Libxc_PZ81 ------------------------
   SUBROUTINE Libxc_lda_set(nps,rhoS,vxcS,exc)
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN)  :: rhoS(nps,nspin)
      REAL(DP),INTENT(OUT) :: vxcS(nps,nspin)
      REAL(DP),OPTIONAL,INTENT(OUT) :: exc
!LOCAL
      TYPE(xc_f90_pointer_t) :: x_func,c_func
      TYPE(xc_f90_pointer_t) :: x_info,c_info
      REAL(DP) :: vxt(NSPIN),vct(NSPIN),rhot(NSPIN) !temp array
      REAL(DP) :: ext,ect

      REAL(DP) :: exc_loc

      INTEGER(I4B) :: Ip
!
      LOGICAL :: lexc
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lexc=PRESENT(exc)
!initial exchange and correlation
      CALL xc_f90_func_init(x_func,x_info,XC_LDA_X,NSPIN)
      CALL xc_f90_func_init(c_func,c_info,XC_LDA_C_PZ,NSPIN)
!set zero
      IF(lexc)THEN 
         exc=0._DP

         exc_loc=0._DP

      ENDIF
!start
      DO Ip=1,nps
!linear/non-linear core-valance correlation
         IF(lcore_val)THEN
         
            IF(nspin==1)THEN 
               rhot(1)=rhoS(Ip,1)+sph%rhoc(Ip) !/nspin
            ELSEIF(nspin==2)THEN
               rhot(1)=rhoS(Ip,1)+sph%rhoc(Ip)*Snlcc
               rhot(2)=rhoS(Ip,2)+sph%rhoc(Ip)*(1._DP-Snlcc)
            ELSE
               STOP 'Error: Check nspin'
            ENDIF

         ELSE
            rhot(:)=rhoS(Ip,:)
         ENDIF
!Potential
         CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
         CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
!store
         vxcS(Ip,:)=vxt(:)+vct(:)
!energy
         IF(lexc)THEN
            CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)
            CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)

            exc_loc=exc_loc+(ext+ect)*SUM(rhot(:))

         ENDIF
      ENDDO
!total xc energy
      IF(lexc)THEN 

         CALL MPI_ALLREDUCE(exc_loc,exc,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

         exc=exc*dvol
      ENDIF
!exit xc call
      CALL xc_f90_func_end(x_func)
      CALL xc_f90_func_end(c_func)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Libxc_lda_set
!------------------------------------------------------------
!Slater-x
!------------------------------------------------------------
   SUBROUTINE Libxc_lda_x(nps,rhoS,ex,vxS)
!Slater-x
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN)  :: rhoS(nps,nspin)
      REAL(DP),INTENT(OUT) :: ex
      REAL(DP),OPTIONAL,INTENT(OUT) :: vxS(nps,nspin)
!LOCAL
      TYPE(xc_f90_pointer_t) :: x_func
      TYPE(xc_f90_pointer_t) :: x_info
      REAL(DP) :: vxt(nspin),rhot(nspin) !temp array
      REAL(DP) :: ext

      REAL(DP) :: ex_loc

      INTEGER(I4B) :: Ip
!
      LOGICAL :: lvx
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lvx=PRESENT(vxS)
!initial exchange and correlation
      CALL xc_f90_func_init(x_func,x_info,XC_LDA_X,nspin)
!set zero
      ex=0._DP

      ex_loc=0._DP

!start
      DO Ip=1,nps
!set rho
         rhot(:)=rhoS(Ip,:)
!Potential
         IF(lvx)THEN
            CALL xc_f90_lda_vxc(x_func,1,rhot(1),vxt(1))
!store
            vxS(Ip,:)=vxt(:)
         ENDIF
!energy
         CALL xc_f90_lda_exc(x_func,1,rhot(1),ext)

         ex_loc=ex_loc+(ext)*SUM(rhot(:))

      ENDDO
!total xc energy

      CALL MPI_ALLREDUCE(ex_loc,ex,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

      ex=ex*dvol
!exit xc call
      CALL xc_f90_func_end(x_func)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Libxc_lda_x
!------------------------------------------------------------
!VWN1RPA-c
!------------------------------------------------------------
   SUBROUTINE Libxc_vwn1rpa_c(nps,rhoS,ec,vcS)
!VWN1RPA-c
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN)  :: rhoS(nps,nspin)
      REAL(DP),INTENT(OUT) :: ec
      REAL(DP),OPTIONAL,INTENT(OUT) :: vcS(nps,nspin)
!LOCAL
      TYPE(xc_f90_pointer_t) :: c_func
      TYPE(xc_f90_pointer_t) :: c_info
      REAL(DP) :: vct(nspin),rhot(nspin) !temp array
      REAL(DP) :: ect

      REAL(DP) :: ec_loc

      INTEGER(I4B) :: Ip
!
      LOGICAL :: lvc
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lvc=PRESENT(vcS)
!initial exchange and correlation
      CALL xc_f90_func_init(c_func,c_info,XC_LDA_C_VWN_RPA,nspin)
!set zero
      ec=0._DP

      ec_loc=0._DP

!start
      DO Ip=1,nps
!set rhot
         rhot(:)=rhoS(Ip,:)
!Potential
         IF(lvc)THEN
            CALL xc_f90_lda_vxc(c_func,1,rhot(1),vct(1))
!store
            vcS(Ip,:)=vct(:)
         ENDIF
!energy
         CALL xc_f90_lda_exc(c_func,1,rhot(1),ect)

         ec_loc=ec_loc+(ect)*SUM(rhot(:))

      ENDDO
!total xc energy

      CALL MPI_ALLREDUCE(ec_loc,ec,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

      ec=ec*dvol
!exit xc call
      CALL xc_f90_func_end(c_func)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Libxc_vwn1rpa_c
!------------------------------------------------------------
!##########################GGA-level########################!
!------------------------------------------------------------
   SUBROUTINE Libxc_gga_set(nps,rhoS,vxcS,exc)
      USE finite_module, ONLY : real_nabla1_1d
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN)  :: rhoS(nps,nspin)
      REAL(DP),INTENT(OUT) :: vxcS(nps,nspin)
      REAL(DP),OPTIONAL,INTENT(OUT) :: exc
!LOCAL
      TYPE(xc_f90_pointer_t) :: x_func,c_func
      TYPE(xc_f90_pointer_t) :: x_info,c_info
      REAL(DP) :: vxt(NSPIN),vct(NSPIN),rhot(NSPIN) !temp array
      REAL(DP) :: ext,ect

      REAL(DP) :: exc_loc

      INTEGER(I4B) :: Ip
!
      LOGICAL :: lexc
      REAL(DP),ALLOCATABLE :: sigma(:,:) &
                     &, grad(:,:,:)  
      REAL(DP) :: rhoSc(nps,nspin)
      REAL(DP),DIMENSION(2*nspin-1) :: sigmat,vxsigmat,vcsigmat
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lexc=PRESENT(exc)
!initial exchange and correlation
      SELECT CASE(iXCDF)
      CASE(-1) !PW91
         CALL xc_f90_func_init(x_func,x_info,XC_GGA_X_PW91,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_GGA_C_PW91,NSPIN)
      CASE(-2) !PBE
         CALL xc_f90_func_init(x_func,x_info,XC_GGA_X_PBE,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_GGA_C_PBE,NSPIN)
      CASE(-3) !BLYP
         CALL xc_f90_func_init(x_func,x_info,XC_GGA_X_B88,NSPIN)
         CALL xc_f90_func_init(c_func,c_info,XC_GGA_C_LYP,NSPIN)
      CASE default
         WRITE(*,*) "XC Functional: now haven't consider this case"
         STOP
      ENDSELECT
!linear/non-linear core-valance correlation
      IF(lcore_val)THEN

         IF(nspin==1)THEN 
            rhoSc(:,1)=rhoS(:,1)+sph%rhoc(:) !/nspin
         ELSEIF(nspin==2)THEN
            rhoSc(:,1)=rhoS(:,1)+sph%rhoc(:)*Snlcc
            rhoSc(:,2)=rhoS(:,2)+sph%rhoc(:)*(1._DP-Snlcc)
         ELSE
            STOP 'Error: Check nspin'
         ENDIF
         
      ELSE
         rhoSc(:,:)=rhoS(:,:)
      ENDIF
!sigma
      IF(nspin==1)THEN !spin-unpolarized
         ALLOCATE(sigma(nps,1))
         ALLOCATE(grad(nps,3,1))
!The gradient of density
         CALL real_nabla1_1d(rhoSc(:,1),grad(:,:,1),sigma(:,1))
!up*up
         sigma(:,1)=sigma(:,1)**2
      ELSE !spin-polarized
         ALLOCATE(sigma(nps,3))
         ALLOCATE(grad(nps,3,2))
!The gradient of density
         CALL real_nabla1_1d(rhoSc(:,1),grad(:,:,1),sigma(:,1))
         CALL real_nabla1_1d(rhoSc(:,2),grad(:,:,2),sigma(:,3))
!up*up
         sigma(:,1)=sigma(:,1)**2
!up*down
         sigma(:,2) = grad(:,1,1)*grad(:,1,2) + &
                   &  grad(:,2,1)*grad(:,2,2) + &
                   &  grad(:,3,1)*grad(:,3,2)
!down*down
         sigma(:,3)=sigma(:,3)**2
      ENDIF
!set zero
      IF(lexc)THEN
         exc=0._DP

         exc_loc=0._DP

      ENDIF
!start
      DO Ip=1,nps
!set rho
         rhot(:)=rhoSc(Ip,:)
!set sigma
         sigmat(:)=sigma(Ip,:)
!Potential
         CALL xc_f90_gga_vxc(x_func,1,rhot(1),sigmat(1),vxt(1),vxsigmat(1))
         CALL xc_f90_gga_vxc(c_func,1,rhot(1),sigmat(1),vct(1),vcsigmat(1))
!store
         vxcS(Ip,:)=vxt(:)+vct(:)

         IF(lexc)THEN
!energy
            CALL xc_f90_gga_exc(x_func,1,rhot(1),sigmat(1),ext)
            CALL xc_f90_gga_exc(c_func,1,rhot(1),sigmat(1),ect)

            exc_loc=exc_loc+(ext+ect)*SUM(rhot(:))

         ENDIF

      ENDDO

      IF(lexc)THEN
!total xc energy

         CALL MPI_ALLREDUCE(exc_loc,exc,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

         exc=exc*dvol
      ENDIF
!exit xc call
      CALL xc_f90_func_end(x_func)
      CALL xc_f90_func_end(c_func)
      DEALLOCATE(sigma,grad)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Libxc_gga_set
!------------------------------------------------------------
!Becke-gga-x (1988)
!------------------------------------------------------------
   SUBROUTINE Libxc_b88_x(nps,rhoS,sigma,ex,vxS)
!Slater-x
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN)  :: rhoS(nps,nspin)
      REAL(DP),INTENT(IN)  :: sigma(:,:)
      REAL(DP),INTENT(OUT) :: ex
      REAL(DP),OPTIONAL,INTENT(OUT) :: vxS(nps,nspin)
!LOCAL
      TYPE(xc_f90_pointer_t) :: x_func
      TYPE(xc_f90_pointer_t) :: x_info
      REAL(DP) :: vxt(nspin),rhot(nspin) !temp array
      REAL(DP) :: ext

      REAL(DP) :: ex_loc

      INTEGER(I4B) :: Ip
      REAL(DP),DIMENSION(SIZE(sigma,2)) :: sigmat,vsigmat
!
      LOGICAL :: lvx
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lvx=PRESENT(vxS)
!initial exchange and correlation
      CALL xc_f90_func_init(x_func,x_info,XC_GGA_X_B88,nspin)
!set zero
      ex=0._DP

      ex_loc=0._DP

!start
      DO Ip=1,nps
!set rho
         rhot(:)=rhoS(Ip,:)
!set sigma
         sigmat(:)=sigma(Ip,:)
!Potential
         IF(lvx)THEN
            CALL xc_f90_gga_vxc(x_func,1,rhot(1),sigmat(1),vxt(1),vsigmat(1))
!store
            vxS(Ip,:)=vxt(:)
         ENDIF
!energy
         CALL xc_f90_gga_exc(x_func,1,rhot(1),sigmat(1),ext)

         ex_loc=ex_loc+(ext)*SUM(rhot(:))

      ENDDO
!total xc energy

      CALL MPI_ALLREDUCE(ex_loc,ex,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

      ex=ex*dvol
!exit xc call
      CALL xc_f90_func_end(x_func)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Libxc_b88_x
!------------------------------------------------------------
!Lee-Yang-Parr-gga-c (1988)
!------------------------------------------------------------
   SUBROUTINE Libxc_lyp_c(nps,rhoS,sigma,ec,vcS)
!LYP-c
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(IN)  :: rhoS(nps,nspin)
      REAL(DP),INTENT(IN)  :: sigma(:,:)
      REAL(DP),INTENT(OUT) :: ec
      REAL(DP),OPTIONAL,INTENT(OUT) :: vcS(nps,nspin)
!LOCAL
      TYPE(xc_f90_pointer_t) :: c_func
      TYPE(xc_f90_pointer_t) :: c_info
      REAL(DP) :: vct(nspin),rhot(nspin) !temp array
      REAL(DP) :: ect

      REAL(DP) :: ec_loc

      INTEGER(I4B) :: Ip
      REAL(DP),DIMENSION(SIZE(sigma,2)) :: sigmat,vsigmat
!
      LOGICAL :: lvc
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      lvc=PRESENT(vcS)
!initial exchange and correlation
      CALL xc_f90_func_init(c_func,c_info,XC_GGA_C_LYP,nspin)
!CALL xc_f90_func_init(c_func,c_info,XC_GGA_XC_XLYP,nspin)
!set zero
      ec=0._DP

      ec_loc=0._DP

!start
      DO Ip=1,nps
!set rho
         rhot(:)=rhoS(Ip,:)
!set sigma
         sigmat(:)=sigma(Ip,:)
!Potential
         IF(lvc)THEN
            CALL xc_f90_gga_vxc(c_func,1,rhot(1),sigmat(1),vct(1),vsigmat(1))
!store
            vcS(Ip,:)=vct(:)
         ENDIF
!energy
         CALL xc_f90_gga_exc(c_func,1,rhot(1),sigmat(1),ect)

         ec_loc=ec_loc+(ect)*SUM(rhot(:))

      ENDDO
!total xc energy

      CALL MPI_ALLREDUCE(ec_loc,ec,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

      ec=ec*dvol
!exit xc call
      CALL xc_f90_func_end(c_func)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Libxc_lyp_c
!------------------------------------------------------------
!B3LYP functional
!------------------------------------------------------------
!------------------------------------------------------------
ENDMODULE xc_module
