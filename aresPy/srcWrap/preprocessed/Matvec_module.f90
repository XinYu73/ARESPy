MODULE matvec_module
  USE constants

  USE smpi_math_module

CONTAINS
!-----------------------for cmplex matver------------------------
  SUBROUTINE cmplx_matvec(Ik,veff1d,p,q,dimen)
     USE finite_module, ONLY : cmplx_keop 
     IMPLICIT NONE
!
     INTEGER(I4B),INTENT(IN) :: Ik,dimen
     REAL(DP),INTENT(IN) :: veff1d(dimen)
     COMPLEX(DCP),INTENT(IN) :: p(dimen)
     COMPLEX(DCP),INTENT(OUT) :: q(dimen)
!> local
     COMPLEX(DCP)     :: kp(dimen) !kinetic operator phi
     INTEGER(I4B) :: I,Ix,Iy,Iz
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!nonlocal part
     CALL cmplx_nlocmatvec(Ik,p,q)
!calculate the kinetic part(see Bloch Theorm.)
     CALL cmplx_keop(p,Ik,kp)
!collection
     q(:)=q(:)+kp(:)+veff1d(:)*p(:)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE cmplx_matvec
!------------------------for real matver------------------------
  SUBROUTINE cmplx_nlocmatvec(Ik,p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot
     USE grid_module ,ONLY : dvol
     IMPLICIT NONE
!IN/OUT
     INTEGER(I4B),INTENT(IN) :: Ik
     COMPLEX(DCP),INTENT(IN)     :: p(:)
     COMPLEX(DCP),INTENT(OUT)    :: q(:)
!LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     COMPLEX(DCP) :: dots(max_nproj,natom),tmp0

     COMPLEX(DCP) :: dots_loc(max_nproj,natom)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
     q(:)=0.d0
     IF(max_nproj==0)RETURN
!
     dots(:,:)=0.d0

     dots_loc(:,:)=0._DP

!all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0) CYCLE
!all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1

           IF(nlpot(Ia)%npts==0) CYCLE

!all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0
!all points
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
!<proj_lm|p>
                 tmp0=tmp0+CONJG(nlpot(Ia)%proj_phs(Ip,Ipj,Ik))*p(Id)
              ENDDO
!save dots
! dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)

              dots_loc(1:psp(Ity)%nproj,Ia)=dots_loc(1:psp(Ity)%nproj,Ia)  &
                       & + tmp0*psp(Ity)%Dij(:,Ipj)

           ENDDO
        ENDDO
     ENDDO

     CALL MPI_ALLREDUCE(dots_loc,dots,SIZE(dots),MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

!scale
     dots(:,:)=dots(:,:)*dvol
!Multiply the nonlocal vectors
     DO Ity=1,naty
!all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1

           IF(nlpot(Ia)%npts==0) CYCLE

!all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
!all points
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
!q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj_phs(Ip,Ipj,Ik)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE cmplx_nlocmatvec
!------------------------for real matver------------------------
  SUBROUTINE real_matvec(veff1d,p,q,dimen)
     USE finite_module, ONLY : real_nabla2_1d
     use grid_module,only:grid
     IMPLICIT NONE
!
     INTEGER(I4B),INTENT(IN) :: dimen
     REAL(DP),INTENT(IN) :: veff1d(dimen)
     REAL(DP),INTENT(IN) :: p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
!> local
     REAL(DP)     :: kp(dimen) !kinetic operator phi
     INTEGER(I4B) :: I,Ix,Iy,Iz
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!nonlocal part
     CALL real_nlocmatvec(p,q)
!calculate the kinetic part(kp=-0.5*nabla2*p)
     CALL real_nabla2_1d(p,kp)
     kp(:)=-0.5_DP*kp(:)
!collection
     q(:)=q(:)+kp(:) +veff1d(:)*p(:)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_matvec
!-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE real_nlocmatvec(p,q)
     USE pspot_module , ONLY : psp ,max_nproj
     USE struct_module , ONLY : naty,natom,struct
     USE nlpot_module , ONLY : nlpot
     USE grid_module ,ONLY : dvol
     IMPLICIT NONE
!IN/OUT
     REAL(DP),INTENT(IN)     :: p(:)
     REAL(DP),INTENT(OUT)    :: q(:)
!LOCAL
     INTEGER(I4B) :: Ity,Ia,Ipj,Ip,Id
     REAL(DP) :: dots(max_nproj,natom),tmp0

     REAL(DP) :: dots_loc(max_nproj,natom)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
     q(:)=0.d0
!return
     IF(max_nproj==0)RETURN
!
     dots(:,:)=0.d0

     dots_loc(:,:)=0._DP

!all type
     DO Ity=1,naty
        IF(psp(Ity)%nproj==0) CYCLE
!all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1

           IF(nlpot(Ia)%npts==0) CYCLE

!all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=0.d0
!all points
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
!<proj_lm|p>
                 tmp0=tmp0+nlpot(Ia)%proj(Ip,Ipj)*p(Id)
              ENDDO
!save dots
! dots(Ipj,Ia)=tmp0*psp(Ity)%D0(Ipj,Ipj)

              dots_loc(1:psp(Ity)%nproj,Ia)=dots_loc(1:psp(Ity)%nproj,Ia)  &
                       & + tmp0*psp(Ity)%Dij(:,Ipj)

           ENDDO
        ENDDO
     ENDDO

     CALL MPI_ALLREDUCE(dots_loc,dots,SIZE(dots),MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)

!scale
     dots(:,:)=dots(:,:)*dvol
!Multiply the nonlocal vectors
     DO Ity=1,naty
!all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1

           IF(nlpot(Ia)%npts==0) CYCLE

!all projectors
           DO Ipj=1,psp(Ity)%nproj
              tmp0=dots(Ipj,Ia)
!all points
              DO Ip=1,nlpot(Ia)%npts
                 Id=nlpot(Ia)%Id(Ip)
!q=\SUM_lm{dots_lm*|proj_lm>}
                 q(Id)=q(Id)+tmp0*nlpot(Ia)%proj(Ip,Ipj)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_nlocmatvec
!-----------------------PARTING-LINE--------------------------
  SUBROUTINE real_matvec_m(mat,p,q,dimen)
     IMPLICIT NONE
!
     INTEGER(I4B),INTENT(IN) :: dimen
     REAL(DP),INTENT(IN) :: mat(dimen,dimen),p(dimen)
     REAL(DP),INTENT(OUT) :: q(dimen)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     q(:)=MATMUL(mat,p)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_matvec_m
!-----------------------PARTING-LINE--------------------------
ENDMODULE matvec_module
