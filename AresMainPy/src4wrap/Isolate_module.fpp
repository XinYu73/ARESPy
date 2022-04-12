# 1 "Isolate_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Isolate_module.f90"
MODULE IsolateSet
  USE constants
  USE parameters,     ONLY: RadiusMax, NVC, norder=>ISOnorder, NFCD, &
       & TOLCG, Lmax=>ISOLmax, iprec_fmm, hartree_method & !IFMM, hartree_direct&
                                !& ,MCM_flag&
       & ,IsoRmaxbohr,Lcellbohr    !not necessary
  USE grid_module
  USE struct_module!,   ONLY:struct.SLmax =>Lmax
  USE m_time_evaluate

  USE poisson_isf , ONLY : poisson_isf_method

  IMPLICIT none
  REAL(DP)              :: Center(3)
  INTEGER(I4B)          :: CellLengthN(3),Thickness,CellLeft,CellRight(3),BandRight(3)
  REAL(DP),ALLOCATABLE  ::BoundaryrhoS(:,:,:),BoundaryVCoulomb(:,:,:),ap(:,:,:),&
       aphi(:,:,:),res(:,:,:),res1(:,:,:)
  !INTEGER(I4B)          ::NVC=6,norder=4,NFCD=4,Lmax=9
  !REAL(DP)              ::TOLCG=1.d-2
  REAL(DP),allocatable   :: VCoulomb_old(:,:,:)
CONTAINS
  !-----------------------DIVIDER-LINE--------------------------



  SUBROUTINE IsolateVCoulomb_a(rho_s,V_s)
    USE smpi_math_module , ONLY : parallel

    USE m_time_evaluate, ONLY: memory_sum,memory_free
    IMPLICIT NONE




    REAL(DP),INTENT(IN)    :: rho_s(:)
    REAL(DP),INTENT(INOUT) :: V_s(:)

    INTEGER(I4B)           :: iter=1,niter_poisson=1000,i
    ! REAL(DP)               :: VCoulomb_old(n1,n2,n3)
    INTEGER(I4B),save :: counter=0
    !> flag for time evaluating
    LOGICAL           :: t_p=.false.  !> time of different poisson solver
    !

    REAL(DP) :: rhoS(n1,n2,n3)
    REAL(DP) :: VCoulomb(n1,n2,n3)

    !===========================================================

    call memory_sum('IsolateVCoulomb',real(size(rhoS),DP)*DP+size(VCoulomb)*DP)
    CALL parallel_s2g(rho_s,rhoS)

    !##TEST FOR LFMM3D LIBRARY
    ! IF(IFMM)THEN
    SELECT CASE(hartree_method)
    CASE(1)
       CALL time_start('FMM time consuming',t_p)
       ! print*,"USE FAST MULTIPOLE METHOD"
       CALL Vhartree_FMM(rhoS,VCoulomb)
       CALL time_end('FMM time consuming',t_p)
       CALL time_output('FMM time consuming',t_p)
       ! ENDIF
       ! IF(hartree_direct)THEN
    CASE(2)
       CALL time_start('hartree direct time consuming',t_p)
       ! print*,"USE FAST MULTIPOLE METHOD"
       CALL Vhartree_direct(rhoS,VCoulomb)
       CALL time_end('hartree direct time consuming',t_p)
       CALL time_output('hartree direct time consuming',t_p)
       ! ENDIF
       ! IF(MCM_flag)THEN
    CASE(3)
       CALL time_start('MCM time consuming',t_p)
       ! print*,"USE FAST MULTIPOLE METHOD"
       CALL poisson_MCM(rhoS,VCoulomb)
       CALL time_end('MCM time consuming',t_p)
       CALL time_output('MCM time consuming',t_p)
       ! ENDIF
       ! IF(MCM_flag)THEN
    CASE(4)
       CALL time_start('Cutoff-Method time consuming',t_p)
       ! print*,"USE FAST MULTIPOLE METHOD"
       CALL cutoff_method(rhoS,VCoulomb)
       CALL time_end('Cutoff-Method time consuming',t_p)
       CALL time_output('Cutoff-Method time consuming',t_p)
    CASE(5)
       CALL time_start('MG-CG-Method time consuming',t_p)
       ! print*,"USE FAST MULTIPOLE METHOD"
       CALL GMG_CG_hart(rhoS,VCoulomb)
       CALL time_end('MG-CG-Method time consuming',t_p)
       CALL time_output('MG-CG-Method time consuming',t_p)
       ! ENDIF
    CASE(6)
       CALL time_start('ISF-Method time consuming',t_p)
       ! print*,"USE FAST MULTIPOLE METHOD"
       CALL poisson_isf_method(rhoS,VCoulomb)
       CALL time_end('ISF-Method time consuming',t_p)
       CALL time_output('ISF-Method time consuming',t_p)

       CALL parallel_g2s(V_s,VCoulomb)

    CASE default
       print*,"poisson-solver parameters err: PLEASE CHECK the input.dat"
       stop
    ENDSELECT


    call memory_free('IsolateVCoulomb',real(size(rhoS),DP)*DP+size(VCoulomb)*DP)

  END SUBROUTINE IsolateVCoulomb_a
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE EvaluateBoundryGrid_b(rhoS)
    implicit none
    REAL(DP),INTENT(IN)  ::rhoS(:,:,:)
    !core size
    CellLengthN(1)=size(rhoS,1)
    CellLengthN(2)=size(rhoS,2)
    CellLengthN(3)=size(rhoS,3)
    !centural position,arbtirary in theory, there is grid center for read_module set
    Center(1)=(CellLengthN(1)+2)*0.5d0*gap(1)
    Center(2)=(CellLengthN(2)+2)*0.5d0*gap(2)
    Center(3)=(CellLengthN(3)+2)*0.5d0*gap(3)
    !shell thickness and size
    Thickness=1
    CellLeft=Thickness  !+1-1  !start at 0
    CellRight(:)=CellLeft+CellLengthN(:)-1   !n+1 dot per grid,so gridright in fact
    BandRight(:)=CellRight(:)+Thickness+1
    allocate(BoundaryrhoS(0:BandRight(1),0:BandRight(2),0:BandRight(3)))
    allocate(BoundaryVCoulomb(0:BandRight(1),0:BandRight(2),0:BandRight(3)))
    BoundaryrhoS=0
    BoundaryrhoS(CellLeft:CellRight(1),CellLeft:CellRight(2),CellLeft:CellRight(3))=rhoS
    !!allocate(dr(CellLengthN(1),CellLengthN(2),CellLengthN(3)))
    !!allocate(cost(CellLengthN(1),CellLengthN(2),CellLengthN(3)))
    !!allocate(sint(CellLengthN(1),CellLengthN(2),CellLengthN(3)))
    !!allocate(cns(CellLengthN(1),CellLengthN(2),CellLengthN(3)))
    !!allocate(indx(CellLengthN(1),CellLengthN(2),CellLengthN(3)))
    !!CALL rcs(CellLengthN(1),CellLengthN(2),CellLengthN(3),dr,cost,sint,cns,indx)
  END SUBROUTINE EvaluateBoundryGrid_b
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE EndEvaluate_b()
    IMPLICIT NONE
    !deallocate the expansion shell
    deallocate(BoundaryrhoS)
    deallocate(BoundaryVCoulomb)
    !------------
    DEALLOCATE(ap,aphi,res,res1)
    !DEALLOCATE(dr,cost,sint,cns,indx)
    !CALL system_clock(it2)
    !       print*,'tot iter time',(it2-it1)/10000.d0
    !       print*,'boundary time',(it4-it3)/10000.d0
    !>>>>>>>>>>>>>>>>>>>END   CGPHI<<<<<<<<<<<<<<<<<<<<<<<
  END SUBROUTINE EndEvaluate_b
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE CalculateBoundryPotential_b(rho,Lmax)!rhoS
    implicit none
    INTEGER(I4B)          :: Lmax
    REAL(DP),INTENT(IN)   :: rho(:,:,:)
    REAL(DP)              :: Q_l0(0:Lmax)
    COMPLEX(DCP)          :: Q_lm(Lmax,Lmax)
    !
    !calculate Qlm
    CALL cal_qlm(rho,Lmax,Q_l0,Q_lm)
    ! calculate the fininate coefficen
    ALLOCATE(ap(BandRight(1)-1,BandRight(2)-1,BandRight(3)-1))  !CorePotential
    !ALLOCATE(bp(0:dnndn+1,0:dnndn+1,0:dnndn+1))
    ALLOCATE(aphi(BandRight(1)-1,BandRight(2)-1,BandRight(3)-1))!CoreVCoulomb
    !ALLOCATE(bphi(0:dnndn+1,0:dnndn+1,0:dnndn+1))
    ALLOCATE(res(BandRight(1)-1,BandRight(2)-1,BandRight(3)-1))
    ALLOCATE(res1(BandRight(1)-1,BandRight(2)-1,BandRight(3)-1))
    !boundary_0(rhos,bf),book bp as bf
    res(:,:,:)=0.d0
    res(CellLeft:CellRight(1),CellLeft:CellRight(2),CellLeft:CellRight(3))=-4*pi*rho(:,:,:)
    !calculate  right hands
    !CALL laplb(NFCD,res,nrhs)
    CALL laplb(NFCD,res,ap)
    !CALL system_clock(it3)
    !center of origini & potential boundary by multi-pole method
    !orig(:)=(dnndn+1)*grid%gaps(1)*0.5d0
    !> compare the direct accumulate and multipole expansion
    ! CALL apply_boundary00(rho,BoundaryVCoulomb)
    ! OPEN (1111,file='boundaryDir')
    ! write(1111,*)BoundaryVCoulomb
    ! CLOSE(1111)
    CALL apply_boundary(Lmax,Center,Q_l0,Q_lm,BoundaryVCoulomb)
    ! OPEN (1111,file='boundaryQlj')
    ! write(1111,*)BoundaryVCoulomb
    ! CLOSE(1111)
    ! stop
    !CALL system_clock(it4)
    !print*,'boundary time',(it4-it3)/10000.d0
  END SUBROUTINE CalculateBoundryPotential_b
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE LaplaceEquationSlover_b(iter)
    USE array_io , ONLY: output
    IMPLICIT NONE
    INTEGER(I4B)     :: iter
    REAL(DP)         :: rr, rrnew
    REAL(DP)         :: err,TOL
    REAL(DP)         :: beta, alpha
    REAL(DP)         :: pap
    INTEGER(I4B)     :: niter_poisson = 1000
    !
    !finite left hand
    CALL lapla(NFCD,gap,BoundaryVCoulomb,aphi)!NFCD
    ! CALL output(size(BoundaryVCoulomb),reshape(BoundaryVCoulomb,(/size(BoundaryVCoulomb)/)),'aphi_n')
    !resduial
    !res(:,:,:)=nrhs(:,:,:)-aphi(:,:,:)
    res(:,:,:)=ap(:,:,:)-aphi(:,:,:)
    err= SQRT(sum(res*res))
    !the tolerance = relative + absolate tolerance
    !TOL=TOLCG*err+1e-4*TOLCG
    TOL=TOLCG
    !---------------------------------------------------------
    IF(err < TOL)THEN
       !DEALLOCATE(ap,BoundaryrhoS,aphi,BoundaryVCoulomb,res,res1)
       RETURN
    ENDIF
    !---------------------------------------------------------
    res1(:,:,:)=0.d0
    CALL v_cycle(NFCD,res1,res,gap,NVC)!^NFCD !dxyz,NVC$
    !res1=res
    !perpare first cg step:p=r
    BoundaryrhoS(:,:,:)=0.d0
    BoundaryrhoS(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1)=res1(:,:,:)
    rr = sum(res*res1)
    ! start cg iter
    DO iter=1,niter_poisson
       CALL lapla(NFCD,gap,BoundaryrhoS,ap)!NFCD,gap
       !step length
       pap = SUM(BoundaryrhoS(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1)*ap(:,:,:))
       alpha = rr / pap
       !update new potential in big area
       BoundaryVCoulomb(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1) = BoundaryVCoulomb(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1) + alpha*BoundaryrhoS(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1)
       !update new resduial
       !if(iter==1)then
       !	 print*,alpha,rr,pap
       !	 open(111,file="ooooo", status="UNKNOWN", action="write")
       !	 write(111,*)BoundaryrhoS(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1)
       !	 close(111)
       !endif
       res(:,:,:) = res(:,:,:) - alpha * ap(:,:,:)
       ! CALL output(size(ap),reshape(ap,(/size(ap)/)),'ap')
       !error evulate
       err = SQRT(sum(res(:,:,:)*res(:,:,:)))
       !print*,err
       if (err < TOL) exit
       res1(:,:,:)=0.d0
       CALL v_cycle(NFCD,res1,res,gap,nvc)!^NFCD,!dxyz,nvc$
       !res1=res
       !beta = (new_r, new_r) / (r, r) & new conjugate
       rrnew= sum(res1(:,:,:) * res(:,:,:))
       beta = rrnew/ rr
       !p = r + beta*p
       !--------------------
       BoundaryrhoS(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1) = res1(:,:,:) + beta*BoundaryrhoS(CellLeft:CellRight(1)+1,CellLeft:CellRight(2)+1,CellLeft:CellRight(3)+1)
       !old rr
       rr = rrnew
    ENDDO
    !===================
    !print*,"ITER",iter
    !===================
  END SUBROUTINE LaplaceEquationSlover_b
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE VCoulombAssignment_b(iter,niter_poisson,phi)
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: iter,niter_poisson
    REAL(DP),INTENT(INOUT)  :: phi(:,:,:)
    !-------------------------------------------------
    IF(iter==niter_poisson)THEN
       WRITE(6,*)'CGPHI:CG FILED'
       STOP
    ENDIF
    ! print *,'CG_steps : ',iter
    phi(:,:,:)=BoundaryVCoulomb(CellLeft:CellRight(1),CellLeft:CellRight(2),CellLeft:CellRight(3))
  END SUBROUTINE VCoulombAssignment_b
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE cal_qlm(rho,Lmax,Q_l0,Q_lm)
    !  USE math , ONLY: cal_plm
    USE math, ONLY: plgndr
    REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: rho
    INTEGER(I4B),INTENT(IN) :: Lmax
    REAL(DP),INTENT(OUT)  :: Q_l0(0:)
    COMPLEX(DCP),INTENT(OUT)  :: Q_lm(:,:)
    INTEGER(I4B) :: ix,iy,iz,n1,n2,n3,l,m
    REAL(DP) :: qrl,hcub,plm(Lmax,0:lmax),clm(Lmax,0:Lmax)
    !Calculate the Qlm
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n1=SIZE(rho,1)
    n2=SIZE(rho,2)
    n3=SIZE(rho,3)
    hcub=gap(1)*gap(2)*gap(3)
    !CALL cal_clm(Lmax,clm)
    !CALL calclm(Lmax,clm)
    CAll  calclm(Lmax,clm)
    Q_l0(:)=0.d0
    !Q_l0(0)=0.d0
    Q_lm(:,:)=(0.d0,0.0)
    !all over the internal grids
    DO iz=1,n3
       DO iy=1,n2
          DO ix=1,n1
             IF(indx(ix,iy,iz)==0) CYCLE
             Q_l0(0)=Q_l0(0)+rho(ix,iy,iz)
             CALL  cal_plm(Lmax,cost(ix,iy,iz),sint(ix,iy,iz),plm)
             !CALL  calplm(Lmax,grid%cost(ix,iy,iz),grid%sint(ix,iy,iz),plm)
             DO l=1,Lmax
                qrl=rho(ix,iy,iz)*hcub*dr(ix,iy,iz)**l
                !Q_l0=\sum_{r'} q'*r'^l*P_l^0(cos\theta')
                Q_l0(l)=Q_l0(l)+qrl*plm(l,0)
                ! Q_l0(l)=Q_l0(l)+qrl*plgndr(l,0,cost(ix,iy,iz))
                !Q_lm=\sum_{r'} c(l,m)*q'*r'^l*P_l^m(cos\theta')*exp(-i*m*\phi)
                DO m=1,l
                   Q_lm(l,m)=Q_lm(l,m)+clm(l,m)*qrl*plm(l,m)*cns(ix,iy,iz)**m
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    Q_l0(0)=Q_l0(0)*hcub
    ! print*, "Ql0",Q_l0
  ENDSUBROUTINE cal_qlm
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE car2spe(ORIG,x,y,z,r,cost,sint,cosp,sinp)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: x,y,z
    REAL(8),INTENT(IN) :: ORIG(3)
    REAL(8),INTENT(OUT) :: r,cost,sint,cosp,sinp
    REAL(8) ::rr,rx,ry,rz
    INTEGER :: l
    !>>>>>>>>>>>>>>>>>>>>>>>start CAR2SPE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    rx=gap(1)*REAL(x,8)-ORIG(1)
    ry=gap(2)*REAL(y,8)-ORIG(2)
    rz=gap(3)*REAL(z,8)-ORIG(3)
    r=SQRT(rx**2+ry**2+rz**2)
    rr=SQRT(rx**2+ry**2)
    IF(r>0.d0)THEN
       cost=rz/r
       sint=rr/r
    ELSE
       !cost=1.d0
       !sint=0.d0
       WRITE(6,*) 'Multi-pole:the origin on grid'
       STOP
    ENDIF

    IF(rr>0.d0)THEN
       cosp=rx/rr
       sinp=ry/rr
    ELSE
       !cosp=1.d0
       !sinp=0.d0
       WRITE(6,*) 'Multi-pole:the origin on grid'
       STOP
    ENDIF
    !<<<<<<<<<<<<<<<<<<<<<<<<<end CAR2SPE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE car2spe
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE cal_plm(Lmax,x,z,plm)
    INTEGER(I4B),INTENT(IN) :: Lmax
    REAL(DP),INTENT(IN)     :: x,z
    REAL(DP),INTENT(OUT)    :: plm(Lmax,0:Lmax)
    REAL(DP)    :: p(0:Lmax,0:Lmax)
    INTEGER(I4B) :: l,m
    REAL(DP)     :: fact,pmm,y

    p=0.d0
    pmm=1.d0
    fact=1.d0
    p(0,0)=1.d0
    y=abs(z)
    do m=1,Lmax
       pmm=-pmm*fact*z
       p(m,m)=pmm
       fact=fact+2.d0
    enddo

    !if(l.eq.m) then
    !   p=pmm
    !else
    DO l=0,Lmax-1
       p(l+1,l)=x*(2*l+1)*p(l,l)
    ENDDO

    DO m=0,Lmax-2
       do l=m+2,Lmax
          p(l,m)=(x*(2*l-1)*p(l-1,m)-(l+m-1)*p(l-2,m))/(l-m)
       enddo
    ENDDO
    plm(:,:)=p(1:Lmax,:)
  ENDSUBROUTINE cal_plm
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE cal1pot(Iex,Iey,Iez,orig,Lmax,Q_l0,Q_lm,phi)
    !USE math , ONLY: cal_plm
    INTEGER(I4B),INTENT(IN) :: Iex,Iey,Iez,Lmax
    REAL(DP),INTENT(IN) :: Q_l0(0:),orig(:)
    COMPLEX(DCP),INTENT(IN)  :: Q_lm(:,:)
    REAL(DP),INTENT(OUT) :: phi
    INTEGER(I4B) :: l,m
    REAL(DP) :: invral1,ra,plm(Lmax,0:Lmax)
    REAL(DP) :: COSthe,SINthe,COSphi,SINphi
    REAL(DP) :: b00,bp0,bp1
    COMPLEX(DCP)  :: bpm
    !
    !b00=0.d0
    bp0=0.d0
    bp1=0.d0
    CALL car2spe(ORIG,Iex,Iey,Iez,ra,COSthe,SINthe,COSphi,SINphi)
    CALL cal_plm(Lmax,COSthe,SINthe,plm)
    !CALL calplm(Lmax,COSthe,SINthe,plm)
    DO l=1,Lmax
       invral1=1.d0/(ra**(l+1))
       !bp0
       bp0=bp0+invral1*plm(l,0)*Q_l0(l)
       !Bp1
       bpm=(0.d0,0.d0)
       DO m=1,l
          bpm=bpm+plm(l,m)*cmplx(COSphi,-SINphi)**m*Q_lm(l,m)
       ENDDO
       bp1=bp1+invral1*2.d0*REAL(bpm,8)
    ENDDO
    b00=Q_l0(0)/ra
    phi=b00+bp0+bp1
  ENDSUBROUTINE cal1pot
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE apply_boundary(Lmax,orig,Q_l0,Q_lm,bphi)
    INTEGER(I4B),INTENT(IN) :: Lmax
    REAL(DP),INTENT(IN) :: Q_l0(0:),orig(3)
    COMPLEX(DCP),INTENT(IN)  :: Q_lm(:,:)
    REAL(DP),DIMENSION(0:,0:,0:),INTENT(OUT) :: bphi
    INTEGER(I4B) :: Iex,Iey,Iez,n1,n2,n3
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bphi(:,:,:)=0.d0
    !-100-
    n1=SIZE(bphi,1)-2
    n2=SIZE(bphi,2)-2
    n3=SIZE(bphi,3)-2
    DO Iez=0,n3+1
       DO Iey=0,n2+1
          CALL cal1pot(0,Iey,Iez,orig,Lmax,Q_l0,Q_lm,bphi(0,Iey,Iez))
          CALL cal1pot(n1+1,Iey,Iez,orig,Lmax,Q_l0,Q_lm,bphi(n1+1,Iey,Iez))
       ENDDO
    ENDDO
    !-010-
    DO Iez=0,n3+1
       DO Iex=1,n1
          CALL cal1pot(Iex,0,Iez,orig,Lmax,Q_l0,Q_lm,bphi(Iex,0,Iez))
          CALL cal1pot(Iex,n2+1,Iez,orig,Lmax,Q_l0,Q_lm,bphi(Iex,n2+1,Iez))
       ENDDO
    ENDDO
    !-001-
    DO Iey=1,n2
       DO Iex=1,n1
          CALL cal1pot(Iex,Iey,0,orig,Lmax,Q_l0,Q_lm,bphi(Iex,Iey,0))
          CALL cal1pot(Iex,Iey,n3+1,orig,Lmax,Q_l0,Q_lm,bphi(Iex,Iey,n3+1))
       ENDDO
    ENDDO
  ENDSUBROUTINE apply_boundary
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE apply_boundary00(rho,bphi)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)              :: rho(:,:,:)
    REAL(DP),DIMENSION(0:,0:,0:),INTENT(OUT) :: bphi
    INTEGER(I4B) :: Iex,Iey,Iez,n1,n2,n3

    bphi(:,:,:)=0.d0
    !-100-
    n1=SIZE(bphi,1)-2
    n2=SIZE(bphi,2)-2
    n3=SIZE(bphi,3)-2
    DO Iez=0,n3+1
       DO Iey=0,n2+1
          CALL cal_Vsrcpot(size(rho),rho,grid%rVec(1:3,:),1,bphi(0:0,Iey,Iez),(/0.d0*gap(1),Iey*gap(1),Iez*gap(1)/))
          CALL cal_Vsrcpot(size(rho),rho,grid%rVec(1:3,:),1,bphi(n1+1:n1+1,Iey,Iez),(/(n1+1)*gap(1),Iey*gap(1),Iez*gap(1)/))
          ! CALL cal1pot(n1+1,Iey,Iez,orig,Lmax,Q_l0,Q_lm,bphi(n1+1,Iey,Iez))
       ENDDO
    ENDDO
    !-010-
    DO Iez=0,n3+1
       DO Iex=1,n1
          CALL cal_Vsrcpot(size(rho),rho,grid%rVec(1:3,:),1,bphi(Iex,0:0,Iez),(/Iex*gap(1),0*gap(1),Iez*gap(1)/))
          CALL cal_Vsrcpot(size(rho),rho,grid%rVec(1:3,:),1,bphi(Iex,n2+1:n2+1,Iez),(/Iex*gap(1),(n2+1)*gap(1),Iez*gap(1)/))
          ! CALL cal1pot(Iex,0,Iez,orig,Lmax,Q_l0,Q_lm,bphi(Iex,0,Iez))
          ! CALL cal1pot(Iex,n2+1,Iez,orig,Lmax,Q_l0,Q_lm,bphi(Iex,n2+1,Iez))
       ENDDO
    ENDDO
    !-001-
    DO Iey=1,n2
       DO Iex=1,n1
          CALL cal_Vsrcpot(size(rho),rho,grid%rVec(1:3,:),1,bphi(Iex,Iey,0:0),(/Iex*gap(1),Iey*gap(1),0*gap(1)/))
          CALL cal_Vsrcpot(size(rho),rho,grid%rVec(1:3,:),1,bphi(Iex,Iey,n3+1:n3+1),(/Iex*gap(1),Iey*gap(1),(n3+1)*gap(1)/))
          ! CALL cal1pot(Iex,Iey,0,orig,Lmax,Q_l0,Q_lm,bphi(Iex,Iey,0))
          ! CALL cal1pot(Iex,Iey,n3+1,orig,Lmax,Q_l0,Q_lm,bphi(Iex,Iey,n3+1))
       ENDDO
    ENDDO

  ENDSUBROUTINE apply_boundary00
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE laplb(Ford,bf,f)
    USE constants
    !USE finite_module, ONLY:laplb_csix,laplb_cfour
    INTEGER(I4B),INTENT(IN) :: Ford
    REAL(DP),INTENT(IN)     :: bf(:,:,:)
    REAL(DP),INTENT(OUT)    :: f(:,:,:)
    IF(Ford==6)THEN
       CALL laplb_csix(bf,f)
    ELSEIF(Ford==4)THEN
       CALL laplb_cfour(bf,f)
    ELSE
       WRITE(6,*) "STOP: Ford=4 or 6 for compact finite difference"
       STOP
    ENDIF
  ENDSUBROUTINE laplb
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE lapla(Ford,dxyz,f,af)
    USE constants
    !USE finite_module,ONLY:csix_coe,cfour_coe,lapla_csix,lapla_cfour
    INTEGER(I4B),INTENT(IN) :: Ford
    REAL(DP),INTENT(IN) :: dxyz(3)
    REAL(DP),INTENT(IN)     :: f(0:,0:,0:)
    REAL(DP),INTENT(OUT)    :: af(:,:,:)
    REAL(DP)                :: coe(0:3)
    IF(Ford==6)THEN
       CALL csix_coe(dxyz,coe)
       CALL lapla_csix(f,af,coe)
    ELSEIF(Ford==4)THEN
       CALL cfour_coe(dxyz,coe)
       CALL lapla_cfour(f,af,coe)
    ELSE
       WRITE(6,*) "STOP: Ford=4 or 6 for compact finite difference"
       STOP
    ENDIF
  ENDSUBROUTINE lapla
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE lapla_coe(Ford,dxyz,coe)
    USE constants
    !USE finite_module
    INTEGER(I4B),INTENT(IN) :: Ford
    REAL(DP),INTENT(IN) :: dxyz(3)
    REAL(DP) :: coe(0:3)
    IF(Ford==6)THEN
       CALL csix_coe(dxyz,coe)
    ELSEIF(Ford==4)THEN
       CALL cfour_coe(dxyz,coe)
    ELSE
       WRITE(6,*) "STOP: Ford=4 or 6 for compact finite difference"
       STOP
    ENDIF
  ENDSUBROUTINE lapla_coe
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine damped_jacobi_iterate(Ford,dxyz,func,rhs)!{{{
    !
    !  One Gauss-Seidel iteration sweep for the Laplace equation
    !
    use constants
    !
    implicit none
    !
    INTEGER(I4B),INTENT(IN) :: Ford
    real(dp),intent(inout)  :: func(:,:,:)
    real(dp),intent(in) :: rhs(:,:,:)
    REAL(DP)             :: dxyz(3)
    !
    REAL(DP),ALLOCATABLE :: vec(:,:,:)
    !
    INTEGER(I4B)         :: i,ish,ix,iy,iz
    real(dp)             :: sumf, normf
    INTEGER(I4B)                  :: n1,n2,n3
    REAL(DP)                      :: coe(0:3)
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    n1= size(func,1)
    n2= size(func,2)
    n3= size(func,3)
    ALLOCATE(vec(0:n1+1,0:n2+1, 0:n3+1) )
    !
    !call apply_pbc_19(func, vec)
    vec(:,:,:)=0.d0
    vec(1:n1,1:n2,1:n3)=func(:,:,:)
    !
    CALL lapla_coe(Ford,dxyz, coe)
    !
    normf = coe(0)
    !$OMP PARALLEL DO PRIVATE(iz,iy,ix,ish)
    do iz=1,n3
       do iy=1,n2
          do ix=1,n1
             sumf = -rhs(ix,iy,iz) + &
                  coe(1) *( (vec(ix-1,iy,iz)+vec(ix+1,iy,iz)) + &
                  (vec(ix,iy-1,iz)+vec(ix,iy+1,iz)) + &
                  (vec(ix,iy,iz-1)+vec(ix,iy,iz+1))) + &
                  coe(2) *((vec(ix-1,iy-1,iz)+vec(ix+1,iy-1,iz)+vec(ix-1,iy+1,iz)+vec(ix+1,iy+1,iz)) + &
                  (vec(ix,iy-1,iz-1)+vec(ix,iy+1,iz-1)+vec(ix,iy-1,iz+1)+vec(ix,iy+1,iz+1)) + &
                  (vec(ix-1,iy,iz-1)+vec(ix+1,iy,iz-1)+vec(ix-1,iy,iz+1)+vec(ix+1,iy,iz+1)))+ &
                  coe(3) *((vec(ix-1,iy-1,iz-1)+vec(ix-1,iy-1,iz+1)+vec(ix-1,iy+1,iz-1)+vec(ix-1,iy+1,iz+1))+&
                  vec(ix+1,iy-1,iz-1)+vec(ix+1,iy-1,iz+1)+vec(ix+1,iy+1,iz-1)+vec(ix+1,iy+1,iz+1) )
             !
             func(ix,iy,iz) = sumf / normf
             !func(ix,iy,iz) = (func(ix,iy,iz) - 2.d0/3.d0*sumf) / normf
             !func(ix,iy,iz) = 2.d0/3.d0*sumf / normf
             !func(ix,iy,iz) = (2.d0/3.d0 * sumf - rhs(ix,iy,iz))/ normf
          enddo
       enddo
    ENDDO
    !$OMP END PARALLEL DO
    IF (ALLOCATED(vec)) DEALLOCATE(vec)
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
  END SUBROUTINE damped_jacobi_iterate !}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine gauss_seidel_iterate(Ford,dxyz,func,rhs)!{{{
    !
    !  One Gauss-Seidel iteration sweep for the Laplace equation
    !
    use constants
    !
    implicit none
    !
    INTEGER(I4B),INTENT(IN) :: Ford
    real(dp),intent(inout)  :: func(:,:,:)
    real(dp),intent(in) :: rhs(:,:,:)
    REAL(DP)             :: dxyz(3)
    !
    REAL(DP),ALLOCATABLE :: vec(:,:,:)
    !
    INTEGER(I4B)         :: i,ish,ix,iy,iz
    real(dp)             :: sumf, normf
    INTEGER(I4B)                  :: n1,n2,n3
    REAL(DP)                      :: coe(0:3)
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    n1= size(func,1)
    n2= size(func,2)
    n3= size(func,3)
    ALLOCATE(vec(0:n1+1,0:n2+1, 0:n3+1) )
    !
    !call apply_pbc_19(func, vec)
    vec(:,:,:)=0.d0
    vec(1:n1,1:n2,1:n3)=func(:,:,:)
    !
    CALL lapla_coe(Ford,dxyz, coe)
    !
    normf = coe(0)
    !$OMP PARALLEL DO PRIVATE(iz,iy,ix,ish)
    do iz=1,n3
       do iy=1,n2
          do ix=1,n1
             sumf = -rhs(ix,iy,iz) + &
                  coe(1) *((vec(ix-1,iy,iz)+vec(ix+1,iy,iz)) + &
                  (vec(ix,iy-1,iz)+vec(ix,iy+1,iz)) + &
                  (vec(ix,iy,iz-1)+vec(ix,iy,iz+1))) + &
                  coe(2) *((vec(ix-1,iy-1,iz)+vec(ix+1,iy-1,iz)+vec(ix-1,iy+1,iz)+vec(ix+1,iy+1,iz)) + &
                  (vec(ix,iy-1,iz-1)+vec(ix,iy+1,iz-1)+vec(ix,iy-1,iz+1)+vec(ix,iy+1,iz+1)) + &
                  (vec(ix-1,iy,iz-1)+vec(ix+1,iy,iz-1)+vec(ix-1,iy,iz+1)+vec(ix+1,iy,iz+1)))+ &
                  coe(3) *((vec(ix-1,iy-1,iz-1)+vec(ix-1,iy-1,iz+1)+vec(ix-1,iy+1,iz-1)+vec(ix-1,iy+1,iz+1))+&
                  vec(ix+1,iy-1,iz-1)+vec(ix+1,iy-1,iz+1)+vec(ix+1,iy+1,iz-1)+vec(ix+1,iy+1,iz+1) )
             !
             vec(ix,iy,iz) = sumf / normf
             !func(ix,iy,iz) = sumf / normf
             !func(ix,iy,iz) = (func(ix,iy,iz) - 2.d0/3.d0*sumf) / normf
             !func(ix,iy,iz) = 2.d0/3.d0*sumf / normf
             !func(ix,iy,iz) = (2.d0/3.d0 * sumf - rhs(ix,iy,iz))/ normf
          enddo
       enddo
    ENDDO
    !$OMP END PARALLEL DO
    func(1:n1,1:n2,1:n3) = vec(1:n1,1:n2,1:n3)
    IF (ALLOCATED(vec)) DEALLOCATE(vec)
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
  END SUBROUTINE gauss_seidel_iterate !}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine refine(var, nxyz, nxyz_fine, var_fine)!{{{
    !
    !  Refine (i.e. fine-grain by interpolating) f from original to finer grid
    !
    real(DP),dimension(:,:,:) :: var, var_fine
    integer(I4B), dimension(3)  :: nxyz, nxyz_fine
    !
    intent(in)  :: var, nxyz, nxyz_fine
    intent(out) :: var_fine
    !
    !
    ! We simply use trilinear interpolation, which works fine
    !
    call interpolate(var, nxyz, nxyz_fine, var_fine)
    !
  END SUBROUTINE refine!}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine interpolate(var1, nxyz1, nxyz2, var2)!{{{
    !
    !  Trilinear interpolation
    !
    real(DP),dimension(:,:,:) :: var1, var2
    integer(I4B), dimension(3)  :: nxyz1, nxyz2
    real(DP)                    :: idx,idy,idz, rx,ry,rz
    integer(I4B)               :: ix,iy,iz
    integer(I4B)               :: jx,jy,jz
    !
    intent(in)  :: var1, nxyz1, nxyz2
    intent(out) :: var2
    REAL(DP)                      :: scalen(3)
    REAL(DP),ALLOCATABLE          :: vec(:,:,:)
    !
    INTEGER(I4B)                  :: nbor =1 ! only need one one the boundary
    REAL(DP)                      :: safety_factor
    !
    ! wd: TO DO: Precalculate idxx(:), rxx1(:), etc.
    !
    IF (ALLOCATED(vec)) DEALLOCATE(vec)
    ALLOCATE(vec(-nbor+1:nxyz1(1)+nbor, -nbor+1:nxyz1(2)+nbor,-nbor+1:nxyz1(3)+nbor))

    safety_factor = 1 + 2*epsilon(1.d0) ! ensure floor(1.d0) = 1

    DO ix =1, 3
       scalen(ix) = real(nxyz1(ix)) / real(nxyz2(ix))
       !scalen(ix) = real(nxyz1(ix)-1) / real(nxyz2(ix)-1)
    ENDDO
    !PRINT *,"scalen",scalen
    !
    !call apply_pbc_27(var1, vec)
    vec(:,:,:)=0.d0
    vec(1:nxyz1(1),1:nxyz1(2),1:nxyz1(3))=var1(:,:,:)
    !
    do jz=1,nxyz2(3)
       idz = (jz-1) * scalen(3) + safety_factor
       iz = floor(idz);  rz = idz - iz
       do jy=1,nxyz2(2)
          idy = (jy-1) * scalen(2) + safety_factor
          iy = floor(idy);  ry = idy - iy
          do jx=1,nxyz2(1)
             idx = (jx-1) * scalen(1) + safety_factor
             !
             ix = floor(idx);  rx = idx - ix
             !
             ! Interpolation formula
             !
             var2(jx,jy,jz) &
                  =   vec(ix  , iy  , iz  ) * (1-rx) * (1-ry) * (1-rz) &
                                !
                  + vec(ix+1, iy  , iz  ) * rx     * (1-ry) * (1-rz) &
                  + vec(ix  , iy+1, iz  ) * (1-rx) * ry     * (1-rz) &
                  + vec(ix  , iy  , iz+1) * (1-rx) * (1-ry) * rz     &
                                !
                  + vec(ix+1, iy+1, iz  ) * rx     * ry     * (1-rz) &
                  + vec(ix  , iy+1, iz+1) * (1-rx) * ry     * rz     &
                  + vec(ix+1, iy  , iz+1) * rx     * (1-ry) * rz     &
                                !
                  + vec(ix+1, iy+1, iz+1) * rx     * ry     * rz
             !-----------------------------------------------------------------------
          enddo
       enddo
    ENDDO
    IF (ALLOCATED(vec)) DEALLOCATE(vec)
    !
  END SUBROUTINE interpolate !}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine smooth(f)!{{{
    !
    !  Apply `full weighting' smoothing ([1,2,1]/4 weighting for each
    !  direction) to var.
    !
    real(DP),dimension(:,:,:) :: f
    integer(I4B)               :: n1, n2, n3
    !
    intent(inout) :: f
    !
    n1 = size(f,1)
    n2 = size(f,2)
    n3 = size(f,3)
    !
    f(2:n1-1, :, :)  &
         = 0.25 * (    f(1:n1-2, :, :) &
         + 2*f(2:n1-1, :, :) &
         +   f(3:n1, :, :) )
    f(1, :, :)  &
         = 0.25 * (    f(1, :, :) &
         + 2*f(2, :, :) &
         + 0.d0  )!f(n1, :, :) )
    f(n1, :, :)  &
         = 0.25 * ( 0.d0 &   !f(1, :, :) &
         + 2*f(n1-1, :, :) &
         +   f(n1, :, :) )
    !
    f(:, 2:n2-1, :) &
         = 0.25 * (    f(: ,1:n2-2, :) &
         + 2*f(: ,2:n2-1, :) &
         +    f(: ,3:n2, :) )
    f(:, 1, :) &
         = 0.25 * (    f(: ,1, :) &
         + 2*f(: ,2, :) &
         + 0.d0  )! f(: ,n2, :) )
    f(:, n2, :) &
         = 0.25 * ( 0.d0 &  ! f(: ,1, :) &
         + 2*f(: ,n2-1, :) &
         +   f(: ,n2, :) )
    !
    f(:, :, 2:n3-1) &
         = 0.25 * (    f(:, :, 1:n3-2) &
         + 2*f(:, :, 2:n3-1) &
         +   f(:, :, 3:n3) )
    f(:, :, 1) &
         = 0.25 * (  f(:, :, 1) &
         + 2*f(:, :, 2) &
         + 0.d0 )  ! f(:, :, n3) )
    f(:, :, n3) &
         = 0.25 * ( 0.d0 &  ! f(:, :, 1) &
         + 2*f(:, :, n3-1) &
         +   f(:, :, n3) )
    !
  END SUBROUTINE smooth!}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine restrict(var, nxyz, nxyz_coarse, var_coarse)!{{{
    !
    !  Restrict, (i.e. coarse-grain) var from the original to coarser grid.
    !  We use `full weighting' on the original grid, followed by trilinear
    !  interpolation to the (cell-centered) coarse grid we really want.
    !
    !  This could be `chunked' (i.e. applied in blocks) to save some memory.
    !
    real(DP),dimension(:,:,:)          :: var, var_coarse
    integer(I4B), dimension(3)           :: nxyz, nxyz_coarse
    real(DP),dimension(size(var,1),size(var,2),size(var,3)) :: tmp
    !
    intent(in)  :: var, nxyz, nxyz_coarse
    intent(out) :: var_coarse
    !
    ! Smooth, then interpolate
    !
    tmp = var
    call smooth(tmp)
    call interpolate(tmp, nxyz, nxyz_coarse, var_coarse)
    !
  END SUBROUTINE restrict!}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  subroutine residual(Ford,dxyz,f,rhs, res)!{{{
    !
    !  Calculate the residual r = L*f - rhs
    !
    INTEGER(I4B),INTENT(IN)             :: Ford
    real(DP),dimension(:,:,:)           :: f
    real(DP),dimension(:,:,:)           :: rhs,res
    real(DP),dimension(3)               :: dxyz
    !
    intent(in)  :: f, rhs, dxyz
    intent(out) :: res
    INTEGER(I4B)                  :: n1,n2,n3
    real(dp),DIMENSION(0:SIZE(f,1)+1,0:SIZE(f,2)+1,0:SIZE(f,3)+1)  :: bf
    !
    n1=SIZE(f,1)
    !CALL llapa_cfour_pbc(f,res,dxyz,n1,n2,n3)
    bf(:,:,:)=0.d0
    bf(1:n1,1:n1,1:n1)=f(:,:,:)
    CALL lapla(Ford,dxyz,bf,res)
    res = rhs - res
  END SUBROUTINE residual!}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  RECURSIVE SUBROUTINE v_cycle(Ford,f,rhs, dxyz, nvc)!{{{
    !
    !  Do one full multigrid V cycle, i.e. go from finest to coarsest grid
    !  and back:
    !
    !    h            h
    !     \          /
    !     2h        2h
    !       \      /
    !        4h   4h
    !         \  /
    !          8h
    !
    !  to get an approximate solution f of (\Laplace ) f = rhs.
    !
    !  This subroutine calls itself recursively, so from the viewpoint of the
    !  numerical grid the V cycle in question may as well be
    !
    !     2h        2h
    !       \      /
    !        4h   4h
    !         \  /
    !          8h
    !
    !  or
    !
    !        4h   4h
    !         \  /
    !          8h
    !
    !
    !  18-may-2007/wolf: coded
    !
    !USE parameters,        only : norder=>finite_order
    INTEGER(I4B),INTENT(IN)             :: Ford
    real(DP),dimension(:,:,:)           :: f,rhs
    real(DP),dimension(size(f,1),size(f,2),size(f,3)) :: res,df
    real(DP),allocatable, dimension(:,:,:) :: res_coarse, c_coarse, df_coarse
    real(DP),dimension(3)                  :: dxyz,dxyz_coarse
    integer(I4B), dimension(3)               :: nxyz, nxyz_coarse
    INTEGER(I4B)                  :: nvc
    INTEGER(I4B)                  :: nvc1
    !
    intent(in)    :: rhs, dxyz
    intent(inout) :: f
    !
    nxyz = (/ size(f,1), size(f,2), size(f,3) /)
    !
    ! Residual
    !
    call residual(Ford,dxyz,f, rhs, res)
    !WRITE(6,'(A,I4,": ",3(", ",1pG10.3)," ",L1)') &
    !'nxyz(1), residual (min,max,rms) : ', &
    !nxyz(1), minval(res), maxval(res), sqrt(sum(res**2)/size(res))
    ! Set up coarse grid
    !
    nxyz_coarse = (nxyz+1)/2  ! so 3 -> 2, not 3 -> 1
    nxyz_coarse = max(nxyz_coarse, norder)!norder
    dxyz_coarse = dxyz * nxyz / nxyz_coarse
    !
    allocate(res_coarse(nxyz_coarse(1), nxyz_coarse(2), nxyz_coarse(3)))
    allocate(c_coarse(  nxyz_coarse(1), nxyz_coarse(2), nxyz_coarse(3)))
    allocate(df_coarse( nxyz_coarse(1), nxyz_coarse(2), nxyz_coarse(3)))
    !
    ! ..and restrict to coarse grid
    !
    !PRINT *,"rest"
    !OPEN(UNIT=121, FILE="dens")
    !WRITE(121,*)nxyz
    !WRITE(121,*) res
    !close(121)
    call restrict(res, nxyz, nxyz_coarse, res_coarse)
    !PRINT *,"coarse"
    !OPEN(UNIT=121, FILE="coarse")
    !WRITE(121,*)nxyz_coarse
    !WRITE(121,*) res_coarse
    !close(121)
    !pause
    !
    ! Do one relaxation step (i.e. damp larger-scale errors on coarse grid)
    !
    df_coarse = 0.d0            ! initial value
    nvc1 = nvc -1
    !PRINT *,"nvc",nvc

    if (nvc > 0) then
       !if (any(nxyz > 10)) then
       call v_cycle(Ford,df_coarse, res_coarse, dxyz_coarse, nvc1)
    endif
    !
    ! Refine (interpolate) coarse-grid correction to finer grid, and
    ! apply
    !
    ! TODO:
    !   Can be memory-optimized to immediately add df(i,j,k) to f(i,j,k)
    !   -> no df needed
    call refine(df_coarse, nxyz_coarse, nxyz, df)
    f = f - df
    !PRINT *,"df_coarse"
    !OPEN(UNIT=121, FILE="dfcoarse")
    !WRITE(121,*)nxyz_coarse
    !WRITE(121,*) df_coarse
    !close(121)
    !!PRINT *,"df"
    !OPEN(UNIT=121, FILE="df")
    !WRITE(121,*)nxyz
    !WRITE(121,*) df
    !close(121)
    !pause
    !
    ! Do one smoothing iteration (i.e. damp small-scale errors)
    !
    !call gauss_seidel_iterate(f, rhs, dxyz)
    call damped_jacobi_iterate(Ford,dxyz,f, rhs)
    !
    deallocate(res_coarse)
    deallocate(c_coarse)
    deallocate(df_coarse)
    !
  END SUBROUTINE v_cycle!}}}
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE csix_coe(dxyz,coe)
    !
    USE constants
    IMPLICIT NONE
    REAL(DP),INTENT(IN)  :: dxyz(0:)
    REAL(DP),INTENT(OUT) :: coe(0:3)
    REAL(DP)             :: temp
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    temp=1/(30.d0*dxyz(1)*dxyz(1))
    coe(0)=-128.d0*temp
    coe(1)=14.d0*temp
    coe(2)=3.d0*temp
    coe(3)=temp
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE csix_coe
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE lapla_csix(f,af,coe)
    !
    USE constants
    IMPLICIT NONE
    !
    REAL(DP),INTENT(IN)  :: f(0:,0:,0:)
    REAL(DP),INTENT(OUT) :: af(:,:,:)
    REAL(DP),INTENT(IN)  :: coe(0:)
    INTEGER(I4B)         :: nn,ix,iy,iz
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nn = SIZE(af,1)
    DO iz=1,nn
       DO iy=1,nn
          DO ix=1,nn
             af(ix,iy,iz)= coe(0)*f(ix,iy,iz) +                       &
                  &    coe(1)*(                                    &
                  &            f(ix-1,iy,iz)+f(ix+1,iy,iz)       &
                  &        +   f(ix,iy-1,iz)+f(ix,iy+1,iz)       &
                  &        +   f(ix,iy,iz-1)+f(ix,iy,iz+1)       &
                  &           )            +                       &
                  &   coe(2)*(                                    &
                  & f(ix-1,iy-1,iz)+f(ix-1,iy+1,iz)+f(ix+1,iy-1,iz)+f(ix+1,iy+1,iz)&
                  &+f(ix-1,iy,iz-1)+f(ix-1,iy,iz+1)+f(ix+1,iy,iz-1)+f(ix+1,iy,iz+1)&
                  &+f(ix,iy-1,iz-1)+f(ix,iy-1,iz+1)+f(ix,iy+1,iz-1)+f(ix,iy+1,iz+1)&
                  &          )            +                       &
                  &   coe(3)*(                                    &
                  &           f(ix-1,iy-1,iz-1)+f(ix-1,iy-1,iz+1) &
                  &       +   f(ix-1,iy+1,iz-1)+f(ix-1,iy+1,iz+1) &
                  &       +   f(ix+1,iy-1,iz-1)+f(ix+1,iy-1,iz+1) &
                  &       +   f(ix+1,iy+1,iz-1)+f(ix+1,iy+1,iz+1) &
                  &          )
          ENDDO
       ENDDO
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapla_csix
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE laplb_csix(bff,f)
    !
    USE constants
    IMPLICIT NONE
    !
    REAL(DP),INTENT(IN)  :: bff(:,:,:)
    REAL(DP),INTENT(OUT) :: f(:,:,:)
    REAL(DP)             :: coe(0:3)
    INTEGER(I4B)         :: nn,ix,iy,iz
    REAL(DP),ALLOCATABLE :: bf(:,:,:)
    !
    nn=SIZE(f,1)
    ALLOCATE(bf(-1:nn+2,-1:nn+2,-1:nn+2))
    bf(:,:,:)=0.d0
    bf(1:nn,1:nn,1:nn)=bff(:,:,:)
    coe(0)=41.d0/60.d0
    coe(1)=1.d0/36.d0
    coe(2)=1.d0/360.d0
    coe(3)=1.d0/90.d0
    DO iz=1,nn
       DO iy=1,nn
          DO ix=1,nn
             f(ix,iy,iz)=coe(0)*bf(ix,iy,iz)+                           &
                  &     coe(1)*(                                       &
                  &             bf(ix+1,iy,iz)+bf(ix-1,iy,iz)          &
                  &  +          bf(ix,iy+1,iz)+bf(ix,iy-1,iz)          &
                  &  +          bf(ix,iy,iz+1)+bf(ix,iy,iz-1)          &
                  &             )          +                           &
                  &     coe(2)*(&
                  &             bf(ix+2,iy,iz)+bf(ix-2,iy,iz)          &
                  &  +          bf(ix,iy+2,iz)+bf(ix,iy-2,iz)          &
                  &  +          bf(ix,iy,iz+2)+bf(ix,iy,iz-2)          &
                  &             )          +                           &
                  &    coe(3)*(&
                  & bf(ix+1,iy+1,iz)+bf(ix+1,iy-1,iz)+bf(ix-1,iy+1,iz)+bf(ix-1,iy-1,iz)&
                  &+bf(ix,iy+1,iz+1)+bf(ix,iy+1,iz-1)+bf(ix,iy-1,iz+1)+bf(ix,iy-1,iz-1)&
                  &+bf(ix+1,iy,iz+1)+bf(ix+1,iy,iz-1)+bf(ix-1,iy,iz+1)+bf(ix-1,iy,iz-1)&
                  &           )
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(bf)
  ENDSUBROUTINE laplb_csix
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE cfour_coe(gaps,coe)
    !
    REAL(DP)                      :: coe(0:3)
    REAL(DP),INTENT(IN)            :: gaps(3)
    INTEGER(I4B)                  :: i
    REAL(DP)                      :: temp
    !>>>>>>>>>>>>>>>>>>>BEGIN CFOURE_COE>>>>>>>>>>>>>>>>>>>>>>
    temp=1.d0/(6.d0*gaps(1)*gaps(1))
    coe(0) = -24.d0*temp
    coe(1) = 2.d0*temp
    coe(2) = temp
    coe(3) = 0.d0
    !<<<<<<<<<<<<<<<<<<<<END CFOURE_COE<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE cfour_coe
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE lapla_cfour(u,au,coe)
    REAL(DP),INTENT(IN)     :: u(0:,0:,0:)
    REAL(DP),INTENT(IN)     :: coe(0:3)

    REAL(DP),INTENT(OUT)    :: au(:,:,:)
    INTEGER(I4B)            :: nn,ix,iy,iz
    !>>>>>>>>>>>>>>>>>>>>>>>>BEGIN LAPLA_CFOUR>>>>>>>>>>>>>>>>>>>>>>>>
    nn=SIZE(au,1)
    DO iz=1,nn
       DO iy=1,nn
          DO ix=1,nn
             au(ix,iy,iz) = coe(0) *u(ix,iy,iz) + &
                  coe(1) *((u(ix-1,iy,iz)+u(ix+1,iy,iz)) + &
                  (u(ix,iy-1,iz)+u(ix,iy+1,iz)) + &
                  (u(ix,iy,iz-1)+u(ix,iy,iz+1)) ) + &
                  coe(2)*((u(ix-1,iy-1,iz)+u(ix+1,iy-1,iz)+u(ix-1,iy+1,iz)+u(ix+1,iy+1,iz)) + &
                  (u(ix,iy-1,iz-1)+u(ix,iy+1,iz-1)+u(ix,iy-1,iz+1)+u(ix,iy+1,iz+1)) + &
                  (u(ix-1,iy,iz-1)+u(ix+1,iy,iz-1)+u(ix-1,iy,iz+1)+u(ix+1,iy,iz+1)))
          ENDDO
       ENDDO
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<END   LAPLA_CFOUR<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE lapla_cfour
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE laplb_cfour(rhs,nrhs)
    REAL(DP),INTENT(IN)  :: rhs(:,:,:)
    REAL(DP),INTENT(OUT) :: nrhs(:,:,:)

    INTEGER(I4B) :: nn,ix,iy,iz
    REAL(DP) :: b0,b1
    REAL(DP),ALLOCATABLE :: rhsf(:,:,:)
    !>>>>>>>>>>>>>>>>>>>>>>>>>BEGIN LAPLB_CFOUR>>>>>>>>>>>>>>>>>>>>>>>>
    nn=SIZE(nrhs,1)
    ALLOCATE(rhsf(0:nn+1,0:nn+1,0:nn+1))
    rhsf(:,:,:)=0.d0
    rhsf(1:nn,1:nn,1:nn)=rhs(:,:,:)
    b0=0.5d0
    b1=1.d0/12.d0
    DO iz=1,nn
       DO iy=1,nn
          DO ix=1,nn
             nrhs(ix,iy,iz)=b0*rhsf(ix,iy,iz) + &
                  &b1*(rhsf(ix-1,iy,iz)+rhsf(ix+1,iy,iz)+rhsf(ix,iy-1,iz)+rhsf(ix,iy+1,iz)+rhsf(ix,iy,iz-1)+rhsf(ix,iy,iz+1))
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(rhsf)
    !<<<<<<<<<<<<<<<<<<<<<<<<<END   LAPLB_CFOUR<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE laplb_cfour
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE calclm(Lmax,clm)
    !USE Math, ONLY : c
    INTEGER(I4B) :: Lmax
    REAL(DP)     :: clm(Lmax,0:Lmax)
    INTEGER(I4B) :: l,m
    DO l=1,Lmax
       DO m=0,l
          clm(l,m)=c(l,m)
       ENDDO
    ENDDO
  ENDSUBROUTINE calclm
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  FUNCTION c(l,m)
    INTEGER(I4B):: l,m,i
    REAL(DP)    :: c
    if(m==0) then
       c=1.0
       return
    else
       c=1.0
       do i=l-m+1,l+m
          c=i*c
       end do
       c=1.0d0/c
    end if
    return
  ENDFUNCTION c
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE rcs(n1,n2,n3,dr,cost,sint,cns,indx)
    !USE struct_info, ONLY: Rermax
    USE struct_module
    INTEGER(I4B),INTENT(IN) :: n1,n2,n3
    REAL(DP),DIMENSION(n1,n2,n3),INTENT(OUT)    :: dr,cost,sint
    COMPLEX(DCP),DIMENSION(n1,n2,n3),INTENT(OUT)  :: cns
    INTEGER(I4B),DIMENSION(n1,n2,n3),INTENT(OUT)  :: indx
    INTEGER(I4B) :: ix,iy,iz
    REAL(DP) :: orig(3),cosp,sinp
    !>>>>>>>>>>>>>>>>>>>>>>>>>START F_DR>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !cubic cell?
    IF(n1/=n2.OR.n2/=n3.OR.n1/=n3) THEN
       WRITE(6,*) 'CGPHI:STOP!!! need cubic cell'
       STOP
    ENDIF
    orig(:)=0.5d0*gap(1)*(n1+2)

    DO iz=1,n3
       DO iy=1,n2
          DO ix=1,n1
             CALL car2spe(ORIG,ix,iy,iz,dr(ix,iy,iz),cost(ix,iy,iz),sint(ix,iy,iz),cosp,sinp)
             cns(ix,iy,iz)=CMPLX(cosp,sinp)

             IF(dr(ix,iy,iz) <= RadiusMax)THEN
                indx(ix,iy,iz) = 1
             ELSE
                indx(ix,iy,iz) = 0
             ENDIF

          ENDDO
       ENDDO
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<< END  F_DR<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE rcs
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE Vhartree_FMM(rhoS,VCoulomb)
    !===============================================
    !##CALCULATE VHARTREE VIA FMM CALLED FROM LFMM3D
    !===============================================
    IMPLICIT NONE
    REAL(DP),INTENT(INOUT) :: VCoulomb(:,:,:) !Vhartree in grid
    REAL(DP),INTENT(IN)    :: rhoS(:,:,:)
    !=======================================================================
    !##LOCAL
    INTEGER(I4B)  :: i,j,k,m
    !...FOR CALL LFMM3D...
    !###INPUT:
    INTEGER(I4B)  :: iprec=1
    !## precision flag. Allowed values are
    !##     iprec = −2 for least squares errors < 0.5 10 0 ,
    !##		 iprec = −1 for least squares errors < 0.5 10 −1 ,
    !##		 iprec = 0 for least squares errors < 0.5 10 −2 ,
    !##		 iprec = 1 for least squares errors < 0.5 10 −3 .
    !##		 iprec = 2 for least squares errors < 0.5 10 −6 .
    !##		 iprec = 3 for least squares errors < 0.5 10 −9 .
    !##		 iprec = 4 for least squares errors < 0.5 10 −12 .
    !##		 iprec = 5 for least squares errors < 0.5 10 −14 .
    INTEGER(I4B)  :: nsource
    !## number of sources
    REAL(DP)      :: source(3,rho_calc%OneDLength)
    ! REAL(DP)      :: source(3,n)
    !## sources(k,j) is the kth component of the jth source in R 3 .
    INTEGER(I4B)  :: ifcharge=1
    !## charge flag. If icharge = 1, then include the effect of the
    !## charge sources. Otherwise,omit.
    COMPLEX(DP)   :: charge(rho_calc%OneDLength)
    ! COMPLEX(DP)   :: charge(n)
    !## charge(j) is the strength of the jth charge
    INTEGER(I4B)  :: ifdipole=0
    !## dipole flag. If idipole = 1, then include the effect of the
    !## dipole sources. Otherwise,omit.
    COMPLEX(DP)   :: dipstr!(rho_calc%OneDLength)
    !## dipstr(j) is the strength of the jth dipole (p j in the formula (1)).
    REAL(DP)      :: dipvec(3)!,rho_calc%OneDLength)
    !## dipvec(k,j) is the kth component of the orientation vector of the
    !## jth dipole
    INTEGER(I4B)  :: ifpot=1
    !## potential flag. If ifpot = 1, the potential is computed. Otherwise, it is not.
    INTEGER(I4B)  :: iffld=0
    !## field (gradient) flag. If iffld = 1 the gradient of the potential is
    !## computed. Otherwise, it is not.
    !###OUPUT:
    INTEGER(I4B)  :: ier
    !## Error return codes.
    !##     ier = 0: Successful completion of code.
    !##     ier = 4: failure to allocate memory for oct-tree
    !##     ier = 8: failure to allocate memory for FMM workspaces
    !##     ier = 16: failure to allocate meory for multipole/local expansions
    COMPLEX(DP)   :: pot(rho_calc%OneDLength)
    ! COMPLEX(DP)   :: pot(n)
    !COMPLEX(DP)   :: pot
    !## pot(i) is the potential at the ith source
    COMPLEX(DP)   :: fld(3)!,rho_calc%OneDLength)
    !## fld(k,i) is the kth component of the field (-gradient of the potential)
    !## at the ith source
    !##TEST BY lfmm3dparttarg
    INTEGER(I4B)  :: ntarget,ifpottarg=0,iffldtarg=0
    REAL(DP)      :: target_fmm!(3,rho_calc%OneDLength)
    COMPLEX(DP)   :: pottarg!(rho_calc%OneDLength),
    COMPLEX(DP)   :: fldtarg(3)
    !===============================================
    ! iprec=iprec_fmm
    ! !## ARRAY ASSIGNMENT
    ! nsource=rho_calc%OneDLength
    ! ! nsource=n
    ! source(1,:)=rho_calc%x(:)*gap(1)-gap(1)
    ! source(2,:)=rho_calc%y(:)*gap(2)-gap(2)
    ! source(3,:)=rho_calc%z(:)*gap(3)-gap(3)
    ! ! source(:,:)=grid%rvec(1:3,:)!-gap(1)
    ! !source=0.d0
    ! DO i=1,rho_calc%OneDLength
    !    ! charge(i)=cmplx(rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i)),rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i)))*dvol
    !    charge(i)=cmplx(rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i)),0.d0)*dvol
    !    !##a failier test for one point
    !    !charge=cmplx(1.d0)
    ! ENDDO
    ! ! charge(:)=cmplx(reshape(rhoS*dvol,(/n/)))
    ! ! m=0
    ! ! do k=1,n3,1
    ! ! do j=1,n2,1
    ! ! do i=1,n1,1
    ! !    m=m+1
    ! !    charge(m)=cmplx(rhoS(i,j,k)*dvol,0.d0)
    ! ! enddo
    ! ! enddo
    ! ! enddo
    ! ! print*,"what",sum(rhoS*dvol),"m",m,"n",n
    ! !================================
    ! !##OUTPUT CHARGE
    ! !open(1111,file="charge")
    ! !write(1111,*)charge
    ! !close(1111)
    ! !================================
    ! dipstr=cmplx(0.d0,0.d0)
    ! dipvec=0.d0
    ! !=========================
    ! !## lfmm3dparttarg

    ! ntarget=0!rho_calc%OneDLength
    ! target_fmm=0.d0
    ! !=========================
    ! !## lfmm3dpartself(ier,iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)
    ! !============
    ! !##FOR DIRECTLY GET Vhart
    ! !goto 5073
    ! !============
    ! !CALL lfmm3dpartself(ier, iprec, nsource, source, ifcharge, charge, &
    ! !                   & ifdipole, dipstr, dipvec,ifpot, pot, iffld, fld)
    ! !##a failier test for one point
    ! CALL lfmm3dparttarg(ier, iprec, nsource, source, ifcharge, charge, &
    !      & ifdipole, dipstr, dipvec,ifpot, pot, iffld, fld,&
    !      & ntarget,target_fmm,ifpottarg,pottarg,iffldtarg,fldtarg)
    ! ! CALL l3dpartdirect(nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, &
    ! ! &     iffld, fld, ntarget, target_fmm, ifpottarg, pottarg, iffldtarg, fldtarg)
    ! !=======================================================================
    ! !##FOR JUDGE OUTPUT
    ! SELECT CASE(ier)
    ! CASE(0)
    !    CONTINUE
    ! CASE(4)
    !    print*,"Vhartree_FMM:failure to allocate memory for oct-tree"
    !    stop
    ! CASE(8)
    !    print*,"Vhartree_FMM:failure to allocate memory for oct-tree"
    !    stop
    ! CASE(16)
    !    print*,"Vhartree_FMM:failure to allocate memory for oct-tree"
    !    stop
    ! CASE DEFAULT
    !    print*,"ERROR IN Vhartree_FMM"
    !    stop
    ! END SELECT
    ! !## FORMAT TRANSLATION
    ! !=======
    ! !5073 !FOR GET Vhart DIRECTLY
    ! !======
    ! VCoulomb=0.d0
    ! !================================================
    ! !##READ Vhart
    ! !open(1111,file="vhart")
    ! !read(1111,*)pot
    ! !close(1111)
    ! !================================================
    ! DO i=1,rho_calc%OneDLength
    !    VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))=REAL(pot(i),DP)
    !    !##a failier test for one point
    !    !VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))=pottarg
    ! ENDDO
    ! ! VCoulomb=reshape(REAL(pot(:),DP),(/n1,n2,n3/))
    ! print*,"sum pot",sum(pot)
    ! !##a failier test for one point
    ! !print*,"sum pot",pottarg
    ! !========================================
    ! !open(unit=2223, file="VCoulomb_1", status="UNKNOWN", action="write")
    ! !write(unit=2223, fmt=*) VCoulomb
    ! !print*,"VCoulomb writed"
    ! !close(2223)
    ! !========================================
    ! !open(123,file="fr_fmm")
    ! !do i=1,n1
    ! !  write(123,*)VCoulomb(i,i,i)
    ! !enddo
    ! !close(123)
    ! !open(123,file="r_fmm")
    ! !do i=1,n1
    ! !  write(123,*)(3*(i-1)**2*gap(1)**2)**0.5d0
    ! !enddo
    ! !close(123)
    ! !========================================
  END SUBROUTINE Vhartree_FMM
  !-----------------------------PARTING-LINE-------------------------------
  SUBROUTINE Vhartree_direct(rhoS,VCoulomb)
    !===============================================
    !##CALCULATE VHARTREE VIA FMM CALLED FROM LFMM3D
    !===============================================
    USE Math, ONLY: NORM
    IMPLICIT NONE
    REAL(DP),INTENT(INOUT) :: VCoulomb(:,:,:) !Vhartree in grid
    REAL(DP),INTENT(IN)    :: rhoS(:,:,:)
    ! INTEGER(I4B)           :: i1,i2,i3,j1,j2,j3
    INTEGER(I4B)           :: i,j
    INTEGER(I4B)           :: na,nb
    ! INTEGER(I4B)           :: x1,y1,z1
    ! INTEGER(I4B)          /:: x2,y2,z2
    ! REAL(DP)               :: dx2,dy2,dz2
    REAL(DP)               :: coeff1
    REAL(DP)               :: coeff2
    REAL(DP)           :: VCoulomb_correct(n1,n2,n3)

    CALL cal_Vsrcpot(size(rhoS),rhoS,grid%rVec(1:3,:),size(VCoulomb),VCoulomb,grid%rVec(1:3,:))
    RETURN
    ! coeff1=dvol/gap(1)
    ! VCoulomb=0.d0
    ! do i=1,rho_calc%OneDLength,1
    !    x1=rho_calc%x(i)
    !    y1=rho_calc%y(i)
    !    z1=rho_calc%z(i)
    !    do j=1,rho_calc%OneDLength,1
    !       IF(i==j)cycle
    !       x2=rho_calc%x(j)
    !       y2=rho_calc%y(j)
    !       z2=rho_calc%z(j)
    !       dx2=real(x2-x1,DP)**2
    !       dy2=real(y2-y1,DP)**2
    !       dz2=real(z2-z1,DP)**2
    !       VCoulomb(x1,y1,z1)=VCoulomb(x1,y1,z1) + rhoS(x2,y2,z2)/dsqrt(dx2+dy2+dz2)*coeff1
    !    enddo
    ! enddo
    coeff1=dvol/gap(1)
    ! coeff2=0.2d0*(gap(1)**2*(36.d0*pi**2)**(1.d0/3))
    ! VCoulomb_correct(:,:,:)=pi/3.d0*rhoS(:,:,:)
    VCoulomb=0.d0
    ! print *, "HAhAhaAH" ,sum(pi/3.d0*rhoS)*dvol
    ! print *, "HAhA  rhos" ,sum(rhoS)*dvol
    ! print *, "HAhAhaAH" ,pi/3.d0
    do i=1,rho_calc%OneDLength,1
       ! na=rho_calc%x(i)+(rho_calc%y(i)-1)*n1+(rho_calc%z(i)-1)*n1*n2
       do j=1,rho_calc%OneDLength,1
          IF(i==j)THEN
             cycle
             ! VCoulomb(/ho_sphere%x(i),rho_calc%y(i),rho_calc%z(i))=&
             !      & VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))&
             !      & +rhoS(rho_calc%x(j),rho_calc%y(j),rho_calc%z(j))*pi/3.d0/coeff1
          ELSE
             VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))=&
                  & VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))&
                  & +rhoS(rho_calc%x(i),rho_calc%y(j),rho_calc%z(j))&
                  &/dsqrt(real((rho_calc%x(j)-rho_calc%x(i))**2&
                  &           +(rho_calc%y(j)-rho_calc%y(i))**2&
                  &           +(rho_calc%z(j)-rho_calc%z(i))**2,DP))
             !    nb=rho_calc%x(j)+(rho_calc%y(j)-1)*n1+(rho_calc%z(j)-1)*n1*n2
             ! VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))=&
             !      & VCoulomb(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))&
             !      & + rhoS(rho_calc%x(j),rho_calc%y(j),rho_calc%z(j))&
             !      & /norm(grid%rVec(1:3,na)-grid%rVec(1:3,nb))
          ENDIF
       enddo
    enddo
    ! do i3=1,n3,1
    ! do i2=1,n2,1
    ! do i1=1,n1,1
    !    do j3=1,n3,1
    !    do j2=1,n2,1
    !    do j1=1,n1,1
    !       IF(i1==j1.and.i2==j2.and.i3==j3)cycle
    !       VCoulomb(i1,i2,i3)=&
    !            & VCoulomb(i1,i2,i3)&
    !            & + rhoS(j1,j2,j3)&
    !            &/dsqrt(real((j1-i1)**2&
    !            &           +(j2-i2)**2&
    !            &           +(j3-i3)**2,DP))!*coeff1
    !    enddo
    !    enddo
    !    enddo
    ! enddo
    ! enddo
    ! enddo
    VCoulomb=VCoulomb*coeff1 !+VCoulomb_correct
  ENDSUBROUTINE Vhartree_direct
  !=======================================================================
  SUBROUTINE cal_Vsrcpot(size_src,src,pos_src,size_tar,pot_tar,pos_tar)
    USE math,ONLY: NORM
    !> src: the source charge,transport a array(n1,n2,n3) is support
    !> when size_src is the product of n1,n2,n3
    !> pos_src: source charge position
    !> pot_tar : the potential calculated
    !> pos_tar: target charge position
    !>> dvol -- grid_module
    IMPLICIT NONE
    !> in/out
    INTEGER(I4B) :: size_src,size_tar
    REAL(DP) :: src(size_src),pos_src(3,size_src)
    REAL(DP) :: pot_tar(size_tar),pos_tar(3,size_tar)
    !> local variable
    INTEGER(I4B) :: i_src,j_tar,i_ignore
    REAL(DP)     :: r_module

    pot_tar(:)=0.d0
    i_ignore=0

    do j_tar=1,size_tar,1
       do i_src=1,size_src,1
          r_module=NORM(pos_src(:,i_src)-pos_tar(:,j_tar))
          if(r_module <= 0.1d-12)then
             i_ignore=i_ignore+1
             cycle
          endif
          pot_tar(j_tar)=pot_tar(j_tar)+src(i_src)/r_module*dvol
       enddo
    enddo

    if(i_ignore/=0.and.i_ignore/=size_src)then
       print*,"direct product the electron interaction, ignore [",i_ignore,"] grids"
    endif

  ENDSUBROUTINE cal_Vsrcpot
  SUBROUTINE vhart_iso(rho,vhart)
    USE grid_module , ONLY : ng1,ng2,ng3,n1,n2,n3,grid
    USE FOURIER
    IMPLICIT NONE
    !
    REAL(DP),INTENT(IN) :: rho(:,:,:)
    REAL(DP),INTENT(OUT) :: vhart(:,:,:)
    !
    COMPLEX(DP),DIMENSION(ng1,ng2,ng3) :: rho_recip,vh_recip
    INTEGER(I4B) :: Ig,I1,I2,I3
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    rho_recip(:,:,:)=FFT(rho)
    !q=0
    grid%gVec(4,1)=1.d0
    !v(q)=4*pi*rho(q)/(q^2)
    Ig=0
    DO I3=1,ng3
       DO I2=1,ng2
          DO I1=1,ng1
             Ig=Ig+1
             vh_recip(I1,I2,I3)=(4.d0*pi)*rho_recip(I1,I2,I3)/(grid%gVec(4,Ig))**2
          ENDDO
       ENDDO
    ENDDO
    !
    grid%gVec(4,1)=0.d0
    !remove q=0 because of netural system
    vh_recip(1,1,1)=(0.d0,0.d0)
    !
    vhart(:,:,:)=FFT(vh_recip)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE vhart_iso
  SUBROUTINE poisson_MCM(rho,pot)
    !> Multipole corrections method
    !> Can. J. Phys. 81: 1151–1164 (2003) doi: 10.1139/P03-078
    USE parameters ,ONLY: Lcellbohr
    USE struct_module, ONLY: ncharge
    USE array_io,ONLY: output
    IMPLICIT NONE
    REAL(DP) :: rho(:,:,:)
    REAL(DP) :: pot(:,:,:)
    !> local
    REAL(DP)              :: Q_l0(0:Lmax)
    COMPLEX(DCP)          :: Q_lm(Lmax,Lmax)
    COMPLEX(DCP),DIMENSION(ng1,ng2,ng3) :: rho_recip,vh_recip
    REAL(DP) :: Vhart_fft(n1,n2,n3),V_corr(n1,n2,n3),V_corr2(n1,n2,n3),V0
    INTEGER(I4B)  :: x1,x2,x3
    REAL(DP) :: alpha
    REAL(DP) :: rhoaux(n)
    REAL(DP) :: coeff_charge
    !> time evaluating flag
    LOGICAL  :: t_p=.false. !time of MCM poisson solver

    !> calculate Qlm
    CALL time_start('Qlm_calculate_consuming',t_p)
    CALL cal_qlm(rho,Lmax,Q_l0,Q_lm)
    CALL time_end('Qlm_calculate_consuming',t_p)
    CALL time_output('Qlm_calculate_consuming',t_p)
    !> calculate auxiliary charge density

    !> n_aux(r)=\sum_l{\sum_m{n_lm(r)}}
    !> n_lm(r)=Q_lm*{ 2**(l+2) / (a**(2l+3)*sqrt(pi)*(2l+1)!!) * r**l * exp(-(r/a)**2) * Y_lm(r)}
    !> V_H=V_Hfft-V_corr+V0
    !> Vcorr=\sum_lm{Q_lm*\psi_lm}

    ! CALL time_start('rho-aux_calculate_consuming',t_p)
    ! CALL MCM_cal_naux(Q_l0,Q_lm,rhoaux)
    ! CALL time_end('rho-aux_calculate_consuming',t_p)
    ! CALL time_output('rho-aux_calculate_consuming',t_p)

    CALL time_start('FFT_calculate_consuming',t_p)
    CALL Vhart_iso(rho,Vhart_fft)
    CALL time_end('FFT_calculate_consuming',t_p)
    CALL time_output('FFT_calculate_consuming',t_p)

    ! CALL MCM_calVcorr(Q_l0,Q_lm,V_corr)
    ! print*,"fft_vcorr1",0.5d0*sum(V_corr*rho)*dvol*hart2eV
    ! coeff_charge=ncharge/(sum(rhoaux)*dvol)
    ! rhoaux=rhoaux*ncharge/(sum(rhoaux)*dvol)
    ! print*,"ssss",sum(rhoaux)*dvol
    ! stop
    ! CALL Vhart_iso(reshape(rhoaux,(/n1,n2,n3/)),V_corr)
    ! CALL GMG_CG_hart(reshape(rhoaux,(/n1,n2,n3/)),V_corr2)
    CALL MCM_calVcorr(Q_l0,Q_lm,V_corr2)
    ! pot=V_corr2
    ! V_corr=V_corr-V_corr2
    ! print*,0.5d0*sum(Vhart_fft*rho)*dvol*hart2eV
    ! print*,"fft_vcorr",0.5d0*sum(V_corr*rho)*dvol*hart2eV
    ! print*,"fft_vcorr2",0.5d0*sum(V_corr2*rho)*dvol*hart2eV*coeff_charge
    print*,"fft_vcorr2-no-normalize",0.5d0*sum(V_corr2*rho)*dvol*hart2eV
    print*,"fft_v",0.5d0*sum(Vhart_fft*rho)*dvol*hart2eV
    pot=Vhart_fft-V_corr2
    ! pot=Vhart_fft-V_corr+V_corr2*coeff_charge
    ! CALL output(n,reshape(V_corr2*coeff_charge,(/n/)),"Vcorr2")
    CALL output(n,reshape(V_corr2,(/n/)),"Vcorr2")
    V0=0.d0
    alpha=1.586718d0
    ! V0=-(sum(pot(:,:,1:2))+sum(pot(:,:,n3-1:n3))+sum(pot(:,1:2,3:n3-2))+sum(pot(:,n2-1:n2,3:n3-2))+sum(pot(1:2,3:n2-2,3:n3-2))+&
    ! & sum(pot(n1-1:n1,3:n2-2,3:n3-2)))/real((n1**2*4+(n1-4)*n1*4+(n1-4)**2*4),DP)+alpha*ncharge/Lcellbohr
    V0=-(sum(pot(:,:,1))+sum(pot(:,:,n3))+sum(pot(:,1,2:n3-1))+sum(pot(:,n2,2:n3-1))+sum(pot(1,2:n2-1,2:n3-1))+&
         & sum(pot(n1,2:n2-1,2:n3-1)))/real((n1**2*2+(n1-2)*n1*2+(n1-2)**2*2),DP)+alpha*ncharge/Lcellbohr
    print*,"V0_corr",0.5d0*sum(V0*rho)*dvol*hart2eV
    pot=pot+V0
  ENDSUBROUTINE poisson_MCM
  SUBROUTINE MCM_calVcorr(Q_l0,Q_lm,V_corr)
    USE struct_module, ONLY: volume
    USE array_io,ONLY: output
    USE FOURIER
    USE math, ONLY:  gammp,integral
    IMPLICIT NONE
    !> Vcorr=\sum_lm{Q_lm*\psi_lm}
    !> \psi_lm=4*pi**2/\Omega*( i**l/(2l+1)!! )*FFT{G**(l-2)*exp(-a**2*G**2/4)*Y_lm(G)}-sqrt(pi)*2**(l+3)/(2*l+1)!!/r**(l+1)*I_l(r/a)*Y_lm(r)
    !> in/out
    REAL(DP)  :: Q_l0(0:),V_corr(:,:,:)
    COMPLEX(DCP) :: Q_lm(:,:)
    !> local
    COMPLEX(DCP)  :: rhoaux_fft(ng1,ng2,ng3,0:Lmax,0:Lmax)
    REAL(DP)  :: rhoaux_anyl(n1,n2,n3,0:Lmax,0:Lmax)
    INTEGER(I4B) :: ig1,ig2,ig3,Ig,in1,in2,in3,in
    COMPLEX(DCP)  :: coeff1(0:lmax)
    REAL(DP)      :: coeff2(0:Lmax)
    INTEGER(I4B)  :: l_factorial
    REAL(DP),DIMENSION(ng1,ng2,ng3)        :: dg,gcost,gsint
    COMPLEX(DCP),DIMENSION(ng1,ng2,ng3)    :: gcns
    ! INTEGER(I4B),DIMENSION(n1,n2,n3)    :: gindx
    COMPLEX(DCP)  :: ylm(0:Lmax,0:Lmax)
    REAL(DP),DIMENSION(n1,n2,n3,0:Lmax,0:Lmax)        :: phi1,phi2
    INTEGER(I4B) :: il,im
    REAL(DP)     :: gauss_a
    REAL(DP),parameter :: four_pi=4.d0*pi
    !> flag of time evaluating
    LOGICAL :: t_pV=.False.  !> t_poisson_Vcorr
    LOGICAL :: t_pV1=.true.  !> t_poisson_Vcorr1
    LOGICAL :: t_cycle=.false.  !> t_flag for big cycle

    ! gauss_a=4.124d0*gap(1)
    gauss_a=0.8448
    print*,"gauss_a",gauss_a
    Ig=0
    !>
    CALL time_start('Vcorr1_coeff_calculate_consuming',t_pV1)
    l_factorial=1.d0
    do il=0,Lmax,1
       l_factorial=l_factorial*(2*il+1)
       coeff1(il)=(4.d0*pi)**2/volume*cmplx(0.d0,1.d0)**il/l_factorial*(2.d0*il+1)/four_pi
       ! coeff1(il)=(4.d0*pi)**2/l_factorial
       ! coeff1(il)=cmplx(0.d0,1.d0)**il/l_factorial
       ! print*,coeff1(il)
    enddo
    CALL time_end('Vcorr1_coeff_calculate_consuming',t_pV1)
    CALL time_output('Vcorr1_coeff_calculate_consuming',t_pV1)
    ! print*,coeff1
    !> calcoordinate in sphere
    CALL time_start('Vcorr1_recipgrid_calculate_consuming',t_pV1)
    CALL rcs_rec(ng1,ng2,ng3,gcost,gsint,gcns)
    CALL time_end('Vcorr1_recipgrid_calculate_consuming',t_pV1)
    CALL time_output('Vcorr1_recipgrid_calculate_consuming',t_pV1)
    !CALL output(n,reshape(dg,(/size(dg)/)),'dg')
    !CALL output(n,grid%gvec(4,:),'grid%gvec')
    !> rhoaux_recip
    CALL time_start('Vcorr1-recip-potential_calculate_consuming',t_pV1)
    do ig3=1,ng3,1
       do ig2=1,ng2,1
          do ig1=1,ng1,1
             Ig=Ig+1
             if(Ig1==1.and.ig2==1.and.ig3==1)cycle
             t_cycle=.false.
             ! if(Ig==5)t_cycle=.true.
             CALL time_start('Vcorr1_calylm',t_cycle)
             CALL cal_ylm(gcost(ig1,ig2,ig3),gsint(ig1,ig2,ig3),gcns(ig1,ig2,ig3),ylm)
             CALL time_end('Vcorr1_calylm',t_cycle)
             CALL time_output('Vcorr1_calylm',t_cycle)
             CALL time_start('Vcorr1_calp',t_cycle)
             do il=0,Lmax,1
                do im=1,il,1
                   rhoaux_fft(ig1,ig2,ig3,il,im)=grid%gVec(4,Ig)**(il-2)*exp(-gauss_a**2*grid%gVec(4,Ig)**2/4.d0)*Ylm(il,im)*coeff1(il)&
                        & *exp(cmplx(0.d0,1.d0)*sum(grid%gVec(1:3,Ig)*(0.5d0*(n1)*gap(1:3))))
                   ! print*,rhoaux_fft(ig1,ig2,ig3,il,im),
                enddo
                rhoaux_fft(ig1,ig2,ig3,il,0)=grid%gVec(4,Ig)**(il-2)*exp(-gauss_a**2*grid%gVec(4,Ig)**2/4.d0)*Ylm(il,0)*coeff1(il)&
                     & *exp(cmplx(0.d0,1.d0)*sum(grid%gVec(1:3,Ig)*(0.5d0*(n1)*gap(1:3))))
             enddo
             CALL time_end('Vcorr1_calp',t_cycle)
             CALL time_output('Vcorr1_calp',t_cycle)
          enddo
       enddo
    enddo
    CALL time_end('Vcorr1-recip-potential_calculate_consuming',t_pV1)
    CALL time_output('Vcorr1-recip-potential_calculate_consuming',t_pV1)
    rhoaux_fft(1,1,1,:,:)=(0.d0,0.d0)
    !> FFT(rhoaux_fft)
    CALL time_start('Vcorr1-real-potential_calculate_consuming',t_pV1)
    do il=0,Lmax,1
       do im=0,il,1
          phi1(:,:,:,il,im)=FFT(rhoaux_fft(:,:,:,il,im))
       enddo
    enddo
    CALL time_end('Vcorr1-real-potential_calculate_consuming',t_pV1)
    CALL time_output('Vcorr1-real-potential_calculate_consuming',t_pV1)
    !> ==============================================
    !> cal phi2
    !> sqrt(pi)*2**(l+3)/(2*l+1)!!/r**(l+1)*I_l(r/a)*Y_lm(r)
    !> ==============================================
    CALL time_start('Vcorr2_calculate_consuming',t_pV)
    l_factorial=1.d0
    do il=0,Lmax,1
       l_factorial=l_factorial*(2*il+1)
       coeff2(il)=sqrt(pi)*2**(il+3)/real(l_factorial)*(2*il+1)/(four_pi)
    enddo
    in=0
    do in3=1,n3,1
       do in2=1,n2,1
          do in1=1,n1,1
             CALL cal_ylm(cost(in1,in2,in3),sint(in1,in2,in3),cns(in1,in2,in3),ylm)
             ! print*,ylm
             ! print*,"-------"
             in=in+1
             do il=0,Lmax,1
                do im=0,il,1
                   rhoaux_anyl(in1,in2,in3,il,im)=1.d0/dr(in1,in2,in3)**(il+1)*Ylm(il,im)*coeff2(il)*0.5d0*gamma(real(il,DP)+0.5d0)*gammp(real(il,DP)+0.5d0,dr(in1,in2,in3)/gauss_a)
                   ! rhoaux_anyl(in1,in2,in3,il,im)=1.d0/dr(in1,in2,in3)**(il+1)*Ylm(il,im)*coeff2(il)*integrade_i(il,dr(in1,in2,in3)/gauss_a)
                   ! rhoaux_anyl(in1,in2,in3,il,im)=1.d0/dr(in1,in2,in3)**(il+1)*Ylm(il,im)*coeff2(il)*integral(il,0.001d0,dr(in1,in2,in3)/gauss_a)
                   ! print*,integrade_i(il,dr(in1,in2,in3)/gauss_a)
                enddo
             enddo
          enddo
       enddo
    enddo
    ! print*,"gauss_a",gauss_a
    ! print*,"ddddddd",gammp(0.5d0,5.d0)
    ! CALL output(n,reshape(V_corr,(/n/)),'V_corr')
    CALL time_end('Vcorr2_calculate_consuming',t_pV)
    CALL time_output('Vcorr2_calculate_consuming',t_pV)

    CALL time_start('Vcorr3_calculate_consuming',t_pV)
    V_corr=0.d0
    do il=1,Lmax,1
       do im=1,il,1
          V_corr=V_corr+2.d0*real(Q_lm(il,im)*(phi1(:,:,:,il,im)-rhoaux_anyl(:,:,:,il,im)),DP)
          ! V_corr=V_corr+2.d0*real(Q_lm(il,im)*(phi1(:,:,:,il,im)),DP)!-rhoaux_anyl(:,:,:,il,im)),DP)
          ! V_corr=V_corr+2.d0*real(Q_lm(il,im)*(rhoaux_anyl(:,:,:,il,im)),DP)
       enddo
    enddo
    do il=0,Lmax,1
       V_corr=V_corr+real(Q_l0(il)*(phi1(:,:,:,il,0)-rhoaux_anyl(:,:,:,il,0)),DP)
       ! V_corr=V_corr+real(Q_l0(il)*(phi1(:,:,:,il,0)),DP)!-rhoaux_anyl(:,:,:,il,0)),DP)
       ! V_corr=V_corr+real(Q_l0(il)*(rhoaux_anyl(:,:,:,il,0)),DP)
    enddo
    CALL time_end('Vcorr3_calculate_consuming',t_pV)
    CALL time_output('Vcorr3_calculate_consuming',t_pV)
    ! print *,'Ql0',Q_l0
    ! print *,'Qlm',Q_lm
    ! CALL output(n,reshape(V_corr,(/n/)),'V_corr')
    ! CALL output(n,reshape(V_corr,(/n/)),'V_anyl')
  ENDSUBROUTINE MCM_calVcorr

  SUBROUTINE MCM_cal_naux(Q_l0,Q_lm,rhoaux)
    !> n_lm(r)=Q_lm*{ 2**(l+2) / (a**(2l+3)*sqrt(pi)*(2l+1)!!) * r**l * exp(-(r/a)**2) * Y_lm(r)}
    USE array_io , only: output
    IMPLICIT NONE
    !> in/out
    !   INTEGER(I4B) :: l,m
    REAL(DP)     :: Q_l0(0:)
    ! REAL(DP),INTENT(IN)     :: Q_l0(0:)
    COMPLEX(DCP) :: Q_lm(:,:)
    ! COMPLEX(DCP),INTENT(IN) :: Q_lm(:,:)
    REAL(DP)            :: rhoaux(n)
    !   REAL(DP)     :: n_lm(n)
    !   !> local
    !   INTEGER(I4B) :: ix,iy,iz,in
    REAL(DP)     :: a
    REAL(DP) :: l_factorial,dual_factorial(0:Lmax)
    COMPLEX(DCP) :: ylm(0:Lmax,0:Lmax)
    REAL(DP) :: plm(Lmax,0:Lmax)
    !REAL(DP) :: n_lm(n)
    INTEGER(I4B) :: in,in1,in2,in3,il,im

    a=3.9d0*gap(1)
    l_factorial=1.d0
    do il=0,Lmax,1
       l_factorial=l_factorial*(2*il+1)
       dual_factorial(il)=real(l_factorial)
    enddo
    in=0
    rhoaux=0.d0
    ! CALL output(n,reshape(cost,(/n/)),"cos")
    ! CALL output(n,reshape(sint,(/n/)),"sin")
    do in3=1,n3,1
       do in2=1,n2,1
          do in1=1,n1,1
             CALL cal_ylm(cost(in1,in2,in3),sint(in1,in2,in3),cns(in1,in2,in3),ylm)
             ! CALL cal_plm(Lmax,cost(in1,in2,in3),sint(in1,in2,in3),plm)
             ! print*,ylm
             ! print*,"-------"
             in=in+1
             do il=1,Lmax,1
                do im=1,il,1
                   IF(indx(in1,in2,in3)==0) CYCLE
                   ! rhoaux_anyl(in1,in2,in3,il,im)=1.d0/dr(in1,in2,in3)**(il+1)*Ylm(il,im)*coeff2(il)*0.5d0*gamma(real(il,DP)+0.5d0)*gammp(real(il,DP)+0.5d0,dr(in1,in2,in3)/gauss_a)
                   rhoaux(in)=rhoaux(in)+real( Q_lm(il,im)*conjg(Ylm(il,im)),DP)*2.d0*2**(il+2) / (a**(2*il+3)*sqrt(pi)*dual_factorial(il)) * dr(in1,in2,in3)**il * exp(-(dr(in1,in2,in3)/a)**2)*(2*il+1)/(4.d0*pi)
                   ! rhoaux(in)=rhoaux(in)+real( Q_lm(il,im) *(Ylm(il,im)),DP)*2.d0 * dr(in1,in2,in3)**il*exp(-(dr(in1,in2,in3)/a)**2)/ (a*sqrt(pi))
                   ! if(il==1.and.im==1)then !rhoaux(in)=Ylm(il,im)
                   !    rhoaux(in)=c(il,im)*real(Ylm(il,im) *conjg(Ylm(il,im)))*(2*il+1)/(4*pi)
                   ! endif
                enddo
             enddo
             do il=0,Lmax,1
                rhoaux(in)=rhoaux(in)+ real(Q_l0(il) * conjg(Ylm(il,0)),DP) / (a**(2*il+3)*sqrt(pi)*dual_factorial(il)) * dr(in1,in2,in3)**il * exp(-(dr(in1,in2,in3)/a)**2)*2**(il+2)*(2*il+1)/(4.d0*pi)
                ! rhoaux(in)=rhoaux(in)+real( Q_l0(il) *Ylm(il,0),DP) * dr(in1,in2,in3)**il  !* exp(-(dr(in1,in2,in3)/a)**2)/ (a*sqrt(pi))
             enddo
          enddo
       enddo
    enddo
    print*,"ll",n*dvol,volume,Lcellbohr
    !> 1/3*r**3=SSS(Y_lm*Y_lm^*)dxdydz
    ! print*,"ww",IsoRmaxbohr**3/3
    ! CALL output(n,reshape(rhoaux,(/n/)),"rhoaux")
  ENDSUBROUTINE MCM_cal_naux
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE cal_ylm(mycost,mysint,e_phi,ylm)
    IMPLICIT NONE
    !INTEGER(I4B) :: ix,iy,iz
    COMPLEX(DCP) :: ylm(0:Lmax,0:Lmax),e_phi
    REAL(DP)     :: plm(Lmax,0:Lmax)
    REAL(DP)     :: mycost,mysint
    !> local
    INTEGER(I4B) :: l,m

    CALL cal_plm(Lmax,mycost,mysint,plm)
    DO l=1,Lmax
       ! qrl=rho(ix,iy,iz)*hcub*dr(ix,iy,iz)**l
       !> Q_l0=\sum_{r'} q'*r'^l*P_l^0(cos\theta')
       ! Q_l0(l)=Q_l0(l)+qrl*plm(l,0)
       ylm(l,0)=plm(l,0)
       !> Q_lm=\sum_{r'} c(l,m)*q'*r'^l*P_l^m(cos\theta')*exp(-i*m*\phi)
       DO m=1,l
          ylm(l,m)=plm(l,m)*e_phi**m
       ENDDO
    ENDDO
    ylm(0,0)=CMPLX(1.d0,0.d0)
  ENDSUBROUTINE cal_ylm

  SUBROUTINE rcs_rec(ng1,ng2,ng3,cost,sint,cns)
    !USE struct_info, ONLY: Rermax
    USE struct_module
    INTEGER(I4B),INTENT(IN) :: ng1,ng2,ng3
    REAL(DP),DIMENSION(ng1,ng2,ng3),INTENT(OUT)    :: cost,sint
    COMPLEX(DCP),DIMENSION(ng1,ng2,ng3),INTENT(OUT)  :: cns
    ! INTEGER(I4B),DIMENSION(n1,n2,n3),INTENT(OUT)  :: indx
    INTEGER(I4B) :: ix,iy,iz,ig
    REAL(DP) :: orig(3),cosp,sinp
    !>>>>>>>>>>>>>>>>>>>>>>>>>START F_DR>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! !cubic cell?
    ! IF(n1*2/=n2.OR.n2/=n3.OR.n1*2/=n3) THEN
    !    WRITE(6,*) 'CGPHI:STOP!!! need cubic cell'
    !    STOP
    ! ENDIF
    orig(:)=(/0.d0,0.d0,0.d0/)!MATMUL(recip_lat,(/1,1,1/))
    ! orig(:)=MATMUL(recip_lat,(/1,1,1/))

    ig=0
    DO iz=1,ng3
       DO iy=1,ng2
          DO ix=1,ng1
             ig=ig+1
             IF(ix==1.and.iy==1.and.iz==1)cycle
             CALL car2spe_rec(ORIG,grid%gVec(1,ig),grid%gVec(2,ig),grid%gVec(3,ig),cost(ix,iy,iz),sint(ix,iy,iz),cosp,sinp)
             cns(ix,iy,iz)=CMPLX(cosp,sinp)

             ! IF(dr(ix,iy,iz) <= RadiusMax)THEN
             !    indx(ix,iy,iz) = 1
             ! ELSE
             !    indx(ix,iy,iz) = 0
             ! ENDIF

          ENDDO
       ENDDO
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<< END  F_DR<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE rcs_rec

  SUBROUTINE car2spe_rec(ORIG,x,y,z,cost,sint,cosp,sinp)
    IMPLICIT NONE
    REAL(DP),INTENT(IN) :: x,y,z
    REAL(DP),INTENT(IN) :: ORIG(3)
    REAL(DP),INTENT(OUT) :: cost,sint,cosp,sinp
    REAL(DP) ::rr,rx,ry,rz,r(3),r_out
    INTEGER :: l
    !>>>>>>>>>>>>>>>>>>>>>>>start CAR2SPE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !r(:)=MATMUL(recip_lat,(/x,y,z/))
    ! rx=gap(1)*REAL(x,8)-ORIG(1)
    ! ry=gap(2)*REAL(y,8)-ORIG(2)
    ! rz=gap(3)*REAL(z,8)-ORIG(3)
    !rx=r(1)-ORIG(1)
    !ry=r(2)-ORIG(2)
    !rz=r(3)-ORIG(3)
    !r_out=SQRT(rx**2+ry**2+rz**2)
    !rr=SQRT(rx**2+ry**2)
    r_out=sqrt(x**2+y**2+z**2)
    rr=sqrt(x**2+y**2)
    IF(r_out>0.d0)THEN
       cost=z/r_out
       sint=rr/r_out
    ELSE
       !cost=1.d0
       !sint=0.d0
       WRITE(6,*) 'Multi-pole:the origin on grid'
       STOP
    ENDIF

    IF(rr>0.d0)THEN
       cosp=x/rr
       sinp=y/rr
    ELSE
       cosp=1.d0
       sinp=0.d0
       ! WRITE(6,*) 'Multi-pole:the origin on grid'
       ! WRITE(6,*) 'r',r_out,""
       ! STOP
    ENDIF
    !<<<<<<<<<<<<<<<<<<<<<<<<<end CAR2SPE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE car2spe_rec
  FUNCTION Integrade_i(l,x)
    IMPLICIT NONE
    INTEGER(I4B)  :: l
    REAL(DP),INTENT(IN)  :: x
    REAL(DP) :: integrade_i
    !> local
    REAL(DP) :: integral
    REAL(DP)     :: sqrtPi
    REAL(DP)     :: coe0,coe2,coe4,coe6,coe8,coe10,coe12,coe14,coe16,coe18,coe20,coe22,coe24,coe26,coe28,coe30,coe32,coe34,coe36,coe38,coe40

    sqrtPi=sqrt(pi)
    select case(l)
    case(0)
       integral=0.5d0*sqrtPi*erf(x)
    case(1)
       integral=0.25d0*(-2.d0*exp(-x**2)*x+sqrtPi*erf(x))
    case(2)
       integral=-0.25d0*exp(-x**2)*x*(3.d0+2.d0*x**2)+3.d0/8.d0*sqrtPi*erf(x)
    case(3)
       integral=-0.125d0*exp(-x**2)*x*(15.d0+10.d0*x**2+4*x**4)+15.d0/16.d0*sqrtPi*erf(x)
    case(4)
       integral=1.d0/32.d0*(-2.d0*exp(-x**2)*x*(105.d0+70.d0*x**2+28*x**4+8*x**6)+105.d0*sqrtPi*erf(x))
    case(5)
       integral=1.d0/64.d0*(-2.d0*exp(-x**2)*(945.d0+630.d0*x**2+252.d0*x**4+72*x**6+16*x**8)+945.d0*sqrtPi*erf(x))
    case(6)
       integral=1.d0/128.d0*(-2.d0*exp(-x**2)*x*(10395+6930*x**2+2772*x**4+792*x**6+176*x**8+32*x**10)-10395*sqrtPi*erf(x))
    case(7)
       integral=1.d0/256.d0*(-2.d0*exp(-x**2)*x*(135135.d0+90090.d0*x**2+36036*x**4+10296*x**6+2288*x**8+416*x**10+64*x**12)+135135*sqrtPi*erf(x))
    case(8)
       integral=1.d0/512.d0*(-2*exp(-x**2)*x*(2027025+1351350*x**2+540540*x**4+154440*x**6+34320*x**8+6340*x**10+960*x**12+128*x**14)+2027025*sqrtPi*erf(x))
    case(9)
       integral=1.d0/1024.d0*(-2*exp(-x**2)*x*(34459425+22972950*x**2+9189180*x**4+2625480*x**6+583440*x**8+106080*x**10+16320*x**12+2176*x**14+256*x**16)+34459425*sqrtPi*erf(x))
    case(10)
       integral=1.d0/2048.d0*(-2*exp(-x**2)*x*(654729075+436486050*x**2+174594420*x**4+49884120*x**6+11085360*x**6+2015520*x**10+310080*x**12+41344*x**14+4864*x**16+512*x**18)&
            & +654729075*sqrtPi*erf(x))
    case(11)
       coe0 =13749310575.d0
       coe2 =9166207050.d0
       coe4 =3666482820.d0
       coe6 =1047566520.d0
       coe8 =232792560.d0
       coe10=42325920.d0
       coe12=6511680.d0
       coe14=868224.d0
       coe16=102144.d0
       coe18=10752.d0
       coe20=1024.d0
       integral=1.d0/4096.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20)&
            & +13749310575.d0*erf(x)*sqrtPi)
    case(12)
       coe0 =316234143225.d0
       coe2 =210822762150.d0
       coe4 =84329104860.d0
       coe6 =24094029960.d0
       coe8 =5354228880.d0
       coe10=973496160.d0
       coe12=149768640.d0
       coe14=19969152.d0
       coe16=2349312.d0
       coe18=247296.d0
       coe20=23552.d0
       coe22=2048.d0
       integral=1.d0/8192.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22)&
            & +316234143225.d0*erf(x)*sqrtPi)
    case(13)
       coe0 =7905853580625.d0
       coe2 =5270569053750.d0
       coe4 =2108227621500.d0
       coe6 =602350749000.d0
       coe8 =133855722000.d0
       coe10=24337404000.d0
       coe12=3744216000.d0
       coe14=499228800.d0
       coe16=58732800.d0
       coe18=6182400.d0
       coe20=588800.d0
       coe22=51200.d0
       coe24=4096.d0
       integral=1.d0/16384.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24)&
            & +7905853580625.d0*erf(x)*sqrtPi)
    case(14)
       coe0 =213458046676875.d0
       coe2 =142305364451250.d0
       coe4 =56922145780500.d0
       coe6 =16263470223000.d0
       coe8 =3614104494000.d0
       coe10=657109908000.d0
       coe12=101093832000.d0
       coe14=13479177600.d0
       coe16=1585785600.d0
       coe18=166924800.d0
       coe20=15897600.d0
       coe22=1382400.d0
       coe24=110592.d0
       coe26=8192.d0
       integral=1.d0/32768.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26)&
            & +213458046676875.d0*erf(x)*sqrtPi)
    case(15)
       coe0 =6190283353629375.d0
       coe2 =4126855569086250.d0
       coe4 =1650742227634500.d0
       coe6 =471640636467000.d0
       coe8 =104809030326000.d0
       coe10=19056187332000.d0
       coe12=2931721128000.d0
       coe14=390896150400.d0
       coe16=45987782400.d0
       coe18=4840819200.d0
       coe20=461030400.d0
       coe22=40089600.d0
       coe24=3207168.d0
       coe26=237568.d0
       coe28=16384.d0
       integral=1.d0/65536.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26+coe28*x**28)&
            & +6190283353629375.d0*erf(x)*sqrtPi)
    case(16)
       coe0 =191898783962510625.d0
       coe2 =127932522641673750.d0
       coe4 =51173009056669500.d0
       coe6 =14620859730477000.d0
       coe8 =3249079940106000.d0
       coe10=590741807292000.d0
       coe12=90883354968000.d0
       coe14=12117780662400.d0
       coe16=1425621254400.d0
       coe18=150065395200.d0
       coe20=14291942400.d0
       coe22=1242777600.d0
       coe24=99422208.d0
       coe26=7364608.d0
       coe28=507904.d0
       coe30=32768.d0
       integral=1.d0/4096.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26+coe28*x**28+coe30*x**30)&
            & +191898783962510625.d0*erf(x)*sqrtPi)
    case(17)
       coe0 =6332659870762850625.d0
       coe2 =4221773247175233750.d0
       coe4 =1688709298870093500.d0
       coe6 =482488371105741000.d0
       coe8 =107219638023498000.d0
       coe10=19494479640636000.d0
       coe12=2999150713944000.d0
       coe14=399886761859200.d0
       coe16=47045501395200.d0
       coe18=4952158041600.d0
       coe20=471634099200.d0
       coe22=41011660800.d0
       coe24=3280932864.d0
       coe26=243032064.d0
       coe28=16760832.d0
       coe30=1081344.d0
       coe32=65536.d0
       integral=1.d0/262144.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26+coe28*x**28+coe30*x**30+coe32*x**32)&
            & +6332659870762850625.d0*erf(x)*sqrtPi)
    case(18)
       coe0 =221643095476699771875.d0
       coe2 =147762063651133181250.d0
       coe4 =59104825460453272500.d0
       coe6 =16887092988700935000.d0
       coe8 =3752687330822430000.d0
       coe10=682306787422260000.d0
       coe12=104970274988040000.d0
       coe14=13996036665072000.d0
       coe16=1646592548832000.d0
       coe18=173325531456000.d0
       coe20=16507193472000.d0
       coe22=1435408128000.d0
       coe24=114832650240.d0
       coe26=8506122240.d0
       coe28=586629120.d0
       coe30=37847040.d0
       coe32=2293760.d0
       coe34=131072.d0
       integral=1.d0/524288.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26+coe28*x**28+coe30*x**30+coe32*x**32+coe34*x**34)&
            & +221643095476699771875.d0*erf(x)*sqrtPi)
    case(19)
       coe0 =8200794532637891559375.d0
       coe2 =5467196355091927706250.d0
       coe4 =2186878542036771082500.d0
       coe6 =624822440581934595000.d0
       coe8 =138849431240429910000.d0
       coe10=25245351134623620000.d0
       coe12=3883900174557480000.d0
       coe14=517853356607664000.d0
       coe16=60923924306784000.d0
       coe18=6413044663872000.d0
       coe20=610766158464000.d0
       coe22=53110100736000.d0
       coe24=4248808058880.d0
       coe26=314726522880.d0
       coe28=21705277440.d0
       coe30=1400340480.d0
       coe32=84869120.d0
       coe34=4849664.d0
       coe36=262144.d0
       integral=1.d0/1048576.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26+coe28*x**28+coe30*x**30+coe32*x**32+coe34*x**34+coe36*x**36)&
            & +8200794532637891559375.d0*erf(x)*sqrtPi)
    case(20)
       coe0 =319830986772877770815625.d0
       coe2 =213220657848585180543750.d0
       coe4 =85288263139434072217500.d0
       coe6 =24368075182695449205000.d0
       coe8 =5415127818376766490000.d0
       coe10=984568694250321180000.d0
       coe12=151472106807741720000.d0
       coe14=20196280907698896000.d0
       coe16=2376033047964576000.d0
       coe18=250108741891008000.d0
       coe20=23819880180096000.d0
       coe22=2071293928704000.d0
       coe24=165703514296320.d0
       coe26=12274334392320.d0
       coe28=846505820160.d0
       coe30=54613278720.d0
       coe32=3309895680.d0
       coe34=189136896.d0
       coe36=10223616.d0
       coe38=524288.d0
       integral=1.d0/2097152.d0*(-2.d0*exp(-x**2)*x*(coe0+coe2*x**2+coe4*x**4+coe6*x**6+coe8*x**8+coe10*x**10+coe12*x**12+coe14*x**14+coe16*x**16+coe18*x**18+coe20*x**20+coe22*x**22&
            & +coe24*x**24+coe26*x**26+coe28*x**28+coe30*x**30+coe32*x**32+coe34*x**34+coe36*x**36+coe38*x**38)&
            & +319830986772877770815625.d0*erf(x)*sqrtPi)
    case default
       print *,"Multipole order(.gt.20) is not support"
    end select
    integrade_i=integral
  END FUNCTION Integrade_i
  SUBROUTINE GMG_CG_hart(rhoS,VCoulomb)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)    :: rhoS(:,:,:)
    REAL(DP),INTENT(INOUT) :: VCoulomb(:,:,:)
    INTEGER(I4B)           :: iter=1,niter_poisson=1000,i
    ! REAL(DP)               :: VCoulomb_old(n1,n2,n3)
    INTEGER(I4B),save :: counter=0
    !
    !===========================================================
    IF(counter==0)THEN
       allocate(VCoulomb_old(n1,n2,n3))
       VCoulomb_old=0.d0
    ENDIF
    counter=counter+1
    !establish the boundary
    CALL EvaluateBoundryGrid_b(rhoS)
    !calculate the boundary potential
    CALL CalculateBoundryPotential_b(rhoS,Lmax)
    !Laplace slover
    ! CALL output(size(VCoulomb_old),reshape(VCoulomb_old,(/size(VCoulomb_old)/)),'VCoulomb')
    ! CALL output(size(VCoulomb),reshape(VCoulomb,(/size(VCoulomb)/)),'VCoulomb')
    ! BoundaryVCoulomb(CellLeft:CellRight(1),CellLeft:CellRight(2),CellLeft:CellRight(3))=VCoulomb_old
    CALL LaplaceEquationSlover_b(iter)
    !VCoulomb assignment
    CALL VCoulombAssignment_b(iter,niter_poisson,VCoulomb)
    ! VCoulomb_old=VCoulomb
    !delete the useless matrix among the subroutine
    !===================================================
    !open(123,file="fr_GMG")
    !do i=1,n1
    !	write(123,*)VCoulomb(i,i,i)
    !enddo
    !close(123)
    !open(123,file="r_GMG")
    !do i=1,n1
    !	write(123,*)(3*(i-1)**2*gap(1)**2)**0.5d0
    !enddo
    !close(123)
    !===================================================
    CALL EndEvaluate_b()
  END SUBROUTINE GMG_CG_HART
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE cutoff_method(rho,Vh)
    USE FOURIER
    USE parameters, ONLY: Lcellbohr
    IMPLICIT NONE
    REAL(DP),intent(in)  :: rho(n1,n2,n3)
    REAL(DP),intent(out) :: Vh(n1,n2,n3)
    !> local
    INTEGER(I4B)  :: ig,ig1,ig2,ig3
    REAL(DP)      :: rho_recip(ng1,ng2,ng3)
    REAL(DP)      :: four_pi,sqrt3Lc

    !> initialize parameters
    four_pi=4.d0*pi
    sqrt3Lc=dsqrt(3.d0)*Lcellbohr

    !> FFT(n(r'))
    rho_recip=FFT(rho)
    !> n(g)*FFT(f(1/(r'-r)))
    ig=0
    do ig3=1,ng3,1
       do ig2=1,ng2,1
          do ig1=1,ng1,1
             ig=ig+1
             if(ig==1)cycle
             rho_recip(ig1,ig2,ig3)=rho_recip(ig1,ig2,ig3)*four_pi*(1.d0-dcos(grid%gvec(4,ig)*sqrt3Lc))/grid%gvec(4,ig)**2
          enddo
       enddo
    enddo
    print*,"ssssssasadad",grid%gvec(4,1)
    rho_recip(1,1,1)=rho_recip(1,1,1)*four_pi*1.5d0*Lcellbohr**2
    print*,"ssssssasadad",rho_recip(1:2,1,1)
    !> FFT(n(g)*f(g))
    Vh=FFT(rho_recip)
  ENDSUBROUTINE cutoff_method
END MODULE IsolateSet
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!------------------------MODULE-DI V IDER-LINE---------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!MODULE constants
!	IMPLICIT NONE
!	INTEGER,  PARAMETER      :: I4B = SELECTED_INT_KIND(9)
!	INTEGER,  PARAMETER      :: DP = KIND(1.0D0)
!	INTEGER,  PARAMETER      :: DCP = KIND((1.0D0,1.0D0))
!	REAL(DP), PARAMETER      :: pi=3.141592653589793238462643383279502884197_DP
!END MODULE constants
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!------------------------MODULE-DI V IDER-LINE---------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!MODULE struct_module
!	USE constants
!	IMPLICIT NONE
!	REAL(DP)             :: RadiusMax
!END MODULE struct_module
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!------------------------MODULE-DI V IDER-LINE---------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!MODULE grid_module
!	USE constants
!	USE struct_module
!	IMPLICIT NONE
!	REAL(DP)                 :: gap(3)
!	REAL(DP),pointer     :: dr(:,:,:)
!	REAL(DP),pointer     :: cost(:,:,:)
!	REAL(DP),pointer     :: sint(:,:,:)
!	COMPLEX(DCP),pointer :: cns(:,:,:)
!	INTEGER(I4B),pointer :: indx(:,:,:)
!CONTAINS
!	SUBROUTINE ReadGapAndRho(rho,VCoulomb)
!		IMPLICIT NONE
!		REAL(DP),ALLOCATABLE  :: rho(:,:,:),VCoulomb(:,:,:)
!		INTEGER(I4B)          :: n1,n2,n3
!		open(unit=2223, file="GapAndRho", status="UNKNOWN", action="read")
!		read(unit=2223, fmt=*) gap
!		read(unit=2223, fmt=*) n1
!		read(unit=2223, fmt=*) n2
!		read(unit=2223, fmt=*) n3
!		read(unit=2223, fmt=*) RadiusMax
!		allocate(rho(n1,n2,n3))
!		allocate(VCoulomb(n1,n2,n3))
!		read(unit=2223, fmt=*) rho
!		close(unit=2223)
!		!gap=0.320292564885347
!	END SUBROUTINE ReadGapAndRho
!END MODULE grid_module
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!------------------------MODULE-DI V IDER-LINE---------------------------
!-----------------------MODULE-  DIVID  ER-LINE--------------------------
!---------------------MODUL    E-DIVIDER    -LINE------------------------
!PROGRAM TestIso
!	USE constants
!	USE grid_module
!	USE IsolateSet
!	IMPLICIT NONE
!	REAL(DP),allocatable :: rho(:,:,:)
!	REAL(DP),allocatable :: VCoulomb(:,:,:)
!	CALL ReadGapAndRho(rho,VCoulomb)
!	CALL IsolateBoundry_a(rho,VCoulomb)
!	open(unit=2223, file="VCoulomb_1", status="UNKNOWN", action="write")
!	write(unit=2223, fmt=*) VCoulomb
!	print*,"VCoulomb writed"
!END PROGRAM TestIso
!BUGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!XXXXXXXXXXX        XXXXXX  XXXXXX  XXXXXX           XXXXXXXXXXX!
!XXXXXXXXXXX  XXXXX   XXXX  XXXXXX  XXXXXX  XXXXXXX  XXXXXXXXXXX!
!XXXXXXXXXXX        XXXXXX  XXXXXX  XXXXXX  XXXXXXXXXXXXXXXXXXXX!
!XXXXXXXXXXX  XXXX   XXXXX  XXXXXX  XXXXXX  XXXXX    XXXXXXXXXXX!
!XXXXXXXXXXX  XXXXX   XXXX  XXXXXX  XXXXXX  XXXXX    XXXXXXXXXXX!
!XXXXXXXXXXX         XXXXXX        XXXXXXX           XXXXXXXXXXX!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!BUG1: Text file expend some accuracy in read and write. Or, Reading and writing text files caused a certain amount of error.
!BUG2: An isolated unit cell is divided into n segments and n+1 grid points are required
!BUG3: a minus transport in a subroutine, calculating time increase large, oh,a large array operation in cycle!!!
