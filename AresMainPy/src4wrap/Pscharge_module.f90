MODULE Pscharge_module
  USE constants
  USE Grid_module , ONLY : n,nsp,n1,n2,n3,dvol
  IMPLICIT NONE
CONTAINS
  SUBROUTINE Initial_Pscharge()
    USE math , ONLY : Norm,CubicHermiteInterp
    USE grid_module , ONLY : grid,n1,n2,n3,nb,Cubic2Sphere,Sphere2Cubic,gap
    USE struct_module , ONLY : volume,struct,naty,Eself,Eself_ref,Ec
    USE pspot_module , ONLY : psp
    USE compact , ONLY : PoissonCG 
    IMPLICIT NONE
    !IN/OUT
    !LOCAL
    INTEGER(I4B) :: ix,iy,iz,id,Ity,Ia
    REAL(DP) :: rvec(3),rR(3),rRnorm,dent,dena
    REAL(DP),DIMENSION(n1,n2,n3) :: vc &
         &, brho,brho_ref
    REAL(DP),DIMENSION(nb,nb,nb) :: Meshb,Meshv &
         &,Meshb_ref,Meshv_ref
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !initial
    Eself=0._DP
    Eself_ref=0._DP
    brho(:,:,:)=0._DP
    brho_ref(:,:,:)=0._DP
    vc(:,:,:)=0._DP
      !--------
    DO Ity=1,naty
       
       DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
          !get rhobI and vaI           
          CALL Pscharge(Ity,Ia,Meshb,Meshv,brho,struct%selfE(Ia),&
               &  Meshb_ref,Meshv_ref,brho_ref,struct%selfE_ref(Ia),vc)
          Eself=Eself+struct%selfE(Ia)
          Eself_ref=Eself_ref+struct%selfE_ref(Ia)
       ENDDO
    ENDDO
    !OPEN(123,FILE='nurho.dat')
    !    WRITE(123,*) n1-1,n2-1,n3-1
    !    WRITE(123,*) grid%brho(1:n1-1,1:n2-1,1:n3-1)
    !CLOSE(123)
    !STOP
    Ec=0.5_DP*SUM(vc*(brho+brho_ref))*dvol+Eself-Eself_ref
    !IF(ABS(Ec)<1e-10)  Ec=0._DP
    
    CALL Cubic2Sphere(brho,grid%brho1d)
    PRINT*,'Total PS-Charge (Z,E-self)',-SUM(brho)*dvol,Eself*hart2ev
    PRINT*,'reference PS-Charge (Z,E-self)',-SUM(brho_ref)*dvol,Eself_ref*hart2ev
    PRINT*,'Correction of repulsive energy',Ec*hart2ev
    !      CALL PoissonCG(grid%brho,vlp)
    !      CALL Cubic2Sphere(vlp,vlp1d)
!OPEN(123,FILE='vlpp.dat')
    !    WRITE(123,*) vlp1d(:)
    !    !READ(123,*) vlp1d(:)
    !CLOSE(123)
    !      PRINT*,'Total PS-Charge (Z,E-self)',0.5_DP*SUM(grid%brho*vlp)*dvol*hart2ev
    !      PRINT*,'Total PS-Charge (Z,E-self)',0.5_DP*SUM(grid%brho1d*vlp1d)*dvol !*hart2ev
    !STOP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_Pscharge
  !-----------------------------------------------------------------
  SUBROUTINE Pscharge(Ity,Ia,rhoa,vapp,brho,es,rhoa_ref,vapp_ref,brho_ref,es_ref,vc)
    USE math , ONLY : Norm,CubicHermiteInterp,CubicSplineInterp
    USE grid_module , ONLY : grid,n1,n2,n3,nb,gap
    USE struct_module , ONLY : volume,struct,naty,Eself,Eionion
    USE pspot_module , ONLY : psp
    USE compact , ONLY : lapla
    USE parameters , ONLY : NFCD,nord=>finite_order
    USE finite_module  , ONLY : real_nabla2_wb
    IMPLICIT NONE
    !IN/OUT
    INTEGER(I4B) :: Ity,Ia
    REAL(DP),INTENT(OUT) :: rhoa(nb,nb,nb),rhoa_ref(nb,nb,nb)
    REAL(DP),INTENT(OUT) :: vapp(nb,nb,nb),vapp_ref(nb,nb,nb)
    REAL(DP),INTENT(INOUT) :: brho(n1,n2,n3),brho_ref(n1,n2,n3),vc(n1,n2,n3)
    REAL(DP),INTENT(OUT) :: es,es_ref !self energy
    !LOCAL
    INTEGER(I4B) :: ix,iy,iz,i1,i2,i3
    INTEGER(I4B) :: is,ie,js,je,ks,ke
    REAL(DP) :: rvec(3),rR(3),rRnorm,vat
    REAL(DP) :: ztmp,zres
    REAL(DP) :: Ztol=1e-8
    REAL(DP),DIMENSION(-nord+1:nb+nord,-nord+1:nb+nord,-nord+1:nb+nord):: &
         & va,varef
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !initlize
    !x
    is=grid%Maindx(1,1,Ia)
    ie=grid%Maindx(1,2,Ia)
    !y
    js=grid%Maindx(2,1,Ia)
    je=grid%Maindx(2,2,Ia)
    !z
    ks=grid%Maindx(3,1,Ia)
    ke=grid%Maindx(3,2,Ia)
    !now get the potential
    DO i3=-nord+1,nb+nord
       iz=ks+i3-1
       DO i2=-nord+1,nb+nord
          iy=js+i2-1
          DO i1=-nord+1,nb+nord
             ix=is+i1-1
             !--------------------------------------------
             rvec(:)=grid%rVec(1:3,ix,iy,iz)
             rR(:)=rvec(:)-struct%poscar(:,Ia)
             rRnorm=Norm(rR)
             !vat=CubicHermiteInterp(psp(Ity)%vr,psp(Ity)%dvdr,psp(Ity)%rmax, &
             !     & psp(Ity)%rspacing, rRnorm,struct%Zion(Ity))
             vat=CubicSplineInterp(psp(Ity)%vr,psp(Ity)%ddvr,psp(Ity)%rmax, &
                  & psp(Ity)%rspacing, rRnorm,struct%Zion(Ity))
             !vatom
             va(i1,i2,i3)=vat 
             !vatom ref
             varef(i1,i2,i3)=refPot(rRNorm,psp(Ity)%rnoverlap,struct%Zion(Ity))
             !--------------------------------------------
          ENDDO
       ENDDO
    ENDDO
    !now get the pscharge brho
    !CALL lapla(NFCD,gap,va,rhoa) 
    !PS-charge density
    CALL real_nabla2_wb(va,nord,rhoa)
    rhoa=-rhoa/(4*pi)
    !PS-reference charge density
    CALL real_nabla2_wb(varef,nord,rhoa_ref)
    rhoa_ref=-rhoa_ref/(4*pi)
    !store potential
    vapp(:,:,:)=va(1:nb,1:nb,1:nb)
    vapp_ref(:,:,:)=varef(1:nb,1:nb,1:nb)
    !self energy
    es=0.5_DP*SUM(vapp*rhoa)*dvol
    es_ref=0.5_DP*SUM(vapp_ref*rhoa_ref)*dvol
    !to 3D mesh
    brho(is:ie,js:je,ks:ke)=brho(is:ie,js:je,ks:ke)+rhoa(1:nb,1:nb,1:nb)
    brho_ref(is:ie,js:je,ks:ke)=brho_ref(is:ie,js:je,ks:ke)+rhoa_ref(1:nb,1:nb,1:nb)
    !calculate the correct potential
    vc(is:ie,js:je,ks:ke)=vc(is:ie,js:je,ks:ke)+vapp_ref(:,:,:)-vapp(:,:,:)
    !PRINT*,'PsCharge test:',SUM(rhoa)*dvol,es*hart2ev
    ztmp=-SUM(rhoa)*dvol
    zres=ABS(ztmp-psp(Ity)%Zion)/psp(Ity)%Zion
    IF( zres>Ztol )THEN
       PRINT*,'PS-Charge have more round of err'
       PRINT*,'The Ps-Charge res:', zres
       PRINT*,'ztmp,Zion',ztmp,psp(Ity)%Zion
       !STOP
    ENDIF
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Pscharge
END MODULE Pscharge_module