module ewald
  !##########################################################!{{{
  !* CREATED_TIME  : 2013-4-11
  !* AUTHOR        : Vincent Su
  !* CHANGE        :
  !* ADD           :
  !* DESCRIPTION   :
  !     This module caculates the ewald energy, ewald forces and ewald stress.
  !     Units: a.u.
  ! function ewald_energy(latticev,ionpositions,eachionnumbers,ioncharges)
  !   real(dp) latticev(3,3): The lattice vector. Latticev(:,1) is vector a.
  !   real(dp) ionpositions(3,:): Positions of each atom(fractional coordinates).
  !   INTEGER(I4B) eachionnumbers(:): Numbers of each ion in the cell.
  !   INTEGER(I4B) ioncharges(:): Charges of each ion. It equals the number of
  !     valance electron number(e.g., Si->4, Ga->2). You could get this number
  !     from pseudopotential file. Size(ioncharges) equals size(eachionnumbers).
  !   real(dp) ewald_energy: The answer.
  ! function ewald_forces(latticev,ionpositions,eachionnumbers,ioncharges)
  !   real(dp) ewald_forces(3,:) : The ewald forces of each ion.
  ! function ewald_stress(latticev,ionpositions,eachionnumbers,ioncharges)
  !   real(dp) ewald_stress(3,3): ewald_stress(1,1) is the stress in direction
  !     XX. ewald_stress(2,2)->YY, (3,3)->ZZ, (1,2)-XY, (2,3)->YZ, (1,3)->ZX.
  !     MADE IN CHINA, Changchun, Jilin University
  !* REFERENCES    :
  !
  !* LOG           :
  !     2015-05-08 :
  !* LAST MODIFIED : 2015-05-11 09:00:24 AM
  !##########################################################
  use constants,only : i4b,dp,pi,bohr2ang,rydberg

  implicit none

  type :: ion
     real(dp) :: charge
     real(dp) :: fcd(3)
  end type ion
  real(dp),parameter :: bohr=bohr2ang,hartreetoev=rydberg*2
  !real(dp),parameter :: bohr=0.529177249_dp,hartreetoev=27.21165200_dp  !vasp
  real(dp),private :: &
       sqrtpi=1.77245385090552_dp,&
       ert=1.0e-10_dp/hartreetoev,&
       maxrc=12._dp/bohr,&
       etainc=.01_dp*bohr,&
       recipinc=.01_dp*bohr
  real(dp),private :: &
       cellvol,&
       sumc,sumcs,&
       cellrecip(3,3),&
       eta,realct,recipct,realer,reciper

  INTEGER(I4B),private :: &
       nar,nbr,ncr,narp,nbrp,ncrp,&
       numion,numiontype

contains
  !}}}

  function ewald_energy(latticev,ionpositions,iontpid,ioncharges)!{{{
    real(dp),intent(in) :: latticev(3,3),ionpositions(:,:)
    INTEGER(I4B),intent(in)  :: iontpid(:)
    REAL(dp),intent(in)  :: ioncharges(:)
    real(dp) :: ewald_energy
    type(ion),allocatable :: iontb(:)
    INTEGER(I4B) :: i,k,natom

    natom=iontpid(size(iontpid))-1
    if(allocated(iontb))   deallocate(iontb)
    allocate(iontb(natom))
    do i=1,size(iontpid)-1
       iontb(iontpid(i):iontpid(i+1)-1)%charge=ioncharges(i)
    enddo
    k=1
    do i=1,natom
       iontb(k)%fcd=ionpositions(:,k)
       k=k+1
    end do

    ewald_energy=ewaldenergy(latticev,iontb,iontpid)
    if(allocated(iontb))      deallocate(iontb)
  end function ewald_energy!}}}
  !====================================================================================
  function ISO_ewald_energy(latticev,ionpositions,iontpid,ioncharges)!{{{
    real(dp),intent(in) :: latticev(3,3),ionpositions(:,:) !Cartesian Position
    INTEGER(I4B),intent(in)  :: iontpid(:)
    REAL(DP),intent(in)  :: ioncharges(:)
    real(dp) :: ISO_ewald_energy,rik,iienergy
    type(ion),allocatable :: iontb(:)
    INTEGER(I4B) :: i,k,natom

    natom=iontpid(size(iontpid))-1
    if(allocated(iontb))   deallocate(iontb)
    allocate(iontb(natom))
    do i=1,size(iontpid)-1
       iontb(iontpid(i):iontpid(i+1)-1)%charge=ioncharges(i)
    enddo
    k=1
    do i=1,natom
       iontb(k)%fcd=ionpositions(:,k)
       k=k+1
    end do
    !calculation
    iienergy=0.0d0
    DO i=1,natom
       DO k=1,natom
          IF(i/=k)THEN
             rik=sqrt(sum((iontb(i)%fcd-iontb(k)%fcd)**2))
             iienergy=iienergy+iontb(i)%charge*iontb(k)%charge/rik
          ENDIF
       ENDDO
    ENDDO
    ISO_ewald_energy=iienergy*0.5d0
    if(allocated(iontb))      deallocate(iontb)
  end function ISO_ewald_energy !}}}
  !
  function ISO_ewald_forces(latticev,ionpositions,iontpid,ioncharges)!{{{
#ifdef MPI
    USE smpi_math_module , ONLY : atom_split,parallel,smpi_reduce_sum_real_2d
#endif
    USE m_time_evaluate, ONLY: memory_sum,memory_free
    !> IN/OUT
    real(dp),intent(in) :: latticev(3,3),ionpositions(:,:) !Cartesian Position
    INTEGER(I4B),intent(in)  :: iontpid(:)
    REAL(DP),intent(in)  :: ioncharges(:)
    !> local
    real(dp) :: ISO_ewald_forces(3,iontpid(size(iontpid))-1),rik,force(3,iontpid(size(iontpid))-1)
    type(ion),allocatable :: iontb(:)
    !Euler angle
    real(DP) :: costheta, sintheta, cosphi, sinphi, temp_force, rxy
    INTEGER(I4B) :: i,k,natom
#ifdef MPI
    integer(I4B) :: mysize,atom_id=0,id_core
    INTEGER(I4B),allocatable :: atom_index(:)
#endif
    !==========================================
    call memory_sum('iso_ewald_local',real(size(force),DP)*2*DP)
    natom=iontpid(size(iontpid))-1
    if(allocated(iontb))then
       call memory_free('iso_ewald_iontb',real(size(iontb),DP)*DP)
       deallocate(iontb)
    endif
    allocate(iontb(natom))
    call memory_sum('iso_ewald_iontb',real(size(iontb),DP)*DP)
    do i=1,size(iontpid)-1
       iontb(iontpid(i):iontpid(i+1)-1)%charge=ioncharges(i)
    enddo
    k=1
    do i=1,natom
       iontb(k)%fcd=ionpositions(:,k)
       k=k+1
    end do
#ifdef MPI
    !> initial parallel config
    IF(.not.allocated(atom_index))THEN
       ALLOCATE(atom_index(natom/parallel%numprocs+1))
       call memory_sum('iso_ewald_atom_index',real(size(atom_index),DP)*DP)
    ENDIF
    ! CALL start_time('init_density_0')
    CALL atom_split(mysize,natom,atom_index)
    ! print *,'atom_index',atom_index,'id',parallel%myid
    ! CALL end_time('init_density_0')
    ! CALL write_time('init_density_0')
    id_core=1
#endif
    !> calculation
    force=0.0d0
    DO i=1,natom
#ifdef MPI
       If(i==atom_index(id_core))THEN
          if(id_core<mysize)id_core=id_core+1
#endif
          temp_force=0.0d0
          DO k=1,natom
             IF(i/=k)THEN
                rik=sqrt(sum((iontb(i)%fcd-iontb(k)%fcd)**2))
                rxy=sqrt((iontb(i)%fcd(1)-iontb(k)%fcd(1))**2+(iontb(i)%fcd(2)-iontb(k)%fcd(2))**2)
                IF(rxy.lt.0.000001d0)THEN
                   sinphi=0.d0
                   cosphi=0.d0
                   rxy=0.d0
                ELSE
                   sinphi=(iontb(k)%fcd(2)-iontb(i)%fcd(2))/rxy
                   cosphi=(iontb(k)%fcd(1)-iontb(i)%fcd(1))/rxy
                ENDIF
                IF(rik.lt.0.000001d0)stop "ATOM POSITION ERROR"
                sintheta=rxy/rik
                costheta=(iontb(k)%fcd(3)-iontb(i)%fcd(3))/rik
                temp_force=-iontb(k)%charge*iontb(i)%charge/(rik**2)*2
                force(1,i)=force(1,i)+temp_force*sintheta*cosphi
                force(2,i)=force(2,i)+temp_force*sintheta*sinphi
                force(3,i)=force(3,i)+temp_force*costheta
             ENDIF
          ENDDO !> k
#ifdef MPI
       endif
#endif
    ENDDO !> i
#ifdef MPI
    CALL smpi_reduce_sum_real_2d(force)
#endif
    ISO_ewald_forces=force
    if(allocated(iontb))then
       call memory_free('iso_ewald_iontb',real(size(iontb),DP)*DP)
       deallocate(iontb)
    endif
#ifdef MPI
    if(allocated(atom_index))then
       call memory_free('iso_ewald_atom_index',real(size(atom_index),DP)*DP)
       deallocate(atom_index)
    endif
#endif
    call memory_free('iso_ewald_local',real(size(force),DP)*2*DP)
  end function ISO_ewald_forces!}}}

  function ewald_forces(latticev,ionpositions,iontpid,ioncharges)!{{{
    real(dp),intent(in)  :: latticev(3,3),ionpositions(:,:)
    INTEGER(I4B),intent(in)   :: iontpid(:)
    REAL(DP),intent(in)   :: ioncharges(:)
    real(dp)             :: ewald_forces(3,iontpid(size(iontpid))-1)
    real(dp),allocatable :: ewald_f(:,:)
    type(ion),allocatable:: iontb(:)
    INTEGER(I4B) :: i,k,ionn

    ionn=iontpid(size(iontpid))-1
    if(allocated(iontb))   deallocate(iontb)
    if(allocated(ewald_f)) deallocate(ewald_f)
    allocate(iontb(ionn),ewald_f(3,ionn))
    do i=1,size(iontpid)-1
       iontb(iontpid(i):iontpid(i+1)-1)%charge=ioncharges(i)
    enddo
    k=1
    do i=1,ionn
       iontb(k)%fcd=ionpositions(:,k)
       k=k+1
    end do

    !call cewaldfcs(latticev,iontb,ewald_f)
    !ewald_forces=ewald_f
    call cewaldfcs(latticev,iontb,ewald_forces)
    if(allocated(iontb))   deallocate(iontb)
    if(allocated(ewald_f)) deallocate(ewald_f)
  end function ewald_forces!}}}

  function ewald_stress(latticev,ionpositions,iontpid,ioncharges)!{{{
    real(dp),intent(in) :: latticev(3,3),ionpositions(:,:)
    INTEGER(I4B),intent(in) :: iontpid(:)
    REAL(DP),intent(in) :: ioncharges(:)
    real(dp) :: ewald_stress(3,3)
    type(ion),allocatable :: iontb(:)
    INTEGER(I4B) :: i,j,k

    j=iontpid(size(iontpid))-1
    if(allocated(iontb))   deallocate(iontb)
    allocate(iontb(j))
    do i=1,size(iontpid)-1
       iontb(iontpid(i):iontpid(i+1)-1)%charge=ioncharges(i)
    enddo
    k=1
    do i=1,j
       iontb(k)%fcd=ionpositions(:,k)
       k=k+1
    end do

    ewald_stress=ewaldstr(latticev,iontb,iontpid)
    if(allocated(iontb))   deallocate(iontb)

  end function ewald_stress!}}}

  function ewaldenergy(cellr,iontb,iontpid)!{{{

    implicit none
    real(dp),intent(in) :: cellr(3,3)
    type(ion),intent(in) :: iontb(:)
    INTEGER(I4B),intent(in) :: iontpid(:)
    real(dp) :: ewaldenergy,ewaldarray(4)
    INTEGER(I4B) :: i,itp,&
         aind,nta1,nta2,bind,numiontemp,minla(1)
    real(dp) :: averagecs,realclen(3),recipclen(3)
    real(dp) :: temp_array(3)

    numion=size(iontb)
    numiontype=size(iontpid)-1
    !sqrtpi=sqrt(pi)
    realct=maxrc
    cellvol=volume(cellr)
    cellrecip=recipvector(cellr)
    sumc=0._dp
    sumcs=0._dp

    do itp=1,numiontype
       numiontemp=iontpid(itp+1)-iontpid(itp)
       sumc=sumc+numiontemp*iontb(iontpid(itp))%charge
       sumcs=sumcs+numiontemp*iontb(iontpid(itp))%charge**2
    end do
    averagecs=sumcs/numion
    !WRITE(6,*)sumc,sumcs,averagecs
    !WRITE(6,*)dp

    do i=1,3
       realclen(i)=vectorlength(cellr(:,i))
       recipclen(i)=vectorlength(cellrecip(:,i))
    end do
    !WRITE(6,*)realclen,recipclen

    minla=minloc(realclen)
    aind=minla(1)
    select case (aind)
    case(1)
       nta1=2
       nta2=3
    case(2)
       nta1=1
       nta2=3
    case(3)
       nta1=1
       nta2=2
    end select

    if(recipclen(nta1)<recipclen(nta2)) then
       bind=nta1
    else
       bind=nta2
    endif

    realct=min(maxrc,sqrt(sum(realclen**2))/2)
    nar=ceiling(realct*vectorlength(cellrecip(:,1))/pi/2)
    nbr=ceiling(realct*vectorlength(cellrecip(:,2))/pi/2)
    ncr=ceiling(realct*vectorlength(cellrecip(:,3))/pi/2)
    realct=min(nar/vectorlength(cellrecip(:,1)),&
         nbr/vectorlength(cellrecip(:,2)),&
         ncr/vectorlength(cellrecip(:,3)))*2*pi*.8239

    !nar=ceiling(realct*vectorlength(cellrecip(:,1))/pi/2)
    !nbr=ceiling(realct*vectorlength(cellrecip(:,2))/pi/2)
    !ncr=ceiling(realct*vectorlength(cellrecip(:,3))/pi/2)

    temp_array=crossp(cellr(:,aind),cellr(:,bind))
    eta=sqrt(pi*vectorlength(temp_array)&
         /realclen(aind)/cellvol)
    !eta=eta*9
    !WRITE(6,*)eta/bohr
    do
       realer=pi*numion**2*averagecs/cellvol/eta**2*erfc(realct*eta)
       if (realer>ert) exit
       eta=eta-etainc
    enddo
    do
       realer=pi*numion**2*averagecs/cellvol/eta**2*erfc(realct*eta)
       if (realer<ert) exit
       eta=eta+etainc
    enddo
    !WRITE(6,*)eta/bohr
    !WRITE(6,*)realer*hartreetoev

    recipct=0._dp
    do
       reciper=numion**2*averagecs/sqrtpi*eta*erfc(recipct/2/eta)
       if(reciper<ert) exit
       recipct=recipct+recipinc
    enddo
    !WRITE(6,*)reciper*hartreetoev
    !WRITE(6,*)eta,realct,recipct


    narp=ceiling(recipct*vectorlength(cellr(:,1))/pi/2)
    nbrp=ceiling(recipct*vectorlength(cellr(:,2))/pi/2)
    ncrp=ceiling(recipct*vectorlength(cellr(:,3))/pi/2)
    !write(51,*)nar,nbr,ncr,(2*nar+1)*(2*nbr+1)*(2*ncr+1)
    !write(51,*)narp,nbrp,ncrp,(2*narp+1)*(2*nbrp+1)*(2*ncrp+1)

    call ewaldrs
    ewaldarray(2)=-eta/sqrtpi*sumcs !selfterm
    call ewaldrps
    ewaldarray(4)=-pi/2/eta**2/cellvol*sumc**2
    ewaldenergy=sum(ewaldarray)
    !write(51,*)eta/bohr,realct*bohr,recipct/bohr
    !write(51,*)realer*hartreetoev,reciper*hartreetoev
    !write(51,*)ewaldarray(1:4)*hartreetoev,ewaldenergy*hartreetoev

  contains

    subroutine ewaldrs !real space

      implicit none
      INTEGER(I4B) :: i,ia,ib,ic,j,jtp
      real(dp) :: chargeij,distanceij,icod(3),jcod(3)
      real(dp) :: ewaldarr1

      !  WRITE(6,*)nar,nbr,ncr
      ewaldarray(1)=0
      ewaldarr1 = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:ewaldarr1) &
      !$OMP PRIVATE(i, icod,jtp,j,ia,ib,ic,jcod,distanceij)
      do i=1,numion
         icod=matmul(cellr,iontb(i)%fcd)
         do jtp=1,numiontype
            chargeij=iontb(i)%charge*iontb(iontpid(jtp))%charge
            !  WRITE(6,*)chargeij
            do j=iontpid(jtp),iontpid(jtp+1)-1
               do ia=-nar,nar
                  do ib=-nbr,nbr
                     do ic=-ncr,ncr
                        jcod=(ia+iontb(j)%fcd(1))*cellr(:,1)+&
                             (ib+iontb(j)%fcd(2))*cellr(:,2)+&
                             (ic+iontb(j)%fcd(3))*cellr(:,3)
                        temp_array=icod-jcod
                        distanceij=vectorlength(temp_array)
                        !  WRITE(6,*)distanceij,realct
                        if(distanceij<realct .and. &
                             (i/=j.or.ia/=0.or.ib/=0.or.ic/=0)) then
                           !WRITE(6,*)ewaldarray(1)
                           ewaldarr1=ewaldarr1+&
                                chargeij*erfc(eta*distanceij)/distanceij
                           !if (eta*distanceij>4) WRITE(6,*)eta*distanceij
                        endif
                     end do
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ewaldarray(1)=ewaldarr1/2._dp

    end subroutine ewaldrs

    subroutine ewaldrps

      INTEGER(I4B) :: ia,ib,ic,iblow,iclow,k
      real(dp) :: strfs,strfc,msqa,m(3),realcod(3)

      ewaldarray(3)=0._dp
      do ia=0,narp
         if(ia==0) then
            iblow=0
         else
            iblow=-nbrp
         end if
         do ib=iblow,nbrp
            if(ia==0.and.ib==0) then
               iclow=1
            else
               iclow=-ncrp
            endif
            do ic=iclow,ncrp
               m=ia*cellrecip(:,1)+ib*cellrecip(:,2)+ic*cellrecip(:,3)
               if (vectorlength(m)<=recipct) then
                  strfs=0._dp
                  strfc=0._dp
                  do k=1,numion
                     realcod=matmul(cellr,iontb(k)%fcd)
                     strfs=strfs+iontb(k)%charge*sin(dot_product(m,realcod))
                     strfc=strfc+iontb(k)%charge*cos(dot_product(m,realcod))
                  enddo
                  msqa=dot_product(m,m)
                  ewaldarray(3)=ewaldarray(3)+exp(-(msqa/eta**2/4._dp))&
                       /msqa*(strfs**2+strfc**2)
               endif
            end do
         enddo
      enddo
      ewaldarray(3)=ewaldarray(3)*4*pi/cellvol
    end subroutine ewaldrps

  end function ewaldenergy!}}}

  subroutine cewaldfcs(cellr,iontb,ewaldfcs)!{{{

    real(dp),intent(in) :: cellr(3,3)
    type(ion),intent(in) :: iontb(:)
    real(dp),intent(out) :: ewaldfcs(:,:)
    INTEGER(I4B) :: i
    real(dp) :: icod(3)
    real(dp) :: ewaldrsf(3)
    INTEGER(I4B) :: charge,ia,ib,ic,j
    real(dp) :: distanceij,etadistance,jcod(3),temp_array(3)

    ewaldfcs=0
    do i=1,numion
       icod=matmul(cellr,iontb(i)%fcd)
       ewaldrsf=0
       do j=1,numion
          do ia=-nar,nar
             do ib=-nbr,nbr
                do ic=-ncr,ncr
                   jcod=(ia+iontb(j)%fcd(1))*cellr(:,1)+&
                        (ib+iontb(j)%fcd(2))*cellr(:,2)+&
                        (ic+iontb(j)%fcd(3))*cellr(:,3)
                   temp_array=icod-jcod
                   distanceij=vectorlength(temp_array)
                   !WRITE(6,*) distanceij
                   if(distanceij<=realct.and.&
                        (i/=j.or.ia/=0.or.ib/=0.or.ic/=0)) then
                      etadistance=eta*distanceij
                      ewaldrsf=ewaldrsf+(icod-jcod)*iontb(j)%charge*&
                           (erfc(etadistance)/etadistance**3+2/sqrtpi &
                           *exp(-(etadistance**2))/etadistance**2)
                   endif
                enddo
             enddo
          enddo
       enddo
       ewaldrsf=ewaldrsf*iontb(i)%charge*eta**3
       ewaldfcs(:,i)=ewaldrsf+ewaldrpf(icod,iontb(i)%charge)
    enddo
  contains

    !function ewaldrsf(icod,charge)
    !real(dp) :: icod(3),ewaldrsf(3)
    !INTEGER(I4B) :: charge,ia,ib,ic,j
    !real(dp) :: distanceij,etadistance,jcod(3)

    !ewaldrsf=0
    !do j=1,numion
    !do ia=-nar,nar
    !do ib=-nbr,nbr
    !do ic=-ncr,ncr
    !jcod=(ia+iontb(j)%fcd(1))*cellr(:,1)+&
    !(ib+iontb(j)%fcd(2))*cellr(:,2)+&
    !(ic+iontb(j)%fcd(3))*cellr(:,3)
    !distanceij=vectorlength(icod-jcod)
    !!WRITE(6,*) distanceij
    !if(distanceij<=realct.and.&
    !(i/=j.or.ia/=0.or.ib/=0.or.ic/=0)) then
    !etadistance=eta*distanceij
    !ewaldrsf=ewaldrsf+(icod-jcod)*iontb(j)%charge*&
    !(erfc(etadistance)/etadistance**3+2/sqrtpi &
    !*exp(-(etadistance**2))/etadistance**2)
    !endif
    !enddo
    !enddo
    !enddo
    !enddo
    !ewaldrsf=ewaldrsf*charge*eta**3
    !WRITE(6,*)"eeee",ewaldrsf
    !end function

    function ewaldrpf(icod,charge)
      real(dp) :: icod(3),ewaldrpf(3)
      INTEGER(I4B) :: i,ia,ib,ic,iblow,iclow,k
      real(dp) :: strfs,strfc,msqa,m(3),realcod(3),charge

      ewaldrpf=0
      do ia=0,narp
         if(ia==0) then
            iblow=0
         else
            iblow=-nbrp
         end if
         do ib=iblow,nbrp
            if(ia==0.and.ib==0) then
               iclow=1
            else
               iclow=-ncrp
            endif
            do ic=iclow,ncrp
               !    WRITE(6,*)ic
               m=ia*cellrecip(:,1)+ib*cellrecip(:,2)+ic*cellrecip(:,3)
               if (vectorlength(m)<=recipct) then
                  strfs=0
                  strfc=0
                  do k=1,numion
                     realcod=matmul(cellr,iontb(k)%fcd)
                     strfs=strfs+iontb(k)%charge*sin(dot_product(m,realcod))
                     strfc=strfc+iontb(k)%charge*cos(dot_product(m,realcod))
                  enddo
                  msqa=dot_product(m,m)
                  ewaldrpf=ewaldrpf+exp(-(msqa/eta**2/4._dp))&
                       /msqa*(strfc*sin(dot_product(m,icod))-&
                       strfs*cos(dot_product(m,icod)))*m
                  !  WRITE(6,*)ewaldrpf
               endif
            end do
         enddo
      enddo
      ewaldrpf=ewaldrpf*8*pi/cellvol*charge
      !  WRITE(6,*)ewaldrpf The error is about 1.e-17(a.u.)

    end function ewaldrpf

  end subroutine cewaldfcs!}}}

  function ewaldstr(cellr, iontb, iontpid)!{{{
   real(dp) :: cellr(3, 3), ewaldstr(3, 3)
   type(ion) :: iontb(:)
   INTEGER(I4B) :: iontpid(:)
   real(dp) ::ewaldrstr(3, 3)
   INTEGER(I4B) :: a, b, i, ia, ib, ic, j, jtp
   real(dp) :: chargeij, distanceij, icod(3), jcod(3), tep&
        &, temp_array(3)
   real(dp) ::  ewaldrpstr(3, 3)
   INTEGER(I4B) :: iblow, iclow, k
   real(dp) :: strfs, strfc, msqa, m(3), realcod(3), tp1, tp2
   real(dp) :: ewaldavstr(3, 3)
   ewaldrstr = 0
   !!!!!!!!!!!!!!!!!!!!!!!!
   do i = 1, numion
       icod = matmul(cellr, iontb(i)%fcd)
       do jtp = 1, numiontype
           chargeij = iontb(i)%charge*iontb(iontpid(jtp))%charge
           !    WRITE(6,*)chargeij
           do j = iontpid(jtp), iontpid(jtp + 1) - 1
               do ia = -nar, nar
                   do ib = -nbr, nbr
                       do ic = -ncr, ncr
                           jcod = (ia + iontb(j)%fcd(1))*cellr(:, 1) + &
                                  (ib + iontb(j)%fcd(2))*cellr(:, 2) + &
                                  (ic + iontb(j)%fcd(3))*cellr(:, 3)
                           temp_array = icod - jcod
                           distanceij = vectorlength(temp_array)
                           !    WRITE(6,*)distanceij,realct
                           if (distanceij < realct .and. &
                               (i /= j .or. ia /= 0 .or. ib /= 0 .or. ic /= 0)) then
                               tep = chargeij*(erfc(eta*distanceij)/distanceij &
                                               + 2*eta*exp(-(eta*distanceij)**2)/sqrtpi)/distanceij**2
                               do a = 1, 3
                                   do b = a, 3
                                       ewaldrstr(a, b) = ewaldrstr(a, b) + tep &
                                                         *(icod(a) - jcod(a))*(icod(b) - jcod(b))
                                   end do
                               end do
                           end if
                       end do
                   end do
               end do
           end do
       end do
   end do
   ewaldrstr = -ewaldrstr/2/cellvol
   ewaldrpstr = 0
   do ia = 0, narp
       if (ia == 0) then
           iblow = 0
       else
           iblow = -nbrp
       end if
       do ib = iblow, nbrp
           if (ia == 0 .and. ib == 0) then
               iclow = 1
           else
               iclow = -ncrp
           end if
           do ic = iclow, ncrp
               m = ia*cellrecip(:, 1) + ib*cellrecip(:, 2) + ic*cellrecip(:, 3)
               if (vectorlength(m) <= recipct) then
                   strfs = 0._dp
                   strfc = 0._dp
                   do k = 1, numion
                       realcod = matmul(cellr, iontb(k)%fcd)
                       strfs = strfs + iontb(k)%charge*sin(dot_product(m, realcod))
                       strfc = strfc + iontb(k)%charge*cos(dot_product(m, realcod))
                   end do
                   msqa = dot_product(m, m)
                   tp1 = exp(-(msqa/eta**2/4))/msqa*(strfs**2 + strfc**2)
                   tp2 = tp1*2/msqa*(1 + msqa/eta**2/4)
                   do a = 1, 3
                       do b = a, 3
                           if (a == b) then
                               ewaldrpstr(a, b) = ewaldrpstr(a, b) + tp2*m(a)*m(b) - tp1
                           else
                               ewaldrpstr(a, b) = ewaldrpstr(a, b) + tp2*m(a)*m(b)
                           end if
                       end do
                   end do
               end if
           end do
       end do
   end do
   ewaldrpstr = ewaldrpstr*4*pi/cellvol**2
   ewaldavstr = 0
   do a = 1, 3
       ewaldavstr(a, a) = pi/2/eta**2/cellvol**2*sumc**2
   end do
   ewaldstr = ewaldrstr + ewaldrpstr + ewaldavstr
   !contains
   ! function ewaldrstr()
   !     real(dp) ::ewaldrstr(3, 3)
   !     INTEGER(I4B) :: a, b, i, ia, ib, ic, j, jtp
   !     real(dp) :: chargeij, distanceij, icod(3), jcod(3), tep&
   !          &, temp_array(3)

   !     ewaldrstr = 0
   !     do i = 1, numion
   !         icod = matmul(cellr, iontb(i)%fcd)
   !         do jtp = 1, numiontype
   !             chargeij = iontb(i)%charge*iontb(iontpid(jtp))%charge
   !             !    WRITE(6,*)chargeij
   !             do j = iontpid(jtp), iontpid(jtp + 1) - 1
   !                 do ia = -nar, nar
   !                     do ib = -nbr, nbr
   !                         do ic = -ncr, ncr
   !                             jcod = (ia + iontb(j)%fcd(1))*cellr(:, 1) + &
   !                                    (ib + iontb(j)%fcd(2))*cellr(:, 2) + &
   !                                    (ic + iontb(j)%fcd(3))*cellr(:, 3)
   !                             temp_array = icod - jcod
   !                             distanceij = vectorlength(temp_array)
   !                             !    WRITE(6,*)distanceij,realct
   !                             if (distanceij < realct .and. &
   !                                 (i /= j .or. ia /= 0 .or. ib /= 0 .or. ic /= 0)) then
   !                                 tep = chargeij*(erfc(eta*distanceij)/distanceij &
   !                                                 + 2*eta*exp(-(eta*distanceij)**2)/sqrtpi)/distanceij**2
   !                                 do a = 1, 3
   !                                     do b = a, 3
   !                                         ewaldrstr(a, b) = ewaldrstr(a, b) + tep &
   !                                                           *(icod(a) - jcod(a))*(icod(b) - jcod(b))
   !                                     end do
   !                                 end do
   !                             end if
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     ewaldrstr = -ewaldrstr/2/cellvol
   !     !WRITE(6,*)ewaldrstr
   ! end function ewaldrstr

   ! function ewaldrpstr()
   !     real(dp) ::  ewaldrpstr(3, 3)
   !     INTEGER(I4B) :: a, b, ia, ib, ic, iblow, iclow, k
   !     real(dp) :: strfs, strfc, msqa, m(3), realcod(3), tp1, tp2

   !     ewaldrpstr = 0
   !     do ia = 0, narp
   !         if (ia == 0) then
   !             iblow = 0
   !         else
   !             iblow = -nbrp
   !         end if
   !         do ib = iblow, nbrp
   !             if (ia == 0 .and. ib == 0) then
   !                 iclow = 1
   !             else
   !                 iclow = -ncrp
   !             end if
   !             do ic = iclow, ncrp
   !                 m = ia*cellrecip(:, 1) + ib*cellrecip(:, 2) + ic*cellrecip(:, 3)
   !                 if (vectorlength(m) <= recipct) then
   !                     strfs = 0._dp
   !                     strfc = 0._dp
   !                     do k = 1, numion
   !                         realcod = matmul(cellr, iontb(k)%fcd)
   !                         strfs = strfs + iontb(k)%charge*sin(dot_product(m, realcod))
   !                         strfc = strfc + iontb(k)%charge*cos(dot_product(m, realcod))
   !                     end do
   !                     msqa = dot_product(m, m)
   !                     tp1 = exp(-(msqa/eta**2/4))/msqa*(strfs**2 + strfc**2)
   !                     tp2 = tp1*2/msqa*(1 + msqa/eta**2/4)
   !                     do a = 1, 3
   !                         do b = a, 3
   !                             if (a == b) then
   !                                 ewaldrpstr(a, b) = ewaldrpstr(a, b) + tp2*m(a)*m(b) - tp1
   !                             else
   !                                 ewaldrpstr(a, b) = ewaldrpstr(a, b) + tp2*m(a)*m(b)
   !                             end if
   !                         end do
   !                     end do
   !                 end if
   !             end do
   !         end do
   !     end do
   !     ewaldrpstr = ewaldrpstr*4*pi/cellvol**2
   !     !WRITE(6,*)ewaldrpstr

   ! end function ewaldrpstr

   ! function ewaldavstr()

   !     real(dp) :: ewaldavstr(3, 3)
   !     INTEGER(I4B) :: a

   !     ewaldavstr = 0
   !     do a = 1, 3
   !         ewaldavstr(a, a) = pi/2/eta**2/cellvol**2*sumc**2
   !     end do
   !     !WRITE(6,*)ewaldavstr
   ! end function ewaldavstr
   end function ewaldstr!}}}

  function vectorlength(vc)!{{{
    real(dp) :: vc(3),vectorlength
    vectorlength=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
  end function vectorlength!}}}

  function recipvector(lat)!{{{
    real(dp),intent(in) :: lat(:,:)
    real(dp) :: recipvector(3,3)

    recipvector(:,1)=crossp(lat(:,2),lat(:,3))
    recipvector(:,2)=crossp(lat(:,3),lat(:,1))
    recipvector(:,3)=crossp(lat(:,1),lat(:,2))
    recipvector=recipvector/volume(lat)*pi*2._dp

  end function recipvector!}}}

  function volume(lat)!{{{
    real(dp),intent(in) :: lat(:,:)
    real(dp) :: volume

    volume=abs(sum(lat(:,1)*crossp(lat(:,2),lat(:,3))))

  end function volume!}}}

  function crossp(va,vb)!{{{
    real(dp),intent(in) :: va(3),vb(3)
    real(dp) :: crossp(3)

    crossp(1)=va(2)*vb(3)-va(3)*vb(2)
    crossp(2)=va(3)*vb(1)-va(1)*vb(3)
    crossp(3)=va(1)*vb(2)-va(2)*vb(1)
  end function crossp!}}}

  function erfc(x)!{{{
    !When -.2<x, the max absolute error is about 4.44e-15 (compared to MATLAB),
    !   for big x, the error is little.
    !   Reference:
    !     W.J. Cody.  Rational Chebyshev Approximations for the Error Function.
    !     Mathematics of Computation, 23(107):7, 1969.
    real(dp) :: x,erfc,temp,x2
    INTEGER(I4B) :: i
    real(dp),parameter :: &
         top1(5)= &
         (/ 1.857777061846031526730E-1_dp,&
         3.161123743870565596947E00_dp,&
         1.138641541510501556495E02_dp,&
         3.774852376853020208137E02_dp,&
         3.209377589138469472562E03_dp  /), &
         bot1(5)=&
         (/ 1.000000000000000000000E00_dp,&
         2.360129095234412093499E01_dp,&
         2.440246379344441733056E02_dp,&
         1.282616526077372275645E03_dp,&
         2.844236833439170622273E03_dp /)  ,&
         top2(9)=&
         (/ 2.15311535474403846343E-8_dp,&
         5.64188496988670089180E-1_dp,&
         8.88314979438837594118E00_dp,&
         6.61191906371416294775E01_dp,&
         2.98635138197400131132E02_dp,&
         8.81952221241769090411E02_dp,&
         1.71204761263407058314E03_dp,&
         2.05107837782607146532E03_dp,&
         1.23033935479799725272E03_dp /),   &
         bot2(9)=&
         (/ 1.00000000000000000000E00_dp,&
         1.57449261107098347253E01_dp,&
         1.17693950891312499305E02_dp,&
         5.37181101862009857509E02_dp,&
         1.62138957456669018874E03_dp,&
         3.29079923573345962678E03_dp,&
         4.36261909014324715820E03_dp,&
         3.43936767414372163696E03_dp,&
         1.23033935480374942043E03_dp /),&
         coef(15)=&
         (/ 13028445231.7429_dp,&
         -965070017.166138_dp,&
         77205601.3732910_dp,&
         -6713530.55419922_dp,&
         639383.862304688_dp,&
         -67303.5644531250_dp,&
         7918.06640625000_dp,&
         -1055.74218750000_dp,&
         162.421875000000_dp,&
         -29.5312500000000_dp,&
         6.56250000000000_dp,&
         -1.87500000000000_dp,&
         0.750000000000000_dp,&
         -0.500000000000000_dp,&
         1.00000000000000_dp /)

    if(x<.47_dp) then
       x2=x*x
       temp=top1(1)
       do i=2,5
          temp=temp*x2+top1(i)
       enddo
       erfc=x*temp
       temp=bot1(1)
       do i=2,5
          temp=temp*x2+bot1(i)
       enddo
       erfc=1-erfc/temp
    else if(x<4.3_dp) then
       temp=top2(1)
       do i=2,9
          temp=temp*x+top2(i)
       enddo
       erfc=exp(-x*x)*temp
       temp=bot2(1)
       do i=2,9
          temp=temp*x+bot2(i)
       enddo
       erfc=erfc/temp
    else
       erfc=coef(1)
       x2=x*x
       do i=2,15
          erfc=coef(i)+erfc/x2
       end do
       erfc=erfc*exp(-x2)/sqrtpi/x
    endif
  end function erfc!}}}

end module ewald
