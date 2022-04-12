MODULE struct_module
   !!####################################################!!
   !!*************  author: Qiang Xu          ***********!!
   !!*************  date  : 2017-07-11        ***********!!
   !!*************  Save Structure Information        ***!!
   !!####################################################!!
   USE constants
   IMPLICIT NONE
   TYPE struct_type
      REAL(DP),ALLOCATABLE      :: Zion (:)
      INTEGER(I4B),ALLOCATABLE  :: nati(:)  !number of atoms for each atom type
      INTEGER(I4B),ALLOCATABLE  :: eleid(:)   !atoms for each atom type id
      REAL(DP),ALLOCATABLE      :: pos(:,:)  !direct
      REAL(DP),ALLOCATABLE      :: poscar(:,:) ! car
      REAL(DP)                  :: stress(3,3)
      REAL(DP),ALLOCATABLE      :: forces(:,:)       ! total forces
      REAL(DP),ALLOCATABLE      :: mass(:)       ! atomic mass
      !for STO
      REAL(DP),ALLOCATABLE      :: zeta(:,:)    ! zeta for STO
      INTEGER(I4B),ALLOCATABLE  :: prinq(:,:)   !atom's principle quantum
      INTEGER(I4B),ALLOCATABLE  :: Lmax(:)   !max L need to consider
      CHARACTER(len=3),ALLOCATABLE :: elements(:)
      !for MO coefficient { coeff(atom_id,l,m,MO_id) }
      REAL(DP),ALLOCATABLE      :: coeff(:,:,:,:)
      REAL(DP),ALLOCATABLE      :: occupy(:)
      INTEGER(I4B)                  :: Noccupy
   END TYPE struct_type
   !-------------------------------------------------------------
   INTEGER(I4B)              :: natom! total number of atoms
   INTEGER(I4B)              :: naty ! number of atom types
   INTEGER(I4B)              :: ncharge    ! total change of the system
   REAL(DP)                  :: charge_ave ! ave change of the system
   REAL(DP)                  :: volume
   REAL(DP)                  :: lat_mat(3,3)
   REAL(DP)                  :: lat_para(6)
   REAL(DP)                  :: recip_lat(3,3) !!! the lattice matrix in reciprocal space
   REAL(DP)                  :: reclat_para(6)
   REAL(DP)                  :: energy(10)      ! energy
   TYPE(struct_type) :: struct
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !--------------------------------------------------------------
   SUBROUTINE creat_struct(numtyp,numatom)
     USE m_time_evaluate, only: memory_sum
      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: numtyp
      INTEGER(I4B),INTENT(IN) :: numatom
      !TYPE(struct_type),INTENT(OUT) :: struct
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      naty=numtyp
      natom=numatom
      CALL destroy_struct()
      ALLOCATE(struct%Zion(numtyp))
      ALLOCATE(struct%nati(numtyp))
      ALLOCATE(struct%eleid(numtyp+1))
      ALLOCATE(struct%pos(3,natom))
      ALLOCATE(struct%poscar(3,natom))
      ALLOCATE(struct%forces(3,natom))
      ALLOCATE(struct%mass(numtyp))
      !
      ALLOCATE(struct%zeta(4,numtyp))
      ALLOCATE(struct%prinq(4,numtyp))
      ALLOCATE(struct%Lmax(numtyp))
      ALLOCATE(struct%elements(numtyp))
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_sum("creat_struct",real(size(struct%Zion),DP)*DP&
           &+size(struct%nati)*I4B+size(struct%eleid)*I4B&
           &+size(struct%pos)*DP+size(struct%poscar)*DP&
           &+size(struct%forces)*DP+size(struct%mass)*DP&
           &+size(struct%zeta)*DP+size(struct%prinq)*I4B&
           &+size(struct%Lmax)*I4B+size(struct%elements)*3)
   ENDSUBROUTINE creat_struct
   !--------------------------------------------------------------
   SUBROUTINE destroy_struct()
     USE m_time_evaluate, ONLY: memory_free
      IMPLICIT NONE
      !
      !TYPE(struct_type),INTENT(INOUT) :: struct
      call memory_free("creat_struct",real(size(struct%Zion),DP)*DP&
           &+size(struct%nati)*I4B+size(struct%eleid)*I4B&
           &+size(struct%pos)*DP+size(struct%poscar)*DP&
           &+size(struct%forces)*DP+size(struct%mass)*DP&
           &+size(struct%zeta)*DP+size(struct%prinq)*I4B&
           &+size(struct%Lmax)*I4B+size(struct%elements)*3&
           &+size(struct%coeff)*DP)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ALLOCATED(struct%Zion))      DEALLOCATE(struct%Zion)
      IF(ALLOCATED(struct%nati))      DEALLOCATE(struct%nati)
      IF(ALLOCATED(struct%eleid))     DEALLOCATE(struct%eleid)
      IF(ALLOCATED(struct%pos))       DEALLOCATE(struct%pos)
      IF(ALLOCATED(struct%poscar))    DEALLOCATE(struct%poscar)
      IF(ALLOCATED(struct%forces))    DEALLOCATE(struct%forces)
      IF(ALLOCATED(struct%mass))    DEALLOCATE(struct%mass)
      IF(ALLOCATED(struct%zeta))    DEALLOCATE(struct%zeta)
      IF(ALLOCATED(struct%prinq))    DEALLOCATE(struct%prinq)
      IF(ALLOCATED(struct%Lmax))    DEALLOCATE(struct%Lmax)
      IF(ALLOCATED(struct%elements))  DEALLOCATE(struct%elements)
      IF(ALLOCATED(struct%coeff))  DEALLOCATE(struct%coeff)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_struct
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE struct_module
!---------------------------PARTING LINE--------------------------
