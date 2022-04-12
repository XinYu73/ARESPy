# 1 "Write_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Write_module.f90"
MODULE write_module
   !##########################################################!{{{
   !* CREATED_TIME  : 2015-07-20 20:09:18
   !* AUTHOR        : Xuecheng Shao
   !* CHANGE        : Xuecheng Shao
   !* Mail          : sxc@calypso.cn
   !* DESCRIPTION   :
   !     Output of something
   !* REFERENCES    :
   !     ------
   !* LOG           :
   !     2015-07-20
   !* LAST MODIFIED : 2015-07-30 05:38:54 PM
   !##########################################################!}}}
   USE constants

   use smpi_math_module

CONTAINS
SUBROUTINE write_POSCAR(filename, lat_mat, eleid, pos, atom_symbol,fixpos)!{{{
   !
   CHARACTER(LEN=*)                       :: filename
   REAL(DP)                               :: lat_mat(3,3)
   INTEGER(I4B)                           :: eleid(:)
   REAL(DP),DIMENSION(:,:)                :: pos
   CHARACTER(LEN=*),OPTIONAL,DIMENSION(:) :: atom_symbol
   LOGICAL,OPTIONAL,DIMENSION(:,:) :: fixpos
   !
   INTEGER(I4B)                           :: naty
   CHARACTER(LEN=20)             :: FMT
   !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>

   if (parallel%isroot) then

      naty = SIZE( eleid ) -1
      OPEN(UNIT=121, FILE=filename,STATUS="REPLACE", ACCESS="SEQUENTIAL",ACTION="WRITE")
      write(121,'(A)') filename
      write(121,'(A)') "1.0"
      write(121,'(3F20.12)') (lat_mat(:,i)*angs, i = 1, 3)
      !PRINT *,'kkkkkkkk'
      !PRINT *,lat_mat
      !print *, 'naty',naty
      IF ( PRESENT(atom_symbol) ) THEN
         write(FMT,*) naty
         WRITE(121,'('//trim(adjustl(FMT))//'A6)') (atom_symbol(i), i = 1, naty)
         !WRITE(121,'(<naty>A6)') (atom_symbol(i), i = 1, naty)
      endif
      write(FMT,*) naty
      write(121,'('//trim(adjustl(FMT))//'I10)') (eleid(i+1)-eleid(i), i = 1, naty)
      !write(121,'(<naty>I10)') (eleid(i+1)-eleid(i), i = 1, naty)
      WRITE(121,'(A)') 'Direct'
      DO i = 1, naty
         DO j = eleid(i), eleid(i+1)-1
            IF ( PRESENT(fixpos) ) THEN
               write(121,'(3F22.16,2X,3L3)') pos(:,j),fixpos(:,j)
            else
               write(121,'(3F22.16)') pos(:,j)
            endif
         enddo
      enddo
      close(121)

   endif

   !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
END SUBROUTINE!}}}
!---------------------------- DIVIDER LINE -----------------------------
SUBROUTINE write_cif(filename, lat_para, nati, pos, atom_symbol)!{{{
   !
   CHARACTER(LEN=*)                       :: filename
   REAL(DP),DIMENSION(6)                  :: lat_para
   INTEGER(I4B)                           :: nati(:)
   REAL(DP),DIMENSION(:,:)                :: pos
   CHARACTER(LEN=*),OPTIONAL,DIMENSION(:) :: atom_symbol
   !
   INTEGER(I4B)                           :: naty
   INTEGER(I4B)                           :: i,j,k,dig
   CHARACTER(LEN=20)             :: FMT
   !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
   naty = SIZE( nati )
   OPEN(UNIT=121, FILE=filename,STATUS="REPLACE", ACCESS="SEQUENTIAL",ACTION="WRITE")
   WRITE(121,*)
   write(121,'(A)') 'data_created_by_atlas'
   WRITE(121,*)
   write(121,"(A50)") "_pd_phase_name","'CIF files'"                                             
   write(121,'(A14,F36.8)') '_cell_length_a',    lat_para(1)
   write(121,'(A14,F36.8)') '_cell_length_b',    lat_para(2)
   write(121,'(A14,F36.8)') '_cell_length_c',    lat_para(3)
   write(121,'(A17,F33.4)') '_cell_angle_alpha', lat_para(4)
   write(121,'(A16,F34.4)') '_cell_angle_beta',  lat_para(5)
   write(121,'(A17,F33.4)') '_cell_angle_gamma', lat_para(6)
   write(121,'(A30,15X,A5)') "_symmetry_space_group_name_H-M", "'P 1'"
   write(121,'(A27,22X,A1)') "_symmetry_Int_Tables_number", "1"
   WRITE(121,*)
   write(121,'(A)')"loop_"
   write(121,'(A)')"_symmetry_equiv_pos_as_xyz"
   write(121,'(3X,A)')"'x, y, z'" 
   WRITE(121,*)
   write(121,'(A)')"loop_"
   write(121,'(3X,A)')"_atom_site_label"
   write(121,'(3X,A)')"_atom_site_occupancy"
   write(121,'(3X,A)')"_atom_site_fract_x"
   write(121,'(3X,A)')"_atom_site_fract_y"
   write(121,'(3X,A)')"_atom_site_fract_z"
   write(121,'(3X,A)')"_atom_site_thermal_displace_type"
   write(121,'(3X,A)')"_atom_site_U_iso_or_equiv"
   write(121,'(3X,A)')"_atom_site_type_symbol"
   k=1
   do i=1,naty
      do j=1, nati(i)
         IF ( j < 10 ) then
            dig = 1
         ELSE IF ( j < 100 ) then
            dig = 2
         ELSE IF ( j < 1000 ) then
            dig = 3
         ELSE IF ( j < 10000 ) then
            dig = 4
         ELSE IF ( j < 100000 ) then
            dig = 5
         ELSE
            dig = 10
         ENDIF
         write(FMT,*) diag
         !write(121, '(3x, A,I<dig>, 1x, F5.1, 1x, F12.6, 1x,  &
         write(121, '(3x, A,I'//trim(adjustl(FMT))//', 1x, F5.1, 1x, F12.6, 1x,  &
            & F12.6, 1x, F12.6, 1x, A6, 1x, F7.2, 1x, A)') &
            trim(adjustl(atom_symbol(i))), j, 1.0, pos(1,k), pos(2,k), pos(3,k), &
            "Uiso", 1.0, trim( adjustl(atom_symbol(i)) )
         k=k+1
      end do
   end do
   close(121)
   !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
END SUBROUTINE!}}}
!---------------------------- DIVIDER LINE -----------------------------

SUBROUTINE Write3DAt(filename, rho)!{{{
   !
   CHARACTER(LEN=*)                       :: filename
   REAL(DP), dimension(:,:,:) :: rho
   INTEGER(I4B)                  :: nxyz
   integer(kind=MPI_OFFSET_KIND) :: offset
   INTEGER(I4B)                  :: n1,n2,n3
   INTEGER(I4B)                  :: head
   INTEGER(I4B)                  :: gridn(3)
   !
   !open(unit=111, access="stream", action="write",file="den.bin", form="unformatted")       
   !!WRITE(111) grid%n1,grid%n2,grid%n3
   !WRITE(111)grid%rho
   !endfile(111)
   !close(111)
   !-----------------------------------------------------------------------
   gridn = shape(rho)
   nxyz = gridn(1)*gridn(2)*gridn(3)
   call MPI_FILE_OPEN(parallel%comm,trim(filename),MPI_MODE_WRONLY + MPI_MODE_CREATE, &
      MPI_INFO_NULL, funit,mpinfo)
   if (parallel%isroot) then
      offset=0
      gridn(3) = gridn(3) * parallel%numprocs
      call MPI_FILE_WRITE_AT(funit,offset,gridn, 3, MPI_INTEGER4, smpi_status, mpinfo) 
   endif
   head = 3 * I4B
   offset=parallel%myid * nxyz * DP + head
   call MPI_FILE_WRITE_AT(funit,offset,rho, nxyz, MPI_REAL8, smpi_status, mpinfo) 
   call MPI_FILE_CLOSE(funit,mpinfo) 
   print *, rho(1,1,1:3)
   !-----------------------------------------------------------------------
   END SUBROUTINE!}}}
!---------------------------- DIVIDER LINE -----------------------------
SUBROUTINE Write3D(filename, rho)!{{{
   !
   CHARACTER(LEN=*)                       :: filename
   REAL(DP), dimension(:,:,:) :: rho
   INTEGER(I4B)                  :: nxyz
   integer(kind=MPI_OFFSET_KIND) :: offset
   INTEGER(I4B)                  :: n1,n2,n3
   INTEGER(I4B)                  :: head
   INTEGER(I4B)                  :: gridn(3)
   !-----------------------------------------------------------------------
   gridn = shape(rho)
   nxyz = gridn(1)*gridn(2)*gridn(3)
   call MPI_FILE_OPEN(parallel%comm,trim(filename),MPI_MODE_WRONLY + MPI_MODE_CREATE, &
      MPI_INFO_NULL, funit,mpinfo)
   gridn(3) = SmpiSum(gridn(3))
   if (parallel%isroot) then
      offset=0
      call MPI_FILE_WRITE_AT(funit,offset,gridn, 3, MPI_INTEGER4, smpi_status, mpinfo) 
   endif
   head = 3 * I4B
   offset=parallel%myid * nxyz * DP + head
   call MPI_FILE_SET_VIEW(funit, offset, MPI_REAL8, MPI_REAL8, 'native', & 
      MPI_INFO_NULL, mpinfo) 
   call MPI_FILE_WRITE(funit,rho, nxyz, MPI_REAL8, MPI_STATUS_IGNORE, mpinfo) 
   call MPI_FILE_CLOSE(funit,mpinfo) 
   !-----------------------------------------------------------------------
END SUBROUTINE!}}}
!---------------------------- DIVIDER LINE -----------------------------
# 237 "Write_module.f90"
!---------------------------- DIVIDER LINE -----------------------------
SUBROUTINE WriteDensity(filename, rho)!{{{
   !
   CHARACTER(LEN=*)                       :: filename
   REAL(DP), dimension(:,:,:) :: rho



   call Write3D(filename,rho)

END SUBROUTINE WriteDensity!}}}

!---------------------------- DIVIDER LINE -----------------------------
END MODULE
