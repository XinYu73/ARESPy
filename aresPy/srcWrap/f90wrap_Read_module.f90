! Module read_module defined in file Read_module.f90

subroutine f90wrap_read_file(infile)
    use read_module, only: read_file
    implicit none
    
    character(7), intent(in) :: infile
    call read_file(infile=infile)
end subroutine f90wrap_read_file

subroutine f90wrap_read_poscar(nty, filename)
    use read_module, only: read_poscar
    implicit none
    
    integer(8), intent(in) :: nty
    character*(*), intent(in) :: filename
    call read_poscar(nty=nty, filename=filename)
end subroutine f90wrap_read_poscar

subroutine f90wrap_resetpos(natom, lat, pos, poscar, n0, n1)
    use read_module, only: resetpos
    implicit none
    
    integer(8), intent(in) :: natom
    real(8), intent(inout), dimension(3,3) :: lat
    real(8), intent(inout), dimension(3,n0) :: pos
    real(8), intent(inout), dimension(3,n1) :: poscar
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,1)
    integer :: n1
    !f2py intent(hide), depend(poscar) :: n1 = shape(poscar,1)
    call resetpos(natom=natom, lat=lat, pos=pos, poscar=poscar)
end subroutine f90wrap_resetpos

! End of module read_module defined in file Read_module.f90

