! Module write_module defined in file Write_module.fpp

subroutine f90wrap_write_poscar(filename, lat_mat, eleid, pos, atom_symbol, fixpos, n0, n1, n2, n3, n4, n5)
    use write_module, only: write_poscar
    implicit none
    
    character*(*) :: filename
    real(8), dimension(3,3) :: lat_mat
    integer(4), dimension(n0) :: eleid
    real(8), dimension(n1,n2) :: pos
    character*(*), optional, dimension(n3) :: atom_symbol
    logical, optional, dimension(n4,n5) :: fixpos
    integer :: n0
    !f2py intent(hide), depend(eleid) :: n0 = shape(eleid,0)
    integer :: n1
    !f2py intent(hide), depend(pos) :: n1 = shape(pos,0)
    integer :: n2
    !f2py intent(hide), depend(pos) :: n2 = shape(pos,1)
    integer :: n3
    !f2py intent(hide), depend(atom_symbol) :: n3 = shape(atom_symbol,0)
    integer :: n4
    !f2py intent(hide), depend(fixpos) :: n4 = shape(fixpos,0)
    integer :: n5
    !f2py intent(hide), depend(fixpos) :: n5 = shape(fixpos,1)
    call write_poscar(filename=filename, lat_mat=lat_mat, eleid=eleid, pos=pos, atom_symbol=atom_symbol, fixpos=fixpos)
end subroutine f90wrap_write_poscar

subroutine f90wrap_write_cif(filename, lat_para, nati, pos, atom_symbol, n0, n1, n2, n3)
    use write_module, only: write_cif
    implicit none
    
    character*(*) :: filename
    real(8), dimension(6) :: lat_para
    integer(4), dimension(n0) :: nati
    real(8), dimension(n1,n2) :: pos
    character*(*), optional, dimension(n3) :: atom_symbol
    integer :: n0
    !f2py intent(hide), depend(nati) :: n0 = shape(nati,0)
    integer :: n1
    !f2py intent(hide), depend(pos) :: n1 = shape(pos,0)
    integer :: n2
    !f2py intent(hide), depend(pos) :: n2 = shape(pos,1)
    integer :: n3
    !f2py intent(hide), depend(atom_symbol) :: n3 = shape(atom_symbol,0)
    call write_cif(filename=filename, lat_para=lat_para, nati=nati, pos=pos, atom_symbol=atom_symbol)
end subroutine f90wrap_write_cif

subroutine f90wrap_write3dat(filename, rho, n0, n1, n2)
    use write_module, only: write3dat
    implicit none
    
    character*(*) :: filename
    real(8), dimension(n0,n1,n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    call write3dat(filename=filename, rho=rho)
end subroutine f90wrap_write3dat

subroutine f90wrap_write3d(filename, rho, n0, n1, n2)
    use write_module, only: write3d
    implicit none
    
    character*(*) :: filename
    real(8), dimension(n0,n1,n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    call write3d(filename=filename, rho=rho)
end subroutine f90wrap_write3d

subroutine f90wrap_writedensity(filename, rho, n0, n1, n2)
    use write_module, only: writedensity
    implicit none
    
    character*(*) :: filename
    real(8), dimension(n0,n1,n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    call writedensity(filename=filename, rho=rho)
end subroutine f90wrap_writedensity

! End of module write_module defined in file Write_module.fpp

