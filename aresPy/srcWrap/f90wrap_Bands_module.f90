! Module band_structure defined in file Bands_module.f90

subroutine f90wrap_init_bandstruct(numk, kvec, n0, n1)
    use band_structure, only: init_bandstruct
    implicit none
    
    integer(8), intent(in) :: numk
    real(8), intent(in), dimension(n0,n1) :: kvec
    integer :: n0
    !f2py intent(hide), depend(kvec) :: n0 = shape(kvec,0)
    integer :: n1
    !f2py intent(hide), depend(kvec) :: n1 = shape(kvec,1)
    call init_bandstruct(numk=numk, kvec=kvec)
end subroutine f90wrap_init_bandstruct

subroutine f90wrap_band_begin(numk, kvec, n0, n1)
    use band_structure, only: band_begin
    implicit none
    
    integer(8), intent(in) :: numk
    real(8), intent(in), dimension(n0,n1) :: kvec
    integer :: n0
    !f2py intent(hide), depend(kvec) :: n0 = shape(kvec,0)
    integer :: n1
    !f2py intent(hide), depend(kvec) :: n1 = shape(kvec,1)
    call band_begin(numk=numk, kvec=kvec)
end subroutine f90wrap_band_begin

subroutine f90wrap_build_bands
    use band_structure, only: build_bands
    implicit none
    
    call build_bands()
end subroutine f90wrap_build_bands

subroutine f90wrap_read_bandpath(infile)
    use band_structure, only: read_bandpath
    implicit none
    
    character(11), intent(in) :: infile
    call read_bandpath(infile=infile)
end subroutine f90wrap_read_bandpath

subroutine f90wrap_band_pathread(infile)
    use band_structure, only: band_pathread
    implicit none
    
    character(11), intent(in) :: infile
    call band_pathread(infile=infile)
end subroutine f90wrap_band_pathread

subroutine f90wrap_read_density(infile, rho, n0, n1, n2)
    use band_structure, only: read_density
    implicit none
    
    character(8), intent(in) :: infile
    real(8), intent(inout), dimension(n0,n1,n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    call read_density(infile=infile, rho=rho)
end subroutine f90wrap_read_density

subroutine f90wrap_cal_band(rhos, nev, eigval, n0, n1, n2, n3, n4, n5)
    use band_structure, only: cal_band
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    integer(8), intent(in) :: nev
    real(8), intent(inout), dimension(n4,n5) :: eigval
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(eigval) :: n4 = shape(eigval,0)
    integer :: n5
    !f2py intent(hide), depend(eigval) :: n5 = shape(eigval,1)
    call cal_band(rhoS=rhos, nev=nev, EIGVAL=eigval)
end subroutine f90wrap_cal_band

! End of module band_structure defined in file Bands_module.f90

