! Module read_module defined in file Read_module.fpp

subroutine f90wrap_attribute__get__value(this, f90wrap_value)
    use read_module, only: attribute
    implicit none
    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    integer, intent(in)   :: this(2)
    type(attribute_ptr_type) :: this_ptr
    character(120), intent(out) :: f90wrap_value
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_value = this_ptr%p%value
end subroutine f90wrap_attribute__get__value

subroutine f90wrap_attribute__set__value(this, f90wrap_value)
    use read_module, only: attribute
    implicit none
    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    integer, intent(in)   :: this(2)
    type(attribute_ptr_type) :: this_ptr
    character(120), intent(in) :: f90wrap_value
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%value = f90wrap_value
end subroutine f90wrap_attribute__set__value

subroutine f90wrap_attribute_initialise(this)
    use read_module, only: attribute
    implicit none
    
    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    type(attribute_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_attribute_initialise

subroutine f90wrap_attribute_finalise(this)
    use read_module, only: attribute
    implicit none
    
    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    type(attribute_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_attribute_finalise

subroutine f90wrap_read_file(infile)
    use read_module, only: read_file
    implicit none
    
    character*(*), intent(in) :: infile
    call read_file(infile=infile)
end subroutine f90wrap_read_file

subroutine f90wrap_read_pos(nty, filename)
    use read_module, only: read_pos
    implicit none
    
    integer(4), intent(in) :: nty
    character*(*), intent(in) :: filename
    call read_pos(nty=nty, filename=filename)
end subroutine f90wrap_read_pos

subroutine f90wrap_resetlattice
    use read_module, only: resetlattice
    implicit none
    
    call resetlattice()
end subroutine f90wrap_resetlattice

subroutine f90wrap_read_pspot_atom(ity, filename, ps)
    use pspot_module, only: pspot
    use read_module, only: read_pspot_atom
    implicit none
    
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer(4), intent(in) :: ity
    character(30), intent(in) :: filename
    type(pspot_ptr_type) :: ps_ptr
    integer, intent(out), dimension(2) :: ps
    allocate(ps_ptr%p)
    call read_pspot_atom(Ity=ity, filename=filename, ps=ps_ptr%p)
    ps = transfer(ps_ptr, ps)
end subroutine f90wrap_read_pspot_atom

subroutine f90wrap_read_pspot(nty, filenames, n0)
    use read_module, only: read_pspot
    implicit none
    
    integer(4), intent(in) :: nty
    character(30), intent(in), dimension(n0) :: filenames
    integer :: n0
    !f2py intent(hide), depend(filenames) :: n0 = shape(filenames,0)
    call read_pspot(nty=nty, filenames=filenames)
end subroutine f90wrap_read_pspot

subroutine f90wrap_read_realpot_atom(ity, filename, ps)
    use pspot_module, only: pspot
    use read_module, only: read_realpot_atom
    implicit none
    
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer(4), intent(in) :: ity
    character(30), intent(in) :: filename
    type(pspot_ptr_type) :: ps_ptr
    integer, intent(out), dimension(2) :: ps
    allocate(ps_ptr%p)
    call read_realpot_atom(Ity=ity, filename=filename, ps=ps_ptr%p)
    ps = transfer(ps_ptr, ps)
end subroutine f90wrap_read_realpot_atom

subroutine f90wrap_exist_in(string1, ret_exist_in, string2)
    use read_module, only: exist_in
    implicit none
    
    character*(*) :: string1
    logical, intent(out) :: ret_exist_in
    character*(*) :: string2
    ret_exist_in = exist_in(string1=string1, string2=string2)
end subroutine f90wrap_exist_in

subroutine f90wrap_exist_ibegin(string1, ret_exist_ibegin, string2)
    use read_module, only: exist_ibegin
    implicit none
    
    character*(*) :: string1
    integer(4), intent(out) :: ret_exist_ibegin
    character*(*) :: string2
    ret_exist_ibegin = exist_ibegin(string1=string1, string2=string2)
end subroutine f90wrap_exist_ibegin

subroutine f90wrap_scan_head(file_unit, title, start_old)
    use read_module, only: scan_head
    implicit none
    
    integer(4) :: file_unit
    character*(*) :: title
    logical :: start_old
    call scan_head(file_unit=file_unit, title=title, start_old=start_old)
end subroutine f90wrap_scan_head

subroutine f90wrap_scan_tail(file_unit, title)
    use read_module, only: scan_tail
    implicit none
    
    integer(4) :: file_unit
    character*(*) :: title
    call scan_tail(file_unit=file_unit, title=title)
end subroutine f90wrap_scan_tail

subroutine f90wrap_read_upf(ity, filename, ps)
    use pspot_module, only: pspot
    use read_module, only: read_upf
    implicit none
    
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer(4), intent(in) :: ity
    character(30), intent(in) :: filename
    type(pspot_ptr_type) :: ps_ptr
    integer, intent(in), dimension(2) :: ps
    ps_ptr = transfer(ps, ps_ptr)
    call read_upf(Ity=ity, filename=filename, ps=ps_ptr%p)
end subroutine f90wrap_read_upf

subroutine f90wrap_read_pseudo_header(zion, mesh_size, nproj)
    use read_module, only: read_pseudo_header
    implicit none
    
    real(8) :: zion
    integer(4) :: mesh_size
    integer(4) :: nproj
    call read_pseudo_header(Zion=zion, mesh_size=mesh_size, nproj=nproj)
end subroutine f90wrap_read_pseudo_header

subroutine f90wrap_read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l, n0, n1, n2, n3, n4)
    use read_module, only: read_pseudo_nonlocal
    implicit none
    
    integer(4) :: unit_upf
    integer(4) :: nl
    real(8), dimension(n0,n1) :: beta_r
    real(8), dimension(n2,n3) :: d0
    real(8) :: rcut
    integer(4), dimension(n4) :: proj_l
    integer :: n0
    !f2py intent(hide), depend(beta_r) :: n0 = shape(beta_r,0)
    integer :: n1
    !f2py intent(hide), depend(beta_r) :: n1 = shape(beta_r,1)
    integer :: n2
    !f2py intent(hide), depend(d0) :: n2 = shape(d0,0)
    integer :: n3
    !f2py intent(hide), depend(d0) :: n3 = shape(d0,1)
    integer :: n4
    !f2py intent(hide), depend(proj_l) :: n4 = shape(proj_l,0)
    call read_pseudo_nonlocal(unit_UPF=unit_upf, nl=nl, beta_r=beta_r, D0=d0, rcut=rcut, proj_l=proj_l)
end subroutine f90wrap_read_pseudo_nonlocal

subroutine f90wrap_get_mo_coefficient
    use read_module, only: get_mo_coefficient
    implicit none
    
    call get_mo_coefficient()
end subroutine f90wrap_get_mo_coefficient

subroutine f90wrap_parse_headline(str, atom_sign, atom_id, atom_orbital)
    use read_module, only: parse_headline
    implicit none
    
    character(100) :: str
    character(3), dimension(3) :: atom_sign
    integer(4), dimension(3) :: atom_id
    character(6), dimension(3) :: atom_orbital
    call parse_headline(str=str, atom_sign=atom_sign, atom_id=atom_id, atom_orbital=atom_orbital)
end subroutine f90wrap_parse_headline

subroutine f90wrap_sign2lm(atom_orbital, l, m)
    use read_module, only: sign2lm
    implicit none
    
    character(6), intent(in), dimension(3) :: atom_orbital
    integer(4), dimension(3) :: l
    integer(4), dimension(3) :: m
    call sign2lm(atom_orbital=atom_orbital, l=l, m=m)
end subroutine f90wrap_sign2lm

subroutine f90wrap_destroy_moinit
    use read_module, only: destroy_moinit
    implicit none
    
    call destroy_moinit()
end subroutine f90wrap_destroy_moinit

subroutine f90wrap_out_concar(filename)
    use read_module, only: out_concar
    implicit none
    
    character*(*) :: filename
    call out_concar(filename=filename)
end subroutine f90wrap_out_concar

subroutine f90wrap_read_chgcar
    use read_module, only: read_chgcar
    implicit none
    
    call read_chgcar()
end subroutine f90wrap_read_chgcar

subroutine f90wrap_get_value_int(char_in, char_find, variable, find_flag)
    use read_module, only: get_value
    implicit none
    
    character(120), intent(in) :: char_in
    character*(*), intent(in) :: char_find
    integer(4) :: variable
    logical :: find_flag
    call get_value(char_in=char_in, char_find=char_find, variable=variable, find_flag=find_flag)
end subroutine f90wrap_get_value_int

subroutine f90wrap_get_value_real(char_in, char_find, variable, find_flag)
    use read_module, only: get_value
    implicit none
    
    character(120), intent(in) :: char_in
    character*(*), intent(in) :: char_find
    real(8) :: variable
    logical :: find_flag
    call get_value(char_in=char_in, char_find=char_find, variable=variable, find_flag=find_flag)
end subroutine f90wrap_get_value_real

! End of module read_module defined in file Read_module.fpp

