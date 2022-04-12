! Module ewald defined in file Ewald.fpp

subroutine f90wrap_ion__get__charge(this, f90wrap_charge)
    use ewald, only: ion
    implicit none
    type ion_ptr_type
        type(ion), pointer :: p => NULL()
    end type ion_ptr_type
    integer, intent(in)   :: this(2)
    type(ion_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_charge
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_charge = this_ptr%p%charge
end subroutine f90wrap_ion__get__charge

subroutine f90wrap_ion__set__charge(this, f90wrap_charge)
    use ewald, only: ion
    implicit none
    type ion_ptr_type
        type(ion), pointer :: p => NULL()
    end type ion_ptr_type
    integer, intent(in)   :: this(2)
    type(ion_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_charge
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%charge = f90wrap_charge
end subroutine f90wrap_ion__set__charge

subroutine f90wrap_ion__array__fcd(this, nd, dtype, dshape, dloc)
    use ewald, only: ion
    implicit none
    type ion_ptr_type
        type(ion), pointer :: p => NULL()
    end type ion_ptr_type
    integer, intent(in) :: this(2)
    type(ion_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%fcd)
    dloc = loc(this_ptr%p%fcd)
end subroutine f90wrap_ion__array__fcd

subroutine f90wrap_ion_initialise(this)
    use ewald, only: ion
    implicit none
    
    type ion_ptr_type
        type(ion), pointer :: p => NULL()
    end type ion_ptr_type
    type(ion_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_ion_initialise

subroutine f90wrap_ion_finalise(this)
    use ewald, only: ion
    implicit none
    
    type ion_ptr_type
        type(ion), pointer :: p => NULL()
    end type ion_ptr_type
    type(ion_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_ion_finalise

subroutine f90wrap_ewald_energy(latticev, ionpositions, iontpid, ret_ewald_energy, ioncharges, n0, n1, n2, n3)
    use ewald, only: ewald_energy
    implicit none
    
    real(8), intent(in), dimension(3,3) :: latticev
    real(8), intent(in), dimension(n0,n1) :: ionpositions
    integer(4), intent(in), dimension(n2) :: iontpid
    real(8), intent(out) :: ret_ewald_energy
    real(8), intent(in), dimension(n3) :: ioncharges
    integer :: n0
    !f2py intent(hide), depend(ionpositions) :: n0 = shape(ionpositions,0)
    integer :: n1
    !f2py intent(hide), depend(ionpositions) :: n1 = shape(ionpositions,1)
    integer :: n2
    !f2py intent(hide), depend(iontpid) :: n2 = shape(iontpid,0)
    integer :: n3
    !f2py intent(hide), depend(ioncharges) :: n3 = shape(ioncharges,0)
    ret_ewald_energy = ewald_energy(latticev=latticev, ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
end subroutine f90wrap_ewald_energy

subroutine f90wrap_iso_ewald_energy(latticev, ionpositions, iontpid, ret_iso_ewald_energy, ioncharges, n0, n1, n2, n3)
    use ewald, only: iso_ewald_energy
    implicit none
    
    real(8), intent(in), dimension(3,3) :: latticev
    real(8), intent(in), dimension(n0,n1) :: ionpositions
    integer(4), intent(in), dimension(n2) :: iontpid
    real(8), intent(out) :: ret_iso_ewald_energy
    real(8), intent(in), dimension(n3) :: ioncharges
    integer :: n0
    !f2py intent(hide), depend(ionpositions) :: n0 = shape(ionpositions,0)
    integer :: n1
    !f2py intent(hide), depend(ionpositions) :: n1 = shape(ionpositions,1)
    integer :: n2
    !f2py intent(hide), depend(iontpid) :: n2 = shape(iontpid,0)
    integer :: n3
    !f2py intent(hide), depend(ioncharges) :: n3 = shape(ioncharges,0)
    ret_iso_ewald_energy = iso_ewald_energy(latticev=latticev, ionpositions=ionpositions, iontpid=iontpid, &
        ioncharges=ioncharges)
end subroutine f90wrap_iso_ewald_energy

subroutine f90wrap_iso_ewald_forces(latticev, ionpositions, iontpid, ret_iso_ewald_forces, ioncharges, n0, n1, n2, n3, &
    n4)
    use ewald, only: iso_ewald_forces
    implicit none
    
    real(8), intent(in), dimension(3,3) :: latticev
    real(8), intent(in), dimension(n0,n1) :: ionpositions
    integer(4), intent(in), dimension(n2) :: iontpid
    real(8), intent(out), dimension(3,n3) :: ret_iso_ewald_forces
    real(8), intent(in), dimension(n4) :: ioncharges
    integer :: n0
    !f2py intent(hide), depend(ionpositions) :: n0 = shape(ionpositions,0)
    integer :: n1
    !f2py intent(hide), depend(ionpositions) :: n1 = shape(ionpositions,1)
    integer :: n2
    !f2py intent(hide), depend(iontpid) :: n2 = shape(iontpid,0)
    integer :: n3
    integer :: n4
    !f2py intent(hide), depend(ioncharges) :: n4 = shape(ioncharges,0)
    ret_iso_ewald_forces = iso_ewald_forces(latticev=latticev, ionpositions=ionpositions, iontpid=iontpid, &
        ioncharges=ioncharges)
end subroutine f90wrap_iso_ewald_forces

subroutine f90wrap_ewald_forces(latticev, ionpositions, iontpid, ret_ewald_forces, ioncharges, n0, n1, n2, n3, n4)
    use ewald, only: ewald_forces
    implicit none
    
    real(8), intent(in), dimension(3,3) :: latticev
    real(8), intent(in), dimension(n0,n1) :: ionpositions
    integer(4), intent(in), dimension(n2) :: iontpid
    real(8), intent(out), dimension(3,n3) :: ret_ewald_forces
    real(8), intent(in), dimension(n4) :: ioncharges
    integer :: n0
    !f2py intent(hide), depend(ionpositions) :: n0 = shape(ionpositions,0)
    integer :: n1
    !f2py intent(hide), depend(ionpositions) :: n1 = shape(ionpositions,1)
    integer :: n2
    !f2py intent(hide), depend(iontpid) :: n2 = shape(iontpid,0)
    integer :: n3
    integer :: n4
    !f2py intent(hide), depend(ioncharges) :: n4 = shape(ioncharges,0)
    ret_ewald_forces = ewald_forces(latticev=latticev, ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
end subroutine f90wrap_ewald_forces

subroutine f90wrap_ewald_stress(latticev, ionpositions, iontpid, ret_ewald_stress, ioncharges, n0, n1, n2, n3)
    use ewald, only: ewald_stress
    implicit none
    
    real(8), intent(in), dimension(3,3) :: latticev
    real(8), intent(in), dimension(n0,n1) :: ionpositions
    integer(4), intent(in), dimension(n2) :: iontpid
    real(8), dimension(3,3), intent(out) :: ret_ewald_stress
    real(8), intent(in), dimension(n3) :: ioncharges
    integer :: n0
    !f2py intent(hide), depend(ionpositions) :: n0 = shape(ionpositions,0)
    integer :: n1
    !f2py intent(hide), depend(ionpositions) :: n1 = shape(ionpositions,1)
    integer :: n2
    !f2py intent(hide), depend(iontpid) :: n2 = shape(iontpid,0)
    integer :: n3
    !f2py intent(hide), depend(ioncharges) :: n3 = shape(ioncharges,0)
    ret_ewald_stress = ewald_stress(latticev=latticev, ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
end subroutine f90wrap_ewald_stress

subroutine f90wrap_ewaldrpstr(ret_ewaldrpstr, eta)
    use ewald, only: ewaldrpstr
    implicit none
    
    real(8), dimension(3,3), intent(out) :: ret_ewaldrpstr
    real(8) :: eta
    ret_ewaldrpstr = ewaldrpstr(eta=eta)
end subroutine f90wrap_ewaldrpstr

subroutine f90wrap_ewaldavstr(ret_ewaldavstr, eta)
    use ewald, only: ewaldavstr
    implicit none
    
    real(8), dimension(3,3), intent(out) :: ret_ewaldavstr
    real(8) :: eta
    ret_ewaldavstr = ewaldavstr(eta=eta)
end subroutine f90wrap_ewaldavstr

subroutine f90wrap_vectorlength(ret_vectorlength, vc)
    use ewald, only: vectorlength
    implicit none
    
    real(8), intent(out) :: ret_vectorlength
    real(8), dimension(3) :: vc
    ret_vectorlength = vectorlength(vc=vc)
end subroutine f90wrap_vectorlength

subroutine f90wrap_recipvector(ret_recipvector, lat, n0, n1)
    use ewald, only: recipvector
    implicit none
    
    real(8), dimension(3,3), intent(out) :: ret_recipvector
    real(8), intent(in), dimension(n0,n1) :: lat
    integer :: n0
    !f2py intent(hide), depend(lat) :: n0 = shape(lat,0)
    integer :: n1
    !f2py intent(hide), depend(lat) :: n1 = shape(lat,1)
    ret_recipvector = recipvector(lat=lat)
end subroutine f90wrap_recipvector

subroutine f90wrap_volume(ret_volume, lat, n0, n1)
    use ewald, only: volume
    implicit none
    
    real(8), intent(out) :: ret_volume
    real(8), intent(in), dimension(n0,n1) :: lat
    integer :: n0
    !f2py intent(hide), depend(lat) :: n0 = shape(lat,0)
    integer :: n1
    !f2py intent(hide), depend(lat) :: n1 = shape(lat,1)
    ret_volume = volume(lat=lat)
end subroutine f90wrap_volume

subroutine f90wrap_crossp(va, ret_crossp, vb)
    use ewald, only: crossp
    implicit none
    
    real(8), intent(in), dimension(3) :: va
    real(8), dimension(3), intent(out) :: ret_crossp
    real(8), intent(in), dimension(3) :: vb
    ret_crossp = crossp(va=va, vb=vb)
end subroutine f90wrap_crossp

subroutine f90wrap_erfc(ret_erfc, x)
    use ewald, only: erfc
    implicit none
    
    real(8), intent(out) :: ret_erfc
    real(8) :: x
    ret_erfc = erfc(x=x)
end subroutine f90wrap_erfc

subroutine f90wrap_ewald__get__bohr(f90wrap_bohr)
    use ewald, only: ewald_bohr => bohr
    implicit none
    real(8), intent(out) :: f90wrap_bohr
    
    f90wrap_bohr = ewald_bohr
end subroutine f90wrap_ewald__get__bohr

subroutine f90wrap_ewald__get__hartreetoev(f90wrap_hartreetoev)
    use ewald, only: ewald_hartreetoev => hartreetoev
    implicit none
    real(8), intent(out) :: f90wrap_hartreetoev
    
    f90wrap_hartreetoev = ewald_hartreetoev
end subroutine f90wrap_ewald__get__hartreetoev

! End of module ewald defined in file Ewald.fpp

