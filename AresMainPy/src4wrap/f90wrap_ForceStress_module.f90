! Module forcestress_module defined in file ForceStress_module.fpp

subroutine f90wrap_cal_force_c(rhos, uik, force, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    use forcestress_module, only: cal_force_c
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    complex(8), intent(in), dimension(n4,n5,n6,n7) :: uik
    real(8), intent(inout), dimension(n8,n9) :: force
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(uik) :: n4 = shape(uik,0)
    integer :: n5
    !f2py intent(hide), depend(uik) :: n5 = shape(uik,1)
    integer :: n6
    !f2py intent(hide), depend(uik) :: n6 = shape(uik,2)
    integer :: n7
    !f2py intent(hide), depend(uik) :: n7 = shape(uik,3)
    integer :: n8
    !f2py intent(hide), depend(force) :: n8 = shape(force,0)
    integer :: n9
    !f2py intent(hide), depend(force) :: n9 = shape(force,1)
    call cal_force_c(rhoS=rhos, Uik=uik, force=force)
end subroutine f90wrap_cal_force_c

subroutine f90wrap_cal_force_gamma(rhos, uik, force, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    use forcestress_module, only: cal_force_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(in), dimension(n4,n5,n6,n7) :: uik
    real(8), intent(inout), dimension(n8,n9) :: force
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(uik) :: n4 = shape(uik,0)
    integer :: n5
    !f2py intent(hide), depend(uik) :: n5 = shape(uik,1)
    integer :: n6
    !f2py intent(hide), depend(uik) :: n6 = shape(uik,2)
    integer :: n7
    !f2py intent(hide), depend(uik) :: n7 = shape(uik,3)
    integer :: n8
    !f2py intent(hide), depend(force) :: n8 = shape(force,0)
    integer :: n9
    !f2py intent(hide), depend(force) :: n9 = shape(force,1)
    call cal_force_gamma(rhoS=rhos, Uik=uik, force=force)
end subroutine f90wrap_cal_force_gamma

subroutine f90wrap_cal_force_r(rhos, uik, force, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use forcestress_module, only: cal_force_r
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(in), dimension(n4,n5,n6) :: uik
    real(8), intent(inout), dimension(n7,n8) :: force
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(uik) :: n4 = shape(uik,0)
    integer :: n5
    !f2py intent(hide), depend(uik) :: n5 = shape(uik,1)
    integer :: n6
    !f2py intent(hide), depend(uik) :: n6 = shape(uik,2)
    integer :: n7
    !f2py intent(hide), depend(force) :: n7 = shape(force,0)
    integer :: n8
    !f2py intent(hide), depend(force) :: n8 = shape(force,1)
    call cal_force_r(rhoS=rhos, Uik=uik, force=force)
end subroutine f90wrap_cal_force_r

subroutine f90wrap_locforce(rhos, lforce, n0, n1, n2, n3, n4, n5)
    use forcestress_module, only: locforce
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5) :: lforce
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(lforce) :: n4 = shape(lforce,0)
    integer :: n5
    !f2py intent(hide), depend(lforce) :: n5 = shape(lforce,1)
    call locforce(rhoS=rhos, lforce=lforce)
end subroutine f90wrap_locforce

subroutine f90wrap_locforce_r(rhos, lforce, n0, n1, n2, n3, n4, n5)
    use forcestress_module, only: locforce_r
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5) :: lforce
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(lforce) :: n4 = shape(lforce,0)
    integer :: n5
    !f2py intent(hide), depend(lforce) :: n5 = shape(lforce,1)
    call locforce_r(rhoS=rhos, lforce=lforce)
end subroutine f90wrap_locforce_r

subroutine f90wrap_nonlforce(uik, nlforce, n0, n1, n2, n3, n4, n5)
    use forcestress_module, only: nonlforce
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3) :: uik
    real(8), intent(inout), dimension(n4,n5) :: nlforce
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(uik) :: n3 = shape(uik,3)
    integer :: n4
    !f2py intent(hide), depend(nlforce) :: n4 = shape(nlforce,0)
    integer :: n5
    !f2py intent(hide), depend(nlforce) :: n5 = shape(nlforce,1)
    call nonlforce(Uik=uik, nlforce=nlforce)
end subroutine f90wrap_nonlforce

subroutine f90wrap_nonlforce_gamma(uik, nlforce, n0, n1, n2, n3, n4, n5)
    use forcestress_module, only: nonlforce_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: uik
    real(8), intent(inout), dimension(n4,n5) :: nlforce
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(uik) :: n3 = shape(uik,3)
    integer :: n4
    !f2py intent(hide), depend(nlforce) :: n4 = shape(nlforce,0)
    integer :: n5
    !f2py intent(hide), depend(nlforce) :: n5 = shape(nlforce,1)
    call nonlforce_gamma(Uik=uik, nlforce=nlforce)
end subroutine f90wrap_nonlforce_gamma

subroutine f90wrap_nonlforce_r_dg(uik, nlforce, n0, n1, n2, n3, n4)
    use forcestress_module, only: nonlforce_r_dg
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: uik
    real(8), intent(inout), dimension(n3,n4) :: nlforce
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(nlforce) :: n3 = shape(nlforce,0)
    integer :: n4
    !f2py intent(hide), depend(nlforce) :: n4 = shape(nlforce,1)
    call nonlforce_r_dg(Uik=uik, nlforce=nlforce)
end subroutine f90wrap_nonlforce_r_dg

subroutine f90wrap_nonlforce_r(uik, nlforce, n0, n1, n2, n3, n4)
    use forcestress_module, only: nonlforce_r
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: uik
    real(8), intent(inout), dimension(n3,n4) :: nlforce
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(nlforce) :: n3 = shape(nlforce,0)
    integer :: n4
    !f2py intent(hide), depend(nlforce) :: n4 = shape(nlforce,1)
    call nonlforce_r(Uik=uik, nlforce=nlforce)
end subroutine f90wrap_nonlforce_r

subroutine f90wrap_cal_stress(rhos, uik, stress, n0, n1, n2, n3, n4, n5, n6, n7)
    use forcestress_module, only: cal_stress
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    complex(8), intent(in), dimension(n4,n5,n6,n7) :: uik
    real(8), dimension(3,3), intent(inout) :: stress
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(uik) :: n4 = shape(uik,0)
    integer :: n5
    !f2py intent(hide), depend(uik) :: n5 = shape(uik,1)
    integer :: n6
    !f2py intent(hide), depend(uik) :: n6 = shape(uik,2)
    integer :: n7
    !f2py intent(hide), depend(uik) :: n7 = shape(uik,3)
    call cal_stress(rhoS=rhos, Uik=uik, stress=stress)
end subroutine f90wrap_cal_stress

subroutine f90wrap_cal_stress_gamma(rhos, uik, stress, n0, n1, n2, n3, n4, n5, n6, n7)
    use forcestress_module, only: cal_stress_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(in), dimension(n4,n5,n6,n7) :: uik
    real(8), dimension(3,3), intent(inout) :: stress
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(uik) :: n4 = shape(uik,0)
    integer :: n5
    !f2py intent(hide), depend(uik) :: n5 = shape(uik,1)
    integer :: n6
    !f2py intent(hide), depend(uik) :: n6 = shape(uik,2)
    integer :: n7
    !f2py intent(hide), depend(uik) :: n7 = shape(uik,3)
    call cal_stress_gamma(rhoS=rhos, Uik=uik, stress=stress)
end subroutine f90wrap_cal_stress_gamma

subroutine f90wrap_lda_stress(vxc, rho, ret_lda_stress, elda, n0, n1, n2, n3, n4, n5)
    use forcestress_module, only: lda_stress
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: vxc
    real(8), intent(in), dimension(n3,n4,n5) :: rho
    real(8), dimension(3,3), intent(out) :: ret_lda_stress
    real(8), intent(in) :: elda
    integer :: n0
    !f2py intent(hide), depend(vxc) :: n0 = shape(vxc,0)
    integer :: n1
    !f2py intent(hide), depend(vxc) :: n1 = shape(vxc,1)
    integer :: n2
    !f2py intent(hide), depend(vxc) :: n2 = shape(vxc,2)
    integer :: n3
    !f2py intent(hide), depend(rho) :: n3 = shape(rho,0)
    integer :: n4
    !f2py intent(hide), depend(rho) :: n4 = shape(rho,1)
    integer :: n5
    !f2py intent(hide), depend(rho) :: n5 = shape(rho,2)
    ret_lda_stress = lda_stress(vxc=vxc, rho=rho, Elda=elda)
end subroutine f90wrap_lda_stress

subroutine f90wrap_hart_stress(rhorecip, ret_hart_stress, eh, n0, n1, n2)
    use forcestress_module, only: hart_stress
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: rhorecip
    real(8), dimension(3,3), intent(out) :: ret_hart_stress
    real(8), intent(in) :: eh
    integer :: n0
    !f2py intent(hide), depend(rhorecip) :: n0 = shape(rhorecip,0)
    integer :: n1
    !f2py intent(hide), depend(rhorecip) :: n1 = shape(rhorecip,1)
    integer :: n2
    !f2py intent(hide), depend(rhorecip) :: n2 = shape(rhorecip,2)
    ret_hart_stress = hart_stress(rhoRecip=rhorecip, Eh=eh)
end subroutine f90wrap_hart_stress

subroutine f90wrap_kin_stress(uik, kinstress, n0, n1, n2, n3)
    use forcestress_module, only: kin_stress
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3) :: uik
    real(8), dimension(3,3), intent(inout) :: kinstress
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(uik) :: n3 = shape(uik,3)
    call kin_stress(Uik=uik, Kinstress=kinstress)
end subroutine f90wrap_kin_stress

subroutine f90wrap_kin_stress_gamma(uik, kinstress, n0, n1, n2, n3)
    use forcestress_module, only: kin_stress_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: uik
    real(8), dimension(3,3), intent(inout) :: kinstress
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(uik) :: n3 = shape(uik,3)
    call kin_stress_gamma(Uik=uik, Kinstress=kinstress)
end subroutine f90wrap_kin_stress_gamma

subroutine f90wrap_ion_nl_stress(uik, nlstress, n0, n1, n2, n3)
    use forcestress_module, only: ion_nl_stress
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3) :: uik
    real(8), dimension(3,3), intent(inout) :: nlstress
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(uik) :: n3 = shape(uik,3)
    call ion_nl_stress(Uik=uik, nlstress=nlstress)
end subroutine f90wrap_ion_nl_stress

subroutine f90wrap_ion_nl_stress_gamma(uik, nlstress, n0, n1, n2, n3)
    use forcestress_module, only: ion_nl_stress_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: uik
    real(8), dimension(3,3), intent(inout) :: nlstress
    integer :: n0
    !f2py intent(hide), depend(uik) :: n0 = shape(uik,0)
    integer :: n1
    !f2py intent(hide), depend(uik) :: n1 = shape(uik,1)
    integer :: n2
    !f2py intent(hide), depend(uik) :: n2 = shape(uik,2)
    integer :: n3
    !f2py intent(hide), depend(uik) :: n3 = shape(uik,3)
    call ion_nl_stress_gamma(Uik=uik, nlstress=nlstress)
end subroutine f90wrap_ion_nl_stress_gamma

subroutine f90wrap_ionele_stress(rhorecip, ret_ionele_stress, energy, n0, n1, n2)
    use forcestress_module, only: ionele_stress
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: rhorecip
    real(8), dimension(3,3), intent(out) :: ret_ionele_stress
    real(8), intent(in) :: energy
    integer :: n0
    !f2py intent(hide), depend(rhorecip) :: n0 = shape(rhorecip,0)
    integer :: n1
    !f2py intent(hide), depend(rhorecip) :: n1 = shape(rhorecip,1)
    integer :: n2
    !f2py intent(hide), depend(rhorecip) :: n2 = shape(rhorecip,2)
    ret_ionele_stress = ionele_stress(rhoRecip=rhorecip, energy=energy)
end subroutine f90wrap_ionele_stress

subroutine f90wrap_pseudopotdifflookup(ity, ret_pseudopotdifflookup, qnorm)
    use forcestress_module, only: pseudopotdifflookup
    implicit none
    
    integer(4), intent(in) :: ity
    real(8), intent(out) :: ret_pseudopotdifflookup
    real(8), intent(in) :: qnorm
    ret_pseudopotdifflookup = pseudopotdifflookup(Ity=ity, qNorm=qnorm)
end subroutine f90wrap_pseudopotdifflookup

subroutine f90wrap_cal_force_stress
    use forcestress_module, only: cal_force_stress
    implicit none
    
    call cal_force_stress()
end subroutine f90wrap_cal_force_stress

subroutine f90wrap_forcestress_module__get__cellpress(f90wrap_cellpress)
    use forcestress_module, only: forcestress_module_cellpress => cellpress
    implicit none
    real(8), intent(out) :: f90wrap_cellpress
    
    f90wrap_cellpress = forcestress_module_cellpress
end subroutine f90wrap_forcestress_module__get__cellpress

subroutine f90wrap_forcestress_module__set__cellpress(f90wrap_cellpress)
    use forcestress_module, only: forcestress_module_cellpress => cellpress
    implicit none
    real(8), intent(in) :: f90wrap_cellpress
    
    forcestress_module_cellpress = f90wrap_cellpress
end subroutine f90wrap_forcestress_module__set__cellpress

! End of module forcestress_module defined in file ForceStress_module.fpp

