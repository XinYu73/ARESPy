! Module energy_module defined in file Energy_module.fpp

subroutine f90wrap_totalenergy(psi, rhos, eval, fmenergy, ets, llast, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use energy_module, only: totalenergy
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2,n3) :: psi
    real(8), intent(in), dimension(n4,n5,n6,n7) :: rhos
    real(8), intent(in), dimension(n8,n9,n10) :: eval
    real(8), intent(in) :: fmenergy
    real(8), intent(in) :: ets
    logical, intent(in) :: llast
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,3)
    integer :: n4
    !f2py intent(hide), depend(rhos) :: n4 = shape(rhos,0)
    integer :: n5
    !f2py intent(hide), depend(rhos) :: n5 = shape(rhos,1)
    integer :: n6
    !f2py intent(hide), depend(rhos) :: n6 = shape(rhos,2)
    integer :: n7
    !f2py intent(hide), depend(rhos) :: n7 = shape(rhos,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call totalenergy(psi=psi, rhoS=rhos, eval=eval, fmenergy=fmenergy, ets=ets, Llast=llast)
end subroutine f90wrap_totalenergy

subroutine f90wrap_out_ke_potential(efm, vie, vh, vxc, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use energy_module, only: out_ke_potential
    implicit none
    
    real(8) :: efm
    real(8), dimension(n0,n1,n2) :: vie
    real(8), dimension(n3,n4,n5) :: vh
    real(8), dimension(n6,n7,n8) :: vxc
    integer :: n0
    !f2py intent(hide), depend(vie) :: n0 = shape(vie,0)
    integer :: n1
    !f2py intent(hide), depend(vie) :: n1 = shape(vie,1)
    integer :: n2
    !f2py intent(hide), depend(vie) :: n2 = shape(vie,2)
    integer :: n3
    !f2py intent(hide), depend(vh) :: n3 = shape(vh,0)
    integer :: n4
    !f2py intent(hide), depend(vh) :: n4 = shape(vh,1)
    integer :: n5
    !f2py intent(hide), depend(vh) :: n5 = shape(vh,2)
    integer :: n6
    !f2py intent(hide), depend(vxc) :: n6 = shape(vxc,0)
    integer :: n7
    !f2py intent(hide), depend(vxc) :: n7 = shape(vxc,1)
    integer :: n8
    !f2py intent(hide), depend(vxc) :: n8 = shape(vxc,2)
    call out_ke_potential(Efm=efm, Vie=vie, Vh=vh, Vxc=vxc)
end subroutine f90wrap_out_ke_potential

subroutine f90wrap_totalenergy_gamma(psi, rhos, eval, fmenergy, ets, llast, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)
    use energy_module, only: totalenergy_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: psi
    real(8), intent(in), dimension(n4,n5,n6,n7) :: rhos
    real(8), intent(in), dimension(n8,n9,n10) :: eval
    real(8), intent(in) :: fmenergy
    real(8), intent(in) :: ets
    logical, intent(in) :: llast
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,3)
    integer :: n4
    !f2py intent(hide), depend(rhos) :: n4 = shape(rhos,0)
    integer :: n5
    !f2py intent(hide), depend(rhos) :: n5 = shape(rhos,1)
    integer :: n6
    !f2py intent(hide), depend(rhos) :: n6 = shape(rhos,2)
    integer :: n7
    !f2py intent(hide), depend(rhos) :: n7 = shape(rhos,3)
    integer :: n8
    !f2py intent(hide), depend(eval) :: n8 = shape(eval,0)
    integer :: n9
    !f2py intent(hide), depend(eval) :: n9 = shape(eval,1)
    integer :: n10
    !f2py intent(hide), depend(eval) :: n10 = shape(eval,2)
    call totalenergy_gamma(psi=psi, rhoS=rhos, eval=eval, fmenergy=fmenergy, ets=ets, Llast=llast)
end subroutine f90wrap_totalenergy_gamma

subroutine f90wrap_ehartree(rho, vhart, ehart, n0, n1, n2, n3, n4, n5)
    use energy_module, only: ehartree
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(in), dimension(n3,n4,n5) :: vhart
    real(8), intent(out) :: ehart
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(vhart) :: n3 = shape(vhart,0)
    integer :: n4
    !f2py intent(hide), depend(vhart) :: n4 = shape(vhart,1)
    integer :: n5
    !f2py intent(hide), depend(vhart) :: n5 = shape(vhart,2)
    call ehartree(rho=rho, vhart=vhart, Ehart=ehart)
end subroutine f90wrap_ehartree

subroutine f90wrap_ehartree_iso(rho, vhart, ehart, n0, n1)
    use energy_module, only: ehartree_iso
    implicit none
    
    real(8), intent(in), dimension(n0) :: rho
    real(8), intent(in), dimension(n1) :: vhart
    real(8), intent(out) :: ehart
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(vhart) :: n1 = shape(vhart,0)
    call ehartree_iso(rho=rho, vhart=vhart, Ehart=ehart)
end subroutine f90wrap_ehartree_iso

subroutine f90wrap_e_ie(rho, vion, eie, n0, n1, n2, n3, n4, n5)
    use energy_module, only: e_ie
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(in), dimension(n3,n4,n5) :: vion
    real(8), intent(out) :: eie
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(vion) :: n3 = shape(vion,0)
    integer :: n4
    !f2py intent(hide), depend(vion) :: n4 = shape(vion,1)
    integer :: n5
    !f2py intent(hide), depend(vion) :: n5 = shape(vion,2)
    call e_ie(rho=rho, vion=vion, Eie=eie)
end subroutine f90wrap_e_ie

subroutine f90wrap_ebands(eval, wke, eband, n0, n1, n2, n3, n4, n5)
    use energy_module, only: ebands
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: eval
    real(8), intent(in), dimension(n3,n4,n5) :: wke
    real(8), intent(out) :: eband
    integer :: n0
    !f2py intent(hide), depend(eval) :: n0 = shape(eval,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,1)
    integer :: n2
    !f2py intent(hide), depend(eval) :: n2 = shape(eval,2)
    integer :: n3
    !f2py intent(hide), depend(wke) :: n3 = shape(wke,0)
    integer :: n4
    !f2py intent(hide), depend(wke) :: n4 = shape(wke,1)
    integer :: n5
    !f2py intent(hide), depend(wke) :: n5 = shape(wke,2)
    call ebands(eval=eval, wke=wke, Eband=eband)
end subroutine f90wrap_ebands

subroutine f90wrap_bvk_totalenergy(rhos, eval, fmenergy, ets, llast, n0, n1, n2, n3, n4, n5)
    use energy_module, only: bvk_totalenergy
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(in), dimension(n4,n5) :: eval
    real(8), intent(in) :: fmenergy
    real(8), intent(in) :: ets
    logical, intent(in) :: llast
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,1)
    call bvk_totalenergy(rhoS=rhos, eval=eval, fmenergy=fmenergy, ets=ets, Llast=llast)
end subroutine f90wrap_bvk_totalenergy

subroutine f90wrap_bvk_ebands(eval, focc, eband, n0, n1, n2, n3)
    use energy_module, only: bvk_ebands
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: eval
    real(8), intent(in), dimension(n2,n3) :: focc
    real(8), intent(out) :: eband
    integer :: n0
    !f2py intent(hide), depend(eval) :: n0 = shape(eval,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,1)
    integer :: n2
    !f2py intent(hide), depend(focc) :: n2 = shape(focc,0)
    integer :: n3
    !f2py intent(hide), depend(focc) :: n3 = shape(focc,1)
    call bvk_ebands(eval=eval, focc=focc, Eband=eband)
end subroutine f90wrap_bvk_ebands

subroutine f90wrap_prr_ebands(veff, phi, pbar, eband, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    use energy_module, only: prr_ebands
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: veff
    real(8), intent(in), dimension(n4,n5,n6) :: phi
    real(8), intent(in), dimension(n7,n8,n9) :: pbar
    real(8), intent(out) :: eband
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(veff) :: n3 = shape(veff,3)
    integer :: n4
    !f2py intent(hide), depend(phi) :: n4 = shape(phi,0)
    integer :: n5
    !f2py intent(hide), depend(phi) :: n5 = shape(phi,1)
    integer :: n6
    !f2py intent(hide), depend(phi) :: n6 = shape(phi,2)
    integer :: n7
    !f2py intent(hide), depend(pbar) :: n7 = shape(pbar,0)
    integer :: n8
    !f2py intent(hide), depend(pbar) :: n8 = shape(pbar,1)
    integer :: n9
    !f2py intent(hide), depend(pbar) :: n9 = shape(pbar,2)
    call prr_ebands(veff=veff, Phi=phi, Pbar=pbar, Eband=eband)
end subroutine f90wrap_prr_ebands

subroutine f90wrap_prr_totalenergy(rhos, phi, pbar, fmenergy, ets, llast, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    use energy_module, only: prr_totalenergy
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(in), dimension(n4,n5,n6) :: phi
    real(8), intent(in), dimension(n7,n8,n9) :: pbar
    real(8), intent(in) :: fmenergy
    real(8), intent(in) :: ets
    logical, intent(in) :: llast
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(phi) :: n4 = shape(phi,0)
    integer :: n5
    !f2py intent(hide), depend(phi) :: n5 = shape(phi,1)
    integer :: n6
    !f2py intent(hide), depend(phi) :: n6 = shape(phi,2)
    integer :: n7
    !f2py intent(hide), depend(pbar) :: n7 = shape(pbar,0)
    integer :: n8
    !f2py intent(hide), depend(pbar) :: n8 = shape(pbar,1)
    integer :: n9
    !f2py intent(hide), depend(pbar) :: n9 = shape(pbar,2)
    call prr_totalenergy(rhoS=rhos, Phi=phi, Pbar=pbar, fmenergy=fmenergy, ets=ets, Llast=llast)
end subroutine f90wrap_prr_totalenergy

subroutine f90wrap_iso_totalenergy(rhos, eval, fmenergy, ets, llast, n0, n1, n2, n3)
    use energy_module, only: iso_totalenergy
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(in), dimension(n2,n3) :: eval
    real(8), intent(in) :: fmenergy
    real(8), intent(in) :: ets
    logical, intent(in) :: llast
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(eval) :: n2 = shape(eval,0)
    integer :: n3
    !f2py intent(hide), depend(eval) :: n3 = shape(eval,1)
    call iso_totalenergy(rhoS=rhos, eval=eval, fmenergy=fmenergy, ets=ets, Llast=llast)
end subroutine f90wrap_iso_totalenergy

subroutine f90wrap_iso_ebands(eval, focc, eband, n0, n1, n2, n3)
    use energy_module, only: iso_ebands
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: eval
    real(8), intent(in), dimension(n2,n3) :: focc
    real(8), intent(out) :: eband
    integer :: n0
    !f2py intent(hide), depend(eval) :: n0 = shape(eval,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,1)
    integer :: n2
    !f2py intent(hide), depend(focc) :: n2 = shape(focc,0)
    integer :: n3
    !f2py intent(hide), depend(focc) :: n3 = shape(focc,1)
    call iso_ebands(eval=eval, focc=focc, Eband=eband)
end subroutine f90wrap_iso_ebands

subroutine f90wrap_iso_e_ie(rho, vionlpp, eie, n0, n1)
    use energy_module, only: iso_e_ie
    implicit none
    
    real(8), intent(in), dimension(n0) :: rho
    real(8), intent(in), dimension(n1) :: vionlpp
    real(8), intent(out) :: eie
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(vionlpp) :: n1 = shape(vionlpp,0)
    call iso_e_ie(rho=rho, vionlpp=vionlpp, Eie=eie)
end subroutine f90wrap_iso_e_ie

subroutine f90wrap_lda_energy(ret_lda_energy, rhoreal, n0, n1, n2)
    use energy_module, only: lda_energy
    implicit none
    
    real(8), intent(out) :: ret_lda_energy
    real(8), intent(in), dimension(n0,n1,n2) :: rhoreal
    integer :: n0
    !f2py intent(hide), depend(rhoreal) :: n0 = shape(rhoreal,0)
    integer :: n1
    !f2py intent(hide), depend(rhoreal) :: n1 = shape(rhoreal,1)
    integer :: n2
    !f2py intent(hide), depend(rhoreal) :: n2 = shape(rhoreal,2)
    ret_lda_energy = lda_energy(rhoReal=rhoreal)
end subroutine f90wrap_lda_energy

subroutine f90wrap_energy_module__get__Etot(f90wrap_Etot)
    use energy_module, only: energy_module_Etot => Etot
    implicit none
    real(8), intent(out) :: f90wrap_Etot
    
    f90wrap_Etot = energy_module_Etot
end subroutine f90wrap_energy_module__get__Etot

subroutine f90wrap_energy_module__set__Etot(f90wrap_Etot)
    use energy_module, only: energy_module_Etot => Etot
    implicit none
    real(8), intent(in) :: f90wrap_Etot
    
    energy_module_Etot = f90wrap_Etot
end subroutine f90wrap_energy_module__set__Etot

subroutine f90wrap_energy_module__get__Ekine(f90wrap_Ekine)
    use energy_module, only: energy_module_Ekine => Ekine
    implicit none
    real(8), intent(out) :: f90wrap_Ekine
    
    f90wrap_Ekine = energy_module_Ekine
end subroutine f90wrap_energy_module__get__Ekine

subroutine f90wrap_energy_module__set__Ekine(f90wrap_Ekine)
    use energy_module, only: energy_module_Ekine => Ekine
    implicit none
    real(8), intent(in) :: f90wrap_Ekine
    
    energy_module_Ekine = f90wrap_Ekine
end subroutine f90wrap_energy_module__set__Ekine

subroutine f90wrap_energy_module__get__Ehart(f90wrap_Ehart)
    use energy_module, only: energy_module_Ehart => Ehart
    implicit none
    real(8), intent(out) :: f90wrap_Ehart
    
    f90wrap_Ehart = energy_module_Ehart
end subroutine f90wrap_energy_module__get__Ehart

subroutine f90wrap_energy_module__set__Ehart(f90wrap_Ehart)
    use energy_module, only: energy_module_Ehart => Ehart
    implicit none
    real(8), intent(in) :: f90wrap_Ehart
    
    energy_module_Ehart = f90wrap_Ehart
end subroutine f90wrap_energy_module__set__Ehart

subroutine f90wrap_energy_module__get__Exc(f90wrap_Exc)
    use energy_module, only: energy_module_Exc => Exc
    implicit none
    real(8), intent(out) :: f90wrap_Exc
    
    f90wrap_Exc = energy_module_Exc
end subroutine f90wrap_energy_module__get__Exc

subroutine f90wrap_energy_module__set__Exc(f90wrap_Exc)
    use energy_module, only: energy_module_Exc => Exc
    implicit none
    real(8), intent(in) :: f90wrap_Exc
    
    energy_module_Exc = f90wrap_Exc
end subroutine f90wrap_energy_module__set__Exc

subroutine f90wrap_energy_module__get__Eband(f90wrap_Eband)
    use energy_module, only: energy_module_Eband => Eband
    implicit none
    real(8), intent(out) :: f90wrap_Eband
    
    f90wrap_Eband = energy_module_Eband
end subroutine f90wrap_energy_module__get__Eband

subroutine f90wrap_energy_module__set__Eband(f90wrap_Eband)
    use energy_module, only: energy_module_Eband => Eband
    implicit none
    real(8), intent(in) :: f90wrap_Eband
    
    energy_module_Eband = f90wrap_Eband
end subroutine f90wrap_energy_module__set__Eband

subroutine f90wrap_energy_module__get__Eie(f90wrap_Eie)
    use energy_module, only: energy_module_Eie => Eie
    implicit none
    real(8), intent(out) :: f90wrap_Eie
    
    f90wrap_Eie = energy_module_Eie
end subroutine f90wrap_energy_module__get__Eie

subroutine f90wrap_energy_module__set__Eie(f90wrap_Eie)
    use energy_module, only: energy_module_Eie => Eie
    implicit none
    real(8), intent(in) :: f90wrap_Eie
    
    energy_module_Eie = f90wrap_Eie
end subroutine f90wrap_energy_module__set__Eie

subroutine f90wrap_energy_module__get__Eienl(f90wrap_Eienl)
    use energy_module, only: energy_module_Eienl => Eienl
    implicit none
    real(8), intent(out) :: f90wrap_Eienl
    
    f90wrap_Eienl = energy_module_Eienl
end subroutine f90wrap_energy_module__get__Eienl

subroutine f90wrap_energy_module__set__Eienl(f90wrap_Eienl)
    use energy_module, only: energy_module_Eienl => Eienl
    implicit none
    real(8), intent(in) :: f90wrap_Eienl
    
    energy_module_Eienl = f90wrap_Eienl
end subroutine f90wrap_energy_module__set__Eienl

subroutine f90wrap_energy_module__get__Eewald(f90wrap_Eewald)
    use energy_module, only: energy_module_Eewald => Eewald
    implicit none
    real(8), intent(out) :: f90wrap_Eewald
    
    f90wrap_Eewald = energy_module_Eewald
end subroutine f90wrap_energy_module__get__Eewald

subroutine f90wrap_energy_module__set__Eewald(f90wrap_Eewald)
    use energy_module, only: energy_module_Eewald => Eewald
    implicit none
    real(8), intent(in) :: f90wrap_Eewald
    
    energy_module_Eewald = f90wrap_Eewald
end subroutine f90wrap_energy_module__set__Eewald

subroutine f90wrap_energy_module__get__Efm(f90wrap_Efm)
    use energy_module, only: energy_module_Efm => Efm
    implicit none
    real(8), intent(out) :: f90wrap_Efm
    
    f90wrap_Efm = energy_module_Efm
end subroutine f90wrap_energy_module__get__Efm

subroutine f90wrap_energy_module__set__Efm(f90wrap_Efm)
    use energy_module, only: energy_module_Efm => Efm
    implicit none
    real(8), intent(in) :: f90wrap_Efm
    
    energy_module_Efm = f90wrap_Efm
end subroutine f90wrap_energy_module__set__Efm

subroutine f90wrap_energy_module__get__Tnad(f90wrap_Tnad)
    use energy_module, only: energy_module_Tnad => Tnad
    implicit none
    real(8), intent(out) :: f90wrap_Tnad
    
    f90wrap_Tnad = energy_module_Tnad
end subroutine f90wrap_energy_module__get__Tnad

subroutine f90wrap_energy_module__set__Tnad(f90wrap_Tnad)
    use energy_module, only: energy_module_Tnad => Tnad
    implicit none
    real(8), intent(in) :: f90wrap_Tnad
    
    energy_module_Tnad = f90wrap_Tnad
end subroutine f90wrap_energy_module__set__Tnad

subroutine f90wrap_energy_module__get__FE(f90wrap_FE)
    use energy_module, only: energy_module_FE => FE
    implicit none
    real(8), intent(out) :: f90wrap_FE
    
    f90wrap_FE = energy_module_FE
end subroutine f90wrap_energy_module__get__FE

subroutine f90wrap_energy_module__set__FE(f90wrap_FE)
    use energy_module, only: energy_module_FE => FE
    implicit none
    real(8), intent(in) :: f90wrap_FE
    
    energy_module_FE = f90wrap_FE
end subroutine f90wrap_energy_module__set__FE

subroutine f90wrap_energy_module__get__FE0(f90wrap_FE0)
    use energy_module, only: energy_module_FE0 => FE0
    implicit none
    real(8), intent(out) :: f90wrap_FE0
    
    f90wrap_FE0 = energy_module_FE0
end subroutine f90wrap_energy_module__get__FE0

subroutine f90wrap_energy_module__set__FE0(f90wrap_FE0)
    use energy_module, only: energy_module_FE0 => FE0
    implicit none
    real(8), intent(in) :: f90wrap_FE0
    
    energy_module_FE0 = f90wrap_FE0
end subroutine f90wrap_energy_module__set__FE0

subroutine f90wrap_energy_module__get__WmaxL(f90wrap_WmaxL)
    use energy_module, only: energy_module_WmaxL => WmaxL
    implicit none
    real(8), intent(out) :: f90wrap_WmaxL
    
    f90wrap_WmaxL = energy_module_WmaxL
end subroutine f90wrap_energy_module__get__WmaxL

subroutine f90wrap_energy_module__set__WmaxL(f90wrap_WmaxL)
    use energy_module, only: energy_module_WmaxL => WmaxL
    implicit none
    real(8), intent(in) :: f90wrap_WmaxL
    
    energy_module_WmaxL = f90wrap_WmaxL
end subroutine f90wrap_energy_module__set__WmaxL

subroutine f90wrap_energy_module__get__WminL(f90wrap_WminL)
    use energy_module, only: energy_module_WminL => WminL
    implicit none
    real(8), intent(out) :: f90wrap_WminL
    
    f90wrap_WminL = energy_module_WminL
end subroutine f90wrap_energy_module__get__WminL

subroutine f90wrap_energy_module__set__WminL(f90wrap_WminL)
    use energy_module, only: energy_module_WminL => WminL
    implicit none
    real(8), intent(in) :: f90wrap_WminL
    
    energy_module_WminL = f90wrap_WminL
end subroutine f90wrap_energy_module__set__WminL

! End of module energy_module defined in file Energy_module.fpp

