! Module constants defined in file Constants.f90

subroutine f90wrap_constants__get__I4B(f90wrap_I4B)
    use constants, only: constants_I4B => I4B
    implicit none
    integer, intent(out) :: f90wrap_I4B
    
    f90wrap_I4B = constants_I4B
end subroutine f90wrap_constants__get__I4B

subroutine f90wrap_constants__get__I2B(f90wrap_I2B)
    use constants, only: constants_I2B => I2B
    implicit none
    integer, intent(out) :: f90wrap_I2B
    
    f90wrap_I2B = constants_I2B
end subroutine f90wrap_constants__get__I2B

subroutine f90wrap_constants__get__I1B(f90wrap_I1B)
    use constants, only: constants_I1B => I1B
    implicit none
    integer, intent(out) :: f90wrap_I1B
    
    f90wrap_I1B = constants_I1B
end subroutine f90wrap_constants__get__I1B

subroutine f90wrap_constants__get__SP(f90wrap_SP)
    use constants, only: constants_SP => SP
    implicit none
    integer, intent(out) :: f90wrap_SP
    
    f90wrap_SP = constants_SP
end subroutine f90wrap_constants__get__SP

subroutine f90wrap_constants__get__DP(f90wrap_DP)
    use constants, only: constants_DP => DP
    implicit none
    integer, intent(out) :: f90wrap_DP
    
    f90wrap_DP = constants_DP
end subroutine f90wrap_constants__get__DP

subroutine f90wrap_constants__get__SCP(f90wrap_SCP)
    use constants, only: constants_SCP => SCP
    implicit none
    integer, intent(out) :: f90wrap_SCP
    
    f90wrap_SCP = constants_SCP
end subroutine f90wrap_constants__get__SCP

subroutine f90wrap_constants__get__DCP(f90wrap_DCP)
    use constants, only: constants_DCP => DCP
    implicit none
    integer, intent(out) :: f90wrap_DCP
    
    f90wrap_DCP = constants_DCP
end subroutine f90wrap_constants__get__DCP

subroutine f90wrap_constants__get__LGT(f90wrap_LGT)
    use constants, only: constants_LGT => LGT
    implicit none
    integer, intent(out) :: f90wrap_LGT
    
    f90wrap_LGT = constants_LGT
end subroutine f90wrap_constants__get__LGT

subroutine f90wrap_constants__get__pi(f90wrap_pi)
    use constants, only: constants_pi => pi
    implicit none
    real(8), intent(out) :: f90wrap_pi
    
    f90wrap_pi = constants_pi
end subroutine f90wrap_constants__get__pi

subroutine f90wrap_constants__get__pio2(f90wrap_pio2)
    use constants, only: constants_pio2 => pio2
    implicit none
    real(8), intent(out) :: f90wrap_pio2
    
    f90wrap_pio2 = constants_pio2
end subroutine f90wrap_constants__get__pio2

subroutine f90wrap_constants__get__twopi(f90wrap_twopi)
    use constants, only: constants_twopi => twopi
    implicit none
    real(8), intent(out) :: f90wrap_twopi
    
    f90wrap_twopi = constants_twopi
end subroutine f90wrap_constants__get__twopi

subroutine f90wrap_constants__get__imag(f90wrap_imag)
    use constants, only: constants_imag => imag
    implicit none
    complex(8), intent(out) :: f90wrap_imag
    
    f90wrap_imag = constants_imag
end subroutine f90wrap_constants__get__imag

subroutine f90wrap_constants__get__inputunit(f90wrap_inputunit)
    use constants, only: constants_inputunit => inputunit
    implicit none
    integer(8), intent(out) :: f90wrap_inputunit
    
    f90wrap_inputunit = constants_inputunit
end subroutine f90wrap_constants__get__inputunit

subroutine f90wrap_constants__get__errunit(f90wrap_errunit)
    use constants, only: constants_errunit => errunit
    implicit none
    integer(8), intent(out) :: f90wrap_errunit
    
    f90wrap_errunit = constants_errunit
end subroutine f90wrap_constants__get__errunit

subroutine f90wrap_constants__get__outputunit(f90wrap_outputunit)
    use constants, only: constants_outputunit => outputunit
    implicit none
    integer(8), intent(out) :: f90wrap_outputunit
    
    f90wrap_outputunit = constants_outputunit
end subroutine f90wrap_constants__get__outputunit

subroutine f90wrap_constants__get__mdposunit(f90wrap_mdposunit)
    use constants, only: constants_mdposunit => mdposunit
    implicit none
    integer(8), intent(out) :: f90wrap_mdposunit
    
    f90wrap_mdposunit = constants_mdposunit
end subroutine f90wrap_constants__get__mdposunit

subroutine f90wrap_constants__get__rydberg(f90wrap_rydberg)
    use constants, only: constants_rydberg => rydberg
    implicit none
    real(8), intent(out) :: f90wrap_rydberg
    
    f90wrap_rydberg = constants_rydberg
end subroutine f90wrap_constants__get__rydberg

subroutine f90wrap_constants__get__bohr2ang(f90wrap_bohr2ang)
    use constants, only: constants_bohr2ang => bohr2ang
    implicit none
    real(8), intent(out) :: f90wrap_bohr2ang
    
    f90wrap_bohr2ang = constants_bohr2ang
end subroutine f90wrap_constants__get__bohr2ang

subroutine f90wrap_constants__get__ang2bohr(f90wrap_ang2bohr)
    use constants, only: constants_ang2bohr => ang2bohr
    implicit none
    real(8), intent(out) :: f90wrap_ang2bohr
    
    f90wrap_ang2bohr = constants_ang2bohr
end subroutine f90wrap_constants__get__ang2bohr

subroutine f90wrap_constants__get__hart2ev(f90wrap_hart2ev)
    use constants, only: constants_hart2ev => hart2ev
    implicit none
    real(8), intent(out) :: f90wrap_hart2ev
    
    f90wrap_hart2ev = constants_hart2ev
end subroutine f90wrap_constants__get__hart2ev

subroutine f90wrap_constants__get__ev2hart(f90wrap_ev2hart)
    use constants, only: constants_ev2hart => ev2hart
    implicit none
    real(8), intent(out) :: f90wrap_ev2hart
    
    f90wrap_ev2hart = constants_ev2hart
end subroutine f90wrap_constants__get__ev2hart

subroutine f90wrap_constants__get__scf_tol(f90wrap_scf_tol)
    use constants, only: constants_scf_tol => scf_tol
    implicit none
    real(8), intent(out) :: f90wrap_scf_tol
    
    f90wrap_scf_tol = constants_scf_tol
end subroutine f90wrap_constants__get__scf_tol

subroutine f90wrap_constants__get__au2ev_force(f90wrap_au2ev_force)
    use constants, only: constants_au2ev_force => au2ev_force
    implicit none
    real(8), intent(out) :: f90wrap_au2ev_force
    
    f90wrap_au2ev_force = constants_au2ev_force
end subroutine f90wrap_constants__get__au2ev_force

subroutine f90wrap_constants__get__au2gpa(f90wrap_au2gpa)
    use constants, only: constants_au2gpa => au2gpa
    implicit none
    real(8), intent(out) :: f90wrap_au2gpa
    
    f90wrap_au2gpa = constants_au2gpa
end subroutine f90wrap_constants__get__au2gpa

subroutine f90wrap_constants__get__golden(f90wrap_golden)
    use constants, only: constants_golden => golden
    implicit none
    real(8), intent(out) :: f90wrap_golden
    
    f90wrap_golden = constants_golden
end subroutine f90wrap_constants__get__golden

subroutine f90wrap_constants__get__vlight(f90wrap_vlight)
    use constants, only: constants_vlight => vlight
    implicit none
    real(8), intent(out) :: f90wrap_vlight
    
    f90wrap_vlight = constants_vlight
end subroutine f90wrap_constants__get__vlight

subroutine f90wrap_constants__get__CONST_ME_AU(f90wrap_CONST_ME_AU)
    use constants, only: constants_CONST_ME_AU => CONST_ME_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_ME_AU
    
    f90wrap_CONST_ME_AU = constants_CONST_ME_AU
end subroutine f90wrap_constants__get__CONST_ME_AU

subroutine f90wrap_constants__get__CONST_E_AU(f90wrap_CONST_E_AU)
    use constants, only: constants_CONST_E_AU => CONST_E_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_E_AU
    
    f90wrap_CONST_E_AU = constants_CONST_E_AU
end subroutine f90wrap_constants__get__CONST_E_AU

subroutine f90wrap_constants__get__CONST_EH_AU(f90wrap_CONST_EH_AU)
    use constants, only: constants_CONST_EH_AU => CONST_EH_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_EH_AU
    
    f90wrap_CONST_EH_AU = constants_CONST_EH_AU
end subroutine f90wrap_constants__get__CONST_EH_AU

subroutine f90wrap_constants__get__CONST_LEN_AU(f90wrap_CONST_LEN_AU)
    use constants, only: constants_CONST_LEN_AU => CONST_LEN_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_LEN_AU
    
    f90wrap_CONST_LEN_AU = constants_CONST_LEN_AU
end subroutine f90wrap_constants__get__CONST_LEN_AU

subroutine f90wrap_constants__get__CONST_HBAR_AU(f90wrap_CONST_HBAR_AU)
    use constants, only: constants_CONST_HBAR_AU => CONST_HBAR_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_HBAR_AU
    
    f90wrap_CONST_HBAR_AU = constants_CONST_HBAR_AU
end subroutine f90wrap_constants__get__CONST_HBAR_AU

subroutine f90wrap_constants__get__CONST_EP_AU(f90wrap_CONST_EP_AU)
    use constants, only: constants_CONST_EP_AU => CONST_EP_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_EP_AU
    
    f90wrap_CONST_EP_AU = constants_CONST_EP_AU
end subroutine f90wrap_constants__get__CONST_EP_AU

subroutine f90wrap_constants__get__CONST_F_AU(f90wrap_CONST_F_AU)
    use constants, only: constants_CONST_F_AU => CONST_F_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_F_AU
    
    f90wrap_CONST_F_AU = constants_CONST_F_AU
end subroutine f90wrap_constants__get__CONST_F_AU

subroutine f90wrap_constants__get__CONST_MT_AU(f90wrap_CONST_MT_AU)
    use constants, only: constants_CONST_MT_AU => CONST_MT_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_MT_AU
    
    f90wrap_CONST_MT_AU = constants_CONST_MT_AU
end subroutine f90wrap_constants__get__CONST_MT_AU

subroutine f90wrap_constants__get__CONST_TIME_AU(f90wrap_CONST_TIME_AU)
    use constants, only: constants_CONST_TIME_AU => CONST_TIME_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_TIME_AU
    
    f90wrap_CONST_TIME_AU = constants_CONST_TIME_AU
end subroutine f90wrap_constants__get__CONST_TIME_AU

subroutine f90wrap_constants__get__CONST_I_AU(f90wrap_CONST_I_AU)
    use constants, only: constants_CONST_I_AU => CONST_I_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_I_AU
    
    f90wrap_CONST_I_AU = constants_CONST_I_AU
end subroutine f90wrap_constants__get__CONST_I_AU

subroutine f90wrap_constants__get__CONST_TEMP_AU(f90wrap_CONST_TEMP_AU)
    use constants, only: constants_CONST_TEMP_AU => CONST_TEMP_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_TEMP_AU
    
    f90wrap_CONST_TEMP_AU = constants_CONST_TEMP_AU
end subroutine f90wrap_constants__get__CONST_TEMP_AU

subroutine f90wrap_constants__get__CONST_P_AU(f90wrap_CONST_P_AU)
    use constants, only: constants_CONST_P_AU => CONST_P_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_P_AU
    
    f90wrap_CONST_P_AU = constants_CONST_P_AU
end subroutine f90wrap_constants__get__CONST_P_AU

subroutine f90wrap_constants__get__CONST_V_AU(f90wrap_CONST_V_AU)
    use constants, only: constants_CONST_V_AU => CONST_V_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_V_AU
    
    f90wrap_CONST_V_AU = constants_CONST_V_AU
end subroutine f90wrap_constants__get__CONST_V_AU

subroutine f90wrap_constants__get__CONST_KE_AU(f90wrap_CONST_KE_AU)
    use constants, only: constants_CONST_KE_AU => CONST_KE_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_KE_AU
    
    f90wrap_CONST_KE_AU = constants_CONST_KE_AU
end subroutine f90wrap_constants__get__CONST_KE_AU

subroutine f90wrap_constants__get__CONST_MU_AU(f90wrap_CONST_MU_AU)
    use constants, only: constants_CONST_MU_AU => CONST_MU_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_MU_AU
    
    f90wrap_CONST_MU_AU = constants_CONST_MU_AU
end subroutine f90wrap_constants__get__CONST_MU_AU

subroutine f90wrap_constants__get__CONST_MA_AU(f90wrap_CONST_MA_AU)
    use constants, only: constants_CONST_MA_AU => CONST_MA_AU
    implicit none
    real(8), intent(out) :: f90wrap_CONST_MA_AU
    
    f90wrap_CONST_MA_AU = constants_CONST_MA_AU
end subroutine f90wrap_constants__get__CONST_MA_AU

subroutine f90wrap_constants__get__CONST_KB_SI(f90wrap_CONST_KB_SI)
    use constants, only: constants_CONST_KB_SI => CONST_KB_SI
    implicit none
    real(8), intent(out) :: f90wrap_CONST_KB_SI
    
    f90wrap_CONST_KB_SI = constants_CONST_KB_SI
end subroutine f90wrap_constants__get__CONST_KB_SI

subroutine f90wrap_constants__get__angs(f90wrap_angs)
    use constants, only: constants_angs => angs
    implicit none
    real(8), intent(out) :: f90wrap_angs
    
    f90wrap_angs = constants_angs
end subroutine f90wrap_constants__get__angs

subroutine f90wrap_constants__get__force2ev(f90wrap_force2ev)
    use constants, only: constants_force2ev => force2ev
    implicit none
    real(8), intent(out) :: f90wrap_force2ev
    
    f90wrap_force2ev = constants_force2ev
end subroutine f90wrap_constants__get__force2ev

subroutine f90wrap_constants__get__CLen(f90wrap_CLen)
    use constants, only: constants_CLen => CLen
    implicit none
    integer, intent(out) :: f90wrap_CLen
    
    f90wrap_CLen = constants_CLen
end subroutine f90wrap_constants__get__CLen

subroutine f90wrap_constants__get__xtiny(f90wrap_xtiny)
    use constants, only: constants_xtiny => xtiny
    implicit none
    real(8), intent(out) :: f90wrap_xtiny
    
    f90wrap_xtiny = constants_xtiny
end subroutine f90wrap_constants__get__xtiny

subroutine f90wrap_constants__get__autogpa(f90wrap_autogpa)
    use constants, only: constants_autogpa => autogpa
    implicit none
    real(8), intent(out) :: f90wrap_autogpa
    
    f90wrap_autogpa = constants_autogpa
end subroutine f90wrap_constants__get__autogpa

subroutine f90wrap_constants__get__hartree2ev(f90wrap_hartree2ev)
    use constants, only: constants_hartree2ev => hartree2ev
    implicit none
    real(8), intent(out) :: f90wrap_hartree2ev
    
    f90wrap_hartree2ev = constants_hartree2ev
end subroutine f90wrap_constants__get__hartree2ev

! End of module constants defined in file Constants.f90

