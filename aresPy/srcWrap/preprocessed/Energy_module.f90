!#############################################################!
!---               Calculate the Energies                  ---!
!---Author:              Qiang Xu                          ---!
!---Date:               2018/01/10                         ---!
!#############################################################!
!-----------------------------------------------------------
MODULE Energy_module
    USE constants
    USE struct_module, ONLY: Eii => Eionion

    USE smpi_math_module

    IMPLICIT NONE
    REAL(DP) :: Etot   & !total energy of system
           &   , Eband  & !Bands' energy
           &, Eh      & !Hartree
           &, Eext    & !external
           &, Exc     & !Exc
           &, Eele      & !electrstatic energy
           &, FE        & !Free-energy
           &, FE0         !Free-energy (T->0)
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    SUBROUTINE TotalEnergy(nps, eig, rhoS, rho)
        USE parameters, ONLY: nspin
        USE grid_module, ONLY: sph => grid, dvol, eigen_type
        USE xc_module, ONLY: xc_functional
        USE smearing_module, ONLY: wke, ets
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: nps !number of grid points
        REAL(DP), INTENT(IN) :: rhoS(nps, nspin), rho(nps)
        TYPE(eigen_type), INTENT(IN) :: eig
!LOCAL
        REAL(DP) :: Evxc, Evele

        REAL(DP) :: E_loc

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!Band Energy
        CALL Ebands(eig%val, wke, Eband)
!xc
        CALL xc_functional(nps, rhoS, sph%vxcS, Exc)

        E_loc = SUM(sph%vxcS*rhoS)*dvol
        CALL MPI_ALLREDUCE(E_loc, Evxc, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

        E_loc = 0.5_DP*SUM(rho*sph%vh)*dvol
        CALL MPI_ALLREDUCE(E_loc, Eh, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

!Eext

        E_loc = SUM(rho*sph%vlpp)*dvol
        CALL MPI_ALLREDUCE(E_loc, Eext, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

!Etot
        Etot = Eband - Evxc - Eh + Exc + Eii
!Free energy
        FE = Etot - ets
        FE0 = Etot - 0.5d0*ets
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE TotalEnergy
!---------------------For  band energy-----------------------
    SUBROUTINE Ebands(eval, wke, Eband)
        USE parameters, ONLY: nev => Nstates, nev_tot => Nstates_global, LMOM
!
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: eval(:, :, :)
        REAL(DP), INTENT(IN) :: wke(:, :, :)
        REAL(DP), INTENT(OUT) :: Eband
!LOCAL
        REAL(DP) :: WmaxL, WminL
        LOGICAL  :: Lcheck
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        Eband = SUM(eval(:, :, :)*wke(:, :, :))
!yhoo, Check out!!!

        WmaxL = MAXVAL(wke(nev_tot, :, :))
        WminL = MINVAL(wke(nev_tot, :, :))

        Lcheck = (ABS(WmaxL) > 1e-14) .OR. (ABS(WminL) > 1e-14)
        IF (Lcheck .AND. (.NOT. LMOM)) THEN
            WRITE (*, *) 'STOP: Please Add More Empty Bands.'
            STOP
        END IF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE Ebands
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END MODULE Energy_module
!-----------------------------------------------------------
!#############################################################!
!---           Calculate the force and stress              ---!
!---Author:              Qiang Xu                          ---!
!---Date:               2018/01/10                         ---!
!#############################################################!
!-----------------------------------------------------------
!#############################################################!
!---                    Output data                        ---!
!---Author:              Qiang Xu                          ---!
!---Date:               2018/01/10                         ---!
!#############################################################!
!-----------------------------------------------------------
MODULE Output_module
    USE constants
    USE Energy_module
    USE struct_module, ONLY: Eshift_ps, Eshift_tot

    USE smpi_math_module

    IMPLICIT NONE
    INTEGER(I4B) :: time_total0, time_total1 &
                 &, time_scf0, time_scf1
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    SUBROUTINE Output()
        USE parameters, ONLY: LCHARGE, LWAVE
        USE grid_module, ONLY: grid, Sphere2Cubic, n1, n2, n3
        USE struct_module, ONLY: natom, lat_mat, struct
        IMPLICIT NONE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        IF (parallel%isroot) THEN

!ground data
            PRINT *, '*******************OUTPUT DATA***********************'

            CALL write_time('[SCF Time]', .TRUE.)

            PRINT *, '-----------------------------------------------------'
            WRITE (*, 10) '*Band energy     ->', Eband*hart2ev, 'eV'
            WRITE (*, 10) '*XC energy       ->', Exc*hart2ev, 'eV'
            WRITE (*, 10) '*Hartree energy  ->', Eh*hart2ev, 'eV'
            WRITE (*, 10) '*I-E energy      ->', Eext*hart2ev, 'eV'
            WRITE (*, 10) '*Ion-Ion energy  ->', Eii*hart2ev, 'eV'
            PRINT *, '-----------------------------------------------------'
            WRITE (*, 10) '*Total-energy    ->', Etot*hart2eV, 'eV'
            WRITE (*, 10) '*Free-energy     ->', FE*hart2eV, 'eV'
            WRITE (*, 10) '*0K Free-energy  ->', FE0*hart2eV, 'eV'
            IF (ABS(Eshift_ps) > 1e-6) THEN
                PRINT *, '-----------------------------------------------------'
                WRITE (*, 10) '*Cohesive energy ->', (Etot + Eshift_ps)*hart2ev, 'eV'
                WRITE (*, 10) '*AE energy (test)->', (Etot + Eshift_tot)*hart2ev, 'eV'
            END IF

        END IF

!Write File
!Output charge density
        IF (LCHARGE) CALL write_density()
!Output wave functions
!IF(LWAVE) CALL write_wvf()
!Output band structures
        CALL write_band()

        IF (parallel%isroot) THEN

            PRINT *, '*******************OUTPUT DATA***********************'

        END IF

10      FORMAT(A20, 1X, ' ', 1X, SP, F23.7, A3)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE Output
!------------------------------------------------------------------------
    SUBROUTINE write_density()
        USE parameters, ONLY: nspin
        USE grid_module, ONLY: ni1 => global_n1, ni2 => global_n2, ni3 => global_n3 &
                    & , Sphere2Cubic, sph => grid, nps => n
        IMPLICIT NONE
        REAL(DP) :: rhot(ni1, ni2, ni3, 2)
        INTEGER(I4B) :: isp
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (nspin == 1) THEN
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1), rhot(:, :, :, 1))
            rhot(:, :, :, 2) = 0._DP
        ELSE
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1) + sph%rhoS(:, 2), rhot(:, :, :, 1))
            CALL Sphere2Cubic(nps, sph%rhoS(:, 1) - sph%rhoS(:, 2), rhot(:, :, :, 2))
        END IF

        IF (parallel%isroot) THEN

            OPEN (298, FILE='ares.chg')
            DO isp = 1, 2
                WRITE (298, *) ni1, ni2, ni3
                WRITE (298, *) rhot(1:ni1, 1:ni2, 1:ni3, isp)
            END DO
            CLOSE (298)

        END IF

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE write_density
!------------------------------------------------------------------------

!------------------------------------------------------------------------
    SUBROUTINE write_band()
        USE parameters, ONLY: nev => Nstates, nspin, MOMsigma
        USE smearing_module, ONLY: wke, ets
        USE grid_module, ONLY: eigen
        IMPLICIT NONE
        INTEGER(I4B) :: Ii, Is
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        IF (parallel%isroot) THEN

            OPEN (300, FILE='ares.band')
            WRITE (300, "(1X,I6,2X,I2,2X,F8.5,4X,A44)") nev, nspin, MOMsigma, '# of states, nspin, Width of MOM smearing'
            WRITE (300, *) '=i-state==i-spin==Occupation==Orbital energy (eV)='
            DO Ii = 1, nev
                DO Is = 1, nspin
                    WRITE (300, "(1X,I6,4X,I2,6X,F10.8,2X,F15.9)") &
                     &  Ii, Is, wke(Ii, 1, Is), eigen%val(Ii, 1, Is)*hart2eV
                END DO
            END DO
            CLOSE (300)

        END IF

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE write_band
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END MODULE Output_module
