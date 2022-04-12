# 1 "IonLocalPotentialAssignment.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "IonLocalPotentialAssignment.f90"
MODULE GetVLocalPseudoPotential
    USE constants
    USE parameters, ONLY: NDG, LDG, inpol
    USE grid_module, ONLY: grid, n1, n2, n3, gap
    USE pspot_module
    USE nlpot_module, ONLY: nlp_wijk_dg

    USE smpi_math_module, ONLY: parallel, mpinfo

    IMPLICIT NONE
    !==================================
    !##DO NOT NEED, JUST FOR COMPILED
    INTEGER(I4B)            :: RadiusN
    REAL(DP), ALLOCATABLE    :: fr(:)
    !==================================
    !weight !! for DoubleGrid
    REAL(DP), ALLOCATABLE :: wijk(:), drijk(:, :)
    INTEGER(I4B)         :: Ilft, Irit, numdg
    
    private CubicSplineInterp
CONTAINS
    !-----------------------DIVIDER-LINE--------------------------
    SUBROUTINE CalVlpp()
        USE struct_module, ONLY: struct, naty, lat_mat
        !##xlt
        USE grid_module, ONLY: rho_calc, n
        USE m_time_evaluate, ONLY: memory_sum, memory_free
        IMPLICIT NONE
        INTEGER(I4B)         :: i, j




        REAL(DP)             :: temp(parallel%mygrid_range(3)), Vsphere(parallel%mygrid_range(3), 1) &
             &, Vtemp(parallel%mygrid_range(3), 1)

        call memory_sum('Vlpp', (real(size(temp), DP) + size(Vsphere) + size(Vtemp))*DP)
        !> Double Grid
        ! REAL(DP),allocatable :: Vcomp(:)
        ! INTEGER(I4B) :: n_Vcomp
        !==========================================================
        IF (LDG) THEN
            !> up bound and down bound
            Ilft = -(ndg + 1)*((inpol + 2)/2) !1-2*(ndg+1)
            Irit = (ndg + 1)*((inpol + 1)/2) !2*(ndg+1)-1
            ! Ilft=-ndg-1
            ! Irit=ndg+1
            !> number of dense-grid
            ! numdg=(n_near*2+1)**3
            numdg = ((inpol + 1)*(ndg + 1) + 1)**3
            ! numdg=(n_near*2+1+2*n_near)**3
            !numdg=(Irit-Ilft+1)**3
            ALLOCATE (wijk(numdg), drijk(3, numdg))
            call memory_sum('double grid', (real(size(wijk), DP) + size(drijk))*DP)
            CALL nlp_wijk_dg(NDG, inpol, numdg, Ilft, Irit, wijk, drijk)
        END IF



        rho_calc%vlpp = 0.d0

        DO i = 1, naty, 1
            !=============================
            !##OLD PSEDUPOTENTIAL CALCULATE
            !CALL CoordinateTranlate(i)
            DO j = 1, struct%nati(i), 1
                IF (LDG) THEN
                    ! n_Vcomp=size(psp(i)%r_real)
                    ! allocate(Vcomp(n_Vcomp))
                    ! CALL set_vcomp(n_Vcomp,psp(i)%Zion,0.5d0,psp(i)%r_real,psp(i)%V_loc,Vcomp)
                    ! CALL IonPotentialAssignment_Dg(i,n_Vcomp,Vcomp,psp(i)%Zion,struct%poscar(:,struct%eleid(i)+j-1),temp)
                    CALL IonPotentialAssignment_Dg(i, psp(i)%Zion, struct%poscar(:, struct%eleid(i) + j - 1), temp)
                    ! deallocate(Vcomp)
                else
                    CALL IonPotentialAssignment(i, psp(i)%Zion, struct%poscar(:, struct%eleid(i) + j - 1), temp)
                end if



                rho_calc%vlpp = rho_calc%vlpp + temp

                !========================
                !##PLOT SAMPLE
                !open(unit=111, file="fr")
                !write(unit=111, fmt=*)fr
                !close(111)
                !========================
                !##PRINT INFO
                !print*,"poscar",struct%poscar(:,struct%eleid(i)+j-1)
                !print*,"atomid",struct%eleid(i)+j-1
                !print*,"gap",gap
                !print*,"lat_mat",lat_mat
                !========================
            END DO
        END DO
        !> clean useless array
        IF (LDG) THEN
            call memory_free('double grid', (real(size(wijk), DP) + size(drijk))*DP)
            deallocate (wijk, drijk)
        END IF
        call memory_free('Vlpp', (real(size(temp), DP) + size(Vsphere) + size(Vtemp))*DP)

        !#ifdef DEBUG
        ! #ifdef 1
        !     if(parallel%isroot)then
        !     open(1212,file='Vl_debug_p')
        !     write(1212,*)rho_calc%vlpp(1:1000)
        !     close(1212)
        !     endif
        ! #else
        !     DO i=1,rho_calc%OneDLength
        !        out_sphere(i)=grid%vlpp(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))
        !     ENDDO
        !     open(1212,file='Vl_debug_s')
        !     write(1212,*)out_sphere(1:1000)
        !     close(1212)
        ! #endif
        !#endif
    END SUBROUTINE CalVlpp
    !-----------------------PARTING-LINE--------------------------
    SUBROUTINE IonPotentialAssignment(Ity, zion, poscar, temp)
        USE constants
        USE struct_module, ONLY: struct, lat_mat
        USE MathSplines, ONLY: polynom
        USE grid_module, ONLY: rho_calc
        USE pspot_module, ONLY: max_rcut
        IMPLICIT NONE



        REAL(DP), INTENT(OUT)       :: temp(parallel%mygrid_range(3))

        INTEGER(I4B), INTENT(IN)    :: Ity  !element number
        REAL(DP), INTENT(IN)        :: Zion
        REAL(DP)                   :: rr, poscar(3), dfr(RadiusN + 1), CarPos(3), DirPos(3)
        INTEGER(I4B)               :: i, j, k, NInOneD, fdim, m, n = 0
        REAL(DP)                   :: c(3)
        INTEGER(I4B)               :: max_m
        !========================================================================
        REAL(DP), ALLOCATABLE :: fv(:), r_fv(:)!,dfdx_fv(fdim),dr=0.d0,r_fv(fdim)!,fv_i(40001),df2dx2_fv(10001)
        !==========================================================xlt test
        REAL(DP) :: out_sphere(rho_calc%OneDLength)
        INTEGER(I4B) :: aa

        !>
        c = 0.d0
        !##INITIALIZE CALCULATE
        ! open(998,file='Vloc_r')
        ! open(999,file='Vloc_f')



        temp = 0.d0

        max_m = size(psp(Ity)%r_real)
# 166 "IonLocalPotentialAssignment.f90"
                    DO i = 1, parallel%mygrid_range(3), 1
                        ! rr=grid%rvec(4,rho_calc%n(i))
                        rr = sqrt((grid%rvec(1, rho_calc%n(i)) - poscar(1))**2 + &
                             & (grid%rvec(2, rho_calc%n(i)) - poscar(2))**2 + &
                             & (grid%rvec(3, rho_calc%n(i)) - poscar(3))**2)

                        !## OLDINTERPOLATING
                        !##INTERPOLATING BY THREE POINT
                        !##FIND THE THREE POINT
                        m = 1
                        DO WHILE ((psp(Ity)%r_real(m) .lt. rr) .and. (m .le. max_m))
                            !DO WHILE ( r_fv(m+1).lt.rr )
                            m = m + 1
                        END DO
                        !##CALCULATE POTENTIAL IN "RR"







                        if (m .le. max_m) then
                            temp(i) = polynom(0, 3, psp(Ity)%r_real(m - 1:m + 1), psp(Ity)%V_loc(m - 1:m + 1), c, rr)
                        else
                            temp(i) = -psp(Ity)%Zion/rr
                        end if

                        !temp(i,j,k)=polynom(0,3,r_fv(m-1:m+1),fv(m-1:m+1),c,rr)
                        !##PLOT
                        ! IF(rr< max_rcut)THEN
                        !  write(998,"(F20.15)")rr
                        !  write(999,"(F20.15)")temp(i)
                        ! ENDIF
                        !=================================================================





        end do     !> i=1,parallel%mygrid_range(3),1

        !=======================================================================
        ! close(998)
        ! close(999)
        !============================
        !##TEST
        ! #ifdef 1
        !     if(parallel%isroot)then
        !     open(1212,file='Vl_debug_p')
        !     write(1212,*)temp(1:1000)
        !     ! write(1212,*)psp(Ity)%r_real
        !     close(1212)
        !     endif
        ! #else
        !     DO i=1,rho_calc%OneDLength
        !        out_sphere(i)=temp(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))
        !     ENDDO
        !     open(1212,file='Vl_debug_s')
        !     write(1212,*)out_sphere(1:1000)
        !     ! write(1212,*)psp(Ity)%r_real
        !     close(1212)
        ! #endif
        !      stop "ion step"
        !============================
    END SUBROUTINE IonPotentialAssignment
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    FUNCTION sphbess(l, x)
        ! Arguments
        integer, intent(in) :: l
        real(DP), intent(in) :: x
        REAL(DP) :: sphbess
        ! Local variables
        integer :: j
        real(kind=DP), parameter :: third = 1.0_DP/3.0_DP
        real(kind=DP), parameter :: ftnth = 1.0_DP/14.0_DP
        real(kind=DP) :: x2, sb0, sb1, by, bym, byp, ux
        if (abs(x) > 0.001_DP) then
            sb0 = sin(x)/x
        else
            x2 = 0.5_DP*x*x
            sb0 = 1.0_DP - third*x2*(1.0_DP - 0.1_DP*x2)
        end if
        if (l == 0) then
            sphbess = sb0
        else
            if (abs(x) > 0.001_DP) then
                sb1 = (sb0 - cos(x))/x
            else
                sb1 = third*x*(1.0_DP - (0.2_DP*x2)*(1.0_DP - ftnth*x2))
            end if
            if (l == 1) then
                sphbess = sb1
            else if (x == 0.0_DP) then
                sphbess = 0.0_DP
            else
                by = sb1
                bym = sb0
                ux = 1.0_DP/x
                do j = 1, l - 1
                    byp = real(2*J + 1, DP)*ux*by - bym
                    bym = by
                    by = byp
                end do
                sphbess = by
            end if
        end if
    end function sphbess
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    SUBROUTINE FourBess_gr(g, fg, r, fr)
        IMPLICIT NONE
        REAL(DP), INTENT(IN)  :: r(:), g(:), fg(:)
        REAL(DP), INTENT(OUT) :: fr(:)
        INTEGER(I4B) :: Ng, Nr
        INTEGER(I4B) :: Ig, Ir
        REAL(DP) :: y, dg
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        Nr = SIZE(r, 1)
        Ng = SIZE(g, 1)
        DO Ir = 1, Nr
            dg = g(2) - g(1)
            y = fg(1)*g(1)**2*sphbess(0, g(1)*r(Ir))*dg
            DO Ig = 2, Ng - 1
                dg = 0.5d0*(g(Ig + 1) - g(Ig - 1))
                y = y + fg(Ig)*g(Ig)**2*sphbess(0, g(Ig)*r(Ir))*dg
            END DO
            dg = g(Ng) - g(Ng - 1)
            y = y + fg(Ng)*g(Ng)**2*sphbess(0, g(Ng)*r(Ir))*dg
            !store
            fr(Ir) = y
        END DO
        fr(:) = fr(:)/(2.d0*pi**2)
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE FourBess_gr
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    FUNCTION CubicSplineInterp(fun, ddfdx2, xmax, dx, x, Zion)
        !just for uniform grid now , x>0
        !IN/OUT
        REAL(DP), INTENT(IN) :: fun(:)  &
             &, ddfdx2(:) &
             &, xmax, dx, x
        INTEGER(I4B), OPTIONAL, INTENT(IN)  :: Zion
        REAL(DP) :: CubicSplineInterp
        !LOCAL
        REAL(DP) :: yval  &
             &, ypp_left &
             &, ypp_right &
             &, y_left  &
             &, y_right &
             &, pos &
             &, dt
        INTEGER(I4B) :: left  &
             &, right
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (x < 0.d0) THEN
            WRITE (6, *) 'CubicSplineInterp: r must > 0'
            STOP
        END IF
        !Some special case
        IF (x == 0.d0) THEN
            CubicSplineInterp = fun(1)
            RETURN
        END IF
        !For r > rmax , fun=0.d0
        IF (x > xmax) THEN
            CubicSplineInterp = -Zion/x
            RETURN
        END IF
        !vlpp
        IF (PRESENT(Zion)) THEN
            ! read the following words to see why I do this check here
            ! since NaN is not equal to anyting in FORTRAN
            ! we use the following to see if 4*pi/qNorm is too big to case a NaN
            IF (-4._DP*pi/x**2_DP .NE. -4_DP*pi/x**2_DP) THEN
                WRITE (6, *) &
                    ' There is another case which should be considered: very large bulk. ', &
                    '  because larger system will give denser q points and will put q points ', &
                    ' very closer to q=0 and this will make -4*pi*Z/q**2 to -infinity', &
                    ' If computer find -4*pi*Z/q**2 is too big, it will generate NaN', &
                    ' But currently, in our group, we have not meet a system large enough to ', &
                    ' cause CPU to generate NaN, if you meet such kind of problem, do something!, code STOP!!'
                STOP
            END IF
        END IF
        !interpolation:
        !1.Find the left and right point , dt
        pos = x/dx !+ 1.d0
        !left
        left = FLOOR(pos)
        !right
        right = left + 1
        !dt
        dt = pos - REAL(left, DP)
        !2.Evaulate the polynomial
        y_left = fun(left)
        y_right = fun(right)
        ypp_left = ddfdx2(left)
        ypp_right = ddfdx2(right)
        !just do it
        yval = y_left + dt*(y_right - y_left &
                            - (ypp_right/6.0_DP + ypp_left/3.0_DP) &
                            + dt*(0.5_DP*ypp_left &
                                  + dt*((ypp_right - ypp_left)/6.0_DP)))
        !output
        !#IF(PRESENT(Zion))THEN
        !#    CubicSplineInterp=yval - 4.d0*pi*REAL(Zion,DP) / x**2._DP
        !#ELSE
        CubicSplineInterp = yval
        !#ENDIF
        RETURN
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END FUNCTION CubicSplineInterp
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    FUNCTION CubicHermiteInterp(fun, dfdx, xmax, h, x, zion)
        !##############################################!
        !y=A0*y0+A1*y1+B0*dy0+B1*dy1                   !
        !defined : dx0=x-x0 and dx1=x-x1 and h=x1-x0   !
        !A0=( 1+2*dx0/h )*(dx1/h)**2                   !
        !A1=( 1-2*dx1/h )*(dx0/h)**2                   !
        !B0=dx0*(dx1/h)**2                             !
        !B1=dx1*(dx0/h)**2                             !
        !##############################################!
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: fun(:)   &  !in fun
             &, dfdx(:)  &  !first derivertive
             &, xmax     &  !max x
             &, h        &  !grid size
             &, x           !in x
        REAL(DP) :: CubicHermiteInterp
        !LOCAL
        REAL(DP) :: pos, y0, y1, dy0, dy1
        REAL(DP) :: dx0, dx1, dx0_h, dx1_h, dx0_h2, dx1_h2
        INTEGER(I4B) :: left, right, zion
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (x < 0.d0) THEN
            WRITE (6, *) 'CubicSplineInterp: r must > 0'
            STOP
        END IF
        !Some special case
        IF (x == 0.d0) THEN
            CubicHermiteInterp = fun(1)
            RETURN
        END IF
        !For r > rmax , fun=0.d0
        IF (x > xmax) THEN
            CubicHermiteInterp = -zion/x
            RETURN
        END IF
        IF (-4._DP*pi/x**2_DP .NE. -4_DP*pi/x**2_DP) THEN
            WRITE (6, *) &
                ' There is another case which should be considered: very large bulk. ', &
                '  because larger system will give denser q points and will put q points ', &
                ' very closer to q=0 and this will make -4*pi*Z/q**2 to -infinity', &
                ' If computer find -4*pi*Z/q**2 is too big, it will generate NaN', &
                ' But currently, in our group, we have not meet a system large enough to ', &
                ' cause CPU to generate NaN, if you meet such kind of problem, do something!, code STOP!!'
            STOP
        END IF
        !interpolation:
        !1.Find the left and right point , dt
        pos = x/h + 1.d0
        !left
        left = FLOOR(pos)
        !right
        right = left + 1
        !dxi
        dx0 = x - (left - 1)*h  !x-x0
        dx1 = dx0 - h              !x-x1
        !dxi/h
        dx0_h = dx0/h
        dx1_h = dx1/h
        !(dxi/h)**2
        dx0_h2 = dx0_h**2
        dx1_h2 = dx1_h**2
        !2.interpolate
        y0 = fun(left)
        y1 = fun(right)
        dy0 = dfdx(left)
        dy1 = dfdx(right)
        !
        CubicHermiteInterp = (y0*(1.d0 + 2.d0*dx0_h) + dy0*dx0)*dx1_h2 &
             & + (y1*(1.d0 - 2.d0*dx1_h) + dy1*dx1)*dx0_h2
        !IF(x>1.56.AND.X<1.60)print*,"right,y1,dy1",right,y1,dy1  DEBUG array overstep
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END FUNCTION CubicHermiteInterp
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    SUBROUTINE dfdr(np, h, f, zion, df)
        !just for first order of the don-dimension fun,r : uniform grid
        USE parameters, ONLY: finite_order
        IMPLICIT NONE
        !INOUT
        INTEGER(I4B), INTENT(IN) :: np, zion
        REAL(DP), INTENT(IN)  :: f(np)  &   !f
             & , h  !grid size
        REAL(DP), INTENT(OUT) :: df(np)  !df_dr
        !LOCAL
        REAL(DP) :: coe(-finite_order:finite_order)
        REAL(DP) :: ft(-finite_order:np + finite_order), tmp
        INTEGER(I4B) :: i, ish
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ft(1:np) = f(1:np)
        !get coe
        CALL finite_factor(1, finite_order, coe)
        coe(:) = coe(:)/h
        !boundary set
        DO i = -finite_order, 0
            ft(i) = f(1)
        END DO
        ft(np + 1:np + finite_order) = zion/((/(i, i=np + 1, np + finite_order, 1)/)*h)
        !finite difference
        DO i = 1, np
            tmp = 0.d0
            DO ish = -finite_order, finite_order, 1
                tmp = tmp + coe(ish)*ft(ish + i)
            END DO
            df(i) = tmp
        END DO
        !first few point
        !df(1:finite_order-1)=0.d0
        !Last few points
        !df(np-finite_order+1:np)=0.d0
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE dfdr
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    subroutine finite_factor(fnor, norder, coe)!{{{
        !##########################################################
        !* CREATED_TIME  : 2013-03-12
        !* AUTHOR        : Yanchao Wang
        !* CHANGE        : Xuecheng Shao
        !* ADD           :
        !* DESCRIPTION   :
        !     ------
        !* REFERENCES    :
        !     ------
        !* LOG           :
        !     2015-05-08 :
        !* LAST MODIFIED : 2015-07-24 05:31:05 PM
        !##########################################################
        use constants
        !
        implicit none
        INTEGER(I4B), intent(in)   :: norder
        INTEGER(I4B), intent(in)   :: fnor
        real(dp), intent(out) :: coe(-norder:norder)
        INTEGER(I4B)               :: i
        !
        !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
        if (fnor .eq. 1) then
            select case (norder)
            case (1)
                coe(1) = 0.50000000000000D+00
            case (2)
                coe(1:2) = [ &
                           0.666666666666667, -8.333333333333333E-002]
            case (3)
                coe(1:3) = [ &
                           0.750000000000000, -0.150000000000000, &
                           1.666666666666667E-002]
            case (4)
                coe(1:4) = [ &
                           0.800000000000000, -0.200000000000000, &
                           3.809523809523809E-002, -3.571428571428571E-003]
            case (5)
                coe(1:5) = [ &
                           0.833333333333333, -0.238095238095238, &
                           5.952380952380952E-002, -9.920634920634920E-003, &
                           7.936507936507937E-004]
            case (6)
                coe(1:6) = [ &
                           0.857142857142857, -0.267857142857143, &
                           7.936507936507936E-002, -1.785714285714286E-002, &
                           2.597402597402598E-003, -1.803751803751804E-004]
            case (7)
                coe(1:7) = [ &
                           0.875000000000000, -0.291666666666667, &
                           9.722222222222224E-002, -2.651515151515151E-002, &
                           5.303030303030303E-003, -6.798756798756799E-004, &
                           4.162504162504163E-005]
            case (8)
                coe(1:8) = [ &
                           0.888888888888889, -0.311111111111111, &
                           0.113131313131313, -3.535353535353535E-002, &
                           8.702408702408702E-003, -1.554001554001554E-003, &
                           1.776001776001776E-004, -9.712509712509713E-006]
            case (9)
                coe(1:9) = [ &
                           0.900000000000000, -0.327272727272727, &
                           0.127272727272727, -4.405594405594405E-002, &
                           1.258741258741259E-002, -2.797202797202797E-003, &
                           4.495504495504496E-004, -4.627725215960510E-005, &
                           2.285296402943462E-006]
            case (10)
                coe(1:10) = [ &
                            0.909090909090909, -0.340909090909091, &
                            0.139860139860140, -5.244755244755244E-002, &
                            1.678321678321678E-002, -4.370629370629371E-003, &
                            8.814714697067639E-004, -1.285479226655697E-004, &
                            1.202787580496559E-005, -5.412544112234515E-007]
            end select
            coe(0) = 0.d0
            do i = 1, norder
                coe(-i) = -coe(i)
            end do
        else if (fnor .eq. 2) then

            select case (norder)
            case (1)
                coe(0:1) = [ &
                           -2.00000000000000, 1.00000000000000]
            case (2)
                coe(0:2) = [ &
                           -2.50000000000000, 1.33333333333333, &
                           -8.333333333333333E-002]
            case (3)
                coe(0:3) = [ &
                           -2.72222222222222, 1.50000000000000, &
                           -0.150000000000000, 1.111111111111111E-002]
            case (4)
                coe(0:4) = [ &
                           -2.84722222222222, 1.60000000000000, &
                           -0.200000000000000, 2.539682539682540E-002, &
                           -1.785714285714286E-003]
            case (5)
                coe(0:5) = [ &
                           -2.92722222222222, 1.66666666666667, &
                           -0.238095238095238, 3.968253968253968E-002, &
                           -4.960317460317460E-003, 3.174603174603175E-004]
            case (6)
                coe(0:6) = [ &
                           -2.98277777777778, 1.71428571428571, &
                           -0.267857142857143, 5.291005291005291E-002, &
                           -8.928571428571428E-003, 1.038961038961039E-003, &
                           -6.012506012506013E-005]
            case (7)
                coe(0:7) = [ &
                           -3.02359410430839, 1.75000000000000, &
                           -0.291666666666667, 6.481481481481481E-002, &
                           -1.325757575757576E-002, 2.121212121212121E-003, &
                           -2.266252266252267E-004, 1.189286903572618E-005]
            case (8)
                coe(0:8) = [ &
                           -3.05484410430839, 1.77777777777778, &
                           -0.311111111111111, 7.542087542087542E-002, &
                           -1.767676767676768E-002, 3.480963480963481E-003, &
                           -5.180005180005181E-004, 5.074290788576503E-005, &
                           -2.428127428127428E-006]
            case (9)
                coe(0:9) = [ &
                           -3.07953546233308, 1.80000000000000, &
                           -0.327272727272727, 8.484848484848484E-002, &
                           -2.202797202797203E-002, 5.034965034965034E-003, &
                           -9.324009324009329E-004, 1.284429855858427E-004, &
                           -1.156931303990128E-005, 5.078436450985472E-007]
            case (10)
                coe(0:10) = [ &
                            -3.09953546233308, 1.81818181818182, &
                            -0.340909090909091, 9.324009324009326E-002, &
                            -2.622377622377622E-002, 6.713286713286712E-003, &
                            -1.456876456876457E-003, 2.518489913447896E-004, &
                            -3.213698066639244E-005, 2.672861289992354E-006, &
                            -1.082508822446903E-007]
            end select
            do i = 1, norder
                coe(-i) = coe(i)
            end do
        end if
        !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
    end subroutine finite_factor
    !-----------------------PARTING-LINE--------------------------
    !-----------------------DIVIDER-LINE--------------------------
    SUBROUTINE dir2car_single(cry_coo, ort_coo, lat)
        IMPLICIT NONE
        INTEGER :: i, IDEM
        REAL(DP), DIMENSION(:) :: cry_coo, ort_coo
        REAL(DP) :: lat(3, 3)
        !>>>>>>>>>>>>>>>>>>>>>Main body>>>>>>>>>>>>>>>>>>>>>>>
        ort_coo(:) = matmul(lat, cry_coo(:))
        !<<<<<<<<<<<<<<<<<<<<<End body<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE dir2car_single
    !---------------------MODUL    E-DIVIDER    -LINE------------------------
    !-----------------------MODULE-  DIVID  ER-LINE--------------------------
    !------------------------MODULE-DI V IDER-LINE---------------------------
    !-----------------------MODULE-  DIVID  ER-LINE--------------------------
    !---------------------MODUL    E-DIVIDER    -LINE------------------------

    ! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
    ! This file is distributed under the terms of the GNU Lesser General Public
    ! License. See the file COPYING for license details.

    !BOP
    ! !ROUTINE: polynom
    ! !INTERFACE:
    function polynom(m, np, xa, ya, c, x)
        ! !INPUT/OUTPUT PARAMETERS:
        !   m  : order of derivative (in,integer)
        !   np : number of points to fit (in,integer)
        !   xa : abscissa array (in,real(np))
        !   ya : ordinate array (in,real(np))
        !   c  : work array (out,real(np))
        !   x  : evaluation abscissa (in,real)
        ! !DESCRIPTION:
        !   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
        !   function returns the $m$th derviative of the polynomial at $x$, while for
        !   $m<0$ the integral of the polynomial from the first point in the array to
        !   $x$ is returned.
        !
        ! !REVISION HISTORY:
        !   Created October 2002 (JKD)
        !EOP
        !BOC
        implicit none
        ! argmuments
        REAL(DP)                        :: polynom
        integer(i4b), intent(in) :: m
        integer(i4b), intent(in) :: np
        real(DP), intent(in) :: xa(np), ya(np)
        real(DP), intent(out) :: c(np)
        real(DP), intent(in) :: x
        ! local variables
        integer(i4b) :: i, j, k
        real(DP) :: x0, x1, x2, x3, y0, y1, y2, y3
        real(DP) :: t0, t1, t2, t3, t4, t5, t6
        real(DP) :: c1, c2, c3, sum
        ! fast evaluations for small np
        select case (np)
        case (1)
            select case (m)
            case (:-1)
                polynom = ya(1)*(x - xa(1))
            case (0)
                polynom = ya(1)
            case default
                polynom = 0.d0
            end select
            return
        case (2)
            c1 = (ya(2) - ya(1))/(xa(2) - xa(1))
            t1 = x - xa(1)
            select case (m)
            case (:-1)
                polynom = t1*(ya(1) + 0.5d0*c1*t1)
            case (0)
                polynom = c1*t1 + ya(1)
            case (1)
                polynom = c1
            case default
                polynom = 0.d0
            end select
            return
        case (3)
            x0 = xa(1)
            x1 = xa(2) - x0
            x2 = xa(3) - x0
            y0 = ya(1)
            y1 = ya(2) - y0
            y2 = ya(3) - y0
            t0 = 1.d0/(x1*x2*(x2 - x1))
            t1 = x1*y2
            t2 = x2*y1
            c1 = x2*t2 - x1*t1
            c2 = t1 - t2
            t1 = x - x0
            select case (m)
            case (:-1)
                polynom = t1*(y0 + t0*t1*(0.5d0*c1 + 0.3333333333333333333d0*c2*t1))
            case (0)
                polynom = y0 + t0*t1*(c1 + c2*t1)
            case (1)
                polynom = t0*(2.d0*c2*t1 + c1)
            case (2)
                polynom = t0*2.d0*c2
            case default
                polynom = 0.d0
            end select
            return
        case (4)
            x0 = xa(1)
            x1 = xa(2) - x0
            x2 = xa(3) - x0
            x3 = xa(4) - x0
            y0 = ya(1)
            y1 = ya(2) - y0
            y2 = ya(3) - y0
            y3 = ya(4) - y0
            t0 = 1.d0/(x1*x2*x3*(x1 - x2)*(x1 - x3)*(x2 - x3))
            t1 = x1*x2*y3
            t2 = x2*x3*y1
            t3 = x3*x1*y2
            c3 = t1*(x1 - x2) + t2*(x2 - x3) + t3*(x3 - x1)
            t6 = x3**2
            t5 = x2**2
            t4 = x1**2
            c2 = t1*(t5 - t4) + t2*(t6 - t5) + t3*(t4 - t6)
            c1 = t1*(x2*t4 - x1*t5) + t2*(x3*t5 - x2*t6) + t3*(x1*t6 - x3*t4)
            t1 = x - x0
            select case (m)
            case (:-1)
                polynom = t1*(y0 + t0*t1*(0.5d0*c1 + t1*(0.3333333333333333333d0*c2 &
                                                         + 0.25d0*c3*t1)))
            case (0)
                polynom = y0 + t0*t1*(c1 + t1*(c2 + c3*t1))
            case (1)
                polynom = t0*(c1 + t1*(2.d0*c2 + 3.d0*c3*t1))
            case (2)
                polynom = t0*(6.d0*c3*t1 + 2.d0*c2)
            case (3)
                polynom = t0*6.d0*c3
            case default
                polynom = 0.d0
            end select
            return
        end select
        if (np .le. 0) then
            WRITE (6, *)
            WRITE (6, '("Error(polynom): np <= 0 : ",I8)') np
            WRITE (6, *)
            stop
        end if
        if (m .ge. np) then
            polynom = 0.d0
            return
        end if
        ! find the polynomial coefficients in divided differences form
        c(:) = ya(:)
        do i = 2, np
            do j = np, i, -1
                c(j) = (c(j) - c(j - 1))/(xa(j) - xa(j + 1 - i))
            end do
        end do
        ! special case m=0
        if (m .eq. 0) then
            sum = c(1)
            t1 = 1.d0
            do i = 2, np
                t1 = t1*(x - xa(i - 1))
                sum = sum + c(i)*t1
            end do
            polynom = sum
            return
        end if
        x0 = xa(1)
        ! convert to standard form
        do j = 1, np - 1
            do i = 1, np - j
                k = np - i
                c(k) = c(k) + (x0 - xa(k - j + 1))*c(k + 1)
            end do
        end do
        if (m .gt. 0) then
            ! take the m th derivative
            do j = 1, m
                do i = m + 1, np
                    c(i) = c(i)*dble(i - j)
                end do
            end do
            t1 = c(np)
            t2 = x - x0
            do i = np - 1, m + 1, -1
                t1 = t1*t2 + c(i)
            end do
            polynom = t1
        else
            ! find the integral
            t1 = c(np)/dble(np)
            t2 = x - x0
            do i = np - 1, 1, -1
                t1 = t1*t2 + c(i)/dble(i)
            end do
            polynom = t1*t2
        end if
        return
    end function polynom
    !EOC
    ! SUBROUTINE IonPotentialAssignment_dg(Ity,n_Vcomp,Vcomp,zion,poscar,temp)
    SUBROUTINE IonPotentialAssignment_dg(Ity, zion, poscar, temp)
        USE constants
        USE struct_module, ONLY: struct, lat_mat
        USE MathSplines, ONLY: polynom
        USE grid_module, ONLY: rho_calc
        USE pspot_module, ONLY: max_rcut
        IMPLICIT NONE



        REAL(DP), INTENT(OUT)       :: temp(parallel%mygrid_range(3))

        INTEGER(I4B), INTENT(IN)    :: Ity !,n_Vcomp  !element number
        REAL(DP), intent(in)        :: Zion!Vcomp(n_Vcomp)
        REAL(DP)                   :: rr, poscar(3), dfr(RadiusN + 1), CarPos(3), DirPos(3)
        INTEGER(I4B)               :: i, j, k, NInOneD, fdim, m, n = 0
        REAL(DP)                   :: c(3)
        REAL(DP)                   :: rr_3d(3) !> relative position
        !========================================================================
        !==========================================================xlt test
        c = 0.d0
        ! open(998,file='Vloc_r')
        ! open(999,file='Vloc_f')
        !##INITIALIZE CALCULATE



        temp = 0.d0

        !##CALCULATE Vion GENERATE BY ATOM IN "struct%poscar(:,struct%eleid(i)+j-1)"
# 888 "IonLocalPotentialAssignment.f90"
                    !> debug
                    ! open(1180,file='Vloc')
                    ! open(1181,file='Vcomp')
                    j = 0
                    DO i = 1, parallel%mygrid_range(3), 1
                        rr = sqrt((grid%rvec(1, rho_calc%n(i)) - poscar(1))**2 + &
                             & (grid%rvec(2, rho_calc%n(i)) - poscar(2))**2 + &
                             & (grid%rvec(3, rho_calc%n(i)) - poscar(3))**2)

                        if (rr < max_rcut) then
                            rr_3d = grid%rvec(:, rho_calc%n(i)) - poscar

                            ! CALL set_Vloc_dg(Ity,n_Vcomp,Vcomp,rr_3d,temp(i))
                            CALL set_Vloc_dg(Ity, rr_3d, temp(i))
                            ! temp(i)=temp(i)-zion*erf(rr/0.5d0)/rr
                            !> debug
                            ! WRITE(1180,*)rr,temp(i)
                            ! WRITE(1181,*)rr,temp(i)+zion*erf(rr/1.0d0)/rr



                        else
                            m = 1
                            DO WHILE (psp(Ity)%r_real(m) .lt. rr)
                                !DO WHILE ( r_fv(m+1).lt.rr )
                                m = m + 1
                            END DO
                            !##CALCULATE POTENTIAL IN "RR"



                            temp(i) = polynom(0, 3, psp(Ity)%r_real(m - 1:m + 1), psp(Ity)%V_loc(m - 1:m + 1), c, rr)

                        end if
                        !temp(i,j,k)=polynom(0,3,r_fv(m-1:m+1),fv(m-1:m+1),c,rr)
                        !##PLOT
                        ! IF(rr< max_rcut)THEN
                        !  write(998,"(F20.15)")rr
                        !  write(999,"(F20.15)")temp(i)
                        ! ENDIF
                        !=================================================================





        end do
        !> debug
        ! print *,'zion',zion
        ! close(1180)
        ! close(1181)

        !=======================================================================
        ! close(998)
        ! close(999)
        !============================
        !##TEST
        ! stop "ion step"
        !============================
    END SUBROUTINE IonPotentialAssignment_Dg

    ! SUBROUTINE set_Vloc_dg(Ity,n_Vcomp,Vcomp,xyz,Vloc)
    SUBROUTINE set_Vloc_dg(Ity, xyz, Vloc)
        USE math, ONLY: Norm
        USE MathSplines, ONLY: polynom
        USE pspot_module, ONLY: max_rcut
        IMPLICIT NONE
        !>
        INTEGER(I4B), intent(IN) :: Ity !> I_type of atoms
        ! INTEGER(I4B),intent(IN) :: n_Vcomp
        ! REAL(DP),intent(IN) :: Vcomp(n_Vcomp)
        REAL(DP), intent(IN) :: xyz(3)  !> corse grid position relative atom position
        REAL(DP), intent(OUT):: Vloc    !> sum of interplote point
        !LOCAL
        REAL(DP) :: d_grid(3, numdg), d_Vloc(numdg)
        REAL(DP) :: dense_r(numdg)
        INTEGER(I4B) :: I, ix, iy, iz, m
        REAL(DP)     :: c(3)
        ! LOGICAL      :: L_isBoundary




        !>
        ! L_isBoundary=.true.
        I = 0
        c = 0.d0
        dense_r = 0.d0
        !> judge the boundary grid
        ! if(norm(xyz+gap)>max_rcut)L_isBoundary=.true.
        !> insure the dense-grid out the Boundary is not calculated
        d_Vloc = 0.d0
        do I = 1, numdg
            !do iz=Ilft,Irit
            !do iy=Ilft,Irit
            !do ix=Ilft,Irit
            !   I=I+1
            d_grid(:, I) = xyz + drijk(:, I)
            dense_r(I) = Norm(d_grid(:, I))
            ! if(L_isBoundary.and.(dense_r(I)>Norm(xyz)))cycle
            m = 1
            DO WHILE (psp(Ity)%r_real(m) .lt. dense_r(I))
                m = m + 1
            END DO
            !##CALCULATE POTENTIAL IN "dense_r"
            d_Vloc(I) = polynom(0, 3, psp(Ity)%r_real(m - 1:m + 1), psp(Ity)%V_loc(m - 1:m + 1), c, dense_r(I))
            ! d_Vloc(I)=polynom(0,3,psp(Ity)%r_real(m-1:m+1),Vcomp(m-1:m+1),c,dense_r(I))
            !enddo
            !enddo
            !enddo
        end do
        Vloc = sum(d_Vloc*wijk)
# 1019 "IonLocalPotentialAssignment.f90"
    END SUBROUTINE set_Vloc_dg
    !-----------------------PARTING-LINE--------------------------
    SUBROUTINE set_vcomp(n_inp, zion, rgauss, r_inp, V_inp, Vcomp)
        IMPLICIT NONE
        !> in/out
        INTEGER(I4B), intent(in) ::n_inp, zion
        REAL(DP), intent(in)  :: rgauss, r_inp(n_inp), V_inp(n_inp)
        REAL(DP), intent(out) :: Vcomp(n_inp)
        !> local
        INTEGER(I4B) :: i

        ! open(1180,file='V_loc')
        ! open(1181,file='V_comp')
        !> cal
        Vcomp = 0.d0
        do i = 1, n_inp
            Vcomp(i) = V_inp(i) + zion*erf(r_inp(i)/rgauss)/r_inp(i)
            ! write(1180,*)r_inp(i),V_inp(i)
            ! write(1181,*)r_inp(i),Vcomp(i)
        end do
        !> debug
        ! close(1180)
        ! close(1181)
        ! stop 'set vcomp'
    END SUBROUTINE set_vcomp
    !-----------------------PARTING-LINE--------------------------
END MODULE GetVLocalPseudoPotential
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!XXXXXXXXXXX        XXXXXX  XXXXXX  XXXXXX           XXXXXXXXXXX!
!XXXXXXXXXXX  XXXXX   XXXX  XXXXXX  XXXXXX  XXXXXXX  XXXXXXXXXXX!
!XXXXXXXXXXX        XXXXXX  XXXXXX  XXXXXX  XXXXXXXXXXXXXXXXXXXX!
!XXXXXXXXXXX  XXXX   XXXXX  XXXXXX  XXXXXX  XXXXX    XXXXXXXXXXX!
!XXXXXXXXXXX  XXXXX   XXXX  XXXXXX  XXXXXX  XXXXX    XXXXXXXXXXX!
!XXXXXXXXXXX         XXXXXX        XXXXXXX           XXXXXXXXXXX!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!BUG1: dir2car not be used in
!BUG2: dir2car used wrong
!BUG3: temp(:,:,:)or say temp write to temp(n1,n2,n3)
!BUG4: array cross the boundary and appear a non-trivial magnitude value
