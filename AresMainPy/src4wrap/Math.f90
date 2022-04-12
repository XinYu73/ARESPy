Module Math
  !##########################################################!{{{
  !* CREATED_TIME  : 2013-03-12
  !* AUTHOR        : Yanchao Wang
  !* CHANGE        : Xuecheng Shao
  !* ADD           : Qiang Xu
  !* DESCRIPTION   :
  !     ------
  !* REFERENCES    :
  !     ------
  !* LOG           :
  !     2017-07-24 :
  !* LAST MODIFIED : 2017-07-24 11:24:35 AM
  !##########################################################!}}}
  use constants
  implicit none
  INTERFACE gasdev!{{{
     MODULE PROCEDURE gasdev_s_sp, gasdev_s_dp, &
          gasdev_v_sp, gasdev_v_dp
  END INTERFACE gasdev!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  interface norm!{{{
     module procedure norm_real
     module procedure norm_complex
  end interface norm!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  interface cross!{{{
     module procedure cross_real
     module procedure cross_complex
  end interface cross!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  interface sort_id!{{{
     module procedure integer_index
     module procedure real_index
  end interface sort_id!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  !> wangsheng, add by Lantian Xue
    INTERFACE  Three2one_dim  !{{{
      MODULE PROCEDURE Three2One_d_cplx
      MODULE PROCEDURE Three2One_d_real
    END INTERFACE

    INTERFACE One2Three_dim
      MODULE PROCEDURE One2Three_d_cplx
      MODULE PROCEDURE One2Three_d_real
    END INTERFACE  !}}}
   !---------------------------- DIVIDER LINE -----------------------------
CONTAINS
  SUBROUTINE change_case(instr,str,fun)!{{{
    !##########################################################
    !* CREATED_TIME  : 2015-05-12
    !* AUTHOR        : Xuecheng Shao
    !* CHANGE        : Xuecheng Shao
    !* ADD           : Xuecheng Shao
    !* DESCRIPTION   :
    !     ------
    !* REFERENCES    :
    !     ------
    !* LOG           :
    !     2015-05-08 :
    !* LAST MODIFIED : 2015-05-12 04:50:12 PM
    !##########################################################
    implicit none
    INTEGER(I4B)              :: fun
    character(len=*)     :: instr
    character(len=*)     :: str
    !
    INTEGER(I4B)              :: ia ! 32
    INTEGER(I4B)              :: i,l_str
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    l_str=len_trim(instr)
    str = instr
    ia = ichar('a') - ichar('A') ! 32
    DO i=1,l_str
       IF (fun == 1 ) then ! lower to upper
          !if(str(i:i) >= 'a' .and. str(i:i) <= 'z') then
          if( lge(instr(i:i),'a') .and. lle(instr(i:i),'z') ) then
             str(i:i)=char( ichar(instr(i:i)) - ia )
          else
             str(i:i)=instr(i:i)
          ENDIF
       ELSE ! upper to lowe
          !if(str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
          if( lge(instr(i:i),'A') .and. lle(instr(i:i),'Z') ) then
             str(i:i)=char( ichar(instr(i:i)) + ia )
          ELSE
             str(i:i)=instr(i:i)
          ENDIF
       endif
    ENDDO
  END subroutine change_case!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine find_keywords(str,ch_mark,id_key,id_value)!{{{
    !##########################################################
    !* CREATED_TIME  : 2015-05-12
    !* AUTHOR        : Xuecheng Shao
    !* CHANGE        : Xuecheng Shao
    !* ADD           : Xuecheng Shao
    !* DESCRIPTION   :
    !     ------
    !* REFERENCES    :
    !     ------
    !* LOG           :
    !     2015-05-08 :
    !* LAST MODIFIED : 2015-05-12 04:50:12 PM
    !##########################################################
    implicit none
    character(len=*)          :: str
    CHARACTER                 :: ch_mark
    INTEGER(I4B)              :: id_key
    INTEGER(I4B)              :: id_value
    INTEGER(I4B)              :: i,l_str
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    l_str=len_trim(str)
    DO i=1,l_str
       IF ( str(i:i) == " " .or. str(i:i) == ch_mark ) exit
    ENDDO
    id_key=i
    DO i=id_key,l_str
       IF ( str(i:i) /= " " .and. str(i:i) /= ch_mark ) exit
    ENDDO
    id_value=i
    id_key=id_key-1
  END subroutine find_keywords!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine find_nword(str,ch_comma,nword)!{{{
    !##########################################################
    !* CREATED_TIME  : 2015-05-12
    !* AUTHOR        : Xuecheng Shao
    !* CHANGE        : Xuecheng Shao
    !* ADD           : Xuecheng Shao
    !* DESCRIPTION   :
    !     ------
    !* REFERENCES    :
    !     ------
    !* LOG           :
    !     2015-05-08 :
    !* LAST MODIFIED : 2015-05-12 04:50:12 PM
    !##########################################################
    implicit none
    character(len=*)          :: str
    CHARACTER                 :: ch_comma
    INTEGER(I4B)              :: nword
    INTEGER(I4B)              :: i,l_str
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    l_str=len_trim(str)
    nword = 1
    DO i = 1, l_str -1
       IF ( (str(i:i) == " " .or. str(i:i) == ch_comma) .and. &
            str(i+1:i+1) /= " " ) then
          nword = nword + 1
       ENDIF
    ENDDO
  END subroutine find_nword !}}}
  !---------------------------- DIVIDER LINE -----------------------------
  function norm_real(a)!{{{
    real(dp), intent(in) :: a(3)
    real(dp)             :: norm_real
    norm_real= sqrt(DOT_PRODUCT(a,a))
  end function norm_real!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  function norm_complex(a)!{{{
    complex(dp), intent(in) :: a(3)
    REAL(dp)             :: norm_complex
    norm_complex= sqrt(REAL(DOT_PRODUCT(a,a),DP))
  end function norm_complex!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  function cross_complex(a,b)!{{{
    complex(dp), intent(in) :: a(3),b(3)
    complex(dp)             :: cross_complex(3)

    cross_complex(1) = a(2)*b(3) - a(3)*b(2)
    cross_complex(2) = a(3)*b(1) - a(1)*b(3)
    cross_complex(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_complex!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  function cross_real(a,b)!{{{
    real(dp), intent(in) :: a(3),b(3)
    real(dp)             :: cross_real(3)

    cross_real(1) = a(2)*b(3) - a(3)*b(2)
    cross_real(2) = a(3)*b(1) - a(1)*b(3)
    cross_real(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_real!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  FUNCTION Det(matrix)!{{{
    real(DP),  intent(in) :: Matrix(3,3)
    real(DP)              :: Det


    Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
         -Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
         +Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

  END FUNCTION Det!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  function inv_33(M)!{{{
    !##########################################################!{{{
    !* CREATED_TIME  : 2015-05-19
    !* AUTHOR        : Xuecheng Shao
    !* CHANGE        : Xuecheng Shao
    !* ADD           : Xuecheng Shao
    !* DESCRIPTION   :
    !         a b c               1   ei-hf  -(bi-hc)  bf-ce
    !    mat=[d e f]  inv_33(mat)=-- [fg-id  -(cg-ia)  cd-af]
    !         g h i          det(mat) dh-ge  -(ah-gb)  ae-bd
    !* REFERENCES    :
    !     ------
    !* LOG           :
    !     2015-05-08 :
    !* LAST MODIFIED : 2015-05-19 09:15:08 AM
    !##########################################################!}}}
    real(dp),intent(in) :: M(3,3)
    real(dp)            :: Inv(3,3),inv_33(3,3)
    real(dp)            ::  d
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    d = Det(M)
    inv(1,1)=M(2,2)*M(3,3) - M(2,3)*M(3,2)
    inv(2,1)=M(2,3)*M(3,1) - M(2,1)*M(3,3)
    inv(3,1)=M(2,1)*M(3,2) - M(2,2)*M(3,1)
    inv(1,2)=M(1,3)*M(3,2) - M(1,2)*M(3,3)
    inv(2,2)=M(1,1)*M(3,3) - M(1,3)*M(3,1)
    inv(3,2)=M(1,2)*M(3,1) - M(1,1)*M(3,2)
    inv(1,3)=M(1,2)*M(2,3) - M(1,3)*M(2,2)
    inv(2,3)=M(1,3)*M(2,1) - M(1,1)*M(2,3)
    inv(3,3)=M(1,1)*M(2,2) - M(1,2)*M(2,1)
    inv_33=inv/d
    !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
  end function inv_33!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  FUNCTION LindG(eta,lambda,mu)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!This function from PROFESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(kind=DP), INTENT(IN)      :: &
         eta, &                 ! the point at which the function gets evaluated.
         lambda, &              ! the TF multiplier for compensating.
         mu                     ! the vW multiplier

    REAL(kind=DP)                  :: &
         LindG                  ! The Lindhard G function as described in [1]

    !>> INTERNAL VARIABLES <<!
    REAL(kind=DP)                  :: &
         eta2, &                ! eta ** 2
         invEta2                ! 1/eta**2
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>

    IF (eta<0._DP) THEN
       LindG = 0._DP

       ! Limit for small eta
    ELSE IF (eta < 1E-10_DP) THEN
       LindG = 1._DP - lambda + eta**2 * (1._DP / 3._DP - 3._DP * mu)

       ! Around the singularity
    ELSE IF (ABS(eta-1._DP) < 1E-10_DP) THEN
       LindG = 2._DP - lambda - 3._DP * mu + 20._DP * (eta-1._DP)

       ! Taylor expansion for high eta
    ELSE IF (eta > 3.65_DP) THEN ! we determined empircally that 3.65 was a
       ! good crossover point to the taylor expansion
       eta2 = eta**2
       invEta2 = 1._DP/eta2
       LindG = 3._DP*(1._DP-mu)*eta2 &
            -lambda-0.6_DP &
            + invEta2 * (-0.13714285714285712_DP &
            + invEta2 * (-6.39999999999999875E-2_DP &
            + invEta2 * (-3.77825602968460128E-2_DP &
            + invEta2 * (-2.51824061652633074E-2_DP &
            + invEta2 * (-1.80879839616166146E-2_DP &
            + invEta2 * (-1.36715733124818332E-2_DP &
            + invEta2 * (-1.07236045520990083E-2_DP &
            + invEta2 * (-8.65192783339199453E-3_DP &
            + invEta2 * (-7.1372762502456763E-3_DP &
            + invEta2 * (-5.9945117538835746E-3_DP &
            + invEta2 * (-5.10997527675418131E-3_DP &
            + invEta2 * (-4.41060829979912465E-3_DP &
            + invEta2 * (-3.84763737842981233E-3_DP &
            + invEta2 * (-3.38745061493813488E-3_DP &
            + invEta2 * (-3.00624946457977689E-3_DP)))))))))))))))

    ELSE
       LindG = 1._DP / (0.5_DP + 0.25_DP * (1._DP-eta**2) * LOG((1._DP + eta)/&
            abs(1._DP-eta))/eta) - 3._DP * mu * eta**2 - lambda
    END IF
    !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
  END FUNCTION LindG!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  Function int_to_char( int )!{{{
    !-----------------------------------------------------------------------
    !
    INTEGER(I4B), intent(in) :: int
    character (len=6)   :: int_to_char
    if ( int <= -10 ) then
       write( unit = int_to_char , fmt = "(i3)" ) int
    else if(int == 0 ) then
       write(unit = int_to_char,fmt= "(i1)") int
    else if(int< 0 ) then
       write( unit = int_to_char , fmt = "(i2)" ) int
    else if ( int < 10 ) then
       write( unit = int_to_char , fmt = "(i1)" ) int
    else if ( int < 100 ) then
       write( unit = int_to_char , fmt = "(i2)" ) int
    else if ( int < 1000 ) then
       write( unit = int_to_char , fmt = "(i3)" ) int
    else if ( int < 10000 ) then
       write( unit = int_to_char , fmt = "(i4)" ) int
    else
       write( unit = int_to_char , fmt = "(i5)" ) int
    end if
    return
  end function int_to_char!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine lat2matrix(lat_para,lat_mat,flag)!{{{
    !##########################################################!{{{
    !* CREATED_TIME  : 2015-05-11 10:58:47
    !* AUTHOR        : Xuecheng Shao
    !* CHANGE        : Xuecheng Shao
    !* ADD           : Xuecheng Shao
    !* DESCRIPTION   :
    !     ------
    !    [ A   B   C  ]
    !      M11 M12 M13  M13 = lA * cos(b)
    !      0   M22 M23  M23 = ( lB*lC*cos(x) - M13*M12 ) / M22
    !      0   0   M33
    !* REFERENCES    :
    !     ------
    !* LOG           :
    !     2015-05-11
    !* LAST MODIFIED : 2015-05-14 10:27:21 AM
    !##########################################################!}}}
    !
    USE constants,     only : DP
    !
    implicit none
    INTEGER(I4B)              :: flag
    real(DP)             :: lat_para(6)
    real(DP)             :: lat_mat(3,3)
    real(DP)             :: angle(3)
    if (flag==1) then
       lat_mat=0.0
       lat_mat(1,1) = lat_para(1)
       lat_mat(1,2) = lat_para(2)*cos(lat_para(6))
       lat_mat(2,2) = lat_para(2)*sin(lat_para(6))
       lat_mat(1,3) = lat_para(3)*cos(lat_para(5))
       lat_mat(2,3) = ( lat_para(2)*lat_para(3)*cos(lat_para(4)) - lat_mat(1,3)*lat_mat(1,2) )/ lat_mat(2,2)
       lat_mat(3,3) = sqrt(lat_para(3)**2 - lat_mat(1,3)**2 - lat_mat(2,3)**2)
    else
       lat_para(1)=dsqrt(sum(lat_mat(:,1)**2))
       lat_para(2)=dsqrt(sum(lat_mat(:,2)**2))
       lat_para(3)=dsqrt(sum(lat_mat(:,3)**2))
       angle(1)=(dot_product(lat_mat(:,2),lat_mat(:,3))) / (lat_para(2)*lat_para(3))
       angle(2)=(dot_product(lat_mat(:,1),lat_mat(:,3))) / (lat_para(1)*lat_para(3))
       angle(3)=(dot_product(lat_mat(:,1),lat_mat(:,2))) / (lat_para(1)*lat_para(2))
       lat_para(4)=acos(angle(1))
       lat_para(5)=acos(angle(2))
       lat_para(6)=acos(angle(3))
    endif
  end subroutine lat2matrix!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine one2three(id,n_dens,pos)!{{{
    implicit none
    INTEGER(I4B)                      :: id,n_dens(3),pos(3)
    INTEGER(I4B)                      :: ixx,iyz,iyy,izz
    ixx=mod(id,n_dens(1))
    iyz=id/n_dens(1)
    iyy=mod(iyz,n_dens(2))
    izz=iyz/n_dens(2)
    if (ixx==0) then
       pos(1)=n_dens(1)
    else
       pos(1)=ixx
    endif
    if (ixx==0 .and. iyy==0) then
       pos(2)=n_dens(2)
       pos(3)=izz
    else
       pos(3)=izz+1
       if (ixx==0) then
          pos(2)=iyy
       else
          pos(2)=iyy+1
       endif
    endif
  END subroutine one2three!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine real_index(array,n,id) !{{{
    USE constants,         only : DP
    implicit none
    INTEGER(I4B)                     :: n
    real(DP)                    :: array(n)
    INTEGER(I4B)                     :: id(n)
    INTEGER(I4B)                     :: i,j,k,l,m
    real(DP),parameter          :: eps=1.d-14
    do i=1,n
       id(i)=i
    enddo
    if (n==1) return
    l=n/2+1
    k=n
    do while(.true.)
       if (l>1) then
          l=l-1
          m=id(l)
       else
          m=id(k)
          id(k)=id(1)
          k=k-1
          if (k==1) then
             id(1)=m
             return
          endif
       endif
       i=l
       j=l+1
       do while(j<=k)
          if (j<k) then
             if (array(id(j))<array(id(j+1))+eps) j=j+1
          endif
          if (array(m)<array(id(j))+eps) then
             id(i)=id(j)
             i=j
             j=j+j
          else
             j=k+1
          endif
       enddo
       id(i)=m
    enddo
  end subroutine real_index!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine integer_index(array,n,id) !{{{
    USE constants,    only : DP
    implicit none
    INTEGER(I4B)                     :: n
    INTEGER(I4B)                     :: array(n)
    INTEGER(I4B)                     :: id(n)
    INTEGER(I4B)                     :: i,j,k,l,m
    real(DP),parameter          :: eps=1.d-14
    do i=1,n
       id(i)=i
    enddo
    if (n==1) return
    l=n/2+1
    k=n
    do while(.true.)
       if (l>1) then
          l=l-1
          m=id(l)
       else
          m=id(k)
          id(k)=id(1)
          k=k-1
          if (k==1) then
             id(1)=m
             return
          endif
       endif
       i=l
       j=l+1
       do while(j<=k)
          if (j<k) then
             if (array(id(j))<array(id(j+1))+eps) j=j+1
          endif
          if (array(m)<array(id(j))+eps) then
             id(i)=id(j)
             i=j
             j=j+j
          else
             j=k+1
          endif
       enddo
       id(i)=m
    enddo
  end subroutine integer_index!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine init_random_seed()!{{{
    !##########################################################
    !* CREATED_TIME  : 2015-07-20 19:38:02
    !* AUTHOR        : Xuecheng Shao
    !* CHANGE        : Xuecheng Shao
    !* Mail          : sxc@calypso.cn
    !* DESCRIPTION   :
    !     ------
    !* REFERENCES    :
    !   https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
    !* LOG           :
    !     2015-07-20
    !* LAST MODIFIED : 2015-07-20 08:16:23 PM
    !##########################################################
    use iso_fortran_env, only: int64
    ! In order to use the function GETPID
    ! ifort need add ifport module, gfortran no need
    ! use ifport
    implicit none
    INTEGER(I4B), allocatable :: seed(:)
    INTEGER(I4B) :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       ! pid = getpid()
       pid = 1024
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer(i4b) :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  SUBROUTINE gasdev_s_sp(harvest, inmu, insigma)!{{{
    IMPLICIT NONE
    REAL(SP), INTENT(OUT) :: harvest
    REAL(SP), OPTIONAL, INTENT(IN) :: inmu
    REAL(SP), OPTIONAL, INTENT(IN) :: insigma
    REAL(SP)                       :: mu
    REAL(SP)                       :: sigma
    REAL(SP) :: rsq,v1,v2
    REAL(SP), SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    IF ( .not. PRESENT(inmu) ) mu = 0_sp
    IF ( .not. PRESENT(insigma) ) sigma = 1_sp

    if (gaus_stored) then
       harvest = g * sigma + mu
       gaus_stored=.false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1=2.0_sp*v1-1.0_sp
          v2=2.0_sp*v2-1.0_sp
          rsq=v1**2+v2**2
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       rsq=sqrt(-2.0_sp*log(rsq)/rsq)
       harvest=v1*rsq
       g=v2*rsq
       harvest = harvest * sigma + mu
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev_s_sp !}}}
  !---------------------------- DIVIDER LINE -----------------------------
  SUBROUTINE gasdev_s_dp(harvest, inmu, insigma)!{{{
    IMPLICIT NONE
    REAL(DP), OPTIONAL, INTENT(IN) :: inmu
    REAL(DP), OPTIONAL, INTENT(IN) :: insigma
    REAL(DP)                       :: mu
    REAL(DP)                       :: sigma
    REAL(DP), INTENT(OUT) :: harvest
    REAL(DP) :: rsq,v1,v2
    REAL(DP), SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    IF ( .not. PRESENT(inmu) ) mu = 0_Dp
    IF ( .not. PRESENT(insigma) ) sigma = 1_Dp
    if (gaus_stored) then
       harvest = g * sigma + mu
       gaus_stored=.false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1=2.0_sp*v1-1.0_sp
          v2=2.0_sp*v2-1.0_sp
          rsq=v1**2+v2**2
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       rsq=sqrt(-2.0_sp*log(rsq)/rsq)
       harvest=v1*rsq
       g=v2*rsq
       harvest = harvest * sigma + mu
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev_s_dp !}}}
  !---------------------------- DIVIDER LINE -----------------------------
  SUBROUTINE gasdev_v_sp(harvest, inmu, insigma)!{{{
    IMPLICIT NONE
    REAL(SP), OPTIONAL, INTENT(IN) :: inmu
    REAL(SP), OPTIONAL, INTENT(IN) :: insigma
    REAL(SP)                       :: mu
    REAL(SP)                       :: sigma
    REAL(SP), INTENT(OUT)                   :: harvest(:)
    REAL(SP)                                :: rsq
    REAL(SP), dimension(size(harvest))      :: v1,v2
    REAL(SP), allocatable,dimension(:),SAVE :: g
    INTEGER(I4B), SAVE                      :: last_allocated=0
    LOGICAL, SAVE                           :: gaus_stored=.false.
    INTEGER(I4B)                            :: i,n
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    IF ( .not. PRESENT(inmu) ) mu = 0_sp
    IF ( .not. PRESENT(insigma) ) sigma = 1_sp
    n=size(harvest)
    if (n /= last_allocated) then
       if (last_allocated /= 0) deallocate(g)
       allocate(g(n))
       last_allocated=n
       gaus_stored=.false.
    end if
    if (gaus_stored) then
       harvest(:) = g(:) * sigma + mu
       gaus_stored=.false.
    else
       do i = 1, size(harvest)
          do
             call random_number(v1(i))
             call random_number(v2(i))
             v1(i) = 2.0_sp*v1(i) -1.0_sp
             v2(i) = 2.0_sp*v2(i) -1.0_sp
             rsq = v1(i)**2 + v2(i)**2
             if (rsq > 0.0 .and. rsq < 1.0) exit
          end do
          rsq=sqrt(-2.0_sp*log(rsq)/rsq)
          harvest(i) = v1(i) * rsq
          g(i) = v2(i) * rsq
          harvest(i) = harvest(i) * sigma + mu
       enddo
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev_v_sp!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  SUBROUTINE gasdev_v_dp(harvest, inmu, insigma)!{{{
    IMPLICIT NONE
    REAL(DP), OPTIONAL, INTENT(IN) :: inmu
    REAL(DP), OPTIONAL, INTENT(IN) :: insigma
    REAL(DP)                       :: mu
    REAL(DP)                       :: sigma
    REAL(DP), INTENT(OUT)                   :: harvest(:)
    REAL(DP)                                :: rsq
    REAL(DP), dimension(size(harvest))      :: v1,v2
    REAL(DP), allocatable,dimension(:),SAVE :: g
    INTEGER(I4B), SAVE                      :: last_allocated=0
    LOGICAL, SAVE                           :: gaus_stored=.false.
    INTEGER(I4B)                            :: i,n
    !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
    IF ( .not. PRESENT(inmu) ) mu = 0_dp
    IF ( .not. PRESENT(insigma) ) sigma = 1_dp
    n=size(harvest)
    if (n /= last_allocated) then
       if (last_allocated /= 0) deallocate(g)
       allocate(g(n))
       last_allocated=n
       gaus_stored=.false.
    end if
    if (gaus_stored) then
       harvest(:) = g(:) * sigma + mu
       gaus_stored=.false.
    else
       do i = 1, size(harvest)
          do
             call random_number(v1(i))
             call random_number(v2(i))
             v1(i) = 2.0_sp*v1(i) -1.0_sp
             v2(i) = 2.0_sp*v2(i) -1.0_sp
             rsq = v1(i)**2 + v2(i)**2
             if (rsq > 0.0 .and. rsq < 1.0) exit
          end do
          rsq=sqrt(-2.0_sp*log(rsq)/rsq)
          harvest(i) = v1(i) * rsq
          g(i) = v2(i) * rsq
          harvest(i) = harvest(i) * sigma + mu
       enddo
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev_v_dp!}}}
  !---------------------------- DIVIDER LINE -----------------------------
  SUBROUTINE atom_mass(atom_name,mass)!{{{
    use constants,  only : DP, CONST_MA_AU
    implicit none
    !--------------------------------------------
    character(len=*)    :: atom_name
    integer             :: i
    real(DP)            :: mass
    !------------------------------------------------
    type :: elements
       integer          ::atom_number
       character(len=3) :: symbols
       real(DP)         :: mass
    end type elements
    type(elements) ::atom(104)
    !-------------------------------------------
    !print *,"name",atom_name
    do i=1,104
       if(adjustl(atom_name)==atom(i)%symbols) then
          mass=atom(i)%mass*CONST_MA_AU
          !print *,"masssss",mass
          !print *,"name",atom_name
       endif
    enddo
    data atom(1:104) /&
         & elements(  1, "H",  1.0079470000d0),&
         & elements(  2, "He", 4.0026022000d0),&
         & elements(  3, "Li", 6.9412000000d0),&
         & elements(  4, "Be", 9.0121823000d0),&
         & elements(  5, "B",  10.811700000d0),&
         & elements(  6, "C",  12.010780000d0),&
         & elements(  7, "N",  14.006720000d0),&
         & elements(  8, "O",  15.999400000d0),&
         & elements(  9, "F",  18.998403250d0),&
         & elements( 10, "Ne", 20.179760000d0),&
         & elements( 11, "Na", 22.989769282d0),&
         & elements( 12, "Mg", 24.305060000d0),&
         & elements( 13, "Al", 26.981538680d0),&
         & elements( 14, "Si", 28.085530000d0),&
         & elements( 15, "P",  30.973762200d0),&
         & elements( 16, "S",  32.065500000d0),&
         & elements( 17, "Cl", 35.453200000d0),&
         & elements( 18, "Ar", 39.948100000d0),&
         & elements( 19, "K",  39.098310000d0),&
         & elements( 20, "Ca", 40.078400000d0),&
         & elements( 21, "Sc", 44.955912600d0),&
         & elements( 22, "Ti", 47.867100000d0),&
         & elements( 23, "V",  50.941510000d0),&
         & elements( 24, "Cr", 51.996160000d0),&
         & elements( 25, "Mn", 54.938045500d0),&
         & elements( 26, "Fe", 55.845200000d0),&
         & elements( 27, "Co", 58.933195220d0),&
         & elements( 28, "Ni", 58.693442000d0),&
         & elements( 29, "Cu", 63.546300000d0),&
         & elements( 30, "Zn", 65.382000000d0),&
         & elements( 31, "Ga", 69.723500000d0),&
         & elements( 32, "Ge", 72.641000000d0),&
         & elements( 33, "As", 74.921595600d0),&
         & elements( 34, "Se", 78.963000000d0),&
         & elements( 35, "Br", 79.904100000d0),&
         & elements( 36, "Kr", 83.798200000d0),&
         & elements( 37, "Rb", 85.467830000d0),&
         & elements( 38, "Sr", 87.621000000d0),&
         & elements( 39, "Y",  88.905852000d0),&
         & elements( 40, "Zr", 91.224000000d0),&
         & elements( 41, "Nb", 92.906382000d0),&
         & elements( 42, "Mo", 95.962000000d0),&
         & elements( 43, "Tc", 98.000000000d0),&
         & elements( 44, "Ru", 101.07200000d0),&
         & elements( 45, "Rh", 102.90550300d0),&
         & elements( 46, "Pd", 106.42100000d0),&
         & elements( 47, "Ag", 107.86822000d0),&
         & elements( 48, "Cd", 112.41180000d0),&
         & elements( 49, "In", 114.81830000d0),&
         & elements( 50, "Sn", 118.71070000d0),&
         & elements( 51, "Sb", 121.76010000d0),&
         & elements( 52, "Te", 127.60300000d0),&
         & elements( 53, "I",  126.90447300d0),&
         & elements( 54, "Xe", 131.29360000d0),&
         & elements( 55, "Cs", 132.90545192d0),&
         & elements( 56, "Ba", 137.32770000d0),&
         & elements( 57, "La", 138.90547700d0),&
         & elements( 58, "Ce", 140.11610000d0),&
         & elements( 59, "Pr", 140.90765200d0),&
         & elements( 60, "Nd", 144.24230000d0),&
         & elements( 61, "Pm", 145.00000000d0),&
         & elements( 62, "Sm", 150.36200000d0),&
         & elements( 63, "Eu", 151.96410000d0),&
         & elements( 64, "Gd", 157.25300000d0),&
         & elements( 65, "Tb", 158.92535200d0),&
         & elements( 66, "Dy", 162.50010000d0),&
         & elements( 67, "Ho", 164.93032300d0),&
         & elements( 68, "Er", 167.25930000d0),&
         & elements( 69, "Tm", 168.93421200d0),&
         & elements( 70, "Yb", 173.05450000d0),&
         & elements( 71, "Lu", 174.96681000d0),&
         & elements( 72, "Hf", 178.49200000d0),&
         & elements( 73, "Ta", 180.94788200d0),&
         & elements( 74, "W",  183.84100000d0),&
         & elements( 75, "Re", 186.20710000d0),&
         & elements( 76, "Os", 190.23000000d0),&
         & elements( 77, "Ir", 192.21730000d0),&
         & elements( 78, "Pt", 195.08490000d0),&
         & elements( 79, "Au", 196.96656940d0),&
         & elements( 80, "Hg", 200.59200000d0),&
         & elements( 81, "Ti", 204.38332000d0),&
         & elements( 82, "Pb", 207.21000000d0),&
         & elements( 83, "Bi", 208.98040100d0),&
         & elements( 84, "Po", 209.00000000d0),&
         & elements( 85, "At", 210.00000000d0),&
         & elements( 86, "Rn", 222.00000000d0),&
         & elements( 87, "Fr", 223.00000000d0),&
         & elements( 88, "Ra", 226.00000000d0),&
         & elements( 89, "Ac", 227.00000000d0),&
         & elements( 90, "Th", 232.03806200d0),&
         & elements( 91, "Pa", 231.03588200d0),&
         & elements( 92, "U",  238.02891300d0),&
         & elements( 93, "Np", 237.00000000d0),&
         & elements( 94, "Pu", 244.00000000d0),&
         & elements( 95, "Am", 243.00000000d0),&
         & elements( 96, "Cm", 247.00000000d0),&
         & elements( 97, "Bk", 247.00000000d0),&
         & elements( 98, "Cf", 251.00000000d0),&
         & elements( 99, "Es", 252.00000000d0),&
         & elements(100, "Fm", 257.00000000d0),&
         & elements(101, "Md", 258.00000000d0),&
         & elements(102, "No", 259.00000000d0),&
         & elements(103, "Lr", 262.00000000d0),&
         & elements(104, "Rf", 265.00000000d0)/
    ! & elements(105, "Db",268.000000)
    ! & elements(106, "Sg",271.000000)
    ! & elements(107, "Bh",272.000000)
    ! & elements(108, "Hs",277.000000)
    ! & elements(109, "Mt",276.000000)
    ! & elements(110, "Ds",281.000000)
    ! & elements(111, "Rg",280.000000)
    ! & elements(112, "Cn",285.000000)
    ! & elements(113, "Uut",284.000000)
    ! & elements(114, "Fl",289.000000)
    ! & elements(115, "Uup",288.000000)
    ! & elements(116, "Lv",293.000000)
    ! & elements(117, "Uus",294.000000)
    ! & elements(118, "Uuo",294.000000) /
    !--------------------------------------------------
  END SUBROUTINE atom_mass !}}}
  !---------------------------- DIVIDER LINE -----------------------------
  subroutine newton_inter(n,x,y,m,tx,ty)!{{{
    !       《数值计算方法》
    !       林成森  科学出版社
    !----------------------------------------------------
    INTEGER(i4b) :: n
    INTEGER :: m
    REAL(DP)  :: tx(m),ty(m)
    REAL(DP)  :: x(0:n-1),y(0:n-1)
    !
    INTEGER(i4b) :: i,j,k
    REAL(DP)  :: qx(n)
    REAL(DP)  :: Q(0:n-1,0:n-1)
    q(:,0)=y(:)
    do i=1,n-1
       do j=1,i
          Q(i,j)=(Q(i,j-1)-Q(i-1,j-1))/(x(i)-x(i-j))
       end do
    end do
    qx(n)=Q(n-1,n-1)
    DO i = 1, m
       do k=n-1,1,-1
          qx(k)=Q(k-1,k-1)+qx(k+1)*(tx(i)-x(k-1))
       end do
       ty(i)=qx(1)
    ENDDO
  end subroutine newton_inter !}}}
  !---------------------------- DIVIDER LINE -----------------------------
  SUBROUTINE diag(N, inA, W, Q)!{{{
    ! ----------------------------------------------------------------------------
    ! Numerical diagonalization of nxn matrcies
    ! Copyright (C) 2006  Joachim Kopp
    ! ----------------------------------------------------------------------------
    ! This library is free software; you can redistribute it and/or
    ! modify it under the terms of the GNU Lesser General Public
    ! License as published by the Free Software Foundation; either
    ! version 2.1 of the License, or (at your option) any later version.

    ! This library is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    ! Lesser General Public License for more details.

    ! You should have received a copy of the GNU Lesser General Public
    ! License along with this library; if not, write to the Free Software
    ! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues and normalized eigenvectors of a symmetric nxn
    ! matrix A using the Jacobi algorithm.
    ! The upper triangular part of A is destroyed during the calculation,
    ! the diagonal elements are read but not destroyed, and the lower
    ! triangular elements are not referenced at all.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   Q: Storage buffer for eigenvectors
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    !     .. Arguments ..
    implicit none
    !     .. Parameters ..
    INTEGER(I4B)                     :: N
    real(DP) :: inA(N,N)
    real(DP) :: A(N,N)
    real(DP) :: Q(N,N)
    real(DP) :: W(N)
    !     .. Local Variables ..
    real(DP) :: SD, SO
    real(DP) :: S, C, T
    real(DP) :: G, H, Z, THETA
    real(DP) :: THRESH
    INTEGER(I4B)                     :: I, X, Y, R

    !     Initialize Q to the identitity matrix
    !     --- This loop can be omitted if only the eigenvalues are desired ---
    A = inA
    !WRITE(6,*)"AA"
    !WRITE(6,*)A
    DO X = 1, N
       Q(X,X) = 1.0D0
       DO Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
       END DO
    END DO

    !     Initialize W to diag(A)
    DO X = 1, N
       W(X) = A(X, X)
    END DO

    !     Calculate SQR(tr(A))
    SD = 0.0D0
    DO X = 1, N
       SD = SD + ABS(W(X))
    END DO
    SD = SD**2

    !     Main iteration loop
    DO I = 1, 50
       !       Test for convergence
       SO = 0.0D0
       DO X = 1, N
          DO Y = X+1, N
             SO = SO + ABS(A(X, Y))
          END DO
       END DO
       IF (SO == 0.0D0) THEN
          RETURN
       END IF

       IF (I < 4) THEN
          THRESH = 0.2D0 * SO / N**2
       ELSE
          THRESH = 0.0D0
       END IF

       !       Do sweep
       DO X = 1, N
          DO Y = X+1, N
             G = 100.0D0 * ( ABS(A(X, Y)) )
             IF ( I > 4 .AND. ABS(W(X)) + G == ABS(W(X)) &
                  .AND. ABS(W(Y)) + G == ABS(W(Y)) ) THEN
                A(X, Y) = 0.0D0
             ELSEIF (ABS(A(X, Y)) > THRESH) THEN
                !             Calculate Jacobi transformation
                H = W(Y) - W(X)
                IF ( ABS(H) + G == ABS(H) ) THEN
                   T = A(X, Y) / H
                ELSE
                   THETA = 0.5D0 * H / A(X, Y)
                   IF (THETA < 0.0D0) THEN
                      T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                   ELSE
                      T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                   END IF
                END IF

                C = 1.0D0 / SQRT( 1.0D0 + T**2 )
                S = T * C
                Z = T * A(X, Y)

                !             Apply Jacobi transformation
                A(X, Y) = 0.0D0
                W(X)    = W(X) - Z
                W(Y)    = W(Y) + Z
                DO R = 1, X-1
                   T       = A(R, X)
                   A(R, X) = C * T - S * A(R, Y)
                   A(R, Y) = S * T + C * A(R, Y)
                END DO
                DO R = X+1, Y-1
                   T       = A(X, R)
                   A(X, R) = C * T - S * A(R, Y)
                   A(R, Y) = S * T + C * A(R, Y)
                END DO
                DO R = Y+1, N
                   T       = A(X, R)
                   A(X, R) = C * T - S * A(Y, R)
                   A(Y, R) = S * T + C * A(Y, R)
                END DO
                !Update eigenvectors
                !--- This loop can be omitted if only the eigenvalues are desired ---
                DO R = 1, N
                   T       = Q(R, X)
                   Q(R, X) = C * T - S * Q(R, Y)
                   Q(R, Y) = S * T + C * Q(R, Y)
                END DO
             ENDIF
          END DO
       END DO
    END DO

    PRINT *, "No convergence."

  END SUBROUTINE diag!}}}
  !-----------------------------------------------------------------------
  function boltzmann_distribution(rnull,width)!{{{
    real(DP) :: num1,num2
    real(DP) :: boltzmann_distribution,width,rnull
    REAL(DP),parameter       :: tpi = 8.d0*atan(1.d0)
    !real(q),parameter :: twopi = 6.283185307179586_q

    call random_number(num1)
    num2=0.d0
    do
       call random_number(num2)
       num2=abs(num2)
       if (num2 .gt. 1e-08) exit
    enddo
    !WRITE(6,*) 'num1,num2',num1,num2
    boltzmann_distribution= cos( tpi*num1 ) * sqrt( 2.d0 *abs(log(num2)) )
    !WRITE(6,*) 'boltzmann_distribution_',boltzmann_distribution,width,rnull
    boltzmann_distribution = width * boltzmann_distribution  +  rnull
    !WRITE(6,*) 'boltzmann_distribution',boltzmann_distribution
  END FUNCTION boltzmann_distribution !}}}
  !>>>>ADD by Qiang Xu
  !#####################################################!
  !cry_cyyoordinates to ort_coordinates
  !#####################################################!
  SUBROUTINE dir2car(cry_coo,ort_coo,lat)
    IMPLICIT NONE
    INTEGER :: i,IDEM
    REAL(DP),DIMENSION(:,:) :: cry_coo,ort_coo
    REAL(DP) :: lat(3,3)
    !>>>>>>>>>>>>>>>>>>>>>Main body>>>>>>>>>>>>>>>>>>>>>>>
    IDEM = size(cry_coo,2)
    DO i = 1,IDEM
       ort_coo(:,i) =  matmul(lat,cry_coo(:,i))
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<End body<<<<<<<<<<<<<<<<<<<<<<<<
  END SUBROUTINE dir2car
  !#####################################################!
  !ort_coordinates to cry_coordinates
  !#####################################################!
  !----------------------ort to cry------------------------
  SUBROUTINE car2dir(ort_coo,cry_coo,lat)
    IMPLICIT NONE
    INTEGER :: i,IDEM
    REAL(DP),DIMENSION(:,:) :: cry_coo,ort_coo
    REAL(DP):: inv_lattice(3,3),lat(3,3)
    !>>>>>>>>>>>>>>>>>>>>>Main body>>>>>>>>>>>>>>>>>>>>>>>>
    IDEM = size(ort_coo,2)
    inv_lattice =inv_33(lat)
    DO i=1,IDEM
       cry_coo(:,i)= matmul(inv_lattice,ort_coo(:,i))
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<End body<<<<<<<<<<<<<<<<<<<<<<<<<
  END SUBROUTINE car2dir
  !----------------------3D-to-matrix--------------------------
  SUBROUTINE thr2mat(n1,n2,n3,I,J,K,dimnu)
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: I,J,K,n1,n2,n3
    INTEGER(I4B),INTENT(OUT) :: dimnu
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    dimnu=(K-1)*n2*n1+(J-1)*n1+I
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE thr2mat
  !---------------------matrix-to-3D---------------------------
  SUBROUTINE mat2thr(n1,n2,n3,I,Ix,Iy,Iz)
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: n1,n2,n3,I
    INTEGER(I4B),INTENT(OUT) :: Ix,Iy,Iz
    !
    INTEGER(I4B) :: timz,timy,delta,deltad
    INTEGER(I4B) :: n21
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    n21=n2*n1
    !Iz
    timz=I/n21
    !delta=I-timz*n21
    delta=MOD(I,n21)
    IF(delta==0)THEN
       Iz=timz
       Iy=n2
       Ix=n1
    ELSE
       Iz=timz+1
       !------------------------
       timy=delta/n1
       !deltad=delta-timy*n1
       deltad=MOD(delta,n1)
       IF(deltad==0)THEN
          Iy=timy
          Ix=n1
       ELSE
          Iy=timy+1
          Ix=deltad
       ENDIF
    ENDIF
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE mat2thr
  !------------------------------------------------------------
  SUBROUTINE SOPO(A,LDA,N,B,LDB,M,W)
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: &
         N, M, LDA, LDB
    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: &
         A
    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: &
         B
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: &
         W
    !C---------------------------------------------------------------
    !C      SOLUTION OF A SET OF LINEAR EQUATIONS WITH A POSITIVE
    !C      DEFINITE REAL SYMMETRIC MATRIX AND WITH MULTIPLE R.H.S.
    !C---------------------------------------------------------------

    INTEGER :: &
         I, J, L, MQ, Ntemp

    REAL(KIND=8) :: &
         DUM
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Ntemp=LDA+LDB
    !C                    ROZKLAD
    !> make the positive define real symmetric matrix A to be a upper triangular matrix
    DO L=1,N-1
       A(L,L)=SQRT(A(L,L))
       DUM=1.d0/A(L,L)
       DO I=L+1,N
          W(I)=A(I,L)*DUM
          A(I,L)=W(I)
       ENDDO
       DO J=L+1,N
          DO I=J,N
             A(I,J)=A(I,J)-W(I)*W(J)
          ENDDO
       ENDDO
    ENDDO
    A(N,N)=SQRT(A(N,N))

    !C                   CYKLUS PRES PRAVE STRANY
    DO MQ=1,M
       !C                   INVERZE DOLNI TROJUH. MATICE
       DO L=1,N-1
          W(L)=B(L,MQ)/A(L,L)
          DO I=L+1,N
             B(I,MQ)=B(I,MQ)-A(I,L)*W(L)
          ENDDO
       ENDDO
       W(N)=B(N,MQ)/A(N,N)
       !C                   INVERZE HORNI TROJUH. MATICE
       B(N,MQ)=W(N)/A(N,N)
       DO L=N-1,1,-1
          DO I=N,L+1,-1
             W(L)=W(L)-B(I,MQ)*A(I,L)
          ENDDO
          B(L,MQ)=W(L)/A(L,L)
       ENDDO

    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE SOPO
  !------------------------------------------------------------
  SUBROUTINE csort_eigen(nev,arr,brr)
    !#################################################################!
    ! For a given set of eigenvalues (for all representations), sort  !
    ! them in ascending order. Use straight insertion, no need for    !
    ! fancy algorithms since the arrays are always short.             !
    ! ISBN 7-03-010217-7                                              !
    !#################################################################!
    !
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: nev
    REAL(DP),INTENT(INOUT) :: arr(:)
    COMPLEX(DP),INTENT(INOUT) :: brr(:,:)
    !
    REAL(DP) :: tmpa
    COMPLEX(DP) :: tmpb(SIZE(brr,1))
    INTEGER(I4B) :: I,J
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    DO J=2,nev
       tmpa=arr(J)
       tmpb(:)=brr(:,J)
       DO I=J-1,1,-1
          IF(arr(I)<=tmpa) EXIT
          arr(I+1)=arr(I)
          brr(:,I+1)=brr(:,I)
       ENDDO

       IF(arr(I)<=tmpa)THEN
          arr(I+1)=tmpa
       ELSE
          I=0
          arr(I+1)=tmpa
       ENDIF
       brr(:,I+1)=tmpb(:)
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE csort_eigen
  !--------------------------------------------------------------------------
  SUBROUTINE rsort_eigen(nev,arr,brr)
    !#################################################################!
    ! For a given set of eigenvalues (for all representations), sort  !
    ! them in ascending order. Use straight insertion, no need for    !
    ! fancy algorithms since the arrays are always short.             !
    ! ISBN 7-03-010217-7                                              !
    !#################################################################!
    !
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: nev
    REAL(DP),INTENT(INOUT) :: arr(:)
    REAL(DP),INTENT(INOUT) :: brr(:,:)
    !
    REAL(DP) :: tmpa
    REAL(DP) :: tmpb(SIZE(brr,1))
    INTEGER(I4B) :: I,J
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    DO J=2,nev
       tmpa=arr(J)
       tmpb(:)=brr(:,J)
       DO I=J-1,1,-1
          IF(arr(I)<=tmpa) EXIT
          arr(I+1)=arr(I)
          brr(:,I+1)=brr(:,I)
       ENDDO

       IF(arr(I)<=tmpa)THEN
          arr(I+1)=tmpa
       ELSE
          I=0
          arr(I+1)=tmpa
       ENDIF
       brr(:,I+1)=tmpb(:)
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE rsort_eigen
  !-----------------------------------------------------------------------
  SUBROUTINE sort_eigval(n,arr)
    !
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: n
    REAL(DP),INTENT(INOUT)  :: arr(n)
    !LOCAL
    INTEGER(I4B) :: i,j
    REAL(DP) :: a
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    DO j=2,n

       a=arr(j)

       DO i=j-1,1,-1
          IF(arr(i)<=a) EXIT
          arr(i+1)=arr(i)
       ENDDO

       IF(arr(i)<=a)THEN
          arr(i+1)=a
       ELSE
          i=0
          arr(i+1)=a
       ENDIF

    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE sort_eigval
  !--------------------------------------------------------------------------
  SUBROUTINE realInt_sort(nev,arr,brr,crr)!{{{
    !#################################################################!
    ! For a given set of eigenvalues (for all representations), sort  !
    ! them in ascending order. Use straight insertion, no need for    !
    ! fancy algorithms since the arrays are always short.             !
    ! ISBN 7-03-010217-7                                              !
    !#################################################################!
    !    
    IMPLICIT NONE 
    INTEGER(I4B),INTENT(IN) :: nev
    REAL(DP),INTENT(INOUT) :: arr(:)
    INTEGER(I4B),INTENT(INOUT) :: brr(:,:)
    REAL(DP),OPTIONAL,INTENT(INOUT) :: crr(3,nev)
    !    
    REAL(DP) :: tmpa 
    INTEGER(I4B) :: tmpb(SIZE(brr,1))
    REAL(DP) :: tmpc(3)
    INTEGER(I4B) :: I,J
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    DO J=2,nev
       tmpa=arr(J)
       tmpb(:)=brr(:,J)
       IF(PRESENT(crr)) tmpc(:)=crr(:,J)
       DO I=J-1,1,-1
          IF(arr(I)<=tmpa) EXIT 
          arr(I+1)=arr(I)
          brr(:,I+1)=brr(:,I)
          IF(PRESENT(crr)) crr(:,I+1)=crr(:,I)
          IF(I==1) EXIT 
       ENDDO

       IF(arr(I)<=tmpa)THEN
          arr(I+1)=tmpa
       ELSE 
          I=0  
          arr(I+1)=tmpa
       ENDIF
       brr(:,I+1)=tmpb(:)
       IF(PRESENT(crr)) crr(:,I+1)=tmpc(:)
    ENDDO
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE realInt_sort!}}}

  !-------------------point group function operator-----------------------
  SUBROUTINE PGFO(Omat,n1,n2,n3,Ix,Iy,Iz,Ex,Ey,Ez)
    !
    IMPLICIT NONE
    !INOUT
    REAL(DP),INTENT(IN) :: Omat(3,3)
    INTEGER(I4B),INTENT(IN)  :: n1,n2,n3  & !number of grids in real space
         &    ,         Ix,Iy,Iz    !The point we like
    INTEGER(I4B),INTENT(OUT) :: Ex,Ey,Ez    !the point to read in origin fun
    !LOCAL
    REAL(DP) :: inv_Omat(3,3)
    REAL(DP) :: rin(3),rout(3)
    INTEGER(I4B) :: temp(3)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !real for matrix product
    rin(1)=REAL(Ix,DP)
    rin(2)=REAL(Ix,DP)
    rin(3)=REAL(Ix,DP)
    !inver Omat
    inv_Omat(:,:)=inv_33( Omat(:,:) )
    !operator
    rout=MATMUL(inv_Omat,rin)
    !to grid mesh
    temp(:)=NINT(rout(:))
    !regularlizetion
    Ex=MODULO(temp(1),n1)
    IF(Ex==0) Ex=1

    Ey=MODULO(temp(2),n2)
    IF(Ey==0) Ey=1

    Ez=MODULO(temp(3),n3)
    IF(Ez==0) Ez=1
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE PGFO
  !-----------------------CubicsplineInterpolation-------------------------
  FUNCTION CubicSplineInterp(fun,ddfdx2,xmax,dx,x,Zion)
    !just for uniform grid now , x>0
    !IN/OUT
    REAL(DP),INTENT(IN) :: fun(:)  &
         &, ddfdx2(:) &
         &,xmax , dx , x
    REAL(DP),OPTIONAL,INTENT(IN)  :: Zion
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
         &,  right
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !
    IF(x<0.d0)THEN
       WRITE(6,*) 'CubicSplineInterp: r must > 0'
       STOP
    ENDIF
    !Some special case
    IF(x==0.d0)THEN
       CubicSplineInterp=fun(1)
       RETURN
    ENDIF
    !For r > rmax , fun=0.d0
    IF(x>xmax)THEN
       CubicSplineInterp=0.d0
       RETURN
    ENDIF
    !vlpp
    IF(PRESENT(Zion))THEN
       ! read the following words to see why I do this check here
       ! since NaN is not equal to anyting in FORTRAN
       ! we use the following to see if 4*pi/qNorm is too big to case a NaN
       IF (-4._DP*pi/x**2_DP .NE. -4_DP*pi/x**2_DP ) THEN
          WRITE(6,*)&
               ' There is another case which should be considered: very large bulk. ', &
               '  because larger system will give denser q points and will put q points ', &
               ' very closer to q=0 and this will make -4*pi*Z/q**2 to -infinity',&
               ' If computer find -4*pi*Z/q**2 is too big, it will generate NaN',&
               ' But currently, in our group, we have not meet a system large enough to ',&
               ' cause CPU to generate NaN, if you meet such kind of problem, do something!, code STOP!!'
          STOP
       ENDIF
    ENDIF
    !interpolation:
    !1.Find the left and right point , dt
    pos=x / dx + 1.d0
    !left
    left=FLOOR(pos)
    !right
    right=left+1
    !dt
    dt=pos-REAL(left,DP)
    !2.Evaulate the polynomial
    y_left=fun(left)
    y_right=fun(right)
    ypp_left=ddfdx2(left)
    ypp_right=ddfdx2(right)
    !just do it
    yval = y_left + dt * ( y_right - y_left &
         - ( ypp_right / 6.0_DP + ypp_left / 3.0_DP ) &
         + dt * ( 0.5_DP * ypp_left &
         + dt * ( ( ypp_right - ypp_left ) / 6.0_DP ) ) )
    !output
    IF(PRESENT(Zion))THEN
       CubicSplineInterp=yval - 4.d0*pi*Zion / x**2._DP
    ELSE
       CubicSplineInterp=yval
    ENDIF

    RETURN
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION CubicSplineInterp
  !------------------------------------------------------------------------
  subroutine finite_factor(fnor,norder,coe)!{{{
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
    if(fnor.eq.1) then
       select case (norder)
       case (1)
          coe(1) =  0.50000000000000D+00
       case (2)
          coe(1:2) = [ &
               0.666666666666667      , -8.333333333333333E-002 ]
       case (3)
          coe(1:3) = [ &
               0.750000000000000      , -0.150000000000000      , &
               1.666666666666667E-002 ]
       case (4)
          coe(1:4) = [ &
               0.800000000000000      , -0.200000000000000      , &
               3.809523809523809E-002 , -3.571428571428571E-003 ]
       case (5)
          coe(1:5) = [ &
               0.833333333333333      , -0.238095238095238      , &
               5.952380952380952E-002 , -9.920634920634920E-003 , &
               7.936507936507937E-004 ]
       case (6)
          coe(1:6) = [ &
               0.857142857142857      , -0.267857142857143      , &
               7.936507936507936E-002 , -1.785714285714286E-002 , &
               2.597402597402598E-003 , -1.803751803751804E-004 ]
       case (7)
          coe(1:7) = [ &
               0.875000000000000      , -0.291666666666667      , &
               9.722222222222224E-002 , -2.651515151515151E-002 , &
               5.303030303030303E-003 , -6.798756798756799E-004 , &
               4.162504162504163E-005 ]
       case (8)
          coe(1:8) = [ &
               0.888888888888889      , -0.311111111111111      , &
               0.113131313131313      , -3.535353535353535E-002 , &
               8.702408702408702E-003 , -1.554001554001554E-003 , &
               1.776001776001776E-004 , -9.712509712509713E-006]
       case (9)
          coe(1:9) = [ &
               0.900000000000000      , -0.327272727272727      , &
               0.127272727272727      , -4.405594405594405E-002 , &
               1.258741258741259E-002 , -2.797202797202797E-003 , &
               4.495504495504496E-004 , -4.627725215960510E-005 , &
               2.285296402943462E-006 ]
       case (10)
          coe(1:10) = [ &
               0.909090909090909      , -0.340909090909091      , &
               0.139860139860140      , -5.244755244755244E-002 , &
               1.678321678321678E-002 , -4.370629370629371E-003 , &
               8.814714697067639E-004 , -1.285479226655697E-004 , &
               1.202787580496559E-005 , -5.412544112234515E-007 ]
       end select
       coe(0) = 0.d0
       do i=1,norder
          coe(-i) =-coe(i)
       enddo
    else if (fnor.eq.2) then

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
               -2.72222222222222  , 1.50000000000000       , &
               -0.150000000000000 , 1.111111111111111E-002 ]
       case (4)
          coe(0:4) = [&
               -2.84722222222222  , 1.60000000000000       , &
               -0.200000000000000 , 2.539682539682540E-002 , &
               -1.785714285714286E-003 ]
       case (5)
          coe(0:5) = [ &
               -2.92722222222222       , 1.66666666666667       , &
               -0.238095238095238      , 3.968253968253968E-002 , &
               -4.960317460317460E-003 , 3.174603174603175E-004 ]
       case (6)
          coe(0:6) = [ &
               -2.98277777777778       , 1.71428571428571       , &
               -0.267857142857143      , 5.291005291005291E-002 , &
               -8.928571428571428E-003 , 1.038961038961039E-003 , &
               -6.012506012506013E-005 ]
       case (7)
          coe(0:7) = [ &
               -3.02359410430839       , 1.75000000000000       , &
               -0.291666666666667      , 6.481481481481481E-002 , &
               -1.325757575757576E-002 , 2.121212121212121E-003 , &
               -2.266252266252267E-004 , 1.189286903572618E-005 ]
       case (8)
          coe(0:8) = [ &
               -3.05484410430839       , 1.77777777777778       , &
               -0.311111111111111      , 7.542087542087542E-002 , &
               -1.767676767676768E-002 , 3.480963480963481E-003 , &
               -5.180005180005181E-004 , 5.074290788576503E-005 , &
               -2.428127428127428E-006 ]
       case (9)
          coe(0:9) = [ &
               -3.07953546233308       , 1.80000000000000       , &
               -0.327272727272727      , 8.484848484848484E-002 , &
               -2.202797202797203E-002 , 5.034965034965034E-003 , &
               -9.324009324009329E-004 , 1.284429855858427E-004 , &
               -1.156931303990128E-005 , 5.078436450985472E-007 ]
       case (10)
          coe(0:10) = [ &
               -3.09953546233308       , 1.81818181818182       , &
               -0.340909090909091      , 9.324009324009326E-002 , &
               -2.622377622377622E-002 , 6.713286713286712E-003 , &
               -1.456876456876457E-003 , 2.518489913447896E-004 , &
               -3.213698066639244E-005 , 2.672861289992354E-006 , &
               -1.082508822446903E-007 ]
       end select
       do i=1,norder
          coe(-i) = coe(i)
       enddo
    endif
    !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
  end subroutine finite_factor
  SUBROUTINE finite_factor_new(fnor,norder,coe)!{{{
    !##########################################################
    !* CREATED_TIME  : 2018-09-03
    !* AUTHOR        : Qiang Xu
    !* Ref.          : 'Real-space grid representation of
    ! momentum and kinetic energy operators for energy operators
    ! for electronic structure calculations'
    !*               : Physics meaning
    !##########################################################
    use constants
    !    
    implicit none 
    INTEGER(I4B), intent(in)   :: norder
    INTEGER(I4B), intent(in)   :: fnor 
    real(dp), intent(out) :: coe(-norder:norder)
    INTEGER(I4B)               :: i 
    REAL(DP) :: Omega(norder)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    CALL OmegaMm(Omega)
    IF(fnor==1)THEN
       !Gradient operator
       coe(0)=0._DP
       DO i=1,norder
          coe(i)=1._DP/(2._DP*i*Omega(i))
          coe(-i)=-coe(i)
       ENDDO
    ELSEIF(fnor==2)THEN
       !Laplace
       !coe(i==0)
       coe(0)=0._DP
       DO i=1,norder
          coe(0)=coe(0)-2._DP/(i**2*Omega(i))
       ENDDO
       !coe(i/=0)
       DO i=1,norder
          coe(i)=1._DP/(i**2*Omega(i))
          coe(-i)=coe(i)
       ENDDO
    ENDIF
    !    
  CONTAINS
    !Omega function
    SUBROUTINE OmegaMm(Omeg)
      IMPLICIT NONE 
      REAL(DP),INTENT(OUT) :: Omeg(norder)
      INTEGER(I4B) :: m,l
      REAL(DP) :: tmp
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO m=1,norder
         tmp=1._DP
         DO l=1,norder
            IF(l/=m) tmp=tmp*(1._DP-(REAL(m,DP)/l)**2)
         ENDDO
         Omeg(m)=tmp
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE OmegaMm
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE finite_factor_new!}}}
  !--------------------------dfdr--------------------------------
  SUBROUTINE dfdr(np,h,f,df)
    !just for first order of the don-dimension fun,r : uniform grid
    USE parameters , ONLY : finite_order
    IMPLICIT NONE
    !INOUT
    INTEGER(I4B),INTENT(IN) :: np
    REAL(DP),INTENT(IN)  :: f(np)  &   !f
         & , h  !grid size
    REAL(DP),INTENT(OUT) :: df(np)  !df_dr
    !LOCAL
    REAL(DP) :: coe(-finite_order:finite_order)
    REAL(DP) :: ft(-finite_order:np),tmp
    INTEGER(I4B) :: i,ish
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ft(1:np)=f(1:np)
    !get coe
    CALL finite_factor(1,finite_order,coe)
    coe(:)=coe(:)/h
    !center symmetry of origin set now
    DO i=-finite_order,0
       ft(i)=f(2-i)
    ENDDO
    !finite difference
    DO i=finite_order,np-finite_order
       tmp=0.d0
       DO ish=-finite_order,finite_order,1
          tmp=tmp + coe(ish)*ft(ish+i)
       ENDDO
       df(i) = tmp
    ENDDO
    !first few point
    df(1:finite_order-1)=0.d0
    !Last few points
    df(np-finite_order+1:np)=0.d0
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE dfdr
  !--------------------interploate value-------------------------
  FUNCTION CubicHermiteInterp(fun,dfdx,xmax,h,x)
    !##############################################!
    !y=A0*y0+A1*y1+B0*dy0+B1*dy1                   !
    !defined : dx0=x-x0 and dx1=x-x1 and h=x1-x0   !
    !A0=( 1+2*dx0/h )*(dx1/h)**2                   !
    !A1=( 1-2*dx1/h )*(dx0/h)**2                   !
    !B0=dx0*(dx1/h)**2                             !
    !B1=dx1*(dx0/h)**2                             !
    !##############################################!
    IMPLICIT NONE
    REAL(DP),INTENT(IN) :: fun(:)   &  !in fun
         &,  dfdx(:)  &  !first derivertive
         &,  xmax     &  !max x
         &,  h        &  !grid size
         &,  x           !in x
    REAL(DP) :: CubicHermiteInterp
    !LOCAL
    REAL(DP) :: pos,y0,y1,dy0,dy1
    REAL(DP) :: dx0,dx1,dx0_h,dx1_h,dx0_h2,dx1_h2
    INTEGER(I4B) :: left,right
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IF(x<0.d0)THEN
       WRITE(6,*) 'CubicSplineInterp: r must > 0'
       STOP
    ENDIF
    !Some special case
    IF(x==0.d0)THEN
       CubicHermiteInterp=fun(1)
       RETURN
    ENDIF
    !For r > rmax , fun=0.d0
    IF(x>xmax)THEN
       CubicHermiteInterp=0.d0
       RETURN
    ENDIF
    !interpolation:
    !1.Find the left and right point , dt
    pos=x / h + 1.d0
    !left
    left=FLOOR(pos)
    !right
    right=left+1
    !dxi
    dx0 = x-( left - 1 )*h  !x-x0
    dx1=  dx0-h              !x-x1
    !dxi/h
    dx0_h=dx0/h
    dx1_h=dx1/h
    !(dxi/h)**2
    dx0_h2=dx0_h**2
    dx1_h2=dx1_h**2
    !2.interpolate
    y0=fun(left)
    y1=fun(right)
    dy0=dfdx(left)
    dy1=dfdx(right)
    !
    CubicHermiteInterp= ( y0*(1.d0+2.d0*dx0_h)+dy0*dx0 )*dx1_h2 &
         & +   ( y1*(1.d0-2.d0*dx1_h)+dy1*dx1 )*dx0_h2
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION CubicHermiteInterp
  !--------------------simple interploate------------------------
  FUNCTION SimpleInterp(fun,xmax,h,x)
    IMPLICIT NONE
    REAL(DP),INTENT(IN) :: fun(:)   &  !in fun
         &,  xmax     &  !max x
         &,  h        &  !grid size
         &,  x           !in x
    REAL(DP) :: SimpleInterp
    !LOCAL
    INTEGER(I4B) :: left,right
    REAL(DP) :: y0,y1,dx,pos
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IF(x<0.d0)THEN
       WRITE(6,*) 'CubicSplineInterp: r must > 0'
       STOP
    ENDIF
    !Some special case
    IF(x==0.d0)THEN
       SimpleInterp=fun(1)
       RETURN
    ENDIF
    !For r > rmax , fun=0.d0
    IF(x>xmax)THEN
       SimpleInterp=0.d0
       RETURN
    ENDIF
    !interpolation:
    !1.Find the left and right point , dt
    pos=x / h + 1.d0
    !left
    left=FLOOR(pos)
    !right
    right=left+1
    !dx
    dx= x - (left-1)*h
    !
    y0=fun(left)
    y1=fun(right)
    !2.interpolate
    SimpleInterp=y0+dx*(y1-y0)/h
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION SimpleInterp
  !-----------------------dYlm-------------------------------
  SUBROUTINE r_dYlm(l,m,x,y,z,rmod,f)
    IMPLICIT NONE
    !IN/OUT
    INTEGER(I4B),INTENT(IN) :: l
    INTEGER(I4B),INTENT(IN) :: m
    REAL(DP),INTENT(IN) :: x,y,z,rmod
    REAL(DP),INTENT(INOUT) :: f(0:3)
    !LOCAL
    REAL(DP) :: scal
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IF(.NOT.(rmod>0.d0).AND.l/=0)THEN
       f(:)=0.d0
       RETURN
    ENDIF
    !
    scal= 1.0_dp /SQRT(4.0_dp*pi) !*(-IMAG)**l
    !Apply Y_lm* (spherical harmonic)
    select case(l)
    case(0)
       select case(m)
       case(0)
          f(0) = scal*(1)
          f(1:3)=0
       case default
          WRITE(6,*) 'Apply Ylm:abs(m)>l'
          STOP
       end select
    case(1)
       select case(m)
       case(-1)
          !f = f*scal*(sqrt(3.0_dp)*x)
          f(0) = scal*(sqrt(3.0_dp)*x)
          f(1) = scal*(sqrt(3.0_dp)*(1.d0-x**2.d0)/rmod)
          f(2) = scal*(sqrt(3.0_dp)*(-x*y/rmod))
          f(3) = scal*(sqrt(3.0_dp)*(-x*z/rmod))
       case(0)
          !f = f*scal*(sqrt(3.0_dp)*z)
          f(0) = scal*(sqrt(3.0_dp)*z)
          f(1) = scal*(sqrt(3.0_dp)*(-z*x/rmod))
          f(2) = scal*(sqrt(3.0_dp)*(-z*y/rmod))
          f(3) = scal*(sqrt(3.0_dp)*(1.d0-z**2.d0)/rmod)
       case(1)
          !f = f*scal*(-(sqrt(3.0_dp)*y))
          f(0) = scal*(sqrt(3.0_dp)*y)
          f(1) = scal*(sqrt(3.0_dp)*(-y*x/rmod))
          f(2) = scal*(sqrt(3.0_dp)*(1.d0-y**2.d0)/rmod)
          f(3) = scal*(sqrt(3.0_dp)*(-y*z/rmod))
       case default
          WRITE(6,*) 'Apply Ylm:abs(m)>l'
          STOP
       end select
    case(2)
       select case(m)
       case(-2)
          !f = f*scal*(-(sqrt(15.0_dp)*x*y)) !xx
          f(0) = scal*(-(sqrt(15.0_dp)*x*y))
          f(1) = scal*(-(sqrt(15.0_dp)*y*(1.d0-2.d0*x**2.d0)/rmod))
          f(2) = scal*(-(sqrt(15.0_dp)*x*(1.d0-2.d0*y**2.d0)/rmod))
          f(3) = scal*(-(sqrt(15.0_dp)*(-2.d0*x*y*z)/rmod))
       case(-1)
          !f = f*scal*(sqrt(15.0_dp)*x*z)
          f(0) = scal*((sqrt(15.0_dp)*x*z))
          f(1) = scal*((sqrt(15.0_dp)*z*(1.d0-2.d0*x**2.d0)/rmod))
          f(2) = scal*((sqrt(15.0_dp)*(-2.d0*x*y*z)/rmod))
          f(3) = scal*((sqrt(15.0_dp)*x*(1.d0-2.d0*z**2.d0)/rmod))
       case(0)
          !f = f*scal*((sqrt(5.0_dp)*(-1 + 3*z**2))/2)
          f(0) = scal*((sqrt(5.0_dp)*(-1.d0 + 3.d0*z**2.d0))/2.d0)
          f(1) = scal*sqrt(5.0_dp)*3.d0*(-x*z**2.d0/rmod)
          f(2) = scal*(sqrt(5.0_dp)*(3.d0*(-y*z**2.d0/rmod)))
          f(3) = scal*((sqrt(5.0_dp)*(3.d0*z*(1.d0-z**2.d0)/rmod)))
       case(1)
          !f = f*scal*(-(sqrt(15.0_dp)*y*z)) !xx
          f(0) = scal*(-(sqrt(15.0_dp)*y*z))
          f(1) = scal*(-(sqrt(15.0_dp)*(-2.d0*x*y*z)/rmod))
          f(2) = scal*(-(sqrt(15.0_dp)*z*(1.d0-2.d0*y**2.d0)/rmod))
          f(3) = scal*(-(sqrt(15.0_dp)*y*(1.d0-2.d0*z**2.d0)/rmod))
       case(2)
          !f = f*scal*((sqrt(15.0_dp)*(-x**2 + y**2))/2)
          f(0) = scal*(sqrt(15.0_dp)*(-x**2.d0 + y**2.d0)/2.d0)
          f(1) = scal*(sqrt(15.0_dp)*(-x*(1.d0-x**2.d0+y**2.d0))/rmod)
          f(2) = scal*(sqrt(15.0_dp)*y*(1.d0+x**2.d0-y**2.d0)/rmod)
          f(3) = -scal*(sqrt(15.0_dp)*z*(y**2.d0-x**2.d0)/rmod)
       case default
          WRITE(6,*) 'Apply Ylm:abs(m)>l'
          STOP
       end select
    case(3)
       select case(m)
       case(-3)
          !f = f*scal*((sqrt(17.5_dp)*x*(-1 + 4*y**2 + z**2))/2)
          f(0) = scal*((sqrt(17.5_dp)*x*(-1 + 4*y**2 + z**2))/2)
          f(1) = scal*(sqrt(17.5_dp)*((y**2+z**2)/rmod*(-1 + 4*y**2 + z**2)+x*(4*(-2*x*y**2)+(-2*x*z**2))/rmod)/2)
          f(2) = scal*((sqrt(17.5_dp)*((-x*y/rmod)*(-1 + 4*y**2 + z**2)+x*(4*2*y*(x**2+z**2)+(-2*y*z**2))/rmod))/2)
          f(3) = scal*((sqrt(17.5_dp)*((-x*z/rmod)*(-1 + 4*y**2 + z**2)+x*(4*2*z*(-2*z*y**2)+(2*z*(x**2+y**2)))/rmod))/2)
       case(-2)
          !f = f*scal*(-(sqrt(105.0_dp)*x*y*z))
          f(0) = scal*(-(sqrt(105.0_dp)*x*y*z))
          f(1) = scal*(-(sqrt(105.0_dp)*(y*z*(y**2+z**2)-2*x**2*y*z))/rmod)
          f(2) = scal*(-(sqrt(105.0_dp)*(x*z*(x**2+z**2)-2*x*y**2*z))/rmod)
          f(3) = scal*(-(sqrt(105.0_dp)*(x*y*(x**2+y**2)-2*x*y*z**2))/rmod)
       case(-1)
          !f = f*scal*((sqrt(10.5_dp)*x*(-1 + 5*z**2))/2)
          f(0) = scal*((sqrt(10.5_dp)*x*(-1 + 5*z**2))/2)
          f(1) = scal*((sqrt(10.5_dp)*((y**2+z**2)/rmod*(-1 + 5*z**2)+x*5*(-2*x*z**2)/rmod))/2)
          f(2) = scal*((sqrt(10.5_dp)*((-x*y/rmod)*(-1 + 5*z**2)+x*5*(-2*y*z**2)/rmod))/2)
          f(3) = scal*((sqrt(10.5_dp)*((-x*z/rmod)*(-1 + 5*z**2)+x*5*(2*z*(x**2+y**2))/rmod))/2)
       case(0)
          !f = f*scal*((sqrt(7.0_dp)*z*(-3 + 5*z**2))/2)
          f(0) = scal*((sqrt(7.0_dp)*z*(-3 + 5*z**2))/2)
          f(1) = scal*((sqrt(7.0_dp)*((-z*x/rmod)*(-3 + 5*z**2)+z*5*(-2*x*z**2)/rmod))/2)
          f(2) = scal*((sqrt(7.0_dp)*((-z*y/rmod)*(-3 + 5*z**2)+z*5*(-2*y*z**2)/rmod))/2)
          f(3) = scal*((sqrt(7.0_dp)*((x**2+y**2)/rmod*(-3 + 5*z**2)+z*5*(2*z*(x**2+y**2))/rmod))/2)
       case(1)
          !f = f*scal*((sqrt(10.5_dp)*y*(1 - 5*z**2))/2)
          f(0) = scal*((sqrt(10.5_dp)*y*(1 - 5*z**2))/2)
          f(1) = scal*((sqrt(10.5_dp)*((-y*x/rmod)*(-3 + 5*z**2)+y*5*(-2*x*z**2)/rmod))/2)
          f(2) = scal*((sqrt(10.5_dp)*((x**2+z**2)/rmod*(-3 + 5*z**2)+y*5*(-2*y*z**2)/rmod))/2)
          f(3) = scal*((sqrt(10.5_dp)*((-y*z/rmod)*(-3 + 5*z**2)+y*5*(2*z*(x**2+y**2))/rmod))/2)
       case(2)
          !f = f*scal*((sqrt(105.0_dp)*(-x**2 + y**2)*z)/2)
          f(0) = scal*((sqrt(105.0_dp)*(-x**2 + y**2)*z)/2)
          f(1) = scal*((sqrt(105.0_dp)*((-z*x/rmod)*(-x**2 + y**2)+z*(-2*x*(y**2+z**2)+(-2*x*y**2))/rmod))/2)
          f(2) = scal*((sqrt(105.0_dp)*((-z*y/rmod)*(-x**2 + y**2)+z*(-(-2*y*x**2)+2*y*(x**2+z**2))/rmod))/2)
          f(3) = scal*((sqrt(105.0_dp)*((x**2+z**2)/rmod*(-x**2 + y**2)+z*(-(-2*z*x**2)+(-2*z*y**2))/rmod))/2)
       case(3)
          !f = f*scal*(-(sqrt(17.5_dp)*y*(-3*x**2 + y**2))/2)
          f(0) = scal*(-(sqrt(17.5_dp)*y*(-3*x**2 + y**2))/2)
          f(1) = scal*(-(sqrt(17.5_dp)*((-y*x/rmod)*(-3*x**2 + y**2)+y*(-3*2*x*(1-x**2)+(-2*x*y**2))/rmod))/2)
          f(2) = scal*(-(sqrt(17.5_dp)*((1-y**2)/rmod*(-3*x**2 + y**2)+y*(-3*(-2*y*x**2)+(2*y*(1-y**2)))/rmod))/2)
          f(3) = scal*(-(sqrt(17.5_dp)*((-y*z/rmod)*(-3*x**2 +  y**2)+y*(-3*(-2*z*x**2)+(-2*z*y**2))/rmod))/2)
       case default
          WRITE(6,*) 'Apply Ylm:abs(m)>l'
          STOP
       end select
    case default
       WRITE(6,*) 'Apply Ylm:l>3 not programmed'
       STOP
    end select
    return
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE r_dYlm
  !-------------------------STO------------------------------
  ! SUBROUTINE atom_STO(p,l,m,zeta,r,f)
  ! USE parameters , ONLY : weight=>Wexict
  ! !for calculate the slate-obtial
  ! IMPLICIT NONE
  ! INTEGER(I4B),INTENT(IN) :: p,l,m !quantum number
  ! REAL(DP),INTENT(IN) :: zeta & !
  !                     &, r(:)  !rvec
  ! REAL(DP),INTENT(OUT) :: f
  ! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! IF ( p==0 ) THEN
  !     f=0.d0
  !     RETURN
  ! ENDIF

  ! select case(p)
  ! case (1)
  !     select case (l)
  !     case (0)
  !             IF ( m==0 ) THEN
  !                     f=sqrt(zeta**3/pi)*exp(-zeta*r(4))
  !             ELSE
  !                     stop "m is greater than l or less than -l"
  !             ENDIF
  !     case default
  !             stop "l is greater than n-1 or less than 0"
  !     end select
  ! case (2)
  !     select case (l)
  !     case (0)
  !             IF ( m==0 ) THEN
  !                     f=sqrt(zeta**5/(3*pi))*r(4)*exp(-zeta*r(4))
  !             ELSE
  !                     stop "m is greater than l or less than -l"
  !             ENDIF
  !     case (1)
  !             select case (m)
  !             case (1)
  !                     f=sqrt(zeta**5/pi)*r(1)*exp(-zeta*r(4))
  !             case (0)
  !                     f=sqrt(zeta**5/pi)*r(2)*exp(-zeta*r(4))
  !             case (-1)
  !                     f=sqrt(zeta**5/pi)*r(3)*exp(-zeta*r(4))
  !             case default
  !                     stop "m is greater than l or less than -l"
  !             end select
  !     case default
  !             stop "l is greater than n-1 or less than 0"
  !     end select
  ! case (3)
  !     select case (l)
  !     case (0)
  !             IF ( m==0 ) THEN
  !                     f=sqrt((2*zeta**7)/(45*pi))*r(4)**2*exp(-zeta*r(4))
  !             ELSE
  !                     stop "m is greater than l or less than -l"
  !             ENDIF
  !     case (1)
  !             select case (m)
  !             case (1)
  !                     f=sqrt(2*zeta**7/(15*pi))*r(4)*r(1)*exp(-zeta*r(4))
  !             case (0)
  !                     f=sqrt(2*zeta**7/(15*pi))*r(4)*r(2)*exp(-zeta*r(4))
  !             case (-1)
  !                     f=sqrt(2*zeta**7/(15*pi))*r(4)*r(3)*exp(-zeta*r(4))
  !             case default
  !                     stop "m is greater than l or less than -l"
  !             end select
  !     case (2)
  !             select case (m)
  !             case (2)
  !                     f=(1.d0/3)*sqrt(zeta**7/(2*pi))*(3*r(3)**2-r(4)**2)*exp(-zeta*r(4))
  !             case (1)
  !                     f=sqrt(2*zeta**7/(3*pi))*r(1)*r(3)*exp(-zeta*r(4))
  !             case (0)
  !                     f=sqrt(2*zeta**7/(3*pi))*r(2)*r(3)*exp(-zeta*r(4))
  !             case (-1)
  !                     f=sqrt(zeta**7/(6*pi))*(r(1)**2-r(2)**2)*exp(-zeta*r(4))
  !             case (-2)
  !                     f=sqrt(2*zeta**7/(3*pi))*r(1)*r(2)*exp(-zeta*r(4))
  !             case default
  !                     stop "m is greater than l or less than -l"
  !             end select
  !     case default
  !             stop "l is greater than n-1 or less than 0"
  !     end select
  ! case default
  !     stop "principl quantum number is greater than 3 or less than 1"
  ! end select

  ! !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! ENDSUBROUTINE atom_STO
  !----------------------zeta data---------------------------
  ! SUBROUTINE atom_effcharge(atom_name,Lmax,nquan,zeta)!{{{
  !    implicit none
  !    !--------------------------------------------
  !    character(len=*)     :: atom_name
  !    integer              :: i
  !    integer,INTENT(OUT)  :: Lmax  &
  !                       &,nquan(3)        !only three dimensions
  !    REAL(DP),INTENT(OUT) :: zeta(3)              !only three dimensions
  !    !------------------------------------------------
  !    type :: elements
  !       integer          ::atom_number
  !       character(len=3) :: symbols
  !       integer(I4B)     ::Lmax,nspd(3)
  !       real(DP)         ::zeta(3)
  !    end type elements
  !    !
  !    type(elements) ::atom(50)
  !    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !    !print *,"name",atom_name
  !    do i=1,50
  !       if(adjustl(atom_name)==atom(i)%symbols) then
  !          Lmax=atom(i)%Lmax
  !          nquan=atom(i)%nspd
  !          zeta=atom(i)%zeta
  !       endif
  !    enddo
  !    data atom(1:50) /&
  !       & elements(  1, "H",  0,(/1,0,0/),(/1.000000000000000d0,0.000000000000000d0,0.0d0/)),&
  !       & elements(  2, "He", 0,(/1,0,0/),(/1.700000000000000d0,0.000000000000000d0,0.0d0/)),&
  !       & elements(  3, "Li", 1,(/2,2,0/),(/0.650000000000000d0,0.650000000000000d0,0.0d0/)),&
  !       & elements(  4, "Be", 1,(/2,2,0/),(/0.975000000000000d0,0.975000000000000d0,0.0d0/)),&
  !       & elements(  5, "B",  1,(/2,2,0/),(/1.300000000000000d0,1.300000000000000d0,0.0d0/)),&
  !       & elements(  6, "C",  1,(/2,2,0/),(/1.625000000000000d0,1.625000000000000d0,0.0d0/)),&
  !       & elements(  7, "N",  1,(/2,2,0/),(/1.950000000000000d0,1.950000000000000d0,0.0d0/)),&
  !       & elements(  8, "O",  1,(/2,2,0/),(/2.275000000000000d0,2.275000000000000d0,0.0d0/)),&
  !       & elements(  9, "F",  1,(/2,2,0/),(/2.600000000000000d0,2.600000000000000d0,0.0d0/)),&
  !       & elements( 10, "Ne", 1,(/2,2,0/),(/2.925000000000000d0,2.925000000000000d0,0.0d0/)),&
  !       & elements( 11, "Na", 1,(/3,3,0/),(/0.733333333333333d0,0.733333333333333d0,0.0d0/)),&
  !       & elements( 12, "Mg", 1,(/3,3,0/),(/0.950000000000000d0,0.950000000000000d0,0.0d0/)),&
  !       & elements( 13, "Al", 1,(/3,3,0/),(/1.166666666666666d0,1.166666666666666d0,0.0d0/)),&
  !       & elements( 14, "Si", 1,(/3,3,0/),(/1.383333333333333d0,1.383333333333333d0,0.0d0/)),&
  !       & elements( 15, "P",  2,(/3,3,3/),(/1.600000000000000d0,1.600000000000000d0,0.533333333333333d0/)),&
  !       & elements( 16, "S",  2,(/3,3,3/),(/1.816666666666666d0,1.816666666666666d0,0.583333333333333d0/)),&
  !       & elements( 17, "Cl", 2,(/3,3,3/),(/2.033333333333333d0,2.033333333333333d0,0.633333333333333d0/)),&
  !       & elements( 18, "Ar", 2,(/3,3,3/),(/2.250000000000000d0,2.250000000000000d0,0.683333333333333d0/)),&
  !       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ADD>NEW>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !       & elements( 19, "K",  1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 20, "Ca", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 21, "Sc", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 22, "Ti", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 23, "V",  1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 24, "Cr", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 25, "Mn", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 26, "Fe", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 27, "Co", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 28, "Ni", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 29, "Cu", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 30, "Zn", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 31, "Ga", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 32, "Ge", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 33, "As", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 34, "Se", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 35, "Br", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 36, "Kr", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 37, "Rb", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 38, "Sr", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 39, "Y",  1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 40, "Zr", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 41, "Nb", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 42, "Mo", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 43, "Tc", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 44, "Ru", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 45, "Rh", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 46, "Pd", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 47, "Ag", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 48, "Cd", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 49, "In", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/)),&
  !       & elements( 50, "Sn", 1,(/3,3,0/),(/1.816666666666666d0,0.0d0,0.0d0/))/
  !    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! ENDSUBROUTINE
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  SUBROUTINE atom_effcharge(atom_name,Lmax,nquan,zeta)!{{{
    implicit none
    !--------------------------------------------
    character(len=*)     :: atom_name
    integer              :: i
    integer,INTENT(OUT)  :: Lmax  &
         &,nquan(4)        !only three dimensions
    REAL(DP),INTENT(OUT) :: zeta(4)              !only three dimensions
    !------------------------------------------------
    type :: elements
       integer          ::atom_number
       character(len=3) :: symbols
       integer(I4B)     ::Lmax,nspd(4)
       real(DP)         ::zeta(4)
    end type elements
    !
    type(elements) ::atom(55)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !print *,"name",atom_name
    data atom(1:55) /&
         & elements(  1, "H",0,(/1,0,0,0/),(/1.000000000000000d0,0.000000000000000d0,0.0d0,0.0d0/)),&
         & elements(  2, "He",0,(/1,0,0,0/),(/1.700000000000000d0,0.000000000000000d0,0.0d0,0.0d0/)),&
         & elements(  3, "Li",1,(/2,2,0,0/),(/0.650000000000000d0,0.650000000000000d0,0.0d0,0.0d0/)),&
         & elements(  4, "Be",1,(/2,2,0,0/),(/0.975000000000000d0,0.975000000000000d0,0.0d0,0.0d0/)),&
         & elements(  5, "B",1,(/2,2,0,0/),(/1.300000000000000d0,1.300000000000000d0,0.0d0,0.0d0/)),&
         & elements(  6, "C",1,(/2,2,0,0/),(/1.625000000000000d0,1.625000000000000d0,0.0d0,0.0d0/)),&
         & elements(  7, "N",1,(/2,2,0,0/),(/1.950000000000000d0,1.950000000000000d0,0.0d0,0.0d0/)),&
         & elements(  8, "O",1,(/2,2,0,0/),(/2.275000000000000d0,2.275000000000000d0,0.0d0,0.0d0/)),&
         & elements(  9, "F",1,(/2,2,0,0/),(/2.600000000000000d0,2.600000000000000d0,0.0d0,0.0d0/)),&
         & elements( 10, "Ne",1,(/2,2,0,0/),(/2.925000000000000d0,2.925000000000000d0,0.0d0,0.0d0/)),&
         & elements( 11, "Na",1,(/3,3,0,0/),(/0.733333333333333d0,0.733333333333333d0,0.0d0,0.0d0/)),&
         & elements( 12, "Mg",1,(/3,3,0,0/),(/0.950000000000000d0,0.950000000000000d0,0.0d0,0.0d0/)),&
         & elements( 13, "Al",1,(/3,3,0,0/),(/1.166666666666666d0,1.166666666666666d0,0.0d0,0.0d0/)),&
         & elements( 14, "Si",1,(/3,3,0,0/),(/1.383333333333333d0,1.383333333333333d0,0.0d0,0.0d0/)),&
         & elements( 15, "P",2,(/3,3,3,0/),(/1.600000000000000d0,1.600000000000000d0,0.533333333333333d0,0.0d0/)),&
         & elements( 16, "S",2,(/3,3,3,0/),(/1.816666666666666d0,1.816666666666666d0,0.583333333333333d0,0.0d0/)),&
         & elements( 17, "Cl",2,(/3,3,3,0/),(/2.033333333333333d0,2.033333333333333d0,0.633333333333333d0,0.0d0/)),&
         & elements( 18, "Ar",2,(/3,3,3,0/),(/2.250000000000000d0,2.250000000000000d0,0.683333333333333d0,0.0d0/)),&
                                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ADD>NEW>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         & elements( 19, "K",1,(/4,4,0,0/),(/0.594594594594595d0,0.594594594594595d0,0.000000000000000d0,0.000000000000000d0/)),&
         & elements( 20, "Ca",1,(/4,4,0,0/),(/0.770270270270270d0,0.770270270270270d0,0.000000000000000d0,0.000000000000000d0/)),&
         & elements( 21, "Sc",2,(/4,4,3,0/),(/0.486486486486486d0,0.486486486486486d0,1.400000000000000d0,0.000000000000000d0/)),&
                                !>orbital 3d under 4s and upper
         & elements( 22, "Ti",2,(/4,4,3,0/),(/0.527027027027027d0,0.527027027027027d0,1.616666666666670d0,0.000000000000000d0/)),&
                                !>than 3s3p
         & elements( 23, "V",2,(/4,4,3,0/),(/0.567567567567568d0,0.567567567567568d0,1.833333333333330d0,0.000000000000000d0/)),&
         & elements( 24, "Cr",2,(/4,4,3,0/),(/0.472972972972973d0,0.472972972972973d0,2.266666666666670d0,0.000000000000000d0/)),&
         & elements( 25, "Mn",2,(/4,4,3,0/),(/0.648648648648649d0,0.648648648648649d0,2.266666666666670d0,0.000000000000000d0/)),&
         & elements( 26, "Fe",2,(/4,4,3,0/),(/0.689189189189189d0,0.689189189189189d0,2.483333333333330d0,0.000000000000000d0/)),&
         & elements( 27, "Co",2,(/4,4,3,0/),(/0.729729729729730d0,0.729729729729730d0,2.700000000000000d0,0.000000000000000d0/)),&
         & elements( 28, "Ni",2,(/4,4,3,0/),(/0.770270270270270d0,0.770270270270270d0,2.916666666666670d0,0.000000000000000d0/)),&
         & elements( 29, "Cu",2,(/4,4,3,0/),(/0.675675675675676d0,0.675675675675676d0,3.350000000000000d0,0.000000000000000d0/)),&
         & elements( 30, "Zn",2,(/4,4,3,0/),(/0.851351351351351d0,0.851351351351351d0,3.566666666666670d0,0.000000000000000d0/)),&
         & elements( 31, "Ga",2,(/4,4,3,0/),(/1.027027027027030d0,1.027027027027030d0,3.900000000000000d0,0.000000000000000d0/)),&
         & elements( 32, "Ge",2,(/4,4,3,0/),(/1.202702702702700d0,1.202702702702700d0,4.233333333333330d0,0.000000000000000d0/)),&
         & elements( 33, "As",2,(/4,4,4,0/),(/1.378378378378380d0,1.378378378378380d0,0.432432432432432d0,0.000000000000000d0/)),&
                                !>due to the 4p has been occupied
         & elements( 34, "Se",2,(/4,4,4,0/),(/1.554054054054050d0,1.554054054054050d0,0.472972972972973d0,0.000000000000000d0/)),&
                                !>3d replaced by 4d
         & elements( 35, "Br",2,(/4,4,4,0/),(/1.729729729729730d0,1.729729729729730d0,0.513513513513513d0,0.000000000000000d0/)),&
                                !!
         & elements( 36, "Kr",2,(/4,4,4,0/),(/1.905405405405410d0,1.905405405405410d0,0.554054054054054d0,0.000000000000000d0/)),&
                                !!
         & elements( 37, "Rb",1,(/5,5,0,0/),(/0.550000000000000d0,0.550000000000000d0,0.000000000000000d0,0.000000000000000d0/)),&
         & elements( 38, "Sr",1,(/5,5,0,0/),(/0.712500000000000d0,0.712500000000000d0,0.000000000000000d0,0.000000000000000d0/)),&
         & elements( 39, "Y",2,(/5,5,4,0/),(/0.450000000000000d0,0.450000000000000d0,1.135135135135140d0,0.000000000000000d0/)),&
         & elements( 40, "Zr",2,(/5,5,4,0/),(/0.487500000000000d0,0.487500000000000d0,1.310810810810810d0,0.000000000000000d0/)),&
         & elements( 41, "Nb",2,(/5,5,4,0/),(/0.400000000000000d0,0.400000000000000d0,1.391891891891890d0,0.000000000000000d0/)),&
         & elements( 42, "Mo",3,(/5,5,4,4/),(/0.437500000000000d0,0.437500000000000d0,1.567567567567570d0,0.650000000000000d0/)),&
         & elements( 43, "Tc",3,(/5,5,4,4/),(/0.600000000000000d0,0.600000000000000d0,1.837837837837840d0,0.900000000000000d0/)),&
         & elements( 44, "Ru",3,(/5,5,4,4/),(/0.512500000000000d0,0.512500000000000d0,1.918918918918920d0,0.725000000000000d0/)),&
         & elements( 45, "Rh",3,(/5,5,4,4/),(/0.550000000000000d0,0.550000000000000d0,2.094594594594590d0,0.762500000000000d0/)),&
         & elements( 46, "Pd",3,(/5,5,4,4/),(/0.587500000000000d0,0.587500000000000d0,2.175675675675680d0,0.587500000000000d0/)),&
                                !>valence electrons is 4d10,transition occured for 5s5p such as format is
                                !uniform
         & elements( 47, "Ag",3,(/5,5,4,4/),(/0.625000000000000d0,0.625000000000000d0,2.445945945945950d0,0.837500000000000d0/)),&
         & elements( 48, "Cd",3,(/5,5,4,4/),(/0.787500000000000d0,0.787500000000000d0,2.716216216216220d0,1.087500000000000d0/)),&
         & elements( 49, "In",3,(/5,5,4,4/),(/0.950000000000000d0,0.950000000000000d0,2.986486486486490d0,1.337500000000000d0/)),&
         & elements( 50, "Sn",3,(/5,5,4,4/),(/1.112500000000000d0,1.112500000000000d0,3.256756756756760d0,1.587500000000000d0/)),&
         & elements( 51, "Sb",3,(/5,5,5,4/),(/1.275000000000000d0,1.275000000000000d0,0.400000000000000d0,1.837500000000000d0/)),&
                                !>for the same reason the
         & elements( 52, "Te",3,(/5,5,5,4/),(/1.437500000000000d0,1.437500000000000d0,0.437500000000000d0,2.087500000000000d0/)),&
                                !>4d is replaced by 5d,but
         & elements( 53, "I",3,(/5,5,5,4/),(/1.600000000000000d0,1.600000000000000d0,0.475000000000000d0,2.337500000000000d0/)),&
                                !>4f residualed reasonally
         & elements( 79, "Au",3,(/5,5,4,4/),(/0.625000000000000d0,0.625000000000000d0,2.445945945945950d0,0.837500000000000d0/)),&
         & elements( 54, "Xe",3,(/5,5,5,4/),(/1.762500000000000d0,1.762500000000000d0,0.512500000000000d0,2.587500000000000d0/))/
    !print *,"name",atom_name
    do i=1,55
       if(adjustl(atom_name)==atom(i)%symbols) then
          Lmax=atom(i)%Lmax
          nquan=atom(i)%nspd
          zeta=atom(i)%zeta
       endif
    enddo
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  end subroutine atom_effcharge !}}}

  ! SUBROUTINE atom_effcharge_cal()(atom_name,Lmax,nquan,zeta)!{{{
  ! implicit none
  ! !> IN/OUT
  ! character(len=*)     :: atom_name
  ! integer              :: i
  ! integer,INTENT(OUT)  :: Lmax  &
  !                     &,nquan(4)        !only three dimensions
  ! REAL(DP),INTENT(OUT) :: zeta(4)              !only three dimensions
  ! !> local
  ! type :: elements
  !     integer          ::atom_number
  !     character(len=3) :: symbols
  ! end type elements
  ! !
  ! type(elements) ::atom(104)
  ! INTEGER(I4B)   :: nu_atom,nu_valence
  ! REAL(DP)   :: zeta_s,zeta_p,zeta_d,zeta_f
  ! !> ----------------------------------
  ! data atom(1:104)/&
  !     & elements(  1, "H" ),&
  !     & elements(  2, "He"),&
  !     & elements(  3, "Li"),&
  !     & elements(  4, "Be"),&
  !     & elements(  5, "B" ),&
  !     & elements(  6, "C" ),&
  !     & elements(  7, "N" ),&
  !     & elements(  8, "O" ),&
  !     & elements(  9, "F" ),&
  !     & elements( 10, "Ne"),&
  !     & elements( 11, "Na"),&
  !     & elements( 12, "Mg"),&
  !     & elements( 13, "Al"),&
  !     & elements( 14, "Si"),&
  !     & elements( 15, "P" ),&
  !     & elements( 16, "S" ),&
  !     & elements( 17, "Cl"),&
  !     & elements( 18, "Ar"),&
  !     & elements( 19, "K" ),&
  !     & elements( 20, "Ca"),&
  !     & elements( 21, "Sc"),&
  !     & elements( 22, "Ti"),&
  !     & elements( 23, "V" ),&
  !     & elements( 24, "Cr"),&
  !     & elements( 25, "Mn"),&
  !     & elements( 26, "Fe"),&
  !     & elements( 27, "Co"),&
  !     & elements( 28, "Ni"),&
  !     & elements( 29, "Cu"),&
  !     & elements( 30, "Zn"),&
  !     & elements( 31, "Ga"),&
  !     & elements( 32, "Ge"),&
  !     & elements( 33, "As"),&
  !     & elements( 34, "Se"),&
  !     & elements( 35, "Br"),&
  !     & elements( 36, "Kr"),&
  !     & elements( 37, "Rb"),&
  !     & elements( 38, "Sr"),&
  !     & elements( 39, "Y" ),&
  !     & elements( 40, "Zr"),&
  !     & elements( 41, "Nb"),&
  !     & elements( 42, "Mo"),&
  !     & elements( 43, "Tc"),&
  !     & elements( 44, "Ru"),&
  !     & elements( 45, "Rh"),&
  !     & elements( 46, "Pd"),&
  !     & elements( 47, "Ag"),&
  !     & elements( 48, "Cd"),&
  !     & elements( 49, "In"),&
  !     & elements( 50, "Sn"),&
  !     & elements( 51, "Sb"),&
  !     & elements( 52, "Te"),&
  !     & elements( 53, "I"),&
  !     & elements( 54, "Xe"),&
  !     & elements( 55, "Cs"),&
  !     & elements( 56, "Ba"),&
  !     & elements( 57, "La"),&
  !     & elements( 58, "Ce"),&
  !     & elements( 59, "Pr"),&
  !     & elements( 60, "Nd"),&
  !     & elements( 61, "Pm"),&
  !     & elements( 62, "Sm"),&
  !     & elements( 63, "Eu"),&
  !     & elements( 64, "Gd"),&
  !     & elements( 65, "Tb"),&
  !     & elements( 66, "Dy"),&
  !     & elements( 67, "Ho"),&
  !     & elements( 68, "Er"),&
  !     & elements( 69, "Tm"),&
  !     & elements( 70, "Yb"),&
  !     & elements( 71, "Lu"),&
  !     & elements( 72, "Hf"),&
  !     & elements( 73, "Ta"),&
  !     & elements( 74, "W"),&
  !     & elements( 75, "Re"),&
  !     & elements( 76, "Os"),&
  !     & elements( 77, "Ir"),&
  !     & elements( 78, "Pt"),&
  !     & elements( 79, "Au"),&
  !     & elements( 80, "Hg"),&
  !     & elements( 81, "Ti"),&
  !     & elements( 82, "Pb"),&
  !     & elements( 83, "Bi"),&
  !     & elements( 84, "Po"),&
  !     & elements( 85, "At"),&
  !     & elements( 86, "Rn"),&
  !     & elements( 87, "Fr"),&
  !     & elements( 88, "Ra"),&
  !     & elements( 89, "Ac"),&
  !     & elements( 90, "Th"),&
  !     & elements( 91, "Pa"),&
  !     & elements( 92, "U"),&
  !     & elements( 93, "Np"),&
  !     & elements( 94, "Pu"),&
  !     & elements( 95, "Am"),&
  !     & elements( 96, "Cm"),&
  !     & elements( 97, "Bk"),&
  !     & elements( 98, "Cf"),&
  !     & elements( 99, "Es"),&
  !     & elements(100, "Fm"),&
  !     & elements(101, "Md"),&
  !     & elements(102, "No"),&
  !     & elements(103, "Lr"),&
  !     & elements(104, "Rf")/
  ! !> get the atom number and effcharge
  ! do i=1,104
  !    if(adjustl(atom_name)==atom(i)%symbols) then
  !       nu_atom=atom(i)%atom_number
  !    endif
  ! enddo

  ! select case(nu_atom)
  ! case(1:2)
  !    !> n=1
  !    !> number of valence electrons
  !    nu_valence=nu_atom-0
  !    Lmax=0
  !    nquan=(/1,0,0,0/)
  !    zeta=(/real(nu_atom,DP)-(nu_valence-1)*0.3d0,0,0,0/)/1.d0
  ! case(3:10)
  !    !> n=2
  !    !> number of valence electrons
  !    nu_valence=nu_atom-2
  !    select case(nu_valence)
  !    case(:2) !> 2s
  !       Lmax=0
  !       nquan=(/2,0,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-1)*0.35d0-2*0.85
  !       zeta_p=0.d0
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    case default !> 2p
  !       Lmax=1
  !       nquan=(/2,2,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-1)*0.35d0-2*0.85
  !       zeta_p=zeta_s
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    end select
  ! case(11:18)
  !    !> n=3
  !    !> number of valence electrons
  !    nu_valence=nu_atom-10
  !    select case(nu_valence)
  !    case(:2) !> 3s
  !       Lmax=0
  !       nquan=(/3,0,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-1)*0.35d0-8*0.85-2*1
  !       zeta_p=0.d0
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    case default !> 3p
  !       Lmax=1
  !       nquan=(/3,3,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-1)*0.35d0-8*0.85-2*1
  !       zeta_p=zeta_s
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    end select
  ! case(19:36)
  !    !> n=4
  !    !> number of valence electrons
  !    nu_valence=nu_atom-18
  !    select case(nu_valence)
  !    case(1:2) !> 4s
  !       Lmax=0
  !       nquan=(/4,0,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-1)*0.35d0-8*0.85-(2+8)*1
  !       zeta_p=0.d0
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    case (3:12) !> 3d
  !       Lmax=2
  !       nquan=(/4,4,3,0/)
  !       zeta_s=real(nu_atom,DP)-(2-1)*0.35d0-(8+nu_valence-2)*0.85-(2+8)*1
  !       zeta_p=zeta_s
  !       zeta_d=real(nu_atom,DP)-(nu_valence-2)*0.35d0-(2+8+8)*1
  !       zeta_f=0.d0
  !       if(nu_atom==29)then !> 3d10 4s1
  !          zeta_s=real(nu_atom,DP)-(1-1)*0.35d0-18*0.85-(2+8)*1
  !          zeta_p=zeta_s
  !          zeta_d=real(nu_atom,DP)-9*0.35d0-(2+8+8)*1
  !          zeta_f=0.d0
  !       endif
  !       if(nu_atom==24)then !> 3d5 4s1
  !          zeta_s=real(nu_atom,DP)-(1-1)*0.35d0-13*0.85-(2+8)*1
  !          zeta_p=zeta_s
  !          zeta_d=real(nu_atom,DP)-4*0.35d0-(2+8+8)*1
  !          zeta_f=0.d0
  !       endif
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    case(13:18)
  !       Lmax=1
  !       nquan=(/4,4,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-10)*0.35d0-(8+10)*0.85-(2+8)*1
  !       zeta_p=zeta_s
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/3.7d0
  !    end select
  ! case(37:54)
  !    !> n=5
  !    !> number of valence electrons
  !    nu_valence=nu_atom-36
  !    select case(nu_valence)
  !    case(1:2) !> 5s
  !       Lmax=0
  !       nquan=(/5,0,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-1)*0.35d0-8*0.85-(2+8+18)*1
  !       zeta_p=0.d0
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    case (3:12) !> 4d
  !       Lmax=2
  !       nquan=(/5,5,3,0/)
  !       zeta_s=real(nu_atom,DP)-(2-1)*0.35-(8+nu_valence-2)*0.85-(2+8+18)*1
  !       zeta_p=zeta_s
  !       zeta_d=real(nu_atom,DP)-(nu_valence-2)*0.35d0-(2+8+18+8)*1
  !       zeta_f=0.d0
  !       if(nu_atom==41)then !> 3d10 4s1
  !          zeta_s=real(nu_atom,DP)-(1-1)*0.35d0-18*0.85-(2+8)*1
  !          zeta_p=zeta_s
  !          zeta_d=real(nu_atom,DP)-9*0.35d0-(2+8+8)*1
  !          zeta_f=0.d0
  !       endif
  !       if(nu_atom==42)then !> 3d10 4s1
  !          zeta_s=real(nu_atom,DP)-(1-1)*0.35d0-18*0.85-(2+8)*1
  !          zeta_p=zeta_s
  !          zeta_d=real(nu_atom,DP)-9*0.35d0-(2+8+8)*1
  !          zeta_f=0.d0
  !       endif
  !       if(nu_atom==46)then !> 3d10 4s1
  !          zeta_s=real(nu_atom,DP)-(1-1)*0.35d0-18*0.85-(2+8)*1
  !          zeta_p=zeta_s
  !          zeta_d=real(nu_atom,DP)-9*0.35d0-(2+8+8)*1
  !          zeta_f=0.d0
  !       endif
  !       if(nu_atom==47)then !> 3d10 4s1
  !          zeta_s=real(nu_atom,DP)-(1-1)*0.35d0-18*0.85-(2+8)*1
  !          zeta_p=zeta_s
  !          zeta_d=real(nu_atom,DP)-9*0.35d0-(2+8+8)*1
  !          zeta_f=0.d0
  !       endif
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/2
  !    case(13:18)
  !       Lmax=1
  !       nquan=(/4,4,0,0/)
  !       zeta_s=real(nu_atom,DP)-(nu_valence-10)*0.35d0-(8+10)*0.85-(2+8)*1
  !       zeta_p=zeta_s
  !       zeta_d=0.d0
  !       zeta_f=0.d0
  !       zeta=(/zeta_s,zeta_p,zeta_d,zeta_f/)/3.7d0
  !    end select
  ! case(55:86)
  !    !> n=6
  ! case(87:104)
  !    !> n=7
  ! ENDSUBROUTINE atom_effcharge_cal

  SUBROUTINE atom_STO(p,l,m,zeta,r,f)!{{{
    USE parameters , ONLY : weight=>Wexict
    !for calculate the slate-obtial
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: p,l,m !quantum number
    REAL(DP),INTENT(IN) :: zeta & !
         &, r(:)  !rvec
    REAL(DP),INTENT(OUT) :: f
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IF ( p==0 .OR. r(4)> 10.d0 ) THEN
    ! IF ( p==0  ) THEN
       f=0.d0
       RETURN
    ENDIF

    select case(p)
    case (1)
       select case (l)
       case (0)
          IF ( m==0 ) THEN
             f=sqrt(zeta**3/pi)*exp(-zeta*r(4))
          ELSE
             stop "m is greater than l or less than -l"
          ENDIF
       case default
          stop "l is greater than n-1 or less than 0"
       end select
    case (2)
       select case (l)
       case (0)
          IF ( m==0 ) THEN
             f=sqrt(zeta**5/(3*pi))*r(4)*exp(-zeta*r(4))
          ELSE
             stop "m is greater than l or less than -l"
          ENDIF
       case (1)
          select case (m)
          case (1)
             f=sqrt(zeta**5/pi)*r(1)*exp(-zeta*r(4))
          case (0)
             f=sqrt(zeta**5/pi)*r(2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(zeta**5/pi)*r(3)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case default
          stop "l is greater than n-1 or less than 0"
       end select
    case (3)
       select case (l)
       case (0)
          IF ( m==0 ) THEN
             f=sqrt((2*zeta**7)/(45*pi))*r(4)**2*exp(-zeta*r(4))
          ELSE
             stop "m is greater than l or less than -l"
          ENDIF
       case (1)
          select case (m)
          case (1)
             f=sqrt(2*zeta**7/(15*pi))*r(4)*r(1)*exp(-zeta*r(4))
          case (0)
             f=sqrt(2*zeta**7/(15*pi))*r(4)*r(2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(2*zeta**7/(15*pi))*r(4)*r(3)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case (2)
          select case (m)
          case (0)
             f=(1.d0/3)*sqrt(zeta**7/(2*pi))*(3*r(3)**2-r(4)**2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(2*zeta**7/(3*pi))*r(1)*r(3)*exp(-zeta*r(4))
          case (1)
             f=-sqrt(2*zeta**7/(3*pi))*r(2)*r(3)*exp(-zeta*r(4))
          case (2)
             f=sqrt(zeta**7/(6*pi))*(r(1)**2-r(2)**2)*exp(-zeta*r(4))
          case (-2)
             f=-sqrt(2*zeta**7/(3*pi))*r(1)*r(2)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case default
          stop "l is greater than n-1 or less than 0"
       end select
    case (4)
       select case (l)
       case (0)
          IF ( m==0 ) THEN
             f=sqrt((zeta**9)/(315*pi))*r(4)**3*exp(-zeta*r(4))
          ELSE
             stop "m is greater than l or less than -l"
          ENDIF
       case (1)
          select case (m)
          case (1)
             f=sqrt(zeta**9/(105*pi))*r(4)**2*r(1)*exp(-zeta*r(4))
          case (0)
             f=sqrt(zeta**9/(105*pi))*r(4)**2*r(2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(zeta**9/(105*pi))*r(4)**2*r(3)*exp(-zeta*r(4))
             !f=sqrt(2*zeta**7/(45*pi))*r(4)*(r(3)+r(2)+r(1))*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case (2)
          select case (m)
          case (0)
             f=(1.d0/6)*sqrt(zeta**9/(7*pi))*r(4)*(3*r(3)**2-r(4)**2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(zeta**9/(21*pi))*r(4)*r(1)*r(3)*exp(-zeta*r(4))
          case (1)
             f=-sqrt(zeta**9/(21*pi))*r(4)*r(2)*r(3)*exp(-zeta*r(4))
          case (2)
             f=sqrt(zeta**9/(84*pi))*r(4)*(r(1)**2-r(2)**2)*exp(-zeta*r(4))
          case (-2)
             f=-sqrt(zeta**9/(21*pi))*r(4)*r(1)*r(2)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case (3)
          select case (m)
          case (3)
             f=-sqrt(zeta**9/(72*pi))*r(2)*(-3*r(1)**2+r(2)**2)*exp(-zeta*r(4))
          case (2)
             f=sqrt(zeta**9/(12*pi))*r(3)*(-r(1)**2+r(2)**2)*exp(-zeta*r(4))
          case (1)
             f=sqrt(zeta**9/(120*pi))*r(2)*(r(4)**2-5*r(3)**2)*exp(-zeta*r(4))
          case (0)
             f=sqrt(zeta**9/(180*pi))*r(3)*(-3*r(4)**2+5*r(3)**2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(zeta**9/(120*pi))*r(1)*(-r(4)**2+5*r(3)**2)*exp(-zeta*r(4))
          case (-2)
             f=-sqrt(zeta**9/(3*pi))*r(1)*r(2)*r(3)*exp(-zeta*r(4))
          case (-3)
             f=sqrt(zeta**9/(72*pi))*r(1)*(-r(4)**2+4*r(2)**2+r(3)**2)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case default
          stop "l is greater than n-1 or less than 0"
       end select
    case (5)
       select case (l)
       case (0)
          IF ( m==0 ) THEN
             f=sqrt(2.d0*(zeta**11)/(14175*pi))*r(4)**4*exp(-zeta*r(4))
          ELSE
             stop "m is greater than l or less than -l"
          ENDIF
       case (1)
          select case (m)
          case (1)
             f=sqrt(2.d0*zeta**11/(4725*pi))*r(4)**3*r(1)*exp(-zeta*r(4))
          case (0)
             f=sqrt(2.d0*zeta**11/(4725*pi))*r(4)**3*r(2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(2.d0*zeta**11/(4725*pi))*r(4)**3*r(3)*exp(-zeta*r(4))
             !f=sqrt(2*zeta**7/(45*pi))*r(4)*(r(3)+r(2)+r(1))*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case (2)
          select case (m)
          case (0)
             f=sqrt(zeta**11/(2835*pi))*r(4)**2*(3*r(3)**2-r(4)**2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(2.d0*zeta**11/(945*pi))*r(4)**2*r(1)*r(3)*exp(-zeta*r(4))
          case (1)
             f=-sqrt(2.d0*zeta**11/(945*pi))*r(4)**2*r(2)*r(3)*exp(-zeta*r(4))
          case (2)
             f=sqrt(zeta**11/(945*pi))*r(4)**2*(r(2)**2-r(1)**2)*exp(-zeta*r(4))
          case (-2)
             f=-sqrt(2.d0*zeta**11/(945*pi))*r(4)**2*r(1)*r(2)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case (3)
          select case (m)
          case (3)
             f=-sqrt(zeta**11/(1620*pi))*r(4)*r(2)*(-3*r(1)**2+r(2)**2)*exp(-zeta*r(4))
          case (2)
             f=sqrt(zeta**11/(270*pi))*r(4)*r(3)*(-r(1)**2+r(2)**2)*exp(-zeta*r(4))
          case (1)
             f=sqrt(zeta**11/(2700*pi))*r(4)*r(2)*(r(4)**2-5*r(3)**2)*exp(-zeta*r(4))
          case (0)
             f=sqrt(zeta**11/(4050*pi))*r(4)*r(3)*(-3*r(4)**2+5*r(3)**2)*exp(-zeta*r(4))
          case (-1)
             f=sqrt(zeta**11/(2700*pi))*r(4)*r(1)*(-r(4)**2+5*r(3)**2)*exp(-zeta*r(4))
          case (-2)
             f=-sqrt(2.d0*zeta**11/(135*pi))*r(4)*r(1)*r(2)*r(3)*exp(-zeta*r(4))
          case (-3)
             f=sqrt(zeta**11/(1620*pi))*r(4)*r(1)*(-r(4)**2+4*r(2)**2+r(3)**2)*exp(-zeta*r(4))
          case default
             stop "m is greater than l or less than -l"
          end select
       case default
          stop "l is greater than n-1 or less than 0"
       end select
    case default
       stop "principl quantum number is greater than 4 or less than 1"
    end select
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE atom_STO!}}}
  ! FUNCTION inc_gamma(alpha,x)
  !   IMPLICIT NONE
  !   REAL(DP) :: inc_gamma
  !   real(DP) ::  alpha
  !   COMPLEX(DCP) :: x
  !   REAL(DP) :: re,one,p,q
  !   REAL(DP) ::  xlim,zero
  !   INTEGER(I4B) :: ibuf,i,ilim

  !   re=0.36787944117144232
  !   one=1
  !   xlim=1.d0
  !   zero=0.d0
  !   ibuf=34
  !   if(dnrm(x).lt.xlim.or.real(x).lt.zero.and.abs(dimag(x)).lt.xlim)then
  !      inc_gamma=re/cdh(alpha,one)
  !      ilim=real(x/re)
  !      do i=0,ibuf-ilim
  !         call term(alpha,x,i,p,q)
  !         inc_gamma=inc_gamma+p*q
  !      enddo
  !   else
  !      inc_gamma=cdexp(-x+alpha*cdlog(x))/cdh(alpha,x)
  !   endif
  !   return
  ! ENDFUNCTION inc_gamma
  ! !>-----------------------The Dividing Line-----------------------------
  ! SUBROUTINE term(alpha,x,i,p,q)
  !   IMPLICIT NONE
  !   REAL(DP) :: alpha,p,q,ci,alphai
  !   COMPLEX(DCP) :: x
  !   REAL(DP) :: zero,one,two,cdlx
  !   REAL(DP) :: tol,xlim,dnrm
  !   INTEGER(I4B) :: i

  !   zero=0.d0
  !   one=1.d0
  !   two=2.d0
  !   tol=3.d-7
  !   xlim=39.d0

  !   if(i.eq.0)q=one
  !   ci=i
  !   alphai=alpha+ci
  !   if(x.eq.zero)then
  !      p=one/alphai
  !      if(i.ne.0)q=-q/ci
  !      return
  !   endif
  !   cdlx=cdlog(x)
  !   ! --- If(1-x**alphai)=-x**alphai on the computer
  !   ! --- then change the inductive scheme to avoid overflow
  !   if(real(alphai*cdlx).gt.xlim.and.i.ne.0)then
  !      p=p*(alphai-one)/alphai
  !      q=-q*x/ci
  !      return
  !   endif
  !   if(dnrm(alphai).gt.tol)then
  !      p=(one-x**alphai)/alphai
  !   else
  !      p=-cdlx*(one+cdlx*alphai/two)
  !   endif
  !   if(i.ne.0)q=-q/ci
  !   return
  ! END SUBROUTINE term
  ! !>-----------------------The Dividing Line-----------------------------
  ! FUNCTION cdh(alpha,x)
  !   ! -- Written By Eric Kostian & Dmitry Gokhman
  !   ! -- March 1986
  !   IMPLICIT NONE
  !   REAL(DP) :: cdh
  !   REAL(DP) :: alpha,cdhs
  !   COMPLEX(DCP) :: x
  !   REAL(DP) :: one,term,sum,cn,alpha1
  !   REAL(DP) :: buf
  !   INTEGER(I4B) :: n,i

  !   one=1.d0
  !   buf=0.d0

  !   ! -- If Re(alpha-x) is too big, shift alpha.
  !   n=real(alpha-x)-buf
  !   if(n.gt.0)then
  !      cn=n+1
  !      alpha1=alpha-cn
  !      term=one/x
  !      sum=term
  !      do i=1,n
  !         cn=n-i+1
  !         term=term*(alpha1+cn)/x
  !         sum=term+sum
  !      enddo
  !      sum=sum+term*alpha1/cdhs(alpha1,x)
  !      cdh=one/sum
  !   else
  !      cdh=cdhs(alpha,x)
  !   endif
  !   return
  ! END FUNCTION cdh
  ! !>-----------------------The Dividing Line-----------------------------
  ! FUNCTION cdhs(alpha,x)
  !   IMPLICIT NONE
  !   ! -- Written By Eric Kostian & Dmitry Gokhman
  !   ! -- March 1986
  !   REAL(DP) :: cdhs
  !   REAL(DP) :: zero,half,one,alpha,x
  !   REAL(DP) :: p0,q0,p1,q1,r0,r1,ci,factor
  !   REAL(DP) :: tol1,tol2,error,dnrm
  !   INTEGER(I4B) :: ilim,i

  !   zero=0.d0
  !   half=0.5d0
  !   one=1.d0
  !   tol1=1.d10
  !   tol2=1.d-10
  !   error=5.d-18
  !   ilim=100000
  !   q0=one
  !   q1=one
  !   p0=x
  !   p1=x+one-alpha
  !   do i=1,ilim
  !      ci=i
  !      if(p0.ne.zero.and.q0.ne.zero.and.q1.ne.zero)then
  !         r0=p0/q0
  !         r1=p1/q1
  !         if(dnrm(r0-r1).le.dnrm(r1)*error)then
  !            cdhs=r1
  !            return
  !         endif
  !         ! -- Occasionally renormalize the sequences to avoid over(under)flow.
  !         if(dnrm(p0).gt.tol1.or.dnrm(p0).lt.tol2.or.&
  !              & dnrm(q0).gt.tol1.or.dnrm(q0).lt.tol2)then
  !            factor=p0*q0
  !            p0=p0/factor
  !            q0=q0/factor
  !            p1=p1/factor
  !            q1=q1/factor
  !         endif
  !      endif
  !      p0=x*p1+ci*p0
  !      q0=x*q1+ci*q0
  !      p1=p0+(ci+one-alpha)*p1
  !      q1=q0+(ci+one-alpha)*q1
  !   enddo
  !   ! -- If the peripheral routines are written correctly
  !   ! -- the following four statements should never be executed.
  !   print*,'cdhs: *** Warning: i >',ilim
  !   print*,'cdhs: *** r0,r1=',r0,r1
  !   cdhs=half*(r0+r1)
  !   return
  ! END FUNCTION cdhs
  ! !>-----------------------The Dividing Line-----------------------------
  ! FUNCTION dnrm(z)
  !   IMPLICIT NONE
  !   REAL(DP) :: dnrm
  !   COMPLEX(DCP) :: z
  !   dnrm=abs(real(z))+abs(aimag(z))
  !   return
  ! END FUNCTION dnrm
  !>-----------------------The Dividing Line-----------------------------
  FUNCTION gammp(a,x)
    IMPLICIT NONE
    REAL(DP) :: a,gammp,x
    !> USES gcf,gser
    !> Returns the incomplete gamma function P(a,x)
    REAL(DP) :: gammcf,gamser,gln

    if(x.lt.0..or.a.le.0.)stop 'bad arguments in gammp'
    if(x.lt.a+1.)then
       !> Use  the  series  representation
       call gser(gamser,a,x,gln)
       gammp=gamser
    else
       !> Use  the  continued  fraction  representation
       !> and  take  its  complement.
       call gcf(gammcf,a,x,gln)
       gammp=1.d0-gammcf
    endif
    return
  ENDFUNCTION gammp
  !>-----------------------The Dividing Line-----------------------------
  SUBROUTINE gser(gamser,a,x,gln)
    IMPLICIT NONE
    REAL(DP) ::     a,gamser,gln,x
    INTEGER(I4B),parameter :: ITMAX=100
    REAL(DP),parameter :: EPS=3.d-7
    !> USES gammln
    !> Returns the incomplete gamma function P(a,x) evaluated by its series
    !> representation as gamser. Also returns ln{\gamma(a)} as gln.
    INTEGER(I4B) :: n
    REAL(DP) :: ap,del,sum

    gln=gammln(a)
    if(x.le.0.)then
       if(x.lt.0.)stop 'x < 0 in gser'
       gamser=0.d0
       return
    endif
    ap=a
    sum=1.d0/a
    del=sum
    do n=1,ITMAX
       ap=ap+1
       del=del*x/ap
       sum=sum+del
       if(abs(del).lt.abs(sum)*EPS)goto 1
    enddo
    stop 'a too large, ITMAX too small in gser'
1   gamser=sum*exp(-x+a*log(x)-gln)
    return
  END SUBROUTINE gser
  !>-----------------------The Dividing Line-----------------------------
  SUBROUTINE gcf(gammcf,a,x,gln)
    IMPLICIT NONE
    REAL(DP) :: a,gammcf,gln,x
    INTEGER(I4B),parameter :: ITMAX=100
    REAL(DP),parameter :: EPS=3.d-7,FPMIN=1.d-30
    !> USES gammln
    !> Returns the incomplete gamma function Q(a,x) evaluates by its continued
    !> fraction representation as gammcf. Also returns ln{\gamma(a)} as gln.
    !> Parameters: ITMAX is the maxium allowed number of iterations; EPS is the
    !> relative accuracy; FPMIN is a number near the smallest representable
    !> floating-point number
    INTEGER(I4B) :: i
    REAL(DP) :: an,b,c,d,del,h

    gln=gammln(a)

    !> Set up for evaluating continued  fraction by modified
    !> Lentz’s method (§5.2) with b0=0
    b=x+1.d0-a
    c=1.d0/FPMIN
    d=1.d0/b
    h=d

    !> Iterate  to  convergence.
    do i=1,ITMAX
       an=-i*(i-a)
       b=b+2.d0
       d=an*d+b
       if(abs(d).lt.FPMIN)d=FPMIN
       c=b+an/c
       if(abs(c).lt.FPMIN)c=FPMIN
       d=1.d0/d
       del=d*c
       h=h*del
       if(abs(del-1.d0).lt.EPS)goto 1
    enddo
    stop 'a too large, ITMAX too small in gcf'

    !> Put factors in front.
1   gammcf=exp(-x+a*log(x)-gln)*h
    return
  END SUBROUTINE gcf
  !>-----------------------The Dividing Line-----------------------------
  FUNCTION gammln(xx)
    REAL(DP) :: gammln,xx
    !> Return the value ln{\gamma(xx)} for xx > 0.
    INTEGER(I4B) :: j
    REAL(DP) :: ser,tmp,x,y
    !> Internal arithmetic will be done in double precision, a nicety that you
    !> can omit if five-figure accuracy is good enough.
    REAL(DP),save :: cof(6)=(/76.18009172947146d0,-86.50532032941677d0,&
         & 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
         & -.5395239384953d-5/),stp=2.5066282746310005d0

    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    return
  END FUNCTION gammln
  !>-----------------------The Dividing Line-----------------------------
  FUNCTION integral(l,dx,x)
    IMPLICIT NONE
    REAL(DP) :: integral
    INTEGER(I4B) :: l
    REAL(DP)  :: dx,x
    !> integral myfun in circulation
    INTEGER(I4B) :: i,j
    REAL(DP) :: t

    t=0.d0
    integral=0.d0
    do i=1,int(x/dx,I4B),1
       t=t+dx
       integral=integral+myfun(l,t)*dx
    enddo
  ENDFUNCTION integral
  !>-----------------------The Dividing Line-----------------------------
  FUNCTION myfun(l,t)
    IMPLICIT NONE
    INTEGER(I4B) :: l
    REAL(DP) :: t
    REAL(DP) :: myfun

    myfun=t**(2*l)*exp(-t**2)
  ENDFUNCTION myfun
  !>-----------------------The Dividing Line-----------------------------
  FUNCTION plgndr(l,m,x)
    IMPLICIT NONE
    INTEGER(I4B) :: l,m
    REAL(DP)     :: plgndr,x
    !> Computes the associated Legendre polynomial P l m (x). Here m and l
    !> are integers satisfying 0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1.
    INTEGER(I4B) :: i,ll
    REAL(DP)     :: fact,pll,pmm,pmmp1,somx2
    if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0)stop 'bad arguments in plgndr'
    pmm=1.d0 !> Compute P_m^m
    if(m.gt.0) then
       somx2=sqrt((1.d0-x)*(1.d0+x))
       fact=1.d0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
       enddo
    endif
    if(l.eq.m)then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.m+1)then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          enddo
          plgndr=pll
       endif
    endif
    return
  ENDFUNCTION plgndr
  !>-----------------------The Dividing Line--------------------------
  SUBROUTINE Lagrange_interpolation_coe(npoint,scatter_x,coe)
    implicit none
    !> in/out
    INTEGER(I4B) :: npoint
    REAL(DP) :: scatter_x(npoint)
    REAL(DP) :: coe(npoint)
    !> local
    INTEGER(I4B) :: i,j
    REAL(DP) :: x_store(npoint-1),product_x

    !>
    do i=1,npoint
       product_x=1
       do j=1,npoint
          if(i==j)cycle
          product_x=product_x*(scatter_x(i)-scatter_x(j))
       enddo
       coe(i)=1.d0/product_x
    enddo
  ENDSUBROUTINE Lagrange_interpolation_coe
  function lagrange_interpolation_x(npoint,x_sample,x_in)
    implicit none
    INTEGER(I4B),intent(in) :: npoint
    REAL(DP),intent(in)     :: x_sample(npoint),x_in
    REAL(DP)         :: lagrange_interpolation_x(npoint)
    !> local
    INTEGER(I4B) :: i,j
    REAL(DP)     :: product_x
    do i=1,npoint
       product_x=1.d0
       do j=1,npoint
          if(i==j)cycle
          product_x=product_x*(x_in-x_sample(j))
       enddo
       lagrange_interpolation_x(i)=product_x
    enddo
  endfunction lagrange_interpolation_x
  SUBROUTINE interpolation_test()
    implicit none
    !> local
    REAL(DP) :: x(9)
    REAL(DP) :: y(9)
    REAL(DP) :: coe(9)
    REAL(DP) :: decompose_x(9)
    REAL(DP) :: x_p,y_p
    INTEGER(I4B) :: i

    !>
    x=(/2,4,6,8,10,12,14,16,18/)
    y=(/2,4,6,8,12,15,17,21,32/)
    CALL lagrange_interpolation_coe(9,x,coe)
    ! print*,'coe',coe
    x_p=1
    open(1173,file='interpolate')
    do i=1,100
       x_p=x_p+0.2
       decompose_x=lagrange_interpolation_x(9,x,x_p)
       y_p=sum(decompose_x*coe*y)
       ! print*,'x_d',decompose_x
       write(1173,*)x_p,y_p
    enddo
    close(1173)
    stop 'interpolate'
  ENDSUBROUTINE interpolation_test

  SUBROUTINE direct_productlm(nll,nml,index_ll,index_ml,mat_in,mat_out)
    !>
    !>
    !>
    IMPLICIT NONE
    INTEGER(I4B),intent(in) :: nll,nml
    INTEGER(I4B),intent(in) :: index_ll(nll),index_ml(nml)
    REAL(DP),intent(in)     :: mat_in(nll,nll)
    REAL(DP),intent(out),target    :: mat_out(nml,nml)
    !> local
    INTEGER(I4B) :: il,i,j
    INTEGER(I4B) :: ia,ib,ima,imb !> boundary for direct product arrays
    INTEGER(I4B) :: imap   !> map from l to l in lm
    INTEGER(I4b) :: big_li,big_lj,two_l
    INTEGER(I4B) :: itest
    LOGICAL :: flag=.false.

    !> calculate the boundary
    !> small matrix
    ia= 0
    ib= 0
    !> big matrix
    ima= 0
    imb= 0
    imap= 0
    ! open(1231,file='test')
    ! write(1231,*)'nll',nll,'index_ll',index_ll
    do il= 1,nll,1
       imap=imap+index_ll(il)*2+1
       flag=.false.
       if(il==nll)then
          flag=.true.
       else
          if( index_ll(il)/=index_ll(il+1) )flag= .true.
       endif
       if(flag)then
          ia= ib+1
          ib= il
          ima= imb+1
          imb= imap
          !> set array value in turns
          do i= ia,ib,1
             do j= ia,ib,1
                two_l= index_ll(il)*2
                big_li= (i-ia)*(two_l+1)+ima
                big_lj= (j-ia)*(two_l+1)+ima
                mat_out(big_li:big_li+two_l,big_lj:big_lj+two_l)= &
                     & unity_matrix(mat_in(i,j),index_ll(il)*2+1)
             enddo
          enddo
          !> test
          ! write(1231,*)'shape small_l',ia,ib
          ! do itest=ia,ib,1
          !    write(1231,'(1000(F6.2,2X))')mat_in(ia:ib,itest)
          ! enddo
          ! write(1231,*)'shape mat_out',ima,imb
          ! do itest=ima,imb,1
          !    write(1231,'(1000(F6.2,2X))')mat_out(ima:imb,itest)
          ! enddo
       endif
    enddo
    ! close(1231)

    CONTAINS
      FUNCTION unity_matrix(val,nsize)
        IMPLICIT NONE
        REAL(DP)     :: val
        INTEGER(I4B) :: nsize
        REAL(DP)     :: unity_matrix(nsize,nsize)
        !> local
        INTEGER(I4B) :: i

        unity_matrix= 0.d0
        do i= 1,nsize,1
           unity_matrix(i,i)= val
        enddo
        ! write(1231,*)'unity_matrix'
        ! do i=1,nsize,1
        !    write(1231,'(1000(F6.2,2X))')unity_matrix(:,i)
        ! enddo
      ENDFUNCTION unity_matrix
  ENDSUBROUTINE direct_productlm
  !--------------------------------------------------------------
  FUNCTION sphbess(l,x)
    ! Arguments
    integer, intent(in) :: l
    real(DP), intent(in) :: x
    REAL(DP) :: sphbess
    ! Local variables
    integer :: j
    real(kind=DP), parameter :: third = 1.0_DP / 3.0_DP
    real(kind=DP), parameter :: ftnth = 1.0_DP / 14.0_DP
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
          sb1 = (sb0 - cos(x)) / x
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
          ux = 1.0_DP / x
          do j=1,l-1
             byp = real(2*J+1,DP)*ux*by - bym
             bym = by
             by = byp
          end do
          sphbess = by
       end if
    end if
  end function sphbess
  !--------------------------------------------------------
  subroutine fourier_1d(nr,rr,rab,vr,ll,nql,yp,vql,vt)

    !use constants
    implicit none
    !     Input/Output variables:
    !     number of grid points
    integer, intent(in) :: nr
    !     coordinates of grid points and differential, dr
    real(dp), intent(in) :: rr(nr), rab(nr)
    !     input function in real space
    real(dp), intent(in) :: vr(nr)
    !     order of Bessel function
    integer, intent(in) :: ll
    !     number of reciprocal space points
    integer, intent(in) :: nql
    !     coordinates of reciprocal space points
    real(dp), intent(in) :: yp(nql)
    !     output function
    real(dp), intent(out) :: vql(nql)
    !     compensating charge
    real(dp), intent(in) :: vt
    !     Work variables:
    !     counters
    integer :: ii, j
    real(dp), dimension(nr+4) :: &
         &     vtmp,     &           ! input function minus compensating function
         &     y                    ! integrand in 1-D Fourier transform
    !
    !     External function:
    real(dp), external :: besselj
    y = 0._DP
    vtmp = 0._DP
    vql = 0._DP
    do ii = 2, nr
       vtmp(ii) = rr(ii)*(rr(ii)*vr(ii)+vt)
    enddo
    !open(1111,FILE='vtmp.dat')
    !   write(1111,*) vtmp(:)
    !close(1111)
    !g=0
    CALL integ_new(rab,vtmp,vql(1))
    vql(1)=vql(1)*4*pi
    !g>0
    do j = 2, nql
       do ii = 2, nr
          y(ii) = vtmp(ii)*sphbess(ll,yp(j)*rr(ii))
       enddo
       CALL integ_new(rab,y,vql(j))
       !do ii = 1, nr, 4
       !   vql(j) = vql(j) + 7.0*y(ii) + 32.0*y(ii+1) +  &
       !   &    12.0*y(ii+2) + 32.0*y(ii+3) + 7.0*y(ii+4)
       !enddo
       !vql(j) = (2._DP*yp(j)*yp(j)*vql(j)/45.0 - vt)*2._DP/pi
       vql(j) = (yp(j)*yp(j)*vql(j) - vt)*2._DP/pi
       vql(j) = vql(j)*2._DP*pi**2/yp(j)**2
    enddo
    !vql(1)=-2._DP*vt/pi
    !open(1112,FILE='vql.dat')
    !   write(1112,*) vql(:)
    !close(1112)
  end subroutine fourier_1d
  !--------------------------------------------------------
  SUBROUTINE invfourier_1d(g,fg,ll,r,fr)
    IMPLICIT NONE
    REAL(DP),INTENT(IN) :: g(:) , fg(:) ,r(:)
    REAL(DP),INTENT(OUT) :: fr(:)
    INTEGER(I4B),INTENT(IN) :: ll
    !LOCAL
    REAL(DP) :: dg(SIZE(g))
    INTEGER(I4B) :: ng,nr,i,j
    REAL(DP),DIMENSION(SIZE(g)+4) :: yp,ypt
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nr=SIZE(r,1)
    ng=SIZE(g,1)
    fr(:)=0._DP
    yp(:)=0._DP
    ypt(:)=0._DP
    dg(1)=g(2)-g(1)
    dg(ng)=g(ng)-g(ng-1)
    DO i=1,ng
       dg(i)=0.5_DP*(g(i+1)-g(i-1))
    ENDDO

    DO i=1,ng
       yp(i)=fg(i)*g(i)**2
    ENDDO
    CALL integ_new(dg,yp,fr(1))
    DO j=2,nr
       DO i=2,ng
          ypt(i)=yp(i)*sphbess(ll,r(j)*g(i))
       ENDDO
       CALL integ_new(dg,ypt,fr(j))
    ENDDO
    fr(:)=fr(:)/(2*pi**2)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE invfourier_1d
  !--------------------------------------------------------
  SUBROUTINE integ_new(rab,y,f)
    IMPLICIT NONE
    REAL(DP),INTENT(IN) :: rab(:),y(:)
    REAL(DP),INTENT(OUT) :: f
    !LOCAL
    INTEGER(I4B) :: nr,i
    REAL(DP) :: yp(SIZE(y))
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nr=SIZE(rab)
    f=0._DP
    yp=0._DP
    DO i=2,nr
       yp(i)=y(i)*rab(i)
    ENDDO
    !
    DO i=1,nr,4
       f=f+ 7._DP*yp(i) + 32._DP*yp(i+1) + 12._DP*yp(i+2) &
            & + 32._DP*yp(i+3) + 7._DP*yp(i+4)
    ENDDO
    f=f*2._DP/45._DP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE integ_new
  !--------------------------------------------------------
  FUNCTION interp(np,f,r,rnorm,Z)
    USE MathSplines, only: polynom
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN) :: np
    REAL(DP),INTENT(IN) :: f(:),r(:),rnorm
    REAL(DP),OPTIONAL :: Z
    REAL(DP) :: interp
    !LOCAL
    REAL(DP) :: rmax,rmin
    INTEGER(I4B) :: I
    REAL(DP) :: c(3)
    LOGICAL  :: lz
    LOGICAL  :: lfind=.FALSE.
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    lz=PRESENT(Z)
    c(:)=0._DP
    rmax=r(np)
    rmin=r(1)
    IF(rnorm>=rmax)THEN
       IF(lz)THEN
          interp=-Z/rnorm !for local part potential
          RETURN
       ELSE
          interp=0._DP
          RETURN
       ENDIF
    ELSEIF(rnorm<rmin)THEN
       interp=f(1)
       RETURN
    ELSE
       lfind=.FALSE.
       DO I=1,np-1
          IF((rnorm>=r(I)).AND.(rnorm<r(I+1)))THEN
             lfind=.TRUE.
             EXIT
          ENDIF
       ENDDO
       !Check
       IF(.NOT.lfind)THEN
          WRITE(*,*) 'STOP: Check the interp()',rnorm,rmax
          STOP
       ENDIF
       !interpole it 
       IF(I>=np-2)THEN
          interp=polynom(0,3,r(np-2:np),f(np-2:np),c,rnorm)  
          RETURN
       ELSEIF(I<=2)THEN
          interp=polynom(0,3,r(1:3),f(1:3),c,rnorm)  
          RETURN
       ELSE
          interp=polynom(0,3,r(I-1:I+1),f(I-1:I+1),c,rnorm)  
          RETURN
       ENDIF
    ENDIF
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION interp
  !> -----------------------------------------------------
  !> parallel code from wangsheng, add by xlt
  Subroutine Three2one_d_cplx(amat,bmat)  !{{{
    implicit none
    complex(dp),dimension(:,:,:),intent(in)  :: amat
    complex(dp),dimension(:),intent(inout)  :: bmat
    integer(i4b)  :: Ii,iz,iy,ix
    !-------------------------------------------------------------
    !allocate(bmat(grid%n))
    Ii = 0
    do iz = 1, size(amat,3)
       do iy = 1, size(amat,2)
          do ix = 1, size(amat,1)
             Ii = Ii + 1
             bmat(Ii) = amat(ix,iy,iz) 
          end do
       end do
    end do

  End Subroutine three2one_d_cplx !}}}
  !
  Subroutine Three2one_d_real(amat,bmat)  !{{{
    implicit none
    real(dp),dimension(:,:,:),intent(in)  :: amat
    real(dp),dimension(:),intent(inout)  :: bmat
    integer(i4b)  :: Ii,iz,iy,ix
    !-------------------------------------------------------------
    !allocate(bmat(grid%n))
    Ii = 0
    do iz = 1,size(amat,3) 
       do iy = 1,size(amat,2) 
          do ix = 1,size(amat,1) 
             Ii = Ii + 1
             bmat(Ii) = amat(ix,iy,iz) 
          end do
       end do
    end do

  End Subroutine three2one_d_real  !}}}
  !
  Subroutine One2three_d_cplx(amat,bmat)  !{{{
    implicit none
    complex(dp),dimension(:,:,:),intent(inout) :: bmat
    complex(dp),dimension(:),intent(in)      :: amat
    integer(i4b)  :: Ii ,iz,iy,ix
    !-------------------------------------------------------------
    !allocate(bmat(grid%n1,grid%n2,grid%n3))
    Ii = 0
    do iz = 1,size(bmat,3) 
       do iy = 1,size(bmat,2)
          do ix = 1,size(bmat,1)
             Ii = Ii + 1
             bmat(ix,iy,iz) = amat(Ii)
          end do
       end do
    end do

  End Subroutine One2three_d_cplx  !}}}
  !
  Subroutine One2three_d_real(amat,bmat)  !{{{
    implicit none
    real(dp),dimension(:,:,:),intent(inout) :: bmat
    real(dp),dimension(:),intent(in)      :: amat
    integer(i4b)  :: Ii ,iz,iy,ix
    !-------------------------------------------------------------
    !allocate(bmat(grid%n1,grid%n2,grid%n3))
    Ii = 0
    do iz = 1,size(bmat,3)  
       do iy = 1,size(bmat,2)
          do ix = 1,size(bmat,1)
             Ii = Ii + 1
             bmat(ix,iy,iz) = amat(Ii)
          end do
       end do
    end do

  End Subroutine One2three_d_real  !}}}
  Function GetFileUnit()!{{{
    implicit none 
    INTEGER(I4B),save :: funit=5090
    INTEGER(I4B)      :: GetFileUnit
    LOGICAL           :: lopen
    funit = funit + 1
    do   
       inquire( UNIT=funit, opened=lopen )
       if ( lopen ) then 
          funit = funit + 1
       else 
          exit 
       endif
    enddo
    GetFileUnit=funit
  end function GetFileUnit!}}}

  function kahan_sum(n,array)
    implicit none
    REAL(DP) :: kahan_sum
    INTEGER(I4B) :: n
    REAL(DP) :: array(n)
    !> local
    INTEGER(I4B) :: i
    REAL(DP) :: eps,y,t

    eps=0.d0
    kahan_sum=0.d0
    do i=1,n
       y=array(i)-eps
       t=kahan_sum+y
       eps=(t-kahan_sum)-y
       kahan_sum=t
    enddo

  end function kahan_sum

ENDMODULE Math
