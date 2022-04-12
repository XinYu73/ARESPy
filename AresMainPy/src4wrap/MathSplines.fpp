# 1 "MathSplines.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "MathSplines.f90"
MODULE MathSplines
!
! 		   GNU LESSER GENERAL PUBLIC LICENSE
!                        Version 3, 29 June 2007
!
!  Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
!  Everyone is permitted to copy and distribute verbatim copies
!  of this license document, but changing it is not allowed.
!
!
!   This version of the GNU Lesser General Public License incorporates
! the terms and conditions of version 3 of the GNU General Public
! License, supplemented by the additional permissions listed below.
!
!   0. Additional Definitions.
!
!   As used herein, "this License" refers to version 3 of the GNU Lesser
! General Public License, and the "GNU GPL" refers to version 3 of the GNU
! General Public License.
!
!   "The Library" refers to a covered work governed by this License,
! other than an Application or a Combined Work as defined below.
!
!   An "Application" is any work that makes use of an interface provided
! by the Library, but which is not otherwise based on the Library.
! Defining a subclass of a class defined by the Library is deemed a mode
! of using an interface provided by the Library.
!
!   A "Combined Work" is a work produced by combining or linking an
! Application with the Library.  The particular version of the Library
! with which the Combined Work was made is also called the "Linked
! Version".
!
!   The "Minimal Corresponding Source" for a Combined Work means the
! Corresponding Source for the Combined Work, excluding any source code
! for portions of the Combined Work that, considered in isolation, are
! based on the Application, and not on the Linked Version.
!
!   The "Corresponding Application Code" for a Combined Work means the
! object code and/or source code for the Application, including any data
! and utility programs needed for reproducing the Combined Work from the
! Application, but excluding the System Libraries of the Combined Work.
!
!   1. Exception to Section 3 of the GNU GPL.
!
!   You may convey a covered work under sections 3 and 4 of this License
! without being bound by section 3 of the GNU GPL.
!
!   2. Conveying Modified Versions.
!
!   If you modify a copy of the Library, and, in your modifications, a
! facility refers to a function or data to be supplied by an Application
! that uses the facility (other than as an argument passed when the
! facility is invoked), then you may convey a copy of the modified
! version:
!
!    a) under this License, provided that you make a good faith effort to
!    ensure that, in the event an Application does not supply the
!    function or data, the facility still operates, and performs
!    whatever part of its purpose remains meaningful, or
!
!    b) under the GNU GPL, with none of the additional permissions of
!    this License applicable to that copy.
!
!   3. Object Code Incorporating Material from Library Header Files.
!
!   The object code form of an Application may incorporate material from
! a header file that is part of the Library.  You may convey such object
! code under terms of your choice, provided that, if the incorporated
! material is not limited to numerical parameters, data structure
! layouts and accessors, or small macros, inline functions and templates
! (ten or fewer lines in length), you do both of the following:
!
!    a) Give prominent notice with each copy of the object code that the
!    Library is used in it and that the Library and its use are
!    covered by this License.
!
!    b) Accompany the object code with a copy of the GNU GPL and this license
!    document.
!
!   4. Combined Works.
!
!   You may convey a Combined Work under terms of your choice that,
! taken together, effectively do not restrict modification of the
! portions of the Library contained in the Combined Work and reverse
! engineering for debugging such modifications, if you also do each of
! the following:
!
!    a) Give prominent notice with each copy of the Combined Work that
!    the Library is used in it and that the Library and its use are
!    covered by this License.
!
!    b) Accompany the Combined Work with a copy of the GNU GPL and this license
!    document.
!
!    c) For a Combined Work that displays copyright notices during
!    execution, include the copyright notice for the Library among
!    these notices, as well as a reference directing the user to the
!    copies of the GNU GPL and this license document.
!
!    d) Do one of the following:
!
!        0) Convey the Minimal Corresponding Source under the terms of this
!        License, and the Corresponding Application Code in a form
!        suitable for, and under terms that permit, the user to
!        recombine or relink the Application with a modified version of
!        the Linked Version to produce a modified Combined Work, in the
!        manner specified by section 6 of the GNU GPL for conveying
!        Corresponding Source.
!
!        1) Use a suitable shared library mechanism for linking with the
!        Library.  A suitable mechanism is one that (a) uses at run time
!        a copy of the Library already present on the user's computer
!        system, and (b) will operate properly with a modified version
!        of the Library that is interface-compatible with the Linked
!        Version.
!
!    e) Provide Installation Information, but only if you would otherwise
!    be required to provide such information under section 6 of the
!    GNU GPL, and only to the extent that such information is
!    necessary to install and execute a modified version of the
!    Combined Work produced by recombining or relinking the
!    Application with a modified version of the Linked Version. (If
!    you use option 4d0, the Installation Information must accompany
!    the Minimal Corresponding Source and Corresponding Application
!    Code. If you use option 4d1, you must provide the Installation
!    Information in the manner specified by section 6 of the GNU GPL
!    for conveying Corresponding Source.)
!
!   5. Combined Libraries.
!
!   You may place library facilities that are a work based on the
! Library side by side in a single library together with other library
! facilities that are not Applications and are not covered by this
! License, and convey such a combined library under terms of your
! choice, if you do both of the following:
!
!    a) Accompany the combined library with a copy of the same work based
!    on the Library, uncombined with any other library facilities,
!    conveyed under the terms of this License.
!
!    b) Give prominent notice with the combined library that part of it
!    is a work based on the Library, and explaining where to find the
!    accompanying uncombined form of the same work.
!
!   6. Revised Versions of the GNU Lesser General Public License.
!
!   The Free Software Foundation may publish revised and/or new versions
! of the GNU Lesser General Public License from time to time. Such new
! versions will be similar in spirit to the present version, but may
! differ in detail to address new problems or concerns.
!
!   Each version is given a distinguishing version number. If the
! Library as you received it specifies that a certain numbered version
! of the GNU Lesser General Public License "or any later version"
! applies to it, you have the option of following the terms and
! conditions either of that published version or of any later version
! published by the Free Software Foundation. If the Library as you
! received it does not specify a version number of the GNU Lesser
! General Public License, you may choose any version of the GNU Lesser
! General Public License ever published by the Free Software Foundation.
!
!   If the Library as you received it specifies that a proxy can decide
! whether future versions of the GNU Lesser General Public License shall
! apply, that proxy's public statement of acceptance of any version is
! permanent authorization for you to choose that version for the
! Library.
!
! THIS SOFTWARE WAS OBTAINED FROM http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.f90
! IT HAS BEEN SHORTENED IN ORDER TO BE BETTER USABLE, THE ABOVE LICENSING APPLIES SOLELY
! TO THESE FUNCTIONS!


  USE CONSTANTS, ONLY : DP, I4B

  CONTAINS

subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
!*******************************************************************************
!
!! SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
!
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, INTEGER(I4B) N, the number of data points; N must be at least 2.
!
!    Input, real T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real Y(N), the data values to be interpolated.
!
!    Input, INTEGER(I4B) IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real YBCBEG, the left boundary value, if needed.
!
!    Input, INTEGER(I4B) IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real YBCEND, the right boundary value, if needed.
!
!    Output, real YPP(N), the second derivatives of the cubic spline.
!
  implicit none
!
  INTEGER(I4B) :: n
!
  real(kind=dp) ::  diag(n)
  INTEGER(I4B) :: i
  INTEGER(I4B) :: ibcbeg
  INTEGER(I4B) :: ibcend
  real(kind=dp) ::  sub(2:n)
  real(kind=dp) ::  sup(1:n-1)
  real(kind=dp) ::  t(n)
  real(kind=dp) ::  y(n)
  real(kind=dp) ::  ybcbeg
  real(kind=dp) ::  ybcend
  real(kind=dp) ::  ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i6)' ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.0E+00
    diag(1) = 1.0E+00
    sup(1) = -1.0E+00
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    diag(1) = ( t(2) - t(1) ) / 3.0E+00
    sup(1) = ( t(2) - t(1) ) / 6.0E+00
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    diag(1) = 1.0E+00
    sup(1) = 0.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
    stop
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    sub(i) = ( t(i) - t(i-1) ) / 6.0E+00
    diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00
    sup(i) = ( t(i+1) - t(i) ) / 6.0E+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0.0E+00
    sub(n) = -1.0E+00
    diag(n) = 1.0E+00
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    sub(n) = ( t(n) - t(n-1) ) / 6.0E+00
    diag(n) = ( t(n) - t(n-1) ) / 3.0E+00
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    sub(n) = 0.0E+00
    diag(n) = 1.0E+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  The input value is IBCEND = ', ibcend
    stop
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0E+00
    ypp(2) = 0.0E+00
!
!  Solve the linear system.
!
  else

    call s3_fs ( sub, diag, sup, n, ypp, ypp )

  end if

  return
end  SUBROUTINE

subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
!
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Modified:
!
!    20 November 2000
!    30 January 2014 (intents and optimizations)
!
!  Author:
!
!    John Burkardt
!    Johannes M Dieterich
!
!  Parameters:
!
!    Input, INTEGER(I4B) N, the number of data values.
!
!    Input, real T(N), the knot values.
!
!    Input, real Y(N), the data values at the knots.
!
!    Input, real YPP(N), the second derivatives of the spline at the knots.
!
!    Input, real TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none
!
  INTEGER(I4B),intent(in):: n
  real(kind=dp),intent(in)::  t(n)
  real(kind=dp),intent(in)::  tval
  real(kind=dp),intent(in)::  y(n)
  real(kind=dp),intent(in)::  ypp(n)
  real(kind=dp),intent(out)::  yppval
  real(kind=dp),intent(out)::  ypval
  real(kind=dp),intent(out)::  yval
  !
  real(kind=dp) :: dt
  real(kind=dp) :: h
  real(kind=dp) :: tmp1, tmp2, tmp3, tmp4
  INTEGER(I4B) :: left
  INTEGER(I4B) :: right
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  tmp1 = ( y(right) - y(left) ) / h
  tmp2 = - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h
  tmp3 = ( ypp(right) - ypp(left) )
  tmp4 = tmp3/h

  yval = y(left) &
       + dt * ( ( tmp1 ) &
       + tmp2 &
       + dt * ( 0.5E+00 * ypp(left) &
       + dt * ( tmp3 / ( 6.0E+00 * h ) ) ) )

  ypval = ( tmp1 ) &
       + tmp2 &
       + dt * ( ypp(left) &
       + dt * ( 0.5E+00 * tmp4 ) )

  yppval = ypp(left) + dt * tmp4


  return
end SUBROUTINE

  SUBROUTINE rvec_bracket ( n, x, xval, left, right )
  !  this subroutine is used for spline by spline_cubic_set() and
  !  spline_cubic_val() above.
  !*******************************************************************************
  !
  !! RVEC_BRACKET searches a sorted array for successive brackets of a value.
  !
  !
  !  Discussion:
  !
  !    If the values in the vector are thought of as defining intervals
  !    on the real line, then this routine searches for the interval
  !    nearest to or containing the given value.
  !
  !  Modified:
  !
  !    06 April 1999
  !    30 January 2014 (intents and opts)
  !
  !  Author:
  !
  !    John Burkardt
  !    Johannes M Dieterich
  !
  !  Parameters:
  !
  !    Input, INTEGER(I4B) N, length of input array.
  !
  !    Input, real X(N), an array sorted into ascending order.
  !
  !    Input, real XVAL, a value to be bracketed.
  !
  !    Output, INTEGER(I4B) LEFT, RIGHT, the results of the search.
  !    Either:
  !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
  !      XVAL > X(N), when LEFT = N-1, RIGHT = N;
  !    or
  !      X(LEFT) <= XVAL <= X(RIGHT).
  !
    implicit none
    !
    INTEGER(I4B),intent(in)::n
    INTEGER(I4B),intent(out)::left
    INTEGER(I4B),intent(out)::right
    real(kind=dp),intent(in)::x(n)
    real(kind=dp),intent(in)::xval
    !
    INTEGER(I4B) i1, i2, i3!, i
!goto 20
    !
    ! a faster version, has been tested for Al fcc bulk,
    ! JMD: no idea who was the original author here?!
    !      but it certainly wasn't as fast as it could
    !      have been and now is
    !
    i1 = 1
    i3 = n
    i2 = (i3+i1)/2 !floor(real(i3+i1,kind=dp)*0.5_DP)

  do while (i2/=i1)
    if (xval == x(i2)) then
      left = i2
      right = i2+1
      return
    else if (xval < x(i2)) then
      i3 = i2
    else if (xval > x(i2)) then
      i1 = i2
    end if
    i2 = (i3+i1)/2 !floor(real(i3+i1,kind=dp)*0.5_DP)
  end do

  if (xval == x(i3) ) then
    left = i3
    right = i3 + 1
  else if (xval==x(i1)) then
    left = i1
    right = i1 + 1
  else
    left = i1
    right = i3
  endif
  return

!
! old version
!
!20 continue
!    do i = 2, n - 1
!       if ( xval < x(i) ) then
!           left = i - 1
!           right = i
!           return
!       end if
!    end do
!
!    left = n - 1
!    right = n
!
!    return
  END SUBROUTINE rvec_bracket



  subroutine s3_fs ( a1, a2, a3, n, b, x )
  !
  !*******************************************************************************
  !
  !! S3_FS factors and solves a tridiagonal linear system.
  !
  !
  !  Note:
  !
  !    This algorithm requires that each diagonal entry be nonzero.
  !
  !  Modified:
  !
  !    05 December 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
  !    On input, the nonzero diagonals of the linear system.
  !    On output, the data in these vectors has been overwritten
  !    by factorization information.
  !
  !    Input, INTEGER(I4B) N, the order of the linear system.
  !
  !    Input/output, real B(N).
  !    On input, B contains the right hand side of the linear system.
  !    On output, B has been overwritten by factorization information.
  !
  !    Output, real X(N), the solution of the linear system.
  !
    implicit none
    !
      INTEGER(I4B) :: n
      !
        real(kind=dp) :: a1(2:n)
        real(kind=dp) :: a2(1:n)
        real(kind=dp) :: a3(1:n-1)
        real(kind=dp) :: b(n)
        INTEGER(I4B) :: i
        real(kind=dp) :: x(n)
        real(kind=dp) :: xmult
        !
                    !  The diagonal entries can't be zero.
                    !
         do i = 1, n
             if ( a2(i) == 0.0E+00 ) then
                 write ( *, '(a)' ) ' '
                 write ( *, '(a)' ) 'S3_FS - Fatal error!'
                 write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
                 return
             end if
         end do

         do i = 2, n-1
             xmult = a1(i) / a2(i-1)
             a2(i) = a2(i) - xmult * a3(i-1)

             b(i) = b(i) - xmult * b(i-1)

         end do


         xmult = a1(n) / a2(n-1)
         a2(n) = a2(n) - xmult * a3(n-1)

         x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
         do i = n-1, 1, -1
             x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
         end do

         return
END SUBROUTINE
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: polynom
! !INTERFACE:
function polynom(m,np,xa,ya,c,x)
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
real(DP), intent(in) :: xa(np),ya(np)
real(DP), intent(out) :: c(np)
real(DP), intent(in) :: x
! local variables
integer(i4b) :: i,j,k
real(DP) :: x0,x1,x2,x3,y0,y1,y2,y3
real(DP) :: t0,t1,t2,t3,t4,t5,t6
real(DP) :: c1,c2,c3,sum
! fast evaluations for small np
select case(np)
case(1)
  select case(m)
  case(:-1)
    polynom=ya(1)*(x-xa(1))
  case(0)
    polynom=ya(1)
  case default
    polynom=0.d0
  end select
  return
case(2)
  c1=(ya(2)-ya(1))/(xa(2)-xa(1))
  t1=x-xa(1)
  select case(m)
  case(:-1)
    polynom=t1*(ya(1)+0.5d0*c1*t1)
  case(0)
    polynom=c1*t1+ya(1)
  case(1)
    polynom=c1
  case default
    polynom=0.d0
  end select
  return
case(3)
  x0=xa(1)
  x1=xa(2)-x0
  x2=xa(3)-x0
  y0=ya(1)
  y1=ya(2)-y0
  y2=ya(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2
  t2=x2*y1
  c1=x2*t2-x1*t1
  c2=t1-t2
  t1=x-x0
  select case(m)
  case(:-1)
    polynom=t1*(y0+t0*t1*(0.5d0*c1+0.3333333333333333333d0*c2*t1))
  case(0)
    polynom=y0+t0*t1*(c1+c2*t1)
  case(1)
    polynom=t0*(2.d0*c2*t1+c1)
  case(2)
    polynom=t0*2.d0*c2
  case default
    polynom=0.d0
  end select
  return
case(4)
  x0=xa(1)
  x1=xa(2)-x0
  x2=xa(3)-x0
  x3=xa(4)-x0
  y0=ya(1)
  y1=ya(2)-y0
  y2=ya(3)-y0
  y3=ya(4)-y0
  t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
  t1=x1*x2*y3
  t2=x2*x3*y1
  t3=x3*x1*y2
  c3=t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1)
  t6=x3**2
  t5=x2**2
  t4=x1**2
  c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
  c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
  t1=x-x0
  select case(m)
  case(:-1)
    polynom=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.3333333333333333333d0*c2 &
     +0.25d0*c3*t1)))
  case(0)
    polynom=y0+t0*t1*(c1+t1*(c2+c3*t1))
  case(1)
    polynom=t0*(c1+t1*(2.d0*c2+3.d0*c3*t1))
  case(2)
    polynom=t0*(6.d0*c3*t1+2.d0*c2)
  case(3)
    polynom=t0*6.d0*c3
  case default
    polynom=0.d0
  end select
  return
end select
if (np.le.0) then
  WRITE(6,*)
  WRITE(6,'("Error(polynom): np <= 0 : ",I8)') np
  WRITE(6,*)
  stop
end if
if (m.ge.np) then
  polynom=0.d0
  return
end if
! find the polynomial coefficients in divided differences form
c(:)=ya(:)
do i=2,np
  do j=np,i,-1
    c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
  end do
end do
! special case m=0
if (m.eq.0) then
  sum=c(1)
  t1=1.d0
  do i=2,np
    t1=t1*(x-xa(i-1))
    sum=sum+c(i)*t1
  end do
  polynom=sum
  return
end if
x0=xa(1)
! convert to standard form
do j=1,np-1
  do i=1,np-j
    k=np-i
    c(k)=c(k)+(x0-xa(k-j+1))*c(k+1)
  end do
end do
if (m.gt.0) then
! take the m th derivative
  do j=1,m
    do i=m+1,np
      c(i)=c(i)*dble(i-j)
    end do
  end do
  t1=c(np)
  t2=x-x0
  do i=np-1,m+1,-1
    t1=t1*t2+c(i)
  end do
  polynom=t1
else
! find the integral
  t1=c(np)/dble(np)
  t2=x-x0
  do i=np-1,1,-1
    t1=t1*t2+c(i)/dble(i)
  end do
  polynom=t1*t2
end if
return
end function

end module MathSplines
