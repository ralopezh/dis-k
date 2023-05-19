!
! A couple of modified routines from special_functions.f90
!
subroutine ik01ae ( x, bi0, di0, bi1, di1)

!*****************************************************************************80
!
!! IK01A compute Bessel function I0(x), I1(x), Exponentially scalled *Exp(-|x|)
!
!  Discussion:
!
!    This procedure computes modified Bessel functions I0(x), I1(x),
!    and their derivatives Exponentially scaled.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!    14 July 2021 R.A.L
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modified by R. A. Lopez: the result is now exponentially scaled
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BI0, DI0, BI1, DI1, the
!    values of I0(x), I0'(x), I1(x), I1'(x).
!
  implicit none

  real ( kind = 8 ), save, dimension ( 12 ) :: a = (/ &
    0.125D+00, 7.03125D-02, &
    7.32421875D-02, 1.1215209960938D-01, &
    2.2710800170898D-01, 5.7250142097473D-01, &
    1.7277275025845D+00, 6.0740420012735D+00, &
    2.4380529699556D+01, 1.1001714026925D+02, &
    5.5133589612202D+02, 3.0380905109224D+03 /)
  real ( kind = 8 ), save, dimension ( 12 ) :: b = (/ &
    -0.375D+00, -1.171875D-01, &
    -1.025390625D-01, -1.4419555664063D-01, &
    -2.7757644653320D-01, -6.7659258842468D-01, &
    -1.9935317337513D+00, -6.8839142681099D+00, &
    -2.7248827311269D+01, -1.2159789187654D+02, &
    -6.0384407670507D+02, -3.3022722944809D+03 /)
  real ( kind = 8 ) bi0
  real ( kind = 8 ) bi1
  real ( kind = 8 ) ca
  real ( kind = 8 ) di0
  real ( kind = 8 ) di1
  real ( kind = 8 ) el
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) xr

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  x2 = x * x

  if ( x == 0.0D+00 ) then

    bi0 = 1.0D+00
    bi1 = 0.0D+00
    di0 = 0.0D+00
    di1 = 0.5D+00
    return

  else if ( x <= 18.0D+00 ) then

    bi0 = 1.0D+00
    r = 1.0D+00
    do k = 1, 50
      r = 0.25D+00 * r * x2 / ( k * k )
      bi0 = bi0 + r
      if ( abs ( r / bi0 ) < 1.0D-15 ) then
        exit
      end if
    end do
    bi0 = bi0*exp( -ABS(x) )

    bi1 = 1.0D+00
    r = 1.0D+00
    do k = 1, 50
      r = 0.25D+00 * r * x2 / ( k * ( k + 1 ) )
      bi1 = bi1 + r
      if ( abs ( r / bi1 ) < 1.0D-15 ) then
        exit
      end if
    end do

    bi1 = 0.5D+00 * x * bi1 * exp( -ABS(x) )

  else

    if ( x < 35.0D+00 ) then
      k0 = 12
    else if ( x < 50.0D+00 ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = 1.0D+00 / sqrt ( 2.0D+00 * pi * x )
    bi0 = 1.0D+00
    xr = 1.0D+00 / x
    do k = 1, k0
      bi0 = bi0 + a(k) * xr ** k
    end do
    bi0 = ca * bi0
    bi1 = 1.0D+00
    do k = 1, k0
      bi1 = bi1 + b(k) * xr ** k
    end do
    bi1 = ca * bi1

  end if

  di0 = bi1
  di1 = bi0 - bi1 / x

  return
end
subroutine iknae ( n, x, nm, bi, di )

!*****************************************************************************80
!
!! IKNA compute Bessel function In(x) and derivatives.
!  Exponentially scaled *EXP(-|x|)
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!    14 July 2021 RAL
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modified by R. A. Lopez: the result is now exponentially scaled 
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of In(x) and Kn(x).
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) BI(0:N), DI(0:N), 
!    the values of In(x), In'(x).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bi(0:n)
  real ( kind = 8 ) bi0
  real ( kind = 8 ) bi1
  real ( kind = 8 ) di(0:n)
  real ( kind = 8 ) di0
  real ( kind = 8 ) di1
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) h
  real ( kind = 8 ) h0
  real ( kind = 8 ) h1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) nm
  real ( kind = 8 ) s0
  real ( kind = 8 ) x

  nm = n

  if ( x <= 1.0D-100 ) then
    do k = 0, n
      bi(k) = 0.0D+00
      di(k) = 0.0D+00
    end do
    bi(0) = 1.0D+00
    di(1) = 0.5D+00
    return
  end if

  call ik01ae ( x, bi0, di0, bi1, di1 )
  bi(0) = bi0
  bi(1) = bi1
  di(0) = di0
  di(1) = di1

  if ( n <= 1 ) then
    return
  end if

  if ( 40.0D+00 < x .and. n < int ( 0.25D+00 * x ) ) then

    h0 = bi0
    h1 = bi1
    do k = 2, n
      h = -2.0D+00 * ( k - 1.0D+00 ) / x * h1 + h0
      bi(k) = h
      h0 = h1
      h1 = h
    end do

  else

    m = msta1 ( x, 200 )

    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f0 = 0.0D+00
    f1 = 1.0D-100
    do k = m, 0, -1
      f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x + f0
      if ( k <= nm ) then
        bi(k) = f
      end if
      f0 = f1
      f1 = f
    end do
    s0 = bi0 / f
    do k = 0, nm
      bi(k) = s0 * bi(k)
    end do
  end if
  

  do k = 2, nm
    di(k) = bi(k-1) - k / x * bi(k)
  end do

  return
end
function msta1 ( x, mp )

!*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward  
!    recurrence such that the magnitude of    
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) MP, the negative logarithm of the 
!    desired magnitude.
!
!    Output, integer ( kind = 4 ) MSTA1, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) x

  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
    nn = n1 - int ( real ( n1 - n0, kind = 8 ) / ( 1.0D+00 - f0 / f1 ) )               
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta1 = nn

  return
end
function msta2 ( x, n, mp )

!*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!    Jianming Jin supplied a modification to this code on 12 January 2016.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Jn(x).
!
!    Input, integer ( kind = 4 ) N, the order of Jn(x).
!
!    Input, integer ( kind = 4 ) MP, the number of significant digits.
!
!    Output, integer ( kind = 4 ) MSTA2, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) ejn
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) hmp
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) obj
  real ( kind = 8 ) x

  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1D+00 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1D+00 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
    nn = n1 - int ( real ( n1 - n0, kind = 8 ) / ( 1.0D+00 - f0 / f1 ) )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta2 = nn + 10

  return
end
function envj ( n, x )

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Discussion:
!
!    ENVJ estimates -log(Jn(x)) from the estimate
!    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modifications suggested by Vincent Lafage, 11 January 2016.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bessel function.
!
!    Input, real ( kind = 8 ) X, the absolute value of the argument.
!
!    Output, real ( kind = 8 ) ENVJ, the value.
!
  implicit none

  real ( kind = 8 ) envj
  real ( kind = 8 ) logten
  integer ( kind = 4 ) n
  real ( kind = 8 ) n_r8
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x
!
!  Original code
!
  if ( .true. ) then

    envj = 0.5D+00 * log10 ( 6.28D+00 * n ) &
      - n * log10 ( 1.36D+00 * x / n )
!
!  Modification suggested by Vincent Lafage.
!
  else

    n_r8 = real ( n, kind = 8 )
    logten = log ( 10.0D+00 )
    envj = r8_gamma_log ( n_r8 + 1.0D+00 ) / logten - n_r8 * log10 ( x )

  end if

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2.  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0 are similar to those in
!    reference 1, but are unpublished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ) :: d1 = -5.772156649015328605195174D-01
  real ( kind = 8 ) :: d2 = 4.227843350984671393993777D-01
  real ( kind = 8 ) :: d4 = 1.791759469228055000094023D+00
  real ( kind = 8 ), parameter :: frtbig = 2.25D+76
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 2.55D+305
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.79D+308
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq

  y = x

  if ( 0.0D+00 < y .and. y <= xbig ) then

    if ( y <= epsilon ( y ) ) then

      res = - log ( y )
!
!  EPS < X <= 1.5.
!
    else if ( y <= 1.5D+00 ) then

      if ( y < 0.6796875D+00 ) then
        corr = -log ( y )
        xm1 = y
      else
        corr = 0.0D+00
        xm1 = ( y - 0.5D+00 ) - 0.5D+00
      end if

      if ( y <= 0.5D+00 .or. 0.6796875D+00 <= y ) then

        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm1 + p1(i)
          xden = xden * xm1 + q1(i)
        end do

        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

      else

        xm2 = ( y - 0.5D+00 ) - 0.5D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

      end if
!
!  1.5 < X <= 4.0.
!
    else if ( y <= 4.0D+00 ) then

      xm2 = y - 2.0D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
!
!  4.0 < X <= 12.0.
!
    else if ( y <= 12.0D+00 ) then

      xm4 = y - 4.0D+00
      xden = -1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm4 + p4(i)
        xden = xden * xm4 + q4(i)
      end do

      res = d4 + xm4 * ( xnum / xden )
!
!  Evaluate for 12 <= argument.
!
    else

      res = 0.0D+00

      if ( y <= frtbig ) then

        res = c(7)
        ysq = y * y

        do i = 1, 6
          res = res / ysq + c(i)
        end do

      end if

      res = res / y
      corr = log ( y )
      res = res + sqrtpi - 0.5D+00 * corr
      res = res + y * ( corr - 1.0D+00 )

    end if
!
!  Return for bad arguments.
!
  else

    res = xinf

  end if
!
!  Final adjustments and return.
!
  r8_gamma_log = res

  return
end
