!********************************************************************************!
!                    DIS-K
!     Plasma dispersion function adapted from J.A. Araneda IDL code
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
FUNCTION zpad16(z0)
  !*************************************************************************
  !*                                                                       *
  !*                     PLASMA DISPERSION FUNCTION                        *
  !*                                                                       *
  !*************************************************************************
  !*                >> This Version: zpad16m.pro  <<                       *
  !*                                                                       *
  !*          Author: Jaime Araneda                                        *
  !*                  Max-Planck-Institut fuer Aeronomie                   *
  !*                  37191 Katlenburg-Lindau, Germany                     *
  !*                  and                                                  *
  !*                  Universidad de Concepcion                            *
  !*                  Facultad de Ciencias Fisicas y Matematicas           *
  !*                  Cas. 160-C, Concepcion, Chile.                       *
  !*                                                                       *
  !*************************************************************************
  IMPLICIT NONE
  COMPLEX(KIND=8) :: zpad16, z, zup, zdw, s, z0
  COMPLEX(KIND=8) :: a0, a1, a2, a3, a4, a5, a6, a7, a8
  COMPLEX(KIND=8) :: a9, a10, a11, a12, a13, a14, a15
  COMPLEX(KIND=8) :: b0, b1, b2, b3, b4, b5, b6, b7, b8
  COMPLEX(KIND=8) :: b9, b10, b11, b12, b13, b14
  REAL(KIND=8) :: x0, y0, yi, sigma
  REAL(KIND=8), PARAMETER :: pi = 3.141592653589793D+00
!
  x0 = REAL(z0)
  yi = AIMAG(z0)
  y0 = ABS(yi)
  s = DCMPLX(x0,y0)
!
  a0 = DCMPLX(4.51110739069151373311155837155d7,0.0d0)
  a1 = DCMPLX(0.0d0, -1.96466871099982095478328658013d8)
  a2 = DCMPLX(-4.09230286151779071237033774845d8, 0.0d0)
  a3 = DCMPLX(0.0d0, 5.4164772918858688654664636532d8)
  a4 = DCMPLX(5.1038859260812639084721686937d8, 0.0d0)
  a5 = DCMPLX(0.0d0, -3.63453930635165917551035026585d8)
  a6 = DCMPLX(-2.02579486447453453725874063827d8, 0.0d0)
  a7 = DCMPLX(0.0d0, 9.0274050783922234545623581774d7)
  a8 = DCMPLX(3.25543933925238207202571127297d7, 0.0d0)
  a9 = DCMPLX(0.0d0, -9.5485558117058112533420437956d6)
  a10 = DCMPLX(-2.27494516742535953559592125592d6, 0.0d0)
  a11 = DCMPLX(0.0d0, 436641.4936529503807308097427410d0)
  a12 = DCMPLX(66366.6985055704504521842782210d0, 0.0d0)
  a13 = DCMPLX(0.0d0, -7747.34312765784084062545831270d0)
  a14 = DCMPLX(-657.651665257707244027477678990d0, 0.0d0)
  a15 = DCMPLX(0.0d0, 36.44936596083223830326901197760d0)
!
  b0 = DCMPLX(0.0d0,7.9957296664795077217158359511d7)
  b1 = DCMPLX(2.58006314442690626013556993414d8, 0.0d0)
  b2 = DCMPLX(0.0d0,-4.12365351061717968390257088387d8)
  b3 = DCMPLX(-4.29665394844948144809272592526d8, 0.0d0)
  b4 = DCMPLX(0.0d0,3.24709384845321600507838711434d8)
  b5 = DCMPLX(1.87888347374865872126633533205d8, 0.0d0)
  b6 = DCMPLX(0.0d0,-8.5812966928908801214488359652d7)
  b7 = DCMPLX(-3.14654692983179605892157387897d7, 0.0d0)
  b8 = DCMPLX(0.0d0, 9.3359772296639028831602893886d6)
  b9 = DCMPLX(2.24225318192151759080284972507d6,0.0d0)
  b10= DCMPLX(0.0d0,-432795.1591135920844892244653430d0)
  b11= DCMPLX(-66038.6226729415968301705393810d0,0.0d0)
  b12= DCMPLX(0.0d0,7729.11844467742472147382380670d0)
  b13= DCMPLX(657.15166525770724402747767899d0,0.0d0)
  b14= DCMPLX(0.0d0,-36.44936596083223830326901197760d0)
!
  zup = b0+(b1+b2*s+b3*s**2+b4*s**3+b5*s**4+b6*s**5+b7*s**6)*s &
        +b8*s**8+b9*s**9+b10*s**10+b11*s**11+b12*s**12+b13*s**13+b14*s**14 - s**15
!
  zdw = a0+(a1+a2*s+a3*s**2+a4*s**3+a5*s**4+a6*s**5+a7*s**6)*s &
       +(a8+a9*s+a10*s**2+a11*s**3+a12*s**4+a13*s**5+a14*s**6+a15*s**7)*s**8 + s**16
!
  z = zup/zdw
!  IF (yi < 0.0) THEN
!    sigma = 2.0
!  ELSE IF (yi == 0.0) THEN
!    sigma = 1.0
!  ELSE
!    sigma = 0.0
!  ENDIF
  !
!  IF (ABS(z0) < 10) THEN
    IF (yi < 0.0) THEN
      zpad16 = CONJG(z) + EXP(-z0*z0)*DCMPLX(0.0d0,2.0d0*SQRT(pi))
    ELSE
      zpad16 = z
    ENDIF
!  ELSE
!    zpad16 = -1.0/z0 - 0.5/(z0*z0*z0) - 0.75/(z0*z0*z0*z0) + EXP(-z0*z0)*DCMPLX(0.0d0,sigma*SQRT(pi))
!  ENDIF
!  
END FUNCTION zpad16
!********************************************************************************!
!********************************************************************************!
FUNCTION zplas(z0)
  !*************************************************************************
  !    Plasma dispersion function using
  !    Imaginary Error Function ErfI
  !*************************************************************************
  IMPLICIT NONE
  COMPLEX(KIND=8) :: z0, erf, erfp, erfi, ii, zplas
  REAL(KIND=8), PARAMETER :: pi = 3.141592653589793D+00
  !
  ii = DCMPLX(0.0d0,1.0d0)
  !CALL cerf(ii*z0,erf,erfp)
  !erfi = -ii*erf
  zplas = 0.0!SQRT(pi)*EXP(-z0*z0)*(ii-erfi)
  !
END FUNCTION zplas
!********************************************************************************!


