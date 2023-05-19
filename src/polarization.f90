!********************************************************************************!
!                    DIS-K
!     Routine to compute the polarization and fields ratios
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
!
SUBROUTINE polarization(w,ki,thi,pol,rat)
  USE params
  USE tensor_elements
  !
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(IN) :: w
  REAL(KIND=8), INTENT(IN) :: ki,thi
  COMPLEX(KIND=8), INTENT(OUT) :: pol
  REAL(KIND=8), INTENT(OUT), DIMENSION(8) :: rat
  COMPLEX(KIND=8) :: ii, disp, xd, R, T
  REAL(KIND=8) :: ex2, ey2, ez2, bx2, by2, bz2
  REAL(KIND=8) :: wr, costh, sinth
  REAL(KIND=8), PARAMETER :: pi = 3.141592653589793D+00
  !
  ii = DCMPLX(0.0d0,1.0d0)
  costh = cos(2.0*pi*thi/360.0d0)
  sinth = sin(2.0*pi*thi/360.0d0)
  !
  ! call disp
  xd = disp(w,ki,thi)
  !
  ! Compute the ratios R=Ex/Ey, T=Ez/Ey
  !
  R = ii*(d12*d33+d13*d23)/(d13*d13-d11*d33)
  T = ii*(d11*d23+d12*d13)/(d11*d33-d13*d13)
  wr = REAL(w)
  pol = ii*R*wr/ABS(wr)
!
  ex2 = ABS(R)**2/(1.0d0+ABS(R)**2+ABS(T)**2) ! |Ex/E|^2
  ey2 = 1.0d0/(1.0d0+ABS(R)**2+ABS(T)**2)     ! |Ey/E|^2
  ez2 = ABS(T)**2/(1.0d0+ABS(R)**2+ABS(T)**2) ! |Ez/E|^2
  !
  bx2 = 1.0/(1./costh**2+ABS(R)**2+(sinth**2/costh**2)*ABS(T)**2 &
&   -(sinth/costh)*2.0*REAL(CONJG(R)*T))  ! |Bx/B|^2
  by2 = (costh**2*ABS(R)**2+sinth**2*ABS(T)**2-costh*sinth*2.0*REAL(CONJG(R)*T))&
&   /(1.0+costh**2*ABS(R)**2+sinth**2*ABS(T)**2-costh*sinth*2.0*REAL(CONJG(R)*T)) ! |By/B|^2
  bz2 = 1.0/(1./sinth**2+(costh**2/sinth**2)*ABS(R)**2+ABS(T)**2 &
&   -(costh/sinth)*2.0*REAL(CONJG(R)*T))  ! |Bz/B|^2
  !
  rat(1) = ex2
  rat(2) = ey2
  rat(3) = ez2
  rat(4) = bx2
  rat(5) = by2
  rat(6) = bz2
  rat(7) = R
  rat(8) = T
!
END SUBROUTINE polarization
!********************************************************************************!
