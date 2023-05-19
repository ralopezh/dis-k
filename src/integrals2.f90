!********************************************************************************!
!                    DIS-K
!     Function compute the integrals I11, I12, I13, I22, I23, I33
!     for the kappa case
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
!
SUBROUTINE i11(result)
  ! Integral I_11
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(INOUT) :: result
  REAL(KIND=8), EXTERNAL :: i11r, i11i
  REAL(KIND=8) :: intr, inti
  REAL(KIND=8), PARAMETER :: epsabs = 0.0d0
  REAL(KIND=8), PARAMETER :: epsrel = 1.0d-06
  REAL(KIND=8) :: err
  INTEGER :: neval, ier, last
  INTEGER, PARAMETER :: limit = 100
  INTEGER, PARAMETER :: lenw = limit*4
  INTEGER, DIMENSION(limit) :: iwork
  REAL(KIND=8), DIMENSION(lenw) :: work
  !
  CALL dqagi(i11r,0.0d0,1,epsabs,epsrel,intr,err,neval,ier,limit,lenw,last,iwork,work)
  CALL dqagi(i11i,0.0d0,1,epsabs,epsrel,inti,err,neval,ier,limit,lenw,last,iwork,work)
  !
  result = DCMPLX(intr,inti)
END SUBROUTINE i11
!********************************************************************************!
SUBROUTINE i12(result)
  ! Integral I_12
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(INOUT) :: result
  REAL(KIND=8), EXTERNAL :: i12r, i12i
  REAL(KIND=8) :: intr, inti
  REAL(KIND=8), PARAMETER :: epsabs = 0.0d0
  REAL(KIND=8), PARAMETER :: epsrel = 1.0d-06
  REAL(KIND=8) :: err
  INTEGER :: neval, ier, last
  INTEGER, PARAMETER :: limit = 100
  INTEGER, PARAMETER :: lenw = limit*4
  INTEGER, DIMENSION(limit) :: iwork
  REAL(KIND=8), DIMENSION(lenw) :: work
  !
  CALL dqagi(i12r,0.0d0,1,epsabs,epsrel,intr,err,neval,ier,limit,lenw,last,iwork,work)
  CALL dqagi(i12i,0.0d0,1,epsabs,epsrel,inti,err,neval,ier,limit,lenw,last,iwork,work)
  !
  result = DCMPLX(intr,inti)
END SUBROUTINE i12
!********************************************************************************!
SUBROUTINE i13(result)
  ! Integral I_13
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(INOUT) :: result
  REAL(KIND=8), EXTERNAL :: i13r, i13i
  REAL(KIND=8) :: intr, inti
  REAL(KIND=8), PARAMETER :: epsabs = 0.0d0
  REAL(KIND=8), PARAMETER :: epsrel = 1.0d-06  
  REAL(KIND=8) :: err
  INTEGER :: neval, ier, last
  INTEGER, PARAMETER :: limit = 100
  INTEGER, PARAMETER :: lenw = limit*4
  INTEGER, DIMENSION(limit) :: iwork
  REAL(KIND=8), DIMENSION(lenw) :: work
  !
  CALL dqagi(i13r,0.0d0,1,epsabs,epsrel,intr,err,neval,ier,limit,lenw,last,iwork,work)
  CALL dqagi(i13i,0.0d0,1,epsabs,epsrel,inti,err,neval,ier,limit,lenw,last,iwork,work)
  !
  result = DCMPLX(intr,inti)
END SUBROUTINE i13
!********************************************************************************!
SUBROUTINE i22(result)
  ! Integral I_22
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(INOUT) :: result
  REAL(KIND=8), EXTERNAL :: i22r, i22i
  REAL(KIND=8) :: intr, inti
  REAL(KIND=8), PARAMETER :: epsabs = 0.0d0
  REAL(KIND=8), PARAMETER :: epsrel = 1.0d-06  
  REAL(KIND=8) :: err
  INTEGER :: neval, ier, last
  INTEGER, PARAMETER :: limit = 100
  INTEGER, PARAMETER :: lenw = limit*4
  INTEGER, DIMENSION(limit) :: iwork
  REAL(KIND=8), DIMENSION(lenw) :: work
  !
  CALL dqagi(i22r,0.0d0,1,epsabs,epsrel,intr,err,neval,ier,limit,lenw,last,iwork,work)
  CALL dqagi(i22i,0.0d0,1,epsabs,epsrel,inti,err,neval,ier,limit,lenw,last,iwork,work)
  !
  result = DCMPLX(intr,inti)
END SUBROUTINE i22
!********************************************************************************!
SUBROUTINE i23(result)
  ! Integral I_23
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(INOUT) :: result
  REAL(KIND=8), EXTERNAL :: i23r, i23i
  REAL(KIND=8) :: intr, inti
  REAL(KIND=8), PARAMETER :: epsabs = 0.0d0
  REAL(KIND=8), PARAMETER :: epsrel = 1.0d-06  
  REAL(KIND=8) :: err
  INTEGER :: neval, ier, last
  INTEGER, PARAMETER :: limit = 100
  INTEGER, PARAMETER :: lenw = limit*4
  INTEGER, DIMENSION(limit) :: iwork
  REAL(KIND=8), DIMENSION(lenw) :: work
  !
  CALL dqagi(i23r,0.0d0,1,epsabs,epsrel,intr,err,neval,ier,limit,lenw,last,iwork,work)
  CALL dqagi(i23i,0.0d0,1,epsabs,epsrel,inti,err,neval,ier,limit,lenw,last,iwork,work)
  !
  result = DCMPLX(intr,inti)
END SUBROUTINE i23
!********************************************************************************!
SUBROUTINE i33(result)
  ! Integral I_33
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(INOUT) :: result
  REAL(KIND=8), EXTERNAL :: i33r, i33i
  REAL(KIND=8) :: intr, inti
  REAL(KIND=8), PARAMETER :: epsabs = 0.0d0
  REAL(KIND=8), PARAMETER :: epsrel = 1.0d-06  
  REAL(KIND=8) :: err
  INTEGER :: neval, ier, last
  INTEGER, PARAMETER :: limit = 100
  INTEGER, PARAMETER :: lenw = limit*4
  INTEGER, DIMENSION(limit) :: iwork
  REAL(KIND=8), DIMENSION(lenw) :: work
  !
  CALL dqagi(i33r,0.0d0,1,epsabs,epsrel,intr,err,neval,ier,limit,lenw,last,iwork,work)
  CALL dqagi(i33i,0.0d0,1,epsabs,epsrel,inti,err,neval,ier,limit,lenw,last,iwork,work)
  !
  result = DCMPLX(intr,inti)
END SUBROUTINE i33
!********************************************************************************!
!********************************************************************************!
!********************************************************************************!
!********************************************************************************!
!********************************************************************************!
!********************************************************************************!
!                    DIS-K
!     Function evaluate the argument in the
!     integrals I11, I12, I13, I22, I23, I33
!     for the kappa case
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
FUNCTION i11r(x)
  ! REAL Argument of integral I_11
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg
  REAL(KIND=8) i11r, den, sqden, sql
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  arg = 2.0*x*BESSEL_JN(n1,x*sql)**2/(den**(k1+1.5))*Han
  i11r = REAL(arg)
END FUNCTION i11r
FUNCTION i11i(x)
  ! IMAGINARY Argument of integral I_11
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg
  REAL(KIND=8) i11i, den, sqden, sql
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  arg = 2.0*x*BESSEL_JN(n1,x*sql)**2/(den**(k1+1.5))*Han
  i11i = AIMAG(arg)
END FUNCTION i11i
!********************************************************************************!
FUNCTION i12r(x)
  ! REAL Argument of integral I_12
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg
  REAL(KIND=8) i12r, den, sqden, sql, hss
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  hss = x**2*BESSEL_JN(n1,x*sql)*(BESSEL_JN(n1-1,x*sql)-BESSEL_JN(n1+1,x*sql))
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  arg = (2.0/sql)*(hss/den**(k1+1.5))*Han
  i12r = REAL(arg)
END FUNCTION i12r
FUNCTION i12i(x)
  ! IMAGINARY Argument of integral I_12
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg
  REAL(KIND=8) i12i, den, sqden, sql, hss
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  hss = x**2*BESSEL_JN(n1,x*sql)*(BESSEL_JN(n1-1,x*sql)-BESSEL_JN(n1+1,x*sql))
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  arg = (2.0/sql)*(hss/den**(k1+1.5))*Han
  i12i = AIMAG(arg)
END FUNCTION i12i
!********************************************************************************!
FUNCTION i13r(x)
  ! REAL Argument of integral I_13
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg, Kan
  REAL(KIND=8) i13r, den, sqden, sql
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  Kan = (1.0-0.25/(k1*k1))*xia+(ztna+vdrft)*Han
  arg = 4.0*x*BESSEL_JN(n1,x*sql)**2/(den**(k1+1.5))*Kan
  i13r = REAL(arg)
END FUNCTION i13r
FUNCTION i13i(x)
  ! IMAGINARY Argument of integral I_13
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg, Kan
  REAL(KIND=8) i13i, den, sqden, sql
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  Kan = (1.0-0.25/(k1*k1))*xia+(ztna+vdrft)*Han
  arg = 4.0*x*BESSEL_JN(n1,x*sql)**2/(den**(k1+1.5))*Kan
  i13i = AIMAG(arg)
END FUNCTION i13i
!********************************************************************************!
FUNCTION i23r(x)
  ! REAL Argument of integral I_23
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg, Kan
  REAL(KIND=8) i23r, den, sqden, sql, hss
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  hss = x**2*BESSEL_JN(n1,x*sql)*(BESSEL_JN(n1-1,x*sql)-BESSEL_JN(n1+1,x*sql))
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  Kan = (1.0-0.25/(k1*k1))*xia+(ztna+vdrft)*Han
  arg = -(4.0/sql)*hss/(den**(k1+1.5))*Kan
  i23r = REAL(arg)
END FUNCTION i23r
FUNCTION i23i(x)
  ! IMAGINARY Argument of integral I_23
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg, Kan
  REAL(KIND=8) i23i, den, sqden, sql, hss
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  hss = x**2*BESSEL_JN(n1,x*sql)*(BESSEL_JN(n1-1,x*sql)-BESSEL_JN(n1+1,x*sql))
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  Kan = (1.0-0.25/(k1*k1))*xia+(ztna+vdrft)*Han
  arg = -(4.0/sql)*hss/(den**(k1+1.5))*Kan
  i23i = AIMAG(arg)
END FUNCTION i23i
!********************************************************************************!
FUNCTION i22r(x)
  ! REAL Argument of integral I_22
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg
  REAL(KIND=8) i22r, den, sqden, sql, hss
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  hss = (n1*n1/la)*BESSEL_JN(n1,x*sql)**2&
&   -2.0*x*x*BESSEL_JN(n1-1,x*sql)*BESSEL_JN(n1+1,x*sql)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  arg = 2.0*x*(hss/den**(k1+1.5))*Han
  i22r = REAL(arg)
END FUNCTION i22r
FUNCTION i22i(x)
  ! IMAGINARY Argument of integral I_22
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg
  REAL(KIND=8) i22i, den, sqden, sql, hss
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  hss = (n1*n1/la)*BESSEL_JN(n1,x*sql)**2&
&   -2.0*x*x*BESSEL_JN(n1-1,x*sql)*BESSEL_JN(n1+1,x*sql)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  arg = 2.0*x*(hss/den**(k1+1.5))*Han
  i22i = AIMAG(arg)
END FUNCTION i22i
!********************************************************************************!
FUNCTION i33r(x)
  ! REAL Argument of integral I_33
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg, Qan
  REAL(KIND=8) i33r, den, sqden, sql
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  Qan = xia*(1.0-0.25/(k1*k1))*(ztna+2.0*vdrft)+(ztna+vdrft)**2*Han
  arg = 2.0*x*BESSEL_JN(n1,x*sql)**2/(den**(k1+1.5))*Qan
  i33r = REAL(arg)
END FUNCTION i33r
FUNCTION i33i(x)
  ! IMAGINARY Argument of integral I_33
  USE integral_params
  !
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: x
  COMPLEX(KIND=8) :: Han, zk12, arg, Qan
  REAL(KIND=8) i33i, den, sqden, sql
  !
  den = 1.0+x*x/k1
  sqden = SQRT(den)
  sql = SQRT(2.0*la)
  Han = (xia-aniso*ztna)*zk12(ztna/sqden,k1)/sqden-(1.0-0.25/(k1*k1))*aniso
  Qan = xia*(1.0-0.25/(k1*k1))*(ztna+2.0*vdrft)+(ztna+vdrft)**2*Han
  arg = 2.0*x*BESSEL_JN(n1,x*sql)**2/(den**(k1+1.5))*Qan
  i33i = AIMAG(arg)
END FUNCTION i33i
!********************************************************************************!

