!********************************************************************************!
!                    DIS-K
!     Function evaluate the modified plasma dispersion function
!     zk12 for integer kappa
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
!
FUNCTION zk12(z,kp)
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(IN) :: z
  INTEGER, INTENT(IN) :: kp
  COMPLEX(KIND=8) :: zk12, ii, suma, fl, fact, fact2
  REAL(KIND=8) :: sqk
  INTEGER :: l
  REAL(KIND=8), PARAMETER :: pi=3.141592653589793D+00
  !
  ii = DCMPLX(0.0,1.0)
  sqk = SQRT(REAL(kp))
  suma = DCMPLX(0.0,0.0)
  !
  IF ((REAL(z) == 0) .AND. (ABS(AIMAG(z)) == sqk)) THEN
    zk12 = (AIMAG(z)/ABS(AIMAG(z)))*ii*(GAMMA(kp+2.5))/(kp**2.5*(kp+2)*GAMMA(kp-0.5))
  ELSE IF ((REAL(z) == 0) .AND. (AIMAG(z) == 0)) THEN
    zk12 = ii*SQRT(pi)*(GAMMA(kp+2.0))/(kp**2.5*GAMMA(kp-0.5))
  ELSE
    fl = (0.5d0-z/(2*ii*sqk))
    DO l=0,kp+1
      suma = suma+FACTORIAL(l+kp+1)/FACTORIAL(l)*fl**l
    ENDDO
    fact2 = ((0.5+z/(2.0*ii*sqk))**(kp+2))*suma/FACTORIAL(kp+1)
    fact = (2.0*ii*SQRT(pi)/(kp**2.5))*FACTORIAL(kp+1)/GAMMA(kp-0.5)&
    *(1.0+z**2/kp)**(-kp-2)
    !
    zk12 = fact*(1.0d0-fact2)
  ENDIF
!
  CONTAINS
  RECURSIVE INTEGER FUNCTION FACTORIAL (k) RESULT (valor)
    INTEGER :: k
    IF (k<=1) THEN
      valor = 1
    ELSE
      valor = k*FACTORIAL(k-1)
    ENDIF
  END FUNCTION FACTORIAL
END FUNCTION zk12
!********************************************************************************!
