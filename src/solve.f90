!********************************************************************************!
!                    DIS-K
!     Routine to solve the dispersion relation for all elements in the vector k
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
SUBROUTINE solve(kvec,k0,thi,z0,np,nroots,wsol,pols,ratios)
!
! Solve the dispersion relation ABS(\Lambda) = 0
! The complex root is found using an optimal
! Muller’s method by FX_ROOT
!
  USE solve_params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: np, nroots
  REAL(KIND=8), INTENT(IN), DIMENSION(np) :: kvec
  COMPLEX(KIND=8), INTENT(IN), DIMENSION(nroots) :: z0
  COMPLEX(KIND=8), DIMENSION(nroots) :: zc
  COMPLEX(KIND=8), INTENT(OUT), DIMENSION(np,nroots) :: wsol, pols
  REAL(KIND=8), INTENT(OUT), DIMENSION(np,nroots,8) :: ratios
  REAL(KIND=8), DIMENSION(8) :: rats
  REAL(KIND=8), DIMENSION(nroots), INTENT(IN) :: k0, thi
  REAL(KIND=8) :: ki, cubic_fit, poly_real, poly_imag
  COMPLEX(KIND=8), DIMENSION(3) :: seed
  COMPLEX(KIND=8) :: root_finder, muller, newseed, polwi
  INTEGER ::  ncount, i, i0, nex
!
  wsol = 0.0
  pols = 0.0
  ratios = 0.0
  !
  DO nex=1,nroots
    IF (np==1) THEN
      i0 = 1
    ELSE
      i0 = MINLOC(ABS(kvec-k0(nex)),DIM=1)
    ENDIF
    !
    ncount = 0
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(i,ki,seed,newseed,zc,polwi,rats,poly_real,poly_imag)
!$OMP SECTION
    zc=z0
    !DO i=i0,1,-1
    i = i0
    DO WHILE (i>=1)
      ki = kvec(i)
      IF ((ncount > 4) .AND. (interpolate == 1)) THEN
        poly_real = cubic_fit(kvec(i+1:i+4),REAL(wsol(i+1:i+4,nex)),ki) !use the previous 4 solutions to extrapolate
        poly_imag = cubic_fit(kvec(i+1:i+4),AIMAG(wsol(i+1:i+4,nex)),ki)
        zc(nex) = DCMPLX(poly_real,poly_imag)
      ENDIF
      ncount = ncount+1
      seed(1) = zc(nex) - DCMPLX(epsi,epsi)
      seed(2) = zc(nex)
      seed(3) = zc(nex) + DCMPLX(epsi,epsi)
      !newseed = root_finder(seed,ki,thi(nex),itmax,tol)
      newseed = muller(seed,ki,thi(nex),itmax,tol)
      !print*,ki,thi(nex),newseed
      !IF (ZABS(newseed) >= 5.0*ZABS(zc(nex))) newseed = 0.0*newseed
      IF (AIMAG(newseed)/AIMAG(zc(nex)) >= 10.0) newseed = 0.0*newseed
      !IF (ABS(REAL(newseed)) <= 1.0e-4) newseed = 0.0*newseed
      CALL polarization(newseed,ki,thi(nex),polwi,rats)
      wsol(i,nex) = newseed
      pols(i,nex) = polwi
      IF (ISNAN(REAL(polwi)) .OR. ISNAN(AIMAG(polwi))) pols(i,nex) = 0.0
      ratios(i,nex,:) = rats
      zc(nex) = newseed
      IF ((ISNAN(AIMAG(newseed))) .OR. (damped == 0 .AND. AIMAG(newseed) < gammatol) ) THEN
        zc(nex) = DCMPLX(0.0,0.0)
        wsol(i,nex) = DCMPLX(0.0,0.0)
        !EXIT
        i = 0
      ENDIF
      !PRINT*,ki,REAL(newseed),AIMAG(newseed)
      i = i-1
    ENDDO
!$OMP SECTION
    zc = z0
    !DO i=i0+1,np
    i = i0+1
    DO WHILE (i<=np)
      ki = kvec(i)
      IF ((ncount > 4) .AND. (interpolate == 1)) THEN
        poly_real = cubic_fit(kvec(i-4:i-1),REAL(wsol(i-4:i-1,nex)),ki) !use the previous 4 solutions to extrapolate
        poly_imag = cubic_fit(kvec(i-4:i-1),AIMAG(wsol(i-4:i-1,nex)),ki)
        zc(nex) = DCMPLX(poly_real,poly_imag)
      ENDIF
      ncount = ncount+1
      seed(1) = zc(nex) - DCMPLX(epsi,epsi)
      seed(2) = zc(nex)
      seed(3) = zc(nex) + DCMPLX(epsi,epsi)
      !newseed = root_finder(seed,ki,thi(nex),itmax,tol)
      newseed = muller(seed,ki,thi(nex),itmax,tol)
      !print*,ki,thi(nex),newseed
      !IF (ZABS(newseed) >= 5.0*ZABS(zc(nex))) newseed = 0.0*newseed
      IF (AIMAG(newseed)/AIMAG(zc(nex)) >= 10.0) newseed = 0.0*newseed
      !IF (ABS(REAL(newseed)) <= 1.0e-4) newseed = 0.0*newseed
      CALL polarization(newseed,ki,thi(nex),polwi,rats)
      wsol(i,nex) = newseed
      pols(i,nex) = polwi
      IF (ISNAN(REAL(polwi)) .OR. ISNAN(AIMAG(polwi))) pols(i,nex) = 0.0
      ratios(i,nex,:) = rats
      zc(nex) = newseed
      IF ((ISNAN(AIMAG(newseed))) .OR. (damped == 0 .AND. AIMAG(newseed) < gammatol)) THEN
        zc(nex) = DCMPLX(0.0,0.0)
        wsol(i,nex) = DCMPLX(0.0,0.0)
        !EXIT
        i = np+1
      ENDIF
      !PRINT*,ki,REAL(newseed),AIMAG(newseed)
      i = i+1
    ENDDO
!$OMP END PARALLEL SECTIONS
  ENDDO
!  
END SUBROUTINE solve
!********************************************************************************!
!********************************************************************************!
SUBROUTINE solve1p(ki,thi,z0)
!
! Solve the dispersion relation ABS(\Lambda) = 0
! For a single point 
! The complex root is found using an optimal
! Muller’s method by FX_ROOT
!
  USE solve_params
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: ki, thi
  COMPLEX(KIND=8), INTENT(INOUT) :: z0
  COMPLEX(KIND=8), DIMENSION(3) :: seed
  COMPLEX(KIND=8) :: root_finder, muller, newseed
!
  seed(1) = z0 - DCMPLX(epsi,epsi)
  seed(2) = z0
  seed(3) = z0 + DCMPLX(epsi,epsi)
  !newseed = root_finder(seed,ki,thi,itmax,tol)
  newseed = muller(seed,ki,thi,itmax,tol)
  z0 = newseed
!
END SUBROUTINE solve1p
!********************************************************************************!
