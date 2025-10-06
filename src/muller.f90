!********************************************************************************!
!                    DIS-K
!     Nonlinear root_solver based on Muller method
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
FUNCTION muller(xi,yi,thi,itmax,tol)
!
! Nonlinear Root solver using Muller method
! From Barrodale and Wilson, J. Comput. Appl. Math. 4 (1978)
! 
  IMPLICIT NONE
  COMPLEX(KIND=8), DIMENSION(3) :: xi, x 
  COMPLEX(KIND=8) :: muller, disp, f1, f2, f3
  COMPLEX(KIND=8) :: evalf
  COMPLEX(KIND=8) :: w, fx3x2, fx3x1, fx2x1, fx3x2x1, sq, den, den1, den2
  REAL(KIND=8) :: yi, thi, tol
  INTEGER :: cond, it, itmax
  !
  x = xi+0.0d0
  IF (SIZE(x) /= 3) PRINT*,'ROOTS: x must be a 3-element initial guess vector.'
  !
  !Initialize stopping criterion and iteration count.
  cond = 0
  it = 0
  !Begin to iteratively compute a root of the nonlinear function.
  f1 = disp(x(1),yi,thi)
  f2 = disp(x(2),yi,thi)
  f3 = disp(x(3),yi,thi)
  DO WHILE (it < itmax .AND. cond /= 1)
    !
    fx3x2 = (f2-f3)/(x(2)-x(3))
    fx3x1 = (f1-f3)/(x(1)-x(3))
    fx2x1 = (f1-f2)/(x(1)-x(2))
    w = fx3x2+fx3x1-fx2x1
    fx3x2x1 = (fx2x1-fx3x2)/(x(1)-x(3))
    !
    sq = ZSQRT(w*w-4.0*f3*fx3x2x1)
    den1 = w+sq
    den2 = w-sq
    !
    IF (ZABS(den1) >= ZABS(den2)) THEN
      den = den1
    ELSE 
      den = den2
    ENDIF
    !
    IF (ZABS(den) == 0.0) THEN
      muller = 2.0*x(3)-x(2)
    ELSE
      muller = x(3) - 2.0*f3/den
    ENDIF
    !
    IF (ZABS(muller-x(3)) < tol) THEN !Absolute error tolerance.
      cond = 1
    ELSE
      evalf = disp(muller,yi,thi) !Functional error tolerance.
      IF (ZABS(evalf) < tol) THEN
        cond = 1
      ELSE IF (evalf == 0.0) THEN
        cond = 1
      ENDIF
    ENDIF
    x(1) = x(2)
    f1 = f2
    x(2) = x(3)
    f2 = f3
    x(3) = muller
    f3 = disp(x(3),yi,thi)
    it = it+1
  ENDDO
  IF (it > itmax .AND. cond==0) PRINT*,'ROOT: Algorithm failed to converge within given parameters.'  
END FUNCTION muller
!********************************************************************************!

