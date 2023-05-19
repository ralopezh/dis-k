!********************************************************************************!
!                    DIS-K
!     Function compute a cubic fit using 4 points
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
!
FUNCTION cubic_fit(xi,yi,xj)
  !
  ! Compute the coeficients of the polinomia
  ! p(x) = c1 + c2*xi + c3*xi**2 + c4*xi**4
  ! using 4 points
  !
  ! xi : abscissa array of dimension 4
  ! yi : ordinate array of dimension 4
  ! xj : abscissa value to be extrapolated
  !
  IMPLICIT NONE
  !
  REAL(KIND=8), INTENT(IN), DIMENSION(4) :: xi, yi
  REAL(KIND=8), INTENT(IN) :: xj
  REAL(KIND=8), DIMENSION(4,4) :: A, Ainv
  REAL(KIND=8), DIMENSION(4) :: coef
  REAL(KIND=8) :: cubic_fit, det, detinv
  INTEGER :: i, j
  !
  DO i=1,4
    A(i,1) = 1.0
    A(i,2) = xi(i)
    A(i,3) = xi(i)*xi(i)
    A(i,4) = xi(i)*xi(i)*xi(i)
  ENDDO
  !
  det = (xi(4)-xi(1))*(xi(3)-xi(1))*(xi(2)-xi(1))*(xi(4)-xi(2))&
&   *(xi(3)-xi(2))*(xi(4)-xi(3))
  !
  IF (det == 0) THEN
    cubic_fit = yi(4)
  ELSE
    detinv = 1.0/det
    ! Calculate the inverse of the matrix
    Ainv(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
&     +A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    Ainv(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
&     +A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    Ainv(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))&
&     +A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    Ainv(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))&
&     +A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    Ainv(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
&     +A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    Ainv(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
&     +A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    Ainv(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))&
&     +A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    Ainv(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))&
&     +A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    Ainv(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))&
&     +A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    Ainv(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))&
&     +A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    Ainv(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))&
&     +A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    Ainv(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))&
&     +A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    Ainv(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))&
&     +A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    Ainv(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))&
&     +A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    Ainv(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))&
&     +A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    Ainv(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
&     +A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    !
    coef(1) = yi(1)*Ainv(1,1)+yi(2)*Ainv(1,2)+yi(3)*Ainv(1,3)+yi(4)*Ainv(1,4)
    coef(2) = yi(1)*Ainv(2,1)+yi(2)*Ainv(2,2)+yi(3)*Ainv(2,3)+yi(4)*Ainv(2,4)
    coef(3) = yi(1)*Ainv(3,1)+yi(2)*Ainv(3,2)+yi(3)*Ainv(3,3)+yi(4)*Ainv(3,4)
    coef(4) = yi(1)*Ainv(4,1)+yi(2)*Ainv(4,2)+yi(3)*Ainv(4,3)+yi(4)*Ainv(4,4)
    !
    cubic_fit = coef(1)+coef(2)*xj+coef(3)*xj**2+coef(4)*xj**3
  ENDIF
    !
END FUNCTION cubic_fit
