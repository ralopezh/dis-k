!********************************************************************************!
!                    DIS-K
!     Function compute the determinant of the dispersion tensor
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
!
FUNCTION disp(w,k,th)
  USE params
  USE tensor_elements
  USE integral_params
  !
  IMPLICIT NONE
  COMPLEX(KIND=8) :: w, disp, z0, zpad16, xi, zeta, zn, zplas, dispfunct
  COMPLEX(KIND=8) :: An, Bn, Cn, wx2
  REAL(KIND=8) :: k, th, theta, lambda, yac, fact, fact3, Aa, vd
  REAL(KIND=8) :: sinth, costh, Lambdan, Lambdanp1, dLambdan
  INTEGER, ALLOCATABLE, DIMENSION(:) :: narr
  INTEGER :: i, is, n, kp, nm, ni, nmax, msta1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bInex, di, bk, dkn
  COMPLEX(KIND=8) :: sum11, sum22, sum33, sum12, sum13, sum23, zk12
  COMPLEX(KIND=8) :: fi11, fi12, fi13, fi22, fi23, fi33
  REAL(KIND=8), PARAMETER :: pi = 3.141592653589793D+00
  !
  theta = 2.0d0*pi*th/360.0d0
  sinth = SIN(theta)
  costh = COS(theta)
!
  sum11 = DCMPLX(0.0d0,0.0d0)
  sum22 = DCMPLX(0.0d0,0.0d0)
  sum33 = DCMPLX(0.0d0,0.0d0)
  sum12 = DCMPLX(0.0d0,0.0d0)
  sum13 = DCMPLX(0.0d0,0.0d0)
  sum23 = DCMPLX(0.0d0,0.0d0)
!
  DO is=1,nsp
    kp = INT(kappa(is))
    xi = (w-wpwc*k*costh*Ua(is))/(k*alphapal(is)*costh*wpwc)
    lambda = 0.5*(wpwc*rm(is)*k*alphaper(is)*sinth/rq(is))**2
    yac = k*alphapal(is)*costh*wpwc
    fact = ns(is)*rq(is)**2/rm(is)
    fact3 = (rq(is)/ABS(rq(is)))*fact/SQRT(anis(is))
    Aa = 1.0-anis(is)
    vd = Ua(is)/alphapal(is)
    !
    IF (ISNAN(kappa(is))) THEN
      ! Maxwellian case
      nmax = nsummax
      nm = 0
      ALLOCATE(bInex(0:nmax+1),di(0:nmax+1))
      CALL iknae(nmax,lambda,nm,bInex,di) !In(lambda)*Exp(-lambda)
      IF (ABS(lambda) <= 708) THEN
        nmax = 0
        DO WHILE (bInex(nmax) >= Bessel_min .AND. nmax < nsummax)
          nmax = nmax+1
        ENDDO
      ENDIF
      ! summation order [nmax,-nmax, nmax-1, -(nmax-1),...,2,-2,1,-1,0]
      ALLOCATE(narr(2*nmax+1))
      DO i=1,2*nmax,2
        narr(i) = (2*nmax-(i-1))/2
        narr(i+1) = -(2*nmax-(i-1))/2
      ENDDO
      narr(2*nmax+1) = 0
      !
      DO ni=1,2*nmax+1
        n = narr(ni)
        zeta = xi-n*rq(is)/(rm(is)*yac)
        zn = zpad16(zeta)
        !zn = zplas(zeta)
        Lambdan = bInex(ABS(n))
        Lambdanp1 = bInex(ABS(n+1))
        dLambdan = (n/lambda-1.0)*Lambdan + Lambdanp1
        An = anis(is)-1.0 + (xi+(anis(is)-1.0)*zeta)*zn
        Bn = xi+(zeta+vd)*An
        Cn = xi*zeta+An*(zeta+vd)**2
        sum11 = sum11 + fact*n*n*(Lambdan/lambda)*An
        sum22 = sum22 + fact*(n*n*Lambdan/lambda -2.0*lambda*dLambdan)*An
        sum12 = sum12 + fact*n*dLambdan*An
        sum13 = sum13 + fact3*n*SQRT(2.0)*(Lambdan/SQRT(lambda))*Bn
        sum23 = sum23 - fact3*SQRT(2.0*lambda)*dLambdan*Bn
        sum33 = sum33 + (fact/anis(is))*2.0*Lambdan*Cn
      ENDDO
      sum33 = sum33 + (fact/anis(is))*2.0*vd*(vd+2*xi)
      DEALLOCATE(narr,bInex,di)
    ELSE
      ! KAPPA case
      nmax = nsum
      ! summation order [nmax,-nmax, nmax-1, -(nmax-1),...,2,-2,1,-1,0]
      ALLOCATE(narr(2*nmax+1))
      DO i=1,2*nmax,2
        narr(i) = (2*nmax-(i-1))/2
        narr(i+1) = -(2*nmax-(i-1))/2
      ENDDO
      narr(2*nmax+1) = 0
    !
      DO ni=1,2*nmax+1
        n = narr(ni)
        zeta = xi-n*rq(is)/(rm(is)*yac)
        xia = xi
        ztna = zeta
        la = lambda
        aniso = Aa
        vdrft = vd
        n1 = n
        k1 = kp
        CALL i11(fi11)
        CALL i12(fi12)
        CALL i13(fi13)
        CALL i22(fi22)
        CALL i23(fi23)
        CALL i33(fi33)
        sum11 = sum11 + fact*(n*n/lambda)*fi11
        sum22 = sum22 + fact*fi22
        sum12 = sum12 + fact*n*fi12
        sum13 = sum13 + fact3*(n/SQRT(2.0*lambda))*fi13
        sum23 = sum23 + fact3*SQRT(0.5*lambda)*fi23
        sum33 = sum33 + 2.0*(fact/anis(is))*fi33
      ENDDO
      sum33 = sum33 + 2.0*(fact/anis(is))*(1.0-0.5/kp)*vd*vd
    DEALLOCATE(narr)
    ENDIF
  ENDDO
!
  !*** Dispersion ****;
  wx2=wpwc**2/w**2
  d11 = 1.0+wx2*(-k**2*costh**2+sum11)
  d22 = 1.0+wx2*(-k**2+sum22)
  d12 = wx2*sum12
  d13 = wx2*(k**2*sinth*costh+sum13)
  d23 = wx2*sum23
  d33 = 1.0+wx2*(-k**2*sinth**2+sum33)
!
  disp = d11*d22*d33-d11*d23**2-d33*d12**2-d22*d13**2-2.0*d12*d13*d23
!  
END FUNCTION disp
!********************************************************************************!
