!********************************************************************************!
!                    DIS-K
!     Modules
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
MODULE params
  IMPLICIT NONE
  INTEGER :: nsp, np, nth, nroots, nsum, nsummax, nw, logscale_k, logscale_th
  REAL(KIND=8) :: rpem, wpewce, angle, wpwc
  CHARACTER(80) :: version, sourcefile, subindex, home
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: rm, ns, rq, betapal, anis, kappa, Ua
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: alphapal, alphaper
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: w0r, w0i, k0, th0
  REAL(KIND=8) :: kmin, kmax, thmin, thmax, Bessel_min
  REAL(KIND=8), DIMENSION(2) :: wr_range, wi_range, k_range, th_range
END MODULE params
!
MODULE solve_params
  IMPLICIT NONE
  INTEGER :: damped, itmax, interpolate, prevmax
  REAL(KIND=8) :: gammatol, tol, epsi
END MODULE solve_params
!
MODULE tensor_elements
  IMPLICIT NONE
  COMPLEX(KIND=8) :: d11, d12, d13, d22, d23, d33
END MODULE tensor_elements
!
MODULE integral_params
  COMPLEX(KIND=8) :: xia, ztna
  REAL(KIND=8) :: la, aniso, vdrft
  INTEGER :: n1, k1
!$OMP THREADPRIVATE(xia,ztna,la,aniso,vdrft,n1,k1)
END MODULE
