!********************************************************************************!
!                    DIS-K
!     Routine to read the input data
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
SUBROUTINE load
  USE params
  USE solve_params
  IMPLICIT NONE
  LOGICAL :: existence
!
  NAMELIST /inputdat/version,home,nsp,rpem,wpewce
  NAMELIST /particles/subindex,wpwc,rm,ns,rq,betapal,anis,kappa,Ua
  NAMELIST /wavenumber/np,nth,k_range,th_range,logscale_k,logscale_th
  NAMELIST /roots/nroots
  NAMELIST /seeds/k0,th0,w0r,w0i
  NAMELIST /solver/damped,gammatol,tol,itmax,nsum,nsummax,Bessel_min, epsi, interpolate, prevmax
  NAMELIST /nyq/nw,wr_range,wi_range
!
  
  ! Check whether file exists.
  INQUIRE (file=sourcefile, EXIST=existence)
  IF (existence .EQV. .FALSE.) THEN
    PRINT*,'Error: input file does not exist'
    CALL EXIT
    RETURN
  END IF
  !      
  OPEN(2, file=sourcefile)
  READ(2, inputdat)
  !
  ALLOCATE(rm(nsp),ns(nsp),rq(nsp),betapal(nsp),anis(nsp),kappa(nsp),Ua(nsp))
  ALLOCATE(alphapal(nsp),alphaper(nsp))
  !
  READ(2, particles)
  READ(2, wavenumber)
  READ(2, roots)
  !
  ALLOCATE(k0(nroots),th0(nroots),w0r(nroots),w0i(nroots))
  !
  READ(2, seeds)
  READ(2, solver)
  READ(2, nyq)
  CLOSE(2)
  WRITE(*, inputdat)
  WRITE(*, particles)
  WRITE(*, wavenumber)
  WRITE(*, roots)
  WRITE(*, seeds)
  WRITE(*, solver)
  WRITE(*, nyq)
!
  alphapal = SQRT(betapal/(ns*rm))/wpwc !alpha_{\parallel a}/c
  alphaper = SQRT(anis)*alphapal
  PRINT*, 'ALPHA_PAL/c = ',alphapal

  kmin = k_range(1)
  kmax = k_range(2)
  thmin = th_range(1)
  thmax = th_range(2)
END SUBROUTINE load
!********************************************************************************!

