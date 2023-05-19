!********************************************************************************!
!                    DIS-K
!     Pogram solve the dispersion relation in the theta vs. k map
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
SUBROUTINE solve2d(kvec,thvec,z0,wsol,pols2d,ratios2d,xwmaxs,kmaxs)
!
! 
!
  USE params
  USE solve_params
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN), DIMENSION(np) :: kvec
  REAL(KIND=8), INTENT(IN), DIMENSION(nth) :: thvec
  COMPLEX(KIND=8), INTENT(IN), DIMENSION(nroots) :: z0
  COMPLEX(KIND=8), INTENT(OUT), DIMENSION(np,nth,nroots) :: wsol, pols2d
  REAL(KIND=8), INTENT(OUT), DIMENSION(np,nth,nroots,8) :: ratios2d
  COMPLEX(KIND=8), DIMENSION(np,1) :: xw, poli
  REAL(KIND=8), DIMENSION(np,1,8) :: ratsi
  COMPLEX(KIND=8), DIMENSION(1) :: zc
  REAL(KIND=8) :: kc, thc, gammamax, local1D
  INTEGER ::  nex, ic, ith, ith0, igm, ncount,l
  COMPLEX(KIND=8), INTENT(OUT), DIMENSION(nth,nroots) :: xwmaxs
  REAL(KIND=8), INTENT(OUT), DIMENSION(nth,nroots) :: kmaxs
!
  wsol = 0.0
  pols2d = 0.0
  ratios2d = 0.0
  xwmaxs = 0.0
  kmaxs = 0.0
  !
  DO nex=1,nroots
    PRINT*,'Solving Root No. ',nex
    ! Locate k0 and th0 in the arrays kvec and thvec
    ic = MINLOC(ABS(kvec-k0(nex)),DIM=1)
    ith0 = MINLOC(ABS(thvec-th0(nex)),DIM=1)
    igm = ic
    !
!!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(thc,zc,kc,igm,ith,xw,poli,ratsi,gammamax,ncount,l)
!!$OMP SECTION
  !
    zc = z0(nex)
    kc = k0(nex)
    !
    IF(ISNAN(AIMAG(zc(1)))) THEN
      !zc = 0.0*zc
      EXIT
    ENDIF
    !
    ! Start the Loop from
    ! theta = th0 ----> theta = thmax
    !
    ncount = 0
    DO ith=ith0,nth
      ncount = ncount+1
      FLUSH(0)
      WRITE(0,'(F5.1,A2,A1)',ADVANCE='NO') ncount*100.0/REAL(nth),'%',ACHAR(13)
      FLUSH(0)
      !
      thc = thvec(ith)
      CALL solve(kvec,kc,thc,zc,np,1,xw,poli,ratsi)
      wsol(:,ith,nex) = xw(:,1)
      pols2d(:,ith,nex) = poli(:,1)
      ratios2d(:,ith,nex,:) = ratsi(:,1,:)
      !
      !Next seed
      IF (damped == 0 .AND. prevmax == 0) THEN
        !for unstable solution
        !igm = MAXLOC(AIMAG(xw(:,1)),DIM=1)
        gammamax = local1D(AIMAG(xw(:,1)),np,AIMAG(xw(igm,1)),igm)
        !print*,igm,AIMAG(xw(igm,1))
        !print*,AIMAG(xw(:,1))
        !PRINT*,''
        !PRINT*,thc,gammamax,igm,MAXVAL(AIMAG(xw(:,1)))
        !pause
        kc = kvec(igm)
        zc = xw(igm,1)
        xwmaxs(ith,nex) = xw(igm,1)
        kmaxs(ith,nex) = kvec(igm)
        IF (AIMAG(zc(1)) < gammatol .OR. ISNAN(AIMAG(zc(1))) ) THEN
          EXIT
        ENDIF
      ELSE
        kc = kvec(ic)
        zc = xw(ic,1)
        IF (damped == 0 .AND. (AIMAG(zc(1)) < gammatol .OR. ISNAN(AIMAG(zc(1)))) ) THEN
          EXIT
        ENDIF        
      ENDIF
      !print*,ith,thc,kc,AIMAG(zc)
      CALL save_temp(kvec,thvec,wsol)
    ENDDO
!!$OMP SECTION
    !
    ! Now from theta = th0-1 ----> theta = thmin
    !
    zc = z0(nex)
    kc = k0(nex)
    igm = ic
    DO ith=ith0-1,1,-1
      ncount = ncount+1
      FLUSH(0)
      WRITE(0,'(F5.1,A2,A1)',ADVANCE='NO') ncount*100.0/REAL(nth),'%',ACHAR(13)
      FLUSH(0)
      !
      thc = thvec(ith)
      CALL solve(kvec,kc,thc,zc,np,1,xw,poli,ratsi)
      wsol(:,ith,nex) = xw(:,1)
      pols2d(:,ith,nex) = poli(:,1)
      ratios2d(:,ith,nex,:) = ratsi(:,1,:)
      !
      !Next seed
      IF (damped == 0 .AND. prevmax == 0) THEN
        !for unstable solution
        !igm = MAXLOC(AIMAG(xw(:,1)),DIM=1)
        gammamax = local1D(AIMAG(xw(:,1)),np,AIMAG(xw(igm,1)),igm)
        !PRINT*,thc,gammamax,igm
        kc = kvec(igm)
        zc = xw(igm,1)
        xwmaxs(ith,nex) = xw(igm,1)
        kmaxs(ith,nex) = kvec(igm)
        IF (AIMAG(zc(1)) < gammatol .OR. ISNAN(AIMAG(zc(1))) ) THEN
          EXIT
        ENDIF
      ELSE
        kc = kvec(ic)
        zc = xw(ic,1)
        IF (damped == 0 .AND. (AIMAG(zc(1)) < gammatol .OR. ISNAN(AIMAG(zc(1)))) ) THEN
          EXIT
        ENDIF        
      ENDIF
      !print*,ith,thc,kc,AIMAG(zc)
    ENDDO
!!$OMP END PARALLEL SECTIONS
    CALL save_temp(kvec,thvec,wsol)
   ENDDO
!
END SUBROUTINE solve2d
!********************************************************************************!
!********************************************************************************!
FUNCTION local1D(array,n,seed,pos)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND=8), DIMENSION(n), INTENT(IN) :: array
  REAL(KIND=8), INTENT(IN) :: seed
  INTEGER, INTENT(INOUT) :: pos
  REAL(KIND=8) :: local1D, df
  INTEGER :: indx, i, dir
  !
  indx = pos!MINLOC(ABS(array-seed),DIM=1)!position of seed in array
  !print*,indx
  !
  IF(indx == 1) THEN
    indx = indx+1
  ELSE IF(indx == n) THEN
    indx = indx-1
  ENDIF
  DO WHILE((indx /= 1) .OR. (indx /= n) .OR. array(indx) <= 0.0)
    IF((array(indx-1) < array(indx)) .AND. (array(indx) > array(indx+1))) THEN
      pos = indx
      local1D = array(indx)
      EXIT
    ELSE IF((array(indx-1) < array(indx)) .AND. array(indx) < array(indx+1)) THEN
      indx = indx+1
    ELSE IF((array(indx-1) > array(indx)) .AND. array(indx) > array(indx+1)) THEN
      indx = indx-1
    ELSE
      indx = MAXLOC(array,DIM=1)
      EXIT
    ENDIF
  ENDDO
!
  !pause
  pos = indx
  local1D = array(indx)
END FUNCTION local1D
!********************************************************************************!
!********************************************************************************!
SUBROUTINE save_temp(kvec,thvec,wsol)
  USE params
  USE solve_params
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN), DIMENSION(np) :: kvec
  REAL(KIND=8), INTENT(IN), DIMENSION(nth) :: thvec
  COMPLEX(KIND=8), INTENT(IN), DIMENSION(np,nth,nroots) :: wsol

  OPEN(28,file=TRIM(HOME)//'temp.save',form='unformatted')
  WRITE(28) np,nth,nroots
  WRITE(28) kvec
  WRITE(28) thvec
  WRITE(28) wsol
  CLOSE(28)
!
END SUBROUTINE save_temp
