!********************************************************************************!
!                    DIS-K
!     Main Pogram Dispersion Solver Drifting bi-Kappa
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
PROGRAM main
!
!
!
!
  USE params
  USE solve_params
!
  IMPLICIT NONE
  CHARACTER(100) :: instruction
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: kvec, thvec
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: z0
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: wsol1d, pols1d, wkth1d
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: wsol, pols, wkth
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: ratios1d
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ratios
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xwmaxs
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: kmaxs
  INTEGER, DIMENSION(2) :: igm
  REAL(KIND=8) :: dk, dth
  INTEGER :: i
!
  !$ CALL OMP_SET_NUM_THREADS(2,2)
  !$ CALL OMP_SET_NESTED(.TRUE.) 
!
  !Read the command line argument
  IF(COMMAND_ARGUMENT_COUNT().NE.1) THEN
    WRITE(*,*)'ERROR, INCORRECT NUMBER OF ARGUMENTS'
    STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,sourcefile)
  ! load the parameters
  CALL load
!
  !******************************************
  ! initial seed
  ALLOCATE(z0(nroots))
  DO i=1,nroots
    z0(i) = DCMPLX(w0r(i),w0i(i))
  ENDDO
!
  ! create wavenumber vector
  ALLOCATE(kvec(np))
  IF (logscale_k == 0) THEN
    dk = (kmax-kmin)/REAL(np-1)
    kvec = [( dk*(i-1) + kmin, i=1, np)]
  ELSE
    dk = (DLOG10(kmax)-DLOG10(kmin))/REAL(np-1)
    kvec = [( 10**(dk*(i-1) + DLOG10(kmin)), i=1, np)]
  ENDIF
!
  !********************************************
  IF (nth == 1) THEN ! FIXED ANGLE
    !
    ALLOCATE(wsol1d(np,nroots),pols1d(np,nroots),ratios1d(np,nroots,8))
    !
    PRINT*,'******************************************************'
    PRINT*,'Solving the dispersion relation for theta = ',th0
    !
    ! Fix Initial seed
    !
    DO i=1,nroots
      CALL solve1p(k0(i),th0(i),z0(i))
      PRINT*,'Initial seed',i,z0(i)
    ENDDO
    !
    CALL solve(kvec,k0,th0,z0,np,nroots,wsol1d,pols1d,ratios1d)
    PRINT*,'Done'
    PRINT*,''
    !
    ! write unformatted data
    !
    OPEN(20,file=TRIM(HOME)//'output.save',form='unformatted')
    WRITE(20) np, nth, nsp, nroots, damped
    WRITE(20) rpem, wpewce, wpwc
    WRITE(20) k0, th0
    WRITE(20) rm, ns, rq, betapal, anis, kappa, Ua
    WRITE(20) kvec
    WRITE(20) wsol1d
    WRITE(20) pols1d
    WRITE(20) ratios1d
    WRITE(20) ICHAR(TRIM(subindex))
    CALL FLUSH(20)
    CLOSE(20)
    !
    DEALLOCATE(wsol1d,pols1d,ratios1d)
    PRINT*,'Plotting'
    instruction = 'cd '//TRIM(home)//' && python3 plot_disp.py'
    CALL SYSTEM(instruction)
    !
  ELSE    ! CREATE MAP k vs. th
    !
    PRINT*,'******************************************************'
    PRINT*,'Creating theta vs. k map'
    ALLOCATE(thvec(np),wkth(np,nth,nroots),pols(np,nth,nroots),ratios(np,nth,nroots,8))
    ALLOCATE(xwmaxs(nth,nroots),kmaxs(nth,nroots))
    IF (logscale_th == 0) THEN
      dth = (thmax-thmin)/REAL(nth-1)
      thvec = [( dth*(i-1) + thmin, i=1, nth)]
    ELSE
      dth = (DLOG10(thmax)-DLOG10(thmin))/REAL(nth-1)
      thvec = [( 10**(dth*(i-1) + DLOG10(thmin)), i=1, nth)]
    ENDIF
    !
    ! Fix Initial seed
    !
    DO i=1,nroots
      CALL solve1p(k0(i),th0(i),z0(i))
      PRINT*,'Initial seed',i,z0(i)
    ENDDO
    !
    CALL solve2d(kvec,thvec,z0,wkth,pols,ratios,xwmaxs,kmaxs)
    PRINT*,'Done'
    PRINT*,''
    DO i=1,nroots
      igm = MAXLOC(AIMAG(wkth(:,:,i)))
      PRINT*,'gmax',kvec(igm(1)),thvec(igm(2)),wkth(igm(1),igm(2),i)
    ENDDO
    !
    !
    ! write unformatted data
    !
    OPEN(22,file=TRIM(HOME)//'output_kth.save',form='unformatted')
    WRITE(22) np, nth, nsp, nroots, damped
    WRITE(22) rpem, wpewce, wpwc
    WRITE(22) k0, th0
    WRITE(22) rm, ns, rq, betapal, anis, kappa, Ua
    WRITE(22) kvec
    WRITE(22) thvec
    WRITE(22) wkth
    WRITE(22) pols
    WRITE(22) ratios
    WRITE(22) ICHAR(TRIM(subindex))
    CALL FLUSH(22)
    CLOSE(22)
    !
    DEALLOCATE(thvec,wkth)
    DEALLOCATE(xwmaxs,kmaxs)
    DEALLOCATE(rm,ns,rq,betapal,anis,kappa,Ua,alphapal,alphaper)
    DEALLOCATE(k0,th0,w0r,w0i)
    PRINT*,'Plotting'
    instruction = 'cd '//TRIM(home)//' && python3 plot_map.py'
    CALL SYSTEM(instruction)
  ENDIF
  !
  DEALLOCATE(z0,kvec)
END PROGRAM main
!********************************************************************************!





