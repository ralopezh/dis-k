!********************************************************************************!
!                    DIS-K
!     Program: Nyquist diagram
!     Author: Rodrigo A. Lopez
!             Universidad de Santiago de Chile, Usach
!             rlopez186@gmail.com
!             ORCID: 0000-0003-3223-1498
!
!********************************************************************************!
PROGRAM nyquist
!
!
  USE params
  USE integral_params
  IMPLICIT NONE
  CHARACTER(100) :: instruction
  REAL(KIND=8) :: wrmin, wrmax, wimin, wimax
  REAL(KIND=8) :: dwr, dwi
  INTEGER :: i, j, nex
  COMPLEX(KIND=8) :: disp
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: wrvec, wivec
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: z
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: zr, zi
!
  !$ CALL OMP_SET_NUM_THREADS(8)
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
! Run nyquist to set the initial seed
!
  wrmin = wr_range(1)
  wrmax = wr_range(2)
  wimin = wi_range(1)
  wimax = wi_range(2)
!
  PRINT*,'Calculating Nyquist Contour'
! Vectors wr and wi
  ALLOCATE(wrvec(nw),wivec(nw),z(nw,nw,nroots),zr(nw,nw,nroots),zi(nw,nw,nroots))
  dwr = (wrmax-wrmin)/REAL(nw-1)
  dwi = (wimax-wimin)/REAL(nw-1)
  wrvec = [( dwr*(i-1), i=1, nw)] + wrmin
  wivec = [( dwi*(i-1), i=1, nw)] + wimin
!
! EVALUATION z = DISP(wr+i*wi, k0)
  DO nex=1,nroots
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
!$OMP DO
    DO j=1,nw
      WRITE(*,'(F5.1,A2,A1)',ADVANCE='NO') j*100.0/REAL(nw),'%',ACHAR(13)
      DO i=1,nw
        z(i,j,nex) = disp(DCMPLX(wrvec(i),wivec(j)),k0(nex),th0(nex))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO
!
  zr = REAL(z)
  zi = AIMAG(z)
  !
  ! write unformatted data
  !
  OPEN(23,file=TRIM(HOME)//'diagram.save',form='unformatted')
  WRITE(23) nroots
  WRITE(23) nw
  WRITE(23) k0, th0
  WRITE(23) wrvec, wivec
  WRITE(23) zr
  WRITE(23) zi
  WRITE(23) ICHAR(TRIM(subindex))
  CALL FLUSH(23)
  CLOSE(23)
  !
  DEALLOCATE(wrvec,wivec,z,zr,zi)
  PRINT*,'Nyquist DONE'
  instruction = 'cd '//TRIM(home)//' && python3 nyquist.py'
  CALL SYSTEM(instruction)
END PROGRAM nyquist
!********************************************************************************!
