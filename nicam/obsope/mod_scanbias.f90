MODULE mod_scanbias
  IMPLICIT NONE
  PUBLIC

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE vbc_scan_read(filename, vbc_scan, maxfoot, maxtvsch, ninstrument, &
                         tvsinst, tvsch, ntvsch)
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN)      :: maxfoot, maxtvsch, ninstrument
  INTEGER, INTENT(IN)      :: tvsinst(3,ninstrument)
  INTEGER, INTENT(IN)      :: tvsch(maxtvsch,ninstrument)
  INTEGER, INTENT(IN)      :: ntvsch(ninstrument)
  REAL(8), INTENT(OUT)     :: vbc_scan(maxfoot,maxtvsch,ninstrument)
  
  CHARACTER(len=8):: platname,instname
  INTEGER(4):: isat,ksat,plat,satn,inst,ierr,inum,iflag
  INTEGER(4):: kchan,idxchan,ichan
  INTEGER(4):: iyy,imm,idd,ihr,imn,betatime,nrec
  REAL(8)::  vbc_in(maxfoot)

  vbc_scan = 0.0d0
  PRINT '(2A)','vbc_read reading a file: ',filename
  OPEN(99,file=filename)

  nrec=0
  DO
    READ(99,'(X,A8,I2,2X,A8,I4,120E16.8)',ERR=99,END=99) &
        & platname,satn,instname,ichan,vbc_in
    write(*,*) platname,satn,instname,ichan
    plat=0
    DO plat=1,nplatforms
      IF(platname==platform_name(plat)) EXIT
    END DO
!    IF(instname=='imager  ') THEN
!      SELECT CASE (instname)
!        CASE('goes    '); inst=rttv_inst_goesi
!        CASE('mtsat   '); inst=rttv_inst_mtsati
!      END SELECT
!    ELSE
      DO inst=0,ninst-1
        IF(instname==inst_name(inst)) EXIT
      END DO
!    END IF

    ksat=-1
    DO isat=1,ninstrument
      IF(tvsinst(1,isat)==plat.AND. &
      &  tvsinst(2,isat)==satn.AND. &
      &  tvsinst(3,isat)==inst) THEN
        ksat=isat
      END IF
    END DO
    write(*,*) 'plat, satn, inst', plat, satn, inst
    write(*,*) 'tvsinst(1:3,1)', tvsinst(:,1)
    write(*,*) 'ksat', ksat
    IF(ksat>0) THEN
      idxchan=-1
      DO kchan=1,ntvsch(ksat)
        IF(tvsch(kchan,ksat)==ichan) idxchan=kchan
      END DO
      IF(idxchan>0) THEN
        vbc_scan(:,idxchan,ksat)=vbc_in(:)
        nrec = nrec+1
      ELSE
        PRINT '(A,A8,I2,2X,A8,I4)','  !!! record skipped.. ',&
          & platname,satn,instname,ichan
      END IF
    ELSE
      PRINT '(A,A8,I2,2X,A8,I4)','  !!! record skipped.. ',&
        & platname,satn,instname,ichan
    END IF
  END DO
  99 CONTINUE
  CLOSE(99)
  PRINT '(A,I4)','  varbc_scan records read: ',nrec
  RETURN
END SUBROUTINE vbc_scan_read
!------------------------------------------------------------------------------
SUBROUTINE vbc_scan_write(filename, vbc_scan, maxfoot, maxtvsch, &
                     tvsinst, tvsch, ntvsch)
  USE mod_adm
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN)      :: maxfoot,maxtvsch
  INTEGER, INTENT(IN)      :: tvsinst(3,1)
  INTEGER, INTENT(IN)      :: ntvsch(1)
  INTEGER, INTENT(IN)      :: tvsch(ntvsch(1),1)
  REAL(8), INTENT(IN)      :: vbc_scan(maxfoot,maxtvsch,1)

  CHARACTER(len=8):: platname,instname
  INTEGER(4):: isat,ksat,plat,satn,inst,ierr,inum,iflag
  INTEGER(4):: kchan,idxchan,ichan
  INTEGER(4):: iyy,imm,idd,ihr,imn,betatime,nrec

  OPEN(99,FILE=filename)
  nrec=0
  ksat=1
  DO idxchan=1,ntvsch(ksat)
    iflag=0
    DO inum=1,maxfoot
      IF(vbc_scan(inum,idxchan,ksat).NE.0.d0) iflag=1
    END DO
    IF(iflag==1) THEN
      WRITE(99,'(X,A8,I2,2X,A8,I4,120E16.8)') &
 &      platform_name(tvsinst(1,ksat)), tvsinst(2,ksat), &
 &      inst_name(tvsinst(3,ksat)), tvsch(idxchan,ksat), &
 &      vbc_scan(1:maxfoot,idxchan,ksat)
      nrec = nrec+1
    END IF
  END DO
  CLOSE(99)

END SUBROUTINE vbc_scan_write
!------------------------------------------------------------------------------
SUBROUTINE tvs_ominusb_output(nn, imem, maxtvsch, nobs, ominusb, ifov, lsql, &
                              lwp, fname, tvsname, islot)
  USE mod_adm
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: nn
  INTEGER, INTENT(IN)        :: imem
  INTEGER, INTENT(IN)        :: maxtvsch
  INTEGER, INTENT(IN)        :: nobs
  INTEGER, INTENT(IN)        :: ifov(nobs)
  INTEGER, INTENT(IN)        :: lsql(nobs)
  REAL(4), INTENT(INOUT)     :: lwp(nobs)
  REAL(8), INTENT(IN)        :: ominusb(maxtvsch, nobs)
  CHARACTER(*), INTENT(IN)   :: fname
  CHARACTER(4), INTENT(IN)   :: tvsname
  INTEGER, INTENT(IN)        :: islot

  CHARACTER(256) :: filename
  CHARACTER(2)   :: cslot
  CHARACTER(6)   :: cmem
  INTEGER        :: i, ifoot
    
  WRITE(cslot,'(i2.2)') islot
  WRITE(cmem,'(i6.6)') imem
  WRITE(ADM_LOG_FID,*) 'TRIM(fname)'
  WRITE(ADM_LOG_FID,*) TRIM(fname)
  WRITE(ADM_LOG_FID,*) 'TRIM(tvsname)'
  WRITE(ADM_LOG_FID,*) TRIM(tvsname)
  WRITE(ADM_LOG_FID,*) 'TRIM(cslot)'
  WRITE(ADM_LOG_FID,*) TRIM(cslot)
  WRITE(ADM_LOG_FID,*) 'TRIM(cmem)'
  WRITE(ADM_LOG_FID,*) TRIM(cmem)
  WRITE(ADM_LOG_FID,*) 'nobs =', nobs 
  FLUSH(ADM_LOG_FID)
  filename=TRIM(fname)//'_'//TRIM(tvsname)//TRIM(cslot)//TRIM(cmem)//'.txt'
  OPEN(201,FILE=TRIM(filename))
  DO i = 1, nobs
    ifoot=ifov(i)
    IF( lsql(i)==1 .AND. lwp(i) < 0.02 ) THEN
      WRITE(201,'(A,i4,10F12.5)') tvsname, ifoot, lwp(i), REAL(lsql(i)), ominusb(1:maxtvsch,i)
    END IF 
  END DO
  CLOSE(201)

END SUBROUTINE tvs_ominusb_output
!------------------------------------------------------------------------------
END MODULE mod_scanbias
