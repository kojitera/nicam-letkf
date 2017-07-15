MODULE mod_vbc
  IMPLICIT NONE
  PUBLIC

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE vbc_read(filename, vbc, maxvbc, maxtvsch, ninstrument, &
                    tvsinst, tvsch, ntvsch)
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN)      :: maxvbc, maxtvsch, ninstrument
  INTEGER, INTENT(IN)      :: tvsinst(3,ninstrument)
  INTEGER, INTENT(IN)      :: tvsch(maxtvsch,ninstrument)
  INTEGER, INTENT(IN)      :: ntvsch(ninstrument)
  REAL(8), INTENT(OUT)     :: vbc(maxvbc,maxtvsch,ninstrument)
  
  CHARACTER(len=8):: platname,instname
  INTEGER(4):: isat,ksat,plat,satn,inst,ierr,inum,iflag
  INTEGER(4):: kchan,idxchan,ichan
  INTEGER(4):: iyy,imm,idd,ihr,imn,betatime,nrec
  REAL(8)::  vbc_in(maxvbc)

  vbc = 0.0d0
  PRINT '(2A)','vbc_read reading a file: ',filename
  OPEN(91,file=filename)

  nrec=0
  DO
    READ(91,'(X,A8,I2,2X,A8,I4,8E16.8)',ERR=99,END=99) &
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
!    write(*,*) 'plat, satn, inst', plat, satn, inst
!    write(*,*) 'tvsinst(1:3,1)', tvsinst(:,1)
!    write(*,*) 'ksat', ksat
    IF(ksat>0) THEN
      idxchan=-1
      DO kchan=1,ntvsch(ksat)
        IF(tvsch(kchan,ksat)==ichan) idxchan=kchan
      END DO
      IF(idxchan>0) THEN
        vbc(:,idxchan,ksat)=vbc_in(:)
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
  CLOSE(91)
  PRINT '(A,I4)','  varbc records read: ',nrec
  RETURN
END SUBROUTINE vbc_read
!------------------------------------------------------------------------------
SUBROUTINE vbc_write(filename, vbc, maxvbc, maxtvsch, &
                     tvsinst, tvsch, ntvsch, ntvschan)
  USE mod_adm
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN)      :: maxvbc, maxtvsch
  INTEGER, INTENT(IN)      :: tvsinst(3,1)
  INTEGER, INTENT(IN)      :: ntvsch(1)
  INTEGER, INTENT(IN)      :: tvsch(ntvsch(1),1)
  REAL(8), INTENT(IN)     :: vbc(maxvbc,maxtvsch,1)
  INTEGER, INTENT(IN)      :: ntvschan(ntvsch(1),1)

  CHARACTER(len=8):: platname,instname
  INTEGER(4):: isat,ksat,plat,satn,inst,ierr,inum,iflag
  INTEGER(4):: kchan,idxchan,ichan
  INTEGER(4):: iyy,imm,idd,ihr,imn,betatime,nrec
  REAL(8)::  vbc_in(maxvbc)

  OPEN(92,file=filename)
  nrec = 0
  DO ksat=1,1
    DO idxchan=1,ntvsch(ksat)
      iflag=0
      DO inum=1,maxvbc
        IF(vbc(inum,idxchan,ksat).NE.0.d0) iflag=1
      END DO
      WRITE(ADM_LOG_FID,*) idxchan, tvsch(idxchan, 1), iflag
      IF(iflag==1) THEN
        WRITE(92,'(X,A8,I2,2X,A8,I4,8E16.8,I5)') &
   &      platform_name(tvsinst(1,ksat)), tvsinst(2,ksat), &
   &      inst_name(tvsinst(3,ksat)), tvsch(idxchan,ksat), &
   &      vbc(1:8,idxchan,ksat), ntvschan(idxchan,ksat)
        nrec = nrec+1
      END IF
    END DO
  END DO
  CLOSE(92)
  PRINT '(A,I4)','  varbc_coef records written: ',nrec
  RETURN

END SUBROUTINE vbc_write
!------------------------------------------------------------------------------
SUBROUTINE das_vbc( maxtvsprof, maxvbc, maxtvsch, &
                   tvsname, tvsch, ntvsch,   &
                   tvsdat, hx, pred, vbcf, vbca, qc, err, ntvschan )
  USE mod_adm
  USE common_mtx
  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: maxtvsprof, maxvbc, maxtvsch
  CHARACTER(4), INTENT(IN) :: tvsname
  INTEGER, INTENT(IN)      :: tvsch(maxtvsch)
  INTEGER, INTENT(IN)      :: ntvsch
  REAL(4), INTENT(IN)      :: tvsdat(maxtvsch,maxtvsprof)
  REAL(8), INTENT(IN)      :: hx(maxtvsch,maxtvsprof)
  REAL(8), INTENT(IN)      :: pred(maxvbc,maxtvsch,maxtvsprof)
  REAL(8), INTENT(INOUT)   :: vbcf(maxvbc,maxtvsch)
  REAL(8), INTENT(OUT)     :: vbca(maxvbc,maxtvsch)
  INTEGER, INTENT(IN)      :: qc(maxtvsch,maxtvsprof)
  REAL(4), INTENT(IN)      :: err(maxtvsch,maxtvsprof)
  INTEGER, INTENT(OUT) :: ntvschan(ntvsch)

  REAL(8) :: a(maxvbc,maxvbc)
  REAL(8) :: b(maxvbc)
  REAL(8) :: ainv(maxvbc,maxvbc)
  INTEGER:: i,j,ic,n,islot,nn
  LOGICAL, PARAMETER :: msw_vbc = .TRUE.

  !WRITE(ADM_LOG_FID,*) 'pred'
  !DO n = 1, maxtvsprof
  !DO ic = 1, maxtvsch
  !  WRITE(ADM_LOG_FID,'(2i6,8f10.6)') n, ic, (pred(i,ic,n),i=1,maxvbc)
  !END DO
  !END DO
  !WRITE(ADM_LOG_FID,*)

  vbca = 0.0d0
  DO ic=1,ntvsch
    ntvschan(ic) = SUM(qc(ic,:))
    IF(msw_vbc .AND. ntvschan(ic) /= 0 ) THEN
      WRITE(ADM_LOG_FID,'(3A,I3,A,I6)') &
            ' >> VBC executed for instrument,channel,ntvsl:',&
            & tvsname,',',tvsch(ic),',',ntvschan(ic)
      CALL vbc_local(ic,maxvbc,maxtvsprof,ntvsch,ntvschan(ic),&
                     tvsdat,hx,pred,qc,err,a,b)
      !WRITE(ADM_LOG_FID,*) 'a'
      !DO j = 1, maxvbc
      !  WRITE(ADM_LOG_FID,'(8f12.1)') (a(i,j),i=1,maxvbc)
      !END DO
      CALL mtx_inv(maxvbc,a,ainv)
      !WRITE(ADM_LOG_FID,*) 'ainv'
      !DO j = 1, maxvbc
      !  WRITE(ADM_LOG_FID,'(8E12.5)') (ainv(i,j),i=1,maxvbc)
      !END DO
      !WRITE(ADM_LOG_FID,*) 'b'
      !WRITE(ADM_LOG_FID,'(8f12.4)') (b(i),i=1,maxvbc)
    
      vbca(:,ic) = vbcf(:,ic)
      DO n=1,maxvbc
        vbca(:,ic) = vbca(:,ic) - ainv(:,n)*b(n) !ATTN: sign for beta
      END DO
    ELSE
      WRITE(ADM_LOG_FID,'(3A,I3,A,I6)') &
            ' !! NO VBC executed for instrument,channel,ntvsl:',&
            & tvsname,',',tvsch(ic),',',ntvschan(ic)
      vbca(:,ic) = vbcf(:,ic)
    END IF
  END DO
  vbcf = vbca

END SUBROUTINE das_vbc
!-----------------------------------------------------------------------
!  (out) a = B_beta^-1 + p R^-1 p^T
!  (out) b = p R^-1 d
!-----------------------------------------------------------------------
SUBROUTINE vbc_local(ic, maxvbc, maxtvsprof, maxtvsch, ntvsl, &
                     tvsdat, hx, pred, qc, err, a, b)
  USE mod_adm
  IMPLICIT NONE
  INTEGER,PARAMETER   :: msw=1
  INTEGER,PARAMETER   :: nmin=400
  INTEGER,INTENT(IN)  :: ic,maxvbc,maxtvsch,maxtvsprof,ntvsl
  REAL(4), INTENT(IN) :: tvsdat(maxtvsch,maxtvsprof)
  REAL(8),INTENT(IN)  :: hx(maxtvsch,maxtvsprof)
  REAL(8),INTENT(IN)  :: pred(maxvbc,maxtvsch,maxtvsprof)
  INTEGER, INTENT(IN) :: qc(maxtvsch,maxtvsprof)
  REAL(4), INTENT(IN) :: err(maxtvsch,maxtvsprof)
  REAL(8),INTENT(OUT) :: a(maxvbc,maxvbc)
  REAL(8),INTENT(OUT) :: b(maxvbc)
  REAL(8) :: dep,dep0
  REAL(8) :: bias,bias0
  REAL(8) :: r,tmp
  INTEGER:: islot, iprof, i,j,n

  a = 0.0d0
  b = 0.0d0
  dep = 0.0d0
  dep0 = 0.0d0
  bias = 0.0d0
  bias0 = 0.0d0
  n = 0
  DO iprof=1,maxtvsprof
    IF(qc(ic,iprof)/=1) CYCLE
    r = err(ic,iprof)**2
    DO j=1,maxvbc
      DO i=1,maxvbc
        a(i,j) = a(i,j) &
             & + pred(i,ic,iprof) &
             & * pred(j,ic,iprof) / r
      END DO
    END DO
    IF(msw == 1) THEN ! Y.Sato
      IF(ntvsl < nmin) THEN
        tmp = REAL(nmin,8) / r
      ELSE
        tmp = (REAL(ntvsl,8) &
          & / (LOG10(REAL(ntvsl,8)/REAL(nmin,kind=8))+1.0d0)) / r
      END IF
    ELSE IF(msw == 2) THEN ! D.Dee
      tmp = REAL(ntvsl,8) / r
    ELSE ! Constant
      tmp = 100.0d0
    END IF
    DO i=1,maxvbc
      a(i,i) = a(i,i) + tmp
    END DO
    b(:) = b(:) + pred(:,ic,iprof) / r &
              & *(tvsdat(ic,iprof)-hx(ic,iprof))
    write(*,*) ic, iprof, tvsdat(ic,iprof), hx(ic,iprof), r
    WRITE(ADM_LOG_FID,'(i6,8f12.1)') iprof, (a(i,i),i=1,maxvbc)
    bias = bias+tvsdat(ic,iprof)-hx(ic,iprof)
    dep = dep+(tvsdat(ic,iprof)-hx(ic,iprof))**2
    n = n+1
  END DO

  dep = SQRT(dep / REAL(n,kind=8))
  dep0 = SQRT(dep0 / REAL(n,kind=8))
  bias = bias / REAL(n,kind=8)
  bias0 = bias0 / REAL(n,kind=8)
  !PRINT '(2A,I3,4F12.4)',' >> D monit: ',tvsname(iinst),tvsch(ic,iinst),bias0,bias,dep0,dep

  RETURN
END SUBROUTINE vbc_local


END MODULE mod_vbc
