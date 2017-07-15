MODULE common_tvs_nicam
!=======================================================================
!
! [PURPOSE:] ATOVS observational procedures for NICAM
!
! [HISTORY:]
!   01/22/2015 Koji TERASAKI  created for NICAM
!
!=======================================================================
  USE common
  USE common_nicam
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: ninstrument=8
  INTEGER,PARAMETER :: maxtvsch=4
  INTEGER,PARAMETER :: maxvbc=8
  INTEGER,SAVE :: maxtvsprof
  INTEGER,SAVE :: maxtvsfoot
  INTEGER,SAVE :: ntvs
  CHARACTER(4),SAVE :: tvsname(ninstrument)
  INTEGER,SAVE :: tvsinst(3,ninstrument)
  INTEGER,SAVE :: tvsch(maxtvsch,ninstrument)
  INTEGER,SAVE :: ntvsch(ninstrument)
  INTEGER,SAVE :: nfootp(ninstrument)
  !INTEGER,SAVE :: ntvschan(maxtvsch,ninstrument)   ! 2016.07.06 Koji
  INTEGER,SAVE :: ntvsprof(ninstrument)
  INTEGER,ALLOCATABLE,SAVE :: ntvsprofslots(:,:)

!  REAL(r_size),ALLOCATABLE,SAVE :: tvslon     (:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvslat     (:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvszenith  (:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvsskin    (:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvsstmp    (:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvsclw     (:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvsemis  (:,:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvsdat   (:,:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvserr   (:,:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvsdep   (:,:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvshdxf(:,:,:,:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: tvswgt (:,:,:,:,:)
!  INTEGER,ALLOCATABLE,SAVE :: tvsqc(:,:,:,:)

  integer,parameter:: &
     &    rttv_plat_noaa =1 ,rttv_plat_dmsp =2 &
     &   ,rttv_plat_meteo=3 ,rttv_plat_goes =4 &
     &   ,rttv_plat_gms  =5 ,rttv_plat_fy2  =6 &
     &   ,rttv_plat_trmm =7 ,rttv_plat_ers  =8 &
     &   ,rttv_plat_eos  =9 ,rttv_plat_metop=10 &
     &   ,rttv_plat_envi =11,rttv_plat_msg  =12 &
     &   ,rttv_plat_fy1  =13,rttv_plat_adeos=14 &
     &   ,rttv_plat_mtsat=15,rttv_plat_cori =16 &
     &   ,rttv_plat_jpss =17,rttv_plat_gifts=18 &
     &   ,rttv_plat_metop2=19

  integer,parameter:: &
     &    rttv_inst_hirs  =0, rttv_chan_hirs  =19   & ! HIRS
     &   ,rttv_inst_amsua =3, rttv_chan_amsua =15   & ! AMSU-A
     &   ,rttv_inst_amsub =4, rttv_chan_amsub =5    & ! AMSU-B
     &   ,rttv_inst_ssmi  =6, rttv_chan_ssmi  =7    & ! SSMI
     &   ,rttv_inst_tmi   =9, rttv_chan_tmi   =9    & ! TMI
     &   ,rttv_inst_ssmis =10,rttv_chan_ssmis =24   & ! SSMIS
     &   ,rttv_inst_airs  =11,rttv_chan_airs  =616  & ! AIRS
     &   ,rttv_inst_iasi  =16,rttv_chan_iasi  =8461 & ! IASI
     &   ,rttv_inst_amsr  =17,rttv_chan_amsr  =14   & ! AMSR
     &   ,rttv_inst_atms  =19,rttv_chan_atms  =22   & ! ATMS
     &   ,rttv_inst_mviri =20,rttv_chan_mviri =2    & ! METEOSAT
     &   ,rttv_inst_seviri=21,rttv_chan_seviri=8    & ! MSG
     &   ,rttv_inst_goesi =22,rttv_chan_goesi =4    & ! GOES-IMAGER(IR)
     &   ,rttv_inst_mtsati=24,rttv_chan_mtsati=4      ! MTSAT

CONTAINS
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_ntvs_mpi(cfile,dir)
  IMPLICIT NONE
  CHARACTER(*),INTENT(INOUT) :: cfile
  CHARACTER(*),optional :: dir
  CHARACTER(256) :: filename
  REAL(r_sngl) :: wk(10000)
  INTEGER :: ios
  INTEGER :: n,nrec
  INTEGER :: iunit=92
  LOGICAL :: ex

  if(present(dir)) then

  else
    dir='./'
  end if

  ntvsprof = 0
  DO n=1,ninstrument
    cfile(1:4)=tvsname(n)
    filename=trim(dir)//trim(cfile)
    WRITE(ADM_LOG_FID,'(A,I7.7,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename)
    !WRITE(cfile(1:4),'(A4)') tvsname(n)
    !INQUIRE(FILE=cfile,EXIST=ex)
    INQUIRE(FILE=filename,EXIST=ex)
    IF(ex) THEN
      !OPEN(iunit,file=cfile,form='unformatted',access='sequential')
      OPEN(iunit,file=filename,form='unformatted',access='sequential')
      nrec=8+5*ntvsch(n)
      !nrec=6+5*ntvsch(n)+(nlev-2)*ntvsch(n)
      !nrec=6+5*ntvsch(n)
      DO
        READ(iunit,IOSTAT=ios) wk(1:nrec)
        IF(ios /= 0) EXIT
        ntvsprof(n) = ntvsprof(n) + 1
      END DO
      CLOSE(iunit)
    ELSE
      !write(ADM_LOG_FID,*) cfile,' does not exist -- skipped'
      write(ADM_LOG_FID,*) filename,' does not exist -- skipped'
      !PRINT '(2A)',cfile,' does not exist -- skipped'
      CYCLE
    END IF
  END DO

  WRITE(ADM_LOG_FID,'(I10,A,I6.6)') SUM(ntvsprof(:)),' ATOVS OBSERVATION RECORDS INPUT in MYRANK',myrank 
  WRITE(ADM_LOG_FID,'(8A11)') 'N15/AMSUA','N16/AMSUA','N18/AMSUA','N19/AMSUA','M02/AMSUA', 'M01/IASI', 'M02/IASI', 'J00/ATMS'
  WRITE(ADM_LOG_FID,'(8I11)') (ntvsprof(n),n=1,8)
  !WRITE(ADM_LOG_FID,'(3A12)') 'M01/IASI', 'M02/IASI', 'J00/ATMS'
  !WRITE(ADM_LOG_FID,'(3I12)') (ntvsprof(n),n=6,8)
  WRITE(ADM_LOG_FID,'(A)') '=========================================================='
  FLUSH(ADM_LOG_FID)

  RETURN
END SUBROUTINE get_ntvs_mpi

SUBROUTINE read_tvs_mpi(cfile,elem,rlon,rlat,zenith,skin,stmp,clw,elev,odat,oerr,hdxf,oqc,foot,dir)
  IMPLICIT NONE
  CHARACTER(*),INTENT(INOUT) :: cfile
  CHARACTER(*),optional :: dir
  CHARACTER(256) :: filename
  REAL(r_size),INTENT(OUT) :: elem(maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: rlon(maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: rlat(maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: zenith(maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: skin(maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: stmp(maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: clw(maxtvsprof,ninstrument)
  !REAL(r_size),INTENT(OUT) :: emis(maxtvsch,maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: elev(maxtvsch,maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: odat(maxtvsch,maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: oerr(maxtvsch,maxtvsprof,ninstrument)
  REAL(r_size),INTENT(OUT) :: hdxf(maxtvsch,maxtvsprof,ninstrument)
  INTEGER,INTENT(OUT) ::  oqc(maxtvsch,maxtvsprof,ninstrument)
  INTEGER,INTENT(OUT) ::  foot(maxtvsprof,ninstrument)
  !REAL(r_size),INTENT(OUT) :: wgt(nlev,maxtvsch,maxtvsprof,ninstrument)

  !REAL(r_size),INTENT(OUT) ::  oqc(maxtvsch,maxtvsprof,ninstrument)
  REAL(r_sngl) :: wk(1000)
  INTEGER :: n,nn,nrec,i,irec,k
  INTEGER :: iunit=93

  if(present(dir)) then

  else
    dir='./'
  end if

  DO nn=1,ninstrument
    IF(ntvsprof(nn) == 0) CYCLE
    nrec=8+5*ntvsch(nn)
    !nrec=6+5*ntvsch(nn)+(nlev-2)*ntvsch(nn)
    WRITE(cfile(1:4),'(A4)') tvsname(nn)
    filename=trim(dir)//trim(cfile)
    WRITE(ADM_LOG_FID,'(A,I7.7,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename)
    !WRITE(ADM_LOG_FID,'(A,I7)') 'NUMBER OF RECORDS =  ', nrec
    OPEN(iunit,file=filename,form='unformatted',access='sequential')
    n = 1
    DO
      READ(iunit) wk(1:nrec)
      elem(n,nn)   = REAL(wk(1),r_size)
      skin(n,nn)   = REAL(wk(2),r_size)
      rlon(n,nn)   = REAL(wk(3),r_size)
      rlat(n,nn)   = REAL(wk(4),r_size)
      zenith(n,nn) = REAL(wk(5),r_size)
      stmp(n,nn)   = REAL(wk(6),r_size)
      clw(n,nn)    = REAL(wk(7),r_size)
      foot(n,nn)   =  INT(wk(8),r_size)
      irec = 9
      DO i=1,ntvsch(nn)
        elev(i,n,nn) = REAL(wk(irec),r_size) * 100.0d0
        irec = irec + 1
      END DO
      DO i=1,ntvsch(nn)
        odat(i,n,nn) = REAL(wk(irec),r_size)
        irec = irec + 1
      END DO
      DO i=1,ntvsch(nn)
        oerr(i,n,nn) = REAL(wk(irec),r_size)
        irec = irec + 1
      END DO
      DO i=1,ntvsch(nn)
        hdxf(i,n,nn) = REAL(wk(irec),r_size)
        irec = irec + 1
      END DO
      DO i=1,ntvsch(nn)
        oqc(i,n,nn) = INT(wk(irec))
        irec = irec + 1
      END DO
      !DO i=1,ntvsch(nn)
      !  DO k = 2, nlev-1
      !    wgt(k,i,n,nn) = REAL(wk(irec))
      !    irec = irec + 1
      !  END DO
      !  wgt(1,i,n,nn)=wgt(2,i,n,nn)
      !  wgt(40,i,n,nn)=wgt(39,i,n,nn)
      !END DO

      !IF(skin(n,nn)/=1) THEN
      !  oqc(:,n,nn)=0
      !END IF

      n = n + 1
      IF(n > ntvsprof(nn)) EXIT
    END DO
    CLOSE(iunit)
  END DO

  RETURN
END SUBROUTINE read_tvs_mpi

!-----------------------------------------------------------------------
! LIST ASSIMILATED INSTRUMENT CHANNELS
!-----------------------------------------------------------------------
SUBROUTINE set_instrument
  integer :: ic
  tvsch = 0
  !
  ! NOAA-15 AMSU-A
  !
  tvsname(1) = 'AA15'
  tvsinst(1,1) = rttv_plat_noaa
  tvsinst(2,1) = 15
  tvsinst(3,1) = rttv_inst_amsua
  !tvsch( 1,1)  =  5
  tvsch( 1,1)  =  7
  tvsch( 2,1)  =  8
  !tvsch( 4,1)  =  9
  !tvsch( 5,1)  = 10
  ntvsch(1)=2
  nfootp(1)=30
  !
  ! NOAA-16 AMSU-A
  !
  tvsname(2) = 'AA16'
  tvsinst(1,2) = rttv_plat_noaa
  tvsinst(2,2) = 16
  tvsinst(3,2) = rttv_inst_amsua
  tvsch( 1,2)  =  6
  ntvsch(2)=1
  nfootp(2)=30
  !
  ! NOAA-18 AMSU-A
  !
  tvsname(3) = 'AA18'
  tvsinst(1,3) = rttv_plat_noaa
  tvsinst(2,3) = 18
  tvsinst(3,3) = rttv_inst_amsua
  !tvsch( 1,3)  =  5
  tvsch( 1,3)  =  6
  tvsch( 2,3)  =  7
  tvsch( 3,3)  =  8
  ntvsch(3)=3
  nfootp(3)=30
  !
  ! NOAA-19 AMSU-A
  !
  tvsname(4) = 'AA19'
  tvsinst(1,4) = rttv_plat_noaa
  tvsinst(2,4) = 19
  tvsinst(3,4) = rttv_inst_amsua
  !tvsch( 1,4)  =  5
  tvsch( 1,4)  =  6
  tvsch( 2,4)  =  7
  !tvsch( 4,4)  =  9
  !tvsch( 5,4)  = 10
  ntvsch(4)=2
  nfootp(4)=30
  !
  ! METOP-2 AMSU-A
  !
  tvsname(5) = 'MA02'
  tvsinst(1,5) = rttv_plat_metop
  tvsinst(2,5) = 2
  tvsinst(3,5) = rttv_inst_amsua
  !tvsch(1,5)  = 5
  tvsch(1,5)  = 6
  tvsch(2,5)  = 8
  ntvsch(5)=2
  nfootp(5)=30
  !
  ! METOP-1 IASI
  !
  tvsname(6)   = 'MI01'
  tvsinst(1,6) = rttv_plat_metop
  tvsinst(2,6) = 01
  tvsinst(3,6) = rttv_inst_iasi
  tvsch(1,6)   = 212
  tvsch(2,6)   = 246
  tvsch(3,6)   = 262
  tvsch(4,6)   = 275
  ntvsch(6)    =   4
  nfootp(6)    = 120
  !
  ! METOP-2 IASI
  !
  tvsname(7)   = 'MI02'
  tvsinst(1,7) = rttv_plat_metop
  tvsinst(2,7) = 02
  tvsinst(3,7) = rttv_inst_iasi
  tvsch(1,7)   = 212
  tvsch(2,7)   = 246
  tvsch(3,7)   = 262
  tvsch(4,7)   = 275
  ntvsch(7)    = 4
  nfootp(7)    = 120
  !
  !JPSS-0 (Suomi-NPP)
  !
  tvsname(8)   = 'JP00'
  tvsinst(1,8) = rttv_plat_jpss
  tvsinst(2,8) = 0
  tvsinst(3,8) = rttv_inst_atms
  tvsch( 1,8)  =  7
  tvsch( 2,8)  =  8
  tvsch( 3,8)  =  9
  ntvsch(8)    =  3
  nfootp(8)    = 32

  RETURN
END SUBROUTINE set_instrument
!-----------------------------------------------------------------------
END MODULE common_tvs_nicam

