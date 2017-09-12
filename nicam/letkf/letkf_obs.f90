MODULE letkf_obs
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   04/03/2013 Takemasa MIYOSHI  separating obs operator
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_nicam
  USE common_obs_nicam
  USE common_tvs_nicam ! 2015.01.22 [added] koji
  USE common_mpi_nicam
  USE common_letkf
  USE mod_adm

  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs
  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov
  REAL(r_size),ALLOCATABLE,SAVE :: dlon_zero(:)
  REAL(r_size),SAVE :: dlat_zero
  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obstyp(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
  INTEGER,SAVE :: nobsgrd(nlon,nlat)

  REAL(r_size),ALLOCATABLE,SAVE :: tvselm     (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvslon     (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvslat     (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvslev     (:,:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvszenith  (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvsskin    (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvsstmp    (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvsclw     (:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvsemis  (:,:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvsdat   (:,:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvserr   (:,:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvsdep   (:,:,:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: tvshdxf(:,:,:,:,:)
  !REAL(r_size),ALLOCATABLE,SAVE :: tvswgt (:,:,:,:,:)
  INTEGER,ALLOCATABLE,SAVE :: tvsqc(:,:,:,:)
  INTEGER,ALLOCATABLE,SAVE :: tvsfoot(:,:,:)
  INTEGER,ALLOCATABLE,SAVE :: ntvsgrd(:,:,:,:)
  !INTEGER,SAVE :: ntvsgrd(nlon,nlat,ninstrument,nslots)

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: gross_error=10.0d0
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: wk2d(:,:)
  INTEGER,ALLOCATABLE :: iwk2d(:,:)
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  !REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmptyp(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:)
  REAL(r_size),ALLOCATABLE :: tmp2elm(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:)
!  REAL(r_size),ALLOCATABLE :: tmp2k(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:)
  REAL(r_size),ALLOCATABLE :: tmp2typ(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
  INTEGER,ALLOCATABLE :: nobslots(:)
  !INTEGER :: nobslots(nslots)
  INTEGER :: n,i,j,ierr,islot,nn,l,im, im0
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(21) :: obsfile='prepbufr_TTNNNNNN.dat'
  CHARACTER(256) :: obsfname=''
  LOGICAL :: readobs_dummy = .true.

  WRITE(ADM_LOG_FID,'(A)') 'Hello from set_letkf_obs'

  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
  dlat_zero = dist_zero / pi / re * 180.0d0
  ALLOCATE(dlon_zero(nij1))
  DO i=1,nij1
    dlon_zero(i) = dlat_zero / COS(pi*lat1(i)/180.0d0)
  END DO

  ALLOCATE( nobslots(nslots) )
  IF(myrank == 0) THEN !Assuming all members have the identical obs records
    DO islot=1,nslots
      im = myrank+1
      IF(start_mem_zero) im=myrank
      WRITE(obsfile(10:17),'(I2.2,I6.6)') islot,im
      obsfname=trim(obsprep_basedir)//'/'//trim(obsfile)
      WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ', TRIM(obsfname)
      CALL get_nobs(obsfname,8,nobslots(islot))
    END DO
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nobslots,nslots,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nobs = SUM(nobslots)
  WRITE(ADM_LOG_FID,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'

  IF(nobs == 0) THEN
    WRITE(ADM_LOG_FID,'(A)') 'No observation assimilated'
    RETURN
  END IF
!
! INITIALIZE GLOBAL VARIABLES
!
  ALLOCATE( tmpelm(nobs) )
  ALLOCATE( tmplon(nobs) )
  ALLOCATE( tmplat(nobs) )
  ALLOCATE( tmplev(nobs) )
  ALLOCATE( tmpdat(nobs) )
  ALLOCATE( tmperr(nobs) )
  ALLOCATE( tmpdep(nobs) )
  ALLOCATE( tmptyp(nobs) )
  ALLOCATE( tmphdxf(nobs,nbv) )
  ALLOCATE( tmpqc0(nobs,nbv) )
  ALLOCATE( tmpqc(nobs) )
  tmpqc0 = 0
  tmphdxf = 0.0d0
!
! reading observation data
!
  nn=0
  DO islot=1,nslots
    IF(nobslots(islot) == 0) CYCLE
    l=0
    DO
      IF(start_mem_zero) THEN
        im  = myrank + nprocs * l
        im0 = im + 1
        IF(im > nbv-1) EXIT
        readobs_dummy = .false.
      ELSE
        im  = myrank+1 + nprocs * l
        im0 = im
        IF(im > nbv) EXIT
        readobs_dummy = .false.
      END IF
      WRITE(obsfile(10:17),'(I2.2,I6.6)') islot,im
      obsfname=trim(obsprep_basedir)//'/'//trim(obsfile) !!! 
      WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK ', myrank,' is reading a file ', TRIM(obsfname)
      CALL read_obs2(obsfname,nobslots(islot),&
       & tmpelm(nn+1:nn+nobslots(islot)),      tmplon(nn+1:nn+nobslots(islot)),     &
       & tmplat(nn+1:nn+nobslots(islot)),      tmplev(nn+1:nn+nobslots(islot)),     &
       & tmpdat(nn+1:nn+nobslots(islot)),      tmperr(nn+1:nn+nobslots(islot)),     &
       & tmphdxf(nn+1:nn+nobslots(islot),im0), tmpqc0(nn+1:nn+nobslots(islot),im0), &
       & tmptyp(nn+1:nn+nobslots(islot)))
      l = l+1
    END DO
    nn = nn + nobslots(islot)
  END DO

  IF(readobs_dummy) THEN
    nn=0
    DO islot=1,nslots
      IF(nobslots(islot) == 0) CYCLE
      WRITE(obsfile(10:17),'(I2.2,I6.6)') islot,1
      obsfname=trim(obsprep_basedir)//'/'//trim(obsfile)
      WRITE(ADM_LOG_FID,'(A,I3.3,3A)') 'MYRANK ', myrank, ' is reading a file ', TRIM(obsfname),' [do not use h(x) and qc]'
      CALL read_obs2(obsfname,nobslots(islot),&
       & tmpelm(nn+1:nn+nobslots(islot)), tmplon(nn+1:nn+nobslots(islot)), &
       & tmplat(nn+1:nn+nobslots(islot)), tmplev(nn+1:nn+nobslots(islot)), &
       & tmpdat(nn+1:nn+nobslots(islot)), tmperr(nn+1:nn+nobslots(islot)), &
       & tmpdep(nn+1:nn+nobslots(islot)), tmpqc(nn+1:nn+nobslots(islot)),  & ! tmpdep and tmpqc here are just dummy variables
       & tmptyp(nn+1:nn+nobslots(islot))) 
      nn = nn + nobslots(islot)
    END DO
  END IF
  FLUSH(ADM_LOG_FID)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ALLOCATE(wk2d(nobs,nbv))
  wk2d = tmphdxf
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(wk2d,tmphdxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  DEALLOCATE(wk2d)
  ALLOCATE(iwk2d(nobs,nbv))
  iwk2d = tmpqc0
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(iwk2d,tmpqc0,nobs*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  DEALLOCATE(iwk2d)

  DO n=1,nobs
    tmpqc(n) = MINVAL(tmpqc0(n,:))
    IF(tmpqc(n) /= 1) CYCLE
    tmpdep(n) = tmphdxf(n,1)
    DO i=2,nbv
      tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
    END DO
    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)
    DO i=1,nbv
      tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
    END DO
    tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx
    IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
      tmpqc(n) = -1
    END IF
  END DO
  DEALLOCATE(tmpqc0)

  if( myrank == 0 ) then
    open(100,file='monit_obs_prepbufr.txt')
    do n = 1, nobs
      write(100,'(9F10.2)') tmpelm(n), tmplon(n), tmplat(n), &
                            tmplev(n), tmpdat(n), tmperr(n), &
                            tmpdep(n), real(tmpqc(n)),  tmptyp(n)
    end do
    close(100)
  end if

  do n = 1, nobs
    if(tmpqc(n)/=1) tmpqc(n)=0
  end do
  WRITE(ADM_LOG_FID,'(I10,A)') SUM(tmpqc),' OBSERVATIONS TO BE ASSIMILATED'

  CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)
!
! temporal observation localization
!
  nn = 0
  DO islot=1,nslots
    tmperr(nn+1:nn+nobslots(islot)) = tmperr(nn+1:nn+nobslots(islot)) &
      & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**3)
    nn = nn + nobslots(islot)
  END DO
!
! SELECT OBS IN THE NODE
!
  nn = 0
  DO n=1,nobs
    IF(tmpqc(n) /= 1) CYCLE
    nn = nn+1
    tmpelm(nn) = tmpelm(n)
    tmplon(nn) = tmplon(n)
    tmplat(nn) = tmplat(n)
    tmplev(nn) = tmplev(n)
    tmpdat(nn) = tmpdat(n)
    tmperr(nn) = tmperr(n)
    tmpdep(nn) = tmpdep(n)
    tmptyp(nn) = tmptyp(n)
    tmphdxf(nn,:) = tmphdxf(n,:)
    tmpqc(nn) = tmpqc(n)
  END DO
  nobs = nn
  WRITE(ADM_LOG_FID,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank

!
! SORT
!
  ALLOCATE( tmp2elm(nobs) )
  ALLOCATE( tmp2lon(nobs) )
  ALLOCATE( tmp2lat(nobs) )
  ALLOCATE( tmp2lev(nobs) )
  ALLOCATE( tmp2dat(nobs) )
  ALLOCATE( tmp2err(nobs) )
  ALLOCATE( tmp2dep(nobs) )
  ALLOCATE( tmp2typ(nobs) )
  ALLOCATE( tmp2hdxf(nobs,nbv) )
  ALLOCATE( obselm(nobs) )
  ALLOCATE( obslon(nobs) )
  ALLOCATE( obslat(nobs) )
  ALLOCATE( obslev(nobs) )
  ALLOCATE( obsdat(nobs) )
  ALLOCATE( obserr(nobs) )
  ALLOCATE( obsdep(nobs) )
  ALLOCATE( obstyp(nobs) )
  ALLOCATE( obshdxf(nobs,nbv) )
  nobsgrd = 0
  nj = 0


!!$OMP PARALLEL PRIVATE(i,j,n,nn)
!!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
  DO n=1,nobs
    IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
    nj(j) = nj(j) + 1
  END DO
  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    njs(j) = SUM(nj(0:j-1))
  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    nn = 0
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
      nn = nn + 1
      tmp2elm(njs(j)+nn) = tmpelm(n)
      tmp2lon(njs(j)+nn) = tmplon(n)
      tmp2lat(njs(j)+nn) = tmplat(n)
      tmp2lev(njs(j)+nn) = tmplev(n)
      tmp2dat(njs(j)+nn) = tmpdat(n)
      tmp2err(njs(j)+nn) = tmperr(n)
!      tmp2k(njs(j)+nn) = tmpk(n)
      tmp2dep(njs(j)+nn) = tmpdep(n)
      tmp2typ(njs(j)+nn) = tmptyp(n)
      tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
    END DO
  END DO
!!$OMP END DO
!!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    IF(nj(j) == 0) THEN
      nobsgrd(:,j) = njs(j)
      CYCLE
    END IF
    nn = 0
    DO i=1,nlon
    DO n=njs(j)+1,njs(j)+nj(j)
      IF(i < nlon) THEN
        IF(tmp2lon(n) < lon(i) .OR. lon(i+1) <= tmp2lon(n)) CYCLE
      ELSE
        IF(tmp2lon(n) < lon(nlon) .OR. 360.0d0 <= tmp2lon(n)) CYCLE
      END IF
      nn = nn + 1
      obselm(njs(j)+nn) = tmp2elm(n)
      obslon(njs(j)+nn) = tmp2lon(n)
      obslat(njs(j)+nn) = tmp2lat(n)
      obslev(njs(j)+nn) = tmp2lev(n)
      obsdat(njs(j)+nn) = tmp2dat(n)
      obserr(njs(j)+nn) = tmp2err(n)
!      obsk(njs(j)+nn) = tmp2k(n)
      obsdep(njs(j)+nn) = tmp2dep(n)
      obstyp(njs(j)+nn) = tmp2typ(n)
      obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
    END DO
    nobsgrd(i,j) = njs(j) + nn
    END DO
    IF(nn /= nj(j)) THEN
!!$OMP CRITICAL
      WRITE(ADM_LOG_FID,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
      WRITE(ADM_LOG_FID,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
      WRITE(ADM_LOG_FID,'(F6.2,A,F6.2)') MINVAL(tmp2lat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(tmp2lat(njs(j)+1:njs(j)+nj(j)))
!!$OMP END CRITICAL
    END IF
  END DO
!!$OMP END DO
!!$OMP END PARALLEL
  DEALLOCATE( tmp2elm )
  DEALLOCATE( tmp2lon )
  DEALLOCATE( tmp2lat )
  DEALLOCATE( tmp2lev )
  DEALLOCATE( tmp2dat )
  DEALLOCATE( tmp2err )
  DEALLOCATE( tmp2dep )
  DEALLOCATE( tmp2typ )
  DEALLOCATE( tmp2hdxf )
  DEALLOCATE( tmpelm )
  DEALLOCATE( tmplon )
  DEALLOCATE( tmplat )
  DEALLOCATE( tmplev )
  DEALLOCATE( tmpdat )
  DEALLOCATE( tmperr )
  DEALLOCATE( tmpdep )
  DEALLOCATE( tmptyp )
  DEALLOCATE( tmphdxf )
  DEALLOCATE( tmpqc )

  RETURN
END SUBROUTINE set_letkf_obs
!-----------------------------------------------------------------------
! Satelite Data
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_tvs
  IMPLICIT NONE
  INTEGER,PARAMETER :: err_unit=6
  INTEGER,PARAMETER :: verbosity_level=1
  REAL(r_size),PARAMETER :: gross_error=3.0d0
  CHARACTER(16) :: cfile='inst00000000.dat'

  REAL(r_size),ALLOCATABLE :: tmpelm(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmplon(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmplat(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpzenith(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpskin(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpstmp(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpclw(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpemis(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmplev(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmplev_tmp(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpdat_tmp(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmperr(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:,:,:,:)
  !REAL(r_size),ALLOCATABLE :: tmpwgt(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpdum(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmppred(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpqc_out(:,:,:,:)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:,:,:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:,:,:,:)
  INTEGER,ALLOCATABLE :: tmpfoot(:,:,:)

  REAL(r_size),ALLOCATABLE :: tmp2elm(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2zenith(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2skin(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2stmp(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2clw(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2emis(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: depstat(:,:,:)
  INTEGER,ALLOCATABLE :: num_depstat(:,:,:)
  INTEGER,ALLOCATABLE :: tmp2qc(:,:,:,:)
  INTEGER,ALLOCATABLE :: tmp2foot(:,:,:)

  INTEGER,ALLOCATABLE      :: wk5i(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: wk5d(:,:,:,:,:)
  REAL(r_size) :: wgtmax

  INTEGER :: n,i,j,ierr,islot,nn,ic,nnn,itmp,im,l,ifoot,im0
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  INTEGER :: ntvsinput
  
  LOGICAL :: readobs_dummy = .true.

  ALLOCATE(ntvsprofslots(ninstrument,nslots))
  ALLOCATE(ntvsgrd(nlon,nlat,ninstrument,nslots))
  ntvsgrd = 0
  ntvsprofslots=0
  !ntvschan=0
  WRITE(ADM_LOG_FID,'(A)') 'Hello from set_letkf_tvs'
  CALL set_instrument

  IF( myrank==0 ) THEN
    im = myrank + 1
    IF(start_mem_zero) im=myrank
    DO islot=1,nslots
      WRITE(cfile(5:12),'(I2.2,I6.6)') islot, im
      write(ADM_LOG_FID,*) cfile
      CALL get_ntvs_mpi(cfile,dir=trim(obssate_basedir))
      ntvsprofslots(:,islot) = ntvsprof(:)
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ntvsprofslots,ninstrument*nslots,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  maxtvsprof=MAXVAL(ntvsprofslots)
  maxtvsfoot=MAXVAL(nfootp)

  ntvsinput=0
  DO n = 1,ninstrument
    ntvsinput=ntvsinput+sum(ntvsprofslots(n,:))*ntvsch(n)
  END DO
  WRITE(ADM_LOG_FID,'(I10,A)') ntvsinput, ' TOTAL Satellite Radiances INPUT'
  WRITE(ADM_LOG_FID,*) "maxtvsprof", maxtvsprof
 
  IF(maxtvsprof /= 0) THEN
    ALLOCATE( tmpelm(   maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmplon(   maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmplat(   maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpzenith(maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpskin(  maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpstmp(  maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpclw (  maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpemis(  maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmplev(   maxtvsch,maxtvsprof,ninstrument,nslots,nbv) )
    ALLOCATE( tmplev_tmp(   maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpdat(   maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpdat_tmp(   maxtvsch,maxtvsprof,ninstrument,nslots,nbv) )
    ALLOCATE( tmperr(   maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpdep(   maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmphdxf(  maxtvsch,maxtvsprof,ninstrument,nslots,nbv) )
    ALLOCATE( tmpdum(   nlev,maxtvsch,maxtvsprof,ninstrument,nslots) )
    !ALLOCATE( tmpwgt(   nlev,maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmppred( maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpqc0(    maxtvsch,maxtvsprof,ninstrument,nslots,nbv) )
    ALLOCATE( tmpqc(    maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpqc_out(    maxtvsch,maxtvsprof,ninstrument,nslots) )
    ALLOCATE( tmpfoot(  maxtvsprof,ninstrument,nslots) )
    ALLOCATE( depstat(  maxtvsch,maxtvsfoot,ninstrument) )
    ALLOCATE( num_depstat(  maxtvsch,maxtvsfoot,ninstrument) )

    tmpqc0 =0
    tmphdxf=0.d0
    tmplev=0.d0
    !tmpwgt =0.d0
    tmperr = 0.0d0
    tmpdat = 0.0d0
    tmpdep = 0.0d0
    tmpdat_tmp = 0.0d0

    DO islot=1,nslots
      ntvsprof(:) = ntvsprofslots(:,islot)
      IF(MAXVAL(ntvsprof) /= 0) THEN
        l=0
        DO
          IF(start_mem_zero) THEN
            im  = myrank + nprocs * l
            im0 = im + 1 
            IF(im > nbv-1) EXIT
            readobs_dummy = .false.
          ELSE
            im  = myrank+1 + nprocs * l
            im0 = im
            IF(im > nbv) EXIT
            readobs_dummy = .false.
          END IF
          WRITE(cfile(5:12),'(I2.2,I6.6)') islot, im
          WRITE(ADM_LOG_FID,*) cfile
          FLUSH(ADM_LOG_FID)
          CALL read_tvs_mpi(cfile, &
            &    tmpelm   (1,1,islot), tmplon   (1,1,islot), tmplat (1,1,islot), &
            &    tmpzenith(1,1,islot), tmpskin(1,1,islot), &
            &    tmpstmp  (1,1,islot), tmpclw (1,1,islot), &
            &    tmplev (1,1,1,islot,im0), tmpdat_tmp(1,1,1,islot,im0), tmperr (1,1,1,islot),&
            &    tmphdxf(1,1,1,islot,im0), tmpqc0 (1,1,1,islot,im0), tmpfoot(1,1,islot), dir=trim(obssate_basedir)) 
          l=l+1
        END DO
      END IF
    END DO

    IF(readobs_dummy) THEN
      nn=0
      DO islot=1,nslots
        ntvsprof(:) = ntvsprofslots(:,islot)
        IF(MAXVAL(ntvsprof) /= 0) THEN
          WRITE(cfile(5:12),'(I2.2,I6.6)') islot, 1
          WRITE(ADM_LOG_FID,'(A,I3.3,3A)') 'MYRANK ',myrank,' is reading a file ', cfile,' [do not use h(x) and qc]'
          FLUSH(ADM_LOG_FID)
          CALL read_tvs_mpi(cfile, &
            &    tmpelm   (1,1,islot), tmplon   (1,1,islot), tmplat (1,1,islot), &
            &    tmpzenith(1,1,islot), tmpskin(1,1,islot), &
            &    tmpstmp  (1,1,islot), tmpclw (1,1,islot), &
            &    tmpdum (1,1,1,1,islot), tmpdat (1,1,1,islot), tmperr (1,1,1,islot),&
            &    tmpdep(1,1,1,islot), tmpqc (1,1,1,islot), tmpfoot(1,1,islot), dir=trim(obssate_basedir))
        END IF
      END DO
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
    ALLOCATE(wk5d(maxtvsch,maxtvsprof,ninstrument,nslots,nbv))
    wk5d=tmphdxf
    nn=maxtvsch*maxtvsprof*ninstrument*nslots*nbv
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(wk5d,tmphdxf,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    DEALLOCATE(wk5d)
  
    ALLOCATE(wk5d(maxtvsch,maxtvsprof,ninstrument,nslots,nbv))
    wk5d=tmplev
    nn=maxtvsch*maxtvsprof*ninstrument*nslots*nbv
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(wk5d,tmplev,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    DEALLOCATE(wk5d)

    ALLOCATE(wk5i(maxtvsch,maxtvsprof,ninstrument,nslots,nbv))
    wk5i=tmpqc0
    nn=maxtvsch*maxtvsprof*ninstrument*nslots*nbv
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(wk5i,tmpqc0,nn,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    DEALLOCATE(wk5i)
  
    ALLOCATE(wk5d(maxtvsch,maxtvsprof,ninstrument,nslots,nbv))
    wk5d=tmpdat_tmp
    nn=maxtvsch*maxtvsprof*ninstrument*nslots*nbv
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(wk5d,tmpdat_tmp,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    DEALLOCATE(wk5d)

    WHERE(tmpqc0<0) tmpqc0=0 !! Koji 2016.05.25

    DO islot=1,nslots
    DO nn = 1,ninstrument
    DO i = 1,ntvsprofslots(nn,islot)
    DO ic= 1,ntvsch(nn)
      tmpdat(ic,i,nn,islot)=sum(tmpdat_tmp(ic,i,nn,islot,1:nbv))
    END DO
    END DO
    END DO
    END DO

    tmpdat(:,:,:,:)=tmpdat(:,:,:,:)/dble(nbv)
    
    write(ADM_LOG_FID,*)  'CHECK Averaging observations'
    DO nn = 1,ninstrument
    DO ic= 1,ntvsch(nn)
      write(ADM_LOG_FID,'(2i6,3f12.7)') &
                  nn, ic, minval(tmpdat_tmp(ic,1,nn,4,1:nbv)), &
                          maxval(tmpdat_tmp(ic,1,nn,4,1:nbv)), &
                          tmpdat(ic,1,nn,4)
    END DO
    END DO

    DO islot=1,nslots
    DO nn = 1,ninstrument
    DO i = 1,ntvsprofslots(nn,islot)
    DO ic= 1,ntvsch(nn)
      tmplev_tmp(ic,i,nn,islot)=sum(log(tmplev(ic,i,nn,islot,1:nbv)))
    END DO
    END DO
    END DO
    END DO

    DO islot=1,nslots
    DO nn = 1,ninstrument
    DO i = 1,ntvsprofslots(nn,islot)
    DO ic= 1,ntvsch(nn)
      tmplev(ic,i,nn,islot,1)=exp( tmplev_tmp(ic,i,nn,islot)/real(nbv) )
    END DO
    END DO
    END DO
    END DO

    !islot=1
    !nn=1
    !DO i = 1,ntvsprofslots(nn,islot)
    !  write(ADM_LOG_FID,*) 'average, exponential'
    !  write(ADM_LOG_FID,*) (tmplev(ic,i,nn,islot,1),ic=1,ntvsch(nn))
    !END DO
    !FLUSH(ADM_LOG_FID)
    !ALLOCATE(wk5d(nlev,maxtvsch,maxtvsprof,ninstrument,nslots))
    !wk5d=tmpwgt
    !nn=nlev*maxtvsch*maxtvsprof*ninstrument*nslots
    !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !CALL MPI_ALLREDUCE(wk5d,tmpwgt,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    !DEALLOCATE(wk5d)
    !
    ! Hdx & y-Hx
    !
    tmpqc = 0
    DO islot=1,nslots
    DO nn=1,ninstrument
    DO n=1,ntvsprofslots(nn,islot)
    DO ic=1,ntvsch(nn)
      tmpqc(ic,n,nn,islot) = SUM(tmpqc0(ic,n,nn,islot,:))
      IF(REAL(tmpqc(ic,n,nn,islot)) < 0.6*REAL(nbv)) THEN
        tmpqc_out(ic,n,nn,islot)=-REAL(tmpqc(ic,n,nn,islot))
        tmpqc(ic,n,nn,islot)=0
      ELSE
        tmpqc_out(ic,n,nn,islot)=REAL(tmpqc(ic,n,nn,islot))
        tmpqc(ic,n,nn,islot)=1
      END IF 

      tmpdep(ic,n,nn,islot) = tmphdxf(ic,n,nn,islot,1)
      DO i=2,nbv
        tmpdep(ic,n,nn,islot) = tmpdep(ic,n,nn,islot) &
                            & + tmphdxf(ic,n,nn,islot,i)
      END DO
      tmpdep(ic,n,nn,islot) = tmpdep(ic,n,nn,islot) / REAL(nbv,r_size)
      DO i=1,nbv
        tmphdxf(ic,n,nn,islot,i) = tmphdxf(ic,n,nn,islot,i) &
                               & - tmpdep(ic,n,nn,islot) ! Hdx
      END DO
      tmpdep(ic,n,nn,islot) = tmpdat(ic,n,nn,islot) &
                          & - tmpdep(ic,n,nn,islot) ! y-Hx
      IF(ABS(tmpdep(ic,n,nn,islot)) > gross_error*tmperr(ic,n,nn,islot)) THEN ! Gross error check
        tmpqc_out(ic,n,nn,islot)=-2.0*real(nbv)
        tmpqc (ic,n,nn,islot) = 0
      ELSE IF(tmpdat(ic,n,nn,islot)/=tmpdat(ic,n,nn,islot)) then
        tmpqc_out(ic,n,nn,islot)=-3.0*real(nbv)
        tmpqc (ic,n,nn,islot) = 0
      END IF
    END DO
    END DO
    END DO
    END DO

    if( myrank == 0 ) then
      depstat(:,:,:)=0.0d0
      num_depstat(:,:,:)=0
      open(100,file='monit_obs_amsua.txt')
      open(101,file='stat_obs_amsua.txt')
      do islot = 1, nslots
      do nn = 1, ninstrument
      do n = 1, ntvsprofslots(nn,islot)
      do ic = 1, ntvsch(nn)
        if(tmpqc(ic,n,nn,islot)==1) then
        !if(tvsch(ic,nn)>=5 .and. tmpqc(ic,n,nn,islot)==1) then
          depstat(ic,tmpfoot(n,nn,islot),nn)=depstat(ic,tmpfoot(n,nn,islot),nn)+&
                                          tmpdep(ic,n,nn,islot)
          num_depstat(ic,tmpfoot(n,nn,islot),nn)=num_depstat(ic,tmpfoot(n,nn,islot),nn)+1
          !write(ADM_LOG_FID,'(5i5,F9.3)') islot, nn, n, ic,
          !tmpfoot(n,nn,islot), tmpdep(ic,n,nn,islot) 
        !else if(tvsch(ic,nn)<=4) then
        !  depstat(ic,tmpfoot(n,nn,islot),nn)=depstat(ic,tmpfoot(n,nn,islot),nn)+&
        !                                  tmpdep(ic,n,nn,islot)
        !  num_depstat(ic,tmpfoot(n,nn,islot),nn)=num_depstat(ic,tmpfoot(n,nn,islot),nn)+1
        end if
      end do
      end do
      end do
      end do

      do nn = 1, ninstrument
      do ifoot = 1, maxtvsfoot
      do ic = 1, ntvsch(nn)
        if(num_depstat(ic,ifoot,nn)==0) then
          depstat(ic,ifoot,nn)=0.0d0
          !depstat(ic,ifoot,nn)=-999.d0
        else
          depstat(ic,ifoot,nn)=depstat(ic,ifoot,nn)/real(num_depstat(ic,ifoot,nn),kind=8)
        end if
      end do
      end do
      end do

      do nn = 1, ninstrument
      do ifoot = 1, maxtvsfoot
        write(101,'(A,i5,20F10.3)') tvsname(nn), ifoot, &
             (depstat(ic,ifoot,nn),ic=1,ntvsch(nn)), &
             (real(num_depstat(ic,ifoot,nn)),ic=1,ntvsch(nn))
      end do
      end do
      close(101)

      do islot = 1, nslots
      do nn = 1, ninstrument
      do n = 1, ntvsprofslots(nn,islot)
        write(100,'(A,50F10.3)') tvsname(nn), tmplon(n,nn,islot),          &
                                tmplat(n,nn,islot),                       &
                                real(tmpfoot(n,nn,islot)),                &
                                tmplev(1:ntvsch(nn),n,nn,islot,1)*0.01d0, &
                                tmpdat(1:ntvsch(nn),n,nn,islot),          &
                                tmpdep(1:ntvsch(nn),n,nn,islot),          &
                                tmpqc_out(1:ntvsch(nn),n,nn,islot),       &
                                tmpskin(n,nn,islot)
      end do
      end do
      end do
      close(100)
    end if

  END IF
!
! temporal observation localization
!
    DO islot=1,nslots
    DO nn=1,ninstrument
    DO n=1,ntvsprofslots(nn,islot)
    DO ic=1,ntvsch(nn)
      IF( tmpqc(ic,n,nn,islot) == 1 ) THEN
        tmperr(ic,n,nn,islot) = tmperr(ic,n,nn,islot) &
         & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
      END IF
    END DO
    END DO
    END DO
    END DO
!
! SELECT OBS IN THE NODE
!
    ntvs=0
    DO islot=1,nslots
    DO nn=1,ninstrument
      nnn = 0
      DO n=1,ntvsprofslots(nn,islot)
        IF( MAXVAL(tmpqc(:,n,nn,islot)) /= 1 ) CYCLE
        nnn = nnn + 1
        tmpelm    (nnn,nn,islot)  =tmpelm    (n,nn,islot)
        tmplon    (nnn,nn,islot)  =tmplon    (n,nn,islot)
        tmplat    (nnn,nn,islot)  =tmplat    (n,nn,islot)
        tmplev  (:,nnn,nn,islot,1)=tmplev    (:,n,nn,islot,1)
        tmpzenith (nnn,nn,islot)  =tmpzenith (n,nn,islot)
        tmpskin   (nnn,nn,islot)  =tmpskin   (n,nn,islot)
        tmpstmp   (nnn,nn,islot)  =tmpstmp   (n,nn,islot)
        tmpclw    (nnn,nn,islot)  =tmpclw    (n,nn,islot)
        tmpemis (:,nnn,nn,islot)  =tmpemis (:,n,nn,islot)
        tmpdat  (:,nnn,nn,islot)  =tmpdat  (:,n,nn,islot)
        tmperr  (:,nnn,nn,islot)  =tmperr  (:,n,nn,islot)
        tmpdep  (:,nnn,nn,islot)  =tmpdep  (:,n,nn,islot)
        tmpqc   (:,nnn,nn,islot)  =tmpqc   (:,n,nn,islot)
        tmpfoot(nnn,nn,islot)   =tmpfoot(n,nn,islot)
        !tmpwgt(:,:,nnn,nn,islot)  =tmpwgt(:,:,n,nn,islot)
        DO i=1,nbv
          tmphdxf (:,nnn,nn,islot,i)=tmphdxf (:,n,nn,islot,i)
        END DO
        ntvs=ntvs+SUM(tmpqc(:,nnn,nn,islot))
      END DO
      ntvsprofslots(nn,islot) = nnn
    END DO
    END DO
    maxtvsprof=MAXVAL(ntvsprofslots)
    WRITE(ADM_LOG_FID,'(I10,A,I3.3)') ntvs,' ATOVS CHANNELS TO BE ASSIMILATED IN MYRANK ',myrank
    FLUSH(ADM_LOG_FID)
!
! SORT
!
  IF(maxtvsprof /= 0) THEN
    ALLOCATE( tmp2elm   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2lon   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2lat   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2lev   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2zenith(maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2skin  (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2stmp  (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2clw   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2emis  (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2dat   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2err   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2dep   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2hdxf  (nbv, maxtvsch,maxtvsprof,ninstrument,nslots))
    !ALLOCATE( tmp2wgt   (nlev,maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2qc    (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tmp2foot  (maxtvsprof,ninstrument,nslots))
!
    ALLOCATE( tvselm  (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsskin  (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvslon   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvslat   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvslev   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvszenith(maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsstmp  (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsclw   (maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsemis  (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsdat   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvserr   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsdep   (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvshdxf  (nbv, maxtvsch,maxtvsprof,ninstrument,nslots))
    !ALLOCATE( tvswgt   (nlev,maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsqc    (maxtvsch,maxtvsprof,ninstrument,nslots))
    ALLOCATE( tvsfoot  (maxtvsprof,ninstrument,nslots))
    ntvsgrd = 0
    tmp2qc = 0
    tvsqc = 0
    DO itmp = 1, ninstrument * nslots
      islot =     (itmp-1)/ninstrument +1
      nn    = mod((itmp-1),ninstrument)+1
      IF( ntvsprofslots(nn,islot) == 0) CYCLE
      nj = 0
      DO j=1,nlat-1
      DO n=1,ntvsprofslots(nn,islot)
        IF(tmplat(n,nn,islot) < lat(j) .OR. &
         & tmplat(n,nn,islot) >=lat(j+1)) CYCLE
        nj(j) = nj(j) + 1
      END DO
      END DO
      DO j=1,nlat-1
        njs(j) = SUM(nj(0:j-1))
      END DO
      DO j=1,nlat-1
        nnn = 0
        DO n=1,ntvsprofslots(nn,islot)
          IF(tmplat(n,nn,islot) < lat(j) .OR. &
           & tmplat(n,nn,islot) >=lat(j+1)) CYCLE
          nnn = nnn + 1
          tmp2elm    (njs(j)+nnn,nn,islot)  =tmpelm    (n,nn,islot)
          tmp2lon    (njs(j)+nnn,nn,islot)  =tmplon    (n,nn,islot)
          tmp2lat    (njs(j)+nnn,nn,islot)  =tmplat    (n,nn,islot)
          tmp2zenith (njs(j)+nnn,nn,islot)  =tmpzenith (n,nn,islot)
          tmp2skin   (njs(j)+nnn,nn,islot)  =tmpskin   (n,nn,islot)
          tmp2stmp   (njs(j)+nnn,nn,islot)  =tmpstmp   (n,nn,islot)
          tmp2clw    (njs(j)+nnn,nn,islot)  =tmpclw    (n,nn,islot)
          tmp2foot   (njs(j)+nnn,nn,islot)  =tmpfoot   (n,nn,islot)
          DO ic=1,ntvsch(nn)
            tmp2lev  (ic,njs(j)+nnn,nn,islot)  =tmplev  (ic,n,nn,islot,1)
            tmp2emis (ic,njs(j)+nnn,nn,islot)  =tmpemis (ic,n,nn,islot)
            tmp2dat  (ic,njs(j)+nnn,nn,islot)  =tmpdat  (ic,n,nn,islot)
            tmp2err  (ic,njs(j)+nnn,nn,islot)  =tmperr  (ic,n,nn,islot)
            tmp2dep  (ic,njs(j)+nnn,nn,islot)  =tmpdep  (ic,n,nn,islot)
            tmp2qc   (ic,njs(j)+nnn,nn,islot)  =tmpqc   (ic,n,nn,islot)
            DO i=1,nbv
              tmp2hdxf (i,ic,njs(j)+nnn,nn,islot)=tmphdxf (ic,n,nn,islot,i)
            END DO
            !tmp2wgt(:,ic,njs(j)+nnn,nn,islot)  =tmpwgt(:,ic,n,nn,islot)
          END DO
        END DO
      END DO
      DO j=1,nlat-1
        IF(nj(j) == 0) THEN
          ntvsgrd(:,j,nn,islot) = njs(j)
          CYCLE
        END IF
        nnn = 0
        DO i=1,nlon
        DO n=njs(j)+1,njs(j)+nj(j)
          IF(i < nlon) THEN
            IF(tmp2lon(n,nn,islot) < lon(i) .OR. &
             & tmp2lon(n,nn,islot) >=lon(i+1)) CYCLE
          ELSE
            IF(tmp2lon(n,nn,islot) < lon(nlon) .OR. &
             & tmp2lon(n,nn,islot) >=360.0d0) CYCLE
          END IF
          nnn = nnn + 1
          tvselm     (njs(j)+nnn,nn,islot) = tmp2elm     (n,nn,islot)
          tvslon     (njs(j)+nnn,nn,islot) = tmp2lon     (n,nn,islot)
          tvslat     (njs(j)+nnn,nn,islot) = tmp2lat     (n,nn,islot)
          tvslev   (:,njs(j)+nnn,nn,islot) = tmp2lev   (:,n,nn,islot)
          tvszenith  (njs(j)+nnn,nn,islot) = tmp2zenith  (n,nn,islot)
          tvsskin    (njs(j)+nnn,nn,islot) = tmp2skin    (n,nn,islot)
          tvsstmp    (njs(j)+nnn,nn,islot) = tmp2stmp    (n,nn,islot)
          tvsclw     (njs(j)+nnn,nn,islot) = tmp2clw     (n,nn,islot)
          tvsemis  (:,njs(j)+nnn,nn,islot) = tmp2emis  (:,n,nn,islot)
          tvsdat   (:,njs(j)+nnn,nn,islot) = tmp2dat   (:,n,nn,islot)
          tvserr   (:,njs(j)+nnn,nn,islot) = tmp2err   (:,n,nn,islot)
          tvsdep   (:,njs(j)+nnn,nn,islot) = tmp2dep   (:,n,nn,islot)
          tvshdxf(:,:,njs(j)+nnn,nn,islot) = tmp2hdxf(:,:,n,nn,islot)
          !tvswgt (:,:,njs(j)+nnn,nn,islot) = tmp2wgt (:,:,n,nn,islot)
          tvsqc    (:,njs(j)+nnn,nn,islot) = tmp2qc    (:,n,nn,islot)
          tvsfoot  (njs(j)+nnn,nn,islot)   = tmp2foot  (n,nn,islot)
        END DO
        ntvsgrd(i,j,nn,islot) = njs(j) + nnn
        END DO
        IF(nnn /= nj(j)) THEN
          WRITE(ADM_LOG_FID,'(A,2I6)') 'ATOVS DATA SORT ERROR: ',nn,nj(j)
          WRITE(ADM_LOG_FID,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
          WRITE(ADM_LOG_FID,'(F6.2,A,F6.2)') MINVAL(tmp2lat(njs(j)+1:njs(j)+nj(j),nn,islot)),' &
                  & < TVSLAT <',MAXVAL(tmp2lat(njs(j)+1:njs(j)+nj(j),nn,islot))
          FLUSH(ADM_LOG_FID)
        END IF
      END DO
    END DO
!
    DEALLOCATE( tmp2elm   )
    DEALLOCATE( tmp2lon   )
    DEALLOCATE( tmp2lat   )
    DEALLOCATE( tmp2lev   )
    DEALLOCATE( tmp2zenith)
    DEALLOCATE( tmp2skin  )
    DEALLOCATE( tmp2stmp  )
    DEALLOCATE( tmp2clw  )
    DEALLOCATE( tmp2emis  )
    DEALLOCATE( tmp2dat   )
    DEALLOCATE( tmp2err   )
    DEALLOCATE( tmp2dep   )
    DEALLOCATE( tmp2hdxf  )
    !DEALLOCATE( tmp2wgt   )
    DEALLOCATE( tmp2qc   )
    DEALLOCATE( tmp2foot )
  END IF

  IF( ALLOCATED(tmplon) ) THEN
    DEALLOCATE( tmpelm    )
    DEALLOCATE( tmplon    )
    DEALLOCATE( tmplat    )
    DEALLOCATE( tmplev    )
    DEALLOCATE( tmpzenith )
    DEALLOCATE( tmpskin   )
    DEALLOCATE( tmpstmp   )
    DEALLOCATE( tmpclw   )
    DEALLOCATE( tmpemis   )
    DEALLOCATE( tmpdat    )
    DEALLOCATE( tmperr    )
    DEALLOCATE( tmpdep    )
    DEALLOCATE( tmphdxf   )
    !DEALLOCATE( tmpwgt    )
    DEALLOCATE( tmpqc     )
    DEALLOCATE( tmpfoot   )
  END IF

END SUBROUTINE set_letkf_tvs

END MODULE letkf_obs
