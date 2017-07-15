MODULE common_nicam
!=======================================================================
!
! [PURPOSE:] Common Information for SPEEDY
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!
!=======================================================================
!$USE OMP_LIB
  USE mod_adm
  USE common
  USE common_mpi
  use mod_fio, only : &  ! [add] H.Yashiro 20110826
    datainfo

  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER :: nbv=0
  INTEGER,PARAMETER :: nlon=160
  INTEGER,PARAMETER :: nlat=nlon/2
  INTEGER,PARAMETER :: nlev=40
  !INTEGER,PARAMETER :: nlev=40
  INTEGER,PARAMETER :: nv3d=8 ! u,v,t,q
  INTEGER,PARAMETER :: nv2d=29 ! ps,rain,slp
  ! 3D VARIABLES (ANALYZED)
  INTEGER,PARAMETER :: iv3d_p  = 1
  INTEGER,PARAMETER :: iv3d_t  = 2
  INTEGER,PARAMETER :: iv3d_u1 = 3
  INTEGER,PARAMETER :: iv3d_u2 = 4
  INTEGER,PARAMETER :: iv3d_u3 = 5
  INTEGER,PARAMETER :: iv3d_w  = 6
  INTEGER,PARAMETER :: iv3d_qv = 7
  INTEGER,PARAMETER :: iv3d_qc = 8
  INTEGER,PARAMETER :: iv3d_u  = 1
  INTEGER,PARAMETER :: iv3d_v  = 2
  ! 2D VARIABLES (ANALYZED)
  INTEGER,PARAMETER :: iv2d_ps   = 1
  INTEGER,PARAMETER :: iv2d_rain = 2
  INTEGER,PARAMETER :: iv2d_slp  = 3
  INTEGER,PARAMETER :: iv2d_ts   = 4
  INTEGER,PARAMETER :: iv2d_q2m  = 5
  INTEGER,PARAMETER :: iv2d_u10m = 6
  INTEGER,PARAMETER :: iv2d_v10m = 7
  INTEGER,PARAMETER :: iv2d_t2m  = 8
  INTEGER,PARAMETER :: iv2d_cldw = 9
  ! 2D VARIABLES (NOT UPDATED)
  INTEGER,PARAMETER :: iv2d_cfrac   = 10
  INTEGER,PARAMETER :: iv2d_lwu_toa = 11
  INTEGER,PARAMETER :: iv2d_lwd_toa = 12
  INTEGER,PARAMETER :: iv2d_swu_toa = 13
  INTEGER,PARAMETER :: iv2d_swd_toa = 14
  INTEGER,PARAMETER :: iv2d_lwu_sfc = 15
  INTEGER,PARAMETER :: iv2d_lwd_sfc = 16
  INTEGER,PARAMETER :: iv2d_swu_sfc = 17
  INTEGER,PARAMETER :: iv2d_swd_sfc = 18
  INTEGER,PARAMETER :: iv2d_sh_sfc  = 19
  INTEGER,PARAMETER :: iv2d_lh_sfc  = 20
  INTEGER,PARAMETER :: iv2d_cldi    = 21
  INTEGER,PARAMETER :: iv2d_ll_tg   = 22
  INTEGER,PARAMETER :: iv2d_ll_ts   = 23
  INTEGER,PARAMETER :: iv2d_ll_wg   = 24
  INTEGER,PARAMETER :: iv2d_sst     = 25
  INTEGER,PARAMETER :: iv2d_ice     = 26
  INTEGER,PARAMETER :: iv2d_snw     = 27
  INTEGER,PARAMETER :: iv2d_ist     = 28
  INTEGER,PARAMETER :: iv2d_icr     = 29

  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,SAVE :: nij0
  INTEGER,SAVE :: ngpv
  REAL(r_size),ALLOCATABLE,SAVE :: ico_lon(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: ico_lat(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: ico_lon_tmp(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: ico_lat_tmp(:,:)
  REAL(r_size),SAVE :: lon(nlon)
  REAL(r_size),SAVE :: lat(nlat)
  REAL(r_size),SAVE :: sig(nlev)
  REAL(r_size),SAVE :: dx(nlat)
  REAL(r_size),SAVE :: dy(nlat)
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)
  REAL(r_size),SAVE :: GRD_z(nlev)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)
  CHARACTER(10),SAVE :: CDATE
  CHARACTER(14),SAVE :: CDATE14='' ! [add] Koji 20141024
  CHARACTER(10),SAVE :: ODATE
  CHARACTER(10),SAVE :: TDATE
  CHARACTER(128),SAVE :: gues_basedir ! [add] Koji 20140207
  CHARACTER(128),SAVE :: infl_basedir ! [add] Koji 20140208
  CHARACTER(128),SAVE :: anal_basedir ! [add] Koji 20140916
  CHARACTER(128),SAVE :: obsprep_basedir ! [add] Koji 20160706
  CHARACTER(128),SAVE :: obssate_basedir ! [add] Koji 20160706
  character(LEN=16), private, save :: restart_layername = ''

  REAL(r_size),SAVE :: dlon
  REAL(r_size),SAVE :: dlat

  INTEGER,SAVE :: nslots=1 ! number of time slots for 3D-LETKF
  INTEGER,SAVE :: nbslot=1 ! basetime slot

  REAL(r_size),SAVE :: sigma_obs=400.0d3 ! [m]
  !REAL(r_size),SAVE :: sigma_obsv=0.4d0 ! [log(hPa)]
  REAL(r_size),SAVE :: sigma_obsv=0.4d0 ! [log(hPa)]
  REAL(r_size),SAVE :: sigma_obst=3.0d0

  REAL(r_size),SAVE :: TIME_CTIME

  type(datainfo),SAVE::   dinfo

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_nicam
  USE mod_adm
  USE mpi
  use mod_misc, only :       &
       MISC_get_latlon,      &
       MISC_make_idstr
  use mod_cnst, only :       &
       CNST_setup
  use mod_comm, only :       &
       COMM_setup
  use mod_fio, only : &  ! [add] H.Yashiro 20110826
       FIO_setup,     &
       FIO_input,     &
       FIO_HLONG,     &
       FIO_INTEG_FILE,&
       FIO_SPLIT_FILE,&
       FIO_ICOSAHEDRON,&
       FIO_BIG_ENDIAN
       
  use mod_grd, only :        &
       GRD_setup,            &
       GRD_x
  use mod_gmtr, only :       &
       GMTR_setup
  use mod_mnginfo_light, only : &
    MNG_mnginfo_input,   &
    MNG_mnginfo_noinput, &
    MNG_PALL,            &
    MNG_prc_rnum,        &
    MNG_prc_tab   !(num_of_rgn,num_of_proc)
  IMPLICIT NONE
  INTEGER :: i,j
  integer :: iolen
  REAL(8) :: llat, llon
  INTEGER :: g, l
  INTEGER :: ierr
  character(128) :: icolatlon_fname = 'icolatlon_gl5_rl0.dat'
  character(128) :: zgrid_fname = ''
  character(128) :: msg = 'msg'
  character(128) :: fname
  character(LEN=FIO_HLONG)  :: mnginfo             = ''
  integer        :: fmode
  logical        :: complete = .false.
  INTEGER :: p

  namelist / current_time / &
    ADM_glevel,             &
    ADM_rlevel,             &
    ADM_vlayer,             &
    CDATE,                  &
    CDATE14,                & ! [add] Koji 20141024
    ODATE,                  &
    TDATE,                  &
    obsprep_basedir,        & ! [add] Koji 20160706
    obssate_basedir,        & ! [add] Koji 20160706
    gues_basedir,           & ! [add] Koji 20140207
    infl_basedir,           & ! [add] Koji 20140207
    anal_basedir,           & ! [add] Koji 20140916
    icolatlon_fname,        & ! [add] Koji 20140220
    zgrid_fname,            & ! [add] Koji 20140501
    nbv,                    & ! [add] Koji 20140916
    nslots,                 & ! [add] Koji 20141128
    nbslot,                 & ! [add] Koji 20141128
    sigma_obs,              & ! [add] Koji 20141128
    sigma_obsv,             & ! [add] Koji 20141128
    sigma_obst,             & ! [add] Koji 20141128
    mnginfo

  open(1,file='time.cnf')
  read(1,nml=current_time)
  if(myrank == 0) then            ! KK
    write(*,*) cdate
    write(*,*) trim(gues_basedir) ! [add] Koji 20140207
    write(*,*) trim(infl_basedir) ! [add] Koji 20140207
    write(*,*) trim(anal_basedir) ! [add] Koji 20140916
    write(*,*) trim(zgrid_fname)  ! [add] Koji 20140501
  end if                          ! KK
  close(1)

  ADM_prc_me=myrank+1
  call MISC_make_idstr(fname, trim(msg), 'pe', myrank+1)
  open(ADM_LOG_FID,file=adjustl(trim(fname)),form='formatted')

  WRITE(ADM_LOG_FID,'(A)') 'Hello from set_common_nicam'

  ADM_gall = (2**(ADM_glevel-ADM_rlevel) + 2)**2
  ADM_rgn_nmax = 10*4**ADM_rlevel
  ADM_lall = 10*4**ADM_rlevel
  !ADM_lall = 10*4**ADM_rlevel/MNG_PALL
  ADM_kall = nlev

  ALLOCATE( ico_lon(ADM_gall, ADM_rgn_nmax) ) 
  ALLOCATE( ico_lat(ADM_gall, ADM_rgn_nmax) ) 

  open(1,file=trim(icolatlon_fname),form='unformatted',access='sequential')
  read(1) ico_lon
  read(1) ico_lat

  nij0=ADM_gall*ADM_rgn_nmax
  ngpv=nij0*ADM_kall

  !
  ! Model Variables
  !
  element(iv3d_u1)        = 'U1  '
  element(iv3d_u2)        = 'U2  '
  element(iv3d_u3)        = 'U3  '
  element(iv3d_w)         = 'W   '
  element(iv3d_t)         = 'T   '
  element(iv3d_p)         = 'P   '
  element(iv3d_qv)        = 'QV  '
  element(iv3d_qc)        = 'QC  '
  element(nv3d+iv2d_ps)   = 'PS  '
  element(nv3d+iv2d_rain) = 'RAIN'
  element(nv3d+iv2d_slp)  = 'SLP'
  element(nv3d+iv2d_ts)   = 'TS'
  element(nv3d+iv2d_q2m)  = 'Q2M'
  element(nv3d+iv2d_u10m) = 'U10M'
  element(nv3d+iv2d_v10m) = 'V10M'
  !
  ! Lon, Lat, Sigma
  !
  dlon = 360.d0/dble(nlon)
  dlat = 180.d0/dble(nlat)
  !
  ! Lon, Lat, Sigma
  !
  DO i=1,nlon
    lon(i) = dlon*dble(i-1)
  END DO
  DO j=1,nlat
    lat(j) = dlat*(0.5d0 + dble(j-1))-90.d0
  END DO

  sig(1) = .95d0
  sig(2) = .835d0
  sig(3) = .685d0
  sig(4) = .51d0
  sig(5) = .34d0
  sig(6) = .2d0
  sig(7) = .08d0
  !
  ! dx and dy
  !
  dx(:) = 2.0d0 * pi * re * cos(lat(:) * pi / 180.0d0) / REAL(nlon,r_size)

  DO i=1,nlat-1
    dy(i) = 2.0d0 * pi * re * (lat(i+1) - lat(i)) / 360.0d0
  END DO
  dy(nlat) = 2.0d0 * pi * re * (90.0d0 - lat(nlat)) / 180.0d0

  DO i=2,nlat
    dy2(i) = (dy(i-1) + dy(i)) * 0.5d0
  END DO
  dy2(1) = (dy(nlat) + dy(1)) * 0.5d0
  !
  ! Corioris parameter
  !
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
  !
  ! Surface geoptential (Read Orography file)
  !
  phi0 = 0.0

  !--- prepare region infomation
  if (complete) then ! all region
    fmode = FIO_INTEG_FILE
    call MNG_mnginfo_noinput( ADM_rlevel )
  else               ! region specified by mnginfo
    fmode = FIO_SPLIT_FILE
    call MNG_mnginfo_input( ADM_rlevel, trim(mnginfo) )
  endif

  ADM_lall = 10*4**ADM_rlevel/MNG_PALL

  RETURN
END SUBROUTINE set_common_nicam
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
SUBROUTINE read_icogrd(filename1,filename2,v3d,v2d)
  USE mod_misc
  use mod_fio, only : &  ! [add] H.Yashiro 20110826
    FIO_setup,     &
    FIO_input,     &
    FIO_HMID,      &
    FIO_HSHORT,    &
    FIO_HLONG,     &
    FIO_INTEG_FILE,&
    FIO_SPLIT_FILE,&
    FIO_FREAD,      &
    FIO_ICOSAHEDRON,&
    FIO_BIG_ENDIAN, &
    headerinfo,     &
    datainfo
  use mod_mnginfo_light, only : &
    MNG_mnginfo_input,   &
    MNG_mnginfo_noinput, &
    MNG_PALL,            &
    MNG_prc_rnum,        &
    MNG_prc_tab   !(num_of_rgn,num_of_proc)
  use mod_cnst, only :       &
     CNST_UNDEF, &
     CNST_UNDEF4

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename1
  CHARACTER(*),INTENT(IN) :: filename2
  character(LEN=FIO_HLONG) :: infname   = " "
  REAL(r_size),INTENT(INOUT) :: v3d(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2d(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: l, nv, ierr, p, did
  CHARACTER(256) :: fname
  integer, allocatable :: ifid(:)
  integer, allocatable :: prc_tab_C(:)
  REAL(r_sngl), ALLOCATABLE :: data4allrgn(:)
  REAL(r_size), ALLOCATABLE :: data8allrgn(:)
  REAL(r_sngl), ALLOCATABLE :: icodata4(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: icodata(:,:,:,:)
  logical                  :: allvar = .true.

  type(headerinfo) hinfo 
  type(datainfo)   dinfo

  character(LEN=FIO_HSHORT) :: varname3d(nv3d)
  character(LEN=FIO_HSHORT) :: varname2d(nv2d)

  integer, parameter :: max_nvar   = 500
  integer, parameter :: max_nstep  = 1500
  integer, parameter :: max_nlayer = 200
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''
  integer :: v, k, kall
  integer :: istep
  logical :: addvar

  varname3d(1)='pre'
  varname3d(2)='tem'
  varname3d(3)='vx'
  varname3d(4)='vy'
  varname3d(5)='vz'
  varname3d(6)='w'
  varname3d(7)='qv'
  varname3d(8)='qc'

  varname2d( 1)='ss_ps'
  varname2d( 2)='sa_tppn'
  varname2d( 3)='ss_slp'
  varname2d( 4)='sa_tem_sfc'
  varname2d( 5)='ss_q2m'
  varname2d( 6)='ss_u10m'
  varname2d( 7)='ss_v10m'
  varname2d( 8)='ss_t2m'
  varname2d( 9)='ss_cldw'
  varname2d(10)='ss_cld_frac'
  varname2d(11)='sa_lwu_toa'
  varname2d(12)='sa_lwd_toa'
  varname2d(13)='sa_swu_toa'
  varname2d(14)='sa_swd_toa'
  varname2d(15)='sa_lwu_sfc'
  varname2d(16)='sa_lwd_sfc'
  varname2d(17)='sa_swu_sfc'
  varname2d(18)='sa_swd_sfc'
  varname2d(19)='sa_sh_sfc'
  varname2d(20)='sa_lh_sfc'
  varname2d(21)='ss_cldi'
  varname2d(22)='ll_tg'
  varname2d(23)='ll_ts'
  varname2d(24)='ll_wg'
  varname2d(25)='ol_sst'
  varname2d(26)='ol_ice'
  varname2d(27)='ol_snw'
  varname2d(28)='ol_ist'
  varname2d(29)='ol_icr'

  !restart_layername='ZSALL80'
  !restart_layername='ZSALL96'
  restart_layername='ZSALL40'
  allocate( ifid(MNG_PALL) )

  !WRITE(ADM_LOG_FID,*) 'ADM_glevel=', ADM_glevel
  !WRITE(ADM_LOG_FID,*) 'ADM_rlevel=', ADM_rlevel
  !WRITE(ADM_LOG_FID,*) 'MNG_prc_rnum=', MNG_prc_rnum
  !WRITE(ADM_LOG_FID,*) 'MNG_prc_tab=', MNG_prc_tab

  call fio_syscheck()

  write(ADM_LOG_FID,*) 'MNG_PALL=', MNG_PALL
  do p = 1, MNG_PALL
    call fio_mk_fname(infname,trim(filename1),'pe',p-1,6)

    allocate( prc_tab_C(MNG_prc_rnum(p))   )
    prc_tab_C(:) = MNG_prc_tab(:,p)-1

    call fio_put_commoninfo(  FIO_SPLIT_FILE,  &
                              FIO_BIG_ENDIAN,  &
                              FIO_ICOSAHEDRON, &
                              ADM_glevel,      &
                              ADM_rlevel,      &
                              MNG_prc_rnum(p), &
                              prc_tab_C        )
    call fio_register_file(ifid(p),trim(infname))
    call fio_fopen(ifid(p),FIO_FREAD)
    call fio_read_allinfo( ifid(p) )

    allocate( data8allrgn(ADM_gall*ADM_kall*MNG_prc_rnum(p)) )
    allocate( icodata    (ADM_gall,ADM_kall,MNG_prc_rnum(p),nv3d) )

    data8allrgn(:)  = CNST_UNDEF

    do nv = 1, nv3d
      !write(ADM_LOG_FID,*) 'READING ', p, trim(varname3d(nv))
      !flush(ADM_LOG_FID)
      call fio_seek_datainfo(did,ifid(p),varname3d(nv),1)
      call fio_get_datainfo(ifid(p),did,dinfo)
      call fio_read_data(ifid(p),did,data8allrgn(:))
      icodata(:,:,:,nv) = reshape( data8allrgn(:), shape(icodata(:,:,:,nv)) )
    end do

    do nv = 1, nv3d
      do l = 1, MNG_prc_rnum(p)
        do k = 1, ADM_kall
          v3d(:,MNG_prc_tab(l,p),k,nv)=icodata(:,k,l,nv)
        end do
      end do
    end do

    deallocate( prc_tab_C  )
    deallocate( data8allrgn  )
    deallocate( icodata  )
    call fio_fclose(ifid(p))

  end do

  TIME_CTIME=dinfo%time_start
  write(ADM_LOG_FID,*) dinfo


  write(ADM_LOG_FID,*) 'guess: pre', maxval(v3d(:,:,2:nlev-1,1)), minval(v3d(:,:,2:nlev-1,1))
  write(ADM_LOG_FID,*) 'guess: tem', maxval(v3d(:,:,2:nlev-1,2)), minval(v3d(:,:,2:nlev-1,2))
  write(ADM_LOG_FID,*) 'guess: vx ', maxval(v3d(:,:,2:nlev-1,3)), minval(v3d(:,:,2:nlev-1,3))
  write(ADM_LOG_FID,*) 'guess: vy ', maxval(v3d(:,:,2:nlev-1,4)), minval(v3d(:,:,2:nlev-1,4))
  write(ADM_LOG_FID,*) 'guess: vz ', maxval(v3d(:,:,2:nlev-1,5)), minval(v3d(:,:,2:nlev-1,5))
  write(ADM_LOG_FID,*) 'guess: w  ', maxval(v3d(:,:,2:nlev-1,6)), minval(v3d(:,:,2:nlev-1,6))
  write(ADM_LOG_FID,*) 'guess: qv ', maxval(v3d(:,:,2:nlev-1,7)), minval(v3d(:,:,2:nlev-1,7))
  write(ADM_LOG_FID,*) 'guess: qc ', maxval(v3d(:,:,2:nlev-1,8)), minval(v3d(:,:,2:nlev-1,8))
  FLUSH(ADM_LOG_FID)

  kall=1
  restart_layername='ZSSFC1'
  do p = 1, MNG_PALL
    call fio_mk_fname(infname,trim(filename2),'pe',p-1,6)

    allocate( prc_tab_C(MNG_prc_rnum(p))   )
    prc_tab_C(:) = MNG_prc_tab(:,p)-1

    call fio_put_commoninfo(  FIO_SPLIT_FILE,  &
                              FIO_BIG_ENDIAN,  &
                              FIO_ICOSAHEDRON, &
                              ADM_glevel,      &
                              ADM_rlevel,      &
                              MNG_prc_rnum(p), &
                              prc_tab_C        )
    call fio_register_file(ifid(p),trim(infname))
    call fio_fopen(ifid(p),FIO_FREAD)
    call fio_read_allinfo( ifid(p) )

    allocate( data4allrgn(ADM_gall*kall*MNG_prc_rnum(p)) )
    allocate( icodata4   (ADM_gall,kall,MNG_prc_rnum(p),nv2d) )

    data4allrgn(:)  = CNST_UNDEF4

    ! THESE VARIABLES ARE UPDATED (hourly data)
    do nv = 1, 9
      istep=6
      if(varname2d(nv)== 'sa_tppn') istep=1
      !write(ADM_LOG_FID,*) 'READING ', trim(varname2d(nv))
      !flush(ADM_LOG_FID)
      call fio_seek_datainfo(did,ifid(p),varname2d(nv),istep)
      call fio_get_datainfo(ifid(p),did,dinfo)
      call fio_read_data(ifid(p),did,data4allrgn(:))
      icodata4(:,:,:,nv) = reshape( data4allrgn(:), shape(icodata4(:,:,:,nv)) )
    end do

    ! THESE VARIABLES ARE UPDATED (6 hourly data)
    do nv = 10, nv2d
      !write(ADM_LOG_FID,*) 'READING ', trim(varname2d(nv))
      !flush(ADM_LOG_FID)
      call fio_seek_datainfo(did,ifid(p),varname2d(nv),1)
      call fio_get_datainfo(ifid(p),did,dinfo)
      call fio_read_data(ifid(p),did,data4allrgn(:))
      icodata4(:,:,:,nv) = reshape( data4allrgn(:), shape(icodata4(:,:,:,nv)) )
    end do

    do nv = 1, nv2d
      do l = 1, MNG_prc_rnum(p)
        v2d(:,MNG_prc_tab(l,p),nv)=icodata4(:,1,l,nv)
      end do
    end do

    deallocate( prc_tab_C  )
    deallocate( data4allrgn  )
    deallocate( icodata4 )
    call fio_fclose(ifid(p))
 
  end do

  deallocate(ifid)

  do nv = 1, nv3d
    write(ADM_LOG_FID,'(2A15,3(A,ES24.16))') &
            '+',      trim(varname3d(nv) ),     &
            ' max=', maxval(v3d(:,:,:,nv)),     &
            ' min=', minval(v3d(:,:,:,nv)),     &
            ' sum=',    sum(v3d(:,:,:,nv))
  end do

  do nv = 1, nv2d
    write(ADM_LOG_FID,'(2A15,3(A,ES24.16))') &
            '+    ',   trim(varname2d(nv) ), &
            ' max=', maxval(v2d(:,:,nv)),    &
            ' min=', minval(v2d(:,:,nv)),    &
            ' sum=',    sum(v2d(:,:,nv))
  end do

  flush(ADM_LOG_FID) 

END SUBROUTINE read_icogrd

SUBROUTINE read_icogrd_legacy(filename,v3d,v2d)
  USE mod_misc
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(INOUT) :: v3d(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2d(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: l, nv, k, irec
  CHARACTER(256) :: fname

!  do l = 1, ADM_rgn_nmax
!    call MISC_make_idstr(fname,trim(filename),'rgn',l)
!    open(1,file=trim(fname),form='unformatted',       &
!                  access='direct',          &
!                  recl=ADM_gall*ADM_kall*8, &
!                  status='old'              )
!    read(1,rec=1) v3d(:,l,:,1)
!    read(1,rec=2) v3d(:,l,:,2)
!    read(1,rec=3) v3d(:,l,:,3)
!    read(1,rec=4) v3d(:,l,:,4)
!    read(1,rec=5) v3d(:,l,:,5)
!    read(1,rec=6) v3d(:,l,:,6)
!    read(1,rec=7) v3d(:,l,:,7)
!    read(1,rec=8) v3d(:,l,:,8)
!    close(1)
!  end do

  do l = 1, ADM_rgn_nmax
    call MISC_make_idstr(fname,trim(filename),'rgn',l)
    open(1,file=trim(fname),form='unformatted',       &
                  access='direct',          &
                  recl=ADM_gall*8)
    irec=1
    do nv = 1, nv3d
      do k = 1, ADM_kall
        read(1,rec=irec) v3d(:,l,k,nv)
        irec=irec+1
      end do
    end do

    do nv = 1,nv2d
      read(1,rec=irec) v2d(:,l,nv)
    end do
    close(1)
  end do


  write(ADM_LOG_FID,*) 'guess: pre', maxval(v3d(:,:,2:nlev-1,1)), minval(v3d(:,:,2:nlev-1,1))
  write(ADM_LOG_FID,*) 'guess: tem', maxval(v3d(:,:,2:nlev-1,2)), minval(v3d(:,:,2:nlev-1,2))
  write(ADM_LOG_FID,*) 'guess: vx ', maxval(v3d(:,:,2:nlev-1,3)), minval(v3d(:,:,2:nlev-1,3))
  write(ADM_LOG_FID,*) 'guess: vy ', maxval(v3d(:,:,2:nlev-1,4)), minval(v3d(:,:,2:nlev-1,4))
  write(ADM_LOG_FID,*) 'guess: vz ', maxval(v3d(:,:,2:nlev-1,5)), minval(v3d(:,:,2:nlev-1,5))
  write(ADM_LOG_FID,*) 'guess: w  ', maxval(v3d(:,:,2:nlev-1,6)), minval(v3d(:,:,2:nlev-1,6))

END SUBROUTINE read_icogrd_legacy
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_icogrd(filename1,filename2,v3d,v2d)
  USE mod_misc
  use mod_fio, only : &  ! [add] H.Yashiro 20110826
    FIO_setup,     &
    FIO_input,     &
    FIO_HMID,      &
    FIO_HSHORT,    &
    FIO_HLONG,     &
    FIO_INTEG_FILE,&
    FIO_SPLIT_FILE,&
    FIO_FWRITE,     &
    FIO_OUTPUT,     &
    FIO_ICOSAHEDRON,&
    FIO_BIG_ENDIAN, &
    headerinfo,     &
    datainfo,       &
    FIO_REAL8
  use mod_mnginfo_light, only : &
    MNG_mnginfo_input,   &
    MNG_mnginfo_noinput, &
    MNG_PALL,            &
    MNG_prc_rnum,        &
    MNG_prc_tab   !(num_of_rgn,num_of_proc)
  use mod_cnst, only :       &
     CNST_UNDEF

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename1
  CHARACTER(*),INTENT(IN) :: filename2
  REAL(r_size),INTENT(IN) :: v3d(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size), ALLOCATABLE :: v3d_tmp(:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d_tmp(:,:)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec, g
  INTEGER :: l, nv
  CHARACTER(256) :: fname
  character(LEN=FIO_HLONG) :: infname   = " "

  CHARACTER(LEN=FIO_HMID) :: desc = 'INITIAL/RESTART DATA of PROGNOSTIC VARIABLES'

  integer, allocatable :: ifid(:)
  integer, allocatable :: prc_tab_C(:)
  REAL(r_size), ALLOCATABLE :: data8allrgn(:)
  REAL(r_size), ALLOCATABLE :: icodata(:,:,:,:)
  logical                  :: allvar = .true.

  type(headerinfo) hinfo

  character(LEN=FIO_HSHORT) :: varname3d(nv3d)
  character(LEN=FIO_HSHORT) :: varname2d(nv2d)
  character(LEN=FIO_HSHORT) :: unit3d(nv3d)
  character(LEN=FIO_HSHORT) :: unit2d(nv2d)

  integer, parameter :: max_nvar   = 500
  integer, parameter :: max_nstep  = 1500
  integer, parameter :: max_nlayer = 200
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''
  character(LEN=FIO_HMID) :: desc = ''
  integer :: v, p, did
  logical :: addvar

  varname3d(1)='pre'
  varname3d(2)='tem'
  varname3d(3)='vx'
  varname3d(4)='vy'
  varname3d(5)='vz'
  varname3d(6)='w'
  varname3d(7)='qv'
  varname3d(8)='qc'

  unit3d(1)='Pa'
  unit3d(2)='K'
  unit3d(3)='m/s'
  unit3d(4)='m/s'
  unit3d(5)='m/s'
  unit3d(6)='m/s'
  unit3d(7)='kg/kg'
  unit3d(8)='kg/kg'

  varname2d( 1)='ss_ps'
  varname2d( 2)='sa_tppn'
  varname2d( 3)='ss_slp'
  varname2d( 4)='sa_tem_sfc'
  varname2d( 5)='ss_q2m'
  varname2d( 6)='ss_u10m'
  varname2d( 7)='ss_v10m'
  varname2d( 8)='ss_t2m'
  varname2d( 9)='ss_cldw'
  varname2d(10)='ss_cld_frac'
  varname2d(11)='sa_lwu_toa'
  varname2d(12)='sa_lwd_toa'
  varname2d(13)='sa_swu_toa'
  varname2d(14)='sa_swd_toa'
  varname2d(15)='sa_lwu_sfc'
  varname2d(16)='sa_lwd_sfc'
  varname2d(17)='sa_swu_sfc'
  varname2d(18)='sa_swd_sfc'
  varname2d(19)='sa_sh_sfc'
  varname2d(20)='sa_lh_sfc'
  varname2d(21)='ss_cldi'
  varname2d(22)='ll_tg'
  varname2d(23)='ll_ts'
  varname2d(24)='ll_wg'
  varname2d(25)='ol_sst'
  varname2d(26)='ol_ice'
  varname2d(27)='ol_snw'
  varname2d(28)='ol_ist'
  varname2d(29)='ol_icr'

  unit2d( 1)='Pa'
  unit2d( 2)='mm/s'
  unit2d( 3)='Pa'
  unit2d( 4)='K'
  unit2d( 5)='kg/kg'
  unit2d( 6)='m/s'
  unit2d( 7)='m/s'
  unit2d( 8)='K'
  unit2d( 9)='kg/m2'
  unit2d(10)='NONE'
  unit2d(11)='J/m2'
  unit2d(12)='J/m2'
  unit2d(13)='J/m2'
  unit2d(14)='J/m2'
  unit2d(15)='J/m2'
  unit2d(16)='J/m2'
  unit2d(17)='J/m2'
  unit2d(18)='J/m2'
  unit2d(19)='J/m2'
  unit2d(20)='J/m2'
  unit2d(21)='kg/kg'
  unit2d(22)='K'
  unit2d(23)='K'
  unit2d(24)='kg'
  unit2d(25)='K'
  unit2d(26)='kg'
  unit2d(27)='kg'
  unit2d(28)='K'
  unit2d(29)='NONE'

  !restart_layername='ZSALL80'
  !restart_layername='ZSALL96'
  restart_layername='ZSALL40'
  allocate( ifid(MNG_PALL) )

  WRITE(ADM_LOG_FID,*) 'MNG_PALL', MNG_PALL
  WRITE(ADM_LOG_FID,*) 'ADM_gall', ADM_gall
  WRITE(ADM_LOG_FID,*) 'ADM_lall', ADM_lall
  WRITE(ADM_LOG_FID,*) 'ADM_kall', ADM_kall

  ALLOCATE( v3d_tmp(ADM_gall,ADM_kall,ADM_lall) )

  !--- append data to the file
  dinfo%description  = desc
  dinfo%layername    = restart_layername
  dinfo%note         = ''
  dinfo%datasize     = ADM_gall * ADM_lall * ADM_kall * 8
  dinfo%datatype     = FIO_REAL8
  dinfo%num_of_layer = ADM_kall
  dinfo%step         = 1
  dinfo%time_start   = int( TIME_CTIME, kind=8 )
  dinfo%time_end     = int( TIME_CTIME, kind=8 )
!
  desc = 'INITIAL/RESTART DATA of PROGNOSTIC VARIABLES'

  !--- output for restart
  do p = 1, MNG_PALL
    call fio_mk_fname(infname,trim(filename1),'pe',p-1,6)
    allocate( prc_tab_C(MNG_prc_rnum(p))   )
    prc_tab_C(:) = MNG_prc_tab(:,p)-1
    call fio_register_file(ifid(p),trim(infname))
    call fio_fopen(ifid(p),FIO_FWRITE)
    call fio_put_write_pkginfo(ifid(p),desc,"")
    do nv = 1, nv3d
      dinfo%varname=varname3d(nv)
      dinfo%unit   =unit3d(nv)
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,:,nv))
    end do
    deallocate( prc_tab_C  )
    call fio_fclose(ifid(p))
  end do

  !--- output for bias correction
  restart_layername='ZSSFC1'
  dinfo%layername    = restart_layername
  dinfo%datasize     = ADM_gall * ADM_lall * 8
  dinfo%num_of_layer = 1
  if(filename2 /= '') then
    do p = 1, MNG_PALL
      call fio_mk_fname(infname,trim(filename2),'pe',p-1,6)
  
      allocate( prc_tab_C(MNG_prc_rnum(p))   )
      prc_tab_C(:) = MNG_prc_tab(:,p)-1
  
      call fio_register_file(ifid(p),trim(infname))
      call fio_fopen(ifid(p),FIO_FWRITE)
      call fio_put_write_pkginfo(ifid(p),desc,"")

      !restart_layername  ='ZSDEF94'
      restart_layername  ='ZSDEF38'
      !restart_layername  ='ZSDEF78'
      dinfo%layername    = restart_layername
      dinfo%datasize     = ADM_gall * ADM_lall * (ADM_kall-2) * 8
      dinfo%num_of_layer = ADM_kall-2
 
      dinfo%varname      ='ms_pres'
      dinfo%unit         ='Pa'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,2:nlev-1,iv3d_p))

      ! THIS IS DUMMY (PRESSURE)
      dinfo%varname      ='ms_u'
      dinfo%unit         ='m/s'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,2:nlev-1,iv3d_p))

      ! THIS IS DUMMY (PRESSURE)
      dinfo%varname      ='ms_v'
      dinfo%unit         ='m/s'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,2:nlev-1,iv3d_p))

      dinfo%varname      ='ms_tem'
      dinfo%unit         ='K'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,2:nlev-1,iv3d_t))

      dinfo%varname      ='ms_qv'
      dinfo%unit         ='kg/kg'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,2:nlev-1,iv3d_qv))
      
      dinfo%varname      ='ms_qc'
      dinfo%unit         ='kg/kg'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v3d(:,p,2:nlev-1,iv3d_qc))

      restart_layername='ZSSFC1'
      dinfo%layername    = restart_layername
      dinfo%datasize     = ADM_gall * ADM_lall * 8
      dinfo%num_of_layer = 1

      dinfo%varname      ='sa_tem_sfc'
      dinfo%unit         ='K'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,4))

      dinfo%varname      ='ss_q2m'
      dinfo%unit         ='kg/kg'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,5))

      dinfo%varname      ='ss_ps'
      dinfo%unit         ='Pa'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,1))

      dinfo%varname      ='ss_u10m'
      dinfo%unit         ='m/s'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,6))

      dinfo%varname      ='ss_v10m'
      dinfo%unit         ='m/s'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,7))

      dinfo%varname      ='ss_t2m'
      dinfo%unit         ='K'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,8))

      dinfo%varname      ='ss_cldw'
      dinfo%unit         ='kg/m2'
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,9))

      do nv = 10,nv2d
        dinfo%varname    = varname2d(nv)
        dinfo%unit       = unit2d(nv)
        call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,nv))
      end do

      nv = 2
      dinfo%varname    = varname2d(nv)
      dinfo%unit       = unit2d(nv)
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,nv))

      nv = 3
      dinfo%varname    = varname2d(nv)
      dinfo%unit       = unit2d(nv)
      call fio_put_write_datainfo_data(did,ifid(p),dinfo,v2d(:,p,nv))

      deallocate( prc_tab_C  )
      call fio_fclose(ifid(p))
    end do
  end if
!
  do nv = 1, nv3d
    write(ADM_LOG_FID,'(2A15,3(A,ES24.16))') &
            '+',      trim(varname3d(nv) ),     &
            ' max=', maxval(v3d(:,:,:,nv)),     &
            ' min=', minval(v3d(:,:,:,nv)),     &
            ' sum=',    sum(v3d(:,:,:,nv))
  end do

  do nv = 1, nv2d
    write(ADM_LOG_FID,'(2A15,3(A,ES24.16))') &
            '+    ',   trim(varname2d(nv) ), &
            ' max=', maxval(v2d(:,:,nv)),    &
            ' min=', minval(v2d(:,:,nv)),    &
            ' sum=',    sum(v2d(:,:,nv))
  end do

  !write(ADM_LOG_FID,*) 'anal: pre', maxval(v3d(:,:,2:nlev-1,1)), minval(v3d(:,:,2:nlev-1,1))
  !write(ADM_LOG_FID,*) 'anal: tem', maxval(v3d(:,:,2:nlev-1,2)), minval(v3d(:,:,2:nlev-1,2))
  !write(ADM_LOG_FID,*) 'anal: vx ', maxval(v3d(:,:,2:nlev-1,3)), minval(v3d(:,:,2:nlev-1,3))
  !write(ADM_LOG_FID,*) 'anal: vy ', maxval(v3d(:,:,2:nlev-1,4)), minval(v3d(:,:,2:nlev-1,4))
  !write(ADM_LOG_FID,*) 'anal: vz ', maxval(v3d(:,:,2:nlev-1,5)), minval(v3d(:,:,2:nlev-1,5))
  !write(ADM_LOG_FID,*) 'anal: w  ', maxval(v3d(:,:,2:nlev-1,6)), minval(v3d(:,:,2:nlev-1,6))
  !write(ADM_LOG_FID,*) 'anal: ps  ', maxval(v2d(:,:,1)), minval(v2d(:,:,1))
  !write(ADM_LOG_FID,*) 'anal: prep', maxval(v2d(:,:,2)), minval(v2d(:,:,2))
  !write(ADM_LOG_FID,*) 'anal: slp ', maxval(v2d(:,:,3)), minval(v2d(:,:,3))
  !write(ADM_LOG_FID,*) 'anal: ts   ', maxval(v2d(:,:,4)), minval(v2d(:,:,4))
  !write(ADM_LOG_FID,*) 'anal: q2m  ', maxval(v2d(:,:,5)), minval(v2d(:,:,5))
  !write(ADM_LOG_FID,*) 'anal: u10m ', maxval(v2d(:,:,6)), minval(v2d(:,:,6))
  !write(ADM_LOG_FID,*) 'anal: v10m ', maxval(v2d(:,:,7)), minval(v2d(:,:,7))
  !write(ADM_LOG_FID,*) 'anal: t2m  ', maxval(v2d(:,:,8)), minval(v2d(:,:,8))
  !write(ADM_LOG_FID,*) 'anal: cldw ', maxval(v2d(:,:,9)), minval(v2d(:,:,9))

  RETURN
END SUBROUTINE write_icogrd
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_icogrd_legacy(filename,v3d,v2d)
  USE mod_misc
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  INTEGER :: l, nv
  CHARACTER(256) :: fname

  do l = 1, ADM_rgn_nmax
    call MISC_make_idstr(fname,trim(filename),'rgn',l)
    open(1,file=trim(fname),form='unformatted',       &
                  access='direct',          &
                  recl=ADM_gall*8)
                  !recl=ADM_gall*ADM_kall*8)
    irec=1
    do nv = 1, nv3d
      do k = 1, ADM_kall
        write(1,rec=irec) dble(v3d(:,l,k,nv))
        irec=irec+1
      end do
    end do

    do nv = 1,nv2d
      write(1,rec=irec) dble(v2d(:,l,nv))
    end do
    close(1)
  end do

  write(ADM_LOG_FID,*) 'anal: pre', maxval(v3d(:,:,2:nlev-1,1)), minval(v3d(:,:,2:nlev-1,1))
  write(ADM_LOG_FID,*) 'anal: tem', maxval(v3d(:,:,2:nlev-1,2)), minval(v3d(:,:,2:nlev-1,2))
  write(ADM_LOG_FID,*) 'anal: vx ', maxval(v3d(:,:,2:nlev-1,3)), minval(v3d(:,:,2:nlev-1,3))
  write(ADM_LOG_FID,*) 'anal: vy ', maxval(v3d(:,:,2:nlev-1,4)), minval(v3d(:,:,2:nlev-1,4))
  write(ADM_LOG_FID,*) 'anal: vz ', maxval(v3d(:,:,2:nlev-1,5)), minval(v3d(:,:,2:nlev-1,5))
  write(ADM_LOG_FID,*) 'anal: w  ', maxval(v3d(:,:,2:nlev-1,6)), minval(v3d(:,:,2:nlev-1,6))

  RETURN
END SUBROUTINE write_icogrd_legacy

SUBROUTINE write_icogrd4(filename,v3d,v2d)
  USE mod_misc
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  INTEGER :: l, nv
  CHARACTER(256) :: fname

  do l = 1, ADM_rgn_nmax
    call MISC_make_idstr(fname,trim(filename),'rgn',l)
    open(1,file=trim(fname),form='unformatted',       &
                  access='direct',          &
                  recl=ADM_gall*ADM_kall*4)
    do nv = 1, nv3d
      write(1,rec=nv) v3d(:,l,:,nv)
    end do
    close(1)
  end do

  RETURN
END SUBROUTINE write_icogrd4
!-----------------------------------------------------------------------
! p_full => ensemble mean
!-----------------------------------------------------------------------
SUBROUTINE calc_pfull(ix,jy,ps,p_full)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ix,jy
  REAL(r_size),INTENT(IN) :: ps(ix,jy)
  REAL(r_size),INTENT(OUT) :: p_full(ix,jy,nlev)
  INTEGER :: i,j,k

!!$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=1,nlev
    DO j=1,jy
      DO i=1,ix
        p_full(i,j,k) = ps(i,j) * sig(k)
      END DO
    END DO
  END DO
!!$OMP END PARALLEL DO

  RETURN
END SUBROUTINE calc_pfull
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(ADM_LOG_FID,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(ADM_LOG_FID,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    WRITE(ADM_LOG_FID,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

  RETURN
END SUBROUTINE monit_grd
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

  DO n=1,nv3d
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd

END MODULE common_nicam
