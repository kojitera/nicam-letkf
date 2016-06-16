!-------------------------------------------------------------------------------
!
!+  Program obsope
!
!-------------------------------------------------------------------------------
program prg_vbc

  use mpi
  use mod_misc, only :         &
       MISC_make_idstr,        &
       MISC_get_available_fid
  use mod_adm, only  :         &
       ADM_CTL_FID,            &
       ADM_LOG_FID,            &
       ADM_MAXFNAME,           &
       ADM_NSYS,               &
       ADM_setup,              &
       ADM_MULTI_PRC,          &
       ADM_proc_init,          &
       ADM_proc_stop,          &
       ADM_proc_finish,        &
       ADM_prc_me,             &
       ADM_prc_run_master,     &
       ADM_prc_run_all,        &
       ADM_prc_tab,            &
       ADM_gall,               &
       ADM_KNONE,              &
       ADM_lall,               &
       ADM_GALL_PL,            &
       ADM_LALL_PL,            &
       ADM_gall_in,            &
       ADM_vlayer,             &
       ADM_prc_all
  use mod_cnst, only :         &
       CNST_setup,             &
       CNST_PI,                &
       CNST_ERADIUS,           &
       CNST_UNDEF,             &
       CNST_UNDEF4
  use mod_comm, only :         &
       COMM_setup
  use mod_grd, only  :         &
       GRD_setup
  use mod_gmtr, only :         &
       GMTR_setup,             &
       GMTR_P_var,             &
       GMTR_P_LAT,             &
       GMTR_P_LON
  use mod_oprt, only :         &
       OPRT_setup
  use mod_vmtr, only :         &
       VMTR_setup
  use mod_gtl, only :          &
       GTL_clip_region,        &
       GTL_clip_region_1layer, &
       GTL_input_var2
  use mod_fio, only : &
    FIO_HSHORT,       &
    FIO_HMID,         &
    FIO_HLONG,        &
    FIO_REAL4,        &
    FIO_REAL8,        &
    FIO_BIG_ENDIAN,   &
    FIO_ICOSAHEDRON,  &
    FIO_IGA_LCP,      &
    FIO_IGA_MLCP,     &
    FIO_INTEG_FILE,   &
    FIO_SPLIT_FILE,   &
!    FIO_MPIIO_NOUSE,  & 
    FIO_FREAD,        &
    headerinfo,       &
    datainfo
  use mod_calendar, only : &
    calendar_ss2yh,   &
    calendar_yh2ss
  use mod_mnginfo_light, only : &
    MNG_mnginfo_input,   &
    MNG_mnginfo_noinput, &
    MNG_PALL,            &
    MNG_prc_rnum,        &
    MNG_prc_tab   !(num_of_rgn,num_of_proc)  
  use common_tvs_nicam, only : &
       set_instrument, &
       tvsname,        &
       tvsinst,        &
       tvsch,          &
       ntvsch,         &
       ntvschan,       &
       nfootp

  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  REAL(8),PARAMETER :: q2ppmv=1.60771704d6

  integer, parameter :: max_nvar   = 500
  integer, parameter :: max_nstep  = 1500
  integer, parameter :: max_nlayer = 200

  integer, parameter :: flim = 1
  integer,      save :: fmax

  real(8), parameter :: UNDEF = -999.9d0
  integer, parameter :: limit_itrmax = 100000
  !
  !+ namelist variables
  character(ADM_MAXFNAME), save :: input_fname = ''
  character(ADM_MAXFNAME), save :: pres_basename = ''
  real(8), save :: alpha = 0.0d0
  real(8), save :: dalpha = 0.0d0
  character(ADM_NSYS), save :: sfc_type = 'RIGID'
  !                                 'ALMOST'
  !                                 'FREE'
  logical, save :: h_intrpl = .false.
  integer, save :: search_itrmax = 3
  character(ADM_MAXFNAME), save :: output_basename = ''
  character(ADM_MAXFNAME), save :: veg_base = ''
  real(8), save :: ocean_value = 0.0d0
  logical, save :: out_gsfile = .false.
  ! :: idate(1:6) year, month, day, hour, min, sec
  integer :: idate(1:6)
  real(8) :: time_tmp
  integer(8), allocatable :: time_obs(:)
  character(ADM_MAXFNAME) :: pres_fname

  !--- NAMELIST
  integer                   :: glevel              = -1
  integer                   :: rlevel              = -1
  character(LEN=FIO_HSHORT) :: grid_topology       = 'ICOSAHEDRON'
                                                   ! 'LCP'
                                                   ! 'MLCP'
  logical                   :: complete            = .false.
  character(LEN=FIO_HLONG)  :: mnginfo             = ''
  character(LEN=FIO_HLONG)  :: layerfile_dir       = ''
  character(LEN=FIO_HLONG)  :: llmap_base          = ''
  character(LEN=FIO_HLONG)  :: infile(flim)        = ''
  integer                   :: step_str            = 1
  integer                   :: step_end            = max_nstep
  character(LEN=FIO_HLONG)  :: outfile_dir         = '.'
  character(LEN=FIO_HLONG) :: outfile_prefix      = 'sim_obs'
  !character(LEN=FIO_HSHORT) :: outfile_prefix      = 'sim_obs'
  integer                   :: outfile_rec         = 1
  logical                   :: output_grads        = .true.
  logical                   :: datainfo_nodep_pe   = .false.   !   <- can be .true. if data header do not depend on pe.
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''
  character(LEN=FIO_HSHORT) :: large_memory_var(max_nvar) = '' ! [add] 13-04-18
  logical                   :: help = .false.
  character(LEN=512) :: obsfile

  character(LEN=FIO_HLONG) :: infname   = " "
  character(LEN=FIO_HLONG) :: outbase   = " "
  character(LEN=FIO_HLONG) :: outbase2  = " "
  character(LEN=FIO_HLONG) :: layerfile = " "
  integer                  :: fmode
  integer                  :: gtopology
  logical                  :: allvar = .true.

  integer :: nsite
  integer, allocatable, save :: obsnum(:)      ! number
  real(8), allocatable, save :: lat(:)         ! latitude
  real(8), allocatable, save :: lon(:)         ! longitude
  !real(8), allocatable, save :: elev(:)        ! pressure [hPa]
  !integer, allocatable, save :: elev(:)        ! elevation [m]
  integer, allocatable, save :: agl(:)         ! above ground level = 1, if not = 0
  !integer, allocatable, save :: land(:)        ! land = 1, ocean = 0
  integer, allocatable, save :: idirec(:,:)    ! 1: primary wind direction
  !                                            ! 2: secondary wind direction
  character(7), parameter :: cdirec_ref(0:8) = (/ &
       'ANY    ',   & !                             
       'N      ',   & !                      1      0: ANY direction
       'NE     ',   & !                  8   N   2
       'E      ',   & !                   NW | NE
       'SE     ',   & !                7 W ----- E 3
       'S      ',   & !                   SW | SE
       'SW     ',   & !                  6   S   4
       'W      ',   & !                      5 
       'NW     '    &
       /)

  character(7), allocatable :: cdirec(:)
  character(7) :: cdrc0
  character(7) :: cdrc(2)
  !
  integer, allocatable, save :: l_index(:)
  integer, allocatable, save :: n1_index(:)
  integer, allocatable, save :: n2_index(:)
  integer, allocatable, save :: n3_index(:)
  real(8), allocatable, save :: w1(:)
  real(8), allocatable, save :: w2(:)
  real(8), allocatable, save :: w3(:)
  integer, allocatable, save :: klev(:,:)
  real(8), allocatable, save :: kfact(:)
  integer, allocatable, save :: date(:,:)
  character(len=128), allocatable, save :: var(:)
  real(8), allocatable, save :: obserr(:,:)
  real(8), allocatable, save :: odat(:,:)
  integer, allocatable, save :: elem(:)
  integer, allocatable, save :: id_obs(:)
  integer, allocatable, save :: oqc(:,:)
  integer, allocatable, save :: oqc_out(:,:)


  real(8), allocatable :: lon_obs(:)
  real(8), allocatable :: lat_obs(:)
  real(8), allocatable :: lev_obs(:)

  INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300) !Standard real type
  !real(kind=JPRB), allocatable :: saza(:), saaz(:), soza(:), soaz(:)
  !real(kind=JPRB), allocatable :: land(:), elev(:)
  !real(kind=JPRB), allocatable :: said(:)
  real(8), allocatable :: saza(:), saaz(:), soza(:), soaz(:)
  real(8), allocatable :: land(:), elev(:)
  real(8), allocatable :: said(:)
  real(8), allocatable, save :: lsql(:,:,:)


  !real(8), allocatable, save :: saza(:)
  !real(8), allocatable, save :: saaz(:)
  !real(8), allocatable, save :: soza(:)
  !real(8), allocatable, save :: soaz(:)


  logical, allocatable, save :: inprc(:)
  real(8), allocatable, save :: veg(:,:,:)
  real(8), allocatable, save :: veg_pl(:,:,:)
  real(8), allocatable, save :: veg_in(:,:)
  integer, allocatable, save :: icoland(:,:)
  real(8), allocatable, save :: icolat(:,:)
  real(8), allocatable, save :: icolon(:,:)
  !
  real(8), allocatable, save :: lat_model(:)
  real(8), allocatable, save :: lon_model(:)
  real(4), allocatable, save :: land_model(:)
  integer, allocatable, save :: idirec_model(:)
  !
  integer, save :: max_num_latlon
  integer :: sum_max_num_latlon
  !
  real(8), allocatable, save :: var0(:,:)
  real(8), allocatable, save :: var1(:,:)
  real(8) :: tmp, tmp2

  integer :: s, n, l, i, nn, j, k
  integer :: n1, n2, n3
  integer :: prc
  integer :: fid, fid2
  integer :: ierr
  integer :: rgnid
  character(5) :: prcnum
  character(5) :: crgn
  character(ADM_MAXFNAME) :: fname2
  !
  INTEGER,PARAMETER :: id_temp_nicam=1
  INTEGER,PARAMETER :: id_qvap_nicam=2
  INTEGER,PARAMETER :: id_qcld_nicam=3
  INTEGER,PARAMETER :: id_pres_nicam=4
  INTEGER,PARAMETER :: id_tsfc_nicam=5
  !INTEGER,PARAMETER :: id_cldf_nicam=6
  INTEGER,PARAMETER :: id_qv2m_nicam=6
  INTEGER,PARAMETER :: id_surp_nicam=7
  INTEGER,PARAMETER :: id_u10m_nicam=8
  INTEGER,PARAMETER :: id_v10m_nicam=9
  INTEGER,PARAMETER :: id_te2m_nicam=10
  INTEGER,PARAMETER :: id_cldw_nicam=11
  !
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_z_obs=2567
  INTEGER,PARAMETER :: id_s_obs=3332
  INTEGER,PARAMETER :: id_rain_obs=9999
  INTEGER,PARAMETER :: id_bt_obs=21023

  INTEGER,PARAMETER :: id_u=1
  INTEGER,PARAMETER :: id_v=2
  INTEGER,PARAMETER :: id_t=3
  INTEGER,PARAMETER :: id_q=4
  INTEGER,PARAMETER :: id_ps=5
  INTEGER,PARAMETER :: id_rain=6

  INTEGER :: id_satellite
  INTEGER :: num_satellite
  INTEGER :: nobs_tmp1
  INTEGER :: nobs_tmp2
  INTEGER :: nobs_noaa15
  INTEGER :: nobs_noaa16
  INTEGER :: nobs_noaa17
  INTEGER :: nobs_noaa18
  INTEGER :: nobs_noaa19
  INTEGER :: nobs_metop2
  INTEGER,PARAMETER :: id_NOAA15=206
  INTEGER,PARAMETER :: id_NOAA16=207
  INTEGER,PARAMETER :: id_NOAA17=208
  INTEGER,PARAMETER :: id_NOAA18=209
  INTEGER,PARAMETER :: id_NOAA19=223
  INTEGER,PARAMETER :: id_METOP2=4

  INTEGER,PARAMETER :: nslots=1
  INTEGER,PARAMETER :: islot=1
  INTEGER,PARAMETER :: ninstrument=5
  INTEGER,PARAMETER :: maxtvsch=15
  INTEGER,PARAMETER :: maxvbc=8
  INTEGER,PARAMETER :: maxfoot=30
  INTEGER,SAVE :: maxtvsprof
  INTEGER,SAVE :: maxtvsfoot
  INTEGER,SAVE :: ntvsprof(ninstrument)
  INTEGER,PARAMETER :: nlev=38

  !CHARACTER(4),SAVE :: tvsname(ninstrument)
  CHARACTER(256),SAVE :: rttovcoef_fname(ninstrument)
  !INTEGER,SAVE :: tvsinst(3,ninstrument)
  !INTEGER,SAVE :: tvsch(maxtvsch,ninstrument)
  !INTEGER,SAVE :: ntvsch(ninstrument)
  INTEGER,SAVE :: ntvs(ninstrument)
  REAL(8),ALLOCATABLE,SAVE :: vbc_pred(:,:,:,:,:)
  REAL(8),SAVE :: vbc(maxvbc,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbca(maxvbc,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbcf_scan(maxfoot,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbca_scan(maxfoot,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbcf_scan_tmp(maxfoot,maxtvsch,ninstrument)
  CHARACTER(4),SAVE :: tvsname_scan

  REAL(8),ALLOCATABLE,SAVE :: tvslon(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvslat(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvselev(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsdat(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsdep(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsqc(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvserr(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsfoot(:,:,:)

  REAL(Kind=8), allocatable :: bt(:,:,:,:)
  !REAL(Kind=8), allocatable :: bt(:,:)
  REAL(Kind=8), allocatable :: bt_tmp(:,:,:,:)
  !REAL(Kind=8), allocatable :: bt_tmp(:,:)
  REAL(Kind=8), allocatable :: tran(:,:,:,:,:)
  REAL(Kind=8), allocatable :: weight(:,:,:,:)
  INTEGER, allocatable :: weight_maxlev(:,:,:)
  !REAL(Kind=8), allocatable :: weight_maxlev(:,:)
  INTEGER :: ichan, ilev, iobs, iobs1, ifoot, ic
  REAL(Kind=8) :: tmp_lev
  REAL(Kind=8) :: iwlr

  ! ico data
  integer              :: PALL_global
  integer              :: LALL_global
  integer              :: LALL_local

  ! for MPI
  integer          :: prc_nall, prc_nlocal
  integer          :: pstr, pend
  integer          :: prc_myrank

  ! ico data information
  integer, allocatable :: ifid(:)
  integer, allocatable :: prc_tab_C(:)
  type(headerinfo) hinfo
  type(datainfo)   dinfo

  integer                                :: num_of_data
  integer                                :: nvar
  character(LEN=FIO_HSHORT), allocatable :: var_name(:)
  character(LEN=FIO_HMID),   allocatable :: var_desc(:)
  character(LEN=FIO_HSHORT), allocatable :: var_unit(:)
  character(LEN=FIO_HSHORT), allocatable :: var_layername(:)
  integer,                   allocatable :: var_datatype(:)
  integer,                   allocatable :: var_nlayer(:)
  integer,                   allocatable :: var_nstep(:)
  integer(8),                allocatable :: var_time_str(:)
  integer(8),                allocatable :: var_dt(:)
  real(8),                   allocatable :: var_zgrid(:,:)
  ! header
  character(LEN=16),         allocatable :: var_gthead(:,:)

  character(LEN=FIO_HLONG) :: fname
  character(LEN=20)        :: tmpl
  character(LEN=16)        :: gthead(64)
  integer(8)               :: nowsec
  integer(8)               :: recsize ! [mod] 12-04-19 H.Yashiro
  integer                  :: kmax, num_of_step, step, date_str(6)

  logical :: addvar
  integer :: did, ofid, ofid2, irec
  integer :: t, v, p, pp  ! loop-index
  real(8) :: fac1, fac2, fac3, fac4, fac5, fac6, fac_sum
  integer :: ks, ke

  ! ico data
  integer              :: GALL
  integer              :: LALL
  real(4), allocatable :: data4allrgn(:)
  real(8), allocatable :: data8allrgn(:)
  real(4), allocatable :: icodata4(:,:,:)
  real(4), allocatable :: icodata(:,:,:,:)

  real(4) :: wk(34)
  real(8), allocatable :: pres(:,:,:)
  real(8), allocatable :: lnpres(:,:,:)
  real(8), allocatable :: tem(:,:,:)
  real(8), allocatable :: qv(:,:,:)
  real(8), allocatable :: ps(:,:)
  real(8), allocatable :: lnps(:,:)
  real(8), allocatable :: ts(:,:)
  real(8), allocatable :: t2m(:,:)
  real(8), allocatable :: qs(:,:)
  real(8), allocatable :: ms_u(:,:,:)
  real(8), allocatable :: ms_v(:,:,:)
  real(8), allocatable :: ms_t(:,:,:)
  real(8), allocatable :: ms_qv(:,:,:)
  !real(kind=JPRB), allocatable :: obsdata(:,:,:)
  !real(kind=JPRB), allocatable :: obsdata_out(:,:,:)
  real(8), allocatable :: obsdata(:,:,:)
  real(8), allocatable :: obsdata_out(:,:,:)
  !real(kind=JPRB), allocatable :: obsdata_jprb(:,:,:)
  real(kind=JPRB), allocatable :: obsdata_jprb(:,:,:,:)
  !real(4), allocatable :: obsdata(:)
  !real(4), allocatable :: obsdata_out(:)
  integer :: sobs, eobs

  LOGICAL,PARAMETER :: msw_vbc = .TRUE.

  logical :: ocheck=.false.
  real(8) :: tmp_time(10)

  character(4) :: cfile='inst'

  integer,parameter:: &
     &    rttv_plat_noaa =1 ,rttv_plat_dmsp =2 &
     &   ,rttv_plat_meteo=3 ,rttv_plat_goes =4 &
     &   ,rttv_plat_gms  =5 ,rttv_plat_fy2  =6 &
     &   ,rttv_plat_trmm =7 ,rttv_plat_ers  =8 &
     &   ,rttv_plat_eos  =9 ,rttv_plat_metop=10 &
     &   ,rttv_plat_envi =11,rttv_plat_msg  =12 &
     &   ,rttv_plat_fy1  =13,rttv_plat_adeos=14 &
     &   ,rttv_plat_mtsat=15,rttv_plat_cori =16 &
     &   ,rttv_plat_npoes=17,rttv_plat_gifts=18 &
     &   ,rttv_plat_metop2=19

  integer,parameter:: &
     &    rttv_inst_hirs  =0, rttv_chan_hirs  =19   & ! HIRS
     &   ,rttv_inst_amsua =3, rttv_chan_amsua =15   & ! AMSU-A
     &   ,rttv_inst_amsub =4, rttv_chan_amsub =5    & ! AMSU-B
     &   ,rttv_inst_ssmi  =6, rttv_chan_ssmi  =7    & ! SSMI
     &   ,rttv_inst_tmi   =9, rttv_chan_tmi   =9    & ! TMI
     &   ,rttv_inst_ssmis =10,rttv_chan_ssmis =24   & ! SSMIS
     &   ,rttv_inst_airs  =11,rttv_chan_airs  =2378 & ! AIRS
     &   ,rttv_inst_amsr  =17,rttv_chan_amsr  =14   & ! AMSR
     &   ,rttv_inst_mviri =20,rttv_chan_mviri =2    & ! METEOSAT
     &   ,rttv_inst_seviri=21,rttv_chan_seviri=8    & ! MSG
     &   ,rttv_inst_goesi =22,rttv_chan_goesi =4    & ! GOES-IMAGER(IR)
     &   ,rttv_inst_mtsati=24,rttv_chan_mtsati=4      ! MTSAT

  namelist /obsope_param/ &
       input_fname,        &
       alpha,              &
       dalpha,             &
       sfc_type,           &
       h_intrpl,           &
       search_itrmax,      &
       output_basename,    &
       veg_base,           &
       ocean_value,        &
       idate,              &
       out_gsfile,         &
       rttovcoef_fname,    &
       pres_basename,      &
       ocheck

  namelist /OPTION/ glevel,            &
                    rlevel,            &
                    grid_topology,     &
                    complete,          &
                    mnginfo,           &
                    layerfile_dir,     &
                    llmap_base,        &
                    infile,            &
                    step_str,          &
                    step_end,          &
                    outfile_dir,       &
                    outfile_prefix,    &
                    outfile_rec,       &
                    output_grads,      &
                    datainfo_nodep_pe, &  ! [add] 13-04-18
                    selectvar,         &
                    large_memory_var,  &  ! [add] 13-04-18
                    obsfile,           &  ! [add] 13-06-06
                    help


  tmp_time(1)=MPI_WTIME()

  call ADM_proc_init(ADM_MULTI_PRC)
  call ADM_setup('vbc.cnf')
  call COMM_setup
  call CNST_setup
  call GRD_setup
  call GMTR_setup
  call OPRT_setup
  call VMTR_setup
  call readoption !! set fmax, infile

  call set_instrument

  rewind(ADM_CTL_FID)
  read(ADM_CTL_FID, nml=obsope_param, iostat=ierr)
  !rttovcoef_fname(1)='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_15_amsua.dat'
  !rttovcoef_fname(2)='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_16_amsua.dat'
  !rttovcoef_fname(3)='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_18_amsua.dat'
  !rttovcoef_fname(4)='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_19_amsua.dat'

  write(ADM_LOG_FID,*) trim(rttovcoef_fname(1))
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(2))
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(3))
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(4))
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(5))


  fid=40
  open(fid, file=trim(input_fname), form='unformatted', access='sequential', &
       status='old', iostat=ierr)
  if(ierr /= 0) then
     write(ADM_LOG_FID,*) 'ERROR for opening the file', trim(input_fname)
     write(ADM_LOG_FID,*) 'STOP!'
     call ADM_proc_stop
  end if

  read(fid) num_satellite
  read(fid) nsite
  read(fid) ntvs(1:ninstrument)
  !read(fid) nobs_noaa15, nobs_noaa16, nobs_noaa18, nobs_noaa19, nobs_metop2
  maxtvsprof=maxval(ntvs)
  maxtvsfoot=maxval(nfootp)
  nsite=sum(ntvs)
  write(ADM_LOG_FID,*) 'num_satellite=', num_satellite
  write(ADM_LOG_FID,*) 'nsite=        ', nsite
  write(ADM_LOG_FID,*) 'ntvs=         ', ntvs(1:ninstrument)
  write(ADM_LOG_FID,*) 'maxtvsprof=   ', maxtvsprof
  write(ADM_LOG_FID,*) 'maxtvsfoot=   ', maxtvsfoot
  FLUSH(ADM_LOG_FID)
  !write(ADM_LOG_FID,*) nobs_noaa15, nobs_noaa16, nobs_noaa18, nobs_noaa19, nobs_metop2

  !call get_nsite(fid,nsite) ! number of observation

  allocate( tvslat(maxtvsprof,ninstrument,nslots) )
  allocate( tvslon(maxtvsprof,ninstrument,nslots) )
  allocate( tvselev(maxtvsprof,ninstrument,nslots) )
  allocate( tvsdat(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvserr(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvsdep(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvsqc(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvsfoot(maxtvsprof,ninstrument,nslots) )

  allocate(veg(ADM_gall,ADM_KNONE,ADM_lall))
  allocate(veg_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL))
  allocate(icoland(ADM_gall,ADM_lall))
  allocate(icolat(ADM_gall,ADM_lall))
  allocate(icolon(ADM_gall,ADM_lall))

  allocate(pres(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(lnpres(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(tem(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(qv(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ps(ADM_gall,ADM_lall))
  allocate(lnps(ADM_gall,ADM_lall))
  allocate(ts(ADM_gall,ADM_lall))
  allocate(t2m(ADM_gall,ADM_lall))
  allocate(qs(ADM_gall,ADM_lall))
  allocate(ms_u(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ms_v(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ms_t(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ms_qv(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(icodata(ADM_gall,ADM_vlayer,ADM_lall,10))
  allocate(obsdata(nlev,nsite,11))
  allocate(obsdata_out(nlev,nsite,11))
  allocate(obsdata_jprb(nlev,maxtvsprof,ninstrument,11))
 ! allocate(obsdata_out(nsite))

  allocate(obsnum(nsite))
  allocate(lon(nsite))
  allocate(lat(nsite))
  allocate(elev(nsite))
  allocate(agl(nsite))
  allocate(land(nsite))
  allocate(cdirec(nsite))
  allocate(idirec(2,nsite))
  allocate(date(6,nsite))
  allocate(var(nsite))
  allocate(obserr(15,nsite))
  allocate(elem(nsite))
  allocate(id_obs(nsite))
  allocate(time_obs(nsite))
  idirec(:,:) = 0
  !
  allocate(l_index(nsite))
  allocate(n1_index(nsite))
  allocate(n2_index(nsite))
  allocate(n3_index(nsite))
  allocate(w1(nsite))
  allocate(w2(nsite))
  allocate(w3(nsite))
  allocate(klev(2,nsite))
  allocate(kfact(nsite))
  allocate(odat(15,nsite))
  allocate(oqc(15,nsite))
  allocate(oqc_out(15,nsite))
  allocate( lon_obs(nsite) )
  allocate( lat_obs(nsite) )
  allocate( lev_obs(nsite) )

  allocate( lsql(maxtvsprof,ninstrument,nslots) )
  allocate( saza(nsite) )
  allocate( saaz(nsite) )
  allocate( soza(nsite) )
  allocate( soaz(nsite) )
  allocate( said(nsite) )

  !
  allocate(inprc(nsite))
  inprc = .false.
  !
  allocate(lat_model(nsite))
  allocate(lon_model(nsite))
  allocate(land_model(nsite))
  allocate(idirec_model(nsite))
  !
  icolat(:,:) = GMTR_P_var(:,ADM_KNONE,:,GMTR_P_LAT)
  icolon(:,:) = GMTR_P_var(:,ADM_KNONE,:,GMTR_P_LON)
  l_index(:) = 0
  n1_index(:) = -1
  n2_index(:) = -1
  n3_index(:) = -1
  w1(:) = 0.0d0
  w2(:) = 0.0d0
  w3(:) = 0.0d0
  klev(:,:) = -1
  kfact(:) = 0.0d0
  land_model(:) = -1.0
  idirec_model(:) = -1
  elem(:)=id_bt_obs
  obsdata(:,:,:)=0.0

  !call set_instrument

  call calendar_yh2ss( time_tmp, idate )
  time_obs(:)=dnint(time_tmp)

  tmp_time(2)=MPI_WTIME()

  !rewind(fid)
  s=1
  do nn = 1, ninstrument
    do n = 1, ntvs(nn)
    !do s=1, nsite
      read(fid) wk
      lat(s) = dble(wk(7))
      tvslat(n,nn,1)  = dble(wk(7))
      lon(s) = dble(wk(8))
      tvslon(n,nn,1)  = dble(wk(8))
      said(s) = dble(wk(9))
      odat(1:15,s) = dble(wk(19:33)) !! 2014.07.17 fix
      tvsfoot(n,nn,1) = dble(wk(11))
      lsql(n,nn,1) = dble(wk(12))
      saza(s) = dble(wk(13))
      soza(s) = dble(wk(14))
      elev(s) = dble(wk(15))
      tvselev(n,nn,1) = dble(wk(15))
      soaz(s) = dble(wk(17))
      saaz(s) = dble(wk(18))
      !obserr(:,s) = 0.3
      !obserr(:,s) = 0.5
      s=s+1
    enddo
  enddo

  tmp_time(3)=MPI_WTIME()

  do s=1, nsite
    if( lat(s) > 90.0 .or. lat(s) < -90.0 ) then
      write(ADM_LOG_FID,*) s, lat(s), lon(s)
    end if
    if( lon(s) > 180.0 .or. lat(s) < -180.0 ) then
      write(ADM_LOG_FID,*) s, lat(s), lon(s)
    end if
  enddo

  lon(:) = lon(:)*CNST_PI/180.0d0
  lat(:) = lat(:)*CNST_PI/180.0d0

  tvsqc(:,:,:,:)=0
  i=1
  write(ADM_LOG_FID,*) 'tvsdat'
  do nn = 1, ninstrument
    do n = 1, ntvs(nn)
      tvsdat(:,n,nn,1)=odat(:,i)
      i=i+1
    end do
  end do
  tvserr(:,:,:,:)=0.5d0
  !tvserr(:,:,:,:)=0.4d0
  !tvserr(:,:,:,:)=0.35d0

  do nn = 1, ninstrument
    do ic = 1, ntvsch(nn)
      ichan=tvsch(ic,nn)
      tvsqc(ichan,1:ntvs(nn),nn,1)=1
    end do
  end do

  do nn = 1, ninstrument
    write(ADM_LOG_FID,'(i5,10f10.3)') nn, (tvsdat(tvsch(ic,nn),1,nn,1),ic=1,ntvsch(nn))
  end do

!  do nn = 1, ninstrument
!    tvserr(1,:,nn,:)=0.3
!    tvserr(2,:,nn,:)=0.2
!    tvserr(3,:,nn,:)=0.25
!    tvserr(4,:,nn,:)=0.28
!    tvserr(5,:,nn,:)=0.3
!    tvserr(6,:,nn,:)=0.4
!  end do

  call GTL_input_var2(trim(veg_base), veg, veg_pl, ADM_KNONE, ADM_KNONE, 8)

  icoland(:,:) = 1
  do l=1, ADM_lall
     do n=1, ADM_gall
        if(veg(n,ADM_KNONE,l) == ocean_value) icoland(n,l) = 0
     end do
  end do
  !
  !--- get the ico. grids corresponding to the lat-lon grids.
  if(.not. h_intrpl) then
     call geticogrid
     w1(:) = 1.0d0
  else
     call setup_ico2latlon_mapping('GET_NUM')
     !
     call MPI_ALLREDUCE(max_num_latlon, sum_max_num_latlon, 1, MPI_INTEGER, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
     write(ADM_LOG_FID,*) 'max_num_latlon: ', ADM_prc_me, max_num_latlon
     !
     if(sum_max_num_latlon /= nsite) then
        call ADM_proc_stop
     end if
     call setup_ico2latlon_mapping('SET_INDEX')
     !
     idirec_model(:) = 0
     !
     do s=1, nsite
        if(l_index(s) /= 0) then
           inprc(s) = .true.
        end if
        !if(elev(s) < 0.1 ) then
        !   inprc(s) = .false. 
        !end if
     end do
     !
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  tmp_time(4)=MPI_WTIME()

  if ( glevel==-1 .or. rlevel==-1 ) then
     write(ADM_LOG_FID,*) "xxx Set glevel, rlevel. STOP"
     stop
  endif
  if ( step_str < 1 .or. step_end < 1 ) then
     write(ADM_LOG_FID,*) "xxx step must be >= 1. STOP"
     stop
  elseif( step_str > step_end ) then
     write(ADM_LOG_FID,*) "xxx step_str must be < step_end. STOP"
     stop
  endif

  if ( grid_topology=="ICOSAHEDRON" ) then
     gtopology = FIO_ICOSAHEDRON
  elseif( grid_topology=="LCP" ) then
     gtopology = FIO_IGA_LCP
  elseif( grid_topology=="MLCP" ) then
     gtopology = FIO_IGA_MLCP
  else
     write(ADM_LOG_FID,*) "Unknown type of Grid toporogy:",grid_topology
     stop
  endif

  if ( trim(selectvar(1)) /= '' ) then
     allvar = .false.
  endif

  !--- prepare region infomation
  if (complete) then ! all region
    fmode = FIO_INTEG_FILE
    call MNG_mnginfo_noinput( rlevel )
  else               ! region specified by mnginfo
    fmode = FIO_SPLIT_FILE
    call MNG_mnginfo_input( rlevel, trim(mnginfo) )
  endif

  call fio_syscheck()

  prc_myrank  = ADM_prc_me - 1
  PALL_global = MNG_PALL
  LALL_global = 10 * (4**rlevel)
  LALL_local  = LALL_global / PALL_global

  if ( mod( PALL_global, ADM_prc_all) /= 0 ) then
     write(ADM_LOG_FID,*) "*** Invalid processor number, STOP:", PALL_global, ADM_prc_all
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
     stop
  endif

  prc_nlocal = PALL_global / ADM_prc_all
  pstr       = prc_myrank*prc_nlocal + 1
  pend       = prc_myrank*prc_nlocal + prc_nlocal

  allocate( ifid(MNG_PALL) )
  allocate( var_nstep    (max_nvar) )
  allocate( var_name     (max_nvar) )
  allocate( var_desc     (max_nvar) )
  allocate( var_unit     (max_nvar) )
  allocate( var_layername(max_nvar) )
  allocate( var_datatype (max_nvar) )
  allocate( var_nlayer   (max_nvar) )
  allocate( var_time_str (max_nvar) )
  allocate( var_dt       (max_nvar) )
  allocate( var_zgrid    (max_nlayer, max_nvar) )
  allocate( var_gthead   (64, max_nvar) )

  do v = 1, max_nvar
     if ( trim(large_memory_var(v)) == '' ) exit
     call fio_register_vname_tmpdata( trim(large_memory_var(v)), &
          len_trim(large_memory_var(v)) )
  end do

  do p = pstr, pend
     pp = p - pstr + 1

     if (complete) then ! all region
        infname = trim(infile(1))//'.rgnall'
     else
        call fio_mk_fname(infname,trim(infile(1)),'pe',p-1,6)
     endif
     allocate( prc_tab_C(MNG_prc_rnum(p))   )
     prc_tab_C(:) = MNG_prc_tab(:,p)-1
     if ( pp == 1 ) then
       call fio_put_commoninfo( fmode,           &
                                FIO_BIG_ENDIAN,  &
                                gtopology,       &
                                glevel,          &
                                rlevel,          &
                                MNG_prc_rnum(p), &
                                prc_tab_C        )
     end if

     call fio_register_file(ifid(pp),trim(infname))
     call fio_fopen(ifid(pp),FIO_FREAD)

     ! <-- [add] C.Kodama 13.04.18
     if( datainfo_nodep_pe .and. pp > 1 ) then
        ! assume that datainfo do not depend on pe.
        call fio_read_pkginfo( ifid(pp) )
        call fio_valid_pkginfo( ifid(pp) )
        call fio_copy_datainfo( ifid(pp), ifid(1) )
     else if ( .not. datainfo_nodep_pe  &
               .and. trim(large_memory_var(1)) /= '' ) then
        ! also read field data and store them in temporary buffer
        call fio_read_allinfo_tmpdata( ifid(pp) )
     else
        ! normal way to read pkginfo and datainfo
        call fio_read_allinfo( ifid(pp) )
     end if
     ! -->

     if ( pp == 1 ) then ! only once
        allocate( hinfo%rgnid(MNG_prc_rnum(pp)) )

        call fio_get_pkginfo(ifid(pp),hinfo)

        num_of_data = hinfo%num_of_data
        !write(*,*) '*** get variable informations'
        !write(*,*) 'num_of_data    : ', num_of_data

        nvar = 0
        do did = 0, num_of_data-1
           call fio_get_datainfo(ifid(pp),did,dinfo)

           if (allvar) then ! output all variables
              addvar = .true.
           else             ! select valiables to output
              addvar = .false.

              do v = 1, max_nvar
                 if ( trim(selectvar(v)) == trim(dinfo%varname) ) then
                    addvar = .true.
                    exit
                 elseif( trim(selectvar(v)) == '' ) then
                    exit
                 endif
              enddo
           endif

           do v = 1, nvar
              if ( trim(var_name(v)) == trim(dinfo%varname) ) then
                 var_nstep(v) = var_nstep(v) + 1

                 if( var_nstep(v) == 2 ) var_dt(v) = dinfo%time_start - var_time_str(v)

                 if( var_nstep(v) == step_str ) var_time_str(v) = dinfo%time_start ! [mod] H.Yashiro 20111003

                 addvar = .false.
                 exit
              endif
           enddo

           if (addvar) then
              nvar = nvar + 1
              var_nstep    (nvar) = 1
              var_name     (nvar) = dinfo%varname
              var_desc     (nvar) = dinfo%description
              var_unit     (nvar) = dinfo%unit
              var_layername(nvar) = dinfo%layername
              var_datatype (nvar) = dinfo%datatype
              var_nlayer   (nvar) = dinfo%num_of_layer
              var_time_str (nvar) = dinfo%time_start
              var_dt       (nvar) = dinfo%time_end - dinfo%time_start

              layerfile = trim(layerfile_dir)//'/'//trim(dinfo%layername)//'.txt'

              fid = MISC_get_available_fid()
              open(fid,file=trim(layerfile),form='formatted',status='old',iostat=ierr)
                 if ( ierr /= 0 ) then
                    write(ADM_LOG_FID,*) 'xxx layerfile doesnt exist!', trim(layerfile)
                    stop
                 endif

                 read(fid,*) kmax
                 do k = 1, kmax
                    read(fid,'(F16.4)') var_zgrid(k,nvar)
                 enddo
              close(fid)
           endif

        enddo !--- did LOOP
     endif !--- PE=000000

     deallocate( prc_tab_C )
     call fio_fclose(ifid(pp)) ! [add] 13-04-18
  end do

  if ( nvar == 0 ) then
     write(ADM_LOG_FID,*) 'No variables to convert. Finish.'
     stop
  endif

  write(ADM_LOG_FID,*) '########## Variable List ########## '
  write(ADM_LOG_FID,*) 'ID |NAME            |STEPS|Layername       |START FROM         |DT [sec]'
  do v = 1, nvar
     call calendar_ss2yh( date_str(:), real(var_time_str(v),kind=8) )
     write(tmpl,'(I4.4,"/",I2.2,"/",I2.2,1x,I2.2,":",I2.2,":",I2.2)') date_str(:)
     write(ADM_LOG_FID,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A19,A1,I8)') &
              v,'|',var_name(v),'|',var_nstep(v),'|',var_layername(v),'|', tmpl,'|', var_dt(v)
  enddo

  write(ADM_LOG_FID,*) '*** convert start : PaNDa format to lat-lon data'
  FLUSH(ADM_LOG_FID)
  !tmp_time(5)=MPI_WTIME()

  !#########################################################
  !--- start weighting summation

  variable_loop: do v = 1, nvar

     kmax = var_nlayer(v)
     write(ADM_LOG_FID,*) '       v=',v
     write(ADM_LOG_FID,*) '    kmax=',kmax
     write(ADM_LOG_FID,*) 'var_name=',trim(var_name(v))

     num_of_step = min(step_end,var_nstep(v)) - step_str + 1  ! [mov] 13-04-18

     step_loop: do t = 1, num_of_step

        nowsec = var_time_str(v) + (t-1)*var_dt(v)
        !print*,'nowsec: ',nowsec

        step = t-1 + step_str

         allocate( data4allrgn(ADM_gall*kmax*LALL_local) )
         allocate( data8allrgn(ADM_gall*kmax*LALL_local) )
         allocate( icodata4   (ADM_gall,kmax,LALL_local) )
         data4allrgn(:)  = CNST_UNDEF4
         data8allrgn(:)  = CNST_UNDEF
        !PE_loop: do p = 1, MNG_PALL
        PE_loop: do p = pstr, pend
           pp = p - pstr + 1

           call fio_fopen(ifid(pp),FIO_FREAD)  ! [add] 13-04-18
           if ( t==1 ) write(ADM_LOG_FID,'(A10)',advance='no') ' ->region:'

           !--- seek data ID and get information
           call fio_seek_datainfo(did,ifid(pp),var_name(v),step)
           call fio_get_datainfo(ifid(pp),did,dinfo)

           !--- verify
           if ( did == -1 ) then
              write(ADM_LOG_FID,*) 'xxx data not found! varname:',trim(var_name(v)),", step : ",step
              call ADM_proc_stop
              stop
           endif

           !--- read from pe000xx file

           if ( dinfo%datatype == FIO_REAL4 ) then
              if ( trim(large_memory_var(1)) /= '' ) then
                 call fio_read_data_tmpdata(ifid(pp),did,data4allrgn(:))
              else
                 call fio_read_data(ifid(pp),did,data4allrgn(:))
              end if

           elseif( dinfo%datatype == FIO_REAL8 ) then
              if ( trim(large_memory_var(1)) /= '' ) then
                 call fio_read_data_tmpdata(ifid(pp),did,data8allrgn(:))
              else
                 call fio_read_data(ifid(pp),did,data8allrgn(:))
              end if

              data4allrgn(:) = real(data8allrgn(:),kind=4)
              where( data8allrgn(:) == CNST_UNDEF )
                 data4allrgn(:) = CNST_UNDEF4
              endwhere

           endif
           icodata4(:,:,:) = reshape( data4allrgn(:), shape(icodata4) )
           if( trim(var_name(v))=='ms_pres' .or. trim(var_name(v))=='ss_ps') then
             icodata4(:,:,:)=log(icodata4(:,:,:)*0.01)
           end if
           write(ADM_LOG_FID,*) minval(icodata4(:,:,:)), maxval(icodata4(:,:,:))

           do l = 1, MNG_prc_rnum(p)
              do i = 1, nsite
                 if( inprc(i) ) then
                   do k = 1, kmax
                      fac1 = w1(i)
                      fac2 = w2(i)
                      fac3 = w3(i)
                      fac_sum = fac1 + fac2 + fac3
                      obsdata(k,i,v) = ( fac1 * icodata4(n1_index(i),k,l) &
                                       + fac2 * icodata4(n2_index(i),k,l) &
                                       + fac3 * icodata4(n3_index(i),k,l) &
                                       ) / fac_sum
                      !obsdata(k,i,v) = exp(obsdata(i))*kfact(i)
                      ! if(trim(var_name(v))=='ms_qv' .or. trim(var_name(v))=='ss_q2m' ) then
                      !   if(obsdata(k,i,v) .lt. 0.0d0 ) then
                      !     obsdata(k,i,v)=1.0d-5
                      !   end if
                      ! end if

                   end do
                   if(trim(var_name(v))=='ms_pres' .or. trim(var_name(v))=='ss_ps') then
                      write(ADM_LOG_FID,*) trim(var_name(v)), minval(obsdata(1:kmax,i,v)), maxval(obsdata(1:kmax,i,v))
                      flush(ADM_LOG_FID)
                      obsdata(1:kmax,i,v)=exp(obsdata(1:kmax,i,v))
                   else if(trim(var_name(v))=='ms_qv' &
                      .or. trim(var_name(v))=='ss_q2m' ) then
                      obsdata(1:kmax,i,v)=obsdata(1:kmax,i,v)*q2ppmv
                   end if
                 end if
              end do
           end do

           call fio_fclose(ifid(pp)) ! [add] 13-04-18

        enddo PE_loop

        deallocate( data4allrgn )
        deallocate( data8allrgn )
        deallocate( icodata4    )

     enddo step_loop ! step LOOP

  enddo variable_loop ! variable LOOP

  do v = 1, nvar
     write(ADM_LOG_FID,'(i6, 5f12.5)') v, (obsdata(1,k,v),k=1,5)
  end do

  oqc(:,:)=0

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(5)=MPI_WTIME()

  call MPI_ALLREDUCE(obsdata, obsdata_out, nlev*nsite*nvar, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)

  tmp_time(6)=MPI_WTIME()
 
  do v = 1, nvar
     write(ADM_LOG_FID,'(i6, 5f12.5)') v, (obsdata_out(1,k,v),k=1,5)
  end do
  FLUSH(ADM_LOG_FID)

  i=1
  do nn = 1, ninstrument
    do n = 1, ntvs(nn)
      obsdata_jprb(:,n,nn,:)=real(obsdata_out(:,i,:),kind=jprb)
      i=i+1
    end do
  end do

  !allocate( bt(15,nsite) )
  allocate( bt(15,maxtvsprof,ninstrument,nslots) )
  !allocate( bt_tmp(15,nsite) )
  allocate( bt_tmp(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tran(nlev,15,maxtvsprof,ninstrument,nslots) )

  lon(:)=lon(:)/CNST_PI*180.0d0
  lat(:)=lat(:)/CNST_PI*180.0d0

  do i = 1, nsite
    if(lon(i) < 0.0d0) lon(i)=lon(i)+360.d0
  end do

  if( ADM_prc_me == 1 ) then
    !if(nobs_noaa15 /= 0 ) then
    if(ntvs(1) /= 0 ) then
      WRITE(ADM_LOG_FID,*) 'NOAA-15', ntvs(1)
      FLUSH(ADM_LOG_FID)
      sobs=1
      eobs=ntvs(1)
      write(ADM_LOG_FID,*) nlev
      write(ADM_LOG_FID,*) ntvs(1)
      write(ADM_LOG_FID,*) trim(rttovcoef_fname(1))
      write(ADM_LOG_FID,*)maxval(obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(1),1,id_pres_nicam))
      write(ADM_LOG_FID,*)minval(obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(1),1,id_pres_nicam))
      write(ADM_LOG_FID,*)maxval(obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(1),1,id_temp_nicam))
      write(ADM_LOG_FID,*)minval(obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(1),1,id_temp_nicam))
      write(ADM_LOG_FID,*)maxval(obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(1),1,id_qvap_nicam))
      write(ADM_LOG_FID,*)minval(obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(1),1,id_qvap_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(1),1,id_tsfc_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(1),1,id_tsfc_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(1),1,id_qv2m_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(1),1,id_qv2m_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(1),1,id_surp_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(1),1,id_surp_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(1),1,id_u10m_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(1),1,id_u10m_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(1),1,id_v10m_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(1),1,id_v10m_nicam))
      write(ADM_LOG_FID,*) maxval( soza(sobs:eobs))
      write(ADM_LOG_FID,*) minval( soza(sobs:eobs))
      write(ADM_LOG_FID,*) maxval( soaz(sobs:eobs))
      write(ADM_LOG_FID,*) minval( soaz(sobs:eobs))
      write(ADM_LOG_FID,*) maxval( saza(sobs:eobs))
      write(ADM_LOG_FID,*) minval( saza(sobs:eobs))
      write(ADM_LOG_FID,*) maxval( saaz(sobs:eobs))
      write(ADM_LOG_FID,*) minval( saaz(sobs:eobs))
      flush(ADM_LOG_FID)

      call AMSUA_fwd(nlev, ntvs(1), rttovcoef_fname(1), &
                     obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(1),1,id_pres_nicam), &
                     obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(1),1,id_temp_nicam), &
                     obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(1),1,id_qvap_nicam), &
                     obsdata_jprb(1,1:ntvs(1),1,id_tsfc_nicam), &
                     obsdata_jprb(1,1:ntvs(1),1,id_qv2m_nicam), &
                     obsdata_jprb(1,1:ntvs(1),1,id_surp_nicam), &
                     obsdata_jprb(1,1:ntvs(1),1,id_u10m_nicam), &
                     obsdata_jprb(1,1:ntvs(1),1,id_v10m_nicam), &
                     soza(sobs:eobs), soaz(sobs:eobs), &
                     saza(sobs:eobs), saaz(sobs:eobs), &
                     elev(sobs:eobs)/1000.0d0, lon(sobs:eobs), &
                     lat(sobs:eobs), land(sobs:eobs), &
                     bt(:,1:ntvs(1),1,1), tran(:,:,1:ntvs(1),1,1) )
      nobs_tmp2=ntvs(1)
    end if

    if(ntvs(2) /= 0 ) then
      WRITE(ADM_LOG_FID,*) 'NOAA-16', ntvs(2)
      FLUSH(ADM_LOG_FID)
      sobs=nobs_tmp2+1
      eobs=nobs_tmp2+ntvs(2)
      write(ADM_LOG_FID,*) nlev
      write(ADM_LOG_FID,*) ntvs(2)
      write(ADM_LOG_FID,*) trim(rttovcoef_fname(2))
      write(ADM_LOG_FID,*) maxval(obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(2),2,id_pres_nicam))
      write(ADM_LOG_FID,*) minval(obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(2),2,id_pres_nicam))
      write(ADM_LOG_FID,*) maxval(obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(2),2,id_temp_nicam))
      write(ADM_LOG_FID,*) minval(obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(2),2,id_temp_nicam))
      write(ADM_LOG_FID,*) maxval(obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(2),2,id_qvap_nicam))
      write(ADM_LOG_FID,*) minval(obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(2),2,id_qvap_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(2),2,id_tsfc_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(2),2,id_tsfc_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(2),2,id_qv2m_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(2),2,id_qv2m_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(2),2,id_surp_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(2),2,id_surp_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(2),2,id_u10m_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(2),2,id_u10m_nicam))
      write(ADM_LOG_FID,*) maxval( obsdata_jprb(1,1:ntvs(2),2,id_v10m_nicam))
      write(ADM_LOG_FID,*) minval( obsdata_jprb(1,1:ntvs(2),2,id_v10m_nicam))
      write(ADM_LOG_FID,*) maxval( soza(sobs:eobs))
      write(ADM_LOG_FID,*) minval( soza(sobs:eobs))
      write(ADM_LOG_FID,*) maxval( soaz(sobs:eobs))
      write(ADM_LOG_FID,*) minval( soaz(sobs:eobs))
      write(ADM_LOG_FID,*) maxval( saza(sobs:eobs))
      write(ADM_LOG_FID,*) minval( saza(sobs:eobs))
      write(ADM_LOG_FID,*) maxval( saaz(sobs:eobs))
      write(ADM_LOG_FID,*) minval( saaz(sobs:eobs))
      flush(ADM_LOG_FID)

      call AMSUA_fwd(nlev, ntvs(2), rttovcoef_fname(2), &
                     obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(2),2,id_pres_nicam),&
                     obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(2),2,id_temp_nicam),&
                     obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(2),2,id_qvap_nicam),&
                     obsdata_jprb(1,1:ntvs(2),2,id_tsfc_nicam), &
                     obsdata_jprb(1,1:ntvs(2),2,id_qv2m_nicam), &
                     obsdata_jprb(1,1:ntvs(2),2,id_surp_nicam), &
                     obsdata_jprb(1,1:ntvs(2),2,id_u10m_nicam), &
                     obsdata_jprb(1,1:ntvs(2),2,id_v10m_nicam), &
                     soza(sobs:eobs), soaz(sobs:eobs), &
                     saza(sobs:eobs), saaz(sobs:eobs), &
                     elev(sobs:eobs)/1000.0d0, lon(sobs:eobs), &
                     lat(sobs:eobs), land(sobs:eobs), &
                     bt(:,1:ntvs(2),2,1), tran(:,:,1:ntvs(2),2,1) )
      nobs_tmp2=nobs_tmp2+ntvs(2)
    end if

    if(ntvs(3) /= 0 ) then
      WRITE(ADM_LOG_FID,*) 'NOAA-18', ntvs(3)
      FLUSH(ADM_LOG_FID)
      sobs=nobs_tmp2+1
      eobs=nobs_tmp2+ntvs(3)
      call AMSUA_fwd(nlev, ntvs(3), rttovcoef_fname(3), &
                     obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(3),3,id_pres_nicam),&
                     obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(3),3,id_temp_nicam),&
                     obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(3),3,id_qvap_nicam),&
                     obsdata_jprb(1,1:ntvs(3),3,id_tsfc_nicam), &
                     obsdata_jprb(1,1:ntvs(3),3,id_qv2m_nicam), &
                     obsdata_jprb(1,1:ntvs(3),3,id_surp_nicam), &
                     obsdata_jprb(1,1:ntvs(3),3,id_u10m_nicam), &
                     obsdata_jprb(1,1:ntvs(3),3,id_v10m_nicam), &
                     soza(sobs:eobs), soaz(sobs:eobs), &
                     saza(sobs:eobs), saaz(sobs:eobs), &
                     elev(sobs:eobs)/1000.0d0, lon(sobs:eobs), &
                     lat(sobs:eobs), land(sobs:eobs), &
                     bt(:,1:ntvs(3),3,1), tran(:,:,1:ntvs(3),3,1) )
      nobs_tmp2=nobs_tmp2+ntvs(3)
    end if

    if(ntvs(4) /= 0 ) then
      WRITE(ADM_LOG_FID,*) 'NOAA-19', ntvs(4)
      FLUSH(ADM_LOG_FID)
      sobs=nobs_tmp2+1
      eobs=nobs_tmp2+ntvs(4)
      call AMSUA_fwd(nlev, ntvs(4), rttovcoef_fname(4), &
                     obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(4),4,id_pres_nicam),&
                     obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(4),4,id_temp_nicam),&
                     obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(4),4,id_qvap_nicam),&
                     obsdata_jprb(1,1:ntvs(4),4,id_tsfc_nicam), &
                     obsdata_jprb(1,1:ntvs(4),4,id_qv2m_nicam), &
                     obsdata_jprb(1,1:ntvs(4),4,id_surp_nicam), &
                     obsdata_jprb(1,1:ntvs(4),4,id_u10m_nicam), &
                     obsdata_jprb(1,1:ntvs(4),4,id_v10m_nicam), &
                     soza(sobs:eobs), soaz(sobs:eobs), &
                     saza(sobs:eobs), saaz(sobs:eobs), &
                     elev(sobs:eobs)/1000.0d0, lon(sobs:eobs), &
                     lat(sobs:eobs), land(sobs:eobs), &
                     bt(:,1:ntvs(4),4,1), tran(:,:,1:ntvs(4),4,1) )
      nobs_tmp2=nobs_tmp2+ntvs(4)
    end if

    if(ntvs(5) /= 0 ) then
      WRITE(ADM_LOG_FID,*) 'METOP-2', ntvs(5)
      FLUSH(ADM_LOG_FID)
      sobs=nobs_tmp2+1
      eobs=nobs_tmp2+ntvs(5)
      call AMSUA_fwd(nlev, ntvs(5), rttovcoef_fname(5), &
                     obsdata_jprb(var_nlayer(id_pres_nicam):1:-1,1:ntvs(5),5,id_pres_nicam),&
                     obsdata_jprb(var_nlayer(id_temp_nicam):1:-1,1:ntvs(5),5,id_temp_nicam),&
                     obsdata_jprb(var_nlayer(id_qvap_nicam):1:-1,1:ntvs(5),5,id_qvap_nicam),&
                     obsdata_jprb(1,1:ntvs(5),5,id_tsfc_nicam), &
                     obsdata_jprb(1,1:ntvs(5),5,id_qv2m_nicam), &
                     obsdata_jprb(1,1:ntvs(5),5,id_surp_nicam), &
                     obsdata_jprb(1,1:ntvs(5),5,id_u10m_nicam), &
                     obsdata_jprb(1,1:ntvs(5),5,id_v10m_nicam), &
                     soza(sobs:eobs), soaz(sobs:eobs), &
                     saza(sobs:eobs), saaz(sobs:eobs), &
                     elev(sobs:eobs)/1000.0d0, lon(sobs:eobs), &
                     lat(sobs:eobs), land(sobs:eobs), &
                     bt(:,1:ntvs(5),5,1), tran(:,:,1:ntvs(5),5,1) )
      nobs_tmp2=nobs_tmp2+ntvs(5)
    end if

    tmp_time(7)=MPI_WTIME()
    write(ADM_LOG_FID,*) 'RTTOV END'
    FLUSH(ADM_LOG_FID)
    maxtvsprof=maxval(ntvs(:))
    allocate( weight(nlev,15,maxtvsprof,ninstrument) )
    allocate( weight_maxlev(15,maxtvsprof,ninstrument) )
  
    iobs=1
    do nn = 1, ninstrument
      do n =  1, ntvs(nn)
        do ichan = tvsch(1,nn), tvsch(ntvsch(nn),nn)
          weight(1,ichan,n,nn)= -(tran(2,ichan,n,nn,1)-tran(1,ichan,n,nn,1))/&
                         (log(obsdata_out(nlev-1,iobs,id_pres_nicam))&
                         -log(obsdata_out(nlev,iobs,id_pres_nicam)))
          do ilev = 2, nlev-1
             weight(ilev,ichan,n,nn)= -(tran(ilev+1,ichan,n,nn,1)-tran(ilev-1,ichan,n,nn,1))/&
                         (log(obsdata_out((nlev-ilev),iobs,id_pres_nicam))&
                         -log(obsdata_out((nlev-ilev+2),iobs,id_pres_nicam)))
          end do
          weight(nlev,ichan,n,nn)= -(tran(nlev,ichan,n,nn,1)-tran(nlev-1,ichan,n,nn,1))/&
                         (log(obsdata_out(1,iobs,id_pres_nicam))&
                         -log(obsdata_out(2,iobs,id_pres_nicam)))
          tmp_lev=-1.0
          do ilev = 1, nlev
             if( weight(ilev,ichan,n,nn) > tmp_lev ) then
                weight_maxlev(ichan,n,nn)=nlev-ilev+1
                tmp_lev=weight(ilev,ichan,n,nn)
             end if
          end do
        end do
        iobs=iobs+1
      end do
    end do
  
    write(*,*) 'test'
    nn=1
    do ichan = tvsch(1,nn), tvsch(ntvsch(nn),nn)
      write(*,*) 'channel', ichan
      do ilev = 1, nlev
        write(*,'(i5,3f10.3,i5)') ilev, tran(ilev,ichan,1,1,1), weight(ilev,ichan,1,1), &
              obsdata_out((nlev-ilev+1),1,id_pres_nicam), weight_maxlev(ichan,1,1)
      end do
    end do
  
    oqc(:,:)=1

  !
  ! BIAS CORRECTION
  !
    maxtvsprof=maxval(ntvs(:))
    ALLOCATE(vbc_pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots))
    vbc_pred(:,:,:,:,:)=0.0d0
  !
    iobs1=0
    do nn = 1, ninstrument
      do n =  1, ntvs(nn)
        iobs1=iobs1+1
        vbc_pred(1,:,n,nn,1)=undef ! IWLR is depending on ch (calc. on part2)
        vbc_pred(2,:,n,nn,1)=undef
        vbc_pred(3,:,n,nn,1)=0.d0
        vbc_pred(4,:,n,nn,1)=(obsdata_out(1,iobs1,id_tsfc_nicam)-273.15d0)/10.d0
        vbc_pred(5,:,n,nn,1)=0.d0
        vbc_pred(6,:,n,nn,1)=obsdata_out(1,iobs1,id_cldw_nicam)/30.0d0
        vbc_pred(7,:,n,nn,1)=1.d0/cos(saza(iobs1)*CNST_PI/180.d0)
        vbc_pred(8,:,n,nn,1)=0.d0
        do ic = 1, ntvsch(nn)
          iwlr=0.0d0
          do ilev = 1, nlev-1
            if(real(obsdata_jprb(ilev,n,nn,id_pres_nicam),kind=8)>200.0d0 .and. &
               real(obsdata_jprb(ilev,n,nn,id_pres_nicam),kind=8)<850.0d0) then
              iwlr = iwlr &
                   +(real(obsdata_jprb(ilev+1,n,nn,id_temp_nicam),kind=8)   &
                    -real(obsdata_jprb(ilev  ,n,nn,id_temp_nicam),kind=8))* &
                    (tran(nlev-ilev  , tvsch(ic,nn),n,nn,1)&
                    -tran(nlev-ilev+1, tvsch(ic,nn),n,nn,1))
            end if
          end do
          vbc_pred(1,ic,n,nn,1)=iwlr
        end do
        do ic = 1, ntvsch(nn)
          iwlr=0.0d0
          do ilev = 1, nlev-1
            if(real(obsdata_jprb(ilev,n,nn,id_pres_nicam),kind=8)>50.0d0   .and.&
               real(obsdata_jprb(ilev,n,nn,id_pres_nicam),kind=8)<200.0d0) then
              iwlr = iwlr &
                   +(real(obsdata_jprb(ilev+1,n,nn,id_temp_nicam),kind=8)   &
                    -real(obsdata_jprb(ilev  ,n,nn,id_temp_nicam),kind=8))* &
                    (tran(nlev-ilev  , tvsch(ic,nn),n,nn,1)&
                    -tran(nlev-ilev+1, tvsch(ic,nn),n,nn,1))
            end if
          end do
          vbc_pred(2,ic,n,nn,1)=iwlr
        end do

      end do
    end do

    CALL vbc_read('vbcf_coef.txt',0,vbc)
  
    CALL vbc_scan_read('vbcf_scanbias_coef.txt',0,vbcf_scan)

!    do nn = 1, ninstrument
!      write(ADM_LOG_FID,*) '### BEFORE normalize', sum(vbcf_scan(1:nfootp(nn),1,nn))/real(nfootp(nn)-6)
!      do ichan=1,ntvsch(nn)
!        vbcf_scan_tmp(1:nfootp(nn),ichan,nn) = vbcf_scan(1:nfootp(nn),ichan,nn) - &
!                                        sum(vbcf_scan(1:nfootp(nn),ichan,nn)) / real(nfootp(nn)-6)
!      end do
!      write(ADM_LOG_FID,*) '###  AFTER normalize',sum(vbcf_scan_tmp(1:nfootp(nn),1,nn))/real(nfootp(nn)-6)
!    end do

    do nn = 1, ninstrument
      do n =  1, ntvs(nn)
        ifoot=tvsfoot(n,nn,1)
        do ic=1,ntvsch(nn)
          ichan=tvsch(ic,nn)
          if(ichan<=9) then
            tvsdat(ichan,n,nn,1)=tvsdat(ichan,n,nn,1)-vbcf_scan(ifoot,ic,nn)
          end if
        end do
      end do
    end do

!    do nn = 1, ninstrument
!      write(ADM_LOG_FID,*) maxval(tvsdat(1:ntvsch(nn),1:ntvs(nn),nn,1))
!      write(ADM_LOG_FID,*) minval(tvsdat(1:ntvsch(nn),1:ntvs(nn),nn,1))
!    end do
    do nn = 1, ninstrument
      do n =  1, ntvs(nn)
        do ic=1,ntvsch(nn)
          ichan=tvsch(ic,nn)
          if(ichan<=9) then
            tvsdat(ichan,n,nn,1)=tvsdat(ichan,n,nn,1)+&
                   sum(vbc_pred(:,ic,n,nn,1)*vbc(:,ic,nn))
          end if
        end do
      end do
    end do

!    do nn = 1, ninstrument
!      write(ADM_LOG_FID,*) maxval(tvsdat(1:ntvsch(nn),1:ntvs(nn),nn,1))
!      write(ADM_LOG_FID,*) minval(tvsdat(1:ntvsch(nn),1:ntvs(nn),nn,1))
!    end do
    bt_tmp(:,:,:,:)=0.0d0
    do nn = 1, ninstrument
      do n =  1, ntvs(nn)
        do ic = 1, ntvsch(nn)
          bt_tmp(ic,n,nn,1)=bt(tvsch(ic,nn),n,nn,1)
        end do
      end do
    end do

!    do nn = 1, ninstrument
!      write(ADM_LOG_FID,*) maxval(tvsdat(1:ntvsch(nn),1:ntvs(nn),nn,1))
!      write(ADM_LOG_FID,*) minval(tvsdat(1:ntvsch(nn),1:ntvs(nn),nn,1))
!    end do
    do nn = 1, ninstrument
      do n =  1, ntvsprof(nn)
        !! [QUALITY CHECK FOR AMSU-A]
        !! Add 2016/02/02 Avoid cloud over land (refer to Bormann, 2010)
        !! Clear sky       : chs <=5 are not assimilated
        !! Clear sky       : ch6 is assimilated if z>1500m
        !! Clear sky       : ch7 is assimilated if z>2500m
        !! Clear sky       : chs >= 8 are assimilated
        !! Cloudy or rainy : no assimilation (chs >= 9 are assimilated in JMA)
        !if( lsql(n,nn,islot)==0 ) then
        !  if( abs( tvsdat(4,n,nn,islot)-bt(4,n,nn,islot) ) < 0.7 ) then
        !    do ic = 1, ntvsch(nn)
        !      ! Clear sky       : chs <=5 are not assimilated
        !      if(tvsch(ic,nn) <= 5) then
        !        tvsqc(ic,n,nn,islot)=0
        !      end if
        !      ! Clear sky       : ch6 is assimilated if z>1500m
        !      if(tvsch(ic,nn) == 6 .and. tvselev(n,nn,islot) > 1500) then
        !        tvsqc(ic,n,nn,islot)=0
        !      end if
        !      ! Clear sky       : ch7 is assimilated if z>2500m
        !      if(tvsch(ic,nn) == 7 .and. tvselev(n,nn,islot) > 2500) then
        !        tvsqc(ic,n,nn,islot)=0
        !      end if
        !    end do
        !  end if
        !  if( abs( tvsdat(4,n,nn,islot)-bt(4,n,nn,islot) ) > 0.7 ) then
        !    do ic = 1, ntvsch(nn)
        !      if(tvsch(ic,nn) <= 8) then
        !        tvsqc(ic,n,nn,islot)=0
        !      end if
        !    end do
        !  end if
        !end if
        !! Add 2016/02/02 Avoid cloud over ocean (refer to Bormann, 2010)
        !! Clear sky       : chs >= 4 are assimilated
        !! Cloudy sky      : chs >= 7 are assimilated
        !!!! rainy           : chs >= 9 are assimilated [not ready]
        !if( lsql(n,nn,islot)==1 ) then
        !  if( abs( tvsdat(3,n,nn,islot)-bt(3,n,nn,islot) ) > 3.0 ) then
        !    do ic = 1, ntvsch(nn)
        !      if( tvsch(ic,nn) < 7 ) then
        !        tvsqc(ic,n,nn,islot)=0
        !      end if
        !    end do
        !  end if
        !end if

        ! Do not assimilate over land
        if( lsql(n,nn,1)==0 ) then
          tvsqc(:,n,nn,1)=0
        end if
        do ic=1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          !if(tvslat(n,nn,1)>=-60.0 .and. lsql(n,nn,1)==0 ) then
          !  tvsqc(ichan,n,nn,1)=0
          !end if
          !if(tvsch(ichan,nn)<9) then
          if(tvslat(n,nn,1)<-60.0) tvsqc(ichan,n,nn,1)=0
          if(tvslat(n,nn,1)> 60.0) tvsqc(ichan,n,nn,1)=0
          if(tvsdat(ichan,n,nn,1)<100.0 .or. tvsdat(ichan,n,nn,1)>400.0) then
            tvsqc(ichan,n,nn,1)=0
          end if
        end do
      end do
    end do

    do nn = 1, ninstrument
      do n =  1, ntvs(nn)
        do ic = 1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          if(abs(tvsdat(ichan,n,nn,1)-bt(ichan,n,nn,1))>tvserr(ichan,n,nn,1)*5.0d0) then
            tvsqc(ichan,n,nn,1)=0
          end if
        end do
      end do
    end do

    call das_vbc(bt_tmp, vbc_pred, vbc, vbca)

    call vbc_write('vbca_coef.txt',0,vbca)
  
  end if
  tmp_time(8)=MPI_WTIME()

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if( ADM_prc_me==1 ) then
    call vbc_scan_read('vbcf_scanbias_coef.txt',0,vbcf_scan)

    vbca_scan(:,:,:)=0.0
    open(100,file='scanbias.txt')
    write(*,*) 'TEST SCANBIAS READING'
    do nn = 1, ninstrument
      do ifoot = 1, maxtvsfoot
        read(100,*) tvsname_scan, n, vbca_scan(ifoot,1:ntvsch(nn),nn)
        write(*,'(2i5,30f9.3)') nn, ifoot, vbca_scan(ifoot,1:ntvsch(nn),nn)
      end do
    end do   

    vbcf_scan(:,:,:)=0.97*vbcf_scan(:,:,:)+0.03*vbca_scan(:,:,:)

    call vbc_scan_write('vbca_scanbias_coef.txt',0,vbcf_scan)

  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize:       ', tmp_time(2)-tmp_time(1)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Read obs:         ', tmp_time(3)-tmp_time(2)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Calc coef:        ', tmp_time(4)-tmp_time(3)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Prepare read gues:', tmp_time(5)-tmp_time(4)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Communication:    ', tmp_time(6)-tmp_time(5)
  write(ADM_LOG_FID,'(A,F15.5)') '##### rttov computation:', tmp_time(7)-tmp_time(6)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Output:           ', tmp_time(8)-tmp_time(7)

  call ADM_proc_finish

contains
  !-----------------------------------------------------------------------------
  !> read option
  !-----------------------------------------------------------------------------
  subroutine readoption
    use mod_misc, only : &
      MISC_get_available_fid
    use mod_tool_option, only: &
      OPT_convert, &
      OPT_fid
    implicit none

    integer :: io
    !---------------------------------------------------------------------------

    ! --- Set option
    OPT_fid = MISC_get_available_fid()
    open(OPT_fid,status='SCRATCH')

      call OPT_convert( fmax )

      read(OPT_fid,nml=OPTION,iostat=io)

    close(OPT_fid)

    if (      io /= 0     &
         .OR. fmax == 0   &
         .OR. fmax > flim &
         .OR. help        ) call helpoption

  end subroutine readoption

  !-----------------------------------------------------------------------------
  !> display help for option and abort
  !-----------------------------------------------------------------------------
  subroutine helpoption
    implicit none
    !---------------------------------------------------------------------------

    write(*,OPTION)

    stop
  end subroutine helpoption


  !S.Iga051226 =>
  !-----------------------------------------------------------------------------
  function sec2initplate(datesec) result(template)
    !-- output grads-like template part  like 01JAN0000
    implicit none

    integer(8)        :: datesec
    ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: template
    character(LEN=20) :: template

    integer :: d(6)

    character(LEN=3) :: nmonth(12)
    data nmonth / 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC' /
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call calendar_ss2yh( d(:), real(datesec,kind=8) )

    write(template,'(I2.2,A1,I2.2,A1,I2.2,A3,I4.4)') &
                              d(4), ':', d(5), 'Z', d(3), nmonth(d(2)), d(1)

  end function sec2initplate

  !-----------------------------------------------------------------------------
  function sec2template(datesec) result(template)
    !-- output grads-like template part  like 2005-12-01-23h50m
    implicit none

    integer(8)        :: datesec
    ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: template
    character(LEN=20) :: template

    integer :: d(6)
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call calendar_ss2yh( d(:), real(datesec,kind=8) )

    write(template,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1)') &
                          d(1), '-', d(2), '-', d(3), '-', d(4), 'h', d(5), 'm'

  end function sec2template

  !-----------------------------------------------------------------------------
  function timeincrement(isec) result(template)
    implicit none

    integer       :: isec
    character(20) :: template

    character(18):: tmp
    !---------------------------------------------------------------------------

    write(tmp,*) max(isec/60, 1)

    template = trim(tmp)//'mn'

  end function timeincrement
  !S.Iga051226 <=

  !-----------------------------------------------------------------------------
  function calendar_ss2cc_gtool(datesec) result(template)
    !--- calendar, sec. -> character (YYYYMMDD HHMMSS)
    implicit none

    integer(8)        :: datesec
    character(LEN=16) :: template

    integer :: d(6), i
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call calendar_ss2yh( d(:), real(datesec,kind=8) )

    write (template,'(i4.4,i2.2,i2.2,1x,i2.2,i2.2,i2.2,1x)') (d(i),i=1,6)

  end function calendar_ss2cc_gtool

  subroutine geticogrid
    !
    use mod_adm, only : &
         ADM_prc_run_master, &
         ADM_IopJop_nmax,   &
         ADM_IopJop,         &
         ADM_GIoJo
    !
    implicit none
    !
    logical :: flag
    real(8) :: dif
    real(8) :: dist, mindist
    real(8) :: r
    logical :: first_hit
    integer :: nsum
    integer :: s, m, l, n, nn
    integer :: cnt
    integer :: sumcnt
    integer :: prc
    integer :: minprc
    logical :: failgetgrd = .false.
    character(ADM_NSYS) :: sfc_type0 ! Y.Niwa 070907
    !
    do s=1, nsite
       !
       sfc_type0 = sfc_type   ! Y.Niwa 070907
       !
       m = 0
       do
          m = m + 1
          r = alpha + dalpha*dble(m-1)
          dist = CNST_PI
          !
          do l=1, ADM_lall
             !
             do n=1, ADM_IopJop_nmax
                !
                nn = ADM_IopJop(n,ADM_GIoJo)
                first_hit = .false.
                !
                call in_region(lat(s), lon(s), r,     &
                     icolat(nn,l), icolon(nn,l),        &
                     flag, dif, idirec(1,s))
                if(flag) then
                   if(dif < dist) then
                      if( (trim(sfc_type0) /= 'FREE')  & !Y.Niwa 070907 sfc_type => sfc_type0
                           .and. (land(s) == 0 .or. land(s) == 1) ) then
                         if( land(s) == icoland(nn,l) ) then
                            n1_index(s) = nn
                            l_index(s)  = l
                            land_model(s) = icoland(nn,l)
                            idirec_model(s) = idirec(1,s)
                            first_hit = .true.
                            dist = dif
                         end if
                      else
                         n1_index(s) = nn
                         l_index(s)  = l
                         land_model(s) = icoland(nn,l)
                         idirec_model(s) = idirec(1,s)
                         first_hit = .true.
                         dist = dif
                      end if
                   end if
                end if
                !
                if(.not. first_hit) then
                   if(idirec(2,s) /= 0) then
                      call in_region(lat(s), lon(s), r,     &
                           icolat(nn,l), icolon(nn,l),          &
                           flag, dif, idirec(2,s))
                      if(flag) then
                         if(dif < dist) then
                            if( (trim(sfc_type0) /= 'FREE') & !Y.Niwa 070907 sfc_type => sfc_type0
                                 .and. (land(s) == 0 .or. land(s) == 1) ) then
                               if( land(s) == icoland(nn,l) ) then
                                  n1_index(s) = nn
                                  l_index(s)  = l
                                  land_model(s) = icoland(nn,l)
                                  idirec_model(s) = idirec(2,s)
                                  inprc(s) = .true.
                                  dist = dif
                               end if
                            end if
                         end if
                      end if
                   end if
                end if
             end do
          end do
          !
          call MPI_ALLREDUCE(l_index(s), nsum, 1, MPI_INTEGER, MPI_SUM, &
               MPI_COMM_WORLD, ierr)
          !
          if(nsum /= 0) exit
          !
          if( m == search_itrmax ) then !Y.Niwa 070907 sfc_type => sfc_type0
             if(trim(sfc_type0) == 'ALMOST') then
                sfc_type0 = 'FREE'
             end if
          else if( m == limit_itrmax ) then
             exit
          end if
          !
       end do
       !
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       !
       call MPI_ALLREDUCE(dist, mindist, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
            MPI_COMM_WORLD, ierr)
       !
       cnt = 0
       if(dist /= mindist) then
          inprc(s) = .false.
       else
          inprc(s) = .true.
          cnt = 1
       end if
       !
       call MPI_ALLREDUCE(cnt, sumcnt, 1, MPI_INTEGER, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
       !
       if(sumcnt /= 1) then
          cnt = 0
          prc = 9999
          if(inprc(s)) prc = ADM_prc_me
          call MPI_ALLREDUCE(prc, minprc, 1, MPI_INTEGER, MPI_MIN, &
               MPI_COMM_WORLD, ierr)
          if(ADM_prc_me /= minprc) then
             inprc(s) = .false.
          else
             inprc(s) = .true.
             cnt = 1
          end if
          !
          call MPI_ALLREDUCE(cnt, sumcnt, 1, MPI_INTEGER, MPI_SUM, &
               MPI_COMM_WORLD, ierr)
          !
          if(sumcnt /= 1) then
             if(ADM_prc_me == ADM_prc_run_master) then
                write(ADM_LOG_FID,*) 'Failure to get grid!  obs. #:', obsnum(s), &
                     ' sumcnt = ', sumcnt
             end if
             failgetgrd = .true.
             !
          end if
       end if
       !
    end do
    !
    if(failgetgrd) call ADM_proc_stop
    !
    return
    !
  end subroutine geticogrid
  !=================================================================================
  subroutine in_region(alat, alon, r0, blat, blon, f, diff, dir)
    !
    implicit none
    !
    real(8), intent(in)  :: alat
    real(8), intent(in)  :: alon
    real(8), intent(in)  :: r0
    real(8), intent(in)  :: blat
    real(8), intent(in)  :: blon
    logical, intent(out) :: f
    real(8), intent(out) :: diff
    integer, intent(in) :: dir
    !
    real(8) :: ap(3)
    real(8) :: bp(3)
    real(8) :: cdott
    logical :: indir
    !
    diff = 999.99d9
    !
    ap(1) = cos(alat)*cos(alon)
    ap(2) = cos(alat)*sin(alon)
    ap(3) = sin(alat)
    !
    bp(1) = cos(blat)*cos(blon)
    bp(2) = cos(blat)*sin(blon)
    bp(3) = sin(blat)
    !
    cdott = ap(1)*bp(1) + ap(2)*bp(2) + ap(3)*bp(3)
    !
    if(cdott >  1.0d0) cdott =  1.0d0
    if(cdott < -1.0d0) cdott = -1.0d0
    !
    if(acos(cdott) <= r0) then
       if(dir /= 0) then
          indir = .false.
          call in_direction(alat, alon, blat, blon, dir, indir)
          if(indir) then
             f = .true.
             diff = abs(acos(cdott))
          else
             f = .false.
          end if
       else
          f = .true.
          diff = abs(acos(cdott))
       end if
    else
       f = .false.
    end if
    !
    return
    !
  end subroutine in_region
  !=================================================================================
  subroutine getklev
    !
    use mod_adm, only : &
         ADM_kmin,      &
         ADM_kmax,      &
         ADM_vlayer
    use mod_grd, only : &
         GRD_vz,   &
         GRD_Z,    &
         GRD_zs,   &
         GRD_ZSFC
    !
    implicit none
    !
    real(8) :: z(ADM_gall,ADM_vlayer,ADM_lall)
    integer :: s, k
    real(8) :: alt
    real(8) :: zs(ADM_gall,ADM_lall)
    character(12) :: name
    real(8) :: zz, zz0
    real(8) :: pp, pp0
    real(8) :: tt, tt0, tv
    real(8) :: qq, qq0
    real(8) :: t1, t2
    real(8) :: tmp_elev
    REAL(8),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]

    !
    z(:,1:ADM_vlayer,:) = GRD_vz(:,ADM_kmin:ADM_kmax,:,GRD_Z)
    zs(:,:) = GRD_zs(:,ADM_KNONE,:,GRD_ZSFC)

    do s=1, nsite
       if(inprc(s)) then
          !
          if( elem(s) == id_ps_obs ) then
             if(.not. h_intrpl) then
                alt = z(n1_index(s),1,l_index(s))
             else
                alt = w1(s)*zs(n1_index(s),l_index(s)) &
                    + w2(s)*zs(n2_index(s),l_index(s)) &
                    + w3(s)*zs(n3_index(s),l_Index(s))
                pp  = w1(s)*lnps(n1_index(s),l_index(s)) &
                    + w2(s)*lnps(n2_index(s),l_index(s)) &
                    + w3(s)*lnps(n3_index(s),l_Index(s))
                tt  = w1(s)*ts(n1_index(s),l_index(s)) &
                    + w2(s)*ts(n2_index(s),l_index(s)) &
                    + w3(s)*ts(n3_index(s),l_Index(s))
                qq  = w1(s)*qs(n1_index(s),l_index(s)) &
                    + w2(s)*qs(n2_index(s),l_index(s)) &
                    + w3(s)*qs(n3_index(s),l_Index(s))
             end if

             if( alt-elev(s) /= 0 .and. abs(alt-elev(s)) < 1000.0d0 ) then
               tv=tt*(1.0d0+0.608d0*qq)
               kfact(s) = (tv/(tv+gamma*(elev(s)-alt)))**(9.81d0/(gamma*287.d0))
             else
               kfact(s) = -999.0d0
             end if
             klev(:,s)=1
          else
             tmp_elev=log(elev(s))
             pp0=undef
             do k = 1, ADM_vlayer
                if(.not. h_intrpl) then
                   pp = lnpres(n1_index(s),k,l_index(s))
                else
                   pp = w1(s)*lnpres(n1_index(s),k,l_index(s)) &
                      + w2(s)*lnpres(n2_index(s),k,l_index(s)) &
                      + w3(s)*lnpres(n3_index(s),k,l_Index(s))
                end if
                if ( pp < tmp_elev ) exit
                pp0=pp
             end do
             write(ADM_LOG_FID,*) k, tmp_elev, pp0
             if( k == 1 ) then
                klev(1,s)=1
                klev(2,s)=1
                kfact(s)=1.0d0
                if( elem(s) == id_t_obs ) then
                   kfact(s) = -999.0d0
                end if
             else
                klev(1,s)=k-1
                klev(2,s)=k
                kfact(s) = ( pp - tmp_elev ) / ( pp - pp0 )
             end if ! k
          end if ! elem
       end if ! inprc
    end do
    !
    return
    !
  end subroutine getklev
  !=================================================================================
  subroutine in_direction(y0, x0, y1, x1, dr, flag)
    !
    implicit none
    !
    real(8), intent(in) :: x0, y0, x1, y1
    integer, intent(in) :: dr
    logical, intent(out) :: flag
    real(8) :: x, y
    real(8) :: r
    real(8) :: c, theta
    real(8) :: pi
    real(8) :: pi1, pi2, pi3, pi4
    !
    flag = .false.
    !
    pi = 2.0d0*asin(1.0d0)
    pi1 = pi/8.0d0
    pi2 = pi*3.0d0/8.0d0
    pi3 = pi*5.0d0/8.0d0
    pi4 = pi*7.0d0/8.0d0
    !
    x = x1 - x0
    if(x > 180.0d0) then
       if(x0 > 0.0d0) then
          x = 360.0d0 - x
       else
          x = x - 360.0d0
       end if
    end if
    !
    y = y1 - y0
    r = dsqrt(x*x + y*y)
    c = x/r
    if(c >  1.0d0) c =  1.0d0
    if(c < -1.0d0) c = -1.0d0
    !
    theta = acos(c)
    !
    if(y >= 0.0d0) then
       select case(dr)
       case(3)
          if(theta >= 0.0d0 .and. theta < pi1) flag = .true.
       case(2)
          if(theta >= pi1 .and. theta < pi2) flag = .true.
       case(1)
          if(theta >= pi2 .and. theta < pi3) flag = .true.
       case(8)
          if(theta >= pi3 .and. theta < pi4) flag = .true.
       case(7)
          if(theta >= pi4 .and. theta <= pi) flag = .true.
       end select
    else
       select case(dr)
       case(3)
          if(theta > 0.0d0 .and. theta < pi1) flag = .true.
       case(4)
          if(theta >= pi1 .and. theta < pi2) flag = .true.
       case(5)
          if(theta >= pi2 .and. theta < pi3) flag = .true.
       case(6)
          if(theta >= pi3 .and. theta < pi4) flag = .true.
       case(7)
          if(theta >= pi4 .and. theta < pi) flag = .true.
       end select
    end if
    !
    return
    !
  end subroutine in_direction

  !--------------------------------------------------------------------------------
  subroutine setup_ico2latlon_mapping( what_is_done )
    !
    !**************************************************
    ! + Imported from mod_latlon.f90
    !***************************************************
    !
    use mod_misc, only :               &
         MISC_triangle_area
    use mod_adm, only :              &
         ADM_LOG_FID,              &
         ADM_NSYS,                 &
         ADM_MAXFNAME,             &
         ADM_IooJoo_nmax,          &
         ADM_IooJoo,               &
         ADM_GIoJo,                &
         ADM_GIpJo,                &
         ADM_GIpJp,                &
         ADM_GIoJp,                &
         ADM_GImJo,                &
         ADM_GIoJm,                &
         ADM_GImJm,                &
         ADM_KNONE,                &
         ADM_TI,ADM_TJ,            &
         ADM_lall
    use mod_grd, only :               &
         GRD_x,GRD_x_pl,           &
         GRD_XDIR,                 &
         GRD_YDIR,                 &
         GRD_ZDIR,                 &
         GRD_rscale
    use mod_gmtr, only :            &
         GMTR_P_var, GMTR_P_var_pl, &
         GMTR_P_LAT,                &
         GMTR_P_LON,                &
         GMTR_polygon_type
    use mod_cnst, only :               &
         CNST_PI
    !
    implicit none
    !
    character(LEN=*), intent(in) :: what_is_done
    !
    real(8) :: r1(3), r2(3), r3(3), r0(3)
    real(8) :: v12(3), v23(3), v31(3), v10(3), v20(3), v30(3)
    real(8) :: nvec(3), v12xv10(3),v23xv20(3),v31xv30(3)
    real(8):: judge12, judge23, judge31, rf, rn
    integer :: n, l, i, t
    !
    real(8) :: latmin_l, latmax_l
    real(8) :: lat1
    real(8) :: lat2
    real(8) :: lat3
    real(8) :: lat4
    real(8) :: lonmin_l, lonmax_l
    real(8) :: lon1, lon2, lon3, lon4
    real(8) :: len, area_total
    !
    real(8) :: coslat(nsite), sinlat(nsite), coslon(nsite), sinlon(nsite)
    logical :: log_j(nsite)
    integer :: max_num_lon
    integer :: max_num_lat_lon2(ADM_lall)
    !
    !    ! S.Iga060210 =>
    !    real(8) :: epsi=1d-8
    !    ! S.Iga060210 >=
    !
    ! S.Iga060212 =>
    !-- infinite simal values to prevent 'missing grid' and 'NaN weight'
    !---(1) if you want to prevent the missing grids on the edge of triangles,
    !       give a value to epsi_inpro. However, you should take a risk to
    !       count-in grids even out of the triangle.
    !       In my experience, zero is OK at least for GLEVEL =< 11
    !---(2) if you want to prevent missing grids on the apex of triangles,
    !       give a value to epsi_grid.
    !       In my experience, this needs a value when epsi_inpro=0.
    !---(3) if you want to prevent 'NaN weight', give a value to epsi_area.
    !       In my experience, zero is OK at least for GLEVEL =< 11
    !
    real(8) :: epsi_inpro =0     ! marginal value for inner products used to
                                  ! calcuate judge
    real(8) :: epsi_grid  =0 ! marginal square near grid points (in radian)
    !real(8) :: epsi_grid  =1d-8 ! marginal square near grid points (in radian)
    real(8) :: epsi_area  =0     ! marginal value for inner products used in
                                 ! misc_triangle_area
    ! S.Iga060212 >=
    !
    do i=1,nsite
       coslon(i)=cos(lon(i))
       sinlon(i)=sin(lon(i))
       coslat(i)=cos(lat(i))
       sinlat(i)=sin(lat(i))
    enddo
    !
    !
    max_num_latlon = 0
    !  max_num_latlon2(:) = 0
    !
    do l=1,ADM_lall
       do n=1, ADM_IooJoo_nmax
          lat1=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GMTR_P_LAT)
          lat2=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GMTR_P_LAT)
          lat3=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GMTR_P_LAT)
          lat4=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GMTR_P_LAT)
          lon1=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GMTR_P_LON)
          lon2=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GMTR_P_LON)
          lon3=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GMTR_P_LON)
          lon4=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GMTR_P_LON)
          !
          do t=ADM_TI,ADM_TJ
             !
             !--- construct triagular vertice by clockwise way from the origin
             if(t==ADM_TI) then
                r1(1) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
                r1(2) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
                r1(3) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
                !
                r2(1)= GRD_x(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
                r2(2)= GRD_x(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
                r2(3)= GRD_x(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
                !
                r3(1) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
                r3(2) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
                r3(3) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
                !
             else !--- ADM_TJ 
                !
                r1(1) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
                r1(2) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
                r1(3) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
                !
                r2(1) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
                r2(2) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
                r2(3) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
                !
                r3(1)= GRD_x(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
                r3(2)= GRD_x(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
                r3(3)= GRD_x(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
                !
             end if
             !
             !--- 
             latmin_l=min(lat1,lat2,lat3,lat4)
             latmax_l=max(lat1,lat2,lat3,lat4)
             !
             lonmin_l=min(lon1,lon2,lon3,lon4)
             lonmax_l=max(lon1,lon2,lon3,lon4)
             if (lonmax_l-lonmin_l > CNST_PI) then
                if (lon1 < 0 ) lon1=lon1+CNST_PI+CNST_PI
                if (lon2 < 0 ) lon2=lon2+CNST_PI+CNST_PI
                if (lon3 < 0 ) lon3=lon3+CNST_PI+CNST_PI
                if (lon4 < 0 ) lon4=lon4+CNST_PI+CNST_PI
                lonmin_l=min(lon1,lon2,lon3,lon4)
                lonmax_l=max(lon1,lon2,lon3,lon4)
             endif
             ! S.Iga060210 =>
             lonmin_l = lonmin_l - epsi_grid
             lonmax_l = lonmax_l + epsi_grid
             latmin_l = latmin_l - epsi_grid
             latmax_l = latmax_l + epsi_grid
             ! S.Iga060210 >=
             !
             do i=1, nsite
                log_j(i)=.not.( (latmin_l<=lat(i)).and.(lat(i)<=latmax_l) )
             end do
             !
             do i=1, nsite
                if (.not.( ((lonmin_l<=lon(i)).and.(lon(i)<=lonmax_l))&
                     .or.((lonmin_l<=lon(i)+CNST_PI+CNST_PI).and.(lon(i)+CNST_PI+CNST_PI<=lonmax_l)) )   ) then
                   cycle
                end if
                !
                ! Somehow when this log_j(j) is used instead,  &
                ! results becomes slightly different in order  &
                ! of 1e-13 at least in Linux. (S.Iga)
                !                   if (log_j(j)) then
                if (.not.( (latmin_l<=lat(i)).and.(lat(i)<=latmax_l) )) then
                   cycle
                end if
                !
                !--- set side-vectors by clockwise way from the origin
                v12(1)=r2(1)-r1(1)
                v12(2)=r2(2)-r1(2)
                v12(3)=r2(3)-r1(3)
                !
                v23(1)=r3(1)-r2(1)
                v23(2)=r3(2)-r2(2)
                v23(3)=r3(3)-r2(3)
                !
                v31(1)=r1(1)-r3(1)
                v31(2)=r1(2)-r3(2)
                v31(3)=r1(3)-r3(3)
                !
                !
                !--- calculate normal unit vector to the plane triagle
                nvec(1)=v12(2)*v23(3)-v12(3)*v23(2)
                nvec(2)=v12(3)*v23(1)-v12(1)*v23(3)
                nvec(3)=v12(1)*v23(2)-v12(2)*v23(1)
                len = sqrt(&
                     nvec(1)*nvec(1)&
                     +nvec(2)*nvec(2)&
                     +nvec(3)*nvec(3))
                nvec(1:3)=nvec(1:3)/len
                !
                !--- target latlon point on the sphere
                r0(1)=coslat(i)*coslon(i)
                r0(2)=coslat(i)*sinlon(i)
                r0(3)=sinlat(i)
                !
                !--- remove the case inner product is negative
                if((r0(1)*r1(1)+r0(2)*r1(2)+r0(3)*r1(3))<0.0D0) then
                   cycle
                end if
                !
                !--- mapping r0 to the plane(r1,r2,r3)
                !
                !------ distance from origin to a plane with r0.
                rf=(r0(1)*nvec(1)+r0(2)*nvec(2)+r0(3)*nvec(3))
                !
                !------ distance from origin to a plane with r1 or (r2,r3).
                rn=(r1(1)*nvec(1)+r1(2)*nvec(2)+r1(3)*nvec(3))
                !
                !------ mapping r0 
                r0(1)=r0(1)*(rn/rf)
                r0(2)=r0(2)*(rn/rf)
                r0(3)=r0(3)*(rn/rf)
                !
                !--- calculate vectors from triangler points 
                !--- to the target point
                v10(1)=r1(1)-r0(1)
                v10(2)=r1(2)-r0(2)
                v10(3)=r1(3)-r0(3)
                !
                v20(1)=r2(1)-r0(1)
                v20(2)=r2(2)-r0(2)
                v20(3)=r2(3)-r0(3)
                v30(1)=r3(1)-r0(1)
                v30(2)=r3(2)-r0(2)
                v30(3)=r3(3)-r0(3)
                !
                !
                v12xv10(1)=v12(2)*v10(3)-v12(3)*v10(2)
                v12xv10(2)=v12(3)*v10(1)-v12(1)*v10(3)
                v12xv10(3)=v12(1)*v10(2)-v12(2)*v10(1)
                !
                v23xv20(1)=v23(2)*v20(3)-v23(3)*v20(2)
                v23xv20(2)=v23(3)*v20(1)-v23(1)*v20(3)
                v23xv20(3)=v23(1)*v20(2)-v23(2)*v20(1)
                !
                v31xv30(1)=v31(2)*v30(3)-v31(3)*v30(2)
                v31xv30(2)=v31(3)*v30(1)-v31(1)*v30(3)
                v31xv30(3)=v31(1)*v30(2)-v31(2)*v30(1)
                !
                judge12 = &
                     nvec(1)*v12xv10(1) &
                     +nvec(2)*v12xv10(2) &
                     +nvec(3)*v12xv10(3)
                judge23 = &
                     nvec(1)*v23xv20(1) &
                     +nvec(2)*v23xv20(2) &
                     +nvec(3)*v23xv20(3)
                judge31 = &
                     nvec(1)*v31xv30(1) &
                     +nvec(2)*v31xv30(2) &
                     +nvec(3)*v31xv30(3)
                !
                if(t==ADM_TI) then
                   if (      (judge12<=0.0D0 + epsi_inpro)& ! S.Iga add epsi_inpro  (060212)
                        .and.(judge23<= 0.0D0 + epsi_inpro)&
                        .and.(judge31<=0.0D0 + epsi_inpro) ) then
                      select case(trim(what_is_done))
                      case('GET_NUM')
                         max_num_latlon = max_num_latlon + 1
                         !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                      case('SET_INDEX')
                         !   max_num_latlon = max_num_latlon + 1
                         !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                         !   lon_index(max_num_latlon) = i
                         !   lat_index(max_num_latlon) = j
                         l_index(i) = l
                         !   t_index(max_num_latlon) = t
                         n1_index(i) = ADM_IooJoo(n,ADM_GIoJo)
                         n2_index(i) = ADM_IooJoo(n,ADM_GIpJo)
                         n3_index(i) = ADM_IooJoo(n,ADM_GIpJp)
                         w1(i) = MISC_triangle_area( &
                              r0(1:3),r2(1:3),r3(1:3),& ! S.Iga add epsi_area  (060212)
                              GMTR_polygon_type, 1.0D0,critical=epsi_area)
                         w2(i) = MISC_triangle_area( &
                              r0(1:3),r3(1:3),r1(1:3),& ! S.Iga add epsi_area  (060212)
                              GMTR_polygon_type, 1.0D0,critical=epsi_area)
                         w3(i) = MISC_triangle_area( &
                              r0(1:3),r1(1:3),r2(1:3),& ! S.Iga add epsi_area  (060212)
                              GMTR_polygon_type, 1.0D0,critical=epsi_area)
                         area_total &
                              = w1(i)&
                              + w2(i)&
                              + w3(i)
                         w1(i)=w1(i)/area_total
                         w2(i)=w2(i)/area_total
                         w3(i)=w3(i)/area_total
                         !
                      case default
                      end select

                      cycle
                      ! S.Iga060210=>
                      ! Deal with the points near 'r1'
                   elseif (abs(v10(1)) <= epsi_grid ) then
                      if (abs(v10(2)) <= epsi_grid ) then
                         if (abs(v10(3)) <= epsi_grid ) then
                            select case(trim(what_is_done))
                            case('GET_NUM')
                               max_num_latlon = max_num_latlon + 1
                               !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                            case('SET_INDEX')
                               !   max_num_latlon = max_num_latlon + 1
                               !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                               !   lon_index(max_num_latlon) = i
                               !   lat_index(max_num_latlon) = j
                               l_index(i) = l
                               !   t_index(max_num_latlon) = t
                               n1_index(i) = ADM_IooJoo(n,ADM_GIoJo)
                               n2_index(i) = ADM_IooJoo(n,ADM_GIpJo)
                               n3_index(i) = ADM_IooJoo(n,ADM_GIpJp)
                               w1(i) = 1d0
                               w2(i) = 0
                               w3(i) = 0
                            end select
                         endif
                      endif
                      ! S.Iga060210>=
                   end if
                else !--- ADM_TJ
                   if (      (judge12< 0.0D0 + epsi_inpro)& ! S.Iga add epsi_inpro  (060212)
                        .and.(judge23< 0.0D0 + epsi_inpro)&
                        .and.(judge31<=0.0D0 + epsi_inpro) ) then
                      select case(trim(what_is_done))
                      case('GET_NUM')
                         max_num_latlon = max_num_latlon + 1
                         !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                      case('SET_INDEX')
                         !   max_num_latlon = max_num_latlon + 1
                         !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                         !   lon_index(max_num_latlon) = i
                         !   lat_index(max_num_latlon) = j
                         l_index(i) = l
                         !   t_index(max_num_latlon) = t
                         n1_index(i) = ADM_IooJoo(n,ADM_GIoJo)
                         n2_index(i) = ADM_IooJoo(n,ADM_GIpJp)
                         n3_index(i) = ADM_IooJoo(n,ADM_GIoJp)
                         w1(i) = MISC_triangle_area( &
                              r0(1:3),r2(1:3),r3(1:3),& ! S.Iga add epsi_area  (060212)
                              GMTR_polygon_type, 1.0D0,critical=epsi_area)
                         w2(i) = MISC_triangle_area( &
                              r0(1:3),r3(1:3),r1(1:3),& ! S.Iga add epsi_area  (060212)
                              GMTR_polygon_type, 1.0D0,critical=epsi_area)
                         w3(i) = MISC_triangle_area( &
                              r0(1:3),r1(1:3),r2(1:3),& ! S.Iga add epsi_area  (060212)
                              GMTR_polygon_type, 1.0D0,critical=epsi_area)
                         area_total &
                              = w1(i)&
                              + w2(i)&
                              + w3(i)
                         w1(i)=w1(i)/area_total
                         w2(i)=w2(i)/area_total
                         w3(i)=w3(i)/area_total
                         !
                      case default
                      end select
                      cycle
                   end if
                end if
             end do
          end do
       end do
    end do
    !
  end subroutine setup_ico2latlon_mapping
  !-----------------------------------------------------------------------------
  subroutine get_nsite(fid,nsite)
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: nsite
    !
    real(4) :: wk(6)
    integer :: istat
    integer :: n

    nloop: do n = 1,1000000
      read(fid,iostat=istat) wk
      if( n == 1 )then
         write(ADM_LOG_FID,*) wk(1:5)
      endif
      if( istat /= 0 )then
        nsite = n - 1
        write(ADM_LOG_FID,*) 'nsite: ', nsite
        return
      endif
    enddo nloop

    return
  end subroutine get_nsite

  !-----------------------------------------------------------------------------
  subroutine get_variable(rid,elev,var)
    implicit none
    real(4),intent(in) :: rid, elev
    character(len=128), intent(out) :: var
    !
    integer :: id
    !
    id = nint(rid)

    if( id == 2819 )then
       var = 'ms_u'
    elseif( id == 2820 )then
       var = 'ms_v'
    !elseif( id == 2567 )then ! z_obs
    !   var = 'ms_pres'
    elseif( id == 3073 )then
       var = 'ms_tem'
    elseif( id == 3330 )then
       var = 'ms_qv'
    elseif( id == 3331 )then
       var = 'ms_rh'
    elseif( id == 3332 )then ! specific humidity
       var = 'ms_qv'
    elseif( id == 9999 )then
       var = 'sa_tppn'
    elseif( id == 14593 )then
       var = 'ss_ps'
    elseif( id == 15619 )then ! gps refractive index
       var = 'ms_gps'
    endif

    return
  end subroutine get_variable
!-----------------------------------------------------------------------
! LIST ASSIMILATED INSTRUMENT CHANNELS
!-----------------------------------------------------------------------
!SUBROUTINE set_instrument
!  tvsch = 0
  !
  ! NOAA-15 AMSU-A
  !
!  tvsname(1) = 'AA15'
!  tvsinst(1,1) = rttv_plat_noaa
!  tvsinst(2,1) = 15
!  tvsinst(3,1) = rttv_inst_amsua
!  tvsch( 1,1)  =  5
!  tvsch( 2,1)  =  6
!  tvsch( 1,1)  =  7
!  tvsch( 2,1)  =  8
!  tvsch( 5,1)  =  9
!  tvsch( 6,1)  = 10
!  tvsch( 7,1)  = 11
!  ntvsch(1)=2
  !
  ! NOAA-16 AMSU-A
  !
!  tvsname(2) = 'AA16'
!  tvsinst(1,2) = rttv_plat_noaa
!  tvsinst(2,2) = 16
!  tvsinst(3,2) = rttv_inst_amsua
!  tvsch( 1,2)  =  5
!  tvsch( 2,2)  =  6
!  tvsch( 1,2)  =  7
!  tvsch( 2,2)  =  8
!  tvsch( 5,2)  =  9
!  tvsch( 6,2)  = 10
!  tvsch( 7,2)  = 11
!  ntvsch(2)=2
  !
  ! NOAA-18 AMSU-A
  !
!  tvsname(3) = 'AA18'
!  tvsinst(1,3) = rttv_plat_noaa
!  tvsinst(2,3) = 18
!  tvsinst(3,3) = rttv_inst_amsua
!  tvsch( 1,3)  =  5
!  tvsch( 2,3)  =  6
!  tvsch( 1,3)  =  7
!  tvsch( 2,3)  =  8
!  tvsch( 5,3)  =  9
!  tvsch( 6,3)  = 10
!  tvsch( 7,3)  = 11
!  ntvsch(3)=2
  !
  ! NOAA-19 AMSU-A
  !
!  tvsname(4) = 'AA19'
!  tvsinst(1,4) = rttv_plat_noaa
!  tvsinst(2,4) = 19
!  tvsinst(3,4) = rttv_inst_amsua
!  tvsch( 1,4)  =  5
!  tvsch( 2,4)  =  6
!  tvsch( 1,4)  =  7
!  tvsch( 2,4)  =  8
!  tvsch( 5,4)  =  9
!  tvsch( 6,4)  = 10
!  tvsch( 7,4)  = 11
!  ntvsch(4)=2
!  !
!  ! METOP-2 AMSU-A
!  !
!  tvsname(5) = 'MA02'
!  tvsinst(1,5) = rttv_plat_metop2
!  tvsinst(2,5) = 19
!  tvsinst(3,5) = rttv_inst_amsua
!  tvsch(1,5)  = 4
!  tvsch(2,5)  = 5
!  tvsch(3,5)  = 6
!  tvsch(4,5)  = 8
!  tvsch(5,5)  = 9
!  tvsch(6,5)  = 10
!  tvsch(7,5)  = 11
!  tvsch(8,5)  = 12
!  tvsch(9,5)  = 13
!  ntvsch(5)=9
!  RETURN
!END SUBROUTINE set_instrument
!-----------------------------------------------------------------------
SUBROUTINE vbc_read(filename,basetime,vbc)
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  INTEGER,PARAMETER::  MD=99
  CHARACTER(*),INTENT(IN):: filename
  INTEGER,INTENT(IN):: basetime
  REAL(8),INTENT(OUT)::  vbc(maxvbc,maxtvsch,ninstrument)

  CHARACTER(len=8):: platname,instname
  INTEGER(4):: isat,ksat,plat,satn,inst,ierr,inum,iflag
  INTEGER(4):: kchan,idxchan,ichan
  INTEGER(4):: iyy,imm,idd,ihr,imn,betatime,nrec
  REAL(8)::  vbc_in(maxvbc)

  vbc = 0.0d0
  PRINT '(2A)','vbc_read reading a file: ',filename
  OPEN(MD,file=filename)
  !
  ! varbc_coef header time
  !
  !READ(MD,'(X,A8,I4,X,I2,X,I2,X,I2,X,I2)',ERR=10) &
  !READ(MD,'(X,A8,I4,X,I2,X,I2,X,I2,X,I2)') &
  !    & platname,iyy,imm,idd,ihr,imn
!  CALL nwp_ymdhm2seq(iyy,imm,idd,ihr,imn,betatime)
!  WRITE(6,'(X,A,I4,4I2.2,A)') &
!       & "varbc_coef time: ",iyy,imm,idd,ihr,imn," ***"
!  GOTO 20
!  10 WRITE(6,*) "*** ATTENTION ! : varbc_coef: NO TIME ***"
!     REWIND(MD)
!  20 CONTINUE
  !
  ! varbc_coef
  !
  nrec=0
  DO
    READ(MD,'(X,A8,I2,2X,A8,I4,8E16.8,I5)',ERR=99,END=99) &
        & platname,satn,instname,ichan,vbc_in,inum
    write(*,*) platname,satn,instname
    plat=0
    DO plat=1,nplatforms
      IF(platname==platform_name(plat)) EXIT
    END DO
    IF(instname=='imager  ') THEN
      SELECT CASE (instname)
        CASE('goes    '); inst=rttv_inst_goesi
        CASE('mtsat   '); inst=rttv_inst_mtsati
      END SELECT
    ELSE
      DO inst=0,ninst-1
        IF(instname==inst_name(inst)) EXIT
      END DO
    END IF
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
  CLOSE(MD)
  PRINT '(A,I4)','  varbc_coef records read: ',nrec
!---------------------------------------------------------------
! ^[$B3NG'MQ%m%0=PNO^[(B
!---------------------------------------------------------------
!    DO ksat=1,ninstrument
!      DO ichan=1,ntvsch(ksat)
!        iflag=0
!        DO inum=1,maxvbc
!          if(vbc(inum,ichan,ksat).ne.0.d0) iflag=1
!        END DO
!        IF(iflag==1) THEN
!          WRITE(6,'(X,A8,I2,2X,A8,I4,8E16.8)') &
!     &      platform_name(tvsinst(1,ksat)), tvsinst(2,isat), &
!     &      inst_name(tvsinst(3,ksat)), tvsch(idxchan,ksat), &
!     &      vbc(1:8,ichan,ksat)
!        END IF
!      END DO
!    END DO
!---------------------------------------------------------------
  RETURN
END SUBROUTINE vbc_read

SUBROUTINE vbc_write(filename,basetime,vbc)
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  INTEGER,PARAMETER::  MD=99
  CHARACTER(*),INTENT(IN):: filename
  INTEGER     ,INTENT(IN):: basetime
  REAL(8),INTENT(IN):: vbc(maxvbc,maxtvsch,ninstrument)

  INTEGER(4):: ksat,iflag,inum,idxchan
  INTEGER(4):: iyy,imm,idd,ihr,imn,nrec

  PRINT '(2A)','vbc_write writing a file: ',filename
  OPEN(MD,file=filename)
  !
  ! varbc_coef header time
  !
  !CALL nwp_seq2ymdhm(iyy,imm,idd,ihr,imn,basetime)
  !WRITE(MD,'(X,A8,I4,X,I2,X,I2,X,I2,X,I2)') &
  !     & '#BTIME: ',iyy,imm,idd,ihr,imn
  !
  ! varbc_coef
  !
  nrec = 0
  DO ksat=1,ninstrument
    DO idxchan=1,ntvsch(ksat)
      iflag=0
      DO inum=1,maxvbc
        IF(vbc(inum,idxchan,ksat).NE.0.d0) iflag=1
      END DO
      IF(iflag==1) THEN
        WRITE(MD,'(X,A8,I2,2X,A8,I4,8E16.8,I5)') &
   &      platform_name(tvsinst(1,ksat)), tvsinst(2,ksat), &
   &      inst_name(tvsinst(3,ksat)), tvsch(idxchan,ksat), &
   &      vbc(1:8,idxchan,ksat), ntvschan(idxchan,ksat)
        nrec = nrec+1
      END IF
    END DO
  END DO
  CLOSE(MD)
  PRINT '(A,I4)','  varbc_coef records written: ',nrec
  RETURN
END SUBROUTINE vbc_write
!---------------------------------------------------------------
SUBROUTINE vbc_scan_read(filename,basetime,vbc_scan)
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  INTEGER,PARAMETER::  MD=99
  CHARACTER(*),INTENT(IN):: filename
  INTEGER,INTENT(IN):: basetime
  REAL(8),INTENT(OUT)::  vbc_scan(maxfoot,maxtvsch,ninstrument)

  CHARACTER(len=8):: platname,instname
  INTEGER(4):: isat,ksat,plat,satn,inst,ierr,inum,iflag
  INTEGER(4):: kchan,idxchan,ichan
  INTEGER(4):: iyy,imm,idd,ihr,imn,betatime,nrec
  REAL(8)::  vbc_in(maxfoot)

  vbc_scan = 0.0d0
  PRINT '(2A)','vbc_read reading a file: ',filename
  OPEN(MD,file=filename)
  !
  ! varbc_coef
  !
  nrec=0
  DO
    READ(MD,'(X,A8,I2,2X,A8,I4,30E16.8)',ERR=99,END=99) &
        & platname,satn,instname,ichan,vbc_in
    write(*,*) platname,satn,instname,ichan
    plat=0
    DO plat=1,nplatforms
      IF(platname==platform_name(plat)) EXIT
    END DO
    IF(instname=='imager  ') THEN
      SELECT CASE (instname)
        CASE('goes    '); inst=rttv_inst_goesi
        CASE('mtsat   '); inst=rttv_inst_mtsati
      END SELECT
    ELSE
      DO inst=0,ninst-1
        IF(instname==inst_name(inst)) EXIT
      END DO
    END IF
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
  CLOSE(MD)
  PRINT '(A,I4)','  varbc_scan records read: ',nrec
  RETURN
END SUBROUTINE vbc_scan_read
!---------------------------------------------------------------
SUBROUTINE vbc_scan_write(filename,basetime,vbc_scan)
  USE rttov_const,ONLY: nplatforms,platform_name,ninst,inst_name
  IMPLICIT NONE
  INTEGER,PARAMETER::  MD=99
  CHARACTER(*),INTENT(IN):: filename
  INTEGER     ,INTENT(IN):: basetime
  !REAL(8),INTENT(IN):: vbc(maxvbc,maxtvsch,ninstrument)
  REAL(8),INTENT(IN)::  vbc_scan(maxfoot,maxtvsch,ninstrument)

  INTEGER(4):: ksat,iflag,inum,idxchan
  INTEGER(4):: iyy,imm,idd,ihr,imn,nrec

  PRINT '(2A)','vbc_write writing a file: ',filename
  OPEN(MD,file=filename)
  nrec = 0
  DO ksat=1,ninstrument
    DO idxchan=1,ntvsch(ksat)
      iflag=0
      DO inum=1,maxfoot
        IF(vbc_scan(inum,idxchan,ksat).NE.0.d0) iflag=1
      END DO
      IF(iflag==1) THEN
        WRITE(MD,'(X,A8,I2,2X,A8,I4,30E16.8)') &
   &      platform_name(tvsinst(1,ksat)), tvsinst(2,ksat), &
   &      inst_name(tvsinst(3,ksat)), tvsch(idxchan,ksat), &
   &      vbc_scan(1:maxfoot,idxchan,ksat)
        nrec = nrec+1
      END IF
    END DO
  END DO
  CLOSE(MD)
  PRINT '(A,I4)','  varbc_coef records written: ',nrec
  RETURN
END SUBROUTINE vbc_scan_write
!-----------------------------------------------------------------------
! Data Assimilation for VARBC
!-----------------------------------------------------------------------
SUBROUTINE das_vbc(hx,pred,vbcf,vbca)
  USE common_mtx
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: hx(:,:,:,:)
  REAL(8),INTENT(IN) :: pred(:,:,:,:,:)
  REAL(8),INTENT(INOUT) :: vbcf(maxvbc,maxtvsch,ninstrument)
  REAL(8),INTENT(OUT)   :: vbca(maxvbc,maxtvsch,ninstrument)
  REAL(8) :: a(maxvbc,maxvbc)
  REAL(8) :: b(maxvbc)
  REAL(8) :: ainv(maxvbc,maxvbc)
  !INTEGER:: ntvschan(maxtvsch,ninstrument)
  INTEGER:: ntvschan1(maxtvsch,ninstrument)
  INTEGER:: i,j,k,n,islot,nn

  PRINT *,'Hello from das_vbc'

!  IF(ntvs == 0) THEN
!    PRINT *,'No radiance data: das_vbc skipped..'
!    vbca = vbcf
!    RETURN
!  END IF

  vbca = 0.0d0
  DO k=1,ninstrument
    !DO j=1,maxtvsch
    DO j=1,ntvsch(k)
      !
      ! Parallel processing
      !
      !IF(MOD(j+maxtvsch*(k-1)-1,nprocs) /= myrank) CYCLE
      !
      ! DATA NUMBER
      !
      ichan=tvsch(j,k)
      ntvschan(j,k) = SUM(tvsqc(ichan,:,k,:))
      !ntvschan(j,k) = SUM(tvsqc(j,:,k,:))
      IF(msw_vbc .AND. ntvschan(j,k) /= 0 ) THEN
        PRINT '(3A,I3,A,I6)',' >> VBC executed for instrument,channel,ntvsl:',&
                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
        CALL vbc_local(j,k,ntvschan(j,k),hx,pred,a,b)
        CALL mtx_inv(maxvbc,a,ainv)
        vbca(:,j,k) = vbcf(:,j,k)
        DO n=1,maxvbc
          vbca(:,j,k) = vbca(:,j,k) - ainv(:,n)*b(n) !ATTN: sign for beta
        END DO
      ELSE
        PRINT '(3A,I3,A,I6)',' !! NO VBC executed for instrument,channel,ntvsl:',&
                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
        vbca(:,j,k) = vbcf(:,j,k)
      END IF
    END DO
  END DO
  vbcf = vbca
!  ntvschan1 = ntvschan
!  DEALLOCATE(hx,pred)
!  n = maxvbc*maxtvsch*ninstrument
!  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!  CALL MPI_ALLREDUCE(vbcf,vbca,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,j)
!  n = maxtvsch*ninstrument
!  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!  CALL MPI_ALLREDUCE(ntvschan1,ntvschan,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,j)

  RETURN
END SUBROUTINE das_vbc
!-----------------------------------------------------------------------
!  (in) ichan: channnel
!  (in) iinst: sensor
!  (out) a = B_beta^-1 + p R^-1 p^T
!  (out) b = p R^-1 d
!-----------------------------------------------------------------------
SUBROUTINE vbc_local(ic,iinst,ntvsl,hx,pred,a,b)
  IMPLICIT NONE
  INTEGER,PARAMETER :: msw=1
  INTEGER,PARAMETER :: nmin=400
  !INTEGER,INTENT(IN) :: ichan,iinst,ntvsl
  INTEGER,INTENT(IN) :: ic,iinst,ntvsl
  REAL(8),INTENT(IN) :: hx(maxtvsch,maxtvsprof,ninstrument,nslots)
  REAL(8),INTENT(IN) :: pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots)
  REAL(8),INTENT(OUT) :: a(maxvbc,maxvbc)
  REAL(8),INTENT(OUT) :: b(maxvbc)
  REAL(8) :: dep,dep0
  REAL(8) :: bias,bias0
  REAL(8) :: r,tmp
  INTEGER:: islot, iprof, i,j,n
  INTEGER:: ichan

  a = 0.0d0
  b = 0.0d0
  dep = 0.0d0
  dep0 = 0.0d0
  bias = 0.0d0
  bias0 = 0.0d0
  n = 0
  ichan=tvsch(ic,iinst)
  DO islot=1,nslots
    DO iprof=1,maxtvsprof
      IF(tvsqc(ichan,iprof,iinst,islot)/=1) CYCLE
      !
      ! R
      !
      r = tvserr(ichan,iprof,iinst,islot)**2
      !
      ! p R^-1 p^T
      !
      DO j=1,maxvbc
        DO i=1,maxvbc
          a(i,j) = a(i,j) &
               & + pred(i,ic,iprof,iinst,islot) &
               & * pred(j,ic,iprof,iinst,islot) / r
               !& + pred(i,ichan,iprof,iinst,islot) &
               !& * pred(j,ichan,iprof,iinst,islot) / r
        END DO
      END DO
      !
      ! B_beta^-1
      !
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
      !
      ! p R^-1 d
      !
      b(:) = b(:) + pred(:,ic,iprof,iinst,islot) / r &
                & *(tvsdat(ichan,iprof,iinst,islot)-hx(ic,iprof,iinst,islot))
      !b(:) = b(:) + pred(:,ichan,iprof,iinst,islot) / r &
      !          & *(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))
      write(*,*) islot, iprof, tvsdat(ichan,iprof,iinst,islot), hx(ic,iprof,iinst,islot)
      bias = bias+tvsdat(ichan,iprof,iinst,islot)-hx(ic,iprof,iinst,islot)
      dep = dep+(tvsdat(ichan,iprof,iinst,islot)-hx(ic,iprof,iinst,islot))**2
      bias0= bias0+tvsdep(ichan,iprof,iinst,islot)
      dep0= dep0+tvsdep(ichan,iprof,iinst,islot)**2
      n = n+1
    END DO
  END DO

  dep = SQRT(dep / REAL(n,kind=8))
  dep0 = SQRT(dep0 / REAL(n,kind=8))
  bias = bias / REAL(n,kind=8)
  bias0 = bias0 / REAL(n,kind=8)
  PRINT '(2A,I3,4F12.4)',' >> D monit: ',tvsname(iinst),tvsch(ichan,iinst),bias0,bias,dep0,dep

  RETURN
END SUBROUTINE vbc_local


end program prg_vbc


