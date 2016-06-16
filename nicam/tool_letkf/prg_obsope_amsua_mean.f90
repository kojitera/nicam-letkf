!-------------------------------------------------------------------------------
!
!+  Program obsope
!
!-------------------------------------------------------------------------------
program prg_obsope_amsua_mean

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
  !use common_tvs_nicam
  use common_tvs_nicam, only : &
       set_instrument, &
       tvsname,        &
       tvsinst,        &
       tvsch,          &
       ntvsch,         &
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
  character(ADM_MAXFNAME), save :: input_fname(100) = ''
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
  character(LEN=FIO_HLONG)  :: infile_header(flim) = ''
  character(LEN=FIO_HLONG)  :: infile(flim)        = ''
  integer                   :: step_str            = 1
  integer                   :: step_end            = max_nstep
  character(LEN=FIO_HLONG)  :: outfile_dir         = '.'
  character(LEN=FIO_HLONG)  :: outfile_prefix      = ''
  integer                   :: outfile_rec         = 1
  logical                   :: output_grads        = .true.
  logical                   :: datainfo_nodep_pe   = .false.   !   <- can be .true. if data header do not depend on pe.
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''
  character(LEN=FIO_HSHORT) :: large_memory_var(max_nvar) = '' ! [add] 13-04-18
  logical                   :: help = .false.

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
  integer, allocatable, save :: agl(:)         ! above ground level = 1, if not = 0
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
  real(8), allocatable :: saza(:,:,:), saaz(:,:,:), soza(:,:,:), soaz(:,:,:)
  real(8), allocatable :: land(:), elev(:)
  real(8), allocatable :: said(:,:,:)
  real(8), allocatable, save :: lsql(:,:,:)

  logical, allocatable, save :: inprc(:)
  real(8), allocatable, save :: veg(:,:,:)
  real(8), allocatable, save :: veg_pl(:,:,:)
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
  integer :: s, n, l, i, nn, k
  integer :: fid
  integer :: ofid, ofid2
  integer :: ierr
  !
  INTEGER,PARAMETER :: id_temp_nicam=1
  INTEGER,PARAMETER :: id_qvap_nicam=2
  INTEGER,PARAMETER :: id_qcld_nicam=3
  INTEGER,PARAMETER :: id_pres_nicam=4
  INTEGER,PARAMETER :: id_tsfc_nicam=5
  INTEGER,PARAMETER :: id_qv2m_nicam=6
  INTEGER,PARAMETER :: id_surp_nicam=7
  INTEGER,PARAMETER :: id_u10m_nicam=8
  INTEGER,PARAMETER :: id_v10m_nicam=9
  INTEGER,PARAMETER :: id_te2m_nicam=10
  INTEGER,PARAMETER :: id_cldw_nicam=11

  !INTEGER,PARAMETER :: id_temp_nicam=1
  !INTEGER,PARAMETER :: id_qvap_nicam=2
  !INTEGER,PARAMETER :: id_qcld_nicam=3
  !INTEGER,PARAMETER :: id_pres_nicam=4
  !INTEGER,PARAMETER :: id_tsfc_nicam=5
  !INTEGER,PARAMETER :: id_u10m_nicam=6
  !INTEGER,PARAMETER :: id_v10m_nicam=7
  !INTEGER,PARAMETER :: id_cldw_nicam=8
  !INTEGER,PARAMETER :: id_surp_nicam=9
  !INTEGER,PARAMETER :: id_qv2m_nicam=10
  !INTEGER,PARAMETER :: id_te2m_nicam=11
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

  INTEGER :: num_satellite
  INTEGER :: nobs
  INTEGER,PARAMETER :: id_NOAA15=206
  INTEGER,PARAMETER :: id_NOAA16=207
  INTEGER,PARAMETER :: id_NOAA17=208
  INTEGER,PARAMETER :: id_NOAA18=209
  INTEGER,PARAMETER :: id_NOAA19=223
  INTEGER,PARAMETER :: id_METOP2=4

  REAL(8),PARAMETER :: gross_error=10.0d0
  INTEGER,PARAMETER :: nslots=1
  INTEGER,PARAMETER :: ninstrument=5
  INTEGER,PARAMETER :: maxtvsch=15
  INTEGER,PARAMETER :: maxvbc=8
  INTEGER,PARAMETER :: maxfoot=30
  INTEGER,SAVE :: maxtvsprof
  INTEGER,SAVE :: maxtvsfoot
  INTEGER,SAVE :: ntvsprof(ninstrument,nslots)
  INTEGER :: islot, ic, ifoot

  !CHARACTER(4),SAVE :: tvsname(ninstrument)
  CHARACTER(256) :: rttovcoef_fname(ninstrument) = ''
  !INTEGER,SAVE :: tvsinst(3,ninstrument)
  !INTEGER,SAVE :: tvsch(maxtvsch,ninstrument)
  !INTEGER,SAVE :: ntvsch(ninstrument)
  INTEGER,SAVE :: ntvs(ninstrument)
  REAL(8),ALLOCATABLE,SAVE :: vbc_pred(:,:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsland(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsdat(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsdep(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsqc(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvserr(:,:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvselev(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvslon(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvslat(:,:,:)
  REAL(8),ALLOCATABLE,SAVE :: tvsfoot(:,:,:)
  INTEGER,ALLOCATABLE,SAVE :: ntvsfoot(:,:)

  REAL(8),SAVE :: vbcf(maxvbc,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbca(maxvbc,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbcf_scan(maxfoot,maxtvsch,ninstrument)
  REAL(8),SAVE :: vbca_scan(maxfoot,maxtvsch,ninstrument)

  REAL(Kind=8), allocatable :: bt(:,:,:,:)
  REAL(Kind=8), allocatable :: hx(:,:,:,:)
  REAL(Kind=8), allocatable :: bt_tmp(:,:,:,:)
  REAL(Kind=8), allocatable :: tran(:,:,:,:,:)
  REAL(Kind=8), allocatable :: tran_tmp(:,:,:,:,:)
  REAL(Kind=8), allocatable :: weight(:,:,:,:,:)
  !REAL(Kind=8), allocatable :: dtau(:)
  INTEGER, allocatable :: weight_maxlev(:,:,:,:)
  INTEGER :: ichan, ilev, iobs, iobs1
  !INTEGER, PARAMETER :: nlev=94
  INTEGER, PARAMETER :: nlev=38
  INTEGER :: nlev_RTTOV=nlev
  !INTEGER :: nlev_RTTOV=nlev+6
  REAL(Kind=8) :: p_tmp(nlev)
  REAL(Kind=8) :: tmp_lev
  REAL(Kind=8) :: iwlr
  INTEGER,SAVE :: ntvschan(maxtvsch,ninstrument)

  ! ico data
  integer              :: PALL_global
  integer              :: LALL_global
  integer              :: LALL_local

  ! for MPI
  integer          :: prc_nlocal
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

  !!!character(LEN=FIO_HLONG) :: fname
  character(LEN=20)        :: tmpl
  !!!character(LEN=16)        :: gthead(64)
  integer(8)               :: nowsec
  !!!integer(8)               :: recsize ! [mod] 12-04-19 H.Yashiro
  integer                  :: kmax, num_of_step, step, date_str(6)

  logical :: addvar
  integer :: did
  integer :: t, v, p, pp  ! loop-index
  real(8) :: fac1, fac2, fac3, fac_sum
  !!!integer :: ks, ke

  ! ico data
  real(4), allocatable :: data4allrgn(:)
  real(8), allocatable :: data8allrgn(:)
  real(4), allocatable :: icodata4(:,:,:)
  real(4), allocatable :: icodata(:,:,:,:)

  real(4) :: wk(34)
  real(4), allocatable :: wk2(:)
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
  real(8), allocatable :: obsdata(:,:,:)
  real(8), allocatable :: obsdata_out(:,:,:)
  real(kind=JPRB), allocatable :: obsdata_jprb(:,:,:,:,:)
  real(kind=JPRB) :: USSTD_jprb(6,2)

  integer :: sobs, eobs, obsint, prc_per_inst

  logical :: ocheck=.true.
  real(8) :: tmp_time(100)
  real(8) :: sum_time(100)

  character(4) :: cfile='inst'

  integer :: smem, emem, imem
  character(6) :: csmem, cemem, cimem

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
       input_fname,   &
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
                    help


  sum_time(:)=0.0d0
  tmp_time(1)=MPI_WTIME()

  tmp_time(21)=MPI_WTIME()
  call ADM_proc_init(ADM_MULTI_PRC)
  tmp_time(22)=MPI_WTIME()
  call ADM_setup('obsope.cnf')
  tmp_time(23)=MPI_WTIME()
  call COMM_setup
  tmp_time(24)=MPI_WTIME()
  call CNST_setup
  tmp_time(25)=MPI_WTIME()
  call GRD_setup
  tmp_time(26)=MPI_WTIME()
  call GMTR_setup
  tmp_time(27)=MPI_WTIME()
  call OPRT_setup
  tmp_time(28)=MPI_WTIME()
  call VMTR_setup
  tmp_time(29)=MPI_WTIME()
  call readoption !! set fmax, infile
  tmp_time(30)=MPI_WTIME()

  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(1))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(1))
  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(2))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(2))
  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(3))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(3))
  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(4))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(4))
  FLUSH(ADM_LOG_FID)

  rewind(ADM_CTL_FID)
  read(ADM_CTL_FID, nml=obsope_param, iostat=ierr)
  !rttovcoef_fname(1) ='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_15_amsua.dat'
  !rttovcoef_fname(2) ='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_16_amsua.dat'
  !rttovcoef_fname(3) ='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_18_amsua.dat'
  !rttovcoef_fname(4) ='/scratch/ra000015/koji/LETKF/run_4d/letkf/rttov/rtcoef_noaa_19_amsua.dat'

  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(1))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(1))
  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(2))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(2))
  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(3))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(3))
  write(ADM_LOG_FID,*) "trim(rttovcoef_fname(4))"
  write(ADM_LOG_FID,*) trim(rttovcoef_fname(4))
  do islot = 1, nslots
    write(ADM_LOG_FID,*) islot, trim(input_fname(islot))
  end do
  FLUSH(ADM_LOG_FID)

  fid=40
  do islot = 1, nslots
    FLUSH(ADM_LOG_FID)
    open(fid, file=trim(input_fname(islot)), form='unformatted', access='sequential', &
         status='old', iostat=ierr)
    if(ierr /= 0) then
       write(ADM_LOG_FID,*) 'ERROR for opening the file', trim(input_fname(islot))
       write(ADM_LOG_FID,*) 'STOP!'
       FLUSH(ADM_LOG_FID)
       call ADM_proc_stop
    end if

    read(fid) num_satellite
    read(fid) nsite
    read(fid) ntvsprof(1:ninstrument,islot)
    write(ADM_LOG_FID,*) islot, ntvsprof(1:ninstrument,islot)
    close(fid)
  end do

  do nn = 1, ninstrument
    ntvs(nn)=sum(ntvsprof(nn,1:nslots))
  end do

  maxtvsprof=maxval(ntvsprof)
  maxtvsfoot=maxval(nfootp)
  nsite=sum(ntvs)

  write(ADM_LOG_FID,*) 'num_satellite', num_satellite
  write(ADM_LOG_FID,*) 'nsite', nsite
  write(ADM_LOG_FID,*) 'ntvs', ntvs(1:ninstrument)
  write(ADM_LOG_FID,*) 'maxtvsprof', maxtvsprof
  FLUSH(ADM_LOG_FID)

  allocate( tvsland(maxtvsprof,ninstrument,nslots) )
  allocate( tvsdat(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvserr(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvsdep(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvsqc(maxtvsch,maxtvsprof,ninstrument,nslots) )
  allocate( tvselev(maxtvsprof,ninstrument,nslots) )
  allocate( tvslat(maxtvsprof,ninstrument,nslots) )
  allocate( tvslon(maxtvsprof,ninstrument,nslots) )
  allocate( tvsfoot(maxtvsprof,ninstrument,nslots) )
  allocate( ntvsfoot(maxtvsfoot,ninstrument) )

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
  allocate(obsdata(nlev,nsite,12))
  allocate(obsdata_out(nlev,nsite,12))
  allocate(obsdata_jprb(nlev+6,maxtvsprof,ninstrument,12,nslots))
 ! allocate(obsdata_out(nsite))

  allocate(obsnum(nsite))
  allocate(lon(nsite))
  allocate(lat(nsite))
  allocate(agl(nsite))
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

  allocate( land(nsite) )
  allocate( elev(nsite) )
  allocate( lsql(maxtvsprof,ninstrument,nslots) )
  allocate( saza(maxtvsprof,ninstrument,nslots) )
  allocate( saaz(maxtvsprof,ninstrument,nslots) )
  allocate( soza(maxtvsprof,ninstrument,nslots) )
  allocate( soaz(maxtvsprof,ninstrument,nslots) )
  allocate( said(maxtvsprof,ninstrument,nslots) )

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

  call set_instrument

  write(ADM_LOG_FID,*) 'ntvsch(n)=', ntvsch(1:ninstrument)
  write(ADM_LOG_FID,*) 'tvsch(1)=', tvsch(1:ntvsch(1),1)
  write(ADM_LOG_FID,*) 'tvsch(2)=', tvsch(1:ntvsch(2),2)
  write(ADM_LOG_FID,*) 'tvsch(3)=', tvsch(1:ntvsch(3),3)
  write(ADM_LOG_FID,*) 'tvsch(4)=', tvsch(1:ntvsch(4),4)
  write(ADM_LOG_FID,*) 'tvsch(5)=', tvsch(1:ntvsch(5),5)

  call calendar_yh2ss( time_tmp, idate )
  time_obs(:)=dnint(time_tmp)

  tmp_time(2)=MPI_WTIME()

  i=1
  do islot = 1, nslots
    open(fid, file=trim(input_fname(islot)), form='unformatted', access='sequential', &
         status='old', iostat=ierr)
    read(fid)
    read(fid) nobs
    read(fid)
    do nn = 1, ninstrument
      do n = 1, ntvsprof(nn,islot)
        read(fid) wk
        lat(i)              = dble(wk(7))
        tvslat(n,nn,islot)  = dble(wk(7))
        lon(i)              = dble(wk(8))
        tvslon(n,nn,islot)  = dble(wk(8))
        said(n,nn,islot)    = dble(wk(9))
        tvsfoot(n,nn,islot) = dble(wk(11))
        lsql(n,nn,islot)    = dble(wk(12))
        saza(n,nn,islot)    = dble(wk(13))
        soza(n,nn,islot)    = dble(wk(14))
        elev(i)             = dble(wk(15))
        tvselev(n,nn,islot) = dble(wk(15))
        soaz(n,nn,islot)    = dble(wk(17))
        saaz(n,nn,islot)    = dble(wk(18))
        odat(1:15,i)        = dble(wk(19:33)) !! 2014.07.17 fix
        obserr(:,i) = 0.35
        i=i+1
      enddo
    enddo
    close(fid)
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

  i=1
  do islot = 1, nslots
    do nn = 1, ninstrument
      do n = 1, ntvsprof(nn,islot)
        tvsland(n,nn,islot)=lsql(n,nn,islot)
        tvsdat(:,n,nn,islot)=odat(:,i)
        i=i+1
      end do
    end do
  end do
  tvserr(:,:,:,:)=0.35d0
  do nn = 1, ninstrument
    write(ADM_LOG_FID,'(i5,10f10.3)') nn, (tvsdat(tvsch(ic,nn),1,nn,1),ic=1,ntvsch(nn))
  end do

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
     !write(ADM_LOG_FID,*) 'max_num_latlon: ', ADM_prc_me, max_num_latlon
     !
     if(sum_max_num_latlon /= nsite) then
        write(ADM_LOG_FID,*) 'sum_max_num_latlon /= nsite', sum_max_num_latlon, nsite
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
!
!
!
  infile_header(1)=trim(infile(1))
  do p = pstr, pend
     pp = p - pstr + 1

     if (complete) then ! all region
        infname = trim(infile(1))//'.rgnall'
     else
        call fio_mk_fname(infname,trim(infile_header(1)),'pe',p-1,6)
        !call fio_mk_fname(infname,trim(infile(1)),'pe',p-1,6)
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

  !write(ADM_LOG_FID,*) '########## Variable List ########## '
  !write(ADM_LOG_FID,*) 'ID |NAME            |STEPS|Layername       |START FROM         |DT [sec]'
  do v = 1, nvar
     call calendar_ss2yh( date_str(:), real(var_time_str(v),kind=8) )
     !write(tmpl,'(I4.4,"/",I2.2,"/",I2.2,1x,I2.2,":",I2.2,":",I2.2)') date_str(:)
     !write(ADM_LOG_FID,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A19,A1,I8)') &
     !         v,'|',var_name(v),'|',var_nstep(v),'|',var_layername(v),'|', tmpl,'|', var_dt(v)
  enddo

  !write(ADM_LOG_FID,*) '*** convert start : PaNDa format to lat-lon data'
  FLUSH(ADM_LOG_FID)
  tmp_time(6)=MPI_WTIME()
  sum_time(1)=sum_time(1)+(tmp_time(6)-tmp_time(5))

  !#########################################################
  !--- start weighting summation

  variable_loop: do v = 1, nvar

     kmax = var_nlayer(v)
     !write(ADM_LOG_FID,*) '       v=',v
     !write(ADM_LOG_FID,*) '    kmax=',kmax
     !write(ADM_LOG_FID,*) 'var_name=',trim(var_name(v))

     num_of_step = min(step_end,var_nstep(v)) - step_str + 1  ! [mov] 13-04-18

     step_loop: do t = 1, num_of_step

        nowsec = var_time_str(v) + (t-1)*var_dt(v)

        step = t-1 + step_str

         allocate( data4allrgn(ADM_gall*kmax*LALL_local) )
         allocate( data8allrgn(ADM_gall*kmax*LALL_local) )
         allocate( icodata4   (ADM_gall,kmax,LALL_local) )
         data4allrgn(:)  = CNST_UNDEF4
         data8allrgn(:)  = CNST_UNDEF
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
           write(ADM_LOG_FID,*) trim(var_name(v)), minval(icodata4(:,:,:)), maxval(icodata4(:,:,:))
           FLUSH(ADM_LOG_FID)

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
                   end do
                   if(trim(var_name(v))=='ms_pres' .or. trim(var_name(v))=='ss_ps') then
                      write(ADM_LOG_FID,*) trim(var_name(v)), minval(obsdata(1:kmax,i,v)), maxval(obsdata(1:kmax,i,v))
                      flush(ADM_LOG_FID)
                      obsdata(1:kmax,i,v)=exp(obsdata(1:kmax,i,v))
                      !obsdata(:,i,v)=exp(obsdata(:,i,v))
                   else if(trim(var_name(v))=='ms_qv' &
                      !.or. trim(var_name(v))=='ms_qc' &
                      .or. trim(var_name(v))=='ss_q2m' ) then
                      obsdata(1:kmax,i,v)=obsdata(1:kmax,i,v)*q2ppmv
                      !obsdata(:,i,v)=obsdata(:,i,v)*q2ppmv
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
     !write(ADM_LOG_FID,'(i6, 5f12.5)') v, (obsdata(1,k,v),k=1,5)
  end do

  oqc(:,:)=0

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(7)=MPI_WTIME()
  sum_time(2)=sum_time(2)+(tmp_time(7)-tmp_time(6))

  call MPI_ALLREDUCE(obsdata, obsdata_out, nlev*nsite*nvar, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)

  tmp_time(8)=MPI_WTIME()
  sum_time(3)=sum_time(3)+(tmp_time(8)-tmp_time(7))
 
  do v = 1, nvar
     write(ADM_LOG_FID,'(i6, 5f12.5)') v, (obsdata_out(1,k,v),k=1,5)
  end do
  FLUSH(ADM_LOG_FID)

  i=1
  do islot = 1, nslots
    do nn = 1, ninstrument
      do n = 1, ntvsprof(nn,islot)
        obsdata_jprb(1:nlev,n,nn,:,islot)=real(obsdata_out(1:nlev,i,:),kind=jprb)
        i=i+1
      end do
    end do
  end do
  !
  ! USE U.S STANDARD ATMOSPHERE
  !
  USSTD_jprb(1,1)=1.491
  USSTD_jprb(2,1)=0.7978
  USSTD_jprb(3,1)=0.2196
  USSTD_jprb(4,1)=0.05221
  USSTD_jprb(5,1)=0.01052
  USSTD_jprb(6,1)=0.001836

  USSTD_jprb(1,2)=264.16
  USSTD_jprb(2,2)=270.65
  USSTD_jprb(3,2)=247.02
  USSTD_jprb(4,2)=219.59
  USSTD_jprb(5,2)=198.64
  USSTD_jprb(6,2)=186.87

  !do islot = 1, nslots
  !  do nn = 1, ninstrument
  !    do n = 1, ntvsprof(nn,islot)
  !      obsdata_jprb(nlev+1:nlev+6,n,nn,id_pres_nicam,islot)=USSTD_jprb(1:6,1)
  !      obsdata_jprb(nlev+1:nlev+6,n,nn,id_temp_nicam,islot)=USSTD_jprb(1:6,2)
  !      obsdata_jprb(nlev+1:nlev+6,n,nn,id_qvap_nicam,islot)=obsdata_jprb(nlev,n,nn,id_qvap_nicam,islot)
  !      obsdata_jprb(nlev+1:nlev+6,n,nn,id_qvap_nicam,islot)=0.0
  !    end do
  !  end do
  !end do  


  allocate(     bt(15,maxtvsprof,ninstrument,nslots) )
  allocate( bt_tmp(15,maxtvsprof,ninstrument,nslots) )
  allocate(     tran(nlev_RTTOV,15,maxtvsprof,ninstrument,nslots) )
  allocate( tran_tmp(nlev_RTTOV,15,maxtvsprof,ninstrument,nslots) )
  allocate( hx(maxtvsch,maxtvsprof,ninstrument,nslots) )
  !allocate( dtau(nlev) )

  bt(:,:,:,:)=0.0d0
  tran(:,:,:,:,:)=0.0d0
  bt_tmp(:,:,:,:)=0.0d0
  tran_tmp(:,:,:,:,:)=0.0d0

  lon(:)=lon(:)/CNST_PI*180.0d0
  lat(:)=lat(:)/CNST_PI*180.0d0

  do islot = 1, nslots
    do nn = 1, ninstrument
      do n =  1, ntvsprof(nn,islot)
        if(tvslon(n,nn,islot) < 0.0d0) then
          tvslon(n,nn,islot)=tvslon(n,nn,islot)+360.d0
        end if
      end do
    end do
  end do

  !nlev=var_nlayer(id_pres_nicam)
  prc_per_inst=ADM_prc_all/ninstrument

  WRITE(ADM_LOG_FID,*) "prc_per_inst", prc_per_inst

  if( ADM_prc_me <= prc_per_inst ) then
    islot=1
    if(ntvsprof(1,islot) /= 0 ) then
      obsint=ntvsprof(1,islot)/prc_per_inst
      if(ntvsprof(1,islot)<prc_per_inst) obsint=1
      sobs=obsint*(ADM_prc_me-1)+1
      eobs=obsint*ADM_prc_me
      if(ADM_prc_me==prc_per_inst) eobs=ntvsprof(1,islot)
      WRITE(ADM_LOG_FID,*) 'NOAA-15', ntvsprof(1,islot), ADM_prc_me, sobs, eobs
      FLUSH(ADM_LOG_FID)
      if(ntvsprof(1,islot)>=prc_per_inst .or. &
         ntvsprof(1,islot)< prc_per_inst .and. ADM_prc_me <= ntvsprof(1,islot)) then
      call AMSUA_fwd(nlev_RTTOV, eobs-sobs+1, rttovcoef_fname(1),&
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 1, id_pres_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 1, id_temp_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 1, id_qvap_nicam, islot), &
           obsdata_jprb(              1, sobs:eobs, 1, id_tsfc_nicam, islot), &
           obsdata_jprb(              1, sobs:eobs, 1, id_qv2m_nicam, islot), &
           obsdata_jprb(              1, sobs:eobs, 1, id_surp_nicam, islot), &
           obsdata_jprb(              1, sobs:eobs, 1, id_u10m_nicam, islot), &
           obsdata_jprb(              1, sobs:eobs, 1, id_v10m_nicam, islot), &
           soza(    sobs:eobs,1,islot), &
           soaz(    sobs:eobs,1,islot), &
           saza(    sobs:eobs,1,islot), &
           saaz(    sobs:eobs,1,islot), &
           tvselev( sobs:eobs,1,islot)/1000.0d0, &
           tvslon(  sobs:eobs,1,islot), &
           tvslat(  sobs:eobs,1,islot), &
           tvsland( sobs:eobs,1,islot), &
           bt(    :,sobs:eobs,1,islot), &
           tran(:,:,sobs:eobs,1,islot) )
      end if
    end if
  end if

  if( ADM_prc_me > prc_per_inst .and. ADM_prc_me <= 2*prc_per_inst ) then
    islot=1
    if(ntvsprof(2,islot) /= 0 ) then
      obsint=ntvsprof(2,islot)/prc_per_inst
      if(ntvsprof(2,islot)<prc_per_inst) obsint=1
      sobs=obsint*(ADM_prc_me-prc_per_inst-1)+1
      eobs=obsint*(ADM_prc_me-prc_per_inst)
      if(ADM_prc_me==2*prc_per_inst) eobs=ntvsprof(2,islot)
      WRITE(ADM_LOG_FID,*) 'NOAA-16', ntvsprof(2,islot), ADM_prc_me, sobs, eobs
      FLUSH(ADM_LOG_FID)
      if(ntvsprof(2,islot)>=prc_per_inst .or. &
         ntvsprof(2,islot)< prc_per_inst .and. ADM_prc_me-prc_per_inst <= ntvsprof(2,islot)) then
      call AMSUA_fwd(nlev_RTTOV, eobs-sobs+1, rttovcoef_fname(2),&
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 2, id_pres_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 2, id_temp_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 2, id_qvap_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 2, id_tsfc_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 2, id_qv2m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 2, id_surp_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 2, id_u10m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 2, id_v10m_nicam, islot), &
           soza(    sobs:eobs,2,islot), &
           soaz(    sobs:eobs,2,islot), &
           saza(    sobs:eobs,2,islot), &
           saaz(    sobs:eobs,2,islot), &
           tvselev( sobs:eobs,2,islot)/1000.0d0, &
           tvslon(  sobs:eobs,2,islot), &
           tvslat(  sobs:eobs,2,islot), &
           tvsland( sobs:eobs,2,islot), &
           bt(    :,sobs:eobs,2,islot), &
           tran(:,:,sobs:eobs,2,islot) )
      end if
    end if
  end if

  if( ADM_prc_me > 2*prc_per_inst .and. ADM_prc_me <= 3*prc_per_inst ) then
    islot=1
    if(ntvsprof(3,islot) /= 0 ) then
      obsint=ntvsprof(3,islot)/prc_per_inst
      if(ntvsprof(3,islot)<prc_per_inst) obsint=1
      sobs=obsint*(ADM_prc_me-2*prc_per_inst-1)+1
      eobs=obsint*(ADM_prc_me-2*prc_per_inst)
      if(ADM_prc_me==3*prc_per_inst) eobs=ntvsprof(3,islot)
      WRITE(ADM_LOG_FID,*) 'NOAA-18', ntvsprof(3,islot), ADM_prc_me, sobs, eobs
      FLUSH(ADM_LOG_FID)
      if(ntvsprof(3,islot)>=prc_per_inst .or. &
         ntvsprof(3,islot)< prc_per_inst .and. ADM_prc_me-2*prc_per_inst <= ntvsprof(3,islot)) then
      call AMSUA_fwd(nlev_RTTOV, eobs-sobs+1, rttovcoef_fname(3),&
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 3, id_pres_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 3, id_temp_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 3, id_qvap_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 3, id_tsfc_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 3, id_qv2m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 3, id_surp_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 3, id_u10m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 3, id_v10m_nicam, islot), &
           soza(    sobs:eobs,3,islot), &
           soaz(    sobs:eobs,3,islot), &
           saza(    sobs:eobs,3,islot), &
           saaz(    sobs:eobs,3,islot), &
           tvselev( sobs:eobs,3,islot)/1000.0d0, &
           tvslon(  sobs:eobs,3,islot), &
           tvslat(  sobs:eobs,3,islot), &
           tvsland( sobs:eobs,3,islot), &
           bt(    :,sobs:eobs,3,islot), &
           tran(:,:,sobs:eobs,3,islot) )
      end if
    end if
  end if

  if( ADM_prc_me > 3*prc_per_inst .and. ADM_prc_me <= 4*prc_per_inst ) then
    islot=1
    if(ntvsprof(4,islot) /= 0 ) then
      obsint=ntvsprof(4,islot)/prc_per_inst
      if(ntvsprof(4,islot)<prc_per_inst) obsint=1
      sobs=obsint*(ADM_prc_me-3*prc_per_inst-1)+1
      eobs=obsint*(ADM_prc_me-3*prc_per_inst)
      if(ADM_prc_me==4*prc_per_inst) eobs=ntvsprof(4,islot)
      WRITE(ADM_LOG_FID,*) 'NOAA-19', ntvsprof(4,islot), ADM_prc_me, sobs, eobs
      FLUSH(ADM_LOG_FID)
      if(ntvsprof(4,islot)>=prc_per_inst .or. &
         ntvsprof(4,islot)< prc_per_inst .and. ADM_prc_me-3*prc_per_inst <= ntvsprof(4,islot)) then
      call AMSUA_fwd(nlev_RTTOV, eobs-sobs+1, rttovcoef_fname(4),&
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 4, id_pres_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 4, id_temp_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 4, id_qvap_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 4, id_tsfc_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 4, id_qv2m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 4, id_surp_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 4, id_u10m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 4, id_v10m_nicam, islot), &
           soza(    sobs:eobs,4,islot), &
           soaz(    sobs:eobs,4,islot), &
           saza(    sobs:eobs,4,islot), &
           saaz(    sobs:eobs,4,islot), &
           tvselev( sobs:eobs,4,islot)/1000.0d0, &
           tvslon(  sobs:eobs,4,islot), &
           tvslat(  sobs:eobs,4,islot), &
           tvsland( sobs:eobs,4,islot), &
           bt(    :,sobs:eobs,4,islot), &
           tran(:,:,sobs:eobs,4,islot) )
      end if
    end if
  end if

  if( ADM_prc_me > 4*prc_per_inst .and. ADM_prc_me <= 5*prc_per_inst ) then
    islot=1
    if(ntvsprof(5,islot) /= 0 ) then
      obsint=ntvsprof(5,islot)/prc_per_inst
      if(ntvsprof(5,islot)<prc_per_inst) obsint=1
      sobs=obsint*(ADM_prc_me-4*prc_per_inst-1)+1
      eobs=obsint*(ADM_prc_me-4*prc_per_inst)
      if(ADM_prc_me==5*prc_per_inst) eobs=ntvsprof(5,islot)
      WRITE(ADM_LOG_FID,*) 'METOP-2', ntvsprof(5,islot), ADM_prc_me, sobs, eobs
      FLUSH(ADM_LOG_FID)
      if(ntvsprof(5,islot)>=prc_per_inst .or. &
         ntvsprof(5,islot)< prc_per_inst .and. ADM_prc_me-4*prc_per_inst <= ntvsprof(5,islot)) then
      call AMSUA_fwd(nlev_RTTOV, eobs-sobs+1, rttovcoef_fname(5),&
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 5, id_pres_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 5, id_temp_nicam, islot), &
           obsdata_jprb(nlev_RTTOV:1:-1, sobs:eobs, 5, id_qvap_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 5, id_tsfc_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 5, id_qv2m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 5, id_surp_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 5, id_u10m_nicam, islot), &
           obsdata_jprb(        1, sobs:eobs, 5, id_v10m_nicam, islot), &
           soza(    sobs:eobs,5,islot), &
           soaz(    sobs:eobs,5,islot), &
           saza(    sobs:eobs,5,islot), &
           saaz(    sobs:eobs,5,islot), &
           tvselev( sobs:eobs,5,islot)/1000.0d0, &
           tvslon(  sobs:eobs,5,islot), &
           tvslat(  sobs:eobs,5,islot), &
           tvsland( sobs:eobs,5,islot), &
           bt(    :,sobs:eobs,5,islot), &
           tran(:,:,sobs:eobs,5,islot) )
      end if
    end if


  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(9)=MPI_WTIME()
  sum_time(4)=sum_time(4)+(tmp_time(9)-tmp_time(8))

  call MPI_ALLREDUCE(bt, bt_tmp, 15*maxtvsprof*ninstrument*nslots,&
                     MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(tran, tran_tmp, nlev_RTTOV*15*maxtvsprof*ninstrument*nslots,&
                     MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  tmp_time(10)=MPI_WTIME()
  sum_time(5)=sum_time(5)+(tmp_time(10)-tmp_time(9))

  !if( ADM_prc_me == 1 ) then
  if( ADM_prc_me <= ninstrument .and. ntvs(ADM_prc_me) /= 0 ) then
    nn = ADM_prc_me
    maxtvsprof=maxval(ntvs(:))
    allocate( weight(nlev_RTTOV,15,maxtvsprof,ninstrument,nslots) )
    allocate( weight_maxlev(15,maxtvsprof,ninstrument,nslots) )
  
    tmp_time(31)=MPI_WTIME()
    iobs=1
    do islot = 1, nslots
      !do nn = 1, ninstrument
        do n =  1, ntvsprof(nn,islot)
          p_tmp(:)=log(real(obsdata_jprb(:,n,nn,id_pres_nicam,islot),kind=8))
          do ichan = tvsch(1,nn), tvsch(ntvsch(nn),nn)
            weight(1,ichan,n,nn,islot)= &
                   -(tran_tmp(2,ichan,n,nn,islot)-tran_tmp(1,ichan,n,nn,islot))/&
                    (p_tmp(nlev_RTTOV-1)-p_tmp(nlev_RTTOV))
            do ilev = 2, nlev_RTTOV-1
               weight(ilev,ichan,n,nn,islot)= &
                   -(tran_tmp(ilev+1,ichan,n,nn,islot)-tran_tmp(ilev-1,ichan,n,nn,islot))/&
                    (p_tmp(nlev_RTTOV-ilev)-p_tmp(nlev_RTTOV-ilev+2))
            end do
            weight(nlev_RTTOV,ichan,n,nn,islot)= &
                   -(tran_tmp(nlev_RTTOV,ichan,n,nn,islot)-tran_tmp(nlev_RTTOV-1,ichan,n,nn,islot))/&
                    (p_tmp(1)-p_tmp(2))
            tmp_lev=-1.0
            do ilev = 1, nlev_RTTOV
               if( weight(ilev,ichan,n,nn,islot) > tmp_lev ) then
                  weight_maxlev(ichan,n,nn,islot)=nlev_RTTOV-ilev+1
                  tmp_lev=weight(ilev,ichan,n,nn,islot)
               end if
            end do
          end do
          iobs=iobs+1
        !end do
      end do
    end do
    tmp_time(32)=MPI_WTIME()
    sum_time(31)=sum_time(31)+(tmp_time(32)-tmp_time(31))
  
    !nn=1
    do ichan = tvsch(1,nn), tvsch(ntvsch(nn),nn)
      write(*,*) 'channel', ichan
      do ilev = 1, nlev_RTTOV
        write(*,'(i5,3f10.3,i5)') ilev, tran_tmp(ilev,ichan,1,nn,1), weight(ilev,ichan,1,nn,1), &
              real(obsdata_jprb((nlev_RTTOV-ilev+1),1,nn,id_pres_nicam,1),kind=8), weight_maxlev(ichan,1,nn,1)
      end do
    end do
    tmp_time(33)=MPI_WTIME()
    sum_time(32)=sum_time(32)+(tmp_time(33)-tmp_time(32))
    !
    ! BIAS CORRECTION
    !
    maxtvsprof=maxval(ntvsprof(:,:))
    ALLOCATE(vbc_pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots))
    !
    iobs1=0
    do islot = 1, nslots
      !do nn = 1, ninstrument
        do n =  1, ntvsprof(nn,islot)
          iobs1=iobs1+1
          vbc_pred(1,:,n,nn,islot)=undef ! IWLR is depending on ch (calc. on part2)
          vbc_pred(2,:,n,nn,islot)=undef ! IWLR is depending on ch (calc. on part2)
          vbc_pred(3,:,n,nn,islot)=0.d0
          vbc_pred(4,:,n,nn,islot)=(obsdata_out(1,iobs1,id_tsfc_nicam)-273.15d0)/10.d0
          vbc_pred(5,:,n,nn,islot)=0.d0
          vbc_pred(6,:,n,nn,islot)=obsdata_out(1,iobs1,id_cldw_nicam)/30.0d0
          vbc_pred(7,:,n,nn,islot)=1.d0/cos(saza(n,nn,islot)*CNST_PI/180.d0)
          vbc_pred(8,:,n,nn,islot)=0.d0
          do ic = 1, ntvsch(nn)
            iwlr=0.0d0
            do ilev = 1, nlev-1
              if(real(obsdata_jprb(ilev,n,nn,id_pres_nicam,islot),kind=8)>200.0d0 .and. &
                 real(obsdata_jprb(ilev,n,nn,id_pres_nicam,islot),kind=8)<850.0d0) then
                iwlr = iwlr &
                     +(real(obsdata_jprb(ilev+1,n,nn,id_temp_nicam,islot),kind=8) &
                      -real(obsdata_jprb(ilev  ,n,nn,id_temp_nicam,islot),kind=8))*&
                       (tran_tmp(nlev-ilev  , tvsch(ic,nn),n,nn,islot)&
                       -tran_tmp(nlev-ilev+1, tvsch(ic,nn),n,nn,islot))
              end if
            end do
            vbc_pred(1,ic,n,nn,islot)=iwlr
          end do
          do ic = 1, ntvsch(nn)
            iwlr=0.0d0
            do ilev = 1, nlev-1
              if(real(obsdata_jprb(ilev,n,nn,id_pres_nicam,islot),kind=8)>50.0d0.and. &
                 real(obsdata_jprb(ilev,n,nn,id_pres_nicam,islot),kind=8)<200.0d0) then
                iwlr = iwlr &
                     +(real(obsdata_jprb(ilev+1,n,nn,id_temp_nicam,islot),kind=8)&
                      -real(obsdata_jprb(ilev,n,nn,id_temp_nicam,islot),kind=8))*&
                       (tran_tmp(nlev-ilev  , tvsch(ic,nn),n,nn,islot)&
                       -tran_tmp(nlev-ilev+1, tvsch(ic,nn),n,nn,islot))
              end if
            end do
            vbc_pred(2,ic,n,nn,islot)=iwlr
          end do
        end do
      !end do
    end do
  
    tmp_time(34)=MPI_WTIME()
    sum_time(33)=sum_time(33)+(tmp_time(34)-tmp_time(33))
    CALL vbc_read('vbcf_coef.txt',0,vbcf)
    CALL vbc_scan_read('vbcf_scanbias_coef.txt',0,vbcf_scan)

    
!    write(ADM_LOG_FID,*) '### READING SCAN BIAS'
!    do ic = 1, ntvsch(nn)
!      write(ADM_LOG_FID,'(2i5,30F6.2)') nn, ic, (vbcf_scan(ifoot,ic,nn),ifoot=1,nfootp(nn))
!    end do 
!
!    write(ADM_LOG_FID,*) '### BEFORE normalize', sum(vbcf_scan(1:nfootp(nn),1,nn))/real(nfootp(nn)-6)
!    do ic=1,ntvsch(nn)
!      vbcf_scan(1:nfootp(nn),ic,nn) = vbcf_scan(1:nfootp(nn),ic,nn) - &
!                                  sum(vbcf_scan(1:nfootp(nn),ic,nn)) / real(nfootp(nn)-6)
!    end do
!    write(ADM_LOG_FID,*) '###  AFTER normalize', sum(vbcf_scan(1:nfootp(nn),1,nn))/real(nfootp(nn)-6)
!  
    tmp_time(35)=MPI_WTIME()
    sum_time(34)=sum_time(34)+(tmp_time(35)-tmp_time(34))

    do islot = 1, nslots
      do n =  1, ntvsprof(nn,islot)
        ifoot=tvsfoot(n,nn,islot)
        do ic=1,ntvsch(nn)
          ichan=tvsch(ic,nn)
          !tvsdat(ic,n,nn,islot)=tvsdat(ic,n,nn,islot)-vbcf_scan(ifoot,ic,nn)
          tvsdat(ichan,n,nn,islot)=tvsdat(ichan,n,nn,islot)-vbcf_scan(ifoot,ic,nn)
        end do
      end do
    end do

    !
    ! Do not correct the airmass bias
    !
    !!vbc_pred(:,:,:,:,:)=0.0d0

    do islot = 1, nslots
      do n =  1, ntvsprof(nn,islot)
        do ic=1,ntvsch(nn)
          ichan=tvsch(ic,nn)
          !tvsdat(ic,n,nn,islot)=tvsdat(ic,n,nn,islot)+&
          tvsdat(ichan,n,nn,islot)=tvsdat(ichan,n,nn,islot)+&
                   sum(vbc_pred(:,ic,n,nn,islot)*vbcf(:,ic,nn))
          write(ADM_LOG_FID,'(3i5,9F8.3)') islot, n, ic, &
                 (vbc_pred(i,ic,n,nn,islot)*vbcf(i,ic,nn),i=1,8), &
                 bt_tmp(tvsch(ic,nn),n,nn,islot)
        end do
      end do
    end do
  
    tmp_time(36)=MPI_WTIME()
    sum_time(35)=sum_time(35)+(tmp_time(36)-tmp_time(35))

    tvsqc(:,:,:,:)=1

    do islot = 1, nslots
      do n =  1, ntvsprof(nn,islot)
        ! [QUALITY CHECK FOR AMSU-A]
        ! Add 2016/02/02 Avoid cloud over land (refer to Bormann, 2010)
        ! Clear sky       : chs <=5 are not assimilated
        ! Clear sky       : ch6 is assimilated if z>1500m
        ! Clear sky       : ch7 is assimilated if z>2500m
        ! Clear sky       : chs >= 8 are assimilated
        ! Cloudy or rainy : no assimilation (chs >= 9 are assimilated in JMA)
        if( lsql(n,nn,islot)==0 ) then
          tvsqc(:,n,nn,islot)=0
!          if( abs( tvsdat(4,n,nn,islot)-bt_tmp(4,n,nn,islot) ) < 0.7 ) then
!            do ic = 1, ntvsch(nn)
!              ! Clear sky       : chs <=5 are not assimilated
!              if(tvsch(ic,nn) <= 5) then
!                tvsqc(ic,n,nn,islot)=0
!              end if
!              ! Clear sky       : ch6 is assimilated if z>1500m
!              if(tvsch(ic,nn) == 6 .and. tvselev(n,nn,islot) > 1500) then
!                tvsqc(ic,n,nn,islot)=0
!              end if
!              ! Clear sky       : ch7 is assimilated if z>2500m
!              if(tvsch(ic,nn) == 7 .and. tvselev(n,nn,islot) > 2500) then
!                tvsqc(ic,n,nn,islot)=0
!              end if
!            end do
!          end if
!          if( abs( tvsdat(4,n,nn,islot)-bt_tmp(4,n,nn,islot) ) > 0.7 ) then
!            do ic = 1, ntvsch(nn)
!              if(tvsch(ic,nn) <= 8) then
!                tvsqc(ic,n,nn,islot)=0
!              end if
!            end do
!          end if
        end if
        ! Add 2016/02/02 Avoid cloud over ocean (refer to Bormann, 2010)
        ! Clear sky       : chs >= 4 are assimilated
        ! Cloudy sky      : chs >= 7 are assimilated
        !!! rainy           : chs >= 9 are assimilated [not ready]
        if( lsql(n,nn,islot)==1 ) then
          if( abs( tvsdat(3,n,nn,islot)-bt_tmp(3,n,nn,islot) ) > 3.0 ) then
            do ic = 1, ntvsch(nn)
              if( tvsch(ic,nn) < 7 ) then
                tvsqc(ic,n,nn,islot)=0
              end if
            end do
          end if
        end if
        do ic=1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          if(tvslat(n,nn,islot)<-60.0) tvsqc(ic,n,nn,islot)=0
          if(tvslat(n,nn,islot)> 60.0) tvsqc(ic,n,nn,islot)=0
          if(tvsdat(ic,n,nn,islot)<100.0 .or. tvsdat(ic,n,nn,islot)>400.0) then
            tvsqc(ic,n,nn,islot)=0
          end if
        end do
      end do
    end do

    i=1
    islot=1
    !do nn = 1, ninstrument
      write(cfile(1:4),'(A4)') tvsname(nn)
      !outbase =trim(outfile_dir)//'/'//trim(cfile)//trim(outfile_prefix)//trim(cimem)
      outbase =trim(outfile_dir)//'/'//trim(cfile)//trim(outfile_prefix)
      ofid = 1
      write(*,*) 'Output: ', trim(outbase)//'.dat'
      open( unit   = ofid,                  &
            file   = trim(outbase)//'.dat', &
            form   = 'unformatted',         &
            access = 'sequential',          &
            status = 'unknown'              )
      do n = 1, ntvsprof(nn,islot)
        write(ofid) real(lsql(n,nn,islot)), &
                    real(tvslon(n,nn,islot)), &
                    real(tvslat(n,nn,islot)), &
                    real(saza(n,nn,islot)), &
                    real(obsdata_out(1,i,id_tsfc_nicam)), &
                    real(0.0), &
                    real(tvsfoot(n,nn,islot)), &
                    real(obsdata_out(weight_maxlev(tvsch(1:ntvsch(nn),nn),n,nn,islot),i,id_pres_nicam)),&
                    real(tvsdat(tvsch(1,nn):tvsch(ntvsch(nn),nn),n,nn,islot)),   &
                    real(tvserr(1:ntvsch(nn),n,nn,islot)),   &
                    real(bt_tmp(tvsch(1:ntvsch(nn),nn),n,nn,islot)), &
                    real( tvsqc(1:ntvsch(nn),n,nn,islot))
        i=i+1
      enddo
      close(ofid)
    !enddo

    close(1)
    tmp_time(37)=MPI_WTIME()
    sum_time(36)=sum_time(36)+(tmp_time(37)-tmp_time(36))
!
    !if(ocheck) then
      i=1
      !do nn = 1, ninstrument
        write(cfile(1:4),'(A4)') tvsname(nn)
        !outbase2 = trim(outfile_dir)//'/'//trim(cfile)//trim(outfile_prefix)//trim(cimem)
        outbase2 = trim(outfile_dir)//'/'//trim(cfile)//trim(outfile_prefix)
        ofid2 = 2
        write(ADM_LOG_FID,*) 'nn=', nn
        write(ADM_LOG_FID,*) 'ntvsch(nn)=', ntvsch(nn)
        write(ADM_LOG_FID,*) 'tvsch=', tvsch(1:ntvsch(nn),nn)
        write(*,*) 'Output: ', trim(outbase2)//'.txt'
  !
        open( unit   = ofid2,                   &
              file   = trim(outbase2)//'.txt', &
              form   = 'formatted',            &
              status = 'unknown'              )
        write(ofid2,*) 'nn=', nn
        write(ofid2,*) 'ntvsch(nn)=', ntvsch(nn)
        write(ofid2,*) 'tvsch=', tvsch(1:ntvsch(nn),nn)
        do n = 1, ntvsprof(nn,islot)
          write(ofid2,'(30F8.2)') &
                    real(lsql(n,nn,islot)), &
                    real(tvslon(n,nn,islot)), &
                    real(tvslat(n,nn,islot)), &
                    real(saza(n,nn,islot)), &
                    real(obsdata_out(1,i,id_tsfc_nicam)), &
                    real(0.0), &
                    real(obsdata_out(weight_maxlev(tvsch(1:ntvsch(nn),nn),n,nn,islot),i,id_pres_nicam)),&
                    !real(tvsdat(1:ntvsch(nn),n,nn,islot)),   &
                    real(tvsdat(tvsch(1,nn):tvsch(ntvsch(nn),nn),n,nn,islot)),   &
                    real(bt_tmp(tvsch(1:ntvsch(nn),nn),n,nn,islot)), &
                    real(tvsfoot(n,nn,islot))
          i=i+1
        end do
      !end do
      close(ofid2)
    !end if ! ocheck

    tmp_time(38)=MPI_WTIME()
    sum_time(37)=sum_time(37)+(tmp_time(38)-tmp_time(37))
    deallocate( weight )
    DEALLOCATE( vbc_pred )
    deallocate( weight_maxlev )
  end if

  deallocate( bt )
  deallocate( bt_tmp) 
  deallocate( tran) 
  deallocate( tran_tmp) 
  deallocate( hx )
  !deallocate( dtau )
  tmp_time(11)=MPI_WTIME()
  sum_time(6)=sum_time(6)+(tmp_time(11)-tmp_time(10))
  tmp_time(12)=MPI_WTIME()

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize:        ', tmp_time(2)-tmp_time(1)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Read obs:          ', tmp_time(3)-tmp_time(2)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Calc coef:         ', tmp_time(4)-tmp_time(3)
  !write(ADM_LOG_FID,'(A,F15.5)') '##### Prepare read gues: ', tmp_time(5)-tmp_time(4)
  !write(ADM_LOG_FID,'(A,F15.5)') '##### Communication:     ', tmp_time(6)-tmp_time(5)
  !write(ADM_LOG_FID,'(A,F15.5)') '##### rttov computation: ', tmp_time(7)-tmp_time(6)
  !write(ADM_LOG_FID,'(A,F15.5)') '##### Communication:     ', tmp_time(8)-tmp_time(7)
  !write(ADM_LOG_FID,'(A,F15.5)') '##### VBC and output:    ', tmp_time(9)-tmp_time(8)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Prepare read gues: ', sum_time(1)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Read Gues:         ', sum_time(2)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Communication:     ', sum_time(3)
  write(ADM_LOG_FID,'(A,F15.5)') '##### rttov computation: ', sum_time(4)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Communication:     ', sum_time(5)
  write(ADM_LOG_FID,'(A,F15.5)') '##### VBC and output:    ', sum_time(6)
  write(ADM_LOG_FID,'(A,F15.5)') '##### TOTAL ELAPSED TIME:', tmp_time(12)-tmp_time(1)

  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize01       ', tmp_time(22)-tmp_time(21)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize02       ', tmp_time(23)-tmp_time(22)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize03       ', tmp_time(24)-tmp_time(23)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize04       ', tmp_time(25)-tmp_time(24)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize05       ', tmp_time(26)-tmp_time(25)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize06       ', tmp_time(27)-tmp_time(26)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize07       ', tmp_time(28)-tmp_time(27)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize08       ', tmp_time(29)-tmp_time(28)
  write(ADM_LOG_FID,'(A,F15.5)') '##### Initialize09       ', tmp_time(30)-tmp_time(29)

  if(ADM_prc_me==1) then
    write(ADM_LOG_FID,'(A,F15.5)') '##### Compute weight:   ', sum_time(31)
    write(ADM_LOG_FID,'(A,F15.5)') '##### Write text:       ', sum_time(32)
    write(ADM_LOG_FID,'(A,F15.5)') '##### Compute predictor:', sum_time(33)
    write(ADM_LOG_FID,'(A,F15.5)') '##### Read vbc:         ', sum_time(34)
    write(ADM_LOG_FID,'(A,F15.5)') '##### Bias correction:  ', sum_time(35)
    write(ADM_LOG_FID,'(A,F15.5)') '##### OUTPUT binary:    ', sum_time(36)
    write(ADM_LOG_FID,'(A,F15.5)') '##### OUTPUT text:      ', sum_time(37)
  end if

  call ADM_proc_stop

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
    !!!character(12) :: name
    !!!real(8) :: zz, zz0
    real(8) :: pp, pp0
    real(8) :: tt, tv
    real(8) :: qq
    !!!real(8) :: t1, t2
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
    !!!integer :: max_num_lon
    !!!integer :: max_num_lat_lon2(ADM_lall)
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
!!  tvsch( 2,3)  =  6
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
!  READ(MD,'(X,A8,I4,X,I2,X,I2,X,I2,X,I2)') &
!      & platname,iyy,imm,idd,ihr,imn
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
  RETURN
END SUBROUTINE vbc_read
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

end program prg_obsope_amsua_mean


