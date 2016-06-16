!-------------------------------------------------------------------------------
!
!+  Program obsope
!
!-------------------------------------------------------------------------------
program prg_obsope

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

  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
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
  character(LEN=FIO_HLONG)  :: infile_header(flim) = ''
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
  real(8), allocatable, save :: lat_obs(:)     ! latitude
  real(8), allocatable, save :: lon_obs(:)     ! longitude
  real(8), allocatable, save :: elev(:)        ! pressure [hPa]
  real(8), allocatable, save :: typ(:)         ! observation type
  !integer, allocatable, save :: elev(:)        ! elevation [m]
  integer, allocatable, save :: agl(:)         ! above ground level = 1, if not = 0
  integer, allocatable, save :: land(:)        ! land = 1, ocean = 0
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
  real(8), allocatable, save :: obserr(:)
  real(8), allocatable, save :: odat(:)
  integer, allocatable, save :: elem(:)
  integer, allocatable, save :: id_obs(:)
  integer, allocatable, save :: oqc(:)
  integer, allocatable, save :: oqc_out(:)

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
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_z_obs=2567
  INTEGER,PARAMETER :: id_s_obs=3332
  INTEGER,PARAMETER :: id_rain_obs=9999

  INTEGER,PARAMETER :: id_u=1
  INTEGER,PARAMETER :: id_v=2
  INTEGER,PARAMETER :: id_t=3
  INTEGER,PARAMETER :: id_q=4
  INTEGER,PARAMETER :: id_ps=5
  INTEGER,PARAMETER :: id_rain=6

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

  real(4) :: wk(7)
  real(8), allocatable :: pres(:,:,:)
  real(8), allocatable :: lnpres(:,:,:)
  real(8), allocatable :: tem(:,:,:)
  real(8), allocatable :: qv(:,:,:)
  real(8), allocatable :: ps(:,:)
  real(8), allocatable :: lnps(:,:)
  real(8), allocatable :: ts(:,:)
  real(8), allocatable :: qs(:,:)
  real(8), allocatable :: ms_u(:,:,:)
  real(8), allocatable :: ms_v(:,:,:)
  real(8), allocatable :: ms_t(:,:,:)
  real(8), allocatable :: ms_qv(:,:,:)
  real(4), allocatable :: obsdata(:)
  real(4), allocatable :: obsdata_out(:)

  logical :: ocheck=.false.
  real(8) :: tmp_time(20)
  real(8) :: time(10)
  real(8) :: time1(10,10)

  integer :: smem, emem, imem
  character(6) :: csmem, cemem, cimem

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
       pres_basename,      &
       ocheck

  namelist /OPTION/ glevel,            &
                    rlevel,            &
                    smem,              &
                    emem,              &
                    grid_topology,     &
                    complete,          &
                    mnginfo,           &
                    layerfile_dir,     &
                    llmap_base,        &
                    infile,     &
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

  time(:)=0.0d0
  tmp_time(1)=MPI_WTIME()

  call ADM_proc_init(ADM_MULTI_PRC)
  call ADM_setup('obsope.cnf')
  call COMM_setup
  call CNST_setup
  call GRD_setup
  call GMTR_setup
  call OPRT_setup
  call VMTR_setup
  call readoption !! set fmax, infile

  write(csmem(1:6),'(I6.6)') smem
  write(cemem(1:6),'(I6.6)') emem
  write(ADM_LOG_FID,*) 'csmem', csmem
  write(ADM_LOG_FID,*) 'cemem', cemem

  rewind(ADM_CTL_FID)
  read(ADM_CTL_FID, nml=obsope_param, iostat=ierr)

  fid=40
  open(fid, file=trim(input_fname), form='unformatted', access='sequential', &
       status='old', iostat=ierr)
  if(ierr /= 0) then
     write(ADM_LOG_FID,*) 'ERROR for opening the file', trim(input_fname)
     write(ADM_LOG_FID,*) 'STOP!'
     call ADM_proc_stop
  end if

  tmp_time(2)=MPI_WTIME()
  time(1)=tmp_time(2)-tmp_time(1)

  call get_nsite(fid,nsite) ! number of observation

  tmp_time(3)=MPI_WTIME()
  time(2)=tmp_time(3)-tmp_time(2)

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
  allocate(qs(ADM_gall,ADM_lall))
  allocate(ms_u(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ms_v(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ms_t(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(ms_qv(ADM_gall,ADM_vlayer,ADM_lall))
  allocate(icodata(ADM_gall,ADM_vlayer,ADM_lall,10))
  allocate(obsdata(nsite))
  allocate(obsdata_out(nsite))

  allocate(obsnum(nsite))
  allocate(lon(nsite))
  allocate(lat(nsite))
  allocate(lon_obs(nsite))
  allocate(lat_obs(nsite))
  allocate(elev(nsite))
  allocate(typ(nsite))
  allocate(agl(nsite))
  allocate(land(nsite))
  allocate(cdirec(nsite))
  allocate(idirec(2,nsite))
  allocate(date(6,nsite))
  allocate(var(nsite))
  allocate(obserr(nsite))
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
  allocate(odat(nsite))
  allocate(oqc(nsite))
  allocate(oqc_out(nsite))
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
  klev(:,:) = -999
  kfact(:) = 0.0d0
  land_model(:) = -1.0
  idirec_model(:) = -1

  tmp_time(4)=MPI_WTIME()
  time(3)=tmp_time(4)-tmp_time(3)

  call calendar_yh2ss( time_tmp, idate )
  time_obs(:)=dnint(time_tmp)

  rewind(fid)
  do s=1, nsite
    read(fid) wk
    call get_variable(wk(1),wk(4),var(s))
    lon_obs(s)=wk(2)
    lat_obs(s)=wk(3)
    if(wk(2) <= 180.0)then
       lon(s) = dble(wk(2))
    else
       lon(s) = dble(wk(2)) - 360.d0
    endif
    lat(s) = dble(wk(3))
    elev(s) = dble(wk(4))
    odat(s) = dble(wk(5))
    obserr(s) = dble(wk(6))
    typ(s) = dble(wk(7))
    elem(s) = wk(1)
    select case(int(wk(1)))
      case(id_u_obs)
        id_obs(s)=id_u
      case(id_v_obs)
        id_obs(s)=id_v
      case(id_t_obs)
        id_obs(s)=id_t
      case(id_q_obs)
        id_obs(s)=id_q
      case(id_ps_obs)
        id_obs(s)=id_ps
      case(id_rain_obs)
        id_obs(s)=id_rain
    end select
  enddo

  do s=1, nsite
    if( lat(s) == -90.0 ) then
      lat(s)=-89.99
    end if
    if( lat(s) ==  90.0 ) then
      lat(s)= 89.99
    end if
    if( lat(s) > 90.0 .or. lat(s) < -90.0 ) then
      write(ADM_LOG_FID,*) s, lat(s), lon(s)
    end if
    if( lon(s) > 180.0 .or. lat(s) < -180.0 ) then
      write(ADM_LOG_FID,*) s, lat(s), lon(s)
    end if
  enddo

  lon(:) = lon(:)*CNST_PI/180.0d0
  lat(:) = lat(:)*CNST_PI/180.0d0

  tmp_time(5)=MPI_WTIME()
  time(4)=tmp_time(5)-tmp_time(4)

  write(ADM_LOG_FID,*) 'ADM_GALL_PL=', ADM_GALL_PL
  write(ADM_LOG_FID,*) 'ADM_LALL_PL=', ADM_LALL_PL
  FLUSH(ADM_LOG_FID)
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
     call MPI_ALLREDUCE(max_num_latlon, sum_max_num_latlon, 1, MPI_INTEGER, MPI_SUM, &
            MPI_COMM_WORLD, ierr)
     !
     if(sum_max_num_latlon /= nsite) then
        !write(ADM_LOG_FID,*) 'ERROR!: max_num_latlon /= nsite'
        !write(ADM_LOG_FID,*) 'max_num_latlon:', sum_max_num_latlon
        !write(ADM_LOG_FID,*) 'STOP!'
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
  tmp_time(6)=MPI_WTIME()
  time(5)=tmp_time(6)-tmp_time(5)

  do imem = smem, emem

    tmp_time(7)=MPI_WTIME()
    write(cimem(1:6),'(I6.6)') imem
    infile_header(1)=trim(infile(1))//'/'//trim(cimem)//'/history'
    do p = pstr, pend
       pp = p - pstr + 1

       if (complete) then ! all region
        infname = trim(infile_header(1))//'.rgnall'
       else
          call fio_mk_fname(infname,trim(infile_header(1)),'pe',p-1,6)
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

           call fio_fclose(ifid(pp)) ! [add] 13-04-18

        enddo PE_loop

           if( trim(var_name(v)) == 'ms_pres' ) then
             pres(:,:,:)=icodata4(:,:,:)*0.01d0
             lnpres(:,:,:)=log(pres(:,:,:))
           else if( trim(var_name(v)) == 'ss_ps' ) then
             ps(:,:)=icodata4(:,1,:)*0.01d0
             lnps(:,:)=log(ps(:,:))
             icodata(:,1,:,5)=lnps(:,:)
           else if( trim(var_name(v)) == 'ss_q2m' ) then
             qs(:,:)=icodata4(:,1,:)
           else if( trim(var_name(v)) == 'ss_t2m' ) then
             ts(:,:)=icodata4(:,1,:)
           else if( trim(var_name(v)) == 'ms_u' ) then
             icodata(:,:,:,1)=icodata4(:,:,:)
           else if( trim(var_name(v)) == 'ms_v' ) then
             icodata(:,:,:,2)=icodata4(:,:,:)
           else if( trim(var_name(v)) == 'ms_tem' ) then
             icodata(:,:,:,3)=icodata4(:,:,:)
           else if( trim(var_name(v)) == 'ms_qv' ) then
             icodata(:,:,:,4)=icodata4(:,:,:)
           end if

        deallocate( data4allrgn )
        deallocate( data8allrgn )
        deallocate( icodata4    )

     enddo step_loop ! step LOOP

  enddo variable_loop ! variable LOOP

  tmp_time(8)=MPI_WTIME()
  time(6)=time(6)+(tmp_time(8)-tmp_time(7))

  call getklev

  tmp_time(9)=MPI_WTIME()
  time(7)=time(7)+(tmp_time(9)-tmp_time(8))

  oqc(:)=0
  obsdata(:)=0.0
  k = 1
  do p = pstr, pend
    do l = 1, ADM_lall
      do i = 1,nsite
        if( l_index(i) == MNG_prc_tab(l,p)-pstr+1 .and. &
            time_obs(i) == nowsec )then
          ks = klev(1,i)
          ke = klev(2,i)
          if( ks==-999 ) then
            obsdata(i)=CNST_UNDEF4
            oqc(i)=0
          else if( icodata(n1_index(i),ks,l,id_obs(i))  == CNST_UNDEF4 &
             .OR.  icodata(n2_index(i),ks,l,id_obs(i))  == CNST_UNDEF4 &
             .OR.  icodata(n3_index(i),ks,l,id_obs(i))  == CNST_UNDEF4 &
             .OR.  kfact(i)==-999 ) then
            obsdata(i)=CNST_UNDEF4
            oqc(i)=0
          else 
            if(id_obs(i)==5) then !! Surface pressure
              fac1 = w1(i)
              fac2 = w2(i)
              fac3 = w3(i)
              fac_sum = fac1 + fac2 + fac3
              obsdata(i) = ( fac1 * icodata(n1_index(i),ks,l,id_obs(i)) &
                           + fac2 * icodata(n2_index(i),ks,l,id_obs(i)) &
                           + fac3 * icodata(n3_index(i),ks,l,id_obs(i)) &
                           ) / fac_sum
              !obsdata(i) = exp(obsdata(i)*kfact(i))
              obsdata(i) = exp(obsdata(i))*kfact(i)
            else
              fac1 = w1(i) * kfact(i)
              fac2 = w2(i) * kfact(i)
              fac3 = w3(i) * kfact(i)
              fac4 = w1(i) * (1.d0-kfact(i))
              fac5 = w2(i) * (1.d0-kfact(i))
              fac6 = w3(i) * (1.d0-kfact(i))
              fac_sum = fac1 + fac2 + fac3 + fac4 + fac5 + fac6
              obsdata(i) = ( fac1 * icodata(n1_index(i),ks,l,id_obs(i)) &
                           + fac2 * icodata(n2_index(i),ks,l,id_obs(i)) &
                           + fac3 * icodata(n3_index(i),ks,l,id_obs(i)) &
                           + fac4 * icodata(n1_index(i),ke,l,id_obs(i)) &
                           + fac5 * icodata(n2_index(i),ke,l,id_obs(i)) &
                           + fac6 * icodata(n3_index(i),ke,l,id_obs(i)) &
                           ) / fac_sum
            end if
            oqc(i)=1
            !exit
          end if
        end if
      end do
    end do
  end do

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(10)=MPI_WTIME()
  time(8)=time(8)+(tmp_time(10)-tmp_time(9))

  call MPI_ALLREDUCE(obsdata, obsdata_out, nsite, MPI_REAL, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)

  call MPI_ALLREDUCE(oqc, oqc_out, nsite, MPI_INTEGER, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)

  tmp_time(11)=MPI_WTIME()
  time(9)=time(9)+(tmp_time(11)-tmp_time(10))

  where( obsdata_out(:) == CNST_UNDEF4 ) obsdata_out(:)=-999.9

!  do i = 1, nsite
!    write(ADM_LOG_FID,'(i5,4f12.5)') i, lon(i), lat(i), obsdata_out(i), obsdata(i)
!    write(20001,'(i5,4f12.5)') i, lon(i), lat(i), obsdata_out(i), real(oqc_out(i))
!  end do

  if(ADM_prc_me==1) then
    outbase = trim(outfile_dir)//'/'//trim(outfile_prefix)//trim(cimem)
    ofid = 1
    !write(*,*) 'Output: ', trim(outbase)//'.dat'
  
    open( unit   = ofid,                  &
          file   = trim(outbase)//'.dat', &
          form   = 'unformatted',         &
          access = 'sequential',          &
          status = 'unknown'              )
  
    do i = 1,nsite
       write(ofid) real(elem(i)),        real(lon_obs(i)), &
                   real(lat_obs(i)),     real(elev(i)),    &
                   real(odat(i)),        real(obserr(i)),  &
                   real(obsdata_out(i)), real(oqc_out(i)), &
                   real(typ(i))
    enddo
  
    close(1)

    if(ocheck) then
      outbase2 = trim(outfile_dir)//'/'//trim(outfile_prefix)//trim(cimem)
      ofid2 = 2
      !write(*,*) 'Output: ', trim(outbase2)//'.txt'
    
      open( unit   = ofid2,                   &
            file   = trim(outbase2)//'.txt', &
            form   = 'formatted',            &
            status = 'unknown'              )

      do i = 1,nsite
         write(ofid2,'(9F10.3)') &
               real(elem(i)),        real(lon_obs(i)), &
               real(lat_obs(i)),     real(elev(i)),    &
               real(odat(i)),        real(obserr(i)),  &
               real(obsdata_out(i)), real(oqc_out(i)), &
               real(typ(i))
      end do

      close(ofid2) 

    end if

  end if

  end do
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(12)=MPI_WTIME()
  time(10)=time(10)+(tmp_time(12)-tmp_time(11))

  call MPI_Allgather(time(1),    10,  MPI_DOUBLE, &
                     time1(1,1), 10,  MPI_DOUBLE, MPI_COMM_WORLD, ierr)

  write(ADM_LOG_FID,'(A,3F15.5)') '##### Initialize:     ', time(1),&
                      maxval(time1(1,:)), sum(time1(1,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Get nobs:       ', time(2),&
                      maxval(time1(2,:)), sum(time1(2,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Allocation:     ', time(3),&
                      maxval(time1(3,:)), sum(time1(3,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Read obs:       ', time(4),&
                      maxval(time1(4,:)), sum(time1(4,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Compute x-coef: ', time(5),&
                      maxval(time1(5,:)), sum(time1(5,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Read gues:      ', time(6),&
                      maxval(time1(6,:)), sum(time1(6,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Compute z-coef: ', time(7),&
                      maxval(time1(7,:)), sum(time1(7,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### intepolation:   ', time(8),&
                      maxval(time1(8,:)), sum(time1(8,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Allreduce:      ', time(9),&
                      maxval(time1(9,:)), sum(time1(9,:))/10.0d0
  write(ADM_LOG_FID,'(A,3F15.5)') '##### Output:         ', time(10),&
                      maxval(time1(10,:)), sum(time1(10,:))/10.0d0

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
             !write(ADM_LOG_FID,*) k, tmp_elev, pp0
             if( k == 1 ) then
                klev(1,s)=1
                klev(2,s)=1
                kfact(s)=1.0d0
                if( elem(s) == id_t_obs ) then
                   kfact(s) = -999.0d0
                end if
             else if ( k >= ADM_vlayer ) then
                kfact(s) = -999.0d0
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

end program prg_obsope

