program main
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
       ADM_proc_finish,        &
       ADM_prc_me,             &
       ADM_prc_run_master,     &
       ADM_prc_run_all,        &
       ADM_prc_tab,            &
       ADM_kall,               &
       ADM_gall,               &
       ADM_KNONE,              &
       ADM_lall,               &
       ADM_GALL_PL,            &
       ADM_LALL_PL,            &
       ADM_gall_in,            &
       ADM_vlayer,             &
       ADM_glevel,             &
       ADM_rlevel,             &
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
    FIO_setup,        &
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
    datainfo,         &
    fio_input,        &
    fio_output,       &
    fio_hmid
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

  integer                   :: glevel              = -1
  integer                   :: rlevel              = -1
  character(14) :: cdate14
  character(6) :: cmem
  character(ADM_MAXFNAME), save :: input_basename = ''
  character(ADM_MAXFNAME), save :: output_basename = ''
  character(ADM_MAXFNAME), save :: input_fname = ''
  character(ADM_MAXFNAME), save :: output_fname = ''
  character(LEN=FIO_HMID) :: desc = 'INITIAL/RESTART DATA of LAND VARIABLES'
  integer :: nmem
  integer :: ierr, iv, imem

  !integer, parameter :: vmax = 11
  integer, parameter :: vmax = 30
  integer, parameter :: kgmax = 5
  integer, parameter :: kwmax = 5
  integer, parameter :: ltsfc = 2      ! max number of surface skin temp.
  integer, parameter :: lwsfc = 2      ! max number of surface water 
  integer, parameter :: ksnmax = 3     ! Max level of snow temperature ! add
  integer, parameter :: nrbnd = 3 
  integer, parameter :: nrdir = 2
  integer :: kmax(vmax)
  integer :: ksta(vmax)
  integer :: kend(vmax)
  integer :: kall

  real(8) :: ctime          ! time
  real(8), allocatable :: tmp(:,:,:,:)
  real(8), allocatable :: mean(:,:,:,:)
  integer, allocatable :: num(:,:,:,:)

  integer :: g, k, l
  integer :: i

  namelist / makemean_param / &
    nmem,                     &
    input_basename,           &
    output_basename,          &
    cdate14

  call ADM_proc_init(ADM_MULTI_PRC)
  call ADM_setup('nhm_driver.cnf')
  call COMM_setup
  call FIO_setup
  call GRD_setup
  call GMTR_setup
  call OPRT_setup
  call VMTR_setup

  glevel = ADM_glevel
  rlevel = ADM_rlevel

  rewind(ADM_CTL_FID)
  read(ADM_CTL_FID, nml=makemean_param, iostat=ierr)

  write(ADM_LOG_FID,*) 'input_basename=', trim(input_basename)
  write(ADM_LOG_FID,*) 'output_basename=', trim(output_basename)

  kmax(1:11)=(/kgmax, kwmax, ltsfc, ltsfc, 1, 1, ksnmax, kwmax, nrbnd, 1, nrbnd*nrdir/)

  kall=0
  !do iv=1,vmax
  do iv=1,11
     ksta(iv)= kall + 1
     kend(iv)= kall + kmax(iv)
     kall    = kall + kmax(iv)
  enddo

  do iv = 12, 17
    ksta(iv)= kall + 1
    kend(iv)= kall + 1
    kall    = kall + 1
  end do

  do iv = 18, 22
    ksta(iv)= kall + 1
    kend(iv)= kall + ADM_kall
    kall    = kall + ADM_kall
  end do

  ksta(23) = kall + 1
  kend(23) = kall + 1
  kall     = kall + 1

  allocate(tmp(ADM_gall,kall+5,ADM_lall,1))
  allocate(mean(ADM_gall,kall+5,ADM_lall,1))
  allocate(num(ADM_gall,kall+5,ADM_lall,1))

  num(:,:,:,:)=0
  mean(:,:,:,:)=0.0d0
  !cdate14='20140601000000'

  do imem = 1, nmem

    write(cmem,'(i6.6)') imem
    write(ADM_LOG_FID,*) cmem
    input_fname=trim(input_basename)//'/'//cmem//'/restart_other'//cdate14
    write(ADM_LOG_FID,*) 'input_fname=', trim(input_fname)

    call FIO_input( tmp(:,ksta(1):kend(1),:,1), &
                    input_fname,'glg','GLKGMAX',1,kgmax,1 )
    call FIO_input( tmp(:,ksta(2):kend(2),:,1), &
                    input_fname,'glw','GLKWMAX',1,kwmax,1 )
    call FIO_input( tmp(:,ksta(3):kend(3),:,1), &
                    input_fname,'glts','GLLTSFC',1,ltsfc,1 )
    call FIO_input( tmp(:,ksta(4):kend(4),:,1), &
                    input_fname,'gltc','GLLTSFC',1,ltsfc,1 )
    call FIO_input( tmp(:,ksta(5):kend(5),:,1), &
                    input_fname,'glwc','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(6):kend(6),:,1), &
                    input_fname,'glsnw','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(7):kend(7),:,1), &
                    input_fname,'gltsn','GLKSNMAX',1,ksnmax,1 )
    call FIO_input( tmp(:,ksta(8):kend(8),:,1), &
                    input_fname,'glfrs','GLKWMAX',1,kwmax,1 )
    call FIO_input( tmp(:,ksta(9):kend(9),:,1), &
                    input_fname,'glasn','GLALB',1,nrbnd,1 )
    call FIO_input( tmp(:,ksta(10):kend(10),:,1), &
                    input_fname,'snrtco','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(11):kend(11),:,1), &
                    input_fname,'albsfc','GRALB',1,nrbnd*nrdir,1)

    call FIO_input( tmp(:,ksta(12):kend(12),:,1), &
                    input_fname,'tem_sfc','ZSSFC1',1,1,1 )

    call FIO_input( tmp(:,ksta(13):kend(13),:,1), &
                    input_fname,'gosst','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(14):kend(14),:,1), &
                    input_fname,'goice','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(15):kend(15),:,1), &
                    input_fname,'goicr','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(16):kend(16),:,1), &
                    input_fname,'gosnw','ZSSFC1',1,1,1 )
    call FIO_input( tmp(:,ksta(17):kend(17),:,1), &
                    input_fname,'goist','ZSSFC1',1,1,1 )

    call FIO_input( tmp(:,ksta(18):kend(18),:,1), &
                    input_fname,'CBMFX','ZSALL40',1,ADM_kall,1 )
    call FIO_input( tmp(:,ksta(19):kend(19),:,1), &
                    input_fname,'qked','ZSALL40',1,ADM_kall,1 )
    call FIO_input( tmp(:,ksta(20):kend(20),:,1), &
                    input_fname,'tsqd','ZSALL40',1,ADM_kall,1 )
    call FIO_input( tmp(:,ksta(21):kend(21),:,1), &
                    input_fname,'qsqd','ZSALL40',1,ADM_kall,1 )
    call FIO_input( tmp(:,ksta(22):kend(22),:,1), &
                    input_fname,'covd','ZSALL40',1,ADM_kall,1 )

    call FIO_input( tmp(:,ksta(23):kend(23),:,1), &
                    input_fname,'ROUGHNESS_SEA','ZSSFC1',1,1,1 )


    do l = 1, ADM_lall
    do k = 1, kall
    do g = 1, ADM_gall
      if(tmp(g,k,l,1).ne.CNST_UNDEF) then
        mean(g,k,l,1) = mean(g,k,l,1) + tmp(g,k,l,1)
        num(g,k,l,1)=num(g,k,l,1)+1
      end if
    end do
    end do
    end do

  end do

  do l = 1, ADM_lall
  do k = 1, kall
  do g = 1, ADM_gall
    if(num(g,k,l,1).eq.0) then
      mean(g,k,l,1) = CNST_UNDEF
    else
      mean(g,k,l,1) = mean(g,k,l,1)/dble(num(g,k,l,1))
    end if
  end do
  end do
  end do

  do iv = 1, 22
  do k = ksta(iv), kend(iv) 
    write(ADM_LOG_FID,*) iv, k, minval(mean(:,k,:,1)), maxval(mean(:,k,:,1))
  end do
  end do

  output_fname = trim(output_basename)//'/mean_land'

  ctime=0

  i=1
  call FIO_output( mean(:,ksta(1):kend(1),:,1),output_fname, desc, '',  &
                  'glg', 'Soil Temperature', '', 'K',                &
                   FIO_REAL8, 'GLKGMAX', 1, kgmax, 1, ctime, ctime   )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(2):kend(2),:,1),output_fname, desc, '',  &
                  'glw', 'Soil Moisture', '', 'm3/m3',               &
                   FIO_REAL8, 'GLKWMAX', 1, kwmax, 1, ctime, ctime   )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(3):kend(3),:,1),output_fname, desc, '',  &
                  'glts', 'Land Skin Temperature', '', 'K',          &
                   FIO_REAL8, 'GLLTSFC', 1, ltsfc, 1, ctime, ctime   )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(4):kend(4),:,1),output_fname, desc, '',  &
                  'gltc', 'Canopy Temperature', '', 'K',             &
                   FIO_REAL8, 'GLLTSFC', 1, ltsfc, 1, ctime, ctime   )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(5):kend(5),:,1),output_fname, desc, '',  &
                  'glwc', 'Canopy Water', '', 'm3/m2',               &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime        )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(6):kend(6),:,1),output_fname, desc, '',  &
                  'glsnw', 'Land Snow Amount', '', 'kg/m2',          &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime        )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(7):kend(7),:,1),output_fname, desc, '',  &
                  'gltsn', 'Land Snow Temperature', '', 'K',         &
                   FIO_REAL8, 'GLKSNMAX', 1, ksnmax, 1, ctime, ctime )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(8):kend(8),:,1),output_fname, desc, '',  &
                  'glfrs', 'Soil Ice', '', 'm3/m3',                  &
                   FIO_REAL8, 'GLKWMAX', 1, kwmax, 1, ctime, ctime   )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(9):kend(9),:,1),output_fname, desc, '',  &
                  'glasn', 'Land Snow Albedo', '', '0-1',            &
                   FIO_REAL8, 'GLALB', 1, nrbnd, 1, ctime, ctime     )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(10):kend(10),:,1),output_fname, desc, '',  &
                  'snrtco', 'Canopy Snow Fraction', '', '0-1',       &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime        )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(11):kend(11),:,1),output_fname, desc, '',  &
                  'albsfc', 'Albedo', '', '0-1',                       &
                   FIO_REAL8, 'GRALB', 1, nrbnd*nrdir, 1, ctime, ctime )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1

  call FIO_output( mean(:,ksta(12):kend(12),:,1),output_fname, desc, '',  &
                  'tem_sfc', 'Surface Temperature', '', 'K', &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1

  call FIO_output( mean(:,ksta(13):kend(13),:,1),output_fname, desc, '',  &
                  'gosst', 'Sea Surface Temperature', '', 'K', &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(14):kend(14),:,1),output_fname, desc, '',  &
                  'goice', 'Sea Ice Mass', '', 'kg/m2',        &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(15):kend(15),:,1),output_fname, desc, '',  &
                  'goicr', 'Sea Ice Fraction', '', '0-1',      &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(16):kend(16),:,1),output_fname, desc, '',  &
                  'gosnw', 'Sea Surface Snow', '', 'kg/m2',    &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(17):kend(17),:,1),output_fname, desc, '',  &
                  'goist', 'Sea Ice Temperature', '', 'K',     &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1

  call FIO_output( mean(:,ksta(18):kend(18),:,1),output_fname, desc, '',  &
                  'CBMFX', 'Cloud-base Mass Flux', '', 'kg/kg', &
                   FIO_REAL8, 'ZSALL40', 1, ADM_kall, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(19):kend(19),:,1),output_fname, desc, '',  &
                  'qked', 'Turbulent Kinetic Energy*2', '', 'J/kg', &
                   FIO_REAL8, 'ZSALL40', 1, ADM_kall, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(20):kend(20),:,1),output_fname, desc, '',  &
                  'tsqd', 'Variance of theta_l', '', 'K', &
                   FIO_REAL8, 'ZSALL40', 1, ADM_kall, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(21):kend(21),:,1),output_fname, desc, '',  &
                  'qsqd', 'Variance of Total Water(qw)', '', '', &
                   FIO_REAL8, 'ZSALL40', 1, ADM_kall, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1
  call FIO_output( mean(:,ksta(22):kend(22),:,1),output_fname, desc, '',  &
                  'covd', 'Covarriance of qw and theta_l', '', 'K', &
                   FIO_REAL8, 'ZSALL40', 1, ADM_kall, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1

  call FIO_output( mean(:,ksta(23):kend(23),:,1),output_fname, desc, '',  &
                   'ROUGHNESS_SEA', 'Sea Roughness Length', '', 'm', &
                   FIO_REAL8, 'ZSSFC1', 1, 1, 1, ctime, ctime  )
  WRITE(ADM_LOG_FID,*) i
  FLUSH(ADM_LOG_FID)
  i=i+1

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call ADM_proc_finish

end program main
