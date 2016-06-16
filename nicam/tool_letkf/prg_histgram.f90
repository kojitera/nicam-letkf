  program main
    use mpi
    use mod_adm, only :         & 
         ADM_proc_init,         &
         ADM_proc_stop,         &
         ADM_proc_finish,       &
         ADM_setup,             &
         ADM_MULTI_PRC,         &
         ADM_gall,              &
         ADM_kall,              &
         ADM_lall,              &
         ADM_gall_pl,           &
         ADM_lall_pl,           &
         ADM_prc_tab,           &
         ADM_prc_me,            &
         ADM_MAXFNAME,          &
         ADM_LOG_FID,           &
         ADM_prc_tab
    use mod_misc
    use mod_fio, only : &  ! [add] H.Yashiro 20110826
      FIO_setup,        &
      FIO_input,        &
      FIO_output,       &
      FIO_finalize,     &
      FIO_HMID,         &
      FIO_HSHORT,       &
      FIO_HLONG,        &
      FIO_INTEG_FILE,   &
      FIO_SPLIT_FILE,   &
      FIO_FREAD,        &
      FIO_ICOSAHEDRON,  &
      FIO_BIG_ENDIAN,   &
      headerinfo,       &
      datainfo,         &
      FIO_REAL8
    use mod_comm, only :           &
         COMM_setup
    use mod_latlon, only :        &
         LATLON_setup,            &
         LATLON_read_outdirname
    use mod_cnst, only :           &
         CNST_setup,               &
         CNST_UNDEF
    use mod_grd, only :            &
         GRD_setup,                &
         GRD_gz,                   &
         GRD_gzh,                  &
         GRD_zs,                   &
         GRD_zs_pl,                &
         GRD_vz,                   &
         GRD_vz_pl
    use mod_gtl, only :            &
         GTL_generate_uv
    use mod_gmtr, only :           &
         GMTR_setup

    implicit none
    integer :: nmem
    integer :: imem
    integer :: g, l, rgnid, n, nv, ip
    integer :: ierr
    integer :: i, j
    character(ADM_MAXFNAME) :: basename
    character(ADM_MAXFNAME) :: fname
    character(ADM_MAXFNAME) :: restart_dirbase
    character(ADM_MAXFNAME) :: restart_basename
    character(ADM_MAXFNAME) :: restart_fname
    character(ADM_MAXFNAME) :: history_basename
    character(ADM_MAXFNAME) :: output_basename
    character(ADM_MAXFNAME) :: restart_layername
    character(ADM_MAXFNAME) :: surface_layername
    character(ADM_MAXFNAME) :: pressure_layername
    character(3) :: ens_num
    real(8), allocatable :: ddd_logp(:,:,:)
    real(8), allocatable :: ddd_3d_org(:,:,:)
    real(8), allocatable :: ddd_3d(:,:,:,:)
    real(8), allocatable :: ddd_2d(:,:,:,:)
    real(8), allocatable :: ddd_uv(:,:,:,:)
    real(8), allocatable :: ddd3_pl(:,:,:,:)
    real(8), allocatable :: output(:)

    !--- 
    real(8), allocatable :: ddd_3d_plev(:,:,:,:)

    real(8) :: plevel(200)
    real(8) :: logp(200)

    real(8), allocatable :: cc(:)
    real(8), allocatable :: ww(:)
    real(8), allocatable :: bx(:)
    real(8), allocatable :: by(:)

    INTEGER :: np
    INTEGER,PARAMETER :: nv3d=8
    INTEGER,PARAMETER :: nv2d=29
    character(LEN=FIO_HSHORT) :: varname3d(nv3d)
    character(LEN=FIO_HSHORT) :: varname2d(nv2d)
    integer, allocatable :: ifid(:)

    !--- ico2ll
    integer :: imax, jmax
    integer,allocatable :: lon_index(:),lat_index(:)
    integer,allocatable :: n1_index(:), n2_index(:), n3_index(:)
    real(8),allocatable :: w1(:), w2(:), w3(:)

    integer,allocatable :: max_num_latlon(:)
    integer,allocatable :: nstart(:), nend(:)
    character(1024) :: fin_llmap_head = '../llmap/llmap'
    character(1024) :: fin_llmap

    real(8),allocatable :: lon(:)  ! -pi - pi
    real(8),allocatable :: lat(:)  ! -pi/2 - pi/2
    real(8),allocatable :: lon_swap(:)

    real(8),allocatable :: var_out_3d(:,:,:,:)
    real(8),allocatable :: var_out_2d(:,:,:,:)
    real(8),allocatable :: var_3d(:,:,:,:)
    real(8),allocatable :: var_2d(:,:,:,:)
    real(8),allocatable :: var_out_swap(:,:)
    character(5) :: c5
    character(6) :: c6
    integer :: imem
    integer :: nmem
 
    logical :: onefile=.false.
    character(256) :: output_dir
    character(10)  :: cdate

    character(LEN=FIO_HMID) :: desc = 'INITIAL/RESTART DATA of PROGNOSTICVARIABLES'

    namelist / RESTARTPARAM /    &
         restart_dirbase,        &
         restart_basename,       &
         output_basename,        &
         varname3d,              &
         varname2d,              &
         plevel,                 &
         np,                     &
         fin_llmap_head,         &
         output_dir,             &
         cdate, onefile, nmem

    call ADM_proc_init(ADM_MULTI_PRC)
    call ADM_setup('nhm_driver.cnf')

    !--- < comm module setup > ---
    !------ setup commnication.
    call COMM_setup
    !
    !--- < cnst module setup > ---
    !------ setup the constants.
    call CNST_setup
    !
    !--- < grid module setup > ---
    !------ setup the icosahedral grid and vertical grid.
    call GRD_setup
    !
    !--- < gmetrics module setup > ---
    !------ setup the horizontal metrics.
    call GMTR_setup
    !allocate( ddd_3d_org(ADM_gall, ADM_kall, ADM_lall, 3) )
    allocate( ddd_3d_org(ADM_gall, ADM_kall, ADM_lall) )
    !allocate( ddd_3d(ADM_gall, ADM_kall, ADM_lall, nv3d) )
    allocate( ddd_logp(ADM_gall, ADM_kall, ADM_lall) )
    allocate( ddd_uv(ADM_gall, ADM_kall, ADM_lall,    2) )
    allocate( ddd_2d(ADM_gall,        1, ADM_lall, nv2d) )
    allocate( ddd3_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl, 5) )

    open(1,file='nhm_driver.cnf')
    read(1,nml=RESTARTPARAM) 
    close(1)
    write(ADM_LOG_FID,*) 'restart_basename,', trim(restart_basename)
    write(ADM_LOG_FID,*) 'nmem,', nmem
    logp(1:np)=log(plevel(1:np))
    allocate( ddd_3d_plev(ADM_gall, np, ADM_lall, nv3d) )
    allocate( ddd_3d    (ADM_gall, ADM_kall, ADM_lall, nmem) )

    restart_layername='ZSALL40'
    surface_layername='ZSSFC1'
    pressure_layername='ZSDEF26'
    !
    !--- < Reading data > ---
    call FIO_setup ! [add] H.Yashiro 20110826

    do imem = 1, nmem
      write(c6,'(I6.6)') imem
      restart_fname=trim(restart_dirbase)//'/'//trim(c6)//'/'//trim(restart_basename)
      FLUSH(ADM_LOG_FID)
      call FIO_input(ddd_3d_org(:,:,:),restart_fname,'tem',restart_layername,1,ADM_kall,1)
      ddd_3d(:,:,:,imem)=ddd_3d_org(:,:,:)
      write(ADM_LOG_FID,*) imem, maxval(ddd_3d_org(:,2,:)), minval(ddd_3d_org(:,2,:))
    end do

    write(ADM_LOG_FID,*) minloc(ddd_3d(:,2,:,91))
    !if(ADM_prc_me==38) then
      do imem = 1, nmem
        write(ADM_LOG_FID,*) imem, ddd_3d(597,2,1,imem)
      end do
    !end if

    call MPI_BARRIER(mpi_comm_world,ierr)
    write(ADM_LOG_FID,*) 'ierr=',ierr
    flush(ADM_LOG_FID)

    call ADM_proc_finish

  end program main
