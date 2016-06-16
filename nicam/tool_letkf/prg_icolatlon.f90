program icolatlon
  use mod_adm
  use mpi
  use mod_misc, only :       &
       MISC_get_latlon,      &
       MISC_make_idstr
  use mod_cnst, only :       &
       CNST_setup
  use mod_comm, only :       &
       COMM_setup
  use mod_fio, only : &  ! [add] H.Yashiro 20110826
       FIO_setup
  use mod_grd, only :        &
       GRD_setup,            &
       GRD_x
  use mod_gmtr, only :       &
       GMTR_setup
  implicit none
 
  REAL(8),ALLOCATABLE,SAVE :: ico_lon(:,:)
  REAL(8),ALLOCATABLE,SAVE :: ico_lat(:,:)
  REAL(8),ALLOCATABLE,SAVE :: ico_lon_tmp(:,:)
  REAL(8),ALLOCATABLE,SAVE :: ico_lat_tmp(:,:)
 
  REAL(8),PARAMETER :: pi=3.1415926535d0

  integer :: g, l, ierr
  integer :: nij0, ngpv

  character(128) :: icolatlon_fname=''

  namelist / icolatlon_cnf / &
    icolatlon_fname

  CALL ADM_proc_init(ADM_MULTI_PRC)

  call ADM_setup('nhm_driver.cnf')
  !--- [add] 09/08/18 .Mitsui 
  ! Extra processes should be used for any other purposes in future.
  if( .not. ADM_myprc_is_run ) then
     call ADM_proc_stop
  end if
  !--- < I/O module setup > ---
  !------ setup file i/o.
  call FIO_setup ! [add] H.Yashiro 20110826
  !--- < comm module setup > ---
  !------ setup commnication.
  call COMM_setup
  !--- < cnst module setup > ---
  !------ setup the constants.
  call CNST_setup
  !--- < grid module setup > ---
  !------ setup the icosahedral grid and vertical grid.
  call GRD_setup
  !--- < gmetrics module setup > ---
  !------ setup the horizontal metrics.
  call GMTR_setup

  ALLOCATE( ico_lon_tmp(ADM_gall, ADM_lall) )
  ALLOCATE( ico_lat_tmp(ADM_gall, ADM_lall) )
  ALLOCATE( ico_lon(ADM_gall, ADM_rgn_nmax) )
  ALLOCATE( ico_lat(ADM_gall, ADM_rgn_nmax) )

  write(ADM_LOG_FID,*) 'ADM_prc_run_all ', ADM_prc_run_all
  write(ADM_LOG_FID,*) 'ADM_myprc_is_run', ADM_myprc_is_run

  nij0=ADM_gall*ADM_rgn_nmax
  ngpv=nij0*ADM_kall

  do l = 1, ADM_lall
    do g = 1, ADM_gall
      call MISC_get_latlon(ico_lat_tmp(g,l), ico_lon_tmp(g,l), GRD_x(g, ADM_KNONE, l, 1), &
               GRD_x(g, ADM_KNONE, l, 2), GRD_x(g, ADM_KNONE, l, 3))
      write(ADM_LOG_FID,'(2i6,2f15.10)') l, g, ico_lat_tmp(g,l)/pi*180.0, &
                                               ico_lon_tmp(g,l)/pi*180.0
    end do
  end do
  ico_lat_tmp(:,:)=ico_lat_tmp(:,:)/pi*180.0
  !ico_lon_tmp(:,:)=ico_lon_tmp(:,:)/pi*180.0+180.0d0
  ico_lon_tmp(:,:)=ico_lon_tmp(:,:)/pi*180.0  ! Fixed [2015.01.21] Koji

  call MPI_GATHER(ico_lon_tmp, ADM_gall*ADM_lall, MPI_REAL8, &
                  ico_lon,     ADM_gall*ADM_lall, MPI_REAL8, &
                  0, MPI_COMM_WORLD, ierr)
  call MPI_GATHER(ico_lat_tmp, ADM_gall*ADM_lall, MPI_REAL8, &
                  ico_lat,     ADM_gall*ADM_lall, MPI_REAL8, &
                  0, MPI_COMM_WORLD, ierr)

  if(ADM_prc_me==1) then
    open(1,file='icolatlon.cnf')
    read(1,nml=icolatlon_cnf)
    close(1)
    open(1,file=trim(icolatlon_fname),form='unformatted',access='sequential')
    write(1) ico_lon
    write(1) ico_lat
    close(1)
    write(*,*) 'minimum value', minval(ico_lon)
    write(*,*) 'maximum value', maxval(ico_lon)
    write(*,*) 'minimum value', minval(ico_lat)
    write(*,*) 'maximum value', maxval(ico_lat)
  end if

  call MPI_Finalize(ierr)
    
end program icolatlon
