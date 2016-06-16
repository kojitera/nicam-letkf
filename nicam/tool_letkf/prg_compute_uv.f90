  program main
    use mod_adm, only :         & 
         ADM_proc_init,         &
         ADM_proc_stop,         &
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
         ADM_LOG_FID
    use mod_misc
    use mod_fio, only : & ! [add] H.Yashiro 20110819
         FIO_input
    use mod_comm, only :           &
         COMM_setup
    use mod_latlon, only :        &
         LATLON_setup,            &
         LATLON_read_outdirname
    use mod_cnst, only :           &
         CNST_setup
    use mod_grd, only :            &
         GRD_setup
    use mod_gtl, only :            &
         GTL_generate_uv
    use mod_gmtr, only :           &
         GMTR_setup

    implicit none
    integer :: nmem
    integer :: imem
    integer :: l, rgnid, n
    character(ADM_MAXFNAME) :: basename
    character(ADM_MAXFNAME) :: fname
    character(ADM_MAXFNAME) :: input_basename
    character(ADM_MAXFNAME) :: output_basename
    character(ADM_MAXFNAME) :: restart_layername
    character(3) :: ens_num
    real(8), allocatable :: ddd1(:,:,:,:)
    real(8), allocatable :: ddd2(:,:,:,:)
    real(8), allocatable :: ddd3_pl(:,:,:,:)

    namelist / RESTARTPARAM /  &
         input_basename,       &
         output_basename

    call ADM_proc_init(ADM_MULTI_PRC)
    call ADM_setup('nhm_driver.cnf')
    !call ADM_setup('compute_uv.cnf')

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
    write(*,*) ADM_gall, ADM_kall, ADM_lall
    allocate( ddd1(ADM_gall, ADM_kall, ADM_lall, 8) )
    allocate( ddd2(ADM_gall, ADM_kall, ADM_lall, 2) )
    allocate( ddd3_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl, 5) )

    open(1,file='nhm_driver.cnf')
    read(1,nml=RESTARTPARAM) 
    close(1)
    write(ADM_LOG_FID,*) trim(input_basename)
    write(ADM_LOG_FID,*) trim(output_basename)
    do l = 1, ADM_lall
      rgnid=ADM_prc_tab(l,ADM_prc_me)
      write(ADM_LOG_FID,*) l, rgnid
      call MISC_make_idstr(fname,Trim(input_basename),'rgn',rgnid)
      write(ADM_LOG_FID,*) trim(fname)
      open(1,file=trim(fname),form='unformatted',access='direct',&
           recl=8*ADM_gall*ADM_kall)
      do n = 1,8
        read(1,rec=n) ddd1(:,:,l,n) 
        write(ADM_LOG_FID,*) minval(ddd1(:,2:ADM_kall-1,l,n)), maxval(ddd1(:,2:ADM_kall-1,l,n))
      end do
      close(1)
    end do

    call GTL_generate_uv(                 &
         ddd2(:,:,:,1), ddd3_pl(:,:,:,1), &
         ddd2(:,:,:,2), ddd3_pl(:,:,:,2), &
         ddd1(:,:,:,3), ddd3_pl(:,:,:,3), &
         ddd1(:,:,:,4), ddd3_pl(:,:,:,4), &
         ddd1(:,:,:,5), ddd3_pl(:,:,:,5), 0)

    write(ADM_LOG_FID,*) 'CHECK COMPUTED U AND V FROM VX, VY, AND VZ'
    write(ADM_LOG_FID,*) minval(ddd2(:,2:ADM_kall-1,:,1)), maxval(ddd2(:,2:ADM_kall-1,:,1))
    write(ADM_LOG_FID,*) minval(ddd2(:,2:ADM_kall-1,:,2)), maxval(ddd2(:,2:ADM_kall-1,:,2))
  
    do l = 1, ADM_lall
      rgnid=ADM_prc_tab(l,ADM_prc_me)
      write(ADM_LOG_FID,*) l, rgnid
      call MISC_make_idstr(fname,Trim(output_basename),'rgn',rgnid)
      write(ADM_LOG_FID,*) trim(fname)
      open(1,file=trim(fname),form='unformatted',access='direct',&
           !recl=8*ADM_gall*(ADM_kall-2))
           recl=8*ADM_gall*ADM_kall)
      do n = 1,2
        write(1,rec=n) ddd1(:,:,l,n)
        !write(1,rec=n) ddd1(:,2:ADM_kall-1,l,n)
      end do
      write(1,rec=3) ddd2(:,:,l,1)
      write(1,rec=4) ddd2(:,:,l,2)
      !write(1,rec=3) ddd2(:,2:ADM_kall-1,l,1)
      !write(1,rec=4) ddd2(:,2:ADM_kall-1,l,2)
      do n = 5,7
        write(1,rec=n) ddd1(:,:,l,n)
        !write(1,rec=n) ddd1(:,2:ADM_kall-1,l,n)
      end do
      close(1)
    end do


    call ADM_proc_stop

  end program main
