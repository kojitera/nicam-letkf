  program main
    use mod_adm, only :         & 
         ADM_proc_init,         &
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
         ADM_LOG_FID
    use mod_misc
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
    integer :: l, rgnid, n, i, g, k
    integer, parameter :: nbv = 128
    character(ADM_MAXFNAME) :: basename
    character(ADM_MAXFNAME) :: fname
    character(ADM_MAXFNAME) :: mean_basename
    character(ADM_MAXFNAME) :: mean_fname
    character(ADM_MAXFNAME) :: member_basename
    character(ADM_MAXFNAME) :: member_fname
    character(ADM_MAXFNAME) :: output_basename
    character(ADM_MAXFNAME) :: restart_layername
    character(3) :: ens_num
    character(6) :: cmem
    real(8), allocatable :: ddd1(:,:,:,:,:)
    real(8), allocatable :: ddd2(:,:,:,:,:)
    real(8), allocatable :: mean(:,:,:,:)
    real(8), allocatable :: sprd(:,:,:,:)
    real(8), allocatable :: sprd_ave(:,:)
    real(8), allocatable :: ddd3_pl(:,:,:,:,:)

    INTEGER,PARAMETER :: nv3d=8 ! u,v,t,q
    character(LEN=FIO_HSHORT) :: varname3d(nv3d)
    integer, allocatable :: ifid(:)

    namelist / RESTARTPARAM /  &
         mean_basename,        &
         member_basename,      &
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

    call FIO_setup ! [add] H.Yashiro 20110826

    write(*,*) ADM_gall, ADM_kall, ADM_lall

    allocate( ddd1(ADM_gall, ADM_kall, ADM_lall, 3, nbv+1) )
    allocate( ddd2(ADM_gall, ADM_kall, ADM_lall, 2, nbv+1) )
    allocate( sprd(ADM_gall, ADM_kall, ADM_lall, 2) )
    allocate( mean(ADM_gall, ADM_kall, ADM_lall, 2) )
    allocate( sprd_ave(ADM_kall, 2) )
    allocate( ddd3_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl, 5, nbv+1) )

    open(1,file='nhm_driver.cnf')
    read(1,nml=RESTARTPARAM) 
    close(1)

    restart_layername='ZSALL40'

    call FIO_input(ddd1(:,:,:,1,nbv+1),mean_basename,'vx', restart_layername,1,ADM_kall,1)
    call FIO_input(ddd1(:,:,:,2,nbv+1),mean_basename,'vy', restart_layername,1,ADM_kall,1)
    call FIO_input(ddd1(:,:,:,3,nbv+1),mean_basename,'vz', restart_layername,1,ADM_kall,1)

    do i = 1, nbv
      write(cmem,'(I6.6)') i
      member_fname=trim(member_basename)//'/'//trim(cmem)//'/restart'
      call FIO_input(ddd1(:,:,:,1,i),member_fname,'vx',restart_layername,1,ADM_kall,1)
      call FIO_input(ddd1(:,:,:,2,i),member_fname,'vy',restart_layername,1,ADM_kall,1)
      call FIO_input(ddd1(:,:,:,3,i),member_fname,'vz',restart_layername,1,ADM_kall,1)
    end do

    do i = 1, nbv+1
      call GTL_generate_uv(                     &
         ddd2(:,:,:,1,i), ddd3_pl(:,:,:,1,i),   &
         ddd2(:,:,:,2,i), ddd3_pl(:,:,:,2,i),   &
         ddd1(:,:,:,1,i), ddd3_pl(:,:,:,3,i),   &
         ddd1(:,:,:,2,i), ddd3_pl(:,:,:,4,i),   &
         ddd1(:,:,:,3,i), ddd3_pl(:,:,:,5,i), icos=0)
    end do

    write(ADM_LOG_FID,*) nbv+1, minval(ddd2(:,2:ADM_kall-1,:,1,nbv+1)), maxval(ddd2(:,2:ADM_kall-1,:,1,nbv+1))
    write(ADM_LOG_FID,*) nbv+1, minval(ddd2(:,2:ADM_kall-1,:,2,nbv+1)), maxval(ddd2(:,2:ADM_kall-1,:,2,nbv+1))
    do i = 1, nbv
      ddd2(:,:,:,1,i) = ( ddd2(:,:,:,1,i) - ddd2(:,:,:,1,nbv+1) )**2
      ddd2(:,:,:,2,i) = ( ddd2(:,:,:,2,i) - ddd2(:,:,:,2,nbv+1) )**2
      !ddd2(:,:,:,1,i) = ddd2(:,:,:,1,i) - ddd2(:,:,:,1,nbv+1)
      !ddd2(:,:,:,2,i) = ddd2(:,:,:,2,i) - ddd2(:,:,:,2,nbv+1)
      !write(ADM_LOG_FID,*) i, minval(ddd2(:,2:ADM_kall-1,:,1,i)), maxval(ddd2(:,2:ADM_kall-1,:,1,i))
      !write(ADM_LOG_FID,*) i, minval(ddd2(:,2:ADM_kall-1,:,2,i)), maxval(ddd2(:,2:ADM_kall-1,:,2,i))
    end do

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
      sprd(g,k,l,1)=sqrt( sum(ddd2(g,k,l,1,1:nbv)) / REAL(nbv-1) )
      sprd(g,k,l,2)=sqrt( sum(ddd2(g,k,l,2,1:nbv)) / REAL(nbv-1) )
    end do
    end do
    end do

    do k = 1, ADM_kall
      sprd_ave(k,1)=sum(sprd(:,k,:,1))/REAL(ADM_gall*ADM_lall)
      sprd_ave(k,2)=sum(sprd(:,k,:,2))/REAL(ADM_gall*ADM_lall)
    end do

    !do k = 1, ADM_kall
    !  write(ADM_LOG_FID,'(i4,3f12.5)') k, minval(sprd(:,k,:,1)), maxval(sprd(:,k,:,1)), sprd_ave(k,1)
    !  write(ADM_LOG_FID,'(i4,3f12.5)') k, minval(sprd(:,k,:,2)), maxval(sprd(:,k,:,2)), sprd_ave(k,2)
    !end do

    do l = 1, ADM_lall
      rgnid=ADM_prc_tab(l,ADM_prc_me)
      write(ADM_LOG_FID,*) l, rgnid
      call MISC_make_idstr(fname,Trim(output_basename),'rgn',rgnid)
      write(ADM_LOG_FID,*) trim(fname)
      open(1,file=trim(fname),form='unformatted',access='direct',&
           recl=4*ADM_gall*ADM_kall)
      do n = 1,2
        write(1,rec=n) sngl(sprd(:,:,l,n))
      end do
      close(1)
    end do


    call ADM_proc_finish

  end program main
