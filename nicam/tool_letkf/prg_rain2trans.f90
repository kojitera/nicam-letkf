!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  prg_rain2trans
!
!  program for precipitation transformation
!
!  created                Apr 2015, Shunji Kotsuki, RIKEN-AICS
!  adopted to NICAM-LETKF May 2015, Shunji Kotsuki, RIKEN-AICS 
!
!-------------------------------------------------------------------------------
program rain2trans


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
  
  ! added by kotsuki       
  use fnc_transgauss     
       
  implicit none
  REAL(8),PARAMETER :: UNDEF8 = -99.9E+33

  integer rain_mem
  
  integer g, l, ibv, ix, iy
  integer gall, lall, xgrd, ygrd
  
  character(255) fname,obsname,gusname,outname,mapname,oname,gname,mskname

  
  ! region independent
  REAL(4), ALLOCATABLE :: obs_rain(:), gus_mean(:),gus_trns_mean(:), prc_mem(:)
  REAL(4), ALLOCATABLE :: gus_rain(:,:)
  REAL(8), ALLOCATABLE :: ico_lat(:,:), ico_lon(:,:)
  REAL(8), ALLOCATABLE :: ppcdf_m(:,:), ppcdf_o(:,:)
  REAL(8), ALLOCATABLE :: ppzero_m(:),  ppzero_o(:)
  REAL(8), ALLOCATABLE :: ppmask_m(:),  ppmask_o(:)
  REAL(8), ALLOCATABLE :: ppsamp_m(:),  ppsamp_o(:)
  REAL(8), ALLOCATABLE :: dfzero_m(:),  dfzero_o(:)
  REAL(8), ALLOCATABLE :: corr_m(:),    corr_o(:)
  REAL(4), ALLOCATABLE :: usage(:)
  REAL(4), PARAMETER :: typ=99.0
  REAL(r_size), ALLOCATABLE :: gus_rain_sub(:,:)

  ! to be gathered
  REAL(4), ALLOCATABLE :: qc(:,:),         qc_tmp(:,:),         qc_rgn(:)     
  REAL(4), ALLOCATABLE :: qc_mem(:),       qc_god(:)     
  REAL(4), ALLOCATABLE :: gus_trns(:,:,:), gus_trns_tmp(:,:,:), gus_trns_rgn(:,:)
  REAL(4), ALLOCATABLE :: obs_trns(:,:),   obs_trns_tmp(:,:),   obs_trns_rgn(:)
  REAL(4), ALLOCATABLE :: obs_err(:,:),    obs_err_tmp(:,:),    obs_err_rgn(:)
  REAL(4), ALLOCATABLE :: elm(:,:),        elm_tmp(:,:),        elm_rgn(:)

  ! to the joint histogram :: 1: from all precip, 2: wo zero precip.
  integer :: icell, hobs, hgus
  integer, save        :: ncell 
  REAL(8), ALLOCATABLE :: hmin(:), hmax(:)
  REAL(4), ALLOCATABLE :: hmsk_obs(:,:),     hmsk_gus(:,:,:)
  REAL(4), ALLOCATABLE :: hmsk_obs_rgn(:),   hmsk_gus_rgn(:,:) 
  REAL(4), ALLOCATABLE :: hmsk_obs_tmp(:,:), hmsk_gus_tmp(:,:,:)
  
  !===> get ratio for treeatment of zero precip.
  REAL(r_size), ALLOCATABLE :: prate_o(:,:),      prate_m(:,:,:)                             
  REAL(4), ALLOCATABLE :: obs_rain_tmp(:,:), gus_rain_tmp(:,:,:) 
  REAL(4), ALLOCATABLE :: obs_map(:,:),      gus_map(:,:,:)
  REAL(4), ALLOCATABLE :: tmp_grd(:,:),      tmp_rgn(:,:)        
  INTEGER, ALLOCATABLE :: vec_rgn(:,:,:),    vec_grd(:,:,:)
  INTEGER, ALLOCATABLE :: cont_m(:)
  integer igrd, ngrd, gg, ll, cont_o

  !===> grid conversion
  integer, allocatable :: g2xgrd(:), g2ygrd(:)

  character(255) tbrname,tbgname

  REAL(4) :: obs_err_org
  REAL(r_size) :: obs_err_p, obs_err_n
  REAL(r_size) :: ym, sigma
  integer :: zero_mem

  !==> nhm_driver.cnf
  integer glevel
  integer rlevel
  integer vlayer
  character(128) rgnmngfname

  !--> paralell computing
  integer :: ierr, iread, nread, iprocs, nprocs, irank, myrank
  integer, allocatable :: proc_l(:,:) 

  !--> random number SFMT
  REAL(r_size) ::  genrand_res53

!========================================================================
  CALL ADM_proc_init(ADM_MULTI_PRC)  
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)  
!========================================================================
  call system_clock(seed)
  call init_gen_rand(seed)

!  ym = genrand_res53()
!  print *, ym
!  stop

  call get_rain2trans()
  call get_gaussvariable()
  call check_transform(myrank)
  call check_opt_ppobserr(myrank)

  if( opt_variloc==1 .and. myrank==0 ) then ! variable localization
    call monit_variloc( extrp_good, 'EXTRATROPIC_GOOD' )
    call monit_variloc( extrp_norm, 'EXTRATROPIC_NORM' )
    call monit_variloc( tropc_good, 'TROPICS_____GOOD' )
    call monit_variloc( tropc_norm, 'TROPICS_____NORM' )
  end if

  namelist / ADMPARAM / &
    glevel,                &
    rlevel,                &
    vlayer,                &
    rgnmngfname
  open(1,file='nhm_driver.cnf')
  read(1,nml=ADMPARAM)
  close(1)
  
  gall   =  (2**(glevel-rlevel) + 2)**2
  lall   =  10*4*rlevel
  xgrd   =  int(  dsqrt(dble(gall))  )
  ygrd   =  xgrd

  allocate ( g2xgrd(gall), g2ygrd(gall) )
  g=0
  do iy=1,ygrd
  do ix=1,xgrd
    g = g + 1
    g2xgrd(g) = ix
    g2ygrd(g) = iy
  end do
  end do
  
  
  !====> Initial Setting
  ALLOCATE ( ico_lat(gall,lall), ico_lon(gall,lall) )
  
  write(mapname,'(a,a,i2.2,a,i2.2,a)') trim(grddir),'/icolatlon_gl',glevel,'_rl',rlevel,'.dat'
  open(1,file=trim(mapname),form='unformatted',access='sequential',action='read')
  read(1) ico_lon  ! [deg]
  read(1) ico_lat  ! [deg]
  close(1)  

  !===> setting for parallel computing
  if( mod(lall,nprocs).ne.0 ) then
    write(6,*) "lall should be devided by the nprocs"
    stop
  endif
  
  nread = lall / nprocs
  allocate ( proc_l(0:nprocs-1,nread) )
  call set_proc_l(lall,nprocs,nread,proc_l)

!  if( myrank ==0 ) then
!    do iprocs=0,nprocs-1
!      write(6,'(i,a,20i)') iprocs," === ",(proc_l(iprocs,iread),iread=1,nread)  
!    end do
!  endif

! region-depending variables
  ALLOCATE ( ppcdf_m(gall,0:ndivid), ppcdf_o(gall,0:ndivid) ) 
  ALLOCATE ( ppzero_m(gall),         ppzero_o(gall)         )
  ALLOCATE ( ppmask_m(gall),         ppmask_o(gall)         )
  ALLOCATE ( ppsamp_m(gall),         ppsamp_o(gall)         )
  ALLOCATE ( dfzero_m(gall),         dfzero_o(gall)         )
  ALLOCATE ( obs_rain(gall),         gus_mean(gall)         )
  ALLOCATE ( gus_rain(gall,nbv),     gus_trns_mean(gall)    )
  ALLOCATE ( gus_rain_sub(gall,nbv)                         )
  ALLOCATE ( qc_mem(gall),           prc_mem(gall)          )
  ALLOCATE ( corr_m(gall),           corr_o(gall)           )
  ALLOCATE ( usage(gall),            qc_god(gall)           )
! variables for parallel computing
  ALLOCATE ( qc_tmp(gall,nprocs),       gus_trns_tmp(gall,nbv,nprocs) ) ! to allgather
  ALLOCATE ( obs_trns_tmp(gall,nprocs), obs_err_tmp(gall,nprocs)      ) ! to allgather
  ALLOCATE ( elm_tmp(gall,nprocs)                                     ) ! to allgather
  ALLOCATE ( qc_rgn(gall),              gus_trns_rgn(gall,nbv)        ) ! to region
  ALLOCATE ( obs_trns_rgn(gall),        obs_err_rgn(gall)             ) ! to region
  ALLOCATE ( elm_rgn(gall)                                            ) ! to region

! variables for all gathered
  ALLOCATE ( qc(gall,lall),          gus_trns(gall,lall,nbv)   )
  ALLOCATE ( obs_err(gall,lall),     obs_trns(gall,lall)       )
  ALLOCATE ( elm(gall,lall)                                    )

!  if( prdhst ) then ! used to the QC 
    ncell = int ( 6.0d0 / hstsiz ) 
    ALLOCATE ( hmin(ncell),              hmax(ncell)                  ) ! to histgram
    ALLOCATE ( hmsk_obs(gall,lall),      hmsk_gus(gall,lall,nbv)      ) ! to histgram
    ALLOCATE ( hmsk_obs_rgn(gall),       hmsk_gus_rgn(gall,nbv)       ) ! to histgram
    ALLOCATE ( hmsk_obs_tmp(gall,nprocs),hmsk_gus_tmp(gall,nbv,nprocs)) ! to histgram

    do icell = 1,ncell
      hmin(icell) = -3.0d0 + hstsiz * dble( icell - 1 )
      hmax(icell) = -3.0d0 + hstsiz * dble( icell - 0 )
    end do    
!  endif 



!!!  !===> get ratio for treeatment of zero precip.
!!!  ALLOCATE ( prate_o(gall,lall)       , prate_m(gall,lall,nbv)        )
!!!  ALLOCATE ( obs_rain_tmp(gall,nprocs), gus_rain_tmp(gall,nbv,nprocs) )
!!!  ALLOCATE ( obs_map(gall,lall),        gus_map(gall,lall,nbv)        )
!!!  ALLOCATE ( tmp_grd(gall,lall),        tmp_rgn(gall,lall)            )
!!!  ALLOCATE ( cont_m(nbv)                                              )
!!!
!!!  write(tbgname,'(a,a,i2.2,a,i2.2,a,i3.3,a)') trim(vecdir),'/vectorgrd_gl',&
!!!          glevel,'_rl',rlevel,'_lrad',int(locrad/1000),'.grd'
!!!  write(tbrname,'(a,a,i2.2,a,i2.2,a,i3.3,a)') trim(vecdir),'/vectorrgn_gl',&
!!!          glevel,'_rl',rlevel,'_lrad',int(locrad/1000),'.grd'
!!!  open(2,file=trim(tbgname),form='unformatted',access='direct',recl=gall*lall*4,action='read')
!!!  open(3,file=trim(tbrname),form='unformatted',access='direct',recl=gall*lall*4,action='read')
!!!  read(2,rec=1) obs_map(:,:)
!!!
!!!  ngrd = int( obs_map(1,1) )
!!!  ALLOCATE ( vec_rgn(gall,lall,ngrd), vec_grd(gall,lall,ngrd) )
!!!  do igrd=1,ngrd 
!!!    read(2,rec=1+igrd) ((tmp_grd(g,l),g=1,gall),l=1,lall)
!!!    read(3,rec=1+igrd) ((tmp_rgn(g,l),g=1,gall),l=1,lall)
!!!    do g=1,gall
!!!    do l=1,lall
!!!      vec_grd(g,l,igrd) = int ( tmp_grd(g,l) ) ! vector for grid
!!!      vec_rgn(g,l,igrd) = int ( tmp_rgn(g,l) ) ! vector for region
!!!    end do
!!!    end do
!!!  end do
!!!  close(2)
!!!  close(3)
!!!
!!!  do iread=1,nread
!!!    l = proc_l(myrank,iread)
!!!   ! read obs
!!!    write(obsname,'(2a,i10.10,a,a)') trim(obsdir),'/',cdate,'/',trim(rainname)
!!!    call MISC_make_idstr(oname,trim(obsname),'rgn',l)
!!!    open(1,file=trim(oname),form='unformatted',access='direct',recl=gall*4,action='read')
!!!      read(1,rec=1) obs_rain(:) ! [mm/6hr]
!!!    close(1)
!!!    do ibv=1,nbv
!!!      write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',cdate,'/',ibv,'/tppn'    
!!!      call MISC_make_idstr(gname,trim(gusname),'rgn',l)
!!!      if(l==1 .and. ibv==1)  write(6,'(2a)') "      Reading Guss. File  ",trim(gname)        
!!!      open(1,file=trim(gname),form='unformatted',access='direct',recl=gall*4,action='read')
!!!        read(1,rec=1) gus_rain(:,ibv) ! [mm/6hr]
!!!      close(1)
!!!    end do
!!!
!!!    CALL MPI_ALLGATHER(obs_rain,     gall,     MPI_REAL,  &
!!!                       obs_rain_tmp, gall,     MPI_REAL, MPI_COMM_WORLD, ierr)
!!!    CALL MPI_ALLGATHER(gus_rain,     gall*nbv, MPI_REAL,  &
!!!                       gus_rain_tmp, gall*nbv, MPI_REAL, MPI_COMM_WORLD, ierr)
!!!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!!
!!!    do iprocs=0,nprocs-1
!!!      l = proc_l(iprocs,iread)
!!!!        if( myrank==0 ) write(6,*) iprocs, " ===> ", l
!!!      obs_map(:,l)       = obs_rain_tmp(:,iprocs+1)     ! [mm/hr] & REAL(4)
!!!      do ibv=1,nbv
!!!        gus_map(:,l,ibv) = gus_rain_tmp(:,ibv,iprocs+1) ! [mm/hr] & REAL(4)
!!!      end do
!!!    enddo ! nprocs
!!!  enddo ! iread
!!!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!!
!!!  prate_o(:,:)   = 0.0d0
!!!  prate_m(:,:,:) = 0.0d0
!!!  do l=1,lall
!!!  do g=1,gall
!!!    cont_o    = 0
!!!    cont_m(:) = 0
!!!    do igrd=1,ngrd
!!!      gg = vec_grd(g,l,igrd)
!!!      ll = vec_rgn(g,l,igrd)
!!!      if( gg==0 .and. ll==0 ) EXIT
!!!      cont_o = cont_o + 1
!!!      if( obs_map(gg,ll)       > ppzero_thres ) prate_o(g,l)     = prate_o(g,l)     + 1.0d0
!!!      do ibv=1,nbv
!!!        cont_m(ibv) = cont_m(ibv) + 1
!!!        if( gus_map(gg,ll,ibv) > ppzero_thres ) prate_m(g,l,ibv) = prate_m(g,l,ibv) + 1.0d0
!!!      end do
!!!    end do ! igrd
!!!    if( cont_o > 4.5 ) then
!!!      prate_o(g,l)  = prate_o(g,l)  / real( cont_o, r_size )
!!!    else
!!!      prate_o(g,l)  = 0.5
!!!    end if 
!!!    do ibv=1,nbv
!!!      if( cont_m(ibv) > 4.5 ) then
!!!        prate_m(g,l,ibv)  = prate_m(g,l,ibv)  / real( cont_m(ibv), r_size )
!!!      else
!!!        prate_m(g,l,ibv)  = 0.5
!!!      end if       
!!!    end do
!!!  end do
!!!  end do
!!!
!!!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!!  !######  get ratio for treeatment of zero precip.

  qc(:,:)         = sngl( undef )
  elm(:,:)        = sngl( undef )
  gus_trns(:,:,:) = sngl( undef )
  obs_err(:,:)    = sngl( undef )
  obs_trns(:,:)   = sngl( undef )

  do iread=1,nread
    l = proc_l(myrank,iread)
    
    if( transform .ne. -1 ) then
      write(mskname,'(a,a,i2.2,a,i2.2,a)') trim(vecdir),'/gridusage_gl',glevel,'_rl',rlevel,'.grd'  
      open(1,file=trim(mskname),form='unformatted',access='direct',recl=gall*4,action='read') 
        read(1,rec=l) usage(:)
      close(1)
    else
      usage(:) = 0.0
    endif

    ppmask_m(:) = 1.0d0
    ppmask_o(:) = 1.0d0
    do g=1,gall
      if( ico_lon(g,l) .lt. 0.0d0 ) ico_lon(g,l) = ico_lon(g,l) + 360.0d0
    end do
    ! read obs
    write(obsname,'(2a,i10.10,a,a)') trim(obsdir),'/',cdate,'/',trim(rainname)
    call MISC_make_idstr(oname,trim(obsname),'rgn',l)
    if(l ==1)  write(6,'(2a)') "      Reading Obs.  File  ",trim(oname)
    open(1,file=trim(oname),form='unformatted',access='direct',recl=gall*4,action='read')
      read(1,rec=1) obs_rain(:) ! [mm/6hr]
    close(1)
    do ibv=1,nbv
      write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',cdate,'/',ibv,'/tppn'    
      call MISC_make_idstr(gname,trim(gusname),'rgn',l)
      !write(*,*) 'gname= ', trim(gname)
      if(l==1 .and. ibv==1)  write(6,'(2a)') "      Reading Guss. File  ",trim(gname)        
      open(1,file=trim(gname),form='unformatted',access='direct',recl=gall*4,action='read')
        read(1,rec=1) gus_rain(:,ibv) ! [mm/6hr]
      close(1)
    end do

    !===> joint probability histogram
    if( prdhst ) then
      hmsk_obs_rgn(:)   = 0.0
      hmsk_gus_rgn(:,:) = 0.0
    endif

    !===> quality control for ensemble member
    qc_rgn(:)  = 0.0    
    qc_mem(:)  = 0.0d0
    qc_god(:)  = 0.0
    elm_rgn(:) = real(id_rain_obs)

    do g=1,gall
      rain_mem = 0
      do ibv=1,nbv
        if( gus_rain(g,ibv) .ge. rain_min_rate ) rain_mem = rain_mem + 1
      end do
      prc_mem(g)  = real(rain_mem) / real(nbv)
      if( prc_mem(g) .ge. rain_min_mem ) then
        qc_mem(g) = 1.0
        if( transform .ne. -1) then ! without precip. assimilation 
          if( (obs_rain(g).ge.0.0d0) .and. (qc_mem(g)>0.5) .and. usage(g)>0.5) qc_rgn(g) = 1.0
        endif
      else
        
      endif
    end do

    !===> reading settlement
    if( transform .ne. -1 ) then ! without spinup-run
      call MISC_make_idstr(gname,trim(mem_cdfname),'rgn',l)
      call MISC_make_idstr(oname,trim(obs_cdfname),'rgn',l) 
      call read_corr(gname,gall,ndivid,corr_m) ! read Corralation (Equal)
      call read_corr(oname,gall,ndivid,corr_o) ! read Corralation (Equal)
      if ( opt_qccorr == 1 ) then
        do g=1,gall
          if( corr_m(g).lt.corr_min_qc ) qc_rgn(g) = 0.0d0
        end do
      endif
      if ( opt_exzero == 1)  then
        do g=1,gall
          if( obs_rain(g) < ppzero_thres ) qc_rgn(g) = 0.0d0
        end do
      endif
    endif

    !===> precip. transformation
    if( transform==2 .or. transform==3 ) then ! gaussian transformation
      ppmask_m(:) = 0.0d0
      ppmask_o(:) = 0.0d0     
      call read_ppcdf(gname,gall,ndivid,ppcdf_m,ppzero_m,ppsamp_m,ppmask_m,dfzero_m) ! read CDF MODEL
      call read_ppcdf(oname,gall,ndivid,ppcdf_o,ppzero_o,ppsamp_o,ppmask_o,dfzero_o) ! read CDF OBS.
      if( l==1 ) write(6,'(2a)') "      Rain2Trans:: Use Model CDF ==> ", trim(gname)
      do g=1,gall
        if( ppmask_m(g) < 0.5 ) qc_rgn(g)   = 0.0d0
        if( ppmask_o(g) < 0.5 ) qc_rgn(g)   = 0.0d0  
        if( ppmask_o(g) < 0.5 ) ppzero_o(g) = undef
        if( (opt_ppguess==1) .and. (ppzero_m(g)<=thres_ppguess) ) qc_rgn(g)   = 0.0d0
      end do
      if( opt_fixzero==1 ) then ! re-compute ensemble member qc
        do g=1,gall
          rain_mem = 0
          do ibv=1,nbv
            if( gus_rain(g,ibv) .ge. dfzero_m(g) ) rain_mem = rain_mem + 1
          end do
          prc_mem(g)  = real(rain_mem) / real(nbv)
          if( prc_mem(g) .ge. rain_min_mem ) then
            qc_mem(g) = 1.0
            if( transform .ne. -1) then ! without precip. assimilation 
              if( (obs_rain(g).ge.0.0d0) .and. (qc_mem(g)>0.5) .and. usage(g)>0.5 ) qc_rgn(g) = 1.0
            endif
          else
            qc_mem(g) = 0.0
            qc_rgn(g) = 0.0
          endif
        end do        
      end if
    end if

    do g=1,gall ! grid loop start
      gus_trns_rgn(g,:) = gus_rain(g,:)                                    ! [mm/6hr]
      obs_trns_rgn(g)   = obs_rain(g)                                      ! [mm/6hr]
      obs_err_rgn(g)    = obs_rain(g) * err_rain                           ! [mm/6hr]
      if( obs_err_rgn(g) .lt. err_min_rain ) obs_err_rgn(g) = err_min_rain ! [mm/6hr]  
      
      obs_err_org = obs_err_rgn(g)  ! original observation error

      if( transform == 1 ) then ! log transformation
        obs_err_p       = dlog( obs_err_org   + obs_trns_rgn(g) + rain_alpha ) &
                        - dlog(                 obs_trns_rgn(g) + rain_alpha ) ! [LOG(mm/6hr)] 
!!!        obs_err_n       = dlog(-obs_err_org   + obs_trns_rgn(g) + rain_alpha ) &
!!!                        + dlog(                 obs_trns_rgn(g) + rain_alpha ) ! [LOG(mm/6hr)]
        obs_err_rgn(g)  = sngl ( obs_err_p ) 

        obs_trns_rgn(g)   = log( obs_rain(g)   + rain_alpha )                 ! [LOG(mm/6hr)]                
        gus_trns_rgn(g,:) = log( gus_rain(g,:) + rain_alpha )                 ! [LOG(mm/6hr)]         
      else if( transform==2 .or. (transform==3 .and. qc_mem(g)<0.5) ) then ! Gaussian transformation with medium zero precipitation
        !---> transformation for ensemble member
        if( ppmask_m(g)>0.5 ) then           
          do ibv=1,nbv               
            gus_trns_rgn(g,ibv) = sngl( pptrans_normal( dble(gus_rain(g,ibv)), ppcdf_m(g,:), ppzero_m(g), dfzero_m(g) ) )
!            gus_trns_rgn(g,ibv) = sngl( pptrans_normal_prate( dble(gus_rain(g,ibv)), ppcdf_m(g,:), ppzero_m(g), dfzero_m(g), prate_m(g,l,ibv) ) )
          end do ! nbv
        else
            gus_trns_rgn(g,:) = sngl( undef )
        end if

        !---> transformation for observation and obs. error
        if (ppmask_o(g)>0.5 .and. obs_rain(g)>-0.000001 ) then
          !-> obs. variable
          obs_trns_rgn(g) = sngl( pptrans_normal( dble(obs_rain(g)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g) ) )
!            obs_trns_rgn(g) = sngl( pptrans_normal_prate( dble(obs_rain(g)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g),prate_o(g,l) ) )

          if( opt_ppobserr==1 ) then ! transformed
            obs_err_rgn(g) = sngl ( ppnormal_obserror( &
              dble(obs_rain(g)), dble(obs_err_org), dble(obs_trns_rgn(g)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g) ) )
          else  ! constant
            obs_err_rgn(g) = err_rain
          end if            
        else
          obs_trns_rgn(g) = sngl( undef )
          obs_err_rgn(g)  = sngl( undef )
        end if
      else if ( (transform==3 .and. qc_mem(g)>0.5) ) then ! Gaussian transformation with modified medium zero precipitation
        !---> transformation for ensemble member
        if( ppmask_m(g)>0.5 ) then           
          gus_rain_sub(g,:) = REAL(gus_rain(g,:),r_size)
          call pptrans_normal_mdzero_def (gus_rain_sub(g,:), ppcdf_m(g,:), ppzero_m(g), dfzero_m(g), zero_mem, ym, sigma)
          gus_trns_rgn(g,:) = sngl ( gus_rain_sub(g,:) )
        else
          gus_trns_rgn(g,:) = sngl( undef )
        end if

        !---> transformation for observation and obs. error
        if (ppmask_o(g)>0.5 .and. obs_rain(g)>-0.000001 ) then
          !-> obs. variable
          if ( ppmask_m(g)>0.5 ) then ! check
            obs_trns_rgn(g) = sngl ( pptrans_normal_mdzero( &
            dble(obs_rain(g)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g), ppzero_m(g), zero_mem, ym, sigma ) )

            if( opt_ppobserr==1 ) then ! transformed obs. error
!              obs_err_p = - dble( obs_trns_rgn(g) ) + pptrans_normal_mdzero( &
!                dble(obs_rain(g)+obs_err_org), ppcdf_o(g,:), ppzero_o(g), ppzero_m(g), zero_mem, ym, sigma )           
!              obs_err_n =   dble( obs_trns_rgn(g) ) - pptrans_normal_mdzero( &
!                dble(obs_rain(g)-obs_err_org), ppcdf_o(g,:), ppzero_o(g), ppzero_m(g), zero_mem, ym, sigma )
              obs_err_p = - dble( obs_trns_rgn(g) ) + pptrans_normal_mdzero( &
                dble(obs_rain(g)*(1.0d0+err_rain)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g), ppzero_m(g), zero_mem, ym, sigma )           
              obs_err_n =   dble( obs_trns_rgn(g) ) - pptrans_normal_mdzero( &
                dble(obs_rain(g)*(1.0d0-err_rain)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g), ppzero_m(g), zero_mem, ym, sigma )
              
              if( obs_err_p < min_ppobserr ) obs_err_p = min_ppobserr
              if( obs_err_n < min_ppobserr ) obs_err_n = min_ppobserr
              obs_err_rgn(g) = 0.5d0 * ( obs_err_p + obs_err_n )            
            else  ! constant
              obs_err_rgn(g) = err_rain
            end if

          else ! do as medium zero precip.
            obs_trns_rgn(g) = sngl( pptrans_normal( dble(obs_rain(g)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g) ) )
            if( opt_ppobserr==1 ) then ! transformed obs. error
              obs_err_rgn(g) = sngl ( ppnormal_obserror( &
                dble(obs_rain(g)), dble(obs_err_org), dble(obs_trns_rgn(g)), ppcdf_o(g,:), ppzero_o(g), dfzero_o(g) ) )
            else  ! constant
              obs_err_rgn(g) = err_rain
            end if 
          endif            
        else
          obs_trns_rgn(g) = sngl( undef )
          obs_err_rgn(g)  = sngl( undef )
        end if
      else
!        if(myrank==0)write(6,*) "not supported option:: ", transform, qc_mem(g)
!        stop            
      end if ! transform

      ! production of joint histgram
      if( transform==1 .or. transform==2 .or. transform==3 ) then ! log & gaussian transformation
      if( qc_rgn(g)>0.5     ) then ! qc includes ppmask_m and ppmask_o
        if( (obs_rain(g)>-0.000001) ) then ! observation
          hmsk_obs_rgn(g) = 1.0
          if( obs_rain(g)      .ge. dfzero_o(g) ) hmsk_obs_rgn(g) = 2.0 ! not zero precip.
          rain_mem = 0
          do ibv=1,nbv
            hmsk_gus_rgn(g,ibv) = 1.0
            if( gus_rain(g,ibv).ge. dfzero_m(g) ) hmsk_gus_rgn(g,ibv) = 2.0 ! not zero precip.
            
            if( hmsk_gus_rgn(g,ibv)>1.5 .and. hmsk_obs_rgn(g)>1.5) then            
              if( ppzero_o(g) .ge. ppzero_m(g) ) then
                hmsk_obs_rgn(g)       = 3.0 ! no bias region
                if( sigma2cdf( dble(gus_trns_rgn(g,ibv)) ) .ge. ppzero_o(g) )  then
                  hmsk_gus_rgn(g,ibv) = 3.0 ! no bias region
                  rain_mem = rain_mem + 1
                endif
              else
                rain_mem = rain_mem + 1
                hmsk_gus_rgn(g,ibv)   = 3.0 ! no bias region
                if( sigma2cdf( dble(obs_trns_rgn(g)) )     .ge. ppzero_m(g) )  &
                  hmsk_obs_rgn(g)     = 3.0 ! no bias region
              end if
            endif
          end do ! ibv

          if ( (real(rain_mem)/real(nbv)) .lt. rain_min_mem ) then
            if ( opt_exzero == 2 ) qc_rgn(g) = 0.0 ! exclude bias region
          else
            qc_god(g) = 1.0 ! good quality 
          end if           
        endif
      endif ! prdhst & qc
      endif ! transform == 1 or 2 or 3
    end do ! g: grid loop end

    !====> variable localization
    do g=1,gall
      if( qc_rgn(g)> 0.5 ) then
        if ( opt_variloc == 1 ) then
          if( -20.0<=ico_lat(g,l) .and. ico_lat(g,l)<=20.0 ) then ! tropic
            if( qc_god(g) > 0.5 ) then 
              call get_varloc( tropc_good , qc_rgn(g), elm_rgn(g) )
            else
              call get_varloc( tropc_norm , qc_rgn(g), elm_rgn(g) )
            end if
          else ! extratropic
            if( qc_god(g) > 0.5 ) then 
              call get_varloc( extrp_good , qc_rgn(g), elm_rgn(g) )
            else
              call get_varloc( extrp_norm , qc_rgn(g), elm_rgn(g) )
            end if            
          end if
        end if ! opt_variloc
        if ( opt_ppreduce == 1 ) then
          ix = g2xgrd(g)
          iy = g2ygrd(g)
          if( ico_lat(g,l) .gt. 20.0 ) then ! NPH
            if( mod(ix,ppred_nph)/=0 .or. mod(iy,ppred_nph)/=0 ) then ! masking out
              qc_rgn(g)  = 0.
!              elm_rgn(g) = real( undef )
            end if
          else if( ico_lat(g,l) .gt. -20.0 ) then ! TRP
            if( mod(ix,ppred_trp)/=0 .or. mod(iy,ppred_trp)/=0 ) then ! masking out
              qc_rgn(g)  = 0.
!              elm_rgn(g) = real( undef )
            end if
          else if( ico_lat(g,l) .ge. -90.0 ) then ! SPH 
            if( mod(ix,ppred_sph)/=0 .or. mod(iy,ppred_sph)/=0 ) then ! masking out
              qc_rgn(g)  = 0.
!              elm_rgn(g) = real( undef )
            end if
          end if
        end if ! opt_ppreduce
      else
        elm_rgn(g) = real( undef )
      end if
    end do    

!    write(6,'(i,5f)') l,(hmsk_obs_rgn(g),g=1,5)
!    write(6,'(i,5f)') l,(hmsk_gus_rgn(g,g),g=1,5)

    !=====> output generation of ensembles
    do ibv=1,nbv
      write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',cdate,'/',ibv,'/tppn_trns'    
      call MISC_make_idstr(gname,trim(gusname),'rgn',l)
!    if(myrank ==1)  write(6,'(a)') trim(gname)        
      open(1,file=trim(gname),form='unformatted',access='direct',recl=gall*4)
        write(1,rec=1) gus_trns_rgn(:,ibv) ! transformed precip.
      close(1)        
    end do      
    
    !===> write ensemble mean for gaussian transformation samples
    call calc_ensmean4(nbv,gall,gus_rain(:,:),gus_mean(:),1)
    write(gusname,'(2a,i10.10,a)') trim(gusdir),'/',cdate,'/tppn_mean'
    call MISC_make_idstr(gname,trim(gusname),'rgn',l)
    open(1,file=trim(gname),form='unformatted',access='direct',recl=gall*4)
      write(1,rec=1) gus_mean(:) ! [mm/6hr]      
    close(1)             

    !=====> output of monit transform
    call calc_ensmean4(nbv,gall,gus_trns_rgn(:,:),gus_trns_mean(:),1)
    write(outname,'(2a)')      trim(datdir),'/monit_transform'
    call MISC_make_idstr(fname,trim(outname),'rgn',l)
    open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4)
      write(1,rec=1)  gus_mean(:)         ! first guess [mm/6hr]
      write(1,rec=2)  obs_rain(:)         ! observation [mm/6hr]
      write(1,rec=3)  gus_trns_mean(:)    ! first guess [transformed unit]
      write(1,rec=4)  obs_trns_rgn(:)     ! observation [transformed unit]
      write(1,rec=5)  obs_err_rgn(:)      ! obs. error  [transformed unit]
      write(1,rec=6)  sngl( ppmask_m(:) ) ! mask for model transform [0 or 1]
      write(1,rec=7)  sngl( ppmask_o(:) ) ! mask for obs.  transform [0 or 1]
      write(1,rec=8)  qc_rgn(:)           ! quality check [0 or 1]
      write(1,rec=9)  prc_mem(:)*100.     ! quality check of precipitating member [%]
      write(1,rec=10) sngl(corr_m(:))     ! correlation between gues and obs.
      write(1,rec=11) qc_god(:)           ! good quality precipitation [0 or 1]
      write(1,rec=12) elm_rgn(:)          ! precipitation element (variable localization)
      write(1,rec=13) sngl( ppzero_m(:) ) ! zero-precip. probavility (model)
      write(1,rec=14) sngl( ppzero_o(:) ) ! zero-precip. probavility (obs.)
    close(1)

    !====> aallgather mpi process
!    write(6,'(2i,3f)') myrank, l, qc_rgn(1), qc_rgn(500) , qc_rgn(gall)
!    write(6,'(2i,3f)') myrank, l, obs_trns_rgn(1), obs_trns_rgn(500) , obs_trns_rgn(gall)
!    write(6,'(2i,3f)') myrank, l, gus_trns_rgn(1,1), gus_trns_rgn(gall,1) , gus_trns_rgn(gall,nbv)

    CALL MPI_ALLGATHER(qc_rgn,        gall,     MPI_REAL,  &
                       qc_tmp,        gall,     MPI_REAL, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(elm_rgn,       gall,     MPI_REAL,  &
                       elm_tmp,       gall,     MPI_REAL, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(obs_err_rgn,   gall,     MPI_REAL,  &
                       obs_err_tmp,   gall,     MPI_REAL, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(obs_trns_rgn,  gall,     MPI_REAL,  &
                       obs_trns_tmp,  gall,     MPI_REAL, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(gus_trns_rgn,  gall*nbv, MPI_REAL,  &
                       gus_trns_tmp,  gall*nbv, MPI_REAL, MPI_COMM_WORLD, ierr)
    if ( prdhst ) then
    CALL MPI_ALLGATHER(hmsk_obs_rgn,  gall,     MPI_REAL,  &
                       hmsk_obs_tmp,  gall,     MPI_REAL, MPI_COMM_WORLD, ierr)     
    CALL MPI_ALLGATHER(hmsk_gus_rgn,  gall*nbv, MPI_REAL,  &
                       hmsk_gus_tmp,  gall*nbv, MPI_REAL, MPI_COMM_WORLD, ierr)  
    endif


    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    do iprocs=0,nprocs-1
      l = proc_l(iprocs,iread)
      qc(:,l)       =     qc_tmp(:,iprocs+1)
      elm(:,l)      =     elm_tmp(:,iprocs+1)
      obs_trns(:,l) =     obs_trns_tmp(:,iprocs+1)
      obs_err(:,l)  =     obs_err_tmp(:,iprocs+1)

      do ibv=1,nbv
        gus_trns(:,l,ibv) = gus_trns_tmp(:,ibv,iprocs+1)
      end do
      if( prdhst ) then
        hmsk_obs(:,l)       = hmsk_obs_tmp(:,iprocs+1)
        do ibv=1,nbv
          hmsk_gus(:,l,ibv) = hmsk_gus_tmp(:,ibv,iprocs+1) 
        end do
!        if( myrank == 0 ) write(6,'(a,i,5f)') "=====> ",l,(hmsk_obs(g,l),g=1,5)
!        if( myrank == 0 ) write(6,'(a,i,5f)') "=====> ",l,(hmsk_gus(g,l,g),g=1,5)
      end if
    enddo ! nprocs
  end do ! iread to paralellize region loop 

!  if( myrank == nprocs-1) then
!    write(6,*) "=================================================================="
!    do l=1,lall
!      write(6,'(2i,3f)') myrank, l, qc(1,l), qc(500,l) , qc(gall,l)
!      write(6,'(2i,3f)') myrank, l, obs_trns(1,l), obs_trns(500,l) , obs_trns(gall,l)
!      write(6,'(2i,3f)') myrank, l, gus_trns(1,l,1), gus_trns(gall,l,1) , gus_trns(gall,l,nbv)
!    end do
!  end if

  ! ???? what to do ?? 20150528
  obs_err_p = 0.0d0
  obs_err_n = 0.0d0
  do l=1,lall
  do g=1,gall
  if( obs_trns(g,l) >undef+1.0d0 ) then 
    obs_err_p = obs_err_p + 1.0d0
    if (qc(g,l)>0.5) obs_err_n = obs_err_n + 1.0d0
  end if
  end do
  end do

  ! finalize output
  if( myrank == 0 ) then

  !===> "opt_modbias" option
    if ( opt_modbias == 1 ) &
    call mod_bias(gall,lall,nbv,hmsk_obs,hmsk_gus,obs_trns,gus_trns) 

    do ibv=1,nbv
      rain_mem=0
      write(gusname,'(2a,i6.6,a)') trim(yhxdir),'/tppn_obs04',ibv,'.dat'
      open(1,file=trim(gusname),form='unformatted',access='sequential')         
      write(gusname,'(2a,i6.6,a)') trim(yhxdir),'/tppn_obs04',ibv,'.txt'
      open(100,file=trim(gusname))         
      do l=1,lall
      do g=1,gall
        if( qc(g,l) > 0.5 ) then
          rain_mem = rain_mem + 1
          write(1) elm(g,l), sngl(ico_lon(g,l)),sngl(ico_lat(g,l)),sngl(base_obsv_rain), &
                    obs_trns(g,l),obs_err(g,l),gus_trns(g,l,ibv),qc(g,l), typ
          write(100,'(8f12.5)') elm(g,l), sngl(ico_lon(g,l)),sngl(ico_lat(g,l)),sngl(base_obsv_rain), &
                    obs_trns(g,l),obs_err(g,l),gus_trns(g,l,ibv),qc(g,l)
!        if( myrank==0 .and. ibv==1) write(6,'(2i,8f)') g,l,  &
!                 elm(g,l), sngl(ico_lon(g,l)),sngl(ico_lat(g,l)),sngl(base_obsv_rain), &
!                  obs_trns(g,l),obs_err(g,l),gus_trns(g,l,ibv),qc(g,l)
                    
        end if
      end do
      end do
      if( ibv==1 ) print *, "     rain2trans ======> num. precip ::",rain_mem
      close(1)
      close(100)
    end do ! ibv  

    if ( prdhst )call get_jhist(gall,lall,nbv,ncell,hmsk_obs,hmsk_gus,hmin,hmax,obs_trns,gus_trns,datdir,opt_modbias)  
  endif  ! myrank



  call MPI_Finalize(ierr)
  if( myrank==0 ) write(6,*) 'Finish Program rain2trans'    
end program

!=========================================================
subroutine calc_ensmean4(nbv,ngrd,ens,mean,nlev)
!=========================================================
  implicit none
  integer nbv,ngrd,nlev,igrd,ilev,ibv
  REAL(4) ens(ngrd,nbv,nlev), mean(ngrd,nlev)

  do igrd=1,ngrd
  do ilev=1,nlev
    mean(igrd,ilev) = ens(igrd,1,ilev)
    do ibv=2,nbv
      mean(igrd,ilev) = mean(igrd,ilev) + ens(igrd,ibv,ilev)
    end do
    mean(igrd,ilev) = mean(igrd,ilev) / real(nbv)
  end do
  end do

end
!=========================================================
subroutine mod_bias(gall,lall,nbv,hmsk_obs,hmsk_gus,obs_trns,gus_trns)
!=========================================================
  implicit none
  integer, intent(in)     :: gall, lall, nbv
  real(4), intent(in)     :: hmsk_obs(gall,lall), hmsk_gus(gall,lall,nbv)
  real(4), intent(inout)  :: obs_trns(gall,lall), gus_trns(gall,lall,nbv)

  integer, parameter :: nhst = 3
  integer g,l,ibv,icell,ihst
  real(4) :: cnt(nhst)
  real(8) :: aveobs(nhst),avegus(nhst),aveyhx(nhst) 

  cnt(:)     = 0.0
  aveobs(:)  = 0.0d0
  avegus(:)  = 0.0d0
  do l=1,lall
  do g=1,gall
    if( hmsk_obs(g,l) > 0.5 ) then
      do ibv=1,nbv
        if( hmsk_gus(g,l,ibv) > 0.5 ) then
          if( hmsk_obs(g,l)>2.5 .and. hmsk_gus(g,l,ibv)>2.5 ) then
            cnt(3)    = cnt(3)    + 1.0
            aveobs(3) = aveobs(3) + obs_trns(g,l)
            avegus(3) = avegus(3) + gus_trns(g,l,ibv)
          else if( hmsk_obs(g,l)>1.5 .and. hmsk_gus(g,l,ibv)>1.5 ) then
            cnt(2)    = cnt(2)    + 1.0
            aveobs(2) = aveobs(2) + obs_trns(g,l)
            avegus(2) = avegus(2) + gus_trns(g,l,ibv)
          else
            cnt(1)    = cnt(1)    + 1.0
            aveobs(1) = aveobs(1) + obs_trns(g,l)
            avegus(1) = avegus(1) + gus_trns(g,l,ibv)
          end if
        end if ! hmsk_gus
      end do ! ibv
    endif ! hmsk_obs
  end do ! g
  end do ! l

  aveobs(:) = aveobs(:) / cnt(:)
  avegus(:) = avegus(:) / cnt(:)
  aveyhx(:) = aveobs(:) - avegus(:)

  do l=1,lall
  do g=1,gall
    if( hmsk_obs(g,l) > 0.5 ) then
      do ibv=1,nbv
        if( hmsk_gus(g,l,ibv) > 0.5 ) then
          if( hmsk_obs(g,l)>2.5 .and. hmsk_gus(g,l,ibv)>2.5 ) then
            gus_trns(g,l,ibv) = gus_trns(g,l,ibv) + aveyhx(3)
          else if( hmsk_obs(g,l)>1.5 .and. hmsk_gus(g,l,ibv)>1.5 ) then
            gus_trns(g,l,ibv) = gus_trns(g,l,ibv) + aveyhx(2)
          else
            gus_trns(g,l,ibv) = gus_trns(g,l,ibv) + aveyhx(1)
          end if
        end if ! hmsk_gus
      end do ! ibv
    endif ! hmsk_obs
  end do ! g
  end do ! l  


end
!=========================================================
subroutine get_jhist(gall,lall,nbv,ncell,hmsk_obs,hmsk_gus,hmin,hmax,obs_trns,gus_trns,datdir)
!=========================================================
  implicit none
  integer, intent(in)  :: gall, lall, nbv, ncell
  real(4), intent(in)  :: hmsk_obs(gall,lall), hmsk_gus(gall,lall,nbv)
  real(4), intent(in)  :: obs_trns(gall,lall), gus_trns(gall,lall,nbv)
  real(8), intent(in)  :: hmin(ncell), hmax(ncell)
  character(128), intent(in) :: datdir

  integer, parameter :: nhst = 5
  integer g,l,ibv,icell,ihst
  integer hobs,hgus
  character(255) :: outname
  real(4) :: cnt(nhst), hst(ncell,ncell,nhst) 
  real(8) :: aveobs(nhst),avegus(nhst),aveyhx(nhst) 
  real(8) :: var_obs(nhst),var_gus(nhst),corr(nhst) 
!!! 1: all precip., 2: zero precip.

!!!!!! start
  cnt(:)     = 0.0
  hst(:,:,:) = 0.0
  aveobs(:)  = 0.0d0
  avegus(:)  = 0.0d0
  do l=1,lall
  do g=1,gall
    if( hmsk_obs(g,l) > 0.5 ) then
          call get_binhist(ncell,hmin,hmax,obs_trns(g,l),    hobs)
      do ibv=1,nbv
        if( hmsk_gus(g,l,ibv) > 0.5 ) then
          call get_binhist(ncell,hmin,hmax,gus_trns(g,l,ibv),hgus)

          if( (1.le.hobs).and.(hobs.le.ncell).and.(1.le.hgus).and.(hgus.le.ncell) ) &
          hst(hobs,hgus,1) = hst(hobs,hgus,1) + 1.0
          cnt(1)           = cnt(1)           + 1.0
          aveobs(1)        = aveobs(1)        + dble( obs_trns(g,l) )
          avegus(1)        = avegus(1)        + dble( gus_trns(g,l,ibv) )

          if( hmsk_obs(g,l)>1.5 .and. hmsk_gus(g,l,ibv)>1.5 ) then ! w/o zero precip.
            if( (1.le.hobs).and.(hobs.le.ncell).and.(1.le.hgus).and.(hgus.le.ncell) ) &
            hst(hobs,hgus,2) = hst(hobs,hgus,2) + 1.0
            cnt(2)           = cnt(2)           + 1.0
            aveobs(2)        = aveobs(2)        + dble( obs_trns(g,l) )
            avegus(2)        = avegus(2)        + dble( gus_trns(g,l,ibv) )
            if( hmsk_obs(g,l)>2.5 .and. hmsk_gus(g,l,ibv)>2.5 ) then ! no bias region
              if( (1.le.hobs).and.(hobs.le.ncell).and.(1.le.hgus).and.(hgus.le.ncell) ) &
              hst(hobs,hgus,3) = hst(hobs,hgus,3) + 1.0
              cnt(3)           = cnt(3)           + 1.0
              aveobs(3)        = aveobs(3)        + dble( obs_trns(g,l) )
              avegus(3)        = avegus(3)        + dble( gus_trns(g,l,ibv) )
            end if
          endif

!          if( hmsk_obs(g,l)>2.5 .and. hmsk_gus(g,l,ibv)>2.5 ) then
!
!          else if( hmsk_obs(g,l)>1.5 .and. hmsk_gus(g,l,ibv)>1.5 ) then
!            if( (1.le.hobs).and.(hobs.le.ncell).and.(1.le.hgus).and.(hgus.le.ncell) ) &
!            hst(hobs,hgus,4) = hst(hobs,hgus,4) + 1.0
!            cnt(4)           = cnt(4)           + 1.0
!            aveobs(4)        = aveobs(4)        + dble( obs_trns(g,l) )
!            avegus(4)        = avegus(4)        + dble( gus_trns(g,l,ibv) )  
!          else 
!            if( (1.le.hobs).and.(hobs.le.ncell).and.(1.le.hgus).and.(hgus.le.ncell) ) &
!            hst(hobs,hgus,5) = hst(hobs,hgus,5) + 1.0
!            cnt(5)           = cnt(5)           + 1.0
!            aveobs(5)        = aveobs(5)        + dble( obs_trns(g,l) )
!            avegus(5)        = avegus(5)        + dble( gus_trns(g,l,ibv) )  
!          end if ! hmsk

        end if ! hmsk_gus
      end do
    endif ! hmsk_obs
  end do
  end do

  write(outname,'(2a,i6.6,a)') trim(datdir),'/joint_histgram.bin'
  open(1,file=trim(outname),access='direct',form='unformatted',recl=ncell*ncell*4)
  do ihst=1,nhst
    hst(:,:,ihst) = hst(:,:,ihst) / cnt(ihst)
    write(1,rec=ihst) hst(:,:,ihst)
    
    aveobs(ihst) = aveobs(ihst) / dble( cnt(ihst) )
    avegus(ihst) = avegus(ihst) / dble( cnt(ihst) )
  end do
  close(1)

  !===> calc correlation
  var_obs(:) = 0.0d0
  var_gus(:) = 0.0d0
  corr(:)    = 0.0d0
  do l=1,lall
  do g=1,gall
    if( hmsk_obs(g,l) > 0.5 ) then
      do ibv=1,nbv
        if( hmsk_gus(g,l,ibv) > 0.5 ) then
          corr(1)    = corr(1)    + ( dble(obs_trns(g,l))     - aveobs(1) )*( dble(gus_trns(g,l,ibv)) - avegus(1) )
          var_obs(1) = var_obs(1) + ( dble(obs_trns(g,l))     - aveobs(1) )**2.0d0
          var_gus(1) = var_gus(1) + ( dble(gus_trns(g,l,ibv)) - avegus(1) )**2.0d0
          if( hmsk_obs(g,l)>1.5 .and. hmsk_gus(g,l,ibv)>1.5 ) then ! w/o zero precip.
            corr(2)    = corr(2)    + ( dble(obs_trns(g,l))     - aveobs(2) )*( dble(gus_trns(g,l,ibv)) - avegus(2) )
            var_obs(2) = var_obs(2) + ( dble(obs_trns(g,l))     - aveobs(2) )**2.0d0
            var_gus(2) = var_gus(2) + ( dble(gus_trns(g,l,ibv)) - avegus(2) )**2.0d0
            if( hmsk_obs(g,l)>2.5 .and. hmsk_gus(g,l,ibv)>2.5 ) then ! no bias region
              corr(3)    = corr(3)    + ( dble(obs_trns(g,l))     - aveobs(3) )*( dble(gus_trns(g,l,ibv)) - avegus(3) )
              var_obs(3) = var_obs(3) + ( dble(obs_trns(g,l))     - aveobs(3) )**2.0d0
              var_gus(3) = var_gus(3) + ( dble(gus_trns(g,l,ibv)) - avegus(3) )**2.0d0          
            endif
          endif
        endif
      end do
    end if
  end do
  end do

  do ihst=1,nhst
  corr(ihst) = corr(ihst) / dsqrt( var_obs(ihst) ) / dsqrt( var_gus(ihst) )
  end do

  print '(a,5I9  )', "      Counter     Num (1) all; (2) w/o zero P; (3) no-bias region ",(int(cnt(ihst)),            ihst=1,nhst)
  print '(a,5f9.3)', "      Correlation R^2 (1) all; (2) w/o zero P; (3) no-bias region ",(corr(ihst)**2.0d0,         ihst=1,nhst)
  print '(a,5f9.3)', "      Average (Y-Hxf) (1) all; (2) w/o zero P; (3) no-bias region ",(aveobs(ihst)-avegus(ihst), ihst=1,nhst)

end
!=========================================================
subroutine get_binhist(ncell,hmin,hmax,vinp,hout)
!=========================================================
  implicit none
  integer, intent(in)  :: ncell
  real(4), intent(in)  :: vinp
  real(8), intent(in)  :: hmin(ncell),hmax(ncell)
  integer, intent(out) :: hout
  integer icell

  hout=0
  do icell=1,ncell
    if( hmin(icell).le.vinp .and. vinp.lt.hmax(icell) ) then
      hout = icell
      goto 10
    end if
  end do
10 continue

end
