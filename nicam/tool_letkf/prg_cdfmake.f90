!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  prg_cdfmake
!
!  program for generating CDF of obs. & first-guess precipitation
!
!  created                Apr 2015, Shunji Kotsuki, RIKEN-AICS
!  adopted to NICAM-LETKF May 2015, Shunji Kotsuki, RIKEN-AICS 
!
!-------------------------------------------------------------------------------

program cdfmake

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
 
  use fnc_transgauss
  implicit none

  integer :: g, l , gg, ll
  integer :: nij0, ngpv 
  integer :: igrd, ngrd

  REAL(4),ALLOCATABLE :: obs_rain_tmp(:,:),gus_rain_tmp(:,:),mem_rain_tmp(:,:,:)
  REAL(4),ALLOCATABLE :: obs_rain_rgn(:),  gus_rain_rgn(:),  mem_rain_rgn(:,:)
  REAL(8),ALLOCATABLE :: obs_rain(:,:,:),  gus_rain(:,:,:),  mem_rain(:,:,:,:)
  REAL(8),ALLOCATABLE :: obs_mask(:),      gus_mask(:),      mem_mask(:)
  REAL(8),ALLOCATABLE :: obs_maskr(:),     gus_maskr(:),     mem_maskr(:)
  REAL(8),ALLOCATABLE :: obs_zero(:),      gus_zero(:),      mem_zero(:)
  REAL(4),ALLOCATABLE :: rgn_tmp(:,:),     grd_tmp(:,:)
  REAL(4),ALLOCATABLE :: usage(:,:)

  REAL(8),ALLOCATABLE :: obs_pdf(:,:),  gus_pdf(:,:),  mem_pdf(:,:)    ! PDF
  REAL(8),ALLOCATABLE :: obs_cdf(:,:),  gus_cdf(:,:),  mem_cdf(:,:)    ! CDF
  REAL(8),ALLOCATABLE :: obs_cdfr(:,:), gus_cdfr(:,:), mem_cdfr(:,:)   ! CDF without rain
  REAL(8),ALLOCATABLE :: obs_norain(:), gus_norain(:), mem_norain(:)   ! Z of no rain

  REAL(8),ALLOCATABLE :: corr_n(:), corr_g(:)
  
  INTEGER,ALLOCATABLE :: obs_smp(:),    gus_smp(:),    mem_smp(:)   
  INTEGER,ALLOCATABLE :: obs_smpr(:),   gus_smpr(:),   mem_smpr(:)
  INTEGER, ALLOCATABLE :: vec_rgn(:,:,:), vec_grd(:,:,:) 
  
  REAL(8), ALLOCATABLE :: obs_org(:),  gus_org(:), mem_org(:)
  REAL(8), ALLOCATABLE :: obs_srt(:),  gus_srt(:), mem_srt(:)
  REAL(8), ALLOCATABLE :: obs_dat(:),  mem_dat(:,:)
  REAL(8), ALLOCATABLE :: obs_trns(:), mem_trns(:,:)
  INTEGER, ALLOCATABLE :: gus_int(:), obs_int(:), mem_int(:)
  
  REAL(8), ALLOCATABLE :: rmin(:),rmax(:)

  integer ibin, icon, idiv, ibv
  integer isamp, nsamp, idata, idata_mem, ndata_loc, ndata_mem
  integer gall, lall
  integer idate,ndate
  character(255) fname, mapname,tbgname,tbrname,mskname
  character(255) obsname, oname
  character(255) gusname, gname
  character(255) memname, mname
  double precision :: dist_zero           ! collection radius to PDF/CDF [m]
  logical exrgn,exgrd,exmsk
  
  !==> nhm_driver.cnf
  integer glevel
  integer rlevel
  integer vlayer
  character(128) rgnmngfname
  
  !==> rain2gauss.cnf
  !&rain2gauss_param

  character(128) :: snpdir=''

  integer :: nmem                    ! ensemble member
  integer :: gsdate                  ! initial date of gaussian transform
  integer :: gedate                  ! ending  date of gaussian transform
  integer :: gcdate                  ! current date of nicam-letkf
  integer :: nhour                   ! assimilation cycle
  integer :: binnum                  ! the number of PDF bins    
  double precision :: binsiz         ! size of the PDF bins [mm/6hr]
  double precision :: minnorain      ! minimum value of the no rain [CDF]
  logical :: prdpdf                  ! production of PDF
  logical :: prdcdf                  ! production of CDF
  logical :: svpcdf                  ! save PDF and CDF 

  double precision :: tmp_min_rate 
  !=> defined by fnc_transgauss.mod
!  character(128) :: obsdir='' 
!  character(128) :: gusdir=''
!  character(128) :: obs_cdfname=''
!  character(128) :: gus_cdfname=''
!  character(128) :: mem_cdfname=''
!  integer :: ndivid                  ! the number of diviation of CDF
!  double precision :: rain_min_rate  ! definition of precipitation [mm/6hr]
!  double precision :: locrad         ! localization radius [m], should be equal to sigma_obs in letkf_obs.f90 

  !--> paralell computing
  integer :: ierr, iread, nread, iprocs, nprocs, irank, myrank
  integer, allocatable :: proc_l(:,:) 
  real :: rtimer00, rtimer

!========================================================================
  CALL ADM_proc_init(ADM_MULTI_PRC)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)  

  if( myrank==0 ) rtimer00=mpi_wtime()
!========================================================================

  namelist / rain2gauss_param / &
    obsdir, gusdir, snpdir, vecdir, &
    obs_cdfname, gus_cdfname, &
    mem_cdfname, rainname,    &
    opt_qccorr, opt_fixzero,  &
    nmem,                     &
    gsdate, gedate, gcdate,   &
    nhour,                    &
    locrad, rain_min_rate,    &
    stp_fixzero, alw_fixzero, &
    minnorain,ndivid,         &
    binnum, binsiz,           &
    prdpdf, prdcdf, svpcdf
    
  open(1,file='rain2gauss.cnf')
  read(1,nml=rain2gauss_param)
  close(1)
    
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

  call get_gaussvariable()
   ! some of the variables cannot set because of originaly for ran2trans & trans2rain
  

  !=====> prep. for parallel computing
  if( mod(lall,nprocs).ne.0 ) then
    write(6,*) "lall should be devided by the nprocs"
    stop
  endif
  
  nread = lall / nprocs
  allocate ( proc_l(0:nprocs-1,nread) )

  iprocs = -1
  iread  = 1
  do l=1,lall
    iprocs = iprocs + 1
    if( iprocs == nprocs ) then
      iprocs = 0
      iread  = iread + 1
    endif
    proc_l(iprocs,iread)    = l
  enddo

   
  !====> Generation/Reading of Reference Table
  ALLOCATE( rgn_tmp(gall,lall), grd_tmp(gall,lall), usage(gall,lall) )

  write(mapname,'(a,a,i2.2,a,i2.2,a)') trim(vecdir),'/icolatlon_gl',glevel,'_rl',rlevel,'.dat'
  write(tbgname,'(a,a,i2.2,a,i2.2,a,i3.3,a)') trim(vecdir),'/vectorgrd_gl',&
          glevel,'_rl',rlevel,'_lrad',int(locrad/1000),'.grd'
  write(tbrname,'(a,a,i2.2,a,i2.2,a,i3.3,a)') trim(vecdir),'/vectorrgn_gl',&
          glevel,'_rl',rlevel,'_lrad',int(locrad/1000),'.grd'  
  write(mskname,'(a,a,i2.2,a,i2.2,a      )') trim(vecdir),'/gridusage_gl',&
          glevel,'_rl',rlevel,'.grd'  

  INQUIRE(FILE=trim(tbgname),EXIST=exgrd)         
  INQUIRE(FILE=trim(tbrname),EXIST=exrgn)
  INQUIRE(FILE=trim(mskname),EXIST=exmsk)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  if( (.not. exgrd) .or. (.not.exrgn) .or. (.not.exmsk)) then ! read the vector file
    open(1,file=trim(mapname),form='unformatted',access='sequential',action='read')
    open(2,file=trim(tbgname),form='unformatted',access='direct',recl=gall*lall*4,action='write')
    open(3,file=trim(tbrname),form='unformatted',access='direct',recl=gall*lall*4,action='write')
    open(4,file=trim(mskname),form='unformatted',access='direct',recl=gall*lall*4,action='write')
    call make_vector(gall,lall,locrad,myrank,nprocs,nread,proc_l)
    close(1)
    close(2)
    close(3)
    close(4)
    if( myrank==0 ) print *, "  Make Smapling Vector   ",trim(tbgname)
    if( myrank==0 ) print *, "  Make Smapling Vector   ",trim(tbrname)
    if( myrank==0 ) print *, "  Make Usage (Overlap)   ",trim(mskname)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  else
    if( myrank==0 ) print *, "  Read Smapling Vector   ",trim(tbgname)
    if( myrank==0 ) print *, "  Read Smapling Vector   ",trim(tbrname)
    if( myrank==0 ) print *, "  Read Usage (Overlap)   ",trim(mskname)
  end if

  open(1,file=trim(mapname),form='unformatted',access='sequential',action='read')
  open(2,file=trim(tbgname),form='unformatted',access='direct',recl=gall*lall*4,action='read')
  open(3,file=trim(tbrname),form='unformatted',access='direct',recl=gall*lall*4,action='read')
  open(4,file=trim(mskname),form='unformatted',access='direct',recl=gall*lall*4,action='read')

  read(2,rec=1) rgn_tmp(:,:)
  ngrd = int( rgn_tmp(1,1) )
  ALLOCATE ( vec_rgn(gall,lall,ngrd), vec_grd(gall,lall,ngrd) )
    do igrd=1,ngrd 
      read(2,rec=1+igrd) ((grd_tmp(g,l),g=1,gall),l=1,lall)
      read(3,rec=1+igrd) ((rgn_tmp(g,l),g=1,gall),l=1,lall)
      do g=1,gall
      do l=1,lall
        vec_grd(g,l,igrd) = int ( grd_tmp(g,l) ) ! vector for grid
        vec_rgn(g,l,igrd) = int ( rgn_tmp(g,l) ) ! vector for region
      end do
      end do
    end do
  read(4,rec=1) usage(:,:)


  if(myrank==0) call check_vector(gall,lall,ngrd,vec_grd,vec_rgn)
!  if(myrank==0) then
!    close(1)
!    open(1,file=trim(mapname),form='unformatted',access='sequential',action='read')
!    call check_usage(gall,lall,usage)  
!  end if
  close(1)
  close(2)
  close(3)
  close(4)

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(END MAKE VECTOR):',rtimer-rtimer00,' // NPROCS ',nprocs
  
  !=====> set up for PDF  
  ALLOCATE ( rmin(binnum), rmax(binnum) )
  
  do ibin=1,binnum
    rmin(ibin) = binsiz * dble( ibin - 1 )
    rmax(ibin) = binsiz * dble( ibin - 0 )
!    if(ADM_prc_me==1) write(6,'(i,2f)') ibin,rmin(ibin),rmax(ibin)
  end do
  
  !=====> get ndate
  idate = gsdate ! initial time
  ndate = 0
  do
    ndate = ndate + 1           
!    if(ADM_prc_me==1) write(6,'(a,3i)') "     CDFMAKE:: CDATE, NCON, EDATE",idate,ndate,gedate
    if( idate == gedate ) goto 10
    call update_date(idate,nhour)
    if ( ndate .ge. 10000 ) then
      write(6,*) 'ndate is too large'
      stop
    end if
  end do
10 continue

  !=====> read obs. and guess

  ALLOCATE( obs_rain(gall, lall, ndate), gus_rain(gall, lall, ndate) )
  ALLOCATE( mem_rain(gall, lall, nmem, ndate) )

  ALLOCATE( obs_rain_rgn(gall),      obs_rain_tmp(gall,nprocs)      )
  ALLOCATE( gus_rain_rgn(gall),      gus_rain_tmp(gall,nprocs)      )
  ALLOCATE( mem_rain_rgn(gall,nmem), mem_rain_tmp(gall,nmem,nprocs) )

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  idate = gsdate ! initial time
  ndate = 0
  obs_rain(:,:,:)   = undef
  gus_rain(:,:,:)   = undef
  mem_rain(:,:,:,:) = undef

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(STS RECORD READ):',rtimer-rtimer00,' // NPROCS ',nprocs

  do
    ndate = ndate + 1           
    if(ADM_prc_me==1) write(6,'(a,3i)') "     Read... CDFMAKE:: CDATE, NCON, EDATE",idate,ndate,gedate    

    do iread=1,nread
      l = proc_l(myrank,iread)
      write(obsname,'(2a,i10.10,a,a)') trim(obsdir),'/',idate,'/',trim(rainname)
      call MISC_make_idstr(oname,trim(obsname),'rgn',l)
      !if( l==1 ) write(6,'(2a)') "reading.. obs.",trim(oname)
      open(1,file=trim(oname),form='unformatted',access='direct',recl=gall*4,action='read')
        read(1,rec=1) obs_rain_rgn(:)
      close(1)
!      call smfilter( gall, obs_rain_rgn )

     
      !write(gusname,'(2a,i10.10,a)') trim(gusdir),'/',idate,'/tppn_mean'
      write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',idate,'/000001/tppn'
      write(*,*) 'gusname= ', trim(gusname)
      call MISC_make_idstr(fname,trim(gusname),'rgn',l)
      open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4,action='read')
        read(1,rec=1) gus_rain_rgn(:)
      close(1)    

      do ibv=1,nmem
        write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',idate,'/',ibv,'/tppn'
        call MISC_make_idstr(fname,trim(gusname),'rgn',l)
!       if(myrank==0)  write(6,'(i,a)') l,trim(fname)
        open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4,action='read')
          read(1,rec=1) mem_rain_rgn(:,ibv)
        close(1)        
      enddo

!      write(6,'(2i,3f)') myrank, l, gus_rain_rgn(1), gus_rain_rgn(500) , gus_rain_rgn(gall)
!      write(6,'(2i,3f)') myrank, l, obs_rain_rgn(1), obs_rain_rgn(500) , obs_rain_rgn(gall) 
!      write(6,'(2i,3f)') myrank, l, mem_rain_rgn(1,1), mem_rain_rgn(gall,1), mem_rain_rgn(gall,nmem)

      CALL MPI_ALLGATHER(obs_rain_rgn, gall,      MPI_REAL,  &
                         obs_rain_tmp, gall,      MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(gus_rain_rgn, gall,      MPI_REAL,  &
                         gus_rain_tmp, gall,      MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(mem_rain_rgn, gall*nmem, MPI_REAL,  &
                         mem_rain_tmp, gall*nmem, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do iprocs=0,nprocs-1
        l = proc_l(iprocs,iread)
!        if( myrank==0 ) write(6,*) iprocs, " ===> ", l
        obs_rain(:,l,ndate) = dble( obs_rain_tmp(:,iprocs+1) )       ! [mm/hr] & REAL(4)
        gus_rain(:,l,ndate) = dble( gus_rain_tmp(:,iprocs+1) )       ! [mm/hr] & REAL(4)
        do ibv=1,nmem
          mem_rain(:,l,ibv,ndate) = dble( mem_rain_tmp(:,ibv,iprocs+1) ) ! [mm/hr] & REAL(4)
        end do
      enddo ! nprocs
    end do ! nread

!    if( myrank==nprocs-1 ) then
!      do l=1,lall
!        write(6,'(2i,3f)') myrank, l, sngl(gus_rain(1,l,ndate)), sngl(gus_rain(500,l,ndate)), sngl(gus_rain(gall,l,ndate))
!        write(6,'(2i,3f)') myrank, l, sngl(mem_rain(1,l,1,ndate)), sngl(mem_rain(gall,l,1,ndate)), sngl(mem_rain(gall,l,nmem,ndate))
!      enddo
!    endif

    if( idate == gedate ) goto 20
    call update_date(idate,nhour)
  end do
20 continue

  !====> production of PDF/CDF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ndata_loc = ndate * ngrd
  ndata_mem = ndate * ngrd * nmem
  ALLOCATE ( gus_org(ndata_loc) , obs_org(ndata_loc), mem_org(ndata_mem) )
  ALLOCATE ( obs_dat(ndata_loc),  mem_dat(ndata_loc,nmem)  )
  ALLOCATE ( obs_trns(ndata_loc), mem_trns(ndata_loc,nmem) )

  if( prdcdf ) then
    if( ndata_loc .lt. ndivid+1 ) then
     write(6,'(a)') 'SAMPLES(ndata_loc) should be larger than NDIVID(PARAM)'
     write(6,'(a,i,a,i)') '  SAMPLES(ndata_loc) ==',ndata_loc,'  NDIVID(PARAM) ==',ndivid      
     stop
    end if 
    ALLOCATE ( obs_cdf(gall,0:ndivid),  gus_cdf(gall,0:ndivid)  ,mem_cdf(gall,0:ndivid)  )
    ALLOCATE ( obs_cdfr(gall,0:ndivid), gus_cdfr(gall,0:ndivid) ,mem_cdfr(gall,0:ndivid) )
  end if
  if(prdpdf) ALLOCATE ( obs_pdf(gall,binnum), gus_pdf(gall,binnum), mem_pdf(gall,binnum) )

  ALLOCATE ( obs_norain(gall), gus_norain(gall), mem_norain(gall) )  
  ALLOCATE ( obs_smp(gall),    gus_smp(gall),    mem_smp(gall)    )
  ALLOCATE ( obs_smpr(gall),   gus_smpr(gall),   mem_smpr(gall)   )
  ALLOCATE ( obs_mask(gall),   gus_mask(gall),   mem_mask(gall)   )
  ALLOCATE ( obs_maskr(gall),  gus_maskr(gall),  mem_maskr(gall)  )
  ALLOCATE ( obs_zero(gall),   gus_zero(gall),   mem_zero(gall)   )
  ALLOCATE ( corr_n(gall),     corr_g(gall) )

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(END RECORD READ):',rtimer-rtimer00,' // NPROCS ',nprocs

  !=======================> get PDF and CDF <=======================!
  do iread=1,nread
  l = proc_l(myrank,iread)

  obs_mask(:)  = 0.0d0
  gus_mask(:)  = 0.0d0
  mem_mask(:)  = 0.0d0
  obs_maskr(:) = 0.0d0
  gus_maskr(:) = 0.0d0
  mem_maskr(:) = 0.0d0
  obs_zero(:)  = rain_min_rate
  gus_zero(:)  = rain_min_rate
  mem_zero(:)  = rain_min_rate

  do g=1,gall
!    if( myrank==0 ) rtimer=mpi_wtime()
!    if( myrank==0 ) write(6,'(A,F13.2,i)') '### TIMER(STS  GATHER):',rtimer-rtimer00,g

    ! get data in the localization circle
    obs_org(:)   = -999.d0
    gus_org(:)   = -999.d0
    mem_org(:)   = -999.d0
    mem_dat(:,:) = -999.d0
    idata     = 0
    idata_mem = 0
    do igrd=1,ngrd
      gg = vec_grd(g,l,igrd)
      ll = vec_rgn(g,l,igrd)
      if( gg==0 .and. ll==0 ) EXIT
      do idate=1,ndate
        idata = idata + 1
        obs_org(idata) = obs_rain(gg,ll,idate) 
        gus_org(idata) = gus_rain(gg,ll,idate)        
        obs_dat(idata) = obs_org(idata)
        if( obs_dat(idata)<-0.00001 ) obs_dat(idata) = undef
        do ibv=1,nmem
          idata_mem = idata_mem + 1
          mem_org(idata_mem) = mem_rain(gg,ll,ibv,idate)
          mem_dat(idata,ibv) = mem_rain(gg,ll,ibv,idate)
          if( mem_dat(idata,ibv)<-0.00001 ) mem_dat(idata,ibv) = undef
        end do
      end do
    end do
                      
!    if( myrank==0 ) rtimer=mpi_wtime()
!    if( myrank==0 ) write(6,'(A,F13.2,i)') '### TIMER(END  GATHER):',rtimer-rtimer00,g
    
    !===> compute CDF & PDF
    call get_cdf(ndata_loc,obs_org,obs_cdfr(g,0:ndivid),obs_smpr(g),ndivid,obs_maskr(g),rain_min_rate) ! CDF without zero precip.
    call get_cdf(ndata_loc,gus_org,gus_cdfr(g,0:ndivid),gus_smpr(g),ndivid,gus_maskr(g),rain_min_rate) ! CDF without zero precip.
    call get_cdf(ndata_mem,mem_org,mem_cdfr(g,0:ndivid),mem_smpr(g),ndivid,mem_maskr(g),rain_min_rate) ! CDF without zero precip.

!    if( myrank ==1 ) &        
    call get_cdf(ndata_loc,obs_org,obs_cdf(g,0:ndivid) ,obs_smp(g), ndivid,obs_mask(g),0.0d0        )  ! CDF with    zero precip.
    call get_cdf(ndata_loc,gus_org,gus_cdf(g,0:ndivid) ,gus_smp(g), ndivid,gus_mask(g),0.0d0        )  ! CDF with    zero precip.
    call get_cdf(ndata_mem,mem_org,mem_cdf(g,0:ndivid) ,mem_smp(g), ndivid,mem_mask(g),0.0d0        )  ! CDF with    zero precip.
    
!    if( myrank ==1 ) &
    call get_norain(ndivid,obs_cdf(g,0:ndivid),obs_norain(g),rain_min_rate,obs_mask(g),minnorain)      ! get No-Rain 
    call get_norain(ndivid,gus_cdf(g,0:ndivid),gus_norain(g),rain_min_rate,gus_mask(g),minnorain)      ! get No-Rain 
    call get_norain(ndivid,mem_cdf(g,0:ndivid),mem_norain(g),rain_min_rate,mem_mask(g),minnorain)      ! get No-Rain

    if( opt_fixzero==1 .and. obs_mask(g)>0.5 ) then
      ! member
      if( mem_mask(g)>0.5 ) then ! member
        tmp_min_rate = rain_min_rate
        if( mem_norain(g) < obs_norain(g) ) then ! nicam has larger precipitation
          do
            if( mem_norain(g) > (obs_norain(g)-alw_fixzero) ) EXIT
            tmp_min_rate = tmp_min_rate + stp_fixzero
            call get_norain(ndivid,mem_cdf(g,0:ndivid),mem_norain(g),tmp_min_rate,mem_mask(g),minnorain)      ! get No-Rain           
          end do
        else if( obs_norain(g) < mem_norain(g) ) then ! nicam has smaller preciptation          
          do
            if( obs_norain(g) > (mem_norain(g)-alw_fixzero) ) EXIT            
            tmp_min_rate = tmp_min_rate - stp_fixzero
            if( tmp_min_rate .le. 0.0d0 ) EXIT
            call get_norain(ndivid,mem_cdf(g,0:ndivid),mem_norain(g),tmp_min_rate,mem_mask(g),minnorain)      ! get No-Rain            
          end do
        end if
        mem_zero(g) = tmp_min_rate
      end if ! mem_mask
    endif

    if ( prdpdf ) then ! PDF production
      call get_pdf(binsiz,binnum,rmin,rmax,ndata_loc,obs_org,obs_pdf(g,1:binnum),rain_min_rate ) ! PDF without zero precip.
      call get_pdf(binsiz,binnum,rmin,rmax,ndata_loc,gus_org,gus_pdf(g,1:binnum),rain_min_rate ) ! PDF without zero precip.
      call get_pdf(binsiz,binnum,rmin,rmax,ndata_mem,mem_org,mem_pdf(g,1:binnum),rain_min_rate ) ! PDF without zero precip.
    end if

!    if( myrank==0 ) rtimer=mpi_wtime()
!    if( myrank==0 ) write(6,'(A,F13.2,i)') '### TIMER(END PDF/CDF):',rtimer-rtimer00,g

    !===> compute gaussian variable for quality control
    if( opt_qccorr .eq. 1 ) then
      do idata = 1,ndata_loc
        obs_trns(idata)   = obs_dat(idata)
        mem_trns(idata,:) = mem_dat(idata,:)

        if( mem_mask(g) > 0.5 ) then
          do ibv=1,nmem
            mem_trns(idata,ibv) = pptrans_normal( mem_dat(idata,ibv), mem_cdf(g,:), mem_norain(g), mem_zero(g)  )
          end do
        else
          mem_trns(idata,:) = undef
        endif

        if( (obs_mask(g)>0.5) .and. obs_dat(idata)>-0.000001 ) then
          obs_trns(idata) = pptrans_normal( obs_dat(idata), obs_cdf(g,:), obs_norain(g), obs_zero(g)  )
        else
          obs_trns(idata) = undef
        endif
      end do

      call calc_corr( ndata_loc, nmem, undef, obs_dat,  mem_dat,  corr_n(g) ) ! normal
      call calc_corr( ndata_loc, nmem, undef, obs_trns, mem_trns, corr_g(g) ) ! gaussina transform       
!      write(6,'(2i,2f)') g,l,corr_n(g),corr_g(g) 
    else
      corr_n(g) = 1.0d0
      corr_g(g) = 1.0d0
    endif
    if( obs_mask(g)<0.5 .or. mem_mask(g)<0.5) corr_n(g) = undef
    if( obs_mask(g)<0.5 .or. mem_mask(g)<0.5) corr_g(g) = undef 

  end do ! g

!  write(6,'(a,i,2f)') 'Region + Max & Min Value of Corr',l,maxval(corr_n(:)),minval(corr_n(:))

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(FWRITE  CDFMAKE):',rtimer-rtimer00,' // NPROCS ',nprocs

  !===> Output Generation
  call MISC_make_idstr(oname,trim(obs_cdfname),'rgn',l)
  call MISC_make_idstr(gname,trim(gus_cdfname),'rgn',l)
  call MISC_make_idstr(mname,trim(mem_cdfname),'rgn',l)

!  if(myrank==1)  write(6,'(i,a,x,a)') l,trim(oname),trim(gname)
  open(1,file=trim(oname),form='unformatted',access='direct',recl=gall*8)
  open(2,file=trim(gname),form='unformatted',access='direct',recl=gall*8)
  open(3,file=trim(mname),form='unformatted',access='direct',recl=gall*8)
    do idiv=0,ndivid
      write(1,rec=idiv+1)          obs_cdf(:,idiv)
      write(2,rec=idiv+1)          gus_cdf(:,idiv)
      write(3,rec=idiv+1)          mem_cdf(:,idiv)
      write(1,rec=ndivid+1+idiv+1) obs_cdfr(:,idiv)
      write(2,rec=ndivid+1+idiv+1) gus_cdfr(:,idiv)
      write(3,rec=ndivid+1+idiv+1) mem_cdfr(:,idiv)
    end do
    write(1,rec=2*(ndivid+1)+1) dble( obs_smp(:) )
    write(2,rec=2*(ndivid+1)+1) dble( gus_smp(:) )
    write(3,rec=2*(ndivid+1)+1) dble( mem_smp(:) )      
    write(1,rec=2*(ndivid+1)+2)  obs_norain(:)
    write(2,rec=2*(ndivid+1)+2)  gus_norain(:)
    write(3,rec=2*(ndivid+1)+2)  mem_norain(:)
    write(1,rec=2*(ndivid+1)+3)  obs_mask(:)
    write(2,rec=2*(ndivid+1)+3)  gus_mask(:)
    write(3,rec=2*(ndivid+1)+3)  mem_mask(:)
    write(1,rec=2*(ndivid+1)+4)  corr_g(:) ! write same correlation
    write(2,rec=2*(ndivid+1)+4)  corr_g(:) ! write same correlation
    write(3,rec=2*(ndivid+1)+4)  corr_g(:) ! write same correlation
    write(1,rec=2*(ndivid+1)+5)  obs_zero(:) 
    write(2,rec=2*(ndivid+1)+5)  gus_zero(:) 
    write(3,rec=2*(ndivid+1)+5)  mem_zero(:)     
  close(1)
  close(2)
  close(3)

  !===> output Generation PDF & CDF
  if ( svpcdf ) then
    if ( prdpdf ) then ! =====> PDF PRODUCTION
      do g=1,1
        write(fname,'(a,a,i5.5,a,i5.5,a,i10.10,a)') trim(snpdir),'/PDF_g',g,'_r',l-1,'_',gcdate,'.txt'
        open(1,file=trim(fname),form='formatted')
          do ibin=1,binnum
            write(1,'(i,f9.4,i,f,i,f,i,f)') ibin,0.5d0*(rmin(ibin)+rmax(ibin)), &
               obs_smp(g), obs_pdf(g,ibin), &
               gus_smp(g), gus_pdf(g,ibin), &
               mem_smp(g), mem_pdf(g,ibin)
          end do
        close(1)
      end do
    end if

    if ( prdcdf ) then ! =====> CDF PRODUCTION
      do g=1,1
        write(fname,'(a,a,i5.5,a,i5.5,a,i10.10,a)') trim(snpdir),'/CDF_g',g,'_r',l-1,'_',gcdate,'.txt'
        open(1,file=trim(fname),form='formatted')
          do ibin=1,ndivid
            write(1,'(i,f7.3,i,3f15.7,i,3f15.7,i,4f15.7)') &
               ibin,sngl( dble(ibin-1) / dble(ndivid-1) ), &
               obs_smp(g), obs_cdf(g,ibin), obs_cdfr(g,ibin), obs_norain(g), &
               gus_smp(g), gus_cdf(g,ibin), gus_cdfr(g,ibin), gus_norain(g), &
               mem_smp(g), mem_cdf(g,ibin), mem_cdfr(g,ibin), mem_norain(g), &
               rain_min_rate
          end do
        close(1)
      end do
    end if          
  
  end if ! svpcdf
  end do ! iread to parallelize "l"

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(END     CDFMAKE):',rtimer-rtimer00,' // NPROCS ',nprocs
  call MPI_Finalize(ierr)

end program

!=========================================================
subroutine calc_corr( ndata, nmem, undef, obs, mem, corr )
!=========================================================
  implicit none
  integer, intent(in) :: ndata, nmem
  double precision, intent(in)  :: obs(ndata), mem(ndata,nmem), undef
  double precision, intent(out) :: corr

  integer idata,ibv,icon
  double precision ave_obs, ave_mem, var_obs, var_mem
  double precision, parameter :: minrain = -10.0d0

  icon   = 0
  ave_obs= 0.0d0
  ave_mem= 0.0d0
  do idata=1,ndata
    do ibv=1,nmem
      if( obs(idata).ge.minrain .and. mem(idata,ibv).ge.minrain ) then
        icon    = icon    + 1
        ave_obs = ave_obs + obs(idata)
        ave_mem = ave_mem + mem(idata,ibv)
      endif
    end do
  end do
  ave_obs = ave_obs / dble( icon )
  ave_mem = ave_mem / dble( icon )

  if( icon==0 .or. (ave_obs.lt.minrain) .or. (ave_obs.lt.minrain) )  then
    corr = undef
  else
    corr    = 0.0d0
    var_obs = 0.0d0
    var_mem = 0.0d0
    do idata=1,ndata
      do ibv=1,nmem
        if( obs(idata).ge.minrain .and. mem(idata,ibv).ge.minrain ) then
          corr    = corr    + (obs(idata)     - ave_obs) * (mem(idata,ibv) - ave_mem)
          var_obs = var_obs + (obs(idata)     - ave_obs)**2.0d0
          var_mem = var_mem + (mem(idata,ibv) - ave_mem)**2.0d0
        endif
      end do
    end do
    corr = corr / dsqrt(var_obs) /dsqrt(var_mem)
  endif

end
!=========================================================
subroutine get_norain(ndivid,cdf,norain,rain_min_rate,mask,minnorain)
!=========================================================
  implicit none
  integer ndivid, idate, idiv
  double precision cdf(0:ndivid), norain, rain_min_rate, mask, minnorain

  norain = 1.0d0
  do idiv = 0,ndivid
    if( cdf(idiv) .ge. rain_min_rate ) then
      norain = dble( idiv ) / dble (ndivid)      
      EXIT
    end if
  end do

  if( norain > minnorain ) mask = 0.0d0

end
!=========================================================
subroutine get_cdf(ndate,org,cdf,nsamp,ndivid,mask,rain_min)
!=========================================================
  implicit none
  integer ndate, nsamp, ndivid, idate, isamp, idiv, ibin
  integer init(ndate)
  double precision mask
  double precision rain_min,org(ndate),srt(ndate),cdf(0:ndivid)

  mask   = 0.0d0
  cdf(:) = 0.0d0
  nsamp  = 0
   
  do idate=1,ndate    
    init(idate) = idate
    srt(idate)  = org(idate)
    if ( srt(idate) .lt. rain_min ) srt(idate) = -999.d0
  end do
      
  call quick_sort(srt,init,1,ndate)
  call get_nsamp(ndate,srt,nsamp) 
  
  if( nsamp .lt. ndivid ) then 
    cdf(:) = 0.0d0 ! no cdf production
  else
    mask = 1.0d0
    do idiv=0,ndivid
      ibin = ( ndate - nsamp ) + nint( dble(nsamp*(idiv+1)/(ndivid+1)) )
      if( idiv==0 ) ibin = ibin + 1
      cdf(idiv) = srt(ibin)    
    end do               
  end if

end
!=========================================================
subroutine get_nsamp(ndate,var,nsamp)
!=========================================================
  implicit none
  integer ndate,nsamp,idate
  double precision var(ndate)
  
  nsamp = 0
  do idate=1,ndate
    if( var(idate).ge.0.0d0 ) nsamp = nsamp + 1
  end do

end
!=========================================================
subroutine get_pdf(binsiz,binnum,rmin,rmax,ndate,org,pdf,rain_min )
!=========================================================
  implicit none
  integer binnum,ndate,idate,ibin,icon
  double precision binsiz, rain_min
  double precision rmin(binnum),rmax(binnum),pdf(binnum)
  double precision org(ndate)
  
  icon   = 0
  pdf(:) = 0.0d0
  do idate=1,ndate
    if ( org(idate) .ge. rain_min ) then ! minimum threshold
    do ibin=1,binnum  
      if( rmin(ibin).lt.org(idate) .and. org(idate).le.rmax(ibin) ) then
        icon = icon + 1
        pdf(ibin) = pdf(ibin) + 1.0d0 
        goto 100
      end if
    end do ! ibin
100 continue
    end if        
  end do ! idate     
  
  do ibin = 1,binnum
    pdf(ibin) = pdf(ibin) / dble(icon)
    if( icon .le. 10 ) pdf(ibin) = 0.0d0
  end do

end
!=========================================================
subroutine make_vector(gall,lall,locrad,myrank,nprocs,nread,proc_l)
!=========================================================
! production of reference table following "obs_local" in letkf_tools.f90 
  use mpi
  implicit none
  !==> parameters
  REAL(8),PARAMETER :: pi=3.1415926535d0
  REAL(8),PARAMETER :: re=6371.3d3
  REAL(8),PARAMETER :: r180=1.0d0/180.0d0

  !==> inner variables
  integer gall,lall, g, l, gg, ll, num_grd, np_rgn, np_grd, sp_rgn, sp_grd
  integer icount,icount_max,icount_min,irec
  REAL(8) locrad  
  REAL(8) ico_lat(gall,lall), ico_lon(gall,lall)
  REAL(8) dist_zero,dlat_zero,dlon_zero(gall,lall)
  
  REAL(8) minlon,maxlon,minlat,maxlat,xlon,xlat,dist
  REAL(8) lon1,lon2,lat1,lat2,cosd
  

  !===> paralell computing
  INTEGER nprocs, nread, iread, ierr, myrank, iprocs
  INTEGER proc_l(0:nprocs-1,nread)
  INTEGER nmax,cnt_pal(gall),cnt_tmp(gall,nprocs),cnt(gall,lall)

  INTEGER, ALLOCATABLE :: vec_grd(:,:,:), vec_grd_pal(:,:), vec_grd_tmp(:,:,:)
  INTEGER, ALLOCATABLE :: vec_rgn(:,:,:), vec_rgn_pal(:,:), vec_rgn_tmp(:,:,:)
  REAL(4), ALLOCATABLE :: msk(:,:)      , msk_pal(:),       msk_tmp(:,:)

  !=====> start
  num_grd = int( sqrt(real(gall)) )    
  read(1) ico_lon  ! [deg]
  read(1) ico_lat  ! [deg]
  
!  dist_zero = locrad * SQRT(10.0d0/3.0d0) * 2.0d0   ! [m]
  dist_zero = locrad  ! [m] kotsuki 20150608
  dlat_zero = dist_zero / pi / re * 180.0d0         ! [deg]
  do g=1,gall
  do l=1,lall
    dlon_zero(g,l) = dlat_zero / DCOS(pi*ico_lat(g,l)/180.0d0)       ! [deg]
    if( ico_lon(g,l) < 0.0d0 ) ico_lon(g,l) = ico_lon(g,l) + 360.0d0 ! [0:360]  
  end do
  end do 
  call set_pole(gall, lall, num_grd, ico_lat, np_rgn, np_grd, sp_rgn, sp_grd)

  ! approxmation of grids in the circle
  call calc_dist_deg(ico_lon(2,1),ico_lon(3,1),ico_lat(2,1),ico_lat(3,1),dist) ! approx. single distance
  nmax = ( int( dist_zero / dist )**2 ) * 10 

  ALLOCATE ( vec_grd(gall,lall,nmax), vec_grd_pal(gall,nmax), vec_grd_tmp(gall,nmax,nprocs) )
  ALLOCATE ( vec_rgn(gall,lall,nmax), vec_rgn_pal(gall,nmax), vec_rgn_tmp(gall,nmax,nprocs) )
  ALLOCATE ( msk(gall,lall)         , msk_pal(gall),          msk_tmp(gall,nprocs)          )

  do iread=1,nread
  l = proc_l(myrank,iread)
  cnt_pal(:)       = 0
  vec_grd_pal(:,:) = 0
  vec_rgn_pal(:,:) = 0
  msk_pal(:)       = 1.0

  do g=1,gall
    minlon = ico_lon(g,l) - dlon_zero(g,l)
    maxlon = ico_lon(g,l) + dlon_zero(g,l)
    minlat = ico_lat(g,l) - dlat_zero
    maxlat = ico_lat(g,l) + dlat_zero
    IF(maxlon - minlon >= 360.0d0) THEN
      minlon = 0.0d0
      maxlon = 360.0d0
    END IF

    if(g >= 1               .and. g <= num_grd        ) msk_pal(g) = 0.0 
    if(g >= gall-num_grd+1  .and. g <= gall           ) msk_pal(g) = 0.0 
    if(mod(g, num_grd) == 0 .or. mod(g, num_grd) == 1 ) msk_pal(g) = 0.0 
    
    icount = 0
    do gg=1,gall
    do ll=1,lall
      xlat = ico_lat(gg,ll)      
      xlon = ico_lon(gg,ll)
      
      if(gg >= 1 .and. gg <= num_grd) cycle                       ! edge exception 
      if(gg >= gall-num_grd+1 .and. gg <= gall) cycle             ! edge exception
      if(mod(gg, num_grd) == 0 .or. mod(gg, num_grd) == 1) cycle  ! edge exception
            
      IF( minlat.le.xlat .and. xlat.le.maxlat ) THEN ! LAT EXCEPTION        
        call calc_dist_deg(xlon,ico_lon(g,l),xlat,ico_lat(g,l),dist)  ! sphere     
!        call calc_length_deg(xlon,ico_lon(g,l),xlat,ico_lat(g,l),dist) ! ellipse
!        write(6,'(7f)') xlon,xlat,xlon-ico_lon(g,l),xlat-ico_lat(g,l),dist,dist_zero
        
        if ( dist < dist_zero ) then
          icount = icount + 1
          if( icount.ge.nmax ) write(6,*) 'ALARM',icount,nmax 
          vec_grd_pal(g,icount) = gg
          vec_rgn_pal(g,icount) = ll
        end if
      END IF
    end do ! ll
    end do ! gg
    
    ! NORTH POLE OR SOUTH POLE
    do gg=1,gall
    do ll=1,lall
      if( (gg==np_grd .and. ll==np_rgn) .or. (gg==sp_grd .and. ll==sp_rgn) ) then   
        xlat = ico_lat(gg,ll)      
        xlon = ico_lon(gg,ll)
        call calc_dist_deg(xlon,ico_lon(g,l),xlat,ico_lat(g,l),dist) ! sphere
!        call calc_length_deg(xlon,ico_lon(g,l),xlat,ico_lat(g,l),dist) ! ellipse       
        if ( dist < dist_zero ) then
          icount = icount + 1
          if( icount.ge.nmax ) write(6,*) 'ALARM',icount,nmax 
          vec_grd_pal(g,icount) = gg
          vec_rgn_pal(g,icount) = ll
!          write(6,'(a,6f)') 'NP POLE',xlon,xlat,ico_lon(g,l),ico_lat(g,l),dist,dist_zero
        end if
      end if
    end do
    end do

    cnt_pal(g) = icount
  end do ! g  
    
    CALL MPI_ALLGATHER(cnt_pal,     gall,      MPI_INTEGER,  &
                       cnt_tmp,     gall,      MPI_INTEGER, MPI_COMM_WORLD, ierr)    
    CALL MPI_ALLGATHER(vec_grd_pal, gall*nmax, MPI_INTEGER,  &
                       vec_grd_tmp, gall*nmax, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(vec_rgn_pal, gall*nmax, MPI_INTEGER,  &
                       vec_rgn_tmp, gall*nmax, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(msk_pal,     gall,      MPI_REAL   ,  &
                       msk_tmp,     gall,      MPI_REAL   , MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
      do iprocs=0,nprocs-1
        l = proc_l(iprocs,iread)
!        if( myrank==0 ) write(6,*) iprocs, " ===> ", l
        do icount=1,nmax
          cnt(:,l)            = cnt_tmp(:,iprocs+1)
          vec_rgn(:,l,icount) = vec_rgn_tmp(:,icount,iprocs+1)
          vec_grd(:,l,icount) = vec_grd_tmp(:,icount,iprocs+1)
          msk(:,l)            = msk_tmp(:,iprocs+1)
        end do
      enddo ! nprocs
  end do ! iread

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if( myrank==0 ) then
    icount_max     = 0
    icount_min     = 1000000000    
    do l=1,lall
    do g=1,gall
      if ( cnt(g,l) .ge. icount_max ) icount_max = cnt(g,l)
      if ( cnt(g,l) .le. icount_min ) icount_min = cnt(g,l)
    end do
    end do
    write(6,'(a,f,2i)') 'LOC & COUNTER MAX/MIN ===>',locrad,icount_max,icount_min
    cnt(:,:) = icount_max
    write(2,rec=1) real( cnt(:,:) )+0.1
    write(3,rec=1) real( cnt(:,:) )+0.1

    irec= 1
    do icount=1,icount_max
      irec = irec + 1
      write(2,rec=irec) vec_grd(:,:,icount)+0.1
      write(3,rec=irec) vec_rgn(:,:,icount)+0.1
    end do
    write(4,rec=1) msk(:,:)
  endif

  DEALLOCATE ( vec_grd, vec_grd_pal, vec_grd_tmp )
  DEALLOCATE ( vec_rgn, vec_rgn_pal, vec_rgn_tmp )  
  DEALLOCATE ( msk,     msk_pal,     msk_tmp     )
end        
!=========================================================
subroutine calc_dist_deg(xlon1,xlon2,xlat1,xlat2,dist)
!=========================================================
  implicit none
  double precision xlon1,xlon2,xlat1,xlat2
  double precision lon1,lon2,lat1,lat2,cosd,dist

  REAL(8),PARAMETER :: pi=3.1415926535d0
  REAL(8),PARAMETER :: re=6371.3d3
  REAL(8),PARAMETER :: r180=1.0d0/180.0d0
  
  lon1 = xlon1 * pi * r180 ! [rad]
  lon2 = xlon2 * pi * r180 ! [rad]
  lat1 = xlat1 * pi * r180 ! [rad]
  lat2 = xlat2 * pi * r180 ! [rad]

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)
  dist = ACOS( cosd ) * re 

end
!--------------------------------------------------------------------!
  subroutine calc_length_deg(xlon1, xlon2, xlat1, xlat2, length) ! modified on Dr. Kondo's program
    implicit none
    real(8), intent(in) :: xlon1, xlon2, xlat1, xlat2
    real(8), intent(out) :: length
    real(8) :: ilat, ilon, jlat, jlon
    real(8) :: A, B
    real(8), parameter :: Rpo = 6334834.0d0, Req = 6377937.0d0
    real(8), parameter :: e = 0.006674d0
    real(8), parameter :: pi=3.1415926535d0
    real(8), parameter :: r180=1.0d0/180.0d0
    real(8) :: P, dP, dR

    ilat = xlat1 * pi * r180 ! [rad]
    ilon = xlon1 * pi * r180 ! [rad]
    jlat = xlat2 * pi * r180 ! [rad]
    jlon = xlon2 * pi * r180 ! [rad]

    P  = (ilat + jlat)*0.5d0
    dP = dabs(ilat - jlat)
    dR = dabs(ilon - jlon)
    if( dR > pi ) then
      dR = 2.d0*pi - dR
    end if

    A = Rpo/(1.d0 - e*dsin(P)**2)**(1.5d0)
    B = Req/(1.d0 - e*dsin(P)**2)**(0.5d0)

    length = dsqrt( (A*dP)**2 + (B*dcos(P)*dR)**2)
  end 
!=========================================================
subroutine set_pole(gall, lall, num_grd, ico_lat, np_rgn, np_grd, sp_rgn, sp_grd) ! modifed on code by Dr. Kondo
!=========================================================
  implicit none

  integer, intent(in) :: gall, lall, num_grd
  real(8), intent(in) :: ico_lat(gall,lall)
  integer, intent(out) :: np_rgn, np_grd, sp_rgn, sp_grd
  integer :: igrd, irgn
  real(8) :: pi

!  pi = 4.d0*datan(1.0d0) ! [rad] 
  pi = 180.0d0 ! [deg] kotsuki 20150422

  np_grd = gall - num_grd + 2
  sp_grd = num_grd * 2

  do irgn = 1, lall
    if( dabs(ico_lat(np_grd,irgn)*2.d0 - pi) .lt. 0.000001 ) then
      np_rgn = irgn
!      write(6,*) 'NP',ico_lat(np_grd,irgn)
      exit
    end if
  end do
  do irgn = 1, lall
    if( dabs(ico_lat(sp_grd,irgn)*2.d0 + pi) .lt. 0.000001 ) then
      sp_rgn = irgn
!      write(6,*) 'SP',ico_lat(sp_grd,irgn)
      exit
    end if
  end do
  
end
!=========================================================
subroutine check_vector(gall,lall,ngrd,vec_grd,vec_rgn)
!=========================================================
  implicit none
  integer gall, lall, ngrd, g, l, igrd, gg, ll, ii
  integer vec_grd(gall,lall,ngrd), vec_rgn(gall,lall,ngrd)
  REAL(8) dif_lon,dif_lat
  REAL(8) ico_lat(gall,lall), ico_lon(gall,lall)
  REAL(8) rem_lat(ngrd), rem_lon(ngrd)

  read(1) ico_lon  ! [deg]
  read(1) ico_lat  ! [deg]

  do g=1,gall
  do l=1,lall
    rem_lat(:) = 0.0d0
    rem_lon(:) = 0.0d0
    do igrd=1,ngrd
      gg = vec_grd(g,l,igrd)
      ll = vec_rgn(g,l,igrd)
      if( gg==0 .and. ll==0 ) exit 
      rem_lon(igrd) = ico_lon(gg,ll)
      rem_lat(igrd) = ico_lat(gg,ll)
      if ( igrd .ne. 1 ) then
        do ii=1,igrd-1
          dif_lon = dabs( rem_lon(ii) - ico_lon(gg,ll) )
          dif_lat = dabs( rem_lat(ii) - ico_lat(gg,ll) )
          if( dif_lon .ge. 180.0d0 ) dif_lon = dabs( dif_lon - 360.0d0 )
!          write(6,'(4i,2f)') gg,ll,igrd,ii,dif_lon,dif_lat
          if( dif_lon<0.000001 .and. dif_lat<0.0000001) then
            write(6,'(a,4f)') ' CAUTION===> SAME LAT LON',ico_lon(gg,ll),ico_lat(gg,ll),rem_lon(ii),rem_lat(ii)
          end if
        end do
      end if
    end do
  end do
  end do
end 
!=========================================================
subroutine check_usage(gall,lall,usage)
!=========================================================
  implicit none
  integer gall, lall, ngrd, g, l, gg, ll
  REAL(4), intent(in) :: usage(gall,lall)
  REAL(8) ico_lat(gall,lall), ico_lon(gall,lall)
  REAL(8) tmp_lon, tmp_lat, dif_lon, dif_lat

  read(1) ico_lon  ! [deg]
  read(1) ico_lat  ! [deg]

  do l=1,lall
  do g=1,gall
    tmp_lon = ico_lon(g,l)
    tmp_lat = ico_lat(g,l)
    if( usage(g,l) > 0.5 ) then
      do ll=1,lall
        if( l == ll ) cycle
        do gg=1,gall
          if( g == gg ) cycle
          if( usage(gg,ll) < 0.5 ) cycle
          dif_lon = dabs( tmp_lon - ico_lon(gg,ll) )
          dif_lat = dabs( tmp_lat - ico_lat(gg,ll) )        

          if( dif_lon<0.000001 .and. dif_lat<0.0000001) &
            print *, " [Err] find same lon lat grid",tmp_lon,tmp_lat
        end do
      end do
    else
      do ll=1,lall
        if( l == ll ) cycle
        do gg=1,gall
          if( g == gg ) cycle
          dif_lon = dabs( tmp_lon - ico_lon(gg,ll) )
          dif_lat = dabs( tmp_lat - ico_lat(gg,ll) )        

          if( dif_lon<0.000001 .and. dif_lat<0.0000001) goto 100
        end do
      end do
      print *, " [Err ???] cannot find same grid",tmp_lon,tmp_lat
100   continue
    end if
  end do
  end do
  
  print *, "fin check usage"
end
!=========================================================
subroutine smfilter( gall, inp )
!=========================================================
  integer, intent(in)    :: gall
  real(4), intent(inout) :: inp(gall)

  integer i,j, ii, jj, icon
  real(4), ALLOCATABLE :: tmp(:,:), flt(:,:)

  mx     =  int( sqrt( real(gall) ) )
  ALLOCATE ( tmp(mx,mx), flt(mx,mx) ) 

  do j=1,mx
  do i=1,mx
    tmp(i,j) = inp( (j-1)*mx + 1 )
  end do
  end do

  flt(:,:) = tmp(:,:)
  do j=2,mx-1
  do i=2,mx-1
    icon = 0
    flt(i,j) = 0.0d0
    do ii=i-1,i+1
    do jj=j-1,j+1
      if( tmp(ii,jj) .ge. -0.0 ) icon     = icon + 1
      if( tmp(ii,jj) .ge. -0.0 ) flt(i,j) = flt(i,j) + tmp(ii,jj)
    enddo
    end do

    if ( icon .ge. 0.5 ) then
      flt(i,j) = flt(i,j) / real(icon)
    else    
      flt(i,j) = tmp(i,j)
    endif
  end do
  end do

  do j=1,mx
  do i=1,mx
    inp( (j-1)*mx + 1 ) = flt(i,j)
  end do
  end do

  DEALLOCATE ( tmp, flt ) 
end
!=========================================================
    recursive subroutine quick_sort(var,init,first,last)
!=========================================================
    integer first,last
    integer init(*)
    double precision var(*)
    double precision x,t
    
    x = var( (first+last) / 2 )
    i = first
    j = last

    do
      do while (var(i) < x)
        i=i+1
      end do
      do while (x < var(j))
        j=j-1
      end do
      if (i >= j) exit
        t       = var(i)
        var(i)  = var(j) 
        var(j)  = t
        it      = init(i)
        init(i) = init(j)
        init(j) = it
        
        i       = i + 1
        j       = j - 1
    end do
    if (first < i - 1 ) call quick_sort(var,init, first, i - 1)
    if (j + 1 < last  ) call quick_sort(var,init, j + 1, last)    
    
    return
    end

