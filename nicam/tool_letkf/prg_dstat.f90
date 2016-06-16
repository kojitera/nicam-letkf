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

program dstat

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

  !==> allocatable variables
  real(r_size), allocatable :: ico_lat(:,:), ico_lon(:,:), landmask(:)
  real(4), allocatable :: tmp(:,:), nph(:,:), trp(:,:), sph(:,:), lnd(:,:), srf(:,:), ocn(:,:)

  real(4), allocatable :: lnd_rgn(:), lnd_tmp(:,:), lnd_all(:,:)
  real(4), allocatable :: gus_rgn(:), gus_tmp(:,:), gus_all(:,:,:) 
  real(4), allocatable :: anl_rgn(:), anl_tmp(:,:), anl_all(:,:,:)
  real(4), allocatable :: obs_rgn(:), obs_tmp(:,:), obs_all(:,:,:) 
  real(4), allocatable :: ddd_rgn(:), ddd_tmp(:,:), ddd_all(:,:,:) ! Y-Hx^f
  real(4), allocatable :: zzz_rgn(:), zzz_tmp(:,:), zzz_all(:,:,:) ! Y-Hx^a
  real(4), allocatable :: msk_rgn(:), msk_tmp(:,:), msk_all(:,:,:)
  real(4), allocatable :: qc_rgn(:), ppmask_o(:), ppmask_m(:)
  real(4), allocatable :: ghst(:), ohst(:), dhst(:)
  real(4), allocatable :: gfnc(:), ofnc(:), dfnc(:)
  real(4) :: gcnt, ocnt, dcnt
  real(4) :: gttl, ottl, dttl, gave, oave, dave, gstd, ostd, dstd

  real(r_size), allocatable :: hmin(:), hmax(:)
  real(r_size) :: hmdl

  !==> global variables
  integer, save :: gall, lall
  integer, save :: idate, ndate
  integer, save :: icell, ncell
  integer, save :: itemp, icase, isurf
  character(255) :: inpname, fname, mapname

  integer, parameter :: ncase = 4 ! (glb,nph,trp,sph)
  integer, parameter :: nsurf = 3 ! (all,lnd,ocn)
  character(3) :: ccase, csurf
  
  !==> dstat.cnf
  integer, save  :: psdate, pedate
  integer, save  :: nhour  ! assimilation cycle
  real(r_size)   :: dsiz, asiz, lngmax
  character(255) :: DSTATDIR, LANDDIR


  !==> nhm_driver.cnf
  integer, save :: glevel
  integer, save :: rlevel
  integer, save :: vlayer
  character(128) rgnmngfname

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
  !==> rain2trans.cnf
  call get_rain2trans()

  !==> nhm_driver.cnf
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

  !==> get ico_lon ico_lat
  allocate( ico_lon(gall,lall), ico_lat(gall,lall) )
  write(mapname,'(a,a,i2.2,a,i2.2,a)') trim(vecdir),'/icolatlon_gl',glevel,'_rl',rlevel,'.dat'
  open(1,file=trim(mapname),form='unformatted',access='sequential',action='read')
  read(1) ico_lon  ! [deg]
  read(1) ico_lat  ! [deg]
  close(1)

  !==> dstat.cnf
  namelist / dstat_param / &
    psdate,  pedate,  &
    LANDDIR, DSTATDIR, nhour, &
    dsiz, asiz, lngmax
  open(1,file='dstat.cnf')
  read(1,nml=dstat_param)
  close(1)


  ncell = int ( (gaussmax-gaussmin) / dsiz ) 
  ALLOCATE ( hmin(ncell),  hmax(ncell) ) ! to histgram
  
  do icell = 1,ncell
    hmin(icell) = gaussmin + dsiz * dble( icell - 1 )
    hmax(icell) = gaussmin + dsiz * dble( icell - 0 )
!    if( myrank==0 ) print *, icell, hmin(icell), hmax(icell)
  end do   

  !=====> get ndate
  idate = psdate ! initial time
  ndate = 0
  do
    ndate = ndate + 1           
!    if(ADM_prc_me==1) write(6,'(a,3i)') "     CDFMAKE:: CDATE, NCON, EDATE",idate,ndate,gedate
    if( idate == pedate ) goto 10
    call update_date(idate,nhour)
    if ( ndate .ge. 10000 ) then
      write(6,*) 'ndate is too large'
      stop
    end if
  end do
10 continue

  !===> setting for parallel computing
  if( mod(lall,nprocs).ne.0 ) then
    write(6,*) "lall should be devided by the nprocs"
    stop
  endif
  
  nread = lall / nprocs
  allocate ( proc_l(0:nprocs-1,nread) )
  call set_proc_l(lall,nprocs,nread,proc_l)

  !===> set allocate
  allocate ( qc_rgn(gall),  ppmask_m(gall), ppmask_o(gall)    )
  allocate ( gus_rgn(gall), gus_tmp(gall,nprocs), gus_all(gall,lall,ndate) )
  allocate ( anl_rgn(gall), anl_tmp(gall,nprocs), anl_all(gall,lall,ndate) )
  allocate ( obs_rgn(gall), obs_tmp(gall,nprocs), obs_all(gall,lall,ndate) )
  allocate ( ddd_rgn(gall), ddd_tmp(gall,nprocs), ddd_all(gall,lall,ndate) )
  allocate ( zzz_rgn(gall), zzz_tmp(gall,nprocs), zzz_all(gall,lall,ndate) )
  allocate ( msk_rgn(gall), msk_tmp(gall,nprocs), msk_all(gall,lall,ndate) )  
  allocate ( lnd_rgn(gall), lnd_tmp(gall,nprocs), lnd_all(gall,lall), landmask(gall) )

  idate = psdate ! initial time
  ndate = 0

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(STS RECORD READ):',rtimer-rtimer00,' // NPROCS ',nprocs
  do
    ndate = ndate + 1           
    if(myrank==0) write(6,'(a,3i)') "     Read... Transform:: CDATE, NCON, EDATE",idate,ndate,pedate

    do iread=1,nread
      l = proc_l(myrank,iread)
      write(inpname,'(2a,i10.10,a)') trim(datdir),'/../',idate,'/monit_transform'      
      call MISC_make_idstr(fname,trim(inpname),'rgn',l)

      open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4,action="read")
        read(1,rec=3) gus_rgn(:)
        read(1,rec=4) obs_rgn(:)
        read(1,rec=6) ppmask_m(:)
        read(1,rec=7) ppmask_o(:)
        read(1,rec=8) qc_rgn(:)
      close(1)

      write(inpname,'(2a,i10.10,a)') trim(datdir),'/../',idate,'/monit_inverse'      
      call MISC_make_idstr(fname,trim(inpname),'rgn',l)
      open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4,action="read")
        read(1,rec=2) anl_rgn(:)
      close(1)

      write(inpname,'(2a)') trim(LANDDIR),'/landmask'      
      call MISC_make_idstr(fname,trim(inpname),'rgn',l)
      open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*8,action="read")
        read(1,rec=1) landmask(:)
        lnd_rgn(:) = sngl( landmask(:) )
      close(1)


      ddd_rgn(:) = real ( undef )
      zzz_rgn(:) = real ( undef )
      msk_rgn(:) = 0.
      do g=1,gall
        if( qc_rgn(g) < 0.5 ) then
          obs_rgn(g) = real( undef )
          gus_rgn(g) = real( undef )
          anl_rgn(g) = real( undef ) 
        else
          ddd_rgn(g)      = obs_rgn(g) - gus_rgn(g)
          zzz_rgn(g)      = obs_rgn(g) - anl_rgn(g)
          msk_rgn(g)      = 1.
        end if  
      end do        

      CALL MPI_ALLGATHER(ddd_rgn, gall, MPI_REAL,  &
                         ddd_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(zzz_rgn, gall, MPI_REAL,  &
                         zzz_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(anl_rgn, gall, MPI_REAL,  &
                         anl_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(gus_rgn, gall, MPI_REAL,  &
                         gus_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(obs_rgn, gall, MPI_REAL,  &
                         obs_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(msk_rgn, gall, MPI_REAL,  &
                         msk_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLGATHER(lnd_rgn, gall, MPI_REAL,  &
                         lnd_tmp, gall, MPI_REAL, MPI_COMM_WORLD, ierr)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do iprocs=0,nprocs-1
        l = proc_l(iprocs,iread)
!        if( myrank==0 ) write(6,*) iprocs, " ===> ", l
        gus_all(:,l,ndate) = gus_tmp(:,iprocs+1)
        anl_all(:,l,ndate) = anl_tmp(:,iprocs+1)
        obs_all(:,l,ndate) = obs_tmp(:,iprocs+1)
        ddd_all(:,l,ndate) = ddd_tmp(:,iprocs+1)
        zzz_all(:,l,ndate) = zzz_tmp(:,iprocs+1)
        msk_all(:,l,ndate) = msk_tmp(:,iprocs+1)
        lnd_all(:,l)       = lnd_tmp(:,iprocs+1)
      enddo ! nprocs
    end do ! iread

    if( idate == pedate ) goto 20
    call update_date(idate,nhour)
  end do
20 continue

  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(END RECORD READ):',rtimer-rtimer00,' // NPROCS ',nprocs

  !==========> start D-stat Hollingsworth and Lonnberg (1986)
  allocate( tmp(gall,lall), nph(gall,lall), trp(gall,lall), sph(gall,lall) )
  allocate( srf(gall,lall), lnd(gall,lall), ocn(gall,lall) )
  nph(:,:) = 0.0
  trp(:,:) = 0.0
  sph(:,:) = 0.0
  lnd(:,:) = 0.0
  ocn(:,:) = 0.0
  do l=1,lall
  do g=1,gall
    !===> hemispheres
    if( ico_lat(g,l).ge.20 ) then
      nph(g,l) = 1.0
    else if( ico_lat(g,l).ge.-20 ) then
      trp(g,l) = 1.0
    else if( ico_lat(g,l).ge.-90 ) then
      sph(g,l) = 1.0
    end if

    !===> land covers
    if( lnd_all(g,l) > 0.5 ) then
      lnd(g,l) = 1.0
    else
      ocn(g,l) = 1.0
    end if
  end do
  end do


  if( myrank==0 ) then
  write(fname,'(2a,i10,a,i10,a)') trim(DSTATDIR),'/monit_',psdate,'-',pedate,'.txt'
  open(11,file=trim(fname),action='write')
  end if

  allocate ( ghst(ncell), ohst(ncell), dhst(ncell) )
  allocate ( gfnc(ncell), ofnc(ncell), dfnc(ncell) )
  do isurf=1,nsurf
    if( isurf.eq.1 ) csurf='all'    
    if( isurf.eq.2 ) csurf='lnd'    
    if( isurf.eq.3 ) csurf='ocn' 
    if( isurf.eq.1 ) srf(:,:) = 1.0   
    if( isurf.eq.2 ) srf(:,:) = lnd(:,:)
    if( isurf.eq.3 ) srf(:,:) = ocn(:,:)

  do icase=1,ncase
    if( icase.eq.1 ) ccase='glb'    
    if( icase.eq.2 ) ccase='nph'    
    if( icase.eq.3 ) ccase='trp'  
    if( icase.eq.4 ) ccase='sph'  
    if( icase.eq.1 ) tmp(:,:) = 1.0       ! mask for region
    if( icase.eq.2 ) tmp(:,:) = nph(:,:)  ! mask for region 
    if( icase.eq.3 ) tmp(:,:) = trp(:,:)  ! mask for region
    if( icase.eq.4 ) tmp(:,:) = sph(:,:)  ! mask for region

    do l=1,lall
    do g=1,gall
      if( tmp(g,l) > 0.5 .and. srf(g,l) < 0.5 ) tmp(g,l) = 0.0
    end do
    end do

    !====> get PDF
    if( csurf == "all" ) then
      call get_avestdhst(gall,lall,ndate,ncell,dsiz,hmin,hmax,msk_all,tmp,gus_all, gcnt,gttl,gave,gstd,ghst,gfnc)
      call get_avestdhst(gall,lall,ndate,ncell,dsiz,hmin,hmax,msk_all,tmp,obs_all, ocnt,ottl,oave,ostd,ohst,ofnc)
      call get_avestdhst(gall,lall,ndate,ncell,dsiz,hmin,hmax,msk_all,tmp,ddd_all, dcnt,dttl,dave,dstd,dhst,dfnc)

      write(fname,'(2a,i10,a,i10,5a)') trim(DSTATDIR),'/dhist_',psdate,'-',pedate,'_',ccase,'_',csurf,'.txt'
      open(10,file=trim(fname),action='write')
      do icell=1,ncell
        hmdl = 0.5d0 * ( hmin(icell) + hmax(icell) )
        write(10,'(7f)') hmdl, ghst(icell),gfnc(icell),ohst(icell),ofnc(icell),dhst(icell),dfnc(icell)
      end do
      close(10)

      write(fname,'(2a,i10,a,i10,5a)') trim(DSTATDIR),'/dstat_',psdate,'-',pedate,'_',ccase,'_',csurf,'.txt'
      open(10,file=trim(fname),action='write')
      write(10,'(9f15.5)') gttl,gave,gstd,ottl,oave,ostd,dttl,dave,dstd
      close(10)
    endif

    !===> make Dstat
    if( myrank==0 ) then
      write(fname,'(2a,i10,a,i10,5a)') trim(DSTATDIR),'/desro_',psdate,'-',pedate,'_',ccase,'_',csurf,'.txt'
      open(10,file=trim(fname),action='write')
      call Desroziers(gall,lall,ndate,msk_all,tmp,ico_lon,ico_lat,ddd_all,zzz_all,asiz,lngmax,myrank,nprocs,rtimer00,10,11,ccase,csurf)
!!!    call Hollingsworth(gall,lall,ndate,msk_all,tmp,ico_lon,ico_lat,ddd_all,dttl,asiz,lngmax)
      close(10)

    end if

!!    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!    stop
  end do ! icase
  end do ! isurf
  if( myrank==0 ) close(11)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(END  READ TRANSFORM):',rtimer-rtimer00,' // NPROCS ',nprocs
  call MPI_Finalize(ierr)

end program
!=========================================================
subroutine Desroziers(gall,lall,ndate,msk,msk_rgn,ico_lon,ico_lat,ddd,zzz,asiz,lngmax,myrank,nprocs,rtimer00,fnum,mnum,ccase,csurf)
!=========================================================
!!! Dstat with Desroziers et al. (2005) R=<(y-Hx^a)(y-Hx^f)T>
!!! Dstat with Hollingsworth (1986)     <ddT> = R + HBHT
  use mpi
  implicit none
  integer, intent(in)  :: gall, lall, ndate, myrank, nprocs, fnum, mnum
  real(8), intent(in)  :: asiz, lngmax, ico_lon(gall,lall), ico_lat(gall,lall)
  real(4), intent(in)  :: ddd(gall,lall,ndate), zzz(gall,lall,ndate), msk(gall,lall,ndate), msk_rgn(gall,lall)
  real(4), intent(in)  :: rtimer00
  character(3), intent(in) :: ccase, csurf
  !===> global variables
  integer :: g, l, idate, icon, nmax, nsamp, ilng, nlng
  real(4) :: rtimer, tmpmsk(gall,lall)
  real(4), allocatable :: d(:), z(:)  
  real(8), allocatable :: dlon(:), dlat(:), acnt(:), ades(:), ahol(:), afnc(:)
  real(4), parameter :: undef = -9999.
  real(8) :: length, angle,   RRR, RHBHT

  !===> cnt
  nlng  = int( lngmax / asiz )  
  allocate ( acnt(0:nlng), ades(0:nlng), ahol(0:nlng), afnc(0:nlng) )

  !===> get maxval
  nmax = 0
  do l=1,lall
  do g=1,gall
    tmpmsk(g,l) = maxval( msk(g,l,1:ndate) )
    if( tmpmsk(g,l)>0.5 ) nmax = nmax + 1
  end do
  end do
  allocate ( d(nmax), z(nmax), dlon(nmax), dlat(nmax) )

  acnt(0:nlng) = 0.0d0 
  ades(0:nlng) = 0.0d0
  ahol(0:nlng) = 0.0d0

  do idate=1,ndate
    nsamp = 0
    d(:) = undef
    z(:) = undef
    do l=1,lall
    do g=1,gall
    if ( msk(g,l,idate)>0.5 .and. msk_rgn(g,l)>0.5 ) then
      nsamp = nsamp + 1
      d(nsamp)    = ddd(g,l,idate)
      z(nsamp)    = zzz(g,l,idate)
      dlon(nsamp) = ico_lon(g,l)
      dlat(nsamp) = ico_lat(g,l) 
    end if
    end do
    end do

    do g = 1,nsamp
    do l = 1,nsamp
      call calc_length_deg(dlon(g), dlon(l), dlat(g), dlat(l), length, angle) !! length [m], angle [deg]
      length = length * 0.001d0 !! [m] --> [km]
      ilng = int( undef )
      if ( length.lt.lngmax+asiz*0.5d0 ) then      
        ilng = int( (length-asiz*0.5d0) / asiz ) + 1
        if( length.lt.0.00001 ) ilng = 0
        acnt(ilng) = acnt(ilng) + 1.0
        ades(ilng) = ades(ilng) + d(g)*z(l) !  R=<(y-Hx^a)(y-Hx^f)T>
        ahol(ilng) = ahol(ilng) + d(g)*d(l) !  <ddT> = R + HBHT
      end if
    end do
    end do
    if( myrank==0 ) rtimer=mpi_wtime()
!    if( myrank==0 ) write(6,'(A,i10,F13.2,A,I)') '### TIMER(DSTAT):',idate,    rtimer-rtimer00,' // NPROCS ',nprocs
  end do 

  do ilng=0,nlng
    length = asiz * dble(ilng)
    if( acnt(ilng)>0.5 ) then
      ades(ilng) = ades(ilng) / acnt(ilng) 
      ahol(ilng) = ahol(ilng) / acnt(ilng)   
    end if
!!!    print '(4f)', length, acnt(ilng), ades(ilng), ahol(ilng)
    write(fnum,'(f8.1,i15,4f15.10)') length, int(acnt(ilng)), &
      ades(ilng), ades(ilng)/ades(0), ahol(ilng), ahol(ilng)/ahol(0)
  end do
  
  RRR    = ades(0)
  RHBHT  = ahol(0)

  print '(a,3f15.10,2a)', "     Estimated R+HBHT, R, R(dev)  ", RHBHT, RRR, dsqrt(RRR), " in ", ccase
  write(mnum,'(3f15.10,3a4)') RHBHT, RRR, dsqrt(RRR),"  ",ccase, csurf

  deallocate ( d, z, dlon, dlat, acnt, ades, ahol, afnc )
return
end  
!=========================================================
subroutine Hollingsworth(gall,lall,ndate,msk,msk_rgn,ico_lon,ico_lat,ddd,ttl,asiz,lngmax)
!=========================================================
  implicit none
  integer, intent(in)  :: gall, lall, ndate
  real(8), intent(in)  :: asiz, lngmax, ico_lon(gall,lall), ico_lat(gall,lall)
  real(4), intent(in)  :: ddd(gall,lall,ndate), msk(gall,lall,ndate),msk_rgn(gall,lall), ttl
  !===> global variables
  integer :: i, j, k, g, l, idate, isamp, nsamp, ilng, nlng, icon
  real(4), allocatable :: d(:)
  real(8), allocatable :: dlon(:), dlat(:), acnt(:), asum(:), afnc(:)
  real(8) :: d2ave, length, angle, angle2, mul_d2ave, std

  !=======> least squares method
  real(8) :: a2, a1, a0, xc, sig2, NNN
  real(8) :: AAAA, BBBB, CCCC, DDDD, NNNN
  ! REFF :: http://www.mk-mode.com/octopress/2014/03/04/fortran-least-squares-method/
  ! REFF :: http://nuclear.phys.tohoku.ac.jp/~ykoba/latex2html/gaussian-fitting/
  


  nsamp = int( ttl )
  allocate ( d(nsamp), dlon(nsamp), dlat(nsamp) )
  nlng  = int( lngmax / asiz )  
  allocate ( acnt(nlng), asum(nlng), afnc(nlng) )

  isamp=0
  d2ave=0.0d0

  do idate=1,ndate
  do l=1,lall
  do g=1,gall
    if ( msk(g,l,idate) > 0.5 .and. msk_rgn(g,l) > 0.5 ) then
      isamp       = isamp + 1
      d2ave       = d2ave + ddd(g,l,idate) ** 2.0d0
      d(isamp)    = ddd(g,l,idate)
      dlon(isamp) = ico_lon(g,l)
      dlat(isamp) = ico_lat(g,l)
    end if
  end do ! gall
  end do ! lall
  end do ! ndate

  d2ave     = d2ave / dble( nsamp )
  mul_d2ave = 1.0d0 / d2ave
  print *, ttl, isamp, nsamp, d2ave, nlng

  acnt(0:nlng) = 0.0d0
  asum(0:nlng) = 0.0d0
  std          = 0.0d0 ! compute with assuming ave = 0.0d0
  icon         = 0
  do g = 1,nsamp
  do l = 1,nsamp
  if( g.ne.l ) then
!!!    call calc_angle_deg(dlon(g), dlon(l), dlat(g), dlat(l), angle)
    call calc_length_deg(dlon(g), dlon(l), dlat(g), dlat(l), length, angle) !! length [m], angle [deg]
    length = length * 0.001d0 !! [m] --> [km]
    if ( asiz*0.5d0.le.length .and. length.lt.lngmax+asiz*0.5d0 ) then      
      ilng = int( (length-asiz*0.5d0) / asiz ) + 1
      acnt(ilng) = acnt(ilng) + 1.0
      asum(ilng) = asum(ilng) + d(g)*d(l) * mul_d2ave
      icon = icon + 1 ! samples for fitting
    end if
  end if
  end do
  end do

  do ilng=1,nlng
    length = asiz * dble(ilng)
    if( acnt(ilng)>0.5 ) asum(ilng) = asum(ilng) / acnt(ilng)    
!    print '(4f)', length, acnt(ilng), asum(ilng), dlog( asum(ilng) )
  end do

!!  !=======> least squares method start
!!  AAAA = 0.0d0
!!  BBBB = 0.0d0
!!  CCCC = 0.0d0
!!  DDDD = 0.0d0
!!  NNNN = 0.0d0
!!  do ilng=1,nlng
!!    length = asiz * dble(ilng)
!!    if( acnt(ilng)>0.5 ) then
!!      AAAA = AAAA + ( length**2.0d0 ) * dlog( asum(ilng) )
!!      BBBB = BBBB + ( length**4.0d0 )
!!      CCCC = CCCC + ( length**2.0d0 )
!!      DDDD = DDDD + dlog( asum(ilng) )
!!      NNNN = NNNN + 1.0d0
!!    end if
!!  end do
!!
!!  a0 = ( BBBB*DDDD - AAAA*CCCC ) / ( BBBB*NNNN - CCCC**2.0d0 )   
!!  a1 = 0.0d0
!!  a2 = ( AAAA*NNNN - CCCC*DDDD ) / ( BBBB*NNNN - CCCC**2.0d0 ) 
!!!  do ilng=1,nlng
!!!    length = asiz * dble(ilng)
!!!    afnc(ilng) = a0 + a2 * ( length**2.0d0 )
!!!    print '(4f)',  length, dlog( asum(ilng)) , afnc(ilng)
!!!  end do
!!  xc   = 0.0d0 - 0.5d0 * a1 / a2
!!  NNN  = dexp( a0 - 0.25d0*(a1**2.0d0)/a2 )
!!  sig2 = 0.0d0 - 0.5d0 * a2

  do ilng=1,nlng   
    length = asiz * dble(ilng)
!    afnc(ilng) = NNN * dexp( - 0.5d0*((length-xc)**2.0d0)/sig2 )
    afnc(ilng) = NNN * dexp( - 0.5d0*((length-xc)**2.0d0)*sig2 )
    print '(4f)',  acnt(ilng), length, asum(ilng), afnc(ilng)
  end do

  stop

    
      
  !==> least squares method

  





return
end
!=========================================================
subroutine get_avestdhst(gall,lall,ndate,ncell,dsiz,hmin,hmax,msk,msk_rgn, var,cnt,ttl,ave,std,hst,fnc)
!=========================================================
  implicit none
  integer, intent(in)  :: gall, lall, ndate, ncell
  real(8), intent(in)  :: dsiz, hmin(ncell), hmax(ncell)
  real(4), intent(in)  :: var(gall,lall,ndate), msk(gall,lall,ndate),msk_rgn(gall,lall)
  real(4), intent(out) :: cnt,ttl,ave,std,hst(ncell),fnc(ncell)
  !==> global variables
  integer :: g, l, idate, icell
  real(8) :: hmdl
  real(8), parameter :: pi = dacos(-1.0d0)

  cnt    = 0.0
  ttl    = 0.0
  ave    = 0.0
  std    = 0.0
  hst(:) = 0.0

  do idate=1,ndate
  do l=1,lall
  do g=1,gall
    if ( msk(g,l,idate) > 0.5 .and. msk_rgn(g,l) > 0.5 ) then
      if( hmin(1).le.var(g,l,idate) .and. var(g,l,idate).lt.hmax(ncell) ) then
        icell = int ( ( var(g,l,idate) - hmin(1) ) / dsiz ) + 1
        cnt        = cnt        + 1.0
        hst(icell) = hst(icell) + 1.0
      end if
      ttl = ttl + 1.
      ave = ave + var(g,l,idate)
    end if
  end do ! gall
  end do ! lall
  end do ! ndate
 
  ave    = ave    / ttl  

  do idate=1,ndate
  do l=1,lall
  do g=1,gall
    if ( msk(g,l,idate) > 0.5 .and. msk_rgn(g,l) > 0.5 ) then
      std = std + ( var(g,l,idate) - ave ) ** 2.0d0
    end if
  end do ! gall
  end do ! lall
  end do ! ndate

  std = sqrt ( std / ttl )

  fnc(:) = 0.0d0
  do icell=1,ncell
    hmdl = 0.5d0 * ( hmin(icell) + hmax(icell) )
    fnc(icell) = sngl( dexp (0.0d0 - 0.5d0*((hmdl-ave)**2.0d0)/std**2.0d0 ) / dsqrt(2.0d0*pi*(std**2.0d0)) )
    hst(icell) = hst(icell) / ttl / sngl( dsiz )
  end do

  if( ttl < 0.5 ) then ! regulation
    ave    = 0.0
    std    = 0.0
    fnc(:) = 0.0
    hst(:) = 0.0
  end if

return
end
!--------------------------------------------------------------------!
subroutine calc_angle_deg(xlon1, xlon2, xlat1, xlat2, angle)
  implicit none
  real(8), intent(in)  :: xlon1, xlon2, xlat1, xlat2
  real(8), intent(out) :: angle
  real(8) :: x1, y1, z1, x2, y2, z2
  real(8), parameter :: pi=3.1415926535d0  
  real(8), parameter :: r180=1.0d0/180.0d0

  !===> using INNER PRODUCT cos0 = A.B / |A||B|, |A,B| is set to be 1
  x1 = dcos(xlon1*pi*r180)*dsin(0.5d0*pi-xlat1*pi*r180)
  y1 = dsin(xlon1*pi*r180)*dsin(0.5d0*pi-xlat1*pi*r180)
  z1 =                     dcos(0.5d0*pi-xlat1*pi*r180)
  x2 = dcos(xlon2*pi*r180)*dsin(0.5d0*pi-xlat2*pi*r180)
  y2 = dsin(xlon2*pi*r180)*dsin(0.5d0*pi-xlat2*pi*r180)
  z2 =                     dcos(0.5d0*pi-xlat2*pi*r180)

  angle = dacos( x1*x2 + y1*y2 + z1*z2 ) *180.0d0 / pi

return
end
!-------------------------------------------------------- ------------!
  subroutine calc_length_deg(xlon1, xlon2, xlat1, xlat2, length, angle) ! modified on Dr. Kondo's program
    implicit none
    real(8), intent(in) :: xlon1, xlon2, xlat1, xlat2
    real(8), intent(out) :: length, angle
    real(8) :: ilat, ilon, jlat, jlon
    real(8) :: A, B
!    real(8), parameter :: Rpo = 6334834.0d0, Req = 6377937.0d0
    real(8), parameter :: Rpo = 6378137.0d0, Req = 6378137.0d0
    real(8), parameter :: RRR = 6378137.0d0

    real(8), parameter :: e = 0.006674d0
    real(8), parameter :: pi=3.1415926535d0
    real(8), parameter :: r180=1.0d0/180.0d0
    real(8) :: P, dP, dR

    ilat = xlat1 * pi * r180 ! [rad]
    ilon = xlon1 * pi * r180 ! [rad]
    jlat = xlat2 * pi * r180 ! [rad]
    jlon = xlon2 * pi * r180 ! [rad]

!    length = dacos ( dsin(ilat)*dsin(jlat) + dcos(jlat)*dcos(ilat)*dcos(ilon-jlon) ) * RRR

    P  = (ilat + jlat)*0.5d0
    dP = dabs(ilat - jlat)
    dR = dabs(ilon - jlon)
    if( dR > pi ) then
      dR = 2.d0*pi - dR
    end if

    A = Rpo/(1.d0 - e*dsin(P)**2)**(1.5d0)
    B = Req/(1.d0 - e*dsin(P)**2)**(0.5d0)
    length = dsqrt( (A*dP)**2 + (B*dcos(P)*dR)**2)    

    angle  = length * 180.0d0 / pi / RRR
  end 