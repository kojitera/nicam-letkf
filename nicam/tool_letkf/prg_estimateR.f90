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
  use common_tvs_nicam
 
  !use fnc_transgauss
  implicit none

  !integer :: g, l , gg, ll
  !integer :: nij0, ngpv 
  !integer :: igrd, ngrd

  !==> allocatable variables
  !real(8), allocatable :: ico_lat(:,:), ico_lon(:,:)
  real(4), allocatable :: tmp(:,:)

  REAL(8),ALLOCATABLE :: tmpelm(:,:,:)
  REAL(8),ALLOCATABLE :: tmplon(:,:,:)
  REAL(8),ALLOCATABLE :: tmplat(:,:,:)
  REAL(8),ALLOCATABLE :: tmpzenith(:,:,:)
  REAL(8),ALLOCATABLE :: tmpskin(:,:,:)
  REAL(8),ALLOCATABLE :: tmpstmp(:,:,:)
  REAL(8),ALLOCATABLE :: tmpclw(:,:,:)
  REAL(8),ALLOCATABLE :: tmplev(:,:,:,:)
  REAL(8),ALLOCATABLE :: tmpdat(:,:,:,:)
  REAL(8),ALLOCATABLE :: tmperr(:,:,:,:)
  REAL(8),ALLOCATABLE :: tmpdep(:,:,:,:)
  REAL(8),ALLOCATABLE :: tmphdxf(:,:,:,:)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:,:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:,:,:,:)
  INTEGER,ALLOCATABLE :: tmpfoot(:,:,:)


  real(4), allocatable :: lnd_rgn(:), lnd_tmp(:,:), lnd_all(:,:)
  real(4), allocatable :: gus_rgn(:), gus_tmp(:,:), gus_all(:,:,:) 
  real(4), allocatable :: anl_rgn(:), anl_tmp(:,:), anl_all(:,:,:)
  real(4), allocatable :: obs_rgn(:), obs_tmp(:,:), obs_all(:,:,:) 
  real(4), allocatable :: ddd_rgn(:), ddd_tmp(:,:), ddd_all(:,:,:) ! Y-Hx^f
  real(4), allocatable :: zzz_rgn(:), zzz_tmp(:,:), zzz_all(:,:,:) ! Y-Hx^a
  real(4), allocatable :: msk_rgn(:), msk_tmp(:,:), msk_all(:,:,:)
  real(4), allocatable :: qc_rgn(:), ppmask_o(:), ppmask_m(:)
  !real(4), allocatable :: ghst(:), ohst(:), dhst(:)
  !real(4), allocatable :: gfnc(:), ofnc(:), dfnc(:)
  !real(4) :: gcnt, ocnt, dcnt
  !real(4) :: gttl, ottl, dttl, gave, oave, dave, gstd, ostd, dstd

  real(8), allocatable :: hmin(:), hmax(:)
  real(8) :: hmdl
  real(8) :: tmpsum(2)

  !==> global variables
  integer, save :: gall, lall
  integer, save :: idate, ndate
  integer, save :: icell, ncell
  integer, save :: itemp, icase, isurf
  integer :: icount
  character(255) :: inpname, fname, mapname

  integer, parameter :: ncase = 4 ! (glb,nph,trp,sph)
  integer, parameter :: nsurf = 3 ! (all,lnd,ocn)
  character(3) :: ccase, csurf
  
  !==> dstat.cnf
  integer, save  :: psdate, pedate
  integer, save  :: nhour  ! assimilation cycle
  real(8)   :: dsiz
  real(8)   :: asiz=100
  real(8)   :: lngmax=1000
  character(255) :: DSTATDIR, LANDDIR

  character(128) :: msg = 'msg'
  character(256) gues_fname
  character(256) anal_fname

  INTEGER,ALLOCATABLE,SAVE :: ntvsgrd(:,:,:,:)

  !==> nhm_driver.cnf
  integer, save :: glevel
  integer, save :: rlevel
  integer, save :: vlayer
  character(128) rgnmngfname
  real(8), parameter :: gaussmin = -7.0d0
  real(8), parameter :: gaussmax =  7.0d0

  !--> paralell computing
  !integer :: ierr, iread, nread, iprocs, nprocs, irank, myrank
  integer :: ierr
  integer :: i, nn, ic
  integer, allocatable :: proc_l(:,:) 
  real :: rtimer00, rtimer

!========================================================================
  CALL ADM_proc_init(ADM_MULTI_PRC)

  if( myrank==0 ) rtimer00=mpi_wtime()
!========================================================================

  !namelist / dstat_param / &
  !  psdate,  pedate, nhour
  !open(1,file='dstat.cnf')
  !read(1,nml=dstat_param)
  !close(1)

  !=====> get ndate
  !idate = psdate ! initial time
  !ndate = 0
  !do
  !  ndate = ndate + 1
! !   if(ADM_prc_me==1) write(6,'(a,3i)') "     CDFMAKE:: CDATE, NCON,
! !   EDATE",idate,ndate,gedate
  !  if( idate == pedate ) goto 10
  !  call update_date(idate,nhour)
  !  if ( ndate .ge. 10000 ) then
  !    write(6,*) 'ndate is too large'
  !    stop
  !  end if
  !end do
!10 continue

  CALL set_instrument
  ALLOCATE(ntvsprofslots(ninstrument,1))
  ALLOCATE(ntvsgrd(nlon,nlat,ninstrument,1))
  ntvsgrd = 0
  ntvsprofslots=0

  gues_fname='AA15_gues.dat'
  anal_fname='AA15_anal.dat'

  ADM_prc_me=myrank+1
  call MISC_make_idstr(fname, trim(msg), 'pe', myrank+1)
  open(ADM_LOG_FID,file=adjustl(trim(fname)),form='formatted')

  CALL get_ntvs_mpi(gues_fname)

  ntvsprofslots(:,1)=ntvsprof(:)
  maxtvsprof=MAXVAL(ntvsprofslots)
  maxtvsfoot=MAXVAL(nfootp)
  write(*,*) ntvsprofslots
  write(*,*) ntvsprof

  ALLOCATE( tmpelm(             maxtvsprof,ninstrument,2) )
  ALLOCATE( tmplon(             maxtvsprof,ninstrument,2) )
  ALLOCATE( tmplat(             maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpzenith(          maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpskin(            maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpstmp(            maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpclw (            maxtvsprof,ninstrument,2) )
  ALLOCATE( tmplev(    maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpdat(    maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmperr(    maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpdep(    maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmphdxf(   maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpqc0(    maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpqc(     maxtvsch,maxtvsprof,ninstrument,2) )
  ALLOCATE( tmpfoot(            maxtvsprof,ninstrument,2) )

  CALL read_tvs_mpi( gues_fname,   tmpelm   (1,1,1), &
            &    tmplon   (1,1,1), tmplat   (1,1,1), &
            &    tmpzenith(1,1,1), tmpskin  (1,1,1), &
            &    tmpstmp  (1,1,1), tmpclw   (1,1,1), &
            &    tmplev (1,1,1,1), tmpdat (1,1,1,1), tmperr(1,1,1,1),&
            &    tmphdxf(1,1,1,1), tmpqc0 (1,1,1,1), tmpfoot (1,1,1))

  CALL read_tvs_mpi( anal_fname,   tmpelm   (1,1,1), &
            &    tmplon   (1,1,2), tmplat   (1,1,2), &
            &    tmpzenith(1,1,2), tmpskin  (1,1,2), &
            &    tmpstmp  (1,1,2), tmpclw   (1,1,2), &
            &    tmplev (1,1,1,2), tmpdat (1,1,1,2), tmperr(1,1,1,2),&
            &    tmphdxf(1,1,1,2), tmpqc0 (1,1,1,2), tmpfoot (1,1,2))


  do nn = 1, ninstrument
  gall=ntvsprof(nn)
  do ic = 1, ntvsch(nn)
    write(*,*) nn, ic, ntvsprof(nn), sum(tmpqc0(ic,1:gall,nn,1))
  end do
  end do

  do nn = 1, ninstrument
  do ic = 1, ntvsch(nn)
    gall=ntvsprof(nn)
    write(*,*) nn, ic, sum(tmpdat (ic,1:gall,nn,1))/gall, &
                       sum(tmphdxf(ic,1:gall,nn,1))/gall, &
                       sum(tmphdxf(ic,1:gall,nn,2))/gall
  end do
  end do

  lall=1
  ndate=1
  allocate( ddd_all(maxtvsprof,maxtvsch,ninstrument) )
  allocate( zzz_all(maxtvsprof,maxtvsch,ninstrument) )
  do nn = 1, ninstrument
  do ic = 1, ntvsch(nn)
    gall=ntvsprof(nn)
    ddd_all(1:gall,ic,nn)=tmpdat(ic,1:gall,nn,2)-tmphdxf(ic,1:gall,nn,2) ! Analysis
    zzz_all(1:gall,ic,nn)=tmpdat(ic,1:gall,nn,1)-tmphdxf(ic,1:gall,nn,1) ! Guess
  end do
  end do
  nprocs=1
  ccase='glb'
  csurf='all'

  do nn = 1, ninstrument
  do ic = 1, ntvsch(nn)
  do i = 1, ntvsprof(nn)
    if(abs(zzz_all(i,ic,nn)) > 0.35*5) then
      !write(*,*) nn, ic, i, zzz_all(i,ic,nn)
      tmpqc0(ic,i,nn,1)=0
    end if
  end do
  end do
  end do

  do nn = 1, ninstrument
  gall=ntvsprof(nn)
  do ic = 1, ntvsch(nn)
    write(*,*) nn, ic, ntvsprof(nn), sum(tmpqc0(ic,1:gall,nn,1))
  end do
  end do

  do nn = 1, ninstrument
  gall=ntvsprof(nn)
  do ic = 1, ntvsch(nn)
    tmpsum(:)=0.0d0
    icount=0
    do i = 1, ntvsprof(nn)
      if(tmpqc0(ic,i,nn,1)==1) then
        tmpsum(1)=tmpsum(1)+ddd_all(i,ic,nn)
        tmpsum(2)=tmpsum(2)+zzz_all(i,ic,nn)
        icount=icount+1
      end if
    end do
    tmpsum(1)=tmpsum(1)/sum(tmpqc0(ic,1:gall,nn,1))
    tmpsum(2)=tmpsum(2)/sum(tmpqc0(ic,1:gall,nn,1))
    ddd_all(1:gall,ic,nn)=ddd_all(1:gall,ic,nn) - tmpsum(1)
    zzz_all(1:gall,ic,nn)=zzz_all(1:gall,ic,nn) - tmpsum(2)
  end do
  end do

  open(10,file='10.txt',action='write')
  open(11,file='11.txt',action='write')
  do nn = 1, ninstrument
  do ic = 1, ntvsch(nn)
    gall=ntvsprof(nn)
    call Desroziers(nn, ic, gall,       &
                    lall,       &
                    ndate,      &
                    real(tmpqc0(ic,1:gall,nn,1)),   &
                    real(tmpqc0(ic,1:gall,nn,1)),   &
                    tmplon(1:gall,nn,1),     &
                    tmplat(1:gall,nn,1),     &
                    ddd_all(1:gall,ic,nn),    &
                    zzz_all(1:gall,ic,nn),    &
                    asiz,       &
                    lngmax,     &
                    myrank,     &
                    nprocs,     &
                    rtimer00,   &
                    10,         &
                    11,         &
                    ccase,      &
                    csurf)
    write(*,*) nn, ic, 'finished'
  end do
  end do
  close(10)
  close(11)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if( myrank==0 ) rtimer=mpi_wtime()
  if( myrank==0 ) write(6,'(A,F13.2,A,I)') '### TIMER(END  READ TRANSFORM):',rtimer-rtimer00,' // NPROCS ',nprocs
  call MPI_Finalize(ierr)

end program
!=========================================================
subroutine Desroziers(inst,ic,gall,lall,ndate,msk,msk_rgn,ico_lon,ico_lat,ddd,zzz,asiz,lngmax,myrank,nprocs,rtimer00,fnum,mnum,ccase,csurf)
!=========================================================
!!! Dstat with Desroziers et al. (2005) R=<(y-Hx^a)(y-Hx^f)T>
!!! Dstat with Hollingsworth (1986)     <ddT> = R + HBHT
  use mpi
  implicit none
  integer, intent(in)  :: gall, lall, ndate, myrank, nprocs, fnum, mnum
  real(8), intent(in)  :: asiz, lngmax, ico_lon(gall,lall), ico_lat(gall,lall)
  real(4), intent(in)  :: ddd(gall,lall,ndate), zzz(gall,lall,ndate), msk(gall,lall,ndate), msk_rgn(gall,lall)
  real(4), intent(in)  :: rtimer00
  integer, intent(in)  :: inst, ic
  character(3), intent(in) :: ccase, csurf
  !===> global variables
  integer :: g, l, idate, icon, nmax, nsamp, ilng, nlng
  real(4) :: rtimer, tmpmsk(gall,lall)
  real(4), allocatable :: d(:), z(:)  
  real(8), allocatable :: dlon(:), dlat(:), acnt(:), ades(:), ahol(:), afnc(:)
  real(8), allocatable :: abes(:)
  real(4), parameter :: undef = -9999.
  real(8) :: length, angle,   RRR, RHBHT

  !===> cnt
  nlng  = int( lngmax / asiz )  
  allocate ( acnt(0:nlng), ades(0:nlng), ahol(0:nlng), afnc(0:nlng), abes(0:nlng) )

  write(*,*) 'lngmax =', lngmax
  write(*,*) 'asiz   =', asiz
  write(*,*) 'nlng   =', nlng

  !===> get maxval
  nmax = 0
  do l=1,lall
  do g=1,gall
    tmpmsk(g,l) = maxval( msk(g,l,1:ndate) )
    if( tmpmsk(g,l)>0.5 ) nmax = nmax + 1
  end do
  end do
  allocate ( d(nmax), z(nmax), dlon(nmax), dlat(nmax) )

  write(*,*) 'nmax   =', nmax

  acnt(0:nlng) = 0.0d0 
  ades(0:nlng) = 0.0d0
  abes(0:nlng) = 0.0d0
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

    write(*,*) 'nsamp  =', nsamp

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
        abes(ilng) = abes(ilng) + z(g)*z(l) !  <(y-Hx^f)(y-Hx^f)T>
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
      abes(ilng) = abes(ilng) / acnt(ilng) 
      ahol(ilng) = ahol(ilng) / acnt(ilng)   
    end if
!!!    print '(4f)', length, acnt(ilng), ades(ilng), ahol(ilng)
    write(fnum,'(2i6,f8.1,i15,6f15.10)') inst, ic, length, int(acnt(ilng)), &
      ades(ilng), ades(ilng)/ades(0), &
      abes(ilng), abes(ilng)/abes(0), &
      !abes(ilng)-ades(ilng), &
      !(abes(ilng)/abes(0)) - (ades(ilng)/ades(0)), &
      ahol(ilng), ahol(ilng)/ahol(0)
  end do
  
  RRR    = ades(0)
  RHBHT  = ahol(0)

  print '(a,3f15.10,2a)', "     Estimated R+HBHT, R, R(dev)  ", RHBHT, RRR, dsqrt(RRR), " in ", ccase
  write(mnum,'(3f15.10,3a4)') RHBHT, RRR, dsqrt(RRR),"  ",ccase, csurf

  FLUSH(fnum)
  FLUSH(mnum)

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

!=========================================================
subroutine update_date(date,nhour)
!=========================================================
  implicit none
  integer date,nhour,date1
  integer year,mon,day,hour,mondays

  date1 = date
  year = int( date / 1000000 )
  date = date - year * 1000000
  mon  = int( date / 10000 )
  date = date - mon  * 10000
  day  = int( date / 100 )
  hour = date - day * 100

  hour = hour + nhour
  if( hour.ge.24 ) then
    hour = hour - 24
    day  = day + 1
    call get_mondays(year,mon,mondays)
    if ( day.gt.mondays ) then
      day = 1
      mon = mon + 1
      if( mon.gt.13 ) then
        mon = 1
        year = year + 1
      end if ! mon
    end if ! day
  end if ! hour

  date = year*1000000 + mon*10000 + day*100 + hour
end subroutine update_date

!=========================================================
subroutine get_mondays(year,mon,mondays)
!=========================================================
  implicit none
  integer year,mon,mondays

  mondays=31
  if( mon ==  4 ) mondays = 30
  if( mon ==  6 ) mondays = 30
  if( mon ==  9 ) mondays = 30
  if( mon == 11 ) mondays = 30
  if( mon ==  2 ) then
    mondays = 28
    if( mod(year,4)==0 ) mondays = 29
  end if

end subroutine get_mondays

