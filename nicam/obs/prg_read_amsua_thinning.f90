PROGRAM MAIN
  use mpi
  use mod_calendar
  IMPLICIT NONE
!
  INTEGER,PARAMETER :: ninstrument=11
  INTEGER,PARAMETER :: maxtvsch=10
  CHARACTER(4),SAVE :: tvsname(ninstrument)
  INTEGER,SAVE :: tvsinst(3,ninstrument)
  INTEGER,SAVE :: tvsch(maxtvsch,ninstrument)
  INTEGER,SAVE :: ntvsch(ninstrument)

  REAL(4),PARAMETER :: oint=250.0d3
  REAL(4),PARAMETER :: radius=40000.0d3

  INTEGER,allocatable :: nx1(:)
  INTEGER :: nx
  INTEGER :: ny
  INTEGER,PARAMETER :: maxlev = 1     !Maximum number of BUFR levels
  INTEGER,PARAMETER :: maxevn = 15      !Maximum number of BUFR event sequences

  CHARACTER(8) :: obtype
  INTEGER,PARAMETER :: nobtype = 20
  INTEGER,PARAMETER :: nslots = 7

  INTEGER :: i, j
  INTEGER :: ix, iy
  INTEGER :: icount
  INTEGER :: IREADNS
  INTEGER :: idate,idummy,nlev,n
  INTEGER :: ilev,ievn
  INTEGER :: iobs(nobtype+1),iobs_out(nobtype+1)
  CHARACTER(6) :: obtypelist(nobtype)
  REAL(8) :: station(20)
  CHARACTER(8) :: cs
  REAL(8) :: obs(15,maxlev,maxevn)
  REAL(4) :: wk(34)
  ! 1: Observation ID
  ! 2: Longitude
  ! 3: Latitude
  ! 4: Elevation =0
  ! 5: Observation
  ! 6: Observation error =0
  ! 7: Satelite zenith angle
  ! 8: Bearing of azimuth
  ! 9: solar zenith angle
  !10: solar azimuth
  INTEGER :: islots, iunit, kslots
  INTEGER :: jslots(nslots)=(/4,3,5,2,6,1,7/)
  INTEGER :: iret
  INTEGER,PARAMETER :: num_satellite=5
  INTEGER :: iobs_noaa15(nslots)
  INTEGER :: iobs_noaa16(nslots)
  INTEGER :: iobs_noaa17(nslots)
  INTEGER :: iobs_noaa18(nslots)
  INTEGER :: iobs_noaa19(nslots)
  INTEGER :: iobs_metop2(nslots)
  INTEGER :: nobs_noaa15(nslots)
  INTEGER :: nobs_noaa16(nslots)
  INTEGER :: nobs_noaa17(nslots)
  INTEGER :: nobs_noaa18(nslots)
  INTEGER :: nobs_noaa19(nslots)
  INTEGER :: nobs_metop2(nslots)
  INTEGER :: nobs_tmp(num_satellite,nslots)
  INTEGER :: nobs_out(num_satellite,nslots)
  INTEGER,PARAMETER :: id_NOAA15=206
  INTEGER,PARAMETER :: id_NOAA16=207
  INTEGER,PARAMETER :: id_NOAA17=208
  INTEGER,PARAMETER :: id_NOAA18=209
  INTEGER,PARAMETER :: id_NOAA19=223
  INTEGER,PARAMETER :: id_METOP2=4
  REAL(4) :: odata(34,100000,num_satellite,nslots)
  REAL(4) :: odata_out(34,100000,num_satellite,nslots)
  INTEGER(1) :: odata_flag(100000,num_satellite,nslots)
  REAL(4) :: odata_tmp(34,1)
  REAL(4) :: UNDEF=-999.
  REAL(8) :: DATE
  REAL(8) :: sdate, edate, cdate
  REAL(8) :: DATE_MAX, DATE_MIN
  INTEGER(4) :: SYEAR, SMON, SDAY, SHR, SMIN, SSEC
  INTEGER(4) :: EYEAR, EMON, EDAY, EHR, EMIN, ESEC
  INTEGER :: cdate_tmp(6)
  INTEGER :: sdate_tmp(6)
  INTEGER :: edate_tmp(6)

  REAL(4) :: dist, dist_limit
  REAL(4) :: dist_min
  character(256) :: icolatlon_fname
  REAL(8), allocatable :: icolat(:,:)
  REAL(8), allocatable :: icolon(:,:)
  REAL(4), allocatable :: lat(:), lon(:,:)
  REAL(4) :: tmplon

  REAL(8) :: time(10)
  integer,parameter:: &
     &    rttv_plat_noaa =1 ,rttv_plat_dmsp =2 &
     &   ,rttv_plat_meteo=3 ,rttv_plat_goes =4 &
     &   ,rttv_plat_gms  =5 ,rttv_plat_fy2  =6 &
     &   ,rttv_plat_trmm =7 ,rttv_plat_ers  =8 &
     &   ,rttv_plat_eos  =9 ,rttv_plat_metop=10 &
     &   ,rttv_plat_envi =11,rttv_plat_msg  =12 &
     &   ,rttv_plat_fy1  =13,rttv_plat_adeos=14 &
     &   ,rttv_plat_mtsat=15,rttv_plat_cori =16 &
     &   ,rttv_plat_npoes=17,rttv_plat_gifts=18 &
     &   ,rttv_plat_metop2=19

  integer,parameter:: &
     &    rttv_inst_hirs  =0, rttv_chan_hirs  =19   & ! HIRS
     &   ,rttv_inst_amsua =3, rttv_chan_amsua =15   & ! AMSU-A
     &   ,rttv_inst_amsub =4, rttv_chan_amsub =5    & ! AMSU-B
     &   ,rttv_inst_ssmi  =6, rttv_chan_ssmi  =7    & ! SSMI
     &   ,rttv_inst_tmi   =9, rttv_chan_tmi   =9    & ! TMI
     &   ,rttv_inst_ssmis =10,rttv_chan_ssmis =24   & ! SSMIS
     &   ,rttv_inst_airs  =11,rttv_chan_airs  =2378 & ! AIRS
     &   ,rttv_inst_amsr  =17,rttv_chan_amsr  =14   & ! AMSR
     &   ,rttv_inst_mviri =20,rttv_chan_mviri =2    & ! METEOSAT
     &   ,rttv_inst_seviri=21,rttv_chan_seviri=8    & ! MSG
     &   ,rttv_inst_goesi =22,rttv_chan_goesi =4    & ! GOES-IMAGER(IR)
     &   ,rttv_inst_mtsati=24,rttv_chan_mtsati=4      ! MTSAT


  integer :: ierr

  real(4) :: dx, dy
  real(4) :: sx, sy
  INTEGER :: tmp_i, tmp_n

  CALL MPI_INIT(ierr)

  ny = floor(radius/2.0/oint)
  allocate( lat(ny) )
  allocate(nx1(ny))
  sy =  -90.0 + 180.0/float(ny)*0.5
  dy = 180.0/float(ny)
  do iy = 1, ny
    lat(iy) = sy + dy * (iy-1)
  end do

  do iy = 1, ny
    nx1(iy)=ceiling(radius*cos(lat(iy)/180.0*3.1415)/oint)
  end do

  allocate( lon(maxval(nx1(:)),ny) )

  do iy = 1, ny
    sx = -180.0 + 360.0/float(nx1(iy))*0.5
    dx = 360.0/float(nx1(iy))
    do ix = 1, nx1(iy)
      lon(ix,iy) = sx + dx * ( ix-1 ) 
    end do
  end do
  !read(5,*) SYEAR, SMON, SDAY, SHR, SMIN, SSEC, EYEAR, EMON, EDAY, EHR, EMIN, ESEC, dist_limit
  read(5,*) cdate_tmp(1:6), dist_limit
  write(*,*) cdate_tmp(1:6)
  OPEN(11,FILE='amsua.in',FORM='unformatted')
  CALL OPENBF(11,'IN',11)
  CALL DATELEN(10)

  iobs = 0
  iobs_out = 0

! OUTPUT BUFR TABLE 
  OPEN(200,FILE='amsua.txt')
  CALL DXDUMP(11,200)

  wk(34)=21023
  icount = 0

  iobs = 1
  iobs_noaa15=1
  iobs_noaa16=1
  iobs_noaa18=1
  iobs_noaa19=1
  iobs_metop2=1

  nobs_noaa15=0
  nobs_noaa16=0
  nobs_noaa18=0
  nobs_noaa19=0
  nobs_metop2=0

  call calendar_yh2ss(cdate, cdate_tmp)

  date_max=-9.9*10**30
  date_min=1.0d33

  time(1)=MPI_WTIME()

  DO
    IF(IREADNS(11,obtype,idate) /= 0) EXIT
    !DO n=1,nobtype+1
      !IF(obtype == obtypelist(n) .OR. n > nobtype) THEN
      !  iobs(n) = iobs(n)+1
      !  EXIT
      !END IF
    !END DO

    CALL UFBINT(11,station(1),20,1,idummy,&
    'YEAR MNTH DAYS HOUR MINU &
     SECO CLAT CLON SAID SIID &
     FOVN LSQL SAZA SOZA HOLS &
     HMSL')
    CALL UFBINT(11,station(17),20,1,idummy,&
    'SOLAZI BEARAZ')
    !'YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA HOLS HMSL SOLAZI BEARAZ')
    !'YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAID SIID' )

! BUFR TABLE FOR AMSU-A
!	YEAR  MNTH  DAYS  HOUR  MINU  SECO  207002  CLAT  CLON  207000 
!	SAID  SIID  FOVN  LSQL  SAZA  SOZA  HOLS  202127  HMSL  202000 
!	SOLAZI  BEARAZ  "BRITCSTC"15                                 

    CALL UFBEVN(11,obs,15,maxlev,maxevn,ilev,'TMBR')

    where( abs(obs(1,1,1:15)) > 1000.0 ) obs(1,1,1:15)=UNDEF 
      
    call calendar_yh2ss(date, int(station(1:6)))

    islots=4+NINT((date-cdate)/3600.d0)
!
      wk(1:18)=sngl(station(1:18))
      wk(19:33)=sngl(obs(1,1,1:15))
!  
      write(*,*) wk(7:9), wk(13), wk(11), wk(15), wk(16)
      if( wk(9) == id_NOAA15 ) then
         odata(1:34,iobs_noaa15(islots),1,islots)=wk(1:34)
         iobs_noaa15(islots)=iobs_noaa15(islots)+1
         nobs_noaa15(islots)=nobs_noaa15(islots)+1
      else if( wk(9) == id_NOAA16 ) then
         odata(1:34,iobs_noaa16(islots),2,islots)=wk(1:34)
         iobs_noaa16(islots)=iobs_noaa16(islots)+1
         nobs_noaa16(islots)=nobs_noaa16(islots)+1
      else if( wk(9) == id_NOAA17 ) then

      else if( wk(9) == id_NOAA18 ) then
         odata(1:34,iobs_noaa18(islots),3,islots)=wk(1:34)
         iobs_noaa18(islots)=iobs_noaa18(islots)+1
         nobs_noaa18(islots)=nobs_noaa18(islots)+1
      else if( wk(9) == id_NOAA19 ) then
         odata(1:34,iobs_noaa19(islots),4,islots)=wk(1:34)
         iobs_noaa19(islots)=iobs_noaa19(islots)+1
         nobs_noaa19(islots)=nobs_noaa19(islots)+1
      else if( wk(9) == id_METOP2 ) then
         odata(1:34,iobs_metop2(islots),5,islots)=wk(1:34)
         iobs_metop2(islots)=iobs_metop2(islots)+1
         nobs_metop2(islots)=nobs_metop2(islots)+1
      end if
  END DO
!
  time(2)=MPI_WTIME()
!
  nobs_tmp(1,:)=nobs_noaa15(:)
  nobs_tmp(2,:)=nobs_noaa16(:)
  nobs_tmp(3,:)=nobs_noaa18(:)
  nobs_tmp(4,:)=nobs_noaa19(:)
  nobs_tmp(5,:)=nobs_metop2(:)

  write(*,*) '############### NUMBER OF INPUT OBSERVATIONS ###############'
  write(*,'(8A7)') 'SLOTS', '-3', '-2', '-1', '0', '+1', '+2', '+3'
  write(*,'(A,7i7)') 'NOAA-15', nobs_noaa15(:)
  write(*,'(A,7i7)') 'NOAA-16', nobs_noaa16(:)
  write(*,'(A,7i7)') 'NOAA-18', nobs_noaa18(:)
  write(*,'(A,7i7)') 'NOAA-19', nobs_noaa19(:)
  write(*,'(A,7i7)') 'METOP-2', nobs_metop2(:)
  write(*,'(A,i10)') 'TOTAL NUMBER OF OBSERVATIONS', sum(nobs_tmp(:,:))
!
  iobs_noaa15=1
  iobs_noaa16=1
  iobs_noaa18=1
  iobs_noaa19=1
  iobs_metop2=1
!
  nobs_noaa15=0
  nobs_noaa16=0
  nobs_noaa18=0
  nobs_noaa19=0
  nobs_metop2=0

!!-----------------------------------------------------------------------
!!  SUPEROBBING THE OBSERATION DATA
!!-----------------------------------------------------------------------
  odata_flag(:,:,:)=0
  do iy = 1, ny
    do ix = 1, nx1(iy)
      do kslots = 1, nslots
        islots=jslots(kslots) 
        dist_min=9999999999.0
        icount=0
        do n =  1, num_satellite
          !odata_tmp(:,:)=0.0
          !$omp parallel
          !$omp do reduction(+:icount)
          do i = 1, nobs_tmp(n,islots)
            if( odata(7,i,n,islots) > lat(iy)-0.5 .and. odata(7,i,n,islots) < lat(iy)+0.5 ) then
              if( odata(8,i,n,islots) > lon(ix,iy)-0.5 .and. odata(8,i,n,islots) < lon(ix,iy)+0.5 ) then
                if( odata_flag(i,n,islots) == 0 ) then
                  call com_distll_1(lon(ix,iy),lat(iy),odata(8,i,n,islots),odata(7,i,n,islots),&
                       & dist*cos(lat(iy)/180.0d0*3.1415d0))
                  if(dist<dist_min) then
                    icount=1
                    dist_min=dist
                    odata_tmp(:,icount)=odata(:,i,n,islots)
                    !odata_flag(tmp_i, tmp_n)=0
                    tmp_i=i
                    tmp_n=n
                  end if
                end if
              end if
            end if
          end do
        end do
        if(icount==1) then
          odata_flag(tmp_i,tmp_n,islots)=1
          exit
        end if
      end do
!        !$omp end do
!        !$omp end parallel
!    
      if(icount>0) then
        if( odata_tmp(9,icount) == id_NOAA15 ) then
           do i = 1, 34
             odata_out(i,iobs_noaa15,1,islots)=odata_tmp(i,1)
           end do 
           iobs_noaa15(islots)=iobs_noaa15(islots)+1
           nobs_noaa15(islots)=nobs_noaa15(islots)+1
        else if( odata_tmp(9,icount) == id_NOAA16 ) then
           do i = 1, 34
             odata_out(i,iobs_noaa16,2,islots)=odata_tmp(i,1)
           end do 
           iobs_noaa16(islots)=iobs_noaa16(islots)+1
           nobs_noaa16(islots)=nobs_noaa16(islots)+1
        else if( odata_tmp(9,icount) == id_NOAA17 ) then
 
        else if( odata_tmp(9,icount) == id_NOAA18 ) then
           do i = 1, 34
             odata_out(i,iobs_noaa18,3,islots)=odata_tmp(i,1)
           end do 
           iobs_noaa18(islots)=iobs_noaa18(islots)+1
           nobs_noaa18(islots)=nobs_noaa18(islots)+1
        else if( odata_tmp(9,icount) == id_NOAA19 ) then
           do i = 1, 34
             odata_out(i,iobs_noaa19,4,islots)=odata_tmp(i,1)
           end do 
           iobs_noaa19(islots)=iobs_noaa19(islots)+1
           nobs_noaa19(islots)=nobs_noaa19(islots)+1
        else if( odata_tmp(9,icount) == id_METOP2 ) then
           do i = 1, 34
             odata_out(i,iobs_metop2,5,islots)=odata_tmp(i,1)
           end do 
           iobs_metop2(islots)=iobs_metop2(islots)+1
           nobs_metop2(islots)=nobs_metop2(islots)+1
        end if
      end if
!    
    end do
  end do

  time(3)=MPI_WTIME()
!      
  nobs_tmp(1,:)=nobs_noaa15(:)
  nobs_tmp(2,:)=nobs_noaa16(:)
  nobs_tmp(3,:)=nobs_noaa18(:)
  nobs_tmp(4,:)=nobs_noaa19(:)
  nobs_tmp(5,:)=nobs_metop2(:)
!  write(*,*) 'nobs_tmp'
!  write(*,*) nobs_tmp(:)
!    
!  OPEN(101,FILE='amsua.out',FORM='unformatted',access='sequential')
!    
  nobs_out(:,:)=0
  do islots = 1, nslots
    do n =  1, num_satellite
      do i = 1, nobs_tmp(n,islots)
        !if(mod(i,1)==0) then
          !nobs_out(n,islots)=nobs_out(n,islots)+1
          if(odata_out(8,i,n,islots) < 0.0) then
            tmplon=odata_out(8,i,n,islots)+360.0d0
          else
            tmplon=odata_out(8,i,n,islots)
          end if
          !if(odata_out(12,i,n,islots) == 1) then
            nobs_out(n,islots)=nobs_out(n,islots)+1
            write(100,'(2i7,3f10.3)') islots, n, odata_out(9,i,n,islots), tmplon, &
                                  odata_out(7,i,n,islots)
          !end if
      end do
    end do
  end do

  do islots = 1, nslots
    iunit=90+islots
    !write(*,*) num_satellite
    !write(*,*) sum(nobs_out(:,islots))
    !write(*,*) nobs_out(1:num_satellite,islots)
    write(iunit) num_satellite
    write(iunit) sum(nobs_out(:,islots))
    write(iunit) nobs_out(1:num_satellite,islots)
    do n =  1, num_satellite
      do i = 1, nobs_tmp(n,islots)
        !if(mod(i,1)==0) then
        !if(odata_out(12,i,n,islots) == 1) then ! output over ocean
          write(iunit) odata_out(1:34,i,n,islots)
        !end if
      end do
    end do
  end do
!
!
  time(4)=MPI_WTIME()

  write(*,*) '############### NUMBER OF OUTPUT OBSERVATIONS ###############'
  write(*,'(8A7)') 'SLOTS', '-3', '-2', '-1', '0', '+1', '+2', '+3'
  write(*,'(A,7i7)') 'NOAA-15', nobs_out(1,:)
  write(*,'(A,7i7)') 'NOAA-16', nobs_out(2,:)
  write(*,'(A,7i7)') 'NOAA-18', nobs_out(3,:)
  write(*,'(A,7i7)') 'NOAA-19', nobs_out(4,:)
  write(*,'(A,7i7)') 'METOP-2', nobs_out(5,:)
  write(*,'(A,i10)') 'TOTAL NUMBER OF OBSERVATIONS', sum(nobs_out(:,:))


  write(*,*) 'READING OBS DATA', time(2)-time(1)
  write(*,*) 'AVERAGE OBS DATA', time(3)-time(2)
  write(*,*) 'OUTPUT OBS DATA ', time(4)-time(3)
!  write(*,*) 'date min        ', date_min
!  write(*,*) 'date max        ', date_max

  CALL MPI_FINALIZE(ierr)

CONTAINS
!-----------------------------------------------------------------------
! LIST ASSIMILATED INSTRUMENT CHANNELS
!-----------------------------------------------------------------------
SUBROUTINE set_instrument
  tvsch = 0
  !
  ! NOAA-15 AMSU-A
  !
  tvsname(1) = 'AA15'
  tvsinst(1,1) = rttv_plat_noaa
  tvsinst(2,1) = 15
  tvsinst(3,1) = rttv_inst_amsua
  tvsch(1,1)  = 4
  tvsch(2,1)  = 5
  tvsch(3,1)  = 6
  tvsch(4,1)  = 7
  tvsch(5,1)  = 8
  tvsch(6,1)  = 9
  tvsch(7,1)  = 10
  tvsch(8,1)  = 12
  tvsch(9,1)  = 13
  ntvsch(1)=9
  !
  ! NOAA-16 AMSU-A
  !
  tvsname(2) = 'AA16'
  tvsinst(1,2) = rttv_plat_noaa
  tvsinst(2,2) = 16
  tvsinst(3,2) = rttv_inst_amsua
  tvsch(1,2)  = 4
  tvsch(2,2)  = 5
  tvsch(3,2)  = 6
  tvsch(4,2)  = 7
  tvsch(5,2)  = 8
  tvsch(6,2)  = 9
  tvsch(7,2)  = 10
  tvsch(8,2)  = 11
  tvsch(9,2)  = 12
  tvsch(10,2) = 13
  ntvsch(2)=10
  !
  ! NOAA-18 AMSU-A
  !
  tvsname(2) = 'AA18'
  tvsinst(1,2) = rttv_plat_noaa
  tvsinst(2,2) = 18
  tvsinst(3,2) = rttv_inst_amsua
  tvsch(1,2)  = 4
  tvsch(2,2)  = 5
  tvsch(3,2)  = 6
  tvsch(4,2)  = 7
  tvsch(5,2)  = 8
  tvsch(6,2)  = 9
  tvsch(7,2)  = 10
  tvsch(8,2)  = 11
  tvsch(9,2)  = 12
  tvsch(10,2) = 13
  ntvsch(2)=10
  !
  ! NOAA-19 AMSU-A
  !
  tvsname(4) = 'AA19'
  tvsinst(1,4) = rttv_plat_noaa
  tvsinst(2,4) = 19
  tvsinst(3,4) = rttv_inst_amsua
  tvsch(1,4)  = 4
  tvsch(2,4)  = 5
  tvsch(3,4)  = 6
  tvsch(4,4)  = 7
  tvsch(5,4)  = 8
  tvsch(6,4)  = 9
  tvsch(7,4)  = 10
  tvsch(8,4)  = 11
  tvsch(9,4)  = 12
  tvsch(10,4) = 13
  ntvsch(4)=10
  !
  ! METOP-2 AMSU-A
  !
  tvsname(4) = 'MA02'
  tvsinst(1,4) = rttv_plat_metop2
  tvsinst(2,4) = 19
  tvsinst(3,4) = rttv_inst_amsua
  tvsch(1,4)  = 4
  tvsch(2,4)  = 5
  tvsch(3,4)  = 6
  tvsch(4,4)  = 8
  tvsch(5,4)  = 9
  tvsch(6,4)  = 10
  tvsch(7,4)  = 11
  tvsch(8,4)  = 12
  tvsch(9,4)  = 13
  ntvsch(4)=9
  RETURN
END SUBROUTINE set_instrument

END PROGRAM MAIN

  subroutine MISC_sec2ymdhms(           &
       ctime_second,                    & !--- IN : target time [sec]
       year, month, day, hour, min, sec & !--- OUT : year, month day hour minite, second
       )
    !
    implicit none
    !
    real(8), intent(in)  :: ctime_second
    integer, intent(out) :: sec
    integer, intent(out) :: min
    integer, intent(out) :: hour
    integer, intent(out) :: day
    integer, intent(out) :: month
    integer, intent(out) :: year
    !
!    integer(4) :: second
    integer :: second      ! S.Iga 060212
    !
    second = idnint(ctime_second)
    year = second/(60*60*24*30*12)
    !
    second = mod(second,60*60*24*30*12)
    month = second/(60*60*24*30)
    !
    second = mod(second,60*60*24*30)
    day = second/(60*60*24)
    !
    second = mod(second,60*60*24)
    hour = second/(60*60)
    !
    second = mod(second,60*60)
    min = second/(60)
    !
    second = mod(second,60)
    sec = second

  end subroutine MISC_sec2ymdhms
  !-----------------------------------------------------------------------------
  subroutine MISC_ymdhms2sec(            &
       ctime_second,                     & !--- OUT : target time [sec]
       year, month, day, hour, min, sec  & !--- IN : year, month day hour minite, second
       )
    !
    implicit none
    !
    real(8), intent(out) :: ctime_second
    integer, intent(in)  :: sec
    integer, intent(in)  :: min
    integer, intent(in)  :: hour
    integer, intent(in)  :: day
    integer, intent(in)  :: month
    integer, intent(in)  :: year
    !
  write(*,*) year, month, day, hour, min, sec

    ctime_second &
         = real( year*(60*60*24*30*12) &
               + month*(60*60*24*30)   &
               + day*(60*60*24)        &
               + hour*(60*60)          &
               + min*60                &
               + sec, kind=8           )
  write(*,*) 'ctime_second', ctime_second

  end subroutine MISC_ymdhms2sec
  !-----------------------------------------------------------------------------
  SUBROUTINE com_distll_1(alon,alat,blon,blat,dist)
    IMPLICIT NONE
    REAL(4),INTENT(IN) :: alon
    REAL(4),INTENT(IN) :: alat
    REAL(4),INTENT(IN) :: blon
    REAL(4),INTENT(IN) :: blat
    REAL(4),INTENT(OUT) :: dist
    REAL(4),PARAMETER :: r180=1.0/180.0
    REAL(4) :: lon1,lon2,lat1,lat2
    REAL(4) :: cosd
    REAL(4) :: pi=3.1415926535
    REAL(4) :: re=6371.3e3
  
    lon1 = alon * pi * r180
    lon2 = blon * pi * r180
    lat1 = alat * pi * r180
    lat2 = blat * pi * r180
  
    cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
    cosd = MIN( 1.0,cosd)
    cosd = MAX(-1.0,cosd)
  
    dist = ACOS( cosd ) * re
  
    RETURN
  END SUBROUTINE com_distll_1


!-----------------------------------------------------------------------
! LIST ASSIMILATED INSTRUMENT CHANNELS
!-----------------------------------------------------------------------
!SUBROUTINE set_instrument
!  tvsch = 0
!  !
!  ! NOAA-15 AMSU-A
!  !
!  tvsname(1) = 'AA15'
!  tvsinst(1,1) = rttv_plat_noaa
!  tvsinst(2,1) = 15
!  tvsinst(3,1) = rttv_inst_amsua
!  tvsch(1,1)  = 4
!  tvsch(2,1)  = 5
!  tvsch(3,1)  = 6
!  tvsch(4,1)  = 7
!  tvsch(5,1)  = 8
!  tvsch(6,1)  = 9
!  tvsch(7,1)  = 10
!  tvsch(8,1)  = 12
!  tvsch(9,1)  = 13
!  ntvsch(1)=9
!  !
!  ! NOAA-16 AMSU-A
!  !
!  tvsname(2) = 'AA16'
!  tvsinst(1,2) = rttv_plat_noaa
!  tvsinst(2,2) = 16
!  tvsinst(3,2) = rttv_inst_amsua
!  tvsch(1,2)  = 4
!  tvsch(2,2)  = 5
!  tvsch(3,2)  = 6
!  tvsch(4,2)  = 7
!  tvsch(5,2)  = 8
!  tvsch(6,2)  = 9
!  tvsch(7,2)  = 10
!  tvsch(8,2)  = 11
!  tvsch(9,2)  = 12
!  tvsch(10,2) = 13
!  ntvsch(2)=10
!
!  !
!  ! NOAA-18 AMSU-A
!  !
!  tvsname(2) = 'AA18'
!  tvsinst(1,2) = rttv_plat_noaa
!  tvsinst(2,2) = 18
!  tvsinst(3,2) = rttv_inst_amsua
!  tvsch(1,2)  = 4
!  tvsch(2,2)  = 5
!  tvsch(3,2)  = 6
!  tvsch(4,2)  = 7
!  tvsch(5,2)  = 8
!  tvsch(6,2)  = 9
!  tvsch(7,2)  = 10
!  tvsch(8,2)  = 11
!  tvsch(9,2)  = 12
!  tvsch(10,2) = 13
!  ntvsch(2)=10
!  !
!  ! NOAA-19 AMSU-A
!  !
!  tvsname(4) = 'AA19'
!  tvsinst(1,4) = rttv_plat_noaa
!  tvsinst(2,4) = 19
!  tvsinst(3,4) = rttv_inst_amsua
!  tvsch(1,4)  = 4
!  tvsch(2,4)  = 5
!  tvsch(3,4)  = 6
!  tvsch(4,4)  = 7
!  tvsch(5,4)  = 8
!  tvsch(6,4)  = 9
!  tvsch(7,4)  = 10
!  tvsch(8,4)  = 11
!  tvsch(9,4)  = 12
!  tvsch(10,4) = 13
!  ntvsch(4)=10
!  !
!  ! METOP-2 AMSU-A
!  !
!  tvsname(4) = 'MA02'
!  tvsinst(1,4) = rttv_plat_metop2
!  tvsinst(2,4) = 19
!  tvsinst(3,4) = rttv_inst_amsua
!  tvsch(1,4)  = 4
!  tvsch(2,4)  = 5
!  tvsch(3,4)  = 6
!  tvsch(4,4)  = 8
!  tvsch(5,4)  = 9
!  tvsch(6,4)  = 10
!  tvsch(7,4)  = 11
!  tvsch(8,4)  = 12
!  tvsch(9,4)  = 13
!  ntvsch(4)=9
!  RETURN
!END SUBROUTINE set_instrument

