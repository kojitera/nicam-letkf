PROGRAM MAIN
  use mpi
  IMPLICIT NONE
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_ps_obs=14593

  INTEGER,PARAMETER :: nx     = 160
  INTEGER,PARAMETER :: ny     = 80

  CHARACTER(8) :: obtype
  INTEGER,PARAMETER :: nobtype = 20

  INTEGER :: i, j, n, k
  INTEGER :: ix, iy
  REAL(4) :: wk(7)

  REAL(4) :: dist, dist_limit
  REAL(4) :: dist_min
  LOGICAL :: out_ps
  character(256) :: icolatlon_fname
  REAL(4) :: lat(ny), lon(nx)

  REAL(4) :: time(10)
  integer :: ierr

  real(4) :: dx, dy
  real(4) :: sx, sy
  INTEGER :: tmp_i, tmp_n

  INTEGER :: iobs_u
  INTEGER :: iobs_v
  INTEGER :: iobs_t
  INTEGER :: iobs_q
  INTEGER :: iobs_p

  INTEGER :: nobs_u
  INTEGER :: nobs_v
  INTEGER :: nobs_t
  INTEGER :: nobs_q
  INTEGER :: nobs_p

  INTEGER :: nobs_tmp(nobtype)
  INTEGER :: iobs(nobtype), iobs_out(nobtype), icount
  INTEGER :: itype, iobstype
  INTEGER(1) :: odata_flag(500000,nobtype)

  integer :: iflag

  REAL(4) :: odata(7,500000,nobtype)
  REAL(4) :: odata_tmp(7)
  REAL(4) :: odata_out(7,500000,nobtype)
 
  CHARACTER(6) :: obtypelist(nobtype)

  INTEGER,PARAMETER :: nz = 19
  INTEGER, PARAMETER :: omt=101
  INTEGER, PARAMETER :: omt1=102
  REAL(4) :: plev(nz)
  plev(1:nz) =(/1000.0, 925.0, 850.0, 700.0, &
       &  500.0, 400.0, 300.0, 250.0, 200.0, &
       &  150.0, 100.0,  70.0,  50.0,  40.0, &
       &   30.0,  20.0,  15.0,  10.0,   5.0/)

  obtypelist = (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',&
               & 'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG',&
               & 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND',&
               & 'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW'/)

  read(5,*) dist_limit, out_ps
  !OPEN(11,FILE='prepbufr.in',FORM='unformatted')

  sx =    0.0 + 360.0/float(nx)*0.5
  sy =  -90.0 + 180.0/float(ny)*0.5
  dx = 360.0/float(nx)
  dy = 180.0/float(ny)

  !write(*,*) sx, sy, dx, dy

  do i = 1, ny
    lat(i) = sy + dy * ( i-1 )
  end do
  do i = 1, nx
    lon(i) = sx + dx * ( i-1 )
  end do

  iobs(:) = 0
  iobs_out(:) = 0
  do
    read(11,END=999) wk
    itype=nint(wk(7))
    iobs(itype)=iobs(itype)+1
    odata(:,iobs(itype),itype)=wk(:)
  end do


999 CONTINUE

  ! ADPUPA   | A48102 | UPPER-AIR (RAOB, PIBAL, RECCO, DROPS) REPORTS  
  itype=1
  do i = 1, iobs(itype)
    do k = 1, nz
      if( abs(odata(4,i,itype)-plev(k) ) <= 1.e-2 ) then 
        write(omt) odata(:,i,itype)
        write(omt1,'(2i6,3f12.5)') itype, i, odata(2:3,i,itype), odata(5,i,itype)
        iobs_out(itype)=iobs_out(itype)+1
      end if
    end do
  end do

  ! AIRCFT   | A48104 | AIREP/PIREP,AMDAR(ASDAR/ACARS), E-ADAS(AMDAR BUFR) ACFT 
  itype=3
  do i = 1, iobs(itype)
    write(omt) odata(:,i,itype)
    write(omt1,'(2i6,3f12.5)') itype, i, odata(2:3,i,itype), odata(5,i,itype)
    iobs_out(itype)=iobs_out(itype)+1
  end do

  ! SATWND   | A48105 | SATELLITE-DERIVED WIND REPORTS
  itype=4
  do iy = 1, ny
    do ix = 1, nx
      icount=0
      dist_min=9999999999.0
      do i = 1, iobs(itype)
        if( odata(3,i,itype) > lat(iy)-0.5 .and. odata(3,i,itype) < lat(iy)+0.5 ) then
          if( odata(2,i,itype) > lon(ix)-0.5 .and. odata(2,i,itype) < lon(ix)+0.5 ) then
            if( odata_flag(i,itype) == 0 ) then
              call com_distll_1(lon(ix),lat(iy),odata(2,i,itype),odata(3,i,itype),dist)
              if(dist<dist_min) then
                icount=1
                dist_min=dist
                odata_tmp(:)=odata(:,i,itype)
                tmp_i=i
                tmp_n=itype
              end if
            end if
          end if
        end if
      end do
      if(icount/=0 .and. dist_min < dist_limit ) then
        write(omt) odata_tmp(1:7)
        write(omt1,'(2i6,3f12.5)') itype, i, odata_tmp(2:3), odata_tmp(5)
        iobs_out(itype)=iobs_out(itype)+1
      end if
    end do
  end do

  ! PROFLR   | A48106 | WIND PROFILER REPORTS
  itype=5
  do i = 1, iobs(itype)
    do k = 1, nz
      if( abs(odata(4,i,itype)-plev(k) ) <= 10.0 ) then 
        write(omt) odata(:,i,itype)
        write(omt1,'(2i6,3f12.5)') itype, i, odata(2:3,i,itype), odata(5,i,itype)
        iobs_out(itype)=iobs_out(itype)+1
      end if
    end do
  end do

  ! VADWND   | A48107 | VAD (NEXRAD) WIND REPORTS 
  itype=6
  do i = 1, iobs(itype)
    do k = 1, nz
      if( abs(odata(4,i,itype)-plev(k) ) <= 10.0 ) then
        write(omt) odata(:,i,itype)
        write(omt1,'(2i6,3f12.5)') itype, i, odata(2:3,i,itype), odata(5,i,itype)
        iobs_out(itype)=iobs_out(itype)+1
      end if
    end do
  end do

  ! ADPSFC   | A48109 | SURFACE LAND (SYNOPTIC, METAR) REPORTS
  itype=8
  if( out_ps ) then
  do iy = 1, ny
    do ix = 1, nx
      icount=0
      dist_min=9999999999.0
      do i = 1, iobs(itype)
        if( odata(3,i,itype) > lat(iy)-0.5 .and. odata(3,i,itype) < lat(iy)+0.5 ) then
          if( odata(2,i,itype) > lon(ix)-0.5 .and. odata(2,i,itype) < lon(ix)+0.5 ) then
            if( odata_flag(i,itype) == 0 ) then
              call com_distll_1(lon(ix),lat(iy),odata(2,i,itype),odata(3,i,itype),dist)
              if(dist<dist_min) then
                icount=1
                dist_min=dist
                odata_tmp(:)=odata(:,i,itype)
                tmp_i=i
                tmp_n=itype
              end if
            end if
          end if
        end if
      end do
      if(icount/=0 .and. dist_min < dist_limit ) then
        !write(*,'(2i6,7f12.4)') itype, tmp_i, odata_tmp(1:7)
        write(omt) odata_tmp(1:7)
        write(omt1,'(2i6,3f12.5)') itype, i, odata_tmp(2:3), odata_tmp(5)
        iobs_out(itype)=iobs_out(itype)+1
      end if
    end do
  end do
  end if

  ! SFCSHP   | A48110 | SURFACE MARINE (SHIP, BUOY, C-MAN PLATFORM) REPORTS
  itype=9
  if( out_ps ) then
  do i = 1, iobs(itype)
    !write(*,'(2i6,7f10.4)') itype, i, odata(:,i,itype)
    write(omt) odata(:,i,itype)
    write(omt1,'(2i5,3f12.5)') itype, i, odata(2:3,i,itype), odata(5,i,itype)
    iobs_out(itype)=iobs_out(itype)+1
  end do
  end if

  ! ASCATW   | A48121 | ASCAT SCATTEROMETER DATA (REPROCESSED) 
  itype=20
  do iy = 1, ny
    do ix = 1, nx
      icount=0
      dist_min=9999999999.0
      do i = 1, iobs(itype)
        if( odata(3,i,itype) > lat(iy)-0.5 .and. odata(3,i,itype) < lat(iy)+0.5 ) then
          if( odata(2,i,itype) > lon(ix)-0.5 .and. odata(2,i,itype) < lon(ix)+0.5 ) then
            if( odata_flag(i,itype) == 0 ) then
              call com_distll_1(lon(ix),lat(iy),odata(2,i,itype),odata(3,i,itype),dist)
              if(dist<dist_min) then
                icount=1
                dist_min=dist
                odata_tmp(:)=odata(:,i,itype)
                tmp_i=i
                tmp_n=itype
              end if
            end if
          end if
        end if
      end do
      if(icount/=0 .and. dist_min < dist_limit ) then
        !write(*,'(2i6,7f12.4)') itype, tmp_i, odata_tmp(1:7)
        write(omt) odata_tmp(1:7)
        write(omt1,'(2i6,3f12.5)') itype, i, odata_tmp(2:3), odata_tmp(5)
        iobs_out(itype)=iobs_out(itype)+1
      end if
    end do
  end do


  write(*,'(A)') '================================================================================'
  write(*,'(A)') '                     SUPEROBBING THE NCEP PREPBUFR DATA'
  write(*,'(A)') '--------------------------------------------------------------------------------'
  write(*,'(A,i10)') ' TOTAL NUMBER OF READ-IN RECORDS:',SUM(iobs)
  write(*,'(A,i10)') ' TOTAL NUMBER OF WRITTEN RECORDS:',SUM(iobs_out)
  write(*,'(A)') '--------------------------------------------------------------------------------'
  write(*,'(10(2X,A))') obtypelist(1:10)
  write(*,'(10i8)') iobs(1:10)
  write(*,'(10i8)') iobs_out(1:10)
  write(*,'(A)') '--------------------------------------------------------------------------------'
  write(*,'(10(2X,A))') obtypelist(11:20)
  write(*,'(10i8)') iobs(11:20)
  write(*,'(10i8)') iobs_out(11:20)
  write(*,'(A)') '================================================================================'


END PROGRAM MAIN
  !-----------------------------------------------------------------------------
  !subroutine temporal_thinning

  !end subroutine temporal_thinning
  !-----------------------------------------------------------------------------
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
    ctime_second &
         = real( year*(60*60*24*30*12) &
               + month*(60*60*24*30)   &
               + day*(60*60*24)        &
               + hour*(60*60)          &
               + min*60                &
               + sec, kind=8           )

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
