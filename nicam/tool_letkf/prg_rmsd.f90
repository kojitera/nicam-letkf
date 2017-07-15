PROGRAM main
  INTEGER            :: sdate, edate
  INTEGER            :: idate
  INTEGER            :: ndate
  CHARACTER(10)      :: cdate

  INTEGER, PARAMETER :: nx = 240
  INTEGER, PARAMETER :: ny = 120
  INTEGER, PARAMETER :: np = 26

  CHARACTER(256)     :: analysis_basename
  CHARACTER(256)     :: analysis_fname
  CHARACTER(256)     :: truth_basename
  CHARACTER(256)     :: truth_fname

  REAL(4) :: ref(nx,ny,np)
  REAL(4) :: anl(nx,ny,np)
  REAL(4) :: diff(nx,ny,np)
  REAL(4) :: rmsd(nx,ny,np)

  namelist / rmsdcnf / &
    nx, ny,            &
    sdate, edate,      &
    analysis_basename, &
    truth_basename

  open(1,file='rmsd.cnf')
  read(1,nml=rmsdcnf)
  close(1)

  idate=sdate
  WRITE(cdate,'(I10.10)') idate
  analysis_fname=TRIM(analysis_basename)//TRIM(cdate)//'.bin'
  truth_fname=TRIM(truth_basename)//TRIM(cdate)//'.dat'

  OPEN(2,FILE=TRIM(truth_fname),FORM='unformatted',ACCESS='direct',RECL=4*nx*ny)
  OPEN(3,FILE=TRIM(analysis_fname),FORM='unformatted',ACCESS='direct',RECL=4*nx*ny)

  READ(2,REC=1) ref(:,:,:)
  READ(3,REC=2) anl(:,:,:)


END PROGRAM main
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
      if( mon.ge.13 ) then
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

