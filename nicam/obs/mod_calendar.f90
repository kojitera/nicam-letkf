!-------------------------------------------------------------------------------
!
!+  calendar module
!
!-------------------------------------------------------------------------------
module mod_calendar
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines w.r.t. the calender.
  !
  !++ Orginal code : CCSR/NIES/AGCM5.4.01
  !                  Modified by M.Satoh
  !       
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.01      04-02-18   Imported from the NICAM-subset model
  !      0.03      04-05-31   Sub[calenar_ss2yh] belongs to "public".
  !                07-03-23   Y.Niwa calendar_daymo, private --> public
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: calendar_ointvl
  public :: calendar_ss2cc
  public :: calendar_ssaft
  public :: calendar_secdy
  public :: calendar_ss2ds
  public :: calendar_ss2yd
  public :: calendar_ym2dd
  public :: calendar_ds2ss
  public :: calendar_dayyr
  public :: calendar_ss2yh
  public :: calendar_xx2ss
  public :: calendar_yh2ss
  public :: calendar_PERPR
  public :: calendar_dd2ym
  public :: calendar_daymo  ! 07/03/23 Y.Niwa [add]

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  integer, private :: isecdy, isechr, jy, jy4, jcent, jcent4
  integer, private :: idays0, idy, ileap, id, m, idayyr, jyear, jmonth
  !
  !--- flag of automatic or not
  logical, private, save ::  oauto = .true.
  !<---                      yr=0-999     : 360day
  !<---                      yr=1000-1899 : 365day
  !<---                      yr=1900-     : gregorian
  !
  !--- flag of the gregorian calendar or not
  logical, private, save :: ogrego = .true.
  !
  !--- number of days in a month
  integer, private, save :: monday ( 12,2 )
  data   monday /                            &
       31,28,31,30,31,30,31,31,30,31,30,31,  &
       31,29,31,30,31,30,31,31,30,31,30,31 /
  !
  !--- flag of ideal calender (n day per month)
  logical, private, save :: oideal = .false.
  !------ 1 month = x days in the ideal case
  integer, private, save :: idaymo = 30
  !------ 1 year = x months in the ideal case
  integer, private, save :: imonyr = 12
  !
  !--- flag of perpetual or not 
  logical, private, save :: operpt = .false.
  !------ perpetual date(year)
  integer, private, save :: iyrpp  = 0
  !------ perpetual date(month)
  integer, private, save :: imonpp = 3
  !------ perpetual date(day)
  integer, private, save :: idaypp = 21
  !
  !--- 1 minute = x sec.
  integer, private, save :: isecmn = 60
  !--- 1 hour = x minutes
  integer, private, save :: iminhr = 60
  !--- 1 day = x hours
  integer, private, save :: ihrday = 24
  !
  logical, private, save :: ooperz = .true.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine calendar_perpt( &   
       iyear ,               & !--- in : year
       imonth,               & !--- in : month
       iday                  & !--- in : day
       )
    !--- fixed date
    implicit none

    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday

    operpt = .true.
    iyrpp  = iyear
    imonpp = imonth
    idaypp = iday

    return
  end subroutine calendar_perpt

  !-----------------------------------------------------------------------------
  subroutine calendar_PERPR( &    !! calendar, refer to fixed date
       IYEAR ,               & !--- OUT
       IMONTH,               & !--- OUT
       IDAY,                 & !--- OUT
       OOPERP                & !--- OUT
       )
    implicit none

    integer, intent(out) :: IYEAR
    integer, intent(out) :: IMONTH
    integer, intent(out) :: IDAY
    logical, intent(out) :: OOPERP

    OOPERP = OPERPT
    IYEAR  = IYRPP
    IMONTH = IMONPP
    IDAY   = IDAYPP

    RETURN
  end subroutine calendar_PERPR

  !-----------------------------------------------------------------------------
  subroutine calendar_perpo( &
       ooperp                & !--- in : flag of perpetual
       )
    implicit  none

    !--- fixed date(on/off)
    logical, intent(in) :: ooperp
    ooperz = ooperp

    return
  end subroutine calendar_perpo

  !-----------------------------------------------------------------------------
  subroutine calendar_daymo(&
       ndaymo,              & !--- OUT : day
       iyear,               & !--- IN : year
       imonth               & !--- IN : month
       )
    !--- calendar, no.of day in a month
    implicit none

    integer, intent(out) :: ndaymo
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndaymo = monday( imonth,2 )
       else
          ndaymo = monday( imonth,1 )            
       endif
    else if ( .not. oideal ) then
       ndaymo = monday( imonth,1 )            
    else
       ndaymo = idaymo
    endif

    return
  end subroutine calendar_daymo

  !-----------------------------------------------------------------------------
  subroutine calendar_dayyr(& 
       ndayyr,              & !-- OUT : day
       iyear                & !---IN : year
       )
    !--- calendar, no.of day in an year
    implicit none

    integer, intent(out) :: ndayyr
    integer, intent(in) :: iyear

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndayyr = 366
       else
          ndayyr = 365
       endif
    else if ( .not. oideal ) then
       ndayyr = 365
    else
       ndayyr = idaymo*imonyr
    endif

    return
  end subroutine calendar_dayyr

  !-----------------------------------------------------------------------------
  subroutine calendar_monyr( &   
       nmonyr,               & !--- OUT : month
       iyear                 & !--- IN : year
       )
    !--- calendar, no.of month in an year
    implicit none

    integer, intent(out) :: nmonyr
    integer, intent(in) :: iyear

    if ( oauto  ) then
       nmonyr = 12
    else
       nmonyr = imonyr
    endif

    return
  end subroutine calendar_monyr

  !-----------------------------------------------------------------------------
  subroutine calendar_secdy( &
       nsecdy                & !--- OUT : sec
       )
    !--- calendar, no.of sec. in a day
    implicit none

    integer, intent(out) :: nsecdy
    nsecdy = isecmn*iminhr*ihrday

    return
  end subroutine calendar_secdy

  !-----------------------------------------------------------------------------
  subroutine calendar_secmi( &
       nsecmi                & !--- OUT
       )
    !--- calendar, no of sec. in a minute
    implicit none

    integer, intent(out) :: nsecmi
    nsecmi = isecmn

    return
  end subroutine calendar_secmi

  !-----------------------------------------------------------------------------
  subroutine calendar_sechr( &   
       nsechr                &
       )
    !--- calendar, no.of sec. in an hour
    implicit none

    integer, intent(out) :: nsechr
    nsechr = isecmn*iminhr

    return
  end subroutine calendar_sechr

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2ds( &
       idays ,               & !--- OUT
       rsec  ,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> ddss
    implicit none

    integer, intent(out) :: idays
    real(8), intent(out) :: rsec
    real(8), intent(in) :: dsec

    isecdy = isecmn*iminhr*ihrday
    idays  = int( dsec/dble (isecdy) ) + 1
    rsec   = dsec - dble (idays-1)*dble (isecdy)
    if ( nint( rsec ) >= isecdy ) then
       idays = idays + 1
       rsec  = rsec - dble(isecdy)
    endif

    return
  end subroutine calendar_ss2ds

  !-----------------------------------------------------------------------------
  subroutine calendar_ds2ss( &
       dsec  ,               & !--- OUT
       idays ,               & !--- IN
       rsec                  & !--- IN
       )
    !--- calendar, ddss -> sec.
    implicit none

    real(8), intent(out) :: dsec
    integer, intent(in) :: idays
    real(8), intent(in) :: rsec

    isecdy = isecmn*iminhr*ihrday
    dsec   = dble (idays-1)*dble (isecdy) + dble (rsec)

    return
  end subroutine calendar_ds2ss

  !-----------------------------------------------------------------------------
  subroutine calendar_rs2hm( &
       ihour ,               & !--- OUT
       imin  ,               & !--- OUT
       isec  ,               & !--- OUT
       rsec                  & !--- IN
       )
    !--- calendar, sec. -> hhmmss
    implicit none

    integer, intent(out) :: ihour
    integer, intent(out) :: imin
    integer, intent(out) :: isec
    real(8), intent(in) :: rsec

    isechr = isecmn*iminhr
    ihour  = int ( rsec / dble(isechr ) )
    imin   = int ( ( rsec - dble(ihour*isechr) )/dble(isecmn) )
    isec   = nint( rsec - dble(ihour*isechr) - dble(imin*isecmn) )

    if ( isec >= isecmn ) then
       imin  = imin + 1
       isec  = isec - isecmn
    endif

    if ( imin == iminhr ) then
       ihour = ihour + 1
       imin  = imin  - iminhr
    endif

    return
  end subroutine calendar_rs2hm

  !-----------------------------------------------------------------------------
  subroutine calendar_hm2rs( &
       rsec  ,               & !--- OUT
       ihour ,               & !--- IN
       imin  ,               & !--- IN
       isec                  & !--- IN
       )
    !--- calendar, hhmmss -> sec.
    implicit none

    real(8), intent(out) :: rsec
    integer, intent(in) :: ihour
    integer, intent(in) :: imin
    integer, intent(in) :: isec

    rsec = real( ihour*isecmn*iminhr + imin*isecmn + isec, kind=8 )

    return
  end subroutine calendar_hm2rs

  !-----------------------------------------------------------------------------
  subroutine calendar_dd2ym( &
       iyear ,               & !--- OUT
       imonth,               & !--- OUT
       iday  ,               & !--- OUT
       idays                 & !--- IN
       )
    !--- calendar, day -> yymmdd
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    integer, intent(in)  :: idays

    if ( oauto  ) then
       if ( idays >= 693961 ) then       !" 1900*365+1900/4-19+5
          ogrego = .true.
       else
          ogrego = .false.
          if ( idays >= 1000*365 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( operpt .and. ooperz ) then
       iyear  = iyrpp
       imonth = imonpp
       iday   = idaypp
          return
       endif
    if ( ogrego ) then
       jy     = int(dble(idays)/365.24d0)
1100   continue 
       jy4    = (jy+3)/4
       jcent  = (jy+99)/100
       jcent4 = (jy+399)/400
       idays0 = jy*365+jy4-jcent+jcent4
       if ( idays <= idays0 ) then
          jy = jy -1 
          if ( jy >= 0 ) goto 1100
       endif
       iyear = jy
       idy   = idays - idays0
       if ( ocleap( iyear ) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    else if ( .not. oideal ) then
       iyear = idays/365
       idy   = idays - iyear*365
       ileap = 1
    endif
    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, 12
          id = id + monday(m,ileap)
          if ( idy <= id ) then
             imonth = m
             iday   = idy-id+monday(m,ileap)
             goto 3190
          endif
       end do
3190   continue 
    else 
       idayyr = idaymo*imonyr
       iyear  = ( idays-1 ) / idayyr
       imonth = ( idays-1 - iyear*idayyr )/idaymo+1
       iday   = idays - iyear*idayyr - (imonth-1)*idaymo
    endif

    return
  end subroutine calendar_dd2ym

  !-----------------------------------------------------------------------------
  subroutine calendar_ym2dd( &
       idays ,               & !--- OUT
       iyear ,               & !--- IN
       imonth,               & !--- IN
       iday                  & !--- IN
       )
    !--- calendar, yymmdd -> day
    implicit none
    integer, intent(out) :: idays
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif
    if ( ogrego .or. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       jy4    = (jyear+3)/4
       jcent  = (jyear+99)/100
       jcent4 = (jyear+399)/400
       idays0 = jyear*365+jy4-jcent+jcent4
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       endif
    else if ( .not. oideal ) then
       idays0 = jyear*365
       ileap  = 1
    endif

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       end do
    else
       idays0 = iyear*idaymo*imonyr
       id     = (imonth-1)*idaymo
    endif
    idays = idays0 + id + iday
    !      
    return
  end subroutine calendar_ym2dd
  !-----------------------------------------------------------------------------
  subroutine calendar_ym2yd( &
       idaysy,               & !--- OUT
       iyear ,               & !--- IN
       imonth,               & !--- IN
       iday                  & !--- IN
       )
    !--- calendar, yymmdd -> yydd
    implicit none
    integer, intent(out) :: idaysy
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif
    if ( ogrego .or. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       endif
    else if ( .not. oideal ) then
       ileap  = 1
    endif

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       end do
    else
       id     = (imonth-1)*idaymo
    endif
    idaysy = id + iday
    !      
    return
  end subroutine calendar_ym2yd
  !-----------------------------------------------------------------------------
  subroutine calendar_ss2yh( &
       idate ,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> date
    implicit none
    integer, intent(out) :: idate(6) !" yymmddhhmmss
    real(8), intent(in) :: dsec      !" time

    integer    idays             !" serial no.of day
    real(8)    rsec              !" no. of sec. in a day
    call calendar_ss2ds &
         ( idays , rsec  , &
         dsec            )
    call calendar_dd2ym &
         ( idate(1), idate(2), idate(3), &
         idays                        )
    call calendar_rs2hm &
         ( idate(4), idate(5), idate(6), &
         rsec                          )

    return
  end subroutine calendar_ss2yh

  !-----------------------------------------------------------------------------
  subroutine calendar_yh2ss( & 
       dsec  ,               & !--- OUT
       idate                 & !--- IN
       )
    !--- calendar, date -> sec.
    implicit none

    real(8), intent(out) :: dsec    !" time
    integer, intent(in) :: idate(6) !" yymmddhhmmss

    integer    idays             !" serial no.of day
    real(8)    rsec              !" no. of sec. in a day

    call calendar_ym2dd &
         ( idays   , &
         idate(1), idate(2), idate(3) )

    call calendar_hm2rs &
         ( rsec    , &
         idate(4), idate(5), idate(6) )

    call calendar_ds2ss &
         ( dsec  , &
         idays , rsec   )

    return
  end subroutine calendar_yh2ss

  !-----------------------------------------------------------------------------
  subroutine calendar_dd2yd( &
       iyear ,               & !--- OUT
       idaysy,               & !--- OUT
       idays                 & !--- IN
       )
    !--- calendar, day -> yydd
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: idaysy
    integer, intent(in) :: idays
    integer :: imonth, iday

    call calendar_dd2ym &
         ( iyear , imonth, iday  , &
         idays                   )

    call calendar_ym2yd &
         ( idaysy, &
         iyear , imonth, iday    )

    return
  end subroutine calendar_dd2yd

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2yd( &
       iyear ,               & !--- OUT
       idaysy,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> yydd
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: idaysy
    real(8), intent(in) :: dsec

    integer :: imonth, iday
    call calendar_ss2ym &
         ( iyear , imonth, iday  , &
         dsec                    )
    call calendar_ym2yd &
         ( idaysy, &
         iyear , imonth, iday    )

    return
  end subroutine calendar_ss2yd

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2ym( &
       iyear ,               & !--- OUT
       imonth,               & !--- OUT
       iday  ,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> yymmdd
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    real(8), intent(in) :: dsec

    integer :: idays
    real(8) :: rsec

    call calendar_ss2ds &
         ( idays , rsec  , &
         dsec            )
    call calendar_dd2ym &
         ( iyear , imonth, iday  , &
         idays                  )

    return
  end subroutine calendar_ss2ym

  !-----------------------------------------------------------------------------
  subroutine calendar_xx2ss( &
       ddsec ,               & !--- OUT
       rtdur ,               & !--- IN
       hunit,                & !--- IN
       dsec                 & !--- IN
       )
    !--- calendar, hour ->sec.
    implicit none

    real(8), intent(out) :: ddsec
    real(8), intent(in) :: rtdur
    character(len=*), intent(in) :: hunit
    real(8), intent(in) :: dsec
    character(len=10) :: hunitx
    integer :: isecmi, isechr, isecdy
    integer :: iyear, imonth, iday, ndaymo, ndayyr

    hunitx = hunit
    !c     call cupper( hunitx )
    if      ( hunitx(1:1) == 's' .or. hunitx(1:1) == 'S' ) then
       ddsec = dble (rtdur)
    else if ( hunitx(1:2) == 'mi' .or. hunitx(1:2) == 'MI' ) then
       call calendar_secmi( isecmi )
       ddsec = dble (rtdur)*dble (isecmi)
    else if ( hunitx(1:1) == 'h' .or. hunitx(1:1) == 'H' ) then
       call calendar_sechr( isechr )
       ddsec = dble (rtdur)*dble (isechr)
    else if ( hunitx(1:1) == 'd' .or. hunitx(1:1) == 'D' ) then
       call calendar_secdy( isecdy )
       ddsec = dble (rtdur)*dble (isecdy)
    else if ( hunitx(1:2) == 'mo' .or. hunitx(1:2) == 'MO' ) then
       call calendar_ss2ym &
            ( iyear , imonth, iday  , &
            dsec                    )
       call calendar_daymo &
            ( ndaymo, &
            iyear , imonth  )
       call calendar_secdy( isecdy )
       ddsec = dble (rtdur)*dble (ndaymo)*dble (isecdy)
    else if ( hunitx(1:1) == 'y' .or. hunitx(1:1) == 'Y' ) then
       call calendar_ss2ym &
            ( iyear , imonth, iday  , &
            dsec                    )
       call calendar_dayyr &
            ( ndayyr, &
            iyear   )
       call calendar_secdy( isecdy )
       ddsec = dble (rtdur)*dble (ndayyr)*dble (isecdy)
    else
       write (6,*) ' ### cxx2ss: invalid unit : ', hunit, &
            ' [sec] assumed'
       ddsec = rtdur
    endif

    return
  end subroutine calendar_xx2ss

  !-----------------------------------------------------------------------------
  subroutine calendar_cc2yh( &
       itime ,               & !--- OUT
       htime                 & !--- IN
       )
    !--- calendar, character -> date
    implicit none

    integer, intent(out) :: itime ( 6 )
    character(len=*), intent(in) :: htime

    integer :: i
    read ( htime, 2600 ) (itime(i),i=1,6)
2600 format( i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2 )

    return
  end subroutine calendar_cc2yh

  !-----------------------------------------------------------------------------
  subroutine calendar_yh2cc( &
       htime ,               & !--- OUT
       itime                 & !--- IN
       )
    !--- calendar, date -> character
    implicit none

    character(len=*), intent(out) :: htime
    integer, intent(in) :: itime ( 6 )

    integer :: i
    write ( htime, 600 ) (itime(i),i=1,6)
600 format( i4.4,'/',i2.2,'/',i2.2,'-',i2.2,':',i2.2,':',i2.2 )

    return
  end subroutine calendar_yh2cc

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2cc( &
       htime ,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> character
    implicit none

    character(len=*), intent(out) :: htime
    real(8), intent(in) :: dsec

    integer :: itime(6), i
    call calendar_ss2yh &
         ( itime , &
         dsec   )
    write ( htime, 600 ) (itime(i),i=1,6)
600 format( i4.4,'/',i2.2,'/',i2.2,'-',i2.2,':',i2.2,':',i2.2 )

    return
  end subroutine calendar_ss2cc

  !-----------------------------------------------------------------------------
  function ocleap( &
       iyear       & !--- IN
       )
    !--- calendar :bissextile or not
    implicit none

    logical :: ocleap
    integer, intent(in) :: iyear

    integer :: iy, iycen, icent

    iy     = mod(iyear,4)
    iycen  = mod(iyear,100)
    icent  = mod(iyear/100,4)
    if ( iy == 0 .and. ( iycen /= 0 .or. icent == 0 ) ) then
       ocleap = .true.
    else
       ocleap = .false.
    endif

    return
  end function ocleap

  !-----------------------------------------------------------------------------
  subroutine calendar_ssaft( &
       dseca ,               & !--- OUT
       dsec  ,               & !--- IN
       raftr ,               & !--- IN
       hunit                 & !--- IN
       )
    !--- calendar, time advancing
    implicit none

    real(8), intent(out) :: dseca
    real(8), intent(in) :: dsec
    real(8), intent(in) :: raftr
    character(len=*), intent(in) :: hunit

    integer :: idays, iyear, imonth, iday
    real(8) :: rsec
    real(8) :: ddtime

    if ( hunit(1:1) == 'y' .or. hunit(1:1) == 'Y' &
         .or. hunit(1:1) == 'mo' .or. hunit(1:2) == 'MO' ) then
       call calendar_ss2ds &
            ( idays , rsec  , &
            dsec            )
       call calendar_dd2ym &
            ( iyear , imonth, iday  , &
            idays                  )
       if ( hunit(1:1) == 'y' .or. hunit(1:1) == 'Y' ) then
          iyear  = iyear  + int(raftr)
       else if ( hunit(1:2) == 'mo' .or. hunit(1:2) == 'MO' ) then
          imonth = imonth + int(raftr)
       endif
       call calendar_ym2dd &
            ( idays , &
            iyear , imonth, iday   )
       call calendar_ds2ss &
            ( dseca , &
            idays , rsec    )
    else
       call calendar_xx2ss( ddtime, raftr, hunit, dsec )
       dseca = dsec + ddtime
    endif

    return
  end subroutine calendar_ssaft

  !-----------------------------------------------------------------------------
  function   calendar_ointvl( &
       dtime ,                & !--- IN
       dtprev ,               & !--- IN
       dtorgn,                & !--- IN
       rintv ,                & !--- IN
       htunit                 & !--- IN
       )
    !--- time step passed ?
    implicit none

    logical :: calendar_ointvl
    real(8), intent(in) :: dtime
    real(8), intent(in) :: dtprev
    real(8), intent(in) :: dtorgn
    real(8), intent(in) :: rintv
    character(len=*), intent(in) :: htunit

    real(8) :: ddtime
    character(len=5) :: hunit
    integer :: iyear, imon, iday, iyearp, imonp, idayp
    integer :: iy, imo
    integer :: nmonyr, ndayyr, ndaymo
    real(8) :: ry, rmo

    hunit = htunit

    if ( dtime == dtprev ) then
       calendar_ointvl = .true.         
       return
    endif

    calendar_ointvl = .false.

    call calendar_ss2ym ( iyear, imon, iday, dtime  )
    call calendar_ss2ym ( iyearp, imonp, idayp, dtprev )
    call calendar_xx2ss ( ddtime, rintv, hunit, dtime  )

    if ( dtime >= dtorgn ) then
       if      ( hunit(1:1) == 'y' .or. hunit(1:1) == 'Y' ) then
          call calendar_monyr( nmonyr, iyear )
          call calendar_dayyr( ndayyr, iyear )
          ry = real( iyear-iyearp, kind=8 )&
             + real( imon-imonp, kind=8 ) / real( nmonyr, kind=8 ) &
             + real( iday-idayp, kind=8 ) / real( ndayyr, kind=8 )
          if ( ry >= rintv ) then
             calendar_ointvl = .true.
          endif
       else if ( hunit(1:2) == 'mo' .or. hunit(1:2) == 'MO' ) then
          imo = 0
          do iy = iyearp, iyear-1
             call calendar_monyr( nmonyr, iy )
             imo = imo + nmonyr
          end do
          call calendar_daymo( ndaymo, iyear, imon )
          rmo = real( imon-imonp+imo, kind=8 ) &
              + real( iday-idayp, kind=8 ) / real( ndaymo, kind=8 )
          if ( rmo >= rintv ) then
             calendar_ointvl = .true.
          endif
       else if (      calendar_dgaus((dtime -dtorgn)/ddtime)  &
            > calendar_dgaus((dtprev-dtorgn)/ddtime) ) then
          calendar_ointvl = .true.
       endif
    endif

    return
  end function calendar_ointvl

  !-----------------------------------------------------------------------------
  function calendar_dgaus( &
        dx                 & !--- IN
        )
    !--- dicard gaussian
    implicit none
    real(8) :: calendar_dgaus
    real(8), intent(in) :: dx

    calendar_dgaus = aint(dx) + aint(dx - aint(dx) + 1.d0) - 1.d0

  end function calendar_dgaus

end module mod_calendar
!-------------------------------------------------------------------------------
