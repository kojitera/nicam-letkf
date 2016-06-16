!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  module common_precip
!
!  module for precipitation assimilation
!
!  created  Jan. 2012, Guo-Yuan Lien, UMD
!  adopted to GFS-LETKF and modified, May 2013, Guo-Yuan Lien, UMD
!  modified, September 2013, Guo-Yuan Lien, UMD
!
!  function dinvnorm(p) modified from Ren-Raw Chen, 
!    Rutgers University in New Brunswick, New Jersey
!    http://home.online.no/~pjacklam/notes/invnorm/
!
!  modified, May 2015 by Shunji Kotsuki, RIKEN-AICS
!    to adopt NICAM-LETKF
!
!-------------------------------------------------------------------------------
!
!  subroutine read_ppcdf     (cdffile_m, cdffile_o, ppcdf_m, ppcdf_o, ppzero_m, ppzero_o)
!  subroutine read_ppmask    (maskfile, ppmask)
!  function   pptrans_normal (pp, ppcdf, ppzero)
!  function   pptrans_log    (pp)
!  subroutine pptrans_normal_mdzero_def (pp_ens, ppcdf, ppzero,           zero_mem, ym, sigma)
!  function   pptrans_normal_mdzero     (pp,     ppcdf, ppzero, ppzero_m, zero_mem, ym, sigma)
!  function   compact_tail   (pos_cdf)
!  function   dinvnorm       (p)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fnc_transgauss
  use mod_sfmt
  use common_obs_nicam, only : &
        id_rain_obs,           &
        id_rain_obs_pp,   id_rain_obs_uv,   id_rain_obs_tt, id_rain_obs_qq, &
        id_rain_obs_uvtt, id_rain_obs_uvqq, id_rain_obs_ttqq

  implicit none
  integer, parameter      :: r_size = 8
  real(r_size), parameter :: undef  = -0.99900d35 
  real(r_size), parameter :: gaussmin = -7.0d0
  real(r_size), parameter :: gaussmax =  7.0d0
  real(r_size), parameter :: pi = dacos(-1.0d0)

  !====> options follow GYL 
  INTEGER, SAVE      :: opt_pptrans      ! equal to be "transform"
                                         ! 0: no transformation
                                         ! 1: log transformation
                                         ! 2: Gaussian transformation with median zero rain
                                         ! 3: Gaussian transformation with modified median zero rai
  INTEGER, SAVE      :: opt_ppobserr     ! equal to be "opt_preciperr"
                                         ! 1: transformed obs. error
                                         ! 2: constant obs. error
                                    
  INTEGER, SAVE      :: ncdf             ! equal to be "ndivid"  
  REAL(r_size), SAVE :: ppzero_thres     ! equal to be "rain_min_rate"
  REAL(r_size), SAVE :: gausstail_thres  ! threshold of the gaussian distribution tail
  REAL(r_size), SAVE :: min_ppobserr     ! minimum observation error for gaussian transformation

!  real(r_size), SAVE :: mask_thres                 ! threshold of assimilation area wrt. the mask file

  !==> rain2trans.cnf
  !&rain2trans
  character(128) :: grddir=''
  character(128) :: gusdir=''
  character(128) :: anldir=''
  character(128) :: obsdir=''
  character(128) :: datdir=''
  character(128) :: yhxdir=''
  character(128) :: vecdir=''
  character(128) :: obs_cdfname=''
  character(128) :: gus_cdfname=''
  character(128) :: mem_cdfname=''
  character(128) :: rainname=''
  integer, save :: nbv                     ! ensemble member
  integer, save :: cdate                   ! current date
  integer, save :: transform               ! defenition of transformation
!   ### optionsfor transformation
!   #-1 :: no assimilation
!   # 0 :: without transformation
!   # 1 :: log transformation
!   # 2 :: gaussian transformation with medium zero precip.
!   # 3 :: gaussian transformation with modified medium zero precip.
  integer, save :: opt_preciperr
!   ### options for precipitation error
!   # 0 :: to assimilate  BOTH prepubufr + precipitation
!   # 1 :: to assimilate  ONLY             precipitation
!   # 2 :: to assimilate  ONLY prepbufr
  integer, save :: opt_qccorr
!   ### options for qc from correlation between observation and background
!     # 0 :: do NOT the qc
!     # 1 :: do     the qc
  integer, save :: opt_exzero
!   ### options for qc for exclude zero precipitation
!     # 0 :: do NOT the qc ( all precip.          )
!     # 1 :: do     the qc ( exclude observed zero precip. )
  integer, save :: opt_fixzero      
!   ###  options for fix def. of zero precip. for guess
!     # 0 :: do NOTHING
!     # 1 :: to fix zero precip. probability by modifying pp_thresh for the model
  integer, save :: opt_modbias
!   ###  options to modify the bias between y and hx
!     # 0 :: do NOTHING
!     # 1 :: to modify the bias between y and hx            
  integer, save :: opt_ppguess      
!   ###  options for exclude large zero precip. rate for guess
!     # 0 :: do NOT the qc 
!     # 1 :: do     the qc (exclude the grids have larger zero precip. rate)
  integer, save :: opt_ppreduce
!   ###  options to reduce precip. observation   
!     # 0 :: do nothing
!     # 1 :: do pp_reduction

  integer, save :: opt_variloc      ! 0: do nothing
                                    ! 1: variable localization
      integer, save  :: extrp_good       ! ==> good precipitation in extra-tropic
      integer, save  :: extrp_norm       ! ==> norm precipitation in extra-tropic
      integer, save  :: tropc_good       ! ==> good precipitation in tropic
      integer, save  :: tropc_norm       ! ==> norm precipitation in tropic 

  REAL(8), SAVE :: rain_min_rate                   ! Rain Defenition [mm/6h]
  REAL(8), SAVE :: base_obsv_rain                  ! Rain level [Pa]  :Following Guo-Yuan (2014)
  REAL(8), SAVE :: err_rain                        ! Rain R  rate of error variance
  REAL(8), SAVE :: err_min_rain                    ! Rain R  threshold: minimum error [mm/6h]
  REAL(8), SAVE :: rain_min_mem                    ! Rain QC threshold: minimum precip member rate
  REAL(8), SAVE :: corr_min_qc                     ! Rain QC threshold: minimum background-obs. correlation
  REAL(8), SAVE :: rain_alpha                      ! used in log transformation [mm/hr] : Y=ln( rain + alpha ) by Mahfouf et al. (2007)
  INTEGER, SAVE :: ndivid                          ! the number of diviation of CDF
  REAL(8), SAVE :: gausstail                       ! threshold of the gaussian distribution tail  
  REAL(8), SAVE :: mingausserr                     ! minimum observation error for gaussian transformation
  LOGICAL, SAVE :: prdhst                          ! production of Joint Histogram
  REAL(8), SAVE :: hstsiz                          ! cell size of the of the joint-histogram
  REAL(8), SAVE :: locrad                          ! localization radius [m], should be equal to sigma_obs in letkf_obs.f90
  REAL(8), SAVE :: thres_ppguess                   ! threshold of the zero precip. rate (minimum)
  INTEGER, SAVE :: ppred_nph, ppred_trp, ppred_sph ! Reduce of Obs. (every XXX grids)

  ! only for cdfmake
  REAL(8), SAVE :: stp_fixzero                     ! delta step to fix zero precip. probability
  REAL(8), SAVE :: alw_fixzero                     ! allowable range to fix zero precip. probability

  ! random number
  INTEGER, SAVE :: seed
  !--> random number SFMT
!  REAL(r_size) ::  genrand_res53

   ! assumuming not used

!!!  real(r_size), parameter :: mask_thres = 0.35d0   ! threshold of assimilation area wrt. the mask file

!!!!-------------------------------------------------------------------------------
!!!
!!!  integer, parameter :: pp_bg_nlev = 2
!!!  integer, parameter :: pp_bg_levs(pp_bg_nlev-1) = &
!!!                        (/24/)
!!!  integer, parameter :: pp_ob_nlev = 2
!!!  real(r_size), parameter :: pp_ob_levs(pp_ob_nlev-1) = &
!!!                        (/ppzero_thres/)
!!!  logical, parameter :: pp_criterion(pp_bg_nlev,pp_ob_nlev) = reshape((/ &
!!!!           bg1   , bg2
!!!           .false.,.true., &  ! ob1
!!!           .false.,.true.  &  ! ob2
!!!           /), (/pp_bg_nlev,pp_ob_nlev/))

!  integer, parameter :: pp_bg_nlev = 4
!  integer, parameter :: pp_bg_levs(pp_bg_nlev-1) = &
!                        (/20,24,28/)
!  integer, parameter :: pp_ob_nlev = 3
!  real(r_size), parameter :: pp_ob_levs(pp_ob_nlev-1) = &
!                        (/ppzero_thres,1.d0/)
!  logical, parameter :: pp_criterion(pp_bg_nlev,pp_ob_nlev) = reshape((/ &
!!           bg1   , bg2,  , bg3   , bg4
!           .false.,.false.,.false.,.true., &  ! ob1
!           .false.,.false.,.true. ,.true., &  ! ob2
!           .false.,.true. ,.true. ,.true.  &  ! ob3
!           /), (/pp_bg_nlev,pp_ob_nlev/))

!-------------------------------------------------------------------------------

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_rain2trans()
  implicit none

  namelist / rain2trans /      &
    grddir, gusdir, anldir,    &
    obsdir, datdir, yhxdir,    &
    vecdir,                    &
    obs_cdfname, gus_cdfname,  &
    mem_cdfname, rainname,     &
    nbv, cdate,                &
    transform,                 & 
    opt_preciperr, opt_qccorr, &
    opt_exzero,  opt_fixzero,  &
    opt_modbias, opt_ppguess,  &
    opt_variloc, opt_ppreduce, &
    base_obsv_rain,            &
    err_rain, err_min_rain,    &
    rain_min_rate,             &
    rain_min_mem,              &
    corr_min_qc,               &
    rain_alpha,                &
    ndivid,                    &
    gausstail, mingausserr,    &
    prdhst, hstsiz, locrad,    &
    thres_ppguess,             &
    extrp_good, extrp_norm,    &
    tropc_good, tropc_norm,    &
    ppred_nph,  ppred_trp,     &
    ppred_sph
    
  open(1,file='rain2trans.cnf')
  read(1,nml=rain2trans)
  close(1)

end subroutine get_rain2trans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_transform(myrank)
  implicit none
  integer, intent (in) :: myrank
  
  if( transform == -1 ) then
    if( myrank==0 ) write(6,*) '     No Precip. Assimilation   : transform =',transform
  else if( transform == 0 ) then
    if( myrank==0 ) write(6,*) '     No Transformation         : transform =',transform  
  else if( transform == 1 ) then
    if( myrank==0 ) write(6,*) '     Log Transformation        : transform =',transform 
  else if( transform == 2 ) then
    if( myrank==0 ) write(6,*) '     Gaussian Transformation NM: transform =',transform 
  else if( transform == 3 ) then
    if( myrank==0 ) write(6,*) '     Gaussian Transformation MD: transform =',transform 
  else
    if( myrank==0 ) write(6,*) '     Error in setting of transform',transform 
    stop
  end if  
end subroutine check_transform  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_opt_ppobserr(myrank)
  implicit none
  integer, intent (in) :: myrank
  
  if( opt_ppobserr == 1 ) then
    if( myrank==0 ) write(6,*) '     Transformed Obs. R        : opt_ppobserr =',opt_ppobserr
  else if( opt_ppobserr  == 2 ) then
    if( myrank==0 ) write(6,*) '     Constant Obs. R           : opt_ppobserr =',opt_ppobserr
  else
    if( myrank==0 ) write(6,*) '     Error in setting of opt_ppobserr',opt_ppobserr
    stop
  end if  
end subroutine check_opt_ppobserr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine monit_variloc (id_tmp,raintype)
  implicit none
  integer, intent(in) :: id_tmp
  character(16), intent(in) :: raintype

  SELECT CASE( id_tmp )
    CASE( id_rain_obs )
      print *, '     UPDATE ALL    <<<===============   ',trim(raintype)
    CASE( id_rain_obs_pp )
      print *, '     UPDATE Rain   <<<===============   ',trim(raintype)
    CASE( id_rain_obs_uv )
      print *, '     UPDATE U&V    <<<===============   ',trim(raintype)
    CASE( id_rain_obs_tt )
      print *, '     UPDATE  T     <<<===============   ',trim(raintype)
    CASE( id_rain_obs_qq )
      print *, '     UPDATE  Q     <<<===============   ',trim(raintype)
    CASE( id_rain_obs_uvtt )
      print *, '     UPDATE UVT    <<<===============   ',trim(raintype)
    CASE( id_rain_obs_uvqq )
      print *, '     UPDATE UVQ    <<<===============   ',trim(raintype)
    CASE( id_rain_obs_ttqq )
      print *, '     UPDATE T&Q    <<<===============   ',trim(raintype)
    CASE default ! added by sawada
      print *, '     NO UPDATES    <<<===============   ',trim(raintype)
  END SELECT

end subroutine monit_variloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_varloc( id_tmp , qc, elm )
  implicit none
  integer, intent(in) :: id_tmp
  real, intent(inout) :: qc, elm

  SELECT CASE( id_tmp )
    CASE( id_rain_obs )
      elm = real( id_rain_obs    )
    CASE( id_rain_obs_pp )
      elm = real( id_rain_obs_pp ) 
    CASE( id_rain_obs_uv )
      elm = real( id_rain_obs_uv ) 
    CASE( id_rain_obs_tt )
      elm = real( id_rain_obs_tt )
    CASE( id_rain_obs_qq )
      elm = real( id_rain_obs_qq )
    CASE( id_rain_obs_uvtt )
      elm = real( id_rain_obs_uvtt )
    CASE( id_rain_obs_uvqq )
      elm = real( id_rain_obs_uvqq )
    CASE( id_rain_obs_ttqq )
      elm = real( id_rain_obs_ttqq )
    CASE default
      elm = real ( undef )
      qc  = 0.0
  END SELECT



end subroutine get_varloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_gaussvariable()
  implicit none
!  REAL(8), intent(in) :: rain_min_rate, gausstail
!  INTEGER, intent(in) :: ndivid, transform

  opt_pptrans     = transform
  ncdf            = ndivid
  ppzero_thres    = rain_min_rate
  gausstail_thres = gausstail
  min_ppobserr    = mingausserr 
  opt_ppobserr    = opt_preciperr
!  write(6,'(2i,4f)') ncdf,opt_pptrans,ppzero_thres,gausstail_thres,min_ppobserr
!  stop

end subroutine get_gaussvariable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_proc_l(lall,nprocs,nread,proc_l)
  implicit none
  integer, intent (in)  :: lall, nprocs, nread
  integer, intent (out) :: proc_l(0:nprocs-1,nread)
  integer iprocs, iread, l 

  proc_l(0:nprocs-1,1:nread) = -999

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
end subroutine set_proc_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_ppcdf(fname,gall,ndivid,ppcdf,ppzero,ppsamp,ppmask,dfzero)
  implicit none
  integer g, idiv
  integer, intent(in)        :: gall, ndivid
  character(255), intent(in) :: fname
  REAL(8), intent(out) :: ppcdf(gall,0:ndivid)
  REAL(8), intent(out) :: ppzero(gall), ppsamp(gall), ppmask(gall), dfzero(gall)
  REAL(8) :: ppcdfr(gall,0:ndivid)

  ppcdf(:,:) = undef
  ppzero(:)  = undef
  ppsamp(:)  = undef
  ppmask(:)  = undef 
  
  open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*8,action='read')
  do idiv=0,ndivid
    read(1,rec=idiv+1)          ppcdf(:,idiv)
  end do
  read(1,rec=2*(ndivid+1)+1)  ppsamp(:)
  read(1,rec=2*(ndivid+1)+2)  ppzero(:)
  read(1,rec=2*(ndivid+1)+3)  ppmask(:)
  read(1,rec=2*(ndivid+1)+5)  dfzero(:)
  close(1)

end subroutine read_ppcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_corr(fname,gall,ndivid,corr)
  implicit none
  integer g, idiv
  integer, intent(in)        :: gall, ndivid
  character(255), intent(in) :: fname
  REAL(8), intent(out) :: corr(gall)

  corr(:)  = undef 
  
  open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*8,action='read')
  read(1,rec=2*(ndivid+1)+4)  corr(:)
  close(1)

end subroutine read_corr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pptrans_normal_prate (pp, ppcdf, ppzero, dfzero, prate)
  implicit none
  real(r_size) :: pptrans_normal_prate
  real(r_size), intent(in) :: pp
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero, dfzero, prate
  real(r_size) :: pos_cdf, rr
  integer :: b
  real(r_size) :: rand


!!!  call random_number(rand)
!!!  print *, rand
!!!  stop

  call check_input_cdf( ppcdf(0), ppzero )

!  if (pp < ppzero_thres) then
  if (pp < dfzero ) then ! kotsuki 20150604
    pos_cdf = ppzero * prate
  else ! [pp >= ppzero_thres]
    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
        if( b==(ncdf+1) ) pos_cdf = 1.0d0
      end do
    end if
  end if

  pptrans_normal_prate = dinvnorm(compact_tail(pos_cdf)) ! get gaussian sigma from CDF

end function pptrans_normal_prate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pptrans_normal (pp, ppcdf, ppzero, dfzero)
  implicit none
  real(r_size) :: pptrans_normal
  real(r_size), intent(in) :: pp
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero, dfzero
  real(r_size) :: pos_cdf, rr
  integer :: b
  real(r_size) :: rand


!!!  call random_number(rand)
!!!  print *, rand
!!!  stop

  call check_input_cdf( ppcdf(0), ppzero )

!  if (pp < ppzero_thres) then
  if (pp < dfzero ) then ! kotsuki 20150604
    pos_cdf = ppzero * 0.5d0
  else ! [pp >= ppzero_thres]
    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
        if( b==(ncdf+1) ) pos_cdf = 1.0d0
      end do
    end if
  end if

  pptrans_normal = dinvnorm(compact_tail(pos_cdf)) ! get gaussian sigma from CDF

end function pptrans_normal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ppinverse_normal (sigma, ppcdf, ppzero)
  implicit none
  real(r_size) :: ppinverse_normal
  real(r_size), intent(in) :: sigma ! gaussian sigma
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size), intent(in) :: ppzero
  real(r_size) :: pos_cdf, rr, pp
  integer :: b
  
  call check_input_cdf( ppcdf(0), ppzero )

  if ( gaussmin.le.sigma .and. sigma.le.gaussmax ) then ! inverse transformation
    pp = sigma2cdf(sigma)
    if( pp .le. ppcdf(0) ) then
      ppinverse_normal = 0.0d0 ! no information
    else
      do b = 0,ncdf
        if (pp < dble(b)/dble(ncdf) ) then
          if ( b.ne.0 ) then
            rr = dble(ncdf) * ( dble(b)/dble(ncdf) - pp )
            ppinverse_normal = rr*ppcdf(b-1) + (1.0d0-rr)*ppcdf(b)
          else
            ppinverse_normal = ppcdf(b)
          endif
          exit
        endif
        if( b==ncdf ) ppinverse_normal = ppcdf(b)  
      end do
    endif  
  else
    ppinverse_normal = undef ! not defined
  end if

end function ppinverse_normal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ppnormal_obserror (obs_rain, obs_err_org, obs_trns, ppcdf, ppzero, dfzero)
  implicit none
  real(r_size) :: ppnormal_obserror
  real(r_size), intent(in) :: obs_rain, obs_err_org, obs_trns, ppzero, dfzero
  real(r_size), intent(in) :: ppcdf(0:ncdf)
  real(r_size) :: obs_err_p, obs_err_n

!  obs_err_p = - obs_trns + pptrans_normal( (obs_rain+obs_err_org), ppcdf(:), ppzero )            
!  obs_err_n =   obs_trns - pptrans_normal( (obs_rain-obs_err_org), ppcdf(:), ppzero )
  obs_err_p = - obs_trns + pptrans_normal( (obs_rain*(1.0d0+err_rain)), ppcdf(:), ppzero, dfzero )            
  obs_err_n =   obs_trns - pptrans_normal( (obs_rain*(1.0d0-err_rain)), ppcdf(:), ppzero, dfzero )


  if( obs_err_p < min_ppobserr ) obs_err_p = min_ppobserr
  if( obs_err_n < min_ppobserr ) obs_err_n = min_ppobserr
  ppnormal_obserror = 0.5d0 * ( obs_err_p + obs_err_n )

!  write(6,'(12f)') sngl(obs_rain),sngl(obs_err_org),sngl(ppnormal_obserror),&
!    sngl(ppzero), sngl(obs_trns),  &
!    sngl(pptrans_normal( obs_rain, ppcdf(:), ppzero )), &
!    sngl( pptrans_normal( (obs_rain+obs_err_org), ppcdf(:), ppzero ) ), &
!    sngl( pptrans_normal( (obs_rain-obs_err_org), ppcdf(:), ppzero ) ), &
!    sngl(pptrans_normal( (obs_rain+obs_err_org), ppcdf(:), ppzero )), &
!    sngl(obs_err_p), sngl(obs_err_n)

end function ppnormal_obserror
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pptrans_normal_mdzero_def (pp_ens, ppcdf, ppzero, dfzero, zero_mem, ym, sigma)
  implicit none
  real(r_size), intent(inout) :: pp_ens(nbv)
  real(r_size), intent(in)    :: ppcdf(0:ncdf)
  real(r_size), intent(in)    :: ppzero, dfzero
  
  integer, intent(out)        :: zero_mem
  real(r_size), intent(out)   :: ym
  real(r_size), intent(out)   :: sigma

  real(r_size) :: pos_cdf
  real(r_size) :: ppzero_b, pprain_b
  real(r_size) :: y_trace, y_trace_b
  real(r_size) :: alpha, beta
  logical :: zero(nbv)
  integer :: n
!
!!------
!!  real(r_size) :: pp_ens_ori(nbv)
!!  pp_ens_ori = pp_ens
!!------
!
  ym    = undef
  sigma = undef

  call check_input_cdf( ppcdf(0), ppzero )
  call check_opt_pptrans( opt_pptrans, 3 )

  beta = 0.0d0
  zero_mem = 0
  zero(:) = .false.
  do n = 1, nbv
!    if (pp_ens(n) < ppzero_thres) then
    if (pp_ens(n) < dfzero ) then ! kotsuki 20150604
      zero_mem = zero_mem + 1
      zero(n) = .true.
    else ! to apply gaussian normal transformation 
      pp_ens(n) = pptrans_normal(pp_ens(n), ppcdf, ppzero, dfzero)
      beta = beta + pp_ens(n)
    end if
  end do

  beta = beta / real(nbv, r_size)
  ppzero_b = real(zero_mem, r_size) / real(nbv, r_size)
  pprain_b = 1.0d0 - ppzero_b

  y_trace   = dinvnorm( compact_tail(ppzero)   )
  y_trace_b = dinvnorm( compact_tail(ppzero_b) )

  alpha = 0.0d0 - exp(0.0d0 - 0.5d0*y_trace_b*y_trace_b) / sqrt(2.0d0*pi)
  ym = (alpha * y_trace + beta * y_trace_b) / (alpha + pprain_b * y_trace_b)
  sigma = (pprain_b * y_trace - beta) / (alpha + pprain_b * y_trace_b)

  do n = 1, nbv
    if (zero(n)) then
      pos_cdf = ppzero_b * 0.5d0
      pp_ens(n) = ym + sigma * dinvnorm( compact_tail(pos_cdf) )
    end if
  end do

!  write(6,'(i,6f)') zero_mem, beta, y_trace, y_trace_b, alpha, ym, sigma

!
!!------
!!  print *, '----'
!!  print *, ppzero, ppzero_b, dinvnorm(pos_cdf)
!!  print *, y_trace, ym, sigma
!!  do n = 1, nbv
!!    print *, pp_ens_ori(n), pp_ens(n), zero(n)
!!  end do
!!------
!
end subroutine pptrans_normal_mdzero_def
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pptrans_normal_mdzero (pp, ppcdf, ppzero, dfzero, ppzero_m, zero_mem, ym, sigma)
  implicit none
  real(r_size) :: pptrans_normal_mdzero

  real(r_size), intent(in) :: pp, ppzero, dfzero, ppzero_m, ppcdf(0:ncdf)
  integer, intent(in)      :: zero_mem
  real(r_size), intent(in) :: ym, sigma

  real(r_size) :: pos_cdf, rr
  integer :: b

  call check_input_cdf( ppcdf(0), ppzero )
  call check_opt_pptrans( opt_pptrans, 3 )

!  if (pp < ppzero_thres) then
  if (pp < dfzero) then ! kotsuki 20150604
    pos_cdf = ppzero * 0.5d0
  else ! [pp >= ppzero_thres]
    if (pp < ppcdf(0)) then
      pos_cdf = 0.0d0
    else
      do b = 1, ncdf+1
        if (pp < ppcdf(b)) then
          rr = (pp - ppcdf(b-1)) / (ppcdf(b) - ppcdf(b-1))
          pos_cdf = ((1.0d0-rr) * real(b-1,r_size) + rr * real(b,r_size)) / real(ncdf,r_size)
          exit
        end if
        if( b==(ncdf+1) ) pos_cdf = 1.0d0
      end do
    end if
  end if

!!------
!!  print *, '---###'
!!  print *, ppzero, ppzero_m, zero_mem, ym, sigma
!!  print *, pos_cdf
!!------

  if (pos_cdf < ppzero_m) then
    pos_cdf = (pos_cdf / ppzero_m) * (real(zero_mem, r_size) / real(nbv, r_size))
    pptrans_normal_mdzero = ym + sigma * dinvnorm(compact_tail(pos_cdf))
  else
    pptrans_normal_mdzero = dinvnorm(compact_tail(pos_cdf))
  end if

!!------
!! print *, pos_cdf
!!  print *, pptrans_normal_mdzero
!!------

end function pptrans_normal_mdzero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_input_cdf( ppcdf_zero, ppzero )
  implicit none
  real(r_size), intent(in) :: ppcdf_zero, ppzero

  if (ppcdf_zero < -1.0d0 .or. ppzero < -1.0d0) then
    write (*, *) '[Error] Wrong input CDF.'
    stop
  end if

end subroutine check_input_cdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_opt_pptrans( opt_pptrans, opt_pptrans_correct )
  implicit none
  integer, intent(in) :: opt_pptrans, opt_pptrans_correct

  if (opt_pptrans /= opt_pptrans_correct) then
    write (*, *) '[Error] Unsupported transformation method.'
    stop
  end if
end subroutine check_opt_pptrans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine read_ppmask (maskfile, ppmask)
!
!  implicit none
!
!  character(len=*), intent(in) :: maskfile
!  real(r_size), intent(out) :: ppmask(nlon,nlat)
!
!  real(r_sngl) :: ppmask_s(nlon,nlat)
!  integer :: i, j, iolen
!  logical :: ex
!
!  inquire (iolength=iolen) iolen
!
!  inquire (file=trim(maskfile), exist=ex)
!  if (ex) then
!    open (92, file=trim(maskfile), status='old', form='unformatted', &
!              access='direct', recl=iolen*nlon*nlat)
!    read (92, rec=1) ((ppmask_s(i,j), i=1,nlon), j=1,nlat)
!    close (92)
!    ppmask = real(ppmask_s, r_size)
!  else
!    write (6,'(3A)') "Mask file ", maskfile, " does not exist -- skipped"
!    ppmask = 1.0e10  ! All data are used.
!  end if
!
!end subroutine read_ppmask
!


!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function sigma2cdf (sigma)
! http://www.snap-tck.com/room04/c01/stat/stat99/stat9901.html
  implicit none
  real(r_size), intent(in)  :: sigma
  real(r_size) :: sigma2cdf
  real(r_size), parameter :: epsiron = 1.0e-10
  real(r_size) :: ai, zet, szet, dev

  szet =  0.0d0
  ai   = -1.0d0
  dev  =  1.0d0
  do 
    ai  = ai + 2
    dev = dev*ai
    zet = (sigma**ai) / dev
    if( dabs(zet) < epsiron ) exit
    szet = szet + zet
  end do
  sigma2cdf = 0.5d0 + dexp( - 0.5d0*(sigma**2.0d0) ) * szet /dsqrt( 2.0d0*pi )

end function sigma2cdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compact_tail (pos_cdf)

  implicit none

  real(r_size) :: compact_tail
  real(r_size), intent(in) :: pos_cdf

  compact_tail = pos_cdf
  if (compact_tail < gausstail_thres        ) compact_tail = gausstail_thres
  if (compact_tail > 1.0d0 - gausstail_thres) compact_tail = 1.0d0 - gausstail_thres

end function compact_tail
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ren-raw chen, rutgers business school
! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real*8 function dinvnorm(p)
      real*8 p,p_low,p_high
      real*8 a1,a2,a3,a4,a5,a6
      real*8 b1,b2,b3,b4,b5
      real*8 c1,c2,c3,c4,c5,c6
      real*8 d1,d2,d3,d4
      real*8 z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
end function dinvnorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


end module fnc_transgauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
