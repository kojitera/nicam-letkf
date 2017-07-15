SUBROUTINE AMSUA_fwd(nlevels, nprof, coef_filename, tmp_p, tmp_t, tmp_q, tmp_t2m, tmp_q2m, &
                     tmp_p2m, tmp_u2m, tmp_v2m, tmp_soze, tmp_soaz, tmp_saze, tmp_saaz,       &
                     tmp_elev, tmp_lon, tmp_lat, tmp_land,                  &
                     bt_out, tran_out  ) 
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2010, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !     TEST PROGRAM FOR RTTOV FORWARD MODEL
  !          RTTOV VERSION 11
  !
  ! To run this program you must have the following files
  ! either resident in the same directory or set up as a
  ! symbolic link:
  !   the file containing input profiles (e.g. prof.dat)
  !   the RTTOV coefficient file
  !
  ! The script run_example_fwd.sh may be used to run this program.
  ! The output is generated in a file called example_fwd_output.dat.
  !
  !
  ! To adapt the code to their own applications users should
  ! edit the code highlighted like this:
  !     !================================
  !     !======Read =====start===========
  !          code to be modified
  !     !======Read ===== end ===========
  !     !================================
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date        Comment
  ! -------   ----        -------
  !  1.0    27/04/2004   orginal (based on tstrad) P. Brunel
  !  1.1    09/08/2004   modified to allow for variable no. channels/per profile
  !                       R. Saunders
  !  1.2    13/04/2007   Modified for RTTOV-90
  !  1.3    31/07/2007   Modified for RTTOV-91 R Saunders
  !  1.4    11/10/2007   Parallel version P.Marguinaud
  !  2.0    25/06/2010   Modified for RTTOV-10 J Hocking
  !  2.1    23/11/2011   Updates for v10.2 J Hocking
  !  3.0    23/11/2012   Updated for v11 J Hocking
  !
  ! Code Description:
  !   Language:          Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & platform_name,       &
       & inst_name

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
       & rttov_options,       &
       & rttov_coefs,         &
       & profile_type,        &
       & transmission_type,   &
       & radiance_type,       &
       & rttov_chanprof,      &
       & rttov_emissivity,    &
       & rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_parallel_direct.interface"
#include "rttov_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !INTEGER, INTENT(IN) :: nprof
  Real(Kind=jprb) :: tmp_p(nlevels, nprof)
  Real(Kind=jprb) :: tmp_t(nlevels, nprof)
  Real(Kind=jprb) :: tmp_q(nlevels, nprof)
  Real(Kind=jprb) :: tmp_t2m(1, nprof)
  Real(Kind=jprb) :: tmp_q2m(1, nprof)
  Real(Kind=jprb) :: tmp_p2m(1, nprof)
  Real(Kind=jprb) :: tmp_u2m(1, nprof)
  Real(Kind=jprb) :: tmp_v2m(1, nprof)
  Real(Kind=jprb) :: tmp_soze(1, nprof)
  Real(Kind=jprb) :: tmp_soaz(1, nprof)
  Real(Kind=jprb) :: tmp_saze(1, nprof)
  Real(Kind=jprb) :: tmp_saaz(1, nprof)
  Real(Kind=jprb) :: tmp_elev(1, nprof)
  Real(Kind=jprb) :: tmp_lon(1, nprof)
  Real(Kind=jprb) :: tmp_lat(1, nprof)
  Real(Kind=jprb) :: tmp_land(1, nprof)
  Real(Kind=8) :: bt_out(15,nprof)
  Real(Kind=8) :: tran_out(nlevels,15,nprof)

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output
  INTEGER(KIND=jpim), PARAMETER :: mxchn = 9000 ! max number of channels

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)                  :: opts           ! Options structure
  TYPE(rttov_coefs)                    :: coefs          ! Coefficients structure
  TYPE(rttov_chanprof),    ALLOCATABLE :: chanprof(:)    ! Input channel/profile list
  TYPE(profile_type),      ALLOCATABLE :: profiles(:)    ! Input profiles
  LOGICAL(KIND=jplm),      ALLOCATABLE :: calcemis(:)    ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  ALLOCATABLE :: emissivity(:)  ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      ALLOCATABLE :: calcrefl(:)    ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), ALLOCATABLE :: reflectance(:) ! Input/output surface BRDF
  TYPE(transmission_type)              :: transmission   ! Output transmittances
  TYPE(radiance_type)                  :: radiance       ! Output radiances

  INTEGER(KIND=jpim)                   :: errorstatus    ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status(10)
  CHARACTER(LEN=11)  :: NameOfRoutine = 'example_fwd'

  ! variables for input
  !====================
  INTEGER(KIND=jpim) :: input_chan(mxchn)
  REAL(KIND=jprb)    :: input_ems(mxchn), input_brdf(mxchn)
  CHARACTER(LEN=*) :: coef_filename
  !CHARACTER(LEN=256) :: coef_filename='/scratch/ra000015/koji/rttov/11.1/rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_15_amsua.dat'
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: ivch, ich
  REAL(KIND=jprb)    :: ems_val, brdf_val
  INTEGER(KIND=jpim) :: asw
  REAL(KIND=jprb), ALLOCATABLE :: emis(:), brdf(:)
  INTEGER(KIND=jpim), ALLOCATABLE :: nchan(:)
  REAL(KIND=jprb)    :: trans_out(10)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Set up the chanprof array with the channels/profiles to simulate
  !   4. Allocate RTTOV profile, radiance and transmittance structures
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  errorstatus     = 0_jpim
  alloc_status(:) = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============
  !WRITE(0,*) 'enter name of coefficient file (in current directory)'
  !READ(*,*) coef_filename
  !WRITE(0,*) 'enter name of file containing profile data (in current directory)'
  !READ(*,*) prof_filename
  !WRITE(0,*) 'enter number of profiles'
  !READ(*,*) nprof
  !WRITE(0,*) 'enter number of profile levels'
  !READ(*,*) nlevels
  !WRITE(0,*) 'turn on solar simulations? (0=no, 1=yes)'
  dosolar=0
  !READ(*,*) dosolar


  ! Read channel list including input surface emissivities and reflectances
  ! If the input emissivities/reflectances are zero or less, calcemis/calcrefl is set to true.

!  nchannels = 0_jpim
!  READ(*,*,iostat=ios) ich, ivch, ems_val, brdf_val ! channel number, validity, emissivity, reflectance
!  DO WHILE (ios == 0 )
!    IF (ivch /= 0) THEN
!      nchannels = nchannels + 1_jpim
!      input_chan(nchannels) = ich
!      input_ems(nchannels)  = ems_val
!      input_brdf(nchannels) = brdf_val
!    ENDIF
!    READ(*,*,iostat=ios) ich, ivch, ems_val, brdf_val
!  ENDDO

  nchannels=15
  do ich = 1, nchannels
    input_chan(ich)=ich
  end do
  input_ems(:)=0.0
  input_brdf(:)=0.0

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.          ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.         ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp  = .TRUE.  ! Allow interpolation of input profile
  opts % rt_all % addrefrac         = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds          = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl          = .FALSE. ! Don't include aerosol effects
  opts % rt_ir % ozone_data         = .FALSE.  ! We have an ozone profile
  !opts % rt_ir % ozone_data         = .TRUE.  ! We have an ozone profile
  opts % rt_ir % co2_data           = .FALSE. ! We do not have profiles
  opts % rt_ir % n2o_data           = .FALSE. !   for any other constituents
  opts % rt_ir % ch4_data           = .FALSE. !
  opts % rt_ir % co_data            = .FALSE. !
  opts % rt_mw % clw_data           = .FALSE. !

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  write(*,*) 'coef_filename', trim(coef_filename)
  CALL rttov_read_coefs(errorstatus, coefs, opts, form_coef='formatted', file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  !CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  !IF (errorstatus /= errorstatus_success) THEN
  !  WRITE(*,*) 'error in rttov options'
  !  CALL rttov_exit(errorstatus)
  !ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  ! In general the number of channels simulated per profile can vary, but for
  ! this example we simulate all specified channels for each profile.

  ALLOCATE(nchan(nprof))     ! Number of channels per profile
  nchan(:) = nchannels
  nchanprof = SUM(nchan(:))  ! Size of chanprof array is total number of channels over all profiles

  ! Pack channels and input emissivity arrays
  ALLOCATE(chanprof(nchanprof))
  ALLOCATE(emis(nchanprof))
  ALLOCATE(brdf(nchanprof))

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchan(j)
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = input_chan(jch)
      emis(nch)          = input_ems(jch)
      brdf(nch)          = input_brdf(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 4. Allocate profiles, radiance and transmittance structures
  ! --------------------------------------------------------------------------

  asw = 1 ! Switch for allocation passed into RTTOV subroutines

  ! Allocate input profile arrays
  ALLOCATE(profiles(nprof), stat=alloc_status(1))
  IF (ANY(alloc_status /= 0)) THEN
    WRITE(*,*) 'mem allocation error for profile array'
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_alloc_prof(         &
      & errorstatus,             &
      & nprof,                   &
      & profiles,                &
      & nlevels,                 &
      & opts,                    &
      & asw,                     &
      & coefs=coefs,             &
      & init=.TRUE._jplm         )
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for profile arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate output radiance arrays
  CALL rttov_alloc_rad( &
      & errorstatus,    &
      & nchanprof,      &
      & radiance,       &
      & nlevels-1_jpim, &
      & asw)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for radiance arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate transmittance structure
  CALL rttov_alloc_transmission( &
      & errorstatus,             &
      & transmission,            &
      & nlevels-1_jpim,          &
      & nchanprof,               &
      & asw,                     &
      & init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for transmission arrays'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Allocate arrays for surface emissivity
  ALLOCATE(calcemis(nchanprof), stat=alloc_status(1))
  ALLOCATE(emissivity(nchanprof), stat=alloc_status(2))
  IF (ANY(alloc_status /= 0)) THEN
    WRITE(*,*) 'mem allocation error for emissivity arrays'
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  ! Allocate arrays for surface reflectance
  ALLOCATE(calcrefl(nchanprof), stat=alloc_status(1))
  ALLOCATE(reflectance(nchanprof), stat=alloc_status(2))
  IF (ANY(alloc_status /= 0)) THEN
    WRITE(*,*) 'memallocation error for reflectance arrays'
    CALL rttov_exit(errorstatus_fatal)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== READ profiles == start =============
  WRITE(*,*) 'START SUBSTITUTE PROFILE'
  DO iprof = 1, nprof
    profiles(iprof)%p(:)=tmp_p(:,iprof)
    profiles(iprof)%t(:)=tmp_t(:,iprof)
    profiles(iprof)%q(:)=tmp_q(:,iprof)
    profiles(iprof)%s2m%t=tmp_t2m(1,iprof)
    profiles(iprof)%s2m%q=tmp_q2m(1,iprof)
    profiles(iprof)%s2m%p=tmp_p2m(1,iprof)
    profiles(iprof)%s2m%u=tmp_u2m(1,iprof)
    profiles(iprof)%s2m%v=tmp_v2m(1,iprof)
    profiles(iprof)%s2m%wfetc= 100000.

    profiles(iprof) % skin % t = tmp_t2m(1,iprof)
    profiles(iprof) % skin % fastem(1) = 3.0
    profiles(iprof) % skin % fastem(2) = 5.0
    profiles(iprof) % skin % fastem(3) =15.0
    profiles(iprof) % skin % fastem(4) = 0.1
    profiles(iprof) % skin % fastem(5) = 0.3

    profiles(iprof) % skin % surftype = int(tmp_land(1,iprof))
    profiles(iprof) % skin % watertype = 1

    profiles(iprof) % elevation = tmp_elev(1, iprof)
    profiles(iprof) % latitude  = tmp_lat(1, iprof)
    profiles(iprof) % longitude = tmp_lon(1, iprof)

    profiles(iprof)%zenangle=tmp_saze(1,iprof)
    profiles(iprof)%azangle=tmp_saaz(1,iprof)
    profiles(iprof)%sunzenangle=tmp_soze(1,iprof)
    profiles(iprof)%sunazangle=tmp_soaz(1,iprof)

    profiles(iprof) % ctp       = 500.
    profiles(iprof) % cfraction = 0.0
  ENDDO

  !WRITE(*,*) 'END SUBSTITUTE PROFILE'
  !WRITE(*,*) 'profiles(iprof)%p'
  !WRITE(*,*) profiles(iprof)%p(:)
  !WRITE(*,*) 'profiles(iprof)%s2m%p'
  !WRITE(*,*) profiles(1)%s2m%p
  !WRITE(*,*) minval(profiles(:)%s2m%p), maxval(profiles(:)%s2m%p)
  !WRITE(*,*) minval(tmp_p2m(:,:)), maxval(tmp_p2m(:,:))


!  OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
!  IF (ios /= 0) THEN
!    WRITE(*,*) 'error opening profile file ios= ', ios
!    CALL rttov_exit(errorstatus_fatal)
!  ENDIF
!  CALL rttov_skipcommentline(iup, errorstatus)
!
!  ! Loop over all profiles and read data for each one
!  DO iprof = 1, nprof
!
!    ! READ pressure (hPa), temp (K), WV (ppmv), O3 (ppmv)
!    READ(iup,*) profiles(iprof) % p(:)
!    CALL rttov_skipcommentline(iup, errorstatus)
!    READ(iup,*) profiles(iprof) % t(:)
!    CALL rttov_skipcommentline(iup, errorstatus)
!    READ(iup,*) profiles(iprof) % q(:)
!    CALL rttov_skipcommentline(iup, errorstatus)
    !READ(iup,*) profiles(iprof) % o3(:)
    !CALL rttov_skipcommentline(iup, errorstatus)

    ! 2 meter air variables
!    READ(iup,*) profiles(iprof) % s2m % t, &
!              & profiles(iprof) % s2m % q, &
!              & profiles(iprof) % s2m % p, &
!              & profiles(iprof) % s2m % u, &
!              & profiles(iprof) % s2m % v, &
!              & profiles(iprof) % s2m % wfetc
!    CALL rttov_skipcommentline(iup, errorstatus)
!
!    ! Skin variables
!    READ(iup,*) profiles(iprof) % skin % t, &
!              & profiles(iprof) % skin % fastem   ! FASTEM only applies to MW instruments
!    CALL rttov_skipcommentline(iup, errorstatus)
!
!    ! Surface type and water type
!    READ(iup,*) profiles(iprof) % skin % surftype, &
!              & profiles(iprof) % skin % watertype
!    CALL rttov_skipcommentline(iup, errorstatus)
!
!    ! Elevation, latitude and longitude
!    READ(iup,*) profiles(iprof) % elevation, &
!              & profiles(iprof) % latitude,  &
!              & profiles(iprof) % longitude
!    CALL rttov_skipcommentline(iup, errorstatus)
!
!    ! Satellite and solar angles
!    READ(iup,*) profiles(iprof) % zenangle,    &
!              & profiles(iprof) % azangle,     &
!              & profiles(iprof) % sunzenangle, &
!              & profiles(iprof) % sunazangle
!    CALL rttov_skipcommentline(iup, errorstatus)
!
!    ! Cloud variables for simple cloud scheme (set cfraction to 0. to turn this off)
!    READ(iup,*) profiles(iprof) % ctp, &
!              & profiles(iprof) % cfraction
!    CALL rttov_skipcommentline(iup, errorstatus)

!  ENDDO
!  CLOSE(iup)

  !========== READ profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! Copy emissivities into RTTOV input emissivity array
  emissivity(:) % emis_in = emis(:)

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less
  calcemis(:) = emis(:) <= 0._jprb

  ! Copy BRDFs into RTTOV input reflectance array
  reflectance(:) % refl_in = brdf(:)

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  calcrefl(:) = brdf(:) <= 0._jprb

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  !CALL rttov_parallel_direct(                &
  CALL rttov_direct(                &
        & errorstatus,              &! out   error flag
        & chanprof,                 &! in    channel and profile index structure
        & opts,                     &! in    options structure
        & profiles,                 &! in    profile array
        & coefs,                    &! in    coefficients strucutre
        & transmission,             &! inout computed transmittances
        & radiance,                 &! inout computed radiances
        & calcemis    = calcemis,   &! in    flag for internal emissivity calcs
        & emissivity  = emissivity, &! inout input/output emissivities per channel
        & calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
        & reflectance = reflectance) ! inout input/output BRDFs per channel

        !& nthreads    = 40,         &
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --- Output the results --------------------------------------------------

  ! Open output file where results are written
  OPEN(ioout, file='output_'//NameOfRoutine//'.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' '
  CALL rttov_print_opts(opts, lu=ioout)

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    !
    !     OUTPUT RESULTS
    !
    nprint = 1 + INT((nchannels-1)/10)
    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)

    bt_out(1:15,iprof)=radiance%bt(1+joff:nchannels+joff)
    DO ilev = 1, nlevels
       tran_out(ilev,1:15,iprof)=transmission % tau_levels(ilev,1+joff:nchannels+joff) 
    ENDDO

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SATELLITE REFLECTANCES (BRF):'
      WRITE(ioout,444) (radiance % refl(j), j = 1+joff, nchannels+joff)
    ENDIF
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
    WRITE(ioout,222) (radiance % total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED OVERCAST RADIANCES:'
    WRITE(ioout,222) (radiance % cloudy(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE TO SPACE TRANSMITTANCE:'
    WRITE(ioout,4444) (transmission % tau_total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SURFACE BRDF:'
      WRITE(ioout,444) (reflectance(j) % refl_out, j = 1+joff, nchannels+joff)
    ENDIF

    IF (nchan(nprof) <= 20) THEN
      DO np = 1, nprint
          WRITE(ioout,*)' '
          WRITE(ioout,*)'Level to space transmittances for channels'
          WRITE(ioout,1115) (chanprof(j+joff) % chan, &
                  & j = 1+(np-1)*10, MIN(np*10, nchannels))
          DO ilev = 1, nlevels
            DO j = 1 + (np-1)*10, MIN(np*10, nchannels)
              ! Select transmittance based on channel type (VIS/NIR or IR)
              IF (coefs % coef % ss_val_chn(chanprof(j+joff) % chan) == 2) THEN
                trans_out(j - (np-1)*10) = transmission % tausun_levels_path1(ilev,j+joff)
              ELSE
                trans_out(j - (np-1)*10) = transmission % tau_levels(ilev,j+joff)
              ENDIF
            ENDDO
            WRITE(ioout,4445) ilev, trans_out(1:j-1-(np-1)*10)
          ENDDO
          WRITE(ioout,1115) (chanprof(j+joff) % chan, &
                  & j = 1+(np-1)*10, MIN(np*10, nchannels))
      ENDDO
    ENDIF
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  ! --- End of output section -----------------------------------------------


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (emis,        stat=alloc_status(1))
  DEALLOCATE (brdf,        stat=alloc_status(2))
  DEALLOCATE (nchan,       stat=alloc_status(3))
  DEALLOCATE (chanprof,    stat=alloc_status(4))
  DEALLOCATE (emissivity,  stat=alloc_status(5))
  DEALLOCATE (calcemis,    stat=alloc_status(6))
  DEALLOCATE (reflectance, stat=alloc_status(7))
  DEALLOCATE (calcrefl,    stat=alloc_status(8))

  IF (ANY(alloc_status /= 0)) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  asw = 0 ! Switch for deallocation passed into RTTOV subroutines

  ! Deallocate radiance arrays
  CALL rttov_alloc_rad(errorstatus, nchannels, radiance, nlevels-1_jpim, asw)
  IF(errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'radiance deallocation error'
  ENDIF

  ! Deallocate transmission arrays
  CALL rttov_alloc_transmission(errorstatus, transmission, nlevels-1_jpim, nchannels, asw)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'transmission deallocation error'
  ENDIF

  ! Deallocate profile arrays
  CALL rttov_alloc_prof(errorstatus, nprof, profiles, nlevels, opts, asw)
  IF (errorstatus /= errorstatus_success .or. alloc_status(1) /= 0) THEN
    WRITE(*,*) 'profile deallocation error'
  ENDIF
  DEALLOCATE(profiles, stat=alloc_status(1))
  IF (alloc_status(1) /= 0) THEN
    WRITE(*,*) 'mem deallocation error for profile array'
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
1115 FORMAT(3X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
4444 FORMAT(1X,10F8.4)
4445 FORMAT(1X,I2,10F8.4)
777  FORMAT(/,A,A8,I3)

END SUBROUTINE AMSUA_fwd
