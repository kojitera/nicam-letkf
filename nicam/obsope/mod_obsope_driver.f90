MODULE mod_obsope_driver
  USE mpi
  USE mod_adm
  IMPLICIT NONE
  PUBLIC

  CHARACTER(ADM_MAXFNAME), SAVE :: veg_base  = ''
  CHARACTER(ADM_MAXFNAME), SAVE :: co2_fname = ''
  INTEGER                       :: co2_step
  REAL(8), ALLOCATABLE, SAVE    :: veg(:,:,:)
  REAL(8), ALLOCATABLE, SAVE    :: veg_pl(:,:,:)

  LOGICAL :: assimilate_prepbufr = .false.
  LOGICAL :: assimilate_amsua    = .false.
  LOGICAL :: assimilate_atms     = .false.
  LOGICAL :: assimilate_gsmap    = .false.
  LOGICAL :: assimilate_iasi     = .false.
  LOGICAL :: assimilate_cris     = .false.
  LOGICAL :: assimilate_airs     = .false.

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE obsope_init
  USE mod_gtl, only: GTL_input_var2_da
  USE mod_obsope_common
  USE mod_obsope_prepbufr
  USE mod_obsope_amsua
  USE mod_obsope_atms
  USE mod_obsope_gsmap
  USE mod_obsope_iasi
  USE mod_obsope_cris
  USE mod_obsope_airs
  IMPLICIT NONE
  INTEGER :: l, n, nn
  REAL(8), ALLOCATABLE :: lon(:)
  REAL(8), ALLOCATABLE :: lat(:)
  REAL(8), ALLOCATABLE :: lev(:)

  INTEGER :: it, ix, iy, iz

  NAMELIST / obsope_param /   &
    istep, islot, smem, emem, &
    nbv,                      &
    veg_base,                 &
    assimilate_prepbufr,      &
    assimilate_amsua,         &
    assimilate_atms,          &
    assimilate_gsmap,         &
    assimilate_iasi,          &
    assimilate_cris,          &
    assimilate_airs,          &
    co2_fname, co2_step,      &
    start_mem_zero,           &  ! [add] Koji 2017.09.08
    flush_text

  OPEN(101,file='obsope.cnf')
  READ(101,NML=obsope_param) 
  CLOSE(101)

  ALLOCATE(veg(ADM_gall,ADM_KNONE,ADM_lall))
  ALLOCATE(veg_pl(ADM_GALL_PL,ADM_KNONE,ADM_LALL_PL))
  ALLOCATE(icoland(ADM_gall,ADM_lall))

  ! INDEX (Land(=1) or Ocean(=0)) => (Land(=0) or Ocean(=1) for RTTOV) 
  !CALL GTL_input_var2(trim(veg_base), veg, veg_pl, ADM_KNONE, ADM_KNONE, 8)
  CALL GTL_input_var2_da(trim(veg_base), veg, ADM_KNONE, ADM_KNONE, 1, 8)
  icoland(:,:) = 0.0d0
  DO l=1, ADM_lall
  DO n=1, ADM_gall
    IF(veg(n,ADM_KNONE,l) == ocean_value) icoland(n,l) = 1.0d0
  END DO
  END DO

  IF( TRIM(co2_fname) /= '' ) THEN
    WRITE(ADM_LOG_FID,*) 'READING CO2'
    OPEN(102,FILE=TRIM(co2_fname),FORM='unformatted',access='sequential')
    READ(102) ak
    READ(102) bk
    READ(102) psurf_co2
    DO it = 1, co2_step
    DO iz = 1, nlev_co2
      !READ(102) co2_org(:,:,nlev_co2-iz+1)
      READ(102) co2_org(:,:,iz)
    END DO
    END DO
    CLOSE(102)

    WRITE(ADM_LOG_FID,*) minval(co2_org), maxval(co2_org)

    DO iz = 1, nlev_co2
    DO iy = 1, nlat_co2
    DO ix = 1, nlon_co2
      !pres_co2(ix,iy,nlev_co2-iz+1)=ak(iz)+bk(iz)*psurf_co2(ix,iy)
      pres_co2(ix,iy,iz)=ak(iz)+bk(iz)*psurf_co2(ix,iy)
    END DO
    END DO
    END DO

    DO ix = 1, nlon_co2
      lon_co2(ix)=2.5*REAL(ix-1)
    END DO

    WHERE( lon_co2(:) > 180.0 ) lon_co2(:)=lon_co2(:)-360.0

    DO iy = 1, nlat_co2
      lat_co2(iy)=-90.0 + 2.5*REAL(iy-1)
    END DO

  END IF

  IF(assimilate_prepbufr) THEN
    CALL obsope_prepbufr_init
    CALL obsope_prepbufr_read
    ALLOCATE( lon(nobs_prepbufr) )
    ALLOCATE( lat(nobs_prepbufr) )
    ALLOCATE( lev(nobs_prepbufr) )
    lon(:)=dble(prep_lon)
    lat(:)=dble(prep_lat)
    lev(:)=dble(prep_lev)
    CALL set_coef_interpolate(nobs_prepbufr, lon, lat, prep_inprc,    &
           prep_l_index, prep_n1_index, prep_n2_index, prep_n3_index, &
           prep_w1, prep_w2, prep_w3)
    DEALLOCATE( lon )
    DEALLOCATE( lat )
    DEALLOCATE( lev )
  END IF

  IF(assimilate_amsua) THEN
    CALL obsope_amsua_init
    CALL obsope_amsua_read
    DO nn = 1, num_satellite_amsua
      IF( nobs_amsua(nn) > 0 ) THEN
        ALLOCATE( lon(nobs_amsua(nn)) )
        ALLOCATE( lat(nobs_amsua(nn)) )
        lon(:) = dble( amsua(nn)%lon )
        lat(:) = dble( amsua(nn)%lat )
        CALL set_coef_interpolate(nobs_amsua(nn), lon, lat, amsua(nn)%inprc,  &
             amsua(nn)%l_index, amsua(nn)%n1_index, amsua(nn)%n2_index,       &
             amsua(nn)%n3_index, amsua(nn)%w1, amsua(nn)%w2, amsua(nn)%w3     )
        DEALLOCATE( lon )
        DEALLOCATE( lat )
      END IF
    END DO
  END IF

  IF(assimilate_atms) THEN
    CALL obsope_atms_init
    CALL obsope_atms_read
    DO nn = 1, num_satellite_atms
      IF( nobs_atms(nn) > 0 ) THEN
        ALLOCATE( lon(nobs_atms(nn)) )
        ALLOCATE( lat(nobs_atms(nn)) )
        lon(:) = dble( atms(nn)%lon )
        lat(:) = dble( atms(nn)%lat )
        CALL set_coef_interpolate(nobs_atms(nn), lon, lat, atms(nn)%inprc,  &
             atms(nn)%l_index, atms(nn)%n1_index, atms(nn)%n2_index,       &
             atms(nn)%n3_index, atms(nn)%w1, atms(nn)%w2, atms(nn)%w3     )
        DEALLOCATE( lon )
        DEALLOCATE( lat )
      END IF
    END DO
  END IF

  IF(assimilate_gsmap) THEN
    CALL obsope_gsmap_init
  END IF

  IF(assimilate_iasi) THEN
    CALL obsope_iasi_init
    CALL obsope_iasi_read
    DO nn = 1, num_satellite_iasi
      IF( nobs_iasi(nn) > 0 ) THEN
        ALLOCATE( lon(nobs_iasi(nn)) )
        ALLOCATE( lat(nobs_iasi(nn)) )
        lon(:) = dble( iasi(nn)%lon )
        lat(:) = dble( iasi(nn)%lat )
        CALL set_coef_interpolate(nobs_iasi(nn), lon, lat, iasi(nn)%inprc,  &
             iasi(nn)%l_index,  iasi(nn)%n1_index, iasi(nn)%n2_index,       &
             iasi(nn)%n3_index, iasi(nn)%w1, iasi(nn)%w2, iasi(nn)%w3     )
        DEALLOCATE( lon )
        DEALLOCATE( lat )
      END IF
    END DO
  END IF
  IF(assimilate_cris) THEN
    CALL obsope_cris_init
  END IF
  IF(assimilate_airs) THEN
    CALL obsope_airs_init
  END IF

  nobs_total=nobs_prepbufr


END SUBROUTINE obsope_init
!------------------------------------------------------------------------------
SUBROUTINE obsope_main(imem)
  USE mod_obsope_common
  USE mod_obsope_prepbufr
  USE mod_obsope_amsua
  USE mod_obsope_atms
  USE mod_obsope_gsmap
  USE mod_obsope_iasi
  USE mod_obsope_cris
  USE mod_obsope_airs
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  INTEGER :: nn, ic
  INTEGER :: tmp
  INTEGER :: ierr
  REAL(8), ALLOCATABLE :: lev(:)

  ! PREPBUFR
  IF(assimilate_prepbufr) THEN
    ALLOCATE( lev(nobs_prepbufr) )
    lev(:)=dble(prep_lev)
    CALL getklev(nobs_prepbufr, prep_elem, lev, prep_inprc,   &
         prep_l_index, prep_n1_index, prep_n2_index, prep_n3_index, &
         prep_w1, prep_w2, prep_w3, prep_klev, prep_kfact)
    DEALLOCATE(lev)

    CALL interpolate_prepbufr

    CALL output_prepbufr(imem)

  END IF

  ! AMSU-A
  IF(assimilate_amsua) THEN
    DO nn = 1, num_satellite_amsua
      IF( nobs_amsua(nn) > 0 ) THEN
        CALL interpolate_amsua(nn)
      END IF
    END DO
    CALL calc_radiance_amsua

  
    IF( ADM_prc_me <= num_satellite_amsua ) then
      nn = ADM_prc_me
      tmp=15
      amsua(nn)%weight(:,:,:)=0.0d0
      amsua(nn)%weight_maxlev(:,:)=0
      amsua(nn)%vbc_pred(:,:,:)=0.0d0
      IF( nobs_amsua(nn) /= 0 ) THEN
        CALL calc_weighting_function(ADM_vlayer, nobs_amsua(nn), tmp,   &
                                     DBLE(amsua(nn)%obsdata_3d(:,:,3)), &
                                     DBLE(amsua(nn)%trans(:,:,:)),      &
                                     DBLE(amsua(nn)%weight(:,:,:)),     &
                                     amsua(nn)%weight_maxlev(:,:)    )
        CALL calc_vbc_amsua(nn)

        CALL quality_control_amsua(nn)

        CALL output_amsua(imem, nn)


      END IF
      IF( imem == -1 ) THEN
        CALL update_vbc_amsua( imem, nn )

        CALL update_scanbias_amsua( nn )
      END IF
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  END IF

  ! ATMS
  IF(assimilate_atms) THEN
    DO nn = 1, num_satellite_atms
      IF( nobs_atms(nn) > 0 ) THEN
        CALL interpolate_atms(nn)
      END IF
    END DO
    CALL calc_radiance_atms

    IF( ADM_prc_me <= num_satellite_atms ) then
      nn = ADM_prc_me
      tmp=22
      atms(nn)%weight(:,:,:)=0.0d0
      atms(nn)%weight_maxlev(:,:)=0
      atms(nn)%vbc_pred(:,:,:)=0.0d0
      IF( nobs_atms(nn) /= 0 ) THEN
        CALL calc_weighting_function(ADM_vlayer, nobs_atms(nn), tmp,   &
                                     DBLE(atms(nn)%obsdata_3d(:,:,3)), &
                                     DBLE(atms(nn)%trans(:,:,:)),      &
                                     DBLE(atms(nn)%weight(:,:,:)),     &
                                     atms(nn)%weight_maxlev(:,:)    )
        WRITE(ADM_LOG_FID,*) 'End calc_weighting_function'
        IF(flush_text) FLUSH(ADM_LOG_FID)

        CALL calc_vbc_atms(nn)
        WRITE(ADM_LOG_FID,*) 'End calc_vbc_atms'
        IF(flush_text) FLUSH(ADM_LOG_FID)

        CALL quality_control_atms(nn)
        WRITE(ADM_LOG_FID,*) 'End quality_control_atms'
        IF(flush_text) FLUSH(ADM_LOG_FID)

        CALL output_atms(imem, nn)
        WRITE(ADM_LOG_FID,*) 'End output_atms'
        IF(flush_text) FLUSH(ADM_LOG_FID)


      END IF
      IF( imem == -1 ) THEN
        CALL update_vbc_atms( imem, nn )

        CALL update_scanbias_atms( nn )
      END IF
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  END IF

  ! GSMaP
  IF(assimilate_gsmap) THEN

  END IF

  ! IASI
  IF(assimilate_iasi) THEN
    DO nn = 1, num_satellite_iasi
      IF( nobs_iasi(nn) > 0 ) THEN
        CALL interpolate_iasi(nn)
      END IF
    END DO
    CALL calc_radiance_iasi

    IF( ADM_prc_me <= num_satellite_iasi ) then
      nn  = ADM_prc_me
      tmp = 4
      iasi(nn)%weight(:,:,:)=0.0d0
      iasi(nn)%weight_maxlev(:,:)=0
      iasi(nn)%vbc_pred(:,:,:)=0.0d0
      IF( nobs_iasi(nn) /= 0 ) THEN
        CALL calc_weighting_function(ADM_vlayer, nobs_iasi(nn), tmp,   &
                                     DBLE(iasi(nn)%obsdata_3d(:,:,3)), &
                                     DBLE(iasi(nn)%trans(:,:,:)),      &
                                     DBLE(iasi(nn)%weight(:,:,:)),     &
                                     iasi(nn)%weight_maxlev(:,:)    )
        CALL calc_vbc_iasi(nn)

        CALL quality_control_iasi(nn)

        CALL output_iasi(imem, nn)

      END IF
      IF( imem == -1 ) THEN
        CALL update_vbc_iasi( imem, nn )

        CALL update_scanbias_iasi( nn )
      END IF
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  END IF

  ! CrIS
  IF(assimilate_cris) THEN

  END IF

  ! AIRS
  IF(assimilate_airs) THEN

  END IF

END SUBROUTINE obsope_main
!------------------------------------------------------------------------------
END MODULE mod_obsope_driver
