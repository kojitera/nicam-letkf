MODULE mod_obsope_iasi
  USE common_tvs_nicam, ONLY: &
    rttv_plat_metop,          &
    rttv_plat_metop2,         &
    rttv_inst_iasi,           &
    rttv_chan_iasi 

  IMPLICIT NONE
  PUBLIC

  INTEGER, PRIVATE, PARAMETER   :: ninstrument=2
  INTEGER, PRIVATE, PARAMETER   :: maxtvsch=4
  INTEGER, PRIVATE, PARAMETER   :: maxvbc=8
  INTEGER, PRIVATE, PARAMETER   :: maxfoot=30
  INTEGER, PRIVATE, SAVE        :: ntvs
  INTEGER, PRIVATE, SAVE        :: tvsinst(3,ninstrument)
  INTEGER, PRIVATE, SAVE        :: tvsch(maxtvsch,ninstrument)
  INTEGER, PRIVATE, SAVE        :: ntvsch(ninstrument)
  INTEGER, PRIVATE, SAVE        :: nfootp(ninstrument)
  CHARACTER(4), PRIVATE, SAVE   :: tvsname(ninstrument)
  CHARACTER(256), PRIVATE, SAVE :: rttovcoef_fname(ninstrument) = ''
  CHARACTER(256), PRIVATE       :: input_fname_iasi = ''
  CHARACTER(256), PRIVATE       :: output_dirname_iasi = ''
  CHARACTER(2), PRIVATE         :: output_basename_iasi = ''
  CHARACTER(256), PRIVATE       :: ominusb_fname_iasi   = ''
  CHARACTER(256), PRIVATE       :: scanbias_fname = ''
  CHARACTER(256), PRIVATE       :: vbc_coef_fname = ''
  CHARACTER(256), PRIVATE       :: vbc_coef_out_fname = ''
  CHARACTER(256), PRIVATE       :: scanbias_est_ifname = ''
  CHARACTER(256), PRIVATE       :: scanbias_est_ofname = ''
  LOGICAL, PRIVATE, SAVE        :: output_text = .false.

  INTEGER, PRIVATE, PARAMETER :: nv3d=3
  INTEGER, PRIVATE, PARAMETER :: nv2d=7

  INTEGER :: num_satellite_iasi

  TYPE iasi_base
    INTEGER, ALLOCATABLE :: datetime(:,:)
    INTEGER, ALLOCATABLE :: said(:)
    REAL(4), ALLOCATABLE :: lat(:)
    REAL(4), ALLOCATABLE :: lon(:)
    REAL(4), ALLOCATABLE :: elev(:)
    REAL(4), ALLOCATABLE :: odat(:,:)
    REAL(4), ALLOCATABLE :: odat_bc(:,:)
    REAL(4), ALLOCATABLE :: err(:,:)
    REAL(4), ALLOCATABLE :: saza(:)
    REAL(4), ALLOCATABLE :: soza(:)
    REAL(4), ALLOCATABLE :: soaz(:)
    REAL(4), ALLOCATABLE :: saaz(:)
    INTEGER, ALLOCATABLE :: fov(:)
    INTEGER, ALLOCATABLE :: lsql(:)
    INTEGER, ALLOCATABLE :: qc_tmp(:,:)
    INTEGER, ALLOCATABLE :: qc(:,:)
    REAL(4), ALLOCATABLE :: co2(:,:)
    REAL(4), ALLOCATABLE :: obsdata_3d_tmp(:,:,:)
    REAL(4), ALLOCATABLE :: obsdata_3d(:,:,:)
    REAL(4), ALLOCATABLE :: obsdata_2d_tmp(:,:)
    REAL(4), ALLOCATABLE :: obsdata_2d(:,:)
    REAL(8), ALLOCATABLE :: bt_tmp(:,:)
    REAL(8), ALLOCATABLE :: bt(:,:)
    REAL(8), ALLOCATABLE :: ominusb(:,:)
    REAL(8), ALLOCATABLE :: trans_tmp(:,:,:)
    REAL(8), ALLOCATABLE :: trans(:,:,:)
    REAL(4), ALLOCATABLE :: wk(:,:)
    REAL(4), ALLOCATABLE :: lwp(:)
    LOGICAL, ALLOCATABLE :: inprc(:)
    INTEGER, ALLOCATABLE :: l_index(:)
    INTEGER, ALLOCATABLE :: n1_index(:)
    INTEGER, ALLOCATABLE :: n2_index(:)
    INTEGER, ALLOCATABLE :: n3_index(:)
    REAL(8), ALLOCATABLE :: w1(:)
    REAL(8), ALLOCATABLE :: w2(:)
    REAL(8), ALLOCATABLE :: w3(:)
    REAL(8), ALLOCATABLE :: vbcf_scan(:,:)
    REAL(8), ALLOCATABLE :: weight(:,:,:)
    INTEGER, ALLOCATABLE :: weight_maxlev(:,:)
    REAL(8), ALLOCATABLE :: vbc_pred(:,:,:)
    REAL(8), ALLOCATABLE :: airmass_bias(:,:,:)
    REAL(8), ALLOCATABLE :: lsql_model_tmp(:)
    REAL(8), ALLOCATABLE :: lsql_model(:)
  END TYPE iasi_base
    
  TYPE(iasi_base), allocatable :: iasi(:)

  INTEGER, ALLOCATABLE :: nobs_iasi(:)

  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbcf_scan(:,:,:)
  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbcf(:,:,:)
  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbca(:,:,:)
  REAL(8), PRIVATE :: time_rttov_iasi(2)

  PRIVATE :: set_instrument

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE obsope_iasi_init
  USE mod_adm
  IMPLICIT NONE
  INTEGER :: ierr

  NAMELIST / iasi_cnf /   &
    input_fname_iasi,     &
    output_dirname_iasi,  &
    output_basename_iasi, &
    rttovcoef_fname,       &
    ominusb_fname_iasi,    &
    scanbias_fname,        &
    vbc_coef_fname,        &
    vbc_coef_out_fname,    &
    scanbias_est_ifname,   &
    scanbias_est_ofname,   &
    output_text

  CALL set_instrument

  OPEN(121, FILE='obsope.cnf')
  READ(121, NML=iasi_cnf)
  CLOSE(121)

END SUBROUTINE obsope_iasi_init
!------------------------------------------------------------------------------
SUBROUTINE obsope_iasi_read
  USE mod_adm
  USE mod_obsope_common
  USE mod_scanbias
  USE mod_vbc
  IMPLICIT NONE
  INTEGER :: ierr
  INTEGER :: n, nn
  INTEGER :: ic
  REAL(4) :: cos_tmp

  OPEN(122,file=trim(input_fname_iasi),  form='unformatted', &
       access='sequential', status='old', iostat=ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in opening the file', trim(input_fname_iasi)
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    FLUSH(ADM_LOG_FID)
    CALL ADM_proc_stop
  END IF
  READ(122) num_satellite_iasi
  ALLOCATE( iasi(num_satellite_iasi) )
  ALLOCATE( nobs_iasi(num_satellite_iasi) )

  READ(122) ntvs
  READ(122) nobs_iasi(:)
  WRITE(ADM_LOG_FID,*) 'nobs_iasi ', nobs_iasi
  FLUSH(ADM_LOG_FID)

  DO nn = 1, num_satellite_iasi
    ALLOCATE( iasi(nn)%datetime      ( 6, nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%said          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%lon           ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%lat           ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%elev          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%odat          ( maxtvsch,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%odat_bc       ( maxtvsch,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%err           ( maxtvsch,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%saza          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%soza          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%soaz          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%saaz          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%fov           ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%lsql          ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%qc_tmp        ( maxtvsch,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%qc            ( maxtvsch,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%co2           ( ADM_vlayer,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%obsdata_3d_tmp( ADM_vlayer,nobs_iasi(nn),nv3d) )
    ALLOCATE( iasi(nn)%obsdata_3d    ( ADM_vlayer,nobs_iasi(nn),nv3d) )
    ALLOCATE( iasi(nn)%obsdata_2d_tmp( nobs_iasi(nn),nv2d) )
    ALLOCATE( iasi(nn)%obsdata_2d    ( nobs_iasi(nn),nv2d) )
    ALLOCATE( iasi(nn)%lwp           ( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%bt_tmp        ( ntvsch(nn),nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%bt            ( ntvsch(nn),nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%ominusb       ( maxtvsch,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%trans_tmp     ( ADM_vlayer, ntvsch(nn),nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%trans         ( ADM_vlayer, ntvsch(nn),nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%wk            ( 28,nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%inprc         ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%l_index       ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%n1_index      ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%n2_index      ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%n3_index      ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%w1            ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%w2            ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%w3            ( nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%vbcf_scan     ( maxfoot, maxtvsch ))
    ALLOCATE( iasi(nn)%weight        ( ADM_vlayer, maxtvsch, nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%weight_maxlev ( maxtvsch, nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%vbc_pred      ( maxvbc, maxtvsch, nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%airmass_bias  ( maxvbc, maxtvsch, nobs_iasi(nn) ))
    ALLOCATE( iasi(nn)%lsql_model_tmp( nobs_iasi(nn)) )
    ALLOCATE( iasi(nn)%lsql_model    ( nobs_iasi(nn)) )
  END DO

  DO nn = 1, num_satellite_iasi
  DO n = 1, nobs_iasi(nn)
    READ(122) iasi(nn)%wk(1:18+ntvsch(nn),n)
  END DO
  END DO

  DO nn = 1, num_satellite_iasi
  !WRITE(ADM_LOG_FID,*) 'Satellite ', nn
  DO n = 1, nobs_iasi(nn)
    iasi(nn)%datetime(:,n) = INT(iasi(nn)%wk(1:6,n))
    iasi(nn)%lat(n)        = iasi(nn)%wk( 7,n)
    iasi(nn)%lon(n)        = iasi(nn)%wk( 8,n)
    iasi(nn)%said(n)       = iasi(nn)%wk( 9,n)
    iasi(nn)%fov(n)        = iasi(nn)%wk(11,n)
    iasi(nn)%lsql(n)       = iasi(nn)%wk(12,n)
    iasi(nn)%saza(n)       = iasi(nn)%wk(13,n)
    iasi(nn)%soza(n)       = iasi(nn)%wk(14,n)
    iasi(nn)%elev(n)       = iasi(nn)%wk(15,n)/1000.0
    iasi(nn)%soaz(n)       = iasi(nn)%wk(17,n)
    iasi(nn)%saaz(n)       = iasi(nn)%wk(18,n)
    iasi(nn)%odat(1:ntvsch(nn),n) = iasi(nn)%wk(19:18+ntvsch(nn),n)
    !WRITE(ADM_LOG_FID,'(6f12.5)') iasi(nn)%lon(n), iasi(nn)%lat(n),   &
    !                              iasi(nn)%saza(n), iasi(nn)%saaz(n), &
    !                              iasi(nn)%soza(n), iasi(nn)%soaz(n) 
  END DO
  !WRITE(ADM_LOG_FID,*) minval(iasi(nn)%lat(:)), maxval(iasi(nn)%lat(:))
  !WRITE(ADM_LOG_FID,*) minval(iasi(nn)%lon(:)), maxval(iasi(nn)%lon(:))
  !WRITE(ADM_LOG_FID,*) minval(iasi(nn)%saza(:)), maxval(iasi(nn)%saza(:))
  END DO


  DO nn = 1, num_satellite_iasi
    !iasi(nn)%co2(:,:)=500.0
    iasi(nn)%lwp(:)=0.0
    iasi(nn)%err(:,:)=0.5
    iasi(nn)%qc(:,:)=1
  END DO

  ! Read scan bias
  ALLOCATE( vbcf_scan(maxfoot,maxtvsch,num_satellite_iasi) )
  IF( scanbias_fname /= '' ) THEN
    CAll vbc_scan_read(trim(scanbias_fname), vbcf_scan, maxfoot, maxtvsch, &
                       num_satellite_iasi, tvsinst, tvsch, ntvsch)
  ELSE
    WRITE(ADM_LOG_FID,*) 'Caution! scan bias file is not appointed'
    vbcf_scan(:,:,:)=0.0d0
  END IF

  DO nn = 1, num_satellite_iasi
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(2i3,120F8.3)') nn, tvsch(ic,nn), vbcf_scan(:,ic,nn)
  END DO
  END DO

  ! Read varBC coefficients
  ALLOCATE( vbcf(maxvbc,maxtvsch,ninstrument) )
  IF( vbc_coef_fname /= '' ) THEN
    CAll vbc_read(trim(vbc_coef_fname), vbcf, maxvbc, maxtvsch, &
                       num_satellite_iasi, tvsinst, tvsch, ntvsch)
  ELSE
    WRITE(ADM_LOG_FID,*) 'Caution! Airmass bias coefficient file is not appointed'
    vbcf(:,:,:)=0.0d0
  END IF
  
  !DO nn = 1, num_satellite_iasi
  !DO ic = 1, ntvsch(nn)
  !  WRITE(ADM_LOG_FID,'(2i3,8F8.3)') nn, tvsch(ic,nn), vbcf(:,ic,nn)
  !END DO
  !END DO
!
end SUBROUTINE obsope_iasi_read
!------------------------------------------------------------------------------
SUBROUTINE interpolate_iasi(nn)
  USE mod_adm
  USE mod_cnst
  USE mod_obsope_common
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: i, l, p, k, n
  INTEGER :: ierr, nv
  INTEGER :: tmp_id(nv3d)
  REAL(8) :: fac1, fac2, fac3
  REAL(8) :: fac_sum

  iasi(nn)%obsdata_3d_tmp=0.0
  iasi(nn)%obsdata_3d=0.0
  iasi(nn)%obsdata_2d_tmp=0.0
  iasi(nn)%obsdata_2d=0.0
  iasi(nn)%lsql_model_tmp=0.0
  iasi(nn)%lsql_model=0.0

  tmp_id(1)=id_temp_nicam
  tmp_id(2)=id_qvap_nicam
  tmp_id(3)=id_pres_nicam
  !WRITE(ADM_LOG_FID,*) 'interpolate_iasi'
  DO i = 1, nobs_iasi(nn)
    !WRITE(ADM_LOG_FID,'(6f12.5)') iasi(nn)%lon(i), iasi(nn)%lat(i),   &
    !                              iasi(nn)%saza(i), iasi(nn)%saaz(i), &
    !                              iasi(nn)%soza(i), iasi(nn)%soaz(i) 
    IF( iasi(nn)%inprc(i) ) THEN
      l    = iasi(nn)%l_index(i)
      fac1 = iasi(nn)%w1(i)
      fac2 = iasi(nn)%w2(i)
      fac3 = iasi(nn)%w3(i)
      fac_sum = fac1 + fac2 + fac3
      write(ADM_LOG_FID,'(i5,4f12.5)') l, fac1, fac2, fac3, fac_sum
      write(ADM_LOG_FID,'(3i5)') iasi(nn)%n1_index(i), iasi(nn)%n2_index(i), iasi(nn)%n3_index(i)
      DO nv = 1, nv3d
        DO k = 1, ADM_vlayer
          iasi(nn)%obsdata_3d_tmp(k,i,nv) = &
             ( fac1 * icodata4_3d(iasi(nn)%n1_index(i),k,l,tmp_id(nv)) &
             + fac2 * icodata4_3d(iasi(nn)%n2_index(i),k,l,tmp_id(nv)) &
             + fac3 * icodata4_3d(iasi(nn)%n3_index(i),k,l,tmp_id(nv)) &
           ) / fac_sum
        END DO ! k
        IF(tmp_id(nv) == id_pres_nicam) THEN
          iasi(nn)%obsdata_3d_tmp(:,i,nv) = EXP(iasi(nn)%obsdata_3d_tmp(:,i,nv))
        END IF
        IF(tmp_id(nv) == id_qvap_nicam) THEN
          iasi(nn)%obsdata_3d_tmp(:,i,nv) = iasi(nn)%obsdata_3d_tmp(:,i,nv) * q2ppmv
        END IF
      END DO
    END IF
  END DO

  DO i = 1, nobs_iasi(nn)
    IF( iasi(nn)%inprc(i) ) THEN
      l    = iasi(nn)%l_index(i)
      fac1 = iasi(nn)%w1(i)
      fac2 = iasi(nn)%w2(i)
      fac3 = iasi(nn)%w3(i)
      fac_sum = fac1 + fac2 + fac3
      DO nv = 1, nv2d
        iasi(nn)%obsdata_2d_tmp(i,nv) = &
                ( fac1 * icodata4_2d(iasi(nn)%n1_index(i),1,l,nv) &
                + fac2 * icodata4_2d(iasi(nn)%n2_index(i),1,l,nv) &
                + fac3 * icodata4_2d(iasi(nn)%n3_index(i),1,l,nv) &
              ) / fac_sum
        IF(nv == id_surp_nicam) THEN
          iasi(nn)%obsdata_2d_tmp(i,nv) = EXP(iasi(nn)%obsdata_2d_tmp(i,nv))
        END IF
        IF(nv == id_qv2m_nicam) THEN
          iasi(nn)%obsdata_2d_tmp(i,nv) = iasi(nn)%obsdata_2d_tmp(i,nv) * q2ppmv
        END IF
      END DO

      iasi(nn)%lsql_model_tmp(i) = &
                ( fac1 * icoland(iasi(nn)%n1_index(i),l) &
                + fac2 * icoland(iasi(nn)%n2_index(i),l) &
                + fac3 * icoland(iasi(nn)%n3_index(i),l) &
              ) / fac_sum

    END IF
  END DO

  CALL MPI_ALLREDUCE( iasi(nn)%obsdata_3d_tmp, iasi(nn)%obsdata_3d, &
                      ADM_vlayer*nobs_iasi(nn)*nv3d,                 &
                      MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (iasi(nn)%obsdata_3d)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

  CALL MPI_ALLREDUCE( iasi(nn)%obsdata_2d_tmp, iasi(nn)%obsdata_2d, &
                      nobs_iasi(nn)*nv2d,                            &
                      MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (iasi(nn)%obsdata_2d)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

  CALL MPI_ALLREDUCE( iasi(nn)%lsql_model_tmp, iasi(nn)%lsql_model, &
                      nobs_iasi(nn),                                 &
                      MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  iasi(nn)%lsql(:)=INT(iasi(nn)%lsql_model(:))

  !DO i = 1, nobs_iasi(nn)
  !  WRITE(ADM_LOG_FID,'(i7,3f12.5)') i, (maxval(iasi(nn)%obsdata_3d(:,i,nv)),nv=1,3)
  !END DO

  !DO i = 1, nobs_iasi(nn)
  !  WRITE(ADM_LOG_FID,'(i7,7f12.5)') i, (iasi(nn)%obsdata_2d(i,nv),nv=1,7)
  !END DO

  !DEALLOCATE( iasi(nn)%obsdata_3d_tmp )
  !DEALLOCATE( iasi(nn)%obsdata_2d_tmp )

  !DEALLOCATE( inprc )
  !DEALLOCATE( l_index )
  !DEALLOCATE( n1_index )
  !DEALLOCATE( n2_index )
  !DEALLOCATE( n3_index )
  !DEALLOCATE( w1 )
  !DEALLOCATE( w2 )
  !DEALLOCATE( w3 )

END SUBROUTINE interpolate_iasi
!------------------------------------------------------------------------------
SUBROUTINE calc_radiance_iasi
  USE mod_adm
  USE mod_obsope_common
  IMPLICIT NONE
  INTEGER :: prc_per_inst
  INTEGER :: nn, n, ic
  INTEGER :: obsint
  INTEGER :: sobs, eobs
  INTEGER :: tmp_id(nv3d)
  INTEGER :: ierr
  INTEGER :: pe_in_satellite

  INTEGER :: iz
  REAL(4) :: ri, rj, rk
  REAL(4) :: tmp_lon

  tmp_id(1)=id_temp_nicam
  tmp_id(2)=id_qvap_nicam
  tmp_id(3)=id_pres_nicam

  !WRITE(ADM_LOG_FID,*) 'calc_radiance'
  !DO nn = 1, num_satellite_iasi
  !DO n = 1, nobs_iasi(nn)
  !  WRITE(ADM_LOG_FID,'(6f12.5)') iasi(nn)%lon(n), iasi(nn)%lat(n),   &
  !                                iasi(nn)%saza(n), iasi(nn)%saaz(n), &
  !                                iasi(nn)%soza(n), iasi(nn)%soaz(n)
  !END DO
  !END DO

  DO nn = 1, num_satellite_iasi
    iasi(nn)%bt_tmp(:,:)=0.0
    iasi(nn)%bt(:,:)=0.0
    iasi(nn)%trans_tmp(:,:,:)=0.0
    iasi(nn)%trans(:,:,:)=0.0
  END DO

  prc_per_inst=ADM_prc_all/num_satellite_iasi
  nn=(ADM_prc_me-1)/prc_per_inst+1
  pe_in_satellite=ADM_prc_me - (nn-1)*prc_per_inst

  IF( nobs_iasi(nn) /= 0 ) THEN
    obsint=nobs_iasi(nn)/prc_per_inst
    IF(nobs_iasi(nn)<prc_per_inst) obsint=1
    sobs=obsint*(pe_in_satellite-1)+1
    eobs=obsint*pe_in_satellite
    IF(ADM_prc_me==prc_per_inst*nn) eobs=nobs_iasi(nn)
    IF(nobs_iasi(nn) >= prc_per_inst .OR. &
       nobs_iasi(nn) <  prc_per_inst .AND. pe_in_satellite <= nobs_iasi(nn) ) THEN

      DO n = sobs, eobs
        tmp_lon=iasi(nn)%lon(n)
        IF(tmp_lon < 0.0) tmp_lon=tmp_lon+360.0
        DO iz = 1, ADM_vlayer
          CALL phys2ijk(nlon_co2, nlat_co2, nlev_co2, pres_co2, 0.0, lon_co2, lat_co2, &
               tmp_lon, iasi(nn)%lat(n), iasi(nn)%obsdata_3d(iz,n,3),                  &
               ri, rj, rk )
          CALL itpl_3d(nlon_co2, nlat_co2, nlev_co2, co2_org, ri, rj, rk, iasi(nn)%co2(iz, n) )
          WRITE(ADM_LOG_FID,*) ri, rj, rk, iasi(nn)%obsdata_3d(iz,n,3), iasi(nn)%co2(iz, n)
          FLUSH(ADM_LOG_FID)
        END DO
      END DO
     
      CALL iasi_fwd( ADM_vlayer, eobs-sobs+1, rttovcoef_fname(nn),              &
          iasi(nn)%datetime(1:3,sobs:eobs),                                     &
          DBLE(iasi(nn)%co2(ADM_vlayer:1:-1, sobs:eobs)),                       &
          DBLE(iasi(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 3)),             &
          DBLE(iasi(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 1)),             &
          DBLE(iasi(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 2)),             &
          DBLE(iasi(nn)%obsdata_2d(                 sobs:eobs, id_tsfc_nicam)), &
          DBLE(iasi(nn)%obsdata_2d(                 sobs:eobs, id_qv2m_nicam)), &
          DBLE(iasi(nn)%obsdata_2d(                 sobs:eobs, id_surp_nicam)), &
          DBLE(iasi(nn)%obsdata_2d(                 sobs:eobs, id_u10m_nicam)), &
          DBLE(iasi(nn)%obsdata_2d(                 sobs:eobs, id_v10m_nicam)), &
          DBLE(iasi(nn)%soza       ( sobs:eobs )),                              &
          DBLE(iasi(nn)%soaz       ( sobs:eobs )),                              &
          DBLE(iasi(nn)%saza       ( sobs:eobs )),                              &
          DBLE(iasi(nn)%saaz       ( sobs:eobs )),                              &
          DBLE(iasi(nn)%elev       ( sobs:eobs )),                              &
          DBLE(iasi(nn)%lon        ( sobs:eobs )),                              &
          DBLE(iasi(nn)%lat        ( sobs:eobs )),                              &
          DBLE(iasi(nn)%lsql       ( sobs:eobs )),                              &
          iasi(nn)%bt_tmp          ( :,sobs:eobs ),                             &
          iasi(nn)%trans_tmp       ( :,:,sobs:eobs ),                           &
          ntvsch(nn), tvsch(:,nn), time_rttov_iasi )

    END IF
  END IF

  DO nn = 1, num_satellite_iasi
    CALL MPI_ALLREDUCE(iasi(nn)%bt_tmp, iasi(nn)%bt, &
                       ntvsch(nn)*nobs_iasi(nn),              &
                       MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(iasi(nn)%trans_tmp, iasi(nn)%trans, &
                       ADM_vlayer*ntvsch(nn)*nobs_iasi(nn),         &
                       MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  END DO 

  DO nn = 1, num_satellite_iasi
  DO n = 1, nobs_iasi(nn)
  DO ic = 1, ntvsch(nn)
    iasi(nn)%ominusb(ic,n)=iasi(nn)%odat(ic,n)-iasi(nn)%bt(ic,n)
  END DO
  END DO
  END DO

END SUBROUTINE calc_radiance_iasi
!------------------------------------------------------------------------------
SUBROUTINE calc_vbc_iasi(nn)
  USE mod_obsope_common
  USE mod_adm
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n, ic, k
  INTEGER :: ifoot, ichan
  REAL(8) :: iwlr

  DO n = 1, nobs_iasi(nn)
    iasi(nn)%vbc_pred(1,:,n)=-999.9 ! IWLR(1000-200 hPa)
    iasi(nn)%vbc_pred(2,:,n)=-999.9 ! IWLR( 200- 50 hPa)
    iasi(nn)%vbc_pred(3,:,n)=-999.9
    iasi(nn)%vbc_pred(4,:,n)=0.0d0
    iasi(nn)%vbc_pred(5,:,n)=0.0d0
    iasi(nn)%vbc_pred(6,:,n)=0.0d0
    iasi(nn)%vbc_pred(7,:,n)=0.0d0
    iasi(nn)%vbc_pred(8,:,n)=0.0d0
    !iasi(nn)%vbc_pred(4,:,n)=(iasi(nn)%obsdata_2d(n,id_tsfc_nicam)-273.15d0)/10.0d0
    !iasi(nn)%vbc_pred(5,:,n)=0.0d0
    !iasi(nn)%vbc_pred(6,:,n)=iasi(nn)%obsdata_2d(n,id_cldw_nicam)/30.0d0
    !iasi(nn)%vbc_pred(7,:,n)=1.0d0/COS(iasi(nn)%saza(n)*deg2rad)
    !iasi(nn)%vbc_pred(8,:,n)=1.0d0

    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(iasi(nn)%obsdata_3d(k,n,3) >  200.0d0 .AND. &
         iasi(nn)%obsdata_3d(k,n,3) < 1000.0d0 ) THEN
         iwlr = iwlr + &
           ( iasi(nn)%obsdata_3d(k+1,n,1) - iasi(nn)%obsdata_3d(k,n,1) ) * &
           ( iasi(nn)%trans(ADM_vlayer-k,   ic, n) -              &
             iasi(nn)%trans(ADM_vlayer-k+1, ic, n))
      END IF
      END DO
      iasi(nn)%vbc_pred(1,ic,n)=iwlr
    END DO
    
    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(iasi(nn)%obsdata_3d(k,n,3) >   50.0d0 .AND. &
         iasi(nn)%obsdata_3d(k,n,3) <  200.0d0 ) THEN
         iwlr = iwlr + &
           ( iasi(nn)%obsdata_3d(k+1,n,1) - iasi(nn)%obsdata_3d(k,n,1) ) * &
           ( iasi(nn)%trans(ADM_vlayer-k,   ic, n) - &
             iasi(nn)%trans(ADM_vlayer-k+1, ic, n))
      END IF
      END DO
      iasi(nn)%vbc_pred(2,ic,n)=iwlr
    END DO

    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(iasi(nn)%obsdata_3d(k,n,3) >    5.0d0 .AND. &
         iasi(nn)%obsdata_3d(k,n,3) <   50.0d0 ) THEN
         iwlr = iwlr + &
           ( iasi(nn)%obsdata_3d(k+1,n,1) - iasi(nn)%obsdata_3d(k,n,1) ) * &
           ( iasi(nn)%trans(ADM_vlayer-k,   ic, n) - &
             iasi(nn)%trans(ADM_vlayer-k+1, ic, n))
      END IF
      END DO
      iasi(nn)%vbc_pred(3,ic,n)=iwlr
    END DO

    WRITE(ADM_LOG_FID,*) 'AIRMASS BIAS'
    DO ic = 1, ntvsch(nn)
      iasi(nn)%airmass_bias(:,ic,n)=iasi(nn)%vbc_pred(:,ic,n)*vbcf(:,ic,nn)
      WRITE(ADM_LOG_FID,'(2i5,8f10.5)') nn, ic, iasi(nn)%airmass_bias(:,ic,n)
    END DO

  END DO

  ! Correct Bias
  WRITE(ADM_LOG_FID,*) 'BIAS CORRECTION'
  iasi(nn)%odat_bc(:,:)=0.0
  DO n = 1, nobs_iasi(nn)
    ifoot=iasi(nn)%fov(n)
    DO ic = 1, ntvsch(nn)
      iasi(nn)%odat_bc(ic,n) = iasi(nn)%odat(ic,n)                & ! y
                             - vbcf_scan(ifoot,ic,nn)             & ! scan bias
                             + sum(iasi(nn)%airmass_bias(:,ic,n))   ! airmass bias
                             !- iasi(nn)%vbcf_scan(ifoot,ic)       & ! scan bias
      !WRITE(ADM_LOG_FID,'(2i5,4f12.7)') &
      !              n, ic, iasi(nn)%odat(ichan,n),                   &
      !              iasi(nn)%odat_bc(ichan,n),                       &
      !              vbcf_scan(ifoot,ic,nn), sum(iasi(nn)%airmass_bias(:,ic,n))
      !              !iasi(nn)%odat(ichan,n)-iasi(nn)%odat_bc(ichan,n)
    END DO
  END DO

END SUBROUTINE calc_vbc_iasi
!------------------------------------------------------------------------------
SUBROUTINE quality_control_iasi(nn)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n

  DO n = 1, nobs_iasi(nn)
    IF( iasi(nn)%lsql(n) == 0 ) THEN ! OVER LAND
      iasi(nn)%qc(:,n)=0
    END IF
  END DO

END SUBROUTINE quality_control_iasi
!------------------------------------------------------------------------------
SUBROUTINE update_vbc_iasi(imem, nn)
  USE mod_adm
  USE mod_vbc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  INTEGER, INTENT(IN) :: nn
  INTEGER :: i, ic, ichan
  INTEGER :: ntvschan( ntvsch(nn) )

  IF( imem > 0 ) RETURN

  DO i = 1, nobs_iasi(nn)
  DO ic = 1, ntvsch(nn)
    IF( ABS( iasi(nn)%odat_bc(ic,i) - iasi(nn)%bt(ic,i) ) > iasi(nn)%err(ic,i)*3.0 ) THEN
      iasi(nn)%qc(ic,i)=0
    END IF
  END DO
  END DO

  ALLOCATE( vbca(maxvbc,maxtvsch,num_satellite_iasi) )
  !WRITE(ADM_LOG_FID,*) 'vbcf', nn
  !DO ic = 1, ntvsch(nn)
  !  WRITE(ADM_LOG_FID,'(8f10.6)') (vbcf(i,ic,nn),i=1,maxvbc)
  !END DO
  !FLUSH(ADM_LOG_FID)

  CALL das_vbc( nobs_iasi(nn), maxvbc, ntvsch(nn),    &
                tvsname(nn), tvsch(:,nn), ntvsch(nn), &
                iasi(nn)%odat_bc(1:ntvsch(nn),:),     &
                iasi(nn)%bt(     1:ntvsch(nn),:),     &
                iasi(nn)%vbc_pred(:,1:ntvsch(nn),:),  &
                vbcf(:,:,nn), vbca(:,:,nn),           &
                iasi(nn)%qc(     1:ntvsch(nn),:),     &
                iasi(nn)%err(    1:ntvsch(nn),:),     &
                ntvschan )

  !WRITE(ADM_LOG_FID,*) 'vbca', nn
  !DO ic = 1, ntvsch(nn)
  !  WRITE(ADM_LOG_FID,'(8f10.6)') (vbca(i,ic,nn),i=1,maxvbc)
  !END DO
  !FLUSH(ADM_LOG_FID)

  vbc_coef_out_fname=TRIM(vbc_coef_out_fname)//'_'//TRIM(tvsname(nn))//'_iasi'
  CALL vbc_write(trim(vbc_coef_out_fname), vbca(:,:,nn), maxvbc, maxtvsch, &
                 tvsinst(:,nn), tvsch(:,nn), ntvsch(nn), ntvschan)

END SUBROUTINE update_vbc_iasi
!------------------------------------------------------------------------------
SUBROUTINE update_scanbias_iasi(nn)
  USE mod_adm
  USE mod_obsope_common, only: nslot, nbv
  USE mod_scanbias

  INTEGER, INTENT(IN) :: nn
  CHARACTER(2)   :: cslot
  CHARACTER(4)   :: cfile
  CHARACTER(6)   :: cimem
  CHARACTER(256) :: fname
  LOGICAL        :: ex

  CHARACTER(4)   :: tvsname_scan
  INTEGER        :: ifoot
  REAL(4)        :: lwp
  REAL(4)        :: dum
  REAL(4)        :: scanread(ntvsch(nn))

  INTEGER :: num_scan     (nfootp(nn),         ntvsch(nn))
  REAL(4) :: vbca_scan_tmp(nfootp(nn), 100000, ntvsch(nn))
  REAL(4) :: vbca_scan_ave(nfootp(nn),         ntvsch(nn))
  REAL(4) :: vbca_scan    (nfootp(nn), 100000, ntvsch(nn))

  INTEGER :: is, imem, ios, ic

  num_scan(:,:)=0
  vbca_scan_tmp(:,:,:)=0.0
  vbca_scan_ave(:,:)=0.0

  DO is = 1, nslot
    WRITE(cslot,'(I2.2)') is
    DO imem = 1, nbv
      WRITE(cimem,'(I6.6)') imem
      fname=TRIM(scanbias_est_ifname)//'_'//TRIM(tvsname(nn))//TRIM(cslot)//TRIM(cimem)//'.txt'
      WRITE(ADM_LOG_FID,*) TRIM(fname)
      INQUIRE(FILE=TRIM(fname), EXIST=ex)
      IF(ex) THEN
        OPEN(301,FILE=TRIM(fname))
        DO
          READ(301,*,IOSTAT=ios) tvsname_scan, ifoot, lwp, dum, scanread(1:ntvsch(nn))
          !WRITE(ADM_LOG_FID,'(A,i4,10F10.3)') tvsname_scan, ifoot, lwp, dum, scanread(1:ntvsch(nn))
          IF(ios /= 0) THEN
            EXIT
          END IF
          DO ic = 1, ntvsch(nn)
            IF( ABS(scanread(ic)) < 10.0 ) THEN
              num_scan(ifoot,ic) = num_scan(ifoot,ic) + 1
              vbca_scan_tmp(ifoot, num_scan(ifoot,ic), ic) = scanread(ic)
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO

  DO ic = 1, ntvsch(nn)
  DO ifoot = 1, nfootp(nn)
    IF( num_scan(ifoot,ic) > 0 ) THEN
      vbca_scan_ave(ifoot,ic) = &
         SUM(vbca_scan_tmp(ifoot, 1:num_scan(ifoot,ic), ic)) / REAL(num_scan(ifoot,ic))
    END IF
  END DO
  END DO

  DO ic = 1, ntvsch(nn)
  DO ifoot = 1, nfootp(nn)
    IF( num_scan(ifoot,ic) > 0 ) THEN
      vbcf_scan(ifoot,ic,nn)=0.97*vbcf_scan(ifoot,ic,nn)+&
                             0.03*(vbca_scan_ave(ifoot,ic)-vbcf_scan(ifoot,ic,nn))
    ELSE
      vbcf_scan(ifoot,ic,nn) = vbcf_scan(ifoot,ic,nn)
    END IF
  END DO
  END DO

  scanbias_est_ofname=TRIM(scanbias_est_ofname)//'_'//TRIM(tvsname(nn))//'_iasi'
  CALL vbc_scan_write(trim(scanbias_est_ofname), vbcf_scan(:,1:ntvsch(nn),nn), &
                 nfootp(nn), ntvsch(nn), &
                 tvsinst(:,nn), tvsch(1:ntvsch(nn),nn), ntvsch(nn))

END SUBROUTINE update_scanbias_iasi
!-------------------------------------------------------------------------------
SUBROUTINE output_iasi(imem, nn)
  USE mod_adm
  USE mod_obsope_common
  USE mod_scanbias
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  INTEGER, INTENT(IN) :: nn
  CHARACTER(6)        :: cimem
  CHARACTER(256)      :: fname
  INTEGER :: i

  IF( imem < 0 ) RETURN

  WHERE( iasi(nn)%lon(:) < 0.0 ) iasi(nn)%lon(:) = iasi(nn)%lon(:) + 360.0

  WRITE(cimem(1:6),'(I6.6)') imem + 1
  fname=TRIM(output_dirname_iasi)//'/'//TRIM(tvsname(nn))//&
        TRIM(output_basename_iasi)//TRIM(cimem)//'.dat'
  OPEN(1, FILE=TRIM(fname),FORM='unformatted',ACCESS='sequential')
  DO i = 1, nobs_iasi(nn)
    WRITE(1) REAL(id_bt_obs),        REAL(iasi(nn)%lsql(i)),         &
           & REAL(iasi(nn)%lon(i)), REAL(iasi(nn)%lat(i)),          &
           & REAL(iasi(nn)%saza(i)),                                 &
           & REAL(iasi(nn)%obsdata_2d(i,id_tsfc_nicam)),             &
           & REAL(0.0), REAL(iasi(nn)%fov(i)),                       &
           & REAL(iasi(nn)%obsdata_3d(iasi(nn)%weight_maxlev(&
                  1:ntvsch(nn),i), i, 3)), &
           & REAL(iasi(nn)%odat_bc(1:ntvsch(nn),i)),       &
           & REAL(iasi(nn)%err    (1:ntvsch(nn),i)),       &
           & REAL(iasi(nn)%bt     (1:ntvsch(nn),i)),       &
           & REAL(iasi(nn)%qc     (1:ntvsch(nn),i))

  END DO
  CLOSE(1)

  IF(output_text) THEN
    fname=TRIM(output_dirname_iasi)//'/'//TRIM(tvsname(nn))//&
          TRIM(output_basename_iasi)//TRIM(cimem)//'.txt'
    OPEN(2, FILE=TRIM(fname),FORM='formatted')
    DO i = 1, nobs_iasi(nn)
      WRITE(2,'(40F8.2)')                                               &
             & REAL(iasi(nn)%lwp(i)), REAL(iasi(nn)%lsql(i)),         &
             & REAL(iasi(nn)%lsql_model(i)),                           &
             & REAL(iasi(nn)%lon(i)), REAL(iasi(nn)%lat(i)),          &
             & REAL(iasi(nn)%saza(i)),                                 &
             & REAL(iasi(nn)%obsdata_2d(i,id_tsfc_nicam)),             &
             & REAL(0.0), REAL(iasi(nn)%fov(i)),                       &
             & REAL(iasi(nn)%obsdata_3d(iasi(nn)%weight_maxlev(       &
                    1:ntvsch(nn),i), i, 3)), &
             & REAL(iasi(nn)%odat   (1:ntvsch(nn),i)),       &
             & REAL(iasi(nn)%odat_bc(1:ntvsch(nn),i)),       &
             & REAL(iasi(nn)%err    (1:ntvsch(nn),i)),       &
             & REAL(iasi(nn)%bt     (1:ntvsch(nn),i)),       &
             & REAL(iasi(nn)%qc     (1:ntvsch(nn),i))
    END DO
    CLOSE(2)
  END IF

  CALL tvs_ominusb_output(nn, imem, ntvsch(nn), nobs_iasi(nn),   &
                          iasi(nn)%ominusb(:,:),                 &
                          iasi(nn)%fov(:),                       &
                          iasi(nn)%lsql(:),                      &
                          iasi(nn)%lwp(:),                       &
                          TRIM(ominusb_fname_iasi), tvsname(nn), &
                          islot)

END SUBROUTINE output_iasi
!------------------------------------------------------------------------------
SUBROUTINE set_instrument
  tvsch = 0
  !
  ! METOP-1 IASI
  !
  tvsname(1)   = 'MI01'
  tvsinst(1,1) = rttv_plat_metop
  tvsinst(2,1) = 01
  tvsinst(3,1) = rttv_inst_iasi
  tvsch(1,1)   = 212
  tvsch(2,1)   = 246
  tvsch(3,1)   = 262
  tvsch(4,1)   = 275
  !tvsch(1,1)   = 148
  !tvsch(2,1)   = 187
  !tvsch(3,1)   = 193
  !tvsch(4,1)   = 212
  ntvsch(1)    =   4
  nfootp(1)    = 30
  !
  ! METOP-2 IASI
  !
  tvsname(2)   = 'MI02'
  tvsinst(1,2) = rttv_plat_metop
  tvsinst(2,2) = 02
  tvsinst(3,2) = rttv_inst_iasi
  tvsch(1,2)   = 212
  tvsch(2,2)   = 246
  tvsch(3,2)   = 262
  tvsch(4,2)   = 275
  ntvsch(2)    = 4
  nfootp(2)    = 30
  RETURN

END SUBROUTINE set_instrument
!-----------------------------------------------------------------------
END MODULE mod_obsope_iasi
