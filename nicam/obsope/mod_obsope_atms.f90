MODULE mod_obsope_atms
  USE mod_obsope_common, ONLY: flush_text
  USE common_tvs_nicam, ONLY: &
    rttv_plat_jpss,           &
    rttv_inst_atms,           &
    rttv_chan_atms 

  IMPLICIT NONE
  PUBLIC

  INTEGER, PRIVATE, PARAMETER   :: ninstrument=1
  INTEGER, PRIVATE, PARAMETER   :: maxtvsch=22
  INTEGER, PRIVATE, PARAMETER   :: maxvbc=8
  INTEGER, PRIVATE, PARAMETER   :: maxfoot=32
  INTEGER, PRIVATE, SAVE        :: ntvs
  INTEGER, PRIVATE, SAVE        :: tvsinst(3,ninstrument)
  INTEGER, PRIVATE, SAVE        :: tvsch(maxtvsch,ninstrument)
  INTEGER, PRIVATE, SAVE        :: ntvsch(ninstrument)
  INTEGER, PRIVATE, SAVE        :: nfootp(ninstrument)
  CHARACTER(4), PRIVATE, SAVE   :: tvsname(ninstrument)
  CHARACTER(256), PRIVATE, SAVE :: rttovcoef_fname(ninstrument) = ''
  CHARACTER(256), PRIVATE       :: input_fname_atms = ''
  CHARACTER(256), PRIVATE       :: output_dirname_atms = ''
  CHARACTER(2), PRIVATE         :: output_basename_atms = ''
  CHARACTER(256), PRIVATE       :: ominusb_fname_atms   = ''
  CHARACTER(256), PRIVATE       :: scanbias_fname = ''
  CHARACTER(256), PRIVATE       :: vbc_coef_fname = ''
  CHARACTER(256), PRIVATE       :: vbc_coef_out_fname = ''
  CHARACTER(256), PRIVATE       :: scanbias_est_ifname = ''
  CHARACTER(256), PRIVATE       :: scanbias_est_ofname = ''
  LOGICAL, PRIVATE, SAVE        :: output_text = .false.

  INTEGER, PRIVATE, PARAMETER :: nv3d=3
  INTEGER, PRIVATE, PARAMETER :: nv2d=7

  INTEGER :: num_satellite_atms

  TYPE atms_base
    INTEGER, ALLOCATABLE :: said(:)
    REAL(4), ALLOCATABLE :: lat(:)
    REAL(4), ALLOCATABLE :: lon(:)
    REAL(4), ALLOCATABLE :: elev(:)
    REAL(4), ALLOCATABLE :: elev_model(:)
    REAL(4), ALLOCATABLE :: elev_model_tmp(:)
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
  END TYPE atms_base
    
  TYPE(atms_base), allocatable :: atms(:)

  INTEGER, ALLOCATABLE :: nobs_atms(:)

  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbcf_scan(:,:,:)
  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbcf(:,:,:)
  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbca(:,:,:)

  PRIVATE :: set_instrument

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE obsope_atms_init
  USE mod_adm
  IMPLICIT NONE
  INTEGER :: ierr

  NAMELIST / atms_cnf /   &
    input_fname_atms,     &
    output_dirname_atms,  &
    output_basename_atms, &
    rttovcoef_fname,       &
    ominusb_fname_atms,   &
    scanbias_fname,        &
    vbc_coef_fname,        &
    vbc_coef_out_fname,    &
    scanbias_est_ifname,   &
    scanbias_est_ofname,   &
    output_text

  CALL set_instrument

  OPEN(121, FILE='obsope.cnf')
  READ(121, NML=atms_cnf)
  CLOSE(121)

END SUBROUTINE obsope_atms_init
!------------------------------------------------------------------------------
SUBROUTINE obsope_atms_read
  USE mod_adm
  USE mod_obsope_common
  USE mod_scanbias
  USE mod_vbc
  IMPLICIT NONE
  INTEGER :: ierr
  INTEGER :: n, nn
  INTEGER :: ic
  REAL(4) :: cos_tmp

  OPEN(122,file=trim(input_fname_atms),  form='unformatted', &
       access='sequential', status='old', iostat=ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in opening the file', trim(input_fname_atms)
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    IF(flush_text) FLUSH(ADM_LOG_FID)
    CALL ADM_proc_stop
  END IF
  READ(122) num_satellite_atms
  ALLOCATE( atms(num_satellite_atms) )
  ALLOCATE( nobs_atms(num_satellite_atms) )

  READ(122) ntvs
  READ(122) nobs_atms(:)
  !WRITE(ADM_LOG_FID,*) nobs_atms

  DO nn = 1, num_satellite_atms
    ALLOCATE( atms(nn)%said          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%lon           ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%lat           ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%elev          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%elev_model    ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%elev_model_tmp( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%odat          ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%odat_bc       ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%err           ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%saza          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%soza          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%soaz          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%saaz          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%fov           ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%lsql          ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%qc_tmp        ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%qc            ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%obsdata_3d_tmp( ADM_vlayer,nobs_atms(nn),nv3d) )
    ALLOCATE( atms(nn)%obsdata_3d    ( ADM_vlayer,nobs_atms(nn),nv3d) )
    ALLOCATE( atms(nn)%obsdata_2d_tmp( nobs_atms(nn),nv2d) )
    ALLOCATE( atms(nn)%obsdata_2d    ( nobs_atms(nn),nv2d) )
    ALLOCATE( atms(nn)%lwp           ( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%bt_tmp        ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%bt            ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%ominusb       ( maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%trans_tmp     ( ADM_vlayer, maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%trans         ( ADM_vlayer, maxtvsch,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%wk            ( 48,nobs_atms(nn)) )
    ALLOCATE( atms(nn)%inprc         ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%l_index       ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%n1_index      ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%n2_index      ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%n3_index      ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%w1            ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%w2            ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%w3            ( nobs_atms(nn) ))
    ALLOCATE( atms(nn)%vbcf_scan     ( maxfoot, maxtvsch ))
    ALLOCATE( atms(nn)%weight        ( ADM_vlayer, maxtvsch, nobs_atms(nn) ))
    ALLOCATE( atms(nn)%weight_maxlev ( maxtvsch, nobs_atms(nn) ))
    ALLOCATE( atms(nn)%vbc_pred      ( maxvbc, maxtvsch, nobs_atms(nn) ))
    ALLOCATE( atms(nn)%airmass_bias  ( maxvbc, maxtvsch, nobs_atms(nn) ))
    ALLOCATE( atms(nn)%lsql_model_tmp( nobs_atms(nn)) )
    ALLOCATE( atms(nn)%lsql_model    ( nobs_atms(nn)) )
  END DO

  DO nn = 1, num_satellite_atms
  DO n = 1, nobs_atms(nn)
    READ(122) atms(nn)%wk(:,n)
  END DO
  END DO

  DO nn = 1, num_satellite_atms
  !WRITE(ADM_LOG_FID,*) 'Satellite ', nn
  DO n = 1, nobs_atms(nn)
    atms(nn)%lat(n)    = atms(nn)%wk(11,n)
    atms(nn)%lon(n)    = atms(nn)%wk(12,n)
    atms(nn)%said(n)   = atms(nn)%wk( 7,n)
    atms(nn)%fov(n)    = atms(nn)%wk(10,n)
    !atms(nn)%lsql(n)   = atms(nn)%wk(12,n)
    atms(nn)%saza(n)   = atms(nn)%wk(21,n)
    atms(nn)%soza(n)   = atms(nn)%wk(23,n)
    !atms(nn)%elev(n)   = atms(nn)%wk(15,n)/1000.0
    atms(nn)%soaz(n)   = atms(nn)%wk(24,n)
    atms(nn)%saaz(n)   = atms(nn)%wk(22,n)
    atms(nn)%odat(:,n) = atms(nn)%wk(26:47,n)
    !WRITE(ADM_LOG_FID,'(6f12.5)') atms(nn)%lon(n), atms(nn)%lat(n),   &
    !                              atms(nn)%saza(n), atms(nn)%saaz(n), &
    !                              atms(nn)%soza(n), atms(nn)%soaz(n) 
  END DO
  !WRITE(ADM_LOG_FID,*) minval(atms(nn)%lat(:)), maxval(atms(nn)%lat(:))
  !WRITE(ADM_LOG_FID,*) minval(atms(nn)%lon(:)), maxval(atms(nn)%lon(:))
  !WRITE(ADM_LOG_FID,*) minval(atms(nn)%saza(:)), maxval(atms(nn)%saza(:))
  END DO

  DO nn = 1, num_satellite_atms
    atms(nn)%elev(n) = 0.0
    atms(nn)%err(:,:)= 0.35 ! Following Bormann et al. (2012)
    atms(nn)%qc(:,:) = 1
  END DO
  ! Liquid Water Path (LWP: g/kg)
  DO nn = 1, num_satellite_atms
  DO n = 1, nobs_atms(nn)
    cos_tmp=COS(atms(nn)%saza(n)*deg2rad)
    IF( atms(nn)%odat(1,n) < 285.0 .AND. &
        atms(nn)%odat(2,n) < 285.0 ) THEN
      atms(nn)%lwp(n) = cos_tmp*(&
                         8.24 - ( 2.539 - 1.744*cos_tmp ) * cos_tmp + &
                         0.754 * log( 285.0 - atms(nn)%odat(1,n)) - &
                         2.265 * log( 285.0 - atms(nn)%odat(2,n)))
    ELSE
      atms(nn)%lwp(n) = 999.0
    END IF
  END DO
  END DO

  ! Read scan bias
  ALLOCATE( vbcf_scan(maxfoot,maxtvsch,num_satellite_atms) )
  IF( scanbias_fname /= '' ) THEN
    CAll vbc_scan_read(trim(scanbias_fname), vbcf_scan, maxfoot, maxtvsch, &
                       num_satellite_atms, tvsinst, tvsch, ntvsch)
  ELSE
    WRITE(ADM_LOG_FID,*) 'Caution! scan bias file is not appointed'
    vbcf_scan(:,:,:)=0.0d0
  END IF

  DO nn = 1, num_satellite_atms
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(2i3,32F8.3)') nn, tvsch(ic,nn), vbcf_scan(:,ic,nn)
  END DO
  END DO

  ! Read varBC coefficients
  ALLOCATE( vbcf(maxvbc,maxtvsch,num_satellite_atms) )
  IF( vbc_coef_fname /= '' ) THEN
    CAll vbc_read(trim(vbc_coef_fname), vbcf, maxvbc, maxtvsch, &
                       num_satellite_atms, tvsinst, tvsch, ntvsch)
  ELSE
    WRITE(ADM_LOG_FID,*) 'Caution! Airmass bias coefficient file is not appointed'
    vbcf(:,:,:)=0.0d0
  END IF

  DO nn = 1, num_satellite_atms
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(2i3,8F8.3)') nn, tvsch(ic,nn), vbcf(:,ic,nn)
  END DO
  END DO
!
end SUBROUTINE obsope_atms_read
!------------------------------------------------------------------------------
SUBROUTINE interpolate_atms(nn)
  USE mod_adm
  USE mod_cnst
  USE mod_obsope_common
  USE mod_grd, ONLY : &
    GRD_zs,           &
    GRD_ZSFC

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: i, l, p, k, n
  INTEGER :: ierr, nv
  INTEGER :: tmp_id(nv3d)
  REAL(8) :: fac1, fac2, fac3
  REAL(8) :: fac_sum

  atms(nn)%obsdata_3d_tmp=0.0
  atms(nn)%obsdata_3d=0.0
  atms(nn)%obsdata_2d_tmp=0.0
  atms(nn)%obsdata_2d=0.0
  atms(nn)%lsql_model_tmp=0.0
  atms(nn)%lsql_model=0.0
  atms(nn)%elev_model_tmp=0.0
  atms(nn)%elev_model=0.0

  tmp_id(1)=id_temp_nicam
  tmp_id(2)=id_qvap_nicam
  tmp_id(3)=id_pres_nicam
  !WRITE(ADM_LOG_FID,*) 'interpolate_atms'
  DO i = 1, nobs_atms(nn)
    !WRITE(ADM_LOG_FID,'(6f12.5)') atms(nn)%lon(i), atms(nn)%lat(i),   &
    !                              atms(nn)%saza(i), atms(nn)%saaz(i), &
    !                              atms(nn)%soza(i), atms(nn)%soaz(i) 
    IF( atms(nn)%inprc(i) ) THEN
      l    = atms(nn)%l_index(i)
      fac1 = atms(nn)%w1(i)
      fac2 = atms(nn)%w2(i)
      fac3 = atms(nn)%w3(i)
      fac_sum = fac1 + fac2 + fac3
      write(ADM_LOG_FID,'(i5,4f12.5)') l, fac1, fac2, fac3, fac_sum
      write(ADM_LOG_FID,'(3i5)') atms(nn)%n1_index(i), atms(nn)%n2_index(i), atms(nn)%n3_index(i)
      DO nv = 1, nv3d
        DO k = 1, ADM_vlayer
          atms(nn)%obsdata_3d_tmp(k,i,nv) = &
             ( fac1 * icodata4_3d(atms(nn)%n1_index(i),k,l,tmp_id(nv)) &
             + fac2 * icodata4_3d(atms(nn)%n2_index(i),k,l,tmp_id(nv)) &
             + fac3 * icodata4_3d(atms(nn)%n3_index(i),k,l,tmp_id(nv)) &
             ) / fac_sum
        END DO ! k
        IF(tmp_id(nv) == id_pres_nicam) THEN
          atms(nn)%obsdata_3d_tmp(:,i,nv) = EXP(atms(nn)%obsdata_3d_tmp(:,i,nv))
        END IF
        IF(tmp_id(nv) == id_qvap_nicam) THEN
          atms(nn)%obsdata_3d_tmp(:,i,nv) = atms(nn)%obsdata_3d_tmp(:,i,nv) * q2ppmv
        END IF
      END DO
    END IF
  END DO

  DO i = 1, nobs_atms(nn)
    IF( atms(nn)%inprc(i) ) THEN
      l    = atms(nn)%l_index(i)
      fac1 = atms(nn)%w1(i)
      fac2 = atms(nn)%w2(i)
      fac3 = atms(nn)%w3(i)
      fac_sum = fac1 + fac2 + fac3
      DO nv = 1, nv2d
        atms(nn)%obsdata_2d_tmp(i,nv) = &
                ( fac1 * icodata4_2d(atms(nn)%n1_index(i),1,l,nv) &
                + fac2 * icodata4_2d(atms(nn)%n2_index(i),1,l,nv) &
                + fac3 * icodata4_2d(atms(nn)%n3_index(i),1,l,nv) &
                ) / fac_sum
        IF(nv == id_surp_nicam) THEN
          atms(nn)%obsdata_2d_tmp(i,nv) = EXP(atms(nn)%obsdata_2d_tmp(i,nv))
        END IF
        IF(nv == id_qv2m_nicam) THEN
          atms(nn)%obsdata_2d_tmp(i,nv) = atms(nn)%obsdata_2d_tmp(i,nv) * q2ppmv
        END IF
      END DO

      atms(nn)%lsql_model_tmp(i) = &
                ( fac1 * icoland(atms(nn)%n1_index(i),l) &
                + fac2 * icoland(atms(nn)%n2_index(i),l) &
                + fac3 * icoland(atms(nn)%n3_index(i),l) &
                ) / fac_sum

      atms(nn)%elev_model_tmp(i) = &
                ( fac1 * GRD_zs(atms(nn)%n1_index(i),ADM_KNONE,l,GRD_ZSFC) &
                + fac2 * GRD_zs(atms(nn)%n2_index(i),ADM_KNONE,l,GRD_ZSFC) &
                + fac3 * GRD_zs(atms(nn)%n3_index(i),ADM_KNONE,l,GRD_ZSFC) &
                ) / fac_sum

    END IF
  END DO

  CALL MPI_ALLREDUCE( atms(nn)%obsdata_3d_tmp, atms(nn)%obsdata_3d, &
                      ADM_vlayer*nobs_atms(nn)*nv3d,                 &
                      MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (atms(nn)%obsdata_3d)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

  CALL MPI_ALLREDUCE( atms(nn)%obsdata_2d_tmp, atms(nn)%obsdata_2d, &
                      nobs_atms(nn)*nv2d,                            &
                      MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (atms(nn)%obsdata_2d)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

  CALL MPI_ALLREDUCE( atms(nn)%lsql_model_tmp, atms(nn)%lsql_model, &
                      nobs_atms(nn),                                 &
                      MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  atms(nn)%lsql(:)=INT(atms(nn)%lsql_model(:))

  CALL MPI_ALLREDUCE( atms(nn)%elev_model_tmp, atms(nn)%elev_model, &
                      nobs_atms(nn),                                 &
                      MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)

  atms(nn)%elev(:)=atms(nn)%elev_model(:)

  DO i = 1, nobs_atms(nn)
    WRITE(ADM_LOG_FID,'(i7,3f12.5)') i, (maxval(atms(nn)%obsdata_3d(:,i,nv)),nv=1,3)
  END DO

  DO i = 1, nobs_atms(nn)
    WRITE(ADM_LOG_FID,'(i7,7f12.5)') i, (atms(nn)%obsdata_2d(i,nv),nv=1,7)
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  !DEALLOCATE( atms(nn)%obsdata_3d_tmp )
  !DEALLOCATE( atms(nn)%obsdata_2d_tmp )

  !DEALLOCATE( inprc )
  !DEALLOCATE( l_index )
  !DEALLOCATE( n1_index )
  !DEALLOCATE( n2_index )
  !DEALLOCATE( n3_index )
  !DEALLOCATE( w1 )
  !DEALLOCATE( w2 )
  !DEALLOCATE( w3 )

END SUBROUTINE interpolate_atms
!------------------------------------------------------------------------------
SUBROUTINE calc_radiance_atms
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

  tmp_id(1)=id_temp_nicam
  tmp_id(2)=id_qvap_nicam
  tmp_id(3)=id_pres_nicam

  WRITE(ADM_LOG_FID,*) 'calc_radiance'
  DO nn = 1, num_satellite_atms
  DO n = 1, nobs_atms(nn)
    WRITE(ADM_LOG_FID,'(6f12.5)') atms(nn)%lon(n), atms(nn)%lat(n),   &
                                  atms(nn)%saza(n), atms(nn)%saaz(n), &
                                  atms(nn)%soza(n), atms(nn)%soaz(n)
  END DO
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  DO nn = 1, num_satellite_atms
    atms(nn)%bt_tmp(:,:)=0.0
    atms(nn)%bt(:,:)=0.0
    atms(nn)%trans_tmp(:,:,:)=0.0
    atms(nn)%trans(:,:,:)=0.0
  END DO

  prc_per_inst=ADM_prc_all/num_satellite_atms
  nn=(ADM_prc_me-1)/prc_per_inst+1
  pe_in_satellite=ADM_prc_me - (nn-1)*prc_per_inst

  WRITE(ADM_LOG_FID,*) 'prc_per_inst =', prc_per_inst
  WRITE(ADM_LOG_FID,*) 'nn =', nn
  WRITE(ADM_LOG_FID,*) 'pe_in_satellite =', pe_in_satellite
  IF(flush_text) FLUSH(ADM_LOG_FID)

  IF( nobs_atms(nn) /= 0 ) THEN
    obsint=nobs_atms(nn)/prc_per_inst
    IF(nobs_atms(nn)<prc_per_inst) obsint=1
    sobs=obsint*(pe_in_satellite-1)+1
    eobs=obsint*pe_in_satellite
    IF(ADM_prc_me==prc_per_inst*nn) eobs=nobs_atms(nn)
    IF(nobs_atms(nn) >= prc_per_inst .OR. &
       nobs_atms(nn) <  prc_per_inst .AND. pe_in_satellite <= nobs_atms(nn) ) THEN

      WRITE(ADM_LOG_FID,*) 'nn ', nn
      WRITE(ADM_LOG_FID,*) 'pe_in_satellite ', pe_in_satellite
      WRITE(ADM_LOG_FID,*) 'obsint ', obsint
      WRITE(ADM_LOG_FID,*) 'sobs ', sobs
      WRITE(ADM_LOG_FID,*) 'eobs ', eobs
      DO n = sobs, eobs
        WRITE(ADM_LOG_FID,*) 'iobs=', n
        WRITE(ADM_LOG_FID,*) maxval(atms(nn)%obsdata_3d(:,n,1)), &
                             minval(atms(nn)%obsdata_3d(:,n,1))
        WRITE(ADM_LOG_FID,*) maxval(atms(nn)%obsdata_3d(:,n,2)), &
                             minval(atms(nn)%obsdata_3d(:,n,2))
        WRITE(ADM_LOG_FID,*) maxval(atms(nn)%obsdata_3d(:,n,3)), &
                             minval(atms(nn)%obsdata_3d(:,n,3))
        WRITE(ADM_LOG_FID,*) 'tsfc= ', atms(nn)%obsdata_2d(n,id_tsfc_nicam)
        WRITE(ADM_LOG_FID,*) 'qv2m= ', atms(nn)%obsdata_2d(n,id_qv2m_nicam)
        WRITE(ADM_LOG_FID,*) 'surp= ', atms(nn)%obsdata_2d(n,id_surp_nicam)
        WRITE(ADM_LOG_FID,*) 'u10m= ', atms(nn)%obsdata_2d(n,id_u10m_nicam)
        WRITE(ADM_LOG_FID,*) 'v10m= ', atms(nn)%obsdata_2d(n,id_v10m_nicam)
        WRITE(ADM_LOG_FID,*) 'soza= ', atms(nn)%soza(n)
        WRITE(ADM_LOG_FID,*) 'soaz= ', atms(nn)%soaz(n)
        WRITE(ADM_LOG_FID,*) 'saza= ', atms(nn)%saza(n)
        WRITE(ADM_LOG_FID,*) 'saaz= ', atms(nn)%saaz(n)
        WRITE(ADM_LOG_FID,*) 'elev= ', atms(nn)%elev(n)
        WRITE(ADM_LOG_FID,*) 'lon = ', atms(nn)%lon(n)
        WRITE(ADM_LOG_FID,*) 'lat = ', atms(nn)%lat(n)
        WRITE(ADM_LOG_FID,*) 'lsql= ', atms(nn)%lsql(n)
      END DO
      IF(flush_text) FLUSH(ADM_LOG_FID)
     
      CALL atms_fwd( ADM_vlayer, eobs-sobs+1, rttovcoef_fname(nn),              &
          DBLE(atms(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 3)),             &
          DBLE(atms(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 1)),             &
          DBLE(atms(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 2)),             &
          DBLE(atms(nn)%obsdata_2d(                 sobs:eobs, id_tsfc_nicam)), &
          DBLE(atms(nn)%obsdata_2d(                 sobs:eobs, id_qv2m_nicam)), &
          DBLE(atms(nn)%obsdata_2d(                 sobs:eobs, id_surp_nicam)), &
          DBLE(atms(nn)%obsdata_2d(                 sobs:eobs, id_u10m_nicam)), &
          DBLE(atms(nn)%obsdata_2d(                 sobs:eobs, id_v10m_nicam)), &
          DBLE(atms(nn)%soza ( sobs:eobs )), &
          DBLE(atms(nn)%soaz ( sobs:eobs )), &
          DBLE(atms(nn)%saza ( sobs:eobs )), &
          DBLE(atms(nn)%saaz ( sobs:eobs )), &
          DBLE(atms(nn)%elev ( sobs:eobs )), &
          DBLE(atms(nn)%lon  ( sobs:eobs )), &
          DBLE(atms(nn)%lat  ( sobs:eobs )), &
          DBLE(atms(nn)%lsql ( sobs:eobs )), &
          atms(nn)%bt_tmp   ( :,sobs:eobs ), &
          atms(nn)%trans_tmp( :,:,sobs:eobs ) )

      WRITE(ADM_LOG_FID,*) 'End RTTOV for ATMS'
      IF(flush_text) FLUSH(ADM_LOG_FID)

    END IF
  END IF

  DO nn = 1, num_satellite_atms
    CALL MPI_ALLREDUCE(atms(nn)%bt_tmp, atms(nn)%bt, &
                       maxtvsch*nobs_atms(nn),              &
                       MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(atms(nn)%trans_tmp, atms(nn)%trans, &
                       ADM_vlayer*maxtvsch*nobs_atms(nn),         &
                       MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  END DO 
  WRITE(ADM_LOG_FID,*) 'End ALLREDUCE of BT for ATMS'
  IF(flush_text) FLUSH(ADM_LOG_FID)

  DO nn = 1, num_satellite_atms
  DO n = 1, nobs_atms(nn)
  DO ic = 1, maxtvsch
    atms(nn)%ominusb(ic,n)=atms(nn)%odat(ic,n)-atms(nn)%bt(ic,n)
  END DO
  END DO
  END DO

END SUBROUTINE calc_radiance_atms
!------------------------------------------------------------------------------
SUBROUTINE calc_vbc_atms(nn)
  USE mod_obsope_common
  USE mod_adm
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n, ic, k
  INTEGER :: ifoot, ichan
  REAL(8) :: iwlr

  atms(nn)%vbc_pred(:,:,:)=0.0

  DO n = 1, nobs_atms(nn)
    atms(nn)%vbc_pred(1,:,n)=-999.9 ! IWLR(1000-200 hPa)
    atms(nn)%vbc_pred(2,:,n)=-999.9 ! IWLR( 200- 50 hPa)
    atms(nn)%vbc_pred(3,:,n)=-999.9
    atms(nn)%vbc_pred(4,:,n)=0.0d0
    atms(nn)%vbc_pred(5,:,n)=0.0d0
    atms(nn)%vbc_pred(6,:,n)=0.0d0
    atms(nn)%vbc_pred(7,:,n)=1.0d0/COS(atms(nn)%saza(n)*deg2rad)
    atms(nn)%vbc_pred(8,:,n)=0.0d0

    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(atms(nn)%obsdata_3d(k,n,3) >  200.0d0 .AND. &
         atms(nn)%obsdata_3d(k,n,3) < 1000.0d0 ) THEN
         iwlr = iwlr + &
           ( atms(nn)%obsdata_3d(k+1,n,1) - atms(nn)%obsdata_3d(k,n,1) ) * &
           ( atms(nn)%trans(ADM_vlayer-k,   tvsch(ic,nn), n) -              &
             atms(nn)%trans(ADM_vlayer-k+1, tvsch(ic,nn), n))
      END IF
      END DO
      atms(nn)%vbc_pred(1,ic,n)=iwlr
    END DO
    
    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(atms(nn)%obsdata_3d(k,n,3) >   50.0d0 .AND. &
         atms(nn)%obsdata_3d(k,n,3) <  200.0d0 ) THEN
         iwlr = iwlr + &
           ( atms(nn)%obsdata_3d(k+1,n,1) - atms(nn)%obsdata_3d(k,n,1) ) * &
           ( atms(nn)%trans(ADM_vlayer-k,   tvsch(ic,nn), n) - &
             atms(nn)%trans(ADM_vlayer-k+1, tvsch(ic,nn), n))
      END IF
      END DO
      atms(nn)%vbc_pred(2,ic,n)=iwlr
    END DO

    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(atms(nn)%obsdata_3d(k,n,3) >    5.0d0 .AND. &
         atms(nn)%obsdata_3d(k,n,3) <   50.0d0 ) THEN
         iwlr = iwlr + &
           ( atms(nn)%obsdata_3d(k+1,n,1) - atms(nn)%obsdata_3d(k,n,1) ) * &
           ( atms(nn)%trans(ADM_vlayer-k,   tvsch(ic,nn), n) - &
             atms(nn)%trans(ADM_vlayer-k+1, tvsch(ic,nn), n))
      END IF
      END DO
      atms(nn)%vbc_pred(3,ic,n)=iwlr
    END DO

    WRITE(ADM_LOG_FID,*) 'AIRMASS BIAS'
    DO ic = 1, ntvsch(nn)
      atms(nn)%airmass_bias(:,ic,n)=atms(nn)%vbc_pred(:,ic,n)*vbcf(:,ic,nn)
      WRITE(ADM_LOG_FID,'(2i5,8f10.5)') nn, ic, atms(nn)%airmass_bias(:,ic,n)
    END DO
    IF(flush_text) FLUSH(ADM_LOG_FID)

  END DO

  ! Correct Bias
  atms(nn)%odat_bc(:,:)=0.0
  DO n = 1, nobs_atms(nn)
    ifoot=atms(nn)%fov(n)
    DO ic = 1, ntvsch(nn)
      ichan=tvsch(ic,nn)
      atms(nn)%odat_bc(ichan,n) = atms(nn)%odat(ichan,n)             & ! y
                                - vbcf_scan(ifoot,ic,nn)             & ! scan bias
                                + sum(atms(nn)%airmass_bias(:,ic,n))   ! airmass bias
    END DO
  END DO

END SUBROUTINE calc_vbc_atms
!------------------------------------------------------------------------------
SUBROUTINE quality_control_atms(nn)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n, ic, ichan

  DO n = 1, nobs_atms(nn)

    IF( atms(nn)%lsql(n) == 0 ) THEN
      DO ic = 1, ntvsch(nn)
        ichan=tvsch(ic,nn)
        IF( ichan <= 8 ) THEN
          atms(nn)%qc(ichan,n) = 0
        END IF
      END DO
    END IF

    IF( atms(nn)%lsql(n) ==1 ) THEN
      IF( atms(nn)%lwp(n) > 0.12 ) THEN
        DO ic = 1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          IF( ichan <= 7 ) THEN
            atms(nn)%qc(ichan,n) = 0
          END IF
        END DO
      END IF
      IF( atms(nn)%lwp(n) > 0.15 ) THEN
        DO ic = 1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          IF( ichan == 8 ) THEN
            atms(nn)%qc(ichan,n) = 0
          END IF
        END DO
      END IF

    END IF

  END DO

END SUBROUTINE quality_control_atms
!------------------------------------------------------------------------------
SUBROUTINE update_vbc_atms(imem, nn)
  USE mod_adm
  USE mod_vbc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  INTEGER, INTENT(IN) :: nn
  INTEGER :: i, ic, ichan
  INTEGER :: ntvschan( ntvsch(nn) )

  IF( imem > 0 ) RETURN

  DO i = 1, nobs_atms(nn)
  DO ic = 1, ntvsch(nn)
    ichan=tvsch(ic,nn)
    IF( ABS( atms(nn)%odat_bc(ichan,i) - atms(nn)%bt(ichan,i) ) > atms(nn)%err(ichan,i)*3.0 ) THEN
      atms(nn)%qc(ichan,i)=0
    END IF    
  END DO
  END DO

  ALLOCATE( vbca(maxvbc,maxtvsch,num_satellite_atms) )
  WRITE(ADM_LOG_FID,*) 'vbcf', nn
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(8f10.6)') (vbcf(i,ic,nn),i=1,maxvbc)
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  CALL das_vbc( nobs_atms(nn), maxvbc, ntvsch(nn),             &
                tvsname(nn), tvsch(:,nn), ntvsch(nn),           &
                atms(nn)%odat_bc(tvsch(1:ntvsch(nn),nn),:),    &
                atms(nn)%bt(     tvsch(1:ntvsch(nn),nn),:),    &
                atms(nn)%vbc_pred(:,1:ntvsch(nn),:),           &
                vbcf(:,:,nn), vbca(:,:,nn),                     &
                atms(nn)%qc(     tvsch(1:ntvsch(nn),nn),:),    &
                atms(nn)%err(    tvsch(1:ntvsch(nn),nn),:),    &
                ntvschan )

  WRITE(ADM_LOG_FID,*) 'vbca', nn
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(8f10.6)') (vbca(i,ic,nn),i=1,maxvbc)
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  vbc_coef_out_fname=TRIM(vbc_coef_out_fname)//'_'//TRIM(tvsname(nn))//'_atms'
  CALL vbc_write(trim(vbc_coef_out_fname), vbca(:,:,nn), maxvbc, maxtvsch, &
                 tvsinst(:,nn), tvsch(:,nn), ntvsch(nn), ntvschan)

END SUBROUTINE update_vbc_atms
!------------------------------------------------------------------------------
SUBROUTINE update_scanbias_atms(nn)
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
          !WRITE(ADM_LOG_FID,'(A,i4,10F10.3)') tvsname_scan, ifoot, lwp, dum,
          !scanread(1:ntvsch(nn))
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

  scanbias_est_ofname=TRIM(scanbias_est_ofname)//'_'//TRIM(tvsname(nn))//'_atms'
  !WRITE(ADM_LOG_FID,*) TRIM(scanbias_est_ofname)
  !WRITE(ADM_LOG_FID,*) nfootp(nn)
  !WRITE(ADM_LOG_FID,*) maxtvsch
  !WRITE(ADM_LOG_FID,*) tvsinst(:,nn)
  !WRITE(ADM_LOG_FID,*) tvsch(:,nn)
  !WRITE(ADM_LOG_FID,*) ntvsch(nn)
  !IF(flush_text) FLUSH(ADM_LOG_FID)
  CALL vbc_scan_write(trim(scanbias_est_ofname), vbcf_scan(:,1:ntvsch(nn),nn), &
                 nfootp(nn), ntvsch(nn), &
                 tvsinst(:,nn), tvsch(1:ntvsch(nn),nn), ntvsch(nn))

END SUBROUTINE update_scanbias_atms
!-------------------------------------------------------------------------------
SUBROUTINE output_atms(imem, nn)
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

  WHERE( atms(nn)%lon(:) < 0.0 ) atms(nn)%lon(:) = atms(nn)%lon(:) + 360.0

  WRITE(cimem(1:6),'(I6.6)') imem
  fname=TRIM(output_dirname_atms)//'/'//TRIM(tvsname(nn))//&
        TRIM(output_basename_atms)//TRIM(cimem)//'.dat'
  WRITE(ADM_LOG_FID,*) TRIM(fname)
  IF(flush_text) FLUSH(ADM_LOG_FID)
  OPEN(1, FILE=TRIM(fname),FORM='unformatted',ACCESS='sequential')
  DO i = 1, nobs_atms(nn)
    WRITE(1) REAL(id_bt_obs),       REAL(atms(nn)%lsql(i)),    &
           & REAL(atms(nn)%lon(i)), REAL(atms(nn)%lat(i)),     &
           & REAL(atms(nn)%saza(i)),                           &
           & REAL(atms(nn)%obsdata_2d(i,id_tsfc_nicam)),       &
           & REAL(0.0), REAL(atms(nn)%fov(i)),                 &
           & REAL(atms(nn)%obsdata_3d(atms(nn)%weight_maxlev(  &
           &      tvsch(1:ntvsch(nn),nn),i), i, 3)),           &
           & REAL(atms(nn)%odat_bc(tvsch(1:ntvsch(nn),nn),i)), &
           & REAL(atms(nn)%err    (tvsch(1:ntvsch(nn),nn),i)), &
           & REAL(atms(nn)%bt     (tvsch(1:ntvsch(nn),nn),i)), &
           & REAL(atms(nn)%qc     (tvsch(1:ntvsch(nn),nn),i))

  END DO
  CLOSE(1)

  IF(output_text) THEN
    fname=TRIM(output_dirname_atms)//'/'//TRIM(tvsname(nn))//&
          TRIM(output_basename_atms)//TRIM(cimem)//'.txt'
    OPEN(2, FILE=TRIM(fname),FORM='formatted')
    DO i = 1, nobs_atms(nn)
      WRITE(2,'(200F8.2)')                                       &
             & REAL(atms(nn)%lwp(i)), REAL(atms(nn)%lsql(i)),    &
             & REAL(atms(nn)%lsql_model(i)),                     &
             & REAL(atms(nn)%lon(i)), REAL(atms(nn)%lat(i)),     &
             & REAL(atms(nn)%saza(i)),                           &
             & REAL(atms(nn)%obsdata_2d(i,id_tsfc_nicam)),       &
             & REAL(0.0), REAL(atms(nn)%fov(i)),                 &
             & REAL(atms(nn)%obsdata_3d(atms(nn)%weight_maxlev(  &
             &      tvsch(1:ntvsch(nn),nn),i), i, 3)),           &
             & REAL(atms(nn)%odat   (tvsch(1:ntvsch(nn),nn),i)), &
             & REAL(atms(nn)%odat_bc(tvsch(1:ntvsch(nn),nn),i)), &
             & REAL(atms(nn)%err    (tvsch(1:ntvsch(nn),nn),i)), &
             & REAL(atms(nn)%bt     (tvsch(1:ntvsch(nn),nn),i)), &
             & REAL(atms(nn)%qc     (tvsch(1:ntvsch(nn),nn),i))
    END DO
    CLOSE(2)
  END IF

  CALL tvs_ominusb_output(nn, imem, ntvsch(nn), nobs_atms(nn),        &
                          atms(nn)%ominusb(tvsch(1:ntvsch(nn),nn),:), &
                          atms(nn)%fov(:),                            &
                          atms(nn)%lsql(:),                           &
                          atms(nn)%lwp(:),                            &
                          TRIM(ominusb_fname_atms), tvsname(nn),      &
                          islot)

END SUBROUTINE output_atms
!------------------------------------------------------------------------------
SUBROUTINE set_instrument
  tvsch = 0
  !
  !JPSS-0 (Suomi-NPP)
  !
  tvsname(1) = 'JP00'
  tvsinst(1,1) = rttv_plat_jpss
  tvsinst(2,1) = 0
  tvsinst(3,1) = rttv_inst_atms
  tvsch( 1,1)  =  7
  tvsch( 2,1)  =  8
  tvsch( 3,1)  =  9
  ntvsch(1)=3
  nfootp(1)=32
  RETURN
END SUBROUTINE set_instrument
!-----------------------------------------------------------------------
END MODULE mod_obsope_atms
