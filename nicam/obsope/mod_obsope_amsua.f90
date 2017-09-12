MODULE mod_obsope_amsua
  USE mod_obsope_common, ONLY: flush_text
  USE common_tvs_nicam, ONLY: &
    rttv_plat_noaa,           &
    rttv_plat_metop,          &
    rttv_plat_metop2,         &
    rttv_inst_amsua,          &
    rttv_chan_amsua 

  IMPLICIT NONE
  PUBLIC

  INTEGER, PRIVATE, PARAMETER   :: ninstrument=5
  INTEGER, PRIVATE, PARAMETER   :: maxtvsch=15
  INTEGER, PRIVATE, PARAMETER   :: maxvbc=8
  INTEGER, PRIVATE, PARAMETER   :: maxfoot=30
  INTEGER, PRIVATE, SAVE        :: ntvs
  INTEGER, PRIVATE, SAVE        :: tvsinst(3,ninstrument)
  INTEGER, PRIVATE, SAVE        :: tvsch(maxtvsch,ninstrument)
  INTEGER, PRIVATE, SAVE        :: ntvsch(ninstrument)
  INTEGER, PRIVATE, SAVE        :: nfootp(ninstrument)
  CHARACTER(4), PRIVATE, SAVE   :: tvsname(ninstrument)
  CHARACTER(256), PRIVATE, SAVE :: rttovcoef_fname(ninstrument) = ''
  CHARACTER(256), PRIVATE       :: input_fname_amsua = ''
  CHARACTER(256), PRIVATE       :: output_dirname_amsua = ''
  CHARACTER(2), PRIVATE         :: output_basename_amsua = ''
  CHARACTER(256), PRIVATE       :: ominusb_fname_amsua   = ''
  CHARACTER(256), PRIVATE       :: scanbias_fname = ''
  CHARACTER(256), PRIVATE       :: vbc_coef_fname = ''
  CHARACTER(256), PRIVATE       :: vbc_coef_out_fname = ''
  CHARACTER(256), PRIVATE       :: scanbias_est_ifname = ''
  CHARACTER(256), PRIVATE       :: scanbias_est_ofname = ''
  LOGICAL, PRIVATE, SAVE        :: output_text = .false.

  INTEGER, PRIVATE, PARAMETER :: nv3d=3
  INTEGER, PRIVATE, PARAMETER :: nv2d=7

  INTEGER :: num_satellite_amsua

  TYPE amsua_base
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
  END TYPE amsua_base
    
  TYPE(amsua_base), allocatable :: amsua(:)

  INTEGER, ALLOCATABLE :: nobs_amsua(:)

  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbcf_scan(:,:,:)
  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbcf(:,:,:)
  REAL(8), ALLOCATABLE, PRIVATE, SAVE :: vbca(:,:,:)

  PRIVATE :: set_instrument

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE obsope_amsua_init
  USE mod_adm
  IMPLICIT NONE
  INTEGER :: ierr

  NAMELIST / amsua_cnf /   &
    input_fname_amsua,     &
    output_dirname_amsua,  &
    output_basename_amsua, &
    rttovcoef_fname,       &
    ominusb_fname_amsua,   &
    scanbias_fname,        &
    vbc_coef_fname,        &
    vbc_coef_out_fname,    &
    scanbias_est_ifname,   &
    scanbias_est_ofname,   &
    output_text

  CALL set_instrument

  OPEN(121, FILE='obsope.cnf')
  READ(121, NML=amsua_cnf)
  CLOSE(121)

END SUBROUTINE obsope_amsua_init
!------------------------------------------------------------------------------
SUBROUTINE obsope_amsua_read
  USE mod_adm
  USE mod_obsope_common
  USE mod_scanbias
  USE mod_vbc
  IMPLICIT NONE
  INTEGER :: ierr
  INTEGER :: n, nn
  INTEGER :: ic
  REAL(4) :: cos_tmp

  OPEN(122,file=trim(input_fname_amsua),  form='unformatted', &
       access='sequential', status='old', iostat=ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in opening the file', trim(input_fname_amsua)
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    IF(flush_text) FLUSH(ADM_LOG_FID)
    CALL ADM_proc_stop
  END IF
  READ(122) num_satellite_amsua
  ALLOCATE( amsua(num_satellite_amsua) )
  ALLOCATE( nobs_amsua(num_satellite_amsua) )

  READ(122) ntvs
  READ(122) nobs_amsua(:)
  !WRITE(ADM_LOG_FID,*) nobs_amsua

  DO nn = 1, num_satellite_amsua
    ALLOCATE( amsua(nn)%said           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%lon            ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%lat            ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%elev           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%odat           ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%odat_bc        ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%err            ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%saza           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%soza           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%soaz           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%saaz           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%fov            ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%lsql           ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%qc_tmp         ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%qc             ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%obsdata_3d_tmp ( ADM_vlayer,nobs_amsua(nn),nv3d) )
    ALLOCATE( amsua(nn)%obsdata_3d     ( ADM_vlayer,nobs_amsua(nn),nv3d) )
    ALLOCATE( amsua(nn)%obsdata_2d_tmp ( nobs_amsua(nn),nv2d) )
    ALLOCATE( amsua(nn)%obsdata_2d     ( nobs_amsua(nn),nv2d) )
    ALLOCATE( amsua(nn)%lwp            ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%bt_tmp         ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%bt             ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%ominusb        ( maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%trans_tmp      ( ADM_vlayer, maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%trans          ( ADM_vlayer, maxtvsch,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%wk             ( 34,nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%inprc          ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%l_index        ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%n1_index       ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%n2_index       ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%n3_index       ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%w1             ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%w2             ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%w3             ( nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%vbcf_scan      ( maxfoot, maxtvsch ))
    ALLOCATE( amsua(nn)%weight         ( ADM_vlayer, maxtvsch, nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%weight_maxlev  ( maxtvsch, nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%vbc_pred       ( maxvbc, maxtvsch, nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%airmass_bias   ( maxvbc, maxtvsch, nobs_amsua(nn) ))
    ALLOCATE( amsua(nn)%lsql_model_tmp ( nobs_amsua(nn)) )
    ALLOCATE( amsua(nn)%lsql_model     ( nobs_amsua(nn)) )
  END DO

  DO nn = 1, num_satellite_amsua
  DO n = 1, nobs_amsua(nn)
    READ(122) amsua(nn)%wk(:,n)
  END DO
  END DO

  DO nn = 1, num_satellite_amsua
  !WRITE(ADM_LOG_FID,*) 'Satellite ', nn
  DO n = 1, nobs_amsua(nn)
    amsua(nn)%lat(n)    = amsua(nn)%wk( 7,n)
    amsua(nn)%lon(n)    = amsua(nn)%wk( 8,n)
    amsua(nn)%said(n)   = amsua(nn)%wk( 9,n)
    amsua(nn)%fov(n)    = amsua(nn)%wk(11,n)
    amsua(nn)%lsql(n)   = amsua(nn)%wk(12,n)
    amsua(nn)%saza(n)   = amsua(nn)%wk(13,n)
    amsua(nn)%soza(n)   = amsua(nn)%wk(14,n)
    amsua(nn)%elev(n)   = amsua(nn)%wk(15,n)/1000.0
    amsua(nn)%soaz(n)   = amsua(nn)%wk(17,n)
    amsua(nn)%saaz(n)   = amsua(nn)%wk(18,n)
    amsua(nn)%odat(:,n) = amsua(nn)%wk(19:33,n)
  END DO
  END DO

  DO nn = 1, num_satellite_amsua
    amsua(nn)%err(1:8,:)=0.3
    amsua(nn)%err(9,:)=0.5
    amsua(nn)%err(10:15,:)=1.0
    amsua(nn)%qc(:,:)=1
  END DO
  ! Liquid Water Path (LWP: g/kg)
  DO nn = 1, num_satellite_amsua
  DO n = 1, nobs_amsua(nn)
    cos_tmp=COS(amsua(nn)%saza(n)*deg2rad)
    IF( amsua(nn)%odat(1,n) < 285.0 .AND. &
        amsua(nn)%odat(2,n) < 285.0 ) THEN
      amsua(nn)%lwp(n) = cos_tmp*(&
                         8.24 - ( 2.539 - 1.744*cos_tmp ) * cos_tmp + &
                         0.754 * log( 285.0 - amsua(nn)%odat(1,n)) - &
                         2.265 * log( 285.0 - amsua(nn)%odat(2,n)))
    ELSE
      amsua(nn)%lwp(n) = 999.0
    END IF
  END DO
  END DO

  ! Read scan bias
  ALLOCATE( vbcf_scan(maxfoot,maxtvsch,num_satellite_amsua) )
  IF( scanbias_fname /= '' ) THEN
    CAll vbc_scan_read(trim(scanbias_fname), vbcf_scan, maxfoot, maxtvsch, &
                       num_satellite_amsua, tvsinst, tvsch, ntvsch)
  ELSE
    WRITE(ADM_LOG_FID,*) 'Caution! scan bias file is not appointed'
    vbcf_scan(:,:,:)=0.0d0
  END IF

  DO nn = 1, num_satellite_amsua
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(2i3,30F8.3)') nn, tvsch(ic,nn), vbcf_scan(:,ic,nn)
  END DO
  END DO

  ! Read varBC coefficients
  ALLOCATE( vbcf(maxvbc,maxtvsch,num_satellite_amsua) )
  IF( vbc_coef_fname /= '' ) THEN
    CAll vbc_read(trim(vbc_coef_fname), vbcf, maxvbc, maxtvsch, &
                       num_satellite_amsua, tvsinst, tvsch, ntvsch)
  ELSE
    WRITE(ADM_LOG_FID,*) 'Caution! Airmass bias coefficient file is not appointed'
    vbcf(:,:,:)=0.0d0
  END IF

  DO nn = 1, num_satellite_amsua
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(2i3,8F8.3)') nn, tvsch(ic,nn), vbcf(:,ic,nn)
  END DO
  END DO
!
end SUBROUTINE obsope_amsua_read
!------------------------------------------------------------------------------
SUBROUTINE interpolate_amsua(nn)
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

  amsua(nn)%obsdata_3d_tmp=0.0
  amsua(nn)%obsdata_3d=0.0
  amsua(nn)%obsdata_2d_tmp=0.0
  amsua(nn)%obsdata_2d=0.0
  amsua(nn)%lsql_model_tmp=0.0
  amsua(nn)%lsql_model=0.0

  tmp_id(1)=id_temp_nicam
  tmp_id(2)=id_qvap_nicam
  tmp_id(3)=id_pres_nicam
  !WRITE(ADM_LOG_FID,*) 'interpolate_amsua'
  DO i = 1, nobs_amsua(nn)
    !WRITE(ADM_LOG_FID,'(6f12.5)') amsua(nn)%lon(i), amsua(nn)%lat(i),   &
    !                              amsua(nn)%saza(i), amsua(nn)%saaz(i), &
    !                              amsua(nn)%soza(i), amsua(nn)%soaz(i) 
    IF( amsua(nn)%inprc(i) ) THEN
      l    = amsua(nn)%l_index(i)
      fac1 = amsua(nn)%w1(i)
      fac2 = amsua(nn)%w2(i)
      fac3 = amsua(nn)%w3(i)
      fac_sum = fac1 + fac2 + fac3
      write(ADM_LOG_FID,'(i5,4f12.5)') l, fac1, fac2, fac3, fac_sum
      write(ADM_LOG_FID,'(3i5)') amsua(nn)%n1_index(i), amsua(nn)%n2_index(i), amsua(nn)%n3_index(i)
      DO nv = 1, nv3d
        DO k = 1, ADM_vlayer
          amsua(nn)%obsdata_3d_tmp(k,i,nv) = &
             ( fac1 * icodata4_3d(amsua(nn)%n1_index(i),k,l,tmp_id(nv)) &
             + fac2 * icodata4_3d(amsua(nn)%n2_index(i),k,l,tmp_id(nv)) &
             + fac3 * icodata4_3d(amsua(nn)%n3_index(i),k,l,tmp_id(nv)) &
           ) / fac_sum
        END DO ! k
        IF(tmp_id(nv) == id_pres_nicam) THEN
          amsua(nn)%obsdata_3d_tmp(:,i,nv) = EXP(amsua(nn)%obsdata_3d_tmp(:,i,nv))
        END IF
        IF(tmp_id(nv) == id_qvap_nicam) THEN
          amsua(nn)%obsdata_3d_tmp(:,i,nv) = amsua(nn)%obsdata_3d_tmp(:,i,nv) * q2ppmv
        END IF
      END DO
    END IF
  END DO

  DO i = 1, nobs_amsua(nn)
    IF( amsua(nn)%inprc(i) ) THEN
      l    = amsua(nn)%l_index(i)
      fac1 = amsua(nn)%w1(i)
      fac2 = amsua(nn)%w2(i)
      fac3 = amsua(nn)%w3(i)
      fac_sum = fac1 + fac2 + fac3
      DO nv = 1, nv2d
        amsua(nn)%obsdata_2d_tmp(i,nv) = &
                ( fac1 * icodata4_2d(amsua(nn)%n1_index(i),1,l,nv) &
                + fac2 * icodata4_2d(amsua(nn)%n2_index(i),1,l,nv) &
                + fac3 * icodata4_2d(amsua(nn)%n3_index(i),1,l,nv) &
              ) / fac_sum
        IF(nv == id_surp_nicam) THEN
          amsua(nn)%obsdata_2d_tmp(i,nv) = EXP(amsua(nn)%obsdata_2d_tmp(i,nv))
        END IF
        IF(nv == id_qv2m_nicam) THEN
          amsua(nn)%obsdata_2d_tmp(i,nv) = amsua(nn)%obsdata_2d_tmp(i,nv) * q2ppmv
        END IF
      END DO

      amsua(nn)%lsql_model_tmp(i) = &
                ( fac1 * icoland(amsua(nn)%n1_index(i),l) &
                + fac2 * icoland(amsua(nn)%n2_index(i),l) &
                + fac3 * icoland(amsua(nn)%n3_index(i),l) &
              ) / fac_sum

    END IF
  END DO

  CALL MPI_ALLREDUCE( amsua(nn)%obsdata_3d_tmp, amsua(nn)%obsdata_3d, &
                      ADM_vlayer*nobs_amsua(nn)*nv3d,                 &
                      MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (amsua(nn)%obsdata_3d)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

  CALL MPI_ALLREDUCE( amsua(nn)%obsdata_2d_tmp, amsua(nn)%obsdata_2d, &
                      nobs_amsua(nn)*nv2d,                            &
                      MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (amsua(nn)%obsdata_2d)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

  CALL MPI_ALLREDUCE( amsua(nn)%lsql_model_tmp, amsua(nn)%lsql_model, &
                      nobs_amsua(nn),                                 &
                      MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

END SUBROUTINE interpolate_amsua
!------------------------------------------------------------------------------
SUBROUTINE calc_radiance_amsua
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

  !WRITE(ADM_LOG_FID,*) 'calc_radiance'
  !DO nn = 1, num_satellite_amsua
  !DO n = 1, nobs_amsua(nn)
  !  WRITE(ADM_LOG_FID,'(6f12.5)') amsua(nn)%lon(n), amsua(nn)%lat(n),   &
  !                                amsua(nn)%saza(n), amsua(nn)%saaz(n), &
  !                                amsua(nn)%soza(n), amsua(nn)%soaz(n)
  !END DO
  !END DO

  DO nn = 1, num_satellite_amsua
    amsua(nn)%bt_tmp(:,:)=0.0
    amsua(nn)%bt(:,:)=0.0
    amsua(nn)%trans_tmp(:,:,:)=0.0
    amsua(nn)%trans(:,:,:)=0.0
  END DO

  prc_per_inst=ADM_prc_all/num_satellite_amsua
  nn=(ADM_prc_me-1)/prc_per_inst+1
  pe_in_satellite=ADM_prc_me - (nn-1)*prc_per_inst

  IF( nobs_amsua(nn) /= 0 ) THEN
    obsint=nobs_amsua(nn)/prc_per_inst
    IF(nobs_amsua(nn)<prc_per_inst) obsint=1
    sobs=obsint*(pe_in_satellite-1)+1
    eobs=obsint*pe_in_satellite
    IF(ADM_prc_me==prc_per_inst*nn) eobs=nobs_amsua(nn)
    IF(nobs_amsua(nn) >= prc_per_inst .OR. &
       nobs_amsua(nn) <  prc_per_inst .AND. pe_in_satellite <= nobs_amsua(nn) ) THEN

      !WRITE(ADM_LOG_FID,*) 'nn ', nn
      !WRITE(ADM_LOG_FID,*) 'pe_in_satellite ', pe_in_satellite
      !WRITE(ADM_LOG_FID,*) 'obsint ', obsint
      !WRITE(ADM_LOG_FID,*) 'sobs ', sobs
      !WRITE(ADM_LOG_FID,*) 'eobs ', eobs
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_3d(10,sobs,1)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_3d(10,sobs,2)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_3d(10,sobs,3)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_2d(sobs,1)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_2d(sobs,2)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_2d(sobs,3)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_2d(sobs,4)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%obsdata_2d(sobs,5)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%soza(sobs)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%soaz(sobs)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%saza(sobs)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%saaz(sobs)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%elev(sobs)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%lon(sobs)
      !WRITE(ADM_LOG_FID,*) amsua(nn)%lat(sobs)
      !IF(flush_text) FLUSH(ADM_LOG_FID)
     
      CALL amsua_fwd( ADM_vlayer, eobs-sobs+1, rttovcoef_fname(nn),            &
          DBLE(amsua(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 3)),     &
          DBLE(amsua(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 1)),     &
          DBLE(amsua(nn)%obsdata_3d(ADM_vlayer:1:-1, sobs:eobs, 2)),     &
          DBLE(amsua(nn)%obsdata_2d(                 sobs:eobs, id_tsfc_nicam)), &
          DBLE(amsua(nn)%obsdata_2d(                 sobs:eobs, id_qv2m_nicam)), &
          DBLE(amsua(nn)%obsdata_2d(                 sobs:eobs, id_surp_nicam)), &
          DBLE(amsua(nn)%obsdata_2d(                 sobs:eobs, id_u10m_nicam)), &
          DBLE(amsua(nn)%obsdata_2d(                 sobs:eobs, id_v10m_nicam)), &
          DBLE(amsua(nn)%soza ( sobs:eobs )), &
          DBLE(amsua(nn)%soaz ( sobs:eobs )), &
          DBLE(amsua(nn)%saza ( sobs:eobs )), &
          DBLE(amsua(nn)%saaz ( sobs:eobs )), &
          DBLE(amsua(nn)%elev ( sobs:eobs )), &
          DBLE(amsua(nn)%lon  ( sobs:eobs )), &
          DBLE(amsua(nn)%lat  ( sobs:eobs )), &
          DBLE(amsua(nn)%lsql ( sobs:eobs )), &
          amsua(nn)%bt_tmp   ( :,sobs:eobs ), &
          amsua(nn)%trans_tmp( :,:,sobs:eobs ) )

    END IF
  END IF

  DO nn = 1, num_satellite_amsua
    CALL MPI_ALLREDUCE(amsua(nn)%bt_tmp, amsua(nn)%bt, &
                       maxtvsch*nobs_amsua(nn),              &
                       MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(amsua(nn)%trans_tmp, amsua(nn)%trans, &
                       ADM_vlayer*maxtvsch*nobs_amsua(nn),         &
                       MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  END DO 

  DO nn = 1, num_satellite_amsua
  DO n = 1, nobs_amsua(nn)
  DO ic = 1, maxtvsch
    amsua(nn)%ominusb(ic,n)=amsua(nn)%odat(ic,n)-amsua(nn)%bt(ic,n)
    IF( ic >= 10 .AND. ABS(amsua(nn)%lat(n)) >= 30.0 ) THEN
      amsua(nn)%ominusb(ic,n)=-99.9
    END IF
  END DO
  END DO
  END DO

END SUBROUTINE calc_radiance_amsua
!------------------------------------------------------------------------------
SUBROUTINE calc_vbc_amsua(nn)
  USE mod_obsope_common
  USE mod_adm
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n, ic, k
  INTEGER :: ifoot, ichan
  REAL(8) :: iwlr

  amsua(nn)%vbc_pred(:,:,:)=0.0

  DO n = 1, nobs_amsua(nn)
    amsua(nn)%vbc_pred(1,:,n)=-999.9 ! IWLR(1000-200 hPa)
    amsua(nn)%vbc_pred(2,:,n)=-999.9 ! IWLR( 200- 50 hPa)
    amsua(nn)%vbc_pred(3,:,n)=0.0d0
    amsua(nn)%vbc_pred(4,:,n)=(amsua(nn)%obsdata_2d(n,id_tsfc_nicam)-273.15d0)/10.0d0
    amsua(nn)%vbc_pred(5,:,n)=0.0d0
    amsua(nn)%vbc_pred(6,:,n)=0.0d0
    amsua(nn)%vbc_pred(7,:,n)=1.0d0/COS(amsua(nn)%saza(n)*deg2rad)
    amsua(nn)%vbc_pred(8,:,n)=0.0d0

    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(amsua(nn)%obsdata_3d(k,n,3) >  200.0d0 .AND. &
         amsua(nn)%obsdata_3d(k,n,3) < 1000.0d0 ) THEN
         iwlr = iwlr + &
           ( amsua(nn)%obsdata_3d(k+1,n,1) - amsua(nn)%obsdata_3d(k,n,1) ) * &
           ( amsua(nn)%trans(ADM_vlayer-k,   tvsch(ic,nn), n) -              &
             amsua(nn)%trans(ADM_vlayer-k+1, tvsch(ic,nn), n))
      END IF
      END DO
      amsua(nn)%vbc_pred(1,ic,n)=iwlr
    END DO
    
    DO ic = 1, ntvsch(nn)
      iwlr=0.0d0
      DO k = 1, ADM_vlayer-1
      IF(amsua(nn)%obsdata_3d(k,n,3) >   50.0d0 .AND. &
         amsua(nn)%obsdata_3d(k,n,3) <  200.0d0 ) THEN
         iwlr = iwlr + &
           ( amsua(nn)%obsdata_3d(k+1,n,1) - amsua(nn)%obsdata_3d(k,n,1) ) * &
           ( amsua(nn)%trans(ADM_vlayer-k,   tvsch(ic,nn), n) - &
             amsua(nn)%trans(ADM_vlayer-k+1, tvsch(ic,nn), n))
      END IF
      END DO
      amsua(nn)%vbc_pred(2,ic,n)=iwlr
    END DO

    WRITE(ADM_LOG_FID,*) 'AIRMASS BIAS'
    DO ic = 1, ntvsch(nn)
      amsua(nn)%airmass_bias(:,ic,n)=amsua(nn)%vbc_pred(:,ic,n)*vbcf(:,ic,nn)
      WRITE(ADM_LOG_FID,'(2i5,8f10.5)') nn, ic, amsua(nn)%airmass_bias(:,ic,n)
    END DO

  END DO

  ! Correct Bias
  amsua(nn)%odat_bc(:,:)=0.0
  DO n = 1, nobs_amsua(nn)
    ifoot=amsua(nn)%fov(n)
    DO ic = 1, ntvsch(nn)
      ichan=tvsch(ic,nn)
      amsua(nn)%odat_bc(ichan,n) = amsua(nn)%odat(ichan,n)             & ! y
                                 - vbcf_scan(ifoot,ic,nn)              & ! scan bias
                                 + sum(amsua(nn)%airmass_bias(:,ic,n))   ! airmass bias
    END DO
  END DO

END SUBROUTINE calc_vbc_amsua
!------------------------------------------------------------------------------
SUBROUTINE quality_control_amsua(nn)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n, ic, ichan

  DO n = 1, nobs_amsua(nn)

    IF( amsua(nn)%lsql(n) == 0 ) THEN
      DO ic = 1, ntvsch(nn)
        ichan=tvsch(ic,nn)
        IF( ichan <= 7 ) THEN
          amsua(nn)%qc(ichan,n) = 0
        END IF
      END DO
    END IF

    IF( amsua(nn)%lsql(n) ==1 ) THEN
      IF( amsua(nn)%lwp(n) > 0.12 ) THEN
        DO ic = 1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          IF( ichan <= 6 ) THEN
            amsua(nn)%qc(ichan,n) = 0
          END IF
        END DO
      END IF
      IF( amsua(nn)%lwp(n) > 0.15 ) THEN
        DO ic = 1, ntvsch(nn)
          ichan=tvsch(ic,nn)
          IF( ichan == 7 ) THEN
            amsua(nn)%qc(ichan,n) = 0
          END IF
        END DO
      END IF

    END IF

  END DO

END SUBROUTINE quality_control_amsua
!------------------------------------------------------------------------------
SUBROUTINE update_vbc_amsua(imem, nn)
  USE mod_adm
  USE mod_vbc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  INTEGER, INTENT(IN) :: nn
  INTEGER :: i, ic, ichan
  INTEGER :: ntvschan( ntvsch(nn) )

  IF( imem > 0 ) RETURN

  DO i = 1, nobs_amsua(nn)
  DO ic = 1, ntvsch(nn)
    ichan=tvsch(ic,nn)
    IF( ichan <= 9 ) THEN
      IF( ABS( amsua(nn)%odat_bc(ichan,i) - amsua(nn)%bt(ichan,i) ) > amsua(nn)%err(ichan,i)*3.0 ) THEN
        amsua(nn)%qc(ichan,i)=0
      END IF    
    ELSE IF( ichan >= 10 .AND. ABS(amsua(nn)%lat(i)) <= 30.0 ) THEN
      IF( ABS( amsua(nn)%odat_bc(ichan,i) - amsua(nn)%bt(ichan,i) ) > amsua(nn)%err(ichan,i)*3.0 ) THEN
        amsua(nn)%qc(ichan,i)=0
      END IF    
    END IF
  END DO
  END DO

  ALLOCATE( vbca(maxvbc,maxtvsch,num_satellite_amsua) )
  WRITE(ADM_LOG_FID,*) 'vbcf', nn
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(8f10.6)') (vbcf(i,ic,nn),i=1,maxvbc)
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  CALL das_vbc( nobs_amsua(nn), maxvbc, ntvsch(nn),             &
                tvsname(nn), tvsch(:,nn), ntvsch(nn),           &
                amsua(nn)%odat_bc(tvsch(1:ntvsch(nn),nn),:),    &
                amsua(nn)%bt(     tvsch(1:ntvsch(nn),nn),:),    &
                amsua(nn)%vbc_pred(:,1:ntvsch(nn),:),           &
                vbcf(:,:,nn), vbca(:,:,nn),                     &
                amsua(nn)%qc(     tvsch(1:ntvsch(nn),nn),:),    &
                amsua(nn)%err(    tvsch(1:ntvsch(nn),nn),:),    &
                ntvschan )

  WRITE(ADM_LOG_FID,*) 'vbca', nn
  DO ic = 1, ntvsch(nn)
    WRITE(ADM_LOG_FID,'(8f10.6)') (vbca(i,ic,nn),i=1,maxvbc)
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  vbc_coef_out_fname=TRIM(vbc_coef_out_fname)//'_'//TRIM(tvsname(nn))//'_amsua'
  CALL vbc_write(trim(vbc_coef_out_fname), vbca(:,:,nn), maxvbc, maxtvsch, &
                 tvsinst(:,nn), tvsch(:,nn), ntvsch(nn), ntvschan)

END SUBROUTINE update_vbc_amsua
!------------------------------------------------------------------------------
SUBROUTINE update_scanbias_amsua(nn)
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

  scanbias_est_ofname=TRIM(scanbias_est_ofname)//'_'//TRIM(tvsname(nn))//'_amsua'
  CALL vbc_scan_write(trim(scanbias_est_ofname), vbcf_scan(:,1:ntvsch(nn),nn), &
                 nfootp(nn), ntvsch(nn), &
                 tvsinst(:,nn), tvsch(1:ntvsch(nn),nn), ntvsch(nn))


END SUBROUTINE update_scanbias_amsua
!-------------------------------------------------------------------------------
SUBROUTINE output_amsua(imem, nn)
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

  WHERE( amsua(nn)%lon(:) < 0.0 ) amsua(nn)%lon(:) = amsua(nn)%lon(:) + 360.0

  WRITE(cimem(1:6),'(I6.6)') imem
  fname=TRIM(output_dirname_amsua)//'/'//TRIM(tvsname(nn))//&
        TRIM(output_basename_amsua)//TRIM(cimem)//'.dat'
  WRITE(ADM_LOG_FID,*) TRIM(fname)
  IF(flush_text) FLUSH(ADM_LOG_FID)
  OPEN(1, FILE=TRIM(fname),FORM='unformatted',ACCESS='sequential')
  DO i = 1, nobs_amsua(nn)
    WRITE(1) REAL(id_bt_obs),        REAL(amsua(nn)%lsql(i)),         &
           & REAL(amsua(nn)%lon(i)), REAL(amsua(nn)%lat(i)),          &
           & REAL(amsua(nn)%saza(i)),                                 &
           & REAL(amsua(nn)%obsdata_2d(i,id_tsfc_nicam)),             &
           & REAL(0.0), REAL(amsua(nn)%fov(i)),                       &
           & REAL(amsua(nn)%obsdata_3d(amsua(nn)%weight_maxlev(&
                  tvsch(1:ntvsch(nn),nn),i), i, 3)), &
           & REAL(amsua(nn)%odat_bc(tvsch(1:ntvsch(nn),nn),i)),       &
           & REAL(amsua(nn)%err    (tvsch(1:ntvsch(nn),nn),i)),       &
           & REAL(amsua(nn)%bt     (tvsch(1:ntvsch(nn),nn),i)),       &
           & REAL(amsua(nn)%qc     (tvsch(1:ntvsch(nn),nn),i))

  END DO
  CLOSE(1)

  IF(output_text) THEN
    fname=TRIM(output_dirname_amsua)//'/'//TRIM(tvsname(nn))//&
          TRIM(output_basename_amsua)//TRIM(cimem)//'.txt'
    OPEN(2, FILE=TRIM(fname),FORM='formatted')
    DO i = 1, nobs_amsua(nn)
      WRITE(2,'(40F8.2)')                                               &
             & REAL(amsua(nn)%lwp(i)), REAL(amsua(nn)%lsql(i)),         &
             & REAL(amsua(nn)%lsql_model(i)),                           &
             & REAL(amsua(nn)%lon(i)), REAL(amsua(nn)%lat(i)),          &
             & REAL(amsua(nn)%saza(i)),                                 &
             & REAL(amsua(nn)%obsdata_2d(i,id_tsfc_nicam)),             &
             & REAL(0.0), REAL(amsua(nn)%fov(i)),                       &
             & REAL(amsua(nn)%obsdata_3d(amsua(nn)%weight_maxlev(       &
                    tvsch(1:ntvsch(nn),nn),i), i, 3)), &
             & REAL(amsua(nn)%odat   (tvsch(1:ntvsch(nn),nn),i)),       &
             & REAL(amsua(nn)%odat_bc(tvsch(1:ntvsch(nn),nn),i)),       &
             & REAL(amsua(nn)%err    (tvsch(1:ntvsch(nn),nn),i)),       &
             & REAL(amsua(nn)%bt     (tvsch(1:ntvsch(nn),nn),i)),       &
             & REAL(amsua(nn)%qc     (tvsch(1:ntvsch(nn),nn),i))
    END DO
    CLOSE(2)
  END IF

  CALL tvs_ominusb_output(nn, imem, ntvsch(nn), nobs_amsua(nn),        &
                          amsua(nn)%ominusb(tvsch(1:ntvsch(nn),nn),:), &
                          amsua(nn)%fov(:),                            &
                          amsua(nn)%lsql(:),                           &
                          amsua(nn)%lwp(:),                            &
                          TRIM(ominusb_fname_amsua), tvsname(nn),      &
                          islot)

END SUBROUTINE output_amsua
!------------------------------------------------------------------------------
SUBROUTINE set_instrument
  tvsch = 0
  !
  ! NOAA-15 AMSU-A
  !
  tvsname(1) = 'AA15'
  tvsinst(1,1) = rttv_plat_noaa
  tvsinst(2,1) = 15
  tvsinst(3,1) = rttv_inst_amsua
  tvsch( 1,1)  =  7
  tvsch( 2,1)  =  8
  ntvsch(1)=2
  nfootp(1)=30
  !
  ! NOAA-16 AMSU-A
  !
  tvsname(2) = 'AA16'
  tvsinst(1,2) = rttv_plat_noaa
  tvsinst(2,2) = 16
  tvsinst(3,2) = rttv_inst_amsua
  tvsch( 1,2)  =  6
  ntvsch(2)=1
  nfootp(2)=30
  !
  ! NOAA-18 AMSU-A
  !
  tvsname(3) = 'AA18'
  tvsinst(1,3) = rttv_plat_noaa
  tvsinst(2,3) = 18
  tvsinst(3,3) = rttv_inst_amsua
  tvsch( 1,3)  =  6
  tvsch( 2,3)  =  7
  tvsch( 3,3)  =  8
  ntvsch(3)=3
  nfootp(3)=30
  !
  ! NOAA-19 AMSU-A
  !
  tvsname(4) = 'AA19'
  tvsinst(1,4) = rttv_plat_noaa
  tvsinst(2,4) = 19
  tvsinst(3,4) = rttv_inst_amsua
  tvsch( 1,4)  =  6
  tvsch( 2,4)  =  7
  ntvsch(4)=2
  nfootp(4)=30
  !
  ! METOP-2 AMSU-A
  !
  tvsname(5) = 'MA02'
  tvsinst(1,5) = rttv_plat_metop
  tvsinst(2,5) = 2
  tvsinst(3,5) = rttv_inst_amsua
  tvsch(1,5)  = 6
  tvsch(2,5)  = 8
  ntvsch(5)=2
  nfootp(5)=30

  RETURN
END SUBROUTINE set_instrument
!-----------------------------------------------------------------------
END MODULE mod_obsope_amsua
