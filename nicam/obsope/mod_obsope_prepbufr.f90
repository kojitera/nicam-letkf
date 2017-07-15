MODULE mod_obsope_prepbufr
  IMPLICIT NONE
  PUBLIC

  INTEGER :: nobs_prepbufr=0
  CHARACTER(128) :: input_fname_prepbufr = ''
  CHARACTER(128) :: output_basename_prepbufr = ''
  INTEGER, ALLOCATABLE :: prep_elem(:)
  REAL(4), ALLOCATABLE :: prep_lat(:)
  REAL(4), ALLOCATABLE :: prep_lon(:)
  REAL(4), ALLOCATABLE :: prep_lev(:)
  REAL(4), ALLOCATABLE :: prep_odat(:)
  REAL(4), ALLOCATABLE :: prep_err(:)
  INTEGER, ALLOCATABLE :: prep_typ(:)
  INTEGER, ALLOCATABLE :: prep_qc_tmp(:)
  INTEGER, ALLOCATABLE :: prep_qc(:)
  REAL(4), ALLOCATABLE :: prep_obsdata_tmp(:)
  REAL(4), ALLOCATABLE :: prep_obsdata(:)

  LOGICAL, ALLOCATABLE, SAVE :: prep_inprc(:)
  INTEGER, ALLOCATABLE, SAVE :: prep_l_index(:)
  INTEGER, ALLOCATABLE, SAVE :: prep_n1_index(:)
  INTEGER, ALLOCATABLE, SAVE :: prep_n2_index(:)
  INTEGER, ALLOCATABLE, SAVE :: prep_n3_index(:)
  REAL(8), ALLOCATABLE, SAVE :: prep_w1(:)
  REAL(8), ALLOCATABLE, SAVE :: prep_w2(:)
  REAL(8), ALLOCATABLE, SAVE :: prep_w3(:)

  REAL(8), ALLOCATABLE, SAVE :: prep_klev(:,:)
  REAL(8), ALLOCATABLE, SAVE :: prep_kfact(:)

  INTEGER, ALLOCATABLE   :: id_obs_prep(:)
  LOGICAL, PRIVATE, SAVE :: output_text = .false.

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE obsope_prepbufr_init
  IMPLICIT NONE
  NAMELIST / prepbufr_cnf /   &
    input_fname_prepbufr,     &
    output_basename_prepbufr, &
    output_text

  OPEN(111, FILE='obsope.cnf')
  READ(111, NML=prepbufr_cnf)
  CLOSE(111)

END SUBROUTINE obsope_prepbufr_init
!------------------------------------------------------------------------------
SUBROUTINE obsope_prepbufr_read
  USE mod_adm
  USE mod_obsope_common
  IMPLICIT NONE
  INTEGER :: i
  INTEGER :: ierr
  REAL(4), ALLOCATABLE :: wk(:,:)

  OPEN(112,FILE=trim(input_fname_prepbufr), FORM='unformatted', &
       ACCESS='sequential', IOSTAT=ierr)
  
  ! Get number of observations
  nobs_prepbufr=0
  DO
    READ(112,IOSTAT=ierr)
    IF(ierr /= 0) THEN
      WRITE(ADM_LOG_FID,*) 'nobs_prepbufr =', nobs_prepbufr 
      EXIT 
    END IF
    nobs_prepbufr=nobs_prepbufr+1
  END DO
  REWIND(112)

  ! Temporary array
  ALLOCATE( wk(7,nobs_prepbufr) )
  ! Observations
  ALLOCATE( prep_elem(nobs_prepbufr) )
  ALLOCATE( prep_lat( nobs_prepbufr) )
  ALLOCATE( prep_lon( nobs_prepbufr) )
  ALLOCATE( prep_lev( nobs_prepbufr) )
  ALLOCATE( prep_odat(nobs_prepbufr) )
  ALLOCATE( prep_err( nobs_prepbufr) )
  ALLOCATE( prep_typ( nobs_prepbufr) )
  ALLOCATE( prep_qc_tmp(  nobs_prepbufr) )
  ALLOCATE( prep_qc(  nobs_prepbufr) )
  ALLOCATE( prep_obsdata_tmp(nobs_prepbufr) )
  ALLOCATE( prep_obsdata(nobs_prepbufr) )
  ALLOCATE( id_obs_prep(nobs_prepbufr) )

  ALLOCATE( prep_inprc (nobs_prepbufr) )
  ALLOCATE( prep_l_index (nobs_prepbufr) )
  ALLOCATE( prep_n1_index(nobs_prepbufr) )
  ALLOCATE( prep_n2_index(nobs_prepbufr) )
  ALLOCATE( prep_n3_index(nobs_prepbufr) )
  ALLOCATE( prep_w1(nobs_prepbufr) )
  ALLOCATE( prep_w2(nobs_prepbufr) )
  ALLOCATE( prep_w3(nobs_prepbufr) )
  ALLOCATE( prep_klev(2,nobs_prepbufr) )
  ALLOCATE( prep_kfact(nobs_prepbufr) )

  prep_inprc(:)=.false.

  ! Read all observations
  DO i = 1, nobs_prepbufr
    READ(112) wk(:,i)
  END DO

  ! Store in array
  DO i = 1, nobs_prepbufr
    prep_elem(i) = INT(wk(1,i))
    IF(wk(2,i) <= 180.0) THEN
      prep_lon(i) = wk(2,i)
    ELSE
      prep_lon(i) = wk(2,i) - 360.0
    END IF
    prep_lat(i)  = wk(3,i)
    prep_lev(i)  = wk(4,i)
    prep_odat(i) = wk(5,i)
    prep_err(i)  = wk(6,i)
    prep_typ(i)  = INT(wk(7,i))
    IF( prep_lat(i) == -90.0 ) THEN
      prep_lat(i) = -89.99
    ELSE IF( prep_lat(i) == 90.0 ) THEN
      prep_lat(i) = 89.99
    END IF
    SELECT CASE(prep_elem(i))
      CASE(id_u_obs)
        id_obs_prep(i)=id_u_nicam
      CASE(id_v_obs)
        id_obs_prep(i)=id_v_nicam
      CASE(id_t_obs)
        id_obs_prep(i)=id_temp_nicam
      CASE(id_q_obs)
        id_obs_prep(i)=id_qvap_nicam
      CASE(id_ps_obs)
        id_obs_prep(i)=id_surp_nicam
      CASE(id_rain_obs)
        id_obs_prep(i)=id_rain_nicam
      CASE DEFAULT
        WRITE(ADM_LOG_FID,*) 'Unknown observation type'
        WRITE(ADM_LOG_FID,*) wk(:,i)
        FLUSH(ADM_LOG_FID)
    END SELECT
  END DO

  DEALLOCATE(wk)

END SUBROUTINE obsope_prepbufr_read
!------------------------------------------------------------------------------
SUBROUTINE interpolate_prepbufr
  USE mod_adm
  USE mod_cnst
  USE mod_obsope_common
  IMPLICIT NONE
  INTEGER :: i, l, p
  INTEGER :: ks, ke
  INTEGER :: ierr
  REAL(8) :: fac1, fac2, fac3
  REAL(8) :: fac4, fac5, fac6
  REAL(8) :: fac_sum

  prep_qc_tmp(:)=0
  prep_qc(:)=0
  prep_obsdata_tmp(:)=0.0
  prep_obsdata(:)=0.0

  DO i = 1, nobs_prepbufr
   IF(prep_inprc(i)) THEN
     l  = prep_l_index(i)
     ks = prep_klev(1,i)
     ke = prep_klev(2,i)
     IF(ks==-999) THEN
       prep_obsdata_tmp(i)=-999.9
       prep_qc_tmp(i)=0
     ELSE
       IF(prep_elem(i)==id_ps_obs) THEN
         fac1 = prep_w1(i)
         fac2 = prep_w2(i)
         fac3 = prep_w3(i)
         fac_sum = fac1 + fac2 + fac3
         prep_obsdata_tmp(i) = &
             ( fac1 * icodata4_2d(prep_n1_index(i),ks,l,id_surp_nicam) &
             + fac2 * icodata4_2d(prep_n2_index(i),ks,l,id_surp_nicam) &
             + fac3 * icodata4_2d(prep_n3_index(i),ks,l,id_surp_nicam) &
             ) / fac_sum
         prep_obsdata_tmp(i) = EXP(prep_obsdata_tmp(i)) * prep_kfact(i)
       ELSE
         fac1 = prep_w1(i) * prep_kfact(i)
         fac2 = prep_w2(i) * prep_kfact(i)
         fac3 = prep_w3(i) * prep_kfact(i)
         fac4 = prep_w1(i) * (1.d0-prep_kfact(i))
         fac5 = prep_w2(i) * (1.d0-prep_kfact(i))
         fac6 = prep_w3(i) * (1.d0-prep_kfact(i))
         fac_sum = fac1 + fac2 + fac3 + fac4 + fac5 + fac6
         !WRITE(ADM_LOG_FID,'(A,7i8)') 'TEST', prep_n1_index(i), prep_n2_index(i), prep_n3_index(i), &
         !                     &  ks, ke, l, id_obs_prep(i)
         !WRITE(ADM_LOG_FID,'(4F12.5)') prep_elem(i), prep_lev(i), prep_odat(i), prep_typ(i)
         !FLUSH(ADM_LOG_FID)
         prep_obsdata_tmp(i) = &
             ( fac1 * icodata4_3d(prep_n1_index(i),ks,l,id_obs_prep(i)) &
             + fac2 * icodata4_3d(prep_n2_index(i),ks,l,id_obs_prep(i)) &
             + fac3 * icodata4_3d(prep_n3_index(i),ks,l,id_obs_prep(i)) &
             + fac4 * icodata4_3d(prep_n1_index(i),ke,l,id_obs_prep(i)) &
             + fac5 * icodata4_3d(prep_n2_index(i),ke,l,id_obs_prep(i)) &
             + fac6 * icodata4_3d(prep_n3_index(i),ke,l,id_obs_prep(i)) &
             ) / fac_sum
       END IF
       prep_qc_tmp(i)=1
     END IF
   END IF
  END DO

  CALL MPI_ALLREDUCE(prep_obsdata_tmp, prep_obsdata, nobs_prepbufr, &
                     MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (prep_obsdata)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF
  

  CALL MPI_ALLREDUCE(prep_qc_tmp, prep_qc, nobs_prepbufr, &
                     MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF(ierr /= 0) THEN
    WRITE(ADM_LOG_FID,*) 'Error in MPI_ALLREDUCE (prep_qc)'
    WRITE(ADM_LOG_FID,*) 'Error code is ', ierr
    CALL ADM_proc_stop
  END IF

!  IF(ADM_prc_me==1) THEN
!    DO i = 1, nobs_prepbufr
!      WRITE(ADM_LOG_FID,'(i7,5f10.3,i4,f10.3,i4)') &
!            prep_elem(i),   prep_lon(i), prep_lat(i), &
!            prep_lev(i),   prep_odat(i), prep_err(i), &
!            prep_qc(i), prep_obsdata(i), prep_typ(i)
!    END DO
!  END IF

END SUBROUTINE interpolate_prepbufr
!------------------------------------------------------------------------------
SUBROUTINE output_prepbufr(imem)
  USE mod_adm
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  CHARACTER(6)        :: cimem
  CHARACTER(256)      :: fname
  INTEGER :: i

  IF( imem < 0 ) RETURN

  WHERE( prep_lon(:) < 0.0 ) prep_lon(:)=prep_lon(:)+360.0

  IF(ADM_prc_me == 1) THEN
    WRITE(cimem(1:6),'(I6.6)') imem + 1
    fname=TRIM(output_basename_prepbufr)//TRIM(cimem)//'.dat'
    OPEN(1, FILE=TRIM(fname),FORM='unformatted',ACCESS='sequential')
    DO i = 1, nobs_prepbufr
      WRITE(1) REAL(prep_elem(i)),     REAL(prep_lon(i)), &
             & REAL(prep_lat(i)),      REAL(prep_lev(i)), &
             & REAL(prep_odat(i)),     REAL(prep_err(i)), &
             & REAL(prep_obsdata(i)),  REAL(prep_qc(i)),  &
             & REAL(prep_typ(i)) 
    END DO
    CLOSE(1)

    IF(output_text) THEN
      fname=TRIM(output_basename_prepbufr)//TRIM(cimem)//'.txt'
      OPEN(2, FILE=TRIM(fname),FORM='formatted')
      DO i = 1, nobs_prepbufr
        WRITE(2,'(9F10.3)')                                  &
               & REAL(prep_elem(i)),     REAL(prep_lon(i)), &
               & REAL(prep_lat(i)),      REAL(prep_lev(i)), &
               & REAL(prep_odat(i)),     REAL(prep_err(i)), &
               & REAL(prep_obsdata(i)),  REAL(prep_qc(i)),  &
               & REAL(prep_typ(i))
      END DO
      CLOSE(2)
    END IF

  END IF

END SUBROUTINE output_prepbufr
!------------------------------------------------------------------------------
END MODULE mod_obsope_prepbufr
