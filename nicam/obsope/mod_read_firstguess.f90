MODULE mod_read_firstguess
  USE mod_fio
  IMPLICIT NONE
  PUBLIC

  INTEGER, SAVE :: pstr, pend

  INTEGER, PRIVATE, PARAMETER :: nv3d=6
  INTEGER, PRIVATE, PARAMETER :: nv2d=7
  INTEGER, PRIVATE, PARAMETER :: flim=1

  CHARACTER(3)              :: dir_prefix          = ''
  CHARACTER(LEN=FIO_HLONG)  :: infile(flim)        = ''

  CHARACTER(LEN=FIO_HSHORT) :: varname3d(nv3d)=&
          &(/'ms_pres', 'ms_u', 'ms_v', 'ms_tem', 'ms_qv', &
          &    'ms_qc' /)
  CHARACTER(LEN=FIO_HSHORT) :: varname2d(nv2d)=&
          &(/'sa_tem_sfc', 'ss_u10m', 'ss_v10m', 'ss_cldw', &
          &       'ss_ps',  'ss_q2m',  'ss_t2m'/)

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE read_firstguess_init
  USE mod_adm
  USE mod_obsope_common
  USE mod_mnginfo_light, ONLY : &
       MNG_mnginfo_input,   &
       MNG_mnginfo_noinput, &
       MNG_PALL,            &
       MNG_prc_rnum,        &
       MNG_prc_tab
  IMPLICIT NONE
  CHARACTER(ADM_MAXFNAME) :: rgnmngfname

  NAMELIST / firstguess_cnf /  &
       rgnmngfname,            &
       infile,                 &
       dir_prefix

  OPEN(90,file='obsope.cnf')
  READ(90,NML=firstguess_cnf)
  CLOSE(90)

  CALL MNG_mnginfo_input( ADM_rlevel, trim(rgnmngfname) )
  CALL fio_syscheck()

  ALLOCATE( icodata4_3d(ADM_gall,ADM_vlayer,MNG_prc_rnum(1),nv3d) )
  ALLOCATE( icodata4_2d(ADM_gall,         1,MNG_prc_rnum(1),nv2d) )

END SUBROUTINE read_firstguess_init
!------------------------------------------------------------------------------
SUBROUTINE read_firstguess(imem)
  USE mod_adm
  USE mod_obsope_common
  USE mod_fio, ONLY : &
       FIO_HSHORT,       &
       FIO_HMID,         &
       FIO_HLONG,        &
       FIO_REAL4,        &
       FIO_REAL8,        &
       FIO_BIG_ENDIAN,   &
       FIO_ICOSAHEDRON,  &
       FIO_IGA_LCP,      &
       FIO_IGA_MLCP,     &
       FIO_INTEG_FILE,   &
       FIO_SPLIT_FILE,   &
       FIO_FREAD,        &
       headerinfo,       &
       datainfo
  USE mod_mnginfo_light, ONLY : &
       MNG_mnginfo_input,   &
       MNG_mnginfo_noinput, &
       MNG_PALL,            &
       MNG_prc_rnum,        &
       MNG_prc_tab 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: imem
  INTEGER :: p, pp, nv
  REAL(4), ALLOCATABLE :: data4allrgn(:)
  CHARACTER(6) :: cimem
  CHARACTER(LEN=FIO_HLONG)  :: infile_header(flim) = ''
  CHARACTER(LEN=FIO_HLONG)  :: infname             = ''
  INTEGER, ALLOCATABLE :: ifid(:)
  INTEGER, ALLOCATABLE :: prc_tab_C(:)
  INTEGER :: did
  TYPE(datainfo)   dinfo

  pstr = (ADM_prc_me - 1) * ADM_lall + 1
  pend = ADM_prc_me * ADM_lall

  ALLOCATE( ifid(1) )

  WRITE(cimem(1:6),'(I6.6)') imem
  infile_header(1)=trim(infile(1))//'/'//trim(dir_prefix)//trim(cimem)//'/history'
  ! READ FIRST GUESS
  DO p = pstr, pend
    pp = p - pstr + 1
    CALL fio_mk_fname(infname,trim(infile_header(1)),'pe',p-1,6)
    ALLOCATE( prc_tab_C(MNG_prc_rnum(p)) )
    prc_tab_C(:) = MNG_prc_tab(:,p)-1
    IF ( pp == 1 ) THEN
      CALL fio_put_commoninfo( FIO_SPLIT_FILE,  &
                               FIO_BIG_ENDIAN,  &
                               FIO_ICOSAHEDRON, &
                               ADM_glevel,      &
                               ADM_rlevel,      &
                               MNG_prc_rnum(p), &
                               prc_tab_C        )
    END IF

    write(*,*) trim(infname)
    CALL fio_register_file(ifid(pp),trim(infname))
    CALL fio_fopen(ifid(pp),FIO_FREAD)
    CALL fio_read_allinfo( ifid(pp) )

    ALLOCATE( data4allrgn(ADM_gall*ADM_vlayer*MNG_prc_rnum(p)) )
    DO nv = 1, nv3d
      CALL fio_seek_datainfo(did,ifid(pp),varname3d(nv),istep)
      CALL fio_get_datainfo(ifid(pp),did,dinfo)
      CALL fio_read_data(ifid(pp),did,data4allrgn(:))
      icodata4_3d(:,:,:,nv) = reshape( data4allrgn(:), shape(icodata4_3d(:,:,:,nv)) )
      IF( TRIM(varname3d(nv))=='ms_pres' ) THEN
        icodata4_3d(:,:,:,nv)=LOG(icodata4_3d(:,:,:,nv)*0.01)
      END IF
    END DO

    DEALLOCATE( data4allrgn )
    ALLOCATE( data4allrgn(ADM_gall*MNG_prc_rnum(p)) )

      DO nv = 1, nv2d
      CALL fio_seek_datainfo(did,ifid(pp),varname2d(nv),istep)
      CALL fio_get_datainfo(ifid(pp),did,dinfo)
      CALL fio_read_data(ifid(pp),did,data4allrgn(:))
      icodata4_2d(:,:,:,nv) = reshape( data4allrgn(:), shape(icodata4_2d(:,:,:,nv)) )
      IF( TRIM(varname2d(nv))=='ss_ps' ) THEN
        icodata4_2d(:,:,:,nv)=LOG(icodata4_2d(:,:,:,nv)*0.01)
      END IF
    END DO

    DEALLOCATE( data4allrgn )
    DEALLOCATE( prc_tab_C )

    CALL fio_fclose(ifid(pp))
  END DO ! p

  DEALLOCATE( ifid )

  !icodata4_2d(:,1,:,2)=icodata4_3d(:,1,:,2)
  !icodata4_2d(:,1,:,3)=icodata4_3d(:,1,:,3)

  DO nv = 1, nv3d
    WRITE(ADM_LOG_FID,'(2A15,3(A,ES24.16))')        &
            '+',       TRIM(varname3d(nv) ),        &
            ' max=', MAXVAL(icodata4_3d(:,:,:,nv)), &
            ' min=', MINVAL(icodata4_3d(:,:,:,nv)), &
            ' sum=',    SUM(icodata4_3d(:,:,:,nv))
  END DO

  DO nv = 1, nv2d
    WRITE(ADM_LOG_FID,'(2A15,3(A,ES24.16))')        &
            '+    ',   TRIM(varname2d(nv) ),        &
            ' max=', MAXVAL(icodata4_2d(:,:,:,nv)), &
            ' min=', MINVAL(icodata4_2d(:,:,:,nv)), &
            ' sum=',    SUM(icodata4_2d(:,:,:,nv))
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

END SUBROUTINE read_firstguess
!------------------------------------------------------------------------------
END MODULE mod_read_firstguess
