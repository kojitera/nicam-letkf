PROGRAM main
  USE mod_adm, ONLY :         & 
       ADM_proc_init,         &
       ADM_proc_finish,       &
       ADM_setup,             &
       ADM_MULTI_PRC,         &
       ADM_gall,              &
       ADM_kall,              &
       ADM_lall,              &
       ADM_gall_pl,           &
       ADM_lall_pl,           &
       ADM_prc_tab,           &
       ADM_prc_me,            &
       ADM_MAXFNAME,          &
       ADM_LOG_FID
  USE mod_misc
  USE mod_fio, ONLY : &  ! [add] H.Yashiro 20110826
    FIO_setup,        &
    FIO_input,        &
    FIO_HMID,         &
    FIO_HSHORT,       &
    FIO_HLONG,        &
    FIO_INTEG_FILE,   &
    FIO_SPLIT_FILE,   &
    FIO_FREAD,        &
    FIO_ICOSAHEDRON,  &
    FIO_BIG_ENDIAN,   &
    headerinfo,       &
    datainfo
  USE mod_comm, ONLY :           &
       COMM_setup
  USE mod_latlon, ONLY :         &
       LATLON_setup,             &
       LATLON_READ_outdirname
  USE mod_cnst, ONLY :           &
       CNST_setup
  USE mod_grd, ONLY :            &
       GRD_setup
  USE mod_gtl, ONLY :            &
       GTL_generate_uv
  USE mod_gmtr, ONLY :           &
       GMTR_setup

  IMPLICIT NONE
  INTEGER                   :: nmem
  INTEGER                   :: imem
  INTEGER                   :: l, rgnid, n, i, g, k
  INTEGER                   :: nbv
  CHARACTER(ADM_MAXFNAME)   :: basename
  CHARACTER(ADM_MAXFNAME)   :: fname
  CHARACTER(ADM_MAXFNAME)   :: mean_basename
  CHARACTER(ADM_MAXFNAME)   :: mean_fname
  CHARACTER(ADM_MAXFNAME)   :: member_basename
  CHARACTER(ADM_MAXFNAME)   :: restart_basename
  CHARACTER(ADM_MAXFNAME)   :: member_fname
  CHARACTER(ADM_MAXFNAME)   :: output_basename
  CHARACTER(ADM_MAXFNAME)   :: restart_layername
  CHARACTER(3)              :: dir_prefix = ''
  CHARACTER(3)              :: ens_num
  CHARACTER(6)              :: cmem
  REAL(8), ALLOCATABLE      :: ddd1(:,:,:,:,:)
  REAL(8), ALLOCATABLE      :: ddd2(:,:,:,:,:)
  REAL(8), ALLOCATABLE      :: mean(:,:,:,:)
  REAL(8), ALLOCATABLE      :: sprd(:,:,:,:)
  REAL(8), ALLOCATABLE      :: sprd_ave(:,:)
  REAL(8), ALLOCATABLE      :: ddd3_pl(:,:,:,:,:)

  INTEGER,PARAMETER         :: nv3d=8 ! u,v,t,q
  CHARACTER(LEN=FIO_HSHORT) :: varname3d(nv3d)
  INTEGER, ALLOCATABLE      :: ifid(:)

  LOGICAL, SAVE  :: start_mem_zero = .true.        ! [add] Koji 2017.09.08  

  NAMELIST / RESTARTPARAM /  &
       mean_basename,        &
       member_basename,      &
       restart_basename,     &
       output_basename,      &
       restart_layername,    &
       start_mem_zero,       &
       dir_prefix,           &
       nbv

  CALL ADM_proc_init(ADM_MULTI_PRC)
  CALL ADM_setup('nhm_driver.cnf')
  CALL FIO_setup ! [add] H.Yashiro 20110826

  !--- < comm module setup > ---
  !------ setup commnication.
  CALL COMM_setup
  !
  !--- < cnst module setup > ---
  !------ setup the constants.
  CALL CNST_setup
  !
  !--- < grid module setup > ---
  !------ setup the icosahedral grid and vertical grid.
  CALL GRD_setup
  !
  !--- < gmetrics module setup > ---
  !------ setup the horizontal metrics.
  CALL GMTR_setup

  OPEN(1,FILE='nhm_driver.cnf')
  READ(1,NML=RESTARTPARAM) 
  CLOSE(1)
  
  WRITE(ADM_LOG_FID,*) 'mean_basename      = ', TRIM(mean_basename)
  WRITE(ADM_LOG_FID,*) 'member_basename    = ', TRIM(member_basename)
  WRITE(ADM_LOG_FID,*) 'restart_basename   = ', TRIM(restart_basename)
  WRITE(ADM_LOG_FID,*) 'restart_layername  = ', TRIM(restart_layername)
  WRITE(ADM_LOG_FID,*) 'start_mem_zero     = ', start_mem_zero
  WRITE(ADM_LOG_FID,*) 'nbv                = ', nbv
  FLUSH(ADM_LOG_FID)

  ALLOCATE( ddd1(ADM_gall, ADM_kall, ADM_lall, 3, nbv+1) )
  ALLOCATE( ddd2(ADM_gall, ADM_kall, ADM_lall, 2, nbv+1) )
  ALLOCATE( sprd(ADM_gall, ADM_kall, ADM_lall, 2) )
  ALLOCATE( mean(ADM_gall, ADM_kall, ADM_lall, 2) )
  ALLOCATE( sprd_ave(ADM_kall, 2) )
  ALLOCATE( ddd3_pl(ADM_gall_pl, ADM_kall, ADM_lall_pl, 5, nbv+1) )
  
  CALL FIO_input(ddd1(:,:,:,1,nbv+1),mean_basename,'vx', restart_layername,1,ADM_kall,1)
  CALL FIO_input(ddd1(:,:,:,2,nbv+1),mean_basename,'vy', restart_layername,1,ADM_kall,1)
  CALL FIO_input(ddd1(:,:,:,3,nbv+1),mean_basename,'vz', restart_layername,1,ADM_kall,1)

  DO i = 1, nbv
    IF(start_mem_zero) THEN
      imem = i-1
    ELSE
      imem = i
    END IF
    WRITE(cmem,'(I6.6)') imem
    member_fname=TRIM(member_basename)//'/'//TRIM(dir_prefix)//TRIM(cmem)//'/'//TRIM(restart_basename)
    CALL FIO_input(ddd1(:,:,:,1,i),member_fname,'vx',restart_layername,1,ADM_kall,1)
    CALL FIO_input(ddd1(:,:,:,2,i),member_fname,'vy',restart_layername,1,ADM_kall,1)
    CALL FIO_input(ddd1(:,:,:,3,i),member_fname,'vz',restart_layername,1,ADM_kall,1)
  END DO

  DO i = 1, nbv+1
    CALL GTL_generate_uv(                     &
       ddd2(:,:,:,1,i), ddd3_pl(:,:,:,1,i),   &
       ddd2(:,:,:,2,i), ddd3_pl(:,:,:,2,i),   &
       ddd1(:,:,:,1,i), ddd3_pl(:,:,:,3,i),   &
       ddd1(:,:,:,2,i), ddd3_pl(:,:,:,4,i),   &
       ddd1(:,:,:,3,i), ddd3_pl(:,:,:,5,i), icos=0)
  END DO

  WRITE(ADM_LOG_FID,*) nbv+1, MINVAL(ddd2(:,2:ADM_kall-1,:,1,nbv+1)), MAXVAL(ddd2(:,2:ADM_kall-1,:,1,nbv+1))
  WRITE(ADM_LOG_FID,*) nbv+1, MINVAL(ddd2(:,2:ADM_kall-1,:,2,nbv+1)), MAXVAL(ddd2(:,2:ADM_kall-1,:,2,nbv+1))
  DO i = 1, nbv
    ddd2(:,:,:,1,i) = ( ddd2(:,:,:,1,i) - ddd2(:,:,:,1,nbv+1) )**2
    ddd2(:,:,:,2,i) = ( ddd2(:,:,:,2,i) - ddd2(:,:,:,2,nbv+1) )**2
    !ddd2(:,:,:,1,i) = ddd2(:,:,:,1,i) - ddd2(:,:,:,1,nbv+1)
    !ddd2(:,:,:,2,i) = ddd2(:,:,:,2,i) - ddd2(:,:,:,2,nbv+1)
    !WRITE(ADM_LOG_FID,*) i, MINVAL(ddd2(:,2:ADM_kall-1,:,1,i)), MAXVAL(ddd2(:,2:ADM_kall-1,:,1,i))
    !WRITE(ADM_LOG_FID,*) i, MINVAL(ddd2(:,2:ADM_kall-1,:,2,i)), MAXVAL(ddd2(:,2:ADM_kall-1,:,2,i))
  END DO

  DO l = 1, ADM_lall
  DO k = 1, ADM_kall
  DO g = 1, ADM_gall
    sprd(g,k,l,1)=SQRT( SUM(ddd2(g,k,l,1,1:nbv)) / REAL(nbv-1) )
    sprd(g,k,l,2)=SQRT( SUM(ddd2(g,k,l,2,1:nbv)) / REAL(nbv-1) )
  END DO
  END DO
  END DO

  DO k = 1, ADM_kall
    sprd_ave(k,1)=SUM(sprd(:,k,:,1))/REAL(ADM_gall*ADM_lall)
    sprd_ave(k,2)=SUM(sprd(:,k,:,2))/REAL(ADM_gall*ADM_lall)
  END DO

  DO l = 1, ADM_lall
    rgnid=ADM_prc_tab(l,ADM_prc_me)
    WRITE(ADM_LOG_FID,*) l, rgnid
    CALL MISC_make_idstr(fname,TRIM(output_basename),'rgn',rgnid)
    WRITE(ADM_LOG_FID,*) TRIM(fname)
    OPEN(1,FILE=TRIM(fname),FORM='unformatted',ACCESS='direct',&
         RECL=4*ADM_gall*ADM_kall)
    DO n = 1,2
      WRITE(1,REC=n) SNGL(sprd(:,:,l,n))
    END DO
    CLOSE(1)
  END DO

  CALL ADM_proc_finish

END PROGRAM main
