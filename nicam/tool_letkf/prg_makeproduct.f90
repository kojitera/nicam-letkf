PROGRAM main
  USE mpi
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
       ADM_LOG_FID,           &
       ADM_prc_tab
  USE mod_misc
  USE mod_fio, ONLY :         &  ! [add] H.Yashiro 20110826
    FIO_setup,                &
    FIO_input,                &
    FIO_output,               &
    FIO_HMID,                 &
    FIO_HSHORT,               &
    FIO_HLONG,                &
    FIO_INTEG_FILE,           &
    FIO_SPLIT_FILE,           &
    FIO_FREAD,                &
    FIO_ICOSAHEDRON,          &
    FIO_BIG_ENDIAN,           &
    headerinfo,               &
    datainfo,                 &
    FIO_REAL8
  USE mod_comm, ONLY :        &
       COMM_setup
  USE mod_latlon, ONLY :      &
       LATLON_setup,          &
       LATLON_READ_outdirname
  USE mod_cnst, ONLY :        &
       CNST_setup,            &
       CNST_UNDEF
  USE mod_grd, ONLY :         &
       GRD_setup,             &
       GRD_gz,                &
       GRD_gzh,               &
       GRD_zs,                &
       GRD_zs_pl,             &
       GRD_vz,                &
       GRD_vz_pl
  USE mod_oprt, ONLY :        &
       OPRT_setup
  USE mod_vmtr, ONLY :        &
       VMTR_setup
  USE mod_gtl, ONLY :         &
       GTL_generate_uv
  USE mod_gmtr, ONLY :        &
       GMTR_setup

  IMPLICIT NONE
  INTEGER                   :: nmem
  INTEGER                   :: imem
  INTEGER                   :: g, l, rgnid, n, nv, ip
  INTEGER                   :: ierr
  INTEGER                   :: i, j
  INTEGER                   :: ix, iy
  CHARACTER(ADM_MAXFNAME)   :: basename
  CHARACTER(ADM_MAXFNAME)   :: fname
  CHARACTER(ADM_MAXFNAME)   :: mean_restart_basename
  CHARACTER(ADM_MAXFNAME)   :: mean_history_basename
  CHARACTER(ADM_MAXFNAME)   :: sprd_restart_basename
  CHARACTER(ADM_MAXFNAME)   :: sprd_history_basename
  CHARACTER(ADM_MAXFNAME)   :: sprd_uv_basename
  
  CHARACTER(ADM_MAXFNAME)   :: output_basename
  CHARACTER(ADM_MAXFNAME)   :: restart_layername
  CHARACTER(ADM_MAXFNAME)   :: surface_layername
  CHARACTER(ADM_MAXFNAME)   :: pressure_layername
  
  CHARACTER(3)              :: ens_num
  REAL(8), ALLOCATABLE      :: ddd_logp(:,:,:,:)
  REAL(8), ALLOCATABLE      :: ddd_3d_org(:,:,:,:)
  REAL(8), ALLOCATABLE      :: ddd_3d(:,:,:,:,:)
  REAL(8), ALLOCATABLE      :: ddd_2d(:,:,:,:,:)
  REAL(8), ALLOCATABLE      :: ddd_uv(:,:,:,:)
  REAL(8), ALLOCATABLE      :: ddd3_pl(:,:,:,:)
  REAL(4), ALLOCATABLE      :: sprd_uv(:,:,:,:)
  
  !---   
  REAL(8), ALLOCATABLE      :: ddd_3d_plev(:,:,:,:,:)
  
  REAL(8)                   :: plevel(200)
  REAL(8)                   :: logp(200)
     
  REAL(8), ALLOCATABLE      :: cc(:)
  REAL(8), ALLOCATABLE      :: ww(:)
  REAL(8), ALLOCATABLE      :: bx(:)
  REAL(8), ALLOCATABLE      :: by(:)
  
  INTEGER                   :: np
  INTEGER,PARAMETER         :: nv3d=8
  !  1: pres -> z, 2: temperature, 3: u-wind, 4: v-wind
  !  5: w-wind,    6: qv,          7: qc
  INTEGER,PARAMETER         :: nv2d=1
  CHARACTER(LEN=FIO_HSHORT) :: varname3d(nv3d)
  CHARACTER(LEN=FIO_HSHORT) :: varname2d(nv2d)
  INTEGER, ALLOCATABLE      :: IFid(:)

  !--- ico2ll
  INTEGER :: imax, jmax
  INTEGER,ALLOCATABLE       :: lon_index(:),lat_index(:)
  INTEGER,ALLOCATABLE       :: n1_index(:), n2_index(:), n3_index(:)
  REAL(8),ALLOCATABLE       :: w1(:), w2(:), w3(:)
      
  INTEGER,ALLOCATABLE       :: max_num_latlon(:)
  INTEGER,ALLOCATABLE       :: nstart(:), nend(:)
  CHARACTER(1024)           :: fin_llmap_head = '../llmap/llmap'
  CHARACTER(1024)           :: fin_llmap
      
  REAL(8),ALLOCATABLE       :: lon(:)  ! -pi - pi
  REAL(8),ALLOCATABLE       :: lat(:)  ! -pi/2 - pi/2
  REAL(8),ALLOCATABLE       :: lon_swap(:)
      
  REAL(8),ALLOCATABLE       :: var_out_3d(:,:,:,:,:)
  REAL(8),ALLOCATABLE       :: var_out_2d(:,:,:,:,:)
  REAL(8),ALLOCATABLE       :: var_3d(:,:,:,:,:)
  REAL(8),ALLOCATABLE       :: var_2d(:,:,:,:,:)
  REAL(8),ALLOCATABLE       :: var_out_swap(:,:)
  CHARACTER(5)              :: c5

  REAL(8),ALLOCATABLE       :: var_3d_rh(:,:,:)
           
  LOGICAL                   :: onefile=.false.
  CHARACTER(256)            :: output_dir
  CHARACTER(10)             :: cdate

  CHARACTER(LEN=FIO_HMID)   :: desc = 'INITIAL/RESTART DATA of PROGNOSTICVARIABLES'

  NAMELIST / RESTARTPARAM /    &
       mean_restart_basename,  &
       mean_history_basename,  &
       sprd_restart_basename,  &
       sprd_history_basename,  &
       sprd_uv_basename,       &
       output_basename,        &
       output_basename,        &
       varname3d,              &
       varname2d,              &
       plevel,                 &
       np,                     &
       fin_llmap_head,         &
       output_dir,             &
       cdate, onefile,         &
       restart_layername,      &
       pressure_layername

  CALL ADM_proc_init(ADM_MULTI_PRC)
  CALL ADM_setup('nhm_driver.cnf')
  CALL FIO_setup ! [add] H.Yashiro 20110826
  CALL CNST_setup
  CALL COMM_setup
  CALL FIO_setup
  CALL GRD_setup
  CALL GMTR_setup
  CALL OPRT_setup
  CALL VMTR_setup

  ALLOCATE( ddd_3d_org(ADM_gall,    ADM_kall,    ADM_lall,    3) ) ! mean vx, vy, vz
  ALLOCATE( ddd_3d    (ADM_gall,    ADM_kall,    ADM_lall, nv3d, 2) )
  ALLOCATE( ddd_uv    (ADM_gall,    ADM_kall,    ADM_lall,    2) )
  ALLOCATE( ddd_2d    (ADM_gall,           1,    ADM_lall, nv2d, 2) )
  ALLOCATE( ddd3_pl   (ADM_gall_pl, ADM_kall, ADM_lall_pl,    5) )
  ALLOCATE( ddd_logp  (ADM_gall,    ADM_kall,    ADM_lall,    2) )
  ALLOCATE(    sprd_uv(ADM_gall,    ADM_kall,    ADM_lall,    2) )

  OPEN(1,file='nhm_driver.cnf')
  READ(1,nml=RESTARTPARAM) 
  CLOSE(1)
  WRITE(ADM_LOG_FID,RESTARTPARAM)
  FLUSH(ADM_LOG_FID)
  WRITE(ADM_LOG_FID,*) TRIM(varname2d(1)), TRIM(varname2d(2))
  WRITE(ADM_LOG_FID,*) 'mean_restart_basename,', TRIM(mean_restart_basename)
  WRITE(ADM_LOG_FID,*) 'mean_history_basename,', TRIM(mean_history_basename)
  WRITE(ADM_LOG_FID,*) 'sprd_restart_basename,', TRIM(sprd_restart_basename)
  WRITE(ADM_LOG_FID,*) 'sprd_history_basename,', TRIM(sprd_history_basename)
  WRITE(ADM_LOG_FID,*) 'sprd_uv_basename     ,', TRIM(sprd_uv_basename)
  FLUSH(ADM_LOG_FID)
  logp(1:np)=log(plevel(1:np))
  ALLOCATE( ddd_3d_plev(ADM_gall, np, ADM_lall, nv3d, 2) )
  surface_layername='ZSSFC1'
  !
  !--- < Reading data > ---

  CALL FIO_input(    ddd_3d(:,:,:,1,1),mean_restart_basename,'pre',restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,2,1),mean_restart_basename,'tem',restart_layername,1,ADM_kall,1)
  CALL FIO_input(  ddd_3d_org(:,:,:,1),mean_restart_basename,'vx', restart_layername,1,ADM_kall,1)
  CALL FIO_input(  ddd_3d_org(:,:,:,2),mean_restart_basename,'vy', restart_layername,1,ADM_kall,1)
  CALL FIO_input(  ddd_3d_org(:,:,:,3),mean_restart_basename,'vz', restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,5,1),mean_restart_basename,'w',  restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,6,1),mean_restart_basename,'qv', restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,7,1),mean_restart_basename,'qc', restart_layername,1,ADM_kall,1)

  DO nv = 1, nv2d
    CALL FIO_input(ddd_2d(:,:,:,nv,1),mean_history_basename,TRIM(varname2d(nv)), surface_layername,1,1,1)
    WRITE(ADM_LOG_FID,'(A15,2F15.5)') TRIM(varname2d(nv)), &
          MINVAL(ddd_2d(:,1,:,nv,1)), MAXVAL(ddd_2d(:,1,:,nv,1))
  END DO

  ! READ SPREAD
  CALL FIO_input(    ddd_3d(:,:,:,1,2),sprd_restart_basename,'pre',restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,2,2),sprd_restart_basename,'tem',restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,5,2),sprd_restart_basename,'w',  restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,6,2),sprd_restart_basename,'qv', restart_layername,1,ADM_kall,1)
  CALL FIO_input(    ddd_3d(:,:,:,7,2),sprd_restart_basename,'qc', restart_layername,1,ADM_kall,1)

  DO nv = 1, nv2d
    CALL FIO_input(ddd_2d(:,:,:,nv,2),sprd_history_basename,TRIM(varname2d(nv)), surface_layername,1,1,1)
    WRITE(ADM_LOG_FID,'(A15,2F15.5)') TRIM(varname2d(nv)), &
          MINVAL(ddd_2d(:,1,:,nv,2)), MAXVAL(ddd_2d(:,1,:,nv,2))
  END DO

  DO l = 1, ADM_lall
    rgnid=ADM_prc_tab(l,ADM_prc_me)
    WRITE(ADM_LOG_FID,*) l, rgnid
    CALL MISC_make_idstr(fname,Trim(sprd_uv_basename),'rgn',rgnid)
    WRITE(ADM_LOG_FID,*) TRIM(fname)
    OPEN(1,file=TRIM(fname),FORM='unformatted',ACCESS='direct',&
         RECL=4*ADM_gall*ADM_kall)
    DO n = 1,2
      READ(1,REC=n) sprd_uv(:,:,l,n)
    END DO
    CLOSE(1)
  END DO
  
  ddd_3d(:,:,:,3,2)=DBLE(sprd_uv(:,:,:,1))
  ddd_3d(:,:,:,4,2)=DBLE(sprd_uv(:,:,:,2))

  CALL GTL_generate_uv(                   &
       ddd_3d(:,:,:,3,1), ddd3_pl(:,:,:,1), &
       ddd_3d(:,:,:,4,1), ddd3_pl(:,:,:,2), &
       ddd_3d_org(:,:,:,1), ddd3_pl(:,:,:,3), &
       ddd_3d_org(:,:,:,2), ddd3_pl(:,:,:,4), &
       ddd_3d_org(:,:,:,3), ddd3_pl(:,:,:,5), icos=0)


  ddd_logp(:,:,:,:)=LOG(ddd_3d(:,:,:,1,:)*0.01)
  ddd_3d(:,:,:,1,:)=ddd_3d(:,:,:,1,:)*0.01d0
  DO n = 1, ADM_kall
    WRITE(ADM_LOG_FID,*) n, ddd_logp(1,n,1,1:2)
  END DO

  ALLOCATE( cc(ADM_kall-2) )
  ALLOCATE( ww(ADM_kall-2) )
  ALLOCATE( bx(np) )
  ALLOCATE( by(np) )

  DO i = 1, 2
    DO l = 1, ADM_lall
    DO g = 1, ADM_gall
      IF( ddd_3d(g,2,l,1,1) >= 1000.0d0 ) THEN
        DO nv = 2, nv3d-1
          CALL spline( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                       ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, ww)
          CALL splint( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                       ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, np, logp(np:1:-1), by(np:1:-1))
          ddd_3d_plev(g,:,l,nv,i)=by(:)
        END DO

        CALL spline( ADM_kall-2, ddd_logp(g, ADM_kall-1:2:-1, l, 1), &
                     GRD_vz(g, ADM_kall-1:2:-1, l, 1), cc, ww)
        CALL splint( ADM_kall-2, ddd_logp(g, ADM_kall-1:2:-1, l, 1), &
                     GRD_vz(g, ADM_kall-1:2:-1, l, 1), cc, np, logp(np:1:-1), by(np:1:-1))
        ddd_3d_plev(g,:,l,1,i)=by(:)
      ELSE
        n = 1
        DO 
          IF(ddd_3d(g,2,l,1,1) >= plevel(n) ) THEN
            EXIT
          END IF
          n = n + 1
        END DO 
        DO nv = 2, nv3d-1
          CALL spline( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                       ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, ww)
          CALL splint( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                       ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, np-n+1, logp(np:n:-1), by(np:n:-1))
          ddd_3d_plev(g,n:np,l,nv,i)=by(n:np)
          ddd_3d_plev(g,1:n-1,l,nv,i)=CNST_UNDEF
        END DO
        CALL spline( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                     GRD_vz(g,ADM_kall-1:2:-1,l,1), cc, ww)
        CALL splint( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                     GRD_vz(g,ADM_kall-1:2:-1,l,1), cc, np-n+1, logp(np:n:-1), by(np:n:-1))
        ddd_3d_plev(g,n:np,l,1,i)=by(n:np)
        ddd_3d_plev(g,1:n-1,l,1,i)=CNST_UNDEF
      END IF
    END DO
    END DO
  END DO

  fin_llmap =  TRIM(fin_llmap_head) // '.info'
  OPEN(10,file=fin_llmap, &
       FORM='unformatted',STATUS='old' ,iostat=ierr)
  IF(ierr/=0) THEN
     WRITE(*,*) 'Cannot OPEN llmap info file!'
     stop
  ENDIF
  READ(10) imax
  ALLOCATE(lon(imax))
  ALLOCATE(lon_swap(imax))
  READ(10) lon(:)
  READ(10) jmax
  ALLOCATE(lat(jmax))
  READ(10) lat(:)
  CLOSE(10)
  WRITE(ADM_LOG_FID,*) 'imax =', imax
  WRITE(ADM_LOG_FID,*) 'jmax =', jmax

  ALLOCATE(max_num_latlon(ADM_lall))
  ALLOCATE(nstart(ADM_lall))
  ALLOCATE(nend(ADM_lall))

  ALLOCATE(lon_index(imax*jmax))
  ALLOCATE(lat_index(imax*jmax))
  ALLOCATE(n1_index(imax*jmax))
  ALLOCATE(n2_index(imax*jmax))
  ALLOCATE(n3_index(imax*jmax))
  ALLOCATE(w1(imax*jmax))
  ALLOCATE(w2(imax*jmax))
  ALLOCATE(w3(imax*jmax))

  ALLOCATE(var_out_3d(imax,jmax,np,nv3d,2))
  ALLOCATE(var_out_2d(imax,jmax,1,nv2d,2))
  ALLOCATE(var_3d(imax,jmax,np,nv3d,2))
  ALLOCATE(var_2d(imax,jmax,1,nv2d,2))
  ALLOCATE(var_out_swap(imax,jmax))
  ALLOCATE(var_3d_rh(imax,jmax,np))

  DO l=1, ADM_lall
     OPEN(20,file=fin_llmap,FORM='unformatted',STATUS='old')
     READ(20) max_num_latlon(l)
     CLOSE(20)
  ENDDO

  DO l=1, ADM_lall
    rgnid = ADM_prc_tab(l,ADM_prc_me)
    WRITE(c5,'(I5)') rgnid-1
    DO n=1, 5
      IF( c5(n:n) == ' ' ) c5(n:n) = '0'
    ENDDO

    fin_llmap = TRIM(fin_llmap_head) // '.rgn' // c5
    OPEN(20,FILE=fin_llmap,FORM='unformatted',STATUS='old')
    READ(20) max_num_latlon(l)
    nend(l) = sum(max_num_latlon(1:l))
    nstart(l) = nend(l) - max_num_latlon(l)+1
    IF(max_num_latlon(l)/=0) THEN
      READ(20) lon_index(nstart(l):nend(l))  ! lon(lon_index(i)) :
      READ(20) lat_index(nstart(l):nend(l))
      READ(20) n1_index(nstart(l):nend(l))
      READ(20) n2_index(nstart(l):nend(l))
      READ(20) n3_index(nstart(l):nend(l))
      READ(20) w1(nstart(l):nend(l))
      READ(20) w2(nstart(l):nend(l))
      READ(20) w3(nstart(l):nend(l))
    ENDIF
    CLOSE(20)

    var_out_3d(:,:,:,:,:)=0.0d0
    var_out_2d(:,:,:,:,:)=0.0d0
    DO i = 1, 2
      DO nv = 1, nv3d-1
      DO ip = 1, np
      DO n=nstart(l),nend(l)
        IF(ddd_3d_plev(n1_index(n),ip,l,nv,i) == CNST_UNDEF ) THEN
           var_out_3d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
        ELSE IF(ddd_3d_plev(n2_index(n),ip,l,nv,i) == CNST_UNDEF ) THEN
           var_out_3d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
        ELSE IF(ddd_3d_plev(n3_index(n),ip,l,nv,i) == CNST_UNDEF) THEN
           var_out_3d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
        ELSE
           var_out_3d(lon_index(n),lat_index(n),ip,nv,i) &
                = w1(n)*ddd_3d_plev(n1_index(n),ip,l,nv,i)      &
                + w2(n)*ddd_3d_plev(n2_index(n),ip,l,nv,i)      &
                + w3(n)*ddd_3d_plev(n3_index(n),ip,l,nv,i) 
        ENDIF
      END DO
      END DO
      END DO

      ip=1
      DO nv = 1, nv2d
      DO n=nstart(l),nend(l)
        IF(     ddd_2d(n1_index(n),ip,l,nv,i) == CNST_UNDEF ) THEN
           var_out_2d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
        ELSE IF(ddd_2d(n2_index(n),ip,l,nv,i) == CNST_UNDEF ) THEN
           var_out_2d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
        ELSE IF(ddd_2d(n3_index(n),ip,l,nv,i) == CNST_UNDEF) THEN
           var_out_2d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
        ELSE
           var_out_2d(lon_index(n),lat_index(n),ip,nv,i) &
                = w1(n)*ddd_2d(n1_index(n),ip,l,nv,i)      &
                + w2(n)*ddd_2d(n2_index(n),ip,l,nv,i)      &
                + w3(n)*ddd_2d(n3_index(n),ip,l,nv,i)
        ENDIF
      END DO
      END DO
    END DO

  END DO

  CALL MPI_ALLREDUCE(var_out_3d, var_3d, imax*jmax*np*nv3d*2, MPI_REAL8, MPI_SUM,&
                     MPI_COMM_WORLD, ierr)

  CALL MPI_ALLREDUCE(var_out_2d, var_2d, imax*jmax*nv2d*2, MPI_REAL8, MPI_SUM,&
                     MPI_COMM_WORLD, ierr)

  !--- Compute Relative humidity
  DO ip = 1, np
  DO iy = 1, jmax
  DO ix = 1, imax
    IF(var_3d(ix,iy,ip,2,1)==CNST_UNDEF .AND. var_3d(ix,iy,ip,6,1)==CNST_UNDEF) THEN
      var_3d_rh(ix,iy,ip)=CNST_UNDEF
    ELSE
      var_3d_rh(ix,iy,ip)=calc_rh(var_3d(ix,iy,ip,2,1), plevel(ip)*100.0, var_3d(ix,iy,ip,6,1) )
      WRITE(ADM_LOG_FID,'(4F12.5)') var_3d(ix,iy,ip,2,1), plevel(ip)*100.0, &
                                    var_3d(ix,iy,ip,6,1), var_3d_rh(ix,iy,ip)
    END IF
  END DO
  END DO
  END DO

  !--- OUTPUT
  IF(ADM_prc_me==1) THEN
    IF(onefile) THEN
      fname=TRIM(output_dir)//'/'//TRIM(cdate)//'.dat'
      OPEN(1,file=TRIM(fname),FORM='unformatted',ACCESS='direct',RECL=4*imax*jmax)
      n=1
      i=1
      DO nv = 1, nv3d-1
        DO ip = 1, np
          WRITE(1,REC=n) SNGL(var_3d(:,:,ip,nv,i))
          WRITE(ADM_LOG_FID,*) ip, MINVAL(var_3d(:,:,ip,nv,i)), MAXVAL(var_3d(:,:,ip,nv,i))
          n=n+1
        END DO
        WRITE(ADM_LOG_FID,*) TRIM(varname3d(nv)), MINVAL(var_3d(:,:,:,nv,i)), MAXVAL(var_3d(:,:,:,nv,i))
      END DO

      DO ip = 1, np
        WRITE(1,REC=n) SNGL(var_3d_rh(:,:,ip))
        WRITE(ADM_LOG_FID,*) ip, MINVAL(var_3d_rh(:,:,ip)), MAXVAL(var_3d_rh(:,:,ip))
        n=n+1
      END DO

      DO nv = 1, nv2d
        WRITE(1,REC=n) SNGL(var_2d(:,:,1,nv,i))
        WRITE(ADM_LOG_FID,*) TRIM(varname2d(nv)), MINVAL(var_2d(:,:,:,nv,i)), MAXVAL(var_2d(:,:,:,nv,i))
        n=n+1
      END DO

      i=2
      DO nv = 2, nv3d-1
        DO ip = 1, np
          WRITE(1,REC=n) SNGL(var_3d(:,:,ip,nv,i))
          WRITE(ADM_LOG_FID,*) ip, MINVAL(var_3d(:,:,ip,nv,i)), MAXVAL(var_3d(:,:,ip,nv,i))
          n=n+1
        END DO
        WRITE(ADM_LOG_FID,*) TRIM(varname3d(nv)), MINVAL(var_3d(:,:,:,nv,i)), MAXVAL(var_3d(:,:,:,nv,i))
      END DO

      DO ip = 1, np
        WRITE(1,REC=n) SNGL(var_3d_rh(:,:,ip))
        WRITE(ADM_LOG_FID,*) ip, MINVAL(var_3d_rh(:,:,ip)), MAXVAL(var_3d_rh(:,:,ip))
        n=n+1
      END DO

      DO nv = 1, nv2d
        WRITE(1,REC=n) SNGL(var_2d(:,:,1,nv,i))
        WRITE(ADM_LOG_FID,*) TRIM(varname2d(nv)), MINVAL(var_2d(:,:,:,nv,i)), MAXVAL(var_2d(:,:,:,nv,i))
        n=n+1
      END DO
      
      CLOSE(1)

    ELSE
      DO nv = 1, nv3d-1
        fname=TRIM(output_dir)//'/'//TRIM(varname3d(nv))//'_'//TRIM(cdate)//'.dat'
        OPEN(1,file=TRIM(fname),FORM='unformatted',ACCESS='direct',RECL=4*imax*jmax*np)
        WRITE(1,REC=1) SNGL(var_3d(:,:,:,nv,i))
        CLOSE(1)
      END DO
  
      DO nv = 1, nv2d
        fname=TRIM(output_dir)//'/'//TRIM(varname2d(nv))//'_'//TRIM(cdate)//'.dat'
        OPEN(1,file=TRIM(fname),FORM='unformatted',ACCESS='direct',RECL=4*imax*jmax)
        WRITE(1,REC=1) SNGL(var_2d(:,:,1,nv,i))
        CLOSE(1)
      END DO
    END IF
    
  END IF

  CALL ADM_proc_finish

CONTAINS

!!$  --- This function calculates relative humidity [ % ].
  function calc_rh( tk, prs, qv )
  use mod_cnst, only :        &
       CP_T0 => CNST_TMELT,   &
       CNST_EPSV
!!$ === Declarations === $!!
    implicit none

    ! [ IN ]
    real(8), intent( in ) :: &
         tk, &                          ! Temperature: [ K ]
         prs, &                         ! Pressure: [ Pa ]
         qv                             ! specific humidity [ kg/kg ]
    ! [ RESULT ]
    real(8) :: calc_rh

    ! [ WORK ]
    real(8) :: psat   ! Saturation vapor pressure: [ Pa ]
    real(8) :: pvap   ! actual vapor pressure: [ Pa ]

!!$ --- End of Header --- $!!

    if ( tk >= CP_T0 ) then
       psat = moist_psat_water( tk )
    else
       psat = moist_psat_ice( tk )
    endif
    pvap = qv * prs / (CNST_EPSV + qv * (1.0d0 - CNST_EPSV))
    calc_rh = pvap / psat * 100.0d0

    return
  end function calc_rh
!!$ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ $!!

function moist_psat_water ( tem ) result(psat)
  ! psat : Clasius-Clapeyron: based on CPV, CPL constant
  use mod_cnst, only :        &
       CNST_RVAP,             &
       CNST_CPV,              &
       CNST_CL,               &
       CNST_LH0,              &
       CNST_LH00,             &
       CNST_PSAT0,            &
       CNST_TEM00
  implicit none

  real(8), intent(in)  :: tem
  real(8)              :: psat

  real(8) :: TEM_MIN = 10.d0
  !---------------------------------------------------------------------------

  psat = CNST_PSAT0 &
       * ( max(tem, TEM_MIN) / CNST_TEM00 ) ** ( ( CNST_CPV - CNST_CL ) / CNST_RVAP )     &
       * exp ( CNST_LH00 / CNST_RVAP * ( 1.0d0 / CNST_TEM00 - 1.0d0 / max(tem, TEM_MIN) ) )

end function moist_psat_water

!-----------------------------------------------------------------------------
function moist_psat_ice ( tem ) result(psat)
  ! psat : Clasius-Clapeyron: based on CPV, CPL constant
  ! for ice
  use mod_cnst, only :        &
       CNST_RVAP,             &
       CNST_CPV,              &
       CNST_CI,               &
       CNST_LHS00,            &
       CNST_LHS0,             &
       CNST_PSAT0,            &
       CNST_TEM00
  implicit none

  real(8), intent(in)  :: tem
  real(8)              :: psat

  real(8) :: TEM_MIN = 10.d0
  !---------------------------------------------------------------------------

  psat = CNST_PSAT0 &
       * ( max(tem, TEM_MIN) / CNST_TEM00 ) ** ( ( CNST_CPV - CNST_CI ) / CNST_RVAP )      &
       * exp ( CNST_LHS00 / CNST_RVAP * ( 1.0d0 / CNST_TEM00 - 1.0d0 / max(tem, TEM_MIN) ) )

end function moist_psat_ice


    SUBROUTINE SPLINE(N,X,Y,CC,WW)
    IMPLICIT REAL*8(A-H),REAL*8(O-Z)
    INTEGER N, I, K
    DIMENSION X(N),Y(N),CC(N),WW(N)

    CC(1)=0.
    WW(1)=0.

    DO 10 I=2,N-1
      SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
      P=SIG*CC(I-1)+2.
      CC(I)=(SIG-1.)/P
      WW(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
          /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*WW(I-1))/P
 10 CONTINUE
    QN=0.
    UN=0.
    CC(N)=(UN-QN*WW(N-1))/(QN*CC(N-1)+1.)
    DO 20 K=N-1,1,-1
      CC(K)=CC(K)*CC(K+1)+WW(K)
 20 CONTINUE
    END SUBROUTINE SPLINE

    SUBROUTINE SPLINT(N,AX,AY,CC,M,BX,BY)
    IMPLICIT REAL*8(A-H),REAL*8(O-Z)
    INTEGER N, M, K, J, KLO, KHI
    REAL*8 AX(N),AY(N),CC(N),BX(M),BY(M)

    DO 10 J=1,M
    X=BX(J)
    KLO=1
    KHI=N
  1 IF (KHI-KLO.GT.1) THEN
      K=(KHI+KLO)/2
      IF(AX(K).GT.X)THEN
        KHI=K
      ELSE
        KLO=K
      ENDIF
    GOTO 1
    ENDIF
    H=AX(KHI)-AX(KLO)
    IF (H.EQ.0.) PAUSE 'Bad AX input.'
    A=(AX(KHI)-X)/H
    B=(X-AX(KLO))/H
    Y=A*AY(KLO)+B*AY(KHI)+ &
          ((A**3-A)*CC(KLO)+(B**3-B)*CC(KHI))*(H**2)/6.
    BY(J)=Y
 10 CONTINUE
    END SUBROUTINE SPLINT
    
END PROGRAM main
