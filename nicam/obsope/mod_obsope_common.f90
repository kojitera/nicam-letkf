MODULE mod_obsope_common
  USE mpi
  USE mod_cnst
  USE mod_adm
  IMPLICIT NONE
  PUBLIC

  INTEGER, SAVE :: nobs_total

  INTEGER, SAVE :: istep
  INTEGER, SAVE :: islot
  INTEGER, SAVE :: nslot=7
  INTEGER, SAVE :: smem, emem
  INTEGER, SAVE :: nbv

  REAL(8), PARAMETER :: q2ppmv=1.60771704d6
  REAL(8), PARAMETER :: DEG2RAD = CNST_PI/180.0d0
  REAL(8), PARAMETER :: RAD2DEG = 180.0d0/CNST_PI

  ! ID for observations
  INTEGER, PARAMETER :: id_u_obs=2819
  INTEGER, PARAMETER :: id_v_obs=2820
  INTEGER, PARAMETER :: id_t_obs=3073
  INTEGER, PARAMETER :: id_q_obs=3330
  INTEGER, PARAMETER :: id_rh_obs=3331
  INTEGER, PARAMETER :: id_ps_obs=14593
  INTEGER, PARAMETER :: id_z_obs=2567
  INTEGER, PARAMETER :: id_s_obs=3332
  INTEGER, PARAMETER :: id_rain_obs=9999
  INTEGER, PARAMETER :: id_bt_obs=21023

  ! ID for NICAM [3D]
  INTEGER, PARAMETER :: id_pres_nicam=1
  INTEGER, PARAMETER :: id_u_nicam=2
  INTEGER, PARAMETER :: id_v_nicam=3
  INTEGER, PARAMETER :: id_temp_nicam=4
  INTEGER, PARAMETER :: id_qvap_nicam=5
  INTEGER, PARAMETER :: id_qcld_nicam=6

  ! ID for NICAM [2D]
  INTEGER, PARAMETER :: id_tsfc_nicam=1
  INTEGER, PARAMETER :: id_u10m_nicam=2
  INTEGER, PARAMETER :: id_v10m_nicam=3
  INTEGER, PARAMETER :: id_cldw_nicam=4
  INTEGER, PARAMETER :: id_surp_nicam=5
  INTEGER, PARAMETER :: id_qv2m_nicam=6
  INTEGER, PARAMETER :: id_te2m_nicam=7
  INTEGER, PARAMETER :: id_rain_nicam=8

  INTEGER, PARAMETER :: id_u=1
  INTEGER, PARAMETER :: id_v=2
  INTEGER, PARAMETER :: id_t=3
  INTEGER, PARAMETER :: id_q=4
  INTEGER, PARAMETER :: id_ps=5
  INTEGER, PARAMETER :: id_rain=6

  REAL(4), ALLOCATABLE :: icodata4_3d(:,:,:,:)
  REAL(4), ALLOCATABLE :: icodata4_2d(:,:,:,:)

  REAL(8), SAVE :: alpha = 0.0d0
  REAL(8), SAVE :: dalpha = 0.0d0
  REAL(8), SAVE :: ocean_value = 0.0d0
  REAL(8), ALLOCATABLE, SAVE :: icolat(:,:)
  REAL(8), ALLOCATABLE, SAVE :: icolon(:,:)

  REAL(8), ALLOCATABLE, SAVE :: lat_model(:)
  REAL(8), ALLOCATABLE, SAVE :: lon_model(:)
  REAL(4), ALLOCATABLE, SAVE :: land_model(:)
  INTEGER, ALLOCATABLE, SAVE :: idirec_model(:)

  INTEGER, ALLOCATABLE, SAVE :: land(:)
  INTEGER, ALLOCATABLE, SAVE :: idirec(:,:)
  REAL(8), ALLOCATABLE, SAVE :: icoland(:,:)

  INTEGER, ALLOCATABLE, SAVE :: obsnum(:)      ! number

  INTEGER, PARAMETER :: nlon_co2=144
  INTEGER, PARAMETER :: nlat_co2=73
  INTEGER, PARAMETER :: nlev_co2=60
  REAL(4), SAVE :: lon_co2(nlon_co2), lat_co2(nlat_co2)
  REAL(4), SAVE :: ak(nlev_co2), bk(nlev_co2)
  REAL(4), SAVE :: psurf_co2(nlon_co2,nlat_co2)
  REAL(4), SAVE :: pres_co2(nlon_co2,nlat_co2,nlev_co2)
  REAL(4), SAVE :: co2_org(nlon_co2,nlat_co2,nlev_co2)

  CHARACTER(ADM_NSYS), SAVE :: sfc_type = 'RIGID'

  LOGICAL, SAVE  :: start_mem_zero = .true.		! [add] Koji 2017.09.08  
  LOGICAL, SAVE  :: flush_text     = .false.	! [add] Koji 2017.09.08  

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE set_coef_interpolate(nsite, lon, lat,                     &
           inprc, l_index, n1_index, n2_index, n3_index, w1, w2, w3 )
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: nsite
  REAL(8), INTENT(INOUT) :: lon(nsite)
  REAL(8), INTENT(INOUT) :: lat(nsite)
  LOGICAL, INTENT(INOUT) :: inprc(nsite)
  INTEGER, INTENT(INOUT) :: l_index(nsite)
  INTEGER, INTENT(INOUT) :: n1_index(nsite)
  INTEGER, INTENT(INOUT) :: n2_index(nsite)
  INTEGER, INTENT(INOUT) :: n3_index(nsite)
  REAL(8), INTENT(INOUT) :: w1(nsite)
  REAL(8), INTENT(INOUT) :: w2(nsite)
  REAL(8), INTENT(INOUT) :: w3(nsite)
  !INTEGER, INTENT(IN), OPTIONAL :: elem(nsite)
  !REAL(8), INTENT(IN), OPTIONAL :: elev(nsite)
  !REAL(8), INTENT(INOUT), OPTIONAL :: klev(2,nsite)
  !REAL(8), INTENT(INOUT), OPTIONAL :: kfact(nsite)

  INTEGER :: ierr, i
  INTEGER :: max_num_latlon
  INTEGER :: sum_max_num_latlon

  WRITE(ADM_LOG_FID,*) "minval(lat(:)), maxval(lat(:))"
  WRITE(ADM_LOG_FID,*) minval(lat(:)), maxval(lat(:))
  WRITE(ADM_LOG_FID,*) minval(lon(:)), maxval(lon(:))

  lon(:)=lon(:)*DEG2RAD
  lat(:)=lat(:)*DEG2RAD

  inprc(:)= .false.
  l_index(:) = 0
  n1_index(:) = -1
  n2_index(:) = -1
  n3_index(:) = -1
  w1(:) = 0.0d0
  w2(:) = 0.0d0
  w3(:) = 0.0d0

  CALL setup_ico2latlon_mapping('GET_NUM', nsite, lat, lon, max_num_latlon, &
                  inprc, l_index, n1_index, n2_index, n3_index, w1, w2, w3 )

  CALL MPI_ALLREDUCE(max_num_latlon, sum_max_num_latlon, 1,      &
                     MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  WRITE(ADM_LOG_FID,*) 'nsite= ', nsite
  WRITE(ADM_LOG_FID,*) 'sum_max_num_latlon= ', sum_max_num_latlon
  IF(flush_text) FLUSH(ADM_LOG_FID)
  IF(sum_max_num_latlon /= nsite) THEN
     WRITE(ADM_LOG_FID,*) 'ERROR!: max_num_latlon /= nsite'
     WRITE(ADM_LOG_FID,*) 'max_num_latlon:', sum_max_num_latlon
     WRITE(ADM_LOG_FID,*) 'STOP!'
     CALL ADM_proc_stop
  END IF
  CALL setup_ico2latlon_mapping('SET_INDEX', nsite, lat, lon, max_num_latlon, &
                inprc, l_index, n1_index, n2_index, n3_index, w1, w2, w3 )

  WRITE(ADM_LOG_FID,*) 'END setup_ico2latlon_mapping'
  IF(flush_text) FLUSH(ADM_LOG_FID)
  !WRITE(ADM_LOG_FID,*) 'l_index', l_index
  !WRITE(ADM_LOG_FID,*) 'w1', w1
  !WRITE(ADM_LOG_FID,*) 'w2', w2
  !WRITE(ADM_LOG_FID,*) 'w3', w3
  !DO i = 1, nsite
  !  IF(l_index(i) /= 0) THEN
  !    WRITE(ADM_LOG_FID,'(i5,3f12.6,3i5)') l_index(i), w1(i), w2(i), w3(i), &
  !                                     n1_index(i), n2_index(i), n3_index(i)
  !  END IF
  !END DO

  DO i = 1, nsite
    IF(l_index(i) /= 0) THEN
      inprc(i) = .true.
    END IF
  END DO
  IF(flush_text) FLUSH(ADM_LOG_FID)

  !IF( PRESENT(elev) ) THEN
  !  CALL getklev(nsite, elem, elev, l_index, n1_index, n2_index, n3_index, &
  !               w1, w2, w3, inprc, klev, kfact)
  !END IF

END SUBROUTINE set_coef_interpolate
!------------------------------------------------------------------------------
SUBROUTINE setup_ico2latlon_mapping( what_is_done, nsite, lat, lon, max_num_latlon,&
                         inprc, l_index, n1_index, n2_index, n3_index, w1, w2, w3 )
  !
  !**************************************************
  ! + Imported from mod_latlon.f90
  !***************************************************
  !
  USE mod_misc, ONLY :               &
       MISC_triangle_area
  USE mod_adm, ONLY :              &
       ADM_LOG_FID,              &
       ADM_NSYS,                 &
       ADM_MAXFNAME,             &
       ADM_IooJoo_nmax,          &
       ADM_IooJoo,               &
       ADM_GIoJo,                &
       ADM_GIpJo,                &
       ADM_GIpJp,                &
       ADM_GIoJp,                &
       ADM_GImJo,                &
       ADM_GIoJm,                &
       ADM_GImJm,                &
       ADM_KNONE,                &
       ADM_TI,ADM_TJ,            &
       ADM_lall
  USE mod_grd, ONLY :               &
       GRD_x,GRD_x_pl,           &
       GRD_XDIR,                 &
       GRD_YDIR,                 &
       GRD_ZDIR,                 &
       GRD_rscale
  USE mod_gmtr, ONLY :            &
       GMTR_P_var, GMTR_P_var_pl, &
       GMTR_P_LAT,                &
       GMTR_P_LON,                &
       GMTR_polygon_type
  USE mod_cnst, ONLY :               &
       CNST_PI
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what_is_DOne
  INTEGER, INTENT(IN) :: nsite
  REAL(8), INTENT(IN) :: lon(nsite) ! in radian 
  REAL(8), INTENT(IN) :: lat(nsite) ! in radian 
  INTEGER, INTENT(OUT) :: max_num_latlon
  LOGICAL, INTENT(INOUT) :: inprc(nsite)
  INTEGER, INTENT(INOUT) :: l_index(nsite)
  INTEGER, INTENT(INOUT) :: n1_index(nsite)
  INTEGER, INTENT(INOUT) :: n2_index(nsite)
  INTEGER, INTENT(INOUT) :: n3_index(nsite)
  REAL(8), INTENT(INOUT) :: w1(nsite)
  REAL(8), INTENT(INOUT) :: w2(nsite)
  REAL(8), INTENT(INOUT) :: w3(nsite)
  !
  REAL(8) :: r1(3), r2(3), r3(3), r0(3)
  REAL(8) :: v12(3), v23(3), v31(3), v10(3), v20(3), v30(3)
  REAL(8) :: nvec(3), v12xv10(3),v23xv20(3),v31xv30(3)
  REAL(8):: judge12, judge23, judge31, rf, rn
  INTEGER :: n, l, i, t
  !
  REAL(8) :: latmin_l, latmax_l
  REAL(8) :: lat1
  REAL(8) :: lat2
  REAL(8) :: lat3
  REAL(8) :: lat4
  REAL(8) :: lonmin_l, lonmax_l
  REAL(8) :: lon1, lon2, lon3, lon4
  REAL(8) :: len, area_total
  !
  REAL(8) :: coslat(nsite), sinlat(nsite), coslon(nsite), sinlon(nsite)
  LOGICAL :: log_j(nsite)
  !
  !    ! S.Iga060210 =>
  !    REAL(8) :: epsi=1d-8
  !    ! S.Iga060210 >=
  !
  ! S.Iga060212 =>
  !-- infinite simal values to prevent 'missing grid' and 'NaN weight'
  !---(1) IF you want to prevent the missing grids on the edge of triangles,
  !       give a value to epsi_inpro. However, you should take a risk to
  !       count-in grids even out of the triangle.
  !       In my experience, zero is OK at least for GLEVEL =< 11
  !---(2) IF you want to prevent missing grids on the apex of triangles,
  !       give a value to epsi_grid.
  !       In my experience, this needs a value when epsi_inpro=0.
  !---(3) IF you want to prevent 'NaN weight', give a value to epsi_area.
  !       In my experience, zero is OK at least for GLEVEL =< 11
  !
  REAL(8) :: epsi_inpro =0     ! marginal value for inner products USEd to
                                ! calcuate judge
  REAL(8) :: epsi_grid  =0 ! marginal square near grid points (in radian)
  !REAL(8) :: epsi_grid  =1d-8 ! marginal square near grid points (in radian)
  REAL(8) :: epsi_area  =0     ! marginal value for inner products used in
                               ! misc_triangle_area
  ! S.Iga060212 >=
  !
  DO i=1,nsite
     coslon(i)=cos(lon(i))
     sinlon(i)=sin(lon(i))
     coslat(i)=cos(lat(i))
     sinlat(i)=sin(lat(i))
  END DO
  !
  !
  max_num_latlon = 0
  !  max_num_latlon2(:) = 0
  !
  DO l=1,ADM_lall
     DO n=1, ADM_IooJoo_nmax
        lat1=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GMTR_P_LAT)
        lat2=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GMTR_P_LAT)
        lat3=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GMTR_P_LAT)
        lat4=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GMTR_P_LAT)
        lon1=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GMTR_P_LON)
        lon2=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GMTR_P_LON)
        lon3=GMTR_P_var(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GMTR_P_LON)
        lon4=GMTR_P_var(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GMTR_P_LON)
        !
        DO t=ADM_TI,ADM_TJ
           !
           !--- construct triagular vertice by clockwise way from the origin
           IF(t==ADM_TI) THEN
              r1(1) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
              r1(2) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
              r1(3) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
              !
              r2(1)= GRD_x(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
              r2(2)= GRD_x(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
              r2(3)= GRD_x(ADM_IooJoo(n,ADM_GIpJo),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
              !
              r3(1) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
              r3(2) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
              r3(3) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
              !
           ELSE !--- ADM_TJ 
              r1(1) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
              r1(2) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
              r1(3) = GRD_x(ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
              !
              r2(1) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
              r2(2) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
              r2(3) = GRD_x(ADM_IooJoo(n,ADM_GIpJp),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
              !
              r3(1)= GRD_x(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GRD_XDIR)/GRD_rscale
              r3(2)= GRD_x(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GRD_YDIR)/GRD_rscale
              r3(3)= GRD_x(ADM_IooJoo(n,ADM_GIoJp),ADM_KNONE,l,GRD_ZDIR)/GRD_rscale
              !
           END IF
           !
           !--- 
           latmin_l=min(lat1,lat2,lat3,lat4)
           latmax_l=max(lat1,lat2,lat3,lat4)
           !
           lonmin_l=min(lon1,lon2,lon3,lon4)
           lonmax_l=max(lon1,lon2,lon3,lon4)
           IF (lonmax_l-lonmin_l > CNST_PI) THEN
              IF (lon1 < 0 ) lon1=lon1+CNST_PI+CNST_PI
              IF (lon2 < 0 ) lon2=lon2+CNST_PI+CNST_PI
              IF (lon3 < 0 ) lon3=lon3+CNST_PI+CNST_PI
              IF (lon4 < 0 ) lon4=lon4+CNST_PI+CNST_PI
              lonmin_l=min(lon1,lon2,lon3,lon4)
              lonmax_l=max(lon1,lon2,lon3,lon4)
           END IF
           ! S.Iga060210 =>
           lonmin_l = lonmin_l - epsi_grid
           lonmax_l = lonmax_l + epsi_grid
           latmin_l = latmin_l - epsi_grid
           latmax_l = latmax_l + epsi_grid
           ! S.Iga060210 >=
           !
           DO i=1, nsite
              log_j(i)=.not.( (latmin_l<=lat(i)).and.(lat(i)<=latmax_l) )
           END DO
           !
           DO i=1, nsite
              IF (.not.( ((lonmin_l<=lon(i)).and.(lon(i)<=lonmax_l))&
                   .or.((lonmin_l<=lon(i)+CNST_PI+CNST_PI).and.(lon(i)+CNST_PI+CNST_PI<=lonmax_l)))   ) THEN
                 CYCLE
              END IF
              !
              ! Somehow when this log_j(j) is USEd instead,  &
              ! results becomes slightly dIFferent in order  &
              ! of 1e-13 at least in Linux. (S.Iga)
              !                   IF (log_j(j)) THEN
              IF (.not.( (latmin_l<=lat(i)).and.(lat(i)<=latmax_l) )) THEN
                 CYCLE
              END IF
              !
              !--- set side-vectors by clockwise way from the origin
              v12(1)=r2(1)-r1(1)
              v12(2)=r2(2)-r1(2)
              v12(3)=r2(3)-r1(3)
              !
              v23(1)=r3(1)-r2(1)
              v23(2)=r3(2)-r2(2)
              v23(3)=r3(3)-r2(3)
              !
              v31(1)=r1(1)-r3(1)
              v31(2)=r1(2)-r3(2)
              v31(3)=r1(3)-r3(3)
              !
              !
              !--- calculate normal unit vector to the plane triagle
              nvec(1)=v12(2)*v23(3)-v12(3)*v23(2)
              nvec(2)=v12(3)*v23(1)-v12(1)*v23(3)
              nvec(3)=v12(1)*v23(2)-v12(2)*v23(1)
              len = sqrt(&
                   nvec(1)*nvec(1)&
                   +nvec(2)*nvec(2)&
                   +nvec(3)*nvec(3))
              nvec(1:3)=nvec(1:3)/len
              !
              !--- target latlon point on the sphere
              r0(1)=coslat(i)*coslon(i)
              r0(2)=coslat(i)*sinlon(i)
              r0(3)=sinlat(i)
              !
              !--- remove the CASE inner product is negative
              IF((r0(1)*r1(1)+r0(2)*r1(2)+r0(3)*r1(3))<0.0D0) THEN
                 CYCLE
              END IF
              !
              !--- mapping r0 to the plane(r1,r2,r3)
              !
              !------ distance from origin to a plane with r0.
              rf=(r0(1)*nvec(1)+r0(2)*nvec(2)+r0(3)*nvec(3))
              !
              !------ distance from origin to a plane with r1 or (r2,r3).
              rn=(r1(1)*nvec(1)+r1(2)*nvec(2)+r1(3)*nvec(3))
              !
              !------ mapping r0 
              r0(1)=r0(1)*(rn/rf)
              r0(2)=r0(2)*(rn/rf)
              r0(3)=r0(3)*(rn/rf)
              !
              !--- calculate vectors from triangler points 
              !--- to the target point
              v10(1)=r1(1)-r0(1)
              v10(2)=r1(2)-r0(2)
              v10(3)=r1(3)-r0(3)
              !
              v20(1)=r2(1)-r0(1)
              v20(2)=r2(2)-r0(2)
              v20(3)=r2(3)-r0(3)
              v30(1)=r3(1)-r0(1)
              v30(2)=r3(2)-r0(2)
              v30(3)=r3(3)-r0(3)
              !
              !
              v12xv10(1)=v12(2)*v10(3)-v12(3)*v10(2)
              v12xv10(2)=v12(3)*v10(1)-v12(1)*v10(3)
              v12xv10(3)=v12(1)*v10(2)-v12(2)*v10(1)
              !
              v23xv20(1)=v23(2)*v20(3)-v23(3)*v20(2)
              v23xv20(2)=v23(3)*v20(1)-v23(1)*v20(3)
              v23xv20(3)=v23(1)*v20(2)-v23(2)*v20(1)
              !
              v31xv30(1)=v31(2)*v30(3)-v31(3)*v30(2)
              v31xv30(2)=v31(3)*v30(1)-v31(1)*v30(3)
              v31xv30(3)=v31(1)*v30(2)-v31(2)*v30(1)
              !
              judge12 = &
                   nvec(1)*v12xv10(1) &
                   +nvec(2)*v12xv10(2) &
                   +nvec(3)*v12xv10(3)
              judge23 = &
                   nvec(1)*v23xv20(1) &
                   +nvec(2)*v23xv20(2) &
                   +nvec(3)*v23xv20(3)
              judge31 = &
                   nvec(1)*v31xv30(1) &
                   +nvec(2)*v31xv30(2) &
                   +nvec(3)*v31xv30(3)
              !
              IF(t==ADM_TI) THEN
                 IF (      (judge12<=0.0D0 + epsi_inpro)& ! S.Iga add epsi_inpro  (060212)
                      .and.(judge23<= 0.0D0 + epsi_inpro)&
                      .and.(judge31<=0.0D0 + epsi_inpro) ) THEN
                    SELECT CASE(TRIM(what_is_DOne))
                    CASE('GET_NUM')
                       max_num_latlon = max_num_latlon + 1
                       !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                    CASE('SET_INDEX')
                       !   max_num_latlon = max_num_latlon + 1
                       !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                       !   lon_index(max_num_latlon) = i
                       !   lat_index(max_num_latlon) = j
                       l_index(i) = l
                       !   t_index(max_num_latlon) = t
                       n1_index(i) = ADM_IooJoo(n,ADM_GIoJo)
                       n2_index(i) = ADM_IooJoo(n,ADM_GIpJo)
                       n3_index(i) = ADM_IooJoo(n,ADM_GIpJp)
                       w1(i) = MISC_triangle_area( &
                            r0(1:3),r2(1:3),r3(1:3),& ! S.Iga add epsi_area (060212)
                            GMTR_polygon_type, 1.0D0,critical=epsi_area)
                       w2(i) = MISC_triangle_area( &
                            r0(1:3),r3(1:3),r1(1:3),& ! S.Iga add epsi_area (060212)
                            GMTR_polygon_type, 1.0D0,critical=epsi_area)
                       w3(i) = MISC_triangle_area( &
                            r0(1:3),r1(1:3),r2(1:3),& ! S.Iga add epsi_area (060212)
                            GMTR_polygon_type, 1.0D0,critical=epsi_area)
                       area_total &
                            = w1(i)&
                            + w2(i)&
                            + w3(i)
                       w1(i)=w1(i)/area_total
                       w2(i)=w2(i)/area_total
                       w3(i)=w3(i)/area_total
                       !
                    CASE DEFAULT
                    END SELECT

                    CYCLE
                    ! S.Iga060210=>
                    ! Deal with the points near 'r1'
                 ELSEIF (ABS(v10(1)) <= epsi_grid ) THEN
                    IF (ABS(v10(2)) <= epsi_grid ) THEN
                       IF (ABS(v10(3)) <= epsi_grid ) THEN
                          SELECT CASE(TRIM(what_is_DOne))
                          CASE('GET_NUM')
                             max_num_latlon = max_num_latlon + 1
                             !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                          CASE('SET_INDEX')
                             !   max_num_latlon = max_num_latlon + 1
                             !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                             !   lon_index(max_num_latlon) = i
                             !   lat_index(max_num_latlon) = j
                             l_index(i) = l
                             !   t_index(max_num_latlon) = t
                             n1_index(i) = ADM_IooJoo(n,ADM_GIoJo)
                             n2_index(i) = ADM_IooJoo(n,ADM_GIpJo)
                             n3_index(i) = ADM_IooJoo(n,ADM_GIpJp)
                             w1(i) = 1d0
                             w2(i) = 0
                             w3(i) = 0
                          END SELECT
                       END IF
                    END IF
                    ! S.Iga060210>=
                 END IF
              ELSE !--- ADM_TJ
                 IF (      (judge12< 0.0D0 + epsi_inpro)& ! S.Iga add epsi_inpro  (060212)
                      .and.(judge23< 0.0D0 + epsi_inpro)&
                      .and.(judge31<=0.0D0 + epsi_inpro) ) THEN
                    SELECT CASE(TRIM(what_is_DOne))
                    CASE('GET_NUM')
                       max_num_latlon = max_num_latlon + 1
                       !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                    CASE('SET_INDEX')
                       !   max_num_latlon = max_num_latlon + 1
                       !   max_num_latlon2(l) = max_num_latlon2(l) + 1
                       !   lon_index(max_num_latlon) = i
                       !   lat_index(max_num_latlon) = j
                       l_index(i) = l
                       !   t_index(max_num_latlon) = t
                       n1_index(i) = ADM_IooJoo(n,ADM_GIoJo)
                       n2_index(i) = ADM_IooJoo(n,ADM_GIpJp)
                       n3_index(i) = ADM_IooJoo(n,ADM_GIoJp)
                       w1(i) = MISC_triangle_area( &
                            r0(1:3),r2(1:3),r3(1:3),& ! S.Iga add epsi_area (060212)
                            GMTR_polygon_type, 1.0D0,critical=epsi_area)
                       w2(i) = MISC_triangle_area( &
                            r0(1:3),r3(1:3),r1(1:3),& ! S.Iga add epsi_area (060212)
                            GMTR_polygon_type, 1.0D0,critical=epsi_area)
                       w3(i) = MISC_triangle_area( &
                            r0(1:3),r1(1:3),r2(1:3),& ! S.Iga add epsi_area (060212)
                            GMTR_polygon_type, 1.0D0,critical=epsi_area)
                       area_total &
                            = w1(i)&
                            + w2(i)&
                            + w3(i)
                       w1(i)=w1(i)/area_total
                       w2(i)=w2(i)/area_total
                       w3(i)=w3(i)/area_total
                       !
                    CASE DEFAULT
                    END SELECT
                    CYCLE
                 END IF
              END IF
           END DO
        END DO
     END DO
  END DO
  !
END subroutine setup_ico2latlon_mapping
!-----------------------------------------------------------------------
SUBROUTINE getklev(nsite, elem, elev,                            &
                   inprc, l_index, n1_index, n2_index, n3_index, &
                   w1, w2, w3, klev, kfact)
  !
  use mod_adm, only : &
       ADM_kmin,      &
       ADM_kmax,      &
       ADM_vlayer
  use mod_grd, only : &
       GRD_vz,   &
       GRD_Z,    &
       GRD_zs,   &
       GRD_ZSFC
  !
  implicit none
  !
  INTEGER :: nsite
  INTEGER :: elem(nsite)
  REAL(8) :: elev(nsite)
  LOGICAL, INTENT(IN) :: inprc(nsite)
  INTEGER, INTENT(IN) :: l_index(nsite)
  INTEGER, INTENT(IN) :: n1_index(nsite)
  INTEGER, INTENT(IN) :: n2_index(nsite)
  INTEGER, INTENT(IN) :: n3_index(nsite)
  REAL(8), INTENT(IN) :: w1(nsite)
  REAL(8), INTENT(IN) :: w2(nsite)
  REAL(8), INTENT(IN) :: w3(nsite)
  REAL(8), INTENT(INOUT) :: klev(2,nsite)
  REAL(8), INTENT(INOUT) :: kfact(nsite)

  REAL(8) :: z(ADM_gall,ADM_vlayer,ADM_lall)
  INTEGER :: s, k
  REAL(8) :: alt
  REAL(8) :: zs(ADM_gall,ADM_lall)
  CHARACTER(12) :: name
  REAL(8) :: zz, zz0
  REAL(8) :: pp, pp0
  REAL(8) :: tt, tt0, tv
  REAL(8) :: qq, qq0
  REAL(8) :: t1, t2
  REAL(8) :: tmp_elev
  REAL(8), PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(8), PARAMETER :: undef = -999.9d0

  z(:,1:ADM_vlayer,:) = GRD_vz(:,ADM_kmin:ADM_kmax,:,GRD_Z)
  zs(:,:) = GRD_zs(:,ADM_KNONE,:,GRD_ZSFC)

  DO s=1, nsite
    IF(inprc(s)) THEN
    !WRITE(ADM_LOG_FID,'(i5,3f12.6,3i5)') l_index(s), w1(s), w2(s), w3(s), &
    !                                   n1_index(s), n2_index(s), n3_index(s)
    !IF(flush_text) FLUSH(ADM_LOG_FID)
      IF( elem(s) == id_ps_obs ) THEN
        alt = w1(s)*zs(n1_index(s),l_index(s)) &
            + w2(s)*zs(n2_index(s),l_index(s)) &
            + w3(s)*zs(n3_index(s),l_index(s))
        pp  = w1(s)*icodata4_2d(n1_index(s),1,l_index(s),id_surp_nicam) &
            + w2(s)*icodata4_2d(n2_index(s),1,l_index(s),id_surp_nicam) &
            + w3(s)*icodata4_2d(n3_index(s),1,l_index(s),id_surp_nicam)
        tt  = w1(s)*icodata4_2d(n1_index(s),1,l_index(s),id_te2m_nicam) &
            + w2(s)*icodata4_2d(n2_index(s),1,l_index(s),id_te2m_nicam) &
            + w3(s)*icodata4_2d(n3_index(s),1,l_index(s),id_te2m_nicam)
        qq  = w1(s)*icodata4_2d(n1_index(s),1,l_index(s),id_qv2m_nicam) &
            + w2(s)*icodata4_2d(n2_index(s),1,l_index(s),id_qv2m_nicam) &
            + w3(s)*icodata4_2d(n3_index(s),1,l_index(s),id_qv2m_nicam)
        IF( alt-elev(s) /= 0 .and. ABS(alt-elev(s)) < 1000.0d0 ) THEN
          tv=tt*(1.0d0+0.608d0*qq)
          kfact(s) = (tv/(tv+gamma*(elev(s)-alt)))**(9.81d0/(gamma*287.d0))
        ELSE
          kfact(s) = -999.0d0
        END IF
        klev(:,s)=1
      ELSE
        tmp_elev=log(elev(s))
        pp0=undef
        DO k = 1, ADM_vlayer
          pp = w1(s)*icodata4_3d(n1_index(s),k,l_index(s),id_pres_nicam) &
             + w2(s)*icodata4_3d(n2_index(s),k,l_index(s),id_pres_nicam) &
             + w3(s)*icodata4_3d(n3_index(s),k,l_index(s),id_pres_nicam)
          IF ( pp < tmp_elev ) EXIT
          pp0=pp
        END DO
        IF( k == 1 ) THEN
          klev(1,s)=1
          klev(2,s)=1
          kfact(s)=1.0d0
          IF( elem(s) == id_t_obs ) THEN
            kfact(s) = -999.0d0
          END IF
        ELSE IF ( k >= ADM_vlayer ) THEN
          klev(1,s)=ADM_vlayer
          klev(2,s)=ADM_vlayer
          kfact(s) = -999.0d0
        ELSE
          klev(1,s)=k-1
          klev(2,s)=k
          kfact(s) = ( pp - tmp_elev ) / ( pp - pp0 )
        END IF ! k
      END IF ! elem
    END IF ! inprc
  END DO
  !
  return
END SUBROUTINE getklev
!------------------------------------------------------------------------------
SUBROUTINE calc_weighting_function(kmax, nobs, maxtvsch, pres, tran, &
                          weight, weight_maxlev)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kmax
  INTEGER, INTENT(IN)  :: nobs
  INTEGER, INTENT(IN)  :: maxtvsch
  REAL(8), INTENT(IN)  :: pres(kmax,nobs)
  REAL(8), INTENT(IN)  :: tran(kmax,maxtvsch,nobs)
  REAL(8)  :: weight(kmax,maxtvsch,nobs)
  INTEGER  :: weight_maxlev(maxtvsch,nobs)
  !REAL(8), INTENT(INOUT) :: weight(kmax,maxtvsch,nobs)
  !INTEGER, INTENT(INOUT) :: weight_maxlev(maxtvsch,nobs)

  INTEGER :: n, ichan, k
  REAL(8) :: tmp_lev
  REAL(8) :: logp(kmax)

  DO n = 1, nobs
    logp(:)=LOG(pres(:,n))
    DO ichan = 1, maxtvsch 
      weight(1,ichan,n)=-(tran(2,ichan,n)-tran(1,ichan,n))/&
                         (logp(kmax-1)-logp(kmax))
      DO k = 2, kmax-1
        weight(k,ichan,n)=-(tran(k,ichan,n)-tran(k-1,ichan,n))/&
                           (logp(kmax-k)-logp(kmax-k+2))
      END DO
      weight(kmax,ichan,n)=-(tran(k,ichan,n)-tran(k-1,ichan,n))/&
                            (logp(1)-logp(2))
      
    END DO
  END DO

  DO n = 1, nobs
  DO ichan = 1, maxtvsch 
    tmp_lev=-1.0
    DO k = 1, kmax
      IF( weight(k,ichan,n) > tmp_lev ) THEN
        weight_maxlev(ichan,n)=kmax-k+1
        tmp_lev=weight(k,ichan,n)
      END IF
    END DO
  END DO
  END DO

END SUBROUTINE calc_weighting_function
!------------------------------------------------------------------------------
SUBROUTINE phys2ijk(nlon,nlat,nlev,p_full,elem,lon,lat,rlon,rlat,rlev,ri,rj,rk)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlon, nlat, nlev
  REAL(4),INTENT(IN) :: p_full(nlon,nlat,nlev)
  REAL(4),INTENT(IN) :: elem
  REAL(4),INTENT(IN) :: lon(nlon)
  REAL(4),INTENT(IN) :: lat(nlat)
  REAL(4),INTENT(IN) :: rlon
  REAL(4),INTENT(IN) :: rlat
  REAL(4),INTENT(IN) :: rlev ! pressure levels
  REAL(4),INTENT(OUT) :: ri
  REAL(4),INTENT(OUT) :: rj
  REAL(4),INTENT(OUT) :: rk
  REAL(4) :: aj,ak
  REAL(4) :: lnps(nlon,nlat)
  REAL(4) :: plev(nlev)
  INTEGER :: i,j,k
!
! rlon -> ri
!
  IF(rlon == 0.0 .OR. rlon == 360.0) THEN
    ri = REAL(nlon+1,4)
  ELSE
    ri = rlon / 360.0d0 * REAL(nlon,4) + 1.0d0
  END IF
  IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) RETURN
!
! rlat -> rj
!
  DO j=1,nlat
    IF(rlat < lat(j)) EXIT
  END DO
  IF(j == 1) THEN
    rj = (rlat + 90.0d0) / (lat(1) + 90.0d0)
  ELSE IF(j == nlat+1) THEN
    aj = (rlat - lat(nlat)) / (90.0d0 - lat(nlat))
    rj = REAL(nlat,4) + aj
  ELSE
    aj = (rlat - lat(j-1)) / (lat(j) - lat(j-1))
    rj = REAL(j-1,4) + aj
  END IF
  IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) RETURN
!
! rlev -> rk
!
  IF(NINT(elem) == id_ps_obs) THEN ! surface pressure observation
    rk = 0.0d0
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    DO k=1,nlev
      IF(i <= nlon) THEN
        lnps(i-1:i,j-1:j) = LOG(p_full(i-1:i,j-1:j,k))
      ELSE
        lnps(i-1,j-1:j) = LOG(p_full(i-1,j-1:j,k))
        lnps(1,j-1:j) = LOG(p_full(1,j-1:j,k))
      END IF
      CALL itpl_2d(nlon,nlat,lnps,ri,rj,plev(k))
    END DO
    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! find rk
    !
    DO k=2,nlev-1
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,4) + ak
    IF(rk<1) rk=1
  END IF

  RETURN
END SUBROUTINE phys2ijk
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(nlon,nlat,var,ri,rj,var5)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlon, nlat
  REAL(4), INTENT(IN) :: var(nlon,nlat)
  REAL(4), INTENT(IN) :: ri
  REAL(4), INTENT(IN) :: rj
  REAL(4), INTENT(OUT) :: var5
  REAL(4) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,4)
  j = CEILING(rj)
  aj = rj - REAL(j-1,4)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_3d(nlon,nlat,nlev,var,ri,rj,rk,var5)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlon, nlat, nlev
  REAL(4),INTENT(IN) :: var(nlon,nlat,nlev)
  REAL(4),INTENT(IN) :: ri
  REAL(4),INTENT(IN) :: rj
  REAL(4),INTENT(IN) :: rk
  REAL(4),INTENT(OUT) :: var5
  REAL(4) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,4)
  j = CEILING(rj)
  aj = rj - REAL(j-1,4)
  k = CEILING(rk)
  ak = rk - REAL(k-1,4)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
  ELSE
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
  END IF

  RETURN
END SUBROUTINE itpl_3d

END MODULE mod_obsope_common
