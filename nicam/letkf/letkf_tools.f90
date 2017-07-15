MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with SPEEDY
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!
!=======================================================================
  USE common
  USE common_mpi
  USE common_nicam
  USE common_mpi_nicam
  USE common_letkf
  USE letkf_obs
  USE mod_adm

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf

  INTEGER,SAVE :: nobstotal

  !REAL(r_size),PARAMETER :: cov_infl_mul = 1.05d0 !multiplicative inflation
  !==> kotsuki 20160311 
  REAL(r_size), SAVE :: cov_infl_mul       = -1.05d0 ! multiplicative inflation
  REAL(r_size), SAVE :: glb_infl_mul       = -1.05d0 ! multiplicative inflation (global constatnt)
  REAL(r_size), SAVE :: RELAX_ALPHA        = 0.0d0   ! relaxation to prior perturbation ; RTPP (Zhang et al. 2005)
  REAL(r_size), SAVE :: RELAX_ALPHA_SPREAD = 0.0d0   ! relaxation to prior spread       ; RTPS (Whitaker and Hamill 2012)
! > 0: globally constant covariance inflation
! < 0: 3D inflation values input from a GPV file "infl_mul.grd"

  REAL(r_size),PARAMETER :: sp_infl_add = 0.d0 !additive inflation
!TVS  LOGICAL,PARAMETER :: msw_vbc = .FALSE.
  REAL(r_size) :: var_local(nv3d+nv2d,nid_obs)
  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
!SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d,atmp3d,atmp2d)
  IMPLICIT NONE
  CHARACTER(256) :: inflfile='' ! [Change] Koji 20140208
  !CHARACTER(12) :: inflfile='infl_mul.grd'
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_size),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size),ALLOCATABLE :: logpfm(:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: rtimer00, rtimer                                       ! KK
  integer,parameter :: slot=20                                           ! KK
  REAL(r_size) :: ltimer_all(slot,nprocs), ltimer(slot)=0.d0, ltimer00   ! KK
  REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d)
  LOGICAL :: ex
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr
  INTEGER :: nobsl_max
  INTEGER :: nobsl_min

  !!RTPS!! kotsuki 20151216
  REAL(r_size),ALLOCATABLE :: work3da(:,:,:)     !GYL ; kotsuki 20151216
  REAL(r_size),ALLOCATABLE :: work2da(:,:)       !GYL ; kotsuki 20151216
  REAL(r_size) :: transm(nbv,nv3d+nv2d)          !GYL ; kotsuki 20151216
  REAL(r_size) :: transrlx(nbv,nbv)              !GYL ; kotsuki 20151216
  REAL(r_size) :: pa(nbv,nbv,nv3d+nv2d)          !GYL ; kotsuki 20151216
  REAL(r_size) :: www(nbv,nbv)

  !==> add by kotsuki 20160317
  LOGICAL :: FWRITE_MULT = .false.
  LOGICAL :: FWRITE_RTPS = .false.
  namelist / infl_param /   &
    cov_infl_mul,           &
    glb_infl_mul,           &
    RELAX_ALPHA,            &
    RELAX_ALPHA_SPREAD

  open(1,file='time.cnf')
  read(1,nml=infl_param)
  close(1)

  if(myrank == 0) then
    write(*,*) "  INFL_MULT (READ)::", cov_infl_mul        ! kotsuki 20160311
    write(*,*) "  INFL_GLBL (READ)::", glb_infl_mul        ! kotsuki 20160311
    write(*,*) "  ALPH_RTPP (READ)::", RELAX_ALPHA         ! kotsuki 20160311
    write(*,*) "  ALPH_RTPS (READ)::", RELAX_ALPHA_SPREAD  ! kotsuki 20160311
  end if

  var_local(:,:)=1.0d0

  nobsl_max=0
  nobsl_min=9999999
!                                                                        ! KK
  CALL CPU_TIME(rtimer00)                                                ! KK
!                                                                        ! KK
  WRITE(ADM_LOG_FID,'(A)') 'Hello from das_letkf'
  nobstotal = nobs + ntvs
  WRITE(ADM_LOG_FID,'(A,I8)') 'Target observation numbers : NOBS=',nobs,', NTVS=',ntvs
  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(ADM_LOG_FID,'(A)') 'No observation assimilated'
    anal3d = gues3d
    anal2d = gues2d
    RETURN
  END IF
  !
  ! Variable localization
  !
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
  write(ADM_LOG_FID,*) "var_local_n2n"
  write(ADM_LOG_FID,*) var_local_n2n
!
  CALL CPU_TIME(rtimer)                                                                             ! KK
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '##### TIMER(Variable localization in das_letkf):',rtimer,rtimer-rtimer00  ! KK
  ltimer(14) = rtimer-rtimer00                                                                      ! KK
  rtimer00=rtimer                                                                                   ! KK
!
  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)
  DO n=1,nv3d
    DO m=1,nbv
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        END DO
      END DO
    END DO
  END DO
  DO n=1,nv2d
    DO m=1,nbv
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      END DO
    END DO
  END DO
!
  CALL CPU_TIME(rtimer)                                                                             ! KK
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '##### TIMER(Calc. f_perturvation in das_letkf):',rtimer,rtimer-rtimer00   ! KK
  ltimer(15) = rtimer-rtimer00                                                                      ! KK
  rtimer00=rtimer                                                                                   ! KK
!
  !
  ! multiplicative inflation
  !
  IF(cov_infl_mul > 0.0d0 .and. glb_infl_mul > 0.0d0) THEN ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = cov_infl_mul
    work2d = cov_infl_mul
    work3d(:,nlev,:) = 1.01d0
  END IF
  IF( glb_infl_mul <= 0.0d0) THEN ! read global constant parameter
    FWRITE_MULT = .true.
    inflfile=trim(infl_basedir)//'/'//trim(TDATE)//'/adaptmult_anal.txt'
    INQUIRE(FILE=inflfile,EXIST=ex)
    IF(ex) THEN
      open(1,file=trim(inflfile),form="formatted",status="old")
        read(1,*) glb_infl_mul
      close(1)
    ELSE
      WRITE(ADM_LOG_FID,'(2A)') '!!WARNING: no such file exist: ',inflfile
      glb_infl_mul = -1.0d0 * glb_infl_mul
    END IF
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = glb_infl_mul
    work2d = glb_infl_mul
    work3d(:,nlev,:) = 1.01d0
  ELSE
    IF(cov_infl_mul <= 0.0d0) THEN ! 3D parameter values are read-in
      !ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
      !ALLOCATE( work2dg(nlon,nlat,nv2d) )
      ALLOCATE( work3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d) )
      ALLOCATE( work2dg(ADM_gall,ADM_rgn_nmax,nv2d) )
      ALLOCATE( work3d(nij1,nlev,nv3d) )
      ALLOCATE( work2d(nij1,nv2d) )
      inflfile=trim(infl_basedir)//'/'//trim(TDATE)//'/infl_mul.rgn00000'      !KK
      INQUIRE(FILE=inflfile,EXIST=ex)
      IF(ex) THEN
        IF(myrank == 0) THEN
          inflfile=trim(infl_basedir)//'/'//trim(TDATE)//'/infl_mul'           !KK
          WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading..',inflfile
          FLUSH(ADM_LOG_FID)
          CALL read_icogrd_legacy(inflfile,work3dg,work2dg)
        END IF
        CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
        where(work3d(:,:,:) > 1.2d0)  work3d(:,:,:)=1.2d0
        where(work2d(:,:)   > 1.2d0)  work2d(:,:)  =1.2d0
        where(work3d(:,:,:) < 0.95d0) work3d(:,:,:)=0.95d0
        where(work2d(:,:)   < 0.95d0) work2d(:,:)  =0.95d0
      ELSE
        WRITE(ADM_LOG_FID,'(2A)') '!!WARNING: no such file exist: ',inflfile
        FLUSH(ADM_LOG_FID)
        work3d = -1.0d0 * cov_infl_mul
        work2d = -1.0d0 * cov_infl_mul
      END IF
    END IF
  END IF
  !
  ! adaptive rtps ; kotsuki 20160326
  !
  IF(RELAX_ALPHA_SPREAD < 0.0d0) THEN ! read global constant parameter
    FWRITE_RTPS = .true.
    inflfile=trim(infl_basedir)//'/'//trim(TDATE)//'/adaptrtps_anal.txt'
    INQUIRE(FILE=inflfile,EXIST=ex)
    IF(ex) THEN
      open(1,file=trim(inflfile),form="formatted",status="old")
        read(1,*) RELAX_ALPHA_SPREAD
      close(1)
    ELSE
      WRITE(ADM_LOG_FID,'(2A)') '!!WARNING: no such file exist: ',inflfile
      RELAX_ALPHA_SPREAD = -1.0d0 * RELAX_ALPHA_SPREAD
    END IF
  END IF
  if(myrank == 0) then
    write(*,*) "  INFL_MULT (USED)::", cov_infl_mul        ! kotsuki 20160311
    write(*,*) "  INFL_GLBL (USED)::", glb_infl_mul        ! kotsuki 20160311
    write(*,*) "  ALPH_RTPP (USED)::", RELAX_ALPHA         ! kotsuki 20160311
    write(*,*) "  ALPH_RTPS (USED)::", RELAX_ALPHA_SPREAD  ! kotsuki 20160311
    write(*,*) "  FWRITE_RTPS     ::", FWRITE_RTPS
    write(*,*) "  FWRITE_MULT     ::", FWRITE_MULT
  end if


!
  CALL CPU_TIME(rtimer)                                                                             ! KK
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '##### TIMER(Read inflation in das_letkf):',rtimer,rtimer-rtimer00         ! KK
  ltimer(16) = rtimer-rtimer00                                                                      ! KK
  rtimer00=rtimer                                                                                   ! KK
!
  !
  ! p_full for background ensemble mean
  !
  ALLOCATE(logpfm(nij1,nlev))
  logpfm(:,:) = DLOG(mean3d(:,:,iv3d_p))
!
  CALL CPU_TIME(rtimer)                                                                             ! KK
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '##### TIMER(Calc. fcst ensemble mean in das_letkf):',rtimer,rtimer-rtimer00 ! KK
  FLUSH(ADM_LOG_FID)
  ltimer(17) = rtimer-rtimer00                                                                      ! KK
  rtimer00=rtimer
!
  !
  ! RTPS relaxation
  !
  IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN
    ALLOCATE( work3da(nij1,nlev,nv3d) )  ! GYL ; kotsuki 20151216
    ALLOCATE( work2da(nij1,nv2d) )       ! GYL ; kotsuki 20151216
    work3da = 1.0d0                      ! GYL ; kotsuki 20151216
    work2da = 1.0d0                      ! GYL ; kotsuki 20151216
  END IF
  !
  ! MAIN ASSIMILATION LOOP
  !
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  WRITE(ADM_LOG_FID,'(A,I3)') '3d ilev = ',ilev
  DO ij=1,nij1
    ilev = 1
    DO n=1,nv3d
      IF(var_local_n2n(n) < n) THEN
        trans(:,:,n) = trans(:,:,var_local_n2n(n))
        transm(:,n) = transm(:,var_local_n2n(n)) !GYL ; kotsuki 20151216
        IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN     !GYL ; kotsuki 20151216
          pa(:,:,n) = pa(:,:,var_local_n2n(n))   !GYL ; kotsuki 20151216
        END IF                                   !GYL ; kotsuki 20151216
        work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
      ELSE
        CALL obs_local(ij,ilev,n,hdxf,rdiag,rloc,dep,nobsl,logpfm)
        parm = work3d(ij,ilev,n)
        IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL ; kotsuki 20151216
          CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL ; kotsuki 20151216
                          trans(:,:,n),ltimer,                          &       !GYL ; kotsuki 20151216 
                          transm=transm(:,n),pao=pa(:,:,n))                     !GYL ; kotsuki 20151216
        ELSE                                                                    !GYL ; kotsuki 20151216
          CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL ; kotsuki 20151216
                          trans(:,:,n),ltimer,                          &       !GYL ; kotsuki 20151216 
                          transm=transm(:,n))                                   !GYL ; kotsuki 20151216
        END IF
        work3d(ij,ilev,n) = parm
        if(nobsl > nobsl_max) nobsl_max=nobsl
        if(nobsl < nobsl_min) nobsl_min=nobsl
      END IF

      IF(RELAX_ALPHA /= 0.0d0) THEN                                                              !GYL - RTPP 
        CALL weight_RTPP(trans(:,:,n),transrlx)                                                  !GYL
      ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                                  !GYL - RTPS
        CALL weight_RTPS(trans(:,:,n),pa(:,:,n),gues3d(ij,ilev,:,n),transrlx,work3da(ij,ilev,n)) !GYL
      ELSE                                                                                       !GYL
        transrlx = trans(:,:,n)                                                                  !GYL
      END IF

      CALL CPU_TIME(ltimer00)                                                               ! KK
      DO m=1,nbv
        anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
        DO k=1,nbv
          anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,k,n) * (transrlx(k,m) + transm(k,n))                       !GYL
        END DO
      END DO
!
      CALL CPU_TIME(rtimer)                                                                 ! KK
      ltimer(11) = ltimer(11) + rtimer-ltimer00                                             ! KK
      ltimer00 = rtimer                                                                     ! KK
!
    END DO

    DO n=1,nv2d
      IF(var_local_n2n(nv3d+n) <= nv3d) THEN
        trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
        transm(:,nv3d+n) = transm(:,var_local_n2n(nv3d+n))      !GYL ; kotsuki 20151216
        IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                    !GYL ; kotsuki 20151216
          pa(:,:,nv3d+n) = pa(:,:,var_local_n2n(nv3d+n))        !GYL ; kotsuki 20151216
        END IF                                                  !GYL ; kotsuki 20151216
        work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n))
      ELSE IF(var_local_n2n(nv3d+n) <= nv3d+n) THEN
        trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
        transm(:,nv3d+n) = transm(:,var_local_n2n(nv3d+n))      !GYL ; kotsuki 20151216
        IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                    !GYL ; kotsuki 20151216
          pa(:,:,nv3d+n) = pa(:,:,var_local_n2n(nv3d+n))        !GYL ; kotsuki 20151216
        END IF                                                  !GYL ; kotsuki 20151216
        work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)
      ELSE
        CALL obs_local(ij,ilev,nv3d+n,hdxf,rdiag,rloc,dep,nobsl,logpfm)
        parm = work2d(ij,n)
        IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL ; kotsuki 20151216
          CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL ; kotsuki 20151216
                          trans(:,:,n),ltimer,                          &       !GYL ; kotsuki 20151216 
                          transm=transm(:,n),pao=pa(:,:,n))                     !GYL ; kotsuki 20151216
        ELSE                                                                    !GYL ; kotsuki 20151216
          CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL ; kotsuki 20151216
                          trans(:,:,n),ltimer,                          &       !GYL ; kotsuki 20151216 
                          transm=transm(:,n))                                   !GYL ; kotsuki 20151216
        END IF
        work2d(ij,n) = parm
      END IF

      IF(RELAX_ALPHA /= 0.0d0) THEN                                                              !GYL - RTPP
        CALL weight_RTPP(trans(:,:,nv3d+n),transrlx)                                             !GYL
      ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                                  !GYL - RTPS
        CALL weight_RTPS(trans(:,:,nv3d+n),pa(:,:,nv3d+n),gues2d(ij,:,n),transrlx,work2da(ij,n)) !GYL
      ELSE                                                                                       !GYL
        transrlx = trans(:,:,nv3d+n)                                                             !GYL
      END IF       

      DO m=1,nbv
        anal2d(ij,m,n)  = mean2d(ij,n)
        DO k=1,nbv
           anal2d(ij,m,n) = anal2d(ij,m,n) &                          !GYL - sum trans and transm here ; kotsuki 20151216
              & + gues2d(ij,k,n) * (transrlx(k,m) + transm(k,nv3d+n)) !GYL                             ; kotsuki 20151216
        END DO
      END DO
      !!!!!!!!! Not to update
      IF(n>=10) THEN
        DO m=1,nbv
          anal2d(ij,m,n) = mean2d(ij,n) + gues2d(ij,m,n)
        END DO
      END IF
!
      CALL CPU_TIME(rtimer)                                                                 ! KK
      ltimer(11) = ltimer(11) + rtimer-ltimer00                                             ! KK
      ltimer00 = rtimer                                                                     ! KK
    END DO
  END DO ! ij loop                                                                         

 DO ij=1,nij1                                                                             
    DO ilev=2,nlev-1
!     if(ij==1) WRITE(ADM_LOG_FID,'(A,I3,A,I6)') 'ilev = ',ilev, ', ij = ', ij
      !WRITE(ADM_LOG_FID,'(A,I3,A,I6)') 'ilev = ',ilev, ', ij = ', ij
      !FLUSH(ADM_LOG_FID)
      DO n=1,nv3d
        IF(var_local_n2n(n) < n) THEN
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          transm(:,n) = transm(:,var_local_n2n(n)) !GYL ; kotsuki 20151216
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN     !GYL ; kotsuki 20151216
            pa(:,:,n) = pa(:,:,var_local_n2n(n))   !GYL ; kotsuki 20151216
          END IF                                   !GYL ; kotsuki 20151216
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
        ELSE
          CALL obs_local(ij,ilev,n,hdxf,rdiag,rloc,dep,nobsl,logpfm)
          parm = work3d(ij,ilev,n)
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL ; kotsuki 20151216
            CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL ; kotsuki 20151216
                            trans(:,:,n),ltimer,                          &       !GYL ; kotsuki 20151216
                            transm=transm(:,n),pao=pa(:,:,n))                     !GYL ; kotsuki 20151216
          ELSE                                                                    !GYL ; kotsuki 20151216
            CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL ; kotsuki 20151216
                            trans(:,:,n),ltimer,                          &       !GYL ; kotsuki 20151216
                            transm=transm(:,n))                                   !GYL ; kotsuki 20151216
          END IF
          work3d(ij,ilev,n) = parm
          if(ilev<20) then
            if(nobsl > nobsl_max) nobsl_max=nobsl
            if(nobsl < nobsl_min) nobsl_min=nobsl
          end if
        END IF

        IF(RELAX_ALPHA /= 0.0d0) THEN                                                              !GYL - RTPP 
          CALL weight_RTPP(trans(:,:,n),transrlx)                                                  !GYL
        ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                                  !GYL - RTPS
          CALL weight_RTPS(trans(:,:,n),pa(:,:,n),gues3d(ij,ilev,:,n),transrlx,work3da(ij,ilev,n)) !GYL
        ELSE                                                                                       !GYL
          transrlx = trans(:,:,n)                                                                  !GYL
        END IF

        CALL CPU_TIME(ltimer00)                                                               ! KK
        DO m=1,nbv
          anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
          DO k=1,nbv
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                & + gues3d(ij,ilev,k,n) * (transrlx(k,m) + transm(k,n))                       !GYL
          END DO
        END DO

        IF(n==iv3d_qv) THEN
          IF(ilev >= 25) THEN ! no analysis for upper-level Q [Koji 2916.01.18]
            DO m = 1,nbv
              anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n) &
                                  + gues3d(ij,ilev,m,n)
            END DO
          END IF
        END IF
        IF(n==iv3d_qc) THEN
          DO m = 1,nbv
            anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n) &
                                + gues3d(ij,ilev,m,n)
          END DO
        END IF
!
        CALL CPU_TIME(rtimer)                                                                 ! KK
        ltimer(11) = ltimer(11) + rtimer-ltimer00                                             ! KK
        ltimer00 = rtimer                                                                     ! KK
!
      END DO
    END DO
  END DO

  WRITE(ADM_LOG_FID,*) '### END LETKF ###'
  FLUSH(30)

  DO m=1,nbv
    anal3d(:,nlev,m,:)=anal3d(:,nlev-1,m,:)
    anal3d(:,1,m,:)=anal3d(:,2,m,:)
  END DO
  CALL CPU_TIME(rtimer)                                                                       ! KK
  ltimer(11) = ltimer(11) + rtimer-ltimer00                                                   ! KK
  ltimer00 = rtimer                                                                           ! KK
!

  DEALLOCATE(hdxf,rdiag,rloc,dep)
  IF(cov_infl_mul < 0.0d0) THEN
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN ! [Change] Koji 2014.02.08
      inflfile=trim(infl_basedir)//'/'//trim(CDATE)//'/infl_mul'                              ! KK
      WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',inflfile
      FLUSH(ADM_LOG_FID)
      CALL write_icogrd_legacy(inflfile,work3dg,work2dg)
    END IF
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  END IF
  CALL CPU_TIME(rtimer)
  ltimer(12) = ltimer(12) + rtimer-ltimer00
  ltimer00 = rtimer
  !
  ! Additive inflation
  !
  IF(sp_infl_add > 0.0d0) THEN
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    DO n=1,nv3d
      DO m=1,nbv
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        END DO
      END DO
    END DO

    DEALLOCATE(work3d,work2d)
    WRITE(ADM_LOG_FID,'(A)') '===== Additive covariance inflation ====='
    WRITE(ADM_LOG_FID,'(A,F10.4)') '  parameter:',sp_infl_add
    WRITE(ADM_LOG_FID,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    DO n=1,nv3d
      DO m=1,nbv
        DO ilev=1,nlev
          DO ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_add
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
        END DO
      END DO
    END DO
  END IF

  DO m = 1, nbv
    DO ilev = 1, nlev
      DO ij = 1, nij1
        IF(anal3d(ij,ilev,m,iv3d_qv)<0.0d0) anal3d(ij,ilev,m,iv3d_qv)=0.0d0
        IF(anal3d(ij,ilev,m,iv3d_qc)<0.0d0) anal3d(ij,ilev,m,iv3d_qc)=0.0d0
      END DO
    END DO
  END DO

  ! Added 2017.05.07 Koji
  DO m = 1, nbv
    DO ij = 1, nij1
      IF(anal2d(ij,m,iv2d_q2m)<0.0d0) anal2d(ij,m,iv2d_q2m)=0.0d0
    END DO
  END DO
!
  CALL CPU_TIME(rtimer)
  ltimer(13) = ltimer(13) + rtimer-ltimer00
  ltimer00 = rtimer
!
  CALL CPU_TIME(rtimer)
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '##### TIMER(das_letkf_main):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

  WRITE(ADM_LOG_FID,*) 'gues'
  WRITE(ADM_LOG_FID,*) '3 dimensional variables'
  do n = 1, nv3d
    WRITE(ADM_LOG_FID,'(i5,2f15.5)') n, maxval(mean3d(:,:,n)), minval(mean3d(:,:,n))
  end do
  WRITE(ADM_LOG_FID,*) '2 dimensional variables'
  do n = 1, nv2d
    WRITE(ADM_LOG_FID,'(i5,2f15.5)') n, maxval(mean2d(:,n)), minval(mean2d(:,n))
  end do

  WRITE(ADM_LOG_FID,*) 'Analysis'
  WRITE(ADM_LOG_FID,*) '3 dimensional variables'
  do n = 1, nv3d
    WRITE(ADM_LOG_FID,'(i5,2f15.5)') n, maxval(anal3d(:,:,:,n)), minval(anal3d(:,:,:,n))
  end do
  WRITE(ADM_LOG_FID,*) '2 dimensional variables'
  do n = 1, nv2d
    WRITE(ADM_LOG_FID,'(i5,2f15.5)') n, maxval(anal2d(:,:,n)), minval(anal2d(:,:,n))
  end do

!
  call MPI_Gather(ltimer, slot, MPI_REAL8, ltimer_all, slot, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_01):', ltimer(1), ltimer(1)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_02):', ltimer(2), ltimer(2)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_03):', ltimer(3), ltimer(3)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_04):', ltimer(4), ltimer(4)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_05):', ltimer(5), ltimer(5)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_06):', ltimer(6), ltimer(6)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_07):', ltimer(7), ltimer(7)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_08):', ltimer(8), ltimer(8)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_09):', ltimer(9), ltimer(9)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_10):', ltimer(10), ltimer(10)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_11):', ltimer(11), ltimer(11)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_12):', ltimer(12), ltimer(12)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_13):', ltimer(13), ltimer(13)/sum(ltimer(:))*100
  WRITE(ADM_LOG_FID,'(A,2F10.2)') '####### TIMER(letkf_total):', sum(ltimer(:)), sum(ltimer(:))/sum(ltimer(:))*100
!
  WRITE(ADM_LOG_FID,'(A,I10)')    '####### nobsl_min=', nobsl_min
  WRITE(ADM_LOG_FID,'(A,I10)')    '####### nobsl_max=', nobsl_max
  FLUSH(30)
!  if(myrank==0) then
!    open(9, file='./ltime_log.bin', form='unformatted', access='direct', recl=slot*nprocs*8)
!    write(9,rec=1) ltimer_all
!    close(9)
!  end if
!
  DEALLOCATE(logpfm,mean3d,mean2d)
  RETURN
END SUBROUTINE das_letkf
!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local(ij,ilev,nvar,hdxf,rdiag,rloc,dep,nobsl,logpfm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,nvar
  REAL(r_size),INTENT(IN) :: logpfm(nij1,nlev)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER,ALLOCATABLE:: ntvs_use_prof(:),ntvs_use_inst(:),ntvs_use_slot(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,tvnn,iobs
!
! INITIALIZE
!
  IF( nobs > 0 ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF
  IF( ntvs > 0 ) THEN
    ALLOCATE(ntvs_use_prof(ntvs))
    ALLOCATE(ntvs_use_inst(ntvs))
    ALLOCATE(ntvs_use_slot(ntvs))
  END IF
!
! data search
!
  minlon = lon1(ij) - dlon_zero(ij)
  maxlon = lon1(ij) + dlon_zero(ij)
  minlat = lat1(ij) - dlat_zero
  maxlat = lat1(ij) + dlat_zero
  IF(maxlon - minlon >= 360.0d0) THEN
    minlon = 0.0d0
    maxlon = 360.0d0
  END IF

  DO jmin=1,nlat-2
    IF(minlat < lat(jmin+1)) EXIT
  END DO
  DO jmax=1,nlat-2
    IF(maxlat < lat(jmax+1)) EXIT
  END DO
  nn = 1
  tvnn = 1
  IF(minlon >= 0 .AND. maxlon <= 360.0) THEN
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    IF( nobs > 0 ) &
    & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    IF( ntvs > 0 ) &
    & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
    &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
  ELSE IF(minlon >= 0 .AND. maxlon > 360.0) THEN
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    maxlon = maxlon - 360.0d0
    IF(maxlon > 360.0d0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      IF( ntvs > 0 ) &
      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imax=1,nlon-1
        IF(maxlon < lon(imax+1)) EXIT
      END DO
      IF(imax > imin) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        imin = 1
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
        DO imin=1,nlon-1
          IF(minlon < lon(imin+1)) EXIT
        END DO
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  ELSE IF(minlon < 0 .AND. maxlon <= 360.0d0) THEN
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    minlon = minlon + 360.0d0
    IF(minlon < 0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      IF( ntvs > 0 ) &
      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imin=1,nlon-1
        IF(minlon < lon(imin+1)) EXIT
      END DO
      IF(imin < imax) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        imin = 1
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
        DO imin=1,nlon-1
          IF(minlon < lon(imin+1)) EXIT
        END DO
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  ELSE
    maxlon = maxlon - 360.0d0
    minlon = minlon + 360.0d0
    IF(maxlon > 360.0 .OR. minlon < 0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      IF( ntvs > 0 ) &
      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imin=1,nlon-1
        IF(minlon < lon(imin+1)) EXIT
      END DO
      DO imax=1,nlon-1
        IF(maxlon < lon(imax+1)) EXIT
      END DO
      IF(imin > imax) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        IF( ntvs > 0 ) &
        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  END IF
  nn = nn-1
  tvnn = tvnn -1
  IF( nn < 1 .AND. tvnn < 1 ) THEN
  !IF(nn < 1) THEN
    nobsl = 0
    RETURN
  END IF
!
! CONVENTIONAL
!
  nobsl = 0
  IF(nn > 0) THEN
    DO n=1,nn
      !
      ! vertical localization
      !
      IF(NINT(obselm(nobs_use(n))) == id_ps_obs .AND. ilev > 1) THEN
        dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))
      ELSE IF(NINT(obselm(nobs_use(n))) /= id_ps_obs) THEN
        dlev = ABS(LOG(obslev(nobs_use(n))) - logpfm(ij,ilev))
      ELSE
        dlev = 0.0d0
      END IF
      IF(dlev > dist_zerov) CYCLE
      !
      ! horizontal localization
      !
      tmplon=obslon(nobs_use(n))
      tmplat=obslat(nobs_use(n))
      CALL com_distll_1( tmplon, tmplat,lon1(ij), lat1(ij), dist)
      !write(*,'(6F15.4)') tmplon, tmplat, lon1(ij), lat1(ij), dist, dist_zero

      IF(dist > dist_zero ) CYCLE
      !
      ! variable localization
      !
      SELECT CASE(NINT(obselm(nobs_use(n))))
      CASE(id_u_obs)
        iobs=1
      CASE(id_v_obs)
        iobs=2
      CASE(id_t_obs)
        iobs=3
      CASE(id_qv_obs)
        iobs=4
      CASE(id_rh_obs)
        iobs=5
      CASE(id_ps_obs)
        iobs=6
      CASE(id_rain_obs)
        iobs=7
      CASE default ! added by sawada
        cycle
      END SELECT
      IF(var_local(nvar,iobs) < TINY(var_local)) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
      dep(nobsl)    = obsdep(nobs_use(n))
      !
      ! Observational localization
      !
      tmperr=obserr(nobs_use(n))
      rdiag(nobsl) = tmperr * tmperr
      rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &  ! KK
                  & * var_local(nvar,iobs)                                       ! KK
      !rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2) ) * var_local(nvar,iobs)   ! KK add(20160610)
    END DO
  END IF
!
! ATOVS
!
  !WRITE(ADM_LOG_FID,*) 'CHECK IN THE LETKF_TOOL'
  IF(tvnn > 0) THEN
    DO n=1,tvnn

      tmplon=tvslon(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
      tmplat=tvslat(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
      CALL com_distll_1( tmplon, tmplat, lon1(ij), lat1(ij), dist)
      IF( dist > dist_zero) CYCLE

      DO ichan=1,ntvsch(ntvs_use_inst(n))
        !
        ! vertical localization
        !
        dlev = &
        ABS(LOG(tvslev(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n)))-logpfm(ij,ilev))
        !write(ADM_LOG_FID,'(2i6,4f10.5)') n, ichan, &
        !  LOG(tvslev(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))),&
        !  logpfm(ij,ilev), dlev, dist_zerov
        IF(dlev > dist_zerov) CYCLE
        SELECT CASE(NINT(tvselm(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))))
        CASE(id_bt_obs)
          iobs=8
        END SELECT

        tmperr=tvserr(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
        tmpqc=tvsqc(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
        !tmpwgt(:)=tvswgt(:,ichan, &
        !                 & ntvs_use_prof(n), &
        !                 & ntvs_use_inst(n), &
        !                 & ntvs_use_slot(n))
        !IF( tmpqc == 1 .AND. tmpwgt(ilev) > 0.05D0 ) THEN
        IF( tmpqc == 1 ) THEN
          nobsl = nobsl + 1
          DO im = 1, nbv
            hdxf(nobsl,im) = tvshdxf(im,ichan, &
                              & ntvs_use_prof(n), &
                              & ntvs_use_inst(n), &
                              & ntvs_use_slot(n))
          END DO

          dep(nobsl)    = tvsdep(ichan, &
                              & ntvs_use_prof(n), &
                              & ntvs_use_inst(n), &
                              & ntvs_use_slot(n))
          rdiag(nobsl)  = tmperr * tmperr
          !rloc(nobsl)   = exp(-0.5d0 * (dist/sigma_obs)**2) &
          !              & * (tmpwgt(ilev) * tmpwgt(ilev))
          !write(ADM_LOG_FID,'(3i10,2f20.5,2f12.5)') n,nobstotal, nobsl, dist, dlev, sigma_obs, sigma_obsv
          !flush(ADM_LOG_FID)
          rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &  ! KK
                    & * var_local(nvar,iobs)                                         ! KK
          !rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2) ) * var_local(nvar,iobs)   ! KK add(20160610)
          
        END IF
      END DO
    END DO
  END IF
!
! DEBUG
! IF( ILEV == 1 .AND. ILON == 1 ) &
! & WRITE(ADM_LOG_FID,*) 'ILEV,ILON,ILAT,NN,TVNN,NOBSL=',ilev,ij,nn,tvnn,nobsl
!
  IF( nobsl > nobstotal ) THEN
    WRITE(ADM_LOG_FID,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    !WRITE(ADM_LOG_FID,*) 'IJ,NN,TVNN=', ij, nn, tvnn
    WRITE(ADM_LOG_FID,*) 'IJ,NN=', ij, nn
    STOP 99
  END IF
!
  IF( nobs > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF
  IF( ntvs > 0 ) THEN
    DEALLOCATE(ntvs_use_prof)
    DEALLOCATE(ntvs_use_inst)
    DEALLOCATE(ntvs_use_slot)
  END IF
!
  RETURN
END SUBROUTINE obs_local

SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  DO j=jmin,jmax
    IF(imin > 1) THEN
      ib = nobsgrd(imin-1,j)+1
    ELSE
      IF(j > 1) THEN
        ib = nobsgrd(nlon,j-1)+1
      ELSE
        ib = 1
      END IF
    END IF
    ie = nobsgrd(imax,j)
    n = ie - ib + 1
    IF(n == 0) CYCLE
    DO ip=ib,ie
      IF(nn > nobs) THEN
        WRITE(ADM_LOG_FID,*) 'FATALERROR, NN > NOBS', NN, NOBS
      END IF
      nobs_use(nn) = ip
      nn = nn + 1
    END DO
  END DO

  RETURN
END SUBROUTINE obs_local_sub

SUBROUTINE tvs_local_sub(imin,imax,jmin,jmax,nn,ntvs_prof,ntvs_inst,ntvs_slot)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, ntvs_prof(ntvs), ntvs_inst(ntvs), ntvs_slot(ntvs)
  INTEGER :: j,n,ib,ie,ip
  INTEGER :: islot, iinst

  DO j=jmin,jmax
    DO islot=1,nslots
      DO iinst=1,ninstrument
        IF(imin > 1) THEN
          ib = ntvsgrd(imin-1,j,iinst,islot)+1
        ELSE
          IF(j > 1) THEN
            ib = ntvsgrd(nlon,j-1,iinst,islot)+1
          ELSE
            ib = 1
          END IF
        END IF
        ie = ntvsgrd(imax,j,iinst,islot)
        n = ie - ib + 1
        IF(n == 0) CYCLE
        DO ip=ib,ie
          IF(nn > ntvs) THEN
          !IF(nn > nobs) THEN
            !WRITE(ADM_LOG_FID,*) 'FATALERROR, NN > NTVS', NN, NTVS
          END IF
          ntvs_prof(nn)=ip
          ntvs_inst(nn)=iinst
          ntvs_slot(nn)=islot
          nn = nn + 1
        END DO
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE tvs_local_sub
!TVS!-----------------------------------------------------------------------
!TVS! Data Assimilation for VARBC
!TVS!-----------------------------------------------------------------------
!TVSSUBROUTINE das_vbc(um,vm,tm,qm,qlm,psm,vbcf,vbca)
!TVS  USE common_mtx
!TVS  IMPLICIT NONE
!TVS  REAL(r_size),INTENT(IN) :: um(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: vm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: tm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: qm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: qlm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: psm(nij1)
!TVS  REAL(r_size),INTENT(INOUT) :: vbcf(maxvbc,maxtvsch,ninstrument)
!TVS  REAL(r_size),INTENT(OUT)   :: vbca(maxvbc,maxtvsch,ninstrument)
!TVS  REAL(r_sngl) :: u4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: v4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: t4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: q4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: ql4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: ps4(nlon,nlat)
!TVS  REAL(r_size) :: u(nlon,nlat,nlev)
!TVS  REAL(r_size) :: v(nlon,nlat,nlev)
!TVS  REAL(r_size) :: t(nlon,nlat,nlev)
!TVS  REAL(r_size) :: q(nlon,nlat,nlev)
!TVS  REAL(r_size) :: ql(nlon,nlat,nlev)
!TVS  REAL(r_size) :: ps(nlon,nlat)
!TVS  REAL(r_size) :: p_full(nlon,nlat,nlev)
!TVS  REAL(r_size),ALLOCATABLE :: hx(:,:,:,:)
!TVS  REAL(r_size),ALLOCATABLE :: pred(:,:,:,:,:)
!TVS  INTEGER,ALLOCATABLE :: tmpqc(:,:,:)
!TVS  REAL(r_size),ALLOCATABLE :: tmpwgt(:,:,:,:)
!TVS  REAL(r_size) :: a(maxvbc,maxvbc)
!TVS  REAL(r_size) :: b(maxvbc)
!TVS  REAL(r_size) :: ainv(maxvbc,maxvbc)
!TVS  INTEGER:: ntvschan1(maxtvsch,ninstrument)
!TVS  INTEGER:: i,j,k,n,islot,nn
!TVS
!TVS  PRINT *,'Hello from das_vbc'
!TVS
!TVS  IF(ntvs == 0) THEN
!TVS    PRINT *,'No radiance data: das_vbc skipped..'
!TVS!$OMP PARALLEL WORKSHARE
!TVS    vbca = vbcf
!TVS!$OMP END PARALLEL WORKSHARE
!TVS    RETURN
!TVS  END IF
!TVS
!TVS  CALL gather_grd_mpi(0,um,vm,tm,qm,qlm,psm,u4,v4,t4,q4,ql4,ps4)
!TVS  n = nlon*nlat*nlev
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(u4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(v4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(t4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(q4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(ql4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  n = nlon*nlat
!TVS  CALL MPI_BCAST(ps4(1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS!$OMP PARALLEL WORKSHARE
!TVS  u = REAL(u4,r_size)
!TVS  v = REAL(v4,r_size)
!TVS  t = REAL(t4,r_size)
!TVS  q = REAL(q4,r_size)
!TVS  ql = REAL(ql4,r_size)
!TVS  ps = REAL(ps4,r_size)
!TVS!$OMP END PARALLEL WORKSHARE
!TVS  CALL calc_pfull(ps,p_full)
!TVS
!TVS  ALLOCATE( hx(maxtvsch,maxtvsprof,ninstrument,nslots) )
!TVS  ALLOCATE( pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots) )
!TVS  ALLOCATE( tmpqc(maxtvsch,maxtvsprof,ninstrument) )
!TVS  ALLOCATE( tmpwgt(nlev,maxtvsch,maxtvsprof,ninstrument) )
!TVS  DO islot=1,nslots
!TVS!    IF(SUM(ntvsprofslots(:,islot)) == 0) CYCLE
!TVS    ntvsprof(:) = ntvsprofslots(:,islot)
!TVS    CALL Trans_XtoY_tvs(u,v,t,q,ql,ps,p_full, &
!TVS      & tvslon(:,:,islot),tvslat(:,:,islot),tvszenith(:,:,islot),&
!TVS      & tvsskin(:,:,islot),tvsstmp(:,:,islot),tvsclw(:,:,islot),&
!TVS      & tvsemis(:,:,:,islot),tmpqc,hx(:,:,:,islot),tmpwgt,pred(:,:,:,:,islot))
!TVS  END DO
!TVS  DEALLOCATE(tmpqc,tmpwgt)
!TVS
!TVS!$OMP PARALLEL PRIVATE(j,k,n,a,b,ainv)
!TVS!$OMP WORKSHARE
!TVS  vbca = 0.0d0
!TVS!$OMP END WORKSHARE
!TVS!$OMP DO SCHEDULE(DYNAMIC)
!TVS  DO k=1,ninstrument
!TVS    DO j=1,maxtvsch
!TVS      !
!TVS      ! Parallel processing
!TVS      !
!TVS      IF(MOD(j+maxtvsch*(k-1)-1,nprocs) /= myrank) CYCLE
!TVS      !
!TVS      ! DATA NUMBER
!TVS      !
!TVS      ntvschan(j,k) = SUM(tvsqc(j,:,k,:))
!TVS      IF(msw_vbc .AND. ntvschan(j,k) /= 0 ) THEN
!TVS        PRINT '(3A,I3,A,I6)',' >> VBC executed for instrument,channel,ntvsl: ',&
!TVS                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
!TVS        CALL vbc_local(j,k,ntvschan(j,k),hx,pred,a,b)
!TVS        CALL mtx_inv(maxvbc,a,ainv)
!TVS        vbca(:,j,k) = vbcf(:,j,k)
!TVS        DO n=1,maxvbc
!TVS          vbca(:,j,k) = vbca(:,j,k) - ainv(:,n)*b(n) !ATTN: sign for beta
!TVS        END DO
!TVS      ELSE
!TVS        PRINT '(3A,I3,A,I6)',' !! NO VBC executed for instrument,channel,ntvsl: ',&
!TVS                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
!TVS        vbca(:,j,k) = vbcf(:,j,k)
!TVS      END IF
!TVS    END DO
!TVS  END DO
!TVS!$OMP END DO
!TVS!$OMP WORKSHARE
!TVS  vbcf = vbca
!TVS  ntvschan1 = ntvschan
!TVS!$OMP END WORKSHARE
!TVS!$OMP END PARALLEL
!TVS  DEALLOCATE(hx,pred)
!TVS  n = maxvbc*maxtvsch*ninstrument
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!TVS  CALL MPI_ALLREDUCE(vbcf,vbca,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,j)
!TVS  n = maxtvsch*ninstrument
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!TVS  CALL MPI_ALLREDUCE(ntvschan1,ntvschan,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,j)
!TVS
!TVS  RETURN
!TVSEND SUBROUTINE das_vbc
!TVS!-----------------------------------------------------------------------
!TVS!  (in) ichan: channnel
!TVS!  (in) iinst: sensor
!TVS!  (out) a = B_beta^-1 + p R^-1 p^T
!TVS!  (out) b = p R^-1 d
!TVS!-----------------------------------------------------------------------
!TVSSUBROUTINE vbc_local(ichan,iinst,ntvsl,hx,pred,a,b)
!TVS  IMPLICIT NONE
!TVS  INTEGER,PARAMETER :: msw=1
!TVS  INTEGER,PARAMETER :: nmin=400
!TVS  INTEGER,INTENT(IN) :: ichan,iinst,ntvsl
!TVS  REAL(r_size),INTENT(IN) :: hx(maxtvsch,maxtvsprof,ninstrument,nslots)
!TVS  REAL(r_size),INTENT(IN) :: pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots)
!TVS  REAL(r_size),INTENT(OUT) :: a(maxvbc,maxvbc)
!TVS  REAL(r_size),INTENT(OUT) :: b(maxvbc)
!TVS  REAL(r_size) :: dep,dep0
!TVS  REAL(r_size) :: bias,bias0
!TVS  REAL(r_size) :: r,tmp
!TVS  INTEGER:: islot, iprof, i,j,n
!TVS
!TVS  a = 0.0d0
!TVS  b = 0.0d0
!TVS  dep = 0.0d0
!TVS  dep0 = 0.0d0
!TVS  bias = 0.0d0
!TVS  bias0 = 0.0d0
!TVS  n = 0
!TVS  DO islot=1,nslots
!TVS    DO iprof=1,maxtvsprof
!TVS      IF(tvsqc(ichan,iprof,iinst,islot)/=1) CYCLE
!TVS      !
!TVS      ! R
!TVS      !
!TVS      r = tvserr(ichan,iprof,iinst,islot)**2
!TVS      !
!TVS      ! p R^-1 p^T
!TVS      !
!TVS      DO j=1,maxvbc
!TVS        DO i=1,maxvbc
!TVS          a(i,j) = a(i,j) &
!TVS               & + pred(i,ichan,iprof,iinst,islot) &
!TVS               & * pred(j,ichan,iprof,iinst,islot) / r
!TVS        END DO
!TVS      END DO
!TVS      !
!TVS      ! B_beta^-1
!TVS      !
!TVS      IF(msw == 1) THEN ! Y.Sato
!TVS        IF(ntvsl < nmin) THEN
!TVS          tmp = REAL(nmin,r_size) / r
!TVS
!TVS        ELSE
!TVS          tmp = (REAL(ntvsl,r_size) &
!TVS            & / (LOG10(REAL(ntvsl,r_size)/REAL(nmin,r_size))+1.0d0)) / r
!TVS        END IF
!TVS      ELSE IF(msw == 2) THEN ! D.Dee
!TVS        tmp = REAL(ntvsl,r_size) / r
!TVS      ELSE ! Constant
!TVS        tmp = 100.0d0
!TVS      END IF
!TVS      DO i=1,maxvbc
!TVS        a(i,i) = a(i,i) + tmp
!TVS      END DO
!TVS      !
!TVS      ! p R^-1 d
!TVS      !
!TVS      b(:) = b(:) + pred(:,ichan,iprof,iinst,islot) / r &
!TVS                & *(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))
!TVS      bias = bias+tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot)
!TVS      dep = dep+(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))**2
!TVS      bias0= bias0+tvsdep(ichan,iprof,iinst,islot)
!TVS      dep0= dep0+tvsdep(ichan,iprof,iinst,islot)**2
!TVS      n = n+1
!TVS    END DO
!TVS  END DO
!TVS
!TVS  dep = SQRT(dep / REAL(n,r_size))
!TVS  dep0 = SQRT(dep0 / REAL(n,r_size))
!TVS  bias = bias / REAL(n,r_size)
!TVS  bias0 = bias0 / REAL(n,r_size)
!TVS  PRINT '(2A,I3,4F12.4)',' >> D monit: ',tvsname(iinst),tvsch(ichan,iinst),bias0,bias,dep0,dep
!TVS
!TVS  RETURN
!TVSEND SUBROUTINE vbc_local
!-----------------------------------------------------------------------
! Relaxation via LETKF weight - RTPP method
!-----------------------------------------------------------------------
subroutine weight_RTPP(w, wrlx)
  implicit none
  real(r_size), intent(in) :: w(nbv,nbv)
  real(r_size), intent(out) :: wrlx(nbv,nbv)
  integer :: m

  wrlx = (1.0d0 - RELAX_ALPHA) * w
  do m = 1, nbv
    wrlx(m,m) = wrlx(m,m) + RELAX_ALPHA
  end do

  return
end subroutine weight_RTPP
!-----------------------------------------------------------------------
! Relaxation via LETKF weight - RTPS method
!-----------------------------------------------------------------------
subroutine weight_RTPS(w, pa, xb, wrlx, infl)
  implicit none
  real(r_size), intent(in) :: w(nbv,nbv)
  real(r_size), intent(in) :: pa(nbv,nbv)
  real(r_size), intent(in) :: xb(nbv)
  real(r_size), intent(out) :: wrlx(nbv,nbv)
  real(r_size), intent(out) :: infl
  real(r_size) :: var_g, var_a
  integer :: m, k

  var_g = 0.0d0
  var_a = 0.0d0
  do m = 1, nbv
    var_g = var_g + xb(m) * xb(m)
    do k = 1, nbv
      var_a = var_a + xb(k) * pa(k,m) * xb(m)
    end do
  end do
  if (var_g > 0.0d0 .and. var_a > 0.0d0) then
    infl = RELAX_ALPHA_SPREAD * sqrt(var_g / (var_a * real(nbv-1,r_size))) - RELAX_ALPHA_SPREAD + 1.0d0   ! Whitaker and Hamill 2012
!    infl = sqrt(RELAX_ALPHA_SPREAD * (var_g / (var_a * real(nbv-1,r_size))) - RELAX_ALPHA_SPREAD + 1.0d0) ! Hamrud et al. 2015 (slightly modified)
    wrlx = w * infl
  else
    wrlx = w
    infl = 1.0d0
  end if

  return
end subroutine weight_RTPS

END MODULE letkf_tools
