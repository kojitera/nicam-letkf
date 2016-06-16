MODULE common_letkf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mtx

  IMPLICIT NONE

  PUBLIC
!=======================================================================
!  LEKF Model Independent Parameters
!=======================================================================
!!  INTEGER,PARAMETER :: nbv=100    ! ensemble size

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     nbv              : ensemble size
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,nbv)   : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!   OUTPUT
!     trans(nbv,nbv) : transformation matrix
!=======================================================================
SUBROUTINE letkf_core(nbv,nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,trans,ltimer,transm,pao) ! KK
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nbv
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  REAL(r_size),INTENT(IN) :: hdxb(1:nobs,1:nbv)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobs)
  REAL(r_size),INTENT(IN) :: rloc(1:nobs)
  REAL(r_size),INTENT(IN) :: dep(1:nobs)
  REAL(r_size),INTENT(INOUT) :: parm_infl
  REAL(r_size),INTENT(OUT) :: trans(nbv,nbv)
  REAL(r_size) :: hdxb_rinv(nobsl,nbv)
  REAL(r_size) :: eivec(nbv,nbv)
  REAL(r_size) :: eival(nbv)
  REAL(r_size) :: pa(nbv,nbv)
  REAL(r_size) :: work1(nbv,nbv)
  REAL(r_size) :: work2(nbv,nobsl)
  REAL(r_size) :: work3(nbv)
  REAL(r_size) :: rho
  REAL(r_size) :: parm(4),sigma_o,gain
  REAL(r_size),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl
  REAL(r_size) :: rtimer00, rtimer                                           ! KK
  REAL(r_size),INTENT(INOUT) :: ltimer(20)                                   ! KK
  INTEGER :: i,j,k
  !! EPES (Ruiz et al. 2013); kotsuki 20151203
  LOGICAL, PARAMETER :: EPES = .false.  ! assumes always false            
  REAL(r_size) :: epeslamda
  REAL(r_size) :: wvar(nbv), pachk(nbv,nbv),transchk(nbv,nbv)
  !!RTPP & RTPS!! kotsuki 20151216
  REAL(r_size),INTENT(OUT),OPTIONAL :: transm(nbv)
  REAL(r_size),INTENT(OUT),OPTIONAL :: pao(nbv,nbv)
!
  CALL CPU_TIME(rtimer00)                                                    ! KK
!
  IF(nobsl == 0) THEN
    trans = 0.0d0
    DO i=1,nbv
      trans(i,i) = SQRT(parm_infl)
    END DO
    IF (PRESENT(transm)) THEN                     !GYL ; kotsuki 20151216
      transm = 0.0d0                              !GYL ; kotsuki 20151216
    END IF                                        !GYL ; kotsuki 20151216
    IF (PRESENT(pao)) THEN                        !GYL ; kotsuki 20151216
      pao = 0.0d0                                 !GYL ; kotsuki 20151216
      DO i=1,nbv                                  !GYL ; kotsuki 20151216
        pao(i,i) = parm_infl / REAL(nbv-1,r_size) !GYL ; kotsuki 20151216
      END DO                                      !GYL ; kotsuki 20151216
    END IF                                        !GYL ; kotsuki 20151216
    RETURN
  ELSE
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
    END DO
  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(1) = ltimer(1) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  CALL dgemm('t','n',nbv,nbv,nobsl,1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),&
    & nobsl,0.0d0,work1,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
!      END DO
!    END DO
!  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(2) = ltimer(2) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  rho = 1.0d0 / parm_infl
  DO i=1,nbv
    work1(i,i) = work1(i,i) + REAL(nbv-1,r_size) * rho
  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(3) = ltimer(3) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(4) = ltimer(4) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,&
    & nbv,0.0d0,pa,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      pa(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(5) = ltimer(5) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
  CALL dgemm('n','t',nbv,nobsl,nbv,1.0d0,pa,nbv,hdxb_rinv,&
    & nobsl,0.0d0,work2,nbv)
!  DO j=1,nobsl
!    DO i=1,nbv
!      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
!      DO k=2,nbv
!        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
!      END DO
!    END DO
!  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(6) = ltimer(6) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,nbv
    work3(i) = work2(i,1) * dep(1)
    DO j=2,nobsl
      work3(i) = work3(i) + work2(i,j) * dep(j)
    END DO
  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(7) = ltimer(7) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
  DO j=1,nbv
    rho = SQRT( REAL(nbv-1,r_size) / eival(j) )
    DO i=1,nbv
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,&
    & nbv,0.0d0,trans,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      trans(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(8) = ltimer(8) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!
!-----------------------------------------------------------------------
!  EPES :: kotsuki 20151203
!-----------------------------------------------------------------------
  IF ( EPES ) THEN
    epeslamda=0.0d0
    !!
    !!DO j=1,nbv
    !!  epeslamda = epeslamda + pa(j,j)
    !!END DO
    !!epeslamda = dsqrt( REAL(nbv,r_size) / REAL(nbv-1,r_size) / epeslamda )

    !! following code of Juan
    DO j=1,nbv
      CALL com_covar( nbv, trans(:,j), trans(:,j), wvar(j) )
    END DO
    epeslamda = dsqrt( REAL(nbv,r_size) / REAL(nbv-1,r_size) / sum( wvar(:) ) )

    DO j=1,nbv
      DO i=1,nbv-1
        trans(i,j) = trans(i,j) * epeslamda
      END DO
    END DO
  END IF
!                                                                            !
!                                                                            KK
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  IF (PRESENT(transm)) THEN                !GYL - if transm is present,  
    transm = work3                         !GYL - return both trans and transm without adding them
  ELSE                                     !GYL
    DO j=1,nbv
      DO i=1,nbv
        trans(i,j) = trans(i,j) + work3(i)
      END DO
    END DO
  END IF                                   !GYL
  IF (PRESENT(pao)) pao = pa               !GYL
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(9) = ltimer(9) + rtimer-rtimer00                                    ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
!-----------------------------------------------------------------------
!  Inflation estimation
!-----------------------------------------------------------------------
  parm = 0.0d0
  DO i=1,nobsl
    parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
  END DO
  DO j=1,nbv
    DO i=1,nobsl
      parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
    END DO
  END DO
  parm(2) = parm(2) / REAL(nbv-1,r_size)
  parm(3) = SUM(rloc(1:nobsl))
  parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
!  sigma_o = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rloc(1:nobsl))
  sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
  gain = sigma_b**2 / (sigma_o + sigma_b**2)
  parm_infl = parm_infl + gain * parm(4)
!                                                                            ! KK
  CALL CPU_TIME(rtimer)                                                      ! KK
  ltimer(10) = ltimer(10) + rtimer-rtimer00                                  ! KK
  rtimer00=rtimer                                                            ! KK
!                                                                            ! KK
  RETURN
  END IF
END SUBROUTINE letkf_core

END MODULE common_letkf
