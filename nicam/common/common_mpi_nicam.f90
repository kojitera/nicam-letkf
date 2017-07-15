MODULE common_mpi_nicam
!=======================================================================
!
! [PURPOSE:] MPI procedures
!
! [ATTENTION:]
!   DO NOT COMPILE WITH BOTH INLINE EXPANSION AND OMP OPTIONS TOGETHER
!    (Use ONE if you want, but DON'T USE BOTH AT THE SAME TIME)
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE mod_adm
  USE common
  USE common_mpi
  USE common_nicam
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: mpibufsize=1000
  !INTEGER,PARAMETER :: mpibufsize=10240
  INTEGER,SAVE :: nij1
  INTEGER,SAVE :: nij1max
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: dx1(:),dy1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),SAVE :: time_IO(10)

CONTAINS
SUBROUTINE set_common_mpi_nicam
  REAL(r_size) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  !REAL(r_size) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  !REAL(r_size) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,n

  time_IO(:)=0.0d0
  WRITE(ADM_LOG_FID,'(A)') 'Hello from set_common_mpi_nicam'
  i = MOD(ADM_gall*ADM_rgn_nmax,nprocs)
  nij1max = (ADM_gall*ADM_rgn_nmax - i)/nprocs + 1
  IF(myrank < i) THEN
    nij1 = nij1max
  ELSE
    nij1 = nij1max - 1
  END IF
  WRITE(ADM_LOG_FID,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs))
  nij1node(:)=nij1max
  DO n=1,nprocs
    IF(n-1 < i) THEN
      nij1node(n) = nij1max
    ELSE
      nij1node(n) = nij1max - 1
    END IF
  END DO

  ALLOCATE(phi1(nij1))
  ALLOCATE(dx1(nij1))
  ALLOCATE(dy1(nij1))
  ALLOCATE(lon1(nij1))
  ALLOCATE(lat1(nij1))

  ALLOCATE(v3d(nij1,nlev,nv3d))
  ALLOCATE(v2d(nij1,nv2d))
  ! --- initialization ---
  v3d = 0.0d0
  v2d = 0.0d0
  v3dg = 0.0d0
  v2dg = 0.0d0

  DO i=1,ADM_rgn_nmax
    v3dg(:,i,1,3) = ico_lon(:,i)
    v3dg(:,i,1,4) = ico_lat(:,i)
  END DO
  CALL scatter_grd_mpi(0,v3dg,v2dg,v3d,v2d)
  lon1(:) = v3d(:,1,3)
  lat1(:) = v3d(:,1,4)

  RETURN
END SUBROUTINE set_common_mpi_nicam
!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_alltoall(member,mstart,mend,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: member, mstart, mend
  REAL(r_size),INTENT(IN) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2d(nij1,member,nv2d)
  REAL(r_size) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_size) :: bufs_tmp(nij1max,nprocs)
  REAL(r_size),ALLOCATABLE :: bufr(:,:,:)
  INTEGER :: j,k,m,n,ierr,ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs) ! KK
  INTEGER :: mcount                                                     ! KK

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) THEN
    print *, "CHECK mcount! mcount > nprocs or <= 0. mcount =", mcount
    STOP
  END IF
  ALLOCATE(bufr(nij1max,nlevall,mcount))

  DO n = 1, nv3d
    !j=0
    IF(myrank < mcount) THEN
      !DO n = 1, nv3d
        DO k = 1, nlev
          !j=j+1
          CALL grd_to_buf(v3dg(:,:,k,n), bufs(:,k,:))
          !CALL grd_to_buf(v3dg(:,:,k,n), bufs_tmp(:,:))
          !bufs(:,j,:)=bufs_tmp(:,:)
        END DO
      !END DO
      !DO n = 1, nv2d
      !  j=j+1
      !  !CALL grd_to_buf(v2dg(:,:,n), bufs(:,j,:))
      !  CALL grd_to_buf(v2dg(:,:,n), bufs_tmp(:,:))
      !  bufs(:,j,:)=bufs_tmp(:,:)
      !END DO
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_REAL8, &
                      & bufr, nij1max*nlevall, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ELSE
      CALL set_alltoallv_counts(mcount, nij1max*nlevall, nr, nrt, ns, nst)
      CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL8, &
                      &  bufr, nr, nrt, MPI_REAL8, MPI_COMM_WORLD, ierr)
    END IF

    DO m = mstart, mend
      DO k = 1, nlev
        v3d(:,k,m,n) = REAL(bufr(1:nij1,k,m-mstart+1), r_size)
      END DO
    END DO
  END DO

  DO n = 1, nv2d
    IF(myrank < mcount) THEN
      !DO k = 1, nlev
      CALL grd_to_buf(v2dg(:,:,n), bufs(:,1,:))
      !END DO
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_REAL8, &
                      & bufr, nij1max*nlevall, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ELSE
      CALL set_alltoallv_counts(mcount, nij1max*nlevall, nr, nrt, ns, nst)
      CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL8, &
                      &  bufr, nr, nrt, MPI_REAL8, MPI_COMM_WORLD, ierr)
    END IF

    DO m = mstart, mend
      v2d(:,m,n) = REAL(bufr(1:nij1,1,m-mstart+1), r_size)
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr)

  RETURN
END SUBROUTINE scatter_grd_mpi_alltoall

SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)

  IF(mpibufsize > nij1max) THEN
    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  ELSE
    CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
  END IF

  RETURN
END SUBROUTINE scatter_grd_mpi

SUBROUTINE scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_size) :: tmp(nij1max,nprocs)
  REAL(r_size) :: bufs(mpibufsize,nprocs)
  REAL(r_size) :: bufr(mpibufsize)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ns = mpibufsize
  nr = ns
  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

  DO n=1,nv3d
    DO k=1,nlev
      IF(myrank == nrank) CALL grd_to_buf(v3dg(:,:,k,n),tmp)
      DO iter=1,niter
        IF(myrank == nrank) THEN
          i = mpibufsize * (iter-1)
          DO j=1,mpibufsize
            i=i+1
            IF(i > nij1max) EXIT
            bufs(j,:) = tmp(i,:)
          END DO
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_SCATTER(bufs,ns,MPI_REAL8,&
                       & bufr,nr,MPI_REAL8,nrank,MPI_COMM_WORLD,ierr)
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          v3d(i,k,n) = REAL(bufr(j),r_size)
        END DO
      END DO
    END DO
  END DO

  DO n=1,nv2d
    IF(myrank == nrank) CALL grd_to_buf(v2dg(:,:,n),tmp)
    DO iter=1,niter
      IF(myrank == nrank) THEN
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1max) EXIT
          bufs(j,:) = tmp(i,:)
        END DO
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_SCATTER(bufs,ns,MPI_REAL8,&
                     & bufr,nr,MPI_REAL8,nrank,MPI_COMM_WORLD,ierr)
      i = mpibufsize * (iter-1)
      DO j=1,mpibufsize
        i=i+1
        IF(i > nij1) EXIT
        v2d(i,n) = REAL(bufr(j),r_size)
      END DO
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi_safe

SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_size) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_size) :: bufr(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL8,&
                 & bufr,nr,MPI_REAL8,nrank,MPI_COMM_WORLD,ierr)

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi_fast
!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi_alltoall(member,mstart,mend,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: member, mstart, mend
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size),ALLOCATABLE :: bufs(:,:,:)
  REAL(r_size) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: k, m, n, ierr
  INTEGER :: ns(nprocs), nst(nprocs), nr(nprocs), nrt(nprocs) ! KK
  INTEGER :: mcount                                           ! KK

! ns = nij1max * nlevall

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) THEN
    print *, "CHECK mcount! mcount > nprocs or <= 0. mcount =", mcount
    STOP
  END IF
  ALLOCATE(bufs(nij1max,nlevall,mcount))

  DO n = 1, nv3d
    DO m = mstart, mend
      DO k = 1, nlev
        bufs(1:nij1,k,m-mstart+1) = REAL(v3d(:,k,m,n),r_size)
      END DO
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_REAL8, &
                      & bufr, nij1max*nlevall, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ELSE
      CALL set_alltoallv_counts(mcount, nij1max*nlevall, ns, nst, nr, nrt)
      CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL8, &
                       & bufr, nr, nrt, MPI_REAL8, MPI_COMM_WORLD, ierr)
    END IF

    IF(myrank < mcount) THEN
      DO k = 1, nlev
        CALL buf_to_grd(bufr(:,k,:),v3dg(:,:,k,n))
      END DO
    END IF

  END DO

  DO n = 1, nv2d
    DO m = mstart, mend
        bufs(1:nij1,1,m-mstart+1) = REAL(v2d(:,m,n),r_size)
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_REAL8, &
                      & bufr, nij1max*nlevall, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ELSE
      CALL set_alltoallv_counts(mcount, nij1max*nlevall, ns, nst, nr, nrt)
      CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL8, &
                       & bufr, nr, nrt, MPI_REAL8, MPI_COMM_WORLD, ierr)
    END IF

    IF(myrank < mcount) THEN
      CALL buf_to_grd(bufr(:,1,:),v2dg(:,:,n))
    END IF
  END DO


  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufs)

  RETURN
END SUBROUTINE gather_grd_mpi_alltoall

SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)

  IF(mpibufsize > nij1max) THEN
    CALL gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  ELSE
    CALL gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
  END IF

  RETURN
END SUBROUTINE gather_grd_mpi

SUBROUTINE gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size) :: tmp(nij1max,nprocs)
  REAL(r_size) :: bufs(mpibufsize)
  REAL(r_size) :: bufr(mpibufsize,nprocs)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ns = mpibufsize
  nr = ns
  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

  DO n=1,nv3d
    DO k=1,nlev
      DO iter=1,niter
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          bufs(j) = REAL(v3d(i,k,n),r_size)
        END DO
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_GATHER(bufs,ns,MPI_REAL8,&
                      & bufr,nr,MPI_REAL8,nrank,MPI_COMM_WORLD,ierr)
        IF(myrank == nrank) THEN
          i = mpibufsize * (iter-1)
          DO j=1,mpibufsize
            i=i+1
            IF(i > nij1max) EXIT
            tmp(i,:) = bufr(j,:)
          END DO
        END IF
      END DO
      IF(myrank == nrank) CALL buf_to_grd(tmp,v3dg(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    DO iter=1,niter
      i = mpibufsize * (iter-1)
      DO j=1,mpibufsize
        i=i+1
        IF(i > nij1) EXIT
        bufs(j) = REAL(v2d(i,n),r_size)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_GATHER(bufs,ns,MPI_REAL8,&
                    & bufr,nr,MPI_REAL8,nrank,MPI_COMM_WORLD,ierr)
      IF(myrank == nrank) THEN
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1max) EXIT
          tmp(i,:) = bufr(j,:)
        END DO
      END IF
    END DO
    IF(myrank == nrank) CALL buf_to_grd(tmp,v2dg(:,:,n))
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi_safe

SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  REAL(r_size) :: bufs(nij1max,nlevall)
  REAL(r_size) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL8,&
                & bufr,nr,MPI_REAL8,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi_fast
!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_ens_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_size) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(256) :: filename1=''
  CHARACTER(256) :: filename2=''
  CHARACTER(6) :: cmem
  REAL(r_size) :: tmp_time(3)
  INTEGER :: mstart, mend

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    tmp_time(1)=MPI_WTIME()
    !im = myrank+1 + (l-1)*nprocs
    im = myrank + (l-1)*nprocs ! 2017.05.07 Koji
    IF(im <= member-1) THEN
      WRITE(cmem,'(I6.6)') im
      filename1=trim(gues_basedir)//'/'//CDATE//'/mem'//cmem//'/restart_da'  ! KK
      filename2=trim(gues_basedir)//'/'//CDATE//'/mem'//cmem//'/history'
      WRITE(ADM_LOG_FID,'(A,I3.3,2A)') &
            'MYRANK ',myrank,' is reading a file ',trim(filename1)
      CALL read_icogrd(filename1,filename2,v3dg,v2dg)
    END IF
    tmp_time(2)=MPI_WTIME()

    time_IO(1)=time_IO(1)+tmp_time(2)-tmp_time(1)

!   DO n=0,nprocs-1                                                          ! KK
!     im = n+1 + (l-1)*nprocs                                                ! KK
!     IF(im <= member) THEN                                                  ! KK
!       CALL scatter_grd_mpi_alltoall(n,v3dg,v2dg,v3d(:,:,im,:),v2d(:,im,:)) ! KK
!     END IF                                                                 ! KK
!   END DO                                                                   ! KK
    mstart = 1 + (l-1)*nprocs                                                ! KK (add 20141017)
    mend   = MIN(l*nprocs, member)                                           ! KK (add 20141017)
    CALL scatter_grd_mpi_alltoall(member,mstart,mend,v3dg,v2dg,v3d,v2d)      ! KK (add 20141017)
    tmp_time(3)=MPI_WTIME()
    time_IO(2)=time_IO(2)+tmp_time(3)-tmp_time(2)
  END DO

  RETURN
END SUBROUTINE read_ens_mpi
!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: l,n,ll,im, mstart, mend
  CHARACTER(256) :: filename1
  CHARACTER(256) :: filename2=''
  CHARACTER(6) :: cmem
  REAL(r_size) :: tmp_time(3)

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    tmp_time(1)=MPI_WTIME()
!   DO n=0,nprocs-1                                                    ! KK
!     im = n+1 + (l-1)*nprocs                                          ! KK
!     IF(im <= member) THEN                                            ! KK
!       CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)     ! KK
!     END IF                                                           ! KK
!   END DO                                                             ! KK
    mstart = 1 + (l-1)*nprocs                                          ! KK (add 20141019)
    mend = MIN(l*nprocs, member)                                       ! KK (add 20141019)
    CALL gather_grd_mpi_alltoall(member,mstart,mend,v3d,v2d,v3dg,v2dg) ! KK (add 20141019)
    tmp_time(2)=MPI_WTIME()
    time_IO(3)=time_IO(3)+tmp_time(2)-tmp_time(1)

    !im = myrank+1 + (l-1)*nprocs
    im = myrank + (l-1)*nprocs ! 2017.05.07 Koji
    IF(im <= member-1) THEN
      WRITE(cmem,'(I6.6)') im
      filename1=trim(anal_basedir)//'/'//CDATE//'/mem'//cmem//'/restart_da' 
      CALL write_icogrd(filename1,filename2,v3dg,v2dg)
    END IF
    tmp_time(3)=MPI_WTIME()
    time_IO(4)=time_IO(4)+tmp_time(3)-tmp_time(2)

  END DO

  RETURN
END SUBROUTINE write_ens_mpi
!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  REAL(r_size),INTENT(IN) :: grd(ADM_gall,ADM_rgn_nmax)
  REAL(r_size),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,ADM_gall) + 1
      ilat = (j-ilon+1) / ADM_gall + 1
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  DO m=1,nprocs
    IF(nij1node(m) < nij1max) buf(nij1max,m) = undef
  END DO

  RETURN
END SUBROUTINE grd_to_buf
!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(buf,grd)
  REAL(r_size),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_size),INTENT(OUT) :: grd(ADM_gall,ADM_rgn_nmax)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,ADM_gall) + 1
      ilat = (j-ilon+1) / ADM_gall + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd
!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size) :: v3dm(nij1,nlev,nv3d)
  REAL(r_size) :: v2dm(nij1,nv2d)
  REAL(r_size) :: v3ds(nij1,nlev,nv3d)
  REAL(r_size) :: v2ds(nij1,nv2d)
  REAL(r_size) :: v3dg(ADM_gall,ADM_rgn_nmax,nlev,nv3d)
  REAL(r_size) :: v2dg(ADM_gall,ADM_rgn_nmax,nv2d)
  INTEGER :: i,k,m,n
  CHARACTER(256) :: filename1=''
  CHARACTER(256) :: filename2=''
  REAL(r_size) :: tmp_time(7)

    tmp_time(1)= MPI_WTIME()
  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)
    tmp_time(2)= MPI_WTIME()
    time_IO(5)=time_IO(5)+tmp_time(2)-tmp_time(1)


  CALL gather_grd_mpi(0,v3dm,v2dm,v3dg,v2dg)
    tmp_time(3)= MPI_WTIME()
    time_IO(6)=time_IO(6)+tmp_time(3)-tmp_time(2)

  IF(myrank == 0) THEN
    IF(trim(file)=='anal') THEN
      filename1=trim(anal_basedir)//'/'//CDATE//'/'//trim(file)//'_me'
      filename2=trim(anal_basedir)//'/'//CDATE//'/'//'history_'//file//'_me'
    ELSE IF(trim(file)=='gues') THEN
      filename1=trim(gues_basedir)//'/'//CDATE//'/'//trim(file)//'_me'
      filename2=trim(gues_basedir)//'/'//CDATE//'/'//'history_'//file//'_me'
    END IF
    !WRITE(filename1(1:7),'(A4,A3)') file,'_me'
    !filename2='history_'//file//'_me'
    WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename1
    CALL write_icogrd(filename1,filename2,v3dg,v2dg)
  END IF
    tmp_time(4)=MPI_WTIME()
    time_IO(7)=time_IO(7)+tmp_time(4)-tmp_time(3)

  DO n=1,nv3d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,member
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
      END DO
    END DO
!!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,member
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
    END DO
!!$OMP END PARALLEL DO
  END DO
    tmp_time(5)=MPI_WTIME()
    time_IO(8)=time_IO(8)+tmp_time(5)-tmp_time(4)

  CALL gather_grd_mpi(0,v3ds,v2ds,v3dg,v2dg)
    tmp_time(6)=MPI_WTIME()
    time_IO(9)=time_IO(9)+tmp_time(6)-tmp_time(5)

  IF(myrank == 0) THEN
    IF(trim(file)=='anal') THEN
      filename1=trim(anal_basedir)//'/'//CDATE//'/'//trim(file)//'_sp'
      filename2=trim(anal_basedir)//'/'//CDATE//'/'//'history_'//file//'_sp'
    ELSE IF(trim(file)=='gues') THEN
      filename1=trim(gues_basedir)//'/'//CDATE//'/'//trim(file)//'_sp'
      filename2=trim(gues_basedir)//'/'//CDATE//'/'//'history_'//file//'_sp'
    END IF
    !WRITE(filename1(1:7),'(A4,A3)') file,'_sp'
    !filename2='history_'//file//'_sp'
    WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename1
    CALL write_icogrd(filename1,filename2,v3dg,v2dg)
  END IF
    tmp_time(7)=MPI_WTIME()
    time_IO(10)=time_IO(10)+tmp_time(7)-tmp_time(6)

  RETURN
END SUBROUTINE write_ensmspr_mpi
!-----------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-----------------------------------------------------------------------
SUBROUTINE set_alltoallv_counts(mcount,ngpblock,n_ens,nt_ens,n_mem,nt_mem)
  INTEGER,INTENT(IN) :: mcount,ngpblock
  INTEGER,INTENT(OUT) :: n_ens(nprocs),nt_ens(nprocs),n_mem(nprocs),nt_mem(nprocs)
  INTEGER :: p

  n_ens  = 0
  nt_ens = 0
  n_mem  = 0
  nt_mem = 0
  DO p = 1, mcount
    n_ens(p) = ngpblock
    IF(myrank+1 == p) THEN
      n_mem(:) = ngpblock
    END IF
  END DO
  DO p = 2, nprocs
    nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
    nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
  END DO

  RETURN
END SUBROUTINE set_alltoallv_counts

END MODULE common_mpi_nicam
