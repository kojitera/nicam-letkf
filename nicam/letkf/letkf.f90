PROGRAM letkf
!=======================================================================
!
! [PURPOSE:] Main program of LETKF
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_nicam
  USE common_mpi_nicam
  USE common_letkf
  USE letkf_obs
  USE letkf_tools

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: time_main(:)
  REAL(r_size),ALLOCATABLE :: time(:,:)
  REAL(r_size),ALLOCATABLE :: IO_time(:,:)
  REAL(r_size) :: tmp_time(100)
  REAL(r_size) :: rtimer_tmp1, rtimer_tmp2
  integer, parameter :: slot=20                                ! KK
  REAL(r_size),allocatable :: rtimer_all(:,:), rtimer_tmp(:)   ! KK
  INTEGER :: ierr, i
  CHARACTER(8) :: stdoutf='NOUT-000'
  CHARACTER(4) :: guesf='gs00'
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  tmp_time(1)=MPI_WTIME()
  CALL initialize_mpi
  tmp_time(2)=MPI_WTIME()

  ALLOCATE( time_main(9) )
  ALLOCATE( time(9,nprocs) )
  ALLOCATE( IO_time(10,nprocs) )
!
  CALL set_common_nicam
  WRITE(stdoutf(6:8), '(I3.3)') myrank
  WRITE(ADM_LOG_FID,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  WRITE(ADM_LOG_FID,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(ADM_LOG_FID,'(A)') '============================================='
  WRITE(ADM_LOG_FID,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  WRITE(ADM_LOG_FID,'(A)') '                                             '
  WRITE(ADM_LOG_FID,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(ADM_LOG_FID,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(ADM_LOG_FID,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(ADM_LOG_FID,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(ADM_LOG_FID,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(ADM_LOG_FID,'(A)') '                                             '
  WRITE(ADM_LOG_FID,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(ADM_LOG_FID,'(A)') '                                             '
  WRITE(ADM_LOG_FID,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(ADM_LOG_FID,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(ADM_LOG_FID,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(ADM_LOG_FID,'(A)') '============================================='
  WRITE(ADM_LOG_FID,'(A)') '              LETKF PARAMETERS               '
  WRITE(ADM_LOG_FID,'(A)') ' ------------------------------------------- '
  WRITE(ADM_LOG_FID,'(A,I15)')   '   nbv        :',nbv
  WRITE(ADM_LOG_FID,'(A,I15)')   '   nslots     :',nslots
  WRITE(ADM_LOG_FID,'(A,I15)')   '   nbslot     :',nbslot
  WRITE(ADM_LOG_FID,'(A,F15.2)') '   sigma_obs  :',sigma_obs
  WRITE(ADM_LOG_FID,'(A,F15.2)') '   sigma_obsv :',sigma_obsv
  WRITE(ADM_LOG_FID,'(A,F15.2)') '   sigma_obst :',sigma_obst
  WRITE(ADM_LOG_FID,'(A)') '============================================='
  CALL set_common_mpi_nicam
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))
  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(anal2d(nij1,nbv,nv2d))
  tmp_time(3)=MPI_WTIME()
!!
!-----------------------------------------------------------------------
! Observations
!-----------------------------------------------------------------------
  !
  ! CONVENTIONAL OBS
  !
  CALL set_letkf_obs
  tmp_time(4)=MPI_WTIME()
!!

  CALL set_letkf_tvs
  tmp_time(5)=MPI_WTIME()
!!
!
!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------
  !
  ! READ GUES
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(guesf(3:4),'(I2.2)') nbslot
  CALL read_ens_mpi(guesf,nbv,gues3d,gues2d)
  tmp_time(6)=MPI_WTIME()
!!
!  !
!  ! WRITE ENS MEAN and SPRD
!  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('gues',nbv,gues3d,gues2d)
  tmp_time(7)=MPI_WTIME()
!!
!!
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
  !
  ! LETKF
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL das_letkf(gues3d,gues2d,anal3d,anal2d)
  tmp_time(8)=MPI_WTIME()
!
!
!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------
  !
  ! WRITE ANAL
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ens_mpi('anal',nbv,anal3d,anal2d)
  tmp_time(9)=MPI_WTIME()
!!
  !
  ! WRITE ENS MEAN and SPRD
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('anal',nbv,anal3d,anal2d)
  tmp_time(10)=MPI_WTIME()
!!
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
!  CALL monit_mean('gues')
!  CALL monit_mean('anal')
!

  do i = 1, 9
    time_main(i)=tmp_time(i+1)-tmp_time(i)
  end do

  call MPI_Allgather(time_main(1), 9,  MPI_DOUBLE, &
                     time(1,1),    9,  MPI_DOUBLE, MPI_COMM_WORLD, ierr)
  call MPI_Allgather(time_IO(1),   10, MPI_DOUBLE, &
                     IO_time(1,1), 10, MPI_DOUBLE, MPI_COMM_WORLD, ierr)

  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### Initialize: ', time_main(1), maxval(time(1,:)),&
                                                  sum(time(1,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### setting:    ', time_main(2), maxval(time(2,:)),&
                                                  sum(time(2,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### read PREP:  ', time_main(3), maxval(time(3,:)),&
                                                  sum(time(3,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### read amsua  ', time_main(4), maxval(time(4,:)),&
                                                  sum(time(4,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### read gues:  ', time_main(5), maxval(time(5,:)),&
                                                  sum(time(5,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### write gues: ', time_main(6), maxval(time(6,:)),&
                                                  sum(time(6,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### DAS LETKF:  ', time_main(7), maxval(time(7,:)),&
                                                  sum(time(7,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### write analy:', time_main(8), maxval(time(8,:)),&
                                                  sum(time(8,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### write me sp:', time_main(9), maxval(time(9,:)),&
                                                  sum(time(9,:))/dble(nprocs)
  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### Total time: ', sum(time_main(1:9))

!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### read_ens:   ', time_IO(1), &
!                         maxval(IO_time(1,:)), sum(IO_time(1,:))/64.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### scatter_ens:', time_IO(2), &
!                         maxval(IO_time(2,:)), sum(IO_time(2,:))/640.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### gather_ens: ', time_IO(3), &
!                         maxval(IO_time(3,:)), sum(IO_time(3,:))/640.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### write_ens:  ', time_IO(4), &
!                         maxval(IO_time(4,:)), sum(IO_time(4,:))/64.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### ensem_mean: ', time_IO(5), &
!                         maxval(IO_time(5,:)), sum(IO_time(5,:))/640.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### gather_ens: ', time_IO(6), &
!                         maxval(IO_time(6,:)), sum(IO_time(6,:))/640.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### write_mean: ', time_IO(7), &
!                         maxval(IO_time(7,:)), sum(IO_time(7,:))
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### calc_spread:', time_IO(8), &
!                         maxval(IO_time(8,:)), sum(IO_time(8,:))/640.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### gather_sprd:', time_IO(9), &
!                         maxval(IO_time(9,:)), sum(IO_time(9,:))/640.0d0
!  WRITE(ADM_LOG_FID,'(A,3F15.5)')  '####### write_sprd: ', time_IO(10), &
!                         maxval(IO_time(10,:)), sum(IO_time(10,:))
  FLUSH(ADM_LOG_FID)

!!-----------------------------------------------------------------------
!! Finalize
!!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  STOP
END PROGRAM letkf
