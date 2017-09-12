PROGRAM vbc
  USE mpi
  USE mod_adm
  USE mod_comm
  USE mod_cnst
  USE mod_grd, ONLY : GRD_setup
  USE mod_gmtr, ONLY : &
      GMTR_setup,      &
      GMTR_P_var,      &
      GMTR_P_LAT,      &
      GMTR_P_LON
  USE mod_oprt, ONLY : OPRT_setup
  USE mod_vmtr, ONLY : VMTR_setup
  USE mod_gtl, ONLY :          &
       GTL_clip_region,        &
       GTL_clip_region_1layer, &
       GTL_input_var2

  USE mod_obsope_common
  USE mod_read_history
  USE mod_obsope_driver
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  INTEGER :: imem
  INTEGER :: ierr
  REAL(8) :: tmp_time(10)
  REAL(8) :: time(10)

  CALL ADM_proc_init(ADM_MULTI_PRC)
  CALL ADM_setup('obsope.cnf')
  CALL FIO_setup
  CALL COMM_setup
  CALL CNST_setup
  CALL GRD_setup
  CALL GMTR_setup
  CALL OPRT_setup
  CALL VMTR_setup

  time(:)=0.0d0
  tmp_time(1)=MPI_WTIME()

  CALL obsope_init

  tmp_time(2)=MPI_WTIME()
  time(1)=tmp_time(2)-tmp_time(1)

  CALL read_history_init
  WRITE(ADM_LOG_FID,*) 'read_history_init'
  IF(flush_text) FLUSH(ADM_LOG_FID)

  tmp_time(3)=MPI_WTIME()
  time(2)=tmp_time(3)-tmp_time(2)

  tmp_time(4)=MPI_WTIME()

  CALL read_history(-1)
  WRITE(ADM_LOG_FID,*) 'read_history'
  IF(flush_text) FLUSH(ADM_LOG_FID)

  tmp_time(5)=MPI_WTIME()
  time(3)=time(3)+ (tmp_time(5)-tmp_time(4))

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(6)=MPI_WTIME()
  time(4)=time(4)+ (tmp_time(6)-tmp_time(5))

  CALL obsope_main(-1)
  WRITE(ADM_LOG_FID,*) 'obsope_main'
  IF(flush_text) FLUSH(ADM_LOG_FID)

  tmp_time(7)=MPI_WTIME()
  time(5)=time(5)+ (tmp_time(7)-tmp_time(6))

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tmp_time(8)=MPI_WTIME()
  time(6)=time(6)+ (tmp_time(8)-tmp_time(7))

  WRITE(ADM_LOG_FID,*) 'obsope_init          = ', time(1)
  WRITE(ADM_LOG_FID,*) 'read_history_init    = ', time(2)
  WRITE(ADM_LOG_FID,*) 'read_history         = ', time(3)
  WRITE(ADM_LOG_FID,*) 'MPI_BARRIER          = ', time(4)
  WRITE(ADM_LOG_FID,*) 'obsope_main          = ', time(5)
  WRITE(ADM_LOG_FID,*) 'MPI_BARRIER          = ', time(6)

  CALL ADM_proc_finish  

END PROGRAM vbc
