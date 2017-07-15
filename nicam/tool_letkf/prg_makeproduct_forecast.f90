  program main
    use mpi
    use mod_adm, only :         & 
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
    use mod_misc
    use mod_fio, only : &  ! [add] H.Yashiro 20110826
      FIO_setup,        &
      FIO_input,        &
      FIO_output,       &
      FIO_HMID,         &
      FIO_HSHORT,       &
      FIO_HLONG,        &
      FIO_INTEG_FILE,   &
      FIO_SPLIT_FILE,   &
      FIO_FREAD,        &
      FIO_ICOSAHEDRON,  &
      FIO_BIG_ENDIAN,   &
      headerinfo,       &
      datainfo,         &
      FIO_REAL8
    use mod_comm, only :           &
         COMM_setup
    use mod_latlon, only :        &
         LATLON_setup,            &
         LATLON_read_outdirname
    use mod_cnst, only :           &
         CNST_setup,               &
         CNST_UNDEF
    use mod_grd, only :            &
         GRD_setup,                &
         GRD_gz,                   &
         GRD_gzh,                  &
         GRD_zs,                   &
         GRD_zs_pl,                &
         GRD_vz,                   &
         GRD_vz_pl
    use mod_oprt, only :         &
         OPRT_setup
    use mod_vmtr, only :         &
         VMTR_setup
    use mod_gtl, only :            &
         GTL_generate_uv
    use mod_gmtr, only :           &
         GMTR_setup

    implicit none
    integer :: nmem
    integer :: imem
    integer :: istep
    integer :: g, l, rgnid, n, nv, ip
    integer :: ierr
    integer :: i, j
    character(ADM_MAXFNAME) :: basename
    character(ADM_MAXFNAME) :: fname
    character(ADM_MAXFNAME) :: history_basename

    character(ADM_MAXFNAME) :: output_basename
    character(ADM_MAXFNAME) :: restart_layername
    character(ADM_MAXFNAME) :: surface_layername
    character(ADM_MAXFNAME) :: pressure_layername

    character(3) :: ens_num
    real(8), allocatable :: ddd_logp(:,:,:,:)
    real(8), allocatable :: ddd_3d_org(:,:,:,:)
    real(8), allocatable :: ddd_3d(:,:,:,:,:)
    real(8), allocatable :: ddd_2d(:,:,:,:,:)
    real(8), allocatable :: ddd_uv(:,:,:,:)
    real(8), allocatable :: ddd3_pl(:,:,:,:)
    real(4), allocatable :: sprd_uv(:,:,:,:)

    !--- 
    real(8), allocatable :: ddd_3d_plev(:,:,:,:,:)

    real(8) :: plevel(200)
    real(8) :: logp(200)

    real(8), allocatable :: cc(:)
    real(8), allocatable :: ww(:)
    real(8), allocatable :: bx(:)
    real(8), allocatable :: by(:)

    INTEGER :: np
    INTEGER,PARAMETER :: nv3d=8
    INTEGER,PARAMETER :: nv2d=1
    character(LEN=FIO_HSHORT) :: varname3d(nv3d)
    character(LEN=FIO_HSHORT) :: varname2d(nv2d)
    integer, allocatable :: ifid(:)

    !--- ico2ll
    integer :: imax, jmax
    integer,allocatable :: lon_index(:),lat_index(:)
    integer,allocatable :: n1_index(:), n2_index(:), n3_index(:)
    real(8),allocatable :: w1(:), w2(:), w3(:)

    integer,allocatable :: max_num_latlon(:)
    integer,allocatable :: nstart(:), nend(:)
    character(1024) :: fin_llmap_head = '../llmap/llmap'
    character(1024) :: fin_llmap

    real(8),allocatable :: lon(:)  ! -pi - pi
    real(8),allocatable :: lat(:)  ! -pi/2 - pi/2
    real(8),allocatable :: lon_swap(:)

    real(8),allocatable :: var_out_3d(:,:,:,:,:)
    real(8),allocatable :: var_out_2d(:,:,:,:,:)
    real(8),allocatable :: var_3d(:,:,:,:,:)
    real(8),allocatable :: var_2d(:,:,:,:,:)
    real(8),allocatable :: var_out_swap(:,:)
    character(5) :: c5
 
    logical :: onefile=.false.
    character(256) :: output_dir
    character(10)  :: cdate

    character(LEN=FIO_HMID) :: desc = 'INITIAL/RESTART DATA of PROGNOSTICVARIABLES'

    namelist / RESTARTPARAM /    &
         istep,                  &
         history_basename,       &
         output_basename,        &
         varname3d,              &
         varname2d,              &
         plevel,                 &
         np,                     &
         fin_llmap_head,         &
         output_dir,             &
         cdate, onefile

    call ADM_proc_init(ADM_MULTI_PRC)
    call ADM_setup('nhm_driver.cnf')
    call COMM_setup
    call FIO_setup
    call GRD_setup
    call GMTR_setup
    call OPRT_setup
    call VMTR_setup

    allocate( ddd_3d_org(ADM_gall,    ADM_kall,    ADM_lall,    3) ) ! mean vx, vy, vz
    allocate( ddd_3d    (ADM_gall,    ADM_kall,    ADM_lall, nv3d, 1) )
    allocate( ddd_uv    (ADM_gall,    ADM_kall,    ADM_lall,    2) )
    allocate( ddd_2d    (ADM_gall,           1,    ADM_lall, nv2d, 1) )
    allocate( ddd3_pl   (ADM_gall_pl, ADM_kall, ADM_lall_pl,    5) )
    allocate( ddd_logp  (ADM_gall,    ADM_kall,    ADM_lall,    2) )
    allocate(    sprd_uv(ADM_gall,    ADM_kall,    ADM_lall,    2) )

    open(1,file='nhm_driver.cnf')
    read(1,nml=RESTARTPARAM) 
    close(1)
    write(ADM_LOG_FID,RESTARTPARAM)
    FLUSH(ADM_LOG_FID)
    write(ADM_LOG_FID,*) trim(varname2d(1)), trim(varname2d(2))
    write(ADM_LOG_FID,*) 'history_basename,', trim(history_basename)
    FLUSH(ADM_LOG_FID)
    logp(1:np)=log(plevel(1:np))
    allocate( ddd_3d_plev(ADM_gall, np, ADM_lall, nv3d, 1) )

    restart_layername='ZSDEF38'
    surface_layername='ZSSFC1'
    pressure_layername='ZSDEF26'
    !
    !--- < Reading data > ---
    call FIO_setup ! [add] H.Yashiro 20110826

    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,1,1),history_basename,'ms_pres',restart_layername,1,ADM_kall-2,istep)
    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,2,1),history_basename,'ms_tem', restart_layername,1,ADM_kall-2,istep)
    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,3,1),history_basename,'ms_u',   restart_layername,1,ADM_kall-2,istep)
    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,4,1),history_basename,'ms_v',   restart_layername,1,ADM_kall-2,istep)
    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,5,1),history_basename,'ms_w',   restart_layername,1,ADM_kall-2,istep)
    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,6,1),history_basename,'ms_qv',  restart_layername,1,ADM_kall-2,istep)
    call FIO_input(    ddd_3d(:,2:ADM_kall-1,:,7,1),history_basename,'ms_qc',  restart_layername,1,ADM_kall-2,istep)

    do nv = 1, nv2d
      call FIO_input(ddd_2d(:,:,:,nv,1),history_basename,trim(varname2d(nv)), surface_layername,1,1,istep)
      write(ADM_LOG_FID,'(A15,2F15.5)') trim(varname2d(nv)), &
            minval(ddd_2d(:,1,:,nv,1)), maxval(ddd_2d(:,1,:,nv,1))
    end do

    ddd_logp(:,:,:,:)=log(ddd_3d(:,:,:,1,:)*0.01)
    ddd_3d(:,:,:,1,:)=ddd_3d(:,:,:,1,:)*0.01d0
    do n = 1, ADM_kall
      write(ADM_LOG_FID,*) n, ddd_logp(1,n,1,1)
    end do

    allocate( cc(ADM_kall-2) )
    allocate( ww(ADM_kall-2) )
    allocate( bx(np) )
    allocate( by(np) )


    do i = 1, 1
      do l = 1, ADM_lall
      do g = 1, ADM_gall
        if( ddd_3d(g,2,l,1,1) >= 1000.0d0 ) then
          do nv = 2, nv3d-1
            call spline( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                         ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, ww)
            call splint( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                         ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, np, logp(np:1:-1), by(np:1:-1))
            ddd_3d_plev(g,:,l,nv,i)=by(:)
          end do
  
          call spline( ADM_kall-2, ddd_logp(g, ADM_kall-1:2:-1, l, 1), &
                       GRD_vz(g, ADM_kall-1:2:-1, l, 1), cc, ww)
          call splint( ADM_kall-2, ddd_logp(g, ADM_kall-1:2:-1, l, 1), &
                       GRD_vz(g, ADM_kall-1:2:-1, l, 1), cc, np, logp(np:1:-1), by(np:1:-1))
          ddd_3d_plev(g,:,l,1,i)=by(:)
        else
          n = 1
          do 
            if(ddd_3d(g,2,l,1,1) >= plevel(n) ) then
              exit
            end if
            n = n + 1
          end do 
          do nv = 2, nv3d-1
            call spline( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                         ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, ww)
            call splint( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                         ddd_3d(g,ADM_kall-1:2:-1,l,nv,i), cc, np-n+1, logp(np:n:-1), by(np:n:-1))
            ddd_3d_plev(g,n:np,l,nv,i)=by(n:np)
            ddd_3d_plev(g,1:n-1,l,nv,i)=CNST_UNDEF
          end do
          call spline( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                       GRD_vz(g,ADM_kall-1:2:-1,l,1), cc, ww)
          call splint( ADM_kall-2, ddd_logp(g,ADM_kall-1:2:-1,l,1), &
                       GRD_vz(g,ADM_kall-1:2:-1,l,1), cc, np-n+1, logp(np:n:-1), by(np:n:-1))
          ddd_3d_plev(g,n:np,l,1,i)=by(n:np)
          ddd_3d_plev(g,1:n-1,l,1,i)=CNST_UNDEF
        end if
      end do
      end do
    end do

    fin_llmap =  trim(fin_llmap_head) // '.info'
    open(10,file=fin_llmap, &
         form='unformatted',status='old' ,iostat=ierr)
    if(ierr/=0) then
       write(*,*) 'Cannot open llmap info file!'
       stop
    endif
    read(10) imax
    allocate(lon(imax))
    allocate(lon_swap(imax))
    read(10) lon(:)
    read(10) jmax
    allocate(lat(jmax))
    read(10) lat(:)
    close(10)
    write(ADM_LOG_FID,*) 'imax =', imax
    write(ADM_LOG_FID,*) 'jmax =', jmax

    allocate(max_num_latlon(ADM_lall))
    allocate(nstart(ADM_lall))
    allocate(nend(ADM_lall))
  
    allocate(lon_index(imax*jmax))
    allocate(lat_index(imax*jmax))
    allocate(n1_index(imax*jmax))
    allocate(n2_index(imax*jmax))
    allocate(n3_index(imax*jmax))
    allocate(w1(imax*jmax))
    allocate(w2(imax*jmax))
    allocate(w3(imax*jmax))
  
    allocate(var_out_3d(imax,jmax,np,nv3d,2))
    allocate(var_out_2d(imax,jmax,1,nv2d,2))
    allocate(var_3d(imax,jmax,np,nv3d,2))
    allocate(var_2d(imax,jmax,1,nv2d,2))
    allocate(var_out_swap(imax,jmax))

    do l=1, ADM_lall
       open(20,file=fin_llmap,form='unformatted',status='old')
       read(20) max_num_latlon(l)
       close(20)
    enddo

    do l=1, ADM_lall
      rgnid = ADM_prc_tab(l,ADM_prc_me)
      write(c5,'(I5)') rgnid-1
      do n=1, 5
        if( c5(n:n) == ' ' ) c5(n:n) = '0'
      enddo

      fin_llmap = trim(fin_llmap_head) // '.rgn' // c5
      open(20,file=fin_llmap,form='unformatted',status='old')
      read(20) max_num_latlon(l)
      nend(l) = sum(max_num_latlon(1:l))
      nstart(l) = nend(l) - max_num_latlon(l)+1
      if(max_num_latlon(l)/=0) then
        read(20) lon_index(nstart(l):nend(l))  ! lon(lon_index(i)) :
        read(20) lat_index(nstart(l):nend(l))
        read(20) n1_index(nstart(l):nend(l))
        read(20) n2_index(nstart(l):nend(l))
        read(20) n3_index(nstart(l):nend(l))
        read(20) w1(nstart(l):nend(l))
        read(20) w2(nstart(l):nend(l))
        read(20) w3(nstart(l):nend(l))
      endif
      close(20)

      var_out_3d(:,:,:,:,:)=0.0d0
      var_out_2d(:,:,:,:,:)=0.0d0
      do i = 1, 1
        do nv = 1, nv3d-1
        do ip = 1, np
        do n=nstart(l),nend(l)
          if(ddd_3d_plev(n1_index(n),ip,l,nv,i) == CNST_UNDEF ) then
             var_out_3d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
          else if(ddd_3d_plev(n2_index(n),ip,l,nv,i) == CNST_UNDEF ) then
             var_out_3d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
          else if(ddd_3d_plev(n3_index(n),ip,l,nv,i) == CNST_UNDEF) then
             var_out_3d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
          else
             var_out_3d(lon_index(n),lat_index(n),ip,nv,i) &
                  = w1(n)*ddd_3d_plev(n1_index(n),ip,l,nv,i)      &
                  + w2(n)*ddd_3d_plev(n2_index(n),ip,l,nv,i)      &
                  + w3(n)*ddd_3d_plev(n3_index(n),ip,l,nv,i) 
          endif
        end do
        end do
        end do
  
        ip=1
        do nv = 1, nv2d
        do n=nstart(l),nend(l)
          if(     ddd_2d(n1_index(n),ip,l,nv,i) == CNST_UNDEF ) then
             var_out_2d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
          else if(ddd_2d(n2_index(n),ip,l,nv,i) == CNST_UNDEF ) then
             var_out_2d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
          else if(ddd_2d(n3_index(n),ip,l,nv,i) == CNST_UNDEF) then
             var_out_2d(lon_index(n),lat_index(n),ip,nv,i) = CNST_UNDEF
          else
             var_out_2d(lon_index(n),lat_index(n),ip,nv,i) &
                  = w1(n)*ddd_2d(n1_index(n),ip,l,nv,i)      &
                  + w2(n)*ddd_2d(n2_index(n),ip,l,nv,i)      &
                  + w3(n)*ddd_2d(n3_index(n),ip,l,nv,i)
          endif
        end do
        end do
      end do

    end do

    call MPI_ALLREDUCE(var_out_3d, var_3d, imax*jmax*np*nv3d*2, MPI_REAL8, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE(var_out_2d, var_2d, imax*jmax*nv2d*2, MPI_REAL8, MPI_SUM,&
                       MPI_COMM_WORLD, ierr)


    !--- OUTPUT
    if(ADM_prc_me==1) then
      if(onefile) then
        fname=trim(output_dir)//'/'//trim(cdate)//'.dat'
        open(1,file=trim(fname),form='unformatted',access='direct',recl=4*imax*jmax)
        n=1
        do i = 1, 1
          do nv = 1, nv3d-1
            do ip = 1, np
              write(1,rec=n) sngL(var_3d(:,:,ip,nv,i))
              write(ADM_LOG_FID,*) ip, minval(var_3d(:,:,ip,nv,i)), maxval(var_3d(:,:,ip,nv,i))
              n=n+1
            end do
            write(ADM_LOG_FID,*) trim(varname3d(nv)), minval(var_3d(:,:,:,nv,i)), maxval(var_3d(:,:,:,nv,i))
          end do

          do nv = 1, nv2d
            write(1,rec=n) sngL(var_2d(:,:,1,nv,i))
            write(ADM_LOG_FID,*) trim(varname2d(nv)), minval(var_2d(:,:,:,nv,i)), maxval(var_2d(:,:,:,nv,i))
            n=n+1
          end do
        end do

        close(1)

      else
        do nv = 1, nv3d-1
          fname=trim(output_dir)//'/'//trim(varname3d(nv))//'_'//trim(cdate)//'.dat'
          open(1,file=trim(fname),form='unformatted',access='direct',recl=4*imax*jmax*np)
          write(1,rec=1) sngL(var_3d(:,:,:,nv,i))
          close(1)
        end do
    
        do nv = 1, nv2d
          fname=trim(output_dir)//'/'//trim(varname2d(nv))//'_'//trim(cdate)//'.dat'
          open(1,file=trim(fname),form='unformatted',access='direct',recl=4*imax*jmax)
          write(1,rec=1) sngL(var_2d(:,:,1,nv,i))
          close(1)
        end do
      end if
      
    end if

    call ADM_proc_finish

  end program main

      SUBROUTINE SPLINE(N,X,Y,CC,WW)
      IMPLICIT REAL*8(A-H),REAL*8(O-Z)
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
      RETURN
      END



      SUBROUTINE SPLINT(N,AX,AY,CC,M,BX,BY)
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
      RETURN
      END


