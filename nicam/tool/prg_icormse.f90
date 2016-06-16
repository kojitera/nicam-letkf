  program main
    use common
    implicit none
    integer :: ADM_glevel
    integer :: ADM_rlevel
    integer :: ADM_gall, ADM_lall, ADM_kall
    real(8), allocatable :: true(:,:,:,:)
    real(8), allocatable :: gues(:,:,:,:)
    real(8), allocatable :: anal(:,:,:,:)
    real(8), allocatable :: diff(:,:,:,:)
    real(8), allocatable :: rmse(:,:)
    real(4), allocatable :: zgrid(:)
    character(5) :: rgnid
    character(256) :: fname
    character(256) :: anal_fname
    character(256) :: true_fname
    character(256) :: gues_fname
    character(256) :: zgrid_fname
    integer :: l, k, n
    integer :: ll, it

    namelist / filename / &
      anal_fname,         &
      gues_fname,         &
      true_fname,         &
      zgrid_fname, it

    ADM_glevel=5
    ADM_rlevel=0
    ADM_kall=40

    ADM_gall = (2**(ADM_glevel-ADM_rlevel) + 2)**2
    ADM_lall = 10*4**ADM_rlevel

    allocate( true(ADM_gall, ADM_lall, ADM_kall, 8) ) 
    allocate( gues(ADM_gall, ADM_lall, ADM_kall, 8) ) 
    allocate( anal(ADM_gall, ADM_lall, ADM_kall, 8) ) 
    allocate( diff(ADM_gall, ADM_lall, ADM_kall, 8) ) 
    allocate( rmse(ADM_kall, 8) ) 
    allocate( zgrid(ADM_kall))

    open(1,file='rmse.cnf')
    read(1,nml=filename)
    close(1)

    open(2,file=trim(zgrid_fname))
    read(2,*) zgrid(:)
    close(2)

    !write(*,*) ' READING ANALYSIS DATA'
    do l = 1, ADM_lall
      ll = l-1
      write(rgnid,'(I5.5)') ll
      fname=trim(anal_fname)//rgnid
      open(1,file=trim(fname),form='unformatted',       &
                    access='direct',          &
                    recl=ADM_gall*ADM_kall*8, &
                    status='old'              )
      read(1,rec=1) anal(:,l,:,1)
      read(1,rec=2) anal(:,l,:,2)
      read(1,rec=3) anal(:,l,:,3)
      read(1,rec=4) anal(:,l,:,4)
      read(1,rec=5) anal(:,l,:,5)
      read(1,rec=6) anal(:,l,:,6)
      read(1,rec=7) anal(:,l,:,7)
      read(1,rec=8) anal(:,l,:,8)
      close(1)
    end do
    !write(*,*) ' END READING ANALYSIS DATA'

    !write(*,*) ' READING TRUTH DATA'
    do l = 1, ADM_lall
      ll = l-1
      write(rgnid,'(I5.5)') ll
      fname=trim(true_fname)//rgnid
      open(1,file=trim(fname),form='unformatted',       &
                    access='direct',          &
                    recl=ADM_gall*ADM_kall*8, &
                    status='old'              )
      read(1,rec=1) true(:,l,:,1)
      read(1,rec=2) true(:,l,:,2)
      read(1,rec=3) true(:,l,:,3)
      read(1,rec=4) true(:,l,:,4)
      read(1,rec=5) true(:,l,:,5)
      read(1,rec=6) true(:,l,:,6)
      read(1,rec=7) true(:,l,:,7)
      read(1,rec=8) true(:,l,:,8)
      close(1)
    end do
    !write(*,*) ' END READING TRUTH DATA'

    !write(*,*) ' READING FIRST GUESS DATA'
    do l = 1, ADM_lall
      ll = l-1
      write(rgnid,'(I5.5)') ll
      fname=trim(gues_fname)//rgnid
      open(1,file=trim(fname),form='unformatted',       &
                    access='direct',          &
                    recl=ADM_gall*ADM_kall*8, &
                    status='old'              )
      read(1,rec=1) gues(:,l,:,1)
      read(1,rec=2) gues(:,l,:,2)
      read(1,rec=3) gues(:,l,:,3)
      read(1,rec=4) gues(:,l,:,4)
      read(1,rec=5) gues(:,l,:,5)
      read(1,rec=6) gues(:,l,:,6)
      read(1,rec=7) gues(:,l,:,7)
      read(1,rec=8) gues(:,l,:,8)
      close(1)
    end do
    !write(*,*) ' END READING FIRST GUESS DATA'

    anal(:,:,:,1)=anal(:,:,:,1)*0.01d0
    true(:,:,:,1)=true(:,:,:,1)*0.01d0
    gues(:,:,:,1)=gues(:,:,:,1)*0.01d0

    diff(:,:,:,:)=anal(:,:,:,:)-true(:,:,:,:)

    do n = 1, 8
      do k = 1, ADM_kall
        call com_rms(ADM_gall*ADM_lall, diff(:,:,k,n), rmse(k,n))
        write(1,'(F15.7,2I6,2F15.7)') real(it/4.0), n, k, zgrid(k), rmse(k,n)
      end do
    end do

    diff(:,:,:,:)=gues(:,:,:,:)-true(:,:,:,:)

    do n = 1, 8
      do k = 1, ADM_kall
        call com_rms(ADM_gall*ADM_lall, diff(:,:,k,n), rmse(k,n))
        write(2,'(F15.7,2I6,2F15.7)') real(it/4.0), n, k, zgrid(k), rmse(k,n)
      end do
    end do


  end program main
