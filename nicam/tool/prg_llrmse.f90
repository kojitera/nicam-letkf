  program main
    use common
    implicit none
    integer :: nx, ny ,np
    real(8), allocatable :: true(:,:,:,:)
    real(8), allocatable :: gues(:,:,:,:)
    real(8), allocatable :: anal(:,:,:,:)
    real(8), allocatable :: diff(:,:,:,:)
    real(8), allocatable :: rmse(:,:)
    real(4), allocatable :: zgrid(:)
    real(4), allocatable :: ygrid(:)
    real(4), allocatable :: cs(:)
    character(5) :: rgnid
    character(256) :: fname
    character(256) :: anal_fname
    character(256) :: true_fname
    character(256) :: gues_fname
    character(256) :: ygrid_fname
    character(256) :: zgrid_fname
    integer :: l, k, n
    integer :: ll, it, irec, iy, ix
    integer :: tmp_nx

    namelist / filename / &
      nx, ny, np,         &
      anal_fname,         &
      gues_fname,         &
      true_fname,         &
      ygrid_fname,        &
      zgrid_fname, it

    allocate( true(nx,ny,np,7) ) 
    allocate( gues(nx,ny,np,7) ) 
    allocate( anal(nx,ny,np,7) ) 
    allocate( diff(nx,ny,np,7) ) 
    allocate( rmse(np,7) ) 
    allocate( zgrid(np))
    allocate( ygrid(ny))
    allocate( cs(ny+1))
 
    open(1,file='rmse.cnf')
    read(1,nml=filename)
    close(1)

    open(2,file=trim(zgrid_fname))
    read(2,*) zgrid(:)
    close(2)

    open(3,file=trim(ygrid_fname))
    read(3,*) ygrid(:)
    close(3)

    cs(1:ny)=ygrid(:)/180.0*pi

    !write(*,*) ' READING ANALYSIS DATA'
    irec=4*nx*ny*np*7
    open(11,file=trim(anal_fname),form='unformatted',access='direct',recl=irec)
    read(11,rec=1) anal(:,:,:,:)
    close(11)

    open(12,file=trim(true_fname),form='unformatted',access='direct',recl=irec)
    read(12,rec=1) true(:,:,:,:)
    close(12)

    diff(:,:,:,:)=anal(:,:,:,:)-true(:,:,:,:)

    do n = 1, 7
      do k = 1, np
        do iy = 1, ny
          do ix = 1, nx
            if( abs(diff(ix,iy,k,n)).gt.1.0d10 ) then 
              rmse(k,n)=rmse(k,n)+ diff(ix,iy,k,n)**2
              !rmse(k,n)=rmse(k,n)+sum( sqrt( diff(:,iy,k,n)**2 / float(nx) ) ) * cs(iy)
            end if
          end do
          if(tmp_nx.gt.0) then
            rmse(k,n)=rmse(k,n)/float(tmp_nx)
            rmse(k,n)=sqrt(rmse(k,n))
            rmse(k,n)=rmse(k,n)*cs(iy)
          else
            rmse(k,n)=-999.0
          end if
        end do
        
        rmse(k,n)=rmse(k,n)/cs(ny+1)
        write(1,'(F15.7,2I6,2F15.7)') real(it/4.0), n, k, zgrid(k), rmse(k,n)
      end do
    end do

!    diff(:,:,:,:)=gues(:,:,:,:)-true(:,:,:,:)
!
!    do n = 1, 8
!      do k = 1, ADM_kall
!        call com_rms(ADM_gall*ADM_lall, diff(:,:,k,n), rmse(k,n))
!        write(2,'(F15.7,2I6,2F15.7)') real(it/4.0), n, k, zgrid(k), rmse(k,n)
!      end do
!    end do


  end program main
