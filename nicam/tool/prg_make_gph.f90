program make_gph
  implicit none
  !integer,parameter :: nx = 320
  !integer,parameter :: ny = 160
  integer,parameter :: nz = 38
  integer,parameter :: np = 26
  integer :: nt = 1

  integer :: nx = 240
  integer :: ny = 121

  real(4) :: z(nz)
  real(4) :: p(np)
  real(8) :: lnp(np)
  real(8) :: z1d_full_tmp(nz+2)
  real(8) :: z1d_half_tmp(nz+2)
  real(4) :: z1d_full(nz)
  real(4) :: ztop

  !
  !real(4) :: prs(nx,ny,nz)
  !real(4) :: zs(nx,ny)
  !real(8) :: lnprs(nx,ny,nz)
  !real(4) :: gph(nx,ny,nz)
  real(4),allocatable :: prs(:,:,:)
  real(4),allocatable :: zs(:,:)
  real(8),allocatable :: lnprs(:,:,:)
  real(4),allocatable :: gph(:,:,:)
  
  !
  character(len=256) :: fin
  character(len=256) :: fout = "gph.grd"
  character(len=256) :: fin_vgrid, fin_zs
  integer :: fid = 101 ! vgrid,zs,fin
  integer :: omt = 201 ! fout
  integer :: kall
  integer :: irec_prs = 1
  integer :: ierr

  !
  integer :: i, j, k, it
  integer :: irec

  !
  namelist/param/ nx, ny, &
       & fin, fout, fin_vgrid, fin_zs, irec_prs

  z(1:nz) =(/   80.84,   248.82,   429.88,   625.04,   835.41,  1062.16, &
 &  1306.57,  1570.01,  1853.97,  2160.05,  2489.96,  2845.57,  3228.88, &
 &  3642.04,  4087.38,  4567.41,  5084.82,  5642.53,  6243.68,  6891.64, &
 &  7590.07,  8342.90,  9154.37, 10029.03, 10971.82, 11988.03, 13083.39, &
 & 14264.06, 15536.68, 16908.43, 18387.01, 19980.75, 21698.62, 23550.28, &
 & 25546.15, 28113.21, 31747.21, 36734.94/)

  p(1:np) = (/ 1000.0, 975.0, 950.0, 925.0, 900.0, 850.0, 800.0, &
    &    750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0, &
    &    350.0, 300.0, 250.0, 200.0, 150.0, 100.0,  70.0,  50.0, &
    &     30.0,  20.0,  10.0/)

  p = 100.0 * p
  lnp(1:np) = log(p(1:np))
  !lnp(1:np) = p(1:np)

  read(5,param)
  write(6,param)

  allocate( zs(nx,ny) )
  allocate( prs(nx,ny,nz) )
  allocate( lnprs(nx,ny,nz) )
  allocate( gph(nx,ny,np) )

  open(fid, file=trim(fin_vgrid), form='unformatted', access='sequential', &
       status='old', action='read', iostat=ierr)
  if (ierr /= 0) stop "error in opening file given by vgrid"
    
  read(fid) kall
  read(fid) z1d_full_tmp(:)
  read(fid) z1d_half_tmp(:)
  ztop = z1d_half_tmp(nz+2)
  print'(A,I4,F9.2)','kall, ztop: ',kall, ztop
  !read(fid) z_half(:)
  do k = 1,nz
    z1d_full(k) = z1d_full_tmp(k+1)
    !print'(I4,3F9.2)',k,z1d_full(k)!,z1d_full_tmp(k+1),z1d_half_tmp(k+1)
  enddo
  print'(A,I4,2F9.2)','z1d_full: ',k,z1d_full(1),z1d_full(nz)
  close(fid)

  open(fid, file=trim(fin_zs), form='unformatted', access='direct', &
       status='old', action='read', recl=nx*ny*4,iostat=ierr)
  if (ierr /= 0) stop "error in opening file given by zs(topog)"
  read(fid,rec=1) zs
  print'(A,2F9.2)','zs: ',minval(zs),maxval(zs)
  close(fid)

  open(fid, file=trim(fin), form='unformatted', access='direct', &
       status='old', action='read', recl=nx*ny*nz*4,iostat=ierr)

  open(omt, file=trim(fout), form='unformatted', access='direct', &
       status='unknown', action='write', recl=nx*ny*np*4,iostat=ierr)

  do it = 1,nt
    ! --- read data
    irec = irec_prs
    read(fid,rec=irec) prs
    !lnprs(:,:,:)= log(prs(:,:,:))
    lnprs = prs
    write(*,'(A,2F9.2)') 'prs: ',0.01*minval(prs),0.01*maxval(prs)
    write(*,'(A,2F9.2)') 'lnprs: ',minval(lnprs),maxval(lnprs)

    ! --- convert
    call conv_prs2gph(nx,ny,nz,np,ztop,z,lnp,zs,lnprs,gph)
    
    ! --- output
    write(omt,rec=it) gph

  enddo

end program make_gph

!=============================================================================
subroutine conv_prs2gph(nx,ny,nz,np,ztop,z,lnp,zs,lnprs,gph)
  implicit none
  real(4),parameter :: undef = -9.99e34
  !
  integer,intent(in) :: nx, ny, nz, np
  real(4),intent(in) :: ztop
  real(4),intent(in) :: z(nz)
  real(8),intent(in) :: lnp(np)
  real(4),intent(in) :: zs(nx,ny)
  real(8),intent(in) :: lnprs(nx,ny,nz)
  real(4),intent(out) :: gph(nx,ny,np)

  real(4),allocatable :: zreal(:)
  real(4) :: fac
  real(4) :: a
  integer :: i, j, k, kk, ks

  allocate( zreal(nz) )

  do j = 1, ny
    zreal(1:nz) = zs(1,j) + ((ztop - zs(1,j)) / ztop) * z(1:nz)
    print'(A,3F9.2)','zreal: ',zreal(1),zreal(nz),exp(lnprs(1,j,1))
  end do

  do j = 1,ny
    do i = 1,nx
      zreal(1:nz) = zs(i,j) + ((ztop - zs(i,j)) / ztop) * z(1:nz)
      ks = 1 
      do k = 1,np

        if( lnprs(i,j,1) < lnp(k) )then
          gph(i,j,k) = undef
        elseif( lnprs(i,j,nz) > lnp(k) )then
          gph(i,j,k) = undef
        else
          kk_loop: do kk = ks,nz
            fac = (lnprs(i,j,kk) - lnp(k))*(lnprs(i,j,kk+1) - lnp(k))
            if( fac <= 0.0 )then
              a = (lnprs(i,j,kk)-lnp(k))/(lnprs(i,j,kk)-lnprs(i,j,kk+1))
              gph(i,j,k) = a * zreal(kk+1) + (1.0-a) * zreal(kk)
              ks = kk
              !if( i == 1 .and. j == )then
              if( i == 1 .and. j == 1 )then
                 print'(I4,F10.5,F10.2)',ks,a,gph(i,j,k)
              endif
              exit kk_loop
            endif
          enddo kk_loop
        endif

      enddo
    enddo
  enddo

  deallocate( zreal )

  return
end subroutine conv_prs2gph
