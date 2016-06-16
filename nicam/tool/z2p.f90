  program main

  integer :: nx
  integer :: ny
  integer :: nz
  integer :: np

  real(4), allocatable :: xgrd(:)
  real(4), allocatable :: ygrd(:)
  real(4), allocatable :: zgrd(:)
  real(4), allocatable :: cs(:)
  real(4), allocatable :: pres(:)
  real(4), allocatable :: pgrd(:)
  real(4), allocatable :: datz(:)
  real(4), allocatable :: datp(:)
  real(4), allocatable :: ttt(:,:)
  real(4), allocatable :: ppp(:,:)
  real(4), allocatable :: uuu(:,:)
  real(4), allocatable :: ddd(:,:,:,:)
  real(4), allocatable :: output_data(:,:,:,:)
  real(4), allocatable :: output_data1(:,:,:,:)

  
  character(128) :: ifile_xgrd
  character(128) :: ifile_ygrd
  character(128) :: ifile_zgrd
  character(128) :: ifile_pgrd
  character(128) :: ifile_fname
  character(128) :: ofile_fname
  character(4)   :: yy
  character(2)   :: mm

  real(8), parameter :: pi = 4.0 * atan(1.0)
  real(8), parameter :: rr = 287.04
  real(8), parameter :: gr = 9.80614

  namelist / z2p_cnf / &
    nx,ny,nz,np,nobs,  &
    g_level,           &
    ifile_xgrd,        &
    ifile_ygrd,        &
    ifile_zgrd,        &
    ifile_pgrd,        &
    ifile_fname,       &
    ofile_fname

!  include 'mpif.h'

!  call mpi_init(ierr)
!  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)  

!  write(*,*) 'myrank=',myrank

  open(unit=10,file='z2p.cnf', form='formatted', status='old', iostat=ierr)
  if(ierr /= 0 ) then
    write(*,*) 'Cannot open z2p.cnf!'
    stop
  end if 

  read(10,nml=z2p_cnf)
  write(*,z2p_cnf)
  close(10)
  write(*,*) nx,ny,nz,np

  allocate( xgrd(nx) )
  allocate( ygrd(ny) ) 
  allocate( zgrd(nz) )

  allocate( cs(ny) )
  allocate( pgrd(np) )
  allocate( pres(nz) )
  allocate( datz(nz) )
  allocate( datp(np) )

  allocate( ddd(nx,ny,nz,7) )
  allocate( uuu(nx,nz) )
  allocate( ttt(nx,nz) )
  allocate( ppp(nx,nz) )
  allocate( output_data(nx,0:ny+1,np,7) )
  allocate( output_data1(nx,0:ny+1,np,7) )

  open(unit=11,file=ifile_xgrd)
  open(unit=12,file=ifile_ygrd)
  open(unit=13,file=ifile_zgrd)
  open(unit=14,file=ifile_pgrd)

  read(11,*) xgrd
  read(12,*) ygrd
  read(13,*) zgrd
  read(14,*) pgrd
  close(11)
  close(12)
  close(13)
  close(14)

  write(*,*) 'zgrd'
  write(*,*) zgrd
  write(*,*) 'pgrd'
  write(*,*) pgrd

!  if(myrank == 0) then
!    write(*,*) ' Grid information'
!    write(*,*) ' Zonal direction'
!    write(*,*) xgrd
!    write(*,*) ' Meridional direction'
!    write(*,*) ygrd
!    write(*,*) ' Vertical direction'
!    write(*,*) zgrd
!  end if

  cs(:) = cos( ygrd(:) * pi / 180.0 )  

  irec1=4*nx*ny*nz*7
  open(21,file=trim(ifile_fname),form='unformatted',access='direct',&
       recl=irec1)
  write(*,*) 'recl=',irec1
  read(21,rec=1) ddd
  close(21)
  
  do n = 1, 7
    write(*,*) n, minval(ddd(:,:,:,n)), maxval(ddd(:,:,:,n))
  end do

  do iy=1,ny
!!!!  T
    do ix=1,nx
      do iz=nz,1,-1
        if(abs(ddd(ix,iy,iz,2)).gt.1.0E20)  then
          ddz=zgrd(iz+1)-zgrd(iz)
          ddd(ix,iy,iz,2) = ddd(ix,iy,iz+1,2) + 0.006*ddz
        endif
      end do
    end do
!!!! P Z T
    do ix=1,nx
      do iz=nz,1,-1
        if(abs(ddd(ix,iy,iz,1)).gt.1.0E20) then
          ddz=zgrd(iz+1)-zgrd(iz)
          sch=rr*0.5*(ddd(ix,iy,iz+1,2)+ddd(ix,iy,iz,2))/gr
          ddd(ix,iy,iz,1)=ddd(ix,iy,iz+1,1)*exp(ddz/sch)
        endif
      end do
!      write(*,*) 'Z'
      do iz=1,nz
        datz(iz)=zgrd(iz)
        pres(iz)=ddd(ix,iy,iz,1)/100.0
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
        output_data(ix,iy,ip,5)=datp(ip)
      end do
!
!      write(*,*) 'T'
      do iz=1,nz
        pres(iz)=ddd(ix,iy,iz,1)/100.0
        datz(iz)=ddd(ix,iy,iz,2)
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
          output_data(ix,iy,ip,4)=datp(ip)
      end do
    end do
!
!!!!  u
!      write(*,*) 'U'
    do ix=1,nx
      do iz=1,nz
        if(abs(ddd(ix,iy,iz,3)).gt.1.0E20) ddd(ix,iy,iz,3) = 0.0
        pres(iz)=ddd(ix,iy,iz,1)/100.0
        datz(iz)=ddd(ix,iy,iz,3)
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
        output_data(ix,iy,ip,1)=datp(ip)
      end do
    end do

!!!!  v
!      write(*,*) 'V'
    do ix=1,nx
      do iz=1,nz
        if(abs(ddd(ix,iy,iz,4)).gt.1.0E20) ddd(ix,iy,iz,4) = 0.0
        pres(iz)=ddd(ix,iy,iz,1)/100.0
        datz(iz)=ddd(ix,iy,iz,4)
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
        output_data(ix,iy,ip,2)=datp(ip)
      end do
    end do
!
!!!!  w  
!      write(*,*) 'W'
    do ix=1,nx
      do iz=1,nz
        if(abs(ddd(ix,iy,iz,5)).gt.1.0E20) ddd(ix,iy,iz,5) = 0.0
        pres(iz)=ddd(ix,iy,iz,1)/100.0
        datz(iz)=ddd(ix,iy,iz,5)
        !datz(iz)=-gr*ppp(ix,iz)/(rr*ttt(ix,iz))*uuu(ix,iz)
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
        output_data(ix,iy,ip,3)=datp(ip)
      end do
    end do
!
!!!  qv
!      write(*,*) 'Q'
    do ix=1,nx
      do iz=1,nz
        if(abs(ddd(ix,iy,iz,6)).gt.1.0E20) ddd(ix,iy,iz,6) = 0.0
        pres(iz)=ddd(ix,iy,iz,1)/100.0
        datz(iz)=ddd(ix,iy,iz,6)
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
        output_data(ix,iy,ip,6)=datp(ip)
      end do
    end do
!
!!!  qc
!      write(*,*) 'Q'
    do ix=1,nx
      do iz=1,nz
        if(abs(ddd(ix,iy,iz,7)).gt.1.0E20) ddd(ix,iy,iz,7) = 0.0
        pres(iz)=ddd(ix,iy,iz,1)/100.0
        datz(iz)=ddd(ix,iy,iz,7)
      end do
      call z2p(nz,pres,datz,np,pgrd,datp)
      do ip=1,np
        output_data(ix,iy,ip,7)=datp(ip)
      end do
    end do

  end do

  do n = 1, 7
    write(*,*) n, minval(output_data(:,1:ny,:,n)), maxval(output_data(:,1:ny,:,n))
  end do

  do ix = 1, nx
    output_data(ix,0,:,:)=(output_data(ix,1,:,:)+output_data(ix+nx/2,1,:,:))*0.5
    output_data(ix,ny+1,:,:)=(output_data(ix,ny,:,:)+output_data(ix+nx/2,ny,:,:))*0.5
  end do
  do n = 1, 7
    do ip = 1, np
      output_data(:, 0,ip,n)=sum(output_data(:, 0,ip,n))/float(nx)
      output_data(:,ny,ip,n)=sum(output_data(:,ny,ip,n))/float(nx)
    end do
  end do

  do n = 1, 7
    do ip = 1, np
      do iy = 0, ny+1
        do ix = 2, nx
          output_data1(ix,iy,ip,n)=0.5d0*(output_data(ix-1,iy,ip,n)+output_data(ix,iy,ip,n))
        end do
        output_data1(1,iy,ip,n)=0.5d0*(output_data(nx,iy,ip,n)+output_data(1,iy,ip,n))
      end do
    end do
  end do

  do n = 1, 7
    write(*,*) n, minval(output_data1(:,:,:,n)), maxval(output_data1(:,:,:,n))
  end do

  irec1=4*nx*(ny+2)*np*7
  open(31,file=trim(ofile_fname),form='unformatted',access='direct',&
       recl=irec1)
  write(*,*) 'recl=',irec1
  write(31,rec=1) output_data1
  close(31)

end program main



      subroutine z2p (n,p1,u1,m,p2,u2)
      !parameter (nz=40,np=26)
      dimension p1(n),u1(n),p2(m),u2(m)
      dimension z1(n),z2(m),cc(n),ww(n)
      !dimension z1(nz),z2(np),cc(nz),ww(nz)
      h0=8000.0
      p0=1000.0
      do 10 i=1,n
   10 z1(i)=h0*alog(p0/p1(i))
      do 20 i=1,m
   20 z2(i)=h0*alog(p0/p2(i))
!      write(*,*) "z1"
!      write(*,*) z1(1:n)
!      write(*,*) "u1"
!      write(*,*) u1(1:n)
      call spline (n,z1,u1,cc,ww)
      call splint (n,z1,u1,cc,m,z2,u2)
!      write(*,*) "u2"
!      write(*,*) u2(1:m)
      return
      end

      SUBROUTINE SPLINE(N,X,Y,CC,WW)
      DIMENSION X(N),Y(N),CC(N),WW(N)
      CC(1)=0.
      WW(1)=0.
!
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
      DIMENSION AX(N),AY(N),CC(N),BX(M),BY(M)
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

