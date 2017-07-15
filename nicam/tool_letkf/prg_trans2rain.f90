!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  prg_trans2rain
!
!  program for precipitation inverse transformation
!
!  created                Apr 2015, Shunji Kotsuki, RIKEN-AICS
!  adopted to NICAM-LETKF May 2015, Shunji Kotsuki, RIKEN-AICS 
!
!-------------------------------------------------------------------------------
program trans2rain

  use mod_adm
  use mpi
  use mod_misc, only :       &
       MISC_get_latlon,      &
       MISC_make_idstr
  use mod_cnst, only :       &
       CNST_setup
  use mod_comm, only :       &
       COMM_setup
  use mod_fio, only : &  ! [add] H.Yashiro 20110826
       FIO_setup
  use mod_grd, only :        &
       GRD_setup,            &
       GRD_x
  use mod_gmtr, only :       &
       GMTR_setup
  
  ! added by kotsuki       
  use fnc_transgauss     
       
  implicit none
 
  REAL(8),PARAMETER :: UNDEF8 = -99.9E+33
  REAL(4),PARAMETER :: rain_elem = 9999.d0

  integer rain_mem
  
  integer g, l, ibv
  integer gall, lall
  
  character(255) fname,gusname,anlname,outname,mapname,aname,oname,ename,gname
 
  REAL(8), ALLOCATABLE :: anl_temp(:)
  REAL(4), ALLOCATABLE :: ave_orgn(:)
  REAL(4), ALLOCATABLE :: anl_invm(:), anl_invo(:), anl_trns(:), anl_inve(:)
  REAL(4), ALLOCATABLE :: gus_invm(:), gus_invo(:), gus_trns(:), gus_inve(:)
  REAL(4), ALLOCATABLE :: obs_invm(:), obs_invo(:), obs_trns(:), obs_inve(:)

  REAL(4), ALLOCATABLE :: ave_invm(:), ave_invo(:), ave_trns(:), ave_inve(:)

  REAL(4), ALLOCATABLE :: gen_orig(:,:), gen_trns(:,:), gen_invo(:,:), gen_invm(:,:)
  REAL(4), ALLOCATABLE ::                aen_trns(:,:), aen_invo(:,:), aen_invm(:,:)
  REAL(4), ALLOCATABLE :: ave_gen_invo(:), ave_gen_invm(:)
  REAL(4), ALLOCATABLE :: ave_aen_invo(:), ave_aen_invm(:)
 
  REAL(8), ALLOCATABLE :: ico_lat(:,:), ico_lon(:,:)
  REAL(8), ALLOCATABLE :: ppcdf_m(:,:), ppcdf_o(:,:), ppcdf_e(:,:)
  REAL(8), ALLOCATABLE :: ppzero_m(:),  ppzero_o(:) , ppzero_e(:)
  REAL(8), ALLOCATABLE :: ppmask_m(:),  ppmask_o(:) , ppmask_e(:)
  REAL(8), ALLOCATABLE :: ppsamp_m(:),  ppsamp_o(:) , ppsamp_e(:)
  REAL(8), ALLOCATABLE :: dfzero_m(:),  dfzero_o(:) , dfzero_e(:)

  REAL(r_size) :: anl_rain_tmp

  !==> nhm_driver.cnf
  integer glevel
  integer rlevel
  integer vlayer
  character(128) rgnmngfname
  

  !--> paralell computing
  integer :: ierr, iread, nread, iprocs, nprocs, irank, myrank
  integer, allocatable :: proc_l(:,:) 

  logical :: ex
  logical :: monit_ensemble = .false.

!========================================================================
  CALL ADM_proc_init(ADM_MULTI_PRC)  
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)  
!========================================================================
  call get_rain2trans()
  call get_gaussvariable()
  call check_transform(myrank)
    
  namelist / ADMPARAM / &
    glevel,                &
    rlevel,                &
    vlayer,                &
    rgnmngfname
  open(1,file='nhm_driver.cnf')
  read(1,nml=ADMPARAM)
  close(1)
  
  gall   =  (2**(glevel-rlevel) + 2)**2
  lall   =  10*4**(rlevel)
    
!  !====> Initial Setting
!  ALLOCATE ( ico_lat(gall,lall), ico_lon(gall,lall) )  
!  write(mapname,'(a,a,i2.2,a,i2.2,a)') trim(grddir),'/icolatlon_gl',glevel,'_rl',rlevel,'.dat'
!  open(1,file=trim(mapname),form='unformatted',access='sequential',,action='read')
!  read(1) ico_lon  ! [deg]
!  read(1) ico_lat  ! [deg]
!  close(1)  

  !===> setting for parallel computing
  if( mod(lall,nprocs).ne.0 ) then
    write(6,*) "lall should be devided by the nprocs"
    stop
  endif
  
  nread = lall / nprocs
  allocate ( proc_l(0:nprocs-1,nread) )
  call set_proc_l(lall,nprocs,nread,proc_l)

  ! region-depending variables 
  ALLOCATE ( anl_temp(gall)                                                 )
  ALLOCATE ( ave_orgn(gall)                                                 )
  ALLOCATE ( anl_invm(gall), anl_invo(gall), anl_inve(gall), anl_trns(gall) )
  ALLOCATE ( gus_invm(gall), gus_invo(gall), gus_inve(gall), gus_trns(gall) )  
  ALLOCATE ( obs_invm(gall), obs_invo(gall), obs_inve(gall), obs_trns(gall) )
  ALLOCATE ( ave_invm(gall), ave_invo(gall), ave_inve(gall), ave_trns(gall) )

  ALLOCATE ( gen_orig(gall,nbv), gen_trns(gall,nbv), gen_invo(gall,nbv), gen_invm(gall,nbv) )
  ALLOCATE (                     aen_trns(gall,nbv), aen_invo(gall,nbv), aen_invm(gall,nbv) )
  ALLOCATE ( ave_gen_invo(gall), ave_gen_invm(gall)                         )
  ALLOCATE ( ave_aen_invo(gall), ave_aen_invm(gall)                         )

  gen_orig(:,:)=sngl(undef) ;   gen_trns(:,:)=sngl(undef) ; gen_invo(:,:)=sngl(undef) ; gen_invm(:,:)=sngl(undef)
                                aen_trns(:,:)=sngl(undef) ; aen_invo(:,:)=sngl(undef) ; aen_invm(:,:)=sngl(undef)
  ave_gen_invo(:)=sngl(undef) ; ave_gen_invm(:)=sngl(undef)
  ave_aen_invo(:)=sngl(undef) ; ave_aen_invm(:)=sngl(undef)


  if ( transform==2 .or. transform==3) then   
    ALLOCATE ( ppcdf_m(gall,0:ndivid), ppcdf_o(gall,0:ndivid), ppcdf_e(gall,0:ndivid) ) 
    ALLOCATE ( ppzero_m(gall),         ppzero_o(gall)        , ppzero_e(gall)         )
    ALLOCATE ( ppmask_m(gall),         ppmask_o(gall)        , ppmask_e(gall)         )
    ALLOCATE ( ppsamp_m(gall),         ppsamp_o(gall)        , ppsamp_e(gall)         )
    ALLOCATE ( dfzero_m(gall),         dfzero_o(gall)        , dfzero_e(gall)         )
  end if

  do iread=1,nread
    l = proc_l(myrank,iread)

    ! read first guess and. obs. precip
    !=====> output of monit transform
    write(outname,'(2a,i6.6,a)') trim(datdir),'/monit_transform'
    call MISC_make_idstr(fname,trim(outname),'rgn',l)
    open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4,action='read')
      read(1,rec=1)  ave_orgn(:)   ! gues [original unit    ; mean]
      read(1,rec=3)  gus_trns(:)   ! gues [transformed unit ; mean]
      read(1,rec=4)  obs_trns(:)   ! obss [transformed unit]
    close(1)

    ! read assimilated precip.
    write(anlname,'(2a,i10.10,a)') trim(anldir),'/',cdate,'/anal_land_me'   
    call MISC_make_idstr(aname,trim(anlname),'rgn',l)       
    open(1,file=trim(aname),form='unformatted',access='direct',recl=gall*8,action='read')
      read(1,rec=2) anl_temp(:)    ! analysis    [transformed unit]
      anl_trns(:) = sngl( anl_temp(:) )
    close(1)
!!!    write(anlname,'(2a,i10.10,a)') trim(anldir),'/',cdate,'/tppn_trns_anal_me'   
!!!    call MISC_make_idstr(aname,trim(anlname),'rgn',l)       
!!!    open(1,file=trim(aname),form='unformatted',access='direct',recl=gall*4,action='read')
!!!      read(1,rec=1) anl_trns(:) 
!!!    close(1)

    !===> read ensemble
    do ibv=1,nbv
      !==> anal transform
      write(anlname,'(2a,i10.10,a,i6.6,a)') trim(anldir),'/',cdate,'/',ibv,'/history_land'   
      call MISC_make_idstr(aname,trim(anlname),'rgn',l)
      INQUIRE( file=trim(aname), exist=ex )
      if( ex ) then
        open(1,file=trim(aname),form='unformatted',access='direct',recl=gall*8,action='read')
          read(1,rec=2) anl_temp(:)     ! analysis    [transformed unit]
          aen_trns(:,ibv) = sngl( anl_temp(:) )
        close(1)
        monit_ensemble = .true.
      end if

      !==> gues transform
      write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',cdate,'/',ibv,'/tppn'
      call MISC_make_idstr(gname,trim(gusname),'rgn',l) 
      INQUIRE( file=trim(gname), exist=ex )
      if( ex ) then
        open(1,file=trim(gname),form='unformatted',access='direct',recl=gall*4,action='read')
          read(1,rec=1) gen_orig(:,ibv) ! gues        [original unit]
        close(1)
        monit_ensemble = .true.
      end if

      !==> gues original
      write(gusname,'(2a,i10.10,a,i6.6,a)') trim(gusdir),'/',cdate,'/',ibv,'/tppn_trns'
      call MISC_make_idstr(gname,trim(gusname),'rgn',l) 
      INQUIRE( file=trim(gname), exist=ex )
      if( ex ) then
        open(1,file=trim(gname),form='unformatted',access='direct',recl=gall*4,action='read')
          read(1,rec=1) gen_trns(:,ibv) ! gues        [original unit]
        close(1)
        monit_ensemble = .true.
      end if
    end do


    !===> precip. inverse transformation
    if ( transform == -1 .or. transform == 0 ) then               ! no transform :: [mm/6hr] -> [mm/6hr]
      anl_invm(:) = anl_trns(:) ; anl_inve(:) = anl_trns(:) ; anl_invo(:) = anl_trns(:)                                     
      gus_invm(:) = gus_trns(:) ; gus_inve(:) = gus_trns(:) ; gus_invo(:) = gus_trns(:)       
      obs_invm(:) = obs_trns(:) ; obs_inve(:) = obs_trns(:) ; obs_invo(:) = obs_trns(:)       
      gen_invo(:,:) = gen_trns(:,:)  ;  gen_invm(:,:) = gen_trns(:,:)
      aen_invo(:,:) = aen_trns(:,:)  ;  aen_invm(:,:) = aen_trns(:,:)
    else if ( transform ==1 ) then                                ! logarithm transform :: [LOG(mm/6hr)] -> [mm/6hr]
      anl_invm(:) = exp( anl_trns(:) ) - rain_alpha ; anl_inve(:) = anl_invm(:) ; anl_invo(:) = anl_invm(:)  
      gus_invm(:) = exp( gus_trns(:) ) - rain_alpha ; gus_inve(:) = gus_invm(:) ; gus_invo(:) = gus_invm(:) 
      obs_invm(:) = exp( obs_trns(:) ) - rain_alpha ; obs_inve(:) = obs_invm(:) ; obs_invo(:) = obs_invm(:) 
      gen_invo(:,:) = exp( gen_trns(:,:) ) - rain_alpha ;  gen_invm(:,:) = exp( gen_trns(:,:) ) - rain_alpha
      aen_invo(:,:) = exp( aen_trns(:,:) ) - rain_alpha ;  aen_invm(:,:) = exp( aen_trns(:,:) ) - rain_alpha
    else if( transform==2 .or. transform==3 ) then                ! gaussian transformation
      call MISC_make_idstr(aname,trim(gus_cdfname),'rgn',l)     
      call MISC_make_idstr(oname,trim(obs_cdfname),'rgn',l)
      call MISC_make_idstr(ename,trim(mem_cdfname),'rgn',l)

      call read_ppcdf(aname,gall,ndivid,ppcdf_m,ppzero_m,ppsamp_m,ppmask_m,dfzero_m) ! read CDF m: MODEL-mean
      call read_ppcdf(oname,gall,ndivid,ppcdf_o,ppzero_o,ppsamp_o,ppmask_o,dfzero_o) ! read CDF o: OBS.
      call read_ppcdf(ename,gall,ndivid,ppcdf_e,ppzero_e,ppsamp_e,ppmask_e,dfzero_e) ! read CDF e: MODEL-Ensemble
      if( l==1 ) write(6,'(2a)') "Trans2Rain::     Use Model CDF ==> ", trim(aname)
      if( l==1 ) write(6,'(2a)') "Trans2Rain:: Not Use Model CDF ==> ", trim(ename)

      if(transform==2 .or. transform==3) then ! Gaussian transformation with medium zero precipitation
        do g=1,gall
          !===> pre process (trans for mean gauss ave-orgn=aveCDF=>ave-trns)
          if( ppmask_m(g)>0.5 ) then
            ave_trns(g) = sngl( pptrans_normal( dble(ave_orgn(g)), ppcdf_m(g,:), ppzero_m(g), dfzero_m(g) ) )
          else
            ave_trns(g) = sngl( undef )
          endif          

          if( ppmask_m(g)>0.5 ) then
            anl_invm(g) = sngl ( ppinverse_normal( dble(anl_trns(g)), ppcdf_m(g,:), ppzero_m(g)  ) )
            gus_invm(g) = sngl ( ppinverse_normal( dble(gus_trns(g)), ppcdf_m(g,:), ppzero_m(g)  ) )
            obs_invm(g) = sngl ( ppinverse_normal( dble(obs_trns(g)), ppcdf_m(g,:), ppzero_m(g)  ) )  
            ave_invm(g) = sngl ( ppinverse_normal( dble(ave_trns(g)), ppcdf_m(g,:), ppzero_m(g)  ) )
            anl_inve(g) = sngl ( ppinverse_normal( dble(anl_trns(g)), ppcdf_e(g,:), ppzero_e(g)  ) )
            gus_inve(g) = sngl ( ppinverse_normal( dble(gus_trns(g)), ppcdf_e(g,:), ppzero_e(g)  ) )
            obs_inve(g) = sngl ( ppinverse_normal( dble(obs_trns(g)), ppcdf_e(g,:), ppzero_e(g)  ) )  
            ave_inve(g) = sngl ( ppinverse_normal( dble(ave_trns(g)), ppcdf_e(g,:), ppzero_e(g)  ) )
            do ibv=1,nbv
              gen_invm(g,ibv) = sngl( ppinverse_normal( dble(gen_trns(g,ibv)), ppcdf_e(g,:), ppzero_e(g)  ) )
              aen_invm(g,ibv) = sngl( ppinverse_normal( dble(aen_trns(g,ibv)), ppcdf_e(g,:), ppzero_e(g)  ) )              
            end do
          else
            anl_invm(g)   = sngl ( undef )
            gus_invm(g)   = sngl ( undef )
            obs_invm(g)   = sngl ( undef )
            ave_invm(g)   = sngl ( undef )
            anl_inve(g)   = sngl ( undef )
            gus_inve(g)   = sngl ( undef )
            obs_inve(g)   = sngl ( undef )
            ave_inve(g)   = sngl ( undef )
            gen_invm(g,:) = sngl( undef )
            aen_invm(g,:) = sngl( undef )
          end if

          if( ppmask_o(g)>0.5 ) then
            anl_invo(g) = sngl ( ppinverse_normal( dble(anl_trns(g)), ppcdf_o(g,:), ppzero_o(g)  ) )
            gus_invo(g) = sngl ( ppinverse_normal( dble(gus_trns(g)), ppcdf_o(g,:), ppzero_o(g)  ) )
            obs_invo(g) = sngl ( ppinverse_normal( dble(obs_trns(g)), ppcdf_o(g,:), ppzero_o(g)  ) )
            ave_invo(g) = sngl ( ppinverse_normal( dble(ave_trns(g)), ppcdf_o(g,:), ppzero_o(g)  ) )

            do ibv=1,nbv
              gen_invo(g,ibv) = sngl ( ppinverse_normal( dble(gen_trns(g,ibv)), ppcdf_o(g,:), ppzero_o(g)  ) )
              aen_invo(g,ibv) = sngl ( ppinverse_normal( dble(aen_trns(g,ibv)), ppcdf_o(g,:), ppzero_o(g)  ) )
            end do
          else
            anl_invo(g)   = sngl ( undef )
            gus_invo(g)   = sngl ( undef )
            obs_invo(g)   = sngl ( undef )
            ave_invo(g)   = sngl ( undef )
            gen_invo(g,:) = sngl ( undef )
            aen_invo(g,:) = sngl ( undef )
          end if      
        end do ! gall
      else 
        write(6,*) "Errror, Not Developped"
        stop
      end if
    else    
      write(6,*) "no such option, transform ==", transform
      stop      
    end if                

    do g=1,gall
      ave_gen_invo(g) = sum( gen_invo(g,:) ) / real(nbv)
      ave_gen_invm(g) = sum( gen_invm(g,:) ) / real(nbv)
      ave_aen_invo(g) = sum( aen_invo(g,:) ) / real(nbv)
      ave_aen_invm(g) = sum( aen_invm(g,:) ) / real(nbv)
    end do
    
    !===> write ensemble mean for gaussian transformation
    write(outname,'(2a,i6.6,a)') trim(datdir),'/monit_inverse'
    call MISC_make_idstr(fname,trim(outname),'rgn',l)
    open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4)
      write(1,rec=1 ) gus_trns(:)   ! [transformed unit]     
      write(1,rec=2 ) anl_trns(:)   ! [transformed unit]
      write(1,rec=3 ) obs_trns(:)   ! [transformed unit]
      write(1,rec=4 ) gus_invm(:)   ! [model-like inversed; mm/6hr]      
      write(1,rec=5 ) anl_invm(:)   ! [model-like inversed; mm/6hr]  
      write(1,rec=6 ) obs_invm(:)   ! [model-like inversed; mm/6hr]  
      write(1,rec=7 ) gus_invo(:)   ! [obs.-like  inversed; mm/6hr]       
      write(1,rec=8 ) anl_invo(:)   ! [obs.-like  inversed; mm/6hr]  
      write(1,rec=9 ) obs_invo(:)   ! [obs.-like  inversed; mm/6hr]

      write(1,rec=10) ave_gen_invm (:) ! ens-trns --> ens-invm --> average [mm/6hr]
      write(1,rec=11) ave_gen_invo (:) ! ens-trns --> ens-invo --> average [mm/6hr]
      write(1,rec=12) ave_aen_invm (:) ! ens-trns --> ens-invm --> average [mm/6hr]
      write(1,rec=13) ave_aen_invo (:) ! ens-trns --> ens-invo --> average [mm/6hr]
    close(1)             
    
    !===> write ensemble mean for gaussian transformation
    write(outname,'(2a,i6.6,a)') trim(datdir),'/monit_ensmean'
    call MISC_make_idstr(fname,trim(outname),'rgn',l)
    open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4)
      write(1,rec=1) gus_invm(:)   ! mm/6h 
      write(1,rec=2) gus_inve(:)   ! mm/6h
      write(1,rec=3) gus_invo(:)   ! mm/6h
      write(1,rec=4) ave_invm(:)   ! mm/6h 
      write(1,rec=5) ave_inve(:)   ! mm/6h 
      write(1,rec=6) ave_invo(:)   ! mm/6h
      write(1,rec=7) gus_trns(:)   ! gauss
      write(1,rec=8) ave_trns(:)   ! gauss
    close(1)

    !===> write ensemble
    if( monit_ensemble ) then
      write(outname,'(2a,i6.6,a)') trim(datdir),'/monit_ensemble'
      call MISC_make_idstr(fname,trim(outname),'rgn',l)
      open(1,file=trim(fname),form='unformatted',access='direct',recl=gall*4)
        do ibv=1,nbv
          write(1,rec=(ibv-1)*7+1) gen_trns(:,ibv)  ! gauss
          write(1,rec=(ibv-1)*7+2) aen_trns(:,ibv)  ! gauss
          write(1,rec=(ibv-1)*7+3) gen_orig(:,ibv)  ! mm/6hr
          write(1,rec=(ibv-1)*7+4) gen_invo(:,ibv)  ! mm/6hr
          write(1,rec=(ibv-1)*7+5) gen_invm(:,ibv)  ! mm/6hr
          write(1,rec=(ibv-1)*7+6) aen_invo(:,ibv)  ! mm/6hr
          write(1,rec=(ibv-1)*7+7) aen_invm(:,ibv)  ! mm/6hr
        end do
      close(1)    
    end if

    !==> obs
    ! obs_invm :: obs  =(obsCDF)=> obs-gauss                                         =(aveCDF)=> obs-invGT(mem)
    ! obs_inve :: obs  =(obsCDF)=> obs-gauss                                         =(ensCDF)=> obs-invGT(ens)
    ! obs_invo :: obs  =(obsCDF)=> obs-gauss                                         =(obsCDF)=> obs-invGT(obs)
    !==> analysis 
    ! anl_invm :: gues =(ensCDF)=> mem-gauss =(AVE)=> mean-gauss =(DAS)=> anal-gauss =(aveCDF)=>anal-invGT(mem)
    ! anl_inve :: gues =(ensCDF)=> mem-gauss =(AVE)=> mean-gauss =(DAS)=> anal-gauss =(ensCDF)=>anal-invGT(ens)
    ! anl_invo :: gues =(ensCDF)=> mem-gauss =(AVE)=> mean-gauss =(DAS)=> anal-gauss =(obsCDF)=>anal-invGT(obs)
    !==> guess1
    ! gus_invm :: gues =(ensCDF)=> mem-gauss =(AVE)=> mean-gauss                     =(aveCDF)=>gues-invGT(mem)
    ! gus_inve :: gues =(ensCDF)=> mem-gauss =(AVE)=> mean-gauss                     =(ensCDF)=>gues-invGT(ens)
    ! gus_invo :: gues =(ensCDF)=> mem-gauss =(AVE)=> mean-gauss                     =(obsCDF)=>gues-invGT(obs)
    !==> guess2
    ! ave_invm :: gues =(AVE)=> mean-gues  =(aveCDF)=> "mean-gauss"                  =(aveCDF)=> gave-invGT(mem)     
    ! ave_inve :: gues =(AVE)=> mean-gues  =(aveCDF)=> "mean-gauss"                  =(ensCDF)=> gave-invGT(ens)    
    ! ave_invo :: gues =(AVE)=> mean-gues  =(aveCDF)=> "mean-gauss"                  =(obsCDF)=> gave-invGT(obs)


  end do ! iread to paralellize region loop

  if (myrank==0) write(6,*) 'Finish Program trans2rain'
  call MPI_Finalize(ierr)
end program
