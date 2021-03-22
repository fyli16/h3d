module initialize
  use parameter_mod
  implicit none
  integer*4:: input_error

  contains
  subroutine read_input()
    namelist /datum/ &
    ! global simulation info
    tmax, t_begin, t_end, dtwci, dt, restart, &
    restrt_write, quota, MPI_IO_format, &
    ! simulation domain
    nx, ny, nz, xmax, ymax, zmax, npx, npy, npz, &
    nodey, nodez, &
    xaa, xbb, nax, nbx, yaa, ybb, nay, nby, zaa, zbb, naz, nbz, &
    uniform_loading_in_logical_grid, &
    buffer_zone, moat_zone, profile_power, &
    ! field solver
    n_subcycles, nskipx, nskipy, nskipz, iterb, testorbt, norbskip, &
    ! plasma setup
    nspec, qspec, wspec, frac, denmin, wpiwci, btspec, bete, &
    ieta, resis, netax, netay, netaz, etamin, etamax, eta_par, &
    anisot, gama, ave1, ave2, phib, smoothing, smooth_coef, &
    ! init waves
    dB_B0, num_cycles, &
    ! diagnostic control
    nprint, nwrtdata, nwrtrestart, nwrtparticle, &
    xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r, &
    ! others
    Yee, global, harris, fxsho, nxcel, & 
    rcorr, ishape, teti, setup_mesh, post_process
      
    time_elapsed=0.; time_begin_array=0; time_end_array=0
    buffer_zone=0.  ! set to 0 anyway despite contained in input
    notime=1 ! notime=0 will output detailed timing
    !tracking_binary=.false.

    ! MPI initialization
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

    ! get_simulation_id
    ! if (myid==0) call get_sim_id(sim_id)

    ! time stamp
    !  call date_and_time(values=wall_clock_begin)
    initial_time = MPI_Wtime()

    ! get the i/o directory names from the environment variable
    call get_environment_variable1(data_directory, len(data_directory))
    data_directory = trim(adjustl(data_directory))//'/'
    call get_environment_variable2(restart_directory, len(restart_directory))
    restart_directory = trim(adjustl(restart_directory))//'/'
    restart_index_suffix(1) = '.1'
    restart_index_suffix(2) = '.2'

    my_short_int = myid
    call integer_to_character(myid_char,len(myid_char), my_short_int)
    if (myid_char == '') myid_char='0'

    ! read input deck
    if (myid == 0) then
        open(5, file='input.f90', form='formatted', status='old')
        read(5, nml=datum, iostat=input_error)
        write(6, datum)
    endif
    iwt=0; t_stopped=0. ! set default values

    ! Broadcast info to all ranks
    ! global sim. info
    call MPI_BCAST(tmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(t_begin                ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(t_end                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dtwci                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dt                     ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restart                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restrt_write           ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(quota                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(MPI_IO_format          ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    ! sim. domain
    call MPI_BCAST(nx                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ny                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nz                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ymax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(npx                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(npy                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(npz                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nodey                  ,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nodez                  ,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nax                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nbx                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(yaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nay                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nby                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(naz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nbz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(uniform_loading_in_logical_grid,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(buffer_zone            ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(moat_zone              ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(profile_power          ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    ! field solver
    call MPI_BCAST(n_subcycles            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipx                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipy                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipz                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(iterb                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(testorbt               ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(norbskip               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! plasma setup
    call MPI_BCAST(nspec                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(qspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(frac                   ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(denmin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wpiwci                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(btspec                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(bete                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ieta                   ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(resis                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netax                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netay                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netaz                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamax                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(eta_par                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(anisot                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(gama                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ave1                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ave2                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(phib                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smoothing              ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smooth_coef            ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! init waves
    call MPI_BCAST(dB_B0                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(num_cycles             ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! diagnostic control
    call MPI_BCAST(nprint                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nwrtdata               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nwrtrestart            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nwrtparticle           ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! others 
    call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(global                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(fxsho                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nxcel                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(rcorr                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ishape                 ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(teti                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(setup_mesh             ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(post_process           ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

    ! The unit of dt is 1/wci in input file, but converted to 1/wpi here
    dt = dtwci * wpiwci

    if (myid==0) then
      if (restart) then
        write(6,*) "*** RUN IS RESTARTED FROM "//trim(adjustl(restart_directory))
      else
        write(6,*) "*** NEW RUN "
      endif
    endif 
  end subroutine read_input

end module initialize