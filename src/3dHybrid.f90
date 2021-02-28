!***************************************************************************
!                                                                          *
!                                 VERSION 6.0                              *
!                           YURI'S NONUNIFORM MESH                         *
!                           3D IMPLEMENTATION ONLY                         *
!                      UNIFORM LOADING IN PHYSICAL SPACE                   *
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       *
!                                                                          *
!***************************************************************************

    program H3D 
      use parameter_mod
      use functions_f90
      use MESH2D
      implicit none
      integer*8 i,irecnum,ixe,iye,ize,j,jbt,jet,k,kbt,ket,nplmax6, &
        nwrtparticle,nwrtrestart,nyl,nzl,numvars
      double precision rnorm,pifac
      integer*4:: time_begin(8),time_end(8),input_error,is
      integer*8  itstart, itfinish
      double precision :: clock_time_re1
      double precision, dimension(:,:,:), allocatable:: uniform_mesh      
      !VR : allocating a global mesh can not work on large runs with small
      !VR : amount of memory per rank
      !VR : double precision, dimension(:,:,:), allocatable:: nonuniform_mesh_global
      character (len=240):: filename, filename2
      character(len=128):: tmpfname
      character(len=27) :: sim_id
      integer:: iwrite
      integer:: ierr2, eStrLen
      character(len=1024):: eStr
      external get_environment_variable

      ! namelist to provide from input
      namelist /datum/ max_sim_time, t_begin, t_end, dtwci, dt, restart, &
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
      if (myid==0) call get_sim_id(sim_id)

      ! time stamp  
      ! call date_and_time(values=wall_clock_begin)
      initial_time = MPI_Wtime()
      if (myid == 0) write(6,*) " H3D is starting"

      ! get the i/o data directory name from the environment variable DATA_DIRECTORY
      if (myid==0) then
        call get_environment_variable1(data_directory,len(data_directory))
        data_directory=trim(adjustl(data_directory))//'/'
        call get_environment_variable2(restart_directory,len(restart_directory))
        restart_directory=trim(adjustl(restart_directory))//'/'
        write(6,*)
        write(6,*) "I/O DATA_DIRECTORY = ",trim(adjustl(data_directory))
        write(6,*) "RESTART_DIRECTORY  = ",trim(adjustl(restart_directory))
        write(6,*)
      endif
      call MPI_BCAST(data_directory,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(restart_directory,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      restart_index_suffix(1)='.1'
      restart_index_suffix(2)='.2'

      my_short_int=myid
      call integer_to_character(myid_char,len(myid_char),my_short_int)
      if (myid_char == '') myid_char='0'

      ! set default values
      iwt=0; nskipx=1; nskipy=1; nskipz=1; testorbt=.false.
      pi=acos(-1.d+00); frac=0.d+00; t_stopped=0.
      
      ! read input deck
      if (myid == 0) then
        !  open (5,file=trim(adjustl(data_directory))//'input.f90',form='formatted',status='old')
        open (5,file='./input.f90',form='formatted',status='old')
        read(5,nml=datum,iostat=input_error)
        write(6, datum)
      endif
 
      ! hxv - 12/02/2008 - Automatic restart
      inquire(file=trim(adjustl(restart_directory))//'restart_index.dat',exist=restart)

      ! Broadcast variables 
      ! global sim. info
      call MPI_BCAST(max_sim_time           ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
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

      ! In input.f90, dt is in units of 1/wci, 
      !   now inside the code it is converted to units of 1/wpi
      dt = dtwci * wpiwci

      ! hxv - 12/02/2008 -Automatic restart
      if (restart .and. myid == 0) then
        write(6,*) " "
        write(6,*) " "
        write(6,*) " RUN IS RESTARTED FROM "//trim(adjustl(restart_directory))
        write(6,*) " "
        write(6,*) " "
      else
        if (myid == 0) then 
          call create_id_file(sim_id)
          write(6,*) " "
          write(6,*) " NEW RUN "
          write(6,*) " "
          write(6,*) " "
        endif
      endif

      ! field subcycling
      ! n_subcycles=max(n_subcycles,1_8)

      ! set MPI Cartesian geometry, define stride vector types, obtain new
      ! ID for the processors, perform 2D decomposition of the
      ! computational mesh, and find nearest neighbors (in y and z
      ! directions)
      if (myid==0) then
        write(6,*) " Number of processors available = ",NUMPROCS
      endif

      ! fyli: because of the manual specification below, this selection block is actually not functioning and thus is annotated. 
      ! if (nz == 1.and.ny == 1) then  
      !   ndim=1
      !   dims(1)=1
      !   dims(2)=1
      ! else if (nz == 1) then
      !   ndim=1
      !   dims(2)=1
      ! else
      !   ndim=2
      ! endif
      ! manually specify decomposition (along y, z only; x direction is not decomposed) 
      ndim=2
      dims(1)=nodey
      dims(2)=nodez

      ! if (ndim /= 2) then
      !    if (myid==0) then
      !       print *,"*************************************************************************"
      !       print *," ERROR: FIELD SOLVER HAS NOT BEEN MODIFIED FOR PERIODIC B.C. in 1D and 2D"
      !       print *,"                            H3D TERMINATING                              "
      !       print *,"*************************************************************************"
      !    endif
      !    call MPI_FINALIZE(IERR)
      !    STOP         
      ! endif

      ! Boundary conditions
      PERIODS = .TRUE. 
      REORDER = .TRUE.
      call MPI_DIMS_CREATE(NUMPROCS,NDIM,DIMS,IERR)

      ! npy and npz are changed here!!!!
      ! now npy means number of particles in each core along y
      npy=npy/dims(1)
      npz=npz/dims(2)

      if (myid == 0) then
        do i=1,ndim
          write(6,*) " DIMENSION = ",I," DIMS = ",DIMS(I)
        enddo
      endif
      call MPI_CART_CREATE(MPI_COMM_WORLD,NDIM,DIMS,PERIODS,REORDER,COMM2D,IERR)
      call MPI_COMM_RANK(COMM2D,MYID,IERR)
      call MPI_CART_GET(COMM2D,NDIM,DIMS,PERIODS,COORDS,IERR)
      call MPE_DECOMP1D(NZ,DIMS(2),COORDS(2),KB,KE)
      call MPE_DECOMP1D(NY,DIMS(1),COORDS(1),JB,JE)
      ! print domain decomposition info
      ! write(*,*)'myid=',myid,'jb,je',jb,je,'kb,ke=',kb,ke
 
      nspecm = nspec
      nxmax  = nx+2
      nymax  = ny+2
      nzmax  = nz+2
      nylmax = je-jb+1
      nzlmax = ke-kb+1
      call set_parameters(NUMPROCS)
      myid_stop(myid)=0 
      do is = 1 , nspecm
         qleft(is)=0
         qrite(is)=0
      enddo

      if (myid==0) then
        write(6,*) " LOCAL ARRAY SIZE IN Y-DIRECTION = ",JE-JB+1
        write(6,*) " LOCAL ARRAY SIZE IN Z-DIRECTION = ",KE-KB+1
      endif

      if (nzlmax < ke-kb+1) then
          print *,'myid = ',myid,' nzlmax lt ke-kb+1'
          print *,'myid = ',myid,' nzlmax,ke,kb= ',nzlmax,ke,kb
          myid_stop(myid) = 1
      endif
      if (nylmax < je-jb+1) then
          print *,'myid = ',myid,' nylmax lt je-jb+1'
          print *,'myid = ',myid,' nylmax,je,jb= ',nylmax,je,jb
          myid_stop(myid) = 1
      endif
      do i = 0, npes-1
        ! if (myid.eq.i) then
        !    write(*,*) "Node no:",i,"myid_stop=",MYID_STOP(I)
        ! endif
         i_i = i
         CALL MPI_BCAST(MYID_STOP(i),1,MPI_INTEGER8,i_i,MPI_COMM_WORLD, IERR)
      enddo
      do i = 0,npes-1
         if (myid_stop(i).ne.0) then
            call MPI_FINALIZE(IERR)
            write(*,*)"TEST HERE"
            write(*,*)i, myid_stop(i)
            STOP
         endif
      enddo

      ! Use CART_SHIFT to determine processor to immediate left
      ! (NBRLEFT) and right (NBRRITE) of processor MYID
      ! Since code is aperiodic in z, need to manually set the
      ! left boundary for processor 0 and right boundary for npes-1
      if (ndim == 2) then
        call MPI_CART_SHIFT(COMM2D,0,1,NBRLEFT,NBRRITE,IERR)
        call MPI_CART_SHIFT(COMM2D,1,1,NBRBOT ,NBRTOP ,IERR)
      else if (ndim == 1) then
        call MPI_CART_SHIFT(COMM2D,0,1,NBRLEFT,NBRRITE,IERR)
        NBRTOP=MYID
        NBRBOT=MYID
      else if (ndim == 0) then
        NBRLEFT=MYID
        NBRRITE=MYID
        NBRTOP =MYID
        NBRBOT =MYID
      endif
      call MPI_SENDRECV(NBRTOP    ,1,MPI_INTEGER ,NBRRITE,0,&
                        NBRLEFTTOP,1,MPI_INTEGER ,NBRLEFT,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(NBRTOP    ,1,MPI_INTEGER ,NBRLEFT,0,&
                        NBRRITETOP,1,MPI_INTEGER ,NBRRITE,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(NBRBOT    ,1,MPI_INTEGER ,NBRRITE,0,&
                        NBRLEFTBOT,1,MPI_INTEGER ,NBRLEFT,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(NBRBOT    ,1,MPI_INTEGER ,NBRLEFT,0,&
                        NBRRITEBOT,1,MPI_INTEGER ,NBRRITE,0,&
                        mpi_comm_world,status,ierr) 
 
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        isendid(1)=1
      else
        isendid(1)=0
      endif
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,1)=nbrrite
        irecvid(2,1)=-1
        irecvid(3,1)=nbrleft
        irecvid(4,1)=-1
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,1)=-1
        irecvid(2,1)=nbrtop
        irecvid(3,1)=-1
        irecvid(4,1)=nbrbot
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,1)=nbrritetop
        irecvid(2,1)=nbrlefttop
        irecvid(3,1)=nbrleftbot
        irecvid(4,1)=nbrritebot
      endif
 
 
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        isendid(2)=1
      else
        isendid(2)=0
      endif
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,2)=nbrrite
        irecvid(2,2)=-1
        irecvid(3,2)=nbrleft
        irecvid(4,2)=-1
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,2)=-1
        irecvid(2,2)=nbrtop
        irecvid(3,2)=-1
        irecvid(4,2)=nbrbot
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,2)=nbrritetop
        irecvid(2,2)=nbrlefttop
        irecvid(3,2)=nbrleftbot
        irecvid(4,2)=nbrritebot
      endif
 
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        isendid(3)=1
      else
        isendid(3)=0
      endif
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,3)=nbrrite
        irecvid(2,3)=-1
        irecvid(3,3)=nbrleft
        irecvid(4,3)=-1
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,3)=-1
        irecvid(2,3)=nbrtop
        irecvid(3,3)=-1
        irecvid(4,3)=nbrbot
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,3)=nbrritetop
        irecvid(2,3)=nbrlefttop
        irecvid(3,3)=nbrleftbot
        irecvid(4,3)=nbrritebot
      endif
 
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        isendid(4)=1
      else
        isendid(4)=0
      endif
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,4)=nbrrite
        irecvid(2,4)=-1
        irecvid(3,4)=nbrleft
        irecvid(4,4)=-1
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,4)=-1
        irecvid(2,4)=nbrtop
        irecvid(3,4)=-1
        irecvid(4,4)=nbrbot
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,4)=nbrritetop
        irecvid(2,4)=nbrlefttop
        irecvid(3,4)=nbrleftbot
        irecvid(4,4)=nbrritebot
      endif
      nzl=nzlmax
      nyl=nylmax

      ! if (myid == 0) then
      !    jbglobal(myid)=jb
      !    jeglobal(myid)=je
      !    kbglobal(myid)=kb
      !    keglobal(myid)=ke
      !    do ipe=numprocs-1,1,-1
      !       call MPI_IRECV(jbglobal(ipe),1,MPI_INTEGER8,IPE,0,MPI_COMM_WORLD,req(1),IERR)
      !       call MPI_IRECV(jeglobal(ipe),1,MPI_INTEGER8,IPE,1,MPI_COMM_WORLD,req(2),IERR)
      !       call MPI_IRECV(kbglobal(ipe),1,MPI_INTEGER8,IPE,2,MPI_COMM_WORLD,req(3),IERR)
      !       call MPI_IRECV(keglobal(ipe),1,MPI_INTEGER8,IPE,3,MPI_COMM_WORLD,req(4),IERR)
      !       call MPI_WAITALL(4,req,status_array,IERR)
      !    enddo
      ! else
      !    call MPI_ISEND(jb           ,1,MPI_INTEGER8,0,0,MPI_COMM_WORLD,req(1),IERR)
      !    call MPI_ISEND(je           ,1,MPI_INTEGER8,0,1,MPI_COMM_WORLD,req(2),IERR)
      !    call MPI_ISEND(kb           ,1,MPI_INTEGER8,0,2,MPI_COMM_WORLD,req(3),IERR)
      !    call MPI_ISEND(ke           ,1,MPI_INTEGER8,0,3,MPI_COMM_WORLD,req(4),IERR)
      !    call MPI_WAITALL(4,req,status_array,IERR)
      ! endif
      ! call MPI_BCAST(JBGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      ! call MPI_BCAST(JEGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      ! call MPI_BCAST(KBGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      ! call MPI_BCAST(KEGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

      !VR: this is much simpler
      call MPI_ALLGATHER(jb,1,MPI_INTEGER8,jbglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
      call MPI_ALLGATHER(je,1,MPI_INTEGER8,jeglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
      call MPI_ALLGATHER(kb,1,MPI_INTEGER8,kbglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
      call MPI_ALLGATHER(ke,1,MPI_INTEGER8,keglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)

      ! VR: what is going on here?
      ! if (myid.ne.0) then
      !   do k=kb,ke
      !     do j=1,je-jb+1
      !       jvec(j)=myid
      !     enddo
      !     i_length=je-jb+1
      !     call MPI_ISEND(jvec(1),i_length,MPI_INTEGER8,0,0,MPI_COMM_WORLD,req(1),IERR)
      !     call MPI_WAITALL(1,req,status_array,IERR)
      !   enddo
      ! else
      !   do k=kbglobal(myid),keglobal(myid)
      !       do j=jbglobal(myid),jeglobal(myid)
      !         idmap_yz(j,k)=myid
      !       enddo
      !   enddo
      !   do ipe=1,numprocs-1
      !       jbt=jbglobal(ipe)
      !       jet=jeglobal(ipe)
      !       kbt=kbglobal(ipe)
      !       ket=keglobal(ipe)
      !       do k=kbt,ket
      !         i_length=jet-jbt+1
      !         call MPI_IRECV(idmap_yz(jbt,k),i_length,MPI_INTEGER8,IPE,0,MPI_COMM_WORLD,req(1),IERR)
      !         call MPI_WAITALL(1,req,status_array,IERR)
      !       enddo
      !   enddo
      ! endif

      !VR: again, this is much simpler than the commented block
      do ipe=0,numprocs-1
         do k=kbglobal(ipe),keglobal(ipe)
            do j=jbglobal(ipe),jeglobal(ipe)
               idmap_yz(j,k)=ipe
            enddo
         enddo
      enddo

      ! VR: fill in ghost cells in idmap     
      idmap_yz(1:ny,0)    = idmap_yz(1:ny,nz)
      idmap_yz(1:ny,nz+1) = idmap_yz(1:ny,1)

      idmap_yz(0,1:nz)    = idmap_yz(ny,1:nz)
      idmap_yz(ny+1,1:nz) = idmap_yz(1,1:nz)

      idmap_yz(0,0)       = idmap_yz(ny,nz)
      idmap_yz(0,nz+1)    = idmap_yz(ny,1)
      idmap_yz(ny+1,0)    = idmap_yz(1,nz)
      idmap_yz(ny+1,nz+1) = idmap_yz(1,1)
      ! VR: end fill ghost cells

      ! VR this is not needed since we filled the map locally on each process
      ! call MPI_BCAST(idmap_yz,size(idmap_yz),MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

      ! VR: output neigbor info for each process & idmap
      ! write(tmpfname,"(A,I0,A)") "neighbors.",myid,".dat"
      ! open(unit=512,file=TRIM(tmpfname),status='replace',action='write')
      ! write(512,"(A,I6)") "TOP=",NBRTOP
      ! write(512,"(A,I6)") "TOPLEFT=",NBRLEFTTOP
      ! write(512,"(A,I6)") "LEFT=",NBRLEFT
      ! write(512,"(A,I6)") "BOTLEFT=",NBRLEFTBOT
      ! write(512,"(A,I6)") "BOT=",NBRBOT
      ! write(512,"(A,I6)") "BOTRITE=",NBRRITEBOT
      ! write(512,"(A,I6)") "RITE=",NBRRITE
      ! write(512,"(A,I6)") "RITETOP=",NBRRITETOP
      ! close(512)
     
      ! write(tmpfname,"(A,I0,A)") "idmap_yz.",myid,".dat"
      ! open(unit=512,file=TRIM(tmpfname),status='replace',action='write',access='stream',form='unformatted')
      ! write(512) idmap_yz
      ! close(512)

      ! print *, myid, "size OF id_map is ",size(idmap_yz)

      ! call MPI_TYPE_VECTOR(int(nzl+2,4),int(nx+2,4),int((nx+2)*(nyl+2),4),MPI_DOUBLE_PRECISION,stridery,IERR)
      ! call MPI_TYPE_COMMIT(stridery,IERR)
      ! call MPI_TYPE_VECTOR(int(nyl+2,4),int(nx+2,4),int(nx+2,4)          ,MPI_DOUBLE_PRECISION,STRIDERZ,IERR)
      ! call MPI_TYPE_COMMIT(STRIDERZ,IERR)
 
      nptotp=0
      do is=1,nspec
          nptotp=nptotp+npx(is)*npy(is)*npz(is)
      enddo
      if (nptotp > nplmax) then
          if (myid == 0) then
            write(6,*) ' Increase nplmax in the input file '
            write(6,*) 'nptotp = ',nptotp
          endif
          myid_stop(myid) = 1
      endif
      do i = 0, npes-1
          i_i = i
          CALL MPI_BCAST(MYID_STOP(i),1,MPI_INTEGER8,i_i,MPI_COMM_WORLD, IERR)
        enddo
        do i = 0,npes-1
          if (myid_stop(i).ne.0) then
              call MPI_FINALIZE(IERR)
              write(*,*)"TEST HERE"
              STOP
          endif
        enddo
        if (.not.testorbt) norbskip=1

      call allocate_global_arrays
      ! call pdf_injection 

      ! Initialize nonuniform mesh
      call MESH_INIT(meshX,xaa,xbb,xmax,nax,nbx,nx) ! initialize x-mesh
      call MESH_INIT(meshY,yaa,ybb,ymax,nay,nby,ny) ! initialize y-mesh
      call MESH_INIT(meshZ,zaa,zbb,zmax,naz,nbz,nz) ! initialize z-mesh

      call MESH_INDEX(meshX,CELL,ixv_2_c_map)
      call MESH_INDEX(meshY,CELL,iyv_2_c_map)
      call MESH_INDEX(meshZ,CELL,izv_2_c_map)
      call MESH_INDEX(meshX,NODE,ixv_2_v_map)
      call MESH_INDEX(meshY,NODE,iyv_2_v_map)
      call MESH_INDEX(meshZ,NODE,izv_2_v_map)
      call MESH_INDEX(meshX,CELL,ixc_2_c_map,CELL)
      call MESH_INDEX(meshY,CELL,iyc_2_c_map,CELL)
      call MESH_INDEX(meshZ,CELL,izc_2_c_map,CELL)
      call MESH_INDEX(meshX,NODE,ixc_2_v_map,CELL)
      call MESH_INDEX(meshY,NODE,iyc_2_v_map,CELL)
      call MESH_INDEX(meshZ,NODE,izc_2_v_map,CELL)

      if (myid == 0) then
        ! write(6,*) " nl_x = ",meshX%nl
        ! do i=1,meshX%nl+3 
        !   write(6,*) " i, x;",i,meshX%xn(i),meshX%xc(i)
        ! enddo
        ! write(6,*) " nl_y = ",meshY%nl
        ! do i=1,meshY%nl+3 
        !   write(6,*) " i, y;",i,meshY%xn(i),meshY%xc(i)
        ! enddo
        ! write(6,*) " nl_z = ",meshZ%nl
        ! do i=1,meshZ%nl+3 
        !   write(6,*) " i, z;",i,meshZ%xn(i),meshZ%xc(i)
        ! enddo

        open(unit=11,file='mesh_vertices.dat',status='unknown',form='formatted')

        write(11,*) meshX%nl+1,meshY%nl+1,meshZ%nl+1
        do i=2,meshX%nl+2 
          write(11,*) meshX%xn(i)
        enddo
        do i=2,meshX%nl+2 
          write(11,*) meshX%dxc(i)
        enddo

        do i=2,meshY%nl+2 
          write(11,*) meshY%xn(i)
        enddo
        do i=2,meshY%nl+2 
          write(11,*) meshY%dxc(i)
        enddo

        do i=2,meshZ%nl+2 
          write(11,*) meshZ%xn(i)
        enddo
        do i=2,meshZ%nl+2 
          write(11,*) meshZ%dxc(i)
        enddo
    
        close(unit=11)
      endif

      ! STOP HERE IF SETUP_MESH=.TRUE.

      if (setup_mesh) then
        call MPI_FINALIZE(IERR)
        STOP
      endif

      allocate (uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1))
      ! VR    
      ! allocate (nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1))

      call date_and_time(values=time_begin)
      clock_time_re1=(time_begin(5)*3600.+time_begin(6)*60.+time_begin(7)+time_begin(8)*0.001)
 
      if (myid == 0) then
        if (restart) then
          open(unit=11,file=trim(adjustl(data_directory))//'energy.dat' ,status='old',position='append')
          open(unit=14,file=trim(adjustl(data_directory))//'time.dat' ,status='old',position='append')
          if (.not. tracking_mpi)then
            open(unit=12,file=trim(adjustl(data_directory))//'probes.dat' ,status='old',position='append')
            if (tracking_binary) then
              open(unit=13,file=trim(adjustl(data_directory))//'tracking_b.dat' ,form='unformatted',status='old',position='append')
            else
              open(unit=13,file=trim(adjustl(data_directory))//'tracking.dat' ,status='old',position='append')
            endif
          endif
        else
          open(unit=11,file=trim(adjustl(data_directory))//'energy.dat' ,status='unknown')
          open(unit=14,file=trim(adjustl(data_directory))//'time.dat' ,status='unknown')
          if (.not. tracking_mpi)then
            open(unit=12,file=trim(adjustl(data_directory))//'probes.dat' ,status='unknown')
            if (tracking_binary) then
              open(unit=13,file=trim(adjustl(data_directory))//'tracking_b.dat' ,form='unformatted',status='unknown')
            else
              open(unit=13,file=trim(adjustl(data_directory))//'tracking.dat' ,status='unknown')
            endif
          endif
        endif
      endif

      if (tracking_mpi) then
        ! write(*,*)'filename=',filename
        write(filename,"(a,i4.4,a)")'tracking_',myid,'.dat'
        write(filename2,"(a,i4.4,a)")'probes_',myid,'.dat'
        if (restart) then
          open(unit=12,file=trim(adjustl(data_directory))//filename2, status='old',position='append')
          open(unit=13,file=trim(adjustl(data_directory))//filename,form='unformatted',status='old',access='append')
        else
          open(unit=12,file=trim(adjustl(data_directory))//filename2, status='unknown')
          open(unit=13,file=trim(adjustl(data_directory))//filename,form='unformatted',status='unknown')
        endif
        ! call MPI_File_open(MPI_COMM_WORLD, trim(adjustl(data_directory))//filename, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, tracking_fh, ierr)
        ! if (ierr.ne.MPI_SUCCESS) then
        !   call MPI_Error_string(iErr,eStr,eStrLen,iErr2)
        !   write(0,*)'Error: Could not open file: ',filename
        !   write(0,*)eStr
        !   write(0,*)'Aborted.'
        !   return
        ! endif
      endif
      call opendiagfiles
   
      if (restart) then
        call makelist
        if (myid == 0) then
          ! open(unit=222,file='restart_index.dat' ,status='old')
          open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='old')
          read(222,*) restart_index,itfin
          close(222)
        endif

        call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
        nplmax6 = 6*nplmax
        ! hxv 01/10/2014

        if (myid == 0) then
          write(6,*) " "
          write(6,*) " RESTARTED FROM SET # ",restart_index
          write(6,*) " "
        endif

        ! comment out for timing on LANL machine
        do iwrite = 0,npes_over_60 
          if (mod( int(myid,8) ,npes_over_60 + 1).eq.iwrite) then
              call restrtrw(-1.0,itstart)
              call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
          endif
          ! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        enddo
         
        if (restart_index == 1) then
          restart_index=2
        else
          restart_index=1
        endif

        ! Uniform mesh - Same as is in version 5.0
        yb=(jb-1)*hy
        ye= je   *hy
        zb=(kb-1)*hz
        ze= ke   *hz
         
        ! Nonuniform mesh
        zb=meshZ%xn(kb+1)
        ze=meshZ%xn(ke+2)
        do ipe=0,npes-1
          zbglobal(ipe)=meshZ%xn(kbglobal(ipe)+1)
          zeglobal(ipe)=meshZ%xn(keglobal(ipe)+2)
        enddo
        yb=meshY%xn(jb+1)
        ye=meshY%xn(je+2)
        do ipe=0,npes-1
          ybglobal(ipe)=meshY%xn(jbglobal(ipe)+1)
          yeglobal(ipe)=meshY%xn(jeglobal(ipe)+2)
        enddo
         
         
        volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)
        xb        = 0.
        xe        = xmax
        xb_logical=MESH_UNMAP(meshX,xb)
        xe_logical=MESH_UNMAP(meshX,xe)
        yb_logical=MESH_UNMAP(meshY,yb)
        ye_logical=MESH_UNMAP(meshY,ye)
        zb_logical=MESH_UNMAP(meshZ,zb)
        ze_logical=MESH_UNMAP(meshZ,ze)
         
        do is=1,nspec
          npm=npx(is)*npy(is)*npz(is)*npes
          dfac(is)=real(ny*nz*nx)/real(npm)
          do ixe=1,nx2
              do iye=jb-1,je+1
                do ize=kb-1,ke+1
                    !                 qp_cell(ixe,iye,ize,is) = (meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)/(hx*hy*hz))*dfac(is)*frac(is)
                    qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
                enddo
              enddo
          enddo
        enddo
         
        do i=1,nxmax
          xc_uniform(i) = hx*(i-1.5)
          xv_uniform(i) = hx*(i-2.0)
        enddo
        do j=1,nymax
          yc_uniform(j) = hy*(j-0.5)
          yv_uniform(j) = hy*(j-1.0)
        enddo
        do k=1,nzmax
          zc_uniform(k) = hz*(k-0.5)
          zv_uniform(k) = hz*(k-1.0)
        enddo
         
        if (myid ==0) then
           my_short_int=it
           call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
           cycle_ascii_new=trim(adjustl(cycle_ascii))
           write(6,*) " cycle = ",cycle_ascii_new
        endif
        
        call MPI_BCAST(cycle_ascii    ,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        
      else  !VR: fresh start
        time=0.0
        call makelist
        if (myid == 0) write(6,*) " harris = ",harris
        if (myid == 0) write(6,*) " global = ",global
        if (myid == 0) write(6,*) " Yee    = ",yee
        restart_index=1

        ! call init_IA_wave
        call init_wave
        ! call init_lapd_wave

      endif !VR: end of if (restart)

      ! VR: removed output of the dipole field

      call date_and_time(values=time_end)
      clock_time_init=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
      if (myid == 0) then
        print *,'load time = ',real(clock_time_init-clock_time_re1)  
      endif
      clock_time_old = clock_time_init
      ! write(*,*)'clock time',clock_time_old

      ! hxv 01/10/2014
      if (myid == 0) then
        write(6,*) " START AT CYCLE # ",itfin
      endif

      ! Write particle data to file and exit if post_process=.true.
      ! VR: this is probably done to get particle data out of restart files
      ! VR: simulation restarts, initializes, dumps particle data, and then exits
      if (post_process) then
        if (myid == 0) then
          my_short_int=itfin
          call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
          write(6,*) " calling particle_in_volume_write with cycle_ascii = ",cycle_ascii
        endif
        call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        call particle_in_volume_write
        goto 999
      endif

      ! march forward in time
      itstart = itfin+1

      ! change how itfinish is computed
      itfinish = (max_sim_time-t_stopped)/dtwci+itstart-1

      if (myid == 0) write(6,*) 't_stopped = ',t_stopped
      if (myid == 0) write(6,*) 'itstart, itfinish = ',itstart,' ',itfinish
      if (myid == 0) write(6,*) 'nwrtdata = ',nwrtdata
      it = itstart
      time_elapsed=0.;time_begin_array=0;time_end_array=0

      !------------------------------------------------------!
      !                 Start of main time loop
      !------------------------------------------------------!
      do while(it <= itfinish)

        call get_cleanup_status(len(cleanup_status))

        if ((myid == 0).and.prntinfo) then
          write(6,*) " "
          write(6,*) " "
          WRITE(6,*) "DT = ", DT, ", T_STOPPED = ", T_STOPPED
        endif

        final_time = MPI_Wtime()
        wall_clock_elapsed = final_time-initial_time
        if (wall_clock_elapsed >= quota*3600. .and. myid == 0) then
          cleanup_status = 'CLEANUP_STATUS=TRUE'
          write(6,*) " "
          write(6,*) " "
          write(6,*) " INSUFFICIENT RUNTIME ALLOWED: DUMPING RESTART DATA"
          write(6,*) " MY QUOTA = ",QUOTA
          write(6,*) " "
          write(6,*) " "
        endif
         
        call MPI_BCAST(cleanup_status,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

        if (cleanup_status == 'CLEANUP_STATUS=EXIT') then
          goto 999
        else if (cleanup_status == 'CLEANUP_STATUS=TRUE') then
          if (myid == 0) then
            WRITE(6,*)
            WRITE(6,*) 'WRITING THE RESTART FILE'
            WRITE(6,*)
          endif
           
          itfin = it
          ! comment out for timing on LANL machine
          
          do iwrite = 0,npes_over_60 
            if (mod( int(myid,8) ,npes_over_60 + 1).eq.iwrite) then
                call restrtrw(1.0,itstart)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
          enddo
           
           if (myid == 0) then
              open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='unknown')
              write(222,*) restart_index,itfin
              close(222)
           endif
           goto 998
        endif

        call date_and_time(values=time_begin_array(:,1))

        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        if (notime == 0) then
          write(file_unit_time,"(i4,' begin    ',f15.3)") it,real(clock_time-clock_time_init)
        endif
        if (myid == 0.and.mod(it,10_8) == 0) then
          print *,"it = ",it
          print *,'system time (delta) = ',real(clock_time - clock_time_old)
          print *,'system time (total) = ',real(clock_time - clock_time_init)
          clock_time_old = clock_time
        endif

        ! determine whether or not to print diagnostic information
        if (mod(it,nprint) == 0) then
          prntinfo=.true.
        else
          prntinfo=.false.
        endif

        ! determine whether or not to write data into files
        ! if (mod(it,nwrtdata) == 0 .or. (it> 6000 .and. it<10000 .and. mod(it,25)==0 )) then
        if (mod(it,nwrtdata) == 0 ) then
          wrtdat=.true.
        else
          wrtdat=.false.
        endif

        ! HXV
        call date_and_time(values=time_begin_array(:,3))

        ! VR: we do not need injection from the boundary 
        ! do is=1,nspec
        !   call particle_newinject_linked_list(is)
        ! enddo

        ! if (it == 2 ) then
        !   CALL MPI_FINALIZE(IERR)
        !   STOP
        ! ENDIF
        ! if (notime == 0) then
        !   call date_and_time(values=time_end)
        !   clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        !   write(file_unit_time,"(i4,' injctpar ',f15.3)") it,real(clock_time-clock_time_init)
        ! endif

        call date_and_time(values=time_end_array(:,3))
        call date_and_time(values=time_begin_array(:,2))

        !VR: compute values of resistivity (that could depend on local parameters, such as
        !VR: current
        if (ndim /= 1) then
          call etacalc       ! Dietmar's resistivity
        else
          call etacalc_2d    ! Dietmar's resistivity
        endif
         
        ntot=0 ! for particle tracking

        if (notime == 0) then
          call date_and_time(values=time_end)
          clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
          write(file_unit_time,"(i4,' trans    ',f15.3)") it,real(clock_time-clock_time_init)
        endif
        !VR: trans computes density and v's 
        !VR: note that trans calls parmov, i.e. it also does a particle push
        call trans

        call date_and_time(values=time_end_array(:,2))

        ! write output to energy.dat
        if (myid==0) then
          write(11,*)it, efld, bfld, efluidt, ethermt, eptclt
          write(14,*)it, time_elapsed(1:40)
        endif

        if (notime == 0) then
          call date_and_time(values=time_end)
          clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
          write(file_unit_time,"(i4,' sortit   ',f15.3)") it,real(clock_time-clock_time_init)
        endif

        !VR: sort the particles
        call date_and_time(values=time_begin_array(:,4))
        if (mod(it,10_8) == 0) call sortit    !  sort the particles
        call date_and_time(values=time_end_array(:,4))
 
        if (notime == 0) then
          call date_and_time(values=time_end)
          clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
          write(file_unit_time,"(i4,' field    ',f15.3)") it,real(clock_time-clock_time_init)
        endif

        !VR: call field solver
        call date_and_time(values=time_begin_array(:,5))
        if (.not.testorbt) then
          if (ndim /=1) then 
            call field
          else
            call field_2d
          endif
        endif        
        call date_and_time(values=time_end_array(:,5))

        !if (it == 21000) call inject_wave
        !if (mod(it,100) == 0) call kick

        if (notime == 0) then
          call date_and_time(values=time_end)
          clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
          write(file_unit_time,"(i4,' diagnose ',f15.3)") it,real(clock_time-clock_time_init)
        endif
        !VR: call user diagnostics
        call date_and_time(values=time_begin_array(:,30))
        call user_diagnostics
        call date_and_time(values=time_end_array(:,30))
        call accumulate_time_difference(time_begin_array(1,30),time_end_array(1,30),time_elapsed(30))
        !VR: end user diagnostics
        
        ! VR: this seems to be data output region        
998     if (.not.testorbt.and.(wrtdat.or.cleanup_status == 'CLEANUP_STATUS=TRUE'.or.it == itfinish &
             .or.it == 1)) then
          if (myid ==0) then
            my_short_int=it
            call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
            cycle_ascii_new=trim(adjustl(cycle_ascii))
            write(6,*) " cycle = ",cycle_ascii_new
          endif
          call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
          call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
          
          if (myid == 0 .and. .not. MPI_IO_format) call openfiles
          
          call date_and_time(values=time_begin_array(:,6))
          if (ndim /= 1) then
            call caltemp2_global
          else
            call caltemp2_global_2d
          endif
          call date_and_time(values=time_end_array(:,6))
          call accumulate_time_difference(time_begin_array(1,6),time_end_array(1,6),time_elapsed(6))

          numvars = 18
          irecnum=1
          ! VR: in the old version, "nonuniform_mesh_global" was passed to dataout
          call dataout(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,tperp,                   &
              p_xx,p_xy,p_xz,p_yy,p_yz,p_zz,fox,foy,foz, vxs,vys,vzs,                 &
              nxmax,nymax,nzmax,file_unit,myid,                                       &
              numvars,irecnum,kb,ke,numprocs,wpiwci,jb,je,ny,nz,nylmax,nzlmax,nspecm, &
              eta,eta_times_b_dot_j,eta_par,            &
              uniform_mesh,trim(adjustl(data_directory)), trim(adjustl(cycle_ascii)),      &
              MPI_IO_format)
          if (myid == 0 .and. .not. MPI_IO_format) then
            do j=1,25
                close(file_unit(j))
            enddo
          endif

          ! Write particle data to file and exit if post_process=.true.
          if (myid == 0) then
            write(6,*) " t_begin,time,t_end = ",t_begin*wpiwci,time,t_end*wpiwci
            write(6,*) "it,nwrtparticle = ",it,nwrtparticle
          endif
          if (t_begin*wpiwci <= time .and. time <= t_end*wpiwci .and. mod( int(it,8) ,nwrtparticle)==0 .or. it==1) then
            if (myid == 0) then
              my_short_int=it
              call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
              write(6,*) " calling particle_in_volume_write with cycle_ascii = ",cycle_ascii
            endif
            call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
            call particle_in_volume_write
          endif
           if (cleanup_status == 'CLEANUP_STATUS=TRUE') goto 999
        endif

        ! if (restrt_write == 1 .and.mod(it,nwrtrestart)==0) then
        if (restrt_write == 1.and.(mod( int(it,8) ,nwrtrestart)==0 .or. it == itfinish)) then
          if (myid == 0) then
            WRITE(6,*)
            WRITE(6,*) 'WRITING THE RESTART FILE'
            WRITE(6,*)
          endif

          itfin = it

          do iwrite = 0,npes_over_60  
          if (mod( int(myid,8) ,npes_over_60 + 1).eq.iwrite) then
            call restrtrw(1.0,itstart)
          endif
          call MPI_BARRIER(MPI_COMM_WORLD,IERR)
          enddo

          if (myid == 0) then
          open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='unknown')
            write(222,*) restart_index,itfin
            close(222)
          endif

          if (restart_index == 1) then
            restart_index=2
          else
            restart_index=1
          endif
        endif

        if (notime == 0) then
          call date_and_time(values=time_end)
          clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
          write(file_unit_time,"(i4,' end      ',f15.3)") it,real(clock_time-clock_time_init)
        endif
        ! if (it == 1 ) then
        !   CALL MPI_FINALIZE(IERR)
        !   STOP
        ! endif

        time=time+dt
        it = it + 1


        ! VR: removed orientation of BZ flip

        call date_and_time(values=time_end_array(:,1))

        do j=1,5
          call accumulate_time_difference(time_begin_array(1,j),time_end_array(1,j),time_elapsed(j))
        enddo

        ! if (it == 1 ) then
        !   CALL MPI_FINALIZE(IERR)
        !   STOP
        ! ENDIF

      enddo
      !------------------------------------------------------!
      !                 End of main time loop
      !------------------------------------------------------!

      if (myid == 0) then
        close(unit=11)
        close(unit=12)
        close(unit=13)
        close(unit=14)
      endif
      if (tracking_mpi) then
      ! call MPI_File_close(tracking_fh,ierr)
      close(unit=13)
      endif
      
 999  if (notime == 0) close(file_unit_time)

      if (myid==0) then
        write(6,*) " "
        write(6,*) " "
        write(6,*) " *** RUN COMPLETED *** RUN COMPLETED *** RUN COMPLETED "
        write(6,*) " "
        write(6,*) " "
        write(6,*) " subroutine trans             (s)          =",time_elapsed(2)
        write(6,*) " subroutine injectpar         (s)          =",time_elapsed(3)
        write(6,*) " subroutine sort              (s)          =",time_elapsed(4)
        write(6,*) " subroutine field             (s)          =",time_elapsed(5)
        write(6,*) " subroutine caltemp2          (s)          =",time_elapsed(6)
        write(6,*) " subroutine user_diagnostics  (s)          =",time_elapsed(30)
        write(6,*) " total time                   (s)          =",time_elapsed(1)
        write(6,*) " "
        write(6,*) " "
        write(6,*) " In subroutine caltemp2,"
        write(6,*) "   subroutine xreal           (s)          =",time_elapsed(24)
        write(6,*) "   subroutine xrealbcc        (s)          =",time_elapsed(25)
        write(6,*) "   interpolation              (s)          =",time_elapsed(26)
        write(6,*) "   total caltemp2             (s)          =",time_elapsed(23)
        write(6,*) " "
        write(6,*) " "
        write(6,*) " In subroutine trans," 
        write(6,*) "   subroutine parmov          (s)          =",time_elapsed(7)
        write(6,*) "   subroutine energy          (s)          =",time_elapsed(8)
        write(6,*) "   total trans                (s)          =",time_elapsed(20)
        write(6,*) " "
        write(6,*) " "
        write(6,*) " In subroutine field,"
        write(6,*) "   subroutine pressgrad       (s)          =",time_elapsed(9)
        write(6,*) "   subroutine bcalc           (s)          =",time_elapsed(10)
        write(6,*) "   subroutine ecalc           (s)          =",time_elapsed(11)
        write(6,*) "   subroutine focalc          (s)          =",time_elapsed(12)
        write(6,*) "   total field                (s)          =",time_elapsed(21)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "   In subroutine parmov,"
        write(6,*) "     particle pushing         (s)          =",time_elapsed(13)
        write(6,*) "     particle exchange        (s)          =",time_elapsed(14)
        write(6,*) "     moment calculation       (s)          =",time_elapsed(15)
        write(6,*) "     total parmov             (s)          =",time_elapsed(19)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "   In subroutine bcalc,"
        write(6,*) "     subroutine ecalc         (s)          =",time_elapsed(16)
        write(6,*) "     subroutine xrealbcc      (s)          =",time_elapsed(17)
        write(6,*) "     total bcalc              (s)          =",time_elapsed(22)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "   In subroutine user_diagnostics,"
        write(6,*) "     sub virtual_probes       (s)          =",time_elapsed(31)
        write(6,*) "     sub track_particle       (s)          =",time_elapsed(32)
        write(6,*) "     total user_diagnose      (s)          =",time_elapsed(30)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "     In subroutine ecalc (called from subroutine bcalc and field),"
        write(6,*) "       subroutine xrealbcc    (s)          =",time_elapsed(18)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "communication time/total time (%)          =" &
        ,100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14)+time_elapsed(17) &
        +time_elapsed(18))/time_elapsed(1)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "Further breakdown of communication time "
        write(6,*) "  particle exchage in subroutine parmov (%) =" &
        ,100.*time_elapsed(14)/time_elapsed(1)
        write(6,*) "  subroutine xrealbcc                   (%) =" &
        ,100.*(time_elapsed(25)+time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
        write(6,*) "  subroutine xreal                      (%) =" &
        ,100.*time_elapsed(24)/time_elapsed(1)
      endif

      call MPI_FINALIZE(IERR)
      stop

    end program H3D

    !------------------------------------------------------!
    ! computes velocities?
    ! what is the difference between vxs and vix?
    !------------------------------------------------------!
      subroutine trans        
      use parameter_mod
      use MESH2D
      implicit none

      integer*8:: is,i,j,k,jbmin,jbmax,kbmin,kbmax
      integer*4:: time_begin(8),time_end(8)
      double precision:: dttmp,dns_tmp
 
      call date_and_time(values=time_begin_array(:,20))
 
      do is=1,nspec
        DO K=KB-1,KE+1
          do j=jb-1,je+1
            do i=1,nx2
              dns(i,j,k,is)=1.e-10
              dnsh(i,j,k,is)=1.e-10
              vxs(i,j,k,is)=0.
              vys(i,j,k,is)=0.
              vzs(i,j,k,is)=0.
              if (is == 1) then
                deno(i,j,k)=den(i,j,k)
                vixo(i,j,k)=vix(i,j,k)
                viyo(i,j,k)=viy(i,j,k)
                vizo(i,j,k)=viz(i,j,k)
                den(i,j,k)=0.
                denh(i,j,k)=0.
                vix(i,j,k)=0.
                viy(i,j,k)=0.
                viz(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      enddo
 
      if (notime == 0) then
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        write(file_unit_time,"(i4,' parmov   ',f15.3)") it,real(clock_time-clock_time_init)
      endif
      call date_and_time(values=time_begin_array(:,7))
  
      if (ndim /= 1) then
         call parmov
       else
         call parmov_2d
       endif

      call date_and_time(values=time_end_array(:,7))
      call accumulate_time_difference(time_begin_array(1,7),time_end_array(1,7),time_elapsed(7))
 
      if (testorbt) return

      if (.false.) then
      do is=1,nspec
        do k=kb-1,ke+1
          do j=jb-1,je+1
            do i=1,nx2
              ! Nonuniform mesh
              cell_volume_ratio = hx*hy*hz/(meshX%dxc(i)*meshY%dxc(j+1)*meshZ%dxc(k+1))
              dns_tmp=dnsh(i,j,k,is)
              ! Uniform mesh - Same as is in version 5.0
              ! if (dns_tmp <= denmin) dns_tmp=1.d+10

              ! Nonuniform mesh
              ! if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=1.d+10
              if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=denmin/cell_volume_ratio ! July 21, 2010

              vxs(i,j,k,is)=vxs(i,j,k,is)/dns_tmp
              vys(i,j,k,is)=vys(i,j,k,is)/dns_tmp
              vzs(i,j,k,is)=vzs(i,j,k,is)/dns_tmp

              ! Uniform mesh - Same as is in version 5.0
              ! dns(i,j,k,is)=dns(i,j,k,is)

              ! Nonuniform mesh
              ! dns(i,j,k,is)=dns(i,j,k,is)*cell_volume_ratio
            enddo
          enddo
        enddo
      enddo
      endif
 
      do is=1,nspec
        do k=kb-1,ke+1
          do j=jb-1,je+1
            do i=1,nx2 
            den(i,j,k)=den(i,j,k)+dns(i,j,k,is)*qspec(is) 
            denh(i,j,k)=denh(i,j,k)+dnsh(i,j,k,is)*qspec(is) 
            !vix(i,j,k)=vix(i,j,k)+qspec(is)*dnsh(i,j,k,is)*vxs(i,j,k,is) 
            !viy(i,j,k)=viy(i,j,k)+qspec(is)*dnsh(i,j,k,is)*vys(i,j,k,is) 
            !viz(i,j,k)=viz(i,j,k)+qspec(is)*dnsh(i,j,k,is)*vzs(i,j,k,is)
            vix(i,j,k)=vix(i,j,k)+qspec(is)*vxs(i,j,k,is) 
            viy(i,j,k)=viy(i,j,k)+qspec(is)*vys(i,j,k,is) 
            viz(i,j,k)=viz(i,j,k)+qspec(is)*vzs(i,j,k,is)
            enddo
          enddo
        enddo
      enddo

      ! Apply Boundary Conditions
      if (ndim /= 1) then
      call XREALBCC(DEN,1_8,NX,NY,NZ)
      call XREALBCC(DENH,1_8,NX,NY,NZ)
      call XREALBCC(VIX,1_8,NX,NY,NZ)
      call XREALBCC(VIY,1_8,NX,NY,NZ)
      call XREALBCC(VIZ,1_8,NX,NY,NZ)
      else
      call XREALBCC_2D(DEN,1_8,NX,NY,NZ)
      call XREALBCC_2D(DENH,1_8,NX,NY,NZ)
      call XREALBCC_2D(VIX,1_8,NX,NY,NZ)
      call XREALBCC_2D(VIY,1_8,NX,NY,NZ)
      call XREALBCC_2D(VIZ,1_8,NX,NY,NZ)
      endif

      ! smooth density and velocity
      if (smoothing) then
        if (ndim /=1) then
          call nsmth(DEN)
          call nsmth(DENH)
          call nsmth(VIX)
          call nsmth(VIY)
          call nsmth(VIZ)
        else
          call nsmth_2d(DEN,NX2,NY2,NZ2)
          call nsmth_2d(DENH,NX2,NY2,NZ2)
          call nsmth_2d(VIX,NX2,NY2,NZ2)
          call nsmth_2d(VIY,NX2,NY2,NZ2)
          call nsmth_2d(VIZ,NX2,NY2,NZ2)
        endif
      endif

      do k=kb-1,ke+1
          do j=jb-1,je+1
          do i=1,nx2
            den(i,j,k)=max(denmin,den(i,j,k))
            !pressure calculation moved to pressgrad
            !pe(i,j,k) =te0*den(i,j,k)**gama
            vix(i,j,k)=vix(i,j,k)/denh(i,j,k)
            viy(i,j,k)=viy(i,j,k)/denh(i,j,k)
            viz(i,j,k)=viz(i,j,k)/denh(i,j,k)
          enddo
        enddo
      enddo

      if (it == 0) then
         deno=den;vixo=vix;viyo=viy;vizo=viz
      endif
 
      call date_and_time(values=time_begin_array(:,8))
 
      if (mod(it,10)==0) call energy
 
      call date_and_time(values=time_end_array(:,8))
      call accumulate_time_difference(time_begin_array(1,8),time_end_array(1,8),time_elapsed(8))

      kbmin = kb-1
      kbmax = ke+1

      jbmin = jb-1
      jbmax = je+1
 
      call date_and_time(values=time_end_array(:,20))
      call accumulate_time_difference(time_begin_array(1,20),time_end_array(1,20),time_elapsed(20))
 
      return
    end subroutine trans

    !------------------------------------------------------!
    ! compute perp and par temperature and pressure tensor
    !------------------------------------------------------!
      subroutine caltemp2_global
 
      use parameter_mod
      use MESH2D
      implicit none

      double precision:: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
      integer*8 ix,iy,iz,ixp1,iyp1,izp1,iiy,iiye,iiz,iize,is,l,iix,iixe
      double precision vxa,vya,vza,rfrac,vxavg,vxavg1,vxavg2 &
           ,vyavg,vyavg1,vyavg2,vzavg,vzavg1,vzavg2,wperp2,wpar,wmult
      double precision w1,w2,w3,w4,w5,w6,w7,w8,h,hh,dns1,dns2,bxa,bya,bza,btota,dnst
 
      call date_and_time(values=time_begin_array(:,23))

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

      tpar=0.
      tperp=0.

      p_xx=0.;p_xy=0.;p_xz=0.;p_yy=0.;p_yz=0.;p_zz=0.

      if (nspec >= 2) then
         rfrac = frac(2)/frac(1)
      else
         rfrac = 0.
      endif
 
      call date_and_time(values=time_begin_array(:,26))
      DO IS=1,NSPEC
        wmult=wspec(is)
        h=dt*qspec(is)/wmult
        hh=.5*h
        dpedx = 0.
        DO IIZE = KB-1,KE
          DO IIYE = JB-1,JE
            DO IIXE = 1, NX1
              NP=IPHEAD(IIXE,IIYE,IIZE,IS)

              ! begin advance of particle position and velocity
              ! If dt=0, skip
              DO WHILE (NP.NE.0)
                L=NP

                ! Uniform mesh - Same as is in version 5.0
                ! rx=hxi*x(l)+1.5000000000000001
                ! ry=hyi*y(l)+0.5000000000000001d+00
                ! rz=hzi*z(l)+0.5000000000000001d+00
                ! ix=rx
                ! iy=ry
                ! iz=rz
                ! fx=rx-ix
                ! fy=ry-iy
                ! fz=rz-iz

                ! Nonuniform mesh - using MESH_UNMAP
                rx=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
                ry=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
                rz=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
                ix=rx
                iy=ry
                iz=rz
                fx=rx-ix
                fy=ry-iy
                fz=rz-iz
                iy=iy-1   ! integer index in y direction starts at 0
                iz=iz-1   ! integer index in z direction starts at 0

                ixp1 = ix+1
                iyp1 = iy+1
                izp1 = iz+1

                w1=(1.-fx)*(1.-fy)*(1.-fz)
                w2=    fx *(1.-fy)*(1.-fz)
                w3=(1.-fx)*    fy *(1.-fz)
                w4=    fx*     fy *(1.-fz)
                w5=(1.-fx)*(1.-fy)*    fz
                w6=    fx *(1.-fy)*    fz
                w7=(1.-fx)*    fy*     fz
                w8=    fx*     fy*     fz

                dnst= dnsh(ix  ,iy  ,iz  ,is)*w1+dnsh(ixp1,iy  ,iz  ,is)*w2  &
                +     dnsh(ix  ,iyp1,iz  ,is)*w3+dnsh(ixp1,iyp1,iz  ,is)*w4  &
                +     dnsh(ix  ,iy  ,izp1,is)*w5+dnsh(ixp1,iy  ,izp1,is)*w6  &
                +     dnsh(ix  ,iyp1,izp1,is)*w7+dnsh(ixp1,iyp1,izp1,is)*w8

                vxavg=vxs(ix  ,iy  ,iz  ,is)*w1+vxs(ixp1,iy  ,iz  ,is)*w2  &
                +      vxs(ix  ,iyp1,iz  ,is)*w3+vxs(ixp1,iyp1,iz  ,is)*w4  &
                +      vxs(ix  ,iy  ,izp1,is)*w5+vxs(ixp1,iy  ,izp1,is)*w6  &
                +      vxs(ix  ,iyp1,izp1,is)*w7+vxs(ixp1,iyp1,izp1,is)*w8
                vxavg=vxavg/dnst
                

                vyavg=vys(ix  ,iy  ,iz  ,is)*w1+vys(ixp1,iy  ,iz  ,is)*w2  &
                +      vys(ix  ,iyp1,iz  ,is)*w3+vys(ixp1,iyp1,iz  ,is)*w4  &
                +      vys(ix  ,iy  ,izp1,is)*w5+vys(ixp1,iy  ,izp1,is)*w6  &
                +      vys(ix  ,iyp1,izp1,is)*w7+vys(ixp1,iyp1,izp1,is)*w8
                vyavg=vyavg/dnst
                

                vzavg=vzs(ix  ,iy  ,iz  ,is)*w1+vzs(ixp1,iy  ,iz  ,is)*w2  &
                +      vzs(ix  ,iyp1,iz  ,is)*w3+vzs(ixp1,iyp1,iz  ,is)*w4  &
                +      vzs(ix  ,iy  ,izp1,is)*w5+vzs(ixp1,iy  ,izp1,is)*w6  &
                +      vzs(ix  ,iyp1,izp1,is)*w7+vzs(ixp1,iyp1,izp1,is)*w8
                vzavg=vzavg/dnst


                vxa=vx(l)-vxavg
                vya=vy(l)-vyavg
                vza=vz(l)-vzavg

                bxa  =bx(ix  ,iy  ,iz  )*w1+bx(ixp1,iy  ,iz  )*w2  &
                +     bx(ix  ,iyp1,iz  )*w3+bx(ixp1,iyp1,iz  )*w4  &
                +     bx(ix  ,iy  ,izp1)*w5+bx(ixp1,iy  ,izp1)*w6  &
                +     bx(ix  ,iyp1,izp1)*w7+bx(ixp1,iyp1,izp1)*w8

                bya  =by(ix  ,iy  ,iz  )*w1+by(ixp1,iy  ,iz  )*w2  &
                +     by(ix  ,iyp1,iz  )*w3+by(ixp1,iyp1,iz  )*w4  &
                +     by(ix  ,iy  ,izp1)*w5+by(ixp1,iy  ,izp1)*w6  &
                +     by(ix  ,iyp1,izp1)*w7+by(ixp1,iyp1,izp1)*w8

                bza  =bz(ix  ,iy  ,iz  )*w1+bz(ixp1,iy  ,iz  )*w2  &
                +     bz(ix  ,iyp1,iz  )*w3+bz(ixp1,iyp1,iz  )*w4  &
                +     bz(ix  ,iy  ,izp1)*w5+bz(ixp1,iy  ,izp1)*w6  &
                +     bz(ix  ,iyp1,izp1)*w7+bz(ixp1,iyp1,izp1)*w8

                btota=sqrt(bxa**2+bya**2+bza**2)
                if (btota < 1.e-20) btota=1.e-20
                wpar=(vxa*bxa+vya*bya+vza*bza)/btota
                wperp2=vxa**2+vya**2+vza**2-wpar**2
                xx=vxa*vxa
                xy=vxa*vya
                xz=vxa*vza
                yy=vya*vya
                yz=vya*vza
                zz=vza*vza

                tpar (ix  ,iy  ,iz  ,is)=tpar (ix  ,iy  ,iz  ,is)+qp(np)*w1*wpar*wpar
                tpar (ixp1,iy  ,iz  ,is)=tpar (ixp1,iy  ,iz  ,is)+qp(np)*w2*wpar*wpar 
                tpar (ix  ,iyp1,iz  ,is)=tpar (ix  ,iyp1,iz  ,is)+qp(np)*w3*wpar*wpar 
                tpar (ixp1,iyp1,iz  ,is)=tpar (ixp1,iyp1,iz  ,is)+qp(np)*w4*wpar*wpar 
                tpar (ix  ,iy  ,izp1,is)=tpar (ix  ,iy  ,izp1,is)+qp(np)*w5*wpar*wpar 
                tpar (ixp1,iy  ,izp1,is)=tpar (ixp1,iy  ,izp1,is)+qp(np)*w6*wpar*wpar 
                tpar (ix  ,iyp1,izp1,is)=tpar (ix  ,iyp1,izp1,is)+qp(np)*w7*wpar*wpar 
                tpar (ixp1,iyp1,izp1,is)=tpar (ixp1,iyp1,izp1,is)+qp(np)*w8*wpar*wpar 
                tperp(ix  ,iy  ,iz  ,is)=tperp(ix  ,iy  ,iz  ,is)+qp(np)*w1*wperp2 
                tperp(ixp1,iy  ,iz  ,is)=tperp(ixp1,iy  ,iz  ,is)+qp(np)*w2*wperp2 
                tperp(ix  ,iyp1,iz  ,is)=tperp(ix  ,iyp1,iz  ,is)+qp(np)*w3*wperp2 
                tperp(ixp1,iyp1,iz  ,is)=tperp(ixp1,iyp1,iz  ,is)+qp(np)*w4*wperp2
                tperp(ix  ,iy  ,izp1,is)=tperp(ix  ,iy  ,izp1,is)+qp(np)*w5*wperp2 
                tperp(ixp1,iy  ,izp1,is)=tperp(ixp1,iy  ,izp1,is)+qp(np)*w6*wperp2
                tperp(ix  ,iyp1,izp1,is)=tperp(ix  ,iyp1,izp1,is)+qp(np)*w7*wperp2 
                tperp(ixp1,iyp1,izp1,is)=tperp(ixp1,iyp1,izp1,is)+qp(np)*w8*wperp2
                dpedx(ix  ,iy  ,iz  )=dpedx(ix  ,iy  ,iz  )+qp(np)*w1
                dpedx(ixp1,iy  ,iz  )=dpedx(ixp1,iy  ,iz  )+qp(np)*w2 
                dpedx(ix  ,iyp1,iz  )=dpedx(ix  ,iyp1,iz  )+qp(np)*w3 
                dpedx(ixp1,iyp1,iz  )=dpedx(ixp1,iyp1,iz  )+qp(np)*w4 
                dpedx(ix  ,iy  ,izp1)=dpedx(ix  ,iy  ,izp1)+qp(np)*w5 
                dpedx(ixp1,iy  ,izp1)=dpedx(ixp1,iy  ,izp1)+qp(np)*w6 
                dpedx(ix  ,iyp1,izp1)=dpedx(ix  ,iyp1,izp1)+qp(np)*w7 
                dpedx(ixp1,iyp1,izp1)=dpedx(ixp1,iyp1,izp1)+qp(np)*w8 
 
                p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
                p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
                p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
                p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 
                p_xx (ix  ,iy  ,izp1,is)=p_xx (ix  ,iy  ,izp1,is)+qp(np)*w5*xx 
                p_xx (ixp1,iy  ,izp1,is)=p_xx (ixp1,iy  ,izp1,is)+qp(np)*w6*xx 
                p_xx (ix  ,iyp1,izp1,is)=p_xx (ix  ,iyp1,izp1,is)+qp(np)*w7*xx 
                p_xx (ixp1,iyp1,izp1,is)=p_xx (ixp1,iyp1,izp1,is)+qp(np)*w8*xx 
 
                p_xy (ix  ,iy  ,iz  ,is)=p_xy (ix  ,iy  ,iz  ,is)+qp(np)*w1*xy
                p_xy (ixp1,iy  ,iz  ,is)=p_xy (ixp1,iy  ,iz  ,is)+qp(np)*w2*xy 
                p_xy (ix  ,iyp1,iz  ,is)=p_xy (ix  ,iyp1,iz  ,is)+qp(np)*w3*xy 
                p_xy (ixp1,iyp1,iz  ,is)=p_xy (ixp1,iyp1,iz  ,is)+qp(np)*w4*xy 
                p_xy (ix  ,iy  ,izp1,is)=p_xy (ix  ,iy  ,izp1,is)+qp(np)*w5*xy 
                p_xy (ixp1,iy  ,izp1,is)=p_xy (ixp1,iy  ,izp1,is)+qp(np)*w6*xy 
                p_xy (ix  ,iyp1,izp1,is)=p_xy (ix  ,iyp1,izp1,is)+qp(np)*w7*xy 
                p_xy (ixp1,iyp1,izp1,is)=p_xy (ixp1,iyp1,izp1,is)+qp(np)*w8*xy 
 
                p_xz (ix  ,iy  ,iz  ,is)=p_xz (ix  ,iy  ,iz  ,is)+qp(np)*w1*xz
                p_xz (ixp1,iy  ,iz  ,is)=p_xz (ixp1,iy  ,iz  ,is)+qp(np)*w2*xz 
                p_xz (ix  ,iyp1,iz  ,is)=p_xz (ix  ,iyp1,iz  ,is)+qp(np)*w3*xz 
                p_xz (ixp1,iyp1,iz  ,is)=p_xz (ixp1,iyp1,iz  ,is)+qp(np)*w4*xz 
                p_xz (ix  ,iy  ,izp1,is)=p_xz (ix  ,iy  ,izp1,is)+qp(np)*w5*xz 
                p_xz (ixp1,iy  ,izp1,is)=p_xz (ixp1,iy  ,izp1,is)+qp(np)*w6*xz 
                p_xz (ix  ,iyp1,izp1,is)=p_xz (ix  ,iyp1,izp1,is)+qp(np)*w7*xz 
                p_xz (ixp1,iyp1,izp1,is)=p_xz (ixp1,iyp1,izp1,is)+qp(np)*w8*xz 
 
                p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
                p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
                p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
                p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
                p_yy (ix  ,iy  ,izp1,is)=p_yy (ix  ,iy  ,izp1,is)+qp(np)*w5*yy 
                p_yy (ixp1,iy  ,izp1,is)=p_yy (ixp1,iy  ,izp1,is)+qp(np)*w6*yy 
                p_yy (ix  ,iyp1,izp1,is)=p_yy (ix  ,iyp1,izp1,is)+qp(np)*w7*yy 
                p_yy (ixp1,iyp1,izp1,is)=p_yy (ixp1,iyp1,izp1,is)+qp(np)*w8*yy 
 
                p_yz (ix  ,iy  ,iz  ,is)=p_yz (ix  ,iy  ,iz  ,is)+qp(np)*w1*yz
                p_yz (ixp1,iy  ,iz  ,is)=p_yz (ixp1,iy  ,iz  ,is)+qp(np)*w2*yz 
                p_yz (ix  ,iyp1,iz  ,is)=p_yz (ix  ,iyp1,iz  ,is)+qp(np)*w3*yz 
                p_yz (ixp1,iyp1,iz  ,is)=p_yz (ixp1,iyp1,iz  ,is)+qp(np)*w4*yz 
                p_yz (ix  ,iy  ,izp1,is)=p_yz (ix  ,iy  ,izp1,is)+qp(np)*w5*yz 
                p_yz (ixp1,iy  ,izp1,is)=p_yz (ixp1,iy  ,izp1,is)+qp(np)*w6*yz 
                p_yz (ix  ,iyp1,izp1,is)=p_yz (ix  ,iyp1,izp1,is)+qp(np)*w7*yz 
                p_yz (ixp1,iyp1,izp1,is)=p_yz (ixp1,iyp1,izp1,is)+qp(np)*w8*yz 
 
                p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
                p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
                p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
                p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
                p_zz (ix  ,iy  ,izp1,is)=p_zz (ix  ,iy  ,izp1,is)+qp(np)*w5*zz 
                p_zz (ixp1,iy  ,izp1,is)=p_zz (ixp1,iy  ,izp1,is)+qp(np)*w6*zz 
                p_zz (ix  ,iyp1,izp1,is)=p_zz (ix  ,iyp1,izp1,is)+qp(np)*w7*zz 
                p_zz (ixp1,iyp1,izp1,is)=p_zz (ixp1,iyp1,izp1,is)+qp(np)*w8*zz 
 
                np=link(np)
              ENDDO
            ENDDO
          ENDDO
        ENDDO


        call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(dpedx(1,jb-1,kb-1   ),NX,NY,NZ)

        call XREAL(p_xx (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_xy (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_xz (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_yy (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_yz (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_zz (1,jb-1,kb-1,is),NX,NY,NZ)

        DO IZ = KB-1,KE
          DO IY = JB-1,JE
            DO IX = 1, NX1
              if (dpedx(ix,iy,iz) /= 0.) then
                tpar (ix,iy,iz,is) = tpar (ix,iy,iz,is)/(   tx0(is)*dpedx(ix,iy,iz))
                tperp(ix,iy,iz,is) = tperp(ix,iy,iz,is)/(2.*tx0(is)*dpedx(ix,iy,iz))
              endif
            ENDDO
          ENDDO
        ENDDO

        DO IIZ=KB-1,KE+1
          DO IIY=JB-1,JE+1
            DO IIX=1,NX2
              p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_xy(iix,iiy,iiz,is) = p_xy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_xz(iix,iiy,iiz,is) = p_xz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_yz(iix,iiy,iiz,is) = p_yz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            ENDDO
          ENDDO
        ENDDO

        ! p_xx(:,:,:,is)=p_xx(:,:,:,is)/(tx0(is)*frac(is))
        ! p_xy(:,:,:,is)=p_xy(:,:,:,is)/(tx0(is)*frac(is))
        ! p_xz(:,:,:,is)=p_xz(:,:,:,is)/(tx0(is)*frac(is))
        ! p_yy(:,:,:,is)=p_yy(:,:,:,is)/(tx0(is)*frac(is))
        ! p_yz(:,:,:,is)=p_yz(:,:,:,is)/(tx0(is)*frac(is))
        ! p_zz(:,:,:,is)=p_zz(:,:,:,is)/(tx0(is)*frac(is))

      ENDDO

      call date_and_time(values=time_end_array(:,26))
      call accumulate_time_difference(time_begin_array(1,26) &
                                    ,time_end_array(1,26) &
                                    ,time_elapsed(26))

      ! do is=1,nspec
      !   call date_and_time(values=time_begin_array(:,24))
      !   call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
      !   call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
      !   call date_and_time(values=time_end_array(:,24))
      !   call accumulate_time_difference(time_begin_array(1,24) &
      !                                  ,time_end_array(1,24) &
      !                                  ,time_elapsed(24))

      !   call date_and_time(values=time_begin_array(:,25))
      !   call XREALBCC(tpar (1,jb-1,kb-1,is),1,NX,NY,NZ)
      !   call XREALBCC(tperp(1,jb-1,kb-1,is),1,NX,NY,NZ)
      !   call date_and_time(values=time_end_array(:,25))
      !   call accumulate_time_difference(time_begin_array(1,25) &
      !                                  ,time_end_array(1,25) &
      !                                  ,time_elapsed(25))

      ! enddo


      ! call date_and_time(values=time_begin_array(:,26))
      ! do is=1,nspec
      !   do k=kb-1,ke+1
      !     do j = jb-1,je+1
      !       do i=1,nx2
      !         if (is == 1) then
      !           dns1=dns(i,j,k,1)/(dfac(1)*frac(1))
      !           dns2=0.
      !           denum=dns1+rfrac*dns2
      !         else
      !           denum=dns(i,j,k,is)/(dfac(is)*frac(is))
      !         endif
      !         if (denum < denmin)  then
      !           tpar(i,j,k,is)=1.e-5
      !           tperp(i,j,k,is)=1.e-5
      !         else
      !           denum=denum*tx0(is)
      !           tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
      !           tperp(i,j,k,is)=0.5*tperp(i,j,k,is)*wspec(is)/denum
      !         endif
      !       enddo
      !     enddo
      !   enddo
      ! enddo
      ! call date_and_time(values=time_end_array(:,26))
      ! call accumulate_time_difference(time_begin_array(1,26) &
      ! &                               ,time_end_array(1,26) &
      ! &                               ,time_elapsed(26))

 
      call date_and_time(values=time_end_array(:,23))
      call accumulate_time_difference(time_begin_array(1,23) &
                                    ,time_end_array(1,23) &
                                    ,time_elapsed(23))

      return
    end subroutine caltemp2_global

