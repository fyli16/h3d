!***************************************************************************
!                                                                          *
!                                 VERSION 6.0                              *
!                           YURI'S NONUNIFORM MESH                         *
!                           3D IMPLEMENTATION ONLY                         *
!                      UNIFORM LOADING IN PHYSICAL SPACE                   *
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       *
!                                                                          *
!***************************************************************************

    program h3d 
      use parameter_mod
      use initialize
      use functions_f90
      use mesh2d
      
      implicit none
      integer*8 :: i, irecnum, ixe, iye, ize, j, jbt, jet, k, kbt, ket, nplmax6, numvars
      real*8 :: rnorm, pifac
      integer*4 :: time_begin(8), time_end(8), is
      integer*8 :: itstart, itfinish
      real*8 :: clock_time_re1
      real*8, dimension(:,:,:), allocatable :: uniform_mesh      
      !VR : allocating a global mesh can not work on large runs with small amount of memory per rank
      !VR : double precision, dimension(:,:,:), allocatable:: nonuniform_mesh_global
      character (len=240) :: filename, filename2, tmpfname
      integer :: iwrite, ierr2, eStrLen
      character(len=1024) :: eStr

      call init_mpi()
      call read_input()
      call mpi_decomposition()
      call set_parameters()     

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
      !     do j=jbglobal(myid),jeglobal(myid)
      !         idmap_yz(j,k)=myid
      !     enddo
      !   enddo
      !   do ipe=1,numprocs-1
      !     jbt=jbglobal(ipe)
      !     jet=jeglobal(ipe)
      !     kbt=kbglobal(ipe)
      !     ket=keglobal(ipe)
      !     do k=kbt,ket
      !         i_length=jet-jbt+1
      !         call MPI_IRECV(idmap_yz(jbt,k),i_length,MPI_INTEGER8,IPE,0,MPI_COMM_WORLD,req(1),IERR)
      !         call MPI_WAITALL(1,req,status_array,IERR)
      !     enddo
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

      !VR: fill in ghost cells in idmap     
      idmap_yz(1:ny,0)    = idmap_yz(1:ny,nz)
      idmap_yz(1:ny,nz+1) = idmap_yz(1:ny,1)

      idmap_yz(0,1:nz)    = idmap_yz(ny,1:nz)
      idmap_yz(ny+1,1:nz) = idmap_yz(1,1:nz)

      idmap_yz(0,0)       = idmap_yz(ny,nz)
      idmap_yz(0,nz+1)    = idmap_yz(ny,1)
      idmap_yz(ny+1,0)    = idmap_yz(1,nz)
      idmap_yz(ny+1,nz+1) = idmap_yz(1,1)
      !VR: end fill ghost cells

      !VR this is not needed since we filled the map locally on each process
      ! call MPI_BCAST(idmap_yz,size(idmap_yz),MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

      !VR: output neigbor info for each process & idmap
      
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

      call MPI_TYPE_VECTOR(int(nzlmax+2,4),int(nx+2,4),int((nx+2)*(nylmax+2),4),MPI_DOUBLE_PRECISION,stridery,IERR)
      call MPI_TYPE_COMMIT(stridery,IERR)
      call MPI_TYPE_VECTOR(int(nylmax+2,4),int(nx+2,4),int(nx+2,4)          ,MPI_DOUBLE_PRECISION,STRIDERZ,IERR)
      call MPI_TYPE_COMMIT(STRIDERZ,IERR)
 
     
      nptotp=0
      do is=1,nspec
          nptotp=nptotp+npx(is)*npy(is)*npz(is)
      enddo
      if (nptotp > nplmax) then
          if (myid == 0) then
            write(6,*) 'Increase nplmax in the input file '
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
      !  call pdf_injection

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
      ! VR: allocate (nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1))

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
        !call MPI_File_open(MPI_COMM_WORLD, trim(adjustl(data_directory))//filename, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, tracking_fh, ierr)
        ! if (ierr.ne.MPI_SUCCESS) then
        !   call MPI_Error_string(iErr,eStr,eStrLen,iErr2)
        !   write(0,*)'Error: Could not open file: ',filename
        !   write(0,*) eStr
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
          read(222,*) restart_index, itfin
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
         
        !       Nonuniform mesh
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

        call init_wave
        
      endif !VR: end of if (restart)


      !VR: removed output of the dipole field

      call date_and_time(values=time_end)
      clock_time_init=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
      if (myid == 0) then
        print *,'load time = ',real(clock_time_init-clock_time_re1)  
      endif
      clock_time_old = clock_time_init
      !write(*,*)'clock time',clock_time_old
  
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
      itfinish = (tmax-t_stopped)/dtwci+itstart-1

      if (myid == 0) write(6,*) 't_stopped = ',t_stopped
      if (myid == 0) write(6,*) 'itstart, itfinish = ',itstart,' ',itfinish
      if (myid == 0) write(6,*) 'nwrtdata = ',nwrtdata
      it = itstart
      time_elapsed=0.;time_begin_array=0;time_end_array=0

      !#########################################################
      !>>         start of main time loop
      !#########################################################
      do while(it <= itfinish)
        call get_cleanup_status(len(cleanup_status))

        if ((myid == 0).and.prntinfo) then
          write(6,*) " "
          WRITE(6,*) "DT = ", dt, ", T_STOPPED = ", t_stopped
        endif

        final_time = MPI_Wtime()
        wall_clock_elapsed= final_time-initial_time
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
              write(222,*) restart_index, itfin
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
          print *, 'it=', it, ', delta_time=', real(clock_time - clock_time_old), ', tot_time=', real(clock_time - clock_time_init)
          clock_time_old = clock_time
        endif

        ! determine whether or not to print diagnostic information
        if (mod(it,nprint) == 0) then
          prntinfo=.true.
        else
          prntinfo=.false.
        endif

        ! determine whether or not to write data into files
        !if (mod(it,nwrtdata) == 0 .or. (it> 6000 .and. it<10000 .and. mod(it,25)==0 )) then
        if (mod(it,nwrtdata) == 0 ) then
          wrtdat=.true.
        else
          wrtdat=.false.
        endif

        ! HXV
        call date_and_time(values=time_begin_array(:,3))

        !VR: we do not need injection from the boundary 
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
        
        !VR: this seems to be data output region        
998     if (.not.testorbt.and.(wrtdat.or.cleanup_status == 'CLEANUP_STATUS=TRUE'.or.it == itfinish       &
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
           !VR: in the old version, "nonuniform_mesh_global" was passed to dataout
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
          !  if (myid == 0) then
          !     write(6,*) " t_begin,time,t_end = ",t_begin*wpiwci,time,t_end*wpiwci
          !     write(6,*) "it,nwrtparticle = ",it,nwrtparticle
          !  endif
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

        ! write restart files at specified steps
        if (restrt_write==1 .and. mod(int(it,8), nwrtrestart)==0) then
          if (myid == 0) then
            WRITE(6,*)
            WRITE(6,*) 'WRITING RESTART FILES'
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
            write(222,*) restart_index, itfin
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
        !  if (it == 1 ) then
        !    CALL MPI_FINALIZE(IERR)
        !    STOP
        !  ENDIF

        time = time + dt
        it = it + 1

        ! VR: removed orientation of BZ flip

        call date_and_time(values=time_end_array(:,1))

        do j=1,5
          call accumulate_time_difference(time_begin_array(1,j),time_end_array(1,j),time_elapsed(j))
        enddo

        !  if (it == 1 ) then
        !    CALL MPI_FINALIZE(IERR)
        !    STOP
        !  ENDIF

      enddo  
      !#########################################################
      !>>         end of main time loop
      !#########################################################

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

      if (MYID.EQ.0) then
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
    end program h3d