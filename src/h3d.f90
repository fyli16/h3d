!---------------------------------------------------------------------
!                                 H3D (V6.0)                               
!                           YURI'S NONUNIFORM MESH                         
!                           3D IMPLEMENTATION ONLY                         
!                      UNIFORM LOADING IN PHYSICAL SPACE                   
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       
!---------------------------------------------------------------------

program h3d 

  use m_parameter
  use m_init

  use m_eta
  use m_particle
  use m_field
  use m_diag
  use m_restart

  implicit none

  ! read input deck
  call init_input

  ! init global arrays
  call init_arrays

  ! init mesh 
  call init_mesh

  ! restart or a fresh start
  if (restart) then 
    call init_restart 
  else
    call init_waves
    call init_particles 
  endif 

  ! main simulation loops
  call run_sim

  ! close program and exit
  call close_sim

  contains 

  !---------------------------------------------------------------------
  ! main simulation loops
  !---------------------------------------------------------------------
  subroutine run_sim
    integer :: i

    ! initialize time arrays, and get a time stamp
    ! just before entering the simulation loop
    time_elapsed=0.; time_begin=0; time_end=0
    call get_time(clock_init)
    clock_old = clock_init
    
    ! main simulation loop
    if (myid==0) then 
      print*, " " 
      print*, "Executing main simulation loops:"
      print*, "-------------------------------------------------"
    endif 

    do while(it <= itfinish)

      ! timing a full loop
      call date_and_time(values=time_begin(:,1))

      ! print time & step info
      if (myid==0 .and. mod(it,n_print)==0) then
        call get_time(clock_now)
        write(6,"(A5,I7,A2,I7,A11,F8.3,A14,F8.3,A12,F8.3)") 'it = ', it, '/', itfinish, &
                      ',   time = ', time, &
                      ',   delta_t = ', real(clock_now - clock_old), &
                      ',   tot_t = ', real(clock_now - clock_init)
        clock_old = clock_now
      endif

      ! update resistivity (eta)
      ! (which could depend on local parameters such as current)
      call date_and_time(values=time_begin(:,2))
      call update_eta
      call date_and_time(values=time_end(:,2))
      call add_time(time_begin(1,2),time_end(1,2),time_elapsed(2))

      ! update particles
      ! (calculate density and v's, push & sort particles)
      call date_and_time(values=time_begin(:,3))
      call update_particles
      call date_and_time(values=time_end(:,3))
      call add_time(time_begin(1,3),time_end(1,3),time_elapsed(3))

      ! update fields
      call date_and_time(values=time_begin(:,5))
      call update_fields
      call date_and_time(values=time_end(:,5))
      call add_time(time_begin(1,5),time_end(1,5),time_elapsed(5))

      ! call user diagnostics
      call date_and_time(values=time_begin(:,6))
      call diagnostics
      call date_and_time(values=time_end(:,6))
      call add_time(time_begin(1,6),time_end(1,6),time_elapsed(6))

      ! time/step increment
      time = time + dtwci
      it = it + 1

      call date_and_time(values=time_end(:,1))
      call add_time(time_begin(1,1),time_end(1,1),time_elapsed(1))

    enddo 

  end subroutine run_sim


  !---------------------------------------------------------------------
  ! close simulation and exit
  !---------------------------------------------------------------------
  subroutine close_sim

    if (myid == 0) then
      close(unit=11) ! energy.dat
      close(unit=12) ! probe.dat
      close(unit=13) ! tracking.dat
      ! close(unit=14) ! time.dat
    endif

    if (tracking_mpi) close(unit=13)
        
    if (myid==0) then
      print*, " "
      print*, " "
      print*, "*** Run completed *** "
      print*, " "
      print*, " "
      print*, "total time                        (s)          =",time_elapsed(1)
      print*, "   sub update_eta                 (s)          =",time_elapsed(2)
      print*, "   sub update_particles           (s)          =",time_elapsed(3)
      print*, "   sub sort                       (s)          =",time_elapsed(4)
      print*, "   sub update_fields              (s)          =",time_elapsed(5)
      print*, "   sub diagnostics                (s)          =",time_elapsed(6)
      print*, " "
      print*, " "
      print*, "In update_particles.parmov,"
      print*, "   sub push                       (s)          =",time_elapsed(33)
      print*, "   sub particle_boundary          (s)          =",time_elapsed(34)
      print*, "   sub moment_calculation         (s)          =",time_elapsed(35)
      print*, "   total parmov                   (s)          =",time_elapsed(31)
      print*, " "
      print*, " "
      print*, "In update_fields,"
      print*, "   sub pressgrad                  (s)          =",time_elapsed(51)
      print*, "   sub bcalc                      (s)          =",time_elapsed(52)
      print*, "   sub ecalc                      (s)          =",time_elapsed(53)
      print*, "   sub focalc                     (s)          =",time_elapsed(54)
      print*, "   total update_fields            (s)          =",time_elapsed(5)
      print*, " "
      print*, " "
      print*, "In diagnostics,"
      print*, "   sub diag_mesh                  (s)          =",time_elapsed(61)
      print*, "   sub diag_energy                (s)          =",time_elapsed(62)
      print*, "   sub diag_particle              (s)          =",time_elapsed(63)
      print*, "   sub diag_probe                 (s)          =",time_elapsed(64)
      print*, "   sub diag_tracking              (s)          =",time_elapsed(65)
      print*, "   sub write_restart              (s)          =",time_elapsed(66)
      print*, "   total diagnostics              (s)          =",time_elapsed(6)
    endif

    call MPI_FINALIZE(IERR)
    stop

  end subroutine close_sim

end program h3d