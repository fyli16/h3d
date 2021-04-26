!---------------------------------------------------------------------
!                                 H3D (V6.0)                               
!                           YURI'S NONUNIFORM MESH                         
!                           3D IMPLEMENTATION ONLY                         
!                      UNIFORM LOADING IN PHYSICAL SPACE                   
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       
!---------------------------------------------------------------------

program h3d 
  use m_init
  implicit none
  ! initialize simulation
  call init_sim
  ! execute main loops
  call sim_loops
  ! shutdown the simulation
  call shutdown
end program h3d


!---------------------------------------------------------------------
! main simulation loops
!---------------------------------------------------------------------
subroutine sim_loops
  use m_parameter
  use m_eta
  use m_particle
  use m_field
  use m_diagnostics
  implicit none 

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
    call date_and_time(values=time_begin(:,1)) ! time one-whole-loop

    ! print time & step info
    if (myid==0 .and. mod(it,n_print)==0) then
      call get_time(clock_now)
      write(6,"(A5,I7,A2,I7,A11,F8.3,A14,F8.3,A12,F8.3)") 'it = ', it, '/', itfinish, &
                    ',   time = ', time, &
                    ',   delta_t = ', real(clock_now - clock_old), &
                    ',   tot_t = ', real(clock_now - clock_init)
      clock_old = clock_now
    endif

    ! calculate resistivity (eta)
    ! (which could depend on local parameters such as current)
    call date_and_time(values=time_begin(:,2))
    if (ndim /= 1) then
      call cal_eta       
    else
      call cal_eta_2d    
    endif
    call date_and_time(values=time_end(:,2))
    call accumulate_time(time_begin(1,2),time_end(1,2),time_elapsed(2))

    ! calculate density and v's, and push particles
    call date_and_time(values=time_begin(:,3))
    ntot = 0 ! for particle tracking
    call trans
    call date_and_time(values=time_end(:,3))
    call accumulate_time(time_begin(1,3),time_end(1,3),time_elapsed(3))

    ! sort particles
    call date_and_time(values=time_begin(:,4))
    if (mod(it,n_sort) == 0) call sortit  
    call date_and_time(values=time_end(:,4))
    call accumulate_time(time_begin(1,4),time_end(1,4),time_elapsed(4))

    ! call field solver
    call date_and_time(values=time_begin(:,5))
    if (ndim /=1) then 
      call field
    else
      call field_2d
    endif     
    call date_and_time(values=time_end(:,5))
    call accumulate_time(time_begin(1,5),time_end(1,5),time_elapsed(5))

    ! inject_wave ?
    ! if (it == 21000) call inject_wave
    ! if (mod(it,100) == 0) call kick

    ! call user diagnostics
    call date_and_time(values=time_begin(:,30))
    call diagnostics
    call date_and_time(values=time_end(:,30))
    call accumulate_time(time_begin(1,30),time_end(1,30),time_elapsed(30))

    ! time/step increment
    time = time + dtwci
    it = it + 1

    call date_and_time(values=time_end(:,1))
    call accumulate_time(time_begin(1,1),time_end(1,1),time_elapsed(1))

  enddo 

end subroutine sim_loops


!---------------------------------------------------------------------
! shutdown simulation and exit
!---------------------------------------------------------------------
subroutine shutdown
  use m_parameter
  implicit none 

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
    print*, " *** Run completed *** "
    print*, " "
    print*, " "
    print*, " total time                   (s)          =",time_elapsed(1)
    print*, "   subroutine eta_cal         (s)          =",time_elapsed(2)
    print*, "   subroutine trans           (s)          =",time_elapsed(3)
    print*, "   subroutine sort            (s)          =",time_elapsed(4)
    print*, "   subroutine field           (s)          =",time_elapsed(5)
    print*, "   subroutine caltemp2        (s)          =",time_elapsed(6)
    print*, "   subroutine diagnostics     (s)          =",time_elapsed(30)
    print*, " "
    print*, " "
    print*, " In subroutine caltemp2,"
    print*, "   subroutine xreal           (s)          =",time_elapsed(24)
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(25)
    print*, "   interpolation              (s)          =",time_elapsed(26)
    print*, "   total caltemp2             (s)          =",time_elapsed(23)
    print*, " "
    print*, " "
    print*, " In subroutine trans," 
    print*, "   subroutine parmov          (s)          =",time_elapsed(7)
    print*, "   subroutine energy          (s)          =",time_elapsed(8)
    print*, "   total trans                (s)          =",time_elapsed(3)
    print*, " "
    print*, " "
    print*, " In subroutine field,"
    print*, "   subroutine pressgrad       (s)          =",time_elapsed(9)
    print*, "   subroutine bcalc           (s)          =",time_elapsed(10)
    print*, "   subroutine ecalc           (s)          =",time_elapsed(11)
    print*, "   subroutine focalc          (s)          =",time_elapsed(12)
    print*, "   total field                (s)          =",time_elapsed(5)
    print*, " "
    print*, " "
    print*, " In subroutine parmov,"
    print*, "   particle pushing           (s)          =",time_elapsed(13)
    print*, "   particle exchange          (s)          =",time_elapsed(14)
    print*, "   moment calculation         (s)          =",time_elapsed(15)
    print*, "   total parmov               (s)          =",time_elapsed(7)
    print*, " "
    print*, " "
    print*, " In subroutine bcalc,"
    print*, "   subroutine ecalc           (s)          =",time_elapsed(16)
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(17)
    print*, "   total bcalc                (s)          =",time_elapsed(22)
    print*, " "
    print*, " "
    print*, " In subroutine diagnostics,"
    print*, "   sub virtual_probes         (s)          =",time_elapsed(31)
    print*, "   sub track_particle         (s)          =",time_elapsed(32)
    print*, "   total user_diagnose        (s)          =",time_elapsed(30)
    print*, " "
    print*, " "
    print*, " In subroutine ecalc (called from subroutine bcalc and field),"
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(18)
    print*, " "
    print*, " "
    print*, "communication time / total time (%)          =", &
                100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14) &
                +time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    print*, " "
    print*, " "
    print*, "Further breakdown of communication time "
    print*, "  particle exchage in subroutine parmov (%) =", &
                100.*time_elapsed(14)/time_elapsed(1)
    print*, "  subroutine xrealbcc                   (%) =", &
                100.*(time_elapsed(25)+time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    print*, "  subroutine xreal                      (%) =", &
                100.*time_elapsed(24)/time_elapsed(1)
  endif

  call MPI_FINALIZE(IERR)
  stop

end subroutine shutdown