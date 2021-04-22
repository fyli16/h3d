!---------------------------------------------------------------------
! open history diagnostic files
subroutine open_hist_files
  use m_parameters
  implicit none 

  character (len=240) :: filename1, filename2

  if (myid == 0) then
    if (restart) then
      open(unit=11,file=trim(data_directory)//'energy.dat' ,status='old',position='append')
      ! open(unit=14,file=trim(data_directory)//'time.dat' ,status='old',position='append')
      if (.not. tracking_mpi)then
        open(unit=12,file=trim(data_directory)//'probes.dat' ,status='old',position='append')
        if (tracking_binary) then
          open(unit=13,file=trim(data_directory)//'tracking_b.dat',form='unformatted',status='old',position='append')
        else
          open(unit=13,file=trim(data_directory)//'tracking.dat',status='old',position='append')
        endif
      endif
    else
      open(unit=11,file=trim(data_directory)//'energy.dat' ,status='unknown')
      ! open(unit=14,file=trim(data_directory)//'time.dat' ,status='unknown')
      if (.not. tracking_mpi)then
        open(unit=12,file=trim(data_directory)//'probes.dat' ,status='unknown')
        if (tracking_binary) then
          open(unit=13,file=trim(data_directory)//'tracking_b.dat' ,form='unformatted',status='unknown')
        else
          open(unit=13,file=trim(data_directory)//'tracking.dat' ,status='unknown')
        endif
      endif
    endif
  endif

  if (tracking_mpi) then
    write(filename1,"(a,i4.4,a)") 'probes/probes_', myid, '.dat'
    write(filename2,"(a,i4.4,a)") 'tracking/tracking_', myid, '.dat'
    if (restart) then
      open(unit=12,file=trim(data_directory)//filename1,status='old',position='append')
      open(unit=13,file=trim(data_directory)//filename2,form='unformatted',status='old',access='append')
    else
      open(unit=12,file=trim(data_directory)//filename1,status='unknown')
      open(unit=13,file=trim(data_directory)//filename2,form='unformatted',status='unknown')
    endif
  endif

end subroutine open_hist_files


!---------------------------------------------------------------------
subroutine open_files
  use m_parameters
  implicit none
  
  integer*8 :: i, lenrec
 
  lenrec = (nxmax-2)*recl_for_single
  do i = 1, 25
    file_unit(i) = 250 + i
  enddo

  open (file_unit(1),                                                                         &
        file= trim(trim(data_directory)//'bx/bx_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(2),                                                                         &
        file= trim(trim(data_directory)//'by/by_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(3),                                                                         &
        file= trim(trim(data_directory)//'bz/bz_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(4),                                                                         &
        file= trim(trim(data_directory)//'den/den_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(5),                                                                         &
        file= trim(trim(data_directory)//'ex/ex_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(6),                                                                         &
        file= trim(trim(data_directory)//'ey/ey_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(7),                                                                         &
        file= trim(trim(data_directory)//'ez/ez_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(8),                                                                         &
        file= trim(trim(data_directory)//'vix/vix_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(9),                                                                          &
        file= trim(trim(data_directory)//'viy/viy_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(10),                                                                         &
        file= trim(trim(data_directory)//'viz/viz_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(11),                                                                          &
        file= trim(trim(data_directory)//'tpar/tpar_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(12),                                                                           &
        file= trim(trim(data_directory)//'tperp/tperp_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)

  if (eta_par == 0) then
    open (file_unit(13),                                                                         &
          file= trim(trim(data_directory)//'eta/eta_'//trim(adjustl(cycle_ascii)))//'.gda', &
          form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  else
    open (file_unit(13),                                                                             &
          file= trim(trim(data_directory)//'eta_par/eta_par_'//trim(adjustl(cycle_ascii)))//'.gda', &
          form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  endif

  open (file_unit(14),                                                                      &
        file= trim(trim(data_directory)//'p-xx/p-xx_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(15),                                                                      &
        file= trim(trim(data_directory)//'p-xy/p-xy_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(16),                                                                      &
        file= trim(trim(data_directory)//'p-xz/p-xz_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(17),                                                                      &
        file= trim(trim(data_directory)//'p-yy/p-yy_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(18),                                                                      &
        file= trim(trim(data_directory)//'p-yz/p-yz_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(19),                                                                      &
        file= trim(trim(data_directory)//'p-zz/p-zz_'//trim(adjustl(cycle_ascii)))//'.gda', &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(20),                                                                      &
        file= trim(trim(data_directory)//'fox/fox_'//trim(adjustl(cycle_ascii)))//'.gda',   &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(21),                                                                      &
        file= trim(trim(data_directory)//'foy/foy_'//trim(adjustl(cycle_ascii)))//'.gda',   &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  open (file_unit(22),                                                                      &
        file= trim(trim(data_directory)//'foz/foz_'//trim(adjustl(cycle_ascii)))//'.gda',   &
        form='unformatted',action='write',access='direct', status='unknown',recl=lenrec)
  
  return
end subroutine open_files


!---------------------------------------------------------------------
! write data
!---------------------------------------------------------------------
subroutine write_mesh_data
  use m_parameters
  implicit none  

  integer :: is
  integer*8 :: irec_num, irec_del, irec_start
  real*8 :: rnorm
  character(len=240) :: filename
  character(len=2) :: specname
 
  ! Now determine the starting and ending record number for a given variable and myid
  irec_num = 1
  irec_start = irec_num
  irec_del = nz*ny

  ! do the actual writing
  rnorm = wpiwci
  uniform_mesh = bx 
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'bx/bx_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(1),irec_start,ny,nz)
  endif

  rnorm = wpiwci
  uniform_mesh = by 
  if (MPI_IO_format) then
    filename = trim(trim(data_directory)//'by/by_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh, rnorm, trim(adjustl(filename)), irec_start, ny, nz)
  else
    call write_file_non_mpio(uniform_mesh, rnorm, file_unit(2), irec_start, ny, nz)
  endif

  rnorm = wpiwci
  uniform_mesh=bz
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'bz/bz_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(3),irec_start,ny,nz)
  endif

  rnorm = 1.
  uniform_mesh=den
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'den/den_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(4),irec_start,ny,nz)
  endif

  rnorm = wpiwci**2
  uniform_mesh=ex 
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'ex/ex_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(5),irec_start,ny,nz)
  endif

  rnorm = wpiwci**2
  uniform_mesh=ey 
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'ey/ey_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(6),irec_start,ny,nz)
  endif

  rnorm = wpiwci**2
  uniform_mesh=ez 
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'ez/ez_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(7),irec_start,ny,nz)
  endif

  rnorm = wpiwci
  uniform_mesh=vix
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'vix/vix_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(8),irec_start,ny,nz)
  endif

  rnorm = wpiwci
  uniform_mesh=viy
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'viy/viy_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(9),irec_start,ny,nz)
  endif

  rnorm = wpiwci 
  uniform_mesh=viz
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'viz/viz_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(10),irec_start,ny,nz)
  endif

  do is=1,nspec
    write(specname,'(I1,A)') is, '_'
    rnorm = 1.
    uniform_mesh=tpar(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'tpar/tpar_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(11),irec_start,ny,nz)
    endif

    rnorm = 1.
    uniform_mesh=tperp(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'tperp/tperp_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(12),irec_start,ny,nz)
    endif

    rnorm = wpiwci
    uniform_mesh=vxs(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'vxs/vxs_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(23),irec_start,ny,nz)
    endif

    rnorm = wpiwci
    uniform_mesh=vys(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'vys/vys_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(24),irec_start,ny,nz)
    endif

    rnorm = wpiwci
    uniform_mesh=vzs(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'vzs/vzs_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(25),irec_start,ny,nz)
    endif

    rnorm = one/(tx0(is)*frac(is))
    uniform_mesh=p_xx(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'p-xx/p-xx_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(14),irec_start,ny,nz)
    endif

    rnorm = one/(tx0(is)*frac(is))
    uniform_mesh=p_xy(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'p-xy/p-xy_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(15),irec_start,ny,nz)
    endif

    rnorm = one/(tx0(is)*frac(is))
    uniform_mesh=p_xz(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'p-xz/p-xz_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(16),irec_start,ny,nz)
    endif

    rnorm = one/(tx0(is)*frac(is))
    uniform_mesh=p_yy(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'p-yy/p-yy_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(17),irec_start,ny,nz)
    endif

    rnorm = one/(tx0(is)*frac(is))
    uniform_mesh=p_yz(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'p-yz/p-yz_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(18),irec_start,ny,nz)
    endif

    rnorm = one/(tx0(is)*frac(is))
    uniform_mesh=p_zz(:,:,:,is)
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'p-zz/p-zz_'//specname//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(19),irec_start,ny,nz)
    endif
  enddo

  rnorm=1.
  uniform_mesh=fox
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'fox/fox_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(20),irec_start,ny,nz)
  endif

  rnorm=1.
  uniform_mesh=foy
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'foy/foy_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(21),irec_start,ny,nz)
  endif

  rnorm=1.
  uniform_mesh=foz
  if (MPI_IO_format) then
    filename= trim(trim(data_directory)//'foz/foz_'//trim(adjustl(cycle_ascii)))//'.gda'
    call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
  else
    call write_file_non_mpio(uniform_mesh,rnorm,file_unit(22),irec_start,ny,nz)
  endif

  if (eta_par == 0) then
    rnorm = 1.
    uniform_mesh=eta
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'eta/eta_'//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(13),irec_start,ny,nz)
    endif
  else
    rnorm = 1.
    uniform_mesh=eta_times_b_dot_j
    if (MPI_IO_format) then
      filename= trim(trim(data_directory)//'eta_par/eta_par_'//trim(adjustl(cycle_ascii)))//'.gda'
      call write_file(uniform_mesh,rnorm,trim(adjustl(filename)),irec_start,ny,nz)
    else
      call write_file_non_mpio(uniform_mesh,rnorm,file_unit(13),irec_start,ny,nz)
    endif
  endif

  irec_num = irec_num + irec_del

end subroutine write_mesh_data


!---------------------------------------------------------------------
! read or write restart data
! rw = +1.0: write 
! rw = -1.0: read
!---------------------------------------------------------------------
subroutine write_read_restart_files(rw)
  use m_parameters
  implicit none

  integer*8 :: f_unit, np_count, is, ixe, iye, ize, noresete
  real :: rw
  real*8, dimension(:), allocatable :: particle_tmp_array
  integer, dimension(:), allocatable :: particle_tmp_array2 
 
  if (rw == +1.0) then  ! write restart data 
    t_stopped = t_stopped + (it - itstart + 1) * dtwci
    f_unit = 215 + myid
    open(unit=f_unit, file=trim(restart_directory)//'restfld_'//trim(adjustl(myid_char))  &
            //'.bin'//restart_index_suffix(restart_index), form='unformatted', status='unknown')
  
    ! Old restart format - dump all allocated particle memory
    ! write(f_unit) x,y,z,vx,vy,vz,qp,link,porder

    ! New restart format - dump only non-trivial particles
    write(f_unit) nspec

    do is = 1, nspec
      nptotp=0
      do ize = kb-1, ke
        do iye = jb-1, je
          do ixe = 1, nx1   
              np = IPHEAD(ixe, iye, ize, is)
              do while (np.ne.0)
                nptotp = nptotp + 1
                np = link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) nptotp

      allocate (particle_tmp_array(nptotp))
      allocate (particle_tmp_array2(nptotp))

      ! x
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = x(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! y
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = y(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! z
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = z(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! vx
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = vx(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! vy
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = vy(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! vz
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = vz(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! qp
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array(nptotp) = qp(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array

      ! ptag
      nptotp=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              np=IPHEAD(ixe,iye,ize,is)
              do while (np.ne.0)
                nptotp=nptotp+1
                particle_tmp_array2(nptotp) = ptag(np)
                np=link(np)
              enddo
          enddo
        enddo
      enddo
      write(f_unit) particle_tmp_array2

      deallocate (particle_tmp_array)
      deallocate (particle_tmp_array2)

    enddo

    write(f_unit) ninj,ninj_global,nescape,nescape_global,npart, &
    npart_global,qleft,qrite,t_stopped

    write(f_unit) x0,x1,tx0,vpar,vper,vbal,bbal,rcorr,teti,ishape

    write(f_unit) btspec, qspec, wspec, frac,                    &
    anisot, denmin, resis, wpiwci, bete, fxsho,ave1,             &
    ave2,phib, xmax,ymax,zmax,gamma,                              &
    npx, npy, npz,                                               &
    iterb,nxcel,netax,netay,nspec,   &
    nx,ny,nz,nskipx,nskipy,                                      &
    nskipz, restart, etamin, etamax, ieta, eta_par

    write(f_unit) hx,hy,hz,hxi,hyi,hzi                           &
    ,pi,efld,bfld,efluid,ethermal,eptcl,time,te0                &
    ,itrestart,iwt                                   &
    ,nx1,nx2,ny1,ny2,nz1,nz2,it                                  &
    ! ,ipstore,nptot,npleaving,npentering                &
    ,nptot,npleaving,npentering                        &
    ,iclock_speed,iopen,iseed, file_unit,file_unit_read          &
    ,clock_init,clock_old,clock_now                   &
    ,clock_time1

    write(f_unit) dfac,nskip,ipleft,iprite,ipsendleft,ipsendrite &
    ,iprecv,ipsendtop,ipsendbot,ipsendlefttop,ipsendleftbot      &
    ,ipsendritetop,ipsendritebot,ipsend

    write(f_unit) bx,by,bz,den,pe,eta,ex,ey,ez,fox,foy,foz       &
    ,eta_times_b_dot_j

    write(f_unit) vix, viy, viz, vixo, viyo, vizo

    ! write(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp
    write(f_unit) tpar,tperp,dns,vxs,vys,vzs

    ! save user data
    ! call user_data_write_restart(f_unit)
    
    close(unit=f_unit)

  else if (rw == -1.0) then ! read restart data

    f_unit = 215 + myid
    open(unit=f_unit, file=trim(restart_directory)//'restfld_'//trim(adjustl(myid_char))  &
            //'.bin'//restart_index_suffix(restart_index), form='unformatted', status='unknown')
      
    ! Old restart format - dump all allocated particle memory
    ! read(f_unit) x,y,z,vx,vy,vz,qp,link,porder

    ! New restart format - dump only non-trivial particles
    read(f_unit) nspec

    do is = 1, nspec
      read(f_unit) nptotp
      allocate (particle_tmp_array(nptotp))
      allocate (particle_tmp_array2(nptotp))

      ! x
      read(f_unit) particle_tmp_array
      ixe=2
      iye=jb
      ize=kb
      do np_count = 1, nptotp
        np = ipstore
        x(np)=particle_tmp_array(np_count)
        ipstore=link(np)
        link(np)=iphead(ixe,iye,ize,is)
        iphead(ixe,iye,ize,is)=np
      enddo
      iptemp(ixe,iye,ize,is)=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        iphead(ixe,iye,ize,is)=link(np)
        link(np)=iptemp(ixe,iye,ize,is)
        iptemp(ixe,iye,ize,is)=np
        np=iphead(ixe,iye,ize,is)
      enddo
      iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
      iptemp(ixe,iye,ize,is)=0

      ! y
      read(f_unit) particle_tmp_array
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        y(np)=particle_tmp_array(np_count)
        np=link(np)
      enddo

      ! z
      read(f_unit) particle_tmp_array
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        z(np)=particle_tmp_array(np_count)
        np=link(np)
      enddo

      ! vx
      read(f_unit) particle_tmp_array
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        vx(np)=particle_tmp_array(np_count)
        np=link(np)
      enddo

      ! vy
      read(f_unit) particle_tmp_array
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        vy(np)=particle_tmp_array(np_count)
        np=link(np)
      enddo
   
      ! vz
      read(f_unit) particle_tmp_array
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        vz(np)=particle_tmp_array(np_count)
        np=link(np)
      enddo

      ! qp
      read(f_unit) particle_tmp_array
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        qp(np)=particle_tmp_array(np_count)
        np=link(np)
      enddo

      ! ptag
      read(f_unit) particle_tmp_array2
      np_count=0
      np=iphead(ixe,iye,ize,is)
      do while (np.ne.0)
        np_count=np_count+1
        ptag(np)=particle_tmp_array2(np_count)
        np=link(np)
      enddo

      deallocate (particle_tmp_array)
      deallocate (particle_tmp_array2)

    enddo

    read(f_unit) ninj,ninj_global,nescape,nescape_global,npart,  &
    npart_global,qleft,qrite,t_stopped

    read(f_unit) x0,x1,tx0,vpar,vper,vbal,bbal,rcorr,teti,ishape

    read(f_unit) btspec, qspec, wspec, frac,       &
    anisot, denmin, resis, wpiwci, bete, fxsho,ave1,      &
    ave2,phib, xmax,ymax,zmax,gamma,                  &
    npx, npy, npz,                                               &
    iterb,nxcel,netax,netay,nspec,   &
    nx,ny,nz,nskipx,nskipy,                                      &
    nskipz, restart,etamin,etamax,ieta,eta_par

    read(f_unit) hx,hy,hz,hxi,hyi,hzi                            &
    ,efld,bfld,efluid,ethermal,eptcl,time,te0                                        &
    ,itrestart,iwt                     &
    ,nx1,nx2,ny1,ny2,nz1,nz2,it                                   &
    ! ,ipstore,nptot,npleaving,npentering                 &
    ,nptot,npleaving,npentering                 &
    ,iclock_speed,iopen,iseed, file_unit,file_unit_read           &
    ,clock_init,clock_old,clock_now                    &
    ,clock_time1

    read(f_unit) dfac,nskip,ipleft,iprite,ipsendleft,ipsendrite  &
    ,iprecv,ipsendtop,ipsendbot,ipsendlefttop,ipsendleftbot      &
    ,ipsendritetop,ipsendritebot,ipsend

    read(f_unit) bx,by,bz,den,pe,eta,ex,ey,ez,fox,foy,foz        &
    ,eta_times_b_dot_j

    read(f_unit) vix, viy, viz, vixo, viyo, vizo

    ! read(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp
    read(f_unit) tpar,tperp,dns,vxs,vys,vzs

    ! save user data
    ! call user_diagnostics_restart(f_unit)        

    close(unit=f_unit)

    ! Reset electic field and fox,y,z 
    noresete = 1
    if (noresete == 0) then
      deno=den; vixo=vix; viyo=viy; vizo=viz
      call ecalc( 1 )
      call focalc
    endif

  endif
    
  call sortit  ! see 'utils.f90'

  return
end subroutine write_read_restart_files


!---------------------------------------------------------------------
subroutine write_particle_in_volume
  use m_parameters
  use m_mesh
  implicit none
      
  real*8 :: fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8
  real*8 :: foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8
  real*8 :: foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8
  real*8 :: arb_x,arb_y,arb_z,btota,bxa,bya,bza
  real*8 :: perp1_x,perp1_y,perp1_z,perp2_x,perp2_y,perp2_z
  real*8 :: par_x,par_y,par_z
  integer :: ierror

  integer*8 :: count_kbq
  integer*8 :: nptotp_kbq,npart_kbq(2),np_ijk,Storage_Error_p,Storage_Error
  data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
  data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
  data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
  integer*8 :: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
              ,ixep1,iyep1,izep1,ixp1,iyp1,izp1,is
  real*8 :: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,xpart,ypart,zpart,r_particle
  real*8 :: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
  real*8 :: v_limit,eps2,rx0,ry0,rz0,rrat,sqrr,outer_radius,myranf,fluxran,vxa,vya,vza
  INTEGER*8 :: L, EXIT_CODE_P, EXIT_CODE
  integer*8 :: n_fast_removed,n_fast_removed_local,Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p
  real*8 :: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min,x_disp,y_disp,z_disp          &
                    ,y_disp_max_p,x_disp_max_p,z_disp_max_p,y_disp_max,x_disp_max,z_disp_max
  real*8 :: disp_max_p(3),disp_max(3),tx,ty,tz,v_x,v_y,v_z  
  INTEGER*4 :: nescapearr(8),nescapearr_global(8)
  INTEGER*4 :: ppacket(3),ppacketg(3),dpacket(4),dpacketg(4)
  INTEGER*8 :: epacket(2),epacketg(2),indx,loop
  INTEGER*8, dimension(:), allocatable :: nparr
  real*8, dimension(3,nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bxyz_av
  real*8 :: TEX1,TEX2,TEX3,TEX4,TEX5,TEX6,TEX7,TEX8  
  real*8 :: TEY1,TEY2,TEY3,TEY4,TEY5,TEY6,TEY7,TEY8  
  real*8 :: TEZ1,TEZ2,TEZ3,TEZ4,TEZ5,TEZ6,TEZ7,TEZ8  
  real*8 :: mX_xa,mX_ta,mX_ca1,mX_ca2,mX_xb,mX_dtdx,mX_tb,mX_cb1,mX_cb2
  real*8 :: mY_xa,mY_ta,mY_ca1,mY_ca2,mY_xb,mY_dtdx,mY_tb,mY_cb1,mY_cb2
  real*8 :: mZ_xa,mZ_ta,mZ_ca1,mZ_ca2,mZ_xb,mZ_dtdx,mZ_tb,mZ_cb1,mZ_cb2
  character(len=160) :: filename
  integer :: N_IN_VOLUME
  integer :: iErr1,iErr2,file,eStrLen,subArray,stat(MPI_STATUS_SIZE),mode
  character :: eStr*(1024)
  integer :: N_TOTAL,N_MAX,RECNUM,LENREC,FILENUM,IP
  real*4, dimension(:), allocatable :: xw,yw,zw,vxw,vyw,vzw,qw,vwpar,vwperp1,vwperp2

  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt
 
  N_IN_VOLUME = 0
  do IS=1,NSPEC
    nptotp=0
    do ize=kb-1,ke
      do iye=jb-1,je
        do ixe=1,nx1   
            np=IPHEAD(ixe,iye,ize,is)
            do while (np.ne.0)
              nptotp=nptotp+1
              if ( X(np) <= XBOX_R .and. X(np) >= XBOX_L .and.         &
                    Y(np) <= YBOX_R .and. Y(np) >= YBOX_L .and.         &
                    Z(np) <= ZBOX_R .and. Z(np) >= ZBOX_L                 ) N_IN_VOLUME = N_IN_VOLUME + 1
              np=link(np)
            enddo
        enddo
      enddo
    enddo
  enddo

  call MPI_ALLREDUCE(N_IN_VOLUME,N_TOTAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(N_IN_VOLUME,N_MAX  ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
  allocate ( XW(N_MAX), YW(N_MAX), ZW(N_MAX), VXW(N_MAX), VYW(N_MAX), VZW(N_MAX), QW(N_MAX)       &
            ,VWPAR(N_MAX),VWPERP1(N_MAX),VWPERP2(N_MAX))

  N_IN_VOLUME = 0
  do IS=1,NSPEC
    nptotp=0
    do ize=kb-1,ke
      do iye=jb-1,je
        do ixe=1,nx1   
          np=IPHEAD(ixe,iye,ize,is)
          L=np
          do while (np.ne.0)
            nptotp=nptotp+1
            if ( X(np) <= XBOX_R .and. X(np) >= XBOX_L .and.         &
                  Y(np) <= YBOX_R .and. Y(np) >= YBOX_L .and.         &
                  Z(np) <= ZBOX_R .and. Z(np) >= ZBOX_L                 ) then
              N_IN_VOLUME = N_IN_VOLUME + 1
              XW(N_IN_VOLUME)  = X(np)
              YW(N_IN_VOLUME)  = Y(np)
              ZW(N_IN_VOLUME)  = Z(np)
              VXW(N_IN_VOLUME) = VX(np)
              VYW(N_IN_VOLUME) = VY(np)
              VZW(N_IN_VOLUME) = VZ(np)
              QW(N_IN_VOLUME)  = QP(np)

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
              iy=iy-1             ! integer index in y direction starts at 0
              iz=iz-1             ! integer index in z direction starts at 0

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

              vxa=vx(l)
              vya=vy(l)
              vza=vz(l)

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

              par_x   = bxa / btota
              par_y   = bya / btota
              par_z   = bza / btota
              arb_x   = 0.
              arb_y   = 0.
              arb_z   = 1.
              perp1_x = arb_y*par_z   - arb_z*par_y
              perp1_y = arb_z*par_x   - arb_x*par_z
              perp1_z = arb_x*par_y   - arb_y*par_x
              perp2_x = perp1_y*par_z - perp1_z*par_y
              perp2_y = perp1_z*par_x - perp1_x*par_z
              perp2_z = perp1_x*par_y - perp1_y*par_x

              vwpar  (N_IN_VOLUME)  = vxa*par_x   + vya*par_y   + vza*par_z
              vwperp1(N_IN_VOLUME)  = vxa*perp1_x + vya*perp1_y + vza*perp1_z
              vwperp2(N_IN_VOLUME)  = vxa*perp2_x + vya*perp2_y + vza*perp2_z
            endif
          np=link(np)
          enddo
        enddo
      enddo
    enddo
  enddo
      
  lenrec=10*recl_for_single
  if (myid == 0) then
      filenum = 101
      open (filenum,                                                       &
            file= trim(data_directory)//'particle/particle_'//trim(adjustl(cycle_ascii))//'.bin',                            &
            form='unformatted', action='write',access='direct', status='unknown',recl=lenrec)
      recnum = 1
      write(filenum,rec=recnum) N_TOTAL
      recnum = recnum + 1
      if (N_IN_VOLUME /= 0) then
        write(filenum,rec=recnum) N_IN_VOLUME
        recnum = recnum + 1
        do IP=1,N_IN_VOLUME
          write(filenum,rec=recnum) XW(IP),YW(IP),ZW(IP),VXW(IP),VYW(IP),VZW(IP),QW(IP)                &
                                    ,VWPAR(IP),VWPERP1(IP),VWPERP2(IP)
          recnum = recnum + 1
        enddo
      endif
  endif

  if (myid == 0) then
    do ipe=1,nprocs-1
        call MPI_RECV(N_IN_VOLUME,1,MPI_INTEGER,ipe,ipe,MPI_COMM_WORLD,status,ierror)
        if (N_IN_VOLUME /= 0) then
          call MPI_RECV(XW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(YW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(ZW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(VXW    ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(VYW    ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(VZW    ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(QW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(VWPAR  ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(VWPERP1,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          call MPI_RECV(VWPERP2,N_IN_VOLUME,MPI_REAL4,ipe,ipe+nprocs,MPI_COMM_WORLD,status,ierror)
          write(filenum,rec=recnum) N_IN_VOLUME
          recnum = recnum + 1
          do IP=1,N_IN_VOLUME
            write(filenum,rec=recnum) XW(IP),YW(IP),ZW(IP),VXW(IP),VYW(IP),VZW(IP),QW(IP)                &
                                    ,VWPAR(IP),VWPERP1(IP),VWPERP2(IP)
            recnum = recnum + 1
          enddo
        endif
    enddo
    CLOSE(UNIT=FILENUM)
  else
      call MPI_SEND(N_IN_VOLUME,1,MPI_INTEGER,0,MYID,MPI_COMM_WORLD,ierror)
      if (N_IN_VOLUME /= 0) then
        call MPI_SEND(XW     ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(YW     ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(ZW     ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(VXW    ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(VYW    ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(VZW    ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(QW     ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(VWPAR  ,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(VWPERP1,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
        call MPI_SEND(VWPERP2,N_IN_VOLUME,MPI_REAL4,0,myid+nprocs,MPI_COMM_WORLD,ierror)
      endif
  endif

  deallocate (XW,YW,ZW,VXW,VYW,VZW,QW,VWPAR,VWPERP1,VWPERP2)

return
end subroutine write_particle_in_volume


!---------------------------------------------------------------------
subroutine write_file(dat,rnorm,fileName,irec_start,ny1m,nz1m)
  use m_parameters
  implicit none

  integer :: num_sdat
  integer*8 :: filenum,irec_start,iry1,iry2,irz1,irz2
  real*8, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: dat
  real*4, dimension(1:nxmax-2,jb:je,kb:ke) :: stemp
  integer :: ip, iry, irz, i, j, k, recnum, ii
  integer*8 :: keg, kbg, jeg, jbg, icount,ny1m,nz1m
  real*8 :: rnorm
  integer :: WriteSubArray     ! return 0 upon failure, 1 upon success
  character :: fileName*(*)    ! file to write paricles to
  integer :: dilo,dihi,&       ! bounds of the entire array
          djlo,djhi,&
          dklo,dkhi
  integer :: ailo,aihi,&       ! bounds of our portion of the array
          ajlo,ajhi,&
          aklo,akhi
  integer*8 :: domainCells,subDomainCells8       ! number of cells in the entire array
  integer :: domainDims(3)     ! dimensions of the entire array
  integer :: subDomainCells    ! number of cells in our portion of the array
  integer :: subDomainDims(3)  ! dimensions of our portion of the array
  integer :: subDomainStart(3) ! lo corner of the bounds of our portion
  integer :: iErr1,iErr2,file,eStrLen,subArray,stat(MPI_STATUS_SIZE),mode
  integer(kind=MPI_OFFSET_KIND) DISP1
  character :: eStr*(1024)

  eStrLen=1024
  WriteSubArray=0
  icount=0
  DISP1=0

  ! Open the file, bail out if this fails.
  mode = MPI_MODE_WRONLY + MPI_MODE_CREATE
  call MPI_File_open(MPI_COMM_WORLD,fileName,mode,MPI_INFO_NULL,file,iErr1)
  if (iErr1.ne.MPI_SUCCESS) then
    call MPI_Error_string(iErr1,eStr,eStrLen,iErr2)
    write(0,*)'Error: Could not open file ',fileName
    write(0,*)eStr
    write(0,*)'Write aborted.'
    return
  endif

  !  begin by converting data to REAL*4
  !
  !  do k = kb-1, ke+1
  !    do j = jb-1,je+1
  !      do i = 1, nxmax
  !        sdat(i,j,k) = rnorm*dat(i,j,k)
  !      enddo
  !    enddo
  !  enddo

  ! print *,kb,ke,jb,je
  do k = kb, ke
    do j = jb,je
      do i = 2, nxmax-1
        stemp(i-1,j,k) = rnorm*dat(i,j,k)
      enddo
    enddo
  enddo

  dilo = 1
  dihi = nxmax-2
  ailo = 1
  aihi = nxmax-2
  djlo = jbglobal(0)
  djhi = jeglobal(nprocs-1)
  ajlo = jbglobal(myid)
  ajhi = jeglobal(myid)
  dklo = kbglobal(0)
  dkhi = keglobal(nprocs-1)
  aklo = kbglobal(myid)
  akhi = keglobal(myid)
  call BoundsToDimensions(dilo,dihi,djlo,djhi,dklo,dkhi,domainDims,domainCells)
  call BoundsToDimensions(ailo,aihi,ajlo,ajhi,aklo,akhi,subDomainDims,subDomainCells8)
  subDomainCells=subDomainCells8
  subDomainStart(1)=ailo-dilo !  convert to c-index!
  subDomainStart(2)=ajlo-djlo
  subDomainStart(3)=aklo-dklo
  call MPI_Type_create_subarray(3,domainDims,subDomainDims,subDomainStart, &
      MPI_ORDER_FORTRAN,MPI_REAL,subArray,iErr1)
  call MPI_Type_commit(subArray,iErr1)

  ! Set the file view
  call MPI_File_set_view(file,DISP1,MPI_REAL,subArray,"native",MPI_INFO_NULL,iErr1)

  ! Write
  call MPI_File_write_all(file,stemp,subDomainCells,MPI_REAL,stat,iErr1)
  call MPI_File_close(file,iErr2)
  call MPI_Type_free(subArray,iErr2)

  if (iErr1.ne.MPI_SUCCESS) then
    call MPI_Error_string(iErr1,eStr,eStrLen,iErr2)
    write(0,*)'Error: Could not write to file ',fileName
    write(0,*) eStr(1:eStrLen)
    write(0,*) 'Write aborted by rank',myid
    write(0,*)'Write aborted.'
  endif

  return
end subroutine write_file


!---------------------------------------------------------------------
subroutine write_file_non_mpio(dat,rnorm,filenum,irec_start,ny1m,nz1m) 
  use m_parameters
  implicit none

  integer :: num_sdat
  integer*8 :: filenum

  integer*8 :: irec_start,iry1,iry2,irz1,irz2
  double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: dat
  real*4, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: sdat
  integer :: ip, iry, irz, i, j, k, recnum, ii
  integer*8 :: keg, kbg, jeg, jbg, icount,ny1m,nz1m
  double precision :: rnorm

  icount=0

  ! begin by converting data to REAL*4
  do k = kb-1, ke+1
    do j = jb-1,je+1
      do i = 1, nxmax
        sdat(i,j,k) = rnorm*dat(i,j,k)
      enddo
    enddo
  enddo

  ! send data to myid = 0
  num_sdat = nxmax*(nylmax+2)*(nzlmax + 2)
  if(myid.ne.0) then
    call MPI_ISEND(sdat(1,jb-1,kb-1),num_sdat,MPI_REAL4,&
         0   ,1,MPI_COMM_WORLD,req(1),IERR)
    call MPI_WAITALL(1,req,status_array,IERR)
  else
    do ip = 0, nprocs-1
      keg = keglobal(ip)
      kbg = kbglobal(ip)
      jeg = jeglobal(ip)
      jbg = jbglobal(ip)

      if (jbg.eq.1.and.jeg.eq.ny1m) then
        iry1 = 0        +1
        iry2 = ny1m+1    -1
        goto 10
      endif

      if(jbg.eq.1) then
        iry1 = 0        +1
        iry2 = jeg
        j    = jb - 2
      else if(jeg.eq.ny1m) then
	      iry1 = jbg
	      iry2 = jeg+1    -1
	      j = jb - 1
      else
        iry1 = jbg
        iry2 = jeg
        j = jb - 1
      endif

10    continue

      if (kbg.eq.1.and.keg.eq.nz1m) then
        irz1=0          +1
        irz2=nz1m+1     -1
        goto 20
      endif

      if(kbg.eq.1) then
        irz1 = 0        +1
        irz2 = keg
        k = kb - 2
      else if(keg.eq.nz1m) then
        irz1 = kbg
        irz2 = keg+1    -1
        k = kb - 1
      else
        irz1 = kbg
        irz2 = keg
        k = kb - 1
      endif

20    continue

      if(ip.ne.0) then
        call MPI_IRECV(sdat(1,jb-1,kb-1),num_sdat,MPI_REAL4,&
             IP,1,MPI_COMM_WORLD,req(1),IERR)
        call MPI_WAITALL(1,req,status_array,IERR)
      endif
      
      do irz = irz1, irz2
	      k = k+1
	      do iry = iry1, iry2
          j = j+1 
          icount=icount+1
          recnum = irec_start + iry -1 + (irz-1)*(nymax - 2)
          write(filenum, rec=recnum) (sdat(i,iry-jbg+1,irz-kbg+1),i=1   +1,nxmax    -1)
        enddo
      enddo
    enddo

    write(6,*) " number of records written = ", icount

  endif

  return
end subroutine write_file_non_mpio


