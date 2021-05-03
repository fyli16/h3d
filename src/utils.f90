!---------------------------------------------------------------------
! get current time
!---------------------------------------------------------------------
subroutine get_time(tlabel)
  real*8, intent(out) :: tlabel
  integer :: now(8)

  call date_and_time(values=now)
  tlabel = now(5)*3600. + now(6)*60. + now(7) + now(8)*0.001
  return
end subroutine get_time


!---------------------------------------------------------------------
! accumulate time between two stamps
!---------------------------------------------------------------------
subroutine add_time(tbegin, tend, telapsed)
  integer, dimension(8) :: tbegin, tend
  real*8 :: telapsed
  
  telapsed = telapsed &
      + (tend(3)-tbegin(3))*3600.*24. &
      + (tend(5)-tbegin(5))*3600. &
      + (tend(6)-tbegin(6))*60. &
      + (tend(7)-tbegin(7)) &
      + (tend(8)-tbegin(8))*0.001

  return
end subroutine add_time


!---------------------------------------------------------------------
! bounds to dimensions?
!---------------------------------------------------------------------
subroutine BoundsToDimensions(ilo, ihi, jlo, jhi, klo, khi, dims, nCells)
  integer:: ilo,ihi,&    ! bounds of a 3d array
            jlo,jhi,&
            klo,khi
  integer :: dims(3)      ! dimensions of the array (out)
  integer*8 :: nCells       ! total number of cells (out)

  dims(1) = ihi-ilo+1
  dims(2) = jhi-jlo+1
  dims(3) = khi-klo+1
  nCells = dims(1)*dims(2)*dims(3)
  
  return
end subroutine BoundsToDimensions


!-----------------------------------------------------------------
! domain decompostion util: splits n elements between nprocs processors
! @param n: number of elements
! @param nprocs: number of processors
! @param b, e : begin/end indices
!-----------------------------------------------------------------
subroutine mpe_decomp1d(n, nprocs, myid, b, e)
  integer ::   nprocs, myid
  integer*8 :: n, b, e, nlocal, deficit

  nlocal  = n/nprocs
  b       = myid*nlocal + 1
  deficit = mod(n, int(nprocs,8) )
  b       = b + min( int(myid,8), deficit )
  if (myid<deficit) then
      nlocal = nlocal + 1
  endif
  e = b + nlocal - 1
  if (e>n .or. myid==nprocs-1) e = n

  return
end subroutine mpe_decomp1d


!---------------------------------------------------------------------
! warning 
!---------------------------------------------------------------------
subroutine warning(message)
  character(*), intent(in) :: message
  write(6,*) message
end subroutine warning


!---------------------------------------------------------------------
! error abort
!---------------------------------------------------------------------
subroutine error_abort(message)
  character(*), intent(in) :: message
  call warning(message)
  stop
end subroutine error_abort


