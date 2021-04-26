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
subroutine accumulate_time(tbegin, tend, telapsed)
  integer, dimension(8) :: tbegin, tend
  real*8 :: telapsed
  
  telapsed = telapsed &
      + (tend(3)-tbegin(3))*3600.*24. &
      + (tend(5)-tbegin(5))*3600. &
      + (tend(6)-tbegin(6))*60. &
      + (tend(7)-tbegin(7)) &
      + (tend(8)-tbegin(8))*0.001

  return
end subroutine accumulate_time


!---------------------------------------------------------------------
! smoothing routine for periodic B.C.
! 3D version of 3-point binomial smoothing
!            y(i)=(x(i-1)+2*x(i)+x(i+1))/4
! i.e. 27 points are involved
!---------------------------------------------------------------------
subroutine nsmth (a)
  use m_parameter

  integer*8 :: i,j,k
  real*8, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: temp, a

  ! copy input array "a" to "temp" including ghost cells
  do k=kb-1,ke+1
      do j = jb-1,je+1
        do i=1,nx2
            temp(i,j,k)=a(i,j,k)
        enddo
      enddo
  enddo

  ! smoothing only for inner cells (exclude ghost cells)
  do k = kb, ke
      do j = jb, je
        do i = 2, nx1
          a(i,j,k)=temp(i,j,k)/8. &
            + ( temp(i-1,j,k)+temp(i+1,j,k)+temp(i,j+1,k)+temp(i,j-1,k) &
            + temp(i,j,k+1)+temp(i,j,k-1))/16. &
            + ( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k) &
            + temp(i-1,j-1,k)+temp(i,j+1,k+1)+temp(i,j-1,k+1) &
            + temp(i,j+1,k-1)+temp(i,j-1,k-1) &
            + temp(i+1,j,k+1)+temp(i-1,j,k+1)+temp(i+1,j,k-1) &
            + temp(i-1,j,k-1))/32. &
            +( temp(i+1,j+1,k+1)+temp(i-1,j+1,k+1) &
            + temp(i+1,j-1,k+1)+temp(i-1,j-1,k+1) &
            + temp(i+1,j+1,k-1)+temp(i-1,j+1,k-1) &
            + temp(i+1,j-1,k-1)+temp(i-1,j-1,k-1)) / 64.
        enddo
      enddo
  enddo

  ! apply periodic BCs 
  call XREALBCC(a,0_8,NX,NY,NZ)

  return
end subroutine nsmth


!---------------------------------------------------------------------
subroutine nsmth_2d (a,nx2m,ny2m,nz2m)
  use m_parameter

  integer*8 :: i,j,k
  integer*8 :: nx2m, ny2m, nz2m
  real*8, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: temp, a

  ! smoothing routine--assumes aperiodic in x
  call xrealbcc_2d(a,0_8,NX,NY,NZ)
  temp=a

  do k=kb-1,ke+1
    do j = jb,je
      do i=2,nx1
        a(i,j,k)=temp(i,j,k)/4.&
                  +( temp(i-1,j  ,k)+temp(i+1,j  ,k)+temp(i  ,j+1,k)   &
                  +temp(i  ,j-1,k))/8.&
                  +( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k)   & 
                  +temp(i-1,j-1,k))/16.
      enddo
    enddo
  enddo

  do k=kb-1,ke+1
      do j = jb-1,je+1
        a(1  ,j,k)=a(nx1  ,j,k)
        a(nx2,j,k)=a(2,j,k)
      enddo
  enddo

  call xrealbcc_2d(a,0_8,NX,NY,NZ)

  return
end subroutine nsmth_2d


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


!---------------------------------------------------------------------
! make list?
!---------------------------------------------------------------------
subroutine makelist
  use m_parameter

  integer*8:: ip

  ipstore = 1 ! integer scalar
  ipleft = 0 ! ipleft(nspec)
  iprite = 0
  iprecv = 0
  iphead = 0 ! iphead(nxmax,jb-1:je+1,kb-1:ke+1,nspec)
  iptemp = 0

  do ip = 1, nplmax-1
      link(ip) = ip+1  ! link(nplmax) as defined in m_parameter
  enddo
  link(nplmax)=0 ! why linking in this strange way?

  return
end subroutine makelist


!-----------------------------------------------------------------
! domain decompostion util: splits n elements between nprocs processors
! @param n: number of elements
! @param nprocs: number of processors
! @param b, e : begin/end indices
!-----------------------------------------------------------------
subroutine mpe_decomp1d(n, nprocs, myid, b, e)
  integer ::   nprocs, myid
  integer*8 :: n, b, e, nlocal, deficit

  nlocal  = n / nprocs
  b       = myid * nlocal + 1
  deficit = mod(n, int(nprocs,8) )
  b       = b + min( int(myid,8) ,deficit)
  if (myid  <  deficit) then
      nlocal = nlocal + 1
  endif
  e = b + nlocal - 1
  if (e > n .or. myid == nprocs-1) e = n

  return
end subroutine mpe_decomp1d


!---------------------------------------------------------------------
! error abort
!---------------------------------------------------------------------
subroutine error_abort(message)
  character(*), intent(in) :: message
  write(6,*) message
  stop
end subroutine error_abort


!---------------------------------------------------------------------
! warning 
!---------------------------------------------------------------------
subroutine warning(message)
  character(*), intent(in) :: message
  write(6,*) message
end subroutine warning