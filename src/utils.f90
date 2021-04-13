!---------------------------------------------------------------------
subroutine accumulate_time(time_begin, time_end, time_elapsed)
  implicit none

  integer, dimension(8) :: time_begin, time_end
  real*8 :: time_elapsed
  
  time_elapsed = time_elapsed &
       + (time_end(3)-time_begin(3))*3600.*24. &
       + (time_end(5)-time_begin(5))*3600. &
       + (time_end(6)-time_begin(6))*60. &
       + (time_end(7)-time_begin(7)) &
       + (time_end(8)-time_begin(8))*0.001

  return
end subroutine accumulate_time


!---------------------------------------------------------------------
! XF:  smoothing routine--for periodic B.C.
! 3D version of 3-point binomial smoothing
!            y(i)=(x(i-1)+2*x(i)+x(i+1))/4
! i.e. 27 points are involved
!---------------------------------------------------------------------
subroutine nsmth (a)
  use parameter_mod
  implicit none
  integer*8 i,j,k
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
  do k=kb,ke
      do j = jb,je
        do i=2,nx1
            a(i,j,k)=temp(i,j,k)/8.&
                +( temp(i-1,j,k)+temp(i+1,j,k)+temp(i,j+1,k)+temp(i,j-1,k)&
                +temp(i,j,k+1)+temp(i,j,k-1))/16.&
                +( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k)&
                +temp(i-1,j-1,k)&
                +temp(i,j+1,k+1)+temp(i,j-1,k+1)+temp(i,j+1,k-1)+temp(i,j-1,k-1)&
                +temp(i+1,j,k+1)+temp(i-1,j,k+1)+temp(i+1,j,k-1)&
                +temp(i-1,j,k-1))/32.&
                +( temp(i+1,j+1,k+1)+temp(i-1,j+1,k+1)&
                +temp(i+1,j-1,k+1)+temp(i-1,j-1,k+1)&
                +temp(i+1,j+1,k-1)+temp(i-1,j+1,k-1)&
                +temp(i+1,j-1,k-1)+temp(i-1,j-1,k-1))/64.
        enddo
      enddo
  enddo

  ! apply periodic BCs 
  call XREALBCC(a,0_8,NX,NY,NZ)

  return
end subroutine nsmth


!---------------------------------------------------------------------
subroutine clock_write(iunit,message,i2,i1,is,it)
  implicit none
  integer*8 iunit,i2,i1,is,it
  character*10 message
  write(iunit,"(i4,a10,e20.8)") it, real(i2-i1)/real(is)
  return
end subroutine clock_write


!---------------------------------------------------------------------
subroutine makelist
  use parameter_mod
  implicit none

  integer*8:: ip

  ipstore = 1
  ipleft = 0
  iprite = 0
  iprecv = 0
  iphead = 0
  iptemp = 0

  do ip = 1, nplmax-1
      link(ip) = ip+1
  enddo
  link(nplmax)=0

  return
end subroutine makelist


!-----------------------------------------------------------------
! domain decompostion util: splits n elements between numprocs processors
! @param n: number of elements
! @param numprocs: number of processors
! @param b, e : begin/end indices
!-----------------------------------------------------------------
subroutine mpe_decomp1d(n, numprocs, myid, b, e)
  implicit none  
  integer ::   numprocs, myid
  integer*8 :: n, b, e, nlocal, deficit

  nlocal  = n / numprocs
  b       = myid * nlocal + 1
  deficit = mod(n, int(numprocs,8) )
  b       = b + min( int(myid,8) ,deficit)
  if (myid  <  deficit) then
      nlocal = nlocal + 1
  endif
  e = b + nlocal - 1
  if (e > n .or. myid == numprocs-1) e = n

  return
end subroutine mpe_decomp1d


!---------------------------------------------------------------------
subroutine error_abort(message)
  character(*), intent(in) :: message
  write(6,*) message
  stop
end subroutine error_abort


!---------------------------------------------------------------------
subroutine warning(message)
  character(*), intent(in) :: message
  write(6,*) message
end subroutine warning