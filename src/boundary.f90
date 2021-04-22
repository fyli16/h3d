!---------------------------------------------------------------------
! subroutine 'xreal' exchanges the buffer (ghost) cell information in the planes kb-1 and ke+1.
! VR: note that this routine adds contributions from the ghost cells
! VR: into the boundary cells. this is used for moments (e.g. v's and density)

!  There are up to four exchanges that take place for each processor:
!  1 & 2)  Processors myid and nbrleft exchange contributions in cells kb-1 
!          (myid) and ke+1 (nbrleft)
!  3 & 4)  Processors myid and nbrrite exchange contributions in cells ke+1 
!          (myid) and kb-1 (nbrleft)
! Fewer exchanges occur if myid is either at the left or right boundary of the physical domain

! Blocking is avoided by first sending to the right and then to the left
!---------------------------------------------------------------------
subroutine xreal(a, nx1m, ny1m, nz1m)
  use m_parameters
  implicit none

  integer*8 :: i, j, nx1m, ny1m, nz1m, k
  real*8 :: a(nxmax, jb-1:je+1, kb-1:ke+1), tmp(nxmax, jb-1:je+1, kb-1:ke+1)

  ! VR: update x direction (periodic case)
  a(2   ,:,:)   = a(2   ,:,:)   + a(nx1m+2   ,:,:)
  a(nx1m+1,:,:) = a(nx1m+1,:,:) + a(1,:,:)

  !VR:  Note that stridey and stridez are custom MPI types.
  !VR:  i need to check their definitions, but this probably sends the whole plane
  call MPI_SENDRECV(a    (1    ,je+1,kb-1),1,stridery,nbrrite,0, &
                    tmp  (1    ,jb-1,kb-1),1,stridery,nbrleft,0, &
                    mpi_comm_world,status,ierr)

  a(:,jb,:)=a(:,jb,:)+tmp  (:,jb-1,:)
  call MPI_SENDRECV(a    (1    ,jb-1,kb-1),1,stridery,nbrleft,1, &
                    tmp  (1    ,je+1,kb-1),1,stridery,nbrrite,1, &
                    mpi_comm_world,status,ierr)

  a(:,je,:)=a(:,je,:)+tmp  (:,je+1,:)
  call MPI_SENDRECV(a    (1    ,jb-1,ke+1),1,striderz,nbrtop ,0, &
                    tmp  (1    ,jb-1,kb-1),1,striderz,nbrbot ,0, &
                    mpi_comm_world,status,ierr)

  a(:,:,kb)=a(:,:,kb)+tmp  (:,:,kb-1)
  call MPI_SENDRECV(a    (1    ,jb-1,kb-1),1,striderz,nbrbot ,1, &
                    tmp  (1    ,jb-1,ke+1),1,striderz,nbrtop ,1, &
                    mpi_comm_world,status,ierr)

  a(:,:,ke)=a(:,:,ke)+tmp  (:,:,ke+1)

  return
end subroutine xreal


!---------------------------------------------------------------------
! XREALBCC updates the ghost cells in the planes kb-1 and ke+1 by obtaining
! latest values for neighboring processors (nbrleft and nbrrite)
!
! There are up to four exchanges that take place for each processor:
! 1) cell ke+1 in nbrleft is replaced by values in cell kb  in myid
! 2) cell kb-1 in nbrrite is replaced by values in cell ke  in myid
! 3) cell ke+1 in myid is replaced by values in cell kb  in nbrrite
! 4) cell kb-1 in myid is replaced by values in cell ke  in nbrleft
!
! Fewer exchanges occur if myid is either at the left or right boundary of
! the physical domain
!
! Blocking is avoided by first sending to the right and then to the left
!
! If ibnd = 1, then zero slope boundary conditions for k = 0 and nz1 are 
! set at end of routine
!---------------------------------------------------------------------
subroutine xrealbcc(a, ibnd, nx1m, ny1m,nz1m)
  use m_parameters
  implicit none

  integer*8 :: ibnd,i,j,nx1m,ny1m,nz1m
  real*8 :: a(nxmax,jb-1:je+1,kb-1:ke+1)

  call MPI_SENDRECV(a(1    ,jb-1,ke  ),1,striderz,nbrtop ,0,&
                    a(1    ,jb-1,kb-1),1,striderz,nbrbot ,0,&
                    mpi_comm_world,status,ierr)
  call MPI_SENDRECV(a(1    ,jb-1,kb  ),1,striderz,nbrbot ,1,&
                    a(1    ,jb-1,ke+1),1,striderz,nbrtop ,1,&
                    mpi_comm_world,status,ierr)

  call MPI_SENDRECV(a(1    ,je  ,kb-1),1,stridery,nbrrite,0,&
                    a(1    ,jb-1,kb-1),1,stridery,nbrleft,0,&
                    mpi_comm_world,status,ierr)
  call MPI_SENDRECV(a(1    ,jb  ,kb-1),1,stridery,nbrleft,1,&
                    a(1    ,je+1,kb-1),1,stridery,nbrrite,1,&
                    mpi_comm_world,status,ierr)

  !VR: update ghost cells in x 
  a(1   ,:,:)   = a(nx1m+1   ,:,:)
  a(nx1m+2,:,:) = a(2,:,:)

  return
end subroutine xrealbcc


!---------------------------------------------------------------------
! VR: with perioidc B.C., this is probably equivalent to the 
! xrealbcc_pack_e subroutine
!---------------------------------------------------------------------
subroutine xrealbcc_pack_b(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)
  use m_parameters
  implicit none
  integer*8 ibnd,i,j,nx1m,ny1m,nz1m,k
  double precision a_x(nxmax,jb-1:je+1,kb-1:ke+1)&
                  ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                  ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                  ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                  ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                  ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                  ,packed_data_xy_recv(nxmax,jb-1:je+1,3)

  ! VR: pack the data on ke boundary (right z boundary)
  do j=jb-1,je+1
    do i=1,nxmax
      packed_data_xy_send(i,j,1)=a_x(i,j,ke)
      packed_data_xy_send(i,j,2)=a_y(i,j,ke)
      packed_data_xy_send(i,j,3)=a_z(i,j,ke)
    enddo
  enddo

  ! VR: send it to the top and recieve from the bottom processors
  call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrtop ,0,&
                    packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrbot ,0,&
                    mpi_comm_world,status,ierr)

  ! VR: unpack into the local array
  do j = jb-1, je+1
      do i=1,nxmax
        a_x(i,j,kb-1)=packed_data_xy_recv(i,j,1)
        a_y(i,j,kb-1)=packed_data_xy_recv(i,j,2)
        a_z(i,j,kb-1)=packed_data_xy_recv(i,j,3)
      enddo
  enddo
  
  ! VR: pack data on the kb bondary
  do j=jb-1,je+1
    do i=1,nxmax
      packed_data_xy_send(i,j,1)=a_x(i,j,kb)
      packed_data_xy_send(i,j,2)=a_y(i,j,kb)
      packed_data_xy_send(i,j,3)=a_z(i,j,kb)
    enddo
  enddo
  ! VR: send the kb data to "bottom" and recieve it from the "top" processors
  call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrbot ,1,&
                    packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrtop ,1,&
                    mpi_comm_world,status,ierr)

  ! VR: unpack into the local array
  do j=jb-1,je+1
      do i=1,nxmax
        a_x(i,j,ke+1)=packed_data_xy_recv(i,j,1)
        a_y(i,j,ke+1)=packed_data_xy_recv(i,j,2)
        a_z(i,j,ke+1)=packed_data_xy_recv(i,j,3)
      enddo
  enddo
  
  ! exhange data on the (je) boundary (right in y)
  do k=kb-1,ke+1
    do i=1,nxmax
      packed_data_xz_send(i,k,1)=a_x(i,je,k)
      packed_data_xz_send(i,k,2)=a_y(i,je,k)
      packed_data_xz_send(i,k,3)=a_z(i,je,k)
    enddo
  enddo
  call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrrite,0,&
                    packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrleft,0,&
                    mpi_comm_world,status,ierr)

  ! VR: unpack into the local array
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

  ! VR: exchange the data on the jb boundary
  do k=kb-1,ke+1
    do i=1,nxmax
      packed_data_xz_send(i,k,1)=a_x(i,jb,k)
      packed_data_xz_send(i,k,2)=a_y(i,jb,k)
      packed_data_xz_send(i,k,3)=a_z(i,jb,k)
    enddo
  enddo
  call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrleft,0,&
                    packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrrite,0,&
                    mpi_comm_world,status,ierr)

  ! VR: unpack into the local attay
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

  ! WHY NO X-EXCHANGE HERE?

  ! VR: update ghost cells in x 
  a_x(1   ,:,:)   = a_x(nx1m+1   ,:,:)
  a_x(nx1m+2,:,:) = a_x(2,:,:)
  a_y(1   ,:,:)   = a_y(nx1m+1   ,:,:)
  a_y(nx1m+2,:,:) = a_y(2,:,:)
  a_z(1   ,:,:)   = a_z(nx1m+1   ,:,:)
  a_z(nx1m+2,:,:) = a_z(2,:,:)

  return
end subroutine xrealbcc_pack_b


!---------------------------------------------------------------------
subroutine xrealbcc_pack_e(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)
  use m_parameters
  implicit none

  integer*8 :: ibnd,i,j,nx1m,ny1m,nz1m,k
  real*8 :: a_x(nxmax,jb-1:je+1,kb-1:ke+1) &
            ,a_y(nxmax,jb-1:je+1,kb-1:ke+1) &
            ,a_z(nxmax,jb-1:je+1,kb-1:ke+1) &
            ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
            ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
            ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
            ,packed_data_xy_recv(nxmax,jb-1:je+1,3)

  !VR: pack z right boundary (ke) for all x and y
  do j=jb-1,je+1
    do i=1,nxmax
      packed_data_xy_send(i,j,1)=a_x(i,j,ke)
      packed_data_xy_send(i,j,2)=a_y(i,j,ke)
      packed_data_xy_send(i,j,3)=a_z(i,j,ke)
    enddo
  enddo

  !VR send the right boundary in z to "top" processor and recieve from "bottom" processor
  call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrtop ,0,&
                    packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrbot ,0,&
                    mpi_comm_world,status,ierr)

  !VR : for periodic b.c., simply unpack into the local array
  do j=jb-1,je+1
      do i=1,nxmax
        a_x(i,j,kb-1)=packed_data_xy_recv(i,j,1)
        a_y(i,j,kb-1)=packed_data_xy_recv(i,j,2)
        a_z(i,j,kb-1)=packed_data_xy_recv(i,j,3)
      enddo
  enddo

  !VR: now pack the data on the left z boundary
  do j=jb-1,je+1
    do i=1,nxmax
      packed_data_xy_send(i,j,1)=a_x(i,j,kb)
      packed_data_xy_send(i,j,2)=a_y(i,j,kb)
      packed_data_xy_send(i,j,3)=a_z(i,j,kb)
    enddo
  enddo
  !VR: send the packed left z boundary to the "bottom" proc. and recieve from "top"
  call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrbot ,1,&
                    packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrtop ,1,&
                    mpi_comm_world,status,ierr)

  !VR: for periodic B.C., simply unpack the data
  do j=jb-1,je+1
      do i=1,nxmax
        a_x(i,j,ke+1)=packed_data_xy_recv(i,j,1)
        a_y(i,j,ke+1)=packed_data_xy_recv(i,j,2)
        a_z(i,j,ke+1)=packed_data_xy_recv(i,j,3)
      enddo
  enddo

  !VR: pack the data on the right y boundary (je)
  do k=kb-1,ke+1
    do i=1,nxmax
      packed_data_xz_send(i,k,1)=a_x(i,je,k)
      packed_data_xz_send(i,k,2)=a_y(i,je,k)
      packed_data_xz_send(i,k,3)=a_z(i,je,k)
    enddo
  enddo
  !VR: send the right y boundary to the "right" processor and recieve it from the "left"
  call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrrite,0,&
                    packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrleft,0,&
                    mpi_comm_world,status,ierr)

  !VR: for periodic B.C., simply unpack the data
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

  !VR: pack the data on the left y boundary (jb)
  do k=kb-1,ke+1
    do i=1,nxmax
      packed_data_xz_send(i,k,1)=a_x(i,jb,k)
      packed_data_xz_send(i,k,2)=a_y(i,jb,k)
      packed_data_xz_send(i,k,3)=a_z(i,jb,k)
    enddo
  enddo
  !VR: send it to the "left" and recieve from the "right" processors
  call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrleft,0,&
                    packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrrite,0,&
                    mpi_comm_world,status,ierr)

  !VR: unpack into the local array
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

  ! WHY NO-EXCHANGE?

  !VR: update ghost cells in x 
  a_x(1   ,:,:)   = a_x(nx1m+1   ,:,:)
  a_x(nx1m+2,:,:) = a_x(2,:,:)
  a_y(1   ,:,:)   = a_y(nx1m+1   ,:,:)
  a_y(nx1m+2,:,:) = a_y(2,:,:)
  a_z(1   ,:,:)   = a_z(nx1m+1   ,:,:)
  a_z(nx1m+2,:,:) = a_z(2,:,:)

  return
end subroutine xrealbcc_pack_e