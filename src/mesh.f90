!===============================================================
! By Yuri Omelchenko, Scibernet Inc., 2003
! This module constructs a 1-D mesh by mapping a logical
! node-centered grid, t=[0,1] to a physical node-centered grid x=[0,xl] surrounded by ghost grid points. 

! The cell-centered grid points are also computed.
! The following mapping is used (dt=1/nl):
! 0 <= t <= ta(=dt*na):
! x = xa - xa*[exp(alph1*(ta-t))-1]/[exp(alph1*ta)-1]
! ta <= t <= tb(=dt*nb):
! x = xa + dx*(t-ta)/dt, where dx = (xb-xa)/(nb-na)
! tb <= t <= 1:
! x = xb + (xl-xb)*[exp(alph2*(t-tb))-1]/[exp(alph2*(1-tb))-1]

! We define eps == exp(alph*dt)
! The first steps on both nonuniform grids are matched to be equal to dx.
! We solve nonlinear equations to find eps1 and eps2:
! eps1: (eps1-1)/(eps1**na-1) = dx/xa 
! eps2: (eps2-1)/(eps2**(nl-nb)-1) = dx/(xl-xb)
! Notice that: in order eps1>0 and eps2>0 we require:
! xa/na > dx and (xl-xb)/(nl-nb) > dx
!===============================================================

module mesh_mod
  implicit none

  type mesh
    real*8, dimension(:), pointer :: xn, xc, dxn, dxc
    integer*8 :: na, nb, nl
    real*8 :: xa, xb, xl
    real*8 :: dt, dx, dtdx, ta, tb, epsa, epsb, ca1, ca2, cb1, cb2
  end type mesh

  type meshtype
    integer*8 :: type
  end type meshtype

  type (mesh) :: meshX, meshY, meshZ 
  type (meshtype), parameter :: cell = meshtype(0)
  type (meshtype), parameter :: NODE = meshtype(1)

  interface mesh_index
    module procedure mesh_index_yuri, mesh_index_hxv
  end interface

  contains

  !---------------------------------------------------------------------
  ! initialize mesh attributes 
  subroutine mesh_init(m, xa, xb, xl, na, nb, nl)
    use findexp_mod

    type(mesh), intent(out) :: m
    real*8, intent(in) :: xa, xb, xl
    integer*8, intent(in) :: na, nb, nl
    integer*8 :: i, nbb
    ! real*8 :: findexp
     
    if((xa.ge.xb).or.(na.ge.nb)) then 
      call error_abort('mesh_init(): bad parameters --- stop!')
    endif

    m%na=na ; m%nb=nb ; m%nl=nl
    m%xa=xa ; m%xb=xb ; m%xl=xl
    nbb = nl - nb

    allocate(m%xn(nl+3))  ! -1:nl+1
    allocate(m%xc(nl+2))  ! -1:nl
    allocate(m%dxn(nl+3)) ! -1:nl+1
    allocate(m%dxc(nl+2)) ! -1:nl

    m%dt = 1./real(nl)
    m%dx = (xb-xa)/(nb-na) 
    m%dtdx = m%dt/m%dx
    m%ta = m%dt*na
    m%tb = m%dt*nb

    if(na.gt.0)  then
        m%epsa = findexp(m%dx/xa,na)
        m%ca1 = (m%epsa**na-1.)/xa
        m%ca2 = m%dt/log(m%epsa)
    else
        m%epsa = 1.
        m%ca1 = 0.
        m%ca2 = 0.
    endif

    if(nbb.gt.0) then
        m%epsb = findexp(m%dx/(xl-xb),nbb)
        m%cb1 = (m%epsb**nbb-1.)/(xl-xb)
        m%cb2 = m%dt/log(m%epsb) 
    else
        m%epsb = 1.
        m%cb1 = 0.
        m%cb2 = 0.
    endif

    do i = 0,na-1
        m%xn(i+2) = xa - (m%epsa**(na-i)-1.)/m%ca1 
        m%xc(i+2) = xa - (m%epsa**(na-i-0.5)-1.)/m%ca1 
    enddo
    do i = na,nb
        m%xn(i+2) = xa + m%dx*(i-na) 
        m%xc(i+2) = xa + m%dx*(i+0.5-na) 
    enddo
    do i = nb+1,nl
        m%xn(i+2) = xb + (m%epsb**(i-nb)-1.)/m%cb1
        m%xc(i+2) = xb + (m%epsb**(i+0.5-nb)-1.)/m%cb1
    enddo

    m%xn(2) = 0. ! correct round-off errors

    m%xn(1)    = 2.*m%xn(2)-m%xn(3)
    m%xc(1)    = 2.*m%xn(2)-m%xc(2)
    m%xn(nl+3) = 2.*m%xn(nl+2)-m%xn(nl+1)
    m%xc(nl+2) = 2.*m%xn(nl+2)-m%xc(nl+1)
    ! m%xc(1)    = 0.5*(m%xn(1)+m%xn(2))
    ! m%xc(nl+2) = 0.5*(m%xn(nl+2)+m%xn(nl+3))

    ! compute cell-based mesh sizes
    do i = 1,nl+2
        m%dxc(i) = m%xn(i+1)-m%xn(i)
    enddo

    ! compute node-based mesh sizes
    do i = 2,nl+2
        m%dxn(i) = m%xc(i)-m%xc(i-1)
    enddo
    m%dxn(1) = m%dxn(2) ! fictitous
    m%dxn(nl+3) = m%dxn(nl+2) ! fictitous
      
    return
  end subroutine mesh_init


  !---------------------------------------------------------------------
  ! destroy mesh
  !---------------------------------------------------------------------
  subroutine mesh_destruct(m)

    type(mesh), intent(inout) :: m

    if(associated(m%xn)) then
      deallocate(m%xn)  
      deallocate(m%xc)  
      deallocate(m%dxc)  
      deallocate(m%dxn)  
      nullify(m%xn)
      nullify(m%xc)
      nullify(m%dxc)
      nullify(m%dxn)
    endif

    m%na = 0 ; m%nb = 0 ; m%nl = 0

    return
  end subroutine mesh_destruct


  !---------------------------------------------------------------------
  ! transform physical coordinate to logical space 
  ! MESH_INIT() must be called prior to this call
  double precision function mesh_unmap(m,x) 

    type(mesh), intent(in) :: m
    real*8, intent(in) :: x 
    real*8 :: t

    if(x.lt.m%xa) then
      t = m%ta - m%ca2*log( m%ca1*(m%xa-x)+1. )
    else if(x.le.m%xb) then
      t = m%ta + (x-m%xa)*m%dtdx 
    else 
      t = m%tb + m%cb2*log( m%cb1*(x-m%xb)+1. ) 
    endif
    ! prevent mapping outiside of logical interval [0,1]
    if(t >= 1.) t = 1.-epsilon(real(1)) 
    if(t <= 0.) t = epsilon(real(1)) 
    mesh_unmap = t

    return
  end function mesh_unmap


  !---------------------------------------------------------------------
  ! transform physical coordinate to logical space 
  ! MESH_INIT() must be called prior to this call
  double precision function mesh_map(m,t)

    type(mesh), intent(in) :: m
    real*8, intent(in) :: t 
    real*8 :: x

    if(t.lt.m%ta) then
      x = m%xa-(exp((m%ta-t)/m%ca2)-1.)/m%ca1
    else if(t.le.m%tb) then
      x = m%xa + (t-m%ta)/m%dtdx
    else 
      x = m%xb + (exp((t-m%tb)/m%cb2)-1.)/m%cb1
    endif
    ! prevent mapping outside of physical interval [0,xl]
    if(x >= m%xl) x = m%xl-epsilon(real(1)) 
    if(x <= 0.) x = epsilon(real(1)) 
    mesh_map = x
    return
  end function mesh_map


  !---------------------------------------------------------------------
  ! let xu(i=1:nnx)=[0,xl] be a uniform node-centered grid
  ! then for each xu(i) find ix(i) such as 
  ! m%xc(ix) <= xu(i) < m%xc(ix+1) if inode==cell
  ! m%xn(ix) <= xu(i) < m%xn(ix+1) if inode==NODE
  ! MESH_INIT() must be called prior to this call
  !---------------------------------------------------------------------
  subroutine mesh_index_hxv(m, inode, ix, inode_uniform)

    type(mesh), intent(in) :: m
    type(meshtype), intent(in) :: inode,inode_uniform
    integer*8, intent(out), dimension(:) :: ix
    real*8, dimension(:), pointer :: pp
    real*8 :: x,hx 
    integer*8 i,k,nnx

    ix(:)=0  
    nullify(pp)

    if(inode%type == cell%type) then ! cell-centered
        pp => m%xc 
    else
        pp => m%xn 
    endif

    nnx = size(ix)
    hx = m%xl/(nnx-1) ! nnx == number of nodes
    do i=1,nnx
      if (inode_uniform%type == cell%type) then     ! Interpolate locations of cell uniform mesh
        x=(i-1-0.5)*hx
      else                                      ! Interpolate locations of NODE uniform mesh
        x=(i-1    )*hx
      endif
      do k=2,size(pp)
          if(x.lt.pp(k)) then
            ix(i)=k-1
            exit
          endif
      enddo
    enddo

    nullify(pp)

  end subroutine mesh_index_hxv


  !---------------------------------------------------------------------
  ! let xu(i=1:nnx)=[0,xl] be a uniform node-centered grid
  ! then for each xu(i) find ix(i) such as 
  ! m%xc(ix) <= xu(i) < m%xc(ix+1) if inode==cell
  ! m%xn(ix) <= xu(i) < m%xn(ix+1) if inode==NODE
  ! MESH_INIT() must be called prior to this call
  !---------------------------------------------------------------------
  subroutine mesh_index_yuri(m, inode, ix)

    type(mesh), intent(in) :: m
    type(meshtype), intent(in) :: inode
    integer*8, intent(out), dimension(:) :: ix
    real*8, dimension(:), pointer :: pp
    real*8 :: x,hx 
    integer*8 i,k,nnx

    ix(:)=0  
    nullify(pp)

    if(inode%type == cell%type) then ! cell-centered
        pp => m%xc 
    else
        pp => m%xn 
    endif

    nnx = size(ix)
    hx = m%xl/(nnx-1) ! nnx == number of nodes
    do i=1,nnx
      x=(i-1    )*hx
      do k=2,size(pp)
          if(x.lt.pp(k)) then
            ix(i)=k-1
            exit
          endif
      enddo
    enddo

    nullify(pp)
  end subroutine mesh_index_yuri


  !---------------------------------------------------------------------
  ! interpolate (1:nx2,1:ny2) arrays from a nonuniform cell-centered grid 
  ! to a uniform node-centered grid (1:nnx,1:nny) 
  subroutine mesh_interpolated2d(inode,mx,my,a,ncx,ncy,ap,nnx,nny,lfirst)

    type(meshtype), intent(in) :: inode
    type(mesh), intent(in) :: mx,my 
    integer*8, intent(in) :: ncx,ncy,nnx,nny
    real*8, dimension(ncx,ncy), intent(in) :: a
    real*8, dimension(nnx,nny), intent(out) :: ap
    logical, intent(in) :: lfirst
    ! integer*8, dimension(nnx), save :: ix  ! won't work under Windows
    ! integer*8, dimension(nny), save :: iy  ! won't work under Windows
    integer*8, allocatable, save :: ix(:)
    integer*8, allocatable, save :: iy(:)

    real*8, save :: dx,dy
    real*8, dimension(:), pointer :: ppx,ppy
    integer*8 i,j,ii,jj
    real*8 :: rx,ry,fx,fy,w1,w2,w3,w4

    ap(:,:) = 0 ! zero output array 
    nullify(ppx)
    nullify(ppy)

    if(inode%type == cell%type) then
      ppx => mx%xc
      ppy => my%xc
    else
      ppx => mx%xn
      ppy => my%xn
    endif

    if(lfirst) then
        if(allocated(ix)) then
          deallocate(ix)
          deallocate(iy)
        endif 
        allocate(ix(nnx))
        allocate(iy(nny))
        call mesh_index(mx,inode,ix)
        call mesh_index(my,inode,iy)
        dx = mx%xl/(nnx-1)
        dy = my%xl/(nny-1)
    endif

    do j=1,nny
        ry=(j-1)*dy
        jj=iy(j)
        fy=( ry-ppy(jj) )/( ppy(jj+1)-ppy(jj) )
        do i=1,nnx
          rx=(i-1)*dx
          ii=ix(i)
          fx=( rx-ppx(ii) )/( ppx(ii+1)-ppx(ii) )
          w1=(1.-fx)*(1.-fy)
          w2=fx*(1.-fy)
          w3=(1.-fx)*fy
          w4=fx*fy
          ap(i,j)=w1*a(ii,jj)+w2*a(ii+1,jj)+w3*a(ii,jj+1)+w4*a(ii+1,jj+1)
        enddo
    enddo

    nullify(ppx)
    nullify(ppy)

  end subroutine mesh_interpolated2d


  !---------------------------------------------------------------------
  ! subroutine mesh_interpolated_3d(nonuniform_mesh, uniform_mesh, nonuniform_mesh_global)
  !   use parameter_mod

  !   real*8 :: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,w1,w2,w3,w4,w5,w6,w7,w8
  !   integer*8:: ix,iy,iz,ixp1,iyp1,izp1,i,j,k,jmin,jmax,kmin,kmax
  !   real*8, dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(in) :: nonuniform_mesh
  !   real*8, dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(out) :: uniform_mesh
  !   real*8, dimension(nxmax,0:ny+1,0:nz+1), intent(out):: nonuniform_mesh_global
  !   real*8, dimension(nxmax,0:ny+1,0:nz+1):: nonuniform_mesh_local
  !   real*8 :: xc_uniform_pos,yc_uniform_pos,zc_uniform_pos

  !   dtxi = 1./meshX%dt
  !   dtyi = 1./meshY%dt
  !   dtzi = 1./meshZ%dt

  !   uniform_mesh          = 0.
  !   nonuniform_mesh_local = 0.

  !   if (jb == 1) then
  !     jmin = 0
  !   else 
  !     jmin = jb
  !   endif

  !   if (je == ny) then
  !     jmax = ny+1
  !   else 
  !     jmax = je
  !   endif

  !   if (kb == 1) then
  !     kmin = 0
  !   else 
  !     kmin = kb
  !   endif
    
  !   if (ke == nz) then
  !     kmax = nz+1
  !   else 
  !     kmax = ke
  !   endif

  !   do i=1,nxmax
  !     do j=jmin,jmax
  !       do k=kmin,kmax
  !         nonuniform_mesh_local(i,j,k)=nonuniform_mesh(i,j,k)
  !       enddo
  !     enddo
  !   enddo

  !   call MPI_ALLREDUCE( nonuniform_mesh_local,nonuniform_mesh_global,size(nonuniform_mesh_local)  &
  !                       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

  !   do i=2,nx+1
  !     xc_uniform_pos = (i-1.5)*hx
  !     rx   = dtxi*MESH_UNMAP(meshX,xc_uniform_pos)+1.50000000000d+00
  !     ix   = rx
  !     fx   = rx-ix
  !     ixp1 = ix+1
  !     do j=jb,je
  !       yc_uniform_pos = (j-0.5)*hy
  !       ry   = dtyi*MESH_UNMAP(meshY,yc_uniform_pos)+1.50000000000d+00
  !       iy   = ry
  !       fy   = ry-iy
  !       iy   = iy-1             ! integer index in y direction starts at 0
  !       iyp1 = iy+1
  !       do k=kb,ke
  !         zc_uniform_pos = (k-0.5)*hz
  !         rz   = dtzi*MESH_UNMAP(meshZ,zc_uniform_pos)+1.50000000000d+00
  !         iz   = rz
  !         fz   = rz-iz
  !         iz   = iz-1             ! integer index in z direction starts at 0
  !         izp1 = iz+1
      
  !         w1=(1.-fx)*(1.-fy)*(1.-fz)
  !         w2=fx     *(1.-fy)*(1.-fz)
  !         w3=(1.-fx)*fy     *(1.-fz)
  !         w4=fx     *fy     *(1.-fz)
  !         w5=(1.-fx)*(1.-fy)*fz
  !         w6=fx     *(1.-fy)*fz
  !         w7=(1.-fx)*fy     *fz
  !         w8=fx     *fy     *fz
 
  !         uniform_mesh(i,j,k) =  w1 * nonuniform_mesh_global(ix  ,iy  ,iz  )     &
  !                               + w2 * nonuniform_mesh_global(ixp1,iy  ,iz  )     &
  !                               + w3 * nonuniform_mesh_global(ix  ,iyp1,iz  )     &
  !                               + w4 * nonuniform_mesh_global(ixp1,iyp1,iz  )     &
  !                               + w5 * nonuniform_mesh_global(ix  ,iy  ,izp1)     &
  !                               + w6 * nonuniform_mesh_global(ixp1,iy  ,izp1)     &
  !                               + w7 * nonuniform_mesh_global(ix  ,iyp1,izp1)     &
  !                               + w8 * nonuniform_mesh_global(ixp1,iyp1,izp1)
  !       enddo
  !     enddo
  !   enddo

  !   return
  ! end subroutine mesh_interpolated_3d


  !---------------------------------------------------------------------
  ! set up uniform mesh
  !---------------------------------------------------------------------
  subroutine setup_mesh
    use parameter_mod

    integer :: i

    if (myid==0) then
      print*, " "
      print*, "Setting up mesh ..."
    endif

    ! Initialize uniform mesh
    ! where meshX, meshY, meshZ are declared in 'mesh_mod'
    call mesh_init(meshX,xaa,xbb,xmax,nax,nbx,nx) ! initialize x-mesh
    call mesh_init(meshY,yaa,ybb,ymax,nay,nby,ny) ! initialize y-mesh
    call mesh_init(meshZ,zaa,zbb,zmax,naz,nbz,nz) ! initialize z-mesh

    ! mesh_index_yuri
    call mesh_index(meshX,CELL,ixv_2_c_map)
    call mesh_index(meshY,CELL,iyv_2_c_map)
    call mesh_index(meshZ,CELL,izv_2_c_map)
    call mesh_index(meshX,NODE,ixv_2_v_map)
    call mesh_index(meshY,NODE,iyv_2_v_map)
    call mesh_index(meshZ,NODE,izv_2_v_map)
    ! mesh_index_hxv
    call mesh_index(meshX,CELL,ixc_2_c_map,CELL)
    call mesh_index(meshY,CELL,iyc_2_c_map,CELL)
    call mesh_index(meshZ,CELL,izc_2_c_map,CELL)
    call mesh_index(meshX,NODE,ixc_2_v_map,CELL)
    call mesh_index(meshY,NODE,iyc_2_v_map,CELL)
    call mesh_index(meshZ,NODE,izc_2_v_map,CELL)

    ! write mesh properties into a file
    if (myid == 0) then
      open(unit=10, file='mesh_vertices.dat', status='unknown', form='formatted')
      write(10,*) meshX%nl+1, meshY%nl+1, meshZ%nl+1

      do i = 2, meshX%nl+2 
        write(10,*) meshX%xn(i)
      enddo
      do i = 2, meshX%nl+2 
        write(10,*) meshX%dxc(i)
      enddo

      do i = 2, meshY%nl+2 
        write(10,*) meshY%xn(i)
      enddo
      do i = 2, meshY%nl+2 
        write(10,*) meshY%dxc(i)
      enddo

      do i = 2, meshZ%nl+2 
        write(10,*) meshZ%xn(i)
      enddo
      do i = 2, meshZ%nl+2 
        write(10,*) meshZ%dxc(i)
      enddo
      
      close(unit=10)
    endif

  end subroutine setup_mesh

end module mesh_mod