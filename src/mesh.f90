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

module m_mesh
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
  type (meshtype), parameter :: cell = meshtype(0) ! cell%type=0
  type (meshtype), parameter :: node = meshtype(1) ! node%type=1

  interface mesh_index
    module procedure mesh_index_yuri, mesh_index_hxv
  end interface

  contains

  !---------------------------------------------------------------------
  ! initialize mesh attributes along one dimension
  !---------------------------------------------------------------------
  subroutine mesh_init_1d(m, xa, xb, xl, na, nb, nl)
    use m_functions
    use m_parameter, only : myid

    type(mesh), intent(out) :: m
    real*8, intent(in) :: xa, xb, xl
    integer*8, intent(in) :: na, nb, nl
    integer*8 :: i, nbb
     
    if( (xa.ge.xb).or.(na.ge.nb) ) then 
      call error_abort('mesh_init_1d(): bad parameters --- stop!')
    endif

    m%na=na; m%nb=nb ; m%nl=nl
    m%xa=xa; m%xb=xb ; m%xl=xl
    nbb = nl - nb  ! -na?

    allocate(m%xn(nl+3))  ! -1:nl+1
    allocate(m%xc(nl+2))  ! -1:nl
    allocate(m%dxn(nl+3)) ! -1:nl+1
    allocate(m%dxc(nl+2)) ! -1:nl

    m%dt = 1./real(nl)  ! nb - na = nl
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

    ! for na=0, this block is not executed
    do i = 0,na-1
        m%xn(i+2) = xa - (m%epsa**(na-i)-1.)/m%ca1 
        m%xc(i+2) = xa - (m%epsa**(na-i-0.5)-1.)/m%ca1 
    enddo

    ! for na=0, nb=nl, only this block is executed
    do i = na,nb
        m%xn(i+2) = xa + m%dx*(i-na) ! coordinates at nodes
        m%xc(i+2) = xa + m%dx*(i+0.5-na) ! coordinates at half cells
    enddo

    ! for nb=nl, this block is not executed
    do i = nb+1,nl
        m%xn(i+2) = xb + (m%epsb**(i-nb)-1.)/m%cb1
        m%xc(i+2) = xb + (m%epsb**(i+0.5-nb)-1.)/m%cb1
    enddo

    m%xn(2) = 0. ! correct round-off errors

    m%xn(1)    = 2.*m%xn(2) - m%xn(3)
    m%xc(1)    = 2.*m%xn(2) - m%xc(2)
    m%xn(nl+3) = 2.*m%xn(nl+2) - m%xn(nl+1)
    m%xc(nl+2) = 2.*m%xn(nl+2) - m%xc(nl+1)

    ! compute cell-based mesh cell sizes
    do i = 1,nl+2
        m%dxc(i) = m%xn(i+1)-m%xn(i)
    enddo

    ! compute node-based mesh cell sizes
    do i = 2,nl+2
        m%dxn(i) = m%xc(i)-m%xc(i-1)
    enddo
    m%dxn(1) = m%dxn(2) 
    m%dxn(nl+3) = m%dxn(nl+2) 
      
    return
  end subroutine mesh_init_1d


  !---------------------------------------------------------------------
  ! let xu(i=1:nnx)=[0,xl] be a uniform node-centered grid
  ! then for each xu(i) find ix(i) such as 
  ! m%xc(ix) <= xu(i) < m%xc(ix+1) if inode==cell
  ! m%xn(ix) <= xu(i) < m%xn(ix+1) if inode==node
  ! mesh_init_1d() must be called prior to this call
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
    do i = 1, nnx
      x = (i-1)*hx
      do k = 2, size(pp)
          if(x.lt.pp(k)) then
            ix(i) = k-1
            exit
          endif
      enddo
    enddo

    nullify(pp)
  end subroutine mesh_index_yuri


  !---------------------------------------------------------------------
  ! let xu(i=1:nnx)=[0,xl] be a uniform node-centered grid
  ! then for each xu(i) find ix(i) such as 
  ! m%xc(ix) <= xu(i) < m%xc(ix+1) if inode==cell
  ! m%xn(ix) <= xu(i) < m%xn(ix+1) if inode==node
  ! mesh_init_1d() must be called prior to this call
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
      else                                      ! Interpolate locations of node uniform mesh
        x=(i-1)*hx
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
  ! transform physical coordinate to logical space 
  ! mesh_init_1d() must be called prior to this call
  !---------------------------------------------------------------------
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
  ! transform physical coordinate to logical space 
  ! mesh_init_1d() must be called prior to this call
  !---------------------------------------------------------------------
  double precision function mesh_unmap(m, x) 
    type(mesh), intent(in) :: m
    real*8, intent(in) :: x 
    real*8 :: t

    if(x < m%xa) then
      t = m%ta - m%ca2*log( m%ca1*(m%xa-x)+1. )
    else if(x <= m%xb) then
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
  ! initialize mesh
  !---------------------------------------------------------------------
  subroutine init_mesh
    use m_parameter

    integer :: is, ixe, iye, ize

    if (myid==0) then
      print*, " "
      print*, "Setting up mesh"
      print*, "-------------------------------------------------"
    endif

    ! Initialize uniform mesh
    call mesh_init_1d(meshX,xaa,xbb,xmax,nax,nbx,nx) ! initialize x-mesh
    call mesh_init_1d(meshY,yaa,ybb,ymax,nay,nby,ny) ! initialize y-mesh
    call mesh_init_1d(meshZ,zaa,zbb,zmax,naz,nbz,nz) ! initialize z-mesh

    ! mesh_index_yuri
    call mesh_index(meshX,cell,ixv_2_c_map)
    call mesh_index(meshY,cell,iyv_2_c_map)
    call mesh_index(meshZ,cell,izv_2_c_map)
    call mesh_index(meshX,node,ixv_2_v_map)
    call mesh_index(meshY,node,iyv_2_v_map)
    call mesh_index(meshZ,node,izv_2_v_map)

    ! mesh_index_hxv
    call mesh_index(meshX,cell,ixc_2_c_map,cell)
    call mesh_index(meshY,cell,iyc_2_c_map,cell)
    call mesh_index(meshZ,cell,izc_2_c_map,cell)
    call mesh_index(meshX,node,ixc_2_v_map,cell)
    call mesh_index(meshY,node,iyc_2_v_map,cell)
    call mesh_index(meshZ,node,izc_2_v_map,cell)

    ! some constants of mesh
    dtxi = one/meshX%dt ! dtxi=nx
    dtyi = one/meshY%dt ! dtyi=ny
    dtzi = one/meshZ%dt ! dtzi=nz
    if (myid==0) print*, "dtxi, dtyi, dtzi = ", dtxi, dtyi, dtzi
    nx1 = nx+1; nx2 = nx+2
    ny1 = ny+1; ny2 = ny+2
    nz1 = nz+1; nz2 = nz+2
    hx = xmax/nx; hy = ymax/ny; hz = zmax/nz
    hxi = one/hx; hyi = one/hy; hzi = one/hz

    xb = zero; xe = xmax

    yb = meshY%xn(jb+1); ye = meshY%xn(je+2)
    do ipe = 0, nprocs-1
      ybglobal(ipe)=meshY%xn(jbglobal(ipe)+1)
      yeglobal(ipe)=meshY%xn(jeglobal(ipe)+2)
    enddo

    zb = meshZ%xn(kb+1); ze = meshZ%xn(ke+2)
    do ipe = 0, nprocs-1
      zbglobal(ipe)=meshZ%xn(kbglobal(ipe)+1)
      zeglobal(ipe)=meshZ%xn(keglobal(ipe)+2)
    enddo

    ! fraction of local mesh size to global mesh
    volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)

    ! what are they for?
    xb_logical = mesh_unmap(meshX,xb)
    xe_logical = mesh_unmap(meshX,xe)
    yb_logical = mesh_unmap(meshY,yb)
    ye_logical = mesh_unmap(meshY,ye)
    zb_logical = mesh_unmap(meshZ,zb)
    ze_logical = mesh_unmap(meshZ,ze)
    ! notice here 'xb' is different from 'meshX%xb'
    ! the former refers to the beginning of x
    ! the latter refers to 'xbb' which is actually the end of x
    if (myid==0) then
      print*, "xb, meshX%xa, meshX%xb, meshX%ta = ", xb, meshX%xa, meshX%xb, meshX%ta
      print*, "xb_logical, xe_logical = ", xb_logical, xe_logical
      print*, "yb_logical, ye_logical = ", yb_logical, ye_logical
      print*, "zb_logical, ze_logical = ", zb_logical, ze_logical
    endif 

    ! what is qp_cell
    do is = 1, nspec
      npm = npx(is)*npy(is)*npz(is)*nprocs
      dfac(is)=real(ny*nz*nx)/real(npm)
      do ixe = 1, nx2 
        do iye = jb-1, je+1
          do ize = kb-1, ke+1
            qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe) * meshY%dxc(iye+1) * meshZ%dxc(ize+1) * dfac(is)*frac(is)
          enddo
        enddo
      enddo
    enddo

  end subroutine init_mesh

end module m_mesh