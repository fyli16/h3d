module MESH2D
  use MESH_CLASS
  type(MESH) :: meshX, meshY, meshZ
contains
!=======================================================================!
      subroutine MESH_INTERPOLATED_3D(nonuniform_mesh,uniform_mesh,nonuniform_mesh_global)
      use parameter_mod
!      use MESH2D
      implicit none

      double precision:: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,w1,w2,w3,w4,w5,w6,w7,w8
      integer*8:: ix,iy,iz,ixp1,iyp1,izp1,i,j,k,jmin,jmax,kmin,kmax
      double precision,dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(in) :: nonuniform_mesh
      double precision,dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(out):: uniform_mesh
      double precision,dimension(nxmax,0:ny+1,0:nz+1), intent(out):: nonuniform_mesh_global
      double precision,dimension(nxmax,0:ny+1,0:nz+1):: nonuniform_mesh_local
      double precision:: xc_uniform_pos,yc_uniform_pos,zc_uniform_pos

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

      uniform_mesh          = 0.
      nonuniform_mesh_local = 0.

      if (jb == 1) then
        jmin = 0
      else 
        jmin = jb
      endif
      if (je == ny) then
        jmax = ny+1
      else 
        jmax = je
      endif
      if (kb == 1) then
        kmin = 0
      else 
        kmin = kb
      endif
      if (ke == nz) then
        kmax = nz+1
      else 
        kmax = ke
      endif

      do i=1,nxmax
        do j=jmin,jmax
          do k=kmin,kmax
            nonuniform_mesh_local(i,j,k)=nonuniform_mesh(i,j,k)
          enddo
        enddo
      enddo
      call MPI_ALLREDUCE( nonuniform_mesh_local,nonuniform_mesh_global,size(nonuniform_mesh_local)        &
                         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

      do i=2,nx+1
        xc_uniform_pos = (i-1.5)*hx
        rx   = dtxi*MESH_UNMAP(meshX,xc_uniform_pos)+1.50000000000d+00
        ix   = rx
        fx   = rx-ix
        ixp1 = ix+1
        do j=jb,je
          yc_uniform_pos = (j-0.5)*hy
          ry   = dtyi*MESH_UNMAP(meshY,yc_uniform_pos)+1.50000000000d+00
          iy   = ry
          fy   = ry-iy
          iy   = iy-1             ! integer index in y direction starts at 0
          iyp1 = iy+1
          do k=kb,ke
            zc_uniform_pos = (k-0.5)*hz
            rz   = dtzi*MESH_UNMAP(meshZ,zc_uniform_pos)+1.50000000000d+00
            iz   = rz
            fz   = rz-iz
            iz   = iz-1             ! integer index in z direction starts at 0
            izp1 = iz+1
       
            w1=(1.-fx)*(1.-fy)*(1.-fz)
            w2=fx     *(1.-fy)*(1.-fz)
            w3=(1.-fx)*fy     *(1.-fz)
            w4=fx     *fy     *(1.-fz)
            w5=(1.-fx)*(1.-fy)*fz
            w6=fx     *(1.-fy)*fz
            w7=(1.-fx)*fy     *fz
            w8=fx     *fy     *fz
 
            uniform_mesh(i,j,k) =  w1 * nonuniform_mesh_global(ix  ,iy  ,iz  )     &
                                 + w2 * nonuniform_mesh_global(ixp1,iy  ,iz  )     &
                                 + w3 * nonuniform_mesh_global(ix  ,iyp1,iz  )     &
                                 + w4 * nonuniform_mesh_global(ixp1,iyp1,iz  )     &
                                 + w5 * nonuniform_mesh_global(ix  ,iy  ,izp1)     &
                                 + w6 * nonuniform_mesh_global(ixp1,iy  ,izp1)     &
                                 + w7 * nonuniform_mesh_global(ix  ,iyp1,izp1)     &
                                 + w8 * nonuniform_mesh_global(ixp1,iyp1,izp1)
          enddo
        enddo
      enddo

      return
    end subroutine MESH_INTERPOLATED_3D
end module MESH2D
