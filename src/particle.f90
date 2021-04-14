!---------------------------------------------------------------------
subroutine parmov   ! particle move?    
    use parameter_mod
    use mesh_mod
    implicit none

    real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8, &
              by1,by2,by3,by4,by5,by6,by7,by8, &
              bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8, &
              bxa,bya,bza
    real*8 :: ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8, &
              ey1,ey2,ey3,ey4,ey5,ey6,ey7,ey8, &
              ez1,ez2,ez3,ez4,ez5,ez6,ez7,ez8
    real*8 :: deltime1,deltime2,ff,h,hh
    real*8 :: fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8,foxa
    real*8 :: foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8,foya
    real*8 :: foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8,foza
    real*8 :: w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
    real*8 :: vex,vey,vez,vmag,vx_tmp,vy_tmp,vz_tmp
    real*8 :: p2xs,p2ys,p2zs,q_p,th
    real*8 :: wmult

    integer*8 i,ii,iix,iixe,iiy,iiye,iiz,iize,irepeat,irepeatp,is,itmp
    integer*8 iv,iye_cc,ize_cc,j,jv,k,npleavingp,nprecv,nprecvtmp
    integer*8 Storage_Error_p,Storage_Error
    integer*8:: nsendactual, nsendactualp, nrecvactual, nrecvactualp, &
                jj,kk,ix,iy,iz,ixe,iye,ize, &
                ixep1,iyep1,izep1,ixp1,iyp1,izp1
    real*8 :: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,xpart,ypart,zpart
    real*8 :: rxe, rye, rze, fxe, fye, fze, dtxi, dtyi, dtzi
    real*8 :: x_disp,y_disp,z_disp 
    real*8 :: dth
    real*8 :: myranf,fluxran,vxa,vyz,vza
    integer*8:: L, EXIT_CODE_P, EXIT_CODE
    integer*8:: n_fast_removed, n_fast_removed_local,Field_Diverge,Field_Diverge_p
    real*8 :: tx,ty,tz,v_x,v_y,v_z  
    integer*4 :: nescapearr(8),nescapearr_global(8)
    integer*4 :: ppacket(3), ppacketg(3), dpacket(4), dpacketg(4)
    integer*8 :: epacket(2),epacketg(2),loop
    real*8, dimension(3,nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bxyz_av
    real*8 :: TEX1,TEX2,TEX3,TEX4,TEX5,TEX6,TEX7,TEX8  
    real*8 :: TEY1,TEY2,TEY3,TEY4,TEY5,TEY6,TEY7,TEY8  
    real*8 :: TEZ1,TEZ2,TEZ3,TEZ4,TEZ5,TEZ6,TEZ7,TEZ8  
    real*8 :: mX_xa,mX_ta,mX_ca1,mX_ca2,mX_xb,mX_dtdx,mX_tb,mX_cb1,mX_cb2
    real*8 :: mY_xa,mY_ta,mY_ca1,mY_ca2,mY_xb,mY_dtdx,mY_tb,mY_cb1,mY_cb2
    real*8 :: mZ_xa,mZ_ta,mZ_ca1,mZ_ca2,mZ_xb,mZ_dtdx,mZ_tb,mZ_cb1,mZ_cb2

    integer,dimension(8) :: nsend_to_nbr, nbrs
    integer :: idest, max_nsend, max_nrecv
    real*8,dimension(:,:,:),allocatable,target :: packed_pdata_send
    real*8,dimension(:,:),allocatable,target :: packed_pdata_recv
    real*8,pointer :: pp(:,:)
    integer exchange_send_request(8)
    ! real*8 :: mp_elapsed

    data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
    data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
    data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/

    ! timing 'parmov'
    call date_and_time(values=time_begin_array(:,19))

    dtxi = one/meshX%dt
    dtyi = one/meshY%dt
    dtzi = one/meshZ%dt
    dth = dt/2
 
    ! obtain spatial average of bx, by, bz
    bx_av=0.; by_av=0.; bz_av=0.
    do k = kb-1, ke
      do j = jb-1, je
        do i = 1, nx1
          bx_av(i,j,k)=0.125*( bx(i  ,j  ,k  )             &
                              +bx(i+1,j  ,k  )             &
                              +bx(i  ,j+1,k  )             &
                              +bx(i+1,j+1,k  )             &
                              +bx(i  ,j  ,k+1)             &
                              +bx(i+1,j  ,k+1)             &
                              +bx(i  ,j+1,k+1)             &
                              +bx(i+1,j+1,k+1)             &
                              )
          by_av(i,j,k)=0.125*( by(i  ,j  ,k  )             &
                              +by(i+1,j  ,k  )             &
                              +by(i  ,j+1,k  )             &
                              +by(i+1,j+1,k  )             &
                              +by(i  ,j  ,k+1)             &
                              +by(i+1,j  ,k+1)             &
                              +by(i  ,j+1,k+1)             &
                              +by(i+1,j+1,k+1)             &
                              )
          bz_av(i,j,k)=0.125*( bz(i  ,j  ,k  )             &
                              +bz(i+1,j  ,k  )             &
                              +bz(i  ,j+1,k  )             &
                              +bz(i+1,j+1,k  )             &
                              +bz(i  ,j  ,k+1)             &
                              +bz(i+1,j  ,k+1)             &
                              +bz(i  ,j+1,k+1)             &
                              +bz(i+1,j+1,k+1)             &
                              )
        enddo
      enddo
    enddo
    
    call xrealbcc_pack_b(bx_av,by_av,bz_av,1_8,nx,ny,nz)

    ! if ((myid==0) .and. mod(it,n_print)==0) then
    !   print*, " "
    !   print*, "  Calling parmov ..."
    ! endif 

    ! initalize diagnostic variables that keep track of particle number, injection, and escape
    deltime1 = 0.0
    deltime2 = 0.0
    npleavingp = 0
 
    if (dt .ne. 0) then ! if dt==0, no actual particle push is done
      ! advance particles for a first half step
      call date_and_time(values=time_begin_array(:,13))
      call push
      call date_and_time(values=time_end_array(:,13))
      call accumulate_time(time_begin_array(1,13),time_end_array(1,13),time_elapsed(13))

      ! check particles
      call date_and_time(values=time_begin_array(:,14))
      call particle_boundary
      call date_and_time(values=time_end_array(:,14))
      call accumulate_time(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))

      ! collect Vi, ni at half step
      call date_and_time(values=time_begin_array(:,15))
      do is = 1, nspec
        nptotp = 0
        npart(is) = 0
        do iiz = kb-1, ke
          do iiy = jb-1, je
            do iix = 1, nx1
              np=iphead(iix,iiy,iiz,is)
              do while (np.ne.0)
                nptotp = nptotp+1           !count particles
                npart(is) = npart(is) + 1 !count particles in each species
                L=NP

                q_p = qp(l)

                ! Nonuniform mesh - using mesh_unmap
                rx=dtxi*mesh_unmap(meshX,x(l))+1.50000000000d+00
                ry=dtyi*mesh_unmap(meshY,y(l))+1.50000000000d+00
                rz=dtzi*mesh_unmap(meshZ,z(l))+1.50000000000d+00
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

                w1=q_p*(1.-fx)*(1.-fy)*(1.-fz)
                w2=q_p*fx     *(1.-fy)*(1.-fz)
                w3=q_p*(1.-fx)*fy     *(1.-fz)
                w4=q_p*fx     *fy     *(1.-fz)
                w5=q_p*(1.-fx)*(1.-fy)*fz
                w6=q_p*fx     *(1.-fy)*fz
                w7=q_p*(1.-fx)*fy     *fz
                w8=q_p*fx     *fy     *fz
   
                dnsh(ix  ,iy  ,iz  ,is)=dnsh(ix  ,iy  ,iz  ,is)+w1
                dnsh(ixp1,iy  ,iz  ,is)=dnsh(ixp1,iy  ,iz  ,is)+w2
                dnsh(ix  ,iyp1,iz  ,is)=dnsh(ix  ,iyp1,iz  ,is)+w3
                dnsh(ixp1,iyp1,iz  ,is)=dnsh(ixp1,iyp1,iz  ,is)+w4
                dnsh(ix  ,iy  ,izp1,is)=dnsh(ix  ,iy  ,izp1,is)+w5
                dnsh(ixp1,iy  ,izp1,is)=dnsh(ixp1,iy  ,izp1,is)+w6
                dnsh(ix  ,iyp1,izp1,is)=dnsh(ix  ,iyp1,izp1,is)+w7
                dnsh(ixp1,iyp1,izp1,is)=dnsh(ixp1,iyp1,izp1,is)+w8
  
                vxs(ix  ,iy  ,iz  ,is)=vxs(ix  ,iy  ,iz  ,is)+w1*vx(l)
                vxs(ixp1,iy  ,iz  ,is)=vxs(ixp1,iy  ,iz  ,is)+w2*vx(l)
                vxs(ix  ,iyp1,iz  ,is)=vxs(ix  ,iyp1,iz  ,is)+w3*vx(l)
                vxs(ixp1,iyp1,iz  ,is)=vxs(ixp1,iyp1,iz  ,is)+w4*vx(l)
                vxs(ix  ,iy  ,izp1,is)=vxs(ix  ,iy  ,izp1,is)+w5*vx(l)
                vxs(ixp1,iy  ,izp1,is)=vxs(ixp1,iy  ,izp1,is)+w6*vx(l)
                vxs(ix  ,iyp1,izp1,is)=vxs(ix  ,iyp1,izp1,is)+w7*vx(l)
                vxs(ixp1,iyp1,izp1,is)=vxs(ixp1,iyp1,izp1,is)+w8*vx(l)

                vys(ix  ,iy  ,iz  ,is)=vys(ix  ,iy  ,iz  ,is)+w1*vy(l)
                vys(ixp1,iy  ,iz  ,is)=vys(ixp1,iy  ,iz  ,is)+w2*vy(l)
                vys(ix  ,iyp1,iz  ,is)=vys(ix  ,iyp1,iz  ,is)+w3*vy(l)
                vys(ixp1,iyp1,iz  ,is)=vys(ixp1,iyp1,iz  ,is)+w4*vy(l)
                vys(ix  ,iy  ,izp1,is)=vys(ix  ,iy  ,izp1,is)+w5*vy(l)
                vys(ixp1,iy  ,izp1,is)=vys(ixp1,iy  ,izp1,is)+w6*vy(l)
                vys(ix  ,iyp1,izp1,is)=vys(ix  ,iyp1,izp1,is)+w7*vy(l)
                vys(ixp1,iyp1,izp1,is)=vys(ixp1,iyp1,izp1,is)+w8*vy(l)

                vzs(ix  ,iy  ,iz  ,is)=vzs(ix  ,iy  ,iz  ,is)+w1*vz(l)
                vzs(ixp1,iy  ,iz  ,is)=vzs(ixp1,iy  ,iz  ,is)+w2*vz(l)
                vzs(ix  ,iyp1,iz  ,is)=vzs(ix  ,iyp1,iz  ,is)+w3*vz(l)
                vzs(ixp1,iyp1,iz  ,is)=vzs(ixp1,iyp1,iz  ,is)+w4*vz(l)
                vzs(ix  ,iy  ,izp1,is)=vzs(ix  ,iy  ,izp1,is)+w5*vz(l)
                vzs(ixp1,iy  ,izp1,is)=vzs(ixp1,iy  ,izp1,is)+w6*vz(l)
                vzs(ix  ,iyp1,izp1,is)=vzs(ix  ,iyp1,izp1,is)+w7*vz(l)
                vzs(ixp1,iyp1,izp1,is)=vzs(ixp1,iyp1,izp1,is)+w8*vz(l)
   
                np=link(np)
              enddo ! while
            enddo !for iix
          enddo !for iiy
        enddo !for iiz

        call xreal(DNSH(1,jb-1,kb-1,is),nx,ny,nz)
        call xreal(VXS(1,jb-1,kb-1,is),nx,ny,nz)
        call xreal(VYS(1,jb-1,kb-1,is),nx,ny,nz)
        call xreal(VZS(1,jb-1,kb-1,is),nx,ny,nz)
        call xrealbcc(DNSH(1,jb-1,kb-1,is),1_8,nx,ny,nz)
        call xrealbcc(VXS(1,jb-1,kb-1,is),1_8,nx,ny,nz)
        call xrealbcc(VYS(1,jb-1,kb-1,is),1_8,nx,ny,nz)
        call xrealbcc(VZS(1,jb-1,kb-1,is),1_8,nx,ny,nz)

        do iiz = kb-1, ke+1
          do iiy = jb-1, je+1
            do iix = 1, nx2
              dnsh(iix,iiy,iiz,is)= dnsh(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              vxs(iix,iiy,iiz,is) = vxs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              vys(iix,iiy,iiz,is) = vys(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              vzs(iix,iiy,iiz,is) = vzs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            enddo
          enddo
        enddo

      enddo ! for is
      call date_and_time(values=time_end_array(:,15))
      call accumulate_time(time_begin_array(1,15),time_end_array(1,15),time_elapsed(15))

      ! advance particles for a seconnd half step
      call date_and_time(values=time_begin_array(:,13))
      do is = 1, nspec
        do iize = kb-1, ke
          do iiye = jb-1, je
            do iixe = 1, nx1
              np = iphead(iixe,iiye,iize,is)
              do while (np.ne.0)
                L=NP
                x_disp = dth*vx(l)
                y_disp = dth*vy(l)
                z_disp = dth*vz(l)

                x(l)=x(l)+ x_disp
                y(l)=y(l)+ y_disp
                z(l)=z(l)+ z_disp

                NP=LINK(NP)
              enddo ! while
            enddo ! for iixe
          enddo ! for iiye
        enddo ! for iize
      enddo ! for is
      call date_and_time(values=time_end_array(:,13))
      call accumulate_time(time_begin_array(1,13),time_end_array(1,13),time_elapsed(13))

    endif ! dt>0

    ! check particles
    call date_and_time(values=time_begin_array(:,14))
    call particle_boundary
    call date_and_time(values=time_end_array(:,14))
    call accumulate_time(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))

    ! collect density
    call date_and_time(values=time_begin_array(:,15))
    do is = 1, nspec
      if (testorbt) goto 10
      nptotp = 0
      npart(is) = 0
      do IIZ=KB-1,KE
        do IIY=JB-1,JE
          do IIX=1,NX1
            np=iphead(iix,iiy,iiz,is)
            do while (np.ne.0)
              nptotp=nptotp+1           !count particles
              npart(is) = npart(is) + 1 !count particles in each species
              L = NP

              q_p = qp(l)

              ! Nonuniform mesh - using mesh_unmap
              rx=dtxi*mesh_unmap(meshX,x(l))+1.50000000000d+00
              ry=dtyi*mesh_unmap(meshY,y(l))+1.50000000000d+00
              rz=dtzi*mesh_unmap(meshZ,z(l))+1.50000000000d+00
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
 
              w1=q_p*(1.-fx)*(1.-fy)*(1.-fz)
              w2=q_p*fx     *(1.-fy)*(1.-fz)
              w3=q_p*(1.-fx)*fy     *(1.-fz)
              w4=q_p*fx     *fy     *(1.-fz)
              w5=q_p*(1.-fx)*(1.-fy)*fz
              w6=q_p*fx     *(1.-fy)*fz
              w7=q_p*(1.-fx)*fy     *fz
              w8=q_p*fx     *fy     *fz

              dns(ix  ,iy  ,iz  ,is)=dns(ix  ,iy  ,iz  ,is)+w1
              dns(ixp1,iy  ,iz  ,is)=dns(ixp1,iy  ,iz  ,is)+w2
              dns(ix  ,iyp1,iz  ,is)=dns(ix  ,iyp1,iz  ,is)+w3
              dns(ixp1,iyp1,iz  ,is)=dns(ixp1,iyp1,iz  ,is)+w4
              dns(ix  ,iy  ,izp1,is)=dns(ix  ,iy  ,izp1,is)+w5
              dns(ixp1,iy  ,izp1,is)=dns(ixp1,iy  ,izp1,is)+w6
              dns(ix  ,iyp1,izp1,is)=dns(ix  ,iyp1,izp1,is)+w7
              dns(ixp1,iyp1,izp1,is)=dns(ixp1,iyp1,izp1,is)+w8
 
              np=link(np)
            enddo
          enddo
        enddo
      enddo

 10   continue 

      nescapearr(1) = nescape(is)
      nescapearr(2) = nescape_yz(is)
      nescapearr(3) = nescape_zy(is)
      nescapearr(4) = nescape_xy(is)
      nescapearr(5) = nescape_yx(is)
      nescapearr(6) = nescape_zx(is)
      nescapearr(7) = nescape_xz(is)
      nescapearr(8) = npart(is)

      call MPI_ALLREDUCE(nescapearr,nescapearr_global,8,MPI_INTEGER4,MPI_SUM,COMM2D,IERR)
      nescape_global(is)    = nescapearr_global(1)
      nescape_yz_global(is) = nescapearr_global(2)
      nescape_zy_global(is) = nescapearr_global(3)
      nescape_xy_global(is) = nescapearr_global(4)
      nescape_yx_global(is) = nescapearr_global(5)
      nescape_zx_global(is) = nescapearr_global(6)
      nescape_xz_global(is) = nescapearr_global(7)
      npart_global(is)      = nescapearr_global(8)

      deltime2 = deltime2 + real(clock_time1-clock_time)
 
      call xreal(dns(1,jb-1,kb-1,is),nx,ny,nz)
      call xrealbcc(dns(1,jb-1,kb-1,is),1_8,nx,ny,nz)

      do iiz = kb-1, ke+1
        do iiy = jb-1, je+1
          do iix = 1, nx2
            dns(iix,iiy,iiz,is) = dns(iix,iiy,iiz,is)/(meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          enddo
        enddo
      enddo

    enddo ! for is

    call date_and_time(values=time_end_array(:,15))
    call accumulate_time(time_begin_array(1,15),time_end_array(1,15),time_elapsed(15))

    ! diagnostic info
    epacket(1) = nptotp
    epacket(2) = npleavingp
    call MPI_ALLREDUCE(epacket,epacketg,2,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
    nptot      = epacketg(1)
    npleaving  = epacketg(2)

    ! if (myid == 0.and.mod(it,n_print)==0) then
    !   do is=1,nspec
    !     if (is == 1) then
    !       print*,
    !       print*, "it = ", it
    !       print*, "species #    ninj     nescape     ntot  "
    !     endif
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_global(is), npart_global(is)
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_yz_global(is), npart_global(is)
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_zy_global(is), npart_global(is)
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_xy_global(is), npart_global(is)
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_yx_global(is), npart_global(is)
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_xz_global(is), npart_global(is)
    !     write(6,"(5x,i2,5x,i8,2x,i8,2x,i10)") is, ninj_global(is), nescape_zx_global(is), npart_global(is)
    !   enddo

    !   if (nspec >= 2) then
    !     write(6,"(5x,'sum',4x,i8,2x,i8,2x,i10)") ninj_global(1)+ninj_global(2), &
    !               nescape_global(1)+nescape_global(2), npart_global(1)+npart_global(2)
    !   else
    !     write(6,"(5x,'sum',4x,i8,2x,i8,2x,i10)") ninj_global(1), nescape_global(1), npart_global(1)
    !   endif
    ! endif

    ninj        = 0
    ninj_global = 0

    call date_and_time(values=time_end_array(:,19))
    call accumulate_time(time_begin_array(1,19),time_end_array(1,19),time_elapsed(19))

    return
end subroutine parmov


!---------------------------------------------------------------------
! This subourtine pushes particles for half a step
!---------------------------------------------------------------------
subroutine push
  use parameter_mod
  use mesh_mod
  implicit none

  integer*8 :: is,iixe,iiye,iize,l
  integer*8:: ix,iy,iz,ixe,iye,ize,ixep1,iyep1,izep1,ixp1,iyp1,izp1
  real*8 :: wmult, h, hh, dth
  real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8, &
            by1,by2,by3,by4,by5,by6,by7,by8, &
            bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8, & 
            bxa,bya,bza
  real*8 :: ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8, &
            ey1,ey2,ey3,ey4,ey5,ey6,ey7,ey8, &
            ez1,ez2,ez3,ez4,ez5,ez6,ez7,ez8, &
            exa,eya,eza
  real*8 :: fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8,foxa
  real*8 :: foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8,foya
  real*8 :: foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8,foza
  real*8 :: w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
  real*8 :: vex,vey,vez
  real*8:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
  real*8 :: p2xs,p2ys,p2zs,ff
  real*8:: x_disp,y_disp,z_disp,disp_max_p(3),disp_max(3), &
          y_disp_max_p, x_disp_max_p, z_disp_max_p, &
          y_disp_max, x_disp_max, z_disp_max
  integer*8 :: Courant_Violation,Courant_Violation_p 

  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt

  do is = 1, nspec
    Courant_Violation_p = 0
    x_disp_max_p        = 0
    y_disp_max_p        = 0
    z_disp_max_p        = 0

    wmult = wspec(is)
    h = dt*qspec(is)/wmult
    hh = .5*h
    dth = dt/2
 
    do IIZE = KB-1,KE
      do IIYE = JB-1,JE
        do IIXE = 1, NX1
          NP=IPHEAD(IIXE,IIYE,IIZE,IS)
          do while (NP.ne.0)
            L=NP

            ! Uniform mesh - Same as in version 5.0
            ! rxe=hxi*x(l)+1.5000000000000001d+00
            ! rye=hyi*y(l)+0.5000000000000001d+00
            ! rze=hzi*z(l)+0.5000000000000001d+00
            ! ixe=rxe
            ! iye=rye
            ! ize=rze
            ! ixe=max(1   ,min(ixe,nx1))
            ! iye=max(jb-1,min(iye,je))
            ! ize=max(kb-1,min(ize,ke))
            ! ixep1 = ixe+1
            ! iyep1 = iye+1
            ! izep1 = ize+1
            ! fxe=rxe-ixe
            ! fye=rye-iye
            ! fze=rze-ize

            ! Nonuniform mesh - without using mesh_unmap
            ! rxe=hxi*x(l)+1.500000000000000d+00
            ! rye=hyi*y(l)+1.500000000000000d+00
            ! rze=hzi*z(l)+1.500000000000000d+00
            ! ixe=rxe
            ! iye=rye
            ! ize=rze
            ! ixe=ixc_2_c_map(ixe)
            ! iye=iyc_2_c_map(iye)
            ! ize=izc_2_c_map(ize)
            ! fxe=(x(l)-meshX%xc(ixe))/meshX%dxn(ixep1)
            ! fye=(y(l)-meshY%xc(iye))/meshY%dxn(iyep1)
            ! fze=(z(l)-meshZ%xc(ize))/meshZ%dxn(izep1)

            ! Nonuniform mesh - using mesh_unmap
            rxe=dtxi*mesh_unmap(meshX,x(l))+1.50000000000d+00
            rye=dtyi*mesh_unmap(meshY,y(l))+1.50000000000d+00
            rze=dtzi*mesh_unmap(meshZ,z(l))+1.50000000000d+00
            ixe=rxe
            iye=rye
            ize=rze
            fxe=rxe-ixe
            fye=rye-iye
            fze=rze-ize
            iye=iye-1             ! integer index in y direction starts at 0
            ize=ize-1             ! integer index in z direction starts at 0

            ixep1 = ixe+1
            iyep1 = iye+1
            izep1 = ize+1

            ! Diagnosic test on particle cell index algorithm
            ! if (     fxe < 0. .or. fxe > 1.                &
            !     .or. fye < 0. .or. fye > 1.                &
            !     .or. fze < 0. .or. fze > 1.) then 
            !     print*, " SCATTER LOOP"
            !     print*, fxe,fye,fze
            !     print*, " x;",x(l),meshX%xn(ixe),meshX%xc(ixe)
            !     print*, " y;",y(l),meshY%xn(iye),meshY%xc(iye)
            !     print*, " z;",z(l),meshZ%xn(ize),meshZ%xc(ize)
            !     print*, " r; ",rxe,rye,rze
            !     print*, " ixc_map; ",ixe,iye,ize
            !     print*, " xc     ; ",meshX%xc(ixe),meshY%xc(iye),meshZ%xc(ize)
            ! endif

            w1e=(1.-fxe)*(1.-fye)*(1.-fze)
            w2e=fxe*(1.-fye)*(1.-fze)
            w3e=(1.-fxe)*fye*(1.-fze)
            w4e=fxe*fye*(1.-fze)
            w5e=(1.-fxe)*(1.-fye)*fze
            w6e=fxe*(1.-fye)*fze
            w7e=(1.-fxe)*fye*fze
            w8e=fxe*fye*fze

            ex1=ex(ixe  ,iye  ,ize  )
            ex2=ex(ixep1,iye  ,ize  )
            ex3=ex(ixe  ,iyep1,ize  )
            ex4=ex(ixep1,iyep1,ize  )
            ex5=ex(ixe  ,iye  ,izep1)
            ex6=ex(ixep1,iye  ,izep1)
            ex7=ex(ixe  ,iyep1,izep1)
            ex8=ex(ixep1,iyep1,izep1)
            ey1=ey(ixe  ,iye  ,ize  )
            ey2=ey(ixep1,iye  ,ize  )
            ey3=ey(ixe  ,iyep1,ize  )
            ey4=ey(ixep1,iyep1,ize  )
            ey5=ey(ixe  ,iye  ,izep1)
            ey6=ey(ixep1,iye  ,izep1)
            ey7=ey(ixe  ,iyep1,izep1)
            ey8=ey(ixep1,iyep1,izep1)
            ez1=ez(ixe  ,iye  ,ize  )
            ez2=ez(ixep1,iye  ,ize  )
            ez3=ez(ixe  ,iyep1,ize  )
            ez4=ez(ixep1,iyep1,ize  )
            ez5=ez(ixe  ,iye  ,izep1)
            ez6=ez(ixep1,iye  ,izep1)
            ez7=ez(ixe  ,iyep1,izep1)
            ez8=ez(ixep1,iyep1,izep1)

            bx1=bx_av(ixe  ,iye  ,ize  )
            bx2=bx_av(ixep1,iye  ,ize  )
            bx3=bx_av(ixe  ,iyep1,ize  )
            bx4=bx_av(ixep1,iyep1,ize  )
            bx5=bx_av(ixe  ,iye  ,izep1)
            bx6=bx_av(ixep1,iye  ,izep1)
            bx7=bx_av(ixe  ,iyep1,izep1)
            bx8=bx_av(ixep1,iyep1,izep1)
            by1=by_av(ixe  ,iye  ,ize  )
            by2=by_av(ixep1,iye  ,ize  )
            by3=by_av(ixe  ,iyep1,ize  )
            by4=by_av(ixep1,iyep1,ize  )
            by5=by_av(ixe  ,iye  ,izep1)
            by6=by_av(ixep1,iye  ,izep1)
            by7=by_av(ixe  ,iyep1,izep1)
            by8=by_av(ixep1,iyep1,izep1)
            bz1=bz_av(ixe  ,iye  ,ize  )
            bz2=bz_av(ixep1,iye  ,ize  )
            bz3=bz_av(ixe  ,iyep1,ize  )
            bz4=bz_av(ixep1,iyep1,ize  )
            bz5=bz_av(ixe  ,iye  ,izep1)
            bz6=bz_av(ixep1,iye  ,izep1)
            bz7=bz_av(ixe  ,iyep1,izep1)
            bz8=bz_av(ixep1,iyep1,izep1)

            fox1=fox(ixe  ,iye  ,ize  )
            fox2=fox(ixep1,iye  ,ize  )
            fox3=fox(ixe  ,iyep1,ize  )
            fox4=fox(ixep1,iyep1,ize  )
            fox5=fox(ixe  ,iye  ,izep1)
            fox6=fox(ixep1,iye  ,izep1)
            fox7=fox(ixe  ,iyep1,izep1)
            fox8=fox(ixep1,iyep1,izep1)
            foy1=foy(ixe  ,iye  ,ize  )
            foy2=foy(ixep1,iye  ,ize  )
            foy3=foy(ixe  ,iyep1,ize  )
            foy4=foy(ixep1,iyep1,ize  )
            foy5=foy(ixe  ,iye  ,izep1)
            foy6=foy(ixep1,iye  ,izep1)
            foy7=foy(ixe  ,iyep1,izep1)
            foy8=foy(ixep1,iyep1,izep1)
            foz1=foz(ixe  ,iye  ,ize  )
            foz2=foz(ixep1,iye  ,ize  )
            foz3=foz(ixe  ,iyep1,ize  )
            foz4=foz(ixep1,iyep1,ize  )
            foz5=foz(ixe  ,iye  ,izep1)
            foz6=foz(ixep1,iye  ,izep1)
            foz7=foz(ixe  ,iyep1,izep1)
            foz8=foz(ixep1,iyep1,izep1)

            exa=w1e*ex1+w2e*ex2+w3e*ex3+w4e*ex4      &
                +w5e*ex5+w6e*ex6+w7e*ex7+w8e*ex8      &
                +w1e*fox1+w2e*fox2+w3e*fox3+w4e*fox4  &
                +w5e*fox5+w6e*fox6+w7e*fox7+w8e*fox8
            eya=w1e*ey1+w2e*ey2+w3e*ey3+w4e*ey4      &
                +w5e*ey5+w6e*ey6+w7e*ey7+w8e*ey8      &
                +w1e*foy1+w2e*foy2+w3e*foy3+w4e*foy4  &
                +w5e*foy5+w6e*foy6+w7e*foy7+w8e*foy8
            eza=w1e*ez1+w2e*ez2+w3e*ez3+w4e*ez4      &
                +w5e*ez5+w6e*ez6+w7e*ez7+w8e*ez8      &
                +w1e*foz1+w2e*foz2+w3e*foz3+w4e*foz4  &
                +w5e*foz5+w6e*foz6+w7e*foz7+w8e*foz8

            bxa=w1e*bx1+w2e*bx2+w3e*bx3+w4e*bx4      &
                +w5e*bx5+w6e*bx6+w7e*bx7+w8e*bx8
            bya=w1e*by1+w2e*by2+w3e*by3+w4e*by4      &
                +w5e*by5+w6e*by6+w7e*by7+w8e*by8
            bza=w1e*bz1+w2e*bz2+w3e*bz3+w4e*bz4      &
                +w5e*bz5+w6e*bz6+w7e*bz7+w8e*bz8

            ff=2./(1.+hh*hh*(bxa**2+bya**2+bza**2))
            vex=vx(l)+exa*hh
            vey=vy(l)+eya*hh
            vez=vz(l)+eza*hh
            p2xs=vex+(vey*bza-vez*bya)*hh
            p2ys=vey+(vez*bxa-vex*bza)*hh
            p2zs=vez+(vex*bya-vey*bxa)*hh
            vx(l)=vex+ff*(p2ys*bza-p2zs*bya)*hh+exa*hh
            vy(l)=vey+ff*(p2zs*bxa-p2xs*bza)*hh+eya*hh
            vz(l)=vez+ff*(p2xs*bya-p2ys*bxa)*hh+eza*hh

            ! advance particles for a half step to calcualte Vi
            x_disp = dth*vx(l)
            y_disp = dth*vy(l)
            z_disp = dth*vz(l)

            x(l)=x(l)+ x_disp
            y(l)=y(l)+ y_disp
            z(l)=z(l)+ z_disp

            if ( abs(2*x_disp/meshX%dxn(ixep1)) > 1.0 .or.                                               &
                  abs(2*y_disp/meshY%dxn(iyep1)) > 1.0 .or.                                               &
                  abs(2*z_disp/meshZ%dxn(izep1)) > 1.0) Courant_Violation_p = Courant_Violation_p + 1

            x_disp_max_p = max(x_disp_max_p,abs(x_disp)/meshX%dxn(ixep1))
            y_disp_max_p = max(y_disp_max_p,abs(y_disp)/meshY%dxn(iyep1))
            z_disp_max_p = max(z_disp_max_p,abs(z_disp)/meshZ%dxn(izep1))

            ! particle tracking
            if (ptag(np).ne.0) then
              ntot=ntot+1
              buf_p1(1,ntot)=x(np)
              buf_p1(2,ntot)=y(np)
              buf_p1(3,ntot)=z(np)
              buf_p1(4,ntot)=vx(np)
              buf_p1(5,ntot)=vy(np)
              buf_p1(6,ntot)=vz(np)
              buf_p1(7,ntot)=qp(np)
              buf_p1(8,ntot)=DBLE(ptag(np))
              buf_p1(9,ntot)=exa
              buf_p1(10,ntot)=eya
              buf_p1(11,ntot)=eza
              buf_p1(12,ntot)=bxa
              buf_p1(13,ntot)=bya
              buf_p1(14,ntot)=bza
            endif

            NP=LINK(NP)

          enddo ! while
        enddo
      enddo
    enddo

    disp_max_p(1) = x_disp_max_p
    disp_max_p(2) = y_disp_max_p
    disp_max_p(3) = z_disp_max_p
    call MPI_ALLREDUCE(disp_max_p,disp_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)   !LAURA
    x_disp_max = disp_max(1)
    y_disp_max = disp_max(2)
    z_disp_max = disp_max(3)

    ! if ((myid == 0).and.mod(it,n_print)==0) then
    !   print*, " maximum x-displacement/dx = ",x_disp_max
    !   print*, " maximum y-displacement/dy = ",y_disp_max
    !   print*, " maximum z-displacement/dz = ",z_disp_max
    ! endif

    call MPI_ALLREDUCE(Courant_Violation_p,Courant_Violation,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

    if (Courant_Violation /= 0) then
        if (myid == 0) print*, "Particle displacements exceed cell size ",Courant_Violation," times"
        call MPI_FINALIZE(IERR)
        STOP
    endif
  enddo !is

end subroutine push


!---------------------------------------------------------------------
subroutine particle_boundary
    use parameter_mod
    use mesh_mod
    implicit none
    integer*4 :: ppacket(3),ppacketg(3),dpacket(4),dpacketg(4)
    integer*8 :: epacket(2),epacketg(2),loop
    real*8 :: xpart,ypart,zpart
    real*8 :: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
    integer*8 :: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
    integer*8 :: Storage_Error_p,Storage_Error
    integer*8 :: l
    integer, dimension(8) :: nsend_to_nbr, nbrs
    integer :: idest,max_nsend,max_nrecv
    real*8, dimension(:,:,:),allocatable,target :: packed_pdata_send
    real*8, dimension(:,:),allocatable,target :: packed_pdata_recv
    real*8, pointer :: pp(:,:)
    integer :: exchange_send_request(8)
    integer*8 :: iv,iye_cc,ize_cc,j,jv,k,npleavingp,nprecv,nprecvtmp
    real*8 :: v_limit
    integer*8 :: n_fast_removed,n_fast_removed_local,Field_Diverge,Field_Diverge_p
    integer*8 :: i,ii,iix,iixe,iiy,iiye,iiz,iize,irepeat,irepeatp,is,itmp
    real*8 :: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min
    Storage_Error_p = 0
    Field_Diverge_p = 0
    dtxi = 1./meshX%dt
    dtyi = 1./meshY%dt
    dtzi = 1./meshZ%dt

    ! determine the velocity limit based on CFL condition

    ! Uniform mesh - Same as in version 5.0
    ! if (dt /=0.) then
    !   if (ndim == 1) then
    !   v_limit = min(hx/dt,hy/dt)
    !   else
    !   v_limit = min(hx/dt,hy/dt,hz/dt)
    ! endif
    ! else
    !   v_limit=1.d+10
    ! endif


    ! Nonuniform mesh
    hxmin=meshX%dxc(1)
    hxmax=hxmin
    do i=1,size(meshX%dxc)
      if(meshX%dxc(i) < hxmin) hxmin=meshX%dxc(i)
      if(meshX%dxc(i) > hxmax) hxmax=meshX%dxc(i)
    enddo
    hymin=meshY%dxc(1)
    hymax=hymin
    do i=1,size(meshY%dxc)
      if(meshY%dxc(i) < hymin) hymin=meshY%dxc(i)
      if(meshY%dxc(i) > hymax) hymax=meshY%dxc(i)
    enddo
    hzmin=meshZ%dxc(1)
    hzmax=hzmin
    do i=1,size(meshZ%dxc)
      if(meshZ%dxc(i) < hzmin) hzmin=meshZ%dxc(i)
      if(meshZ%dxc(i) > hzmax) hzmax=meshZ%dxc(i)
    enddo
    cell_size_min = min(hxmin,hymin,hzmin)
    v_limit=(cell_size_min/dtwci)/wpiwci

    do IS=1, nspec    
      irepeatp=0                !  are no fast particles, it will be called 
      nsendp=0
      nrecvp=0
      ipleft (is)=0
      iprite (is)=0
      ipsendleft(is)=0
      ipsendrite(is)=0

      ipsend(is)=0
      ipsendtop(is)=0
      ipsendbot(is)=0
      ipsendlefttop(is)=0
      ipsendleftbot(is)=0
      ipsendritetop(is)=0
      ipsendritebot(is)=0
 
      nescape(is)=0
      nescape_yz(is)=0
      nescape_zy(is)=0
      nescape_xy(is)=0
      nescape_yx(is)=0
      nescape_zx(is)=0
      nescape_xz(is)=0

      n_fast_removed_local = 0 

      do IZE = KB-1,KE
        do IYE = JB-1,JE
          do IXE = 1,NX1
            iptemp(ixe,iye,ize,is)=0

            ! mark particles that need to be sent to other processors
            NP=IPHEAD(IXE,IYE,IZE,IS)  !VR: the first particle in the cell
            do while (NP.ne.0)         !VR: loop over particles in the cell
              xpart=x(np)
              ypart=y(np)
              zpart=z(np)
              l=np

              ! Remove fast particles
              if (abs(vx(l)) >= v_limit .or. abs(vy(l)) >= v_limit .or. abs(vz(l)) >= v_limit) then
                n_fast_removed_local = n_fast_removed_local + 1
                iphead(ixe,iye,ize,is)=link(np)    !VR iphead is now the next particle in the list. Note that iphead
                                                    !VR will be re-assigned after the loop over particles in the cell
                link(np)=ipstore                   !VR current particle points to the old start of empty list
                ipstore=np                         !VR start of empty list is this particle
                np=iphead(ixe,iye,ize,is)          !VR go to the next particle
                goto 15
              endif

              ! VR: removed interaction with a sphere
              ! VR: removed re-injection on y and z global domain boundaries

              ! VR: loop particle through periodic X boundary --------
              ! VR: if the particle is oitside y & z global domain boundaries,
              ! VR: it will be sent to the corresponding processor
              if (x(l) < zero) then
                  x(l) = x(l) + xmax
              endif

              if (x(l) > xmax) then
                  x(l) = x(l) - xmax
              endif

              xpart = x(l)
              ypart = y(l)
              zpart = z(l)
                
              if ((zpart <= ze.and.zpart >= zb).and.(ypart <= ye.and.ypart >= yb)) then
                  ! VR: particle inside local domain (in y and z)
                  ! VR: here we simply re-order the list:
                iphead(ixe,iye,ize,is)=link(np)      !VR: next particle to consider (tmp)
                link(np)=iptemp(ixe,iye,ize,is)      !VR: the next particle in the new list is the prev. part. in the old list
                iptemp(ixe,iye,ize,is)=np            !VR: imptemp becomes this particle (head of the new list)
                np=iphead(ixe,iye,ize,is)            !VR: next particle to consider is from link
                                                      !VR: at the end of the list, the following happens:
                                                      !VR: iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
                                                      !VR: iptemp(ixe,iye,ize,is)=0
              else
                ! VR: particle outside local domain
                if (ypart <= ye.and.ypart >= yb) then
                    ! VR: particle is inside domain in y (but outside in z?)
                  iye_cc=jb    !VR: iye_cc is the index into idmap array. here we simply get the id at the left boundary
                else
                  if (ypart > ye) then !VR: particle is leaving through the right boundary
                    iye_cc=je+1 
                  else                 !VR: particle is leaving through the left boundary
                    iye_cc=jb-1 
                  endif
                endif

                ! VR: repeat the same exercise in z
                if (zpart <= ze.and.zpart >= zb) then
                  ize_cc=kb
                else
                  if (zpart > ze) then
                    ize_cc=ke+1 
                  else
                    ize_cc=kb-1 
                  endif
                endif
                    
                !VR: increase the count of particles that will leave the domain
                nsendp(idmap_yz(iye_cc,ize_cc))=nsendp(idmap_yz(iye_cc,ize_cc))+1
                iphead(ixe,iye,ize,is)=link(np)
                link(np)=ipsend(is)     !VR: add particle to the ipsend list (particles to send)
                ipsend(is)=np           !VR: the new head for ipsend list (this particle)
                np=iphead(ixe,iye,ize,is) !VR: next particle to consider (from the old link list)
              endif ! VR: check particle inside/outside domain

15            CONTINUE
              enddo !VR: loop over particles in the cell
              iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)    !VR: save the new head (re-ordered)
              iptemp(ixe,iye,ize,is)=0                         !VR: reset the temp list
            enddo !VR: grid loop in x
          enddo !VR: grid loop in y
        enddo !VR: grid loop in z

        ! exchange data among processes and compute to see how many
        ! particles each process has to send to, and receive from, other
        ! processes

        ! nsendp(nbrleft) -> how many particles to be sent from myid to 
        !                     left neighbor
        ! nsendp(nbrrite) -> how many particles to be sent from myid to 
        !                     right neighbor
        ! nrecvp(nbrleft) -> how many particles to be sent to myid from 
        !                     left neighbor
        ! nrecvp(nbrrite) -> how many particles to be sent to myid from 
        !                     right neighbor

        ! exchange information about particle numbers in two steps. First, 
        ! send to right and receive from left. 

        ! VR: is this comment still true? I can not see this
        ! Note that the processors on 
        ! the physical boundaries only send (myid == nbrleft) or receive 
        ! (myid == nbrrite)
        call MPI_SENDRECV(nsendp(NBRLEFT   ),1,MPI_INTEGER8,NBRLEFT   ,0,&
                          nrecvp(NBRRITE   ),1,MPI_INTEGER8,NBRRITE   ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRRITE   ),1,MPI_INTEGER8,NBRRITE   ,0,&
                          nrecvp(NBRLEFT   ),1,MPI_INTEGER8,NBRLEFT   ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRTOP    ),1,MPI_INTEGER8,NBRTOP    ,0,&
                          nrecvp(NBRBOT    ),1,MPI_INTEGER8,NBRBOT    ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRBOT    ),1,MPI_INTEGER8,NBRBOT    ,0,&
                          nrecvp(NBRTOP    ),1,MPI_INTEGER8,NBRTOP    ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRLEFTTOP),1,MPI_INTEGER8,NBRLEFTTOP,0,&
                          nrecvp(NBRRITEBOT),1,MPI_INTEGER8,NBRRITEBOT,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRRITEBOT),1,MPI_INTEGER8,NBRRITEBOT,0,&
                          nrecvp(NBRLEFTTOP),1,MPI_INTEGER8,NBRLEFTTOP,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRRITETOP),1,MPI_INTEGER8,NBRRITETOP,0,&
                          nrecvp(NBRLEFTBOT),1,MPI_INTEGER8,NBRLEFTBOT,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRLEFTBOT),1,MPI_INTEGER8,NBRLEFTBOT,0,&
                          nrecvp(NBRRITETOP),1,MPI_INTEGER8,NBRRITETOP,0,&
                          mpi_comm_world,status,ierr)

        nsendtotp=sum(nsendp)
        nrecvtotp=sum(nrecvp)
        
        ppacket(1) = nsendtotp
        ppacket(2) = nrecvtotp
        ppacket(3) = n_fast_removed_local
        call MPI_ALLREDUCE(ppacket,ppacketg,3,MPI_INTEGER4,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
        nsendtot       = ppacketg(1)
        nrecvtot       = ppacketg(2)
        n_fast_removed = ppacketg(3)

        ! if ((myid == 0).and.mod(it,n_print)==0) then
        !   print*,  " Finished compiling lists "
        !   print*, " # of particles to be sent     = ",nsendtot
        !   print*, " # of particles to be received = ",nrecvtot
        !   print*, " # of particles removed because V > Vlimit = ", n_fast_removed
        ! endif

        if (nsendtot.ne.nrecvtot) THEN
          call MPI_FINALIZE(IERR)
          print*, "Error: nsendtot /= nrecvtot. Terminating"
          stop
        endif

        ! VR it's convenient to have neighboring processes in an array
        ! Eventually, we should change all of the code to use this
        nbrs(1) = NBRTOP
        nbrs(2) = NBRLEFTTOP
        nbrs(3) = NBRLEFT
        nbrs(4) = NBRLEFTBOT
        nbrs(5) = NBRBOT
        nbrs(6) = NBRRITEBOT
        nbrs(7) = NBRRITE
        nbrs(8) = NBRRITETOP

        ! VR: maximum numbers of particles to be sent/received to/from a process
        max_nsend = maxval(nsendp)
        max_nrecv = maxval(nrecvp)

        ! VR: allocate tmp arrays for sending/recieving data
        allocate(packed_pdata_send(8,max_nsend,8))
        allocate(packed_pdata_recv(8,max_nrecv))

        ! VR: counters for particles sent/recieved
        nsendactualp = 0
        nrecvactualp = 0

        ! VR: 4 stages of data exhcnage in a simple even->odd, odd->even 2D pattern
        do irepeat=1,4
          if (isendid(irepeat) == 1) then
            nsend_to_nbr = 0
            NP=IPSEND(IS)
            
            ! loop over particles in the ipsend list
            do while (NP.ne.0)           
              nsendactualp = nsendactualp + 1

              ! map this particle to the logical mesh
              !             Nonuniform mesh - using mesh_unmap
              rxe=dtxi*mesh_unmap(meshX,x(np))+1.50000000000d+00
              rye=dtyi*mesh_unmap(meshY,y(np))+1.50000000000d+00
              rze=dtzi*mesh_unmap(meshZ,z(np))+1.50000000000d+00
              ixe=rxe
              iye=rye
              ize=rze
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0
                
              ypart=y(np)
              zpart=z(np)
              if (ypart <= ye.and.ypart >= yb) then
                iye_cc=jb 
              else
                if (ypart > ye) then
                    iye_cc=je+1 
                else
                    iye_cc=jb-1 
                endif
              endif
              if (zpart <= ze.and.zpart >= zb) then
                ize_cc=kb 
              else
                if (zpart > ze) then
                    ize_cc=ke+1 
                else
                    ize_cc=kb-1 
                endif
              endif
                
              ! destination process for this particle
              i_source = idmap_yz(iye_cc,ize_cc)
              
              if (i_source==NBRTOP) then
                idest = 1
              else if (i_source==NBRLEFTTOP) then
                idest = 2
              else if (i_source==NBRLEFT) then
                idest = 3
              else if (i_source==NBRLEFTBOT) then
                idest = 4
              else if (i_source==NBRBOT) then
                idest = 5
              else if (i_source==NBRRITEBOT) then
                idest = 6
              else if (i_source==NBRRITE) then
                idest = 7
              else if (i_source==NBRRITETOP) then
                idest = 8
              else
                print*, 'myid = ', myid, ": trying to send particles to procsees that are not immediate neighbors"
              endif
                
              !VR pack the data
              nsend_to_nbr(idest) = nsend_to_nbr(idest) + 1
              ii = nsend_to_nbr(idest)  
              packed_pdata_send(1,ii,idest)=x(np)
              packed_pdata_send(2,ii,idest)=y(np)
              packed_pdata_send(3,ii,idest)=z(np)
              packed_pdata_send(4,ii,idest)=vx(np)
              packed_pdata_send(5,ii,idest)=vy(np)
              packed_pdata_send(6,ii,idest)=vz(np)
              packed_pdata_send(7,ii,idest)=qp(np)
              packed_pdata_send(8,ii,idest)=DBLE(ptag(np))
          
              np = link(np) !VR next particle in the list
            enddo

            !VR now send all the data
            do idest=1,8 ! 
              if (nsend_to_nbr(idest) > 0) then
                pp => packed_pdata_send(:,:,idest)
                call MPI_ISEND(pp,8*nsend_to_nbr(idest), &
                      MPI_DOUBLE_PRECISION, nbrs(idest),0,MPI_COMM_WORLD,exchange_send_request(idest),IERR)
              else
                exchange_send_request(idest) = MPI_REQUEST_NULL
              endif
            enddo
              
            !VR: wait for all requests to complete
            call MPI_WAITALL(8,exchange_send_request,status_array,ierr)

          else

            !VR this process recieves data in this stage from up to 4 different processes
            do itmp=1,4
              !VR loop over possible neighbors
              ipe=irecvid(itmp,irepeat)
              if (ipe.eq.-1) cycle
              if (nrecvp(ipe) == 0) cycle
              
              !VR get all the data
              call MPI_RECV(packed_pdata_recv,int(8*nrecvp(ipe),4),MPI_DOUBLE_PRECISION,&
                  ipe,0, MPI_COMM_WORLD,STATUS2,IERR)
            
              !VR unpack the data

              do ii=1,nrecvp(ipe)
                nprecv=ipstore
                          
                if (ipstore == 0) then
                    Storage_Error_p = 1
                    exit
                endif
                    
                nrecvactualp=nrecvactualp+1                       
                x (nprecv)=packed_pdata_recv(1,ii)
                y (nprecv)=packed_pdata_recv(2,ii)
                z (nprecv)=packed_pdata_recv(3,ii)
                vx(nprecv)=packed_pdata_recv(4,ii)
                vy(nprecv)=packed_pdata_recv(5,ii)
                vz(nprecv)=packed_pdata_recv(6,ii)
                qp(nprecv)=packed_pdata_recv(7,ii)
                ptag(nprecv)=INT(packed_pdata_recv(8,ii))
                    
                !VR: check y & z global domain boundaries and loop particles in case of periodicity
                if ( y(nprecv) < zero) then
                    y(nprecv) = y(nprecv) + ymax
                endif
                if ( y(nprecv) > ymax) then
                    y(nprecv) = y(nprecv) - ymax
                endif
                !VR: the same for z boundary
                if ( z(nprecv) < zero) then
                    z(nprecv) = z(nprecv) + zmax
                endif
                if ( z(nprecv) > zmax) then
                    z(nprecv) = z(nprecv) - zmax
                endif

                !VR: end of boundary check ----------------------------------
                ! Uniform mesh - Same as in version 5.0
                ! ixe=hxi*x(nprecv)+1.5000000000000001d+00
                ! iye=hyi*y(nprecv)+0.5000000000000001d+00
                ! ize=hzi*z(nprecv)+0.5000000000000001d+00

                ! Nonuniform mesh - using mesh_unmap
                rxe=dtxi*mesh_unmap(meshX,x(nprecv))+1.50000000000d+00
                rye=dtyi*mesh_unmap(meshY,y(nprecv))+1.50000000000d+00
                rze=dtzi*mesh_unmap(meshZ,z(nprecv))+1.50000000000d+00
                ixe=rxe
                iye=rye
                ize=rze
                iye=iye-1             ! integer index in y direction starts at 0
                ize=ize-1             ! integer index in z direction starts at 0
                !VR head into the list of "empty" particles
                ipstore=link(nprecv)
                    
                if ((ixe>nx+1 .or. ixe<1) .or. (iye>je+1 .or. iye<jb-1) .or. (ize>ke+1 .or. ize < kb-1)) then
                    Field_Diverge_p = 1
                    ixe = min(max(iye,1_8 ),nx+1)
                    iye = min(max(iye,jb-1),je+1)
                    ize = min(max(ize,kb-1),ke+1)
                endif
                    
                !VR link this particle into the list of active parrticles
                link(nprecv)=iphead(ixe,iye,ize,is)
                iphead(ixe,iye,ize,is)=nprecv
                    
              enddo
              !VR: prevent a situation where we attempt to receive more than once from the same processo
              !VR  e.g. if NDIM=2
              nrecvp(ipe) = 0  
            enddo
          endif
        enddo  ! end of irepeat loop
        
        ! VR: deallocate packed arrays
        deallocate(packed_pdata_send)
        deallocate(packed_pdata_recv)
        nullify(pp)
        
        ! end of packed exchange        
        dpacket(1) = Storage_Error_p
        dpacket(2) = Field_Diverge_p
        dpacket(3) = nsendactualp
        dpacket(4) = nrecvactualp
        call MPI_ALLREDUCE(dpacket,dpacketg,4,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,IERR)
        Storage_Error = dpacketg(1)
        Field_Diverge = dpacketg(2)
        nsendactual   = dpacketg(3)
        nrecvactual   = dpacketg(4)

        if (Storage_Error /= 0) then
           if (myid == 0) then
             print*, " "
             print*, "Particle storage allocation is exceeded."
             print*,  "H3D is stopped"
             print*, " "
           endif
           call MPI_FINALIZE(IERR)
           stop
        endif

        if (Field_Diverge /= 0) then
           if (myid == 0) then
             print*, " "
             print*, "Field Solver Diverges"
             print*, "3DHybrid is stopped"
             print*," "
           endif
           call MPI_FINALIZE(IERR)
           STOP
        endif

        ! if ((myid == 0).and.mod(it,n_print)==0) then
        !   print*, " Finished exchanging particles"
        !   print*, " # of particles sent     = ", nsendactual
        !   print*, " # of particles received = ", nrecvactual
        ! endif

      enddo ! for IS
 
end subroutine particle_boundary


!---------------------------------------------------------------------
! sort the particles
!---------------------------------------------------------------------
subroutine sortit
  use parameter_mod
  use mesh_mod
  implicit none

  real*8 :: pstore(nplmax)
  integer :: pstore2(nplmax)
  integer*8 :: id, kb1, is, ix, iy, iz, ixe ,iye, ize, l, nttot, nplist
  real*8 :: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
  
  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt

  id = 0
  kb1 = kb-1
  iptemp = 0
  porder = 0
  do is = 1,nspec
    do iz = kb-1,ke
      do iy = jb-1,je
        do ix = 1, nx1
          np = iphead(ix,iy,iz,is)
          do while (np.ne.0)
            ! Uniform mesh - Same as in version 5.0
            ! ixe = hxi*x(np)+1.5000000000000001d+00
            ! iye = hyi*y(np)+0.5000000000000001d+00
            ! ize = hzi*z(np)+0.5000000000000001d+00

            ! Nonuniform mesh - using mesh_unmap
            rxe=dtxi*mesh_unmap(meshX,x(np))+1.50000000000d+00
            rye=dtyi*mesh_unmap(meshY,y(np))+1.50000000000d+00
            rze=dtzi*mesh_unmap(meshZ,z(np))+1.50000000000d+00
            ixe=rxe
            iye=rye
            ize=rze
            iye=iye-1    ! integer index in y direction starts at 0
            ize=ize-1    ! integer index in z direction starts at 0

            porder(np)=iptemp(ixe,iye,ize,is)
            iptemp(ixe,iye,ize,is)=np
            np = link(np)
          enddo
        enddo
      enddo
    enddo
  enddo

  l = 0
  nttot = 0
  do is = 1,nspec
    do iz = kb-1,ke
      do iy = jb-1,je
        do ix = 1, nx1
          np = iptemp(ix,iy,iz,is)
          nplist = 0
          do while (NP.ne.0)
            nplist = nplist+1
            l = l+1
            link(l) = np
            np = porder(np)
          enddo
          nttot = nttot + nplist
          iphead(ix,iy,iz,is) = nplist
        enddo
      enddo
    enddo
  enddo

  id = 0
  kb1 = kb-1
  do l = 1,nttot
    pstore(l) = vx(link(l))
  enddo
  do l = 1,nttot
    vx(l) = pstore(l)
  enddo

  do l = 1,nttot
    pstore(l) = vy(link(l))
  enddo
  do l = 1,nttot
    vy(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = vz(link(l))
  enddo
  do l = 1,nttot
    vz(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = x(link(l))
  enddo
  do l = 1,nttot
    x(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = y(link(l))
  enddo
  do l = 1,nttot
    y(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = z(link(l))
  enddo
  do l = 1,nttot
    z(l) = pstore(l)
  enddo

  do l = 1,nttot
    pstore(l) = qp(link(l))
  enddo
  do l = 1,nttot
    qp(l) = pstore(l)
  enddo

  do l = 1,nttot
    pstore2(l) = ptag(link(l))
  enddo
  do l = 1,nttot
    ptag(l) = pstore2(l)
  enddo

  l=1
  do is = 1,nspec
    do iz = kb-1,ke
      do iy = jb-1,je
        do ix = 1, nx1
          nplist = iphead(ix,iy,iz,is)
          if (nplist.ne.0) then
            iphead(ix,iy,iz,is) = l
            do np = l, l+nplist-1
              link(np) = np+1
            enddo
            link(l+nplist-1) = 0
            l = l + nplist
          else
            iphead(ix,iy,iz,is) = 0
          endif
        enddo
      enddo
    enddo
  enddo

  if (l-1.ne.nttot) then
    print *,'Problem in SORT: l-1 NE NTTOT'
    stop
  endif
  ipstore = nttot + 1
  do l = nttot+1, nplmax-1
    link(l) = l+1
  enddo
  link(nplmax) = 0

  return
end subroutine sortit


!-----------------------------------------------------------------
! computes velocities?
! what is the difference between vxs and vix?
!-----------------------------------------------------------------
subroutine trans        
  use parameter_mod
  use mesh_mod
  implicit none

  integer*8 :: is, i, j, k, jbmin, jbmax, kbmin, kbmax
  integer*4 :: time_begin(8), time_end(8)
  real*8 :: dns_tmp

  ! timing whole 'trans'
  call date_and_time(values=time_begin_array(:,20))

  ! initialize sth
  do is = 1, nspec
    do k = kb-1, ke+1
      do j = jb-1, je+1
        do i = 1, nx2
          dns(i,j,k,is)=1.e-10
          dnsh(i,j,k,is)=1.e-10
          vxs(i,j,k,is)=0.
          vys(i,j,k,is)=0.
          vzs(i,j,k,is)=0.
          if (is==1) then
            deno(i,j,k)=den(i,j,k)
            vixo(i,j,k)=vix(i,j,k)
            viyo(i,j,k)=viy(i,j,k)
            vizo(i,j,k)=viz(i,j,k)
            den(i,j,k)=0.
            denh(i,j,k)=0.
            vix(i,j,k)=0.
            viy(i,j,k)=0.
            viz(i,j,k)=0.
          endif
        enddo
      enddo
    enddo
  enddo

  ! push particles
  call date_and_time(values=time_begin_array(:,7))
  if (ndim /= 1) then
    call parmov
  else
    call parmov_2d
  endif
  call date_and_time(values=time_end_array(:,7))
  call accumulate_time(time_begin_array(1,7),time_end_array(1,7),time_elapsed(7))

  ! if test orbits, just return by here?
  if (testorbt) return

  ! what
  do is = 1, nspec
    do k = kb-1, ke+1
      do j = jb-1, je+1
        do i = 1, nx2 
          den(i,j,k)=den(i,j,k)+dns(i,j,k,is)*qspec(is) 
          denh(i,j,k)=denh(i,j,k)+dnsh(i,j,k,is)*qspec(is) 
          vix(i,j,k)=vix(i,j,k)+qspec(is)*vxs(i,j,k,is) 
          viy(i,j,k)=viy(i,j,k)+qspec(is)*vys(i,j,k,is) 
          viz(i,j,k)=viz(i,j,k)+qspec(is)*vzs(i,j,k,is)
        enddo
      enddo
    enddo
  enddo

  ! apply boundary conditions
  if (ndim /= 1) then
    call xrealbcc(den,1_8,nx,ny,nz)
    call xrealbcc(denh,1_8,nx,ny,nz)
    call xrealbcc(vix,1_8,nx,ny,nz)
    call xrealbcc(viy,1_8,nx,ny,nz)
    call xrealbcc(viz,1_8,nx,ny,nz)
  else
    call xrealbcc_2d(den,1_8,nx,ny,nz)
    call xrealbcc_2d(denh,1_8,nx,ny,nz)
    call xrealbcc_2d(vix,1_8,nx,ny,nz)
    call xrealbcc_2d(viy,1_8,nx,ny,nz)
    call xrealbcc_2d(viz,1_8,nx,ny,nz)
  endif

  ! smooth density and velocity
  if (smoothing) then
    if (ndim /=1) then
      call nsmth(den)
      call nsmth(denh)
      call nsmth(vix)
      call nsmth(viy)
      call nsmth(viz)
    else
      call nsmth_2d(den, nx2, ny2, nz2)
      call nsmth_2d(denh, nx2, ny2, nz2)
      call nsmth_2d(vix, nx2, ny2, nz2)
      call nsmth_2d(viy, nx2, ny2, nz2)
      call nsmth_2d(viz, nx2, ny2, nz2)
    endif
  endif

  do k = kb-1, ke+1
    do j = jb-1, je+1
      do i = 1, nx2
        den(i,j,k)=max(denmin,den(i,j,k))
        ! pe(i,j,k) =te0*den(i,j,k)**gamma ! pressure calculation moved to pressgrad
        vix(i,j,k)=vix(i,j,k)/denh(i,j,k)
        viy(i,j,k)=viy(i,j,k)/denh(i,j,k)
        viz(i,j,k)=viz(i,j,k)/denh(i,j,k)
      enddo
    enddo
  enddo

  ! for 1st step?
  if (it == 0) then
     deno=den; vixo=vix; viyo=viy; vizo=viz
  endif

  ! what does 'energy' do?
  call date_and_time(values=time_begin_array(:,8))
  if (mod(it,10)==0) call energy
  call date_and_time(values=time_end_array(:,8))
  call accumulate_time(time_begin_array(1,8),time_end_array(1,8),time_elapsed(8))

  kbmin = kb-1; kbmax = ke+1
  jbmin = jb-1; jbmax = je+1

  call date_and_time(values=time_end_array(:,20))
  call accumulate_time(time_begin_array(1,20),time_end_array(1,20),time_elapsed(20))

  return
end subroutine trans


!-----------------------------------------------------------------
! compute perp and par temperature and pressure tensor
!-----------------------------------------------------------------
subroutine caltemp2_global
  use parameter_mod
  use mesh_mod
  implicit none

  real*8 :: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
  integer*8 :: ix,iy,iz,ixp1,iyp1,izp1,iiy,iiye,iiz,iize,is,l,iix,iixe
  real*8 :: vxa,vya,vza,rfrac,vxavg,vxavg1,vxavg2 &
            ,vyavg,vyavg1,vyavg2,vzavg,vzavg1,vzavg2,wperp2,wpar,wmult
  real*8 :: w1,w2,w3,w4,w5,w6,w7,w8,h,hh,dns1,dns2,bxa,bya,bza,btota,dnst

  call date_and_time(values=time_begin_array(:,23))

  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt

  tpar=0.; tperp=0.
  p_xx=0.; p_xy=0.; p_xz=0.
  p_yy=0.; p_yz=0.; p_zz=0.

  if (nspec >= 2) then
     rfrac = frac(2)/frac(1)
  else
     rfrac = 0.
  endif

  call date_and_time(values=time_begin_array(:,26))
  do IS=1,NSPEC
    wmult=wspec(is)
    h=dt*qspec(is)/wmult
    hh=.5*h
    dpedx = 0.
    do IIZE = KB-1,KE
      do IIYE = JB-1,JE
        do IIXE = 1, NX1
          NP=IPHEAD(IIXE,IIYE,IIZE,IS)
          !  begin advance of particle position and velocity
          !  If dt=0, skip
          do while (NP.ne.0)
            L=NP
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

            ! Nonuniform mesh - using mesh_unmap
            rx=dtxi*mesh_unmap(meshX,x(l))+1.50000000000d+00
            ry=dtyi*mesh_unmap(meshY,y(l))+1.50000000000d+00
            rz=dtzi*mesh_unmap(meshZ,z(l))+1.50000000000d+00
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

            dnst= dnsh(ix  ,iy  ,iz  ,is)*w1+dnsh(ixp1,iy  ,iz  ,is)*w2  &
            +     dnsh(ix  ,iyp1,iz  ,is)*w3+dnsh(ixp1,iyp1,iz  ,is)*w4  &
            +     dnsh(ix  ,iy  ,izp1,is)*w5+dnsh(ixp1,iy  ,izp1,is)*w6  &
            +     dnsh(ix  ,iyp1,izp1,is)*w7+dnsh(ixp1,iyp1,izp1,is)*w8

            vxavg=vxs(ix  ,iy  ,iz  ,is)*w1+vxs(ixp1,iy  ,iz  ,is)*w2  &
            +      vxs(ix  ,iyp1,iz  ,is)*w3+vxs(ixp1,iyp1,iz  ,is)*w4  &
            +      vxs(ix  ,iy  ,izp1,is)*w5+vxs(ixp1,iy  ,izp1,is)*w6  &
            +      vxs(ix  ,iyp1,izp1,is)*w7+vxs(ixp1,iyp1,izp1,is)*w8
            vxavg=vxavg/dnst
            

            vyavg=vys(ix  ,iy  ,iz  ,is)*w1+vys(ixp1,iy  ,iz  ,is)*w2  &
            +      vys(ix  ,iyp1,iz  ,is)*w3+vys(ixp1,iyp1,iz  ,is)*w4  &
            +      vys(ix  ,iy  ,izp1,is)*w5+vys(ixp1,iy  ,izp1,is)*w6  &
            +      vys(ix  ,iyp1,izp1,is)*w7+vys(ixp1,iyp1,izp1,is)*w8
            vyavg=vyavg/dnst
            

            vzavg=vzs(ix  ,iy  ,iz  ,is)*w1+vzs(ixp1,iy  ,iz  ,is)*w2  &
            +      vzs(ix  ,iyp1,iz  ,is)*w3+vzs(ixp1,iyp1,iz  ,is)*w4  &
            +      vzs(ix  ,iy  ,izp1,is)*w5+vzs(ixp1,iy  ,izp1,is)*w6  &
            +      vzs(ix  ,iyp1,izp1,is)*w7+vzs(ixp1,iyp1,izp1,is)*w8
            vzavg=vzavg/dnst


            vxa=vx(l)-vxavg
            vya=vy(l)-vyavg
            vza=vz(l)-vzavg

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
            wpar=(vxa*bxa+vya*bya+vza*bza)/btota
            wperp2=vxa**2+vya**2+vza**2-wpar**2
            xx=vxa*vxa
            xy=vxa*vya
            xz=vxa*vza
            yy=vya*vya
            yz=vya*vza
            zz=vza*vza

            tpar (ix  ,iy  ,iz  ,is)=tpar (ix  ,iy  ,iz  ,is)+qp(np)*w1*wpar*wpar
            tpar (ixp1,iy  ,iz  ,is)=tpar (ixp1,iy  ,iz  ,is)+qp(np)*w2*wpar*wpar 
            tpar (ix  ,iyp1,iz  ,is)=tpar (ix  ,iyp1,iz  ,is)+qp(np)*w3*wpar*wpar 
            tpar (ixp1,iyp1,iz  ,is)=tpar (ixp1,iyp1,iz  ,is)+qp(np)*w4*wpar*wpar 
            tpar (ix  ,iy  ,izp1,is)=tpar (ix  ,iy  ,izp1,is)+qp(np)*w5*wpar*wpar 
            tpar (ixp1,iy  ,izp1,is)=tpar (ixp1,iy  ,izp1,is)+qp(np)*w6*wpar*wpar 
            tpar (ix  ,iyp1,izp1,is)=tpar (ix  ,iyp1,izp1,is)+qp(np)*w7*wpar*wpar 
            tpar (ixp1,iyp1,izp1,is)=tpar (ixp1,iyp1,izp1,is)+qp(np)*w8*wpar*wpar 
            tperp(ix  ,iy  ,iz  ,is)=tperp(ix  ,iy  ,iz  ,is)+qp(np)*w1*wperp2 
            tperp(ixp1,iy  ,iz  ,is)=tperp(ixp1,iy  ,iz  ,is)+qp(np)*w2*wperp2 
            tperp(ix  ,iyp1,iz  ,is)=tperp(ix  ,iyp1,iz  ,is)+qp(np)*w3*wperp2 
            tperp(ixp1,iyp1,iz  ,is)=tperp(ixp1,iyp1,iz  ,is)+qp(np)*w4*wperp2
            tperp(ix  ,iy  ,izp1,is)=tperp(ix  ,iy  ,izp1,is)+qp(np)*w5*wperp2 
            tperp(ixp1,iy  ,izp1,is)=tperp(ixp1,iy  ,izp1,is)+qp(np)*w6*wperp2
            tperp(ix  ,iyp1,izp1,is)=tperp(ix  ,iyp1,izp1,is)+qp(np)*w7*wperp2 
            tperp(ixp1,iyp1,izp1,is)=tperp(ixp1,iyp1,izp1,is)+qp(np)*w8*wperp2
            dpedx(ix  ,iy  ,iz  )=dpedx(ix  ,iy  ,iz  )+qp(np)*w1
            dpedx(ixp1,iy  ,iz  )=dpedx(ixp1,iy  ,iz  )+qp(np)*w2 
            dpedx(ix  ,iyp1,iz  )=dpedx(ix  ,iyp1,iz  )+qp(np)*w3 
            dpedx(ixp1,iyp1,iz  )=dpedx(ixp1,iyp1,iz  )+qp(np)*w4 
            dpedx(ix  ,iy  ,izp1)=dpedx(ix  ,iy  ,izp1)+qp(np)*w5 
            dpedx(ixp1,iy  ,izp1)=dpedx(ixp1,iy  ,izp1)+qp(np)*w6 
            dpedx(ix  ,iyp1,izp1)=dpedx(ix  ,iyp1,izp1)+qp(np)*w7 
            dpedx(ixp1,iyp1,izp1)=dpedx(ixp1,iyp1,izp1)+qp(np)*w8 

            p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
            p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
            p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
            p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 
            p_xx (ix  ,iy  ,izp1,is)=p_xx (ix  ,iy  ,izp1,is)+qp(np)*w5*xx 
            p_xx (ixp1,iy  ,izp1,is)=p_xx (ixp1,iy  ,izp1,is)+qp(np)*w6*xx 
            p_xx (ix  ,iyp1,izp1,is)=p_xx (ix  ,iyp1,izp1,is)+qp(np)*w7*xx 
            p_xx (ixp1,iyp1,izp1,is)=p_xx (ixp1,iyp1,izp1,is)+qp(np)*w8*xx 

            p_xy (ix  ,iy  ,iz  ,is)=p_xy (ix  ,iy  ,iz  ,is)+qp(np)*w1*xy
            p_xy (ixp1,iy  ,iz  ,is)=p_xy (ixp1,iy  ,iz  ,is)+qp(np)*w2*xy 
            p_xy (ix  ,iyp1,iz  ,is)=p_xy (ix  ,iyp1,iz  ,is)+qp(np)*w3*xy 
            p_xy (ixp1,iyp1,iz  ,is)=p_xy (ixp1,iyp1,iz  ,is)+qp(np)*w4*xy 
            p_xy (ix  ,iy  ,izp1,is)=p_xy (ix  ,iy  ,izp1,is)+qp(np)*w5*xy 
            p_xy (ixp1,iy  ,izp1,is)=p_xy (ixp1,iy  ,izp1,is)+qp(np)*w6*xy 
            p_xy (ix  ,iyp1,izp1,is)=p_xy (ix  ,iyp1,izp1,is)+qp(np)*w7*xy 
            p_xy (ixp1,iyp1,izp1,is)=p_xy (ixp1,iyp1,izp1,is)+qp(np)*w8*xy 

            p_xz (ix  ,iy  ,iz  ,is)=p_xz (ix  ,iy  ,iz  ,is)+qp(np)*w1*xz
            p_xz (ixp1,iy  ,iz  ,is)=p_xz (ixp1,iy  ,iz  ,is)+qp(np)*w2*xz 
            p_xz (ix  ,iyp1,iz  ,is)=p_xz (ix  ,iyp1,iz  ,is)+qp(np)*w3*xz 
            p_xz (ixp1,iyp1,iz  ,is)=p_xz (ixp1,iyp1,iz  ,is)+qp(np)*w4*xz 
            p_xz (ix  ,iy  ,izp1,is)=p_xz (ix  ,iy  ,izp1,is)+qp(np)*w5*xz 
            p_xz (ixp1,iy  ,izp1,is)=p_xz (ixp1,iy  ,izp1,is)+qp(np)*w6*xz 
            p_xz (ix  ,iyp1,izp1,is)=p_xz (ix  ,iyp1,izp1,is)+qp(np)*w7*xz 
            p_xz (ixp1,iyp1,izp1,is)=p_xz (ixp1,iyp1,izp1,is)+qp(np)*w8*xz 

            p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
            p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
            p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
            p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
            p_yy (ix  ,iy  ,izp1,is)=p_yy (ix  ,iy  ,izp1,is)+qp(np)*w5*yy 
            p_yy (ixp1,iy  ,izp1,is)=p_yy (ixp1,iy  ,izp1,is)+qp(np)*w6*yy 
            p_yy (ix  ,iyp1,izp1,is)=p_yy (ix  ,iyp1,izp1,is)+qp(np)*w7*yy 
            p_yy (ixp1,iyp1,izp1,is)=p_yy (ixp1,iyp1,izp1,is)+qp(np)*w8*yy 

            p_yz (ix  ,iy  ,iz  ,is)=p_yz (ix  ,iy  ,iz  ,is)+qp(np)*w1*yz
            p_yz (ixp1,iy  ,iz  ,is)=p_yz (ixp1,iy  ,iz  ,is)+qp(np)*w2*yz 
            p_yz (ix  ,iyp1,iz  ,is)=p_yz (ix  ,iyp1,iz  ,is)+qp(np)*w3*yz 
            p_yz (ixp1,iyp1,iz  ,is)=p_yz (ixp1,iyp1,iz  ,is)+qp(np)*w4*yz 
            p_yz (ix  ,iy  ,izp1,is)=p_yz (ix  ,iy  ,izp1,is)+qp(np)*w5*yz 
            p_yz (ixp1,iy  ,izp1,is)=p_yz (ixp1,iy  ,izp1,is)+qp(np)*w6*yz 
            p_yz (ix  ,iyp1,izp1,is)=p_yz (ix  ,iyp1,izp1,is)+qp(np)*w7*yz 
            p_yz (ixp1,iyp1,izp1,is)=p_yz (ixp1,iyp1,izp1,is)+qp(np)*w8*yz 

            p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
            p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
            p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
            p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
            p_zz (ix  ,iy  ,izp1,is)=p_zz (ix  ,iy  ,izp1,is)+qp(np)*w5*zz 
            p_zz (ixp1,iy  ,izp1,is)=p_zz (ixp1,iy  ,izp1,is)+qp(np)*w6*zz 
            p_zz (ix  ,iyp1,izp1,is)=p_zz (ix  ,iyp1,izp1,is)+qp(np)*w7*zz 
            p_zz (ixp1,iyp1,izp1,is)=p_zz (ixp1,iyp1,izp1,is)+qp(np)*w8*zz 

            np=link(np)
          enddo
        enddo
      enddo
    enddo

    call xreal(tpar (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(tperp(1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(dpedx(1,jb-1,kb-1   ),nx,ny,nz)

    call xreal(p_xx (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_xy (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_xz (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_yy (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_yz (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_zz (1,jb-1,kb-1,is),nx,ny,nz)

    do IZ = KB-1,KE
      do IY = JB-1,JE
        do IX = 1, NX1
          if (dpedx(ix,iy,iz) /= 0.) then
            tpar (ix,iy,iz,is) = tpar (ix,iy,iz,is)/(   tx0(is)*dpedx(ix,iy,iz))
            tperp(ix,iy,iz,is) = tperp(ix,iy,iz,is)/(2.*tx0(is)*dpedx(ix,iy,iz))
          endif
        enddo
      enddo
    enddo

    do iiz = kb-1, ke+1
      do iiy = jb-1, je+1
        do iix = 1, nx2
          p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_xy(iix,iiy,iiz,is) = p_xy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_xz(iix,iiy,iiz,is) = p_xz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yz(iix,iiy,iiz,is) = p_yz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
        enddo
      enddo
    enddo

    ! p_xx(:,:,:,is)=p_xx(:,:,:,is)/(tx0(is)*frac(is))
    ! p_xy(:,:,:,is)=p_xy(:,:,:,is)/(tx0(is)*frac(is))
    ! p_xz(:,:,:,is)=p_xz(:,:,:,is)/(tx0(is)*frac(is))
    ! p_yy(:,:,:,is)=p_yy(:,:,:,is)/(tx0(is)*frac(is))
    ! p_yz(:,:,:,is)=p_yz(:,:,:,is)/(tx0(is)*frac(is))
    ! p_zz(:,:,:,is)=p_zz(:,:,:,is)/(tx0(is)*frac(is))

  enddo
  call date_and_time(values=time_end_array(:,26))
  call accumulate_time(time_begin_array(1,26),time_end_array(1,26),time_elapsed(26))

!  do is=1,nspec
!    call date_and_time(values=time_begin_array(:,24))
!    call xreal(tpar (1,jb-1,kb-1,is),nx,ny,nz)
!    call xreal(tperp(1,jb-1,kb-1,is),nx,ny,nz)
!    call date_and_time(values=time_end_array(:,24))
!    call accumulate_time(time_begin_array(1,24) &
! &                                 ,time_end_array(1,24) &
! &                                 ,time_elapsed(24))

!    call date_and_time(values=time_begin_array(:,25))
!    call xrealbcc(tpar (1,jb-1,kb-1,is),1,nx,ny,nz)
!    call xrealbcc(tperp(1,jb-1,kb-1,is),1,nx,ny,nz)
!    call date_and_time(values=time_end_array(:,25))
!    call accumulate_time(time_begin_array(1,25) &
! &                                 ,time_end_array(1,25) &
! &                                 ,time_elapsed(25))

!  enddo


!  call date_and_time(values=time_begin_array(:,26))
!  do is=1,nspec
!    do k=kb-1,ke+1
!      do j = jb-1,je+1
!        do i=1,nx2
!          if (is == 1) then
!            dns1=dns(i,j,k,1)/(dfac(1)*frac(1))
!            dns2=0.
!            denum=dns1+rfrac*dns2
!          else
!            denum=dns(i,j,k,is)/(dfac(is)*frac(is))
!          endif
!          if (denum < denmin)  then
!           tpar(i,j,k,is)=1.e-5
!           tperp(i,j,k,is)=1.e-5
!          else
!           denum=denum*tx0(is)
!           tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
!           tperp(i,j,k,is)=0.5*tperp(i,j,k,is)*wspec(is)/denum
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!  call date_and_time(values=time_end_array(:,26))
!  call accumulate_time(time_begin_array(1,26) &
! &                               ,time_end_array(1,26) &
! &                               ,time_elapsed(26))

   call date_and_time(values=time_end_array(:,23))
   call accumulate_time(time_begin_array(1,23),time_end_array(1,23),time_elapsed(23))

  return
end subroutine caltemp2_global


!********************************************************
! computes field energy ex^2+ey^2+ez^2 and bx^2+by^2+bz^2
! and particle energies
!********************************************************
subroutine energy
  use parameter_mod
  use mesh_mod
  implicit none
  real*8 :: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
  integer*8 ix,iy,iz,ixp1,iyp1,izp1,iiy,iiye,iiz,iize,is,l,iix,iixe
  real*8 :: vxa,vya,vza,rfrac,vxavg,vxavg1,vxavg2 &
        ,vyavg,vyavg1,vyavg2,vzavg,vzavg1,vzavg2,wperp2,wpar,wmult
  real*8 :: w1,w2,w3,w4,w5,w6,w7,w8,h,hh,dns1,dns2,bxa,bya,bza,btota,dnst
  real*8 :: bfldp,efldp,e_fluid,e_thermal,v2
  integer*8 :: i,j,k
  
  efldp=0.
  bfldp=0.
  e_fluid=0.
  e_thermal=0.
  v2=0.

  ! particle energy calculation -- works for 3D only !!!
  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt
  p_xx=0.; p_yy=0.; p_zz=0.

  do is = 1, 1 ! ??
    do IIZE = KB-1,KE
      do IIYE = JB-1,JE
        do IIXE = 1, NX1
          NP=IPHEAD(IIXE,IIYE,IIZE,IS)
          do while (NP.ne.0)
            L=NP
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

            ! Nonuniform mesh - using mesh_unmap
            rx=dtxi*mesh_unmap(meshX,x(l))+1.50000000000d+00
            ry=dtyi*mesh_unmap(meshY,y(l))+1.50000000000d+00
            rz=dtzi*mesh_unmap(meshZ,z(l))+1.50000000000d+00
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

            dns1= dns(ix  ,iy  ,iz  ,1)*w1+dns(ixp1,iy  ,iz  ,1)*w2  &
            +     dns(ix  ,iyp1,iz  ,1)*w3+dns(ixp1,iyp1,iz  ,1)*w4  &
            +     dns(ix  ,iy  ,izp1,1)*w5+dns(ixp1,iy  ,izp1,1)*w6  &
            +     dns(ix  ,iyp1,izp1,1)*w7+dns(ixp1,iyp1,izp1,1)*w8

            dns2= 0.

            dnst = dns1 + dns2
  
            vxavg1=vxs(ix  ,iy  ,iz  ,1)*w1+vxs(ixp1,iy  ,iz  ,1)*w2  &
            +      vxs(ix  ,iyp1,iz  ,1)*w3+vxs(ixp1,iyp1,iz  ,1)*w4  &
            +      vxs(ix  ,iy  ,izp1,1)*w5+vxs(ixp1,iy  ,izp1,1)*w6  &
            +      vxs(ix  ,iyp1,izp1,1)*w7+vxs(ixp1,iyp1,izp1,1)*w8

            vxavg2= 0.

            vxavg = (dns1*vxavg1 + dns2*vxavg2)/dnst

            vyavg1=vys(ix  ,iy  ,iz  ,1)*w1+vys(ixp1,iy  ,iz  ,1)*w2  &
            +      vys(ix  ,iyp1,iz  ,1)*w3+vys(ixp1,iyp1,iz  ,1)*w4  &
            +      vys(ix  ,iy  ,izp1,1)*w5+vys(ixp1,iy  ,izp1,1)*w6  &
            +      vys(ix  ,iyp1,izp1,1)*w7+vys(ixp1,iyp1,izp1,1)*w8

            vyavg2=0.

            vyavg = (dns1*vyavg1 + dns2*vyavg2)/dnst

            vzavg1=vzs(ix  ,iy  ,iz  ,1)*w1+vzs(ixp1,iy  ,iz  ,1)*w2  &
            +      vzs(ix  ,iyp1,iz  ,1)*w3+vzs(ixp1,iyp1,iz  ,1)*w4  &
            +      vzs(ix  ,iy  ,izp1,1)*w5+vzs(ixp1,iy  ,izp1,1)*w6  &
            +      vzs(ix  ,iyp1,izp1,1)*w7+vzs(ixp1,iyp1,izp1,1)*w8

            vzavg2=0.

            vzavg = (dns1*vzavg1 + dns2*vzavg2)/dnst

            vxa=vx(l)-vxavg
            vya=vy(l)-vyavg
            vza=vz(l)-vzavg

            xx=vxa*vxa
            yy=vya*vya
            zz=vza*vza

            p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
            p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
            p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
            p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 
            p_xx (ix  ,iy  ,izp1,is)=p_xx (ix  ,iy  ,izp1,is)+qp(np)*w5*xx 
            p_xx (ixp1,iy  ,izp1,is)=p_xx (ixp1,iy  ,izp1,is)+qp(np)*w6*xx 
            p_xx (ix  ,iyp1,izp1,is)=p_xx (ix  ,iyp1,izp1,is)+qp(np)*w7*xx 
            p_xx (ixp1,iyp1,izp1,is)=p_xx (ixp1,iyp1,izp1,is)+qp(np)*w8*xx 

            p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
            p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
            p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
            p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
            p_yy (ix  ,iy  ,izp1,is)=p_yy (ix  ,iy  ,izp1,is)+qp(np)*w5*yy 
            p_yy (ixp1,iy  ,izp1,is)=p_yy (ixp1,iy  ,izp1,is)+qp(np)*w6*yy 
            p_yy (ix  ,iyp1,izp1,is)=p_yy (ix  ,iyp1,izp1,is)+qp(np)*w7*yy 
            p_yy (ixp1,iyp1,izp1,is)=p_yy (ixp1,iyp1,izp1,is)+qp(np)*w8*yy 

            p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
            p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
            p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
            p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
            p_zz (ix  ,iy  ,izp1,is)=p_zz (ix  ,iy  ,izp1,is)+qp(np)*w5*zz 
            p_zz (ixp1,iy  ,izp1,is)=p_zz (ixp1,iy  ,izp1,is)+qp(np)*w6*zz 
            p_zz (ix  ,iyp1,izp1,is)=p_zz (ix  ,iyp1,izp1,is)+qp(np)*w7*zz 
            p_zz (ixp1,iyp1,izp1,is)=p_zz (ixp1,iyp1,izp1,is)+qp(np)*w8*zz 

            v2=v2+qp(l)*(vx(l)**2+vy(l)**2+vz(l)**2)

            np=link(np)
          enddo
        enddo
      enddo
    enddo


    call xreal(p_xx (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_yy (1,jb-1,kb-1,is),nx,ny,nz)
    call xreal(p_zz (1,jb-1,kb-1,is),nx,ny,nz)

    do iiz = kb-1, ke+1
      do iiy = jb-1, je+1
        do iix = 1, nx2
          p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
        enddo
      enddo
    enddo

  enddo

  ! field energy calculation
  do k=kb,ke
      do j=jb,je
        do i=2,nx1
            efldp=efldp+ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2
            bfldp=bfldp+bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
            ! particle energy for species 1 only
            e_fluid=e_fluid+dns(i,j,k,1)*(vxs(i,j,k,1)**2+vys(i,j,k,1)**2+vzs(i,j,k,1)**2)
            e_thermal=e_thermal+p_xx(i,j,k,1)+p_yy(i,j,k,1)+p_yy(i,j,k,1)
        enddo
      enddo
  enddo
  efldp=efldp*hx*hy*hz*0.5 ! assuming uniform grids
  bfldp=bfldp*hx*hy*hz*0.5
  e_fluid=e_fluid*hx*hy*hz*0.5
  e_thermal=e_thermal*hx*hy*hz*0.5
  v2=v2*0.5

  ! collect energies
  call MPI_ALLREDUCE(efldp,efld,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(bfldp,bfld,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(e_fluid,efluidt,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(e_thermal,ethermt,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(v2,eptclt,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  return
end subroutine energy