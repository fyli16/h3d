!#######################################################################
!
      subroutine parmov      
      use parameter_mod
      use MESH2D
      implicit none
      double precision bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8,by1,by2,by3,by4,by5,by6,by7,by8,bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8,bxa,bya,bza
      double precision ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ey1,ey2,ey3,ey4,ey5,ey6,ey7,ey8,ez1,ez2,ez3,ez4,ez5,ez6,ez7,ez8

      double precision d_ranf,deltime1,deltime2,epsilon,ff,h,hh
      double precision fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8,foxa
      double precision foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8,foya
      double precision foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8,foza
      double precision w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
      double precision vex,vey,vez,vmag,vx_tmp,vy_tmp,vz_tmp
      double precision p2xs,p2ys,p2zs,q_p,th
      double precision wmult

      integer*8 i,ii,iix,iixe,iiy,iiye,iiz,iize,irepeat,irepeatp,is,itmp
      integer*8 iv,iye_cc,ize_cc,j,jv,k,npleavingp,nprecv,nprecvtmp
      integer*8 Storage_Error_p,Storage_Error
      data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
      data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
      data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
      integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                 ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
      double precision:: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,xpart,ypart,zpart
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
      double precision:: x_disp,y_disp,z_disp 
      double precision:: dth ! dt/2
      double precision:: eps2,myranf,twopi,fluxran,vxa,vyz,vza
      INTEGER*8:: L, EXIT_CODE_P, EXIT_CODE
      integer*8:: n_fast_removed,n_fast_removed_local,nptot_max,Field_Diverge,Field_Diverge_p
      double precision:: tx,ty,tz,v_x,v_y,v_z  
      INTEGER*4 :: nescapearr(8),nescapearr_global(8)
      INTEGER*4 :: ppacket(3),ppacketg(3),dpacket(4),dpacketg(4)
      INTEGER*8 :: epacket(2),epacketg(2),loop
      double precision, dimension(3,nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bxyz_av
      double precision:: TEX1,TEX2,TEX3,TEX4,TEX5,TEX6,TEX7,TEX8  
      double precision:: TEY1,TEY2,TEY3,TEY4,TEY5,TEY6,TEY7,TEY8  
      double precision:: TEZ1,TEZ2,TEZ3,TEZ4,TEZ5,TEZ6,TEZ7,TEZ8  
      double precision:: mX_xa,mX_ta,mX_ca1,mX_ca2,mX_xb,mX_dtdx,mX_tb,mX_cb1,mX_cb2
      double precision:: mY_xa,mY_ta,mY_ca1,mY_ca2,mY_xb,mY_dtdx,mY_tb,mY_cb1,mY_cb2
      double precision:: mZ_xa,mZ_ta,mZ_ca1,mZ_ca2,mZ_xb,mZ_dtdx,mZ_tb,mZ_cb1,mZ_cb2

      integer,dimension(8) :: nsend_to_nbr, nbrs
      integer :: idest,max_nsend,max_nrecv
      double precision,dimension(:,:,:),allocatable,target :: packed_pdata_send
      double precision,dimension(:,:),allocatable,target :: packed_pdata_recv
      double precision, pointer :: pp(:,:)
      integer exchange_send_request(8)

!      double precision :: mp_elapsed


      call date_and_time(values=time_begin_array(:,19))


      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

      dth=dt/2

 
      epsilon= buffer_zone
      d_ranf=1./1001.
      twopi=2.*acos(-1.)

      eps2=1.d-25
 
 

 
      bx_av=0.;by_av=0.;bz_av=0.
      DO K = KB-1,KE
        DO J = JB-1,JE
          DO I = 1, NX1
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
      CALL XREALBCC_PACK_B(BX_AV,BY_AV,BZ_AV,1_8,NX,NY,NZ)

      if ((myid == 0).and.prntinfo) WRITE(6,*) " CALLING PARMOVE, NSPEC = ",NSPEC
!
!=======================================================================
!
!  initalize diagnostic variables that keep track of
!  particle number, injection, and escape
!
      deltime1 = 0.0
      deltime2 = 0.0
!      nptotp=0
      npleavingp=0
 
      if (dt .ne. 0) then
        ! advance particles for a first half step
        call date_and_time(values=time_begin_array(:,13))
        call push
        call date_and_time(values=time_end_array(:,13))
        call accumulate_time_difference(time_begin_array(1,13),time_end_array(1,13),time_elapsed(13))

        ! check particles
        call date_and_time(values=time_begin_array(:,14))
        call particle_boundary
        call date_and_time(values=time_end_array(:,14))
        call accumulate_time_difference(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))

        ! collect Vi, ni at half step
        call date_and_time(values=time_begin_array(:,15))
        DO IS=1,NSPEC
          NPTOTP = 0
          npart(is) = 0
          DO IIZ=KB-1,KE
            DO IIY=JB-1,JE
              DO IIX=1,NX1
                np=iphead(iix,iiy,iiz,is)
                do while (np.ne.0)
                  nptotp=nptotp+1           !count particles
                  npart(is) = npart(is) + 1 !count particles in each species
                  L=NP

                  q_p = qp(l)

  !              Nonuniform mesh - using MESH_UNMAP
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
              ENDDO !for iix
            ENDDO !for iiy
          ENDDO !for iiz

          call XREAL(DNSH(1,jb-1,kb-1,is),NX,NY,NZ)
          call XREAL(VXS(1,jb-1,kb-1,is),NX,NY,NZ)
          call XREAL(VYS(1,jb-1,kb-1,is),NX,NY,NZ)
          call XREAL(VZS(1,jb-1,kb-1,is),NX,NY,NZ)
          call XREALBCC(DNSH(1,jb-1,kb-1,is),1_8,NX,NY,NZ)
          call XREALBCC(VXS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)
          call XREALBCC(VYS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)
          call XREALBCC(VZS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)

          DO IIZ=KB-1,KE+1
            DO IIY=JB-1,JE+1
              DO IIX=1,NX2
                dnsh(iix,iiy,iiz,is) = dnsh(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
                vxs(iix,iiy,iiz,is) = vxs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
                vys(iix,iiy,iiz,is) = vys(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
                vzs(iix,iiy,iiz,is) = vzs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              ENDDO
            ENDDO
          ENDDO

        ENDDO ! for is
        call date_and_time(values=time_end_array(:,15))
        call accumulate_time_difference(time_begin_array(1,15),time_end_array(1,15),time_elapsed(15))

        ! advance particles for a seconnd half step
        call date_and_time(values=time_begin_array(:,13))
        DO IS=1, nspec
          DO IIZE = KB-1,KE
            DO IIYE = JB-1,JE
              DO IIXE = 1, NX1
                NP=IPHEAD(IIXE,IIYE,IIZE,IS)

                DO WHILE (NP.NE.0)
                  L=NP
                  x_disp = dth*vx(l)
                  y_disp = dth*vy(l)
                  z_disp = dth*vz(l)

                  x(l)=x(l)+ x_disp
                  y(l)=y(l)+ y_disp
                  z(l)=z(l)+ z_disp

                  
                  NP=LINK(NP)
                enddo ! while
              ENDDO ! for iiz
            ENDDO ! for iiy
          ENDDO ! for iix
        ENDDO ! for is
        call date_and_time(values=time_end_array(:,13))
        call accumulate_time_difference(time_begin_array(1,13),time_end_array(1,13),time_elapsed(13))

      endif ! dt<>0


      ! check particles
      call date_and_time(values=time_begin_array(:,14))
      call particle_boundary
      call date_and_time(values=time_end_array(:,14))
      call accumulate_time_difference(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))
!
!
! collect density
      call date_and_time(values=time_begin_array(:,15))
      DO IS=1,nspec
        if (testorbt) goto 10
        NPTOTP = 0
        npart(is) = 0
        DO IIZ=KB-1,KE
          DO IIY=JB-1,JE
            DO IIX=1,NX1
              np=iphead(iix,iiy,iiz,is)
              do while (np.ne.0)
                nptotp=nptotp+1           !count particles
                npart(is) = npart(is) + 1 !count particles in each species
                L=NP

                q_p = qp(l)


!               Uniform mesh - Same as in version 5.0
!                rx=hxi*x(l)+1.5000000000000001d+00
!                ry=hyi*y(l)+0.5000000000000001d+00
!                rz=hzi*z(l)+0.5000000000000001d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                ix=max(1,min(ix,nx1))
!                iy=max(jb-1,min(iy,je))
!                iz=max(kb-1,min(iz,ke))
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

!               Nonuniform mesh - without using MESH_UNMAP
!                rx=hxi*x(l)+1.500000000000000d+00
!                ry=hyi*y(l)+1.500000000000000d+00
!                rz=hzi*z(l)+1.500000000000000d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                ix=ixc_2_c_map(ix)
!                iy=iyc_2_c_map(iy)
!                iz=izc_2_c_map(iz)
!                fx=(x(l)-meshX%xc(ix))/meshX%dxn(ix+1)
!                fy=(y(l)-meshY%xc(iy))/meshY%dxn(iy+1)
!                fz=(z(l)-meshZ%xc(iz))/meshZ%dxn(iz+1)

!              Nonuniform mesh - using MESH_UNMAP
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


!                 Diagnosic test on particle cell index algorithm
!                  if (     fx < 0. .or. fx > 1.                &
!                      .or. fy < 0. .or. fy > 1.                &
!                      .or. fz < 0. .or. fz > 1.) then 
!                      write(6,*) " GATHER LOOP "
!                      write(6,*) fx,fy,fz
!                      write(6,*) " x;",x(l),meshX%xn(ix),meshX%xc(ix)
!                      write(6,*) " y;",y(l),meshY%xn(iy),meshY%xc(iy)
!                      write(6,*) " z;",z(l),meshZ%xn(iz),meshZ%xc(iz)
!                      write(6,*) " r; ",rx,ry,rz
!                      write(6,*) " i; ",ixe,iye,ize
!                      write(6,*) " ixc_map; ",ix,iy,iz
!                  endif
 
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
              ENDDO
            ENDDO
          ENDDO
        ENDDO

 10     CONTINUE

        nescapearr(1) = nescape(is)
        nescapearr(2) = nescape_yz(is)
        nescapearr(3) = nescape_zy(is)
        nescapearr(4) = nescape_xy(is)
        nescapearr(5) = nescape_yx(is)
        nescapearr(6) = nescape_zx(is)
        nescapearr(7) = nescape_xz(is)
        nescapearr(8) = npart(is)
        call MPI_ALLREDUCE(nescapearr,nescapearr_global,8,MPI_INTEGER4,&
                         MPI_SUM,COMM2D,IERR)
        nescape_global(is)    = nescapearr_global(1)
        nescape_yz_global(is) = nescapearr_global(2)
        nescape_zy_global(is) = nescapearr_global(3)
        nescape_xy_global(is) = nescapearr_global(4)
        nescape_yx_global(is) = nescapearr_global(5)
        nescape_zx_global(is) = nescapearr_global(6)
        nescape_xz_global(is) = nescapearr_global(7)
        npart_global(is)      = nescapearr_global(8)

        deltime2 = deltime2 + real(clock_time1-clock_time)
 
        call XREAL(DNS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREALBCC(DNS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)


        DO IIZ=KB-1,KE+1
          DO IIY=JB-1,JE+1
            DO IIX=1,NX2
              dns(iix,iiy,iiz,is) = dns(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            ENDDO
          ENDDO
        ENDDO

      ENDDO ! for is
      call date_and_time(values=time_end_array(:,15))
      call accumulate_time_difference(time_begin_array(1,15),time_end_array(1,15),time_elapsed(15))
!
!
      ! diagnostic info

      epacket(1) = nptotp
      epacket(2) = npleavingp
      call MPI_ALLREDUCE(epacket,epacketg,2,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
      nptot      = epacketg(1)
      npleaving  = epacketg(2)

      if (myid == 0.and.prntinfo) then
        do is=1,nspec

          if (is == 1) then
            write(6,*)
            write(6,*) " it = ",it
            write(6,*) "  species #    ninj     nescape     ntot  "
          endif

          write(6,1000) is,ninj_global(is),nescape_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_yz_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_zy_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_xy_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_yx_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_xz_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_zx_global(is)&
                       ,npart_global(is)
 1002     format(3i6,5x,i8,2x,i8,2x,i10)
 1000     format(5x,i2,5x,i8,2x,i8,2x,i10)
        enddo

        if (nspec >= 2) then
          write(6,1005) ninj_global(1)+ninj_global(2),nescape_global(1)+nescape_global(2),npart_global(1)+npart_global(2)
        else
          write(6,1005) ninj_global(1),nescape_global(1),npart_global(1)
        endif

 1005   format(5x,'sum',4x,i8,2x,i8,2x,i10)
 1006   format(2i6,4x,'sum',4x,i8,2x,i8,2x,i10)
      endif


      ninj        = 0
      ninj_global = 0
 
      call date_and_time(values=time_end_array(:,19))
      call accumulate_time_difference(time_begin_array(1,19),time_end_array(1,19),time_elapsed(19))
 

      return
    end subroutine parmov

      subroutine push
      ! This subourtine pushes particles for half a step
      use parameter_mod
      use MESH2D
      implicit none
      integer*8 :: is,iixe,iiye,iize,l
      integer*8:: ix,iy,iz,ixe,iye,ize,ixep1,iyep1,izep1,ixp1,iyp1,izp1
      double precision :: wmult, h, hh, dth
      double precision bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8,by1,by2,by3,by4,by5,by6,by7,by8,bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8,bxa,bya,bza
      double precision ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ey1,ey2,ey3,ey4,ey5,ey6,ey7,ey8,ez1,ez2,ez3,ez4,ez5,ez6,ez7,ez8,exa,eya,eza
      double precision fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8,foxa
      double precision foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8,foya
      double precision foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8,foza
      double precision w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
      double precision vex,vey,vez
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
      double precision p2xs,p2ys,p2zs,ff
      double precision:: x_disp,y_disp,z_disp,disp_max_p(3),disp_max(3)&
                        ,y_disp_max_p,x_disp_max_p,z_disp_max_p,y_disp_max,x_disp_max,z_disp_max
      integer*8 :: Courant_Violation,Courant_Violation_p 

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

      DO IS=1,NSPEC
        Courant_Violation_p = 0
        x_disp_max_p        = 0
        y_disp_max_p        = 0
        z_disp_max_p        = 0
 
        wmult=wspec(is)
        h=dt*qspec(is)/wmult
        hh=.5*h
        dth=dt/2
 
        DO IIZE = KB-1,KE
          DO IIYE = JB-1,JE
            DO IIXE = 1, NX1
              NP=IPHEAD(IIXE,IIYE,IIZE,IS)
              DO WHILE (NP.NE.0)
                L=NP

!                 Uniform mesh - Same as in version 5.0
!                  rxe=hxi*x(l)+1.5000000000000001d+00
!                  rye=hyi*y(l)+0.5000000000000001d+00
!                  rze=hzi*z(l)+0.5000000000000001d+00
!                  ixe=rxe
!                  iye=rye
!                  ize=rze
!                  ixe=max(1   ,min(ixe,nx1))
!                  iye=max(jb-1,min(iye,je))
!                  ize=max(kb-1,min(ize,ke))
!                  ixep1 = ixe+1
!                  iyep1 = iye+1
!                  izep1 = ize+1
!                  fxe=rxe-ixe
!                  fye=rye-iye
!                  fze=rze-ize

!                 Nonuniform mesh - without using MESH_UNMAP
!                  rxe=hxi*x(l)+1.500000000000000d+00
!                  rye=hyi*y(l)+1.500000000000000d+00
!                  rze=hzi*z(l)+1.500000000000000d+00
!                  ixe=rxe
!                  iye=rye
!                  ize=rze
!                  ixe=ixc_2_c_map(ixe)
!                  iye=iyc_2_c_map(iye)
!                  ize=izc_2_c_map(ize)
!                  fxe=(x(l)-meshX%xc(ixe))/meshX%dxn(ixep1)
!                  fye=(y(l)-meshY%xc(iye))/meshY%dxn(iyep1)
!                  fze=(z(l)-meshZ%xc(ize))/meshZ%dxn(izep1)

!                 Nonuniform mesh - using MESH_UNMAP
                rxe=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
                rye=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
                rze=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
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

!                 Diagnosic test on particle cell index algorithm
!                  if (     fxe < 0. .or. fxe > 1.                &
!                      .or. fye < 0. .or. fye > 1.                &
!                      .or. fze < 0. .or. fze > 1.) then 
!                      write(6,*) " SCATTER LOOP"
!                      write(6,*) fxe,fye,fze
!                      write(6,*) " x;",x(l),meshX%xn(ixe),meshX%xc(ixe)
!                      write(6,*) " y;",y(l),meshY%xn(iye),meshY%xc(iye)
!                      write(6,*) " z;",z(l),meshZ%xn(ize),meshZ%xc(ize)
!                      write(6,*) " r; ",rxe,rye,rze
!                      write(6,*) " ixc_map; ",ixe,iye,ize
!                      write(6,*) " xc     ; ",meshX%xc(ixe),meshY%xc(iye),meshZ%xc(ize)
!                  endif


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

              ENDDO ! while
            ENDDO
          ENDDO
        ENDDO
        disp_max_p(1) = x_disp_max_p
        disp_max_p(2) = y_disp_max_p
        disp_max_p(3) = z_disp_max_p
        call MPI_ALLREDUCE(disp_max_p,disp_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)   !LAURA
        x_disp_max = disp_max(1)
        y_disp_max = disp_max(2)
        z_disp_max = disp_max(3)

        if ((myid == 0).and.prntinfo) then
          write(6,*) " maximum x-displacement/dx = ",x_disp_max
          write(6,*) " maximum y-displacement/dy = ",y_disp_max
          write(6,*) " maximum z-displacement/dz = ",z_disp_max
        endif

        call MPI_ALLREDUCE(Courant_Violation_p,Courant_Violation,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

        if (Courant_Violation /= 0) then
           if (myid == 0) write(6,*) "Particle displacements exceed cell size ",Courant_Violation," times"
           call MPI_FINALIZE(IERR)
           STOP
        endif
      ENDDO !is

      end subroutine push

      subroutine particle_boundary
      use parameter_mod
      use MESH2D
      implicit none
      INTEGER*4 :: ppacket(3),ppacketg(3),dpacket(4),dpacketg(4)
      INTEGER*8 :: epacket(2),epacketg(2),loop
      double precision:: xpart,ypart,zpart
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
      integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                 ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
      integer*8 :: Storage_Error_p,Storage_Error
      integer*8 :: l
      integer,dimension(8) :: nsend_to_nbr, nbrs
      integer :: idest,max_nsend,max_nrecv
      double precision,dimension(:,:,:),allocatable,target :: packed_pdata_send
      double precision,dimension(:,:),allocatable,target :: packed_pdata_recv
      double precision, pointer :: pp(:,:)
      integer :: exchange_send_request(8)
      integer*8 :: iv,iye_cc,ize_cc,j,jv,k,npleavingp,nprecv,nprecvtmp
      double precision:: v_limit
      integer*8:: n_fast_removed,n_fast_removed_local,nptot_max,Field_Diverge,Field_Diverge_p
      integer*8 i,ii,iix,iixe,iiy,iiye,iiz,iize,irepeat,irepeatp,is,itmp
      double precision:: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min
      Storage_Error_p = 0
      Field_Diverge_p = 0
      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

      !determine the velocity limit based on CFL condition

!     Uniform mesh - Same as in version 5.0
!      if (dt /=0.) then
!        if (ndim == 1) then
!         v_limit = min(hx/dt,hy/dt)
!        else
!         v_limit = min(hx/dt,hy/dt,hz/dt)
!      endif
!      else
!        v_limit=1.d+10
!      endif
!

!     Nonuniform mesh
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

      DO IS=1, nspec
           
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

        DO IZE = KB-1,KE
          DO IYE = JB-1,JE
            DO IXE = 1,NX1
              iptemp(ixe,iye,ize,is)=0
!
!=======================================================================
!
!     mark particles that need to be sent to other processors
!
              NP=IPHEAD(IXE,IYE,IZE,IS)  !VR: the first particle in the cell
              DO WHILE (NP.NE.0)         !VR: loop over particles in the cell
                xpart=x(np)
                ypart=y(np)
                zpart=z(np)
                l=np
!
!     Remove fast particles
!
                if (abs(vx(l)) >= v_limit .or. abs(vy(l)) >= v_limit .or. abs(vz(l)) >= v_limit) then
                  n_fast_removed_local = n_fast_removed_local + 1
                  iphead(ixe,iye,ize,is)=link(np)    !VR iphead is now the next particle in the list. Note that iphead
                                                     !VR will be re-assigned after the loop over particles in the cell
                  link(np)=ipstore                   !VR current particle points to the old start of empty list
                  ipstore=np                         !VR start of empty list is this particle
                  np=iphead(ixe,iye,ize,is)          !VR go to the next particle
                  goto 15
                endif
!
!
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
                ! VR: --------------------------------------------------

                xpart = x(l)
                ypart = y(l)
                zpart = z(l)

 
               !  if (    xpart < 0..or.xpart > xmax&
               !      .or.ypart < 0..or.ypart > ymax&       
               !      .or.zpart < 0..or.zpart > zmax) then

               !     ! VR: particle outside of the domain: I will get rid of this section
               !    if (xpart < 0.) then
               !      nescape_yz(is)=nescape_yz(is)+1
               !    else if (xpart > xmax) then
               !      nescape_zy(is)=nescape_zy(is)+1
               !    else if (ypart < 0.) then
               !      nescape_xz(is)=nescape_xz(is)+1
               !    else if (ypart > ymax) then
               !      nescape_zx(is)=nescape_zx(is)+1
               !    else if (zpart < 0.) then
               !      nescape_xy(is)=nescape_xy(is)+1
               !    else 
               !      nescape_yx(is)=nescape_yx(is)+1
               !    endif
 
               !    npleavingp=npleavingp+1             !VR: # of particles leaving domain?
               !    nescape(is) = nescape(is) + 1       !VR: # of particles escaping? (what's the difference?)
               !    iphead(ixe,iye,ize,is)=link(np)     !VR: we are removing this particle: iphead points to the next particle
               !    link(np)=ipstore                    !VR: the next particle is in the "empty" list
               !    ipstore=np                          !VR: empty list starts at this particle
               !    np=iphead(ixe,iye,ize,is)           !VR next particle to consider (remember that this is from the link)
 
               ! else  !VR: particle inside the global domain
                
 
                
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
                 ELSE
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
                 ENDIF ! VR: check particle inside/outside domain

!              ENDIF !VR: check particle inside/outside global domain
15            CONTINUE
              ENDDO !VR: loop over particles in the cell
              iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)    !VR: save the new head (re-ordered)
              iptemp(ixe,iye,ize,is)=0                         !VR: reset the temp list
            ENDDO !VR: grid loop in x
          ENDDO !VR: grid loop in y
        ENDDO!VR: grid loop in z
!
!************************************************************************
!
!     exchange data among processes and compute to see how many
!     particles each process has to send to, and receive from, other
!     processes
!
!     nsendp(nbrleft) -> how many particles to be sent from myid to 
!                         left neighbor
!     nsendp(nbrrite) -> how many particles to be sent from myid to 
!                         right neighbor
!     nrecvp(nbrleft) -> how many particles to be sent to myid from 
!                         left neighbor
!     nrecvp(nbrrite) -> how many particles to be sent to myid from 
!                         right neighbor
!
!
!  exchange information about particle numbers in two steps. First, 
!  send to right and receive from left. 

!  VR: is this comment still true? I can not see this
!  Note that the processors on 
!  the physical boundaries only send (myid == nbrleft) or receive 
!  (myid == nbrrite)
!
!

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
!
!***********************************************************************
!
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

        if ((myid == 0).and.prntinfo) then
          write(6,*) " FINISHED COMPILING LISTS "
          write(6,*) " # OF PARTICLES TO BE SENT     = ",NSENDTOT
          write(6,*) " # OF PARTICLES TO BE RECEIVED = ",NRECVTOT
          write(6,*) " # OF PARTICLES REMOVED BECAUSE V > VLIMIT = ",n_fast_removed
        endif
        if (NSENDTOT.NE.NRECVTOT) THEN
          CALL MPI_FINALIZE(IERR)
          WRITE(*,*)"Error: NSENDTOT != NRECVTOT. Terminating"
          STOP
        ENDIF
!
!       Check to see if each processor has enough particle storage
!       to handle incoming particles
!
        ! IF (NPTOTP+NRECVTOTP < NPLMAX) THEN
        !   EXIT_CODE_P = 0
        ! ELSE
        !   EXIT_CODE_P = 1
        !   write(6,*) " PROCESSOR # ",myid," RAN OUT OF PARTICLE STORAGE"
        ! ENDIF
        ! call MPI_ALLREDUCE(EXIT_CODE_P,EXIT_CODE,1,MPI_INTEGER8,MPI_SUM,&
        !                    MPI_COMM_WORLD,IERR)
        ! IF (EXIT_CODE /= 0) THEN
        !   CALL MPI_FINALIZE(IERR)
        !   STOP
        ! ENDIF
!
        !***********************************************************************
        !
        !     exchange particles with other processors 
        !     VR : old version. Here each particle is sent individually
        !
        ! nsendactualp=0
        ! nrecvactualp=0
        ! do irepeat=1,4
        !    if (isendid(irepeat) == 1) then
        !       NP=IPSEND(IS)
        !       DO WHILE (NP.NE.0)
        !          nsendactualp=nsendactualp+1

        !          !             Uniform mesh - Same as in version 5.0
        !          !              ixe=hxi*x(np)   +1.5000000000000001d+00
        !          !              iye=hyi*y(np)   +0.5000000000000001d+00
        !          !              ize=hzi*z(np)   +0.5000000000000001d+00

        !          !             Nonuniform mesh - using MESH_UNMAP
        !          rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
        !          rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
        !          rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
        !          ixe=rxe
        !          iye=rye
        !          ize=rze
        !          iye=iye-1             ! integer index in y direction starts at 0
        !          ize=ize-1             ! integer index in z direction starts at 0

        !          ypart=y(np)
        !          zpart=z(np)
        !          if (ypart <= ye.and.ypart >= yb) then
        !             iye_cc=jb 
        !          else
        !             if (ypart > ye) then
        !                iye_cc=je+1 
        !             else
        !                iye_cc=jb-1 
        !             endif
        !          endif
        !          if (zpart <= ze.and.zpart >= zb) then
        !             ize_cc=kb 
        !          else
        !             if (zpart > ze) then
        !                ize_cc=ke+1 
        !             else
        !                ize_cc=kb-1 
        !             endif
        !          endif
        !          pdata(1)=x(np)
        !          pdata(2)=y(np)
        !          pdata(3)=z(np)
        !          pdata(4)=vx(np)
        !          pdata(5)=vy(np)
        !          pdata(6)=vz(np)
        !          pdata(7)=qp(np)
        !          i_source = idmap_yz(iye_cc,ize_cc)
        !          i_tag    = it
        !          call MPI_SEND(pdata,7,MPI_DOUBLE_PRECISION,&
        !               i_source,i_tag,MPI_COMM_WORLD,IERR)
        !          ipsend(is)=link(np)
        !          link(np)=ipstore
        !          ipstore=np
        !          np=ipsend(is)
        !       ENDDO
        !    else
        !       nprecvtmp=0
        !       do itmp=1,4
        !          ipe=irecvid(itmp,irepeat)
        !          if (ipe.ne.-1) then
        !             nprecvtmp=nprecvtmp+nrecvp(irecvid(itmp,irepeat))
        !             nrecvp(irecvid(itmp,irepeat))=0
        !          endif
        !       enddo
        !       do ii=1,nprecvtmp
        !          nrecvactualp=nrecvactualp+1
        !          nprecv=ipstore

        !          i_tag=it
        !          call MPI_RECV(pdata,7,MPI_DOUBLE_PRECISION,&
        !               MPI_ANY_SOURCE,I_TAG,        &
        !               MPI_COMM_WORLD,STATUS2,IERR)


        !          if (ipstore == 0) then
        !             Storage_Error_p = 1
        !          else

        !             x (nprecv)=pdata(1)
        !             y (nprecv)=pdata(2)
        !             z (nprecv)=pdata(3)
        !             vx(nprecv)=pdata(4)
        !             vy(nprecv)=pdata(5)
        !             vz(nprecv)=pdata(6)
        !             qp(nprecv)=pdata(7)

        !             !VR: check y & z global domain boundaries and loop particles in case of periodicity
        !             if ( y(nprecv) < zero) then
        !                y(nprecv) = y(nprecv) + ymax
        !             endif
        !             if ( y(nprecv) > ymax) then
        !                y(nprecv) = y(nprecv) - ymax
        !             endif
        !             !VR: the same for z boundary
        !             if ( z(nprecv) < zero) then
        !                z(nprecv) = z(nprecv) + zmax
        !             endif
        !             if ( z(nprecv) > zmax) then
        !                z(nprecv) = z(nprecv) - zmax
        !             endif
        !             !VR: end of boundary check ----------------------------------




        !             !             Uniform mesh - Same as in version 5.0
        !             !              ixe=hxi*x(nprecv)+1.5000000000000001d+00
        !             !              iye=hyi*y(nprecv)+0.5000000000000001d+00
        !             !              ize=hzi*z(nprecv)+0.5000000000000001d+00

        !             !             Nonuniform mesh - using MESH_UNMAP
        !             rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000d+00
        !             rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000d+00
        !             rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000d+00
        !             ixe=rxe
        !             iye=rye
        !             ize=rze
        !             iye=iye-1             ! integer index in y direction starts at 0
        !             ize=ize-1             ! integer index in z direction starts at 0

        !             ipstore=link(nprecv)

        !             if ((ixe > nx+1  .or. ixe < 1 ) .or. (iye > je+1    .or. iye < jb-1) .or. (ize > ke+1    .or. ize < kb-1)) then
        !                Field_Diverge_p = 1
        !                ixe = min(max(iye,1   ),nx+1)
        !                iye = min(max(iye,jb-1),je+1)
        !                ize = min(max(ize,kb-1),ke+1)
        !             endif

        !             link(nprecv)=iphead(ixe,iye,ize,is)
        !             iphead(ixe,iye,ize,is)=nprecv

        !          endif

        !       enddo
        !    endif
        ! enddo

!***********************************************************************
!
!     exchange particles with other processors
!     VR: re-written on 04/19/2016 to use 
!     VR: packed exchage that reduces number of MPI requests
!     VR: thus increasing performance (a factor of 5 on my tests on 512 procs.) 
!     VR: stability
!
!***********************************************************************


        ! time diagnostic
        ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ! mp_elapsed = MPI_WTIME()


        !VR it's convenient to have neighboring processes in an array
        !Eventually, we should change all of the code to use this
        nbrs(1) = NBRTOP
        nbrs(2) = NBRLEFTTOP
        nbrs(3) = NBRLEFT
        nbrs(4) = NBRLEFTBOT
        nbrs(5) = NBRBOT
        nbrs(6) = NBRRITEBOT
        nbrs(7) = NBRRITE
        nbrs(8) = NBRRITETOP


        !VR maximum numbers of particles to be sent/received to/from a process
        max_nsend = maxval(nsendp)
        max_nrecv = maxval(nrecvp)

        !VR allocate tmp arrays for sending/recieving data
        allocate(packed_pdata_send(8,max_nsend,8))
        allocate(packed_pdata_recv(8,max_nrecv))

        !VR counters for particles sent/recieved
        nsendactualp = 0
        nrecvactualp = 0

        !VR 4 stages of data exhcnage in a simple even->odd, odd->even 2D pattern
        do irepeat=1,4

           if (isendid(irepeat) == 1) then
              !
              nsend_to_nbr = 0
              NP=IPSEND(IS)
              
              ! loop over particles in the ipsend list
              DO WHILE (NP.NE.0)
                 
                 nsendactualp = nsendactualp + 1

                 ! map this particle to the logical mesh
                 !             Nonuniform mesh - using MESH_UNMAP
                 rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
                 rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
                 rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
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
                 elseif (i_source==NBRLEFTTOP) then
                    idest = 2
                 elseif (i_source==NBRLEFT) then
                    idest = 3
                 elseif (i_source==NBRLEFTBOT) then
                    idest = 4
                 elseif (i_source==NBRBOT) then
                    idest = 5
                 elseif (i_source==NBRRITEBOT) then
                    idest = 6
                 elseif (i_source==NBRRITE) then
                    idest = 7
                 elseif (i_source==NBRRITETOP) then
                    idest = 8
                 else
                    print *,"myid"," Something is wrong:trying to send particles to procsees that are not immediate neighbors"
                 end if
                 
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
              ENDDO

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
                    
                    !             Uniform mesh - Same as in version 5.0
                    !              ixe=hxi*x(nprecv)+1.5000000000000001d+00
                    !              iye=hyi*y(nprecv)+0.5000000000000001d+00
                    !              ize=hzi*z(nprecv)+0.5000000000000001d+00
                    
                    !             Nonuniform mesh - using MESH_UNMAP
                    rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000d+00
                    rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000d+00
                    rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000d+00
                    ixe=rxe
                    iye=rye
                    ize=rze
                    iye=iye-1             ! integer index in y direction starts at 0
                    ize=ize-1             ! integer index in z direction starts at 0
                    !VR head into the list of "empty" particles
                    ipstore=link(nprecv)
                    
                    if ((ixe > nx+1  .or. ixe < 1 ) .or. (iye > je+1    .or. iye < jb-1) .or. (ize > ke+1    .or. ize < kb-1)) then
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

        
        ! VR: for diagnostic purposes: timings
        ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ! mp_elapsed = MPI_WTIME() - mp_elapsed

        ! exchange_time_total = exchange_time_total+mp_elapsed

        ! if (it ==100) then
        !    if (myid==0) print *, "Exchange time average over 100 it:", exchange_time_total/100.0
        ! endif
        
        ! VR: deallocate packed arrays
        deallocate(packed_pdata_send)
        deallocate(packed_pdata_recv)
        nullify(pp)
        
!      ------- end of packed exchange        

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
             write(6,*)" "
             write(6,*)" "
             write(6,*) "Particle storage allocation is exceeded."
             write(6,*) "3DHybrid is stopped"
             write(6,*)" "
             write(6,*)" "
           endif
           call MPI_FINALIZE(IERR)
           STOP
        endif

        if (Field_Diverge /= 0) then
           if (myid == 0) then
             write(6,*)" "
             write(6,*)" "
             write(6,*) "Field Solver Diverges"
             write(6,*) "3DHybrid is stopped"
             write(6,*)" "
             write(6,*)" "
           endif
           call MPI_FINALIZE(IERR)
           STOP
        endif

!
!=======================================================================
!
        if ((myid == 0).and.prntinfo) then
          write(6,*) " FINISHED EXCHANGING PARTICLES "
          write(6,*) " # OF PARTICLES       SENT     = ",NSENDACTUAL
          write(6,*) " # OF PARTICLES       RECEIVED = ",NRECVACTUAL
        endif
!
!=======================================================================
!
        ! NPTOTP=0
        ! do ize=kb-1,ke
        !   do iye=jb-1,je
        !     do ixe=1,nx1
        !           NP=IPHEAD(ixe,iye,ize,is)
        !           DO WHILE (NP.NE.0)
        !             NPTOTP=NPTOTP+1
        !             NP=LINK(NP)
        !           ENDDO
        !     enddo
        !   enddo
        ! enddo
        
        ! call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,&
        !                    MPI_COMM_WORLD,IERR)
        ! IF ((MYID.EQ.0).and.prntinfo) THEN
        !   WRITE(6,*) " IS = ",IS
        !   WRITE(6,*) " TOTAL # OF PARTICLES AFTER  PARMOV = ",NPTOT
        ! ENDIF
 999    CONTINUE

      ENDDO ! for IS
 
      end subroutine particle_boundary

