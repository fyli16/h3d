!---------------------------------------------------------------------
subroutine field_2d
  use parameter_mod
  use mesh_mod
  implicit none

  call date_and_time(values=time_begin_array(:,21))
  call date_and_time(values=time_begin_array(:,9))
  call pressgrad_2d(1)
  call date_and_time(values=time_end_array(:,9))
  call accumulate_time(time_begin_array(1,9),time_end_array(1,9),time_elapsed(9))

  call date_and_time(values=time_begin_array(:,10))
  call bcalc_2d
  call date_and_time(values=time_end_array(:,10))
  call accumulate_time(time_begin_array(1,10),time_end_array(1,10),time_elapsed(10))

  call date_and_time(values=time_begin_array(:,9))
  call pressgrad_2d(0)
  call date_and_time(values=time_end_array(:,9))
  call accumulate_time(time_begin_array(1,9),time_end_array(1,9),time_elapsed(9))

  call date_and_time(values=time_begin_array(:,11))
  call ecalc_2d( 0)
  call date_and_time(values=time_end_array(:,11))
  call accumulate_time(time_begin_array(1,11),time_end_array(1,11),time_elapsed(11))

  call date_and_time(values=time_begin_array(:,12))
  call focalc_2d
  call date_and_time(values=time_end_array(:,12))
  call accumulate_time(time_begin_array(1,12),time_end_array(1,12),time_elapsed(12))

  call date_and_time(values=time_end_array(:,21))
  call accumulate_time(time_begin_array(1,21),time_end_array(1,21),time_elapsed(21))

  return
end subroutine field_2d


!---------------------------------------------------------------------
subroutine pressgrad_2d(iflag)
  use parameter_mod
  use mesh_mod
  implicit none
  integer :: iflag
  integer*8 :: i,j,k
  real*8 :: dena,dxa,dya,dza,a

  do k=kb,ke
    do j = jb,je
      do i=2,nx1
        dena = iflag*0.5*(den(i,j,k)+deno(i,j,k)) + (1.-iflag)*den(i,j,k)
        a=1/dena

        ! Uniform mesh - Same as is in version 5.0
        ! dxa=a/(2.*hx)
        ! dya=a/(2.*hy)
        ! dza=a/(2.*hz)

        ! Nonuniform mesh
        dxa=a/((meshX%dxn(i  )+meshX%dxn(i+1)))
        dya=a/((meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
        dza=a/((meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0

        dpedx(i,j,k)=(  (pe(i+1,j-1,k  )+2.*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4.    &
                      - (pe(i-1,j-1,k  )+2.*pe(i-1,j,k  )+pe(i-1,j+1,k  ))/4. )*dxa
        dpedy(i,j,k)=(  (pe(i-1,j+1,k  )+2.*pe(i,j+1,k  )+pe(i+1,j+1,k  ))/4.    &
                      - (pe(i-1,j-1,k  )+2.*pe(i,j-1,k  )+pe(i+1,j-1,k  ))/4. )*dya
        dpedz(i,j,k)=0.0
      enddo
    enddo
  enddo
  return
end subroutine pressgrad_2d


!---------------------------------------------------------------------
subroutine ecalc_2d( iflag)
  use parameter_mod
  use mesh_mod
  implicit none

  integer :: iflag
  integer*8 :: i,j,k

  real*8 :: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &
                    ,curr_tot
  real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8 
  real*8 :: by1,by2,by3,by4,by5,by6,by7,by8 
  real*8 :: bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8  
  real*8 :: vixa, viya, viza, dena, a, dxa, dya, dza 
  real*8 :: dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy 
  real*8 :: curlbx_scalar,curlby_scalar,curlbz_scalar
  real*8 :: bxav, byav, bzav  
  real*8 :: dexdy, dexdz, deydx, deydz, dezdx,dezdy  

  do k=kb,ke
    do j = jb,je
      do i=2,nx1
        bx1=bx(i+1,j+1,k)  
        bx2=bx(i  ,j+1,k)  
        bx3=bx(i  ,j  ,k)  
        bx4=bx(i+1,j  ,k)  
        by1=by(i+1,j+1,k)  
        by2=by(i  ,j+1,k)  
        by3=by(i  ,j  ,k)  
        by4=by(i+1,j  ,k)  
        bz1=bz(i+1,j+1,k)  
        bz2=bz(i  ,j+1,k)  
        bz3=bz(i  ,j  ,k)  
        bz4=bz(i+1,j  ,k)  
        vixa=(1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k))&
            +iflag*vix(i,j,k)
        viya=(1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k))&
            +iflag*viy(i,j,k)
        viza=(1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k))&
            +iflag*viz(i,j,k)
        dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))&
            +(1.-iflag)*den(i,j,k)
        a=1/dena

        ! Uniform mesh - Same as is in version 5.0
        ! dxa=a/(2.*hx)
        ! dya=a/(2.*hy)
        ! dza=a/(2.*hz)

        ! Nonuniform mesh
        dxa=a/(2.*meshX%dxc(i))
        dya=a/(2.*meshY%dxc(j+1)) ! integer index in y direction starts at 0
        dza=a/(2.*meshZ%dxc(k+1)) ! integer index in z direction starts at 0

        dbxdy=+bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
              -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
        dbxdz= 0.0
        dbydx=+by(i+1,j  ,k  )+by(i+1,j+1,k  )&
              -by(i  ,j  ,k  )-by(i  ,j+1,k  )
        dbydz= 0.0
        dbzdx=+bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
              -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
        dbzdy=+bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
              -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
        curlbx_scalar=dya*dbzdy-dza*dbydz
        curlby_scalar=dza*dbxdz-dxa*dbzdx
        curlbz_scalar=dxa*dbydx-dya*dbxdy
        bxav=0.25*(bx1+bx2+bx3+bx4)
        byav=0.25*(by1+by2+by3+by4)
        bzav=0.25*(bz1+bz2+bz3+bz4)

        ! 6/25/2006 New eta_par option: tensor eta
        xj = curlbx_scalar
        yj = curlby_scalar
        zj = curlbz_scalar

        if (eta_par.eq.0) then
          tenx=eta(i,j,k)*xj
          teny=eta(i,j,k)*yj
          tenz=eta(i,j,k)*zj
        else
          bxx = bxav
          byy = byav
          bzz = bzav
          btot = sqrt(bxx**2 + byy**2 + bzz**2)
          xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
          xtmp2m = 1.d-12
          curr_tot = max(xtmp2m,xtmp1m)

          if (eta_par.eq.1) then
            tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
            tenx = tjdotb*bxx/btot
            teny = tjdotb*byy/btot
            tenz = tjdotb*bzz/btot
          else if (eta_par.eq.2) then
            xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
            xtmp2m = 1.d-12
            curr_tot = max(xtmp2m,xtmp1m)
            tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
            tenx = min(resis,tenx)
            tenx = tenx*xj
            teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
            teny = min(resis,teny)
            teny = teny*yj
            tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
            tenz = min(resis,tenz)
            tenz = tenz*zj
          endif
        endif

        ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)-dpedx(i,j,k)+tenx/a     
        ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)-dpedy(i,j,k)+teny/a     
        ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)-dpedz(i,j,k)+tenz/a                     
      enddo
    enddo
  enddo
 
  ! boundary conditions
  call date_and_time(values=time_begin_array(:,18))
  call XREALBCC_PACK_E_2D(EX,EY,EZ,1_8,NX,NY,NZ)
  call date_and_time(values=time_end_array(:,18))
  call accumulate_time(time_begin_array(1,18),time_end_array(1,18),time_elapsed(18))

  !VR: boundaries in x & z
  ex(nx2,:,:)=ex(2,:,:)
  ey(nx2,:,:)=ey(2,:,:)
  ez(nx2,:,:)=ez(2,:,:)

  ex(1,:,:)=ex(nx1,:,:)
  ey(1,:,:)=ey(nx1,:,:)
  ez(1,:,:)=ez(nx1,:,:)

  ex(:,:,kb-1)=ex(:,:,kb)
  ey(:,:,kb-1)=ey(:,:,kb)
  ez(:,:,kb-1)=ez(:,:,kb)

  ex(:,:,kb+1)=ex(:,:,kb)
  ey(:,:,kb+1)=ey(:,:,kb)
  ez(:,:,kb+1)=ez(:,:,kb)

  ! calculate curl E
  do k=kb,ke+1
    do j = jb,je+1
      do i=2,nx2
        dexdy= ex(i  ,j  ,k)+ex(i-1,j  ,k)&
              -ex(i  ,j-1,k)-ex(i-1,j-1,k)
        dexdz= 0.0
        deydx= ey(i  ,j  ,k)+ey(i  ,j-1,k)&
              -ey(i-1,j  ,k)-ey(i-1,j-1,k)
        deydz= 0.0
        dezdx= ez(i  ,j  ,k)+ez(i  ,j-1,k)&
              -ez(i-1,j  ,k)-ez(i-1,j-1,k)
        dezdy= ez(i  ,j  ,k)+ez(i-1,j  ,k)&
              -ez(i  ,j-1,k)-ez(i-1,j-1,k)

        ! Uniform mesh - Same as is in version 5.0
        ! curlex(i,j,k)=dezdy/(2.*hy)-deydz/(2.*hz)
        ! curley(i,j,k)=dexdz/(2.*hz)-dezdx/(2.*hx)
        ! curlez(i,j,k)=deydx/(2.*hx)-dexdy/(2.*hy)

        ! Nonuniform mesh
        curlex(i,j,k)=dezdy/(2.*meshY%dxn(j+1))-deydz/(2.*meshZ%dxn(k+1)) ! integer index in y and z directions start  at 0
        curley(i,j,k)=dexdz/(2.*meshZ%dxn(k+1))-dezdx/(2.*meshX%dxn(i  )) ! integer index in z       direction  starts at 0
        curlez(i,j,k)=deydx/(2.*meshX%dxn(i  ))-dexdy/(2.*meshY%dxn(j+1)) ! integer index in y       direction  starts at 0
      enddo
    enddo
  enddo
 
  return
end subroutine ecalc_2d


!---------------------------------------------------------------------
subroutine bcalc_2d
  use parameter_mod
  use mesh_mod
  implicit none
  integer*8 :: i, j, k, ii
  real*8 :: dts, dts2, dts6
  real*8 :: tempx1(nxmax,jb-1:je+1,kb-1:ke+1)&
                    ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1)&
                    ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)

  call date_and_time(values=time_begin_array(:,22))

  dts=dt/real(iterb)
  dts2=dts/2.
  dts6=dts/6.

  ! subcycle into iterb interations
  do ii = 1, iterb
    bxs=bx
    bys=by
    bzs=bz

    ! R-K first part
    call date_and_time(values=time_begin_array(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))
 
    ! B = B(n)+dt*K1/2
    bx=bxs-dts2*curlex
    by=bys-dts2*curley
    bz=bzs-dts2*curlez

    ! temp1 = K1
    tempx1=curlex
    tempy1=curley
    tempz1=curlez

    ! R-K part 2
    call date_and_time(values=time_begin_array(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))
 
    ! B = B(n)+dt*K2/2
    bx=bxs-dts2*curlex
    by=bys-dts2*curley
    bz=bzs-dts2*curlez

    ! temp2 = K2
    tempx1=tempx1+2.*curlex
    tempy1=tempy1+2.*curley
    tempz1=tempz1+2.*curlez

    ! R-K  part 3
    call date_and_time(values=time_begin_array(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))

    ! B = B(n)+dt*K3
    bx=bxs-dts*curlex
    by=bys-dts*curley
    bz=bzs-dts*curlez

    ! temp3 = K3
    tempx1=tempx1+2.*curlex
    tempy1=tempy1+2.*curley
    tempz1=tempz1+2.*curlez

    ! R-K  part 4
    call date_and_time(values=time_begin_array(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end_array(:,16)) 
    call accumulate_time(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))
 
    ! B = B(n) + dt*(K1+2K2+2K3+K4)/6
    bx=bxs-dts6*(tempx1+curlex)
    by=bys-dts6*(tempy1+curley)
    bz=bzs-dts6*(tempz1+curlez)

    ! restore boundary values 
    ! bx(1,:   ,:   )=bxs(1,:   ,:   )
    ! bx(:,jb-1,:   )=bxs(:,jb-1,:   )
    ! bx(:,:   ,kb-1)=bxs(:,:   ,kb-1)
    ! by(1,:   ,:   )=bys(1,:   ,:   )
    ! by(:,jb-1,:   )=bys(:,jb-1,:   )
    ! by(:,:   ,kb-1)=bys(:,:   ,kb-1)
    ! bz(1,:   ,:   )=bzs(1,:   ,:   )
    ! bz(:,jb-1,:   )=bzs(:,jb-1,:   )
    ! bz(:,:   ,kb-1)=bzs(:,:   ,kb-1)

    ! VR: why restore bxs values?? I comment this out for now
    ! bx(1,:   ,:   )=bxs(1,:   ,:   )
    ! if (jb == 1) bx(:,jb-1,:   )=bxs(:,jb-1,:   )
    ! if (kb == 1) bx(:,:   ,kb-1)=bxs(:,:   ,kb-1)
    ! by(1,:   ,:   )=bys(1,:   ,:   )
    ! if (jb == 1) by(:,jb-1,:   )=bys(:,jb-1,:   )
    ! if (kb == 1) by(:,:   ,kb-1)=bys(:,:   ,kb-1)
    ! bz(1,:   ,:   )=bzs(1,:   ,:   )
    ! if (jb == 1) bz(:,jb-1,:   )=bzs(:,jb-1,:   )
    ! if (kb == 1) bz(:,:   ,kb-1)=bzs(:,:   ,kb-1)

    ! VR: unlike 3D version, here B is exchanged at every step.
    ! VR: We need to check what is correct

    call XREALBCC_PACK_B_2D(BX,BY,BZ,1_8,NX,NY,NZ)
    !VR: exchange X and Z info
    bx(1  ,:,:)=bx(nx1  ,:,:)
    by(1  ,:,:)=by(nx1  ,:,:)
    bz(1  ,:,:)=bz(nx1  ,:,:)
    
    bx(nx2  ,:,:)=bx(2  ,:,:)
    by(nx2  ,:,:)=by(2  ,:,:)
    bz(nx2  ,:,:)=bz(2  ,:,:)

    bx(:,:,kb-1)=bx(:,:,kb)
    by(:,:,kb-1)=by(:,:,kb)
    bz(:,:,kb-1)=bz(:,:,kb)
    
    bx(:,:,kb+1)=bx(:,:,kb)
    by(:,:,kb+1)=by(:,:,kb)
    bz(:,:,kb+1)=bz(:,:,kb)
  end do

  ! VR not needed if it's included in the loop 
  ! call XREALBCC_PACK_B_2D(BX,BY,BZ,1_8,NX,NY,NZ)

  ! bx(1  ,:,:)=bx(nx1  ,j,k)
  ! by(1  ,:,:)=by(nx1  ,j,k)
  ! bz(1  ,:,:)=bz(nx1  ,j,k)
  
  ! bx(nx2  ,:,:)=bx(2  ,j,k)
  ! by(nx2  ,:,:)=by(2  ,j,k)
  ! bz(nx2  ,:,:)=bz(2  ,j,k)

  ! bx(:,:,kb-1)=bx(:,:,kb)
  ! by(:,:,kb-1)=by(:,:,kb)
  ! bz(:,:,kb-1)=bz(:,:,kb)
  
  ! bx(:,:,kb+1)=bx(:,:,kb)
  ! by(:,:,kb+1)=by(:,:,kb)
  ! bz(:,:,kb+1)=bz(:,:,kb)
 
  call date_and_time(values=time_end_array(:,22))
  call accumulate_time(time_begin_array(1,22),time_end_array(1,22),time_elapsed(22))

  return
end subroutine bcalc_2d


!---------------------------------------------------------------------
subroutine focalc_2d
  use parameter_mod
  use mesh_mod
  implicit none

  real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8
  real*8 :: by1,by2,by3,by4,by5,by6,by7,by8
  real*8 :: bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8
  real*8 :: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &
                    ,curr_tot
  integer*8 :: i,j,k
  real*8 :: dbxdy,dbydx,dbzdx,dbxdz,dbzdy,dbydz
  real*8 :: curlbx_scalar,curlby_scalar,curlbz_scalar,bxav,byav,bzav

  do k = kb, ke 
    do j = jb,je
      do i = 2, nx1
        dbxdy= bx(i+1,j+1,k) + bx(i,j+1,k) - bx(i+1,j,k) - bx(i,j,k)
        dbxdz= 0.0
        dbydx= by(i+1,j+1,k) + by(i+1,j,k) - by(i,j+1,k) - by(i,j,k)
        dbydz= 0.0
        dbzdx= bz(i+1,j+1,k) + bz(i+1,j,k) - bz(i,j+1,k) - bz(i,j,k)
        dbzdy= bz(i+1,j+1,k) + bz(i,j+1,k) - bz(i+1,j,k) - bz(i,j,k)

        ! Uniform mesh - Same as is in version 5.0
        ! curlbx_scalar=dbzdy/(2.*hy)-dbydz/(2.*hz)
        ! curlby_scalar=dbxdz/(2.*hz)-dbzdx/(2.*hx)
        ! curlbz_scalar=dbydx/(2.*hx)-dbxdy/(2.*hy)

        ! Nonuniform mesh
        curlbx_scalar=dbzdy/(2.*meshY%dxc(j+1))-dbydz/(2.*meshZ%dxc(k+1))
        curlby_scalar=dbxdz/(2.*meshZ%dxc(k+1))-dbzdx/(2.*meshX%dxc(i  ))
        curlbz_scalar=dbydx/(2.*meshX%dxc(i  ))-dbxdy/(2.*meshY%dxc(j+1))
 
        ! 6/25/2006 New eta_par option: tensor eta
        bx1=bx(i+1,j+1,k)  
        bx2=bx(i  ,j+1,k)  
        bx3=bx(i  ,j  ,k)  
        bx4=bx(i+1,j  ,k)  
        by1=by(i+1,j+1,k)  
        by2=by(i  ,j+1,k)  
        by3=by(i  ,j  ,k)  
        by4=by(i+1,j  ,k)  
        bz1=bz(i+1,j+1,k)  
        bz2=bz(i  ,j+1,k)  
        bz3=bz(i  ,j  ,k)  
        bz4=bz(i+1,j  ,k)  
        bxav=.25*(bx1+bx2+bx3+bx4)
        byav=.25*(by1+by2+by3+by4)
        bzav=.25*(bz1+bz2+bz3+bz4)
        xj = curlbx_scalar
        yj = curlby_scalar
        zj = curlbz_scalar
        bxx = bxav
        byy = byav
        bzz = bzav
        xtmp1m = sqrt(bxx**2 + byy**2 + bzz**2)
	      btot = max(1.d-12,xtmp1m)
        xtmp1m = sqrt(xj**2 + yj**2 + zj**2) 
        curr_tot = max(1.d-12,xtmp1m)
        
        if (eta_par.eq.0) then
          tenx=eta(i,j,k)*xj
          teny=eta(i,j,k)*yj
          tenz=eta(i,j,k)*zj
	      else
          bxx = bxav
          byy = byav
          bzz = bzav
          xtmp1m = sqrt(bxx**2 + byy**2 + bzz**2)
		      btot = max(1.d-12,xtmp1m)
          xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
	        curr_tot = max(1.d-12,xtmp1m)
		 
          if (eta_par.eq.1) then
            tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
            tenx = tjdotb*bxx/btot
            teny = tjdotb*byy/btot
            tenz = tjdotb*bzz/btot
          else if (eta_par.eq.2) then
            tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
            tenx = min(resis,tenx)
            tenx = tenx*xj
            teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
            teny = min(resis,teny)
            teny = teny*yj
            tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
            tenz = min(resis,tenz)
            tenz = tenz*zj
          endif
          eta_times_b_dot_j(i,j,k) = min(resis,eta(i,j,k)*abs((bxav*xj + byav*yj + bzav*zj))/(curr_tot*btot))
        endif
        
        fox(i,j,k)=-tenx
        foy(i,j,k)=-teny
        foz(i,j,k)=-tenz

      enddo
    enddo
  enddo

  ! boundary conditions
  ! first update internal ghost cells so that all processors
  ! have latest information. Care must be exercised so that
  ! BCs are set wrt proecessor that corresponds to that location. 
  ! To that end, keep z loops set to limits of kb and ke.
  call XREALBCC_PACK_E_2D(fox,foy,foz,1_8,NX,NY,NZ)

  fox(1  ,:,:)=fox(nx1  ,:,:)
  foy(1  ,:,:)=foy(nx1  ,:,:)
  foz(1  ,:,:)=foz(nx1  ,:,:)
  fox(nx2,:,:)=fox(2,:,:)
  foy(nx2,:,:)=foy(2,:,:)
  foz(nx2,:,:)=foz(2,:,:)

  fox(:,:,kb-1)=fox(:,:,kb)
  foy(:,:,kb-1)=foy(:,:,kb)
  foz(:,:,kb-1)=foz(:,:,kb)
  
  fox(:,:,kb+1)=fox(:,:,kb)
  foy(:,:,kb+1)=foy(:,:,kb)
  foz(:,:,kb+1)=foz(:,:,kb)

  return
end subroutine focalc_2d


!---------------------------------------------------------------------
subroutine parmov_2d
    use parameter_mod
    use mesh_mod
    implicit none

    real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8,by1,by2,by3,by4,by5,by6,by7,by8, &
              bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8,bxa,bya,bza
    real*8 :: ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ey1,ey2,ey3,ey4,ey5,ey6,ey7,ey8, &
              ez1,ez2,ez3,ez4,ez5,ez6,ez7,ez8,exa,eya,eza

    real*8 :: d_ranf,deltime1,deltime2,ff,h,hh
    real*8 :: fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8,foxa
    real*8 :: foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8,foya
    real*8 :: foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8,foza
    real*8 :: w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
    real*8 :: vex,vey,vez,vmag,vx_tmp,vy_tmp,vz_tmp
    real*8 :: p2xs,p2ys,p2zs,q_p,th
    real*8 :: wmult

    integer*8 :: i,ii,iix,iixe,iiy,iiye,iiz,iize,irepeat,irepeatp,is,itmp
    integer*8 :: iv,iye_cc,ize_cc,j,jv,k,kspc,npleavingp,nprecv,nprecvtmp
    integer*8 :: icount
    integer*8 :: count_kbq
    integer :: time_begin(8),time_end(8)
    integer*8 nptotp_kbq,npart_kbq(2),np_ijk,Storage_Error_p,Storage_Error
    data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
    data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
    data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
    integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
    real*8 :: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,xpart,ypart,zpart
    real*8 :: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
    real*8 :: v_limit,eps2,myranf,fluxran,vxa,vyz,vza
    INTEGER*8:: L, EXIT_CODE_P, EXIT_CODE
    integer*8:: n_fast_removed,n_fast_removed_local,Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p
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

 
    call date_and_time(values=time_begin_array(:,19))

    Storage_Error_p = 0
    Field_Diverge_p = 0
    Courant_Violation_p = 0
    x_disp_max_p        = 0
    y_disp_max_p        = 0
    z_disp_max_p        = 0

    dtxi = 1./meshX%dt
    dtyi = 1./meshY%dt
    dtzi = 1./meshZ%dt

    d_ranf=1./1001.
    eps2=1.d-25
 
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
 
    bx_av=0.;by_av=0.;bz_av=0.
    do K = KB-1,KE
      do J = JB-1,JE
        do I = 1, NX1
          bx_av(i,j,k)=0.25*( bx(i  ,j  ,k  )             &
                              +bx(i+1,j  ,k  )             &
                              +bx(i  ,j+1,k  )             &
                              +bx(i+1,j+1,k  )             &
                              )
          by_av(i,j,k)=0.25*( by(i  ,j  ,k  )             &
                              +by(i+1,j  ,k  )             &
                              +by(i  ,j+1,k  )             &
                              +by(i+1,j+1,k  )             &
                              )
          bz_av(i,j,k)=0.25*( bz(i  ,j  ,k  )             &
                              +bz(i+1,j  ,k  )             &
                              +bz(i  ,j+1,k  )             &
                              +bz(i+1,j+1,k  )             &
                              )
        enddo
      enddo
    enddo
    call XREALBCC_PACK_B(BX_AV,BY_AV,BZ_AV,1_8,NX,NY,NZ)
    call XREALBCC_PACK_B(BX   ,BY   ,BZ   ,1_8,NX,NY,NZ)

    if ((myid == 0).and.mod(it,n_print)==0) write(6,*) " Calling parmov, nspec = ", nspec

    ! initalize diagnostic variables that keep track of
    ! particle number, injection, and escape
    deltime1 = 0.0
    deltime2 = 0.0
    nptotp=0
    nptotp_kbq=0
    npleavingp=0
 
    ! beginning of main particle loop
    do IS = 1, NSPEC
      call date_and_time(values=time_begin_array(:,13))
      NPTOTP=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1   
              NP=IPHEAD(ixe,iye,ize,is)
              do while (NP.NE.0)
                NPTOTP=NPTOTP+1
                NP=LINK(NP)
              enddo
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
      if ((MYID.EQ.0).and.mod(it,n_print)==0) then
        write(6,*) " IS = ",IS
        write(6,*) " TOTAL # OF PARTICLES BEFORE PARMOV = ",NPTOT
      endif

      wmult=wspec(is)
      h=dt*qspec(is)/wmult
      hh=.5*h
 
      call date_and_time(values=time_end)
      clock_time1=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
      if (DT.NE.0) then
        npart(is) = 0
        npart_kbq(is) = 0
        do IIZE = KB-1,KE
          do IIYE = JB-1,JE
            do IIXE = 1, NX1
              NP=IPHEAD(IIXE,IIYE,IIZE,IS)

              ! begin advance of particle position and velocity
              ! If dt=0, skip
              do while (NP.NE.0)
                L=NP
                npart_kbq(is) = npart_kbq(is)+1
                nptotp_kbq = nptotp_kbq + 1

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

                ! Nonuniform mesh - without using MESH_UNMAP
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

                ! Nonuniform mesh - using MESH_UNMAP
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

                ! Diagnosic test on particle cell index algorithm
                ! if (     fxe < 0. .or. fxe > 1.                &
                !     .or. fye < 0. .or. fye > 1.                &
                !     .or. fze < 0. .or. fze > 1.) then 
                !     write(6,*) " SCATTER LOOP"
                !     write(6,*) fxe,fye,fze
                !     write(6,*) " x;",x(l),meshX%xn(ixe),meshX%xc(ixe)
                !     write(6,*) " y;",y(l),meshY%xn(iye),meshY%xc(iye)
                !     write(6,*) " z;",z(l),meshZ%xn(ize),meshZ%xc(ize)
                !     write(6,*) " r; ",rxe,rye,rze
                !     write(6,*) " ixc_map; ",ixe,iye,ize
                !     write(6,*) " xc     ; ",meshX%xc(ixe),meshY%xc(iye),meshZ%xc(ize)
                ! endif

                w1e=(1.-fxe)*(1.-fye)
                w2e=fxe*(1.-fye)
                w3e=(1.-fxe)*fye
                w4e=fxe*fye*(1.-fze)
 
                ex1=ex(ixe  ,iye  ,ize  )
                ex2=ex(ixep1,iye  ,ize  )
                ex3=ex(ixe  ,iyep1,ize  )
                ex4=ex(ixep1,iyep1,ize  )
                ey1=ey(ixe  ,iye  ,ize  )
                ey2=ey(ixep1,iye  ,ize  )
                ey3=ey(ixe  ,iyep1,ize  )
                ey4=ey(ixep1,iyep1,ize  )
                ez1=ez(ixe  ,iye  ,ize  )
                ez2=ez(ixep1,iye  ,ize  )
                ez3=ez(ixe  ,iyep1,ize  )
                ez4=ez(ixep1,iyep1,ize  )
 
                ! Old code that causes problems with density holes
                ! bx1=bx(ixe  ,iye  ,ize  )+bdipole_x(ixe  ,iye  ,ize  )
                ! bx2=bx(ixep1,iye  ,ize  )+bdipole_x(ixep1,iye  ,ize  )
                ! bx3=bx(ixe  ,iyep1,ize  )+bdipole_x(ixe  ,iyep1,ize  )
                ! bx4=bx(ixep1,iyep1,ize  )+bdipole_x(ixep1,iyep1,ize  )
                ! by1=by(ixe  ,iye  ,ize  )+bdipole_y(ixe  ,iye  ,ize  )
                ! by2=by(ixep1,iye  ,ize  )+bdipole_y(ixep1,iye  ,ize  )
                ! by3=by(ixe  ,iyep1,ize  )+bdipole_y(ixe  ,iyep1,ize  )
                ! by4=by(ixep1,iyep1,ize  )+bdipole_y(ixep1,iyep1,ize  )
                ! bz1=bz(ixe  ,iye  ,ize  )+bdipole_z(ixe  ,iye  ,ize  )
                ! bz2=bz(ixep1,iye  ,ize  )+bdipole_z(ixep1,iye  ,ize  )
                ! bz3=bz(ixe  ,iyep1,ize  )+bdipole_z(ixe  ,iyep1,ize  )
                ! bz4=bz(ixep1,iyep1,ize  )+bdipole_z(ixep1,iyep1,ize  )

                ! New code that fixes problems with density holes - consistent with 3D parmov
                bx1=bx_av(ixe  ,iye  ,ize  )
                bx2=bx_av(ixep1,iye  ,ize  )
                bx3=bx_av(ixe  ,iyep1,ize  )
                bx4=bx_av(ixep1,iyep1,ize  )
                by1=by_av(ixe  ,iye  ,ize  )
                by2=by_av(ixep1,iye  ,ize  )
                by3=by_av(ixe  ,iyep1,ize  )
                by4=by_av(ixep1,iyep1,ize  )
                bz1=bz_av(ixe  ,iye  ,ize  )
                bz2=bz_av(ixep1,iye  ,ize  )
                bz3=bz_av(ixe  ,iyep1,ize  )
                bz4=bz_av(ixep1,iyep1,ize  )

                fox1=fox(ixe  ,iye  ,ize  )
                fox2=fox(ixep1,iye  ,ize  )
                fox3=fox(ixe  ,iyep1,ize  )
                fox4=fox(ixep1,iyep1,ize  )
                foy1=foy(ixe  ,iye  ,ize  )
                foy2=foy(ixep1,iye  ,ize  )
                foy3=foy(ixe  ,iyep1,ize  )
                foy4=foy(ixep1,iyep1,ize  )
                foz1=foz(ixe  ,iye  ,ize  )
                foz2=foz(ixep1,iye  ,ize  )
                foz3=foz(ixe  ,iyep1,ize  )
                foz4=foz(ixep1,iyep1,ize  )
              
                exa=w1e*ex1+w2e*ex2+w3e*ex3+w4e*ex4      &
                    +w1e*fox1+w2e*fox2+w3e*fox3+w4e*fox4  
                eya=w1e*ey1+w2e*ey2+w3e*ey3+w4e*ey4      &
                    +w1e*foy1+w2e*foy2+w3e*foy3+w4e*foy4  
                eza=w1e*ez1+w2e*ez2+w3e*ez3+w4e*ez4      &
                    +w1e*foz1+w2e*foz2+w3e*foz3+w4e*foz4  
                
                bxa=w1e*bx1+w2e*bx2+w3e*bx3+w4e*bx4      
                bya=w1e*by1+w2e*by2+w3e*by3+w4e*by4      
                bza=w1e*bz1+w2e*bz2+w3e*bz3+w4e*bz4      
              
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
              
                x_disp = dt*vx(l)
                y_disp = dt*vy(l)
                
                x(l)=x(l)+ x_disp
                y(l)=y(l)+ y_disp
                  
                if ( abs(x_disp/meshX%dxn(ixep1)) > 1.0 .or.                                               &
                    abs(y_disp/meshY%dxn(iyep1)) > 1.0) Courant_Violation_p = Courant_Violation_p + 1
                x_disp_max_p = max(x_disp_max_p,abs(x_disp)/meshX%dxn(ixep1))
                y_disp_max_p = max(y_disp_max_p,abs(y_disp)/meshY%dxn(iyep1))

                ! periodic B.C. in x
                if (x(l) < zero) x(l) = x(l) + xmax
                if (x(l) > xmax) x(l) = x(l) - xmax
              
                NP=LINK(NP)
              enddo
            enddo
          enddo
        enddo
      endif

      call MPI_ALLREDUCE(x_disp_max_p,x_disp_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(y_disp_max_p,y_disp_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      if ((myid == 0).and.mod(it,n_print)==0) then
        write(6,*) " maximum x-displacement/dx = ",x_disp_max
        write(6,*) " maximum y-displacement/dy = ",y_disp_max
      endif

      call MPI_ALLREDUCE(Courant_Violation_p,Courant_Violation,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

      if (Courant_Violation /= 0) then
          if (myid == 0) write(6,*) "Particle displacements exceed cell size ",Courant_Violation," times"
          call MPI_FINALIZE(IERR)
          STOP
      endif

      call date_and_time(values=time_end_array(:,13))
      call accumulate_time(time_begin_array(1,13),time_end_array(1,13),time_elapsed(13))

      call date_and_time(values=time_begin_array(:,14))
      ICOUNT=0                  !ICOUNT records how many times the particle
11    ICOUNT=ICOUNT+1           !  exchange routine has been run. If there               
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
            NP=IPHEAD(IXE,IYE,IZE,IS)
            do while (NP.NE.0)
              xpart=x(np)
              ypart=y(np)
              zpart=z(np)
              l=np

              ! Remove fast particles
              if (abs(vx(l)) >= v_limit .or. abs(vy(l)) >= v_limit .or. abs(vz(l)) >= v_limit) then
                n_fast_removed_local = n_fast_removed_local + 1
                iphead(ixe,iye,ize,is)=link(np)
                link(np)=ipstore
                ipstore=np
                np=iphead(ixe,iye,ize,is)
                goto 15
              endif
                
              xpart = x(l)
              ypart = y(l)
              zpart = z(l)
              
              if (ypart <= ye.and.ypart >= yb) then
                  !VR: do nothing with particles, re-order the list
                  iphead(ixe,iye,ize,is)=link(np)
                  link(np)=iptemp(ixe,iye,ize,is)
                  iptemp(ixe,iye,ize,is)=np
                  np=iphead(ixe,iye,ize,is)
              else
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

                nsendp(idmap_yz(iye_cc,ize_cc))=nsendp(idmap_yz(iye_cc,ize_cc))+1
                iphead(ixe,iye,ize,is)=link(np)
                link(np)=ipsend(is)
                ipsend(is)=np
                np=iphead(ixe,iye,ize,is)
              endif
      
15            continue
            enddo
            iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
            iptemp(ixe,iye,ize,is)=0
          enddo
        enddo
      enddo

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
      ! send to right and receive from left. Note that the processors on 
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
      call MPI_ALLREDUCE(nsendtotp,nsendtot,1,MPI_INTEGER8,MPI_SUM,&
                          MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(nrecvtotp,nrecvtot,1,MPI_INTEGER8,MPI_SUM,&
                          MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(n_fast_removed_local,n_fast_removed,1,MPI_INTEGER8,MPI_SUM,&
                          MPI_COMM_WORLD,IERR)
      if ((myid == 0).and.mod(it,n_print)==0) then
        write(6,*) " FINISHED COMPILING LISTS "
        write(6,*) " # OF PARTICLES TO BE SENT     = ",NSENDTOT
        write(6,*) " # OF PARTICLES TO BE RECEIVED = ",NRECVTOT
        write(6,*) " # OF PARTICLES REMOVED BECAUSE V > VLIMIT = ",n_fast_removed
      endif

      if (NSENDTOT.NE.NRECVTOT) then
        call MPI_FINALIZE(IERR)
        write(*,*)"HERE TESTSSS"
        STOP
      endif

      ! Check to see if each processor has enough particle storage
      ! to handle incoming particles

      IF (NPTOTP+NRECVTOTP < NPLMAX) THEN
        EXIT_CODE_P = 0
      ELSE
        EXIT_CODE_P = 1
        write(6,*) " PROCESSOR # ",myid," RAN OUT OF PARTICLE STORAGE"
      endif
      call MPI_ALLREDUCE(EXIT_CODE_P,EXIT_CODE,1,MPI_INTEGER8,MPI_SUM,&
                          MPI_COMM_WORLD,IERR)
      IF (EXIT_CODE /= 0) THEN
        call MPI_FINALIZE(IERR)
        STOP
      endif

      ! exchange particles with other processors
      nsendactualp=0
      nrecvactualp=0
      do irepeat=1,4
        if (isendid(irepeat) == 1) then
          NP=IPSEND(IS)
          do while (NP.NE.0)
              nsendactualp=nsendactualp+1
              
            ! Uniform mesh - Same as in version 5.0
            ! ixe=hxi*x(np)   +1.5000000000000001d+00
            ! iye=hyi*y(np)   +0.5000000000000001d+00
            ! ize=hzi*z(np)   +0.5000000000000001d+00

            ! Nonuniform mesh - using MESH_UNMAP
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
            pdata(1)=x(np)
            pdata(2)=y(np)
            pdata(3)=z(np)
            pdata(4)=vx(np)
            pdata(5)=vy(np)
            pdata(6)=vz(np)
            pdata(7)=qp(np)
            i_source = idmap_yz(iye_cc,ize_cc)
            i_tag    = it
            call MPI_SEND(pdata,7,MPI_DOUBLE_PRECISION,&
                          i_source,i_tag,MPI_COMM_WORLD,IERR)
            ipsend(is)=link(np)
            link(np)=ipstore
            ipstore=np
            np=ipsend(is)
          enddo
        else
          nprecvtmp=0
          do itmp=1,4
            ipe=irecvid(itmp,irepeat)
            if (ipe.ne.-1) then
              nprecvtmp=nprecvtmp+nrecvp(irecvid(itmp,irepeat))
              nrecvp(irecvid(itmp,irepeat))=0
            endif
          enddo

          do ii=1,nprecvtmp
            nrecvactualp=nrecvactualp+1
            nprecv=ipstore

            i_tag=it
            call MPI_RECV(pdata,7,MPI_DOUBLE_PRECISION,&
                          MPI_ANY_SOURCE,I_TAG,        &
                          MPI_COMM_WORLD,STATUS2,IERR)
                            
            if (ipstore == 0) then
              Storage_Error_p = 1
            else
              x (nprecv)=pdata(1)
              y (nprecv)=pdata(2)
              z (nprecv)=pdata(3)
              vx(nprecv)=pdata(4)
              vy(nprecv)=pdata(5)
              vz(nprecv)=pdata(6)
              qp(nprecv)=pdata(7)

              !VR: impose periodic B.C.
              if (y (nprecv) < zero) y (nprecv) = y (nprecv) + ymax
              if (y (nprecv) > ymax) y (nprecv) = y (nprecv) - ymax

              ! Uniform mesh - Same as in version 5.0
              ! ixe=hxi*x(nprecv)+1.5000000000000001d+00
              ! iye=hyi*y(nprecv)+0.5000000000000001d+00
              ! ize=hzi*z(nprecv)+0.5000000000000001d+00

              ! Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000d+00
              rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000d+00
              rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000d+00
              ixe=rxe
              iye=rye
              ize=rze
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0

              ipstore=link(nprecv)

              if ((ixe > nx+1  .or. ixe < 1 ) .or. (iye > je+1    .or. iye < jb-1) .or. (ize > ke+1    .or. ize < kb-1)) then
                Field_Diverge_p = 1
                ixe = min(max(iye,1_8 ),nx+1)
                iye = min(max(iye,jb-1),je+1)
                ize = min(max(ize,kb-1),ke+1)
              endif
              link(nprecv)=iphead(ixe,iye,ize,is)
              iphead(ixe,iye,ize,is)=nprecv
            endif
          enddo
        endif
      enddo

      call MPI_ALLREDUCE(Storage_Error_p,Storage_Error,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
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

      call MPI_ALLREDUCE(Field_Diverge_p,Field_Diverge,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
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

      call MPI_ALLREDUCE(nsendactualp,nsendactual,1,MPI_INTEGER8,&
                          MPI_SUM,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(nrecvactualp,nrecvactual,1,MPI_INTEGER8,&
                          MPI_SUM,MPI_COMM_WORLD,IERR)

      if ((myid == 0).and.mod(it,n_print)==0) then
        write(6,*) " FINISHED EXCHANGING PARTICLES "
        write(6,*) " # OF PARTICLES       SENT     = ",NSENDACTUAL
        write(6,*) " # OF PARTICLES       RECEIVED = ",NRECVACTUAL
      endif

      NPTOTP=0
      do ize=kb-1,ke
        do iye=jb-1,je
          do ixe=1,nx1
            NP=IPHEAD(ixe,iye,ize,is)
            do while (NP.NE.0)
              NPTOTP=NPTOTP+1
              NP=LINK(NP)
            enddo
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,&
                          MPI_COMM_WORLD,IERR)
      IF ((myid.eq.0).and.mod(it,n_print)==0) THEN
        write(6,*) " IS = ",IS
        write(6,*) " TOTAL # OF PARTICLES AFTER  PARMOV = ",NPTOT
      endif
999   continue
 
      call date_and_time(values=time_end_array(:,14))
      call accumulate_time(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))

      call date_and_time(values=time_begin_array(:,15))
      NPTOTP = 0
      do IIZ=KB-1,KE
        do IIY=JB-1,JE
          do IIX=1,NX1
            np_ijk=0
            np=iphead(iix,iiy,iiz,is)
            do while (np.ne.0)
              np_ijk=np_ijk+1
              nptotp=nptotp+1           !count particles
              npart(is) = npart(is) + 1 !count particles in each species
              L=NP

              q_p = qp(l)

              ! Uniform mesh - Same as in version 5.0
              ! rx=hxi*x(l)+1.5000000000000001d+00
              ! ry=hyi*y(l)+0.5000000000000001d+00
              ! rz=hzi*z(l)+0.5000000000000001d+00
              ! ix=rx
              ! iy=ry
              ! iz=rz
              ! ix=max(1,min(ix,nx1))
              ! iy=max(jb-1,min(iy,je))
              ! iz=max(kb-1,min(iz,ke))
              ! fx=rx-ix
              ! fy=ry-iy
              ! fz=rz-iz

              ! Nonuniform mesh - without using MESH_UNMAP
              ! rx=hxi*x(l)+1.500000000000000d+00
              ! ry=hyi*y(l)+1.500000000000000d+00
              ! rz=hzi*z(l)+1.500000000000000d+00
              ! ix=rx
              ! iy=ry
              ! iz=rz
              ! ix=ixc_2_c_map(ix)
              ! iy=iyc_2_c_map(iy)
              ! iz=izc_2_c_map(iz)
              ! fx=(x(l)-meshX%xc(ix))/meshX%dxn(ix+1)
              ! fy=(y(l)-meshY%xc(iy))/meshY%dxn(iy+1)
              ! fz=(z(l)-meshZ%xc(iz))/meshZ%dxn(iz+1)

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

              ! Diagnosic test on particle cell index algorithm
              ! if (     fx < 0. .or. fx > 1.                &
              !     .or. fy < 0. .or. fy > 1.                &
              !     .or. fz < 0. .or. fz > 1.) then 
              !     write(6,*) " GATHER LOOP "
              !     write(6,*) fx,fy,fz
              !     write(6,*) " x;",x(l),meshX%xn(ix),meshX%xc(ix)
              !     write(6,*) " y;",y(l),meshY%xn(iy),meshY%xc(iy)
              !     write(6,*) " z;",z(l),meshZ%xn(iz),meshZ%xc(iz)
              !     write(6,*) " r; ",rx,ry,rz
              !     write(6,*) " i; ",ixe,iye,ize
              !     write(6,*) " ixc_map; ",ix,iy,iz
              ! endif

              w1=q_p*(1.-fx)*(1.-fy)
              w2=q_p*fx     *(1.-fy)
              w3=q_p*(1.-fx)*fy     
              w4=q_p*fx     *fy     
 
              dns(ix  ,iy  ,iz  ,is)=dns(ix  ,iy  ,iz  ,is)+w1
              dns(ixp1,iy  ,iz  ,is)=dns(ixp1,iy  ,iz  ,is)+w2
              dns(ix  ,iyp1,iz  ,is)=dns(ix  ,iyp1,iz  ,is)+w3
              dns(ixp1,iyp1,iz  ,is)=dns(ixp1,iyp1,iz  ,is)+w4

              vxs(ix  ,iy  ,iz  ,is)=vxs(ix  ,iy  ,iz  ,is)+w1*vx(l)
              vxs(ixp1,iy  ,iz  ,is)=vxs(ixp1,iy  ,iz  ,is)+w2*vx(l)
              vxs(ix  ,iyp1,iz  ,is)=vxs(ix  ,iyp1,iz  ,is)+w3*vx(l)
              vxs(ixp1,iyp1,iz  ,is)=vxs(ixp1,iyp1,iz  ,is)+w4*vx(l)

              vys(ix  ,iy  ,iz  ,is)=vys(ix  ,iy  ,iz  ,is)+w1*vy(l)
              vys(ixp1,iy  ,iz  ,is)=vys(ixp1,iy  ,iz  ,is)+w2*vy(l)
              vys(ix  ,iyp1,iz  ,is)=vys(ix  ,iyp1,iz  ,is)+w3*vy(l)
              vys(ixp1,iyp1,iz  ,is)=vys(ixp1,iyp1,iz  ,is)+w4*vy(l)

              vzs(ix  ,iy  ,iz  ,is)=vzs(ix  ,iy  ,iz  ,is)+w1*vz(l)
              vzs(ixp1,iy  ,iz  ,is)=vzs(ixp1,iy  ,iz  ,is)+w2*vz(l)
              vzs(ix  ,iyp1,iz  ,is)=vzs(ix  ,iyp1,iz  ,is)+w3*vz(l)
              vzs(ixp1,iyp1,iz  ,is)=vzs(ixp1,iyp1,iz  ,is)+w4*vz(l)

              np=link(np)
            enddo
          enddo
        enddo
      enddo

      kspc=is

      call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape(is),nescape_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape_yz(is),nescape_yz_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape_zy(is),nescape_zy_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape_xy(is),nescape_xy_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape_yx(is),nescape_yx_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape_zx(is),nescape_zx_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)
      call MPI_ALLREDUCE(nescape_xz(is),nescape_xz_global(is),1,MPI_INTEGER8,&
                        MPI_SUM,COMM2D,IERR)

      deltime2 = deltime2 + real(clock_time1-clock_time)

      call xreal_2d(DNS(1,jb-1,kb-1,is),NX,NY,NZ)
      call xreal_2d(VXS(1,jb-1,kb-1,is),NX,NY,NZ)
      call xreal_2d(VYS(1,jb-1,kb-1,is),NX,NY,NZ)
      call xreal_2d(VZS(1,jb-1,kb-1,is),NX,NY,NZ)
      call xrealbcc_2d(DNS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)
      call xrealbcc_2d(VXS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)
      call xrealbcc_2d(VYS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)
      call xrealbcc_2d(VZS(1,jb-1,kb-1,is),1_8,NX,NY,NZ)

      do IIZ=KB-1,KE+1
        do IIY=JB-1,JE+1
          do IIX=1,NX2
            dns(iix,iiy,iiz,is) = dns(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            vxs(iix,iiy,iiz,is) = vxs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            vys(iix,iiy,iiz,is) = vys(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            vzs(iix,iiy,iiz,is) = vzs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          enddo
        enddo
      enddo

      call date_and_time(values=time_end_array(:,15))
      call accumulate_time(time_begin_array(1,15),time_end_array(1,15),time_elapsed(15))
 
    enddo  ! IS DO LOOP

    ! end of main particle loop
    call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(npleavingp,npleaving,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

    if (mod(it,n_print)==0) then
      call date_and_time(values=time_end)
      clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
    endif

    if (myid == 0.and.mod(it,n_print)==0) then
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

1005  format(5x,'sum',4x,i8,2x,i8,2x,i10)
1006  format(2i6,4x,'sum',4x,i8,2x,i8,2x,i10)
    endif

    ninj        = 0
    ninj_global = 0

    call date_and_time(values=time_end_array(:,19))
    call accumulate_time(time_begin_array(1,19),time_end_array(1,19),time_elapsed(19))
 
    return
end subroutine parmov_2d


!---------------------------------------------------------------------
subroutine xreal_2d(a,nx1m,ny1m,nz1m)
  use parameter_mod
  implicit none

  integer*8 :: i, j, nx1m, ny1m, nz1m, k
  real*8 :: a(nxmax, jb-1:je+1, kb-1:ke+1)&
                    ,tmp(nxmax, jb-1:je+1, kb-1:ke+1)

  a(2,:,:) = a(2,:,:) + a(nx2,:,:)
  a(nx1,:,:) = a(nx1,:,:) + a(1,:,:)
  a(:,:,kb) = a(:,:,kb) + a(:,:,kb-1) + a(:,:,kb+1)
  call MPI_SENDRECV(a    (1    ,je+1,kb-1),1,stridery,nbrrite,0, &
                    tmp  (1    ,jb-1,kb-1),1,stridery,nbrleft,0, &
                    mpi_comm_world,status,ierr)
  a(:,jb,:)=a(:,jb,:)+tmp  (:,jb-1,:)
  call MPI_SENDRECV(a    (1    ,jb-1,kb-1),1,stridery,nbrleft,1, &
                    tmp  (1    ,je+1,kb-1),1,stridery,nbrrite,1, &
                    mpi_comm_world,status,ierr)
  a(:,je,:)=a(:,je,:)+tmp  (:,je+1,:)

  a(1   ,:,:   )=a(nx1   ,:,: )
  a(nx2,:,:   ) =a(2,:,: )
  a(:   ,:,kb-1)=a(:   ,:,kb)
  a(:   ,:,kb+1)=a(:   ,:,kb)

  return
end subroutine xreal_2d


!---------------------------------------------------------------------
subroutine xrealbcc_2d(a, ibnd, nx1m, ny1m,nz1m)
  use parameter_mod
  implicit none
  integer*8 :: ibnd,i,j,nx1m,ny1m,nz1m
  real*8 :: a(nxmax,jb-1:je+1,kb-1:ke+1)&
                    ,tmp(nxmax,jb-1:je+1,kb-1:ke+1)

  tmp = a
  a(1   ,:,:   )=a(nx1   ,:,: )
  a(nx2,:,:   ) =a(2,:,: )
  a(:   ,:,kb-1)=a(:   ,:,kb)
  a(:   ,:,kb+1)=a(:   ,:,kb)

  call MPI_SENDRECV(a(1    ,je  ,kb-1),1,stridery,nbrrite,0,&
                  a(1    ,jb-1,kb-1),1,stridery,nbrleft,0,&
                  mpi_comm_world,status,ierr)
  call MPI_SENDRECV(a(1    ,jb  ,kb-1),1,stridery,nbrleft,1,&
                  a(1    ,je+1,kb-1),1,stridery,nbrrite,1,&
                  mpi_comm_world,status,ierr)
  return
end subroutine xrealbcc_2d


!---------------------------------------------------------------------
subroutine xrealbcc_pack_b_2d(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)
  use parameter_mod
  implicit none

  integer*8 :: i,j,nx1m,ny1m,nz1m,k,ibnd
  real*8 :: a_x(nxmax,jb-1:je+1,kb-1:ke+1)&
                    ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                    ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                    ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                    ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                    ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                    ,packed_data_xy_recv(nxmax,jb-1:je+1,3)

  a_x(:,:,kb-1)=a_x(:,:,kb)
  a_y(:,:,kb-1)=a_y(:,:,kb)
  a_z(:,:,kb-1)=a_z(:,:,kb)
  a_x(:,:,ke+1)=a_x(:,:,ke)
  a_y(:,:,ke+1)=a_y(:,:,ke)
  a_z(:,:,ke+1)=a_z(:,:,ke)

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
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

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
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

  return
end subroutine xrealbcc_pack_b_2d

    
!---------------------------------------------------------------------
subroutine caltemp2_global_2d
  use parameter_mod
  use mesh_mod
  implicit none

  real*8 :: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
  integer*8 ix,iy,iz,ixp1,iyp1,izp1,iiy,iiye,iiz,iize,is,l,iix,iixe
  real*8 :: vxa,vya,vza,rfrac,vxavg,vxavg1,vxavg2 &
        ,vyavg,vyavg1,vyavg2,vzavg,vzavg1,vzavg2,wperp2,wpar,wmult
  real*8 :: w1,w2,w3,w4,w5,w6,w7,w8,h,hh,dns1,dns2,bxa,bya,bza,btota,dnst

  call date_and_time(values=time_begin_array(:,23))
 
  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt

  tpar  = 0.
  tperp = 0.
  rfrac = 0.

  p_xx=0.;p_xy=0.;p_xz=0.;p_yy=0.;p_yz=0.;p_zz=0.

  if (nspec >= 2) rfrac = frac(2)/frac(1)

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
          do while (NP.NE.0)
            L=NP

            ! Uniform mesh - Same as is in version 5.0
            ! rx=hxi*x(l)+1.5000000000000001
            ! ry=hyi*y(l)+0.5000000000000001d+00
            ! rz=hzi*z(l)+0.5000000000000001d+00
            ! ix=rx
            ! iy=ry
            ! iz=rz
            ! IZ=1
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
            IZ=1
            fx=rx-ix
            fy=ry-iy
            fz=rz-iz
            iy=iy-1             ! integer index in y direction starts at 0
            iz=iz-1             ! integer index in z direction starts at 0

            ixp1 = ix + 1
            iyp1 = iy + 1

            w1=(1.-fx)*(1.-fy)
            w2=    fx *(1.-fy)
            w3=(1.-fx)*    fy
            w4=    fx*     fy

            dns1= dns(ix,iy  ,iz  ,1)*w1+dns(ix+1,iy  ,iz  ,1)*w2&
            +     dns(ix,iy+1,iz  ,1)*w3+dns(ix+1,iy+1,iz  ,1)*w4

            dns2= 0.

            dnst = dns1 + dns2
    
            vxavg1=vxs(ix,iy  ,iz  ,1)*w1+vxs(ix+1,iy  ,iz  ,1)*w2&
            +      vxs(ix,iy+1,iz  ,1)*w3+vxs(ix+1,iy+1,iz  ,1)*w4

            vxavg2= 0.

            vxavg = (dns1*vxavg1 + dns2*vxavg2)/dnst

            vyavg1=vys(ix,iy  ,iz  ,1)*w1+vys(ix+1,iy  ,iz  ,1)*w2&
            +      vys(ix,iy+1,iz  ,1)*w3+vys(ix+1,iy+1,iz  ,1)*w4

            vyavg2=0.

            vyavg = (dns1*vyavg1 + dns2*vyavg2)/dnst

            vzavg1=vzs(ix,iy  ,iz  ,1)*w1+vzs(ix+1,iy  ,iz  ,1)*w2&
            +      vzs(ix,iy+1,iz  ,1)*w3+vzs(ix+1,iy+1,iz  ,1)*w4

            vzavg2=0.

            vzavg = (dns1*vzavg1 + dns2*vzavg2)/dnst

            vxa=vx(l)-vxavg
            vya=vy(l)-vyavg
            vza=vz(l)-vzavg

            bxa  =bx       (ix,iy  ,iz  )*w1+bx       (ix+1,iy  ,iz  )*w2&
            +     bx       (ix,iy+1,iz  )*w3+bx(       ix+1,iy+1,iz  )*w4

            bya  =by(       ix,iy  ,iz  )*w1+by(       ix+1,iy  ,iz  )*w2&
            +     by(       ix,iy+1,iz  )*w3+by       (ix+1,iy+1,iz  )*w4

            bza  =bz       (ix,iy  ,iz  )*w1+bz       (ix+1,iy  ,iz  )*w2&
            +     bz       (ix,iy+1,iz  )*w3+bz       (ix+1,iy+1,iz  )*w4

            btota=sqrt(bxa**2+bya**2+bza**2)
            if(btota.lt.1.e-20) btota=1.e-20
            wpar=(vxa*bxa+vya*bya+vza*bza)/btota
            wperp2=vxa**2+vya**2+vza**2-wpar**2
            xx=vxa**2
            xy=vxa*vya
            xz=vxa*vza
            yy=vya**2
            yz=vya*vza
            zz=vza**2

            tpar (ix  ,iy  ,iz,is)=tpar (ix  ,iy  ,iz,is)+qp(np)*w1*wpar*wpar
            tpar (ix+1,iy  ,iz,is)=tpar (ix+1,iy  ,iz,is)+qp(np)*w2*wpar*wpar 
            tpar (ix  ,iy+1,iz,is)=tpar (ix  ,iy+1,iz,is)+qp(np)*w3*wpar*wpar 
            tpar (ix+1,iy+1,iz,is)=tpar (ix+1,iy+1,iz,is)+qp(np)*w4*wpar*wpar 
            tperp(ix  ,iy  ,iz,is)=tperp(ix  ,iy  ,iz,is)+qp(np)*w1*wperp2 
            tperp(ix+1,iy  ,iz,is)=tperp(ix+1,iy  ,iz,is)+qp(np)*w2*wperp2 
            tperp(ix  ,iy+1,iz,is)=tperp(ix  ,iy+1,iz,is)+qp(np)*w3*wperp2 
            tperp(ix+1,iy+1,iz,is)=tperp(ix+1,iy+1,iz,is)+qp(np)*w4*wperp2
            dpedx(ix  ,iy  ,iz)=dpedx(ix  ,iy  ,iz)+qp(np)*w1
            dpedx(ix+1,iy  ,iz)=dpedx(ix+1,iy  ,iz)+qp(np)*w2 
            dpedx(ix  ,iy+1,iz)=dpedx(ix  ,iy+1,iz)+qp(np)*w3 
            dpedx(ix+1,iy+1,iz)=dpedx(ix+1,iy+1,iz)+qp(np)*w4
 
            p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
            p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
            p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
            p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 

            p_xy (ix  ,iy  ,iz  ,is)=p_xy (ix  ,iy  ,iz  ,is)+qp(np)*w1*xy
            p_xy (ixp1,iy  ,iz  ,is)=p_xy (ixp1,iy  ,iz  ,is)+qp(np)*w2*xy 
            p_xy (ix  ,iyp1,iz  ,is)=p_xy (ix  ,iyp1,iz  ,is)+qp(np)*w3*xy 
            p_xy (ixp1,iyp1,iz  ,is)=p_xy (ixp1,iyp1,iz  ,is)+qp(np)*w4*xy 

            p_xz (ix  ,iy  ,iz  ,is)=p_xz (ix  ,iy  ,iz  ,is)+qp(np)*w1*xz
            p_xz (ixp1,iy  ,iz  ,is)=p_xz (ixp1,iy  ,iz  ,is)+qp(np)*w2*xz 
            p_xz (ix  ,iyp1,iz  ,is)=p_xz (ix  ,iyp1,iz  ,is)+qp(np)*w3*xz 
            p_xz (ixp1,iyp1,iz  ,is)=p_xz (ixp1,iyp1,iz  ,is)+qp(np)*w4*xz 

            p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
            p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
            p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
            p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
 
            p_yz (ix  ,iy  ,iz  ,is)=p_yz (ix  ,iy  ,iz  ,is)+qp(np)*w1*yz
            p_yz (ixp1,iy  ,iz  ,is)=p_yz (ixp1,iy  ,iz  ,is)+qp(np)*w2*yz 
            p_yz (ix  ,iyp1,iz  ,is)=p_yz (ix  ,iyp1,iz  ,is)+qp(np)*w3*yz 
            p_yz (ixp1,iyp1,iz  ,is)=p_yz (ixp1,iyp1,iz  ,is)+qp(np)*w4*yz 

            p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
            p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
            p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
            p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
 
            np=link(np)
          enddo
        enddo
      enddo
    enddo

    call xreal_2d(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(dpedx(1,jb-1,kb-1   ),NX,NY,NZ)

    call xreal_2d(p_xx (1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(p_xy (1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(p_xz (1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(p_yy (1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(p_yz (1,jb-1,kb-1,is),NX,NY,NZ)
    call xreal_2d(p_zz (1,jb-1,kb-1,is),NX,NY,NZ)

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

    do IIZ=KB-1,KE+1
      do IIY=JB-1,JE+1
        do IIX=1,NX2
          p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_xy(iix,iiy,iiz,is) = p_xy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_xz(iix,iiy,iiz,is) = p_xz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yz(iix,iiy,iiz,is) = p_yz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
        enddo
      enddo
    enddo
    p_xx(:,:,:,is)=p_xx(:,:,:,is)/(tx0(is)*frac(is))
    p_xy(:,:,:,is)=p_xy(:,:,:,is)/(tx0(is)*frac(is))
    p_xz(:,:,:,is)=p_xz(:,:,:,is)/(tx0(is)*frac(is))
    p_yy(:,:,:,is)=p_yy(:,:,:,is)/(tx0(is)*frac(is))
    p_yz(:,:,:,is)=p_yz(:,:,:,is)/(tx0(is)*frac(is))
    p_zz(:,:,:,is)=p_zz(:,:,:,is)/(tx0(is)*frac(is))

  enddo

  call date_and_time(values=time_end_array(:,26))
  call accumulate_time(time_begin_array(1,26),time_end_array(1,26),time_elapsed(26))

  ! do is=1,nspec
  !   call date_and_time(values=time_begin_array(:,24))
  !   call xreal_2d(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
  !   call xreal_2d(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
  !   call date_and_time(values=time_end_array(:,24))
  !   call accumulate_time(time_begin_array(1,24),time_end_array(1,24),time_elapsed(24))


  !   call date_and_time(values=time_begin_array(:,25))
  !   call xrealbcc_2d(tpar (1,jb-1,kb-1,is),1,NX,NY,NZ)
  !   call xrealbcc_2d(tperp(1,jb-1,kb-1,is),1,NX,NY,NZ)
  !   call date_and_time(values=time_end_array(:,25))
  !   call accumulate_time(time_begin_array(1,25),time_end_array(1,25),time_elapsed(25))
  ! enddo
  
  call date_and_time(values=time_begin_array(:,26))
 
  ! GOTO 10  

  ! do is=1,nspec
  !   do k=kb-1,ke+1
  !   do j = jb-1,je+1
  !       do i=1,nx2
  !         if(is.eq.1) then
  !           dns1=dns(i,j,k,1)/(dfac(1)*frac(1))
  !           dns2=0.
  !           denum=dns1+rfrac*dns2
  !         else
  !           denum=dns(i,j,k,is)/(dfac(is)*frac(is))
  !         endif
  !         if(den(i,j,k).le.denmin)  then
  !           tpar(i,j,k,is)=1.e-5
  !           tperp(i,j,k,is)=1.e-5
  !         else
  !           denum=denum*tx0(is)
  !           tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
  !           tperp(i,j,k,is)=0.5*tperp(i,j,k,is)*wspec(is)/denum
  !         endif
  !       enddo
  !     enddo
  !   enddo
  ! enddo

  ! 10   continue

  call date_and_time(values=time_end_array(:,26))
  call accumulate_time(time_begin_array(1,26),time_end_array(1,26),time_elapsed(26))

  call date_and_time(values=time_end_array(:,23))
  call accumulate_time(time_begin_array(1,23),time_end_array(1,23),time_elapsed(23))

  return
end subroutine caltemp2_global_2d


!---------------------------------------------------------------------
subroutine xrealbcc_pack_e_2d(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)
  use parameter_mod
  implicit none

  integer*8 :: ibnd,i,j,nx1m,ny1m,nz1m,k
  real*8 :: a_x(nxmax,jb-1:je+1,kb-1:ke+1)&
                  ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                  ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                  ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                  ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                  ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                  ,packed_data_xy_recv(nxmax,jb-1:je+1,3)

  a_x(:,:,kb-1)=a_x(:,:,kb)
  a_y(:,:,kb-1)=a_y(:,:,kb)
  a_z(:,:,kb-1)=a_z(:,:,kb)
  a_x(:,:,ke+1)=a_x(:,:,ke)
  a_y(:,:,ke+1)=a_y(:,:,ke)
  a_z(:,:,ke+1)=a_z(:,:,ke)

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
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

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
  do k=kb-1,ke+1
      do i=1,nxmax
        a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
        a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
        a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
      enddo
  enddo

  return
end subroutine xrealbcc_pack_e_2d


!---------------------------------------------------------------------
subroutine nsmth_2d (a,nx2m,ny2m,nz2m)
  use parameter_mod
  implicit none

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
real*8 function mahimax(x1, x2)
  implicit none
  real*8 :: x1, x2, x3
  if (x1.gt.x2) then
      x3=x1
  else
      x3=x2
  endif
  mahimax=x3
  return
end function mahimax


!---------------------------------------------------------------------   
real*8 function mahimin(x1, x2)
  implicit none
  real*8 :: x1, x2, x3
  if (x1.gt.x2) then
      x3=x2
  else
      x3=x1
  endif
  mahimin=x3
  return
end function mahimin

