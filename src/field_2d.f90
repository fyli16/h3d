!---------------------------------------------------------------------
subroutine field_2d
  use parameter_mod
  use mesh_mod
  implicit none

  call date_and_time(values=time_begin(:,9))
  call pressgrad_2d(1)
  call date_and_time(values=time_end(:,9))
  call accumulate_time(time_begin(1,9),time_end(1,9),time_elapsed(9))

  call date_and_time(values=time_begin(:,10))
  call bcalc_2d
  call date_and_time(values=time_end(:,10))
  call accumulate_time(time_begin(1,10),time_end(1,10),time_elapsed(10))

  call date_and_time(values=time_begin(:,9))
  call pressgrad_2d(0)
  call date_and_time(values=time_end(:,9))
  call accumulate_time(time_begin(1,9),time_end(1,9),time_elapsed(9))

  call date_and_time(values=time_begin(:,11))
  call ecalc_2d( 0)
  call date_and_time(values=time_end(:,11))
  call accumulate_time(time_begin(1,11),time_end(1,11),time_elapsed(11))

  call date_and_time(values=time_begin(:,12))
  call focalc_2d
  call date_and_time(values=time_end(:,12))
  call accumulate_time(time_begin(1,12),time_end(1,12),time_elapsed(12))

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
subroutine ecalc_2d(iflag)
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
  call date_and_time(values=time_begin(:,18))
  call XREALBCC_PACK_E_2D(EX,EY,EZ,1_8,NX,NY,NZ)
  call date_and_time(values=time_end(:,18))
  call accumulate_time(time_begin(1,18),time_end(1,18),time_elapsed(18))

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

  call date_and_time(values=time_begin(:,22))

  dts=dt/real(iterb)
  dts2=dts/2.
  dts6=dts/6.

  ! subcycle into iterb interations
  do ii = 1, iterb
    bxs=bx
    bys=by
    bzs=bz

    ! R-K first part
    call date_and_time(values=time_begin(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end(:,16))
    call accumulate_time(time_begin(1,16),time_end(1,16),time_elapsed(16))
 
    ! B = B(n)+dt*K1/2
    bx=bxs-dts2*curlex
    by=bys-dts2*curley
    bz=bzs-dts2*curlez

    ! temp1 = K1
    tempx1=curlex
    tempy1=curley
    tempz1=curlez

    ! R-K part 2
    call date_and_time(values=time_begin(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end(:,16))
    call accumulate_time(time_begin(1,16),time_end(1,16),time_elapsed(16))
 
    ! B = B(n)+dt*K2/2
    bx=bxs-dts2*curlex
    by=bys-dts2*curley
    bz=bzs-dts2*curlez

    ! temp2 = K2
    tempx1=tempx1+2.*curlex
    tempy1=tempy1+2.*curley
    tempz1=tempz1+2.*curlez

    ! R-K  part 3
    call date_and_time(values=time_begin(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end(:,16))
    call accumulate_time(time_begin(1,16),time_end(1,16),time_elapsed(16))

    ! B = B(n)+dt*K3
    bx=bxs-dts*curlex
    by=bys-dts*curley
    bz=bzs-dts*curlez

    ! temp3 = K3
    tempx1=tempx1+2.*curlex
    tempy1=tempy1+2.*curley
    tempz1=tempz1+2.*curlez

    ! R-K  part 4
    call date_and_time(values=time_begin(:,16))
    call ecalc_2d( 1 )
    call date_and_time(values=time_end(:,16)) 
    call accumulate_time(time_begin(1,16),time_end(1,16),time_elapsed(16))
 
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
 
  call date_and_time(values=time_end(:,22))
  call accumulate_time(time_begin(1,22),time_end(1,22),time_elapsed(22))

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

  call date_and_time(values=time_begin(:,23))
 
  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt

  tpar  = 0.
  tperp = 0.
  rfrac = 0.

  p_xx=0.;p_xy=0.;p_xz=0.;p_yy=0.;p_yz=0.;p_zz=0.

  if (nspec >= 2) rfrac = frac(2)/frac(1)

  call date_and_time(values=time_begin(:,26))

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

  call date_and_time(values=time_end(:,26))
  call accumulate_time(time_begin(1,26),time_end(1,26),time_elapsed(26))

  ! do is=1,nspec
  !   call date_and_time(values=time_begin(:,24))
  !   call xreal_2d(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
  !   call xreal_2d(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
  !   call date_and_time(values=time_end(:,24))
  !   call accumulate_time(time_begin(1,24),time_end(1,24),time_elapsed(24))


  !   call date_and_time(values=time_begin(:,25))
  !   call xrealbcc_2d(tpar (1,jb-1,kb-1,is),1,NX,NY,NZ)
  !   call xrealbcc_2d(tperp(1,jb-1,kb-1,is),1,NX,NY,NZ)
  !   call date_and_time(values=time_end(:,25))
  !   call accumulate_time(time_begin(1,25),time_end(1,25),time_elapsed(25))
  ! enddo
  
  call date_and_time(values=time_begin(:,26))
 
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

  call date_and_time(values=time_end(:,26))
  call accumulate_time(time_begin(1,26),time_end(1,26),time_elapsed(26))

  call date_and_time(values=time_end(:,23))
  call accumulate_time(time_begin(1,23),time_end(1,23),time_elapsed(23))

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

