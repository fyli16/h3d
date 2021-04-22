!---------------------------------------------------------------------
subroutine field_2d
  use m_parameters
  use m_mesh
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
  use m_parameters
  use m_mesh
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
  use m_parameters
  use m_mesh
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
  use m_parameters
  use m_mesh
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
  use m_parameters
  use m_mesh
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