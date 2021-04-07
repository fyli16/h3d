!---------------------------------------------------------------------
! advances electromagnetic field in time
!---------------------------------------------------------------------
subroutine field
  use parameter_mod
  implicit none

  call date_and_time(values=time_begin_array(:,21))

  call date_and_time(values=time_begin_array(:,9))
  call pressgrad(1)
  call date_and_time(values=time_end_array(:,9))
  call accumulate_time(time_begin_array(1,9), time_end_array(1,9), time_elapsed(9))

  call date_and_time(values=time_begin_array(:,10))
  call bcalc
  call date_and_time(values=time_end_array(:,10))
  call accumulate_time(time_begin_array(1,10), time_end_array(1,10), time_elapsed(10))

  call date_and_time(values=time_begin_array(:,9))
  call pressgrad(0)
  call date_and_time(values=time_end_array(:,9))
  call accumulate_time(time_begin_array(1,9), time_end_array(1,9), time_elapsed(9))

  call date_and_time(values=time_begin_array(:,11))
  call ecalc( 0 )
  call date_and_time(values=time_end_array(:,11))
  call accumulate_time(time_begin_array(1,11), time_end_array(1,11), time_elapsed(11))

  call date_and_time(values=time_begin_array(:,12))
  call focalc
  call date_and_time(values=time_end_array(:,12))
  call accumulate_time(time_begin_array(1,12), time_end_array(1,12), time_elapsed(12))

  call date_and_time(values=time_end_array(:,21))
  call accumulate_time(time_begin_array(1,21), time_end_array(1,21), time_elapsed(21))

  return
end subroutine field


!---------------------------------------------------------------------
! computes electron pressure gradient
! @param iflag: memory flag?
!---------------------------------------------------------------------
subroutine pressgrad(iflag)
  use parameter_mod
  use mesh2d
  implicit none

  integer :: iflag
  integer*8 :: i,j,k
  real*8 :: dena,dxa,dya,dza,a
  
  if (iflag==0) then 
    pe = te0*denh**gama
  else
    pe = te0*den**gama
  endif

  do k=kb,ke
    do j = jb,je
      do i=2,nx1
        ! dena=iflag*0.5*(den(i,j,k) + deno(i,j,k)) + (1.-iflag)*den(i,j,k)
        dena = iflag*den(i,j,k)+(1-iflag)*denh(i,j,k)
        a = one/dena

        ! Uniform mesh - Same as in version 5.0
        ! dxa=a/(4.*hx)
        ! dya=a/(4.*hy)
        ! dza=a/(4.*hz)

        ! Nonuniform mesh
        ! dxa=a/(2.*(meshX%dxn(i  )+meshX%dxn(i+1)))
        ! dya=a/(2.*(meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
        ! dza=a/(2.*(meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0
        dxa=a/(4.*(meshX%dxn(i  )+meshX%dxn(i+1)))
        dya=a/(4.*(meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
        dza=a/(4.*(meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0

        dpedx(i,j,k) = ((pe(i+1,j-1,k+1)+2.*pe(i+1,j,k+1) &
              + pe(i+1,j+1,k+1))/4. &
              + 2.*(pe(i+1,j-1,k  )+2.*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4. &
              + (pe(i+1,j-1,k-1)+2.*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4. &
              - (pe(i-1,j-1,k+1)+2.*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4. &
              - 2.*(pe(i-1,j-1,k  )+2.*pe(i-1,j,k  )+pe(i-1,j+1,k  ))/4. &
              - (pe(i-1,j-1,k-1)+2.*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4.) * dxa
        dpedy(i,j,k)=((pe(i-1,j+1,k+1)+2.*pe(i,j+1,k+1)&
              + pe(i+1,j+1,k+1))/4.&
              + 2.*(pe(i-1,j+1,k  )+2.*pe(i,j+1,k  )+pe(i+1,j+1,k  ))/4. &
              +   (pe(i-1,j+1,k-1)+2.*pe(i,j+1,k-1)+pe(i+1,j+1,k-1))/4. &
              -   (pe(i-1,j-1,k+1)+2.*pe(i,j-1,k+1)+pe(i+1,j-1,k+1))/4. &
              -2.*(pe(i-1,j-1,k  )+2.*pe(i,j-1,k  )+pe(i+1,j-1,k  ))/4. &
              -   (pe(i-1,j-1,k-1)+2.*pe(i,j-1,k-1)+pe(i+1,j-1,k-1))/4.) * dya
        dpedz(i,j,k)=((pe(i+1,j-1,k+1)+2.*pe(i+1,j,k+1) &
              + pe(i+1,j+1,k+1))/4. & 
              + 2.*(pe(i  ,j-1,k+1)+2.*pe(i  ,j,k+1)+pe(i  ,j+1,k+1))/4. &
              + (pe(i-1,j-1,k+1)+2.*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4. &
              - (pe(i+1,j-1,k-1)+2.*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4. &
              - 2.*(pe(i  ,j-1,k-1)+2.*pe(i  ,j,k-1)+pe(i  ,j+1,k-1))/4. &
              - (pe(i-1,j-1,k-1)+2.*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4.) * dza                    
      enddo
    enddo
  enddo

  return
end subroutine pressgrad


!---------------------------------------------------------------------
! computes electric field and curl(E)
!---------------------------------------------------------------------
subroutine ecalc( iflag )
  use parameter_mod
  use mesh2d
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

  if (eta_par == 0) then
    do k = kb,ke
      do j = jb,je
        do i=2,nx1
          bx1=bx(i+1,j+1,k)  
          bx2=bx(i  ,j+1,k)  
          bx3=bx(i  ,j  ,k)  
          bx4=bx(i+1,j  ,k)  
          bx5=bx(i+1,j+1,k+1)
          bx6=bx(i  ,j+1,k+1)
          bx7=bx(i  ,j  ,k+1)
          bx8=bx(i+1,j  ,k+1)
          by1=by(i+1,j+1,k)  
          by2=by(i  ,j+1,k)  
          by3=by(i  ,j  ,k)  
          by4=by(i+1,j  ,k)  
          by5=by(i+1,j+1,k+1)
          by6=by(i  ,j+1,k+1)
          by7=by(i  ,j  ,k+1)
          by8=by(i+1,j  ,k+1)
          bz1=bz(i+1,j+1,k)  
          bz2=bz(i  ,j+1,k)  
          bz3=bz(i  ,j  ,k)  
          bz4=bz(i+1,j  ,k)  
          bz5=bz(i+1,j+1,k+1)
          bz6=bz(i  ,j+1,k+1)
          bz7=bz(i  ,j  ,k+1)
          bz8=bz(i+1,j  ,k+1)

          vixa=(1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k))&
              +iflag*vix(i,j,k)
          viya=(1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k))&
              +iflag*viy(i,j,k)
          viza=(1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k))&
              +iflag*viz(i,j,k)

          ! dena = iflag*0.5*(den(i,j,k) + deno(i,j,k)) + (1.-iflag)*den(i,j,k)
          dena = iflag*den(i,j,k) + (1-iflag)*denh(i,j,k)
          a = one/dena

          ! Uniform mesh - Same as is in version 5.0
          ! dxa=a/(4.*hx)
          ! dya=a/(4.*hy)
          ! dza=a/(4.*hz)

          ! Nonuniform mesh
          dxa=a/(4.*meshX%dxc(i))
          dya=a/(4.*meshY%dxc(j+1))  ! integer index in y direction starts at 0
          dza=a/(4.*meshZ%dxc(k+1))  ! integer index in z direction starts at 0

          dbxdy= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
              +bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
              -bx(i+1,j  ,k+1)-bx(i  ,j  ,k+1)&
              -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
          dbxdz= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
              +bx(i  ,j  ,k+1)+bx(i+1,j  ,k+1)&
              -bx(i+1,j+1,k  )-bx(i  ,j+1,k  )&
              -bx(i  ,j  ,k  )-bx(i+1,j ,k  )
          dbydx= by(i+1,j+1,k+1)+by(i+1,j  ,k+1)&
              +by(i+1,j  ,k  )+by(i+1,j+1,k  )&
              -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)&
              -by(i  ,j  ,k  )-by(i  ,j+1,k  )
          dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)&
              +by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
              -by(i+1,j+1,k  )-by(i  ,j+1,k  )&
              -by(i  ,j  ,k  )-by(i+1,j  ,k  )
          dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)&
              +bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
              -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)&
              -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
          dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)&
              +bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
              -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)&
              -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
          curlbx_scalar=dya*dbzdy-dza*dbydz
          curlby_scalar=dza*dbxdz-dxa*dbzdx
          curlbz_scalar=dxa*dbydx-dya*dbxdy
          bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
          byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
          bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
          xj = curlbx_scalar
          yj = curlby_scalar
          zj = curlbz_scalar

          ! From the eta_par conditional
          tenx = eta(i,j,k)*xj
          teny = eta(i,j,k)*yj
          tenz = eta(i,j,k)*zj
          ! End content from the eta_par conditional

          ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                  -dpedx(i,j,k)+tenx/a
          ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                  -dpedy(i,j,k)+teny/a
          ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                  -dpedz(i,j,k)+tenz/a
        enddo
      enddo
    enddo
  else 
    if (eta_par == 1) then
      do k = kb,ke
        do j = jb,je
          do i=2,nx1
            bx1=bx(i+1,j+1,k)  
            bx2=bx(i  ,j+1,k)  
            bx3=bx(i  ,j  ,k)  
            bx4=bx(i+1,j  ,k)  
            bx5=bx(i+1,j+1,k+1)
            bx6=bx(i  ,j+1,k+1)
            bx7=bx(i  ,j  ,k+1)
            bx8=bx(i+1,j  ,k+1)
            by1=by(i+1,j+1,k)  
            by2=by(i  ,j+1,k)  
            by3=by(i  ,j  ,k)  
            by4=by(i+1,j  ,k)  
            by5=by(i+1,j+1,k+1)
            by6=by(i  ,j+1,k+1)
            by7=by(i  ,j  ,k+1)
            by8=by(i+1,j  ,k+1)
            bz1=bz(i+1,j+1,k)  
            bz2=bz(i  ,j+1,k)  
            bz3=bz(i  ,j  ,k)  
            bz4=bz(i+1,j  ,k)  
            bz5=bz(i+1,j+1,k+1)
            bz6=bz(i  ,j+1,k+1)
            bz7=bz(i  ,j  ,k+1)
            bz8=bz(i+1,j  ,k+1)

            vixa=(1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k))&
                +iflag*vix(i,j,k)
            viya=(1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k))&
                +iflag*viy(i,j,k)
            viza=(1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k))&
                +iflag*viz(i,j,k)

            dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))&
                +(1.-iflag)*den(i,j,k)
            a=one/dena

            ! Uniform mesh - Same as is in version 5.0
            ! dxa=a/(4.*hx)
            ! dya=a/(4.*hy)
            ! dza=a/(4.*hz)

            ! Nonuniform mesh
            dxa=a/(4.*meshX%dxc(i))
            dya=a/(4.*meshY%dxc(j+1))  ! integer index in y direction starts at 0
            dza=a/(4.*meshZ%dxc(k+1))  ! integer index in z direction starts at 0

            dbxdy= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                +bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
                -bx(i+1,j  ,k+1)-bx(i  ,j  ,k+1)&
                -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
            dbxdz= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                +bx(i  ,j  ,k+1)+bx(i+1,j  ,k+1)&
                -bx(i+1,j+1,k  )-bx(i  ,j+1,k  )&
                -bx(i  ,j  ,k  )-bx(i+1,j ,k  )
            dbydx= by(i+1,j+1,k+1)+by(i+1,j  ,k+1)&
                +by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)&
                -by(i  ,j  ,k  )-by(i  ,j+1,k  )
            dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)&
                +by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
                -by(i+1,j+1,k  )-by(i  ,j+1,k  )&
                -by(i  ,j  ,k  )-by(i+1,j  ,k  )
            dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)&
                +bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)&
                -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
            dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)&
                +bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)&
                -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
            curlbx_scalar=dya*dbzdy-dza*dbydz
            curlby_scalar=dza*dbxdz-dxa*dbzdx
            curlbz_scalar=dxa*dbydx-dya*dbxdy
            bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar

            ! From eta_par conditional
            bxx = bxav
            byy = byav
            bzz = bzav
            btot = sqrt(bxx**2 + byy**2 + bzz**2)
            tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
            tenx = tjdotb*bxx/btot
            teny = tjdotb*byy/btot
            tenz = tjdotb*bzz/btot
            ! End content from eta_par conditional

            ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                    -dpedx(i,j,k)+tenx/a
            ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                    -dpedy(i,j,k)+teny/a
            ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                    -dpedz(i,j,k)+tenz/a
          enddo
        enddo
      enddo
    else if (eta_par == 2) then
      do k = kb,ke
        do j = jb,je
          do i=2,nx1
            bx1=bx(i+1,j+1,k)  
            bx2=bx(i  ,j+1,k)  
            bx3=bx(i  ,j  ,k)  
            bx4=bx(i+1,j  ,k)  
            bx5=bx(i+1,j+1,k+1)
            bx6=bx(i  ,j+1,k+1)
            bx7=bx(i  ,j  ,k+1)
            bx8=bx(i+1,j  ,k+1)
            by1=by(i+1,j+1,k)  
            by2=by(i  ,j+1,k)  
            by3=by(i  ,j  ,k)  
            by4=by(i+1,j  ,k)  
            by5=by(i+1,j+1,k+1)
            by6=by(i  ,j+1,k+1)
            by7=by(i  ,j  ,k+1)
            by8=by(i+1,j  ,k+1)
            bz1=bz(i+1,j+1,k)  
            bz2=bz(i  ,j+1,k)  
            bz3=bz(i  ,j  ,k)  
            bz4=bz(i+1,j  ,k)  
            bz5=bz(i+1,j+1,k+1)
            bz6=bz(i  ,j+1,k+1)
            bz7=bz(i  ,j  ,k+1)
            bz8=bz(i+1,j  ,k+1)

            vixa=(1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k))&
                +iflag*vix(i,j,k)
            viya=(1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k))&
                +iflag*viy(i,j,k)
            viza=(1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k))&
                +iflag*viz(i,j,k)

            dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))&
                +(1.-iflag)*den(i,j,k)
            a=one/dena

            ! Uniform mesh - Same as is in version 5.0
            ! dxa=a/(4.*hx)
            ! dya=a/(4.*hy)
            ! dza=a/(4.*hz)

            ! Nonuniform mesh
            dxa=a/(4.*meshX%dxc(i))
            dya=a/(4.*meshY%dxc(j+1))  ! integer index in y direction starts at 0
            dza=a/(4.*meshZ%dxc(k+1))  ! integer index in z direction starts at 0

            dbxdy= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                +bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
                -bx(i+1,j  ,k+1)-bx(i  ,j  ,k+1)&
                -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
            dbxdz= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                +bx(i  ,j  ,k+1)+bx(i+1,j  ,k+1)&
                -bx(i+1,j+1,k  )-bx(i  ,j+1,k  )&
                -bx(i  ,j  ,k  )-bx(i+1,j ,k  )
            dbydx= by(i+1,j+1,k+1)+by(i+1,j  ,k+1)&
                +by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)&
                -by(i  ,j  ,k  )-by(i  ,j+1,k  )
            dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)&
                +by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
                -by(i+1,j+1,k  )-by(i  ,j+1,k  )&
                -by(i  ,j  ,k  )-by(i+1,j  ,k  )
            dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)&
                +bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)&
                -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
            dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)&
                +bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)&
                -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
            curlbx_scalar=dya*dbzdy-dza*dbydz
            curlby_scalar=dza*dbxdz-dxa*dbzdx
            curlbz_scalar=dxa*dbydx-dya*dbxdy
            bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar

            ! From eta_par conditional
            bxx = bxav
            byy = byav
            bzz = bzav
            btot = sqrt(bxx**2 + byy**2 + bzz**2)
            curr_tot = max(1.d-12,sqrt(xj**2 + yj**2 + zj**2))
            tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
            tenx = min(resis,tenx)
            tenx = tenx*xj
            teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
            teny = min(resis,teny)
            teny = teny*yj
            tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
            tenz = min(resis,tenz)
            tenz = tenz*zj
              ! End content from eta_par conditional

            ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                    -dpedx(i,j,k)+tenx/a
            ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                    -dpedy(i,j,k)+teny/a
            ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                    -dpedz(i,j,k)+tenz/a
          enddo
        enddo
      enddo
    endif
  endif

  ! boundary conditions
  call date_and_time(values=time_begin_array(:,18))
  call xrealbcc_pack_e(EX,EY,EZ,1_8,NX,NY,NZ)
  call date_and_time(values=time_end_array(:,18))
  call accumulate_time(time_begin_array(1,18),time_end_array(1,18),time_elapsed(18))

  ! VR: impose periodic B.C. on E (in x)
  ! this should be part of xrealbcc*

  !ex(nx2,:,:)=ex(2,:,:)
  !ey(nx2,:,:)=ey(2,:,:)
  !ez(nx2,:,:)=ez(2,:,:)

  !ex(1,:,:)=ex(nx1,:,:)
  !ey(1,:,:)=ey(nx1,:,:)
  !ez(1,:,:)=ez(nx1,:,:)

  ! do k=kb-1,ke+1
  !   do j = jb-1,je+1
  !     ex(nx2,j,k)=ex(2,j,k)
  !     ey(nx2,j,k)=ey(2,j,k)
  !     ez(nx2,j,k)=ez(2,j,k)

  !     ex(1,j,k)=ex(nx1,j,k)
  !     ey(1,j,k)=ey(nx1,j,k)
  !     ez(1,j,k)=ez(nx1,j,k)
  !   enddo
  ! enddo

  ! calculate curl E
  do k=kb,ke+1
    do j = jb,je+1
      do i=2,nx2
        dexdy=  ex(i  ,j  ,k  ) + ex(i-1,j  ,k  )   &
              + ex(i-1,j  ,k-1) + ex(i  ,j  ,k-1)   &
              - ex(i  ,j-1,k  ) - ex(i-1,j-1,k  )   &
              - ex(i-1,j-1,k-1) - ex(i  ,j-1,k-1)
        dexdz=  ex(i  ,j  ,k  ) + ex(i-1,j  ,k  )   &
              + ex(i-1,j-1,k  ) + ex(i  ,j-1,k  )   &
              - ex(i  ,j  ,k-1) - ex(i-1,j  ,k-1)   &
              - ex(i-1,j-1,k-1) - ex(i  ,j-1,k-1)
        deydx=  ey(i  ,j  ,k  ) + ey(i  ,j-1,k  )   &
              + ey(i  ,j-1,k-1) + ey(i  ,j  ,k-1)   &
              - ey(i-1,j  ,k  ) - ey(i-1,j-1,k  )   &
              - ey(i-1,j-1,k-1) - ey(i-1,j  ,k-1)
        deydz=  ey(i  ,j  ,k  ) + ey(i-1,j  ,k  )   &
              + ey(i-1,j-1,k  ) + ey(i  ,j-1,k  )   &
              - ey(i  ,j  ,k-1) - ey(i-1,j  ,k-1)   &
              - ey(i-1,j-1,k-1) - ey(i  ,j-1,k-1)
        dezdx=  ez(i  ,j  ,k  ) + ez(i  ,j-1,k  )   &
              + ez(i  ,j-1,k-1) + ez(i  ,j  ,k-1)   &
              - ez(i-1,j  ,k  ) - ez(i-1,j-1,k  )   &
              - ez(i-1,j-1,k-1) - ez(i-1,j  ,k-1)
        dezdy=  ez(i  ,j  ,k  ) + ez(i-1,j  ,k  )   &
              + ez(i-1,j  ,k-1) + ez(i  ,j  ,k-1)   &
              - ez(i  ,j-1,k  ) - ez(i-1,j-1,k  )   &
              - ez(i-1,j-1,k-1) - ez(i  ,j-1,k-1)

        ! Uniform mesh - Same as is in version 5.0
        ! curlex(i,j,k)=dezdy/(4.*hy)-deydz/(4.*hz)
        ! curley(i,j,k)=dexdz/(4.*hz)-dezdx/(4.*hx)
        ! curlez(i,j,k)=deydx/(4.*hx)-dexdy/(4.*hy)

        ! Nonuniform mesh
        curlex(i,j,k)=dezdy/(4.*meshY%dxn(j+1))-deydz/(4.*meshZ%dxn(k+1))        ! integer index in y and z directions start  at 0
        curley(i,j,k)=dexdz/(4.*meshZ%dxn(k+1))-dezdx/(4.*meshX%dxn(i  ))        ! integer index in z       direction  starts at 0
        curlez(i,j,k)=deydx/(4.*meshX%dxn(i  ))-dexdy/(4.*meshY%dxn(j+1))        ! integer index in y       direction  starts at 0
      enddo
    enddo
  enddo

  return
end subroutine ecalc


!---------------------------------------------------------------------
! advances magnetic field
!---------------------------------------------------------------------
subroutine bcalc
  use parameter_mod
  use mesh2d
  implicit none

  integer*8 :: i,j,k,ii
  real*8 :: dts,dts2,dts6
  real*8 :: tempx1(nxmax,jb-1:je+1,kb-1:ke+1) &
          ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1) &
          ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)

  call date_and_time(values=time_begin_array(:,22))
  
  dts=dt/real(iterb)
  dts2=dts/2.
  dts6=dts/6.

  ! subcycle into iterb interations
  !VR : E is synchronized between processors at each step of RK.
  !VR : but is it enough to make sure that B is consistent?
  do ii = 1, iterb
    ! save B at start of subcycle
    ! Bs = B(n)
    bxs=bx
    bys=by
    bzs=bz

    ! R-K first part
    call date_and_time(values=time_begin_array(:,16))
    call ecalc( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16) &
                                  ,time_end_array(1,16) &
                                  ,time_elapsed(16))
      
    ! B = B(n)+dt*K1/2
    do k=kb,ke+1
      do j = jb,je+1
          do i=2,nx2
            bx(i,j,k)=bxs(i,j,k)-dts2*curlex(i,j,k)
            by(i,j,k)=bys(i,j,k)-dts2*curley(i,j,k)
            bz(i,j,k)=bzs(i,j,k)-dts2*curlez(i,j,k)
          enddo
      enddo
    enddo
      
    ! temp1 = K1
    do k=kb,ke+1
      do j = jb,je+1
          do i=2,nx2
            tempx1(i,j,k)=curlex(i,j,k)
            tempy1(i,j,k)=curley(i,j,k)
            tempz1(i,j,k)=curlez(i,j,k)
          enddo
      enddo
    enddo
      
    ! R-K part 2
    call date_and_time(values=time_begin_array(:,16))
    call ecalc( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16) &
                                  ,time_end_array(1,16) &
                                  ,time_elapsed(16))
      
    ! B = B(n)+dt*K2/2
    do k=kb,ke+1
      do j = jb,je+1
        do i=2,nx2
          bx(i,j,k)=bxs(i,j,k)-dts2*curlex(i,j,k)
          by(i,j,k)=bys(i,j,k)-dts2*curley(i,j,k)
          bz(i,j,k)=bzs(i,j,k)-dts2*curlez(i,j,k)
        enddo
      enddo
    enddo
      
    ! temp2 = K2
    do k=kb,ke+1
      do j = jb,je+1
        do i=2,nx2
          tempx1(i,j,k)=tempx1(i,j,k)+2.*curlex(i,j,k)
          tempy1(i,j,k)=tempy1(i,j,k)+2.*curley(i,j,k)
          tempz1(i,j,k)=tempz1(i,j,k)+2.*curlez(i,j,k)
        enddo
      enddo
    enddo
      
    ! R-K part 3
    call date_and_time(values=time_begin_array(:,16))
    call ecalc( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16) &
                                  ,time_end_array(1,16) &
                                  ,time_elapsed(16))
      
    ! B = B(n)+dt*K3
    do k=kb,ke+1
      do j = jb,je+1
        do i=2,nx2
          bx(i,j,k)=bxs(i,j,k)-dts*curlex(i,j,k)
          by(i,j,k)=bys(i,j,k)-dts*curley(i,j,k)
          bz(i,j,k)=bzs(i,j,k)-dts*curlez(i,j,k)
        enddo
      enddo
    enddo
      
    ! temp3 = K3
    do k=kb,ke+1
      do j = jb,je+1
        do i=2,nx2
          tempx1(i,j,k)=tempx1(i,j,k)+2.*curlex(i,j,k)
          tempy1(i,j,k)=tempy1(i,j,k)+2.*curley(i,j,k)
          tempz1(i,j,k)=tempz1(i,j,k)+2.*curlez(i,j,k)
        enddo
      enddo
    enddo
      
    ! R-K  part 4
    call date_and_time(values=time_begin_array(:,16))
    call ecalc( 1 )
    call date_and_time(values=time_end_array(:,16))
    call accumulate_time(time_begin_array(1,16) &
                                  ,time_end_array(1,16) &
                                  ,time_elapsed(16))

    ! B = B(n) + dt*(K1+2K2+2K3+K4)/6
    do k=kb,ke+1
      do j = jb,je+1
        do i=2,nx2
          bx(i,j,k)=bxs(i,j,k)-dts6*(tempx1(i,j,k)+curlex(i,j,k))
          by(i,j,k)=bys(i,j,k)-dts6*(tempy1(i,j,k)+curley(i,j,k))
          bz(i,j,k)=bzs(i,j,k)-dts6*(tempz1(i,j,k)+curlez(i,j,k))
        enddo
      enddo
    enddo
      
    ! end of iteration loop
    CALL XREALBCC_PACK_B(BX,BY,BZ,1_8,NX,NY,NZ)  

  end do


  ! this is part of xreal*
  !  set ghost cell values for B: it is only for diagnostics as these B
  !  values are not used anywhere
  ! do k=kb-1,ke+1
  !    do j = jb-1,je+1
  !       bx(1  ,j,k)=bx(nx1  ,j,k)
  !       by(1  ,j,k)=by(nx1  ,j,k)
  !       bz(1  ,j,k)=bz(nx1  ,j,k)
        
  !       bx(nx2  ,j,k)=bx(2  ,j,k)
  !       by(nx2  ,j,k)=by(2  ,j,k)
  !       bz(nx2  ,j,k)=bz(2  ,j,k)
  !    enddo
  ! enddo
  
  
  call date_and_time(values=time_end_array(:,22))
  call accumulate_time(time_begin_array(1,22) &
                                ,time_end_array(1,22) &
                                ,time_elapsed(22))

  return
end subroutine bcalc


!---------------------------------------------------------------------
! ??
!---------------------------------------------------------------------
subroutine focalc
  use parameter_mod
  use mesh2d
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
    do j = jb, je
      do i = 2, nx1
        dbxdy= bx(i,j  ,k)+bx(i-1,j  ,k)&
              +bx(i-1,j  ,k-1)+bx(i,j  ,k-1)&
              -bx(i,j-1,k)-bx(i-1,j-1,k)&
              -bx(i-1,j-1,k-1)-bx(i,j-1,k-1)
        dbxdz= bx(i,j,k  )+bx(i-1,j,k  )&
              +bx(i-1,j-1,k  )+bx(i,j-1,k  )&
              -bx(i,j,k-1)-bx(i-1,j,k-1)&
              -bx(i-1,j-1,k-1)-bx(i,j-1,k-1)
        dbydx= by(i  ,j,k)+by(i  ,j-1,k)&
              +by(i  ,j-1,k-1)+by(i  ,j,k-1)&
              -by(i-1,j,k)-by(i-1,j-1,k)&
              -by(i-1,j-1,k-1)-by(i-1,j,k-1)
        dbydz= by(i,j,k  )+by(i-1,j,k  )&
              +by(i-1,j-1,k  )+by(i,j-1,k  )&
              -by(i,j,k-1)-by(i-1,j,k-1)&
              -by(i-1,j-1,k-1)-by(i,j-1,k-1)
        dbzdx= bz(i  ,j,k)+bz(i  ,j-1,k)&
              +bz(i  ,j-1,k-1)+bz(i  ,j,k-1)&
              -bz(i-1,j,k)-bz(i-1,j-1,k)&
              -bz(i-1,j-1,k-1)-bz(i-1,j,k-1)
        dbzdy= bz(i,j  ,k)+bz(i-1,j  ,k)&
              +bz(i-1,j  ,k-1)+bz(i,j  ,k-1)&
              -bz(i,j-1,k)-bz(i-1,j-1,k)&
              -bz(i-1,j-1,k-1)-bz(i,j-1,k-1)

        ! Uniform mesh - Same as is in version 5.0
        ! curlbx_scalar=dbzdy/(4.*hy)-dbydz/(4.*hz)
        ! curlby_scalar=dbxdz/(4.*hz)-dbzdx/(4.*hx)
        ! curlbz_scalar=dbydx/(4.*hx)-dbxdy/(4.*hy)


        ! Nonuniform mesh
        curlbx_scalar=dbzdy/(4.*meshY%dxc(j+1))-dbydz/(4.*meshZ%dxc(k+1))
        curlby_scalar=dbxdz/(4.*meshZ%dxc(k+1))-dbzdx/(4.*meshX%dxc(i  ))
        curlbz_scalar=dbydx/(4.*meshX%dxc(i  ))-dbxdy/(4.*meshY%dxc(j+1))

        ! 6/25/2006 New eta_par option: tensor eta
        bx1=bx(i+1,j+1,k)  
        bx2=bx(i  ,j+1,k)  
        bx3=bx(i  ,j  ,k)  
        bx4=bx(i+1,j  ,k)  
        bx5=bx(i+1,j+1,k+1)
        bx6=bx(i  ,j+1,k+1)
        bx7=bx(i  ,j  ,k+1)
        bx8=bx(i+1,j  ,k+1)
        by1=by(i+1,j+1,k)  
        by2=by(i  ,j+1,k)  
        by3=by(i  ,j  ,k)  
        by4=by(i+1,j  ,k)  
        by5=by(i+1,j+1,k+1)
        by6=by(i  ,j+1,k+1)
        by7=by(i  ,j  ,k+1)
        by8=by(i+1,j  ,k+1)
        bz1=bz(i+1,j+1,k)  
        bz2=bz(i  ,j+1,k)  
        bz3=bz(i  ,j  ,k)  
        bz4=bz(i+1,j  ,k)  
        bz5=bz(i+1,j+1,k+1)
        bz6=bz(i  ,j+1,k+1)
        bz7=bz(i  ,j  ,k+1)
        bz8=bz(i+1,j  ,k+1)
        bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
        byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
        bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
        xj = curlbx_scalar
        yj = curlby_scalar
        zj = curlbz_scalar
        tenx=eta(i,j,k)*xj
        teny=eta(i,j,k)*yj
        tenz=eta(i,j,k)*zj
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
  call XREALBCC_PACK_E(fox,foy,foz,1_8,NX,NY,NZ)

  !VR : periodic boundary conditions in x 
  ! fox(1  ,:,:)=fox(nx1  ,:,:)
  ! foy(1  ,:,:)=foy(nx1  ,:,:)
  ! foz(1  ,:,:)=foz(nx1  ,:,:)
  ! fox(nx2,:,:)=fox(2,:,:)
  ! foy(nx2,:,:)=foy(2,:,:)
  ! foz(nx2,:,:)=foz(2,:,:)

  !VR: this is not periodic. 
  ! do k=kb-1,ke+1 
  !   do i=1,nx2
  !     if (jb == 1) then
  !       fox(i,jb-1,k)=fox(i,jb,k)
  !       foy(i,jb-1,k)=foy(i,jb,k)
  !       foz(i,jb-1,k)=foz(i,jb,k)
  !     endif
  !     if (je == ny) then
  !       fox(i,je+1,k)=fox(i,je,k)
  !       foy(i,je+1,k)=foy(i,je,k)
  !       foz(i,je+1,k)=foz(i,je,k)
  !     endif
  !   enddo
  ! enddo

  !VR: this is not periodic
  ! if (kb == 1) then
  !   do j = jb-1,je+1
  !     do i=1,nx2
  !       fox(i,j,kb-1)=fox(i,j,kb)
  !       foy(i,j,kb-1)=foy(i,j,kb)
  !       foz(i,j,kb-1)=foz(i,j,kb)
  !     enddo
  !   enddo
  ! endif

  ! if (ke == nz) then
  !   do j = jb-1,je+1
  !     do i=1,nx2
  !       fox(i,j,ke+1)=fox(i,j,ke)
  !       foy(i,j,ke+1)=foy(i,j,ke)
  !       foz(i,j,ke+1)=foz(i,j,ke)
  !     enddo
  !   enddo
  ! endif

  return
end subroutine focalc
