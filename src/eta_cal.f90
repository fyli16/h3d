!---------------------------------------------------------------------
!> calculate a resistivity that depends on physical quantities
! remember, eta and fo{x,y,z} are cell-centered quantities
! ghost cells are only used in diagnostic output
! Here five choices, based on argument ieta:
! (1) constant
! (2) 4th power of B-grad, 2nd power of density
! (3) Arbitrary power of curent (here: 4)
! (4) Homas method: 2nd derivative of current
! (5) (2) and (4) combined, with (4) reduced by a factor (1/5)
!---------------------------------------------------------------------
subroutine eta_calc
  use parameter_mod
  use mesh_mod
  implicit none

  real*8 :: ajl(nxmax, jb-1:jb+nylmax, kb-1:kb+nzlmax) 
  real*8 :: ainv4b, ajg, ajmx, ajmy, ajpx, ajpy, gb2, gb4, expo, eps, &
            dxa, dya, dza, &
            wpiwcigb, cfront, &
            dbxdy, dbydx, dbzdx, dbzdy, &
            ba1, ba2, ba3, ba4, b2
  integer*8 :: i, j, k, l, ietb, ietn, ietj, ietg, itresis

  data eps /1.e-25/
  
  eta = 0. 
  if (ieta == 0) then  ! choice 1
    eta = resis

  else if (ieta == 1) then ! choice 2
    itresis=1000
    do k = kb, ke
      do j = jb, je
        do i = 1, nx2
          ! no z-dependence
          eta(i,j,k) = resis/( cosh((real(i)-0.5*real(nx2+1))/real(netax))  &
                              *cosh((real(j)-0.5*real(ny2+1))/real(netay)) )
        enddo
      enddo
    enddo
    eta = eta*exp(-float(it)/float(itresis)) ! decays with time

  else if ((ieta == 2).or.(ieta == 5) ) then ! choice 3
    ! use gradient of |B|, B, and n
    ! good combination is (grad |B|)**4 / (n**4 * B**2)
    ! set powers, here; and adjust cfront factor accordingly
    ! e.g., cfront = 1e-3 for ietg=4, ietn=4, and ietb=2
    ! this also uses an upper limit on eta, set desired value at top
    cfront = 8.0e-3
    ietg = 4; ietn = 4; ietb = 2
    ainv4b = 1.0/ ( float(4**ietb) )
    wpiwcigb = wpiwci**(ietg -ietb)
    do k = kb, ke
      do j = jb, je
        do i = 2, nx1
          if (den(i,j,k) .le. denmin) then
            eta(i,j,k) = 0.0
          else
            ba1 = ( bx(i+1,j  ,k)      &
                    +bx(i+1,j+1,k) )**2 &
                  +( by(i+1,j  ,k)      &
                    +by(i+1,j+1,k) )**2 &
                  +( bz(i+1,j  ,k)      &
                    +bz(i+1,j+1,k) )**2
            ba2 = ( bx(i  ,j  ,k)      &
                    +bx(i  ,j+1,k) )**2 &
                  +( by(i  ,j  ,k)      &
                    +by(i  ,j+1,k) )**2 &
                  +( bz(i  ,j  ,k)      &
                    +bz(i  ,j+1,k) )**2
            ba3 = ( bx(i  ,j+1,k)      &
                    +bx(i+1,j+1,k) )**2 &
                  +( by(i  ,j+1,k)      &
                    +by(i+1,j+1,k) )**2 &
                  +( bz(i  ,j+1,k)      &
                    +bz(i+1,j+1,k) )**2
            ba4 = ( bx(i  ,j  ,k)      &
                    +bx(i+1,j  ,k) )**2 &
                  +( by(i  ,j  ,k)      &
                    +by(i+1,j  ,k) )**2 &
                  +( bz(i  ,j  ,k)      &
                    +bz(i+1,j  ,k) )**2

            ! Uniform mesh - Same as is in version 5.0
            ! gb2 = 1.0 *(  ((sqrt(ba1) -sqrt(ba2))/ hx)**2    &
            !             +((sqrt(ba3) -sqrt(ba4))/ hy)**2  )

            ! Nonuniform mesh
            gb2 = 1.0 *(  ((sqrt(ba1) -sqrt(ba2))/ meshX%dxc(i  ))**2    &
                          +((sqrt(ba3) -sqrt(ba4))/ meshY%dxc(j+1))**2  )

            gb4 = gb2**(ietg/2)
            b2 = ainv4b * (sqrt(ba1) +sqrt(ba2))**ietb
            eta(i,j,k) = cfront *wpiwcigb *resis *gb4/         &
                          (b2 *den(i,j,k)**ietn)
            eta(i,j,k) = min(eta(i,j,k), etamax)
            if (eta(i,j,k) .lt. etamin) eta(i,j,k) = 0.0
          endif
        enddo
      enddo
    enddo

  else if (ieta == 3) then
    ! use arbitrary power of j
    ! using the sum over all three components
    ! uses arbitrary power ietj
    ! adjust cfront factor accordingly
    ! e.g., cfront = 1.0
    ! this also uses an upper limit on eta, set desired value at top
    ietj = 4
    cfront = 2.0 *1.0e3**ietj

    ! Uniform mesh - Same as is in version 5.0
    ! dxa=1.0/(2.*hx)
    ! dya=1.0/(2.*hy)
    do k = kb, ke
      do j = jb, je
        ! Nonuniform mesh
        dya = 1.0/(2.0*meshY%dxc(j+1))
        do i = 2, nx1
          ! Nonuniform mesh
          dxa = 1.0/(2.0*meshX%dxc(i  ))

          dbzdy= bz(i+1,j+1,k)+bz(i  ,j+1,k) &
                -bz(i+1,j  ,k)-bz(i  ,j  ,k)
          dbzdx= bz(i+1,j+1,k)+bz(i+1,j  ,k) &
                -bz(i  ,j+1,k)-bz(i  ,j  ,k)
          dbydx= by(i+1,j+1,k)+by(i+1,j  ,k) &
                -by(i  ,j+1,k)-by(i  ,j  ,k)
          dbxdy= bx(i+1,j+1,k)+bx(i  ,j+1,k) &
                -bx(i+1,j  ,k)-bx(i  ,j  ,k)
          ajl(i,j,k) = (dya*dbzdy          )**2   &
                      +(         -dxa*dbzdx)**2 &
                      +(dxa*dbydx-dya*dbxdy)**2
        enddo
      enddo

      ! use 1/2 power, because of square root of squared components
      ! i.e., starting from |j| = sqrt(jx**2 +jy**2 +jz**2)
      expo = 0.5 *float(ietj)
      do j = jb, je
        do i = 2, nx1
          eta(i,j,k) = cfront *resis *ajl(i,j,k)**expo
          eta(i,j,k) = min(etamax, eta(i,j,k))
          if (eta(i,j,k) .lt. etamin) eta(i,j,k) = 0.0
        enddo
      enddo

    enddo

  else if (ieta==6) then ! exponential increase at edge of z
    do k = kb, ke
      do j = jb, je
        do i = 1, nx2
          ! if (k<196) then
          !   eta(i,j,k) = 0.
          ! else
          !   eta(i,j,k) = resis*(exp((real(k)-196.)/14.)-1.)
          ! endif 

          ! if (k.le.eta_zs) then
          !   eta(i,j,k) = resis*cos(pi*real(k)/(2.*real(eta_zs)))
          ! else if (k.ge.(nz-eta_zs)) then
          !   eta(i,j,k) = resis*cos(pi*(real(k)-nz)/(2.*real(eta_zs)))
          ! else
          !   eta(i,j,k) = 0.
          ! endif 
          
          if (k.le.eta_zs) then
            eta(i,j,k) = (resis/2)*(1+cos(pi*real(k)/real(eta_zs)))
          else if (k.ge.(nz-eta_zs)) then
            eta(i,j,k) = (resis/2)*(1-cos(pi*(real(k)-nz+eta_zs)/real(eta_zs)))
          else
            eta(i,j,k) = 0.
          endif 
        enddo 
      enddo 
    enddo 

  else if ( (ieta .gt. 6).or.(ieta .lt. 0) ) then
    call ERROR_ABORT('Currently eta_calc only accepts ieta = 0 ~ 6')

  endif

  if (ieta==4 .or. ieta==5) then
    ! use 2nd gradient of j
    ! this routine calculates something related to the
    ! logarithmic derivative of grad j...
    ! seems to be effective at inhibiting whistlers
    ! protect against division by zero with "eps"

    ! using the sum over all three components

    ! uses arbitrary power ietj

    ! adjust cfront factor accordingly
    ! e.g., cfront = 1.0
    ! this also uses an upper limit on eta, set desired value at top
    cfront = 0.5
    ! when adding this resitivity to the gradient method (case 2),
    ! lower the contribution from this method:
    if (ieta == 5) cfront = cfront/ 5.0
    ietj = 4

    ! Uniform mesh - Same as in version 5.0
    ! dxa=1.0/(2.*hx)
    ! dya=1.0/(2.*hy)

    do l = 1, 3
      do k = kb, ke
        do j = jb, je
          ! Nonuniform mesh
          dya=1.0/(2.0*meshY%dxc(j+1))
          do i=2,nx1
            ! Nonuniform mesh
            dxa=1.0/(2.0*meshX%dxc(i  ))
            if (l == 1) then
              dbzdy= bz(i+1,j+1,k)+bz(i  ,j+1,k) &
                    -bz(i+1,j  ,k)-bz(i  ,j  ,k)
              ajl(i,j,k) =  dya*dbzdy
              if (ieta == 4)  eta(i,j,k) = 0.0
            else if (l == 2) then
              dbzdx= bz(i+1,j+1,k)+bz(i+1,j  ,k) &
                    -bz(i  ,j+1,k)-bz(i  ,j  ,k)
              ajl(i,j,k) = -dxa*dbzdx
            else if (l == 3) then
              dbydx= by(i+1,j+1,k)+by(i+1,j  ,k) &
                    -by(i  ,j+1,k)-by(i  ,j  ,k)
              dbxdy= bx(i+1,j+1,k)+bx(i  ,j+1,k) &
                    -bx(i+1,j  ,k)-bx(i  ,j  ,k)
              ajl(i,j,k) =  dxa*dbydx -dya*dbxdy
            endif
          enddo
        enddo
      enddo

      ! set non-periodic boundary conditions
      do k=kb,ke
        do j=jb,je
          ajl(1  ,j,k) = ajl(2  ,j,k)
          ajl(nx2,j,k) = ajl(nx1,j,k)
        enddo
      enddo

      do k=kb,ke
        do i=1,nx2
          if (jb == 1) ajl(i,0   ,k) = ajl(i, 1,k)
          if (je == ny) ajl(i,ny1,k) = ajl(i,ny,k)
        enddo
      enddo

      do k=kb,ke
        do j=jb,je
          do i=2,nx1
            ! Uniform mesh - Same as is in version 5.0
            ajpx = (ajl(i+1,j,k) -ajl(i,  j,k))/ hx
            ajmx = (ajl(i,  j,k) -ajl(i-1,j,k))/ hx
            ajpy = (ajl(i,j+1,k) -ajl(i,  j,k))/ hy
            ajmy = (ajl(i,j  ,k) -ajl(i,j-1,k))/ hy
            ajg = abs( (ajpx - ajmx)/ hx/                  &
                  (abs(ajpx) +abs(ajmx) +eps) )**ietj &
                  *abs( (ajpy - ajmy)/ hy/                  &
                  (abs(ajpy) +abs(ajmy) +eps) )**ietj

            ! Nonuniform mesh
            ! ajpx = (ajl(i+1,j,k) -ajl(i,  j,k))/ meshX%dxn(i+1)
            ! ajmx = (ajl(i,  j,k) -ajl(i-1,j,k))/ meshX%dxn(i+1)
            ! ajpy = (ajl(i,j+1,k) -ajl(i,  j,k))/ meshY%dxn(j+2)
            ! ajmy = (ajl(i,j  ,k) -ajl(i,j-1,k))/ meshY%dxn(j+2)
            ! ajg = abs( (ajpx - ajmx)/ meshX%dxc(i  )/                  &
            !     (abs(ajpx) +abs(ajmx) +eps) )**ietj                   &
            !     *abs( (ajpy - ajmy)/ meshY%dxc(j+1)/                  &
            !     (abs(ajpy) +abs(ajmy) +eps) )**ietj
            eta(i,j,k) = eta(i,j,k) +cfront *resis *ajg
          enddo
        enddo
      enddo

    enddo

    do k = kb, ke
      do j = jb, je
        do i = 2, nx1
          eta(i,j,k) = min(etamax, eta(i,j,k))
          if (eta(i,j,k) .lt. etamin) eta(i,j,k) = 0.0
        enddo
      enddo
    enddo

  endif

end subroutine eta_calc

