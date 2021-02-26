!=======================================================================
!=======================================================================
      subroutine inject_wave
      use parameter_mod
      use MESH2D
      implicit none

      integer*8:: ibp1,ibp2,nptot_max,i,remake,field_subcycle
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi     &
                        ,x_p,y_p,z_p,x_p_logical,y_p_logical        &
                        ,z_p_logical,r_c,q_p,dtsav
 
      integer*8 ip,ipb1,ipb2,is,ixe,iye,ize,j,k,l, iixe,iiye, iize
      double precision vxa,vya,vza,vmag,th,ranval(4)

      real(kind=8) :: x_pos,y_pos,z_pos, B0, VA
      real(kind=8) :: bx_,by_,bz_, ex_,ey_,ez_

      real(kind=8) :: kx,ky,kz,kxmin,kymin,kzmin,dvx_,dvy_,dvz_,sin_factor
      real(kind=8) :: loaded_percentage, print_percentage

      double precision:: w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
      double precision:: vix1,viy1,viz1
      double precision:: vix2,viy2,viz2
      double precision:: vix3,viy3,viz3
      double precision:: vix4,viy4,viz4
      double precision:: vix5,viy5,viz5
      double precision:: vix6,viy6,viz6
      double precision:: vix7,viy7,viz7
      double precision:: vix8,viy8,viz8
      integer ixep1,iyep1,izep1


      dtxi = one/meshX%dt
      dtyi = one/meshY%dt
      dtzi = one/meshZ%dt
 
      !VR: initialize wave parameters

      ! dB_B0 = 0.1                ! RMS amplitude of the pertubation [B0=RMS(B)]
      B0 = one/wpiwci
      !VR Alfven speed
      VA = one/wpiwci      

      kxmin = two*pi/xmax
      kymin = two*pi/ymax
      kzmin = two*pi/zmax
    
      kx = zero
      ky = zero
      kz = num_cycles * kzmin
      
      !VR: end wave parameters

      if (myid==0) write(6,*) "Injecting waves"
      vixo=vix
      viyo=viy
      vizo=viz
      do k=kb-1,ke+1
         z_pos = meshZ%xc(k+1)
         do j=jb-1,je+1  
            y_pos = meshY%xc(j+1)
            do i=1,nx2
               x_pos = meshX%xc(i)   !VR this is not a typo. For some reason, x has different indexing compared to y and z (!!!)               
!! single Alfven wave
               bx_ =  dB_B0*B0*sin(kz*z_pos)
               by_ = -dB_B0*B0*cos(kz*z_pos)
               bz_ = zero
               ex_ = zero
               ey_ = zero
               ez_ = zero
              
               dvx_ = -VA*bx_/B0 
               dvy_ = -VA*by_/B0 
               dvz_ = zero

              !
              bx(i,j,k) =bx(i,j,k)+ bx_
              by(i,j,k) =by(i,j,k)+ by_
              bz(i,j,k) =bz(i,j,k)+ bz_
              ex(i,j,k) =ex(i,j,k)+ ex_
              ey(i,j,k) =ey(i,j,k)+ ey_
              ez(i,j,k) =ez(i,j,k)+ ez_

              ! use vix to temporary store values of V on the grid
              vix(i,j,k) = dvx_
              viy(i,j,k) = dvy_
              viz(i,j,k) = dvz_
              
            enddo
         enddo
      enddo
   
    do is =1, nspec 
      DO iize = kb-1,ke
        DO iiye = jb-1,je
          DO iixe = 1, nx1
            l=iphead(iixe,iiye,iize,is)

            DO WHILE (l.NE.0)
  !         Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
              rye=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
              rze=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
              ixe=rxe
              iye=rye
              ize=rze
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0
  !
              fxe=rxe-ixe
              fye=rye-iye
              fze=rze-ize
              ixep1 = ixe+1
              iyep1 = iye+1
              izep1 = ize+1

              w1e=(1.-fxe)*(1.-fye)*(1.-fze)
              w2e=fxe*(1.-fye)*(1.-fze)
              w3e=(1.-fxe)*fye*(1.-fze)
              w4e=fxe*fye*(1.-fze)
              w5e=(1.-fxe)*(1.-fye)*fze
              w6e=fxe*(1.-fye)*fze
              w7e=(1.-fxe)*fye*fze
              w8e=fxe*fye*fze
            
              vix1=vix(ixe  ,iye  ,ize  )
              vix2=vix(ixep1,iye  ,ize  )
              vix3=vix(ixe  ,iyep1,ize  )
              vix4=vix(ixep1,iyep1,ize  )
              vix5=vix(ixe  ,iye  ,izep1)
              vix6=vix(ixep1,iye  ,izep1)
              vix7=vix(ixe  ,iyep1,izep1)
              vix8=vix(ixep1,iyep1,izep1)
              viy1=viy(ixe  ,iye  ,ize  )
              viy2=viy(ixep1,iye  ,ize  )
              viy3=viy(ixe  ,iyep1,ize  )
              viy4=viy(ixep1,iyep1,ize  )
              viy5=viy(ixe  ,iye  ,izep1)
              viy6=viy(ixep1,iye  ,izep1)
              viy7=viy(ixe  ,iyep1,izep1)
              viy8=viy(ixep1,iyep1,izep1)
              viz1=viz(ixe  ,iye  ,ize  )
              viz2=viz(ixep1,iye  ,ize  )
              viz3=viz(ixe  ,iyep1,ize  )
              viz4=viz(ixep1,iyep1,ize  )
              viz5=viz(ixe  ,iye  ,izep1)
              viz6=viz(ixep1,iye  ,izep1)
              viz7=viz(ixe  ,iyep1,izep1)
              viz8=viz(ixep1,iyep1,izep1)
              
              dvx_=w1e*vix1+w2e*vix2+w3e*vix3+w4e*vix4   &
                   +w5e*vix5+w6e*vix6+w7e*vix7+w8e*vix8  
              dvy_=w1e*viy1+w2e*viy2+w3e*viy3+w4e*viy4   &
                   +w5e*viy5+w6e*viy6+w7e*viy7+w8e*viy8  
              dvz_=w1e*viz1+w2e*viz2+w3e*viz3+w4e*viz4   &
                   +w5e*viz5+w6e*viz6+w7e*viz7+w8e*viz8  
                        
            !interpolate V at the particle position from pre-computed values at the grid


              vx(l)=vx(l)+dvx_
              vy(l)=vy(l)+dvy_
              vz(l)=vz(l)+dvz_

              l=link(l)
   
            enddo ! while
          enddo
        enddo
      enddo
    enddo ! is
      vix=vixo+vix
      viy=viyo+viy
      viz=vizo+viz
     
      
     if (ndim /= 1) then
        call xrealbcc(ex,1_8,nx,ny,nz)
        call xrealbcc(ey,1_8,nx,ny,nz)
        call xrealbcc(ez,1_8,nx,ny,nz)
     else
        call xrealbcc_pack_e_2d(ex,ey,ez,1_8,nx,ny,nz)
     endif
     
 
 
    return
  end subroutine inject_wave
!
!***********************************************************************
      subroutine kick
      use parameter_mod
      use MESH2D
      implicit none

      integer*8:: i,j,k

      real(kind=8) :: x_pos,y_pos,z_pos, B0, VA
      real(kind=8) :: bx_,by_,bz_, ex_,ey_,ez_

      real(kind=8) :: kx,ky,kz,kxmin,kymin,kzmin

!!!!  Alfvenic perturbation with deltaB in the x direction
!!!!  i,j,k are wave numbers in x,y,z
#define DBX_1(k,j,phi) (dB_B0*B0*cos((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define DEY_1(k,j,phi) (-dB_B0*(k/abs(k))*VA*B0*cos((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
!!!! These give velocity & current consistent with Alfven wave
#define DUX_1(k,j,phi) (-dB_B0*(k/abs(k))*VA*cos((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define DJY_1(k,j,phi) (-dB_B0*B0*(k)*kzmin*sin((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define DJZ_1(k,j,phi) (dB_B0*B0*(j)*kymin*sin((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
!
#define BX_PERT_1 DBX_1(1,1,0) + DBX_1(1,2,1.5) + DBX_1(-2,3,3.9)   
#define EY_PERT_1 DEY_1(1,1,0) + DEY_1(1,2,1.5) + DEY_1(-2,3,3.9)  
#define UX_PERT_1 DUX_1(1,1,0) + DUX_1(1,2,1.5) + DUX_1(-2,3,3.9)
#define JY_PERT_1 DJY_1(1,1,0) + DJY_1(1,2,1.5) + DJY_1(-2,3,3.9)
#define JZ_PERT_1 DJZ_1(1,1,0) + DJZ_1(1,2,1.5) + DJZ_1(-2,3,3.9)

!!!!  Alfvenic perturbation with deltaB in the y direction
!!!!  works only for a pair plasma
#define DBY_2(k,i,phi) (dB_B0*B0*cos((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define DEX_2(k,i,phi) (dB_B0*(k/abs(k))*VA*B0*cos((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
!!!! These give velocity & current consistent with Alfven wave
#define DUY_2(k,i,phi) (-dB_B0*(k/abs(k))*VA*cos((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define DJX_2(k,i,phi) (dB_B0*B0*(k)*kzmin*sin((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define DJZ_2(k,i,phi) (-dB_B0*B0*(i)*kxmin*sin((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
!
#define BY_PERT_2 DBY_2(-1,1,0.4) + DBY_2(-1,-2,2.56) + DBY_2(2,-3,4.19)   
#define EX_PERT_2 DEX_2(-1,1,0.4) + DEX_2(-1,-2,2.56) + DEX_2(2,-3,4.19) 
#define UY_PERT_2 DUY_2(-1,1,0.4) + DUY_2(-1,-2,2.56) + DUY_2(2,-3,4.19)
#define JX_PERT_2 DJX_2(-1,1,0.4) + DJX_2(-1,-2,2.56) + DJX_2(2,-3,4.19)
#define JZ_PERT_2 DJZ_2(-1,1,0.4) + DJZ_2(-1,-2,2.56) + DJZ_2(2,-3,4.19)
 
      !VR: initialize wave parameters

      !dB_B0 = 1.2e-3                ! RMS amplitude of the pertubation [B0=RMS(B)]
      ! dB_B0 = 4e-3   
      B0 = one/wpiwci
      !VR Alfven speed
      VA = one/wpiwci      

      kxmin = two*pi/xmax
      kymin = two*pi/ymax
      kzmin = two*pi/zmax
    
      kx = zero
      ky = zero
      kz = num_cycles * kzmin
      
      !VR: end wave parameters

      !if (myid==0) write(6,*) "Kicking"
      do k=kb-1,ke+1
        z_pos = meshZ%xc(k+1)
        do j=jb-1,je+1  
          y_pos = meshY%xc(j+1)
          do i=1,nx2
            x_pos = meshX%xc(i)   !VR this is not a typo. For some reason, x has different indexing compared to y and z (!!!)               
!! single Alfven wave
            !bx_ =  dB_B0*B0*sin(kz*z_pos)
            !by_ = -dB_B0*B0*cos(kz*z_pos)
            !bz_ = zero
            !ex_ = zero
            !ey_ = zero
            !ez_ = zero
! multiple waves
            bx_ = BX_PERT_1
            by_ = BY_PERT_2
              
            bx(i,j,k) =bx(i,j,k)+ bx_
            by(i,j,k) =by(i,j,k)+ by_
            !bz(i,j,k) =bz(i,j,k)+ bz_
            !ex(i,j,k) =ex(i,j,k)+ ex_
            !ey(i,j,k) =ey(i,j,k)+ ey_
            !ez(i,j,k) =ez(i,j,k)+ ez_

          enddo
        enddo
      enddo
    return
  end subroutine kick
