! ***** file = channel_flow_3d.F90*******************************************
! *     turbulant flow in a square cavity                                   *
! *     *  2nd order central FDM in space                                   *
! *     *  2nd order adams-bashforth scheme for convection and diffusion    *
! *     *  1st order backword Euler scheme for continuace equation          *
! *        and pressure gradient                                            *
! *     *  fractional step method                                                      *
! ***************************************************************************

module param
  implicit none
  integer, parameter :: nx = 64,ny = 64, nz = 64
  real(8), parameter :: dxp = 18.d0, dzp = 9.d0                    !non-dimension grid gap
  integer, parameter :: iskip = 1000,istop = 1e7
  real(8), parameter :: Re_t = 300.d0
  real(8), parameter :: dt = 1.d-4                                 !non-dimension time, for wall-unit, multiply Re_t
  real(8), parameter :: T_osc=100.d0                               !Wall oscilation frequency  
end module param

module input
  use param
  implicit none
  real(8), dimension(nx/2,ny,nz/2) :: UPI, VPI, WPI, PPI
end module input

module time_step
  implicit none
  integer :: istep
  real(8) :: tstep
end module time_step

module mesh
  use param
  implicit none

  real(8) :: dx, dxi, dxi2, dz, dzi, dzi2
  real(8) :: Hx,HY,Hz
  integer, dimension(0:nx+1) :: ip,im
  integer, dimension(0:nz+1) :: kp,km
  real(8), dimension(nx) :: xp
  real(8), dimension(nz) :: zp
  real(8), dimension(0:ny) :: dyc,dyci, dyci2
  real(8), dimension(0:ny+1) :: dy,dyi, dyi2
  real(8), dimension(0:ny) :: yv
  real(8), dimension(0:ny+1) :: yp
  real(8), parameter :: ALG = 9.5d-1
  real(8) :: AT, ETA
end module mesh

module velocity
  use param
  implicit none
  real(8), dimension(nx,ny,nz) :: u,w,CX,CZ
  real(8), dimension(nx,0:ny,nz) :: CY,v
  real(8), dimension(nx,ny,nz) :: up,uf,vp,vf,wp,wf
  real(8) :: VX, VY, VZ
  real(8) :: a,b                                         !parameters for explicit method 
end module velocity

module SOR
  use param
  implicit none
  integer :: itrp, itrp_end
  real(8) :: sumq, sumr, resi, B0, errp
  real(8), dimension(nx,ny,nz) :: Q,p

end module SOR

module average_velocity
  use param
  implicit none
  real(8),dimension(ny) :: u_ave, v_ave, w_ave, u_sta,v_sta,w_sta, u_rms, v_rms, w_rms
  real(8),dimension(ny) :: u_ave2, v_ave2, w_ave2, u_sta2,v_sta2, w_sta2
  real(8),dimension(ny) :: sumu, sumv, sumw, sumu2, sumv2, sumw2
end module average_velocity

module lagrange_marker
  use param
  implicit none
  real(8) :: lag_mark_x
  integer :: pass_count
end module lagrange_marker

module vorticity
  use param
  implicit none
  real(8), dimension(nx,ny,nz)::udx,vdx,wdx,udy,vdy,wdy,udz,vdz,wdz,qq
end module vorticity

module save
  use param
  implicit none
  real,dimension(nx,ny,nz) :: xp_real,yp_real,zp_real
  real,dimension(nx,ny,nz) :: u_real, v_real, w_real, p_real, qq_real

end module

!-------------------------------------------------------------------------
!  Main
!-------------------------------------------------------------------------

program main
  use param
  use time_step
  use mesh
  use velocity
  use SOR
  use lagrange_marker
  implicit none
!  real(8) :: etk, energb                                                   !future area

  call Mesh_Generation
  call InitializeFlowField
  write(6,1100)

  pass_count = 0
  lag_mark_x = 0.d0

  do istep=1,istop

     tstep = dt * dble(istep)

     call RHSforVelocity
     call VelocityPrediction
     call PoissonRHS
     call solve_SOR
     call Correction
     call q_value

     if(mod(istep,iskip) == 0) then
        write(6,1000) istep,tstep,itrp,errp,sumq,lag_mark_x, pass_count !,divmax  
        call SaveData
     endif


  enddo
  write(*,*) "end"

1000 format(i9,f8.3,i9,f12.9,f20.2,f8.3,i3)
1100 format(2x,4hstep,7x,1ht,6x,6h  itrp,6x,4herrp,10x, 4hsumq,6x,&
            9hlag-point, 10x, 10hpass_count)

  stop

  return
end program main

!-------------------------------------------------------------------------
!  subroutine. Mesh
!-------------------------------------------------------------------------
subroutine Mesh_Generation
  use mesh

  call mesh_xz
  call mesh_y

  write(6,1200) nx, ny, nz, dx, dy(32),dy(1),dz
1200 format(&
          /' -------------------Mesh condition ----------------------'&
          /' Numbers of grid points             nx     = ',i12        &
          /'                                    ny     = ',i12        &
          /'                                    nz     = ',i12        &
          /' Grid resolution                    dx     = ',f12.3      &
          /'                                    dy_max = ',f12.3      &
          /'                                    dy_min = ',f12.3      &
          /'                                    dz     = ',f12.3      &
          /' --------------------------------------------------------')


  return
end subroutine Mesh_Generation

!-------------------------------------------------------------------------
!  generate mesh for direction X,Z
!-------------------------------------------------------------------------
  subroutine Mesh_xz
    use param
    use mesh
    integer :: i,k

    dx = dxp / Re_t
    dz = dzp / Re_t
    dxi = 1.d0 / dx
    dxi2 = 1.d0 /dx**2
    dzi = 1.d0 / dz
    dzi2 = 1.d0 /dz**2

!   periodic condition
    do i=1,nx
      ip(i)=i+1
      im(i)=i-1
      if(ip(i)>nx) then
        ip(i) = ip(i)-nx
      endif
      if(im(i)< 1) then
        im(i) = im(i)+nx
      endif
    end do
    do k=1, nz
      kp(k)=k+1
      km(k)=k-1
      if(kp(k)>nz) then
        kp(k) = kp(k)-nz
      endif
      if(km(k)<1) then
        km(k) = km(k)+nz
      endif
    end do

    xp(nx) = 0.d0
    zp(nz) = 0.d0

    do i=1,nx
      xp(i) = xp(im(i))+dx
    end do

    do k=1,nz
      zp(k) = zp(km(k))+dz
    end do

    open(10,file='Mesh_Data/xp.dat',status='replace',form='formatted')
    do i=1,nx
      write(10,*) xp(i)
    end do
    close(10)

    open(10,file='Mesh_Data/zp.dat',status='replace',form='formatted')
    do k=1,nz
      write(10,*) zp(k)
    end do
    close(10)

    return
  end subroutine Mesh_xz

!-------------------------------------------------------------------------
!  generate mesh for direction Y
!-------parametes placement image------------
!        yp = center of mesh
!        yv = edge of Mesh
!       _ _ _ _ _ _ _  _ _ _ _ _ _ _
!      ｜ 　         ｜             ｜
!      ｜   yp(j) 　 ｜   yp(j+1)　 ｜
!      ｜     ・←-----------→・     ｜
!      ｜        　 dyc(j)      　  ｜
!      ｜_ _ _ _ _ _ ｜_ _ _ _ _ _ _｜
!    yv(j-1)       yv(j)          yv(j+1)
!      ｜←----------→｜←-----------→｜
!      ｜    dy(j)　 ｜   dy(j+1)   ｜
!--------------------------------------------------------------------------
  subroutine Mesh_y
    use param
    use mesh
    implicit none
    integer :: j

    AT = dlog((1.d0+ALG)/(1.d0-ALG))/2.d0
    yv(0) = 0.d0
    yv(ny) = 1.d0                               !the formula doesn't give exact 1 so 
    do j=1,ny-1
      ETA = AT*(-1.d0+2.d0*dble(j)/dble(ny))
      yv(j) = (dtanh(ETA)/ALG+1.d0)/2.d0
    end do

    do j=1,ny
      ETA = AT*(-1.d0+2.d0*(dble(j)-0.5d0)/dble(ny))
      yp(j) = (dtanh(ETA)/ALG+1.d0)/2.d0
    end do

    open(10,file='Mesh_Data/yp.dat',status='replace',form='formatted')
    do j=1,ny
      write(10,*) yp(j)
    end do
    close(10)

    !half mesh for y direction in outer
    yp(0) = 2.d0*yv(0)-yp(1)
    yp(ny+1) = 2.d0*yv(ny)-yp(ny)

    do j=1, ny
      dy(j) = -yv(j-1)+yv(j)
      dyi(j) = 1.d0/dy(j)
      dyi2(j) = 1.d0/dy(j)**2
    end do

    do j=1,ny-1
      dyc(j)= -yp(j)+yp(j+1)
      dyci(j) = 1.d0/dyc(j)
      dyci2(j) = 1.d0/dyc(j)**2
    end do

!　Noumann B.C at the wall (dp/dy=0) for pressure＠Convection term
    dyci(0) = 0.d0
    dyci(ny) = 0.d0

!  non-source＠Vicous term
    dyi(0) = 0.d0
    dyi(ny+1) = 0.d0

  return
end subroutine Mesh_y

!-------------------------------------------------------------------------
!  subroutine. Initialize FlowField
!-------------------------------------------------------------------------
subroutine InitializeFlowField
  use input
  use param
  use mesh
  use velocity
  use SOR
  implicit none
  integer :: i,j,k

  open(10,file="initia.bin", form = "unformatted", status = "unknown", convert = 'BIG_ENDIAN')
  read(10)UPI
  read(10)VPI
  read(10)WPI
  read(10)PPI
  close(10)

  do k=1,nz/2
    do j=1,nz
      do i=1,nx/2
        u(i,j,k) = UPI(i,j,k)
        v(i,j,k) = VPI(i,j,k)
        w(i,j,k) = WPI(i,j,k)
        p(i,j,k) = PPI(i,j,k)
      end do
    end do
  end do

  do k= nz/2+1,nz
    do j=1,nz
      do i=nx/2+1, nx
        u(i,j,k) = UPI(i-nx/2,j,k-nz/2)
        v(i,j,k) = VPI(i-nx/2,j,k-nz/2)
        w(i,j,k) = WPI(i-nx/2,j,k-nz/2)
        p(i,j,k) = PPI(i-nx/2,j,k-nz/2)
      end do
    end do
  end do

  return
end subroutine InitializeFlowField

subroutine RHSforVelocity
  use param
  use mesh
  use velocity
  implicit none
  integer :: i,j,k

!-------------------------------------------------------------------------
!  Convection term
!-------------------------------------------------------------------------

! components for x
  do k=1,nz
    do j=1,ny
      do i=1,nx
        CX(i,j,k) = dy(j)*dz*(u(i,j,k)+u(im(i),j,k))*(u(i,j,k)-u(im(i),j,k))/4.d0
        CZ(i,j,k) = dx*dy(j)*(w(i,j,k)+w(ip(i),j,k))*(u(i,j,kp(k))-u(i,j,k))/4.d0
      end do
    end do
  end do

  do k=1,nz
    do j=0,ny
      do i=1,nx
        if(j>0 .and. J<ny) then
          CY(i,j,k) = dx*dz*(v(i,j,k)+v(ip(i),j,k))*(-u(i,j,k)+u(i,j+1,k))/4.d0
        !non-slip for wall
        else
          CY(i,j,k) = 0.d0
        end if
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny
      do i=1,nx
        uf(i,j,k) = (CX(i,j,k)+CX(ip(i),j,k)+CY(i,j-1,k)+CY(i,j,k)&
                    +CZ(i,j,km(k))+CZ(i,j,k))*dxi*dyi(j)*dzi
      end do
    end do
  end do

! components for y
  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        CX(i,j,k) = dz*(dy(j)*u(i,j,k)+dy(j+1)*u(i,j+1,k))*(-v(i,j,k)+v(ip(i),j,k))/4.d0
        CZ(i,j,k) = dx*(dy(j)*w(i,j,k)+dy(j+1)*w(i,j+1,k))*(-v(i,j,k)+v(i,j,kp(k)))/4.d0
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny
      do i=1,nx
        CY(i,j,k) = dx*dz*(v(i,j-1,k)+v(i,j,k))*(-v(i,j-1,k)+v(i,j,k))/4.d0
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        vf(i,j,k) = (CX(im(i),j,k)+CX(i,j,k)+CY(i,j+1,k)+CY(i,j,k)&
                    +CZ(i,j,km(k))+CZ(i,j,k))*dxi*dyci(j)*dzi
      end do
    end do
  end do

! components for z
  do k=1,nz
    do j=1,ny
      do i=1,nx
        CX(i,j,k) = dy(j)*dz*(u(i,j,k)+u(i,j,kp(k)))*(w(ip(i),j,k)-w(i,j,k))/4.d0
        CZ(i,j,k) = dx*dy(j)*(w(i,j,k)+w(i,j,km(k)))*(-w(i,j,km(k))+w(i,j,k))/4.d0
      end do
    end do
  end do

  do k=1,nz
    do j=0,ny
      do i=1,nx
        if(j>0 .and. J<ny) then
          CY(i,j,k) = dz*dx*(v(i,j,k)+v(i,j,kp(k)))*(-w(i,j,k)+w(i,j+1,k))/4.d0
        !non-slip on the wall for z
        else
          CY(i,j,k) = 0.d0
        end if
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny
      do i=1,nx
        wf(i,j,k) = (CX(i,j,k)+CX(im(i),j,k)+CY(i,j-1,k)+CY(i,j,k)&
                    +CZ(i,j,kp(k))+CZ(i,j,k))*dxi*dyi(j)*dzi
      end do
    end do
  end do

!-------------------------------------------------------------------------
!  Viscous term
!-------------------------------------------------------------------------

  do k=1,nz
    do j=1,ny
      do i=1,nx
        VX = dxi/Re_t*(-(-u(im(i),j,k)+u(i,j,k))*dxi+(-u(i,j,k)+u(ip(i),j,k))*dxi)
        VZ = dzi/Re_t*(-(-u(i,j,km(k))+u(i,j,k))*dzi+(-u(i,j,k)+u(i,j,kp(k)))*dzi)
        if(j>1 .and. j<ny) then
          VY = dyi(j)/Re_t*(-(-u(i,j-1,k)+u(i,j,k))*dyci(j-1)+(-u(i,j,k)+u(i,j+1,k))*dyci(j)) 
        elseif(j==1) then
          VY = dyi(j)/Re_t*((-u(i,j,k)*yp(2)/yp(1)+u(i,j+1,k)*yp(1)/yp(2))*dyci(j)&
              +(-u(i,j,k)+u(i,j+1,k))*dyci(j))
        else
          VY = dyi(j)/Re_t*(-(-u(i,j-1,k)+u(i,j,k))*dyci(j-1)&
              +(-u(i,j,k)*(1.d0-yp(ny-1))/(1.d0-yp(ny))+u(i,j-1,k)*(1.d0-yp(ny))/(1.d0-yp(ny-1)))*dyci(j-1))
        endif

        uf(i,j,k) = uf(i,j,k)+VX+VY+VZ
      end do
    end do
  end do

! v(i,ny,k)=0 because it's staggered grid
  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        VX = dxi/Re_t*(-(-v(im(i),j,k)+v(i,j,k))*dxi+(-v(i,j,k)+v(ip(i),j,k))*dxi)
        VZ = dzi/Re_t*(-(-v(i,j,km(k))+v(i,j,k))*dzi+(-v(i,j,k)+v(i,j,kp(k)))*dzi)
        VY = dyci(j)/Re_t*(-(-v(i,j-1,k)+v(i,j,k))*dyi(j)+(-v(i,j,k)+v(i,j+1,k))*dyi(j+1))

        vf(i,j,k) = vf(i,j,k)+VX+VY+VZ
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny
      do i=1,nx
        VX = dxi/Re_t*(-(-w(im(i),j,k)+w(i,j,k))*dxi+(-w(i,j,k)+w(ip(i),j,k))*dxi)
        VZ = dzi/Re_t*(-(-w(i,j,km(k))+w(i,j,k))*dzi+(-w(i,j,k)+w(i,j,kp(k)))*dzi)
        if(j>1 .and. j<ny) then
          VY = dyi(j)/Re_t*(-(-w(i,j-1,k)+w(i,j,k))*dyci(j-1)+(-w(i,j,k)+w(i,j+1,k))*dyci(j)) 
        elseif(j==1) then
          VY = dyi(j)/Re_t*((-w(i,j,k)*yp(2)/yp(1)+w(i,j+1,k)*yp(1)/yp(2))*dyci(j)&
              +(-w(i,j,k)+w(i,j+1,k))*dyci(j))
        else
          VY = dyi(j)/Re_t*(-(-w(i,j-1,k)+w(i,j,k))*dyci(j-1)&
              +(-w(i,j,k)*(1.d0-yp(ny-1))/(1.d0-yp(ny))+w(i,j-1,k)*(1.d0-yp(ny))/(1.d0-yp(ny-1)))*dyci(j-1))
        endif

        wf(i,j,k) = wf(i,j,k)+VX+VY+VZ
      end do
    end do
  end do

  return
end subroutine RHSforVelocity

subroutine VelocityPrediction
  use param
  use mesh
  use velocity
  use time_step
  implicit none
  integer :: i,j,k

  if(istep==1) then
    a = 0.d0
    b = 2.d0
  else
    a = -1.d0
    b = 3.d0
  endif

  do k=1,nz
    do j=1,ny
      do i=1,nx
        !ad  pressure gradient
        uf(i,j,k) = uf(i,j,k) + 2.d0

         u(i,j,k) = u(i,j,k)+dt*(a*up(i,j,k)+b*uf(i,j,k))/2.d0
         w(i,j,k) = w(i,j,k)+dt*(a*wp(i,j,k)+b*wf(i,j,k))/2.d0
        up(i,j,k) = uf(i,j,k)
        wp(i,j,k) = wf(i,j,k)
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        v(i,j,k) = v(i,j,k)+dt*(a*vp(i,j,k)+b*vf(i,j,k))/2.d0
        vp(i,j,k) = vf(i,j,k)
      end do
    end do
  end do

  return
end subroutine VelocityPrediction

subroutine PoissonRHS
  use param
  use mesh
  use Velocity
  use SOR
    implicit none
  integer :: i,j,k

! RHS for Poisson equation
  do k=1,nz
    do j=1,ny
      do i=1,nx
          Q(i,j,k) = 1.d0/dt*((-u(im(i),j,k)+u(i,j,k))*dxi&
                     +(-v(i,j-1,k)+v(i,j,k))*dyi(j)&
                     +(-w(i,j,km(k))+w(i,j,k))*dzi)
      end do
    end do
  end do

  return
end subroutine PoissonRHS

subroutine solve_SOR
  use param
  use mesh
  use Velocity
  use SOR
  implicit none
  integer :: i,j,k

  itrp_end = 40
  sumq = 0.d0

  do k=1,nz
    do j=1,ny
      do i=1,nx
        sumq = Q(i,j,k)**2 + sumq
      end do
    end do
  end do

  do itrp=1,itrp_end
    sumr = 0.d0
    do j=1,ny
      B0 = -2.d0*(dxi2+dzi2)-dyi(j)*(dyci(j)+dyci(j-1))
      do k=1,nz
        do i=1,nx
          if(j==1) then
            resi = dxi2*(p(im(i),j,k)+p(ip(i),j,k))&
                  +dyi(j)*dyci(j)*p(i,j+1,k)+B0*p(i,j,k)&
                  +dzi2*(p(i,j,km(k))+p(i,j,kp(k)))&
                  -Q(i,j,k)

          elseif (j == ny) then
            resi = dxi2*(p(im(i),j,k)+p(ip(i),j,k))&
                  +dyi(j)*dyci(j-1)*p(i,j-1,k)+B0*p(i,j,k)&
                  +dzi2*(p(i,j,km(k))+p(i,j,kp(k)))&
                  -Q(i,j,k)

          else
            resi = dxi2*(p(im(i),j,k)+p(ip(i),j,k))&
                  +dyi(j)*(dyci(j-1)*p(i,j-1,k)+dyci(j)*p(i,j+1,k))+B0*p(i,j,k)&
                  +dzi2*(p(i,j,km(k))+p(i,j,kp(k)))&
                  -Q(i,j,k)
          endif

          p(i,j,k) = p(i,j,k)-1.5d0*resi/B0
          sumr = sumr + resi**2

        end do
      end do
    end do

    errp = dsqrt(sumr/sumq)
  !  write(6,*) 'ITRP=',itrp,'   err=',errp
    if(errp < 1.d-5 .or. sumq == 0.d0) exit

  end do

  return
end subroutine solve_SOR

subroutine Correction
  use param
  use time_step
  use mesh
  use velocity
  use SOR
  use average_velocity
  use lagrange_marker
  implicit none
  integer i,j,k

  sumu = 0.d0
  sumv = 0.d0
  sumw = 0.d0
  sumu2 = 0.d0
  sumv2 = 0.d0
  sumw2 = 0.d0

  do k=1,nz
    do j=1,ny
      do i=1,nx
        u(i,j,k) = u(i,j,k) -dt*(-p(i,j,k)+p(ip(i),j,k))*dxi
        w(i,j,k) = w(i,j,k) -dt*(-p(i,j,k)+p(i,j,kp(k)))*dzi
        sumu(j) = sumu(j)+u(i,j,k)
        sumw(j) = sumw(j)+w(i,j,k)
        sumu2(j) = sumu2(j)+u(i,j,k)**2
        sumw2(j) = sumw2(j)+w(i,j,k)**2
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        v(i,j,k) = v(i,j,k) -dt*(-p(i,j,k)*dyci(j)+p(i,j+1,k)*dyci(j))
        sumv(j) = sumv(j)+v(i,j,k)
        sumv2(j) = sumv2(j)+v(i,j,k)**2
      end do
    end do
  end do

  !area average
  do j=1,ny
       u_ave(j) =  sumu(j)/dble(nx*nz)
       v_ave(j) =  sumv(j)/dble(nx*nz)
       w_ave(j) =  sumw(j)/dble(nx*nz)
      u_ave2(j) = sumu2(j)/dble(nx*nz)
      v_ave2(j) = sumv2(j)/dble(nx*nz)
      w_ave2(j) = sumw2(j)/dble(nx*nz)
  end do

 !time average
  if(istep>900000)then
    do j=1,ny
      u_sta(j) = (u_ave(j)+u_sta(j)*dble(istep-900001))/dble(istep-900000)
      v_sta(j) = (v_ave(j)+v_sta(j)*dble(istep-900001))/dble(istep-900000)
      w_sta(j) = (w_ave(j)+w_sta(j)*dble(istep-900001))/dble(istep-900000)
      u_sta2(j) = (u_ave2(j)+u_sta2(j)*dble(istep-900001))/dble(istep-900000)
      v_sta2(j) = (v_ave2(j)+v_sta2(j)*dble(istep-900001))/dble(istep-900000)
      w_sta2(j) = (w_ave2(j)+w_sta2(j)*dble(istep-900001))/dble(istep-900000)

      !RMS
      u_rms(j) = sqrt(u_sta2(j)-u_sta(j)**2)
      v_rms(j) = sqrt(v_sta2(j)-v_sta(j)**2)
      w_rms(j) = sqrt(w_sta2(j)-w_sta(j)**2)
    end do
  endif

 !! set marker and cound how many times the marker goes through the computation region
  lag_mark_x = lag_mark_x + u_ave(ny/2)*dt
  if(lag_mark_x > 3.84)then
    lag_mark_x = 0.d0
    pass_count = pass_count+1
  endif
  !do k=1,nz
  !  do j=1,ny
  !    do i=1,nx
  !      div = (-u(im(i),j,k)+u(i,j,k))*dxi*(-v(i,j-1,k)+v(i,j,k))*dyi(j)*
  !    end do
  !  end do
  !end do

  return
end subroutine Correction

subroutine q_value
  use param
  use mesh
  use velocity
  use vorticity
  implicit none
  integer i,j,k

  do k=1,nz
    do j=1,ny
      do i=1,nx
        udx(i,j,k) = (u(i,j,k)-u(im(i),j,k))/dx
        vdx(i,j,k) = (u(i,j,k)-u(im(i),j,k))/dx
        wdx(i,j,k) = (w(i,j,k)-w(im(i),j,k))/dx
        udz(i,j,k) = (u(i,j,k)-u(i,j,km(k)))/dz
        vdz(i,j,k) = (v(i,j,k)-v(i,j,km(k)))/dz
        wdz(i,j,k) = (w(i,j,k)-w(i,j,km(k)))/dz
      end do
    end do
  end do

  do k=1,nz
    do i=1,nx
      udy(i,1,k) = 0.d0
      vdy(i,1,k) = 0.d0
      wdy(i,1,k) = 0.d0
    end do
  end do

  do k=1,nz
    do j=2,ny
      do i=1,nx
        udy(i,j,k) = (u(i,j,k)-u(i,j-1,k))/dy(j)
        vdy(i,j,k) = (v(i,j,k)-v(i,j-1,k))/dy(j)
        wdy(i,j,k) = (w(i,j,k)-w(i,j-1,k))/dy(j)
      end do
    end do
  end do

  do k=1,nz
    do j=1,ny
      do i=1,nx
        qq(i,j,k) = -0.5*(udx(i,j,k)**2+vdy(i,j,k)**2+wdz(i,j,k)**2)&
                   -(udy(i,j,k)*vdx(i,j,k)&
                   +vdz(i,j,k)*wdy(i,j,k)+&
                   wdx(i,j,k)*udz(i,j,k))
      end do
    end do
  end do

  return 
end subroutine q_value

subroutine SaveData
  use param
  use mesh
  use velocity
  use time_step
  use average_velocity
  use SOR
  use vorticity
  use save
  implicit none
  integer i,j,k
  character :: fnum*7,fname*40
! --------------------------------------------------------------------
  write(fnum,'(i7.7)') istep

 do k=1,nz
   do j=1,ny
     do i=1,nx
       u_real(i,j,k) = real(u(i,j,k))
       v_real(i,j,k) = real(v(i,j,k))
       w_real(i,j,k) = real(w(i,j,k))
       p_real(i,j,k) = real(p(i,j,k))
       qq_real(i,j,k) = real(qq(i,j,k))
       xp_real(i,j,k) = real(xp(i))
       yp_real(i,j,k) = real(yp(j))
       zp_real(i,j,k) = real(zp(k))
     end do
   end do
end do


 !open(10,file='unformatted/xp//xp_'//fnum, status='unknown',form='unformatted')
  !write(10) xp_real
 !close(10)

 !open(10,file='unformatted/yp/yp_'//fnum, status='unknown',form='unformatted')
  !write(10)  yp_real
 !close(10)

 !open(10,file='unformatted/zp/zp_'//fnum, status='unknown',form='unformatted')
  !write(10) zp_real
 !close(10)

 open(10,file='data/unformatted/p/p_'//fnum, status='unknown',form='unformatted')
  write(10) p_real
 close(10)

 open(10,file='data/unformatted/u/u_'//fnum, status='unknown',form='unformatted')
  write(10) u_real
 close(10)

 open(10,file='data/unformatted/v/v_'//fnum, status='unknown',form='unformatted')
  write(10) v_real
 close(10)

 open(10,file='data/unformatted/w/w_'//fnum, status='unknown',form='unformatted')
  write(10) w_real
 close(10)

 open(10,file='data/unformatted/q/qq_'//fnum, status='unknown',form='unformatted')
  write(10) qq_real
 close(10)


  fname='data/Pressure/p_'//fnum//'.dat'
  open(10,file=fname, status='replace',form='formatted')
  do k=1,nz
    do j=1,ny
      do i=1,nx
        write(10,*) xp(i),yp(j),zp(k),p(i,j,k)
      enddo
      write(10,*)' '
    end do
    write(10,*)' '
  end do
  close(10)

  fname='data/Vector/vec_'//fnum//'.dat'
  open(10,file=fname, status='replace',form='formatted')
  do k=1,nz
    do j=1,ny
      do i=1,nx
        write(10,*) xp(i),yp(j),zp(k),u(i,j,k), v(i,j,k), w(i,j,k)
      enddo
      write(10,*)' '
    end do
    write(10,*)' '
  end do
  close(10)

  fname='data/Ave_vel/vel_'//fnum//'.dat'
  open(10,file=fname, status='replace',form='formatted')
  do j=1,ny
    write(10,*) yp(j),u_ave(j), v_ave(j), w_ave(j)
  enddo
  write(10,*)' '
  close(10)

  fname='data/Vel_sta/u_sta_'//fnum//'.dat'
  open(10,file=fname, status='replace',form='formatted')
  do j=1,ny
    write(10,*) yp(j),u_sta(j),v_sta(j),w_sta(j)
  enddo
  write(10,*)' '
  close(10)

  if(istep>200000)then
    fname='data/Vel_RMS/v_rms_'//fnum//'.dat'
    open(10,file=fname, status='replace',form='formatted')
    do j=1,ny
     write(10,*) yp(j),u_rms(j),v_rms(j),w_rms(j)
    enddo
    write(10,*)' '
    close(10)
  endif

  fname='data/Qvalue/qq'//fnum//'.dat'
  open(10,file=fname, status='replace',form='formatted')

  do k=1,nz
    do j=1,ny
      do i=1,nx
        write(10,*) xp(i),yp(j),zp(k),qq(i,j,k)
      end do
    end do
  enddo
  write(10,*)' '
  close(10)

  fname='data/cross_section_v/cs_velo_'//fnum//'.dat'
  open(10,file=fname, status='replace',form='formatted')

  do j=1,ny
    do i=1,nx
      write(10,*) xp(i),yp(j),u(i,j,nz/2)
    end do
  end do
  write(10,*)' '
  close(10)

!  fname='DATA_U/u_'//fnum//'.dat'
!  open(10,file=fname, status='replace',form='formatted')
!  do k=1,nz
!    do j=1,ny
!      do i=1,nx
!        write(10,*) xp(i),yp(j),zp(k),u(i,j,k)
!      enddo
!      write(10,*)' '
!    end do
!    write(10,*)' '
!  end do
!  close(10)

!  fname='DATA_V/v_'//fnum//'.dat'
!  open(10,file=fname, status='replace',form='formatted')
!  do k=1,nz
!    do j=1,ny
!      do i=1,nx
!        write(10,*) xp(i),yp(j),zp(k),v(i,j,k)
!      enddo
!      write(10,*)' '
!    end do
!    write(10,*)' '
!  end do
!  close(10)

!  fname='DATA_W/w_'//fnum//'.dat'
!  open(10,file=fname, status='replace',form='formatted')
!  do k=1,nz
!    do j=1,ny
!      do i=1,nx
!        write(10,*) xp(i),yp(j),zp(k),w(i,j,k)
!      enddo
!      write(10,*)' '
!    end do
!  write(10,*)' '
!  end do
!  close(10)

    return
end subroutine SaveData
