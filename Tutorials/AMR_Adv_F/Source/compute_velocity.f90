module compute_velocity_module

  use ml_layout_module
  use bl_constants_module

  implicit none

  private

  public :: compute_velocity

contains
  
  subroutine compute_velocity(mla,velocity,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: velocity(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)

    ng = velocity(1,1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,nfabs(velocity(n,1))
          dp1 => dataptr(velocity(n,1),i)
          dp2 => dataptr(velocity(n,2),i)
          lo = lwb(get_box(velocity(n,1),i))
          hi = upb(get_box(velocity(n,1),i))
          select case(dm)
          case (2)
             call compute_velocity_2d(dp1(:,:,1,1), dp2(:,:,1,1), ng, &
                                      lo, hi, dx(n), time)
          case (3)
             dp3 => dataptr(velocity(n,3),i)
             call compute_velocity_3d(dp1(:,:,:,1), dp2(:,:,:,1), dp3(:,:,:,1), ng, &
                                      lo, hi, dx(n), time)
          end select
       end do
    end do

  end subroutine compute_velocity

  subroutine compute_velocity_2d(velx, vely, ng, lo, hi, dx, time)

    integer          :: lo(2), hi(2), ng
    double precision :: velx(lo(1)-ng:,lo(2)-ng:)
    double precision :: vely(lo(1)-ng:,lo(2)-ng:)
    double precision :: dx,time
 
    ! local varables
    integer          :: i,j
    double precision :: x,y
    double precision :: psi(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2)

    ! streamfunction psi
      do j = lo(2)-2, hi(2)+2
         y = (dble(j)+0.5d0)*dx
         do i = lo(1)-2, hi(1)+2
            x = (dble(i)+0.5d0)*dx
            psi(i,j) =  sin(M_PI*x)**2 * sin(M_PI*y)**2 * cos (M_PI*time/2.d0) / M_PI
         end do
      end do

    ! x velocity
    do j=lo(2)-1,hi(2)+1
       y = (dble(j)+0.5d0) * dx
       do i=lo(1)-1,hi(1)+2
          x = dble(i) * dx
          velx(i,j) =  -( (psi(i,j+1)+psi(i-1,j+1)) - (psi(i,j-1)+psi(i-1,j-1)) ) / (4.d0*dx)
       end do
    end do

    ! y velocity
    do j=lo(2)-1,hi(2)+2
       y = dble(j) * dx
       do i=lo(1)-1,hi(1)+1
          x = (dble(i)+0.5d0) * dx
          vely(i,j) = ( (psi(i+1,j)+psi(i+1,j-1)) - (psi(i-1,j)+psi(i-1,j-1)) ) / (4.d0*dx)
       end do
    end do

    end subroutine compute_velocity_2d

    subroutine compute_velocity_3d(velx, vely, velz, ng, lo, hi, dx, time)

    integer          :: lo(3), hi(3), ng
    double precision :: velx(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: vely(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: velz(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: dx,time
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z



  end subroutine compute_velocity_3d

end module compute_velocity_module
