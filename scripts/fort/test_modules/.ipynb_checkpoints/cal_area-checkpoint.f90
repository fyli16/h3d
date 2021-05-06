module cal_area
  use constants
  implicit none
  real :: radius, area_val

  contains 
  subroutine area()
    area_val = pi * radius**2.0
  end subroutine area
end module cal_area


