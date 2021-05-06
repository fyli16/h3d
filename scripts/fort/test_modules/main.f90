program module_example     
  use cal_area
  use volume
  implicit none     

  real :: x, ePowerx!, area !radius 
  x = 2.0
  radius = 7.0
  ePowerx = e ** x
  !area = pi * radius**2     

  call show_consts() 
  call area()

  length = x
  call cal_volume
  print*, "Volume of a cubic of length 2.0 = ", vol  

  print*, "e raised to the power of 2.0 = ", ePowerx
  print*, "Area of a circle with radius 7.0 = ", area_val  
   
end program module_example
