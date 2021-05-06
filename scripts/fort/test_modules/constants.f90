module constants  
  implicit none 
  real, parameter :: pi = 3.1415926536  
  real, parameter :: e = 2.7182818285 

  contains      
  subroutine show_consts()          
      print*, "Pi = ", pi          
      print*,  "e = ", e     
  end subroutine show_consts 
end module constants 


module volume
  implicit none
  real :: length, vol 

  contains
  subroutine cal_volume()
      vol = length**3.0
  end subroutine cal_volume
end module volume
