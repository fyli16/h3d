program test_c_interface
  implicit none
  integer:: my_short_int
  character(len=160):: myid_char
  my_short_int=1184
  call integer_to_character(myid_char,len(myid_char),my_short_int)
  write(6,*)myid_char
end program
