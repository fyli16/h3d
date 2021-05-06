integer function myfunction1()
implicit none
integer::i
myfunction1=i
i=i+1
return
end function

integer function myfunction2()
implicit none
integer::i
myfunction2=i
i=i+1
return
end function

program internal1
    implicit none
    print *, myfunction1()
    print *, myfunction1()
    print *, myfunction2()
    print *, myfunction2()
    
end program internal1