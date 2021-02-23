      program endian_test
        implicit none
        integer*4 :: i
        byte :: c
        open(13,file='one.dat',form='unformatted',access='stream',status='unknown')
        write(13)1
        close(13)
        open(13,file='one.dat',form='unformatted',access='stream',status='old')
        do i=1,4
          read(13)c
          write(*,*)c
        enddo
        if (c .eq. 1) then
          write(*,*) 'Big Endian'
        else
          write(*,*) 'Little Endian'
        endif
      end program
