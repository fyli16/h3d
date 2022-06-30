module m_functions
  implicit none 

  contains

  !---------------------------------------------------------------------
  ! F(t) = 0 to solve
  !---------------------------------------------------------------------
  double precision function func(t) 
    real*8, intent(in) :: t 
    real*8 :: rhs
    integer*8 :: n

    common /fparams/ rhs,n

    if( t <=1. ) then
      func = 1./real(n) - rhs
    else
      func = (t-1.)/(t**n-1.) - rhs
    endif
    
    return
  end function func 


  !---------------------------------------------------------------------
  ! simple root finder
  !---------------------------------------------------------------------
  subroutine bisect(f,a,b,tol)
    real*8, intent(in) :: tol
    real*8, intent(inout) :: a,b
    real*8 :: c,u,v,w
    real*8 :: f

    u = f(a) ; v = f(b)      
    if(u*v.gt.0.) then
      call error_abort('bisect(): bad initial interval --- stop!')
    else if(u.eq.0.) then
      b=a
      return
    else if(v.eq.0.) then
      a=b
      return
    endif

    do while((abs(f(a)).gt.tol).or.(abs(f(b)).gt.tol) ) 
      c = 0.5*(a+b) 
      w = f(c)    
      if(w.eq.0.) then 
        a=c ; b=c
        return;
      endif
      if(w*u.lt.0.) then
        b=c ; v=w       
      else
        a=c ; u=w       
      endif
    enddo

    return      
  end subroutine bisect


  !---------------------------------------------------------------------1
  ! Helper functions
  ! Solve for x: (x-1)/(x^N-1) - rhs = 0 
  ! by using a simple bisection method
  !---------------------------------------------------------------------1
  double precision function findexp(rhsi,ni)
    real*8, intent(in) :: rhsi
    integer*8, intent(in) :: ni
    real*8 :: tol,rhs,af,bf
    integer*8 n
    common /fparams/ rhs,n

    tol = 10.*epsilon(real(0))
    ! These are common block parameters
    rhs=rhsi
    n = ni
    ! These are common block parameters
    if(func(1.0D0).le.tol) then
      call error_abort('findexp(): dx_uniform too large --- stop!')
    endif
    af = 1.
    bf = af
    do while(func(bf).ge.0.)
      bf = bf*2.
    enddo
    call bisect(func,af,bf,tol)
    findexp = 0.5*(af+bf)

    return
  end function findexp


end module m_functions