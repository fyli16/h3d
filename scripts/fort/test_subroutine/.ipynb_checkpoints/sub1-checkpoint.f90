subroutine decomp1d( n, numprocs, myid, s, e )
    use parameters
    implicit none  
    integer numprocs, myid
    integer n
    integer s, e
    integer nlocal
    integer deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n, int(numprocs,8) )
    s       = s + min( int(myid,8) ,deficit)
    
    if (myid  <  deficit) nlocal = nlocal + 1
    e = s + nlocal - 1
    if (e  >  n .or. myid  ==  numprocs-1) e = n
    nz=324
    
    !return
end subroutine decomp1d
