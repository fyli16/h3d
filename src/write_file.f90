!---------------------------------------------------------------------
subroutine wrtfile(dat,rnorm,fileName, irec_start,ny1m,nz1m)
  use parameter_mod
  implicit none

  integer:: num_sdat
  integer*8 filenum,irec_start,iry1,iry2,irz1,irz2
  real*8, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: dat
  real*4, dimension(1:nxmax-2,jb:je,kb:ke) :: stemp
  integer :: ip, iry, irz, i, j, k, recnum, ii
  integer*8 :: keg, kbg, jeg, jbg, icount,ny1m,nz1m
  double precision :: rnorm
  integer :: WriteSubArray     ! return 0 upon failure, 1 upon success
  character :: fileName*(*)    ! file to write paricles to
  integer dilo,dihi,&       ! bounds of the entire array
          djlo,djhi,&
          dklo,dkhi
  integer ailo,aihi,&       ! bounds of our portion of the array
          ajlo,ajhi,&
          aklo,akhi
  integer*8 domainCells,subDomainCells8       ! number of cells in the entire array
  integer domainDims(3)     ! dimensions of the entire array
  integer subDomainCells    ! number of cells in our portion of the array
  integer subDomainDims(3)  ! dimensions of our portion of the array
  integer subDomainStart(3) ! lo corner of the bounds of our portion
  integer iErr1,iErr2,file,eStrLen,subArray,stat(MPI_STATUS_SIZE),mode
  INTEGER(KIND=MPI_OFFSET_KIND) DISP1

  character :: eStr*(1024)
  eStrLen=1024
  WriteSubArray=0
  icount=0
  DISP1=0

  ! Open the file, bail out if this fails.
  mode=MPI_MODE_WRONLY+MPI_MODE_CREATE
  call MPI_File_open(&
      MPI_COMM_WORLD,fileName,mode,MPI_INFO_NULL,file,iErr1)
  if (iErr1.ne.MPI_SUCCESS) then
    call MPI_Error_string(iErr1,eStr,eStrLen,iErr2)
    write(0,*)'Error: Could not open file ',fileName
    write(0,*)eStr
    write(0,*)'Write aborted.'
    return
  endif

  !  begin by converting data to REAL*4
  !
  !  do k = kb-1, ke+1
  !    do j = jb-1,je+1
  !      do i = 1, nxmax
  !        sdat(i,j,k) = rnorm*dat(i,j,k)
  !      enddo
  !    enddo
  !  enddo

  ! print *,kb,ke,jb,je
  do k = kb, ke
    do j = jb,je
      do i = 2, nxmax-1
        stemp(i-1,j,k) = rnorm*dat(i,j,k)
      enddo
    enddo
  enddo

  dilo = 1
  dihi = nxmax-2
  ailo = 1
  aihi = nxmax-2
  djlo = jbglobal(0)
  djhi = jeglobal(npes-1)
  ajlo = jbglobal(myid)
  ajhi = jeglobal(myid)
  dklo = kbglobal(0)
  dkhi = keglobal(npes-1)
  aklo = kbglobal(myid)
  akhi = keglobal(myid)
  call BoundsToDimensions(&
      dilo,dihi,djlo,djhi,dklo,dkhi,domainDims,domainCells)
  call BoundsToDimensions(&
      ailo,aihi,ajlo,ajhi,aklo,akhi,subDomainDims,subDomainCells8)
  subDomainCells=subDomainCells8
  subDomainStart(1)=ailo-dilo !  convert to c-index!
  subDomainStart(2)=ajlo-djlo
  subDomainStart(3)=aklo-dklo
  call MPI_Type_create_subarray(&
      3,domainDims,subDomainDims,subDomainStart,&
      MPI_ORDER_FORTRAN,MPI_REAL,subArray,iErr1)
  call MPI_Type_commit(subArray,iErr1)

  ! Set the file view
  call MPI_File_set_view(&
      file,DISP1,MPI_REAL,subArray,"native",MPI_INFO_NULL,iErr1)

  ! Write
  call MPI_File_write_all(file,stemp,subDomainCells,MPI_REAL,stat,iErr1)
  call MPI_File_close(file,iErr2)
  call MPI_Type_free(subArray,iErr2)

  if (iErr1.ne.MPI_SUCCESS) then
    call MPI_Error_string(iErr1,eStr,eStrLen,iErr2)
    write(0,*)'Error: Could not write to file ',fileName
    write(0,*) eStr(1:eStrLen)
    write(0,*) 'Write aborted by rank',myid
    write(0,*)'Write aborted.'
  endif

  return
end subroutine wrtfile


!---------------------------------------------------------------------
subroutine BoundsToDimensions(ilo, ihi, jlo, jhi, klo, khi, dims, nCells)
  implicit none

  integer:: ilo,ihi,&    ! bounds of a 3d array
            jlo,jhi,&
            klo,khi
  integer :: dims(3)      ! dimensions of the array (out)
  integer*8 :: nCells       ! total number of cells (out)

  dims(1) = ihi-ilo+1
  dims(2) = jhi-jlo+1
  dims(3) = khi-klo+1
  nCells = dims(1)*dims(2)*dims(3)
  
  return
end subroutine BoundsToDimensions


!---------------------------------------------------------------------
subroutine wrtfile_non_mpio(dat,rnorm,filenum, irec_start,ny1m,nz1m) 
  use parameter_mod
  implicit none

  integer:: num_sdat,filenum

  integer*8 :: irec_start,iry1,iry2,irz1,irz2
  double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: dat
  real*4, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: sdat
  integer :: ip, iry, irz, i, j, k, recnum, ii
  integer*8 :: keg, kbg, jeg, jbg, icount,ny1m,nz1m
  double precision :: rnorm

  icount=0

  ! begin by converting data to REAL*4
  do k = kb-1, ke+1
    do j = jb-1,je+1
      do i = 1, nxmax
        sdat(i,j,k) = rnorm*dat(i,j,k)
      enddo
    enddo
  enddo

  ! send data to myid = 0
  num_sdat = nxmax*(nylmax+2)*(nzlmax + 2)
  if(myid.ne.0) then
    call MPI_ISEND(sdat(1,jb-1,kb-1),num_sdat,MPI_REAL4,&
         0   ,1,MPI_COMM_WORLD,req(1),IERR)
    call MPI_WAITALL(1,req,status_array,IERR)
  else
    do ip = 0, npes-1
      keg = keglobal(ip)
      kbg = kbglobal(ip)
      jeg = jeglobal(ip)
      jbg = jbglobal(ip)

      if (jbg.eq.1.and.jeg.eq.ny1m) then
        iry1 = 0        +1
        iry2 = ny1m+1    -1
        goto 10
      endif

      if(jbg.eq.1) then
        iry1 = 0        +1
        iry2 = jeg
        j    = jb - 2
      else if(jeg.eq.ny1m) then
	      iry1 = jbg
	      iry2 = jeg+1    -1
	      j = jb - 1
      else
        iry1 = jbg
        iry2 = jeg
        j = jb - 1
      endif

10    continue

      if (kbg.eq.1.and.keg.eq.nz1m) then
        irz1=0          +1
        irz2=nz1m+1     -1
        goto 20
      endif

      if(kbg.eq.1) then
        irz1 = 0        +1
        irz2 = keg
        k = kb - 2
      else if(keg.eq.nz1m) then
        irz1 = kbg
        irz2 = keg+1    -1
        k = kb - 1
      else
        irz1 = kbg
        irz2 = keg
        k = kb - 1
      endif

20    continue

      if(ip.ne.0) then
        call MPI_IRECV(sdat(1,jb-1,kb-1),num_sdat,MPI_REAL4,&
             IP,1,MPI_COMM_WORLD,req(1),IERR)
        call MPI_WAITALL(1,req,status_array,IERR)
      endif
      
      do irz = irz1, irz2
	      k = k+1
	      do iry = iry1, iry2
          j = j+1 
          icount=icount+1
          recnum = irec_start + iry -1 + (irz-1)*(nymax - 2)
          write(filenum, rec=recnum) (sdat(i,iry-jbg+1,irz-kbg+1),i=1   +1,nxmax    -1)
        enddo
      enddo
    enddo

    write(6,*) " number of records written = ", icount

  endif

  return
end subroutine wrtfile_non_mpio


