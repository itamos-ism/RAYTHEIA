module m_readdensity
  use MPI
  use omp_lib
  use m_parameters
  use m_Healpix
  use m_Ray_box
  implicit none

  public::readdensity,initialization
contains

  subroutine readdensity
    implicit none

    ! locals
    integer :: maxpoints_loc,temp,pdr_ptot,GI,GJ,GK
    real(RK) :: xpos,ypos,zpos,denst
    real(RK) :: xc,yc,zc
    real(RK) :: maxdens,mindens

    open(unit=2,file=input,status='old')

    read(2,*); read(2,*)
    pdr_ptot = nxc*nyc*nzc
    allocate(pdr(1:nxnp,1:nynp,1:nznp))
    ! read input file and assign values to local variables
    do II=1,nxc
    do JJ=1,nyc
    do KK=1,nzc
      i=II-IID*nxnp
      j=JJ-JID*nynp
      k=KK-KID*nznp
      read(2,*) xpos,ypos,zpos,denst
      if(i.ge.1.and.i.le.nxnp) then
        if(j.ge.1.and.j.le.nynp) then
          if(k.ge.1.and.k.le.nznp) then
            pdr(i,j,k)%x=xpos
            pdr(i,j,k)%y=ypos
            pdr(i,j,k)%z=zpos
            pdr(i,j,k)%rho=denst
          endif
        endif
      endif
    enddo
    enddo
    enddo

    maxdens = maxval(pdr(:,:,:)%rho)
    mindens = minval(pdr(:,:,:)%rho)
    call MPI_AllReduce(MPI_IN_PLACE, maxdens, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    call MPI_AllReduce(MPI_IN_PLACE, mindens, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)

    if(nrank.eq.0) then
      write(6,*) 'Elements           = ',pdr_ptot
      write(6,*) 'Maximum density    = ',maxdens
      write(6,*) 'Minimum density    = ',mindens
      write(6,*) 'Grid resolution    = ',nxc
      write(6,*) 'Domain size        = ',real(xlx)
    endif

    levels = nint(log(DBLE(nxnp)) / log(2.D0)) + 1
    maxpoints = 500
    if(nrank.eq.0) print*,'Maxpoints          = ',maxpoints

  end subroutine readdensity

  subroutine initialization
    implicit none
    integer :: max_nodes

    allocate(Aveff(nxnp,nynp,nznp))
    do k=1,nznp
    do j=1,nynp
    do i=1,nxnp
      allocate(pdr(i,j,k)%cd(0:nrays-1))
    enddo
    enddo
    enddo

#ifdef AMR
    max_nodes = (8**(levels) - 1) / 7
    allocate(tree(max_nodes))
    ! print*,nrank,levels,max_nodes
#else
    allocate(plength(0:maxpoints))
    allocate(projected(0:maxpoints,3))
#endif

  end subroutine initialization

end module m_readdensity
