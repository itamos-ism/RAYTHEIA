module m_parameters
  use MPI
  implicit none
  private
  integer,parameter,public::RK=KIND(0.D0)
  real(RK),parameter,public:: Pi = 3.141592653589793238462643383279502884_RK
  real(RK),parameter,public:: pc = 3.08568025D+18 !pc in cm

  ! parallel
  integer,public:: nrank  ! local MPI rank 
  integer,public:: nproc, nthreads  ! total number of processors
  integer,public:: IID,JID,KID
  integer,public:: ierror

  ! mesh option
  real(RK),public::xlx,yly,zlz ! domain length
  integer,public::nxp,nyp,nzp ! grid points number
  integer,public::nxc,nyc,nzc ! grid center number
  integer,public::npx,npy,npz ! processor number
  integer,public::nxnp,nynp,nznp ! local grid number
  real(RK),public::dx,dy,dz

  ! variables
  character(len=50),public :: input,indir,outdir
  integer,public :: I,J,K,II,JJ,KK
  integer,public :: level,nrays,nside,ipix
  integer,public :: maxpoints
  integer,public :: epray
  integer,public,allocatable :: projected(:,:)
  real(RK),public,allocatable :: plength(:)
  real(RK),public,allocatable :: Aveff(:,:,:)

  type pdr_node
    real(RK) :: x,y,z
    real(RK) :: rho
    real(RK), allocatable :: cd(:)
  end type pdr_node
  type(pdr_node),public,allocatable :: pdr(:,:,:)

  public::readparams
contains

  subroutine readparams
    implicit none

    open(unit=12,file='params.dat',status='old')
    read(12,*); read(12,*); read(12,*)
    read(12,*) xlx
    read(12,*) yly
    read(12,*) zlz
    read(12,*) nxc
    read(12,*) nyc
    read(12,*) nzc
    read(12,*) npx
    read(12,*) npy
    read(12,*) npz
    read(12,*) level
    read(12,*); read(12,*); read(12,*)
    read(12,*) indir
    read(12,*) input
    read(12,*) outdir
    input = trim(adjustl(indir))//'/'//trim(adjustl(input))

    ! grid parameters
    nxp=nxc+1
    nyp=nyc+1
    nzp=nzc+1
    dx = xlx/real(nxc,kind=RK)
    dy = yly/real(nyc,kind=RK)
    dz = zlz/real(nzc,kind=RK)

    ! mapping global index into local index
    KID = nrank / npx / npy
    JID = (nrank - KID * npx * npy) / npx
    IID =  nrank - KID * npx * npy - JID * npx
    nxnp = nxc/npx ! local grid cells number in x
    nynp = nyc/npy ! local grid cells number in y
    nznp = nzc/npz ! local grid cells number in z

    nside=2**level
    nrays=12*nside**2

    ! parallel parameters check
    if(nproc/=npx*npy*npz) then
      if(nrank==0) print*, 'CPU ACTIVATED:' , nproc, &
                          'CPU REQUIRED:' , npx*npy*npz
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    if(mod(nxc,npx).ne.0) then
      if(nrank==0) print*, 'nxc must be divisible by npx',nxc,npx
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    if (.not. is_power_of_two(nxnp)) then
      if(nrank==0) print*, 'nxnp is not a power of two', nxnp
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    if(mod(nyc,npy).ne.0) then
      if(nrank==0) print*, 'nyc must be divisible by npy',nyc,npy
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    if (.not. is_power_of_two(nynp)) then
      if(nrank==0) print*, 'nynp is not a power of two', nynp
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    if(mod(nzc,npz).ne.0) then
      if(nrank==0) print*, 'nzc must be divisible by npz',nzc,npz
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    if (.not. is_power_of_two(nznp)) then
      if(nrank==0) print*, 'nznp is not a power of two', nznp
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)
      stop
    endif

    ! if (nxnp /= nynp .or. nynp /= nznp) then
    !   if(nrank==0) print*, 'nxnp,nynp and nznp is not equal'
    !   call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !   call MPI_FINALIZE(ierror)
    !   stop
    ! endif

  end subroutine readparams

  function is_power_of_two(n) result(res)
    integer, intent(in) :: n
    logical :: res

    res = (n > 0) .and. (popcnt(n) == 1)
  end function is_power_of_two

end module m_parameters    
