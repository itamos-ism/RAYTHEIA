program main
  use MPI
  use omp_lib
  use m_parameters
  use m_readdensity
  use m_calc_columndens
  use m_outputs
  implicit none

  ! locals
  integer :: nUnit
  real(RK) :: time(20)

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read parameters
  call readparams

  ! read density distribution
  call readdensity

  ! allocate variables
  call initialization

  ! Calculate columndensity
  time(1)=MPI_WTIME()

#ifdef HAMMER
  if(nproc.eq.1) then
    call calc_hammermap(nxc/2,nyc/2,nzc/2)
  else
    stop "Error: hammermap only for OpenMP"
  endif 
#else
  call calc_columndens
#endif

  time(2)=MPI_WTIME()
  if(nrank==0) print*,'Elapse time for variables calculation', time(2)-time(1), 'with nodes number',nproc, 'with threads number',nthreads

  ! write output
#ifdef HAMMER
  open(newunit=nUnit, file=trim(adjustl(outdir))//'/cdMaps.dat', status='replace')
  do ipix=0,nrays-1
    call pix2ang_nest(nside,ipix,thfpix,phfpix)
    write(nUnit,*) thfpix,phfpix,single_cd(ipix)
  enddo
  close(nUnit)
#else
  call writeoutpus
#endif

  deallocate(pdr)
  call MPI_FINALIZE(ierror)
end program main
