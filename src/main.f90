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

  call calc_columndens

  time(2)=MPI_WTIME()
  if(nrank==0) print*,'Elapse time for variables calculation', time(2)-time(1), 'with nodes number',nproc, 'with threads number',nthreads

  ! write output
  call writeoutpus

  deallocate(pdr)
  call MPI_FINALIZE(ierror)
end program main
