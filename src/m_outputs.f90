module m_outputs
    use MPI
    use m_parameters
    use m_readdensity
    use m_calc_columndens
    implicit none
    private

    public::writeoutpus
contains

    subroutine writeoutpus
        implicit none
        integer :: iproc, istart, jstart, kstart, iu, ierr

        real(RK), allocatable,dimension(:,:,:) :: global_rho, global_Aveff

        if(nrank == 0) then
            allocate(global_rho(nxc,nyc,nzc))
            allocate(global_Aveff(nxc,nyc,nzc))
        endif

        do iproc = 0, nproc - 1
            KID = iproc / (npx * npy)
            JID = (iproc - KID * npx * npy) / npx
            IID = iproc - KID * npx * npy - JID * npx

            istart = IID * nxnp + 1
            jstart = JID * nynp + 1
            kstart = KID * nznp + 1

            if (nrank == 0) then
                if (iproc == 0) then
                    global_rho(istart:istart+nxnp-1, jstart:jstart+nynp-1, kstart:kstart+nznp-1) = pdr(1:nxnp, 1:nynp, 1:nznp)%rho
                    global_Aveff(istart:istart+nxnp-1, jstart:jstart+nynp-1, kstart:kstart+nznp-1) = Aveff(1:nxnp, 1:nynp, 1:nznp)
                else
                    call MPI_Recv(global_rho(istart:istart+nxnp-1, jstart:jstart+nynp-1, kstart:kstart+nznp-1), nxnp*nynp*nznp, MPI_DOUBLE_PRECISION, &
                                iproc, iproc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                    call MPI_Recv(global_Aveff(istart:istart+nxnp-1, jstart:jstart+nynp-1, kstart:kstart+nznp-1), nxnp*nynp*nznp, MPI_DOUBLE_PRECISION, &
                                iproc, iproc + nproc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                end if
            else
                if (nrank == iproc) then
                    call MPI_Send(pdr(1:nxnp, 1:nynp, 1:nznp)%rho, nxnp*nynp*nznp, MPI_DOUBLE_PRECISION, 0, nrank, MPI_COMM_WORLD, ierr)
                    call MPI_Send(Aveff(1:nxnp, 1:nynp, 1:nznp), nxnp*nynp*nznp, MPI_DOUBLE_PRECISION, 0, nrank + nproc, MPI_COMM_WORLD, ierr)
                end if
            end if
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do

        ! -------------------------------
        ! Rank 0 写入 ASCII 文件
        ! -------------------------------
        if (nrank == 0) then
            open(newunit=iu, file='Aveff.dat', form='formatted')
            do i = 1, nxc
            do j = 1, nyc
            do k = 1, nzc
                write(iu, *) global_rho(i,j,k), global_Aveff(i,j,k)
            end do
            end do
            end do
            close(iu)
            deallocate(global_rho, global_Aveff)
        end if        
    end subroutine writeoutpus

end module m_outputs