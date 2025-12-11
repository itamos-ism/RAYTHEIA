module m_calc_columndens
  use MPI
  use omp_lib
  use m_parameters
  use m_Healpix
  use m_Ray_box
  implicit none

  public::calc_columndens
contains

  subroutine calc_columndens
    implicit none

    !locals
    real(RK) :: xc,yc,zc,AV_fac,cd,adaptive_step,gamma
    real, allocatable :: remote_rho(:,:,:)
    integer :: cr,sI,sJ,sK,cIID,cJID,cKID,sII,sJJ,sKK,ip,nUnit,ierror

    allocate(remote_rho(nxnp,nynp,nznp))

Aveff = 0.D0
AV_fac = 6.29e-22
gamma = 3.02D0

#ifdef OPENMP
!$OMP PARALLEL DO COLLAPSE(3)
#endif
    do sK = 1,nznp
    do sJ = 1,nynp
    do sI = 1,nxnp
      pdr(sI,sJ,sK)%cd = 0.D0
    enddo
    enddo
    enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

do cr=0,nproc-1
   cKID = cr / npx / npy
   cJID = (cr - cKID * npx * npy) / npx
   cIID =  cr - cKID * npx * npy - cJID * npx

#ifdef OPENMP
!$OMP PARALLEL DO COLLAPSE(3)
#endif
   do sK = 1,nznp
   do sJ = 1,nynp
   do sI = 1,nxnp
      remote_rho(sI,sJ,sK) = REAL(pdr(sI,sJ,sK)%rho)
   enddo
   enddo
   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

   call MPI_BCAST(remote_rho,nxnp*nynp*nznp,MPI_REAL,cr,MPI_COMM_WORLD,ierror)

#ifdef AMR
   box1%min = [DBLE(cIID*nxnp)*dx, DBLE(cJID*nynp)*dy, DBLE(cKID*nznp)*dz]
   box1%max = [DBLE((cIID+1)*nxnp)*dx, DBLE((cJID+1)*nynp)*dy, DBLE((cKID+1)*nznp)*dz]
   call build_tree(box1, cIID, cJID, cKID, levels-1, real(remote_rho))
#endif

#ifdef OPENMP
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) PRIVATE(sII,sJJ,sKK,xc,yc,zc,adaptive_step) &
!$OMP PRIVATE(thfpix,phfpix,ray,box1,epray,projected,plength,cd)
#endif
   do sK=1,nznp
   do sJ=1,nynp
   do sI=1,nxnp
      ! source of rays
      sII=sI+IID*nxnp ! Global index
      sJJ=sJ+JID*nynp
      sKK=sK+KID*nznp
      xc=(real(sII,kind=RK)-0.5D0)*dx
      yc=(real(sJJ,kind=RK)-0.5D0)*dy
      zc=(real(sKK,kind=RK)-0.5D0)*dz
      box1%min = [DBLE(cIID*nxnp)*dx, DBLE(cJID*nynp)*dy, DBLE(cKID*nznp)*dz]
      box1%max = [DBLE((cIID+1)*nxnp)*dx, DBLE((cJID+1)*nynp)*dy, DBLE((cKID+1)*nznp)*dz]
      do ipix=0,nrays-1
         epray = 0
         projected(:,1) = sI
         projected(:,2) = sJ
         projected(:,3) = sK
         plength = 0.D0
         call pix2ang_nest(nside, ipix, thfpix, phfpix)
         ray%origin = [xc, yc, zc]
         ray%angle = [thfpix, phfpix]
#ifdef AMR
         cd = 0.D0
         call raytheia_amr(ray, 1, remote_rho, cd)
         pdr(sI,sJ,sK)%cd(ipix) = pdr(sI,sJ,sK)%cd(ipix) + cd
#else 
         call raytheia_sm(ray, box1, levels-1, cIID, cJID, cKID, epray, projected, plength)
         IF (epray.GT.0) THEN
            do ip=1,epray
               adaptive_step = plength(ip)
               pdr(sI,sJ,sK)%cd(ipix) = pdr(sI,sJ,sK)%cd(ipix) + &
               &DBLE(remote_rho(projected(ip,1),projected(ip,2),projected(ip,3)))*adaptive_step*pc
            enddo
         ENDIF
#endif
      enddo
   enddo
   enddo
   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

enddo

#ifdef OPENMP
!$OMP PARALLEL DO COLLAPSE(3)
#endif
   do sK=1,nznp
   do sJ=1,nynp
   do sI=1,nxnp
      do ipix=0,nrays-1
         Aveff(sI,sJ,sK) = Aveff(sI,sJ,sK) + exp(-gamma*(pdr(sI,sJ,sK)%cd(ipix)*AV_fac))/DBLE(nrays)
      enddo
      Aveff(sI,sJ,sK) = -log(Aveff(sI,sJ,sK))/gamma
   enddo
   enddo
   enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

  end subroutine calc_columndens

end module m_calc_columndens
