module m_readdensity
  use MPI
  use omp_lib
  use m_parameters
  use m_Healpix
  use m_Raytheia
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

    maxpoints = 500
    if(nrank.eq.0) print*,'Maxpoints          = ',maxpoints

    ! 构建树
    call BuildLinearOctree(real(pdr%rho))

    ! ==================================================================
    ! [Verification] Linear octree test (验证平均密度)
    ! ==================================================================
    ! block
    !     integer :: i, fid
    !     integer(i8b) :: code
    !     integer :: lx, ly, lz ! 局部逻辑坐标 (finest unit)
    !     integer :: lvl
    !     real(RK) :: dens
    !     real(RK) :: phys_x, phys_y, phys_z, phys_size
    !     real(RK) :: cell_len_x, cell_len_y, cell_len_z
    !     integer :: max_dim_log
    !     character(len=64) :: filename
    !     logical :: is_padding
        
    !     ! 重新计算逻辑最大维度 (与 Build 过程一致)
    !     max_dim_log = max(nxnp, max(nynp, nznp))

    !     ! 生成文件名: octree_check_000.txt
    !     write(filename, "('octree_check_', I3.3, '.txt')") nrank
    !     fid = 100 + nrank
    !     open(unit=fid, file=trim(filename), status='replace')

    !     ! 写入表头
    !     write(fid, '(A)') "ID, Level, Code, Global_X, Global_Y, Global_Z, Size_Phys, Density, Is_Padding"

    !     print *, "Rank", nrank, ": Verifying", n_linear_leaves, "nodes..."

    !     ! 3. 遍历所有叶子节点
    !     do i = 1, n_linear_leaves
    !         code = LinearCodes(i)
    !         lvl  = LinearLevels(i)
    !         dens = LinearDensity(i)

    !         ! A. 解码 Morton 码 -> 得到局部逻辑坐标 (0-based, based on finest grid)
    !         ! 注意：这里得到的是该 Block 左下角的坐标
    !         call DecodeMorton3D(code, lx, ly, lz)

    !         ! B. 计算物理尺寸
    !         ! Build逻辑：Root是level 0 (size=max_dim), level+1 size减半
    !         ! 所以当前 Grid 单元数 = max_dim / (2^lvl)
    !         ! 注意：Fortran实数除法
    !         phys_size = (real(max_dim_log, RK) / (2.0_RK**lvl)) 
            
    !         ! 各个方向的物理长度 (假设 dx, dy, dz 可能不同)
    !         cell_len_x = phys_size * dx
    !         cell_len_y = phys_size * dy
    !         cell_len_z = phys_size * dz

    !         ! C. 计算全局物理坐标 (Global Physical Min Corner)
    !         ! 局部逻辑坐标 lx + MPI偏移 (IID*nxnp) -> 全局逻辑坐标
    !         ! 乘以 dx -> 物理坐标
    !         phys_x = (real(lx, RK) + real(IID * nxnp, RK)) * dx
    !         phys_y = (real(ly, RK) + real(JID * nynp, RK)) * dy
    !         phys_z = (real(lz, RK) + real(KID * nznp, RK)) * dz

    !         ! D. 检查是否为填充区域 (Padding)
    !         ! 如果局部坐标 lx 超过了物理边界 nxnp，或者 ly > nynp ...
    !         is_padding = .false.
    !         if (lx >= nxnp .or. ly >= nynp .or. lz >= nznp) then
    !             is_padding = .true.
    !         endif

    !         ! E. 写入文件
    !         ! 格式: ID, Lvl, Code, X, Y, Z, Size, Rho, Padding?
    !         write(fid, "(I6, ',', I3, ',', I15, ',', 3(ES14.6, ','), ES14.6, ',', ES14.6, ',', L1)") &
    !             i, lvl, code, phys_x, phys_y, phys_z, cell_len_x, dens, is_padding

    !     end do

    !     close(fid)
    !     print *, "Rank", nrank, ": Octree verification written to ", trim(filename)
    ! end block

    ! ! 停止程序以便查看结果，而不是继续跑光线追踪
    ! call MPI_Barrier(MPI_COMM_WORLD, ierror)
    ! stop "Check complete. Inspect output files."

    ! ==================================================================
    ! [Verification] Diagonality Test (验证对角线积分精度)
    ! ==================================================================
    ! block
    !     type(HEALPix_ray) :: test_ray
    !     real(RK) :: total_len_local, total_len_global, expected_len
    !     real(RK) :: dx_domain, dy_domain, dz_domain, adaptive_step
    !     integer :: i, cr, cIID, cJID, cKID, ip, ierr
        
    !     ! 临时定义 DDA 输出数组 (假设只测试 1 条光线)
    !     integer :: test_epray
    !     integer :: test_projected(0:maxpoints)
    !     real(RK) :: test_plength(0:maxpoints)
    !     integer(i8b) :: temp_codes(n_linear_leaves_max)
    !     integer :: temp_levels(n_linear_leaves_max)

        
    !     test_ray%origin = [0.0_RK, 0.0_RK, 0.0_RK]
    !     test_ray%angle(1) = acos(1.0_RK / sqrt(3.0_RK))
    !     test_ray%angle(2) = 0.785398163397448_RK

    !     total_len_global = 0.0_RK
    !     do cr=0,nproc-1
    !       cKID = cr / npx / npy
    !       cJID = (cr - cKID * npx * npy) / npx
    !       cIID =  cr - cKID * npx * npy - cJID * npx

    !       temp_codes = LinearCodes
    !       temp_levels = LinearLevels

    !       call MPI_Bcast(temp_codes, n_linear_leaves_all(cr), MPI_INTEGER8, &
    !                cr, MPI_COMM_WORLD, ierr)
    !       call MPI_Bcast(temp_levels, n_linear_leaves_all(cr), MPI_INTEGER, &
    !                cr, MPI_COMM_WORLD, ierr)

    !       ! 初始化
    !       box1%min = [DBLE(cIID*nxnp)*dx, DBLE(cJID*nynp)*dy, DBLE(cKID*nznp)*dz]
    !       box1%max = [DBLE((cIID+1)*nxnp)*dx, DBLE((cJID+1)*nynp)*dy, DBLE((cKID+1)*nznp)*dz]
    !       test_epray = 0
    !       test_projected = 0
    !       test_plength = 0.D0
          
    !       ! 2. 调用 DDA 算法
    !       call RayTheia_Linear_DDA(test_ray, box1, n_linear_leaves_all(cr), temp_codes, temp_levels, 0, &
    !                               test_epray, test_projected, test_plength)
    !       IF (test_epray.GT.0) THEN
    !           do ip=1,test_epray
    !             adaptive_step = test_plength(ip)
    !             total_len_global = total_len_global + adaptive_step
    !           enddo
    !       ENDIF           
    !     enddo
        
    !     ! 3. 输出验证结果 (仅 Rank 0)
    !     if (nrank == 0) then
    !         ! 计算理论长度: sqrt(Lx^2 + Ly^2 + Lz^2)
    !         ! 假设 xlx, yly, zlz 是全局物理尺寸变量
    !         dx_domain = real(xlx, RK)
    !         dy_domain = real(yly, RK) ! 如果是立方体可能都是 xlx
    !         dz_domain = real(zlz, RK)
    !         expected_len = sqrt(dx_domain**2 + dy_domain**2 + dz_domain**2)
            
    !         print *
    !         print *, "========== DDA Verification =========="
    !         print *, "Test Ray      : Diagonal (0,0,0) -> (L,L,L)"
    !         print *, "Computed Len  : ", total_len_global
    !         print *, "Expected Len  : ", expected_len
    !         print *, "Abs Error     : ", abs(total_len_global - expected_len)
    !         print *, "Rel Error     : ", abs(total_len_global - expected_len)/expected_len
            
    !         if (abs(total_len_global - expected_len) < 1.0e-3_RK) then
    !              print *, "Result        : [PASS] Geometry & Integration OK"
    !         else
    !              print *, "Result        : [FAIL] Check intersections/EPS/boundary"
    !         endif
    !         print *, "======================================"
    !         print *
    !     endif
        
    !     call MPI_Barrier(MPI_COMM_WORLD, ierror)
    !     stop
        
    ! end block

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

    allocate(plength(0:maxpoints))
    allocate(projected(0:maxpoints))

  end subroutine initialization

end module m_readdensity
