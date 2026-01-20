!Written by Zhengping Zhu
module m_Raytheia
    use MPI
    use m_parameters
    use healpix_types
    implicit none
    private

    type, public :: box
        real(RK) :: min(3), max(3)
    end type box

    type, public :: slab
        real(RK) :: origin(3),  dir(3),  dir_inv(3) 
    end type slab

    type, public :: HEALPix_ray
        integer :: eval
        real(RK) :: length 
        real(RK) :: origin(3), angle(2)
    end type HEALPix_ray

    type(box),public :: box1
    type(HEALPix_ray),public :: ray
    real(kind=dp),public :: thfpix, phfpix, contribution
    real(RK), parameter :: REFINE_CRITERIA = 1.D0 ! Refine criteria：if max number < 1, merge
    real(RK), parameter :: UNIFORM_CRITERIA = 0.01_RK ! Uniform criteria: if (Max - Min) / Max < 1%, merge

    ! Linear octree data structure for AMR
    integer(i8b), public, allocatable :: LinearCodes(:)    ! Store Morton Codes (64-bit)
    integer,      public, allocatable :: LinearLevels(:)   ! Store the level of the grid cell in the octree
    real(RK),     public, allocatable :: LinearDensity(:)  ! Store density
    integer,      public :: n_linear_leaves, n_linear_leaves_max, n_linear_leaves_total ! The number of leave nodes
    integer, allocatable, public :: n_linear_leaves_all(:)
    integer, public :: max_tree_level                      ! Max depth of the Octree

    public:: BuildLinearOctree, RayTheia_Linear_DDA
contains

!###############################
!       Basic functions        !
!###############################

    logical function isfinite(x)
            real(RK), intent(in) :: x
            isfinite = (abs(x) < huge(x) .and. x == x) 
    end function isfinite

    subroutine intersections(ray_xyz, box_in, tmin, tmax, length)
    implicit none
    type(slab), intent(in) :: ray_xyz
    type(box), intent(in) :: box_in
    real(RK), intent(out) :: tmin, tmax
    real(RK), intent(inout) :: length

    ! locals
    integer :: d
    real(RK) :: t1, t2
    real(RK) :: Pmin(3), Pmax(3)
    real(RK) :: distance
    real(RK) :: temp

    ! 初始化 tmin 和 tmax
    tmin = 0.0
    tmax = huge(0.0)  ! 设置 tmax 为无穷大

    ! 遍历 x, y, z 三个维度
    do d = 1, 3
        if (isfinite(ray_xyz%dir_inv(d))) then
            ! 计算 t1 和 t2
            t1 = (box_in%min(d) - ray_xyz%origin(d)) * ray_xyz%dir_inv(d)
            t2 = (box_in%max(d) - ray_xyz%origin(d)) * ray_xyz%dir_inv(d)

            ! 确保 t1 是较小值，t2 是较大值
            if (t1 > t2) then
                temp = t1
                t1 = t2
                t2 = temp
            end if

            ! 更新 tmin 和 tmax
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
        else if (ray_xyz%origin(d) < box_in%min(d) .or. ray_xyz%origin(d) > box_in%max(d)) then
            ! 射线与某个维度的边界盒平行且射线起点不在盒子内
            tmin = huge(0.0)
            exit
        end if
    end do

    ! only can be used when direction is strict unit vector
    if(tmin <= tmax) then
        length = tmax -tmin
    else
        length = 0.D0
    endif

    ! 判断射线是否与盒子相交
    ! if (tmin <= tmax) then
    !     ! 计算交点坐标
    !     Pmin = ray_xyz%origin + tmin / ray_xyz%dir_inv
    !     Pmax = ray_xyz%origin + tmax / ray_xyz%dir_inv
    !     distance = sum((Pmax - Pmin)**2)
    !     length = sqrt(distance)
    ! else
    !     length = 0.0  ! 如果不相交，长度为 0
    ! end if

    end subroutine intersections

    subroutine box_avg_density(box_in, pid, pjd, pkd, rho, avg_rho, min_rho, max_rho)
    implicit none
        type(box),intent(in) :: box_in
        integer, intent(in) :: pid, pjd, pkd
        real, intent(in) :: rho(:,:,:)
        real(RK), intent(out) :: avg_rho
        real(RK), intent(out), optional :: min_rho, max_rho 

        integer :: i1, i2, j1, j2, k1, k2 !local indices
        integer :: ig1, ig2, jg1, jg2, kg1, kg2 ! global indices
        real(RK) :: vol, total_rho
        real(RK) :: local_min, local_max, val
        integer :: i, j, k

        ig1 = nint(box_in%min(1)/dx) + 1
        ig2 = nint(box_in%max(1)/dx)
        jg1 = nint(box_in%min(2)/dy) + 1
        jg2 = nint(box_in%max(2)/dy)
        kg1 = nint(box_in%min(3)/dz) + 1
        kg2 = nint(box_in%max(3)/dz)

        i1 = ig1 - pid * nxnp
        i2 = ig2 - pid * nxnp
        j1 = jg1 - pjd * nynp
        j2 = jg2 - pjd * nynp
        k1 = kg1 - pkd * nznp
        k2 = kg2 - pkd * nznp

        ! boundary protection
        if (i1 < 1) i1 = 1
        if (j1 < 1) j1 = 1
        if (k1 < 1) k1 = 1
        if (i2 > nxnp) i2 = nxnp
        if (j2 > nynp) j2 = nynp
        if (k2 > nznp) k2 = nznp

        total_rho = 0.0_RK
        local_min = huge(0.0_RK)
        local_max = -huge(0.0_RK)
        if (i2 >= i1 .and. j2 >= j1 .and. k2 >= k1) then
            if (present(min_rho) .or. present(max_rho)) then
                do k = k1, k2
                do j = j1, j2
                do i = i1, i2
                    val = rho(i,j,k)
                    total_rho = total_rho + val
                    if (val < local_min) local_min = val
                    if (val > local_max) local_max = val
                enddo
                enddo
                enddo
            else
                total_rho = DBLE(sum(rho(i1:i2, j1:j2, k1:k2)))
            endif
            vol = DBLE((i2-i1+1) * (j2-j1+1) * (k2-k1+1))
            if (vol > 0.0_RK) then
                avg_rho = total_rho / vol
            else
                avg_rho = 0.0_RK
            endif
        else
            ! 越界或无效区域处理：不要 STOP，而是返回真空 (0)
            ! 这在构建八叉树的边界填充（Padding）时非常常见
            avg_rho = 0.0_RK
            local_min = 0.0_RK
            local_max = 0.0_RK          
        endif

        ! 输出可选参数
        if (present(min_rho)) min_rho = local_min
        if (present(max_rho)) max_rho = local_max        
 
    end subroutine box_avg_density

!###############################
!         Morton Codes         !
!###############################

    ! Morton编码函数 (3D Grid Index (ix, iy, iz) -> Morton)
    function EncodeMorton3D(ix, iy, iz) result(code)
        integer, intent(in) :: ix, iy, iz
        integer(i8b) :: code
        integer(i8b) :: x, y, z

        x = int(ix, i8b)
        y = int(iy, i8b)
        z = int(iz, i8b)

        x = SplitBy3(x)
        y = SplitBy3(y)
        z = SplitBy3(z)

        ! bit interleaving: ...zyxzyx
        code = IOR(IOR(x, ishft(y, 1)), ishft(z, 2))
    end function EncodeMorton3D

    ! Bit Interleaving (3D -> 1D)
    function SplitBy3(a) result(r)
        integer(i8b), intent(in) :: a
        integer(i8b) :: r
        integer :: i
        
        r = 0_i8b
        do i = 0, 20
            if (btest(a, i)) r = ibset(r, i*3)
        end do
    end function SplitBy3

    ! Morton解码函数 (Morton -> 3D Grid Index (ix, iy, iz))
    subroutine DecodeMorton3D(code, ix, iy, iz)
        integer(i8b), intent(in) :: code
        integer, intent(out) :: ix, iy, iz
        integer :: i
        
        ix = 0; iy = 0; iz = 0
        do i = 0, 20
            if (btest(code, i*3))   ix = ibset(ix, i)
            if (btest(code, i*3+1)) iy = ibset(iy, i)
            if (btest(code, i*3+2)) iz = ibset(iz, i)
        end do
    end subroutine DecodeMorton3D

!################################
!    Construct Linear Octree    !
!################################

    subroutine BuildLinearOctree(amr_var)
        implicit none
        real, intent(in) :: amr_var(:,:,:)
        integer :: max_possible_nodes, max_dim, ierror
        
        ! 1. 估算最大节点数 (局部网格大小)
        max_possible_nodes = nxnp * nynp * nznp
        
        ! 2. 分配内存 (先按最大分配，后压缩)
        allocate(LinearCodes(max_possible_nodes))
        allocate(LinearDensity(max_possible_nodes))
        allocate(LinearLevels(max_possible_nodes))
        allocate(n_linear_leaves_all(0:nproc-1))
        
        n_linear_leaves = 0
        max_dim = max(nxnp, max(nynp, nznp))
        max_tree_level = nint(log(real(max_dim))/log(2.0))
        
        ! 3. 开始递归构建
        ! 从局部坐标 (1,1,1) 开始，大小为 max_dim (max_dim 为 2 的幂)
        ! level 从 0 开始
        call RecursiveBuild(1, 1, 1, max_dim, 0, amr_var)

        call MPI_AllReduce(n_linear_leaves, n_linear_leaves_max, 1, &
                           MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
        call MPI_AllReduce(n_linear_leaves, n_linear_leaves_total, 1, &
                           MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
        call MPI_AllGather(n_linear_leaves, 1, MPI_INTEGER, &
                           n_linear_leaves_all, 1, MPI_INTEGER, &
                           MPI_COMM_WORLD, ierror)

        if (nrank == 0) then
            print *, "--- Linear Octree Build Stats ---"
            print *, "Total Leaves across all CPUs :", n_linear_leaves_total
            print *, "Max Leaves on single CPU     :", n_linear_leaves_max
            print *, "Compression Ratio (Global)   :", &
                     real(n_linear_leaves_total) / real(nxc*nyc*nzc)
        endif                
        
        ! 4. 内存压缩 (Resize)
        call ResizeArrays()
        
    end subroutine BuildLinearOctree

    recursive subroutine RecursiveBuild(gx, gy, gz, max_dim, level, amr_var)
        integer, intent(in) :: gx, gy, gz, max_dim, level
        real, intent(in) :: amr_var(:,:,:)
        
        real(RK) :: local_min, local_max, local_avg
        logical :: is_leaf
        integer :: half_size
        integer(i8b) :: m_code
        type(box) :: query_box
        integer :: global_ix, global_iy, global_iz
        
        ! 1. 构造当前 Block 的物理包围盒 (Global Physical Coordinates)
        ! IID, JID, KID 是 MPI 全局偏移量 
        ! gx, gy, gz 是 1-based 的局部索引
        
        ! 全局MPI区域起始坐标索引 (0-based 用于计算坐标)
        global_ix = (gx - 1) + IID * nxnp
        global_iy = (gy - 1) + JID * nynp
        global_iz = (gz - 1) + KID * nznp
        
        ! 计算物理坐标
        query_box%min(1) = real(global_ix, RK) * dx
        query_box%min(2) = real(global_iy, RK) * dy
        query_box%min(3) = real(global_iz, RK) * dz
        
        query_box%max(1) = real(global_ix + max_dim, RK) * dx
        query_box%max(2) = real(global_iy + max_dim, RK) * dy
        query_box%max(3) = real(global_iz + max_dim, RK) * dz
        
        ! 2. 调用现有的 box_avg_density 获取统计信息
        ! 注意：传入 IID, JID, KID 以便函数内部能映射回局部索引
        call box_avg_density(query_box, IID, JID, KID, &
                             amr_var, local_avg, local_min, local_max)
        
        ! 3. 决策：是否合并为叶子？
        is_leaf = .false.
        
        if (max_dim == 1) then
            is_leaf = .true.
        else 
            ! 判据 A: 真空合并，不能用avg要用max，因为平均值会被真空稀释
            if (local_max < REFINE_CRITERIA) then
                is_leaf = .true.
            ! 判据 B: 均匀性合并 (如果 Max/Min 差异不大)
            else
                if ((local_max - local_min) / local_max < UNIFORM_CRITERIA) then
                    is_leaf = .true.
                endif
            endif
        endif
        
        ! 4. 存储或分裂
        if (is_leaf) then
            n_linear_leaves = n_linear_leaves + 1
            if (n_linear_leaves > size(LinearCodes)) stop "Octree Overflow"
            
            ! 这里我们将MPI区域内网格转换为 Morton 码
            m_code = EncodeMorton3D(gx-1, gy-1, gz-1)
            
            LinearCodes(n_linear_leaves)   = m_code
            LinearDensity(n_linear_leaves) = local_avg
            LinearLevels(n_linear_leaves)  = level
            
        else
            ! 分裂 (Z-Curve Order 0..7)
            half_size = max_dim / 2
            
            call RecursiveBuild(gx, gy, gz, half_size, level + 1, amr_var)
            call RecursiveBuild(gx + half_size, gy, gz, half_size, level + 1, amr_var)
            call RecursiveBuild(gx, gy + half_size, gz, half_size, level + 1, amr_var)
            call RecursiveBuild(gx + half_size, gy + half_size, gz, half_size, level + 1, amr_var)
            call RecursiveBuild(gx, gy, gz + half_size, half_size, level + 1, amr_var)
            call RecursiveBuild(gx + half_size, gy, gz + half_size, half_size, level + 1, amr_var)
            call RecursiveBuild(gx, gy + half_size, gz + half_size, half_size, level + 1, amr_var)
            call RecursiveBuild(gx + half_size, gy + half_size, gz + half_size, half_size, level + 1, amr_var)
        endif
        
    end subroutine RecursiveBuild

    subroutine ResizeArrays()
        integer(i8b), allocatable :: temp_c(:)
        integer, allocatable :: temp_l(:)
        real(RK), allocatable :: temp_d(:)
        
        allocate(temp_c(n_linear_leaves_max))
        allocate(temp_l(n_linear_leaves_max))
        allocate(temp_d(n_linear_leaves_max))

        ! [关键注意] 只拷贝本地有效的数据 (1 : n_linear_leaves)
        ! 剩下的部分 (n_linear_leaves+1 : n_linear_leaves_max) 是垃圾值或0，不用管
        if (n_linear_leaves > 0) then
            temp_c(1:n_linear_leaves) = LinearCodes(1:n_linear_leaves)
            temp_l(1:n_linear_leaves) = LinearLevels(1:n_linear_leaves)
            temp_d(1:n_linear_leaves) = LinearDensity(1:n_linear_leaves)
        endif

        ! 为了安全，可以将多余部分初始化为 0
        if (n_linear_leaves < n_linear_leaves_max) then
             temp_c(n_linear_leaves+1:) = 0
             temp_l(n_linear_leaves+1:) = -1
             temp_d(n_linear_leaves+1:) = 0.0_RK
        endif
        
        call move_alloc(temp_c, LinearCodes)
        call move_alloc(temp_l, LinearLevels)
        call move_alloc(temp_d, LinearDensity)
    end subroutine ResizeArrays

!#################################################
!    DDA ray traversal based on linear octree    !
!#################################################
subroutine RayTheia_Linear_DDA(ray, current_box, n_tree_nodes, temp_codes, temp_levels, ipix, epray, projected, plength)
        implicit none
        type(HEALPix_ray), intent(in) :: ray
        type(box), intent(in) :: current_box
        integer, intent(in) :: n_tree_nodes
        integer(i8b), intent(in) :: temp_codes(n_linear_leaves_max)
        integer, intent(in) :: temp_levels(n_linear_leaves_max)
        integer, intent(in) :: ipix
        integer, intent(inout) :: epray
        integer, intent(inout) :: projected(0:maxpoints)
        real(RK), intent(inout) :: plength(0:maxpoints)

        ! locals
        type(slab) :: ray_xyz
        type(box) :: node_box
        
        ! 时间/距离变量
        real(RK) :: t_curr, t_next
        real(RK) :: tmin, tmax, length
        real(RK) :: t_node_min, t_node_max, node_len
        
        ! 几何与索引变量
        real(RK) :: pos(3), local_pos(3)
        real(RK) :: d_inv(3)
        integer :: lx, ly, lz, m
        integer(i8b) :: search_code
        integer :: idx, lvl
        
        ! 记录辅助
        integer :: node_lx, node_ly, node_lz, node_grid_size
        integer :: id
        real(RK), parameter :: EPS = 1.0e-5_RK

        ! ---------------------------------------------------
        ! 1. 初始化射线 Slab 几何
        ! ---------------------------------------------------
        ray_xyz%origin = ray%origin
        ray_xyz%dir(1) = sin(ray%angle(1)) * cos(ray%angle(2))
        ray_xyz%dir(2) = sin(ray%angle(1)) * sin(ray%angle(2))
        ray_xyz%dir(3) = cos(ray%angle(1))
        do m = 1, 3
            if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
            else
                ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
            end if
        end do
        
        ! ---------------------------------------------------
        ! 2. [MPI核心] 确定光线在当前进程区域 (current_box) 的进出点
        ! ---------------------------------------------------
        ! 我们只关心光线在 current_box 内部的这一段
        length = 0.0
        call intersections(ray_xyz, current_box, tmin, tmax, length)
        if(length <= 0.D0 .or. tmax < 0.D0) return
        
        ! 初始化追踪起点 t_curr
        ! 如果光线起点在 current_box 内部，则从 t_curr = 0 开始追踪
        ! 如果光线起点在 current_box 外部，则从进入点 t_curr = tmin 开始追踪
        t_curr = max(tmin, 0.D0)
        t_curr = t_curr + EPS

        ! ---------------------------------------------------
        ! 3. DDA 步进循环 (仅在 current_box 范围内)
        ! ---------------------------------------------------
        do while (t_curr < tmax)
            
            ! A. 计算 全局物理位置 P_global
            pos = ray_xyz%origin + t_curr * ray_xyz%dir
            
            ! B. [MPI关键] 转换为 局部逻辑坐标 (Local Logical Coord)
            ! 全局位置 - 本进程左下角坐标 = 本进程内的相对偏移
            local_pos = (pos - current_box%min)
            
            ! 转换为 Grid Index (0-based)
            lx = floor(local_pos(1) / dx)
            ly = floor(local_pos(2) / dy)
            lz = floor(local_pos(3) / dz)
            
            ! 边界截断保护 (Clamp)
            ! 防止 EPS 误差导致坐标变成 -1
            if (lx < 0) lx = 0
            if (ly < 0) ly = 0
            if (lz < 0) lz = 0
            ! 注意：不用检查上限，因为我们限制在 t_domain_exit 内，
            ! 且 EncodeMorton/FindLeafIndex 会处理找不到的情况
            
            ! C. 生成局部 Morton 码并查表
            search_code = EncodeMorton3D(lx, ly, lz)
            call FindLeafIndex(search_code, n_tree_nodes, temp_codes, idx)
            
            if (idx > 0) then
                ! === 命中局部AMR节点 ===
                ! D. 构建该节点的 全局物理包围盒 (node_box)
                lvl = temp_levels(idx)
                call DecodeMorton3D(temp_codes(idx), node_lx, node_ly, node_lz)
                ! 1. 计算节点大小 (物理单位)
                node_grid_size = ishft(1, max_tree_level - lvl) 
                
                ! 2. 映射回全局物理坐标:
                ! Global_Node_Min = Box1_Min (MPI偏移) + Local_Offset * dx
                node_box%min(1) = current_box%min(1) + real(node_lx, RK) * dx
                node_box%min(2) = current_box%min(2) + real(node_ly, RK) * dy
                node_box%min(3) = current_box%min(3) + real(node_lz, RK) * dz
                
                node_box%max(1) = node_box%min(1) + real(node_grid_size, RK) * dx
                node_box%max(2) = node_box%min(2) + real(node_grid_size, RK) * dy
                node_box%max(3) = node_box%min(3) + real(node_grid_size, RK) * dz
                
                ! E. 计算光线穿出该节点的时间 (t_node_max)
                ! 注意：intersections 是纯几何运算，支持全局坐标
                call intersections(ray_xyz, node_box, t_node_min, t_node_max, node_len)
                
                ! 确定下一步位置
                ! 正常情况：t_next = t_node_max
                if (t_node_max > t_curr) then
                    t_next = t_node_max
                else
                    ! 容错：如果计算出的出口比当前还小（浮点误差），强制推进一步
                    t_next = t_curr + min(dx, min(dy, dz))
                end if
                
                ! F. 记录
                epray = epray + 1
                id = epray
                projected(id) = idx
                plength(id) = node_len             
                
                ! G. 推进
                t_curr = t_next
                
            else
                ! === 未找到节点 (异常或Padding外溢) ===
                ! 在 MPI 边界处，八叉树可能为了补齐 2^N 而比 current_box 稍微大一点，
                ! 或者光线刚好处在 current_box 边缘。
                ! 动作：推进一步，尝试重新定位
                t_curr = t_curr + min(dx, min(dy, dz))
            endif
            
            ! H. 跨越边界
            t_curr = t_curr + EPS
            
        end do
        
    end subroutine RayTheia_Linear_DDA
    
    ! Binary search leaf nodes
    subroutine FindLeafIndex(key, n_search_size, temp_codes, idx)
        integer(i8b), intent(in) :: key
        integer, intent(in) :: n_search_size
        integer(i8b), intent(in) :: temp_codes(n_linear_leaves_max)
        integer, intent(out) :: idx

        ! locals
        integer :: left, right, mid
        integer(i8b) :: code_mid
        
        idx = -1
        if (n_search_size == 0) return
        
        left = 1
        right = n_search_size
        
        ! 查找满足 LinearCodes(i) <= key 的最大索引
        do while (left <= right)
            mid = (left + right) / 2
            code_mid = temp_codes(mid)
            if (code_mid <= key) then
                idx = mid
                left = mid + 1
            else
                right = mid - 1
            end if
        end do
    end subroutine FindLeafIndex

end module m_Raytheia

