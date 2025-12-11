!Written by Zhengping Zhu
module m_Ray_box
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

    type, public :: amr_node
        type(box) :: bbox
        real :: avg_rho
        integer :: child(0:7)
    end type amr_node
    type(amr_node), public, allocatable :: tree(:)
    integer, public :: n_nodes

    integer,public :: levels
    integer,public,allocatable :: maxpoints_ray(:)
    type(box),public :: box1
    type(HEALPix_ray),public :: ray
    real(kind=dp),public :: thfpix, phfpix, contribution
    real(kind=dp),public :: corner_min(3), corner_max(3)

    real(RK), parameter :: DENSITY_THRESHOLD = 1.D-8

    public:: intersections, box_avg_density, raytheia_sm, build_tree, raytheia_amr
contains
    logical function isfinite(x)
            real(RK), intent(in) :: x
            isfinite = (abs(x) < huge(x) .and. x == x) 
    end function isfinite

    subroutine split_box(parent, children)
    implicit none
    type(box),intent(in) :: parent
    type(box),intent(out) :: children(0:7)
    real(RK) :: mid(3)
    integer :: i, d

    mid = (parent%min + parent%max) / 2.D0

    do i = 0, 7
        do d = 1, 3
            if (btest(i, d-1)) then
                children(i)%min(d) = mid(d)
                children(i)%max(d) = parent%max(d)
            else
                children(i)%min(d) = parent%min(d)
                children(i)%max(d) = mid(d)
            endif
        enddo
    enddo

    end subroutine split_box

    subroutine intersections(ray_xyz, box_in, length)
    implicit none
    type(slab), intent(in) :: ray_xyz
    type(box), intent(in) :: box_in
    real(RK), intent(inout) :: length
    integer :: i, d
    real(RK) :: tmin, tmax, t1, t2
    real(RK) :: Pmin(3), Pmax(3)
    real(RK) :: distance, diff
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

    ! 判断射线是否与盒子相交
    if (tmin <= tmax) then
        ! 计算交点坐标
        do d = 1, 3
            Pmin(d) = ray_xyz%origin(d) + tmin / ray_xyz%dir_inv(d)
            Pmax(d) = ray_xyz%origin(d) + tmax / ray_xyz%dir_inv(d)
        end do

        ! 计算欧几里得距离
        distance = 0.0
        do d = 1, 3
            diff = Pmax(d) - Pmin(d)
            distance = distance + diff**2
        end do
        length = sqrt(distance)  ! 射线穿过盒子的欧几里得距离
    else
        length = 0.0  ! 如果不相交，长度为 0
    end if

    end subroutine intersections

    recursive subroutine raytheia_sm(ray, parent, level, pid, pjd, pkd, epray, projected, plength) ! save memory
        type(HEALPix_ray), intent(in) :: ray
        type(slab) :: ray_xyz
        type(box), intent(in) :: parent
        type(box) :: children(0:7)
        integer, intent(in) :: level, pid, pjd, pkd
        integer :: epray
        integer :: projected(0:maxpoints,3)
        real(RK) :: plength(0:maxpoints)
        integer :: i, j, k, II, JJ, KK, d, m, parent_index, ir, node_count, start_index, cI, id
        real(RK) :: mid(3), extent(3), center(3)
        logical :: intersect
        real(RK) :: x, y, z, r, xnode, ynode, znode
!        real(RK) :: Aij,c,nu,pc2cm,f,g_i,g_j,thfpix,phfpix,length
        real(RK) :: thfpix,phfpix,length

        if (level > 0) then
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)

            intersect = .false.
            if (length > 0.0) then
                intersect = .true.
            endif

            if(intersect) then
                call split_box(parent, children)

                do i = 0, 7
                    call raytheia_sm(ray, children(i), level-1, pid, pjd, pkd, epray, projected, plength)
                enddo
            endif
        endif

        ! leaf nodes calculations
        if(level == 0) then

            II = nint(parent%max(1)/dx)
            JJ = nint(parent%max(2)/dy)
            KK = nint(parent%max(3)/dz)

            i=II-pid*nxnp
            j=JJ-pjd*nynp
            k=KK-pkd*nznp

            ! penetration length
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)
            if(length.ne.0.D0) then
                epray = epray + 1
                id = epray
                projected(id,1) = i
                projected(id,2) = j
                projected(id,3) = k
                plength(id) = length
            endif

        endif

    end subroutine raytheia_sm

!###############################
!         AMR routines         !
!###############################

    subroutine box_avg_density(box_in, pid, pjd, pkd, rho, avg_rho)
    implicit none
        type(box),intent(in) :: box_in
        integer, intent(in) :: pid, pjd, pkd
        real, intent(in) :: rho(:,:,:)
        real(RK), intent(out) :: avg_rho 

        integer :: i1, i2, j1, j2, k1, k2 !local indices
        integer :: ig1, ig2, jg1, jg2, kg1, kg2 ! global indices
        real(RK) :: vol, total_rho

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

        ! calculte average density
        total_rho = DBLE(sum(rho(i1:i2, j1:j2, k1:k2)))
        vol = DBLE((i2-i1+1) * (j2-j1+1) * (k2-k1+1))
        avg_rho = total_rho / vol
 
    end subroutine box_avg_density

    subroutine build_tree(parent_box, pid, pjd, pkd, level, rho)
        type(box), intent(in) :: parent_box
        real, intent(in) :: rho(:,:,:)
        integer, intent(in) :: pid, pjd, pkd, level

        integer :: root_id

        n_nodes = 0
        call add_node(parent_box, pid, pjd, pkd, level, rho, root_id)
        ! print*,nrank,'n_nodes',n_nodes

    end subroutine build_tree

    recursive subroutine add_node(current_box, pid, pjd, pkd, level, rho, node_id)
        type(box), intent(in) :: current_box
        real, intent(in) :: rho(:,:,:)
        integer, intent(in) :: pid, pjd, pkd, level
        integer, intent(out) :: node_id

        real(RK) :: avg_rho
        type(box) :: children(0:7)
        integer :: i, child_id

        n_nodes = n_nodes + 1
        if(n_nodes > size(tree)) stop 'AMR tree overflow'

        node_id = n_nodes

        tree(node_id)%bbox = current_box
        call box_avg_density(current_box, pid, pjd, pkd, rho, avg_rho)
        tree(node_id)%avg_rho = avg_rho

        if(level > 0 .and. avg_rho >= DENSITY_THRESHOLD) then
            call split_box(current_box, children)
            do i = 0, 7
                call add_node(children(i), pid, pjd, pkd, level - 1, rho, child_id)
                tree(node_id)%child(i) = child_id
            enddo
        else
            tree(node_id)%child(:) = -1 ! no more refine
        endif

    end subroutine add_node

    recursive subroutine raytheia_amr(ray, node_id, rho, integral)
        type(HEALPix_ray), intent(in) :: ray
        integer, intent(in) :: node_id
        real, intent(in) :: rho(:,:,:)
        real(RK), intent(inout) :: integral

        type(slab) :: ray_xyz
        integer :: i, m
        logical :: intersect
        real(RK) :: thfpix,phfpix,length,avg_rho

        ray_xyz%origin = ray%origin
        thfpix = ray%angle(1)
        phfpix = ray%angle(2)
        ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
        ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
        ray_xyz%dir(3) = cos(thfpix)
        do m = 1, 3
            if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
            else
                ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
            end if
        end do

        length = 0.0
        call intersections(ray_xyz, tree(node_id)%bbox, length)
        intersect = (length > 0.D0)

        if (.not. intersect) return ! not intersect

        if (tree(node_id)%child(0) == -1) then
            integral = integral + tree(node_id)%avg_rho * length * pc
        else
            do i = 0, 7
                if (tree(node_id)%child(i) > 0) then
                    call raytheia_amr(ray, tree(node_id)%child(i), &
                                        rho, integral)
                end if
            end do        
        end if

    end subroutine raytheia_amr

end module m_Ray_box

