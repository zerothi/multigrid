module t_mg_interp

  use t_mg

  implicit none

contains

  subroutine grid_restriction(grid)
    type(mg_grid), intent(inout) :: grid

    select case ( grid%RES_PRO_method )
    case ( 1 ) 
       call grid_restriction_full(grid)
    case ( 2 )
       call grid_restriction_half(grid)
    case default
       call grid_restriction_full(grid)
    end select
  end subroutine grid_restriction

  subroutine grid_prolongation(grid)
    type(mg_grid), intent(inout) :: grid

    select case ( grid%RES_PRO_method )
    case ( 1 ) 
       call grid_prolongation_full(grid)
    case ( 2 )
       call grid_prolongation_half(grid)
    case default
       call grid_prolongation_full(grid)
    end select
  end subroutine grid_prolongation

  subroutine grid_restriction_new(grid)
    type(mg_grid), intent(inout) :: grid

    real(grid_p),  pointer :: V(:,:,:), Vc(:,:,:)
    type(mg_grid), pointer :: child

    real(grid_p), parameter :: f1 = 10._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 5._grid_p / 64._grid_p
    real(grid_p), parameter :: f3 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f4 = 1._grid_p / 64._grid_p

    real(grid_p), pointer :: vt
    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%child) ) return

    V  => grid%V
    child => grid%child
    Vc => child%V

!$OMP parallel default(shared) private(x,px,y,py,z,pz,vt)

    ! initialize the child
!$OMP workshare
    Vc = 0._grid_p
!$OMP end workshare

    ! we employ full-weighting

!$OMP do
    do z = 1 , child%n(3)
    pz = z * 2
    if ( pz >= grid%n(3) ) cycle
    do y = 1 , child%n(2)
    py = y * 2
    if ( py >= grid%n(2) ) cycle
    do x = 1 , child%n(1)
       px = x * 2
       if ( px >= grid%n(1) ) cycle

       ! point to the element
       vt => Vc(x,y,z)

       vt =      V(px,py,pz) * f1

       ! neighbours
       vt = vt + V(px-1,py,pz) * f2
       vt = vt + V(px+1,py,pz) * f2
       vt = vt + V(px,py-1,pz) * f2
       vt = vt + V(px,py+1,pz) * f2
       vt = vt + V(px,py,pz-1) * f2
       vt = vt + V(px,py,pz+1) * f2

       ! next-neighbours
       vt = vt + V(px-1,py-1,pz) * f3
       vt = vt + V(px-1,py+1,pz) * f3
       vt = vt + V(px-1,py,pz-1) * f3
       vt = vt + V(px-1,py,pz+1) * f3

       vt = vt + V(px+1,py-1,pz) * f3
       vt = vt + V(px+1,py+1,pz) * f3
       vt = vt + V(px+1,py,pz-1) * f3
       vt = vt + V(px+1,py,pz+1) * f3

       vt = vt + V(px,py-1,pz-1) * f3
       vt = vt + V(px,py-1,pz+1) * f3
       vt = vt + V(px,py+1,pz-1) * f3
       vt = vt + V(px,py+1,pz+1) * f3

    end do
    end do
    end do
!$OMP end do

!$OMP end parallel

    ! we still need the border...

    ! re-instantiate the constant fields
    call grid_setup(child)

  end subroutine grid_restriction_new

  subroutine grid_restriction_full(grid)
    type(mg_grid), intent(inout) :: grid

    real(grid_p),  pointer :: V(:,:,:), Vc(:,:,:), vt
    type(mg_grid), pointer :: child

    real(grid_p), parameter :: f4 = 1._grid_p / 64._grid_p
    real(grid_p), parameter :: f3 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 4._grid_p / 64._grid_p
    real(grid_p), parameter :: f1 = 8._grid_p / 64._grid_p

    integer :: x,y,z, cx,cy,cz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%child) ) return

    V  => grid%V
    child => grid%child
    Vc => child%V

!$OMP parallel default(shared) private(x,cx,y,cy,z,cz,vt)

    ! initialize the child
!$OMP workshare
    Vc = 0._grid_p
!$OMP end workshare

    ! we employ full-weighting

!$OMP do
    do z = 2 , grid%n(3) - 1 , 2
    cz = z / 2

    do y = 2 , grid%n(2) - 1 , 2
    cy = y / 2

    do x = 2 , grid%n(1) - 1 , 2
       cx = x / 2

       vt => Vc(cx,cy,cz)

       ! corners
       vt =      V(x-1,y-1,z-1) * f4
       vt = vt + V(x-1,y-1,z+1) * f4
       vt = vt + V(x-1,y+1,z-1) * f4
       vt = vt + V(x-1,y+1,z+1) * f4 ! 4
       vt = vt + V(x+1,y-1,z-1) * f4
       vt = vt + V(x+1,y-1,z+1) * f4
       vt = vt + V(x+1,y+1,z-1) * f4
       vt = vt + V(x+1,y+1,z+1) * f4 ! 8

       ! middles
       vt = vt + V(x-1,y-1,z) * f3
       vt = vt + V(x-1,y+1,z) * f3
       vt = vt + V(x-1,y,z-1) * f3
       vt = vt + V(x-1,y,z+1) * f3 ! 4
       vt = vt + V(x+1,y-1,z) * f3
       vt = vt + V(x+1,y+1,z) * f3
       vt = vt + V(x+1,y,z-1) * f3
       vt = vt + V(x+1,y,z+1) * f3 ! 8
       vt = vt + V(x,y-1,z-1) * f3
       vt = vt + V(x,y-1,z+1) * f3
       vt = vt + V(x,y+1,z-1) * f3
       vt = vt + V(x,y+1,z+1) * f3 ! 12

       ! neighbours
       vt = vt + V(x-1,y,z) * f2
       vt = vt + V(x+1,y,z) * f2
       vt = vt + V(x,y-1,z) * f2 ! 3
       vt = vt + V(x,y+1,z) * f2
       vt = vt + V(x,y,z-1) * f2
       vt = vt + V(x,y,z+1) * f2 ! 6

       ! center
       vt = vt + V(x,y,z) * f1

    end do
    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    ! re-instantiate the constant fields
    call grid_setup(child)

  end subroutine grid_restriction_full

  subroutine grid_restriction_half(grid)
    type(mg_grid), intent(inout) :: grid

    real(grid_p),  pointer :: V(:,:,:), Vc(:,:,:), vt
    type(mg_grid), pointer :: child

    real(grid_p), parameter :: f3 = 2._grid_p / 56._grid_p
    real(grid_p), parameter :: f2 = 4._grid_p / 56._grid_p
    real(grid_p), parameter :: f1 = 8._grid_p / 56._grid_p

    integer :: x,y,z, cx,cy,cz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%child) ) return

    V  => grid%V
    child => grid%child
    Vc => child%V

!$OMP parallel default(shared) private(x,cx,y,cy,z,cz,vt)

    ! initialize the child
!$OMP workshare
    Vc = 0._grid_p
!$OMP end workshare

    ! we employ full-weighting

!$OMP do
    do z = 2 , grid%n(3) - 1 , 2
    cz = z / 2

    do y = 2 , grid%n(2) - 1 , 2
    cy = y / 2

    do x = 2 , grid%n(1) - 1 , 2
       cx = x / 2

       vt => Vc(cx,cy,cz)

       ! middles
       vt =      V(x-1,y-1,z) * f3
       vt = vt + V(x-1,y+1,z) * f3
       vt = vt + V(x-1,y,z-1) * f3
       vt = vt + V(x-1,y,z+1) * f3 ! 4
       vt = vt + V(x+1,y-1,z) * f3
       vt = vt + V(x+1,y+1,z) * f3
       vt = vt + V(x+1,y,z-1) * f3
       vt = vt + V(x+1,y,z+1) * f3 ! 8
       vt = vt + V(x,y-1,z-1) * f3
       vt = vt + V(x,y-1,z+1) * f3
       vt = vt + V(x,y+1,z-1) * f3
       vt = vt + V(x,y+1,z+1) * f3 ! 12

       ! neighbours
       vt = vt + V(x-1,y,z) * f2
       vt = vt + V(x+1,y,z) * f2
       vt = vt + V(x,y-1,z) * f2 ! 3
       vt = vt + V(x,y+1,z) * f2
       vt = vt + V(x,y,z-1) * f2
       vt = vt + V(x,y,z+1) * f2 ! 6

       ! center
       vt = vt + V(x,y,z) * f1

    end do
    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    ! re-instantiate the constant fields
    call grid_setup(child)

  end subroutine grid_restriction_half


  subroutine grid_prolongation_full(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:), vt
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f4 = 1._grid_p / 64._grid_p
    real(grid_p), parameter :: f3 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 4._grid_p / 64._grid_p
    real(grid_p), parameter :: f1 = 8._grid_p / 64._grid_p

    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! initialize the parent to zero
!$OMP parallel default(shared) private(x,px,y,py,z,pz,vt)

!$OMP workshare
    Vp = 0._grid_p
!$OMP end workshare

    ! do middle loop
!$OMP do
    do pz = 1 , parent%n(3)
    z = max(2, pz / 2)
    if ( z >= grid%n(3) ) cycle

    do py = 1 , parent%n(2)
    y = max(2,py / 2 - py / grid%n(2))
    if ( y >= grid%n(2) ) cycle

    do px = 1 , parent%n(1)
       x = max(2,px / 2 - px / grid%n(1))
       if ( x >= grid%n(1) ) cycle

       ! point to the parent point
       vt => Vp(px,py,pz)

       ! center
       vt = V(x,y,z) * f1

       ! neighbours
       vt = vt + V(x-1,y,z) * f2
       vt = vt + V(x+1,y,z) * f2
       vt = vt + V(x,y-1,z) * f2
       vt = vt + V(x,y+1,z) * f2
       vt = vt + V(x,y,z-1) * f2
       vt = vt + V(x,y,z+1) * f2

       ! next-neighbours
       vt = vt + V(x-1,y-1,z) * f3
       vt = vt + V(x-1,y+1,z) * f3
       vt = vt + V(x-1,y,z-1) * f3
       vt = vt + V(x-1,y,z+1) * f3

       vt = vt + V(x+1,y-1,z) * f3
       vt = vt + V(x+1,y+1,z) * f3
       vt = vt + V(x+1,y,z-1) * f3
       vt = vt + V(x+1,y,z+1) * f3

       vt = vt + V(x,y-1,z-1) * f3
       vt = vt + V(x,y-1,z+1) * f3
       vt = vt + V(x,y+1,z-1) * f3
       vt = vt + V(x,y+1,z+1) * f3

       ! corners
       vt = vt + V(x-1,y-1,z-1) * f4
       vt = vt + V(x-1,y-1,z+1) * f4
       vt = vt + V(x-1,y+1,z-1) * f4
       vt = vt + V(x-1,y+1,z+1) * f4 ! 4
       vt = vt + V(x+1,y-1,z-1) * f4
       vt = vt + V(x+1,y-1,z+1) * f4
       vt = vt + V(x+1,y+1,z-1) * f4
       vt = vt + V(x+1,y+1,z+1) * f4 ! 8
       
    end do
    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    call grid_setup(parent)

  end subroutine grid_prolongation_full

  subroutine grid_prolongation_half(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:), vt
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f3 = 2._grid_p / 56._grid_p
    real(grid_p), parameter :: f2 = 4._grid_p / 56._grid_p
    real(grid_p), parameter :: f1 = 8._grid_p / 56._grid_p

    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! initialize the parent to zero
!$OMP parallel default(shared) private(x,px,y,py,z,pz,vt)

!$OMP workshare
    Vp = 0._grid_p
!$OMP end workshare

    ! do middle loop
!$OMP do
    do pz = 1 , parent%n(3)
    z = max(2, pz / 2)
    if ( z >= grid%n(3) ) cycle

    do py = 1 , parent%n(2)
    y = max(2,py / 2 - py / grid%n(2))
    if ( y >= grid%n(2) ) cycle

    do px = 1 , parent%n(1)
       x = max(2,px / 2 - px / grid%n(1))
       if ( x >= grid%n(1) ) cycle

       ! point to the parent point
       vt => Vp(px,py,pz)

       ! center
       vt = V(x,y,z) * f1

       ! neighbours
       vt = vt + V(x-1,y,z) * f2
       vt = vt + V(x+1,y,z) * f2
       vt = vt + V(x,y-1,z) * f2
       vt = vt + V(x,y+1,z) * f2
       vt = vt + V(x,y,z-1) * f2
       vt = vt + V(x,y,z+1) * f2

       ! next-neighbours
       vt = vt + V(x-1,y-1,z) * f3
       vt = vt + V(x-1,y+1,z) * f3
       vt = vt + V(x-1,y,z-1) * f3
       vt = vt + V(x-1,y,z+1) * f3

       vt = vt + V(x+1,y-1,z) * f3
       vt = vt + V(x+1,y+1,z) * f3
       vt = vt + V(x+1,y,z-1) * f3
       vt = vt + V(x+1,y,z+1) * f3

       vt = vt + V(x,y-1,z-1) * f3
       vt = vt + V(x,y-1,z+1) * f3
       vt = vt + V(x,y+1,z-1) * f3
       vt = vt + V(x,y+1,z+1) * f3

    end do
    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    call grid_setup(parent)

  end subroutine grid_prolongation_half

  subroutine grid_prolongation_old(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:)
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f2 =   .5_grid_p ! 1 / 2
    real(grid_p), parameter :: f4 =  .25_grid_p ! 1 / 4
    real(grid_p), parameter :: f8 = .125_grid_p ! 1 / 8

    real(grid_p) :: v2, v4, v8

    logical :: one_larger(3)
    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! currently this prolongation only accepts
    ! grids of double +1/+2 sizes
    do x = 1 , 3
       select case ( parent%n(x) - grid%n(x) * 2 )
       case ( 1 )
          one_larger(x) = .true.
       case ( 2 )
          one_larger(x) = .false.
       case default
          stop 'Not a functioning grid step...'
       end select
    end do
          
    ! initialize the parent to zero
!$OMP parallel default(shared) private(x,px,y,py,z,pz,v2,v4,v8)

!$OMP workshare
    Vp = 0._grid_p
!$OMP end workshare

    ! ensure to set the top corners correctly

    z  = 1
    pz = 1
!$OMP do 
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( one_larger(2) ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v8
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
!       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px+1,py+1,pz,'x1')
!call print_t(px-1,py-1,pz,'x1')
!call print_t(px-1,py+1,pz,'x1')
!call print_t(px+1,py-1,pz,'x1')
!call print_t(px,py,pz,'x1')

    end do
    end do
!$OMP end do
    if ( grid%n(2)*2 + 1 /= parent%n(2) ) then
    y  = grid%n(2)
    py = parent%n(2)
!$OMP do 
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v8
!       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
!       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px-1,py-1,pz,'x1')
!call print_t(px+1,py-1,pz,'x1')
!call print_t(px,py,pz,'x1')

    end do
!$OMP end do
    end if

    y  = 1
    py = 1
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( one_larger(3) ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v8
!       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px-1,py,pz-1,'y1')
!call print_t(px-1,py,pz+1,'y1')
!call print_t(px+1,py,pz-1,'y1')
!call print_t(px+1,py,pz+1,'y1')
!call print_t(px,py,pz,'y1')

    end do
    end do
!$OMP end do

    x  = 1
    px = 1
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( one_larger(3) ) cycle
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
       py = 2 * y
       if ( one_larger(2) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
!       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
!       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v8
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px,py-1,pz-1,'z1')
!call print_t(px,py-1,pz+1,'z1')
!call print_t(px,py+1,pz-1,'z1')
!call print_t(px,py+1,pz+1,'z1')
!call print_t(px,py,pz,'z1')

    end do
    end do
!$OMP end do
    
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( one_larger(3) ) cycle
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( one_larger(2) ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py-1,pz+1) = Vp(px-1,py-1,pz+1) + v8
       Vp(px-1,py+1,pz-1) = Vp(px-1,py+1,pz-1) + v8
       Vp(px-1,py+1,pz+1) = Vp(px-1,py+1,pz+1) + v8
       Vp(px+1,py-1,pz-1) = Vp(px+1,py-1,pz-1) + v8
       Vp(px+1,py-1,pz+1) = Vp(px+1,py-1,pz+1) + v8
       Vp(px+1,py+1,pz-1) = Vp(px+1,py+1,pz-1) + v8
       Vp(px+1,py+1,pz+1) = Vp(px+1,py+1,pz+1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v4
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v4
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v4
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v4
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v4
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v4
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v4
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v2
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v2
       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v2

       ! center
       Vp(px,py,pz)   = V(x,y,z)

    end do
    end do
    end do
!$OMP end do

    z = grid%n(3)
    pz = parent%n(3)
    if ( .not. one_larger(3) ) then
!$OMP do 
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( parent%n(2) <= py + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py+1,pz-1) = Vp(px-1,py+1,pz-1) + v8
       Vp(px+1,py-1,pz-1) = Vp(px+1,py-1,pz-1) + v8
       Vp(px+1,py+1,pz-1) = Vp(px+1,py+1,pz-1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v4
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v4
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v4
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v2
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v2

    end do
    end do
!$OMP end do
    end if
!$OMP do 
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( parent%n(2) <= py + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v8
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8

       ! middles
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v8
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v8
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v8
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v8

       ! neighbours
       Vp(px,py,pz) = Vp(px,py,pz) + v8

    end do
    end do
!$OMP end do

    y = grid%n(2)
    py = parent%n(2)
    if ( y * 2 /= py ) then
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py-1,pz+1) = Vp(px-1,py-1,pz+1) + v8
       Vp(px+1,py-1,pz-1) = Vp(px+1,py-1,pz-1) + v8
       Vp(px+1,py-1,pz+1) = Vp(px+1,py-1,pz+1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v4
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v4
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v4
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v2
       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v2

    end do
    end do
!$OMP end do
    end if
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v8
!       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

    end do
    end do
!$OMP end do

    x = grid%n(1)
    px = parent%n(1)
    if ( x * 2 /= px ) then
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
       py = 2 * y
       if ( parent%n(2) <= py + 1 ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py-1,pz+1) = Vp(px-1,py-1,pz+1) + v8
       Vp(px-1,py+1,pz-1) = Vp(px-1,py+1,pz-1) + v8
       Vp(px-1,py+1,pz+1) = Vp(px-1,py+1,pz+1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v4
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v4
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v4
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v2
       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v2

    end do
    end do
!$OMP end do
    end if
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
       py = 2 * y
       if ( parent%n(2) <= py + 1 ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
!       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
!       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v8
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    ! we still need the border

    ! re-instantiate the constant fields
    call grid_setup(parent)

  contains

    subroutine print_t(x,y,z,m)
      integer, intent(in) :: x,y,z
      character(len=*), intent(in) :: m
      if ( x == 14 .and. y == 14 .and. z == 1 ) then
         print *,'we are here', m
      end if
    end subroutine print_t

  end subroutine grid_prolongation_old

end module t_mg_interp
