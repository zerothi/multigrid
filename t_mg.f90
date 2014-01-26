module t_mg

  use t_grid_box

  implicit none

  ! a module to sustain a "simple" multi-grid solver using the 

  type :: mg_grid
     ! the grid information
     integer :: n(3) ! size in each direction
     real(grid_p) :: sor ! the SOR value
     real(grid_p) :: a(3) ! the pre-factors for the summation
     real(grid_p) :: tol ! the tolerance of the current grid
                         ! this allows different tolerances for different layer-grids
     integer :: itt ! iterations needed to converge...
     integer :: layer ! the layer that this grid resides in
     real(grid_p), pointer :: V(:) => null() ! black/red update array
     real(grid_p), pointer :: g(:) => null() ! ghost arrays ( one long array for all bounds )
     real(grid_p), pointer :: g_s(:) => null() ! send ghost arrays ( one long array for all bounds )
     type(mg_grid), pointer :: parent => null()
     type(mg_grid), pointer :: child => null()
     integer :: N_box
     type(mg_box), allocatable :: box(:)
  end type mg_grid

contains

  subroutine init_grid(grid, n1, n2, n3, tol, cell, layer)
    type(mg_grid), intent(inout) :: grid
    integer,  intent(in) :: n1, n2, n3
    real(grid_p), intent(in) :: tol
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: layer
    real(dp), intent(in) :: celll(3)
    real(dp) :: tmp
    integer :: i

    ! ensure it is empty
    call delete_grid(grid)

    ! create the ax, ay, az pre-factors for the 3D Poisson solver
    do i = 1 , 3
       celll(i) = &
            cell(1,i) ** 2 + &
            cell(2,i) ** 2 + &
            cell(3,i) ** 2
    end do

    celll(1) = celll(1) / n1
    celll(2) = celll(2) / n2
    celll(3) = celll(3) / n3

    grid%layer = layer
    grid%n1 = n1
    grid%n2 = n2
    grid%n3 = n3
    grid%itt = 0
    if ( present(tol) ) then
       grid%tol = tol
    else if ( associated(grid%parent) ) then
       grid%tol = grid%parent%tol
    else
       stop 
    end if

    ! the SOR parameter
    grid%sor = 2._grid_p / (1._grid_p + 3.141592635_grid_p / max(n1,n2,n3) )

    tmp = 1._dp / ( 2._dp * sum(celll) ) 
    grid%ax = celll(2) * celll(3) * tmp
    grid%ay = celll(1) * celll(3) * tmp
    grid%az = celll(1) * celll(2) * tmp

  end subroutine init_grid

  subroutine init_grid_child(grid)
    type(mg_grid), intent(inout) :: grid
    < do the precontraction >
  end subroutine init_grid_child

  subroutine init_grid_parent(grid)
    type(mg_grid), intent(inout) :: grid
    < do the prolongation >
  end subroutine init_grid_parent

  subroutine grid_hold_back(grid)
    type(mg_grid), intent(inout) :: grid
    
    if ( associated(grid%V)   ) deallocate(grid%V)
    if ( associated(grid%g)   ) deallocate(grid%g)
    if ( associated(grid%g_s) ) deallocate(grid%g_s)

    nullify(grid%V,grid%g,grid%g_s) ! ensure nullification

  end subroutine grid_hold_back

  subroutine grid_bring_back(grid)
    type(mg_grid), intent(inout) :: grid
    
    allocate(grid%V(n1*n2*n3))
    allocate(grid%g(n1*2+n2*2+n3*2))
    allocate(grid%g_s(n1*2+n2*2+n3*2))

  end subroutine grid_bring_back

  recursive subroutine delete_grid(grid)
    type(mg_grid), intent(inout) :: grid

    if ( associated(grid%child) ) then
       call delete_grid(grid%child)
       deallocate(grid%child)
       nullify(grid%child)
    end if

    call grid_hold_back(grid)

  end subroutine delete_grid

  subroutine grid_restriction(grid)
    type(mg_grid), intent(inout) :: grid

    type(mg_grid), intent(inout), target :: grid

    real(grid_p),  pointer :: V(:,:,:), Vc(:,:,:)
    type(mg_grid), pointer :: child

    real(grid_p), parameter :: f1 = 1._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f4 = 4._grid_p / 64._grid_p
    real(grid_p), parameter :: f8 = 8._grid_p / 64._grid_p

    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%child) ) return

    call from1dto3d(grid%n1 ,grid%n2 ,grid%n3 ,grid%V ,V )
    child => grid%child
    call from1dto3d(child%n1,child%n2,child%n3,child%V,Vc)

    ! initialize the child
    Vc = 0._grid_p

    ! we employ full-weighting

    do z = 2 , grid%n3 - 1 , 2
    pz = z / 2 
    do y = 2 , grid%n2 - 1 , 2
    py = y / 2
    do x = 2 , grid%n1 - 1 , 2
       px = x / 2

       ! corners
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y-1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y-1,z+1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y+1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y+1,z+1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y-1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y-1,z+1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y+1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y+1,z+1) * f1

       ! middles
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y-1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y+1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y,z+1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y-1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y+1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y,z+1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y-1,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y-1,z+1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y+1,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y+1,z+1) * f2

       ! neighbours
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y-1,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y+1,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y,z-1) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y,z+1) * f4

       ! center
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y,z) * f8

    end do
    end do
    end do

    ! we still need the border...

  end subroutine grid_restriction

  subroutine grid_prolongation(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:)
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f2 =   .5_grid_p ! 1 / 2
    real(grid_p), parameter :: f4 =  .25_grid_p ! 1 / 4
    real(grid_p), parameter :: f8 = .125_grid_p ! 1 / 8

    real(grid_p) :: v2, v4, v8

    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    call from1dto3d(grid%n1  ,grid%n2  ,grid%n3  ,grid%V  ,V )
    parent => grid%parent
    call from1dto3d(parent%n1,parent%n2,parent%n3,parent%V,Vp)

    ! initialize the parent
    Vp = 0._grid_p

    do z = 1 , grid%n3
    pz = 2 * z
    do y = 1 , grid%n2
    py = 2 * y
    do x = 1 , grid%n1

       px = 2 * x

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
       Vp(px,py,pz) = V(x,y,z)

    end do
    end do
    end do

    ! we still need the border

  end subroutine grid_prolongation

  function is_constant(grid,x,y,z) result(is)
    type(mg_grid), intent(in) :: grid
    integer, intent(in) :: x,y,z
    
    integer :: i
    
    do i = 1 , grid%N_box
       if ( grid%box(i)%box%constant ) then
          if ( in_box(grid%box(i)%box,x,y,z) ) then
             is = .true.
             return
          end if
       end if
    end do
    is = .false.

  end function is_constant

  function val_rho(grid,x,y,z) result(val)
    type(mg_grid), intent(in) :: grid
    integer, intent(in) :: x,y,z
    real(grid_p) :: val
    
    integer :: i
    
    do i = 1 , grid%N_box
       if ( grid%box(i)%box%constant ) then
          if ( in_box(grid%box(i)%box,x,y,z) )then
             val = grid%box(i)%box%val
             return
          end if
       end if
    end do
    ! all values are defaulted to 1
    val = 1._grid_p

  end function val_rho

end module t_mg

subroutine from1dto3d(n1,n2,n3,V1D,V3D)
  use t_mg, only : grid_p
  integer,   intent(in) :: n1, n2, n3
  real(grid_p),  target :: V1D(n1,n2,n3)
  real(grid_p), pointer :: V3D(:,:,:)
  V3D => V1D
end subroutine from1dto3d
