module t_mg

  implicit none

  integer, parameter :: dp = selected_real_kind(p=15)
  integer, parameter :: grid_p = selected_real_kind(p=6)

  ! a module to sustain a "simple" multi-grid solver using the 

  type :: mg_grid
     ! the grid information
     integer :: n1, n2, n3 ! size in each direction
     real(grid_p) :: ax, ay, az ! the pre-factors for the summation
     real(grid_p) :: tol ! the tolerance of the current grid
                         ! this allows different tolerances for different layer-grids
     integer :: itt ! iterations needed to converge...
     integer :: layer ! the layer that this grid resides in
     real(grid_p), pointer :: V(:) => null() ! black/red update array
     real(grid_p), pointer :: g(:) => null() ! ghost arrays ( one long array for all bounds )
     real(grid_p), pointer :: g_s(:) => null() ! send ghost arrays ( one long array for all bounds )
     type(mg_grid), pointer :: parent => null()
     type(mg_grid), pointer :: child => null()
  end type mg_grid

  type :: mg_constant
     type(mg_grid), pointer :: grid => null()
     ! this contains the min/max indices for the constant region
     ! xmin/xmax
     ! ymin/ymax
     ! zmin/zmax
     integer :: place(2,3) = 0
     real(grid_p) :: val = 0._grid_p
  end type mg_constant

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

    tmp = 1._dp / ( 2._dp * sum(celll) ) 
    grid%ax = celll(2) * celll(3) * tmp
    grid%ay = celll(1) * celll(3) * tmp
    grid%az = celll(1) * celll(2) * tmp

    allocate(grid%V(n1*n2*n3))
    allocate(grid%g(n1*2+n2*2+n3*2))
    allocate(grid%g_s(n1*2+n2*2+n3*2))
    
  end subroutine init_grid

  subroutine init_grid_child(grid)
    type(mg_grid), intent(inout) :: grid
    < do the precontraction >
  end subroutine init_grid_child

  subroutine init_grid_parent(grid)
    type(mg_grid), intent(inout) :: grid
    < do the prolongation >
  end subroutine init_grid_parent

  recursive subroutine delete_grid(grid)
    type(mg_grid), intent(inout) :: grid
    if ( associated(grid%child) ) then
       call delete_grid(grid%child)
       deallocate(grid%child)
       nullify(grid%child)
    end if
    if ( associated(grid%V)   ) deallocate(grid%V)
    if ( associated(grid%g)   ) deallocate(grid%g)
    if ( associated(grid%g_s) ) deallocate(grid%g_s)
    nullify(grid%V,grid%g_s) ! ensure nullification
  end subroutine delete_grid

end module t_mg

subroutine from1dto3d(n1,n2,n3,V1D,V3D)
  use t_mg, only : grid_p
  integer,   intent(in) :: n1, n2, n3
  real(grid_p),  target :: V1D(n1,n2,n3)
  real(grid_p), pointer :: V3D(:,:,:)
  V3D => V1D
end subroutine from1dto3d
