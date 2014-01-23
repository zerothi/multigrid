module t_grid_box

  implicit none

  integer, parameter :: dp = selected_real_kind(p=15)
  integer, parameter :: grid_p = selected_real_kind(p=6)

  type :: mg_box
     sequence
     ! this contains the min/max indices for the constant region
     ! xmin/xmax
     ! ymin/ymax
     ! zmin/zmax
     integer :: place(2,3) = 0
     real(grid_p) :: val = 0._grid_p
     logical :: constant = .false.
  end type mg_box

contains

  subroutine mv_xa2uc(na_u,xa,xa_u)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    real(dp), intent(out) :: xa_u(3,na_u)
    
    
    

  subroutine crt_box(na_u,xa,cell,

  subroutine init_box(box, cell, , layer)
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
