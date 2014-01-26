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
     real(grid_p) :: val = 1._grid_p ! *MUST be above 1*
     logical :: constant = .false.
  end type mg_box

contains

  pure function in_box(box,x,y,z) result(in)
    type(mg_box), intent(in) :: box
    integer, intent(in) :: x,y,z
    in = x <= grid%box(i)%box%place(1,1) .and. &
         grid%box(i)%box%place(2,1) <= x .and. &
         y <= grid%box(i)%box%place(1,2) .and. &
         grid%box(i)%box%place(2,2) <= y .and. &
         z <= grid%box(i)%box%place(1,3) .and. &
         grid%box(i)%box%place(2,3) <= z
  end function in_box


  subroutine init_box(grid, box, cell, layer)
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
    
  end subroutine init_box

end module t_mg
