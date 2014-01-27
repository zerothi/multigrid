module t_grid_box

  implicit none

  integer, parameter :: dp = selected_real_kind(p=15)
  integer, parameter :: grid_p = selected_real_kind(p=6)

  type :: mg_box_
     sequence
     ! this contains the min/max indices for the constant region
     ! xmin/xmax
     ! ymin/ymax
     ! zmin/zmax
     integer :: place(2,3) = 0
     real(grid_p) :: val = 1._grid_p ! *MUST be above 1*
     logical :: constant = .false.
  end type mg_box_

  type :: mg_box
     type(mg_box_), pointer :: box => null()
  end type mg_box

contains

  subroutine init_box(grid, box, llc, box_cell, val, constant)
    type(mg_grid), intent(in) :: grid
    type(mg_box), intent(inout) :: box
    real(dp), intent(in) :: llc(3), box_cell(3,3)
    real(grid_p), intent(in) :: val
    logical, intent(in) :: constant

    integer,  intent(in) :: n1, n2, n3
    real(grid_p), intent(in) :: tol
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: layer
    real(dp), intent(in) :: celll(3)
    real(dp) :: tmp
    integer :: i

    ! ensure it is empty
    call delete_box(box)

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

  subroutine delete_box(box)
    type(mg_box), intent(inout) :: box
    if ( associated(box%box) ) then
       deallocate(box%box)
       nullify(box%box)
    end if
  end subroutine delete_box

  pure function in_box(box,x,y,z) result(in)
    type(mg_box), intent(in) :: box
    integer, intent(in) :: x,y,z
    in = x <= box%box%place(1,1) .and. &
         box%box%place(2,1) <= x .and. &
         y <= box%box%place(1,2) .and. &
         box%box%place(2,2) <= y .and. &
         z <= box%box%place(1,3) .and. &
         box%box%place(2,3) <= z
  end function in_box

end module t_grid_box
