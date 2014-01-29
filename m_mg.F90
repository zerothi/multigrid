module m_mg

  use t_mg
  use m_gs_CDS
  !use m_gs_br

  implicit none

  integer :: method = MG_METHOD_GS_TEMPORAL_CDS

  ! the top multi-grid type
  type(mg_grid), target :: top

contains

  subroutine solve_MultiGrid_poison(top)

    type(mg_grid), pointer :: top

    ! local variables
    type(mg_grid), pointer :: grid, tmp_grid
    integer :: N, nn(3)

    !< initialize communication scheme &
    !     each grid has the same neighbours and as such we dont &
    !     need separate schemes for each grid >

    call new_grid_size(top,nn,N)

    ! create all the grids
    grid => top
    do while ( N > 200 ) 
       nullify(tmp_grid)
       allocate(tmp_grid)
       grid%child => tmp_grid
       tmp_grid%parent => grid
       call init_grid(tmp_grid,nn, &
            grid%cell, grid%layer+1, grid%N_box, &
            grid%tol, grid%offset)
       grid => tmp_grid
       call new_grid_size(grid,nn,N)
    end do
    nullify(grid,tmp_grid)

    ! solve the system
    if ( method == MG_METHOD_GS_TEMPORAL_CDS ) then
       call mg_gs_cds(top)
    !else if ( method == GS_BLACK_RED ) then
       !call mg_gs_br(top)
    end if
    
    !< transfer solution to the Vscf array > 

    ! this will recursively delete all grids
    call delete_grid(top)

  contains

    subroutine new_grid_size(grid,n,NN)
      type(mg_grid), intent(in) :: grid
      integer, intent(out) :: n(3), NN
      integer :: i
      
      do i = 1 , 3
         n(i) = (grid%n(i) + 1) / 2
         if ( n(i) < 7 ) then
            n = 0
            NN = 0
            return
         end if
      end do

      NN = n(1) * n(2) * n(3)

    end subroutine new_grid_size

  end subroutine solve_MultiGrid_poison

end module m_mg
