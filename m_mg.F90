module m_mg

  use t_mg
  use m_gs_CDS
  !use m_gs_br

  implicit none

  integer :: method = MG_METHOD_GS_TEMPORAL_CDS

  ! the top multi-grid type
  type(mg_grid), target :: top

contains

  subroutine solve_MultiGrid_poison()

    type(mg_grid), pointer :: grid, tmp_grid

    real(grid_p) :: N_frac
    integer :: nn1, nn2, nn3

    !< initialize communication scheme &
    !     each grid has the same neighbours and as such we dont &
    !     need separate schemes for each grid >

    ! This estimates the needed number of grids
    N_frac = log(2._grid_p) / 2._grid_p

    call init_grid(top,n1,n2,n3)
    call new_grid_size(top,nn1,nn2,nn3,N)
    ! create all the grids
    grid => top
    do while ( N > 200 ) 
       nullify(tmp_grid)
       allocate(tmp_grid)
       grid%child => tmp_grid
       tmp_grid%parent => grid
       grid => tmp_grid
       call init_grid(grid,nn1,nn2,nn3)
       call new_grid_size(grid,nn1,nn2,nn3,N)
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

    subroutine new_grid_size(grid,nn1,nn2,nn3,N)
      type(mg_grid), intent(in) :: grid
      integer, intent(out) :: nn1,nn2,nn3,N
      nn1 = grid%n(1) * N_frac
      nn2 = grid%n(2) * N_frac
      nn3 = grid%n(3) * N_frac
      N = nn1 * nn2 * nn3
    end subroutine new_grid_size

  end subroutine solve_MultiGrid_poison

end module m_mg
