module m_gs_br

  use t_mg
  
  implicit none

  private

  integer, parameter :: MG_METHOD_GS_BR = 2

  public :: mg_gs_br
  public :: MG_METHOD_GS_BR

contains

  subroutine mg_gs_br(top_grid)
    type(mg_grid), pointer :: top_grid
    type(mg_grid), pointer :: grid
    
    grid => top_grid
    do while ( associated(grid%child) ) 
       grid => grid%child
       ! precontraction
       call init_grid_child(grid)
    end do
    
    do while ( associated(grid) ) 
       call grid_solve(grid)
       ! prolongation
       call init_grid_parent(grid)
       grid => grid%parent
    end do

  end subroutine mg_gs_br

  subroutine grid_solve(grid)
    real(grid_p), pointer :: Vb(:), Vr(:)
    real(grid_p) :: tol

    tol = grid%tol + 1._grid_p

    do while ( tol > grid%tol ) 

       ! initialize the tolerance
       tol = 0._dp
       
       < communicate red >

       call gs_black(grid,grid%V, tol)
       
       < inner loop on black >
       
       < wait communicate red >
       
       < calculate black bounds >
       
       < communicate black >

       call gs_red(grid,grid%V, tol)       

       < wait communicate black >

       < calculate red bounds >

    end do

  end subroutine grid_solve

  subroutine gs_black(grid,V,tol)
    integer, intent(in) :: n1, n2, n3
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, a(3)

    a(1) = grid%ax
    a(2) = grid%ay
    a(3) = grid%az
    
    do z = 2 , grid%n3 - 1 , 2
       do y = 2 , grid%n2 - 1 , 2
          do x = 3 - mod(z,2) , grid%n1 - 1 , 2
             ! calculate the current contribution
             vcur = val(a,grid%V,x,y,z)
             ! Calculate the tolerance
             tol = max(abs(V(x,y,z) - vcur),tol)
             V(x,y,z) = vcur
          end do
       end do
    end do
    
  end subroutine gs_black

  subroutine gs_red(grid,V,tol)
    integer, intent(in) :: n1, n2, n3
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(out) :: tol
    real(grid_p) :: vcur, a(3)

    a(1) = grid%ax
    a(2) = grid%ay
    a(3) = grid%az
    
    tol = 0._grid_p
    do z = 2 , grid%n3 - 1 , 2
       do y = 2 , grid%n2 - 1 , 2
          do x = 2 + mod(z,2) , grid%n1 - 1 , 2
             ! calculate the current contribution
             vcur = val(a,grid%V,x,y,z)
             ! Calculate the tolerance
             tol = max(abs(V(x,y,z) - vcur),tol)
             V(x,y,z) = vcur
          end do
       end do
    end do
    
  end subroutine gs_red
  
  pure function val(a,V,x,y,z) 
    real(grid_p), intent(in) :: a(3), V(:,:,:)
    integer, intent(in) :: x,y,z
    real(grid_p) :: val

    val = &
         a(1) * ( V(x-1,y,z) + V(x+1,y,z) ) + &
         a(2) * ( V(x,y-1,z) + V(x,y+1,z) ) + &
         a(3) * ( V(x,y,z+1) + V(x,y,z+1) )

  end function val

end module m_gs_br
