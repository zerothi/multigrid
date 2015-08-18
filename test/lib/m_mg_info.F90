! Module to print information about the grid that we currently have set up
module m_mg_info

  use t_mg

  implicit none

contains

  subroutine grid_save_all(grid,fname)
    use m_mg_save
    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: fname

    ! Save all available formats
    call mg_save(grid,fname,MG_SAVE_CUBE)
#ifdef MG__CDF
    call mg_save(grid,fname,MG_SAVE_CDF)
#endif
    call mg_save(grid,fname,MG_SAVE_BINARY)
    call mg_save(grid,fname,MG_SAVE_ASCII)

  end subroutine grid_save_all

  ! Create a script for python to plot the grid in matplotlib
  subroutine mtpl_grid(top,out)
    type(mg_grid), intent(in), target :: top
    integer, intent(in), optional :: out
    
    type(mg_grid), pointer :: grid
    
    integer :: lu
    
    lu = 6
    if ( present(out) ) lu = out
    
    write(lu,'(a)') "#!/usr/bin/env python"
    write(lu,'(a)') 
    write(lu,'(a)') "import matplotlib.pyplot as plt"
    write(lu,'(a)') "from mpl_toolkits.mplot3d import Axes3D"
    grid => top
    
    do while ( associated(grid) )
       
       write(lu,'(a)') "fig = plt.figure()"
       write(lu,'(a)') "ax = fig.add_subplot(111, projection='3d')"
       
       ! Create full box

       grid => grid%child
    end do
    
  end subroutine mtpl_grid

end module m_mg_info
