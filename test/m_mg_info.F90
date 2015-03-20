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
#ifdef CDF
    call mg_save(grid,fname,MG_SAVE_CDF)
#endif
    call mg_save(grid,fname,MG_SAVE_BINARY)
    call mg_save(grid,fname,MG_SAVE_ASCII)

  end subroutine grid_save_all

  subroutine print_grid(top)
    type(mg_grid), intent(in), target :: top

    type(mg_grid), pointer :: grid

    integer :: i, j
    real(grid_p) :: n
    character(len=10) :: fmt

    i = 1
    fmt = '(t1,'
    grid => top
    write(*,*) '*******************************************'
    do while ( associated(grid) )
       write(*,trim(fmt)//'a,i0,a,3(i4,tr1))')' -- Layer: ',grid%layer,' size: ',grid%n
       write(*,trim(fmt)//'a,e10.3)')' -- tolerance: ',grid_tolerance(grid)
       write(*,trim(fmt)//'a,f6.4)')' -- SOR: ',grid%sor
       write(*,trim(fmt)//'a,3(tr1,e10.4))')' -- a(3): ',grid%a
       if ( associated(grid%child) ) then
          write(*,trim(fmt)//'a)',advance='NO')' -- Child: '
          do j = 1 , 3
             write(*,'(i0,tr1)',advance='NO') grid%n(j)-grid%child%n(j)*2
          end do
          write(*,trim(fmt)//'a)',advance='NO')' -- Child-fac: '
          do j = 1 , 3
             n = real(grid%n(j),dp)/real(grid%child%n(j),dp)
             write(*,'(e10.4,tr1)',advance='NO') n - int(n)
          end do
          write(*,*) ''
       end if
       do j = 1 , grid%N_box
          call print_box(grid%box(j),indent = i)
       end do
       i = i + 1
       write(fmt,'(a,i0,a)') '(t',i,','
       grid => grid%child
    end do
    write(*,*) '*******************************************'
  end subroutine print_grid


  subroutine print_box(box,indent)
    type(mg_box), intent(in) :: box
    integer, intent(in), optional :: indent
    integer :: i
    character(len=10) :: fmt

    i = 1
    if ( present(indent) ) then
       write(fmt,'(a,i0,a)') '(t',indent,','
    else
       fmt = '(t1,'
    end if
    
    write(*,trim(fmt)//'a,e10.3,a,l2,a,6(tr1,i0))')' ++ Box, value: ',&
         box%val,' constant?: ',box%constant, &
         'place: ',box%place
  end subroutine print_box

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
