! Module to print information about the grid that we currently have set up
module m_mg_info

  use t_mg

  implicit none

contains

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
       write(*,trim(fmt)//'a,e10.3)')' -- SOR: ',grid%sor
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

end module m_mg_info
