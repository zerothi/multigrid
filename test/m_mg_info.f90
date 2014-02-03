! Module to print information about the grid that we currently have set up
module m_mg_info

  use t_mg

  implicit none

contains

  subroutine print_grid(top)
    type(mg_grid), intent(in), target :: top

    type(mg_grid), pointer :: grid

    integer :: i, j
    character(len=10) :: fmt

    i = 1
    fmt = '(t1,'
    grid => top
    write(*,*) '*******************************************'
    do while ( associated(grid) )
       write(*,trim(fmt)//'a,i0,a,3(i4,tr1))')' -- Layer: ',grid%layer,' size: ',grid%n
       write(*,trim(fmt)//'a,e10.3)')' -- tolerance: ',grid%tol
       write(*,trim(fmt)//'a,e10.3)')' -- SOR: ',grid%sor
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
    
    write(*,trim(fmt)//'a,e10.3,a,l2)')' ++ Box, value: ',box%val,' constant?: ',box%constant
  end subroutine print_box

end module m_mg_info
