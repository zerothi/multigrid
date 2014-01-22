module m_cube

  use t_mg

  implicit none

contains

  subroutine write_cube(name,grid)

    character(len=*), intent(in) :: name
    type(mg_grid), intent(in) :: grid

    integer :: x, zy, io

    character(len=50) :: fmt
    logical :: exist

    io = 10
    open(io,filename=trim(name)//'.cube',form='formatted')

    write(io,'(i0,3(tr1,e10.5))') 0,0._grid_p,0._grid_p,0._grid_p
    write(io,'(i0,3(tr1,e10.5))') grid%n1,1._grid_p,0._grid_p,0._grid_p
    write(io,'(i0,3(tr1,e10.5))') grid%n2,0._grid_p,1._grid_p,0._grid_p
    if ( grid%n3 == 1 ) then
       write(io,'(i0,3(tr1,e10.5))') 3,0._grid_p,0._grid_p,1._grid_p
    else
       write(io,'(i0,3(tr1,e10.5))') grid%n3,0._grid_p,0._grid_p,1._grid_p
    end if

    write(fmt,'(a,i0,a)') '(',grid%n1,'(e10.5))'
    
    ! write cube data...
    if ( grid%n3 == 1 ) then
       do zy = 0 , grid%n3 * grid%n2 - 1
          write(io,fmt) (0._grid_p,i=1,grid%n1)
       end do
    end if
       
    do zy = 0 , grid%n3 * grid%n2 - 1
       write(io,fmt) V(zy+1:zy+grid%n1)
    end do

    if ( grid%n3 == 1 ) then
       do zy = 0 , grid%n3 * grid%n2 - 1
          write(io,fmt) (0._grid_p,i=1,grid%n1)
       end do
    end if

    close(io)

  end subroutine write_cube

end module m_cube
