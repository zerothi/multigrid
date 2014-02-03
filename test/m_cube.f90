module m_cube

  use t_mg

  implicit none

contains

  subroutine write_cube(name,grid)

    character(len=*), intent(in) :: name
    type(mg_grid), intent(in) :: grid

    integer :: x, zy, io, i

    character(len=50) :: fmt
    logical :: exist

    io = 10
    open(io,file=trim(name)//'.cube',form='formatted')

    write(io,'(i0,3(tr1,e10.5))') 0,grid%offset
    do i = 1 , 3
       if ( i > 1 ) then
          write(io,'(i0,3(tr1,e10.5))') grid%n(i),grid%cell(:,i)
       else
          ! we create a fictional 3D cell
          write(io,'(i0,3(tr1,e10.5))') 3,grid%cell(:,i)
       end if
    end do

    ! create the format string
    write(fmt,'(a,i0,a)') '(',grid%n(1),'(e10.5))'
    
    ! write cube data...
    if ( grid%n(3) == 1 ) then
       do zy = 0 , grid%n(3) * grid%n(2) - 1
          write(io,fmt) (0._grid_p,i=1,grid%n(1))
       end do
    end if
       
    do zy = 0 , grid%n(3) * grid%n(2) - 1
       write(io,fmt) grid%V(zy+1:zy+grid%n(1))
    end do

    if ( grid%n(3) == 1 ) then
       do zy = 0 , grid%n(3) * grid%n(2) - 1
          write(io,fmt) (0._grid_p,i=1,grid%n(1))
       end do
    end if

    close(io)

  end subroutine write_cube
  
end module m_cube
