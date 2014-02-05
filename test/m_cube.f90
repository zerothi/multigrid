module m_cube

  use t_mg

  implicit none

contains

  subroutine write_cube(name,grid)

    character(len=*), intent(in) :: name
    type(mg_grid), intent(in) :: grid

    integer :: x, y, z, io, i

    real(grid_p), pointer :: V(:,:,:)
    character(len=50) :: fmt
    logical :: exist

    V => grid%V

    io = 10
    open(io,file=trim(name)//'.cube',form='formatted')

    ! two comment lines
    write(io,'(a)') 'Comment'
    write(io,'(a)') 'Comment'

    write(io,'(i0,3(tr1,e10.5))') 1,grid%offset
    do i = 1 , 3
       if ( grid%n(i) > 1 ) then
          write(io,'(i5,3(tr1,e12.6))') grid%n(i),grid%cell(:,i)/grid%n(i)
       else
          ! we create a fictional 3D cell
          write(io,'(i5,3(tr1,e12.6))') 3,grid%cell(:,i)/3
       end if
    end do

    ! write out the "default" atom
    write(io,'(i5,4(tr1,f12.6))') 1,(/0._dp,0._dp,0._dp,0._dp/)

    ! create the format string
    write(fmt,'(a,i0,a)') '(',grid%n(1),'(e12.6,tr1))'
    
    fmt = '(e12.6)'

    ! write cube data...
    if ( grid%n(3) == 1 ) then
       do i = 1 , grid%n(2) * grid%n(3)
          write(io,fmt) 0._grid_p
       end do
    end if
       
    do x = 1 , grid%n(1)
       do y = 1 , grid%n(2)
          do z = 1 , grid%n(3)
             write(io,fmt) V(x,y,z)
          end do
       end do
    end do

    if ( grid%n(3) == 1 ) then
       do i = 1 , grid%n(1) * grid%n(2)
          write(io,fmt) 0._grid_p
       end do
    end if

    close(io)

  end subroutine write_cube
  
end module m_cube
