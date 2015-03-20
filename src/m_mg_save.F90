! Module for saving the grid in a few kinds.


! This code has been fully created by 
! Nick Papior Andersen
! nickpapior@gmail.com
! 2014

! The intent for this software is to ease the saving of the
! multi-grid output and save it in various formats.

module m_mg_save

  ! get access to the multi-grid type.
  use t_mg

  private
  public :: mg_save

  integer, parameter, public :: MG_SAVE_CUBE = 1
#ifdef CDF
  integer, parameter, public :: MG_SAVE_CDF = 2
#endif
  integer, parameter, public :: MG_SAVE_BINARY = 3
  integer, parameter, public :: MG_SAVE_ASCII = 4

contains

  ! Routine for saving a grid into file-name.
  !   grid [type(mg_grid)]   grid that we wish to save
  !   filename [character]   name of file that we wish to save in
  !   method [integer]       the method that we use for saving the file
  ! (determines the extension)
  subroutine mg_save(grid,filename,method)

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: method

    integer :: lmethod

    lmethod = MG_SAVE_BINARY
    if ( present(method) ) lmethod = method

    select case( lmethod )
    case ( MG_SAVE_CUBE )
       call mg_cube(grid,filename)
#ifdef CDF
    case ( MG_SAVE_CDF )
       call mg_cdf(grid,filename)
#endif
    case ( MG_SAVE_BINARY )
       call mg_binary(grid,filename)
    case ( MG_SAVE_ASCII )
       call mg_ascii(grid,filename)
    case default
#ifdef CDF
       stop 'Error in file type, must be CUBE/CDF/BIN/ASCII'
#else
       stop 'Error in file type, must be CUBE/BIN/ASCII'
#endif
    end select
    
  end subroutine mg_save

  subroutine mg_cube(grid,name)

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: name

    integer :: x, y, z, io, i

    real(grid_p), pointer :: V(:,:,:)
    character(len=50) :: fmt
    logical :: exist

    V => grid%V

    io = 10
    open(io,file=trim(name)//'.cube',form='formatted')

    ! two comment lines
    write(io,'(a)') 'Created by Nick Papior Andersens multigrid-method'
    write(io,'(a)') 'Line not read'

    ! write out the offset of the grid
    write(io,'(i0,3(tr1,e10.5))') 2,grid%offset
    do i = 1 , 3
       if ( grid%n(i) > 1 ) then
          write(io,'(i5,3(tr1,e12.6))') grid%n(i),grid%cell(:,i)/grid%n(i)
       else
          ! we create a fictional 3D cell
          write(io,'(i5,3(tr1,e12.6))') 3,grid%cell(:,i)/3
       end if
    end do

    ! write out the "default" atom (no atoms)
    write(io,'(i5,4(tr1,f12.6))') 1,(/0._dp,0._dp,0._dp,0._dp/)
    write(io,'(i5,4(tr1,f12.6))') 1, grid%cell(:,1)+ &
         grid%cell(:,2) + grid%cell(:,3)

    ! create the format string
    write(fmt,'(a,i0,a)') '(',grid%n(1),'(e12.6,tr1))'
    
    fmt = '(e12.6)'

    ! write cube data...
    if ( grid%n(3) == 1 ) then
       do i = 1 , grid%n(2) * grid%n(3)
          write(io,fmt) 0._grid_p
       end do
    end if

    ! this is due to the cube file format
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

  end subroutine mg_cube

#ifdef CDF

  subroutine mg_cdf(grid,name)

    use variable
    use dictionary

    use nf_ncdf

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: name

    type(hNCDF) :: ncdf
    type(dict) :: dic

    real(grid_p), pointer :: V(:,:,:)
    character(len=200) :: fname

    V => grid%V

    ! NetCDF file name
    fname = trim(name)//'.nc'

    call ncdf_create(ncdf, fname, &
         mode=ior(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
         overwrite = .true. )

    call ncdf_def_dim(ncdf,'x',grid%n(1))
    call ncdf_def_dim(ncdf,'y',grid%n(2))
    call ncdf_def_dim(ncdf,'z',grid%n(3))
    call ncdf_def_dim(ncdf,'xyz',3)

    ! define the offset, cell (we should probably also add the boxes)
    dic = ('unit'.kv.'Bohr') // ('info'.kv.'Offset of the cell')
    call ncdf_def_var(ncdf,'offset',NF90_DOUBLE,(/'xyz'/), atts=dic)
    dic = dic // ('info'.kv.'Cell dimensions')

    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), atts=dic)

    ! Define the variable. 
    dic = dic//('unit'.kv.'Ry') // ('info'.kv.'Potential') // &
         ('ATT_DELETE'.kv.1)

    if ( dp == grid_p ) then
       call ncdf_def_var(ncdf,'V',NF90_DOUBLE,(/'x','y','z'/), atts=dic)
    else
       call ncdf_def_var(ncdf,'V',NF90_FLOAT,(/'x','y','z'/), atts=dic)
    end if

    dic = ('title'.kv.'Created by Nick Papior Andersen MG')
    call ncdf_put_gatt(ncdf,atts = dic)

    call ncdf_put_var(ncdf,'offset',grid%offset)
    call ncdf_put_var(ncdf,'cell',grid%cell)
    call ncdf_put_var(ncdf,'V',grid%V)

    call ncdf_close(ncdf)
    
  end subroutine mg_cdf

#endif

  subroutine mg_binary(grid,name)

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: name

    integer :: z, io

    real(grid_p), pointer :: V(:,:,:)
    character(len=200) :: fmt
    logical :: exist

    V => grid%V

    io = 10
    open(io,file=trim(name)//'.VMG',form='unformatted')

    ! comment line
    fmt = 'Created by Nick Papior Andersens multigrid-method'
    write(io) fmt

    ! cell number of partitions
    write(io) grid%n
    ! write out the offset of the grid
    write(io) grid%offset
    ! cell size
    write(io) grid%cell
    ! write each x-y plane in their own record (allows finer partitioning 

    do z = 1 , grid%n(3)
       write(io) V(:,:,z)
    end do

    close(io)

  end subroutine mg_binary

  subroutine mg_ascii(grid,name)

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: name

    integer :: io, x, y, z, i

    real(grid_p), pointer :: V(:,:,:)
    character(len=50) :: fmt
    logical :: exist

    V => grid%V

    io = 10
    open(io,file=trim(name)//'.VMGASC',form='formatted')

    ! two comment lines
    write(io,'(a)') 'Created by Nick Papior Andersens multigrid-method'

    ! cell number of partitions
    write(io,'(2(i10,tr1),i10)') grid%n
    ! write out the offset of the grid
    write(io,'(2(e15.10,tr1),e15.10)') grid%offset
    ! cell size
    write(io,'(3(2(e15.10,tr1),e15.10,/))') grid%cell

    ! write each x-y plane in their own record (allows finer partitioning 
    do z = 1 , grid%n(3)
    do y = 1 , grid%n(2)
    do x = 1 , grid%n(1)
       write(io,'(e15.10)') V(x,y,z)
    end do
    end do
    end do

    close(io)
    
  end subroutine mg_ascii

end module m_mg_save
