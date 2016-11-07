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

  implicit none
  private
  public :: mg_save

  integer, parameter, public :: MG_SAVE_CUBE = 1
#ifdef MG__CDF
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

    use m_io, only : lcase

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: method

    integer :: i,lmethod
    character(len=len(filename)+5) :: fname

    ! Base default method on extension
    fname = filename
    i = index(filename,'.')
    if ( i > 0 ) fname = filename(:i-1)
    if ( lcase(filename(i+1:i+4)) == 'cube' ) then
       lmethod = MG_SAVE_CUBE
#ifdef MG__CDF
    else if ( lcase(filename(i+1:i+2)) == 'nc' ) then
       lmethod = MG_SAVE_CDF
#endif
    else if ( lcase(filename(i+1:i+3)) == 'vmg' ) then
       lmethod = MG_SAVE_BINARY
    else if ( lcase(filename(i+1:i+6)) == 'vmgasc' ) then
       lmethod = MG_SAVE_ASCII
    else
       lmethod = MG_SAVE_BINARY
    end if
    if ( present(method) ) lmethod = method

    select case( lmethod )
    case ( MG_SAVE_CUBE )
       call mg_cube(grid,fname)
#ifdef MG__CDF
    case ( MG_SAVE_CDF )
       call mg_cdf(grid,fname)
#endif
    case ( MG_SAVE_BINARY )
       call mg_binary(grid,fname)
    case ( MG_SAVE_ASCII )
       call mg_ascii(grid,fname)
    case default
#ifdef MG__CDF
       stop 'Error in file type, must be CUBE/NC/VMG/VMGASC'
#else
       stop 'Error in file type, must be CUBE/VMG/VMGASC'
#endif
    end select
    
  end subroutine mg_save

  subroutine mg_cube(grid,name)

    use m_io, only : next_unit

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: name

    integer :: x, y, z, io, i

    real(grid_p), pointer :: V(:,:,:)
    character(len=50) :: fmt
    logical :: exist

    V => grid%V

    io = next_unit()
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
    write(io,'(i5,4(tr1,f12.6))') 1, SUM(grid%cell(:,:),DIM=2), 0._dp

    ! create the format string
    write(fmt,'(a,i0,a)') '(',grid%n(1),'(e12.6,tr1))'
    
    fmt = '(5(e12.6,tr1),e12.6)'

    ! write cube data...
    if ( grid%n(3) == 1 ) then
       write(io,fmt) (0._grid_p,i=1,grid%n(2)*grid%n(3))
    end if

    ! this is due to the cube file format
    write(io,fmt) (((V(x,y,z),z=1,grid%n(3)),y=1,grid%n(2)),x=1,grid%n(1))
!!$    do x = 1 , grid%n(1)
!!$       do y = 1 , grid%n(2)
!!$          do z = 1 , grid%n(3)
!!$             write(io,fmt) V(x,y,z)
!!$          end do
!!$       end do
!!$    end do

    if ( grid%n(3) == 1 ) then
       write(io,fmt) (0._grid_p,i=1,grid%n(1)*grid%n(2))
!!$       do i = 1 , grid%n(1) * grid%n(2)
!!$          write(io,fmt) 0._grid_p
!!$       end do
    end if

    close(io)

  end subroutine mg_cube

#ifdef MG__CDF

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
    integer :: i
    real(dp) :: Vminmax(2)

    V => grid%V

    ! NetCDF file name
    fname = trim(name)//'.nc'

    call ncdf_create(ncdf, fname, &
         mode=NF90_64BIT_OFFSEt, &
         overwrite = .true. )

    call ncdf_def_dim(ncdf,'x',grid%n(1))
    call ncdf_def_dim(ncdf,'y',grid%n(2))
    call ncdf_def_dim(ncdf,'z',grid%n(3))
    call ncdf_def_dim(ncdf,'xyz',3)
    call ncdf_def_dim(ncdf,'one',1)


    ! define the offset, cell (we should probably also add the boxes)
    dic = ('unit'.kv.'Bohr') // ('info'.kv.'Offset of the cell')
    call ncdf_def_var(ncdf,'offset',NF90_DOUBLE,(/'xyz'/), atts=dic)
    dic = dic // ('info'.kv.'Cell dimensions')

    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), atts=dic)

    ! Define the variable. 
    dic = dic // ('info'.kv.'Electrostatic potential')
    if ( dp == grid_p ) then
       call ncdf_def_var(ncdf,'V',NF90_DOUBLE,(/'x','y','z'/), atts=dic)
    else
       call ncdf_def_var(ncdf,'V',NF90_FLOAT, (/'x','y','z'/), atts=dic)
    end if
    call delete(dic)

    dic = ('info'.kv.'Maximum and minimum of BC in solution')
    call ncdf_def_var(ncdf,'Vmin',NF90_DOUBLE,(/'one'/),atts=dic)
    call ncdf_def_var(ncdf,'Vmax',NF90_DOUBLE,(/'one'/),atts=dic)
    Vminmax(1) =  huge(1._dp)
    Vminmax(2) = -huge(1._dp)
    do i = 1 , grid%N_box
       Vminmax(1) = min(Vminmax(1),grid%box(i)%val)
       Vminmax(2) = max(Vminmax(2),grid%box(i)%val)
    end do

    dic = ('title'.kv.'Created by Nick R. Papior MG')
    call ncdf_put_gatt(ncdf,atts = dic)

    call ncdf_put_var(ncdf,'offset',grid%offset)
    call ncdf_put_var(ncdf,'cell',grid%cell)
    call ncdf_put_var(ncdf,'V',grid%V(1:grid%n(1),1:grid%n(2),1:grid%n(3)))
    call ncdf_put_var(ncdf,'Vmin',Vminmax(1))
    call ncdf_put_var(ncdf,'Vmax',Vminmax(2))

    call ncdf_close(ncdf)
    
  end subroutine mg_cdf

#endif

  subroutine mg_binary(grid,name)

    use m_io, only : next_unit

    type(mg_grid), intent(in) :: grid
    character(len=*), intent(in) :: name

    integer :: z, io

    real(grid_p), pointer :: V(:,:,:)
    character(len=200) :: fmt
    logical :: exist

    V => grid%V

    io = next_unit()
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
