! Program to create the MULTIGRID method

program mg

  use t_mg

  use m_io
  use m_mg_io
  use m_mg_save

  use m_gs_CDS

  character(len=100) :: filein, fileout
  character(len=300) :: line
  
  integer :: narg

  type(tIO) :: IO
  type(mg_grid) :: grid

  ! solutionmethod
  integer :: method, itmp

  ! Get input file, default to mg.input
  filein = 'mg.input'
  narg = command_argument_count()
  do while ( narg > 0 )
     if (narg == 1) then
        call get_command_argument(narg,filein)
     end if
     narg = narg - 1
  end do

  call io_crt(IO,filein)

  call iomg_read(IO,grid)

  ! Now we have setup the grid, lets compute something
  ! First parse the options of the multigrid-method employed
  call io_open(IO)

  method = CDS_BOTTOM_UP
  line = io_step(IO,'method')
  if ( line(1:1) /= '#' ) then
     line = strip(line)
     if ( has_sub(line,'v-cycle') .or. has_sub(line,'v') .or. &
          has_sub(line,'w') .or. has_sub(line,'w-cycle') ) then
        method = CDS_W_CYCLE
     else if ( has_sub(line,'bottom-up') .or. has_sub(line,'bu') ) then
        method = CDS_BOTTOM_UP
     end if
  end if

  call print_grid(grid)

  ! Whether the user request an initial saving of the grid
  ! Save grid
  line = io_step(IO,'init-save',case='N')
  if ( line(1:1) /= '#' ) then

     ! Create the initial grid
     call grid_bring_back(grid)
     call grid_setup(grid,init=.true.)

     itmp = -100
     do while ( itmp /= IO%il ) 
        if ( itmp == -100 ) itmp = IO%il
        if ( line(1:1) /= '#' ) then
           fileout = strip(line)
           call mg_save(grid,fileout)
        end if
        line = io_step(IO,'init-save',case='N')
     end do
     
     ! Delete the grid information
     call grid_hold_back(grid)

  end if

  call mg_gs_cds(grid,method)

  ! Save grid
  fileout = 'mg.vmg'
  line = io_step(IO,'save',case='N')
  if ( line(1:1) == '#' ) then
     ! We default to saving the vmg
     call mg_save(grid,fileout)
  else
     itmp = -100
     do while ( itmp /= IO%il ) 
        if ( itmp == -100 ) itmp = IO%il
        if ( line(1:1) /= '#' .and. &
             startswith(line,'save') ) then
           fileout = strip(line)
           call mg_save(grid,fileout)
        end if
        line = io_step(IO,'save',case='N')
     end do
  end if

  call io_destroy(IO)

  call delete_grid(grid)
  
end program mg
