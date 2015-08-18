! Module for reading an input file for defining
! structures etc.

! Read in a txt file and parse to create an MG



! This code has been fully created by 
! Nick Papior Andersen
! nickpapior@gmail.com
! 2015

module m_mg_io

  use t_mg
  use m_io

  implicit none

contains

  ! Create a layer stack by this block
  ! If the word is "layers" it defines a default setting for all layers
  ! If an integer, one can customize individual layers
  ! If a negative integer, one counts from the back of the available 
  ! layers.
  ! Star denotes that the option can only be defined
  ! using the "layers" block
  ! offset <val> <val> <val> 
  !   begin cell
  !    a1 a2 a3
  !    b1 b2 b3
  !    c1 c2 c3
  !   end cell
  ! or the orthogonal cell:
  !   cell a1 b2 c3
  ! begin layers | layer [<int>|-<int>]
  !   tolerance <value> ! tolerance for convergence
  !   sor <value> ! over-relaxation parameter
  !   nnn <int> <int> <int> size of grid in each direction
  !   v-steps <int> ! number of steps done for each pass in a V-cycle
  !   *max-layers <int> ! maximum number of grid partitions
      !define the boundary condition:
  !   bc|boundary-condition [abc][|-+01]|all neumann|dirichlet|periodic
  !   restriction half|full
  !   prolongation half|full
  ! end
  subroutine iomg_read(IO,grid)
    use m_unit

    type(tIO), intent(inout) :: IO
    type(mg_grid), intent(inout) :: grid

    ! IO variables
    character(len=300) :: line
    character(len=100) :: opt
    integer :: lBC, i
    type(mg_grid), pointer :: g

    ! Conversion tool
    real(dp) :: conv 

    ! The cell information
    logical :: found(2), mask(2)
    real(dp) :: cell(3,3), offset(3)
    real(dp) :: sor, tol
    integer :: restrict, prolong, steps
    integer :: nnn(3), BC(2,3)
    integer :: N_layers, i_layer

    ! Box information
    integer :: N_boxes, i_box
    real(dp) :: llc(3), value, rho
    logical :: constant
    integer :: w_pref

    ! First destroy the grid
    call delete_grid(grid)

    ! Open the file
    call io_open(IO)

    found = .false.
    offset = 0._dp
    cell = 0._dp
    nnn = 0
    sor = 1.8_dp
    tol = 1.e-4_dp
    BC = MG_BC_DIRICHLET
    restrict = MG_INTERP_FULL
    prolong = MG_INTERP_FULL
    N_layers = 1000
    ! Simply prepare enough boxes
    N_boxes = 10

    ! Read in the options
    line = io_step(IO,'length')
    if ( line(1:1) /= '#' ) then
       line = strip(line)
       if ( startswith(line,'ang') ) then
          conv = Ang
       else if ( startswith(line,'bohr') ) then
          conv = 1._dp
       else if ( startswith(line,'nm') ) then
          conv = 10._dp * Ang
       else
          stop 'Could not decipher the length unit [Ang|Bohr|nm]'
       end if
    else
       ! We default to the Ang input
       conv = Ang
    end if

    line = io_step(IO,'offset')
    if ( line(1:1) /= '#' ) then
       line = strip(line)
       ! Read in the offset
       read(line,*) offset
       offset = offset * conv
    end if

    line = io_step(IO,'cell')
    if ( line(1:1) /= '#' ) then
       found(1) = .true.
       if ( startswith(line,'begin') ) then
          ! we have a block cell
          do i = 1 , 3
             line = io_line(IO)
             read(line,*) cell(:,i)
          end do
       else
          line = strip(line)
          ! Read in the offset
          read(line,*) cell(1,1),cell(2,2),cell(3,3)
       end if
       cell = cell * conv
    end if

    line = io_step(IO,'max-layers')
    if ( line(1:1) /= '#' ) then
       line = strip(line)
       read(line,*) N_layers
    end if

    line = io_step(IO,'boxes')
    if ( line(1:1) /= '#' ) then
       line = strip(line)
       read(line,*) N_boxes
    end if

    ! Cell-based weight, short cell vectors have high impact
    w_pref = 0
    line = io_step(IO,'weight')
    if ( line(1:1) /= '#' ) then
       line = strip(line)
       if ( has_sub(line,'equal',word = .true.) ) then
          w_pref = 0
       else if ( has_sub(line,'short',word = .true.) ) then
          w_pref = 1
       else if ( has_sub(line,'long',word = .true.) ) then
          w_pref = -1
       end if
    end if

    ! Read in layer information
    line = io_step(IO,'layers')
    if ( line(1:1) /= '#' ) then
       found(2) = .true.
       write(*,'(a)') 'Reading generic layer information'
       call populate_layer_info()

       ! Create layers
       call init_grid(grid,nnn,cell,N_boxes, &
            tol=tol,offset=offset,sor=sor,steps=steps)

       call init_grid_children_half(grid,max_layer=N_layers)
       ! Count number of layers
       N_layers = layers(grid)

       ! Create all BC
       call grid_BC(grid,BC(1,1),MG_BC_A0)
       call grid_BC(grid,BC(2,1),MG_BC_A1)
       call grid_BC(grid,BC(1,2),MG_BC_B0)
       call grid_BC(grid,BC(2,2),MG_BC_B1)
       call grid_BC(grid,BC(1,3),MG_BC_C0)
       call grid_BC(grid,BC(2,3),MG_BC_C1)

       do i_layer = 1 , N_layers
          call grid_set(grid,layer=i_layer,weight = w_pref, &
               restrict = restrict, prolong = prolong ) 
       end do

    else
       
       write(*,'(a)') 'Could not read generic layer information'
       write(*,'(a)') 'This information is necessary!'

       stop

    end if

    ! Default to not set the BC
    BC = -1

    line = io_step(IO,'layer')
    i = IO%il
    mask(1) = .true.
    do while ( mask(1) )
       line = io_step(IO,'layer')
       if ( IO%il == i ) mask(1) = .false.

       if ( line(1:1) /= '#' ) then
          ! first retrieve the index
          line = strip(strip(line))
          read(line,*) i_layer
          if ( i_layer < 0 ) then
             i_layer = layers(grid) + i_layer + 1
          else if ( i_layer == 0 ) then
             stop 'The 0th layer is not existing, i < 0 < i, only!'
          end if
          
          if ( 0 < i_layer .and. i_layer <= layers(grid) ) then

             write(*,'(a,i0,a)') 'Reading layer ',i_layer,' information'

             g => grid_layer(grid,i_layer)
             call populate_layer_info()
             
             call grid_set(g,tol=tol,sor=sor,steps=steps, &
                  restrict=restrict,prolong=prolong)

             ! Create all BC
             call grid_BC(grid,BC(1,1),MG_BC_A0)
             call grid_BC(grid,BC(2,1),MG_BC_A1)
             call grid_BC(grid,BC(1,2),MG_BC_B0)
             call grid_BC(grid,BC(2,2),MG_BC_B1)
             call grid_BC(grid,BC(1,3),MG_BC_C0)
             call grid_BC(grid,BC(2,3),MG_BC_C1)
             
          end if
       end if

    end do

    ! Read in all boxes
    i_box = 0
    line = io_step(IO,'box')
    i = IO%il
    mask(1) = .true.
    do while ( mask(1) )
       line = io_step(IO,'box')
       if ( IO%il == i ) mask(1) = .false.

       if ( line(1:1) /= '#' ) then

          ! Increment box-counter
          i_box = i_box + 1
          write(*,'(a,i0,a)') 'Reading box ',i_box,' information'

          call populate_box_info()

          call grid_add_box(grid, llc, cell, value, rho, constant)
          
       end if
       
    end do

    call io_close(IO)
    
  contains

    subroutine populate_layer_info()
      line = io_line(IO)
      do while ( .not. has_sub(line,'end') )
         
         if ( startswith(line,'sor') ) then
            sor = next_real(line)
         else if ( startswith(line,'tol') .or. &
              startswith(line,'tolerance') ) then
            tol = next_real(line)
         else if ( startswith(line,'nnn') .or. &
              startswith(line,'size') .or. &
              startswith(line,'mesh-size') ) then
            nnn(1) = next_int(line)
            nnn(2) = next_int(line)
            nnn(3) = next_int(line)
         else if ( startswith(line,'v-steps') ) then
            steps = next_int(line)
         else if ( startswith(line,'bc') .or. &
              startswith(line,'boundary-condition') ) then
            
            line = strip(line)
            opt = strip(line)
            ! Figure out the BC
            if ( startswith(opt,'periodic') ) then
               lBC = MG_BC_PERIODIC
            else if ( startswith(opt,'dirichlet') ) then
               lBC = MG_BC_DIRICHLET
            else if ( startswith(opt,'neumann') ) then
               lBC = MG_BC_NEUMANN
            end if
            if ( startswith(line,'all') ) then
               BC(:,:) = lBC
            else
               mask = .false.
               ! If it is positive then we remov
               if ( has_sub(line,'+') ) then
                  ! upper bound
                  mask(2) = .true.
               else if ( has_sub(line,'-') ) then
                  ! lower bound
                  mask(1) = .true.
               else
                  ! both
                  mask = .true.
               end if
               if ( startswith(line,'a') ) then
                  if ( mask(1) ) BC(1,1) = lBC
                  if ( mask(2) ) BC(2,1) = lBC
               else if ( startswith(line,'b') ) then
                  if ( mask(1) ) BC(1,2) = lBC
                  if ( mask(2) ) BC(2,2) = lBC
               else if ( startswith(line,'c') ) then
                  if ( mask(1) ) BC(1,3) = lBC
                  if ( mask(2) ) BC(2,3) = lBC
               end if
            end if
         else if ( startswith(line,'interp') ) then
            line = strip(line)
            if ( startswith(line,'half') ) then
               restrict = MG_INTERP_HALF
               prolong = MG_INTERP_HALF
            end if
            if ( startswith(line,'full') ) then
               restrict = MG_INTERP_FULL
               prolong = MG_INTERP_FULL
            end if
         else if ( startswith(line,'restrict') ) then
            line = strip(line)
            if ( startswith(line,'half') ) restrict = MG_INTERP_HALF
            if ( startswith(line,'full') ) restrict = MG_INTERP_FULL
         else if ( startswith(line,'prolong') ) then
            line = strip(line)
            if ( startswith(line,'half') ) prolong = MG_INTERP_HALF
            if ( startswith(line,'full') ) prolong = MG_INTERP_FULL
         end if

         line = io_line(IO)
         
      end do
    end subroutine populate_layer_info

    subroutine populate_box_info()
      integer :: i
      cell = 0._dp
      llc = 0._dp
      value = 0._dp
      ! The density of a point *must* be larger than or equal to 1
      rho = 1._dp
      constant = .true.
      line = io_line(IO)
      do while ( .not. has_sub(line,'end') ) 

         if ( startswith(line,'llc') .or. &
              startswith(line,'lower-left-corner') ) then
            line = strip(line)
            read(line,*) llc
         else if ( startswith(line,'cell') ) then
            line = strip(line)
            read(line,*) cell(1,1),cell(2,2),cell(3,3)
         else if ( has_sub(line,'cell') ) then
            ! We must have begin end cell
            do i = 1 , 3
               line = io_line(IO)
               read(line,*) cell(:,i)
            end do
            ! Pass the end block
            line = io_line(IO)
         else if ( startswith(line,'value') .or. &
              startswith(line,'val') ) then
            value = next_real(line)
         else if ( startswith(line,'density') .or. &
              startswith(line,'rho') ) then
            rho = next_real(line)
         else if ( startswith(line,'constant') ) then
            line = strip(line)
            ! check for true false etc
            if ( startswith(line,'true') .or. startswith(line,'t') ) then
               constant = .true.
            else if ( startswith(line,'false') .or. startswith(line,'f') ) then
               constant = .false.
            else
               write(*,'(a)') 'Warning: Read constant value, could not be true/t/false/f'
            end if
         end if

         line = io_line(IO)
      end do


      llc = llc * conv
      cell = cell * conv

    end subroutine populate_box_info
    
  end subroutine iomg_read

  function next_real(l) result(r)
    character(len=*), intent(inout) :: l
    real(dp) :: r
    l = strip(l)
    read(l,*) r
  end function next_real
  
  function next_int(l) result(i)
    character(len=*), intent(inout) :: l
    integer :: i
    l = strip(l)
    read(l,*) i
  end function next_int

end module m_mg_io
