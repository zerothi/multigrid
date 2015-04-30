module t_mg

  use t_BC

  implicit none

  integer, parameter :: dp = selected_real_kind(p=15)
  integer, parameter :: grid_p = selected_real_kind(p=6)

  ! a module to sustain a "simple" multi-grid solver
  integer, parameter :: MG_CELL_A0 = 1
  integer, parameter :: MG_CELL_A1 = 2
  integer, parameter :: MG_CELL_B0 = 4
  integer, parameter :: MG_CELL_B1 = 8
  integer, parameter :: MG_CELL_C0 = 16
  integer, parameter :: MG_CELL_C1 = 32

  ! Current interpolation methods
  integer, parameter :: MG_INTERP_FULL = 1 
  integer, parameter :: MG_INTERP_HALF = 2

  type :: mg_grid
     ! the grid information
     logical      :: enabled =.true. ! turn on/off the run on this grid
     real(dp)     :: offset(3) ! the offset of the placement of the local grid 
     real(dp)     :: cell(3,3) ! the cell for the local grid
     real(dp)     :: dL(3,3) ! the cell stepping
     real(dp)     :: dLL(3) ! the voxel side lengths
     real(dp)     :: dVol ! volume of voxel
     real(dp)     :: Vol  ! volume of cell
     integer      :: n(3) ! size in each direction
     real(dp)     :: sor  ! the SOR value
     real(grid_p) :: a(3) ! the pre-factors for the summation
     real(dp)     :: tol  ! the tolerance of the current grid
                          ! this allows different tolerances for different layer-grids
     integer :: itt       ! iterations currently processed
     integer :: steps     ! number of steps in each V-cycle
     integer :: layer = 0 ! the layer that this grid resides in
     real(grid_p),  pointer :: V  (:,:,:) => null() ! update array
     real(grid_p),  pointer :: g  (:) => null() ! ghost arrays ( one long array for all bounds )
     real(grid_p),  pointer :: g_s(:) => null() ! send ghost arrays ( one long array for all bounds )
     type(mg_grid), pointer :: parent => null()
     type(mg_grid), pointer :: child => null()

     ! The constant valued boxes in this grid
     integer :: N_box = 0
     type(mg_box), pointer :: box(:) => null()

     ! The max change in the grid
     real(grid_p) :: err

     ! The BC of the cell
     type(tBC) :: BC(2,3)

     ! the prolongation method
     integer :: PRO_method = MG_INTERP_FULL
     ! the restriction method
     integer :: RES_method = MG_INTERP_FULL

  end type mg_grid

  type :: mg_box
     sequence
     ! this contains the min/max indices for the constant region
     ! xmin/xmax
     ! ymin/ymax
     ! zmin/zmax
     integer :: place(2,3) = 0
     real(grid_p) :: val = 0._grid_p ! 
     real(grid_p) :: rho = 1._grid_p ! *MUST be above 1*
     logical :: constant = .false.
  end type mg_box

contains

  subroutine init_grid(grid, n, cell, boxes, tol, offset, sor, steps)
    type(mg_grid), intent(inout) :: grid
    integer,       intent(in)    :: n(3)
    real(dp),      intent(in)    :: cell(3,3)
    integer,       intent(in)    :: boxes
    real(dp),      intent(in), optional :: tol
    real(dp),      intent(in), optional :: offset(3)
    real(dp),      intent(in), optional :: sor
    integer,       intent(in), optional :: steps

    real(dp) :: celll(3), tmp
    integer  :: i

    ! ensure it is empty
    call delete_grid(grid)

    grid%steps = 2
    if ( present(steps) ) grid%steps = steps

    ! create the ax, ay, az pre-factors for the 3D Poisson solver
    ! note that cell is the cell-size for the total cell
    grid%n = n
    do i = 1 , 3

       celll(i) = &
            cell(1,i) ** 2 + &
            cell(2,i) ** 2 + &
            cell(3,i) ** 2
       celll(i) = celll(i) ! / grid%n(i)
       grid%dL(:,i) = cell(:,i) / grid%n(i)
    end do
    do i = 1 , 3
       grid%dLL(i) = sum(grid%dL(i,:))
    end do
    
    grid%Vol = ( cell(2,1)*cell(3,2) - cell(3,1)*cell(2,2) ) * cell(1,3) + &
         ( cell(3,1)*cell(1,2) - cell(1,1)*cell(3,2) ) * cell(2,3) + &
         ( cell(1,1)*cell(2,2) - cell(2,1)*cell(1,2) ) * cell(3,3)

    grid%dVol = ( grid%dL(2,1)*grid%dL(3,2) - grid%dL(3,1)*grid%dL(2,2) ) * grid%dL(1,3) + &
         ( grid%dL(3,1)*grid%dL(1,2) - grid%dL(1,1)*grid%dL(3,2) ) * grid%dL(2,3) + &
         ( grid%dL(1,1)*grid%dL(2,2) - grid%dL(2,1)*grid%dL(1,2) ) * grid%dL(3,3)

    grid%enabled = .true.
    ! set the values for the grid...
    if ( associated(grid%parent) ) then
       grid%layer = grid%parent%layer + 1
    else
       grid%layer = 1
    end if
    grid%offset = 0._grid_p
    if ( present(offset) ) grid%offset = offset

    ! re-construct the actual cell size for this processor
    do i = 1 , 3
       grid%cell(:,i) = grid%dL(:,i) * grid%n(i)
    end do
       
    grid%itt = 0
    if ( present(tol) ) then
       grid%tol = tol
    else if ( associated(grid%parent) ) then
       grid%tol = grid%parent%tol
    else
       stop 'Tolerance not set'
    end if

    ! the SOR parameter
    grid%sor = 2._grid_p / (1._grid_p + 3.1415926535897_grid_p / maxval(grid%n) )
    if ( present(sor) ) grid%sor = sor

    ! Equal weights for each direction
    grid%a(:) = 1._dp / 3._dp

    ! pre-allocate room for the boxes
    grid%N_box = boxes
    if ( boxes > 0 ) then
       allocate(grid%box(boxes))
    end if

  end subroutine init_grid

  subroutine init_grid_children_half(top,max_layer)
    type(mg_grid), intent(inout), target :: top
    integer, intent(in), optional :: max_layer
    integer :: n(3), lmax_layer
    type(mg_grid), pointer :: grid, tmp_grid
    
    ! we will now create all the children according to the
    ! standard algorithms
    lmax_layer = top%layer - 1
    if ( present(max_layer) ) lmax_layer = max_layer

    call new_grid_size(top,n)

    ! create all the grids
    grid => top
    do while ( all(n /= 0) .and. lmax_layer /= grid%layer)
       nullify(tmp_grid)
       allocate(tmp_grid)
       grid%child => tmp_grid
       tmp_grid%parent => grid
       call init_grid(tmp_grid,n, &
            grid%cell, grid%N_box, &
            tol=grid%tol, offset=grid%offset, sor=grid%sor)
       grid => tmp_grid
       call new_grid_size(grid,n)
    end do
    nullify(grid,tmp_grid)

  contains

    subroutine new_grid_size(grid,n)
      type(mg_grid), intent(in) :: grid
      integer, intent(out) :: n(3)
      integer :: i
      
      n(:) = grid%n(:) / 2

      do i = 1 , 3
         if ( n(i) < 20 ) then
            if ( any(n > n(i) * 3) ) then
               n(i) = grid%n(i)
            else
               n(:) = 0
               return
            end if
         end if
      end do

    end subroutine new_grid_size

  end subroutine init_grid_children_half

  subroutine grid_set(grid,sor,tol,layer,weight,offset,steps,restrict,prolong)
    type(mg_grid), intent(inout), target :: grid
    real(dp), intent(in), optional :: sor, tol, offset(3)
    integer, intent(in), optional :: layer, steps, restrict, prolong, weight

    type(mg_grid), pointer :: tmp_grid
    integer :: i
    real(dp) :: celll(3)

    tmp_grid => grid
    if ( present(layer) ) then
       do while ( layer /= tmp_grid%layer ) 
          tmp_grid => tmp_grid%child
          if ( .not. associated(tmp_grid) ) return
       end do
    end if
    if ( present(sor) ) &
         tmp_grid%sor = sor
    if ( present(tol) ) &
         tmp_grid%tol = tol
    if ( present(steps) ) &
         tmp_grid%steps = steps
    if ( present(offset) ) &
         tmp_grid%offset = offset
    if ( present(restrict) ) &
         tmp_grid%RES_method = restrict
    if ( present(prolong) ) &
         tmp_grid%PRO_method = prolong
    if ( present(weight) ) then

       do i = 1 , 3

          celll(i) = &
               tmp_grid%dL(1,i) ** 2 + &
               tmp_grid%dL(2,i) ** 2 + &
               tmp_grid%dL(3,i) ** 2
       end do

       select case ( weight ) 
       case ( 0 )
          ! Equal weights
          tmp_grid%a = 1._dp
       case ( 1 ) 

          ! This additional weighting of each cell
          ! means that different situations might occur
          ! 1. If the "important" direction is along a 
          !    longer voxel direction, the convergence is slower
          ! 2. If the "important" direction is along a 
          !    shorter voxel direction, the convergence is faster
          tmp_grid%a(1) = celll(2) * celll(3)
          tmp_grid%a(2) = celll(1) * celll(3)
          tmp_grid%a(3) = celll(1) * celll(2)

          ! Prefer the short cells (default)
       case ( -1 )

          ! This additional weighting of each cell
          ! means that different situations might occur
          ! 1. If the "important" direction is along a 
          !    short voxel direction, the convergence is slower
          ! 2. If the "important" direction is along a 
          !    longer voxel direction, the convergence is faster
          tmp_grid%a(1) = celll(1) / (celll(2) * celll(3))
          tmp_grid%a(2) = celll(2) / (celll(1) * celll(3))
          tmp_grid%a(3) = celll(3) / (celll(1) * celll(2))
       end select
 
       ! this normalization does not really matter, however,
       ! for the convenience of seeing their appropriate weigths
       ! it could be nice to see...
       tmp_grid%a = tmp_grid%a / sum(tmp_grid%a)

    end if

  end subroutine grid_set

  recursive subroutine grid_BC(grid, BC, plane)
    type(mg_grid), intent(inout) :: grid
    ! Define which plane we want to add the dirichlet boundary
    ! conditions on.
    integer, intent(in) :: BC
    integer, intent(in), optional :: plane

    integer :: lplane

    select case ( BC ) 
    case ( MG_BC_PERIODIC , MG_BC_DIRICHLET , MG_BC_NEUMANN )
       ! correct
    case default
       return
    end select

    lplane = IOR(MG_BC_A0,MG_BC_A1)
    lplane = IOR(lplane,MG_BC_B0)
    lplane = IOR(lplane,MG_BC_B1)
    lplane = IOR(lplane,MG_BC_C0)
    lplane = IOR(lplane,MG_BC_C1)
    if ( present(plane) ) lplane = plane

    ! Create all plane-boxes
    ! Start by creating the boxes starting from origo
    if ( iand(MG_CELL_A0,lplane) == MG_CELL_A0 ) &
         grid%BC(1,1)%method = BC
    if ( iand(MG_CELL_A1,lplane) == MG_CELL_A1 ) &
         grid%BC(2,1)%method = BC
    if ( iand(MG_CELL_B0,lplane) == MG_CELL_B0 ) &
         grid%BC(1,2)%method = BC
    if ( iand(MG_CELL_B1,lplane) == MG_CELL_B1 ) &
         grid%BC(2,2)%method = BC
    if ( iand(MG_CELL_C0,lplane) == MG_CELL_C0 ) &
         grid%BC(1,3)%method = BC
    if ( iand(MG_CELL_C1,lplane) == MG_CELL_C1 ) &
         grid%BC(2,3)%method = BC
    if ( associated(grid%child) ) then
       call grid_BC(grid%child,BC,plane=plane)
    end if
    
  end subroutine grid_BC

  recursive subroutine grid_add_box(grid, llc, box_cell, val, rho, constant, recurse)
    type(mg_grid), intent(inout) :: grid
    real(dp), intent(in) :: llc(3), box_cell(3,3)
    real(dp), intent(in) :: val, rho
    logical, intent(in) :: constant
    logical, intent(in), optional :: recurse
    
    real(dp) :: dz(3), dyz(3), xyz(3), urc(3), offset(3)
    integer :: i,x,y,z
    type(mg_box), pointer :: box
    
    if ( grid%N_box == 0 ) stop 'No boxes'

    nullify(box)
    do i = 1 , grid%N_box
       if ( .not. all(grid%box(i)%place == 0) ) cycle
       box => grid%box(i)
       exit
    end do

    if ( .not. associated(box) ) then
       stop 'No boxes available'
    end if

    ! ensure it is empty
    call delete_box(box)

    box%val = val 
    box%rho = rho 
    if ( rho < 0._grid_p ) stop 'not available'
    box%constant = constant

    ! initialize
    box%place(1,:) = huge(0)
    box%place(2,:) = 0

    offset = grid%offset + (grid%dL(:,1)+grid%dL(:,2)+grid%dL(:,3))*.5_dp
    urc = llc + sum(box_cell(:,:),DIM=2)
    ! TODO currently this does not work with skewed axis

!$OMP parallel do default(shared), private(x,y,z,xyz,dyz,dz)
    do z = 0 , grid%n(3) - 1
    dz  = grid%dL(:,3) * z + offset ! we immediately add the offset
    do y = 0 , grid%n(2) - 1
    dyz = grid%dL(:,2) * y + dz
    do x = 0 , grid%n(1) - 1
       xyz = grid%dL(:,1) * x + dyz
       if ( all(llc <= xyz) ) then
          ! check that the point also
          ! lies inside on the right borders
          if ( all(xyz <= urc) ) then
!$OMP critical
             call insert_point(box%place,x,y,z)
!$OMP end critical
          end if
       else if ( all(urc <= xyz) ) then
          ! check that the point also
          ! lies inside on the right borders
          if ( all(xyz <= llc) ) then
!$OMP critical
             call insert_point(box%place,x,y,z)
!$OMP end critical
          end if
       end if
    end do
    end do
    end do
!$OMP end parallel do
    
    ! ensure at least one point!
    do x = 1, 3
       box%place(2,x) = max(box%place(2,x),box%place(1,x))
    end do

    ! if not present, simply delete it...
    if ( all(box%place(1,:) == huge(0)) .and. &
         all(box%place(2,:) == 0) ) then
       call delete_box(box)
    end if

    ! add the box to the child grid
    if ( associated(grid%child) ) then
       x = 0
       if ( present(recurse) ) then
          if ( recurse ) x = 1
       else
          x = 1
       end if
       if ( x == 1 ) then
          call grid_add_box(grid%child,llc,box_cell,val,rho,constant,recurse = recurse)
       end if
    end if

  contains

    subroutine insert_point(place,x,y,z)
      integer, intent(inout) :: place(2,3)
      integer, intent(in) :: x,y,z
      if ( x < place(1,1) ) &
           place(1,1) = x + 1
      if ( place(2,1) <= x ) &
           place(2,1) = x + 1
      if ( y < place(1,2) ) &
           place(1,2) = y + 1
      if ( place(2,2) <= y ) &
           place(2,2) = y + 1
      if ( z < place(1,3) ) &
           place(1,3) = z + 1
      if ( place(2,3) <= z ) &
           place(2,3) = z + 1
    end subroutine insert_point

  end subroutine grid_add_box

  subroutine grid_add_point(grid, llc, val, rho, constant)
    type(mg_grid), intent(inout) :: grid
    real(dp), intent(in) :: llc(3)
    real(dp), intent(in) :: val, rho
    logical, intent(in) :: constant
    
    real(dp) :: bl(3,3)
    
    bl(:,:) = 0._dp
    call grid_add_box(grid,llc,bl,val,rho,constant)

  end subroutine grid_add_point

  subroutine grid_add_line(grid, llc, dir, l, val, rho, constant)
    type(mg_grid), intent(inout) :: grid
    real(dp), intent(in) :: llc(3), l ! the length of the line
    integer, intent(in) :: dir
    real(dp), intent(in) :: val, rho
    logical, intent(in) :: constant
    
    real(dp) :: bl(3,3)
    
    bl = 0._dp
    bl(dir,dir) = l
    call grid_add_box(grid,llc,bl,val,rho,constant)

  end subroutine grid_add_line

  subroutine grid_setup(grid,init)
    type(mg_grid), intent(inout) :: grid
    logical, intent(in), optional :: init
    integer :: x,y,z
    real(grid_p), pointer :: V(:,:,:)

    V => grid%V
    
    if ( present(init) ) then
       if ( init ) then
!$OMP parallel workshare default(shared)
          V = 0._grid_p
!$OMP end parallel workshare
       end if
    end if

!$OMP parallel default(shared)

    ! set all boxes to their values if constant
!$OMP do private(x,y,z), collapse(3)
    do z = 1 , grid%n(3)
    do y = 1 , grid%n(2)
    do x = 1 , grid%n(1)
       if ( is_constant(grid,x,y,z) ) then
          V(x,y,z) = val_constant(grid,x,y,z)
       end if
    end do
    end do
    end do
!$OMP end do

    select case ( grid%BC(1,1)%method ) 
    case ( MG_BC_PERIODIC )
!$OMP workshare
       V(0,:,:) = V(grid%n(1),:,:)
!$OMP end workshare nowait
    case ( MG_BC_DIRICHLET )
!$OMP workshare
       V(0,:,:) = 0._grid_p
!$OMP end workshare nowait
    case ( MG_BC_NEUMANN ) 
!$OMP workshare
       V(0,:,:) = V(1,:,:)
!$OMP end workshare nowait
    end select

    select case ( grid%BC(2,1)%method ) 
    case ( MG_BC_PERIODIC )
!$OMP workshare
       V(grid%n(1)+1,:,:) = V(1,:,:)
!$OMP end workshare nowait
    case ( MG_BC_DIRICHLET )
!$OMP workshare
       V(grid%n(1)+1,:,:) = 0._grid_p
!$OMP end workshare nowait
    case ( MG_BC_NEUMANN ) 
!$OMP workshare
       V(grid%n(1)+1,:,:) = V(grid%n(1),:,:)
!$OMP end workshare nowait
    end select

    select case ( grid%BC(1,2)%method ) 
    case ( MG_BC_PERIODIC )
!$OMP workshare
       V(:,0,:) = V(:,grid%n(2),:)
!$OMP end workshare nowait
    case ( MG_BC_DIRICHLET )
!$OMP workshare
       V(:,0,:) = 0._grid_p
!$OMP end workshare nowait
    case ( MG_BC_NEUMANN ) 
!$OMP workshare
       V(:,0,:) = V(:,1,:)
!$OMP end workshare nowait
    end select

    select case ( grid%BC(2,2)%method ) 
    case ( MG_BC_PERIODIC )
!$OMP workshare
       V(:,grid%n(2)+1,:) = V(:,1,:)
!$OMP end workshare nowait
    case ( MG_BC_DIRICHLET )
!$OMP workshare
       V(:,grid%n(2)+1,:) = 0._grid_p
!$OMP end workshare nowait
    case ( MG_BC_NEUMANN ) 
!$OMP workshare
       V(:,grid%n(2)+1,:) = V(:,grid%n(2),:)
!$OMP end workshare nowait
    end select

    select case ( grid%BC(1,3)%method ) 
    case ( MG_BC_PERIODIC )
!$OMP workshare
       V(:,:,0) = V(:,:,grid%n(3))
!$OMP end workshare nowait
    case ( MG_BC_DIRICHLET )
!$OMP workshare
       V(:,:,0) = 0._grid_p
!$OMP end workshare nowait
    case ( MG_BC_NEUMANN ) 
!$OMP workshare
       V(:,:,0) = V(:,:,1)
!$OMP end workshare nowait
    end select

    select case ( grid%BC(2,3)%method ) 
    case ( MG_BC_PERIODIC )
!$OMP workshare
       V(:,:,grid%n(3)+1) = V(:,:,1)
!$OMP end workshare nowait
    case ( MG_BC_DIRICHLET )
!$OMP workshare
       V(:,:,grid%n(3)+1) = 0._grid_p
!$OMP end workshare nowait
    case ( MG_BC_NEUMANN ) 
!$OMP workshare
       V(:,:,grid%n(3)+1) = V(:,:,grid%n(3))
!$OMP end workshare nowait
    end select

!$OMP end parallel

  end subroutine grid_setup

  subroutine grid_hold_back(grid)
    type(mg_grid), intent(inout) :: grid
    
    if ( associated(grid%V)   ) deallocate(grid%V)
    if ( associated(grid%g)   ) deallocate(grid%g)
    if ( associated(grid%g_s) ) deallocate(grid%g_s)

    nullify(grid%V,grid%g,grid%g_s) ! ensure nullification

  end subroutine grid_hold_back

  subroutine grid_bring_back(grid)
    type(mg_grid), intent(inout) :: grid

    ! ensure that it is empty
    call grid_hold_back(grid)

    ! Re-allocate
    allocate(grid%V(0:grid%n(1)+1,0:grid%n(2)+1,0:grid%n(3)+1))
    allocate(grid%g(grid%n(1)*2+grid%n(2)*2+grid%n(3)*2))
    allocate(grid%g_s(grid%n(1)*2+grid%n(2)*2+grid%n(3)*2))

  end subroutine grid_bring_back

  recursive subroutine delete_grid(grid)
    type(mg_grid), intent(inout) :: grid

    if ( associated(grid%child) ) then
       call delete_grid(grid%child)
       deallocate(grid%child)
       nullify(grid%child)
    end if

    if ( grid%N_box > 0 ) then
       deallocate(grid%box)
       nullify(grid%box)
       grid%N_box = 0
    end if
    call grid_hold_back(grid)

    grid%layer = 0

  end subroutine delete_grid

  pure function is_constant(grid,x,y,z) result(is)
    type(mg_grid), intent(in) :: grid
    integer, intent(in) :: x,y,z
    logical :: is
    integer :: i
    
    do i = 1 , grid%N_box
       if ( .not. grid%box(i)%constant ) cycle
       if ( in_box(grid%box(i),x,y,z) ) then
          is = .true.
          return
       end if
    end do
    is = .false.

  end function is_constant

  pure function val_constant(grid,x,y,z) result(val)
    type(mg_grid), intent(in) :: grid
    integer, intent(in) :: x,y,z
    real(grid_p) :: val

    integer :: i

    do i = 1 , grid%N_box
       if ( in_box(grid%box(i),x,y,z) )then
          val = grid%box(i)%val
          return
       end if
    end do

    ! all values are defaulted to 0
    val = 0._grid_p

  end function val_constant
  
  pure function val_rho(grid,x,y,z) result(rho)
    type(mg_grid), intent(in) :: grid
    integer, intent(in) :: x,y,z
    real(grid_p) :: rho

    integer :: i

    do i = 1 , grid%N_box
       if ( in_box(grid%box(i),x,y,z) )then
          rho = grid%box(i)%rho
          return
       end if
    end do

    ! the "importance" of a grid is defaulted to 1
    rho = 1._grid_p

  end function val_rho

  ! box-routines
  subroutine delete_box(box)
    type(mg_box), intent(inout) :: box
    box%place = 0
    box%val = 0._grid_p
    box%rho = 1._grid_p
    box%constant = .false.
  end subroutine delete_box

  pure function in_box(box,x,y,z) result(in)
    type(mg_box), intent(in) :: box
    integer, intent(in) :: x,y,z
    logical :: in

    in = box%place(1,1) <= x .and. &
         x <= box%place(2,1) .and. &
         box%place(1,2) <= y .and. &
         y <= box%place(2,2) .and. &
         box%place(1,3) <= z .and. &
         z <= box%place(2,3)

  end function in_box

  function layers(grid,enabled)
    type(mg_grid), intent(in) :: grid
    logical, intent(in), optional :: enabled
    type(mg_grid), pointer :: t
    integer :: layers
    layers = 1
    t => grid%child
    do while ( associated(t) )
       if ( present(enabled) ) then
          if ( enabled .eqv. t%enabled ) then
             layers = layers + 1
          end if
       else
          layers = layers + 1
       end if
       t => t%child
    end do
  end function layers

  subroutine grid_delete_layer(in_grid,layer)
    type(mg_grid), intent(inout), target :: in_grid
    integer, intent(in), optional :: layer
    integer :: llayer
    type(mg_grid), pointer :: grid

    grid => in_grid
    llayer = -1
    if ( present(layer) ) llayer = layer
    if ( llayer < 0 ) then
       llayer = layers(in_grid) + 1 + llayer
    end if

    ! we can't kill the first one...
    if ( llayer == 0 ) return
    
    do while ( llayer > 1 )
       grid => grid%child
       if ( .not. associated(grid) ) return

       llayer = llayer - 1

    end do

    ! delete the grid
    call delete_grid(grid)
    grid => grid%parent
    nullify(grid%child%parent)
    deallocate(grid%child)
    nullify(grid%child)
    
  end subroutine grid_delete_layer

  subroutine grid_onoff_layer(in_grid,onoff,layer)
    type(mg_grid), intent(inout), target :: in_grid
    logical, intent(in) :: onoff
    integer, intent(in), optional :: layer
    integer :: llayer
    type(mg_grid), pointer :: grid

    grid => in_grid
    llayer = -1
    if ( present(layer) ) llayer = layer
    if ( llayer < 0 ) then
       llayer = layers(in_grid) + 1 + llayer
    end if

    ! we can't kill the first one...
    if ( llayer == 0 ) return
    
    do while ( llayer > 1 )
       grid => grid%child
       if ( .not. associated(grid) ) return

       llayer = llayer - 1

    end do

    ! delete the grid
    grid%enabled = onoff
    
  end subroutine grid_onoff_layer

  function grid_sum(grid) result(sum)
    type(mg_grid), intent(in) :: grid
    real(grid_p) :: sum
    real(grid_p), pointer :: V(:,:,:)
    integer :: x,y,z

    V => grid%V

    sum = 0._grid_p
!$OMP parallel do default(shared), private(x,y,z), &
!$OMP&   reduction(+:sum), collapse(3)
    do z = 1 , grid%n(3)
    do y = 1 , grid%n(2)
    do x = 1 , grid%n(1)
       sum = sum + abs(V(x,y,z))
    end do
    end do
    end do
!$OMP end parallel do

  end function grid_sum

  function grid_non_constant_elem(grid) result(elem)
    type(mg_grid), intent(in) :: grid
    integer :: x,y,z,elem

    elem = 0
!$OMP parallel do default(shared), private(x,y,z), &
!$OMP&   reduction(+:elem), collapse(3)
    do z = 1 , grid%n(3)
    do y = 1 , grid%n(2)
    do x = 1 , grid%n(1)
       if ( .not. is_constant(grid,x,y,z) ) then
          elem = elem + 1
       end if
    end do
    end do
    end do
!$OMP end parallel do

  end function grid_non_constant_elem

  function grid_layer(grid,layer) result(lg)
    type(mg_grid), intent(in), target :: grid
    integer, intent(in) :: layer
    type(mg_grid), pointer :: lg
    integer :: ll
    lg => grid
    ll = layer
    if ( ll < 1 ) ll = layers(grid,enabled=.true.) + 1 + ll

    if ( ll == 0 ) return
    
    do while ( ll > 1 )
       if ( .not. associated(lg%child) ) return
       lg => lg%child

       ll = ll - 1

    end do

  end function grid_layer

  function grid_tolerance(grid) result(tol)
    type(mg_grid), intent(in) :: grid
    real(grid_p) :: tol, vmax,vmin
    integer :: i
    vmax = -huge(0._grid_p)
    vmin =  huge(0._grid_p)
    do i = 1 , grid%N_box
       vmin = min(vmin,grid%box(i)%val)
       vmax = max(vmax,grid%box(i)%val)
    end do
    tol = grid%tol * abs(vmax - vmin) !/ maxval(grid%n) !grid_non_constant_elem(grid)
  end function grid_tolerance

  subroutine print_grid(top)
    type(mg_grid), intent(in), target :: top

    type(mg_grid), pointer :: grid

    integer :: i, j, k
    real(grid_p) :: n
    character(len=10) :: fmt, bc(2)

    i = 1
    fmt = '(t1,'
    grid => top
    write(*,*) '*******************************************'
    do while ( associated(grid) )
       write(*,trim(fmt)//'a,i0,a,3(i0,tr1))')' -- Layer: ',grid%layer,' size: ',grid%n
       write(*,trim(fmt)//'a,e10.3)')' -- tolerance: ',grid_tolerance(grid)
       write(*,trim(fmt)//'a,f6.4)')' -- SOR: ',grid%sor
       write(*,trim(fmt)//'a,3(tr1,e10.4))')' -- a(3): ',grid%a
       do j = 1 , 3
          do k = 1 , 2
             select case ( grid%BC(k,j)%method ) 
             case ( MG_BC_PERIODIC )
                BC(k) = 'periodic'
             case ( MG_BC_DIRICHLET )
                BC(k) = 'dirichlet'
             case ( MG_BC_NEUMANN )
                BC(k) = 'neumann'
             end select
          end do
          write(*,trim(fmt)//'a,i0,2a,'' : '',a)')' -- BC(',j,'){ ',trim(BC(1)),trim(BC(2))//' }'
       end do

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
    
    write(*,trim(fmt)//'a,e10.3,a,l1,tr1,a,6(tr1,i0),a)')' ++ Box{ value: ',&
         box%val,' constant: ',box%constant, &
         'place: ',box%place,' }'
  end subroutine print_box

end module t_mg
