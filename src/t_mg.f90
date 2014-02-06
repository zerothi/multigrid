module t_mg

  implicit none

  integer, parameter :: dp = selected_real_kind(p=15)
  integer, parameter :: grid_p = selected_real_kind(p=6)

  ! a module to sustain a "simple" multi-grid solver using the 

  type :: mg_grid
     ! the grid information
     logical      :: enabled =.true. ! turn on/off the run on this grid
     real(dp)     :: offset(3) ! the offset of the placement of the local grid 
     real(dp)     :: cell(3,3) ! the cell for the local grid
     real(dp)     :: dL(3,3) ! the cell stepping
     real(dp)     :: dVol ! volume of voxel
     real(dp)     :: Vol  ! volume of cell
     integer      :: n(3) ! size in each direction
     real(grid_p) :: sor  ! the SOR value
     real(grid_p) :: a(3) ! the pre-factors for the summation
     real(grid_p) :: tol  ! the tolerance of the current grid
                          ! this allows different tolerances for different layer-grids
     integer :: itt       ! iterations currently processed
     integer :: steps     ! number of steps in each V-cycle
     integer :: layer     ! the layer that this grid resides in
     real(grid_p),  pointer :: V  (:,:,:) => null() ! update array
     real(grid_p),  pointer :: g  (:) => null() ! ghost arrays ( one long array for all bounds )
     real(grid_p),  pointer :: g_s(:) => null() ! send ghost arrays ( one long array for all bounds )
     type(mg_grid), pointer :: parent => null()
     type(mg_grid), pointer :: child => null()

     ! The constant valued boxes in this grid
     integer :: N_box = 0
     type(mg_box), pointer :: box(:) => null()
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

  subroutine init_grid(grid, n, cell, layer, boxes, tol, offset, sor)
    type(mg_grid), intent(inout) :: grid
    integer,       intent(in)    :: n(3)
    real(dp),      intent(in)    :: cell(3,3)
    integer,       intent(in)    :: layer
    integer,       intent(in)    :: boxes
    real(grid_p),  intent(in), optional :: tol
    real(dp),      intent(in), optional :: offset(3)
    real(grid_p),  intent(in), optional :: sor

    real(dp) :: celll(3), tmp
    integer  :: i

    ! ensure it is empty
    call delete_grid(grid)

    grid%steps = 3

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

    grid%Vol = ( cell(2,1)*cell(3,2) - cell(3,1)*cell(2,2) ) * cell(1,3) + &
         ( cell(3,1)*cell(1,2) - cell(1,1)*cell(3,2) ) * cell(2,3) + &
         ( cell(1,1)*cell(2,2) - cell(2,1)*cell(1,2) ) * cell(3,3)

    grid%dVol = ( grid%dL(2,1)*grid%dL(3,2) - grid%dL(3,1)*grid%dL(2,2) ) * grid%dL(1,3) + &
         ( grid%dL(3,1)*grid%dL(1,2) - grid%dL(1,1)*grid%dL(3,2) ) * grid%dL(2,3) + &
         ( grid%dL(1,1)*grid%dL(2,2) - grid%dL(2,1)*grid%dL(1,2) ) * grid%dL(3,3)

    grid%enabled = .true.
    ! set the values for the grid...
    grid%layer  = layer
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

    tmp = 1._dp / ( 2._dp * sum(celll) ) 
    grid%a(1) = celll(2) * celll(3) * tmp
    grid%a(2) = celll(1) * celll(3) * tmp
    grid%a(3) = celll(1) * celll(2) * tmp
    grid%a = grid%a / sum(grid%a)
    grid%a = 1._grid_p

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
            grid%cell, grid%layer+1, grid%N_box, &
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

      do i = 1 , 3
         n(i) = (grid%n(i) + mod(grid%n(i),2)) / 2 - 1
         if ( n(i) < 7 ) then
            n(:) = 0
            return
         end if
      end do

    end subroutine new_grid_size

  end subroutine init_grid_children_half

  subroutine grid_set(grid,sor,tol,layer)
    type(mg_grid), intent(inout), target :: grid
    real(grid_p), intent(in), optional :: sor, tol
    integer, intent(in), optional :: layer

    type(mg_grid), pointer :: tmp_grid

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
  end subroutine grid_set

  recursive subroutine grid_add_box(grid, llc, box_cell, val, rho, constant)
    type(mg_grid), intent(inout) :: grid
    real(dp), intent(in) :: llc(3), box_cell(3,3)
    real(grid_p), intent(in) :: val, rho
    logical, intent(in) :: constant
    
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
    urc = llc + box_cell(:,1) + box_cell(:,2) + box_cell(:,3)
    ! TODO currently this does not work with skewed axis


!$OMP parallel do default(shared) private(x,y,z,xyz,dyz,dz)
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
       end if
    end do
    end do
    end do
!$OMP end parallel do

    ! if not present, simply delete it...
    if ( all(box%place(1,:) == huge(0)) .and. &
         all(box%place(2,:) == 0) ) then
       call delete_box(box)
    end if

    ! add the box to the child grid
    if ( associated(grid%child) ) then
       call grid_add_box(grid%child,llc,box_cell,val,rho,constant)
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

  subroutine grid_setup(grid)
    type(mg_grid), intent(inout) :: grid
    integer :: x,y,z
    real(grid_p), pointer :: V(:,:,:)

    V => grid%V

    ! set all boxes to their values if constant
!$OMP parallel do default(shared) private(x,y,z) &
!$OMP    collapse(3)
    do z = 1 , grid%n(3)
    do y = 1 , grid%n(2)
    do x = 1 , grid%n(1)
       if ( is_constant(grid,x,y,z) ) then
          V(x,y,z) = val_constant(grid,x,y,z)
       end if
    end do
    end do
    end do
!$OMP end parallel do

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
    allocate(grid%V(grid%n(1),grid%n(2),grid%n(3)))
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

  end subroutine delete_grid

  subroutine grid_restriction(grid)
    type(mg_grid), intent(inout) :: grid

    real(grid_p),  pointer :: V(:,:,:), Vc(:,:,:)
    type(mg_grid), pointer :: child

    real(grid_p), parameter :: f1 = 1._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f4 = 4._grid_p / 64._grid_p
    real(grid_p), parameter :: f8 = 8._grid_p / 64._grid_p

    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%child) ) return

    V  => grid%V
    child => grid%child
    Vc => child%V

    ! initialize the child
!$OMP parallel workshare
    Vc = 0._grid_p
!$OMP end parallel workshare

    ! we employ full-weighting

!$OMP parallel do default(shared) private(x,px,y,py,z,pz)
    do z = 2 , grid%n(3) - 1 , 2
    pz = z / 2 
    do y = 2 , grid%n(2) - 1 , 2
    py = y / 2
    do x = 2 , grid%n(1) - 1 , 2
       px = x / 2

       ! corners
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y-1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y-1,z+1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y+1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y+1,z+1) * f1 ! 4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y-1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y-1,z+1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y+1,z-1) * f1
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y+1,z+1) * f1 ! 8

       ! middles
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y-1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y+1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y,z+1) * f2 ! 4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y-1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y+1,z) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y,z+1) * f2 ! 8
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y-1,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y-1,z+1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y+1,z-1) * f2
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y+1,z+1) * f2 ! 12

       ! neighbours
       Vc(px,py,pz) = Vc(px,py,pz) + V(x-1,y,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x+1,y,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y-1,z) * f4 ! 3
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y+1,z) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y,z-1) * f4
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y,z+1) * f4 ! 6

       ! center
       Vc(px,py,pz) = Vc(px,py,pz) + V(x,y,z) * f8

    end do
    end do
    end do
!$OMP end parallel do

    ! we still need the border...

    ! re-instantiate the constant fields
    call grid_setup(child)

  end subroutine grid_restriction

  subroutine grid_prolongation(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:)
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f2 =   .5_grid_p ! 1 / 2
    real(grid_p), parameter :: f4 =  .25_grid_p ! 1 / 4
    real(grid_p), parameter :: f8 = .125_grid_p ! 1 / 8

    real(grid_p) :: v2, v4, v8

    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! initialize the parent to zero
!$OMP parallel workshare
    Vp = 0._grid_p
!$OMP end parallel workshare

!$OMP parallel do default(shared) private(x,px,y,py,z,pz,v2,v4,v8)
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( parent%n(2) <= py + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py-1,pz+1) = Vp(px-1,py-1,pz+1) + v8
       Vp(px-1,py+1,pz-1) = Vp(px-1,py+1,pz-1) + v8
       Vp(px-1,py+1,pz+1) = Vp(px-1,py+1,pz+1) + v8
       Vp(px+1,py-1,pz-1) = Vp(px+1,py-1,pz-1) + v8
       Vp(px+1,py-1,pz+1) = Vp(px+1,py-1,pz+1) + v8
       Vp(px+1,py+1,pz-1) = Vp(px+1,py+1,pz-1) + v8
       Vp(px+1,py+1,pz+1) = Vp(px+1,py+1,pz+1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v4
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v4
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v4
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v4
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v4
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v4
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v4
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v2
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v2
       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v2

       ! center
       Vp(px,py,pz)   = V(x,y,z)

    end do
    end do
    end do
!$OMP end parallel do

    ! we still need the border

    ! re-instantiate the constant fields
    call grid_setup(parent)

  end subroutine grid_prolongation

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

    ! all values are defaulted to 1
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

    ! all values are defaulted to 1
    rho = 1._grid_p

  end function val_rho

  ! box-routines
  subroutine delete_box(box)
    type(mg_box), intent(inout) :: box
    box%place = 0
    box%val = 1._grid_p
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
!$OMP parallel do default(shared) private(x,y,z) &
!$OMP   reduction(+:sum) collapse(3)
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
!$OMP parallel do default(shared) private(x,y,z) &
!$OMP   reduction(+:elem) collapse(3)
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
    print *,'returning:',lg%layer

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
    tol = grid%tol * abs(vmax - vmin) / grid_non_constant_elem(grid)
  end function grid_tolerance

end module t_mg
