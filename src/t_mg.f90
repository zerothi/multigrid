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
     real(dp)     :: dLL(3) ! the voxel side lengths
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

  subroutine init_grid(grid, n, cell, boxes, tol, offset, sor)
    type(mg_grid), intent(inout) :: grid
    integer,       intent(in)    :: n(3)
    real(dp),      intent(in)    :: cell(3,3)
    integer,       intent(in)    :: boxes
    real(grid_p),  intent(in), optional :: tol
    real(dp),      intent(in), optional :: offset(3)
    real(grid_p),  intent(in), optional :: sor

    real(dp) :: celll(3), tmp
    integer  :: i

    ! ensure it is empty
    call delete_grid(grid)

    grid%steps = 4

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

    ! This additional weighting of each cell
    ! means that different situations might occur
    ! 1. If the "important" direction is along a 
    !    longer voxel direction, the convergence is slower
    ! 2. If the "important" direction is along a 
    !    shorter voxel direction, the convergence is faster
    tmp = 1._dp / ( 2._dp * sum(celll) ) 
    grid%a(1) = celll(2) * celll(3) * tmp
    grid%a(2) = celll(1) * celll(3) * tmp
    grid%a(3) = celll(1) * celll(2) * tmp

    ! this normalization does not really matter, however,
    ! for the convenience of seeing their appropriate weigths
    ! it could be nice to see...
    grid%a = grid%a / sum(grid%a)

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

      do i = 1 , 3
         n(i) = (grid%n(i) + mod(grid%n(i),2)) / 2 - 1
         if ( n(i) < 7 ) then
            n(:) = 0
            return
         end if
      end do

    end subroutine new_grid_size

  end subroutine init_grid_children_half

  subroutine grid_set(grid,sor,tol,layer,steps)
    type(mg_grid), intent(inout), target :: grid
    real(grid_p), intent(in), optional :: sor, tol
    integer, intent(in), optional :: layer, steps

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
    if ( present(steps) ) &
         tmp_grid%steps = steps
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

  subroutine grid_add_point(grid, llc, val, rho, constant)
    type(mg_grid), intent(inout) :: grid
    real(dp), intent(in) :: llc(3)
    real(grid_p), intent(in) :: val, rho
    logical, intent(in) :: constant
    
    real(dp) :: bl(3,3)
    
    bl(:,:) = 0._dp
    call grid_add_box(grid,llc,bl,val,rho,constant)

  end subroutine grid_add_point

  subroutine grid_add_line(grid, llc, dir, l, val, rho, constant)
    type(mg_grid), intent(inout) :: grid
    real(dp), intent(in) :: llc(3), l ! the length of the line
    integer, intent(in) :: dir
    real(grid_p), intent(in) :: val, rho
    logical, intent(in) :: constant
    
    real(dp) :: bl(3,3)
    
    bl = 0._dp
    bl(dir,dir) = l
    call grid_add_box(grid,llc,bl,val,rho,constant)

  end subroutine grid_add_line

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

!$OMP parallel workshare default(shared)
    grid%V   = 0._grid_p
    grid%g   = 0._grid_p
    grid%g_s = 0._grid_p
!$OMP end parallel workshare

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

  subroutine grid_restriction_new(grid)
    type(mg_grid), intent(inout) :: grid

    real(grid_p),  pointer :: V(:,:,:), Vc(:,:,:)
    type(mg_grid), pointer :: child

    real(grid_p), parameter :: f1 = 10._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 5._grid_p / 64._grid_p
    real(grid_p), parameter :: f3 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f4 = 1._grid_p / 64._grid_p

    real(grid_p), pointer :: vt
    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%child) ) return

    V  => grid%V
    child => grid%child
    Vc => child%V

!$OMP parallel default(shared) private(x,px,y,py,z,pz,vt)

    ! initialize the child
!$OMP workshare
    Vc = 0._grid_p
!$OMP end workshare

    ! we employ full-weighting

!$OMP do
    do z = 1 , child%n(3)
    pz = z * 2
    if ( pz >= grid%n(3) ) cycle
    do y = 1 , child%n(2)
    py = y * 2
    if ( py >= grid%n(2) ) cycle
    do x = 1 , child%n(1)
       px = x * 2
       if ( px >= grid%n(1) ) cycle

       ! point to the element
       vt => Vc(x,y,z)

       vt =      V(px,py,pz) * f1

       ! neighbours
       vt = vt + V(px-1,py,pz) * f2
       vt = vt + V(px+1,py,pz) * f2
       vt = vt + V(px,py-1,pz) * f2
       vt = vt + V(px,py+1,pz) * f2
       vt = vt + V(px,py,pz-1) * f2
       vt = vt + V(px,py,pz+1) * f2

       ! next-neighbours
       vt = vt + V(px-1,py-1,pz) * f3
       vt = vt + V(px-1,py+1,pz) * f3
       vt = vt + V(px-1,py,pz-1) * f3
       vt = vt + V(px-1,py,pz+1) * f3

       vt = vt + V(px+1,py-1,pz) * f3
       vt = vt + V(px+1,py+1,pz) * f3
       vt = vt + V(px+1,py,pz-1) * f3
       vt = vt + V(px+1,py,pz+1) * f3

       vt = vt + V(px,py-1,pz-1) * f3
       vt = vt + V(px,py-1,pz+1) * f3
       vt = vt + V(px,py+1,pz-1) * f3
       vt = vt + V(px,py+1,pz+1) * f3

    end do
    end do
    end do
!$OMP end do

!$OMP end parallel

    ! we still need the border...

    ! re-instantiate the constant fields
    call grid_setup(child)

  end subroutine grid_restriction_new

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

!$OMP parallel default(shared) private(x,px,y,py,z,pz)

    ! initialize the child
!$OMP workshare
    Vc = 0._grid_p
!$OMP end workshare

    ! we employ full-weighting

!$OMP do
    do z = 2 , grid%n(3) - 1 , 2
    pz = z / 2 
!    if ( pz > child%n(3) ) cycle ! these checks are needed for different sizes of different layers
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
!$OMP end do

!$OMP end parallel

    ! we still need the border...

    ! re-instantiate the constant fields
    call grid_setup(child)

  end subroutine grid_restriction


  subroutine grid_prolongation(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:)
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f1 = 10._grid_p / 64._grid_p
    real(grid_p), parameter :: f2 = 5._grid_p / 64._grid_p
    real(grid_p), parameter :: f3 = 2._grid_p / 64._grid_p
    real(grid_p), parameter :: f4 = 1._grid_p / 64._grid_p

    real(grid_p), pointer :: vt
    integer :: x,y,z, px,py,pz


    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! initialize the parent to zero
!$OMP parallel default(shared) private(x,px,y,py,z,pz,vt)

!$OMP workshare
    Vp = 0._grid_p
!$OMP end workshare

    ! do middle loop
!$OMP do
    do pz = 1 , parent%n(3)
    z = max(2, pz / 2 - pz / grid%n(3))
    if ( z >= grid%n(3) ) cycle

    do py = 1 , parent%n(2)
    y = max(2,py / 2 - py / grid%n(2))
    if ( y >= grid%n(2) ) cycle

    do px = 1 , parent%n(1)
       x = max(2,px / 2 - px / grid%n(1))
       if ( x >= grid%n(1) ) cycle

       ! point to the vector
       vt => Vp(px,py,pz)

       vt = V(x,y,z) * f1

       ! neighbours
       vt = vt + V(x-1,y,z) * f2
       vt = vt + V(x+1,y,z) * f2
       vt = vt + V(x,y-1,z) * f2
       vt = vt + V(x,y+1,z) * f2
       vt = vt + V(x,y,z-1) * f2
       vt = vt + V(x,y,z+1) * f2

       ! next-neighbours
       vt = vt + V(x-1,y-1,z) * f3
       vt = vt + V(x-1,y+1,z) * f3
       vt = vt + V(x-1,y,z-1) * f3
       vt = vt + V(x-1,y,z+1) * f3

       vt = vt + V(x+1,y-1,z) * f3
       vt = vt + V(x+1,y+1,z) * f3
       vt = vt + V(x+1,y,z-1) * f3
       vt = vt + V(x+1,y,z+1) * f3

       vt = vt + V(x,y-1,z-1) * f3
       vt = vt + V(x,y-1,z+1) * f3
       vt = vt + V(x,y+1,z-1) * f3
       vt = vt + V(x,y+1,z+1) * f3
       
    end do
    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    call grid_setup(parent)

  end subroutine grid_prolongation


  subroutine grid_prolongation_elaborate(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:)
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f2 =   .5_grid_p ! 1 / 2

    integer :: x,y,z, px,py,pz

    integer :: c(3)
    real(grid_p) :: wp(3), wpf(3)
    real(grid_p) :: wxm1, wx0, wxp1
    real(grid_p) :: wym1, wy0, wyp1
    real(grid_p) :: wzm1, wz0, wzp1


    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! create the fraction graphs
    wp = real(parent%n,dp) / real(grid%n,dp)
    wpf = wp - int(wp)
    wp = 1._grid_p

    ! initialize the parent to zero
!$OMP parallel default(shared) private(x,px,y,py,z,pz,c) &
!$OMP    private(wxm1,wx0,wxp1,wym1,wy0,wyp1,wzm1,wz0,wzp1) &
!$OMP    firstprivate(wpf,wp)

!$OMP workshare
    Vp = 0._grid_p
!$OMP end workshare

    ! its a private variable
    c = 1
!$OMP single
    ! setup the correct values for the stuff
    wzm1 = wp(3)
    wz0  = wp(3)
    wzp1 = wp(3)
    wym1 = wp(2)
    wy0  = wp(2)
    wyp1 = wp(2)
    wxm1 = wp(1)
    wx0  = wp(1)
    wxp1 = wp(1)
!$OMP end single

    ! do middle loop
!$OMP do
    do pz = 2 , parent%n(3) - 1
    z = 2 + pz * 2
    if ( z >= grid%n(3) ) cycle
    call setup_pointer(c(3),wp(3),wpf(3),wzm1,wz0,wzp1)
    y = 2
    c(2) = 1
    do py = 2 , parent%n(2) - 1
    y = y + 2 
    if ( y >= grid%n(2) ) cycle
    call setup_pointer(c(2),wp(2),wpf(2),wym1,wy0,wyp1)
    x = 2
    c(1) = 1
    do px = 2 , parent%n(1) - 1
       x = x + 2
       if ( x >= grid%n(1) ) cycle
       call setup_pointer(c(1),wp(1),wpf(1),wxm1,wx0,wxp1)

       if ( px == 6 .and. py == pz .and. py == 2) then
          print *,Vp(px,py,pz)
       end if

       ! corner
       Vp(px,py,pz) = &
            wxm1 * ( V(x-1,y-1,z-1) * wym1 * wzm1 + &
            V(x-1,y-1,z+1) * wym1 * wzp1 + &
            V(x-1,y+1,z-1) * wyp1 * wzm1 + &
            V(x-1,y+1,z+1) * wyp1 * wzp1 ) + &
            wxp1 * ( V(x+1,y-1,z-1) * wym1 * wzm1 + &
            V(x+1,y-1,z+1) * wym1 * wzp1 + &
            V(x+1,y+1,z-1) * wyp1 * wzm1 + &
            V(x+1,y+1,z+1) * wyp1 * wzp1 )

       if ( px == 6 .and. py == pz .and. py == 2) then
          print *,Vp(px,py,pz) * f2 **3
       end if

       ! middle
       Vp(px,py,pz) = f2 * Vp(px,py,pz) + &
            wxm1 * ( V(x-1,y-1,z)  * wym1 * wz0  + &
            V(x-1,y+1,z) * wyp1 * wz0  + &
            V(x-1,y,z-1) * wy0  * wzm1 + &
            V(x-1,y,z+1) * wy0  * wzp1 ) + &
            wxp1 * ( V(x+1,y-1,z) * wym1 * wz0  + &
            V(x+1,y+1,z) * wyp1 * wz0  + &
            V(x+1,y,z-1) * wy0  * wzm1 + &
            V(x+1,y,z+1) * wy0  * wzp1 ) + &
            wx0 * ( V(x,y-1,z-1) * wym1 * wzm1 + &
            V(x,y-1,z+1) * wym1 * wzp1 + &
            V(x,y+1,z-1) * wyp1 * wzm1 + &
            V(x,y+1,z+1) * wyp1 * wzp1 )

       if ( px == 6 .and. py == pz .and. py == 2) then
          print *,Vp(px,py,pz) * f2 ** 2
       end if

       ! middle
       Vp(px,py,pz) = f2 * Vp(px,py,pz) + &
            V(x-1,y,z) * wxm1 * wy0  * wz0  + &
            V(x+1,y,z) * wxp1 * wy0  * wz0  + &
            wx0 * ( V(x,y-1,z) * wym1 * wz0  + &
            V(x,y+1,z) * wyp1 * wz0  + &
            V(x,y,z-1) * wy0  * wzm1 + &
            V(x,y,z+1) * wy0  * wzp1 )

       if ( px == 6 .and. py == pz .and. py == 2) then
          print *,Vp(px,py,pz) * f2
       end if

       ! center
       Vp(px,py,pz) = f2 * Vp(px,py,pz) + V(x,y,z) * wx0 * wy0 * wz0
       
    end do
    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    call grid_setup(parent)

  contains

    subroutine setup_pointer(c,wp,wpf,wm1,w0,wp1)
      integer, intent(inout) :: c
      real(grid_p), intent(in) :: wp, wpf
      real(grid_p), intent(out) :: wm1, w0, wp1
      real(grid_p) :: w

      if ( c == 1 ) then
         w = wm1
         wm1 = w0
         w0  = wp1
         wp1 = w - wpf
         c = 2
      else if ( c == 2 ) then
         wm1 = wp
         w0  = wp
         wp1 = wpf
         c = 3
      else
         wm1 = wp
         w0  = wp
         wp1 = wpf
         c = 1
      end if

    end subroutine setup_pointer

  end subroutine grid_prolongation_elaborate

  subroutine grid_prolongation_old(grid)
    type(mg_grid), intent(inout), target :: grid

    real(grid_p), pointer :: Vp(:,:,:), V(:,:,:)
    type(mg_grid), pointer :: parent

    real(grid_p), parameter :: f2 =   .5_grid_p ! 1 / 2
    real(grid_p), parameter :: f4 =  .25_grid_p ! 1 / 4
    real(grid_p), parameter :: f8 = .125_grid_p ! 1 / 8

    real(grid_p) :: v2, v4, v8

    logical :: one_larger(3)
    integer :: x,y,z, px,py,pz

    ! if the child does not exist, then return immediately
    if ( .not. associated(grid%parent) ) return

    V  => grid%V
    parent => grid%parent
    Vp => parent%V

    ! currently this prolongation only accepts
    ! grids of double +1/+2 sizes
    do x = 1 , 3
       select case ( parent%n(x) - grid%n(x) * 2 )
       case ( 1 )
          one_larger(x) = .true.
       case ( 2 )
          one_larger(x) = .false.
       case default
          stop 'Not a functioning grid step...'
       end select
    end do
          
    ! initialize the parent to zero
!$OMP parallel default(shared) private(x,px,y,py,z,pz,v2,v4,v8)

!$OMP workshare
    Vp = 0._grid_p
!$OMP end workshare

    ! ensure to set the top corners correctly

    z  = 1
    pz = 1
!$OMP do 
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( one_larger(2) ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v8
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
!       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px+1,py+1,pz,'x1')
!call print_t(px-1,py-1,pz,'x1')
!call print_t(px-1,py+1,pz,'x1')
!call print_t(px+1,py-1,pz,'x1')
!call print_t(px,py,pz,'x1')

    end do
    end do
!$OMP end do
    if ( grid%n(2)*2 + 1 /= parent%n(2) ) then
    y  = grid%n(2)
    py = parent%n(2)
!$OMP do 
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v8
!       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
!       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px-1,py-1,pz,'x1')
!call print_t(px+1,py-1,pz,'x1')
!call print_t(px,py,pz,'x1')

    end do
!$OMP end do
    end if

    y  = 1
    py = 1
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( one_larger(3) ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v8
!       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px-1,py,pz-1,'y1')
!call print_t(px-1,py,pz+1,'y1')
!call print_t(px+1,py,pz-1,'y1')
!call print_t(px+1,py,pz+1,'y1')
!call print_t(px,py,pz,'y1')

    end do
    end do
!$OMP end do

    x  = 1
    px = 1
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( one_larger(3) ) cycle
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
       py = 2 * y
       if ( one_larger(2) ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
!       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
!       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v8
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

!call print_t(px,py-1,pz-1,'z1')
!call print_t(px,py-1,pz+1,'z1')
!call print_t(px,py+1,pz-1,'z1')
!call print_t(px,py+1,pz+1,'z1')
!call print_t(px,py,pz,'z1')

    end do
    end do
!$OMP end do
    
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( one_larger(3) ) cycle
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( one_larger(2) ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( one_larger(1) ) cycle

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
!$OMP end do

    z = grid%n(3)
    pz = parent%n(3)
    if ( .not. one_larger(3) ) then
!$OMP do 
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
       Vp(px-1,py+1,pz-1) = Vp(px-1,py+1,pz-1) + v8
       Vp(px+1,py-1,pz-1) = Vp(px+1,py-1,pz-1) + v8
       Vp(px+1,py+1,pz-1) = Vp(px+1,py+1,pz-1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v4
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v4
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v4
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v2
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v2

    end do
    end do
!$OMP end do
    end if
!$OMP do 
    do y = 1 , grid%n(2)
    py = 2 * y
    if ( parent%n(2) <= py + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v8
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8

       ! middles
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v8
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v8
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v8
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v8

       ! neighbours
       Vp(px,py,pz) = Vp(px,py,pz) + v8

    end do
    end do
!$OMP end do

    y = grid%n(2)
    py = parent%n(2)
    if ( y * 2 /= py ) then
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py-1,pz+1) = Vp(px-1,py-1,pz+1) + v8
       Vp(px+1,py-1,pz-1) = Vp(px+1,py-1,pz-1) + v8
       Vp(px+1,py-1,pz+1) = Vp(px+1,py-1,pz+1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v4
       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v4
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v4
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v2
       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v2

    end do
    end do
!$OMP end do
    end if
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do x = 1 , grid%n(1)
       px = 2 * x
       if ( parent%n(1) <= px + 1 ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v8
!       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v8
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v8
       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
!       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
!       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

    end do
    end do
!$OMP end do

    x = grid%n(1)
    px = parent%n(1)
    if ( x * 2 /= px ) then
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
       py = 2 * y
       if ( parent%n(2) <= py + 1 ) cycle

       v2 = f2 * V(x,y,z)
       v4 = f4 * V(x,y,z)
       v8 = f8 * V(x,y,z)

       ! corners
       Vp(px-1,py-1,pz-1) = Vp(px-1,py-1,pz-1) + v8
       Vp(px-1,py-1,pz+1) = Vp(px-1,py-1,pz+1) + v8
       Vp(px-1,py+1,pz-1) = Vp(px-1,py+1,pz-1) + v8
       Vp(px-1,py+1,pz+1) = Vp(px-1,py+1,pz+1) + v8

       ! middles
       Vp(px-1,py-1,pz) = Vp(px-1,py-1,pz) + v4
       Vp(px-1,py,pz-1) = Vp(px-1,py,pz-1) + v4
       Vp(px-1,py+1,pz) = Vp(px-1,py+1,pz) + v4
       Vp(px-1,py,pz+1) = Vp(px-1,py,pz+1) + v4
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v4
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v4
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v4
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v4

       ! neighbours
       Vp(px-1,py,pz) = Vp(px-1,py,pz) + v2
       Vp(px,py-1,pz) = Vp(px,py-1,pz) + v2
       Vp(px,py,pz-1) = Vp(px,py,pz-1) + v2
       Vp(px,py+1,pz) = Vp(px,py+1,pz) + v2
       Vp(px,py,pz+1) = Vp(px,py,pz+1) + v2

    end do
    end do
!$OMP end do
    end if
!$OMP do 
    do z = 1 , grid%n(3)
    pz = 2 * z
    if ( parent%n(3) <= pz + 1 ) cycle
    do y = 1 , grid%n(2)
       py = 2 * y
       if ( parent%n(2) <= py + 1 ) cycle

       v8 = f8 * V(x,y,z)

       ! middles
!       Vp(px+1,py-1,pz) = Vp(px+1,py-1,pz) + v8
!       Vp(px+1,py,pz-1) = Vp(px+1,py,pz-1) + v8
!       Vp(px+1,py+1,pz) = Vp(px+1,py+1,pz) + v8
!       Vp(px+1,py,pz+1) = Vp(px+1,py,pz+1) + v8
       Vp(px,py-1,pz-1) = Vp(px,py-1,pz-1) + v8
       Vp(px,py-1,pz+1) = Vp(px,py-1,pz+1) + v8
       Vp(px,py+1,pz-1) = Vp(px,py+1,pz-1) + v8
       Vp(px,py+1,pz+1) = Vp(px,py+1,pz+1) + v8

       ! neighbours
!       Vp(px+1,py,pz) = Vp(px+1,py,pz) + v8
       Vp(px,py,pz)   = Vp(px,py,pz) + v8

    end do
    end do
!$OMP end do nowait

!$OMP end parallel

    ! we still need the border

    ! re-instantiate the constant fields
    call grid_setup(parent)

  contains

    subroutine print_t(x,y,z,m)
      integer, intent(in) :: x,y,z
      character(len=*), intent(in) :: m
      if ( x == 14 .and. y == 14 .and. z == 1 ) then
         print *,'we are here', m
      end if
    end subroutine print_t

  end subroutine grid_prolongation_old

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
    tol = grid%tol * abs(vmax - vmin) / maxval(grid%n) !grid_non_constant_elem(grid)
  end function grid_tolerance

  function overlap(tg,tx,ty,tz,bg,bx,by,bz) result(lap)
    type(mg_grid), intent(in) :: tg, bg
    integer, intent(in) :: tx,ty,tz,bx,by,bz
    real(grid_p) :: lap
    real(grid_p) :: tll(3), bll(3)

    tll(:) = tg%offset(:) &
         + tg%dL(:,1) * (tx-1) &
         + tg%dL(:,2) * (ty-1) &
         + tg%dL(:,3) * (tz-1)
    bll(:) = bg%offset(:) &
         + bg%dL(:,1) * (bx-1) &
         + bg%dL(:,2) * (by-1) &
         + bg%dL(:,3) * (bz-1)

    ! TODO currently this overlap only works for rectangles...
!    print '(2(3(3(tr1,e10.3),/)))',tll,tg%dLL,bll,bg%dLL
!    print *, max(0._grid_p, max(tll(1)+tg%dLL(1),bll(1)+bg%dLL(1)) - min(tll(1),bll(1)))
!    print *, max(0._grid_p, max(tll(3)+tg%dLL(3),bll(3)+bg%dLL(3)) - min(tll(3),bll(3)))
    lap = &
         max(0._grid_p, max(tll(1)+tg%dLL(1),bll(1)+bg%dLL(1)) - min(tll(1),bll(1))) * & ! x
         max(0._grid_p, max(tll(2)+tg%dLL(2),bll(2)+bg%dLL(2)) - min(tll(2),bll(2))) * & ! y
         max(0._grid_p, max(tll(3)+tg%dLL(3),bll(3)+bg%dLL(3)) - min(tll(3),bll(3)))     ! z

 !   print *,tg%dVol,bg%dVol,lap
!    lap = max(0._grid_p,min(1._grid_p,lap / ( tg%dVol + bg%dVol - lap )))
    lap = lap / ( tg%dVol + bg%dVol - lap )
  end function overlap

end module t_mg
