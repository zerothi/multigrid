module m_gs_CDS

  use t_mg
  
  implicit none

  private

  integer, parameter :: MG_METHOD_GS_TEMPORAL_CDS = 1

  integer, parameter :: CDS_BOTTOM_UP = 1
  integer, parameter :: CDS_W_CYCLE = 2

  public :: mg_gs_cds
  public :: MG_METHOD_GS_TEMPORAL_CDS
  public :: CDS_BOTTOM_UP, CDS_W_CYCLE

contains


  subroutine mg_gs_cds(top_grid,method)
    type(mg_grid), intent(inout), target :: top_grid
    integer, intent(in), optional :: method
    integer :: lmethod 

    lmethod = CDS_BOTTOM_UP
    if ( present(method) ) lmethod = method

    select case ( lmethod )
    case ( CDS_W_CYCLE )
       
       call gs_W(top_grid)

    case ( CDS_BOTTOM_UP )

       call gs_bottom_up(top_grid)
       
    case default

       stop 'Error on method used for CDS'

    end select
    
  end subroutine mg_gs_cds

  subroutine gs_bottom_up(top_grid)
    type(mg_grid), intent(inout), target :: top_grid
    type(mg_grid), pointer :: grid
    integer :: old_itt, n(3)
    real(grid_p) :: tol, old_sum, new_sum

    grid => top_grid

    do while ( associated(grid%child) ) 

       ! allocate space for the child
       call grid_bring_back(grid%child)
       
       ! restrict data to child
       call grid_restriction(grid)

       ! this will hold-back the grid
       call grid_hold_back(grid)
       
       ! next
       grid => grid%child

    end do

    do while ( associated(grid) ) 

       old_itt = grid%itt

       ! solve the grid
       call grid_solve(grid)

       ! print out number of iterations used on that cycle
       write(*,'(2(a,i0),a)') &
            'Completed (',grid%layer,') cycle in ',grid%itt-old_itt,' cycles'

       ! bring back data
       if ( associated(grid%parent) ) then
          call grid_bring_back(grid%parent)
          call grid_prolongation(grid)
          call grid_hold_back(grid)
       end if

       grid => grid%parent

    end do

  end subroutine gs_bottom_up

  subroutine gs_w(top_grid)
    type(mg_grid), intent(inout), target :: top_grid
    type(mg_grid), pointer :: pg, cg
    integer :: old_itt, n(3)
    real(grid_p) :: nr, old_sum, new_sum
    real(grid_p), target :: tol, itol
    real(grid_p), pointer :: otol

    otol => tol
    pg => top_grid
    cg => top_grid

    ! first we move all the way down to the second lowest...
    do while ( associated(pg%child%child) )
       if ( .not. pg%child%child%enabled ) exit
       call grid_bring_back(pg%child)
       call grid_restriction(pg)
       call grid_hold_back(pg)
       pg => pg%child
       cg => pg%child
    end do

    ! continue until we are under tolerance
    otol = grid_tolerance(pg) + 1._grid_p
    do while ( otol > grid_tolerance(top_grid) )     

       nr = 1._grid_p / grid_non_constant_elem(pg)

       old_sum = grid_sum(pg)

       old_itt = pg%itt

       print '(2(a,i0))','Running between ',pg%layer,' and ',cg%layer

       ! We must have a double loop...
       itol = grid_tolerance(pg) + 1._grid_p
       do while ( itol > grid_tolerance(pg) )

          call gs_V(pg,cg)

          new_sum = grid_sum(pg)
          itol = abs(old_sum - new_sum) * nr
          old_sum = new_sum

          print '(a,3(tr1,e10.3))',' Tol,new_sum',itol,new_sum,grid_tolerance(pg)

       end do

       print '(3(a,i0),a)','Completed: ',pg%layer,':',cg%layer,' in ', &
            pg%itt-old_itt,' itt. per. lvl'

       ! if we are the top layer, limit it to only one layer
       if ( pg%layer == top_grid%layer ) then
          cg   => pg
          itol =  otol
          otol => itol
       end if

       ! step up...
       if ( associated(pg%parent) ) then
          call grid_bring_back(pg%parent)
          call grid_prolongation(pg)
          call grid_hold_back(pg)
          cg => pg
          pg => pg%parent
       end if

    end do

  end subroutine gs_w

  subroutine grid_solve(grid)
    type(mg_grid), intent(inout) :: grid
    real(grid_p) :: tol, old_sum, new_sum
    real(grid_p) :: nr
    integer :: old_itt

    if ( .not. grid%enabled ) return

    nr = 1._grid_p / (grid_non_constant_elem(grid))

    tol = grid_tolerance(grid) + 1._grid_p
    
    old_sum = grid_sum(grid) ! * nr

    old_itt = grid%itt

    do while ( tol > grid_tolerance(grid) ) 

       call gs_step(grid)

       ! the tolerance is the difference between the 
       ! sum of the before iteration and the new iteration
       new_sum = grid_sum(grid) ! * nr
       tol = abs(old_sum - new_sum) * nr
       old_sum = new_sum
       print '(a,3(tr1,e10.3))',' Tol,new_sum',tol,new_sum,grid_tolerance(grid)

       ! print out number of iterations used on that cycle
       !write(*,'(2(a,i0),a)') &
       !     'Completed (',grid%layer,') cycle in ',grid%itt-old_itt,' cycles'

    end do

  end subroutine grid_solve

  subroutine gs_V(pg,cg)
    type(mg_grid), intent(inout), target :: pg, cg
    type(mg_grid), pointer :: grid
    integer :: i,n(3), cur_layer
    
    grid => pg

    ! move down
    do while ( grid%layer /= cg%layer ) 

       ! if the grid is not enabled we immediately exit
       ! the restriction cycle
       if ( .not. grid%child%enabled ) exit

       ! Do two steps
       do i = 1 , grid%steps
          call gs_step(grid)
       end do
       
       ! allocate space for the child
       call grid_bring_back(grid%child)
       
       ! restrict data to child
       call grid_restriction(grid)

       ! this will hold-back the grid
       call grid_hold_back(grid)

       ! next
       grid => grid%child

    end do

    ! go back in cycle...

    do
       
       do i = 1 , grid%steps
          call gs_step(grid)
       end do

       if ( grid%layer == pg%layer ) exit
       
       ! bring back data
       if ( associated(grid%parent) ) then
          call grid_bring_back(grid%parent)
          call grid_prolongation(grid)
          call grid_hold_back(grid)
       end if

       grid => grid%parent
       
    end do
    
  end subroutine gs_V

  subroutine gs_step(grid)
    type(mg_grid), intent(inout) :: grid

    grid%itt = grid%itt + 1

    !< communicate bounds >

    call gs(grid)
    !< wait communicate >       

    call gs_bound(grid)

  end subroutine gs_step

  subroutine gs(grid)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2)
    integer :: x,y,z

    V => grid%V

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

!$OMP parallel do default(shared) collapse(3) &
!$OMP   private(x,y,z) firstprivate(sor)
    do z = 2 , grid%n(3) - 1
       do y = 2 , grid%n(2) - 1
          do x = 2 , grid%n(1) - 1
             ! calculate the current contribution
             V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * val(grid,V,x,y,z)
          end do
       end do
    end do
!$OMP end parallel do
    
  end subroutine gs

  subroutine gs_bound(grid)
    type(mg_grid), intent(inout) :: grid
    
    ! currently we can't handle 2D stuff
    call gs_xb(grid,1)
    call gs_xb(grid,grid%n(1))
    call gs_yb(grid,1)
    call gs_yb(grid,grid%n(2))
    call gs_zb(grid,1)
    call gs_zb(grid,grid%n(3))

  end subroutine gs_bound
  
  subroutine gs_xb(grid,x)
    ! this routine calculates the contribution on the
    ! lower/upper x-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: x
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2), val_r(4)
    integer :: dx,y,z

    V => grid%V
    if ( x == 1 ) then
       dx = 1 
    else
       dx = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)
    
!$OMP parallel default(shared) &
!$OMP   private(y,z,val_r) firstprivate(sor,x)

    ! x, y, z-corners (notice that there are 8 corners)
!$OMP sections
!$OMP section
    call gs_corner(grid,V,sor,x, 1,1, dx, 1, 1)
!$OMP section
    call gs_corner(grid,V,sor,x, grid%n(2),1, dx, -1,  1)
!$OMP section
    call gs_corner(grid,V,sor,x, 1, grid%n(3), dx, 1, -1)
!$OMP section
    call gs_corner(grid,V,sor,x, grid%n(2),grid%n(3), dx, -1, -1)
!$OMP end sections nowait

! consider adding single, and add dynamic scheduling for the "big" loop

    ! take the y plane at (x,z) = (1|grid%n(1),1)

    z = 1
!$OMP do
    do y = 2 , grid%n(2) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z)
          val_r(2) = val_rho(grid,x,y-1,z)
          val_r(3) = val_rho(grid,x,y+1,z)
          val_r(4) = val_rho(grid,x,y,z+1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z+1) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

!$OMP do
    do z = 2 , grid%n(3) - 1
       ! calculate the contribution in y = 1
       ! notice that z is different from each thread,
       ! hence we need not "critical"
       y = 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z)
          val_r(2) = val_rho(grid,x,y+1,z)
          val_r(3) = val_rho(grid,x,y,z-1)
          val_r(4) = val_rho(grid,x,y,z+1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y+1,z) * val_r(2) ) + &
               grid%a(3) * ( V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) ) )
       end if
       do y = 2 , grid%n(2) - 1
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * val_xb(grid,V,x,y,z,dx)
       end do
       y = grid%n(2)
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z)
          val_r(2) = val_rho(grid,x,y-1,z)
          val_r(3) = val_rho(grid,x,y,z-1)
          val_r(4) = val_rho(grid,x,y,z+1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y-1,z) * val_r(2) ) + &
               grid%a(3) * ( V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

    ! take the y plane at (x,z) = (1|grid%n(1),grid%n(3))

    z = grid%n(3)
!$OMP do
    do y = 2 , grid%n(2) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z)
          val_r(2) = val_rho(grid,x,y-1,z)
          val_r(3) = val_rho(grid,x,y+1,z)
          val_r(4) = val_rho(grid,x,y,z-1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z-1) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

!$OMP end parallel

  end subroutine gs_xb

  subroutine gs_yb(grid,y)
    ! this routine calculates the contribution on the
    ! lower/upper y-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: y
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2), val_r(4)
    integer :: x,dy,z

    V => grid%V
    if ( y == 1 ) then
       dy = 1 
    else
       dy = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

!$OMP parallel default(shared) &
!$OMP   private(x,z,val_r) firstprivate(sor,y)

    ! take the x plane at (y,z) = (1|grid%n(2),1)

    z = 1
!$OMP do
    do x = 2 , grid%n(1) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z)
          val_r(2) = val_rho(grid,x+1,y,z)
          val_r(3) = val_rho(grid,x,y+dy,z)
          val_r(4) = val_rho(grid,x,y,z+1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
               grid%a(2) * ( V(x,y+dy,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z+1) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait
    
!$OMP do
    do z = 2 , grid%n(3) - 1
       ! calculate the contribution in y = 1
       ! notice that z is different from each thread,
       ! hence we need not "critical"
       x = 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+1,y,z)
          val_r(2) = val_rho(grid,x,y+dy,z)
          val_r(3) = val_rho(grid,x,y,z-1)
          val_r(4) = val_rho(grid,x,y,z+1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x+1,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y+dy,z) * val_r(2) ) + &
               grid%a(3) * ( V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) ) )
       end if

       do x = 2 , grid%n(1) - 1
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * val_yb(grid,V,x,y,z,dy)
       end do

       x = grid%n(1)
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z)
          val_r(2) = val_rho(grid,x,y+dy,z)
          val_r(3) = val_rho(grid,x,y,z-1)
          val_r(4) = val_rho(grid,x,y,z+1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x-1,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y+dy,z) * val_r(2) ) + &
               grid%a(3) * ( V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

    ! take the x plane at (y,z) = (1|grid%n(2),grid%n(3))

    z = grid%n(3)
!$OMP do
    do x = 2 , grid%n(1) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z)
          val_r(2) = val_rho(grid,x+1,y,z)
          val_r(3) = val_rho(grid,x,y+dy,z)
          val_r(4) = val_rho(grid,x,y,z-1)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
               grid%a(2) * ( V(x,y+dy,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z-1) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

!$OMP end parallel

  end subroutine gs_yb

  subroutine gs_zb(grid,z)
    ! this routine calculates the contribution on the
    ! lower/upper z-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: z
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2), val_r(4)
    integer :: x,y,dz

    V => grid%V
    if ( z == 1 ) then
       dz = 1 
    else
       dz = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

!$OMP parallel default(shared) &
!$OMP   private(x,y,val_r) firstprivate(sor,z)

!$OMP do
    do y = 2 , grid%n(2) - 1
       x = 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+1,y,z)
          val_r(2) = val_rho(grid,x,y-1,z)
          val_r(3) = val_rho(grid,x,y+1,z)
          val_r(4) = val_rho(grid,x,y,z+dz)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x+1,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z+dz) * val_r(4) ) )
       end if

       do x = 2 , grid%n(1) - 1
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * val_zb(grid,V,x,y,z,dz)
       end do

       x = grid%n(1)
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z)
          val_r(2) = val_rho(grid,x,y-1,z)
          val_r(3) = val_rho(grid,x,y+1,z)
          val_r(4) = val_rho(grid,x,y,z+dz)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x-1,y,z) * val_r(1) ) + &
               grid%a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z+dz) * val_r(4) ) )
       end if

    end do
!$OMP end do nowait

!$OMP end parallel

  end subroutine gs_zb

  subroutine gs_corner(grid,V,sor,x,y,z,dx,dy,dz)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(inout) :: V(:,:,:)
    real(grid_p), intent(in) :: sor(2)
    integer, intent(in) :: x,y,z,dx,dy,dz
    real(grid_p) :: vcur, val_r(3)
    if ( is_constant(grid,x,y,z) ) return
    val_r(1) = val_rho(grid,x+dx,y,z)
    val_r(2) = val_rho(grid,x,y+dy,z)
    val_r(3) = val_rho(grid,x,y,z+dz)
    val_r = val_r / sum(val_r)
    vcur = &
         grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
         grid%a(2) * ( V(x,y+dy,z) * val_r(2) ) + &
         grid%a(3) * ( V(x,y,z+dz) * val_r(3) )
    V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
  end subroutine gs_corner

  pure function val(grid,V,x,y,z)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: V(:,:,:)
    integer, intent(in) :: x,y,z
    real(grid_p) :: val, val_r(6)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return
    
    val_r(1) = val_rho(grid,x-1,y,z)
    val_r(2) = val_rho(grid,x+1,y,z)
    val_r(3) = val_rho(grid,x,y-1,z)
    val_r(4) = val_rho(grid,x,y+1,z)
    val_r(5) = val_rho(grid,x,y,z-1)
    val_r(6) = val_rho(grid,x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         grid%a(2) * ( V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) ) + &
         grid%a(3) * ( V(x,y,z-1) * val_r(5) + V(x,y,z+1) * val_r(6) )

  end function val

  pure function val_xb(grid,V,x,y,z,dx) result(val)
    ! calculates the contribution on the lower/upper x-bound
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dx
    real(grid_p) :: val, val_r(5)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return

    val_r(1) = val_rho(grid,x+dx,y,z)
    val_r(2) = val_rho(grid,x,y-1,z)
    val_r(3) = val_rho(grid,x,y+1,z)
    val_r(4) = val_rho(grid,x,y,z-1)
    val_r(5) = val_rho(grid,x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
         grid%a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
         grid%a(3) * ( V(x,y,z-1) * val_r(4) + V(x,y,z+1) * val_r(5) )

  end function val_xb

  pure function val_yb(grid,V,x,y,z,dy) result(val)
    type(mg_grid), intent(in) :: grid
    ! calculates the contribution on the lower/upper x-bound
    real(grid_p), intent(in) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dy
    real(grid_p) :: val, val_r(5)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return

    val_r(1) = val_rho(grid,x-1,y,z)
    val_r(2) = val_rho(grid,x+1,y,z)
    val_r(3) = val_rho(grid,x,y+dy,z)
    val_r(4) = val_rho(grid,x,y,z-1)
    val_r(5) = val_rho(grid,x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         grid%a(2) * ( V(x,y+dy,z) * val_r(3) ) + &
         grid%a(3) * ( V(x,y,z-1) * val_r(4) + V(x,y,z+1) * val_r(5) )
    
  end function val_yb

  pure function val_zb(grid,V,x,y,z,dz) result(val)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dz
    real(grid_p) :: val, val_r(5)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return

    val_r(1) = val_rho(grid,x-1,y,z)
    val_r(2) = val_rho(grid,x+1,y,z)
    val_r(3) = val_rho(grid,x,y-1,z)
    val_r(4) = val_rho(grid,x,y+1,z)
    val_r(5) = val_rho(grid,x,y,z+dz)
    val_r = val_r / sum(val_r)

    val = &
         grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         grid%a(2) * ( V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) ) + &
         grid%a(3) * ( V(x,y,z+dz) * val_r(5) )

  end function val_zb

end module m_gs_CDS
