module m_gs_CDS

  use t_mg
  use t_mg_interp
  
  implicit none

  private

  integer, parameter :: MG_METHOD_GS_TEMPORAL_CDS = 1

  integer, parameter :: CDS_BOTTOM_UP = 1
  integer, parameter :: CDS_W_CYCLE = 2

  public :: mg_gs_cds
  public :: MG_METHOD_GS_TEMPORAL_CDS
  public :: CDS_BOTTOM_UP, CDS_W_CYCLE

contains

  subroutine mg_gs_cds(top_grid,method, init)
    type(mg_grid), intent(inout), target :: top_grid
    integer, intent(in), optional :: method
    logical, intent(in), optional :: init
    integer :: lmethod 

    lmethod = CDS_BOTTOM_UP
    if ( present(method) ) lmethod = method

    ! We default to initialization of the variables and
    ! the arrays.
    ! The user may say "no" to initialization which 
    ! enables the user to change the initial guess
    if ( .not. present(init) ) then
       call grid_bring_back(top_grid)
       call grid_setup(top_grid, init = .true. )
    else if ( init ) then
       call grid_bring_back(top_grid)
       call grid_setup(top_grid, init = .true. )
    end if

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

       ! next
       grid => grid%child

       if ( .not. associated(grid%child) ) then

          call grid_bring_back(grid)

          call grid_setup(grid , init = .false. )

       end if

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

          print '(a,i6,3(tr1,f10.7))',' itt ',pg%itt-old_itt, &
               itol,new_sum*nr,pg%err

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
    
    old_sum = grid_sum(grid)
    print'(a,f10.7)','Initial sum: ',old_sum * nr

    old_itt = grid%itt

    do while ( tol > grid_tolerance(grid) ) 

       call gs_step(grid)

       ! the tolerance is the difference between the 
       ! sum of the before iteration and the new iteration
       new_sum = grid_sum(grid)
       tol = abs(old_sum - new_sum) * nr
       old_sum = new_sum
       print '(a,i6,3(tr1,f10.7))',' itt ',grid%itt-old_itt, &
            tol,new_sum*nr,grid%err

       ! print out number of iterations used on that cycle
       !write(*,'(2(a,i0),a)') &
       !     'Completed (',grid%layer,') cycle in ',grid%itt-old_itt,' cycles'

    end do

  end subroutine grid_solve

  subroutine gs_V(pg,cg)
    type(mg_grid), intent(inout), target :: pg, cg
    type(mg_grid), pointer :: grid
    integer :: i, n(3), cur_layer, ig
    real(grid_p) :: tol, old_sum, new_sum
    real(grid_p) :: nr
    
    grid => pg

    ig = 0
    ! move down
    do while ( grid%layer /= cg%layer ) 

       ! if the grid is not enabled we immediately exit
       ! the restriction cycle
       if ( .not. grid%child%enabled ) exit

       nr = 1._grid_p / (grid_non_constant_elem(grid))
       old_sum = grid_sum(grid)

       do i = 1 , grid%steps
          call gs_step(grid)
       end do

       ig = ig + 1
       new_sum = grid_sum(grid)
       tol = abs(old_sum - new_sum) * nr
       print '(tr6,a,i6,2(tr1,f10.7))',repeat(' ',ig * 2), &
            grid%steps,tol,new_sum*nr
       
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
    ig = ig + 1
    do
       
       nr = 1._grid_p / (grid_non_constant_elem(grid))
       old_sum = grid_sum(grid)

       do i = 1 , grid%steps
          call gs_step(grid)
       end do

       new_sum = grid_sum(grid)
       tol = abs(old_sum - new_sum) * nr
       print '(tr6,a,i6,2(tr1,f10.7))',repeat(' ',ig * 2), &
            grid%steps,tol,new_sum*nr
       ig = ig - 1

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
    real(grid_p) :: err

    ! initialize the max error
    grid%err = 0._grid_p

    grid%itt = grid%itt + 1

    !< communicate bounds >

    call gs(grid)
    !< wait communicate >       

!    call gs_bound(grid)

    ! Take the square root of the abs error
    grid%err = sqrt(grid%err)

  end subroutine gs_step

  subroutine gs(grid)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2), err, vv
    integer :: x,y,z

    V => grid%V

    err = 0._grid_p
    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

!$OMP parallel do default(shared), collapse(3), &
!$OMP&private(x,y,z,vv), firstprivate(sor), reduction(max:err)
    do z = 1 , grid%n(3)
       do y = 1 , grid%n(2)
          do x = 1 , grid%n(1)
             ! calculate the current contribution
             vv = sor(1) * V(x,y,z) + sor(2) * val(grid,V,x,y,z)
             err = max(err,(vv-V(x,y,z))**2)
             V(x,y,z) = vv
          end do
       end do
    end do
!$OMP end parallel do

    grid%err = max(grid%err,err)
    
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
    real(grid_p) :: sor(2), val_r(4), vv, err
    integer :: dx,y,z

    V => grid%V
    if ( x == 1 ) then
       dx = 1 
    else
       dx = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)
    
!$OMP parallel default(shared), private(y,z,val_r,vv,err), &
!$OMP firstprivate(sor,x)

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

    err = 0._grid_p

    z = 1
!$OMP do
    do y = 2 , grid%n(2) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y-1,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y+1,z) * grid%a(2)
          val_r(4) = val_rho(grid,x,y,z+1) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x+dx,y,z) * val_r(1) + &
               V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) + &
               V(x,y,z+1) * val_r(4) )
          
          err = max(err,(vv-V(x,y,z))**2)

          V(x,y,z) = vv
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
          val_r(1) = val_rho(grid,x+dx,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y+1,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y,z-1) * grid%a(3)
          val_r(4) = val_rho(grid,x,y,z+1) * grid%a(3)
          val_r = val_r / sum(val_r)

          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x+dx,y,z) * val_r(1) + &
               V(x,y+1,z) * val_r(2) + &
               V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) )
          
          err = max(err,(vv-V(x,y,z))**2)

          V(x,y,z) = vv
       end if
       do y = 2 , grid%n(2) - 1
          vv = sor(1) * V(x,y,z) + sor(2) * val_xb(grid,V,x,y,z,dx)
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end do
       y = grid%n(2)
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y-1,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y,z-1) * grid%a(3)
          val_r(4) = val_rho(grid,x,y,z+1) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x+dx,y,z) * val_r(1) + &
               V(x,y-1,z) * val_r(2) + &
               V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end if
    end do
!$OMP end do nowait

    ! take the y plane at (x,z) = (1|grid%n(1),grid%n(3))

    z = grid%n(3)
!$OMP do
    do y = 2 , grid%n(2) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+dx,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y-1,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y+1,z) * grid%a(2)
          val_r(4) = val_rho(grid,x,y,z-1) * grid%a(3)
          val_r = val_r / sum(val_r)

          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x+dx,y,z) * val_r(1) + &
               V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) + &
               V(x,y,z-1) * val_r(4) )
          
          err = max(err,(vv-V(x,y,z))**2)
          
          V(x,y,z) = vv
       end if
    end do
!$OMP end do nowait

   ! Update all grid%err, atomically
!$OMP atomic
    grid%err = max(grid%err,err)
!$OMP end atomic

!$OMP end parallel

  end subroutine gs_xb

  subroutine gs_yb(grid,y)
    ! this routine calculates the contribution on the
    ! lower/upper y-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: y
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2), val_r(4), err, vv
    integer :: x,dy,z

    V => grid%V
    if ( y == 1 ) then
       dy = 1 
    else
       dy = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

!$OMP parallel default(shared), private(x,z,val_r,vv,err), &
!$OMP firstprivate(sor,y)

    ! take the x plane at (y,z) = (1|grid%n(2),1)
    err = 0._grid_p

    z = 1
!$OMP do
    do x = 2 , grid%n(1) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x+1,y,z) * grid%a(1)
          val_r(3) = val_rho(grid,x,y+dy,z) * grid%a(2)
          val_r(4) = val_rho(grid,x,y,z+1) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) + &
               V(x,y+dy,z) * val_r(3) + &
               V(x,y,z+1) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
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
          val_r(1) = val_rho(grid,x+1,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y+dy,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y,z-1) * grid%a(3)
          val_r(4) = val_rho(grid,x,y,z+1) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x+1,y,z) * val_r(1) + &
               V(x,y+dy,z) * val_r(2) + &
               V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end if

       do x = 2 , grid%n(1) - 1
          vv = sor(1) * V(x,y,z) + sor(2) * val_yb(grid,V,x,y,z,dy)
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end do

       x = grid%n(1)
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y+dy,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y,z-1) * grid%a(3)
          val_r(4) = val_rho(grid,x,y,z+1) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x-1,y,z) * val_r(1)+ &
               V(x,y+dy,z) * val_r(2) + &
               V(x,y,z-1) * val_r(3) + V(x,y,z+1) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end if
    end do
!$OMP end do nowait

    ! take the x plane at (y,z) = (1|grid%n(2),grid%n(3))

    z = grid%n(3)
!$OMP do
    do x = 2 , grid%n(1) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x+1,y,z) * grid%a(1)
          val_r(3) = val_rho(grid,x,y+dy,z) * grid%a(2)
          val_r(4) = val_rho(grid,x,y,z-1) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) + &
               V(x,y+dy,z) * val_r(3) + &
               V(x,y,z-1) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end if
    end do
!$OMP end do nowait

!$OMP atomic
    grid%err = max(grid%err,err)
!$OMP end atomic

!$OMP end parallel

  end subroutine gs_yb

  subroutine gs_zb(grid,z)
    ! this routine calculates the contribution on the
    ! lower/upper z-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: z
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: sor(2), val_r(4), err, vv
    integer :: x,y,dz

    V => grid%V
    if ( z == 1 ) then
       dz = 1 
    else
       dz = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

!$OMP parallel default(shared), private(x,y,val_r,vv,err), &
!$OMP firstprivate(sor,z)

    err = 0._grid_p

!$OMP do
    do y = 2 , grid%n(2) - 1
       x = 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x+1,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y-1,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y+1,z) * grid%a(2)
          val_r(4) = val_rho(grid,x,y,z+dz) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x+1,y,z) * val_r(1) + &
               V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) + &
               V(x,y,z+dz) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end if

       do x = 2 , grid%n(1) - 1
          vv = sor(1) * V(x,y,z) + sor(2) * val_zb(grid,V,x,y,z,dz)
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end do

       x = grid%n(1)
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
          val_r(2) = val_rho(grid,x,y-1,z) * grid%a(2)
          val_r(3) = val_rho(grid,x,y+1,z) * grid%a(2)
          val_r(4) = val_rho(grid,x,y,z+dz) * grid%a(3)
          val_r = val_r / sum(val_r)
          
          vv = sor(1) * V(x,y,z) + sor(2) * (   &
               V(x-1,y,z) * val_r(1) + &
               V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) + &
               V(x,y,z+dz) * val_r(4) )
          err = max(err,(vv-V(x,y,z))**2)
          V(x,y,z) = vv
       end if

    end do
!$OMP end do nowait

!$OMP atomic
    grid%err = max(grid%err,err)
!$OMP end atomic

!$OMP end parallel

  end subroutine gs_zb

  subroutine gs_corner(grid,V,sor,x,y,z,dx,dy,dz)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), intent(inout) :: V(0:,0:,0:)
    real(grid_p), intent(in) :: sor(2)
    integer, intent(in) :: x,y,z,dx,dy,dz
    real(grid_p) :: vcur, val_r(3)
    if ( is_constant(grid,x,y,z) ) return
    val_r(1) = val_rho(grid,x+dx,y,z) * grid%a(1)
    val_r(2) = val_rho(grid,x,y+dy,z) * grid%a(2)
    val_r(3) = val_rho(grid,x,y,z+dz) * grid%a(3)
    val_r = val_r / sum(val_r)
    vcur = &
         V(x+dx,y,z) * val_r(1) + &
         V(x,y+dy,z) * val_r(2) + &
         V(x,y,z+dz) * val_r(3) 
    val_r(1) = sor(1) * V(x,y,z) + sor(2) * vcur
    ! Currently there could be race-conditions for this
    ! value, hence we skip it
    !grid%err = max(grid%err,(val_r(1)-V(x,y,z))**2)
    V(x,y,z) = val_r(1)
  end subroutine gs_corner

  pure function val(grid,V,x,y,z)
    type(mg_grid), intent(in) :: grid
    ! We need to conform to the grid size padding, start from 0
    real(grid_p), intent(in) :: V(0:,0:,0:)
    integer, intent(in) :: x,y,z
    real(grid_p) :: val, val_r(6)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return
    
    val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
    val_r(2) = val_rho(grid,x+1,y,z) * grid%a(1)
    val_r(3) = val_rho(grid,x,y-1,z) * grid%a(2)
    val_r(4) = val_rho(grid,x,y+1,z) * grid%a(2)
    val_r(5) = val_rho(grid,x,y,z-1) * grid%a(3)
    val_r(6) = val_rho(grid,x,y,z+1) * grid%a(3)
    val_r = val_r / sum(val_r)

    val = &
         V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) + &
         V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) + &
         V(x,y,z-1) * val_r(5) + V(x,y,z+1) * val_r(6) 

  end function val

  pure function val_xb(grid,V,x,y,z,dx) result(val)
    ! calculates the contribution on the lower/upper x-bound
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: V(0:,0:,0:)
    integer, intent(in) :: x,y,z,dx
    real(grid_p) :: val, val_r(5)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return

    val_r(1) = val_rho(grid,x+dx,y,z) * grid%a(1)
    val_r(2) = val_rho(grid,x,y-1,z) * grid%a(2)
    val_r(3) = val_rho(grid,x,y+1,z) * grid%a(2)
    val_r(4) = val_rho(grid,x,y,z-1) * grid%a(3)
    val_r(5) = val_rho(grid,x,y,z+1) * grid%a(3)
    val_r = val_r / sum(val_r)

    val = &
         V(x+dx,y,z) * val_r(1) + &
         V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) + &
         V(x,y,z-1) * val_r(4) + V(x,y,z+1) * val_r(5)

  end function val_xb

  pure function val_yb(grid,V,x,y,z,dy) result(val)
    type(mg_grid), intent(in) :: grid
    ! calculates the contribution on the lower/upper x-bound
    real(grid_p), intent(in) :: V(0:,0:,0:)
    integer, intent(in) :: x,y,z,dy
    real(grid_p) :: val, val_r(5)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return

    val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
    val_r(2) = val_rho(grid,x+1,y,z) * grid%a(1)
    val_r(3) = val_rho(grid,x,y+dy,z) * grid%a(2)
    val_r(4) = val_rho(grid,x,y,z-1) * grid%a(3)
    val_r(5) = val_rho(grid,x,y,z+1) * grid%a(3)
    val_r = val_r / sum(val_r)

    val = &
         V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) + &
         V(x,y+dy,z) * val_r(3) + &
         V(x,y,z-1) * val_r(4) + V(x,y,z+1) * val_r(5)
    
  end function val_yb

  pure function val_zb(grid,V,x,y,z,dz) result(val)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: V(0:,0:,0:)
    integer, intent(in) :: x,y,z,dz
    real(grid_p) :: val, val_r(5)

    ! default value must already be initialized
    val = V(x,y,z)
    if ( is_constant(grid,x,y,z) ) return

    val_r(1) = val_rho(grid,x-1,y,z) * grid%a(1)
    val_r(2) = val_rho(grid,x+1,y,z) * grid%a(1)
    val_r(3) = val_rho(grid,x,y-1,z) * grid%a(2)
    val_r(4) = val_rho(grid,x,y+1,z) * grid%a(2)
    val_r(5) = val_rho(grid,x,y,z+dz) * grid%a(3)
    val_r = val_r / sum(val_r)

    val = &
         V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) + &
         V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) + &
         V(x,y,z+dz) * val_r(5)

  end function val_zb

end module m_gs_CDS
