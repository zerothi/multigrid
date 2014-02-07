module m_gs_CDS

  use t_mg
  
  implicit none

  private

  integer, parameter :: MG_METHOD_GS_TEMPORAL_CDS = 1

  public :: mg_gs_cds
  public :: MG_METHOD_GS_TEMPORAL_CDS

contains

  subroutine mg_gs_cds(top_grid)
    type(mg_grid), intent(inout), target :: top_grid
    type(mg_grid), pointer :: grid
    integer :: old_itt, n(3)
    real(grid_p) :: tol, old_sum, new_sum

    grid => top_grid

!    ! first we move all the way down to the second lowest...
!    do while ( associated(grid%child%child) )
!       call grid_bring_back(grid%child)
!       call grid_restriction(grid)
!       call grid_hold_back(grid)
!       grid => grid%child
!    end do


!    do 
!
!       !grid => grid_layer(top_grid,layer=-2)
!       n = grid%n / 2
!       print *, layers(grid,enabled=.true.)
!       tol = grid_tolerance(grid) + 1._grid_p
!       old_itt = grid%itt
!       old_sum = grid_sum(grid)! * nr
!       do while ( tol > grid_tolerance(grid) ) 
!
!          ! initialize the tolerance and step iteration
!          grid%itt = grid%itt + 1
!
!          call gs_v_c(grid)
!
!          ! the tolerance is the difference between the 
!          ! sum of the before iteration and the new iteration
!          new_sum = grid_sum(grid)! * nr
!
!          tol = abs(old_sum - new_sum)
!          old_sum = new_sum
!          print '(a,3(tr1,e10.3))',' Tol,new_sum',tol,new_sum,grid%V(n(1),n(2),n(3))
!
!          print '(100(tr1,f4.2))',grid%V(n(1),n(2),1:n(3))
!          print '(100(tr1,f4.2))',grid%V(n(1),n(2),n(3)+1:)
!          ! print out number of iterations used on that cycle
!          !write(*,'(2(a,i0),a)') &
!          !     'Completed (',grid%layer,') cycle in ',grid%itt-old_itt,' cycles'
!
!          if ( grid%itt > old_itt + 100 ) then
!             tol = 0.
!          end if
!
!       end do
!
!       if ( layers(grid,enabled=.true.) == 1 ) exit
!       
!       call grid_onoff_layer(grid,.false.,layer=layers(grid,enabled=.true.))
!       
!    end do
!
!    return

    ! < here is the per-grid solver>

    do while ( associated(grid%child) ) 

       ! Do one step
!       call gs_step(grid)

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

  end subroutine mg_gs_cds

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

       ! initialize the tolerance and step iteration
       grid%itt = grid%itt + 1

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

  subroutine gs_V(top_grid)
    type(mg_grid), intent(inout), target :: top_grid
    type(mg_grid), pointer :: grid
    integer :: i,n(3), cur_layer
    
    grid => top_grid
    do while ( associated(grid%child) ) 

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

    n = grid%n / 2
    print '(''a'',i0,100(tr1,f4.2))',grid%layer,grid%V(n(1),n(2),1:n(3))
    print '(i0,100(tr1,f4.2))',grid%layer,grid%V(n(1),n(2),n(3)+1:)

    ! solve the lowest grid...
    call grid_solve(grid)

    do while ( associated(grid) )

       do i = 1 , grid%steps
          call gs_step(grid)
       end do
       n = grid%n / 2
       print '(i0,100(tr1,f4.2))',grid%layer,grid%V(n(1),n(2),1:n(3))
       print '(i0,100(tr1,f4.2))',grid%layer,grid%V(n(1),n(2),n(3)+1:)

       ! bring back data
       if ( associated(grid%parent) ) then
          call grid_bring_back(grid%parent)
          call grid_prolongation(grid)
          call grid_hold_back(grid)
       end if

       grid => grid%parent

    end do

  end subroutine gs_V

  subroutine gs_V_c(top_grid)
    type(mg_grid), intent(inout), target :: top_grid
    type(mg_grid), pointer :: grid
    integer :: i,n(3), cur_layer

    call gs_down(top_grid,grid)

    n = grid%n / 2
    print '(''a'',i0,100(tr1,f4.2))',grid%layer,grid%V(n(1),n(2),1:n(3))
    print '(i0,100(tr1,f4.2))',grid%layer,grid%V(n(1),n(2),n(3)+1:)

    ! solve the lowest grid...
    call grid_solve(grid)

    call gs_up(grid,top_grid)

  end subroutine gs_V_c

  subroutine gs_down(top_grid,bottom_grid)
    type(mg_grid), intent(inout), target :: top_grid
    type(mg_grid), pointer :: bottom_grid
    type(mg_grid), pointer :: grid
    integer :: i
    
    grid => top_grid
    bottom_grid => grid
    do while ( associated(grid%child) ) 

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

       ! update bottom-grid
       bottom_grid => grid

    end do

  end subroutine gs_down

  subroutine gs_up(bottom_grid,top_grid)
    type(mg_grid), intent(inout), target :: bottom_grid
    type(mg_grid), intent(in)            :: top_grid
    type(mg_grid), pointer               :: grid
    integer :: i,n(3)
    
    grid => bottom_grid
    do while ( associated(grid) )
       n = grid%n

       do i = 1 , grid%steps
          call gs_step(grid)
       end do
       n = grid%n
       print '(i0,100(tr1,f4.2))',grid%layer,grid%V(n(1)/2,n(2)/2,1:n(3)/2)
       print '(i0,100(tr1,f4.2))',grid%layer,grid%V(n(1)/2,n(2)/2,n(3)/2+1:)

       ! bring back data
       if ( associated(grid%parent) ) then
          if ( grid%layer == top_grid%layer ) return
          call grid_bring_back(grid%parent)
          call grid_prolongation(grid)
          call grid_hold_back(grid)
       end if

       grid => grid%parent

    end do

  end subroutine gs_up

  subroutine gs_step(grid)
    type(mg_grid), intent(inout) :: grid

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
    
    ! x-corners
!    call gs_corner(a,sor,V,x,      1,1, dx, 1,1,tol)

!$OMP parallel default(shared) &
!$OMP   private(y,z,val_r) firstprivate(sor,x)

! consider adding single, and add dynamic scheduling for the "big" loop

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

    ! x-corners
!    call gs_corner(a,sor,V,x,      1,grid%n(3), dx, 1,-1,tol)
!    call gs_corner(a,sor,V,x,grid%n(2),grid%n(3), dx,-1,-1,tol)

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

!    ! y-corners
!    call gs_corner(grid,sor,V,      1,y,1,  1,dy,1,tol)
!    call gs_corner(grid,sor,V,grid%n(1),y,1, -1,dy,1,tol)

!$OMP parallel default(shared) &
!$OMP   private(x,z,val_r) firstprivate(sor,y)

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

!    ! y-corners
!    call gs_corner(grid,sor,V,      1,y,grid%n(3),  1,dy,-1,tol)
!    call gs_corner(grid,sor,V,grid%n(1),y,grid%n(3), -1,dy,-1,tol)

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

!    ! z-corners
!    call gs_corner(grid,sor,V,      1,1,z,  1,1,dz,tol)
!    call gs_corner(grid,sor,V,grid%n(1),1,z, -1,1,dz,tol)
    
!$OMP parallel default(shared) &
!$OMP   private(x,y,val_r) firstprivate(sor,z)

    y = 1
!$OMP do 
    do x = 2 , grid%n(1) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z)
          val_r(2) = val_rho(grid,x+1,y,z)
          val_r(3) = val_rho(grid,x,y+1,z)
          val_r(4) = val_rho(grid,x,y,z+dz)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
               grid%a(2) * ( V(x,y+1,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z+dz) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

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

    y = grid%n(2)
!$OMP do
    do x = 2 , grid%n(1) - 1
       if ( .not. is_constant(grid,x,y,z) ) then
          val_r(1) = val_rho(grid,x-1,y,z)
          val_r(2) = val_rho(grid,x+1,y,z)
          val_r(3) = val_rho(grid,x,y-1,z)
          val_r(4) = val_rho(grid,x,y,z+dz)
          val_r = val_r / sum(val_r)
          
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * (   &
               grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
               grid%a(2) * ( V(x,y-1,z) * val_r(3) ) + &
               grid%a(3) * ( V(x,y,z+dz) * val_r(4) ) )
       end if
    end do
!$OMP end do nowait

!$OMP end parallel

!    ! z-corners
!    call gs_corner(grid,sor,V,      1,grid%n(2),z,  1,-1,dz,tol)
!    call gs_corner(grid,sor,V,grid%n(1),grid%n(2),z, -1,-1,dz,tol)

  end subroutine gs_zb

  subroutine gs_corner(grid,sor,x,y,z,dx,dy,dz)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: sor(2)
    integer, intent(in) :: x,y,z,dx,dy,dz
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: vcur, val_r(3)
    if ( is_constant(grid,x,y,z) ) return
    V => grid%V
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
