module m_gs_CDS

  use t_mg
  
  implicit none

  private

  integer, parameter :: MG_METHOD_GS_TEMPORAL_CDS = 1

  public :: mg_gs_cds
  public :: MG_METHOD_GS_TEMPORAL_CDS

contains

  subroutine mg_gs_cds(top_grid)
    type(mg_grid), target :: top_grid
    type(mg_grid), pointer :: grid
    integer :: old_itt
    
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

       ! bring back data
       if ( associated(grid%parent) ) then
          call grid_bring_back(grid%parent)
       end if

       call grid_prolongation(grid)

       ! hold back data
       if ( associated(grid%parent) ) then
          call grid_hold_back(grid)
       end if

       ! print out number of iterations used on that cycle
       write(*,'(2(a,i0),a)') &
            'Completed (',grid%layer,') cycle in ',grid%itt-old_itt,' cycles'
       
       grid => grid%parent

    end do

  end subroutine mg_gs_cds

  subroutine grid_solve(grid)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: tol

    tol = grid%tol + 1._grid_p

    V => grid%V

    do while ( tol > grid%tol ) 

       ! initialize the tolerance and step iteration
       tol = 0._dp
       grid%itt = grid%itt + 1
       
       !< communicate bounds >

       call gs(grid, V, tol)
       
       !< wait communicate >
       
       call gs_bound(grid, V, tol)

    end do

  end subroutine grid_solve
  
  subroutine gs(grid,V,tol)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), intent(inout) :: V(grid%n(1),grid%n(2),grid%n(3))
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, a(3), sor(2)
    integer :: x,y,z

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

    do z = 2 , grid%n(3) - 1
       do y = 2 , grid%n(2) - 1
          do x = 2 , grid%n(1) - 1
             ! calculate the current contribution
             vcur = val(grid,V,x,y,z)
             ! Calculate the tolerance
             tol = max(abs(V(x,y,z) - vcur),tol)
             V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
          end do
       end do
    end do
    
  end subroutine gs

  subroutine gs_bound(grid,V,tol)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), intent(inout) :: V(grid%n(1),grid%n(2),grid%n(3))
    real(grid_p), intent(inout) :: tol

    call gs_xb(grid,1,V,tol)
    call gs_xb(grid,grid%n(1),V,tol)
    call gs_yb(grid,1,V,tol)
    call gs_yb(grid,grid%n(2),V,tol)
    call gs_zb(grid,1,V,tol)
    call gs_zb(grid,grid%n(3),V,tol)

  end subroutine gs_bound
  
  subroutine gs_xb(grid,x,V,tol)
    ! this routine calculates the contribution on the
    ! lower/upper x-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: x
    real(grid_p), intent(inout) :: V(grid%n(1),grid%n(2),grid%n(3))
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, sor(2)
    integer :: dx,y,z

    if ( x == 1 ) then
       dx = 1 
    else
       dx = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)
    
    ! x-corners
!    call gs_corner(a,sor,V,x,      1,1, dx, 1,1,tol)
!    z = 1
!    do y = 2 , grid%n(2) - 1
!       if ( is_constant(x,y,z) ) cycle
!       call gs_line(,x,grid%n(2),1, dx,-1,1,tol)
!    end do
!    call gs_corner(a,sor,V,x,grid%n(2),1, dx,-1,1,tol)

    do z = 2 , grid%n(3) - 1
       do y = 2 , grid%n(2) - 1
          ! calculate the current contribution
          vcur = val_xb(grid,V,x,y,z,dx)
          ! Calculate the tolerance
          tol = max(abs(V(x,y,z) - vcur),tol)
          ! we implement the SOR-algorithm
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
       end do
    end do

    ! x-corners
!    call gs_corner(a,sor,V,x,      1,grid%n(3), dx, 1,-1,tol)
!    call gs_corner(a,sor,V,x,grid%n(2),grid%n(3), dx,-1,-1,tol)

!  contains 
!    
!    subroutine gs_line(dy,dz)
!      integer, intent(in) :: dy,dz
!      real(grid_p) :: val_r(4)
!      val_r(1) = val_rho(x-1,y,z)
!      val_r(2) = val_rho(x+1,y,z)
!      val_r(3) = val_rho(x,y+dy,z)
!      val_r(4) = val_rho(x,y,z+dz)
!      val_r = val_r / sum(val_r)
!      vcur = &
!           a(1) * ( V(x-1,y,z)  * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
!           a(2) * ( V(x,y+dy,z) * val_r(3) ) + &
!           a(3) * ( V(x,y,z+dz) * val_r(4) )
!      tol = max(abs(V(x,y,z) - vcur),tol)
!      V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
!    end subroutine gs_line

  end subroutine gs_xb

  subroutine gs_yb(grid,y,V,tol)
    ! this routine calculates the contribution on the
    ! lower/upper y-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: y
    real(grid_p), intent(inout) :: V(grid%n(1),grid%n(2),grid%n(3))
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, sor(2)
    integer :: x,dy,z
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
    
    do z = 2 , grid%n(3) - 1
       do x = 2 , grid%n(1) - 1
          ! calculate the current contribution
          vcur = val_yb(grid,V,x,y,z,dy)
          ! Calculate the tolerance
          tol = max(abs(V(x,y,z) - vcur),tol)
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
       end do
    end do

!    ! y-corners
!    call gs_corner(grid,sor,V,      1,y,grid%n(3),  1,dy,-1,tol)
!    call gs_corner(grid,sor,V,grid%n(1),y,grid%n(3), -1,dy,-1,tol)

  end subroutine gs_yb

  subroutine gs_zb(grid,z,V,tol)
    ! this routine calculates the contribution on the
    ! lower/upper z-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: z
    real(grid_p), intent(inout) :: V(grid%n(1),grid%n(2),grid%n(3))
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, val_r(3), sor(2)
    integer :: x,y,dz

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
    
    do y = 2 , grid%n(2) - 1
       do x = 2 , grid%n(1) - 1
          ! calculate the current contribution
          vcur = val_zb(grid,V,x,y,z,dz)
          ! Calculate the tolerance
          tol = max(abs(V(x,y,z) - vcur),tol)
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
       end do
    end do

!    ! z-corners
!    call gs_corner(grid,sor,V,      1,grid%n(2),z,  1,-1,dz,tol)
!    call gs_corner(grid,sor,V,grid%n(1),grid%n(2),z, -1,-1,dz,tol)

  end subroutine gs_zb

  subroutine gs_corner(grid,sor,V,x,y,z,dx,dy,dz,tol)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: sor(2)
    real(grid_p), intent(inout) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dx,dy,dz
    real(grid_p), intent(inout) :: tol
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
    tol = max(abs(V(x,y,z) - vcur),tol)
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
         grid%a(3) * ( V(x,y,z+1) * val_r(5) + V(x,y,z+1) * val_r(6) )

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
         2._grid_p * grid%a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
         grid%a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
         grid%a(3) * ( V(x,y,z+1) * val_r(4) + V(x,y,z+1) * val_r(5) )

  end function val_xb

  pure function val_yb(grid,V,x,y,z,dy) result(val)
    type(mg_grid), intent(in) :: grid
    ! calculates the contribution on the lower/upper x-bound
    real(grid_p), intent(in) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dy
    real(grid_p) :: val, val_r(5)

    val_r(1) = val_rho(grid,x-1,y,z)
    val_r(2) = val_rho(grid,x+1,y,z)
    val_r(3) = val_rho(grid,x,y+dy,z)
    val_r(4) = val_rho(grid,x,y,z-1)
    val_r(5) = val_rho(grid,x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         2._grid_p * grid%a(2) * ( V(x,y+dy,z) * val_r(3) ) + &
         grid%a(3) * ( V(x,y,z+1) * val_r(4) + V(x,y,z+1) * val_r(5) )
    
  end function val_yb

  pure function val_zb(grid,V,x,y,z,dz) result(val)
    type(mg_grid), intent(in) :: grid
    real(grid_p), intent(in) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dz
    real(grid_p) :: val, val_r(5)

    val_r(1) = val_rho(grid,x-1,y,z)
    val_r(2) = val_rho(grid,x+1,y,z)
    val_r(3) = val_rho(grid,x,y-1,z)
    val_r(4) = val_rho(grid,x,y+1,z)
    val_r(5) = val_rho(grid,x,y,z+dz)
    val_r = val_r / sum(val_r)

    val = &
         grid%a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         grid%a(2) * ( V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) ) + &
         2._grid_p * grid%a(3) * ( V(x,y,z+dz) * val_r(5) )

  end function val_zb

end module m_gs_CDS
