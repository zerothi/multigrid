module m_gs

  use t_mg
  
  implicit none

  private

  integer, parameter :: MG_METHOD_GS_TEMPORAL = 1

  public :: mg_gs
  public :: MG_METHOD_GS_TEMPORAL

contains

  subroutine mg_gs(top_grid)
    type(mg_grid), pointer :: top_grid
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

       call grid_solve(grid)

       ! bring back data
       call grid_bring_back(grid%parent)

       call grid_prolongation(grid)

       ! hold back data
       if ( associated(grid%parent) ) then
          call grid_hold_back(grid)
       end if

       write(*,'(2(a,i0),a)') &
            'Completed (',grid%layer,') cycle in ',grid%itt-old_itt,' cycles'
       
       grid => grid%parent

!       cycle
!
!       ! we decide whether we should move down or up in the V-cycle
!       if ( grid%itt - old_itt > 10 .and. &
!            associated(grid%child) ) then ! move down in V-cycle
!          ! precontraction 
!          call init_grid_child(grid%child)
!          grid => grid%child
!       else if ( grid%itt - old_itt > 10 ) then
!          ! stay-put
!       else 
!          ! move up in V-shape
!          ! prolongation
!          call init_grid_parent(grid)
!          grid => grid%parent
!       end if

    end do

  end subroutine mg_gs

  subroutine grid_solve(grid)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), pointer :: V(:,:,:)
    real(grid_p) :: tol

    tol = grid%tol + 1._grid_p

    call from1Dto3D(grid%n1,grid%n2,grid%n3,grid%V, V)

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
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, a(3), sor(2)
    integer :: x,y,z

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

    a(1) = grid%ax
    a(2) = grid%ay
    a(3) = grid%az
    
    do z = 2 , grid%n3 - 1
       do y = 2 , grid%n2 - 1
          do x = 2 , grid%n1 - 1
             if ( in_constant(x,y,z) ) cycle
             ! calculate the current contribution
             vcur = val(a,V,x,y,z)
             ! Calculate the tolerance
             tol = max(abs(V(x,y,z) - vcur),tol)
             V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
          end do
       end do
    end do
    
  end subroutine gs

  subroutine gs_bound(grid,V,tol)
    type(mg_grid), intent(inout) :: grid
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(inout) :: tol

    call gs_xb(grid,1,V,tol)
    call gs_xb(grid,grid%n1,V,tol)
    call gs_yb(grid,1,V,tol)
    call gs_yb(grid,grid%n2,V,tol)
    call gs_zb(grid,1,V,tol)
    call gs_zb(grid,grid%n3,V,tol)

  end subroutine gs_bound
  
  subroutine gs_xb(grid,x,V,tol)
    ! this routine calculates the contribution on the
    ! lower/upper x-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: x
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, a(3), sor(2)
    integer :: dx,y,z
    if ( x == 1 ) then
       dx = 1 
    else
       dx = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

    a(1) = grid%ax
    a(2) = grid%ay
    a(3) = grid%az
    
    ! x-corners
!    call gs_corner(a,sor,V,x,      1,1, dx, 1,1,tol)
!    z = 1
!    do y = 2 , grid%n2 - 1
!       if ( is_constant(x,y,z) ) cycle
!       call gs_line(,x,grid%n2,1, dx,-1,1,tol)
!    end do
!    call gs_corner(a,sor,V,x,grid%n2,1, dx,-1,1,tol)

    do z = 2 , grid%n3 - 1
       do y = 2 , grid%n2 - 1
          if ( in_constant(x,y,z) ) cycle
          ! calculate the current contribution
          vcur = val_xb(a,V,dx,y,z)
          ! Calculate the tolerance
          tol = max(abs(V(x,y,z) - vcur),tol)
          ! we implement the SOR-algorithm
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
       end do
    end do

    ! x-corners
!    call gs_corner(a,sor,V,x,      1,grid%n3, dx, 1,-1,tol)
!    call gs_corner(a,sor,V,x,grid%n2,grid%n3, dx,-1,-1,tol)

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
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, a(3), sor(2)
    integer :: x,dy,z
    if ( y == 1 ) then
       dy = 1 
    else
       dy = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

    a(1) = grid%ax
    a(2) = grid%ay
    a(3) = grid%az

!    ! y-corners
!    call gs_corner(a,sor,V,      1,y,1,  1,dy,1,tol)
!    call gs_corner(a,sor,V,grid%n1,y,1, -1,dy,1,tol)
    
    do z = 2 , grid%n3 - 1
       do x = 2 , grid%n1 - 1
          if ( in_constant(x,y,z) ) cycle
          ! calculate the current contribution
          vcur = val_yb(a,V,x,dy,z)
          ! Calculate the tolerance
          tol = max(abs(V(x,y,z) - vcur),tol)
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
       end do
    end do

!    ! y-corners
!    call gs_corner(a,sor,V,      1,y,grid%n3,  1,dy,-1,tol)
!    call gs_corner(a,sor,V,grid%n1,y,grid%n3, -1,dy,-1,tol)

  end subroutine gs_yb

  subroutine gs_zb(grid,z,V,tol)
    ! this routine calculates the contribution on the
    ! lower/upper z-bound
    type(mg_grid), intent(inout) :: grid
    integer, intent(in) :: z
    real(grid_p), intent(inout) :: V(grid%n1,grid%n2,grid%n3)
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, a(3), val_r(3), sor(2)
    integer :: x,y,dz
    if ( z == 1 ) then
       dz = 1 
    else
       dz = -1
    end if

    sor(2) = grid%sor
    sor(1) = 1._grid_p - sor(2)

    a(1) = grid%ax
    a(2) = grid%ay
    a(3) = grid%az

!    ! z-corners
!    call gs_corner(a,sor,V,      1,1,z,  1,1,dz,tol)
!    call gs_corner(a,sor,V,grid%n1,1,z, -1,1,dz,tol)
    
    do y = 2 , grid%n2 - 1
       do x = 2 , grid%n1 - 1
          if ( in_constant(x,y,z) ) cycle
          ! calculate the current contribution
          vcur = val_zb(a,V,x,y,dz)
          ! Calculate the tolerance
          tol = max(abs(V(x,y,z) - vcur),tol)
          V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
       end do
    end do

!    ! z-corners
!    call gs_corner(a,sor,V,      1,grid%n2,z,  1,-1,dz,tol)
!    call gs_corner(a,sor,V,grid%n1,grid%n2,z, -1,-1,dz,tol)

  end subroutine gs_zb

  subroutine gs_corner(a,sor,V,x,y,z,dx,dy,dz,tol)
    real(grid_p), intent(in) :: a(3), sor(2)
    real(grid_p), intent(inout) :: V(:,:,:)
    integer, intent(in) :: x,y,z,dx,dy,dz
    real(grid_p), intent(inout) :: tol
    real(grid_p) :: vcur, val_r(3)
    if ( is_constant(x,y,z) ) return
    val_r(1) = val_rho(x+dx,y,z)
    val_r(2) = val_rho(x,y+dy,z)
    val_r(3) = val_rho(x,y,z+dz)
    val_r = val_r / sum(val_r)
    vcur = &
         a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
         a(2) * ( V(x,y+dy,z) * val_r(2) ) + &
         a(3) * ( V(x,y,z+dz) * val_r(3) )
    tol = max(abs(V(x,y,z) - vcur),tol)
    V(x,y,z) = sor(1) * V(x,y,z) + sor(2) * vcur
  end subroutine gs_corner

  pure function val(a,V,x,y,z) 
    real(grid_p), intent(in) :: a(3), V(:,:,:)
    integer, intent(in) :: x,y,z
    real(grid_p) :: val, val_r(6)

    val_r(1) = val_rho(x-1,y,z)
    val_r(2) = val_rho(x+1,y,z)
    val_r(3) = val_rho(x,y-1,z)
    val_r(4) = val_rho(x,y+1,z)
    val_r(5) = val_rho(x,y,z-1)
    val_r(6) = val_rho(x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         a(2) * ( V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) ) + &
         a(3) * ( V(x,y,z+1) * val_r(5) + V(x,y,z+1) * val_r(6) )

  end function val

  pure function val_xb(a,V,dx,y,z) 
    ! calculates the contribution on the lower/upper x-bound
    real(grid_p), intent(in) :: a(3), V(:,:,:)
    integer, intent(in) :: dx,y,z
    real(grid_p) :: val, val_r(5)

    val_r(1) = val_rho(x+dx,y,z)
    val_r(2) = val_rho(x,y-1,z)
    val_r(3) = val_rho(x,y+1,z)
    val_r(4) = val_rho(x,y,z-1)
    val_r(5) = val_rho(x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         2._grid_p * a(1) * ( V(x+dx,y,z) * val_r(1) ) + &
         a(2) * ( V(x,y-1,z) * val_r(2) + V(x,y+1,z) * val_r(3) ) + &
         a(3) * ( V(x,y,z+1) * val_r(4) + V(x,y,z+1) * val_r(5) )

  end function val_xb

  pure function val_yb(a,V,x,dy,z) 
    ! calculates the contribution on the lower/upper x-bound
    real(grid_p), intent(in) :: a(3), V(:,:,:)
    integer, intent(in) :: x,dy,z
    real(grid_p) :: val, val_r(5)

    val_r(1) = val_rho(x-1,y,z)
    val_r(2) = val_rho(x+1,y,z)
    val_r(3) = val_rho(x,y+dx,z)
    val_r(4) = val_rho(x,y,z-1)
    val_r(5) = val_rho(x,y,z+1)
    val_r = val_r / sum(val_r)

    val = &
         a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         2._grid_p * a(2) * ( V(x,y+dx,z) * val_r(3) ) + &
         a(3) * ( V(x,y,z+1) * val_r(4) + V(x,y,z+1) * val_r(5) )
    
  end function val_yb

  pure function val_zb(a,V,x,y,dz) 
    real(grid_p), intent(in) :: a(3), V(:,:,:)
    integer, intent(in) :: x,y,dz
    real(grid_p) :: val, val_r(5)

    val_r(1) = val_rho(x-1,y,z)
    val_r(2) = val_rho(x+1,y,z)
    val_r(3) = val_rho(x,y-1,z)
    val_r(4) = val_rho(x,y+1,z)
    val_r(5) = val_rho(x,y,z+dz)
    val_r = val_r / sum(val_r)

    val = &
         a(1) * ( V(x-1,y,z) * val_r(1) + V(x+1,y,z) * val_r(2) ) + &
         a(2) * ( V(x,y-1,z) * val_r(3) + V(x,y+1,z) * val_r(4) ) + &
         2._grid_p * a(3) * ( V(x,y,z+dz) * val_r(5) )

  end function val

end module m_gs
