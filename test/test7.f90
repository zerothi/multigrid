program test

  use t_mg
  use t_mg_interp
  use m_gs_CDS
  use m_mg_save

  ! test libraries
  use m_mg_info
  use m_time

  implicit none

  ! the top multi-grid type
  type(mg_grid), target :: top
  type(mg_grid), pointer :: g

  integer :: N , nn(3), c1
  real :: time
  real(dp) :: cell(3,3)
  real(dp) :: tol, sor
  integer :: i, p, c, cc

  call init_timing()

  ! tolerance for the convergence
  tol = 1.e-4_grid_p
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/31.75_dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,20._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,31.75_dp/)

  ! we create it to be 100x100x100
  nn = 200
  ! create grid
  call init_grid(top,nn,cell,4,tol=tol)

  write(*,*) '>> Created initial grid...'

  ! create all children
  call init_grid_children_half(top,max_layer=4)

  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)
  call grid_setup(top,init=.true.)

  ! write out the initial cube file
!  call mg_save(top,'test6_initial',MG_SAVE_CUBE)

  call grid_hold_back(top)

  c1 = clock()

  g => top
  do while ( associated(g%child) )

     ! Track parent
     do i = 1 , 3
        N = 0
        c = -1
        do p = 1 , g%n(i)
           cc = g2g(g%n(i),p,g%child%n(i))
           if ( cc > g%child%n(i) ) then
              print *,'TOO BIG',i
           end if
           if ( cc /= c ) then
              N = N + 1
              c = cc
           end if
        end do
        if ( N /= g%child%n(i) ) then
           print *,'Tracking child: ',i,N,g%child%n(i)
        end if
     end do

     g => g%child

     do i = 1 , 3
        N = 0
        c = -1
        do p = 1 , g%n(i)
           cc = g2g(g%n(i),p,g%parent%n(i))
           if ( cc > g%parent%n(i) ) then
              print *,'TOO BIG',i
           end if
           if ( cc /= c ) then
              N = N + 1
              c = cc
           end if
        end do
        if ( N /= g%n(i) ) then
           print *,'Tracking parent: ',i,N,g%n(i)
        end if
     end do

  end do

  time = timing(c1)

  print *,'Timing:', time

  call delete_grid(top)
  
end program test
