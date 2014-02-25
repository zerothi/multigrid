program topbottom

  use t_mg
  use m_gs_CDS
  use m_mg_save

  ! test libraries
  use m_mg_info
  use m_time

  implicit none

  integer :: method = MG_METHOD_GS_TEMPORAL_CDS

  ! the top multi-grid type
  type(mg_grid), target :: top

  integer :: N , nn(3), c1
  real :: time
  real(dp) :: cell(3,3), ll(3), bcell(3,3)
  real(grid_p) :: tol, sor

  call init_timing()

  ! tolerance for the convergence
  tol = 1.e-4_grid_p
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/30._dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,30._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,30._dp/)

  ! we create it to be 100x100x100
  nn = 200
  ! create grid
  call init_grid(top,nn,cell,7,tol=tol)

  write(*,*)'>> Created initial grid...'

  ! create all children
  call init_grid_children_half(top,max_layer=4)

  ! manually set the sor-parameter
  call grid_set(top,layer=1,sor=1.8_grid_p)
  call grid_set(top,layer=2,sor=1.8_grid_p)
  call grid_set(top,layer=3,sor=1.8_grid_p)
  call grid_set(top,layer=4,sor=1.8_grid_p)
  call grid_set(top,layer=5,sor=1.8_grid_p)

  do N = 1 , layers(top)
     call grid_set(top,layer=N,sor=sor,tol=tol)
     !if ( N > 1 ) call grid_onoff_layer(top,.false.,layer=N)
  end do
  call grid_set(top,layer=1,tol=tol)

  write(*,*)' >> Created all children...'

  write(*,*)'  >> Add all the boxes...'
  ! The "electrodes"
  ! the first one is in the middle of the bottom x-cell
  bcell(:,:) = cell(:,:) / 10._dp
  ll = (/ &
       cell(1,1) / 2._dp - bcell(1,1) / 2._dp , &
       0._dp , &
       0._dp /)
  call grid_add_box(top, ll, bcell, 0.5_grid_p, 1._grid_p, .true.)
  ll(2) = cell(2,2) - bcell(2,2)
  ll(3) = cell(3,3) / 10._dp
  call grid_add_box(top, ll, bcell, 1.0_grid_p, 1._grid_p, .true.)
  ll(2) = cell(2,2) / 2._dp - bcell(2,2) / 2._dp
  ll(3) = cell(3,3) - bcell(3,3)
  call grid_add_box(top, ll, bcell, -1._grid_p, 1._grid_p, .true.)

  ! add points to control run-away potentials
  do N = 1 , 3
     bcell(:,N) = cell(:,N) / (nn(N)-.05 * nn(N))
  end do
  ll = (/ 0._dp,0._dp,cell(3,3)/)
  call grid_add_box(top, ll, bcell, 0._grid_p, 1._grid_p, .true.)
  ll = (/ cell(1,1),cell(2,2),cell(3,3)/)
  call grid_add_box(top, ll, bcell, 0._grid_p, 1._grid_p, .true.)
  ll = (/ cell(1,1),cell(2,2),0._dp/)
  call grid_add_box(top, ll, bcell, 0._grid_p, 1._grid_p, .true.)
  ll = (/ 0._dp,cell(2,2),0._dp/)
  call grid_add_box(top, ll, bcell, 0._grid_p, 1._grid_p, .true.)
  
  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)
  call grid_setup(top)

  ! write out the initial cube file
  !call mg_save(top,'initial',MG_SAVE_CUBE)

  c1 = clock()

  call mg_gs_cds(top,CDS_W_CYCLE)

  time = timing(c1)

  print *,'Timing:', time

  c1 = clock()
  call mg_save(top,'test',MG_SAVE_CUBE)
  time = timing(c1)
  print *,'Timing CUBE:', time
  c1 = clock()
  call mg_save(top,'test',MG_SAVE_BINARY)
  time = timing(c1)
  print *,'Timing BINARY:', time

  call delete_grid(top)
  
end program topbottom
