program topbottom

  use t_mg
  use m_gs_CDS

  ! test libraries
  use m_cube
  use m_mg_info

  implicit none

  integer :: method = MG_METHOD_GS_TEMPORAL_CDS

  ! the top multi-grid type
  type(mg_grid), target :: top

  integer :: N , nn(3)
  real :: t1, t2
  real(dp) :: cell(3,3), ll(3)
  real(grid_p) :: tol, sor

  ! tolerance for the convergence
  tol = 1.e-9_grid_p
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/12._dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,12._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,30._dp/)

  ! we create it to be 100x100x100
  nn = 200
  ! create grid
  call init_grid(top,nn,cell,1,2,tol=tol)

  write(*,*)'>> Created initial grid...'

  ! create all children
  call init_grid_children_half(top,max_layer=4)

  ! manually set the sor-parameter
  call grid_set(top,layer=1,sor=1.2_grid_p)
  call grid_set(top,layer=2,sor=1.6_grid_p)
  call grid_set(top,layer=3,sor=1.6_grid_p)
  call grid_set(top,layer=4,sor=1.5_grid_p)
  call grid_set(top,layer=5,sor=1.4_grid_p)

  do N = 1 , layers(top)
     call grid_set(top,layer=N,sor=sor,tol=tol*10**N)
     !if ( N > 1 ) call grid_onoff_layer(top,.false.,layer=N)
  end do
  call grid_set(top,layer=1,tol=tol)

  write(*,*)' >> Created all children...'

  write(*,*)'  >> Add all the boxes...'
  cell(:,3) = cell(:,3) / 10._grid_p
  ll = 0._dp
  call grid_add_box(top, ll, cell, 1._grid_p, 1._grid_p, .true.)
  ll(3) = 3._dp - cell(3,3)
  call grid_add_box(top, ll, cell, -1._grid_p, 1._grid_p, .true.)

  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)
  top%V = 0._grid_p
  call grid_setup(top)

  ! write out the initial cube file
  call write_cube('initial',top)

  CALL CPU_TIME(t1)

  call mg_gs_cds(top)

  CALL CPU_TIME(t2)

  print *,'Timing:', t2-t1

  call write_cube('test',top)
  
  call delete_grid(top)
  
end program topbottom
