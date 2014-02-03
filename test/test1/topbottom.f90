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
  real(dp) :: cell(3,3), ll(3)
  real(grid_p) :: tol, sor

  ! tolerance for the convergence
  tol = 1.e-6_grid_p
  ! over-relaxation parameter
  sor = 1.25_grid_p

  ! initialize the initial grid
  cell(:,1) = (/2._dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,2._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,2._dp/)

  ! we create it to be 100x100x100
  nn = 100
  ! create grid
  call init_grid(top,nn,cell,1,2,tol=tol,sor=sor)

  write(*,*)'>> Created initial grid...'

  ! create all children
  call init_grid_children_half(top)

  write(*,*)' >> Created all children...'

  write(*,*)'  >> Add all the boxes...'
  cell(:,3) = cell(:,3) / 10._grid_p
  ll = 0._dp
  call grid_add_box(top, ll, cell, 1._grid_p, .true.)
  ll(3) = 2._dp - cell(3,3)
  call grid_add_box(top, ll, cell, -1._grid_p, .true.)

  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)
  top%V = 0._grid_p
  call grid_setup(top)
  print*,maxval(top%V)

  ! write out the initial cube file
  call write_cube('initial',top)

  call mg_gs_cds(top)

  call write_cube('test',top)
  
  call delete_grid(top)
  
end program topbottom
