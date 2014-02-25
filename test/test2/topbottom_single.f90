program topbottom_single

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
  real(dp) :: cell(3,3), bcell(3,3), ll(3)
  real(grid_p) :: tol, sor

  call init_timing()

  ! tolerance for the convergence
  tol = 1.e-3_grid_p
  ! over-relaxation parameter
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/12._dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,12._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,30._dp/)

  ! we create it to be 100x100x100
  nn = 200
  ! create grid
  call init_grid(top,nn,cell,2,tol=tol,sor=sor)

  write(*,*)'>> Created initial grid...'

  write(*,*)'  >> Add all the boxes...'
  bcell(:,:) = cell(:,:)
  bcell(:,3) = cell(:,3) / 10._grid_p
  ll = 0._dp
  call grid_add_box(top, ll, bcell, 1._grid_p, 1._grid_p, .true.)
  ll(3) = cell(3,3) - bcell(3,3)
  call grid_add_box(top, ll, bcell, -1._grid_p, 1._grid_p, .true.)

  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)
  call grid_setup(top)

  ! write out the initial cube file
  call mg_save(top,'initial_single',MG_SAVE_CUBE)

  c1 = clock()

  call mg_gs_cds(top)

  time = timing(c1)

  print *,'Timing:', time

  call mg_save(top,'test_single',MG_SAVE_CUBE)
  
  call delete_grid(top)
  
end program topbottom_single
