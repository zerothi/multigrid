program test

  use t_mg
  use m_gs_CDS
  use m_mg_save

  ! test libraries
  use m_mg_info
  use m_time

  implicit none

  ! the top multi-grid type
  type(mg_grid), target :: top

  integer :: N , nn(3), c1
  real :: time
  real(dp) :: cell(3,3), bcell(3,3), ll(3)
  real(dp) :: tol, sor

  call init_timing()

  ! tolerance for the convergence
  tol = 1.e-4_grid_p
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/12._dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,12._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,5._dp/)

  ! we create it to be 100x100x100
  nn = 200
  ! create grid
  call init_grid(top,nn,cell,2,tol=tol)

  write(*,*)'>> Created initial grid...'

  ! create all children
  call init_grid_children_half(top,max_layer=4)

  ! manually set the sor-parameter
  call grid_set(top,layer=1,sor=1.2_dp)
  call grid_set(top,layer=2,sor=1.6_dp)
  call grid_set(top,layer=3,sor=1.6_dp)
  call grid_set(top,layer=4,sor=1.5_dp)
  call grid_set(top,layer=5,sor=1.4_dp)

  do N = 1 , layers(top)
     call grid_set(top,layer=N,sor=sor,tol=tol)
  end do

  write(*,*)' >> Created all children...'

  write(*,*)'  >> Add all the boxes...'
  bcell(:,:) = cell(:,:)
  bcell(:,3) = cell(:,3) / 10._dp
  ll = 0._dp
  call grid_add_box(top, ll, bcell, 1._dp, 1._dp, .true.)
  ll(3) = cell(3,3) - bcell(3,3)
  call grid_add_box(top, ll, bcell, -1._dp, 1._dp, .true.)

  call print_grid(top)

  call grid_bring_back(top)
  call grid_setup(top, init = .true. ) ! only for saving initial

  ! write out the initial cube file (the initial one we only save
  ! in cube)
  call mg_save(top,'initial',MG_SAVE_CUBE)

  call grid_hold_back(top)

  c1 = clock()

  call mg_gs_cds(top,method=CDS_BOTTOM_UP)

  time = timing(c1)

  print *,'Timing - BU:', time

  ! write out the initial cube file
  call grid_save_all(top,'test2_BU')

  call grid_hold_back(top)

  c1 = clock()

  call mg_gs_cds(top,method=CDS_W_CYCLE)

  time = timing(c1)

  print *,'Timing - W: ', time

  ! write out the initial cube file
  call grid_save_all(top,'test2_W')
  
  call delete_grid(top)
  
end program test
