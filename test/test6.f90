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
  real(dp) :: cell(3,3), ll(3), bcell(3,3)
  real(dp) :: tol, sor

  call init_timing()

  ! tolerance for the convergence
  tol = 1.e-4_grid_p
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/31.75_dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,20._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,31.75_dp/)

  ! we create it to be 100x100x100
  nn = 240
  nn(2) = 150
  ! create grid
  call init_grid(top,nn,cell,4,tol=tol)

  write(*,*)'>> Created initial grid...'

  ! create all children
  call init_grid_children_half(top,max_layer=4)

  ! manually set the sor-parameter
  call grid_set(top,layer=1,sor=1.8_dp)
  call grid_set(top,layer=2,sor=1.8_dp)
  call grid_set(top,layer=3,sor=1.8_dp)
  call grid_set(top,layer=4,sor=1.8_dp)
  call grid_set(top,layer=5,sor=1.8_dp)

  do N = 1 , layers(top)
     call grid_set(top,layer=N,sor=sor,tol=tol)
     !if ( N > 1 ) call grid_onoff_layer(top,.false.,layer=N)
  end do
  call grid_set(top,layer=layers(top),tol=tol/100._dp)

  write(*,*)' >> Created all children...'

  write(*,*)'  >> Add all the boxes...'
  ! Four boxes, three with same potential, 1 with different.
  call grid_BC(top,MG_BC_DIRICHLET)
  call grid_BC(top,MG_BC_NEUMANN,plane = MG_BC_C0)

  ll         = (/  0.00000E+00,  0.20100E+01,  0.12340E+02/)
  bcell(:,1) = (/  0.63500E+01,  0.00000E+00,  0.00000E+00/)
  bcell(:,2) = (/  0.00000E+00,  0.70000E+01,  0.00000E+00/)
  bcell(:,3) = (/  0.00000E+00,  0.00000E+00,  0.70000E+01/)
  call grid_add_box(top, ll, bcell, .5_dp, 1._dp, .true.)

  ll         = (/  0.25400E+02,  0.20100E+01,  0.12340E+02/)
  call grid_add_box(top, ll, bcell, -.5_dp, 1._dp, .true.)

  ll         = (/  0.12340E+02,  0.20100E+01,  0.00000E+00/)
  bcell(:,1) = (/  0.70000E+01,  0.00000E+00,  0.00000E+00/)
  bcell(:,2) = (/  0.00000E+00,  0.70000E+01,  0.00000E+00/)
  bcell(:,3) = (/  0.00000E+00,  0.00000E+00,  0.63500E+01/)
  call grid_add_box(top, ll, bcell, .5_dp, 1._dp, .true.)

  ll         = (/  0.12340E+02,  0.20100E+01,  0.25400E+02/)
  call grid_add_box(top, ll, bcell, .5_dp, 1._dp, .true.)

  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)
  call grid_setup(top,init=.true.)

  ! write out the initial cube file
!  call mg_save(top,'test6_initial',MG_SAVE_CUBE)

  call grid_hold_back(top)

  c1 = clock()

  call mg_gs_cds(top)

  time = timing(c1)

  print *,'Timing:', time

  c1 = clock()

  call mg_gs_cds(top,CDS_W_CYCLE)

  time = timing(c1)

  print *,'Timing:', time

  call grid_save_all(top,'test6')

  call delete_grid(top)
  
end program test
