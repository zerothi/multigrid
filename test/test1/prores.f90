program prores_check

  use t_mg

  ! test libraries
  use m_mg_info

  implicit none

  ! the top multi-grid type
  type(mg_grid), target :: top

  integer :: i, N , nn(3)
  real(dp) :: cell(3,3)
  real(grid_p) :: tol, sor

  ! tolerance for the convergence
  tol = 1.e-3_grid_p
  sor = 1.8_grid_p

  ! initialize the initial grid
  cell(:,1) = (/2._dp,0._dp,0._dp/)
  cell(:,2) = (/0._dp,2._dp,0._dp/)
  cell(:,3) = (/0._dp,0._dp,2._dp/)

  i = 7
  do 
  i = i + 1
  if ( i == 16 ) then
     i = 128
  else if ( i > 100 ) then
     exit
  end if
  ! we create it to be 100x100x100
  nn = i
  ! create grid
  call init_grid(top,nn,cell,0,tol=tol)
  allocate(top%child)
  top%child%parent => top
  call init_grid(top%child,(nn+mod(i,2)) / 2 - 1, &
       top%cell, top%N_box, offset=top%offset)

  write(*,*)'>> Created initial grid...'

  do N = 1 , layers(top)
     call grid_set(top,layer=N,sor=sor,tol=tol)
  end do

  write(*,*)' >> Created all children...'

  call print_grid(top)

  ! initialize the grid
  call grid_bring_back(top)

  top%V = 1._grid_p

  ! initialize the grid
  call grid_bring_back(top%child)

  call grid_restriction(top)

  print '(a,e10.3)','Minimum value after restriction: ',minval(top%child%V)
  print '(a,3(i3,tr1))','   -- Minimum position         : ',minloc(top%child%V)
  print '(a,e10.3)','Maximum value after restriction: ',maxval(top%child%V)
  print '(a,3(i3,tr1))','   -- Maximum position         : ',maxloc(top%child%V)

  call grid_prolongation(top%child)

  print '(a,e10.3)','Minimum value after prolongation: ',minval(top%V)
  print '(a,3(i3,tr1))','   -- Minimum position         : ',minloc(top%V)
  print '(a,e10.3)','Maximum value after prolongation: ',maxval(top%V)
  print '(a,3(i3,tr1))','   -- Maximum position         : ',maxloc(top%V)

  call delete_grid(top)

  end do
  
end program prores_check
