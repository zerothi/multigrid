module t_bc

  implicit none

  ! a module to sustain a "simple" multi-grid solver
  integer, parameter :: MG_BC_A0 = 1
  integer, parameter :: MG_BC_A1 = 2
  integer, parameter :: MG_BC_B0 = 4
  integer, parameter :: MG_BC_B1 = 8
  integer, parameter :: MG_BC_C0 = 16
  integer, parameter :: MG_BC_C1 = 32

  integer, parameter :: MG_BC_DIRICHLET = 1
  integer, parameter :: MG_BC_NEUMANN = 2

  type :: tBC
     sequence
     ! Method
     integer :: method = MG_BC_DIRICHLET
  end type tBC
     
end module t_bc
