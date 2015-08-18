
module m_unit

  integer, parameter, private :: dp = selected_real_kind(14,100)

  real(dp), parameter :: Ang    = 1._dp / 0.529177_dp
  real(dp), parameter :: eV     = 1._dp / 13.60580_dp
  
  real(dp), parameter :: Pi = 3.1415926535897932384626433832795028_dp
  real(dp), parameter :: deg = Pi / 180.0_dp
  
end module m_unit
