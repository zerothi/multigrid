module m_time

  implicit none

  private

  integer, save :: cr, cm
  real, save :: rate

  public :: init_timing, clock, timing

  interface timing
     module procedure timing_c1
     module procedure timing_c2
  end interface timing

contains

  subroutine init_timing()

    ! First initialize the system_clock
    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = 1. / real(cr)

  end subroutine init_timing

  function clock()
    integer :: clock
    call system_clock(clock)
  end function clock

  function timing_c1(c1) result(time)
    integer, intent(in) :: c1
    real :: time
    integer :: c2
    c2 = clock()
    time = timing_c2(c1,c2)
  end function timing_c1
  function timing_c2(c1,c2) result(time)
    integer, intent(in) :: c1, c2
    real :: time
    time = ABS(c2 - c1) * rate
  end function timing_c2

end module m_time
