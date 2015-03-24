! Module for reading an input file.

! This code has been fully created by 
! Nick Papior Andersen
! nickpapior@gmail.com
! 2015

module m_io

  implicit none

  integer, parameter :: FILE_LENGTH = 150

  private
  public :: tIO
  public :: io_crt
  public :: io_open
  public :: io_destroy
  public :: io_close
  public :: io_line
  public :: io_step

  ! Generic string routines
  public :: has_sub
  public :: startswith
  public :: strip
  public :: lcase

  type :: tIO
     ! Handle to read in multi-grid information
     ! in a small txt file
     integer :: u = -1 ! unit
     character(FILE_LENGTH) :: file = ' '
     ! Current read line (stepping allows reread until il is seen again)
     integer :: il
  end type tIO

contains

  subroutine io_crt(IO,file)
    type(tIO), intent(inout) :: IO
    character(len=*), intent(in) :: file
    logical :: exists

    ! Ensure that it is clean
    call io_destroy(IO)

    IO%file = file

    inquire(file=IO%file, exist = exists)
    if ( .not. exists ) then
       print '(a)', 'Can not find file: '//trim(IO%file)
       stop
    end if

  end subroutine io_crt

  subroutine io_open(IO)
    type(tIO), intent(inout) :: IO

    call io_close(IO)
    IO%u = next_unit()

    open(unit = IO%u, file = IO%file, &
         form = 'formatted', &
         status = 'old')
    IO%il = 0

  end subroutine io_open

  subroutine io_destroy(IO)
    type(tIO), intent(inout) :: IO
    call io_close(IO)
    IO%file = ' '
  end subroutine io_destroy

  subroutine io_close(IO)
    type(tIO), intent(inout) :: IO
    logical :: opened

    inquire(unit=IO%u,opened=opened)
    if ( opened ) then
       close(IO%u)
    end if

    IO%u = -1
    IO%il = 0

  end subroutine io_close

  function has_sub(str,sub) result(has)
    character(len=*), intent(in) :: str, sub
    logical :: has 
    has = index(str,sub) > 0
  end function has_sub
  
  function strip(l,n) result(ol)
    character(len=*), intent(in) :: l
    integer, intent(in), optional :: n
    character(len=len(l)) :: ol

    integer :: i, ln

    ln = 1
    if ( present(n) ) ln = n

    ol = l
    do while ( ln > 0 )
       ln = ln - 1
       i = index(ol,' ')
       ol = adjustl(ol(i+1:))
    end do

  end function strip

  function startswith(str,sub) result(starts)
    character(len=*), intent(in) :: str, sub
    logical :: starts
    starts = index(str,sub) == 1
  end function startswith

  function io_line(IO) result(line)
    type(tIO), intent(inout) :: IO
    character(len=300) :: line
    integer :: i
    
    line = ' '
    do while ( len_trim(line) == 0 )
       read(IO%u,'(a)',iostat=i) line
       if ( i /= 0 ) then
          ! Signal the end of the file
          line = '#'
          return
       end if
       IO%il = IO%il + 1
       ! Check that it is not a comment
       line = adjustl(line)
       i = scan(line,'#!')
       if ( i > 0 ) then
          ! Clear comments
          line(i:) = ' '
       end if
       line = lcase(line)
    end do

  end function io_line

  function io_step(IO,keyword) result(line)
    type(tIO), intent(inout) :: IO
    character(len=*), intent(in) :: keyword
    character(len=len_trim(keyword)) :: lkeyword
    character(len=300) :: line
    integer :: old_il, in_block
    logical :: reopen

    reopen = .false.
    in_block = 0
    old_il = IO%il
    ! Return to lower-case keyword
    lkeyword = lcase(trim(keyword))

    ! This will pass all comments and will lower-case the line
    line = io_line(IO)

    do while ( .not. has_sub(line,lkeyword) )

       if ( reopen .and. old_il <= IO%il ) then
          ! We have re-read the file and gotten
          ! to the same point again
          line = '#'
          return
       end if
       
       ! If we are in a block, we do not look for keywords
       if ( in_block > 0 ) then
          if ( startswith(line,'end') ) in_block = in_block - 1
       else
          if ( startswith(line,'begin') ) in_block = in_block + 1
          if ( in_block == 1 .and. has_sub(line,lkeyword) ) exit
       end if

       line = io_line(IO)

       if ( line(1:1) == '#' ) then
          call io_close(IO)
          call io_open(IO)
          reopen = .true.
       end if

    end do

  end function io_step
    
  ! Lower-case a string
  pure function lcase(str) 
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lcase
    integer :: ic, i

    character(len=26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'

    ! Capitalize each letter if it is lowecase
    lcase = str
    i = scan(lcase,upper)
    do while ( i > 0 ) 
       ! Get the conversion index
       ic = index(upper,lcase(i:i))
       lcase(i:i) = lower(ic:ic)
       ic = scan(lcase(i+1:),upper)
       if ( ic > 0 ) then
          i = i + ic
       else
          i = 0
       end if
    end do

  end function lcase

  ! Return the next available unit for the system
  function next_unit() result(u)
    integer :: u
    logical :: opened
    u = 9
    opened = .true.
    do while ( opened )
       u = u + 1
       inquire(unit=u,opened=opened)
    end do
  end function next_unit

end module m_io
