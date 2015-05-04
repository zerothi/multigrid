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
  public :: ccase, lcase, ucase
  public :: next_unit

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

    if ( IO%u < 0 ) then
       opened = .false.
    else
       inquire(unit=IO%u,opened=opened)
    end if
    if ( opened ) close(IO%u)

    IO%u = -1
    IO%il = 0

  end subroutine io_close

  function has_sub(str,sub,word) result(has)
    character(len=*), intent(in) :: str, sub
    ! Whether the substring should be a word (with surrounding spaces)
    logical, intent(in), optional :: word
    integer :: i 
    logical :: has 
    i = index(str,sub)
    if ( present(word) ) then
       if ( word ) then
          ! single word somewhere in the line
          i = index(str,' '//trim(sub)//' ')
          if ( i == 0 ) then
             ! single word at the beginning of the line
             i = index(str,trim(sub)//' ')
          end if
       end if
    end if
    has = i > 0
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

  function io_line(IO,case) result(line)
    type(tIO), intent(inout) :: IO
    character(len=1), intent(in), optional :: case
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
       if ( present(case) ) then
          if ( case == 'L' .or. case == 'l' ) &
               line = lcase(line)
          if ( case == 'U' .or. case == 'u' ) &
               line = ucase(line)
       else
          line = lcase(line)
       end if
    end do

  end function io_line

  function io_step(IO,keyword,case) result(line)
    type(tIO), intent(inout) :: IO
    character(len=*), intent(in) :: keyword
    character(len=len_trim(keyword)) :: lkeyword
    character(len=1), intent(in), optional :: case
    character(len=300) :: line
    integer :: old_il, in_block
    logical :: reopen

    reopen = .false.
    in_block = 0
    old_il = IO%il
    ! Return to lower-case keyword
    lkeyword = lcase(trim(keyword))

    ! This will pass all comments and will lower-case the line
    line = io_line(IO,case=case)

    do while ( .not. has_sub(line,lkeyword,word = .true. ) )

       if ( reopen .and. old_il <= IO%il ) then
          ! We have re-read the file and gotten
          ! to the same point again
          line = '#'
          return
       end if
       
       ! If we are in a block, we do not look for keywords
       if ( startswith(line,'begin') ) then
          in_block = in_block + 1
          if ( in_block == 1 .and. has_sub(line,lkeyword, word = .true. ) ) exit
       end if

       line = io_line(IO,case=case)

       if ( startswith(line,'end') ) then
          if ( in_block > 0 ) in_block = in_block - 1
          line = io_line(IO)
       end if

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

  ! Upper-case a string
  pure function ucase(str) 
    character(len=*), intent(in) :: str
    character(len=len(str)) :: ucase
    integer :: ic, i

    character(len=26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'

    ! Capitalize each letter if it is lowecase
    ucase = str
    i = scan(ucase,lower)
    do while ( i > 0 ) 
       ! Get the conversion index
       ic = index(lower,ucase(i:i))
       ucase(i:i) = upper(ic:ic)
       ic = scan(ucase(i+1:),lower)
       if ( ic > 0 ) then
          i = i + ic
       else
          i = 0
       end if
    end do

  end function ucase

  ! Upper-case a string
  pure function ccase(str,case) 
    character(len=*), intent(in) :: str
    character(len=1), intent(in), optional :: case
    character(len=len(str)) :: ccase

    if ( present(case) ) then
       if ( scan(case,'Uu') > 0 ) then
          ccase = ucase(str)
       else if ( scan(case,'Ll') > 0 ) then
          ccase = lcase(str)
       else
          ccase = str
       end if
    else
       ccase = str
    end if
  end function ccase

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
