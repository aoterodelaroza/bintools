module bintools_procedures
  implicit none
  public
contains

  ! convert one or more 64-bit hex values in pat0 to doubles
  subroutine hex2double(pat0)
    implicit none
    character*(*), intent(in) :: pat0

    character*(:), allocatable :: pat
    integer :: i, plen2
    integer*1, allocatable :: arr(:)
    real*8 :: dbl

    ! check the pattern
    pat = trim(pat0)
    call check_hex_pattern(pat)

    ! allocate space for the conversion
    plen2 = len(pat)/2
    allocate(arr(plen2))

    ! do the conversion
    do i = 1, plen2
       read (pat(2*i-1:2*i),'(z2)') arr(i)
    end do

    ! convert and write the doubles
    do i = 1, plen2, 8
       dbl = transfer(arr(i:min(i+7,plen2)),dbl)
       write (*,'(1p,E22.15)') dbl
    end do

  end subroutine hex2double

  ! convert a double to a 64-bit hex value
  subroutine double2hex(dbl)
    implicit none
    real*8, intent(in) :: dbl

    integer*1 :: arr(8)
    integer :: i

    arr = transfer(dbl,arr)

    do i = 1, 8
       write (*,'(z2.2)',advance="no") arr(i)
    end do
    write (*,*)

  end subroutine double2hex

  ! Search for hex pattern pat0 in a binary file (stream access). Print
  ! the locations in hex.
  subroutine grep(pat0,file)
    implicit none
    character*(*), intent(in) :: pat0, file

    integer, parameter :: lu = 10

    integer :: i, ios, plen2
    character*(:), allocatable :: pat
    integer*1 :: byte
    integer*8 :: offset, ishift, nbuf, nbuf2
    integer*1, allocatable :: ipat(:), ibuf(:), ibuf2(:)

    ! open file
    open(file=file,unit=lu,status='old',form='unformatted',access='stream',iostat=ios)
    if (ios /= 0) then
       write (*,'("Error opening file: ",A)') file
       call print_help_and_exit()
    end if

    ! check the pattern
    pat = trim(pat0)
    call check_hex_pattern(pat)
    plen2 = len(pat)/2

    ! convert the pattern to integer*1
    allocate(ipat(0:plen2-1),ibuf(0:plen2-1),ibuf2(0:plen2-1))
    do i = 0, plen2-1
       read (pat(2*i+1:2*i+2),'(z2)') ipat(i)
    end do

    ! process the file
    offset = -1
    ishift = 0
    nbuf = 0
    nbuf2 = 0
    do while (.true.)
       ! read next byte from the primary buffer or the file
       offset = offset + 1
       if (nbuf > 0) then
          nbuf = nbuf - 1
          byte = ibuf(nbuf)
       else
          read(lu,iostat=ios) byte
       end if

       ! check we are not at the end of the file
       if (is_iostat_end(ios)) exit
       if (ios /= 0) then
          write (*,'("Error reading file: ",A)') file
          call print_help_and_exit()
       end if

       ! compare to the pattern
       if (byte == ipat(ishift)) then
          ! advance the pattern pointer and fill the backtrack buffer
          ishift = ishift + 1
          ibuf2(nbuf2) = byte
          nbuf2 = nbuf2 + 1

          ! pattern match - write and discard the backtrack buffer
          if (ishift >= plen2) then
             ishift = 0
             write (*,'(z8)') offset-plen2+1
             nbuf2 = 0
          end if
       else
          ! no match, reset the pattern pointer
          ishift = 0

          ! if the backtrack buffer has any contents, move them to the primary buffer;
          ! discard the first element from the backtrack buffer and add the last byte read;
          ! set the offset to the correct value;
          if (nbuf2 > 0) then
             ibuf(nbuf:nbuf+nbuf2-2) = ibuf2(1:nbuf2-1)
             ibuf(nbuf+nbuf2-1) = byte
             nbuf = nbuf + nbuf2
             offset = offset - nbuf2
          end if

          ! discard the contents of the backtrack buffer
          nbuf2 = 0
       end if
    end do

    close(lu)

  end subroutine grep

  ! check that a hex value is sane
  subroutine check_hex_pattern(pat)
    implicit none
    character*(*), intent(inout) :: pat

    integer :: i
    integer :: plen

    character*(*), parameter :: allowedhex = '0123456789abcdef'

    ! check the pattern
    call lowercase(pat)
    plen = len(pat)
    if (mod(plen,2) /= 0) then
       write (*,'("The pattern must have an even number of hexadecimal digits.")')
       call print_help_and_exit()
    end if
    do i = 1, plen
       if (index(allowedhex,pat(i:i)) == 0) then
          write (*,'("Not a hexadecimal digit (",A1,") in pattern: ",A)') pat(i:i), pat
          call print_help_and_exit()
       end if
    end do

  end subroutine check_hex_pattern

  ! print help message and stop the program
  subroutine print_help_and_exit()
    implicit none

    write (*,'("Usage:")')
    write (*,'("  bintools GREP sequence file")')
    write (*,'("  bintools DOUBLE2HEX double")')
    write (*,'("  bintools HEX2DOUBLE sequence")')
    stop

  end subroutine print_help_and_exit

  ! convert a string to lowercase
  subroutine lowercase(str)
    implicit none
    character*(*), intent(inout) :: str

    character(*), parameter :: lo = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: i, idx

    do i = 1, len(str)
       idx = index(up,str(i:i))
       if (idx > 0) str(i:i) = lo(idx:idx)
    enddo

  end subroutine lowercase

end module bintools_procedures

program bintools
  use bintools_procedures
  implicit none

  integer, parameter :: arglen = 1024

  character(len=arglen) :: argv
  integer :: argc
  character*(:), allocatable :: str, pat, file
  real*8 :: dbl

  argc = command_argument_count()
  if (argc == 0) call print_help_and_exit()

  call getarg(1,argv)
  str = trim(argv)
  call lowercase(str)

  if (str == "grep") then
     ! GREP sequence file
     if (argc /= 3) call print_help_and_exit()
     call getarg(2,argv)
     pat = trim(argv)
     call getarg(3,argv)
     file = trim(argv)
     call grep(pat,file)
  elseif (str == "double2hex" .or. str == "doubletohex") then
     ! DOUBLE2HEX double
     if (argc /= 2) call print_help_and_exit()
     call getarg(2,argv)
     read(argv,*) dbl
     call double2hex(dbl)
  elseif (str == "hex2double" .or. str == "hextodoube") then
     ! DOUBLE2HEX double
     if (argc /= 2) call print_help_and_exit()
     call getarg(2,argv)
     pat = trim(argv)
     call hex2double(pat)
  else
     call print_help_and_exit()
  end if

end program bintools
