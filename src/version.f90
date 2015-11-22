!------------------------------------------------------------------------------!
!  Q v.5.7 makefile                                                            !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Kajsa Ljunjberg, John Marelius, Martin Nervall                              !
!  Maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer, Mauricio Esguerra, !
!  Irek Szeler                                                                 !
!  latest update: November 22, 2015                                                !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! (c) 2015 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! version.f90
! initial date: 2015
! by Ireneusz Szeler
! q version and help print info
!------------------------------------------------------------------------------!

module version
  implicit none

contains


subroutine version_check(q_program, q_version, q_date, q_suffix)
! arguments
  character(*) :: q_program
  character(*) :: q_version
  character(*) :: q_date
  character(*) :: q_suffix
  character(len=32) :: arg

! local
  logical :: fin
  integer :: i
  integer :: count
  integer :: datum(8)
  character(200) :: flag

!  if (nodeid .eq. 0) then
  fin = .false.
  count = command_argument_count()
  if (count .gt. 0) then
    call getarg(1,flag)
    call lowcase(flag)
    if ( flag .eq. '-v' .or. flag .eq. '-version' .or. flag .eq. '--version') fin = .true. 
  end if

! start-of-header
  write (*,'(80a)') ('-',i=1,80)
  if ( q_program .eq. 'qdyn') then
#if defined (dum)
    write(*,'(a,a,a)') 'qdum input checker version ', trim(q_version), ' initializing'
#elif defined(eval)
    write(*,'(a,a,a)') 'qdyn evaluation version ', trim(q_version), ' initializing'
    write(*,'(a)') 'This version is for evaluation purposes only.'
    write(*,'(a)') 'Optimizations are disabled - runs at <20% of maximum speed.'
#endif
  end if

  call date_and_time(values=datum)
  write(*,130) datum(1),datum(2),datum(3)
  write(*,140) datum(5),datum(6),datum(7)
130 format('Program run date = ',i4,'-',i2,'-',i2)
140 format('Program run time = ',i2,':',i2,':',i2)

!  call version_print(q_program, q_version, q_date, q_suffix)
  call lowcase(q_program)
  if ( flag .eq. '-h' .or. flag .eq. '-help' .or. flag .eq. '--help') then
  fin = .true.
  write(*,*)
  write(*,'(a,a)') trim(q_program), ' help information'
  select case (q_program)
    case ('qdyn')
      write(*,*) 
      write(*,'(a)') 'To run calculations use: '
      write(*,'(a,a,a)') '    ',trim(q_program), '5 inputfile.inp > output.file'
      write(*,'(a)') ' or for parallel version'
      write(*,'(a,a,a)') '    mpienvironmment ', trim(q_program), '5p inputfile.inp > output.file'
      write(*,*) 
      write(*,'(a)') 'where:'
      write(*,'(a)') 'mpienvironment - e.q. mpirun -n 4, for more info check cluster information '
      write(*,*) 
    case ('qfep')
      write(*,*) 
      write(*,'(a,a)') 'You are running program: ', trim(q_program) 
      write(*,*) 
    case ('qprep')
      write(*,*) 
      write(*,'(a,a)') 'You are running program: ', trim(q_program) 
      write(*,*) 
    case ('qcalc')
      write(*,*) 
      write(*,'(a,a)') 'You are running program: ', trim(q_program) 
      write(*,*) 
    case default
      write(*,*) 
      write(*,'(a)') 'Is this a new program added to Q?'
      write(*,*) 
    end select
  end if

  if(fin) stop

  select case (q_program)
    case ('qdyn')
      write(*,*)
    case ('qfep')
      write(*,*)
      write(*,'(a,a)') 'Welcome to: ', trim(q_program)
      write(*,'(a,a)') 'Latest changes in the source code: ', trim(q_date) 
      write(*,*)
    case ('qprep')
      write(*,*)
      write(*,'(a,a)') 'Welcome to: ', trim(q_program)
      write(*,'(a,a)') 'Latest changes in the source code: ', trim(q_date)
      write(*,*)
    case ('qcalc')
      write(*,*)
      write(*,'(a,a)') 'Welcome to: ', trim(q_program)
      write(*,'(a,a)') 'Latest changes in the source code: ', trim(q_date) 
      write(*,*) 
    case default
      write(*,*)
      write(*,'(a)') 'Welcome to ...' 
      write(*,*)
  end select
!  end if ! node .eq. 0
end subroutine version_check


elemental subroutine lowcase(word)
  character(*), intent(in out) :: word
  integer :: i, ic, nlen
  nlen = len(word) 
  do i=1, nlen
    ic = ichar(word(i:i))
    if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
  end do
end subroutine lowcase


subroutine version_print(q_program, q_version, q_date, q_suffix)
! arguments
  character(*) :: q_program
  character(*) :: q_version
  character(*) :: q_date
  character(*) :: q_suffix

! local
  integer :: i
  integer :: datum(8)

! start-of-header
  write(*,*)
  write(*,'(a)') 'Build and version information'
  write(*,*)
!  write (*,'(79a)') ('#',i=1,79)
  if ( q_program .eq. 'qdyn') then
#if defined (dum)
    write(*,'(a,a,a)') 'qdum input checker version ', trim(q_version), ' initializing'
#elif defined(eval)
    write(*,'(a,a,a)') 'qdyn evaluation version ', trim(q_version), ' initializing'
    write(*,'(a)') 'This version is for evaluation purposes only.'
    write(*,'(a)') 'Optimizations are disabled - runs at <20% of maximum speed.'
#endif
  end if

  call date_and_time(values=datum)
  write(*,130) datum(1),datum(2),datum(3)
  write(*,140) datum(5),datum(6),datum(7)
130 format('Program run date = ',i4,'-',i2,'-',i2)
140 format('Program run time = ',i2,':',i2,':',i2)

end subroutine version_print


end module version
