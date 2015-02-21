! (c) 2014 uppsala molekylmekaniska hb, uppsala, sweden
! version.f90
! initial date: 2015
! by ireneusz szeler
! q version and help print info

module versions

  implicit none

contains

subroutine version_check(q_program, q_version, q_date, q_suffix)

! arguments
  character(*)  :: q_program
  character(*)  :: q_version
  character(*)  :: q_date
  character(*)  :: q_suffix

! local
  logical  :: fin
  integer  :: num_arg
  character(200) :: flag

!  if (nodeid .eq. 0) then
    
    fin = .false.
    num_arg = command_argument_count()
    if (num_arg .gt. 0) then
      call getarg(1,flag)
      call lowcase(flag)
      if ( flag .eq. '-v' .or. flag .eq. '-version' .or. flag .eq. '--version') fin = .true. 
    end if

    call version_print(q_program, q_version, q_date, q_suffix)
    
    call lowcase(q_program)
    
    if ( flag .eq. '-h' .or. flag .eq. '-help' .or. flag .eq. '--help') then
       fin = .true.
       write(*,*) 
       write(*,'(a,a)') trim(q_program), ' help information'
       select case (q_program)
         case ('qdyn')
             write(*,*) 
             write(*,'(a)') 'to run calculations use: '
             write(*,'(a,a,a)') '    ',trim(q_program), '5 inputfile.inp > output.file'
             write(*,'(a)') ' or for parallel version'
             write(*,'(a,a,a)') '    mpienvirment ', trim(q_program), '5p inputfile.inp > output.file'
             write(*,*) 
             write(*,'(a)') 'where:'
             write(*,'(a)') 'mpienvirment - e.q. mpirun -n 4, for more info check cluster informations '
             write(*,*) 
         case ('qfep')
             write(*,*) 
             write(*,'(a,a)') 'in this moment no info available for ', trim(q_program) 
             write(*,*) 
         case ('qprep')
             write(*,*) 
             write(*,'(a,a)') 'in this moment no info available for ', trim(q_program) 
             write(*,*) 
         case ('qcalc')
             write(*,*) 
             write(*,'(a,a)') 'in this moment no info available for ', trim(q_program) 
             write(*,*) 
         case default
             write(*,*) 
             write(*,'(a)') 'it is some new program added to q?'
             write(*,*) 
         end select
    end if
  
    if(fin) stop

  select case (q_program)
    case ('qdyn')
      write(*,*) 
    case ('qfep')
      write(*,*) 
      write(*,'(a,a,a,a)') 'welcome in ', trim(q_program), ' modification date ', trim(q_date) 
      write(*,*) 
    case ('qprep')
      write(*,*) 
      write(*,'(a,a,a,a)') 'welcome in ', trim(q_program), ' modification date ', trim(q_date) 
      write(*,*) 
    case ('qcalc')
      write(*,*) 
      write(*,'(a,a,a,a)') 'welcome in ', trim(q_program), ' modification date ', trim(q_date) 
      write(*,*) 
    case default
      write(*,*) 
      write(*,'(a)') 'welcome in ...' 
      write(*,*) 
  end select


!  end if ! node .eq. 0

end subroutine version_check

elemental subroutine lowcase(word)

 character(*), intent(in out) :: word
 integer   :: i, ic, nlen

 nlen = len(word) 
 do i=1, nlen
   ic = ichar(word(i:i))
   if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
 end do

end subroutine lowcase

subroutine version_print(q_program, q_version, q_date, q_suffix)

! arguments
  character(*)  :: q_program
  character(*)  :: q_version
  character(*)  :: q_date
  character(*)  :: q_suffix

! local
  integer  :: i
  integer  :: datum(8)

! start-of-header
  write(*,*)
  write(*,'(a)') 'build and version information'
  write(*,*)
!  write (*,'(79a)') ('#',i=1,79)
  if ( q_program .eq. 'qdyn') then
#if defined (dum)
    write(*,'(a,a,a)') 'qdum input checker version ', trim(q_version), ' initialising'
#elif defined(eval)
    write(*,'(a,a,a)') 'qdyn evaluation version ', trim(q_version), ' initialising'
    write(*,'(a)') 'this version is for evaluation purposes only.'
    write(*,'(a)') 'optimisations are disabled - runs at <20% of maximum speed.'
#endif
  end if
#if defined (build_username) && defined (build_hostname) && defined (build_date) && defined (build_source) && defined (build_number) && defined(build_compiler)
  write(*,'(a,a)') 'build number ', build_number
  write(*,'(a,a)') 'build date   ', build_date
  write(*,'(a)')   'built:       '
  write(*,'(a,a)') '      by     ', build_username
  write(*,'(a,a)') '      on     ', build_hostname
  write(*,'(a,a)') '      git id ', build_source
  write(*,'(a,a)') '      with   ', build_compiler
#else
  write(*,'(a,a,a,a,a)')  trim(q_program), ' version ', trim(q_version), trim(q_suffix),' initialising'
#endif
  call date_and_time(values=datum)
  write(*,130) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)
130	format('current date ',i4,'-',i2,'-',i2,' and time ',i2,':',i2,':',i2)

end subroutine version_print

end module versions
