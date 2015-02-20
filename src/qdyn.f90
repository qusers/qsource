! (c) 2015 uppsala molekylmekaniska hb, uppsala, sweden
! qdyn.f90
! initial date: 2000
! by johan åqvist, john marelius, anders kaplan & martin nervall
! qdyn molecular dynamics main program

program qdyn
  use md
  use mpiglob ! use mpi global data

  implicit none
  ! version data
  character(10) :: qdyn_version = '5.7'
  character(12) :: qdyn_date = '2015-02-08'
#if defined (use_mpi)
  character(10) :: qdyn_suffix = '_parallel'
#else
  character(10) :: qdyn_suffix = ''
#endif
  
#if defined (use_mpi)
  ! mpi error code
  integer :: qdyn_ierr
#endif
  
  ! signal handler data and declarations
  integer(4) :: sigret
  integer(4), parameter :: sigint = 2 ! ctrl-c signal
  integer(4), parameter :: sigkill = 9 ! kill/ctrl-break signal
  integer(4), parameter :: sigabrt = 6 ! kill/ctrl-break signal

  external sigint_handler
  external sigkill_handler
  external sigabrt_handler

#if defined (use_mpi)
  ! initialize mpi
  call mpi_init(qdyn_ierr)
  if (qdyn_ierr .ne. mpi_success) call die('failure at mpi init')
  call mpi_comm_rank(mpi_comm_world, nodeid, qdyn_ierr)
  call mpi_comm_size(mpi_comm_world, numnodes, qdyn_ierr)
#else
  nodeid = 0
  numnodes = 1
#endif
  
  ! init signal handlers
  sigret = signal(sigint, sigint_handler, -1_4)
  sigret = signal(sigkill, sigkill_handler, -1_4)
  sigret = signal(sigabrt, sigabrt_handler, -1_4)

  ! initialize static data, display banner etc
  call startup

  if (nodeid .eq. 0) then
    ! master node: read input and initialize
    if (.not.initialize()) call die('invalid data in input file') ! read input data
    call open_files ! open necessary files
    call topology ! read topology
    call prep_coord ! read coords, solvates etc
    if (nstates > 0) call get_fep ! read fep/evb strategy
    !remove things with code 0 and maybe excluded bonded interactions
    call prep_sim ! prepare for simulation (calc. inv. mass, total charge,...)
    call close_input_files ! close input files
    call init_shake
    call make_nbqqlist
    call shrink_topology
    call nbmonitorlist
    call init_trj

    ! generate maxwellian velocities and shake initial x and v if necessary
    if (iseed > 0) then
      call maxwell
      call initial_shaking
      call stop_cm_translation
    end if
  end if

#if defined (use_mpi)
  ! initialize slave nodes
  if (numnodes .gt. 1) call init_nodes
#endif
  
  ! count non-bonded pairs to get the maximum number, then distribute them among the nodes
  call distribute_nonbonds

  ! do the work!
  call md_run

  if (nodeid .eq. 0) then
    ! master node: close output files
    call close_output_files
  end if

  ! deallocate memory etc.
  call shutdown

#if defined (use_mpi)
  ! shut down mpi
  call mpi_finalize(qdyn_ierr)
#endif
  
contains

  !-------------------------------------------------------------------------------
  ! startup code
  !-------------------------------------------------------------------------------
  subroutine startup
    integer :: i
    integer :: datum(8)

    if (nodeid .eq. 0) then
      ! start-of-header
      write (*, '(79a)') ('#', i = 1, 79)
#if defined (dum)
      write(*, '(a,a,a)') 'qdum input checker version ', trim(qdyn_version), ' initializing'
#elif defined(eval)
      write(*, '(a,a,a)') 'qdyn evaluation version ', trim(qdyn_version), ' initializing'
      write(*, '(a)') 'this version is for evaluation purposes only.'
      write(*, '(a)') 'optimisations are disabled - runs at <20% of maximum speed.'
#else
      write(*, '(a,a,a,a)') 'qdyn version ', trim(qdyn_version), trim(qdyn_suffix), ' initializing'
#endif
      call date_and_time(values = datum)
      write(*, 130) datum(1), datum(2), datum(3), datum(5), datum(6), datum(7)
      130 format('current date ', i4, '-', i2, '-', i2, ' and time ', i2, ':', i2, ':', i2)
    end if


    ! initialize used modules
    call md_startup

  end subroutine startup

  !-------------------------------------------------------------------------------
  ! shutdown code
  !-------------------------------------------------------------------------------
  subroutine shutdown
    integer :: i

    if (nodeid .eq. 0) then
#if defined (dum)
      write(*, *) 'qdum input checker version ', qdyn_version, ' terminated normally.'
#else
      write(*, *) 'qdyn version ', trim(qdyn_version), trim(qdyn_suffix), ' terminated normally.'
#endif
      write (*, '(79a)') ('#', i = 1, 79)
    end if

    ! call md's shutdown
    call md_shutdown
  end subroutine shutdown

  !-----------------------------------------------------------------------
  !interface
  integer(4) function signal(signum, proc, flag)
    integer(4) signum, flag
    external proc
  end function
  !end interface

end program qdyn

! signal handlers

integer(4) function sigint_handler(sig_num)
  use md
  implicit none
  integer(4) :: sig_num
  call die('user request (control-c)')
  sigint_handler = 1
end function sigint_handler

integer(4) function sigkill_handler(sig_num)
  use md
  implicit none
  integer(4) :: sig_num
  call die('kill signal')
  sigkill_handler = 1
end function sigkill_handler

integer(4) function sigabrt_handler(sig_num)
  use md
  implicit none
  integer(4) :: sig_num
  call die('kill signal')
  sigabrt_handler = 1
end function sigabrt_handler
