!------------------------------------------------------------------------------!
!  Q v.5.7 makefile                                                            !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Kajsa Ljunjberg, John Marelius, Martin Nervall                              !
!  Maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer, Mauricio Esguerra, !
!  Irek Szeler                                                                 !
!  latest update: february 21, 2015                                            !
!  make [option] [comp=compiler]                                               !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
!  qdyn.f90
!  by Johan Aqvist, John Marelius, Anders Kaplan & Martin Nervall
!  Qdyn molecular dynamics main program
!------------------------------------------------------------------------------!

program qdyn
  use MD                                                                
  use MPIGLOB ! use MPI global data

#if defined (_DF_VERSION_)
  use dfport ! portability lib for signals. Used in windows.
#endif

  implicit none

  ! version data
  character(10)             :: QDYN_VERSION = '5.7'
  character(12)             :: QDYN_DATE = '2015-02-22'
#if defined (USE_MPI)
  character(10)           :: QDYN_SUFFIX = '_parallel'
#else
  character(10)           :: QDYN_SUFFIX = ''
#endif


#if defined (USE_MPI)
  ! MPI error code
  integer                 :: qdyn_ierr
#endif

  ! signal handler data and declarations
  integer(4)                :: sigret
#if defined (_DF_VERSION_)
  ! nothing
#else
  integer(4), parameter   :: SIGINT = 2   ! CTRL-C signal
  integer(4), parameter   :: SIGKILL = 9  ! kill/CTRL-BREAK signal
  integer(4), parameter   :: SIGABRT = 6  ! kill/CTRL-BREAK signal
#endif

  external sigint_handler
  external sigkill_handler
  external sigabrt_handler

#if defined (USE_MPI)
  ! initialize MPI
  call MPI_Init(qdyn_ierr)
  if (qdyn_ierr .ne. MPI_SUCCESS) call die('failure at MPI init')
  call MPI_Comm_rank(MPI_COMM_WORLD, nodeid, qdyn_ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numnodes, qdyn_ierr)
#else
  nodeid = 0
  numnodes = 1
#endif

  ! init signal handlers
  sigret = SIGNAL(SIGINT, sigint_handler, -1_4)
  sigret = SIGNAL(SIGKILL, sigkill_handler, -1_4)
  sigret = SIGNAL(SIGABRT, sigabrt_handler, -1_4)

  ! initialize static data, display banner etc
  call startup

  if (nodeid .eq. 0) then
     ! master node: read input and initialize
     if(.not. initialize()) call die('Invalid data in input file') ! read input data
     call open_files ! open necessary files
     call topology ! read topology
     call prep_coord ! read coords, solvates etc
     if ( nstates > 0 ) call get_fep ! read fep/evb strategy
     !remove things with code 0 and maybe excluded bonded interactions

     call prep_sim ! prepare for simulation (calc. inv. mass, total charge,...)
     call close_input_files ! close input files

     call init_shake
     call make_nbqqlist
     call shrink_topology
     call nbmonitorlist
     call init_trj

     ! generate Maxwellian velocities and shake initial x and v if necessary
     if ( iseed > 0 ) then
        call maxwell
        call initial_shaking
        call stop_cm_translation
     end if
  end if

#if defined (USE_MPI)
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

#if defined (USE_MPI)
  ! shut down MPI
  call MPI_Finalize(qdyn_ierr)
#endif

contains

  !-------------------------------------------------------------------------------
  ! startup subroutine
  !-------------------------------------------------------------------------------
  subroutine startup
    integer :: i
    integer :: datum(8)

    if (nodeid .eq. 0) then
       ! start-of-header
       write (*,'(79a)') ('#',i=1,79)

#if defined (DUM)
       write(*,'(a,a,a)') 'QDum input checker version ', trim(QDYN_VERSION), ' initializing'
#elif defined(EVAL)
       write(*,'(a,a,a)') 'QDyn evaluation version ', trim(QDYN_VERSION), ' initializing'
       write(*,'(a)') 'This version is for evaluation purposes only.'
       write(*,'(a)') 'Optimizations are disabled - runs at <20% of maximum speed.'
#else
       write(*,'(a,a,a,a)') 'QDyn version ', trim(QDYN_VERSION), trim(QDYN_SUFFIX),' initializing'
#endif

       call date_and_time(values=datum)
       write(*,130) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)  
130    format('Current date ',i4,'-',i2,'-',i2,' and time ',i2,':',i2,':',i2)
    end if


    ! initialize used modules
    call md_startup

  end subroutine startup

  !-------------------------------------------------------------------------------
  ! shutdown subroutine
  !-------------------------------------------------------------------------------
  subroutine shutdown
    integer :: i

    if (nodeid .eq. 0) then
#if defined (DUM)
       write(*,*) 'QDum input checker version ', QDYN_VERSION, ' terminated normally.'
#else
       write(*,*) 'QDyn version ', trim(QDYN_VERSION), trim(QDYN_SUFFIX), ' terminated normally.'
#endif
       write (*,'(79a)') ('#',i=1,79)
    end if

    ! call md's shutdown
    call md_shutdown

  end subroutine shutdown

#if ! defined(__PATHCC__)
  integer(4) function signal( signum, proc, flag )
    implicit none
    integer(4) :: signum, flag
    external proc
    signal = 1
  end function signal
#endif


end program qdyn


!-------------------------------------------------------------------------------
! signal handlers
!-------------------------------------------------------------------------------
INTEGER(4) FUNCTION sigint_handler(sig_num)
  use MD
  implicit none
  INTEGER(4) :: sig_num

  call die('user request (control-C)')
  sigint_handler = 1
END FUNCTION sigint_handler

INTEGER(4) FUNCTION sigkill_handler(sig_num)
  use MD
  implicit none
  INTEGER(4) :: sig_num

  call die('kill signal')
  sigkill_handler = 1
END FUNCTION sigkill_handler

INTEGER(4) FUNCTION sigabrt_handler(sig_num)
  use MD
  implicit none
  INTEGER(4) :: sig_num

  call die('kill signal')
  sigabrt_handler = 1
END FUNCTION sigabrt_handler
