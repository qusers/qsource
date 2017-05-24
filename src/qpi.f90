! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! qpi.f90
! based on qdyn.f90
! by Johan Åqvist, John Marelius, Anders Kaplan & Martin Nervall
! main routine for QCP analysis tool
! by Paul Bauer

! this will look like a copy of qdyn,f90, because it needs to do a large part of the same work  done there
! so we also initialize the same variables and stuff, but exclude everything only needed during the actual md

program QPI5
  use QCP
  use VERSIONS
  use MPIGLOB ! use MPI global data
  use SIMPREP
#if defined (_DF_VERSION_)
  use dfport  ! portability lib for signals, used by Compaq Visual Fortran compiler
#endif
!$ use omp_lib

  implicit none
  ! version data
  character(10)					:: QPI_NAME = 'Qpi'
  character(80)					:: QPI_VERSION = ''
  character(80)					:: QPI_DATE    = ''
#if defined (USE_MPI)
  character(10)					:: QPI_SUFFIX = '_parallel'
#else
  character(10)					:: QPI_SUFFIX = ''
#endif


#if defined (USE_MPI)
  ! MPI error code
  integer						:: qpi_ierr
#endif

  ! signal handler data and declarations
  integer(4)					:: sigret
#if defined (_DF_VERSION_)
  ! nothing
#else
  integer(4), parameter			:: SIGINT  = 2 ! CTRL-C signal
  integer(4), parameter			:: SIGABRT = 6 ! kill/CTRL-BREAK signal
  integer(4), parameter			:: SIGKILL = 9 ! kill/CTRL-BREAK signal
#endif
  external sigint_handler
  external sigkill_handler
  external sigabrt_handler

! read in version info
  QPI_VERSION = trim(version_pass())
  QPI_DATE    = trim(date_pass())


#if defined (USE_MPI)
  ! initialize MPI
  call MPI_Init(qpi_ierr)
  if (qpi_ierr .ne. MPI_SUCCESS) call die('failure at MPI init')
  call MPI_Comm_rank(MPI_COMM_WORLD, nodeid, qpi_ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numnodes, qpi_ierr)
#else
	nodeid = 0
	numnodes = 1
#endif

  ! initialize signal handlers
  sigret = qsignal(SIGINT, sigint_handler, -1_4)
  sigret = qsignal(SIGKILL, sigkill_handler, -1_4)
  sigret = qsignal(SIGABRT, sigabrt_handler, -1_4)
#if defined(__INTEL_COMPILER) || defined(__PGI)
      call signal(SIGINT, sigint_handler , -1)
      call signal(SIGABRT,sigkill_handler, -1)
      call signal(SIGKILL,sigabrt_handler, -1)
#elif defined(__GFORTRAN__) || defined(__PATHSCALE__)
      call signal(SIGINT, sigint_handler )
      call signal(SIGABRT,sigkill_handler)
      call signal(SIGKILL,sigabrt_handler)
#else
#error "This code is inteded for use with Intel, PGI, GNU and PATHSCALE compiler. Please add a signalhandler for your compiler."
#endif




  ! initialise static data, display banner etc
  call startup

  if (nodeid .eq. 0) then
#ifdef _OPENMP
!$omp parallel
  threads_num = omp_get_num_threads()
!$omp end parallel
#endif

	! master node: read input and initialise

	if(.not. qcp_initialize()) call die('Invalid data in input file')						! read input data
	call open_files(md=.false.)		! open necessary files
	call topology				! read topology
	call qcp_prep_coord				! read coords, get number of configurations
	call get_fep	! read fep/evb strategy
	!remove things with code 0 and maybe excluded bonded interactions

	call prep_sim						! prepare for simulation (calc. inv. mass, total charge,...)
	call prep_sim_version(version_pass())
	call close_input_files				! close input files

	call init_constraints
!the nb monitor now needs to be initialised after! we got the precomputed interactions
!so it is moved in later 
	call make_nbqqlist
	call shrink_topology

  end if

#if defined (USE_MPI)
  ! initialise data types
  call set_mpi_types
  ! initialise slave nodes
  if (numnodes .gt. 1) call init_nodes
#endif
	! count non-bonded pairs to get the maximum number, then distribute them among the nodes

  call distribute_nonbonds

! before doing actual work, make each node do a precomputation of all possible
! interactions
! so we don't have to do this at every step
! has to be after distribute_nonbonds because we need to know what is actually
! needed on each node
  call precompute_interactions

 ! do the work!
  call qcp_work

  if (nodeid .eq. 0) then
	! master node: close output files
	close (11)
  end if

  ! deallocate memory etc.
  call qpi_shutdown

#if defined (USE_MPI)
  ! shut down MPI
  call MPI_Finalize(qpi_ierr)
#endif

contains

!-----------------------------------------------------------------------

! startup/shutdown code

subroutine startup

  if (nodeid .eq. 0) then
    call version_check(QPI_NAME, QPI_VERSION, QPI_DATE, QPI_SUFFIX) ! print version and chack for flags
  end if

  ! initialise used modules
  call qcp_startup
  call simprep_startup
end subroutine startup

!-----------------------------------------------------------------------

subroutine qpi_shutdown
  integer						:: i

  if (nodeid .eq. 0) then
	write(*,*) 'QPI version ', trim(QPI_VERSION), trim(QPI_SUFFIX), ' terminated normally.'
	write (*,'(79a)') ('#',i=1,79)
  end if

  ! call shutdown
  call qcp_shutdown
end subroutine qpi_shutdown

!-----------------------------------------------------------------------
INTEGER(4) FUNCTION qsignal( signum, proc, sigflag )
!       use MD
       implicit none
       INTEGER(4)                               :: signum, sigflag
       external proc
       qsignal = 1
END FUNCTION qsignal

end program Qpi5

! signal handlers

INTEGER(4) FUNCTION sigint_handler(sig_num)
  use QCP
  implicit none
  INTEGER(4)					:: sig_num

  call die('user request (control-C)')
  sigint_handler = 1
END FUNCTION sigint_handler

INTEGER(4) FUNCTION sigkill_handler(sig_num)
  use QCP
  implicit none
  INTEGER(4)					:: sig_num

  call die('kill signal')
  sigkill_handler = 1
END FUNCTION sigkill_handler

INTEGER(4) FUNCTION sigabrt_handler(sig_num)
  use QCP
  implicit none
  INTEGER(4)					:: sig_num

  call die('kill signal')
  sigabrt_handler = 1
END FUNCTION sigabrt_handler