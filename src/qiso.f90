program Qiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!# Qiso caclulate the H/D KIE baesd on the Path Integral method.
!# the Qfep strategy will be adopted for input section
!# the Qcalc starategy will be used for reading the trajectory files
!# The Qdyne strategy will be used for reading topology and fep files
!  and calculating potenial energy (classic)
!# the free energy calculation has to be done internally since we need
!  the bin averge of centroids.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use MD

  implicit none
	character(*), parameter			::	MODULE_VERSION = '0.10'
	character(*), parameter			::	MODULE_DATE    = '2015-07-01'
	integer					::	h,nfiles,ifile,nskip
	real(8)					::	Temp_static
	character(80)				::	line

	type INPUT_TYPE
		character(80)			::	filnam
		real(8),allocatable		::	lambda(:)
	end type INPUT_TYPE

	type(INPUT_TYPE),allocatable		::	traj(:)


	!defined in MPIGLOB via md.90
	nodeid = 0
	numnodes = 1

	!------------------------------------------------------------------------------------------------
	! INPUT OF PARAMETERS
	! we need # of files, # of states, topology file name, fep file name, trajectory file names, 
	! # of frames to be skiped, tempreture 
	! and the information regarding the free energy calculation which could be done later 
	!the limitation of this approach is the information has to be given in a specific order
	!------------------------------------------------------------------------------------------------

	call prompt ('--> Number of trajectory files: ')
	read (*,*) nfiles	!TODO Read command need control mechanism
	write (*,1) nfiles
1	format('# Number of files                 =',i6)

	call prompt ('--> No. of states ')
	!where the nstates is defined ???
	read (*,*) nstates	!TODO Read command need control mechanism
	write (*,2) nstates
2	format('# Number of states                 =',i6)

	call prompt ('--> Give T & number of trajectories to skip: ')
	read (*,*) Temp_static,nskip	!TODO Read command need control mechanism
	write (*,3) Temp_static,nskip
3	format('# T         =',f7.3,/, '# Number of trajectories to skip   =',i6)

	call prompt ('--> Topology file name ')
	!The fep_file is defined in md.f90
	read (*,*) top_file	!TODO Read command need control mechanism
	write (*,4) top_file
4	format('# Topology file         =',a40)

	call prompt ('--> Fep file name ')
	!The fep_file is defined in md.f90
	read (*,*) fep_file	!TODO Read command need control mechanism
	write (*,5) fep_file
5	format('# Fep file         =',a20)

!itirate over the trajectory files names and put them in a array.
	allocate(traj(nfiles))		!TODO Check allocation status

	do ifile=1,nfiles
		allocate(traj(ifile)%lambda(nstates))	!TODO Check allocation status
		write(line,6) ifile
6		format('--> Name of file number',i4,'       & The lambda values for each state')
		call prompt(line)
		read(*,*) traj(ifile)%filnam,traj(ifile)%lambda(1:nstates)	!TODO do sanity chech for lambda (sum=1) and read
		write (*,7) traj(ifile)%filnam,(traj(ifile)%lambda(h), h=1,nstates)
7		format ('Trajectory file   ',a20,'with lambda ', 7(f8.2))
	end do	!TODO check end of file to stop the memory error

!Read the topology file and assigne
	call topology














Contains

subroutine prompt (outtxt)
	character(*) outtxt
#if defined (__osf__)
	!prompt to STDERR using unit 5=STDIN on OSF/1=DEC UNIX
	integer, parameter			::	f=5
	!write (f,'($,a)') outtxt
#elseif defined (_WIN32)
	!open the ERR file on Win32
	integer, save :: f
	if(f==0) then
		f=17
		open(unit=f,file='ERR')
	end if
	!write (f1,'($,a)') outtxt
#else
	!otherwise prompt to STDOUT
	integer, parameter			::	f=6
	!write (f2,'($,a)') outtxt
#endif
write (f,'(a,$)') outtxt
end subroutine prompt







end program Qiso
