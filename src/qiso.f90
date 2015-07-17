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
  use MPIGLOB ! use MPI global data
  use NRGY
  use PARSE


  implicit none
	character(*), parameter			::	MODULE_VERSION = '0.10'
	character(*), parameter			::	MODULE_DATE    = '2015-07-01'
	integer					::	nfiles

	!defined in MPIGLOB
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
	read (*,*) nfiles
	write (*,1) nfiles
1	format('# Number of files                 =',i6)
	call prompt ('--> No. of states ')
	!where the nstates is defined ???
	read (*,*) nstates
	write (*,2) nstates
2	format('# Number of states                 =',i6)


  
























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
