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
  use MPIGLOB

  implicit none
	character(*), parameter			::	MODULE_VERSION = '0.10'
	character(*), parameter			::	MODULE_DATE    = '2015-07-17'
#if defined (USE_MPI)
	character(10)					:: QDYN_SUFFIX = '_parallel'
	integer						:: qdyn_ierr
#else
  character(10)					:: QDYN_SUFFIX = ''
#endif
	integer					::	h,nfiles,ifile,nskip,framesn,iframe
	real(8)					::	Temp_static
	character(80)				::	line
	character(200)				::	infilename

	type INPUT_TYPE
		character(80)			::	filnam
		real(8),allocatable		::	lambda(:)
	end type INPUT_TYPE

	type(INPUT_TYPE),allocatable		::	traj(:)


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

	! Only initialise static data, The banner is not needed since it is interactive
	call md_startup

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

!	call prompt ('--> No. of states ')
!	!where the nstates is defined ???
!	read (*,*) nstates	!TODO Read command need control mechanism
!	write (*,2) nstates
!2	format('# Number of states                 =',i6)

	call prompt ('--> Give T & number of trajectories to skip: ')
	read (*,*) Temp_static,nskip	!TODO Read command need control mechanism
	write (*,3) Temp_static,nskip
3	format('# T         =',f7.3,/, '# Number of trajectories to skip   =',i6)

!	call prompt ('--> Topology file name ')
!	!The fep_file is defined in md.f90
!	read (*,*) top_file	!TODO Read command need control mechanism
!	write (*,4) top_file
!4	format('# Topology file         =',a40)

!	call prompt ('--> Fep file name ')
!	!The fep_file is defined in md.f90
!	read (*,*) fep_file	!TODO Read command need control mechanism
!	write (*,5) fep_file
!5	format('# Fep file         =',a20)

	call prompt ('--> Sample MD input ') !to read the restraints and other parameters used
	!The fep_file is defined in md.f90
	read (*,*) infilename	!TODO Read command need control mechanism
	write (*,6) infilename
6	format('# Sample MD input         =',a20)

	!itirate over the trajectory files names and put them in a array.
	allocate(traj(nfiles))		!TODO Check allocation status

	do ifile=1,nfiles
		allocate(traj(ifile)%lambda(nstates))	!TODO Check allocation status
		write(line,7) ifile
7		format('--> Name of file number',i4,'       & The lambda values for each state')
		call prompt(line)
		read(*,*) traj(ifile)%filnam,traj(ifile)%lambda(1:nstates)	!TODO do sanity check for lambda (sum=1) and read
		write (*,8) traj(ifile)%filnam,(traj(ifile)%lambda(h), h=1,nstates)
8		format ('Trajectory file   ',a20,'with lambda ', 7(f8.2))
	end do	!TODO check end of file to stop the memory error



	!------------------------------------------------------------------------------------------------
	!Reading various input files and setting up the parameters
	!------------------------------------------------------------------------------------------------
	if (nodeid .eq. 0) then
		! Read input data
		if(.not. initialize_2()) call die_iso('Invalid data in input file')
		! Read the topology file and assigne
		call topology
		! Read coords, solvates etc. This produces extra arrays for force, displacement etc that is not used.
		call prep_coord
		! Read fep/evb strategy
		if ( nstates > 0 ) call get_fep
		! prepare for simulation (calc. inv. mass, total charge,...)
		call prep_sim
		! The list of Q-Q atom nb interactions. taking care of bonded and exclusions
		call make_nbqqlist

	end if

#if defined (USE_MPI)
	! initialise slave nodes
	if (numnodes .gt. 1) call init_nodes
#endif

	! count non-bonded pairs to get the maximum number, then distribute them among the nodes
	call distribute_nonbonds_iso


	!------------------------------------------------------------------------------------------------
	!Reading the trajectory files and Qatom potential calculation
	!------------------------------------------------------------------------------------------------
	if (nodeid .eq. 0) then
		do ifile=1,nfiles
			!open each trajectory file
			if(.not. trj_open(traj(ifile)%filnam)) call die_iso()
			iframe = 0
			!read each trajectory file and put it in the global position array
			do while(trj_read(x))
			iframe = iframe + 1
			!position array is updated
			end do
			framesn = iframe
			call trj_close
		end do
	end if !end of node 0




















Contains
!-----------------------------------------------------------------------

subroutine distribute_nonbonds_iso
!locals
integer					:: npp, npw, nqp, nww, nqw
type(NODE_ASSIGNMENT_TYPE),allocatable  :: node_assignment(:)
real					:: avgload, old_avgload
integer					:: i, last_cgp, last_pair
integer					:: mpitype_pair_assignment, mpitype_node_assignment
integer                                 :: average_pairs,inode,icgp,sum,less_than_sum
integer                                 :: n_bonded, n_nonbonded, master_assign
real                                    :: percent
integer                                 :: master_sum
!!!!Tmp vars för allokering
integer,parameter			:: vars = 5
integer     :: mpi_batch, i_loop
integer     :: blockcnt(vars),type(vars)
integer(8)  :: disp(vars)
!!!


! count the number of nonbonded interactions and distribute them among the nodes

if (nodeid .eq. 0) then
! nice header
call centered_heading('Distribution of charge groups','-')

!Allocate node_assignment
!>>>this includes pp pw ww qw qp<<<
allocate(node_assignment(0:numnodes-1),stat=alloc_status)
call check_alloc('node_assignment')

!Allocate arrays that hold no. pairs per chargegroup.
!>>>this includes pp pw ww qw qp<<<
call allocate_nbxx_per_cgp

!Count stuff for balancing nodes and allocating nonbond arrays nbxx()
nbpp_per_cgp = 0
nbpw_per_cgp = 0
nbqp_per_cgp = 0
nbqw_per_cgp = 0
nbww_per_cgp = 0

call nbpp_count(npp, nbpp_per_cgp) !<<<   !Only for switching atoms!!!?
call nbpw_count(npw, nbpw_per_cgp) !<<<
call nbqp_count(nqp, nbqp_per_cgp)
call nbqw_count(nqw, nbqw_per_cgp) 
call nbww_count(nww, nbww_per_cgp) !<<<

!For keeping track of actual # of nonbonded pairs
totnbpp = npp   !<<<0
totnbpw = npw   !<<<0
totnbww = nww*9 !<<<0
totnbqp = nqp
totnbqw = nqw*3*nqat

if (numnodes .eq. 1) then
! only one node: no load balancing

! make the master node handle everything
calculation_assignment%pp%start = 1
calculation_assignment%pp%end = ncgp_solute
calculation_assignment%pw%start = 1
calculation_assignment%pw%end = ncgp_solute
calculation_assignment%qp%start = 1
calculation_assignment%qp%end = ncgp_solute
calculation_assignment%qw%start = 1
calculation_assignment%qw%end = nwat
calculation_assignment%ww%start = 1
calculation_assignment%ww%end = nwat

#if defined (USE_MPI)
else ! i.e. slave nodes exists


! A simple solution to avoid parallelising the bonded
! Calculate n_bonded and n_nonbonded 
! Approximate time of computing one bonded with one nonbonded
! The number of qq-interactions are neglected
n_bonded = nbonds + nangles + ntors + nimps
n_nonbonded = totnbpp + totnbpw + totnbww + totnbqw + totnbqp !<<<

! Compare to determine how many nonbonded master should get
! A bonded is faster, so this favours an early completion for master
master_assign =  n_nonbonded/numnodes - n_bonded * numnodes

! calculate the assignments ********************
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!Calculate balanced assignment for p-p pairs
icgp=0
sum=0
 !First assign master a small part
node_assignment(0)%pp%start=icgp+1 !<<<0
percent=REAL(totnbpp)/n_nonbonded
less_than_sum = master_assign*percent  ! No. of pp-type to assign master
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpp_per_cgp(icgp)
end do
node_assignment(0)%pp%end=icgp !<<<0
master_sum=sum
 !Now assign slaves
average_pairs=(totnbpp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pp%start=icgp+1 !<<<0
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum) 
     icgp=icgp+1
     sum=sum + nbpp_per_cgp(icgp)
  end do
  node_assignment(inode)%pp%end=icgp !<<<0
end do
node_assignment(numnodes-1)%pp%start=icgp+1 !<<<0
node_assignment(numnodes-1)%pp%end=ncgp_solute !<<<0

!Calculate balanced assignment for p-w pairs
icgp=0
sum=0
node_assignment(0)%pw%start=icgp+1 !<<<0
percent=REAL(totnbpw)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpw_per_cgp(icgp)
end do
node_assignment(0)%pw%end=icgp !<<<0
master_sum=sum
average_pairs=(totnbpw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pw%start=icgp+1 !<<<0
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbpw_per_cgp(icgp)
  end do
  node_assignment(inode)%pw%end=icgp !<<<0
end do
node_assignment(numnodes-1)%pw%start=icgp+1 !<<<0
node_assignment(numnodes-1)%pw%end=ncgp_solute !<<<0
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!Calculate balanced assignment for q-p pairs
icgp=0
sum=0
node_assignment(0)%qp%start=icgp+1
percent=REAL(totnbqp)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbqp_per_cgp(icgp)
end do
node_assignment(0)%qp%end=icgp
master_sum=sum
average_pairs=(totnbqp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%qp%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbqp_per_cgp(icgp)
  end do
  node_assignment(inode)%qp%end=icgp
end do
node_assignment(numnodes-1)%qp%start=icgp+1
node_assignment(numnodes-1)%qp%end=ncgp_solute
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!Calculate balanced assignment for w-w pairs
icgp=0
sum=0
node_assignment(0)%ww%start=icgp+1 !<<<0
percent=REAL(totnbww)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbww_per_cgp(icgp)
end do
node_assignment(0)%ww%end=icgp !<<<0
master_sum=sum
average_pairs=(totnbww-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%ww%start=icgp+1 !<<<0
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbww_per_cgp(icgp)
  end do
  node_assignment(inode)%ww%end=icgp !<<<0
end do
node_assignment(numnodes-1)%ww%start=icgp+1 !<<<0
node_assignment(numnodes-1)%ww%end=nwat !<<<0
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!Calculate balanced assignment for q-w pairs
icgp=0
sum=0
node_assignment(0)%qw%start=icgp+1
percent=REAL(totnbqw)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbqw_per_cgp(icgp)
end do
node_assignment(0)%qw%end=icgp
master_sum=sum
average_pairs=(totnbqw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%qw%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbqw_per_cgp(icgp)
  end do
  node_assignment(inode)%qw%end=icgp
end do
node_assignment(numnodes-1)%qw%start=icgp+1
node_assignment(numnodes-1)%qw%end=nwat

#endif
end if    !if (numnodes .eq. 1)

! deallocate bookkeeping arrays
!deallocate(nppcgp, npwcgp, nqpcgp, nwwmol)

end if   !if (nodeid .eq. 0)

! distribute assignments to the nodes
#if defined (USE_MPI)
if (numnodes .gt. 1) then
    if (nodeid .ne. 0) then
	! Dummy allocation to avoid runtime errors when using pointer checking
	allocate(node_assignment(1),stat=alloc_status)
    endif 
! register data types
call MPI_Type_contiguous(3, MPI_INTEGER, mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')
call MPI_Type_commit(mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')

call MPI_Type_contiguous(5, mpitype_pair_assignment, mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')
call MPI_Type_commit(mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')

! distribute
call MPI_Scatter(node_assignment, 1, mpitype_node_assignment, &
    calculation_assignment, 1, mpitype_node_assignment, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('failure while sending node assignments')

! free data type
call MPI_Type_free(mpitype_node_assignment, ierr)
call MPI_Type_free(mpitype_pair_assignment, ierr)
    if (nodeid .ne. 0) then
	deallocate(node_assignment)
    endif 
end if
#endif

if (nodeid .eq. 0) then
 ! print a status report
 write(*,98) 'solute-solute', 'solute-water', 'water-water', 'Q-solute', 'Q-water'
 write(*,99) 'total', ncgp_solute,ncgp_solute,nwat,ncgp_solute,nwat
if (numnodes .gt. 1) then
 do i=0,numnodes-1
    write(*,100) i, 'assigned cgps', &
         node_assignment(i)%pp%end-node_assignment(i)%pp%start+1, &
         node_assignment(i)%pw%end-node_assignment(i)%pw%start+1, &
         node_assignment(i)%ww%end-node_assignment(i)%ww%start+1, &
         node_assignment(i)%qp%end-node_assignment(i)%qp%start+1, &
         node_assignment(i)%qw%end-node_assignment(i)%qw%start+1
 end do
end if
end if

#if defined (USE_MPI)
blockcnt(:) = 1
type(:) = MPI_INTEGER
call MPI_Bcast(totnbpp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbpw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbqp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbqw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

! allocate
calculation_assignment%pp%max = totnbpp/numnodes + 0.20*totnbpp
allocate(nbpp(calculation_assignment%pp%max), stat=alloc_status)
call check_alloc('solute-solute non-bond list')

calculation_assignment%pw%max = (totnbpw+3)/numnodes + 0.20*totnbpw
allocate(nbpw(calculation_assignment%pw%max), stat=alloc_status)
call check_alloc('solute-solvent non-bond list')

calculation_assignment%ww%max = (totnbww+nwat)/numnodes + 0.20*totnbww
allocate(nbww(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('solvent-solvent non-bond list')

calculation_assignment%qp%max = totnbqp/numnodes + 0.20*totnbqp
allocate(nbqp(calculation_assignment%qp%max), stat=alloc_status)
call check_alloc('Qatom - solute non-bond list')

calculation_assignment%qw%max = nwat
allocate(nbqw(calculation_assignment%qw%max), stat=alloc_status)
call check_alloc('Qatom - water non-bond list')



  if (use_PBC) then
	!allocate array to keep track of chargegroups
	!approximate with one half of the number of atompairs
	allocate(nbpp_cgp(calculation_assignment%pp%max / 2), stat=alloc_status)
	call check_alloc('solute-solute non-bond charge group pair list')
	allocate(nbpw_cgp(calculation_assignment%pw%max / 2), stat=alloc_status)
	call check_alloc('solute-solvent non-bond charge group pair list')
	allocate(nbqp_cgp(calculation_assignment%qp%max / 2), stat=alloc_status)
	call check_alloc('qatom-solute non-bond charge group pair list')
  end if	

!Kanske deallokera nbxx_per_cgp TODO

98 format('node value ',5a13)
99 format(a10,1x,5(1x,i12))
!99 format(a4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)
100 format(i4,1x,a5,1x,5(1x,i12))
!100 format(i4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)

if (nodeid .eq. 0)  call centered_heading('End of distribution', '-')

end subroutine distribute_nonbonds_iso

!-----------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------------------------
subroutine die_iso(cause)
!TODO deallocation of variable defined in this module

	! args
	character(*), optional				::	cause
	! local vars
	integer						::	i
	! flush stuff
	integer(4), parameter				::	stdout_unit = 6
	! external flush disabled for gfortran
	! external flush

if (nodeid .eq. 0) then
        write(*,*)
        call centered_heading('ABNORMAL TERMINATION', '!')
!	write final energies if run has started
!	if (istep > 0) then
!		if ( mod(istep,iout_cycle) .ne. 1 ) call write_out
!	end if
!	if(allocated(v)) then
!		!save restart file for diagnosing coordinate problems
!		write(*,*) 'restart file written at step', istep
!		call write_xfin
!	endif
        write (*,'(79a)') ('!',i=1,79)
        call close_output_files

        ! apologise
        write(*,'(a)') 'ABNORMAL TERMINATION of Qiso'
        if (present(cause)) then
                write(*,'(79a)') 'Terminating due to ', cause
                endif
        write (*,'(79a)') ('!',i=1,79)
#if defined(CRAY)
        !Cray can't flush stdout...
#elif defined(NO_FLUSH)
		!When you can't flush
#else
        ! flush stdout
        call flush(stdout_unit)
#endif
end if	
! clean up
call md_deallocate


#if defined (USE_MPI)
! abort all processes with exit code 255
call MPI_Abort(MPI_COMM_WORLD, 255, ierr)
#else
! stop with a message to stderr
stop 'Qiso terminated abnormally'
#endif

end subroutine die_iso
!----------------------------------------------------------------------------------------------------------------------
logical function initialize_2()
	! local variables
	character					::	text*80
	integer						::	i,j,length
	real(8)						::	stepsize
	real(8)						::	lamda_tmp(max_states)
	integer						::	fu, fstat
	real(8)						::	rjunk
	integer						::	ijunk

! local parameters
	integer						::	num_args
!	character(200)					::	infilename
	logical						::	yes
	logical						::	need_restart
	character(len=80)				::	instring
	logical						::	inlog
	integer						::	mask_rows

! this subroutine will init:
!  nsteps, stepsize, dt
!  Temp0, tau_T, iseed, Tmaxw
!  use_LRF, NBcycle, Rcpp, Rcww, Rcpw, Rcq
!  shake_solute, shake_solvent, shake_hydrogens
! fk_pshell
!  fk_wsphere=-1, wpol_restr, wpol_born
!  fkwpol=-1, Dwmz=-1 (values  ized to -1 will
!    be set in water_sphere, once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle, [trj_file], [ene_file]
!  fep_file
!  nstates, EQ (allocating memory for EQ)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)

! external definition of iargc disabled for gfortran
!integer(4) iargc
!external iargc

! read name of input file from the command line
!num_args = command_argument_count()
!if (num_args .lt. 1) call die('no input file specified on the command line')
!#if defined(CRAY)
!call pxfgetarg(num_args, infilename, 200, i)
!#elif defined(MPICH)
!call getarg(1, infilename)
!#else
!call getarg(num_args, infilename)
!#endif
text = 'Reading input from '//infilename
call centered_heading(trim(text), '-')

initialize_2 = .true. 

if(.not. prm_open_section('PBC', infilename)) then
        box = .false.
        write(*,'(a)') 'Boundary: sphere'
else
        box = .true.
        write(*,'(a)') 'Boundary: periodic box'
        if( .not. prm_get_logical_by_key('rigid_box_centre', rigid_box_centre, .false. ) ) then
                write(*,'(a)') '>>> Error: rigid_box_centre must be on or off'
                initialize_2 = .false.
        end if
        write(*,'(a,a3)') 'Rigid box centre ', onoff(rigid_box_centre)
        if( .not. prm_get_logical_by_key('constant_pressure', constant_pressure, .false.) ) then
                write(*,'(a)') '>>> Error: constant_pressure must be on or off'
                initialize_2 = .false.
        end if

        if( constant_pressure ) then
                write(*,'(a)') 'NPT-ensemble'
                volume_try = 0
                volume_acc = 0
                if( .not. prm_get_real8_by_key('max_volume_displ', max_vol_displ) ) then
                        initialize_2 = .false.
                        write(*,'(a)') '>>> ERROR: maximum volume displacement not specified (section PBC)'
                else
                        write(*,5) max_vol_displ
                end if
5	format ('Maximum volume displacemet = ', f10.3)

                if( .not. prm_get_integer_by_key('pressure_seed', pressure_seed)) then
                        pressure_seed = 3781
                end if

				write(*, '(a, i4 )' ) 'Pressure seed: ', pressure_seed

                if( .not. prm_get_real8_by_key('pressure', pressure) ) then
                        pressure = 1.0
                end if
                write(*,9) pressure
9	format ('Pressure = ',f10.3,'  bar')
                !convert pressure to strange internal unit
                pressure = pressure * 1.43836e-5
                yes = prm_get_logical_by_key('atom_based_scaling', atom_based_scaling, .false.)
                if (atom_based_scaling) then
                        write (*,'(a)') 'Coordinate scaling on volume changes:  Atom based'
                else
                        write (*,'(a)') 'Coordinate scaling on volume changes:  Molecule based'
                end if

        else
                write(*,'(a)') 'NVT-ensemble'
                if( prm_get_line_by_key('control_box', instring) ) then
                        read(instring, *) new_boxl(:)
                        control_box = .true.
                        write(*,'(a, 3f10.3)')'Boxsize will be changed to: ', new_boxl
                else
                        control_box = .false.
                end if


        end if !section constant_pressure

	yes = prm_get_logical_by_key('put_solvent_back_in_box', put_solvent_back_in_box)

	yes = prm_get_logical_by_key('put_solute_back_in_box', put_solute_back_in_box)


	if(put_solute_back_in_box .and. put_solvent_back_in_box) then
           write(*,'(a)') 'Solute and solvent molecules will be put back in box.'
	else
		if (put_solute_back_in_box) then
           	write(*,'(a)') 'Only solute molecules will be put back in box.'
		else
			if (put_solvent_back_in_box) then
				write(*,'(a)') 'Only solvent molecules will be put back in box.'
			else
				write(*,'(a)') 'No molecules will be put back in box.'				
			end if
		end if
	end if



end if !section PBC


if(.not. prm_open_section('md')) then
        call prm_close
        ! open input file
        fu = freefile()
        open(unit=fu, file=infilename, action='read', form='formatted', status='old', iostat=fstat)
        if (fstat .ne. 0) call die('error opening input file '//infilename)
                initialize_2 = old_initialize(fu)
                close(fu)
        return
end if

need_restart = .false. !flag for restart file required
if(.not. prm_get_integer_by_key('steps', nsteps)) then
        write(*,*) '>>> ERROR: steps not specified (section MD)'
        initialize_2 = .false.
end if
if(.not. prm_get_real8_by_key('stepsize', stepsize)) then
        write(*,*) '>>> ERROR: stepsize not specified (section MD)'
        initialize_2 = .false.
end if
write (*,10) nsteps, stepsize
10	format ('Number of MD steps =',i10,'  Stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462*stepsize

! --- Temperature etc.
if(.not. prm_get_real8_by_key('temperature', Temp0)) then
        write(*,*) '>>> ERROR: temperature not specified (section MD)'
        initialize_2 = .false.
end if
if(.not. prm_get_real8_by_key('bath_coupling', tau_T)) then
        write(*,*) 'Temperature bath relaxation time tau_T set to default'
        tau_T = tau_T_default
end if
write (*,15) Temp0,tau_T
tau_T=0.020462*tau_T
if(Temp0 <= 0) then
        write(*,'(a)') &
                '>>> Error: No dynamics at zero temperature!'
        initialize_2 = .false.
end if
if(tau_T < dt) then
        write(*,'(a)') '>>> Error: tau_t must be >= stepsize.'
        initialize_2 = .false.
end if

yes = prm_get_logical_by_key('separate_scaling', separate_scaling, .true.)
if(separate_scaling) then
	   write(*,'(a)') 'Solute and solvent atoms coupled separately to heat bath.'
else 
	   write(*,'(a)') 'Solute and solvent atoms coupled together to heat bath.'
end if

15	format ('Target temperature =',f10.2,'  T-relax time     =',f10.2)

yes = prm_get_integer_by_key('random_seed', iseed, 1) 
if(.not. prm_get_real8_by_key('initial_temperature', Tmaxw)) then
        iseed = 100 !set iseed = 0 if no initial temp
        need_restart = .false. !never check for restart
end if

if (iseed > 0) write (*,16) Tmaxw, iseed
16	format ('Initial velocities will be generated from Maxwell distribution:',&
                /,'Maxwell temperature=',f10.2,' Random number seed=',i10)

! --- shake, LRF
if(.not. prm_get_logical_by_key('shake_solvent', shake_solvent, .true.)) then
        write(*,'(a)') '>>> Error: shake_solvent must be on or off.'
        initialize_2 = .false.
end if
write(*,17) 'all solvent bonds', onoff(shake_solvent)
17	format('Shake constaints on ',a,t42,': ',a3)

if(.not. prm_get_logical_by_key('shake_solute', shake_solute, .false.)) then
        write(*,'(a)') '>>> Error: shake_solute must be on or off.'
        initialize_2 = .false.
end if 
write(*,17) 'all solute bonds', onoff(shake_solute)

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .false.)) then
        write(*,'(a)') '>>> Error: shake_hydrogens must be on or off.'
        initialize_2 = .false.
end if 
write(*,17) 'all bonds to hydrogen', onoff(shake_hydrogens)


       yes = prm_get_logical_by_key('lrf', use_LRF, .true.)
       if(use_LRF) then
               write(*,20) 'LRF Taylor expansion outside cut-off'
       else 
               write(*,20) 'standard cut-off'
       end if

20	format ('Nonbonded method   = ',a)

yes = prm_get_logical_by_key('force_rms', force_rms, .false.)
if(force_rms) then
        write(*,22) 
end if
22	format ('R.M.S. force will be calculated.')


! --- Rcpp, Rcww, Rcpw, Rcq, RcLRF
if(.not. prm_open_section('cut-offs')) then
        write(*,'(a)') 'No cut-offs section, default cut-offs used'
        rcpp = rcpp_default
        rcww = rcww_default
        rcpw = rcpw_default
        rcq = rcq_default
        rcLRF = rcLRF_default
else
        if(.not. prm_get_real8_by_key('solute_solute', rcpp, rcpp_default)) then
                write(*,'(a)') 'solute-solute cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('solvent_solvent', rcww, rcww_default)) then
                write(*,'(a)') 'solvent-solvent cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('solute_solvent', rcpw, rcpw_default)) then
                write(*,'(a)') 'solute-solvent cut-off set to default'
        end if
        if(.not. prm_get_real8_by_key('q_atom', rcq, rcq_default)) then
                write(*,'(a)') 'q-atom cut-off set to default'
        end if
        if(use_LRF) then
                if(.not. prm_get_real8_by_key('lrf', rcLRF, rcLRF_default)) then
                        write(*,'(a)') 'LRF cut-off set to default'
                end if
                if(RcLRF < rcpp .or. RcLRF < rcpw .or. RcLRF < rcww) then
                        write(*,'(a)') &
                                '>>> ERROR; LRF cut-off must not be smaller than solute or solvent cut-offs!'
                        initialize_2 = .false.
                end if
        end if
end if

write (*,25) Rcpp,Rcww,Rcpw,Rcq
if(use_LRF) write(*,26) RcLRF
25	format ('Cut-off radii for non-bonded interactions:',/, &
                'Solute-solute:    ',f6.2,/,&
                'Solvent-solvent:  ',f6.2,/,&
                'Solute-solvent:   ',f6.2,/,&
                'Q-atom-non-Q-atom:',f6.2)
26	format ('LRF:              ',f6.2)

30	format ('>>> WARNING: Ingnoring obsolete keyword ',a,'.')
! --- simulation sphere

if( .not. box ) then
        if(.not. prm_open_section('sphere')) then
			fk_pshell = fk_pshell_default
			print*,'Radius of inner restrained shell set to 85% of exclusion shell radius.'
			rexcl_i = shell_default
			write(*,50) rexcl_i
        else
                if(prm_get_line_by_key('centre', instring)) then
                        write(*,30) 'centre'
                end if
                ! --- rexcl_o, rexcl_i, fk_pshell
                if(prm_get_real8_by_key('radius', rjunk)) then
                        write(*,30) 'radius'
                end if
                if(prm_get_real8_by_key('shell_radius', rexcl_i)) then  !inner radius of restrained shell
                    write(*,50) rexcl_i
                    if(rexcl_i < 0.) then
                      call die('inner radius of restrained shell must be >= 0')
                    end if
                else
                    print*,'Radius of inner restrained shell set to 85% of exclusion shell radius.'
					rexcl_i = shell_default
                    write(*,50) rexcl_i
                end if
50 format('Radius of inner restrained shell       =    ',f8.3) 
                if(.not. prm_get_real8_by_key('shell_force', fk_pshell)) then
                        write(*,'(a)') 'Shell force constant set to default'
                        fk_pshell = fk_pshell_default
                end if
                if(fk_pshell > 0) then
                        write(*,47) fk_pshell
                end if
47		format('Shell restraint force constant         =',f8.2)
                if(.not. prm_get_real8_by_key('excluded_force', fk_fix   )) then
                        write(*,'(a)') 'Excluded atoms force constant set to default'
                        fk_fix = fk_fix_default
                end if
                if(fk_fix > 0) then
                        write(*,48) fk_fix
                else
                        fk_fix = fk_fix_default
                        write(*,'(a)')'Shell restraint force constant can not be less-equal zero, set to default'
                        write(*,48) fk_fix
                end if
48              format('Excluded atom restraint force constant         =',f8.2)

                yes = prm_get_logical_by_key('excluded_freeze', freeze, .false.)
                if(freeze) then
                        write(*,'(a)') &
                                'Excluded atoms will not move.'
                end if

                yes = prm_get_logical_by_key('exclude_bonded', exclude_bonded, .false.)
                if(exclude_bonded) then
                        write(*,'(a)') &
                                'Bonded interactions outside the sphere will be eliminated'
                end if
       end if

        ! --- solvent 
        inlog = prm_open_section('solvent')
        if(.not. inlog) inlog = prm_open_section('water') !try also the old name
        if(.not. inlog) then       !defaults
                fk_wsphere = -1
                Dwmz = -1
                awmz = -1
                wpol_restr = wpol_restr_default
                wpol_born = wpol_restr_default
                fkwpol = -1 
        else
                if(prm_get_real8_by_key('radius', rwat_in)) then
                        write(*,'(a,f8.2)') 'Target solvent radius =',rwat_in
                end if
                if(prm_get_line_by_key('centre', instring)) then
                        write(*,30) 'centre'
                end if
                if(prm_get_real8_by_key('pack', rjunk)) then
                        write(*,30) 'pack'
                end if


          if(.not. prm_get_real8_by_key('radial_force', fk_wsphere)) then
                write(*,'(a)') 'Solvent radial restraint force constant set to default'
                fk_wsphere = -1 ! this will be set in water_sphere, once target radius is known
          end if
          yes=prm_get_logical_by_key('polarisation', wpol_restr, wpol_restr_default)
          !default is on when pol. restr is on, otherwise off
          yes=prm_get_logical_by_key('charge_correction', wpol_born, wpol_restr)
          if(wpol_born .and. .not. wpol_restr) then
                write(*,'(a)') '>>> ERROR: charge_correction on requires polarisation on (section solvent)'
                initialize_2 = .false.
          end if
          if(.not. prm_get_real8_by_key('polarisation_force', fkwpol)) then
                write(*,'(a)') 'Solvent polarisation force constant set to default'
                fkwpol = -1 ! this will be set in water_sphere, once target radius is known
          end if
          yes = prm_get_real8_by_key('morse_depth', Dwmz, -1._8)			
          yes = prm_get_real8_by_key('morse_width', awmz, -1._8)			
          if(prm_get_string_by_key('model', instring)) then
                write(*,30) 'model'
          end if
        end if !if (.not. inlog)
end if !if( .not. box )


if(.not. prm_open_section('intervals')) then
        write(*,'(a)') 'non-bond list update interval set to default.'
        NBcycle = NB_cycle_default
        write(*,'(a)') 'energy summary interval set to default.'
        iout_cycle = iout_cycle_default
        itemp_cycle = iout_cycle_default
        iene_cycle = 0 !no energy
        itrj_cycle = 0 !no trajectory

        ivolume_cycle = ivolume_cycle_default


else
        if(.not. prm_get_integer_by_key('non_bond', NBcycle)) then
                write(*,'(a)') 'non-bond list update interval set to default.'
                NBcycle = NB_cycle_default
        end if
        if(.not. prm_get_integer_by_key('output', iout_cycle)) then
                write(*,'(a)') 'energy summary interval set to default.'
                iout_cycle = iout_cycle_default
        end if
        if(.not. prm_get_integer_by_key('temperature', itemp_cycle)) then
                write(*,'(a)') 'temperature print-out interval set to default.'
                itemp_cycle = iout_cycle_default
        end if
        yes = prm_get_integer_by_key('energy', iene_cycle, 0)
        yes = prm_get_integer_by_key('trajectory', itrj_cycle, 0)

        if( constant_pressure ) then
                if( .not. prm_get_integer_by_key('volume_change', ivolume_cycle) ) then
                        write(*,'(a)') 'volume change intervall set to default'
                        ivolume_cycle = ivolume_cycle_default
                end if
        end if
end if

write(*,84) NBcycle
84	format('Non-bonded pair list update interval   =',i8)
86	format('Energy summary print-out interval      =',i8)
87	format('Temperature print-out interval         =',i8)
88	format('Trajectory write interval              =',i8)
89	format('Energy file write interval             =',i8)
83  format('Volume change interval                 =',i8)

if(iout_cycle > 0) then
        write (*,86) iout_cycle
else
        write(*,'(a)') 'No energy summaries written.'
        iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if
if(itemp_cycle > 0) then
        write (*,87) itemp_cycle
else
        write(*,'(a)') 'No temperatures written.'
        itemp_cycle = -999999999 ! make sure mod(istep, itemp_cycle) never = 0
end if
if(itrj_cycle > 0) then
        write (*,88) itrj_cycle
else
        itrj_cycle = -999999999 !no energy
        write(*,'(a)') 'No trajectory written.'
end if
if(iene_cycle > 0) then
        write (*,89) iene_cycle
else
        iene_cycle = -999999999 !no energy
        write(*,'(a)') 'No energy file written.'
end if
if( constant_pressure ) then
        write(*,83) ivolume_cycle
end if

!read trajectory atom mask
mask_rows = prm_count('trajectory_atoms')
if(itrj_cycle > 0) then
        if(mask_rows == 0) then
                write(*,'(a)') 'All atoms will be included in the trajectory.'
                yes = trj_store_mask('all')
        else
                do i=1,mask_rows
                        yes = prm_get_line(text)
                        yes = trj_store_mask(text)
                end do
        end if
elseif(mask_rows == 0) then
        write(*,'(a)') 'Ignoring section trajectory_atoms.'
end if

if(.not. prm_open_section('files')) then
        write(*,'(a)') '>>> ERROR: files section not found.'
        initialize_2 = .false.
else
        if(.not. prm_get_string_by_key('topology', top_file)) then
                write(*,'(a)') '>>> ERROR: topology not specified (section files)'
                initialize_2 = .false.
        end if
        write (*,60) trim(top_file)
60		format ('Topology file      = ',a)

        if(.not. prm_get_string_by_key('restart', restart_file)) then
                restart = .false.
                if(need_restart) then
                        write(*,'(a)') '>>> ERROR: Restart file required when initial temp. not given.'
                        initialize_2 = .false.
                end if
        else
                restart = .false.   !never check for restart
        end if

        if(restart) then
                write (*,65) trim(restart_file)
        else
                write (*,'(a)') 'Initial coordinates taken from topology.'
                if(iseed == 0) then
                        write(*,'(a)') &
                                '>>> ERROR: Need a random number seed to generate initial velocities, aborting.'
                        initialize_2 = .false.
                end if
        end if
65		format ('Initial coord. file= ',a)

        if(.not. prm_get_string_by_key('final', xfin_file)) then
                write(*,'(a)') '>>> ERROR: final co-ordinate file not specified (section files, keyword final)'
                initialize_2 = .false.
        end if
        write (*,80) trim(xfin_file)
80		format ('Final coord. file  = ',a)

        if(.not. prm_get_string_by_key('trajectory', trj_file)) then
                if(itrj_cycle > 0) then
                        write(*,'(a)') '>>> ERROR: Trajectory file name required to write trajectory!'
                        initialize_2 = .false.
                end if
        else
                if(itrj_cycle < 0) then
                        write(*,*) '>>> Error: Trajectory file given but no output interval'
                        initialize_2 = .false.
                end if
                if(itrj_cycle > 0) write (*,90) trim(trj_file)
        end if
90		format ('Trajectory file    = ',a)

        if(.not. prm_get_string_by_key('energy', ene_file)) then
                if(iene_cycle > 0) then
                        write(*,'(a)') '>>> ERROR: Energy file name required to write energies!'
                        initialize_2 = .false.
                end if
        else

                if(iene_cycle < 0) then

                        write(*,'(a)') '>>> ERROR: Energy file given but no energy interval'

                        initialize_2 = .false.

                end if
                if(iene_cycle > 0) write (*,94) trim(ene_file)
        end if
94		format ('Energy output file = ',a)

        if(.not. prm_get_string_by_key('fep', fep_file)) then
                write(*,'(a)') 'No FEP file.'
                !initialize_2 = .false. !This condition IS OK.
                fep_file = ''
        else
                write (*,95) trim(fep_file)
95			format ('FEP input file     = ',a,/)
        end if
        if(.not. prm_get_string_by_key('restraint', exrstr_file)) then
                implicit_rstr_from_file = 0
        else
                implicit_rstr_from_file = 1
                write (*,104) trim(exrstr_file)
104			format ('External rstr file = ',a,/)
        end if
        if(prm_get_string_by_key('water', instring)) then
                write(*,30) 'water'
        end if
end if			

! --- states, EQ
nstates = 0
if(prm_open_section('lambdas')) then
        do while(prm_get_field(instring))
                nstates = nstates + 1
                read(instring, *, iostat=fstat) lamda_tmp(nstates)
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid lambda value.'
                        initialize_2 = .false.
                        exit
                end if
        end do
end if
if(nstates == 0 .and. fep_file /= '') then
        if(fep_file /= '') then
                write(*,'(a)') 'Defaulting to single FEP state.'
                nstates = 1
                lamda_tmp(1) = 1.
        end if
end if
if(nstates > 0 ) then
        if(fep_file == '') then
                write(*,'(a)') '>>> ERROR: FEP file required to use lambdas!'
                initialize_2 = .false.
        else
                ! allocate memory for EQ
                allocate(EQ(nstates), stat=alloc_status)
                call check_alloc('Q-atom energy array')

                ! init EQ%lambda
                EQ(1:nstates)%lambda = lamda_tmp(1:nstates)
                write (*,98) (EQ(i)%lambda,i=1,nstates)
98			format ('lambda-values      = ',10f8.5)
        end if
end if

!	--- restraints:
write (*,'(/,a)') 'Listing of restraining data:'

! --- nrstr_seq, [rstseq]
nrstr_seq = prm_count('sequence_restraints')
109 format (/,'No. of sequence restraints =',i10)
if ( nrstr_seq .gt. 0 ) then
        ! allocate memory for rstseq
        write (*,109) nrstr_seq
        allocate(rstseq(nrstr_seq), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,110)
110		format ('  atom_i  atom_j      fc  H-flag to_centre')
        do i=1,nrstr_seq
                ! read rstseq(i)
                yes = prm_get_line(text)
                rstseq(i)%to_centre = 0 
                read(text,*, end=111, err=111) rstseq(i)
111			write(*,112) rstseq(i)
112			format (2i8,f8.2,i8,i10)
  end do
end if

! --- nrstr_pos, [rstpos]
nrstr_pos = prm_count('atom_restraints')
115 format (/,'No. of position restratints =',i10)
if ( nrstr_pos .gt. 0 ) then
        write (*,115) nrstr_pos
        ! allocate memory for rstpos
        allocate(rstpos(nrstr_pos), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,120)
120		format ('atom_i      x0      y0      z0     fcx     fcy     fcz   state')
        do i=1,nrstr_pos ! read rstpos(i)
                yes = prm_get_line(text)
                read(text,*, iostat=fstat) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
                        (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid atom restraint data.'
                        initialize_2 = .false.
                        exit
                end if
                write (*,122) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
                        (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
        end do
122		format (i6,6f8.2,i8)
end if

! --- nrstr_dist, [rstdis]
nrstr_dist = prm_count('distance_restraints')
125	format (/,'No. of distance restraints =',i10)
if ( nrstr_dist .gt. 0 ) then
        write (*,125) nrstr_dist
        ! allocate memory for rstdis
        allocate(rstdis(nrstr_dist), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,130)
130		format ('atom_i atom_j   dist1   dist2   fc        state') 
        do i=1,nrstr_dist
                yes=prm_get_line(text)
                ! read rstdis(i)
                if(scan(text, ':') > 0) then !got res:atnr
                  !Store in i&j as res:atnr and assign atom nr after topology is read (prep_coord)
                  read(text,*, iostat=fstat) rstdis(i)%itext,rstdis(i)%jtext,rstdis(i)%d1,& 
                     rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
                else !Plain numbers
                  read(text,*, iostat=fstat) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,&
                        rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
                  rstdis(i)%itext = 'nil'
                  rstdis(i)%jtext = 'nil'
                end if
                if(fstat /= 0) then
                  write(*,'(a)') '>>> ERROR: Invalid distance restraint data.'
                  initialize_2 = .false.
                  exit
                end if
                write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%d2,rstdis(i)%fk, &
                        rstdis(i)%ipsi
        end do
132		format (i6,1x,i6,3f8.2,i8)
end if

! --- nrstr_angl, [rstang]
nrstr_angl = prm_count('angle_restraints')
135     format (/,'No. of angle restraints =',i10)
if ( nrstr_angl .gt. 0 ) then
        write (*,135) nrstr_angl
        ! allocate memory for rstang
        allocate(rstang(nrstr_angl), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,140)
140             format ('atom_i atom_j atom_k   angle   fc        state')
        do i=1,nrstr_angl
                yes=prm_get_line(text)
                ! read rstang(i)
                !if(scan(text, ':') > 0) then !got res:atnr
                  !Store in i&j as res:atnr and assign atom nr after topology is
                  !read (prep_coord)
                !  read(text,*, iostat=fstat) rstang(i)%itext,rstang(i)%jtext,rstang(i)%ktext,&
                !     rstang(i)%ang, rstang(i)%fk, rstang(i)%ipsi
                !else !Plain numbers
                  read(text,*, iostat=fstat) rstang(i)%i,rstang(i)%j,rstang(i)%k,&
                        rstang(i)%ang, rstang(i)%fk, rstang(i)%ipsi
                 ! rstang(i)%itext = 'nil'
                 ! rstang(i)%jtext = 'nil'
                 ! rstang(i)%ktext = 'nil'
                !end if
                if(fstat /= 0) then
                  write(*,'(a)') '>>> ERROR: Invalid angle restraint data.'
                  initialize_2 = .false.
                  exit
                end if
                write (*,142) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang,rstang(i)%fk, &
                        rstang(i)%ipsi
        end do
142             format (i6,1x,i6,1x,i6,2f8.2,i8)
end if


if (.not. box )then
! --- nrstr_wall, [rstwal]
nrstr_wall = prm_count('wall_restraints')
145 format (/,'No. of wall sequence restraints=',i10)
if ( nrstr_wall .gt. 0) then
        write (*,145) nrstr_wall
        ! allocate memory for rstwal
        allocate(rstwal(nrstr_wall), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,150)
150  format ('atom_i atom_j   dist.      fc  aMorse  dMorse  H-flag')
        do i=1,nrstr_wall
                ! read rstwal(:)
                yes = prm_get_line(text)
                read(text,*, iostat=fstat) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
                        rstwal(i)%aMorse, rstwal(i)%dMorse, rstwal(i)%ih
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid wall restraint data.'
                        initialize_2 = .false.
                        exit
                end if
                write (*,152) rstwal(i)     
        end do
152  format (i6,1x,i6,4f8.2,i8)
end if
end if

call prm_close
end function initialize_2
!----------------------------------------------------------------------------------------------------------------------




end program Qiso
