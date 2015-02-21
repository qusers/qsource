!------------------------------------------------------------------------------!
!  Q V.5.7                                                                     !
!  code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Kajsa Ljunjberg, John Marelius, Martin Nervall                              !
!  Maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer,                    !
!  Mauricio Esguerra, Irek Szeler, Masoud Karemi                               !
!  latest update: February 9, 2015                                             !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  md.f90                                                                      !
!  by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg,          !
!  Martin Nervall & Martin Almlöf                                              !
!------------------------------------------------------------------------------!
module md

!used modules
!use profiling
use sizes
use trj
use mpiglob
use qatom
#if defined (_df_version_)
use dfport
use dflib
implicit none
#if defined (use_mpi)
include "mpif.h"
#endif

!-------------------------------------------------------------------------------
!       shared variables
!-------------------------------------------------------------------------------
!constants
character*(*), parameter :: md_version = '5.7'
character*(*), parameter :: md_date = '2015-02-18'
real(8) :: pi, deg2rad !set in sub startup
real, parameter :: rho_wat = 0.0335 ! molecules / a**3
real, parameter :: boltz = 0.001986

!read status
integer :: stat

!print temperature if it changed more than 2% in one time step
real, parameter :: temp_print_threshold=0.02

!memory management
integer                         :: alloc_status

!topology information
! --- atoms
integer                         :: natom
integer                         :: nat3
! --- q-atom number or 0 for non-q
integer(tiny), allocatable      :: iqatom(:)
! --- water topology
! atoms of water, total water molecules, excluded water molecules
integer                         :: nwat
real(4)                         :: crg_ow, crg_hw, mu_w

!-------------------------------------------------------------------
!       periodic box information
!-------------------------------------------------------------------

!******pwadded variable 2001-10-10
logical                         :: box, rigid_box_centre
logical                         :: put_solute_back_in_box, put_solvent_back_in_box

!variables used for constant pressure algorithm
logical                         :: constant_pressure = .false.
integer                         :: volume_try, volume_acc
real(8)                         :: pressure, max_vol_displ
integer                         :: pressure_seed
real(8), allocatable            :: mol_mass(:)
real(8), allocatable            :: mass(:)

!variables used when user changes boxsize in inputfile
logical                         :: control_box
real(8)                         :: new_boxl(3) 

!-----------------------------------------------------------------------
!       dynamics control information     
!-----------------------------------------------------------------------
! --- md parameters
integer                                         ::      nsteps, istep
integer                                         ::      iseed
logical                                         ::      restart
real(8)                                         ::      dt, temp0, tmaxw, tau_t
logical                                         ::      shake_solvent, shake_solute
logical                                         ::      shake_hydrogens
logical                                         ::      separate_scaling

! --- non-bonded strategy
logical                                         ::      use_lrf
integer                                         ::      nbcycle
real(8)                                         ::      rcpp, rcww, rcpw, rcq, rclrf
integer, parameter                              ::      max_atyp = 255
integer                                         ::      ljcod(max_atyp,max_atyp)


! --- output control parameters
integer                                         ::      itrj_cycle, iene_cycle
integer                                         ::      itemp_cycle, iout_cycle
logical                                         ::      force_rms
!******pwadded 2001-10-23
integer                                         ::      ivolume_cycle

! --- protein boundary
logical                                         ::      exclude_bonded
real(8)                                         ::      fk_pshell
real(8), parameter                              ::      fk_fix = 200.0

! --- water sphere
real(8)                                         ::      dwmz, awmz, rwat_in
real(8)                                         ::      fkwpol
logical                                         ::      wpol_restr, wpol_born
real(8)                                         ::      fk_wsphere, crgtot, crgqtot
integer(ai), allocatable                        ::      list_sh(:,:), nsort(:,:)
real(8),  allocatable                           ::      theta(:), theta0(:), tdum(:)
integer                                         ::      nwpolr_shell, n_max_insh


type shell_type
  real                                    ::      rout, dr, cstb
  real                                    ::      avtheta, avn_insh, theta_corr
  integer                                 ::      n_insh
end type shell_type

type(shell_type), allocatable :: wshell(:)

! constants & default values
integer, parameter                      ::      itdis_update            = 100
real,  parameter                        ::      wpolr_layer             = 3.0001
real, parameter                         ::      drout                   = 0.5
real(8), parameter                      ::      tau_t_default           = 10.
real(8), parameter                      ::      rcpp_default            = 10.
real(8), parameter                      ::      rcww_default            = 10.
real(8), parameter                      ::      rcpw_default            = 10.
real(8), parameter                      ::      rcq_default             = 99.
real(8), parameter                      ::      rclrf_default           = 99.
!recxl_i is set to rexcl_o * shell_default
real(8), parameter                      ::      shell_default           = 0.85
real(8), parameter                      ::      fk_pshell_default       = 10.
integer, parameter                      ::      itrj_cycle_default      = 100
integer, parameter                      ::      iene_cycle_default      = 10
integer, parameter                      ::      iout_cycle_default      = 10
integer, parameter                      ::      nb_cycle_default        = 10
real(8), parameter                      ::      fkwpol_default          = 20.
real(8), parameter                      ::      fk_wsphere_default      = 60.
logical, parameter                      ::      wpol_restr_default      = .true.
integer, parameter                      ::      ivolume_cycle_default   = 10
integer, parameter                      ::      ipressure_cycle_default = 100
!yes, and use born corr.
!constants in the sigmoid function giving default dwmz as function of radius

! --- file names
character(len=200)                      ::      top_file
character(len=200)                      ::      restart_file
character(len=200)                      ::      xfin_file
character(len=200)                      ::      trj_file
character(len=200)                      ::      fep_file
character(len=200)                      ::      ene_file
character(len=200)                      ::      exrstr_file
character(len=200)                      ::      xwat_file

! --- restraints
integer                                 ::      implicit_rstr_from_file
integer                                 ::      nrstr_seq, nrstr_pos, nrstr_dist, nrstr_wall

type rstrseq_type
  integer(ai)                           ::      i,j
  real(8)                               ::      fk
  integer(tiny)                         ::      ih
  integer                               ::      to_centre !flag for restraining to geom. or mass centre
end type rstrseq_type

type rstrpos_type
        integer(ai)                     ::      i
        integer(tiny)                   ::      ipsi
        real(8)                         ::      fk(3)
        real(8)                         ::      x(3)
end type rstrpos_type

type rstrdis_type
        integer(ai)                     ::      i,j
        integer(tiny)                   ::      ipsi
        real(8)                         ::      fk
        real(8)                         ::      d1, d2
        character(len=20)               ::  itext,jtext
end type rstrdis_type

type rstrwal_type
        integer(ai)                             ::      i,j
        real(8)                                 ::      d, fk, amorse, dmorse
        integer(tiny)                   ::      ih
end type rstrwal_type

type(rstrseq_type), allocatable::       rstseq(:)
type(rstrpos_type), allocatable::       rstpos(:)
type(rstrdis_type), allocatable::       rstdis(:)
type(rstrwal_type), allocatable::       rstwal(:)


!-----------------------------------------------------------------------
!       coordinates, velocities, forces
!-----------------------------------------------------------------------
real(8), allocatable            ::      d(:)
real(8), allocatable            ::      x(:)
real(8), allocatable            ::      xx(:) !for shake
real(8), allocatable            ::      v(:)
real, allocatable               ::      winv(:)
real(8)                         ::      grms !rms force


!-----------------------------------------------------------------------
!       energies , eq is defined in qatom.f90
!-----------------------------------------------------------------------
type(energies)                          ::      e
real(8)                                         ::      tfree, tfree_solvent, tfree_solute
real(8)                                         ::      temp_solvent, temp_solute, texcl_solute, texcl_solvent


!-----------------------------------------------------------------------
!       nonbonded pair information
!-----------------------------------------------------------------------
type nb_type
        integer(ai)                             ::      i, j
        integer(tiny)                   ::      ljcod
        integer                                 :: cgp_pair ! cgp_pair only used with periodic conditions
end type nb_type

type cgp_pair_type
        integer(ai)                             :: i, j !switching atoms (or equal in case of no switching atoms) of the chargegroups
        real(8)                                 :: x, y, z !periodical shifts
end type cgp_pair_type

type nbqp_type
        integer(ai)                             ::      i, j
        integer(tiny)                   ::      ljcod, qljcod
        integer                                 ::      cgp_pair !this variable only used with periodic conditions
end type nbqp_type


type nbq_type
        integer(ai)                             ::      j       !atom number
        integer(ai)                             ::      iq, jq  !q-atom numbers
        integer(tiny)                   ::      ljcod
        real(8)                 ::  el_scale !scale factor for electostatic interactions in qq-pairs
end type nbq_type

integer                                         ::      nbpp_pair !current no solute-solute pairs
type(nb_type), allocatable, target::nbpp(:)

integer ::      nbww_pair,nbww_true_pair !current no solvent-solvent pairs, implicit and explicit
integer(ai), allocatable, target::nbww(:)

integer                                         ::      nbpw_pair !current no solute-solvent pairs
type(nb_type), allocatable, target::nbpw(:)

integer                                         ::      nbqq_max !max number of q-q pairs in any state
integer(tiny), allocatable      ::      qconn(:,:,:) !q-atom connectivity list
integer                                         ::      nbqq_pair(max_states)
type(nbq_type), allocatable ::nbqq(:,:)

integer                                         ::      nbqp_pair !current no of qatom-solute pairs
type(nbqp_type), allocatable, target::nbqp(:)

integer                                         ::      nbqw_pair !current no of q-atom-water mol. pairs
integer(ai), allocatable        ::      nbqw(:)


!these three used only under periodic conditions
integer                                         :: nbpp_cgp_pair !number of solute-solute chargegroups interacting
type(cgp_pair_type), allocatable:: nbpp_cgp(:)
integer                                         :: nbpw_cgp_pair
type(cgp_pair_type), allocatable :: nbpw_cgp(:)
integer                                         ::      nbqp_cgp_pair
type(cgp_pair_type), allocatable :: nbqp_cgp(:)

!special monitoring of pairs
integer (tiny),allocatable  :: special_ljcod(:,:,:,:)

! lrf related variables
integer(ai), allocatable        ::      iwhich_cgp(:)


type lrf_type
        real(8)                                 ::      cgp_cent(3)
        real(8)                                 ::      phi0
        real(8)                                 ::      phi1(3)
        real(8)                                 ::      phi2(9)
        real(8)                                 ::      phi3(27)
end type lrf_type

type(lrf_type), allocatable     ::      lrf(:)
type(lrf_type), allocatable     ::      old_lrf(:)  !for constant pressure: mc_volume routine


type(node_assignment_type)      :: calculation_assignment


!shake types & variables
!convergence criterion (fraction of distance)
real(8), parameter                      ::      shake_tol = 0.0001
integer, parameter                      ::      shake_max_iter = 1000


type shake_bond_type
        integer(ai)                             ::      i,j
        real(8)                                 ::      dist2
        logical                                 ::      ready
end type shake_bond_type


type shake_mol_type
        integer                                 ::      nconstraints
        type(shake_bond_type), pointer :: bond(:)
end type shake_mol_type


integer                                         ::      shake_constraints, shake_molecules
type(shake_mol_type), allocatable :: shake_mol(:)

!-----------------------------------------------------------------------
!       profiling vars
!-----------------------------------------------------------------------
#if defined (profiling)


integer, parameter              :: num_profiling_times = 11

type profiling_var_type
        character(len=100)      :: name
        real(8)                 :: time = 0.0
end type profiling_var_type

type(profiling_var_type)        :: profile(num_profiling_times)

#if defined (use_mpi)
 !vectors for keeping track of node times
real(8),allocatable                             :: all_node_times(:)
real(8),allocatable                             :: node_times(:)
#endif

#endif

!-----------------------------------------------------------------------
!       temperature calculation variables
!-----------------------------------------------------------------------
integer                                         :: ndegf,ndegfree
integer                                         :: ndegf_solute,ndegfree_solute
integer                                         :: ndegf_solvent,ndegfree_solvent
logical                                         :: detail_temps                 !controls whether or not solute and solvent temps are printed separately (true if solute and solvent degrees of freedom are both not zero)

!----end of shared variables


!----start of public subroutines
contains


!----------------------------------------------------------------------


subroutine md_startup
! initialise used modules
call qatom_startup
call trj_startup
! initialise constants
pi = 4.0*atan(1.0)
deg2rad = pi/180.0
end subroutine md_startup


!----------------------------------------------------------------------

subroutine md_shutdown
! call used modules' shutdown subroutines
call md_deallocate
call topo_deallocate
call qatom_shutdown
end subroutine md_shutdown

!----------------------------------------------------------------------

subroutine die(cause)
! args
character(*), optional          :: cause


!
! exit with an error message
!


! local vars
integer                                         :: i
! flush stuff
integer(4), parameter                   :: stdout_unit = 6

if (nodeid .eq. 0) then
        write(*,*)
        call centered_heading('abnormal termination', '!')
        ! write final energies if run has started
        if (istep > 0) then
                if ( mod(istep,iout_cycle) .ne. 1 ) call write_out
        end if
                if(allocated(v)) then
                !save restart file for diagnosing coordinate problems
                    write(*,*) 'restart file written at step', istep
                        call write_xfin
                endif
        write (*,'(79a)') ('!',i=1,79)
        call close_output_files

        ! apologise
        write(*,'(a)') 'abnormal termination of qdyn5'
        if (present(cause)) then
                write(*,'(79a)') 'terminating due to ', cause
                endif
        write (*,'(79a)') ('!',i=1,79)
#if defined(cray)
        !cray can't flush stdout...
#elif defined(no_flush)
                !when you can't flush
#else
        ! flush stdout
        call flush(stdout_unit)
#endif
end if  
! clean up
call md_deallocate


#if defined (use_mpi)
! abort all processes with exit code 255
call mpi_abort(mpi_comm_world, 255, ierr)
#else
! stop with a message to stderr
stop 'qdyn5 terminated abnormally'
#endif

end subroutine die

!----------------------------------------------------------------------

! --- memory management routines

subroutine allocate_natom_arrays

allocate(x(natom*3), &
                xx(natom*3), &
                v(natom*3), &
                d(natom*3), &
                winv(natom), &
                iqatom(natom), &
                stat=alloc_status)
call check_alloc('atom data arrays')
end subroutine allocate_natom_arrays

!----------------------------------------------------------------------
!allocate arrays that hold no. pairs per chargegroup.
subroutine allocate_nbxx_per_cgp

 allocate(nbpp_per_cgp(ncgp_solute), & 
          nbww_per_cgp(nwat), &
          nbqp_per_cgp(ncgp_solute), &
          nbqw_per_cgp(nwat), &
          nbpw_per_cgp(ncgp_solute), &
          stat=alloc_status)
 call check_alloc('mpi data arrays')

end subroutine allocate_nbxx_per_cgp

!----------------------------------------------------------------------

subroutine allocate_lrf_arrays

if (use_pbc .and. constant_pressure) then
        allocate(iwhich_cgp(natom), lrf(ncgp), old_lrf(ncgp), stat=alloc_status)
else
        allocate(iwhich_cgp(natom), lrf(ncgp), stat=alloc_status)
end if

call check_alloc('lrf arrays')
end subroutine allocate_lrf_arrays

!----------------------------------------------------------------------
#if defined(use_mpi)
subroutine allocate_mpi

if(nodeid .eq. 0) then
 allocate(mpi_status(mpi_status_size,numnodes-1), & 
         request_recv(numnodes-1,3), &
         d_recv(natom*3,numnodes-1), &
         e_recv(numnodes-1), &
         eq_recv(nstates,numnodes-1), &
         stat=alloc_status)
 call check_alloc('mpi data arrays')
else
 allocate(e_send(1), &
          eq_send(nstates), &
          stat=alloc_status)
 call check_alloc('mpi data arrays')
end if
end subroutine allocate_mpi
#endif
!----------------------------------------------------------------------

subroutine allocate_watpol_arrays

allocate(list_sh(n_max_insh, nwpolr_shell), &
nsort(n_max_insh,nwpolr_shell), &
theta(nwat), &
theta0(nwat), &
tdum(nwat), &
stat=alloc_status)
call check_alloc('water polarisation shell arrays')

end subroutine allocate_watpol_arrays

!----------------------------------------------------------------------

subroutine check_alloc(message)

! argument
character(*) message

! local var
character(120) allocmsg

if (alloc_status .ne. 0) then
allocmsg = '>>> out of memory trying to allocate '//message
call die(allocmsg)
end if
end subroutine check_alloc

!----------------------------------------------------------------------

subroutine md_deallocate
! deallocatde this module's own arrays. called by shutdown

! use status to avoid errors if not allocated

! atom arrays
deallocate (x, stat=alloc_status)
deallocate (xx, stat=alloc_status)
deallocate (v, stat=alloc_status)
deallocate (d, stat=alloc_status)
deallocate (winv, stat=alloc_status)
deallocate (iqatom, stat=alloc_status)

! nonbond lists
deallocate(nbpp, nbpw, nbww, nbqq, nbqp, nbqw, qconn, stat=alloc_status)

! watpol arrays
deallocate(wshell, stat=alloc_status)
deallocate(list_sh, stat=alloc_status)
deallocate(nsort, stat=alloc_status)
deallocate(theta, stat=alloc_status)
deallocate(theta0, stat=alloc_status)
deallocate(tdum, stat=alloc_status)

! lrf arrays
deallocate(iwhich_cgp, lrf, stat=alloc_status)

! restraints
deallocate(rstseq, stat=alloc_status)
deallocate(rstpos, stat=alloc_status)
deallocate(rstdis, stat=alloc_status)
deallocate(rstwal, stat=alloc_status)

#if defined (use_mpi)
!mpi arrays
deallocate(nbpp_per_cgp ,stat=alloc_status)
deallocate(nbww_per_cgp ,stat=alloc_status)
deallocate(nbqp_per_cgp ,stat=alloc_status)
deallocate(nbqw_per_cgp ,stat=alloc_status)
deallocate(nbpw_per_cgp ,stat=alloc_status)


if(nodeid .eq. 0) then
   deallocate(d_recv, stat=alloc_status)
   deallocate(e_recv, stat=alloc_status)
   deallocate(eq_recv, stat=alloc_status)
else
   deallocate(e_send, stat=alloc_status)
   deallocate(eq_send, stat=alloc_status)
end if
#endif

end subroutine md_deallocate

!-----------------------------------------------------------------------

subroutine reallocate_nonbondlist_pp

! variables
type(nb_type), allocatable      :: old_nbxx(:)
integer                                         :: old_max

! copy
old_max = calculation_assignment%pp%max
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbpp(1:old_max)

! deallocate and copy back
deallocate(nbpp)
calculation_assignment%pp%max = int(calculation_assignment%pp%max * 1.05) + 200
allocate(nbpp(calculation_assignment%pp%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbpp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%pp%max
100 format('>>> reallocating p-p pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_pp
!----------------------------------------------------------------------
subroutine reallocate_nbpp_cgp

! variables
type(cgp_pair_type), allocatable        :: old_nbxx(:)
integer                                         :: old_max, new_max

! copy
old_max = size(nbpp_cgp, 1)
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbpp_cgp(1:old_max)

! deallocate and copy back
deallocate(nbpp_cgp)
new_max = int( old_max*1.05) + 200
allocate(nbpp_cgp(new_max), stat = alloc_status)
call check_alloc('reallocating non-bonded charge group list')
nbpp_cgp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) new_max
100 format('>>> reallocating p-p charge group pair list, new max is ', i8)

end subroutine reallocate_nbpp_cgp
!----------------------------------------------------------------------
subroutine reallocate_nbpw_cgp

! variables
type(cgp_pair_type), allocatable        :: old_nbxx(:)
integer                                         :: old_max, new_max

! copy
old_max = size(nbpw_cgp, 1)
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbpw_cgp(1:old_max)

! deallocate and copy back
deallocate(nbpw_cgp)
new_max = int( old_max*1.05) + 200
allocate(nbpw_cgp(new_max), stat = alloc_status)
call check_alloc('reallocating non-bonded charge group list')
nbpw_cgp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) new_max
100 format('>>> reallocating p-w charge group pair list, new max is ', i8)

end subroutine reallocate_nbpw_cgp
!----------------------------------------------------------------------
subroutine reallocate_nbqp_cgp

! variables
type(cgp_pair_type), allocatable        :: old_nbxx(:)
integer                                         :: old_max, new_max

! copy
old_max = size(nbqp_cgp, 1)
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded charge group list')
old_nbxx(1:old_max) = nbqp_cgp(1:old_max)

! deallocate and copy back
deallocate(nbqp_cgp)
new_max = int( old_max*1.05) + 200
allocate(nbqp_cgp(new_max), stat = alloc_status)
call check_alloc('reallocating non-bonded charge group list')
nbqp_cgp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) new_max
100 format('>>> reallocating q-p charge group pair list, new max is ', i8)

end subroutine reallocate_nbqp_cgp

!----------------------------------------------------------------------

subroutine reallocate_nonbondlist_pw
! variables
type(nb_type), allocatable      :: old_nbxx(:)
integer                                         :: old_max

! copy
old_max = calculation_assignment%pw%max
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbpw(1:old_max)

! deallocate and copy back
deallocate(nbpw)
calculation_assignment%pw%max = int(calculation_assignment%pw%max * 1.05) + 200
allocate(nbpw(calculation_assignment%pw%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbpw(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%pw%max
100 format('>>> reallocating p-w pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_pw

!----------------------------------------------------------------------

subroutine reallocate_nonbondlist_qp
! variables
type(nbqp_type), allocatable    :: old_nbxx(:)
integer                                         :: old_max

! copy
old_max = calculation_assignment%qp%max
allocate(old_nbxx(old_max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:old_max) = nbqp(1:old_max)

! deallocate and copy back
deallocate(nbqp)
calculation_assignment%qp%max = int(calculation_assignment%qp%max * 1.05) + 200
allocate(nbqp(calculation_assignment%qp%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbqp(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%qp%max
100 format('>>> reallocating q-s pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_qp

!----------------------------------------------------------------------

subroutine reallocate_nonbondlist_ww
! variables
integer(ai), allocatable                :: old_nbxx(:)
integer                                         :: old_max

! copy
allocate(old_nbxx(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('reallocating non-bonded interaction list')
old_nbxx(1:calculation_assignment%ww%max) = nbww(1:calculation_assignment%ww%max)
old_max = calculation_assignment%ww%max

! deallocate and copy back
deallocate(nbww)
calculation_assignment%ww%max = int(calculation_assignment%ww%max * 1.05) + 200 + nwat
allocate(nbww(calculation_assignment%ww%max), stat = alloc_status)
call check_alloc('reallocating non-bonded interaction list')
nbww(1:old_max) = old_nbxx(1:old_max)

! deallocate copy
deallocate(old_nbxx)

! tell the world
write (*,100) calculation_assignment%ww%max
100 format('>>> reallocating w-w pair list, new max is ', i8)

end subroutine reallocate_nonbondlist_ww

!----------------------------------------------------------------------

! --- dynamics subroutines, alphabetically

real(8) function angle(istart, iend)
! *** arguments
integer                                         ::      istart, iend

! *** local variables
integer                                         ::      i,j,k,ia,ic,i3,j3,k3
real(8)                                         ::      bjiinv, bjkinv, bji2inv, bjk2inv
real(8)                                         ::      scp,angv,da,dv,f1
real(8)                                         ::  rji(3),rjk(3),di(3),dk(3) 

! global variables used:
! ang, x, anglib, d

! calculate the total energy of all protein or water angles, depending
! updates d

! reset eangle
angle = 0.

do ia=istart,iend
        ! for each angle in range:

        i  = ang(ia)%i
        j  = ang(ia)%j
        k  = ang(ia)%k
        ic = ang(ia)%cod
        ! calculate rji and rjk
i3=i*3-3
j3=j*3-3
k3=k*3-3
rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)

        ! calculate bjiinv and bjkinv and their squares
bji2inv = 1./(rji(1)**2 + rji(2)**2 + rji(3)**2 )
bjk2inv = 1./(rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
        bjiinv = sqrt(bji2inv)
        bjkinv = sqrt(bjk2inv)

        ! calculate scp and angv
scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
scp = scp * bjiinv*bjkinv
if ( scp .gt.  1.0 ) then
          scp =  1.0
        else if ( scp .lt. -1.0 ) then
          scp = -1.0
        end if
angv = acos(scp)

        ! calculate da and dv
da = angv - anglib(ic)%ang0
angle = angle + 0.5*anglib(ic)%fk*da**2
dv = anglib(ic)%fk*da

        ! calculate f1
f1 = sin ( angv ) 
if ( abs(f1) .lt. 1.e-12 ) then
          ! avoid division by zero
          f1 = -1.e12
        else
  f1 =  -1.0 / f1
        end if

        ! calculate di and dk
di(1) = f1 * ( rjk(1)*bjiinv*bjkinv - scp*rji(1)*bji2inv )
di(2) = f1 * ( rjk(2)*bjiinv*bjkinv - scp*rji(2)*bji2inv )
di(3) = f1 * ( rjk(3)*bjiinv*bjkinv - scp*rji(3)*bji2inv )
dk(1) = f1 * ( rji(1)*bjiinv*bjkinv - scp*rjk(1)*bjk2inv )
dk(2) = f1 * ( rji(2)*bjiinv*bjkinv - scp*rjk(2)*bjk2inv )
dk(3) = f1 * ( rji(3)*bjiinv*bjkinv - scp*rjk(3)*bjk2inv )

        ! update d
d(i3+1) = d(i3+1) + dv*di(1)
d(i3+2) = d(i3+2) + dv*di(2)
d(i3+3) = d(i3+3) + dv*di(3)
d(k3+1) = d(k3+1) + dv*dk(1)
d(k3+2) = d(k3+2) + dv*dk(2)
d(k3+3) = d(k3+3) + dv*dk(3)
d(j3+1) = d(j3+1) - dv*( di(1) + dk(1) )
d(j3+2) = d(j3+2) - dv*( di(2) + dk(2) )
d(j3+3) = d(j3+3) - dv*( di(3) + dk(3) )
end do

end function angle

!-----------------------------------------------------------------------
real(8) function urey_bradley(istart, iend)
! *** arguments
integer                                         ::      istart, iend

! *** local variables
integer                                         ::      i,j,k,ia,ic,i3,j3,k3
real(8)                                         ::      bjiinv, bjkinv, bji2inv, bjk2inv
real(8)                                         ::      scp,angv,da,dv,f1
real(8)                                         ::  rji(3),rjk(3),di(3),dk(3) 
real(8)                                         ::      rik(3), dik, ru, du 
real(8)                                         ::      eurey

! global variables used:
! ang, x, anglib, d

! reset energy
urey_bradley = 0.

do ia=istart,iend
        ! for each angle in range:

        i  = ang(ia)%i
        j  = ang(ia)%j
        k  = ang(ia)%k
        ic = ang(ia)%cod
        ! calculate rji and rjk
i3=i*3-3
j3=j*3-3
k3=k*3-3
        ! 1-3 distance for urey-bradley potential:
    if(anglib(ic)%ureyfk > 0.) then
                rik(1) = x(k3+1) - x(i3+1)
                rik(2) = x(k3+2) - x(i3+2)
                rik(3) = x(k3+3) - x(i3+3)
                dik = sqrt(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                ru = dik - anglib(ic)%ureyr0
                urey_bradley = urey_bradley + anglib(ic)%ureyfk*ru**2
                du = 2*anglib(ic)%ureyfk*ru/dik
                d(k3+1) = d(k3+1) + du*rik(1)
                d(k3+2) = d(k3+2) + du*rik(2)
                d(k3+3) = d(k3+3) + du*rik(3)
                d(i3+1) = d(i3+1) - du*rik(1)
                d(i3+2) = d(i3+2) - du*rik(2)
                d(i3+3) = d(i3+3) - du*rik(3)
        end if
end do

end function urey_bradley

!-----------------------------------------------------------------------

real(8) function bond(istart, iend)
! *** arguments
integer                                         ::      istart, iend

! *** local variables
integer                                         ::      i,j,ib,ic,i3,j3
real(8)                                         ::      b,db,dv
real(8)                                         ::      rij(3)

! global variables used:
! bnd, x, bondlib, d

! reset ebond
bond = 0

do ib=istart,iend
        ! for each bond in range:

        i  = bnd(ib)%i
        j  = bnd(ib)%j
        ic = bnd(ib)%cod
        ! calculate rij
        i3=i*3-3
        j3=j*3-3
        rij(1) = x(j3+1) - x(i3+1)
        rij(2) = x(j3+2) - x(i3+2)
        rij(3) = x(j3+3) - x(i3+3)

        ! calculate b and db, update ebond
        b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
        db = b - bondlib(ic)%bnd0
        bond = bond + 0.5*bondlib(ic)%fk*db**2

        ! calculate dv and update d
        dv = bondlib(ic)%fk*db/b
        d(j3+1) = d(j3+1) + rij(1)*dv
        d(j3+2) = d(j3+2) + rij(2)*dv
        d(j3+3) = d(j3+3) + rij(3)*dv
        d(i3+1) = d(i3+1) - rij(1)*dv
        d(i3+2) = d(i3+2) - rij(2)*dv
        d(i3+3) = d(i3+3) - rij(3)*dv
end do

end function bond
!-----------------------------------------------------------------------
subroutine cgp_centers
! *** local variables
integer                                         ::      ig,i,i3

do ig = 1, ncgp
        lrf(ig)%cgp_cent(:) = 0.
        lrf(ig)%phi0 = 0.
lrf(ig)%phi1(:) = 0.
lrf(ig)%phi2(:) = 0.
lrf(ig)%phi3(:) = 0.

        do i  = cgp(ig)%first, cgp(ig)%last
                lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:) + x(cgpatom(i)*3-2:cgpatom(i)*3)
        end do

lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:)/real(cgp(ig)%last - cgp(ig)%first +1)

end do

end subroutine cgp_centers
!-----------------------------------------------------------------------
subroutine make_nbqqlist
!locals
integer                                         ::      is

call make_qconn
nbqq_max = nbqq_count()
allocate(nbqq(nbqq_max, nstates), stat=alloc_status)
call check_alloc('qatom-qatom non-bond list')
! prepare q-atom nonbond lists that do not need updating
call nbqqlist

do is =1, nstates
        write (*,200) nbqq_pair(is),is
end do
write (*,*)

200 format ('no. of rcq indep. nb pairs involving q-atoms = ',i5, &
' in state :',i3)
end subroutine make_nbqqlist

!-----------------------------------------------------------------------

subroutine distribute_nonbonds
!locals
integer                                 :: npp, npw, nqp, nww, nqw
type(node_assignment_type),allocatable  :: node_assignment(:)
real                                    :: avgload, old_avgload
integer                                 :: i, last_cgp, last_pair
integer                                 :: mpitype_pair_assignment, mpitype_node_assignment
integer                                 :: average_pairs,inode,icgp,sum,less_than_sum
integer                                 :: n_bonded, n_nonbonded, master_assign
real                                   :: percent
integer                                 :: master_sum
!!!!tmp vars f�r allokering
integer,parameter :: vars = 5
integer :: i_loop
!!!


! count the number of nonbonded interactions and distribute them among the nodes

if (nodeid .eq. 0) then
! nice header
call centered_heading('distribution of charge groups','-')

!allocate node_assignment
allocate(node_assignment(0:numnodes-1),stat=alloc_status)
call check_alloc('node_assignment')

!allocate arrays that hold no. pairs per chargegroup.
call allocate_nbxx_per_cgp

!count stuff for balancing nodes and allocating nonbond arrays nbxx()
nbpp_per_cgp = 0
nbpw_per_cgp = 0
nbqp_per_cgp = 0
nbqw_per_cgp = 0
nbww_per_cgp = 0

call nbpp_count(npp, nbpp_per_cgp)    !only for switching atoms!!!?
call nbpw_count(npw, nbpw_per_cgp)
call nbqp_count(nqp, nbqp_per_cgp)
call nbqw_count(nqw, nbqw_per_cgp) 
call nbww_count(nww, nbww_per_cgp) 

!for keeping track of actual # of nonbonded pairs
totnbpp = npp
totnbpw = npw
totnbww = nww*9
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

#if defined (use_mpi)
else ! i.e. slave nodes exists

! a simple solution to avoid parallelising the bonded
! calculate n_bonded and n_nonbonded 
! approximate time of computing one bonded with one nonbonded
! the number of qq-interactions are neglected
n_bonded = nbonds + nangles + ntors + nimps
n_nonbonded = totnbpp + totnbpw + totnbww + totnbqw + totnbqp

! compare to determine how many nonbonded master should get
! a bonded is faster, so this favours an early completion for master
master_assign =  n_nonbonded/numnodes - n_bonded * numnodes

! calculate the assignments ********************

!calculate balanced assignment for p-p pairs
icgp=0
sum=0
 !first assign master a small part
node_assignment(0)%pp%start=icgp+1
percent=real(totnbpp)/n_nonbonded
less_than_sum = master_assign*percent  ! no. of pp-type to assign master
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpp_per_cgp(icgp)
end do
node_assignment(0)%pp%end=icgp
master_sum=sum
 !now assign slaves
average_pairs=(totnbpp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pp%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum) 
     icgp=icgp+1
     sum=sum + nbpp_per_cgp(icgp)
  end do
  node_assignment(inode)%pp%end=icgp
end do
node_assignment(numnodes-1)%pp%start=icgp+1
node_assignment(numnodes-1)%pp%end=ncgp_solute

!calculate balanced assignment for p-w pairs
icgp=0
sum=0
node_assignment(0)%pw%start=icgp+1
percent=real(totnbpw)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpw_per_cgp(icgp)
end do
node_assignment(0)%pw%end=icgp
master_sum=sum
average_pairs=(totnbpw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pw%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbpw_per_cgp(icgp)
  end do
  node_assignment(inode)%pw%end=icgp
end do
node_assignment(numnodes-1)%pw%start=icgp+1
node_assignment(numnodes-1)%pw%end=ncgp_solute

!calculate balanced assignment for q-p pairs
icgp=0
sum=0
node_assignment(0)%qp%start=icgp+1
percent=real(totnbqp)/n_nonbonded
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

!calculate balanced assignment for w-w pairs
icgp=0
sum=0
node_assignment(0)%ww%start=icgp+1
percent=real(totnbww)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbww_per_cgp(icgp)
end do
node_assignment(0)%ww%end=icgp
master_sum=sum
average_pairs=(totnbww-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%ww%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbww_per_cgp(icgp)
  end do
  node_assignment(inode)%ww%end=icgp
end do
node_assignment(numnodes-1)%ww%start=icgp+1
node_assignment(numnodes-1)%ww%end=nwat

!calculate balanced assignment for q-w pairs
icgp=0
sum=0
node_assignment(0)%qw%start=icgp+1
percent=real(totnbqw)/n_nonbonded
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
#if defined (use_mpi)
if (numnodes .gt. 1) then
    if (nodeid .ne. 0) then
        ! dummy allocation to avoid runtime errors when using pointer checking
        allocate(node_assignment(1),stat=alloc_status)
    endif 
! register data types
call mpi_type_contiguous(3, mpi_integer, mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom mpi data type')
call mpi_type_commit(mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom mpi data type')

call mpi_type_contiguous(5, mpitype_pair_assignment, mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom mpi data type')
call mpi_type_commit(mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom mpi data type')

! distribute
call mpi_scatter(node_assignment, 1, mpitype_node_assignment, &
    calculation_assignment, 1, mpitype_node_assignment, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('failure while sending node assignments')

! free data type
call mpi_type_free(mpitype_node_assignment, ierr)
call mpi_type_free(mpitype_pair_assignment, ierr)

    if (nodeid .ne. 0) then
        deallocate(node_assignment)
    endif 
end if
#endif

if (nodeid .eq. 0) then
 ! print a status report
 write(*,98) 'solute-solute', 'solute-water', 'water-water', 'q-solute', 'q-water'
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

#if defined (use_mpi)
call mpi_bcast(totnbpp, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_bcast(totnbpw, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_bcast(totnbqp, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_bcast(totnbqw, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_bcast(totnbww, 1, mpi_integer, 0, mpi_comm_world, ierr)
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
call check_alloc('qatom - solute non-bond list')

calculation_assignment%qw%max = nwat
allocate(nbqw(calculation_assignment%qw%max), stat=alloc_status)
call check_alloc('qatom - water non-bond list')



  if (use_pbc) then
        !allocate array to keep track of chargegroups
        !approximate with one half of the number of atompairs
        allocate(nbpp_cgp(calculation_assignment%pp%max / 2), stat=alloc_status)
        call check_alloc('solute-solute non-bond charge group pair list')
        allocate(nbpw_cgp(calculation_assignment%pw%max / 2), stat=alloc_status)
        call check_alloc('solute-solvent non-bond charge group pair list')
        allocate(nbqp_cgp(calculation_assignment%qp%max / 2), stat=alloc_status)
        call check_alloc('qatom-solute non-bond charge group pair list')
  end if        

!kanske deallokera nbxx_per_cgp todo

98 format('node value ',5a13)
99 format(a10,1x,5(1x,i12))
!99 format(a4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)
100 format(i4,1x,a5,1x,5(1x,i12))
!100 format(i4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)

if (nodeid .eq. 0)  call centered_heading('end of distribution', '-')

end subroutine distribute_nonbonds

!-----------------------------------------------------------------------


subroutine close_input_files
close (1)
if(restart) close (2)
if ( implicit_rstr_from_file .eq. 1 ) close (12)
close (13)

end subroutine close_input_files

!-----------------------------------------------------------------------

subroutine close_output_files
close (3)
if ( itrj_cycle .gt. 0 ) close (10)
if ( iene_cycle .gt. 0 ) close (11)

end subroutine close_output_files

!-----------------------------------------------------------------------

subroutine open_files
! --> restart file (2)
if(restart) then
open (unit=2, file=restart_file, status='old', form='unformatted', action='read', err=2)
end if

! --> final coords (3)
open (unit=3, file=xfin_file, status='unknown', form='unformatted', action='write', err=3)

! --> energy output file (11)
if ( iene_cycle .gt. 0 ) then
open (unit=11, file=ene_file, status='unknown', form='unformatted', action='write', err=11)
end if

! --> external file for implicit position restraints (12)
if ( implicit_rstr_from_file .eq. 1 ) then
open (unit=12, file=exrstr_file, status='old', form='unformatted', action='read', err=12)
end if


return

! crude error handling
2 call die('error opening restart file.')
3 call die('error opening final coordinates file.')
11 call die('error opening energy output file.')
12 call die('error opening position restraints file.')

end subroutine open_files

!-----------------------------------------------------------------------

!restrain all excluded atoms plus heavy solute atoms in the inner shell.
subroutine fix_shell
! local variables
integer                                         ::      i,i3
real(8)                                         ::      fk,r2,erst
real(8)                                         ::      dr(3)

! global variables used:
!  e, nat_pro, excl, shell, heavy, fk_fix, fk_pshell, x, xtop, d

do i = 1, nat_pro
if (excl(i) .or. shell(i)) then
  ! decide which fk to use
if ( excl(i) ) then 
    fk = fk_fix
  else
fk = fk_pshell
  end if
i3 = i*3-3

  ! calculate drift from topology
dr(1)   = x(i3+1) - xtop(i3+1)
dr(2)   = x(i3+2) - xtop(i3+2)
dr(3)   = x(i3+3) - xtop(i3+3)
r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
erst    = 0.5*fk*r2

  ! update restraint energies
if ( excl(i) ) e%restraint%fix   = e%restraint%fix + erst
if ( shell(i) ) e%restraint%shell = e%restraint%shell + erst 

  ! update forces
d(i3+1) = d(i3+1) + fk*dr(1)
d(i3+2) = d(i3+2) + fk*dr(2)
d(i3+3) = d(i3+3) + fk*dr(3)
end if
end do
end subroutine fix_shell

!-----------------------------------------------------------------------

subroutine gauss (am,sd,v,ig)
! arguments
real(8)                                 ::      am,sd,v
integer                                 ::      ig

! local variables
integer                                 ::      i
real(8)                                 ::      a,y

a=0.0
do i=1,12
y=randm(ig)
a=a+y
end do
v=(a-6.0)*sd+am
end subroutine gauss

!-----------------------------------------------------------------------
subroutine get_fep
! local variables
character                                       ::      libtext*80,qaname*2
integer                                 ::      i,j,k,iat
!temp. array for reallocating long-range exclusion list
integer(ai), pointer    ::      tempexlong(:,:)

! --- # states, # q-atoms
if(.not. qatom_load_atoms(fep_file)) then
        call die('failure to load q-atoms from fep file.')
end if

! set flags
do i=1,nqat
        if(iqseq(i) > 0 .and. iqseq(i) <= nat_solute)  then
                iqatom(iqseq(i)) = i
        else if(iqseq(i) == 0) then
                write(*,10) i
        else
                write(*,20) i, iqseq(i)
                call die('invalid q-atom data')
        end if
end do
10      format('>>> warning: q-atom no. ',i2,' is not associated with a topology atom.')
20      format('>>>>> error: q-atom no. ',i2,' has invalid topology number ',i5)
!allocate memory for qatom charges
allocate(qcrg(nqat,nstates), stat=alloc_status)
call check_alloc('qatom charges')

! --- copy topology charges

do i=1,nqat
        do j=1,nstates
                qcrg(i,j)=crg(iqseq(i))
        end do
end do

!initialize softcore lookup array
allocate (sc_lookup(nqat,natyps+nqat,nstates))
sc_lookup(:,:,:)=0.0

!load rest of fep file
if(.not. qatom_load_fep(fep_file)) then
        call die('failure to load fep file.')
end if

!adapt lj parameters to topology
!if arithmetic combination rule take sqrt(epsilon) now
if (qvdw_flag .and. ivdw_rule .eq. 2 ) then
        qbvdw(:,1) = sqrt( qbvdw(:,1) )
        qbvdw(:,3) = sqrt( qbvdw(:,3) )
end if

!remove redefined bonded interactions from topology
if(nqbond > 0 .or. nqangle > 0 .or. nqtor > 0 .or. nqimp > 0 ) then
        write(*,*)
        call centered_heading('removing redefined interactions from topology','-')
230             format('type',t10,' atom1 atom2 atom3 atom4')
        write(*,230)
231             format(a,t10,4i6)
        !remove bonds that were redefined
        do i=1,nbonds
                do j=1,nqbond
                        if ( (bnd(i)%i==qbnd(j)%i .and. bnd(i)%j==qbnd(j)%j) .or. &
                                (bnd(i)%i==qbnd(j)%j .and. bnd(i)%j==qbnd(j)%i) ) then
                                bnd(i)%cod = 0
                                write (*,231) 'bond',bnd(i)%i,bnd(i)%j
                        end if
                end do
        end do

        !remove angles that were redefined
        do i=1,nangles
                do j=1,nqangle
                        if((ang(i)%i.eq.qang(j)%i .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%k)                          .or. &
                                (ang(i)%i.eq.qang(j)%k .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%i) )                         then

                                ang(i)%cod = 0
                                write (*,231) 'angle',ang(i)%i,ang(i)%j,ang(i)%k
                        end if
                end do
        end do

        !remove torsions that were redefined
        do i=1,ntors
                do j=1,nqtor
                                                if(( (tor(i)%i.eq.iqtor(j) .and. tor(i)%j.eq.jqtor(j) .and. &
                                                        tor(i)%k.eq.kqtor(j) .and. tor(i)%l.eq.lqtor(j)) .or. &
                                (tor(i)%i.eq.lqtor(j) .and. tor(i)%j.eq.kqtor(j) .and. &
                                                                tor(i)%k.eq.jqtor(j) .and. tor(i)%l.eq.iqtor(j)) ) .and. &
                                tor(i)%cod /= 0) then
                                tor(i)%cod = 0
                                write (*,231) 'torsion', tor(i)%i,tor(i)%j,tor(i)%k,tor(i)%l
                                                end if
                end do
        end do


        !remove impropers that were redefined
        select case(ff_type)
        case(ff_charmm) !special code for charmm
                do i=1,nimps
                        do j=1,nqimp
                                if(((imp(i)%i .eq. iqimp(j)) .or. &
                                        (imp(i)%i .eq. lqimp(j)) .or. &
                                        (imp(i)%l .eq. iqimp(j)) .or. &
                                        (imp(i)%l .eq. lqimp(j))) .and. &
                                        ((imp(i)%j .eq. iqimp(j)) .or. &
                                        (imp(i)%j .eq. jqimp(j))  .or. &
                                        (imp(i)%j .eq. kqimp(j))  .or. &
                                        (imp(i)%j .eq. lqimp(j))) .and. &
                                        ((imp(i)%k .eq. iqimp(j)) .or. &
                                        (imp(i)%k .eq. jqimp(j)) .or. &
                                        (imp(i)%k .eq. kqimp(j)) .or. &
                                        (imp(i)%k .eq. lqimp(j))) .and. &
                                        imp(i)%cod /= 0) then
                                        imp(i)%cod = 0
                                        write (*,231) &
                                        'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                                end if
                        end do
                end do

        case default
        do i=1,nimps
                do j=1,nqimp
                        if(((imp(i)%j.eq.jqimp(j) .and. imp(i)%k.eq.kqimp(j)) .or. &
                                (imp(i)%j.eq.kqimp(j) .and. imp(i)%k.eq.jqimp(j))) .and. &
                                imp(i)%cod /= 0) then
                                imp(i)%cod = 0
                                write(*,231)'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                        end if
                end do
        end do
        end select
end if

!check special exclusions
!modify exclusion lists to inclue special exclusions between q and non-q
if(nexspec > 0) then
        allocate(tempexlong(2,nexlong+nexspec))
        tempexlong(:, 1:nexlong) = listexlong(:, 1:nexlong)
        deallocate(listexlong)
        listexlong => tempexlong
end if

do k = 1, nexspec
        i = exspec(k)%i
        j = exspec(k)%j
        if(i < 1 .or. i > nat_pro .or. j < 1 .or. j > nat_pro) then
                write(*, 592) k, i, j
                call die('invalid special exclusion data')
        end if
        !if one or more non-q-atoms modify exclusion lists
        if(iqatom(i)==0 .or. iqatom(j)==0) then
                !with non-q-atoms involved only accept all or no states
                if(any(exspec(k)%flag(1:nstates))) then
                        if(.not. all(exspec(k)%flag(1:nstates))) then
                                write(*,594) k
                                call die('invalid special exclusion data')
                        else !exlcude in all states
                                if(abs(j-i) <= max_nbr_range) then
                                        if(i < j) then
                                                listex(j-i,i) = .true.
                                        else
                                                listex(i-j,j) = .true.
                                        end if
                                else
                                        nexlong = nexlong + 1
                                        listexlong(1, nexlong) = i
                                        listexlong(2, nexlong) = j
                                end if
                        end if
                end if
        end if
end do
592     format('>>>>> error: special exclusion pair ',i2,' (',i5,1x,i5,') is invalid')
594     format('>>>>> error: non-q-atom special excl. pair ',i2,' must be on in all or no states')
end subroutine get_fep

!-----------------------------------------------------------------------

subroutine get_fname (text,length,filnam)
! arguments
character                                       ::      text*80,filnam*80
integer                                 ::      length
! local variables


integer                                 ::      i

length=80
do i=1,80
if ( text(i:i) .eq. ' ' ) then
length=i-1
goto 10
end if
end do
10 filnam(1:length)=text(1:length)

end subroutine get_fname

!-----------------------------------------------------------------------

real(8) function improper(istart, iend)
!arguments
integer                                         ::      istart, iend

! evaluate harmonic impropers
! local variables
integer                                         ::      ip
real(8)                                         ::      scp,phi,dv,arg,f1
real(8)                                         ::      bjinv, bkinv, bj2inv, bk2inv
real(8)                                         ::      rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)                                         ::      rki(3),rlj(3),dp(12),di(3),dl(3)
type(tor_type), pointer         ::      t
type(implib_type), pointer      ::      lib

! global variables used:
!  imp, implib, x, pi, d

improper = 0.

do ip = iend, istart,-1
        t => imp(ip)
        lib => implib(t%cod)
rji(1) = x(t%i*3-2) - x(t%j*3-2)
rji(2) = x(t%i*3-1) - x(t%j*3-1)
rji(3) = x(t%i*3-0) - x(t%j*3-0)
rjk(1) = x(t%k*3-2) - x(t%j*3-2)
rjk(2) = x(t%k*3-1) - x(t%j*3-1)
rjk(3) = x(t%k*3-0) - x(t%j*3-0)
rkl(1) = x(t%l*3-2) - x(t%k*3-2)
rkl(2) = x(t%l*3-1) - x(t%k*3-1)
rkl(3) = x(t%l*3-0) - x(t%k*3-0)


rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

bj2inv  = 1./( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
bk2inv  = 1./( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  1.0 ) then
                scp =  1.0
else if ( scp .lt. -1.0 ) then 
                scp = -1.0
        end if
phi = acos ( scp )
if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) &
                 .lt. 0) then 
      phi = -phi
        end if

! ---       energy

arg = phi - lib%imp0
arg = arg - 2.*pi*nint(arg/(2.*pi))
dv  = lib%fk*arg
improper = improper + 0.5*dv*arg

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
di(1) = f1 * ( rnk(1)*bjinv*bkinv - scp*rnj(1)*bj2inv )
di(2) = f1 * ( rnk(2)*bjinv*bkinv - scp*rnj(2)*bj2inv )
di(3) = f1 * ( rnk(3)*bjinv*bkinv - scp*rnj(3)*bj2inv )
dl(1) = f1 * ( rnj(1)*bjinv*bkinv - scp*rnk(1)*bk2inv )
dl(2) = f1 * ( rnj(2)*bjinv*bkinv - scp*rnk(2)*bk2inv )
dl(3) = f1 * ( rnj(3)*bjinv*bkinv - scp*rnk(3)*bk2inv )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)
rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)
dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(t%i*3-2) = d(t%i*3-2) + dv*dp(1)
d(t%i*3-1) = d(t%i*3-1) + dv*dp(2)
d(t%i*3-0) = d(t%i*3-0) + dv*dp(3)
d(t%j*3-2) = d(t%j*3-2) + dv*dp(4)
d(t%j*3-1) = d(t%j*3-1) + dv*dp(5)
d(t%j*3-0) = d(t%j*3-0) + dv*dp(6)
d(t%k*3-2) = d(t%k*3-2) + dv*dp(7)
d(t%k*3-1) = d(t%k*3-1) + dv*dp(8)
d(t%k*3-0) = d(t%k*3-0) + dv*dp(9)
d(t%l*3-2) = d(t%l*3-2) + dv*dp(10)
d(t%l*3-1) = d(t%l*3-1) + dv*dp(11)
d(t%l*3-0) = d(t%l*3-0) + dv*dp(12)
end do
end function improper

!-----------------------------------------------------------------------

real(8) function improper2(istart, iend)
!evaluate periodic impropers
!arguments
integer                                         ::      istart, iend
! local variables
integer                                         ::      ip
real(8)                                         ::      scp,phi,dv,arg,f1
real(8)                                         ::      bjinv, bkinv, bj2inv, bk2inv
real(8)                                         ::      rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)                                         ::      rki(3),rlj(3),dp(12),di(3),dl(3)
type(tor_type), pointer         ::      t
type(implib_type), pointer      ::      lib

! global variables used:
! imp, implib, x, pi, d

improper2 = 0.

do ip = iend, istart,-1
        t => imp(ip)
        lib => implib(t%cod)
rji(1) = x(t%i*3-2) - x(t%j*3-2)
rji(2) = x(t%i*3-1) - x(t%j*3-1)
rji(3) = x(t%i*3-0) - x(t%j*3-0)
rjk(1) = x(t%k*3-2) - x(t%j*3-2)
rjk(2) = x(t%k*3-1) - x(t%j*3-1)
rjk(3) = x(t%k*3-0) - x(t%j*3-0)
rkl(1) = x(t%l*3-2) - x(t%k*3-2)
rkl(2) = x(t%l*3-1) - x(t%k*3-1)
rkl(3) = x(t%l*3-0) - x(t%k*3-0)
rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)


bj2inv  = 1./( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
bk2inv  = 1./( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  1.0 ) then
                scp =  1.0
else if ( scp .lt. -1.0 ) then 
                scp = -1.0
        end if
phi = acos ( scp )
if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) &
                 .lt. 0) then 
      phi = -phi
        end if

! ---       energy

arg = 2*phi - lib%imp0
improper2 = improper2 + lib%fk * (1 + cos(arg))
dv  = -2*lib%fk * sin(arg)

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
di(1) = f1 * ( rnk(1)*bjinv*bkinv - scp*rnj(1)*bj2inv )
di(2) = f1 * ( rnk(2)*bjinv*bkinv - scp*rnj(2)*bj2inv )
di(3) = f1 * ( rnk(3)*bjinv*bkinv - scp*rnj(3)*bj2inv )
dl(1) = f1 * ( rnj(1)*bjinv*bkinv - scp*rnk(1)*bk2inv )
dl(2) = f1 * ( rnj(2)*bjinv*bkinv - scp*rnk(2)*bk2inv )
dl(3) = f1 * ( rnj(3)*bjinv*bkinv - scp*rnk(3)*bk2inv )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)
rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)
dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(t%i*3-2) = d(t%i*3-2) + dv*dp(1)
d(t%i*3-1) = d(t%i*3-1) + dv*dp(2)
d(t%i*3-0) = d(t%i*3-0) + dv*dp(3)
d(t%j*3-2) = d(t%j*3-2) + dv*dp(4)
d(t%j*3-1) = d(t%j*3-1) + dv*dp(5)
d(t%j*3-0) = d(t%j*3-0) + dv*dp(6)
d(t%k*3-2) = d(t%k*3-2) + dv*dp(7)
d(t%k*3-1) = d(t%k*3-1) + dv*dp(8)
d(t%k*3-0) = d(t%k*3-0) + dv*dp(9)
d(t%l*3-2) = d(t%l*3-2) + dv*dp(10)
d(t%l*3-1) = d(t%l*3-1) + dv*dp(11)


d(t%l*3-0) = d(t%l*3-0) + dv*dp(12)
end do
end function improper2

!-----------------------------------------------------------------------

#if defined (use_mpi)
!defines and allocates variables needed in the md-calculations
!the node initiation is written for ai = 4. if changes are made to any size in
! sizes.f90 the mpi-code must be changed accordingly. it is not dynamically
! implemented yet.
subroutine init_nodes
!
! initialise slave nodes, sending to slaves:
!
! variables:
!  natom,nwat,nsteps,use_lrf,nbcycle,crg_ow,crg_hw,rcpp,rcww,rcpw,rcq,xpcent
!  nat_solute,ncgp,ncgp_solute,ivdw_rule,iuse_switch_atom,el14_scale,n14long
!  nexlong,natyps,nljtyp,rexcl_o,nstates,nqat,qvdw_flag,nqlib,rclrf,
!  use_pbc, qswitch, nmol, nat_pro
!
! arrays:
!  x,v,iqatom,ljcod,qconn,iwhich_cgp,lrf,excl,iac,crg,cgpatom,cgp,iaclib
!  list14,listex,list14long,listexlong,iqseq,qiac,qcrg,qavdw,qbvdw,eq(:)%lambda,
!  boxlength, inv_boxl, boxcentre, sc_lookup 
!

integer, parameter                      :: vars = 40    !increment this var when adding data to broadcast in batch 1
integer                         :: blockcnt(vars), ftype(vars) 
integer(kind=mpi_address_kind)                          :: fdisp(vars)
integer                                 :: mpitype_batch,mpitype_batch2
integer                                 :: nat3
real(kind=wp8), allocatable             :: temp_lambda(:)
integer, parameter                      ::maxint=2147483647
real(kind=wp8), parameter                        ::maxreal=1e35
integer  :: mpi_ai_integer, mpi_tiny_integer, i_loop

!external mpi_address
!external mpi_bcast

!**********
!2002-11-28 
!mn-> this will work with new implementations of mpi standard >= 2
!the mpi library at pdc does not support these definitions when i tried to use them.
!using these routines will allow a change made to the sizes in sizes.f90 to
! affect the mpi. without them the variables below marked (ai) and (tiny) will have to 
! be changed manually.
!when using this part make sure the vars marked with comments (ai) and (tiny) are 
! changed to mpi_ai_integer and mpi_tiny_integer.

!external mpi_type_create_f90_integer
!external mpi_sizeof

!define data types
! this is wrong, the 1:st param is "precision, in decimal digits", not bits
!call mpi_type_create_f90_integer((8*ai-1),mpi_ai_integer,ierr)
!call mpi_type_create_f90_integer((8*tiny-1),mpi_tiny_integer,ierr)
!to check the size in bytes of the new types use
!call mpi_sizeof(mpi_ai_integer,size,ierr)
!call mpi_sizeof(mpi_tiny_integer,size,ierr)
!***************************


if (nodeid .eq. 0) call centered_heading('distributing data to slave nodes', '-')

! --- mandatory data, first batch ---

if (nodeid .eq. 0) write (*,'(80a)') 'md data, first batch'

! broadcast some initial variables

! run control constants: natom, nwat, nsteps, nbmethod, nbcycle
call mpi_bcast(natom, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast natom')
call mpi_bcast(nwat, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nwat')
call mpi_bcast(nsteps, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nsteps')
call mpi_bcast(use_lrf, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast use_lrf')
call mpi_bcast(nbcycle, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nbcycle')

! water parameters: crg_ow, crg_hw (used by nonbond_ww)
call mpi_bcast(crg_ow, 1, mpi_real4, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast crg_ow')
call mpi_bcast(crg_hw, 1, mpi_real4, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast crg_hw')

! cutoffs: rcpp, rcww, rcpw, rcq, rclrf (used by pair list generating functions)
call mpi_bcast(rcpp, 1, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rcpp')
call mpi_bcast(rcww, 1, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rcww')
call mpi_bcast(rcpw, 1, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rcpw')
call mpi_bcast(rcq, 1, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rcq')
call mpi_bcast(rclrf, 1, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rclrf')

!periodic boudary condition
call mpi_bcast(use_pbc, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast use_pbc')
call mpi_bcast(boxcentre, 3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast boxcentre')
call mpi_bcast(boxlength, 3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast boxlength')
call mpi_bcast(inv_boxl, 3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast inv_boxl')
call mpi_bcast(qswitch, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qswitch')
call mpi_bcast(constant_pressure, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast constant_pressure')
call mpi_bcast(ivolume_cycle, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast ivolume_cycle')
call mpi_bcast(rigid_box_centre, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rigid_box_centre')
call mpi_bcast(put_solvent_back_in_box, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast put_solvent_back_in_box')
call mpi_bcast(put_solute_back_in_box, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast put_solute_back_in_box')

! xpcent            from topo, needed for listgeneration
call mpi_bcast(xpcent, 3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast xpcent')

!**mn-> beh�vs om shake ska parallelliseras
! shake/temperature parameters
!call mpi_bcast(shake_constraints, 1, mpi_integer, 0, mpi_comm_world, ierr)  !bara i init_shake & md_run
!call mpi_bcast(shake_molecules, 1, mpi_integer, 0, mpi_comm_world, ierr)    !bara i div init_
!call mpi_bcast(ndegf, 1, mpi_integer, 0, mpi_comm_world, ierr)    !bara i div init_
!call mpi_bcast(ndegfree, 1, mpi_integer, 0, mpi_comm_world, ierr)    !bara i div init_

! a bunch of vars from the topo module
call mpi_bcast(nat_solute, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nat_solute')
call mpi_bcast(ncgp, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast ncgp')
call mpi_bcast(ncgp_solute, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast ncgp_solute')
call mpi_bcast(ivdw_rule, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast ivdw_rule')
call mpi_bcast(iuse_switch_atom, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast iuse_switch_atom')
call mpi_bcast(el14_scale, 1, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast el14_scale')
call mpi_bcast(n14long, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast n14long')
call mpi_bcast(nexlong, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nexlong')
call mpi_bcast(natyps, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast natyps')
call mpi_bcast(rexcl_o, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast rexcl')
call mpi_bcast(nmol, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nmol')
call mpi_bcast(nat_pro, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nat_pro')

!vars from qatom
call mpi_bcast(nstates, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nstates')
call mpi_bcast(nqat, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nqat')
call mpi_bcast(qvdw_flag, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qvdw_flag')
call mpi_bcast(nqlib, 1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast nqlib')

!setting all vars not sent to slaves to 2147483647. to avoid hidden bugs.
if (nodeid .ne. 0) then 
shake_constraints=maxint 
shake_molecules=maxint 
ndegf=maxint 
ndegfree=maxint
xwcent(:)=maxreal 
end if

! --- md data, second batch ---

if (nodeid .eq. 0) write (*,'(80a)') 'md data, second batch'

! allocate arrays
if (nodeid .ne. 0) then
call allocate_natom_arrays
end if

! broadcast x, v and winv
nat3 = 3*natom
call mpi_bcast(x, nat3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast x')
call mpi_bcast(v, nat3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast v')

!setting all vars not sent to slaves to 2147483647. to avoid conflicts.
if (nodeid .ne. 0) then 
winv=maxint 
end if

!broadcast iqatom
call mpi_bcast(iqatom, natom, mpi_integer2, 0, mpi_comm_world, ierr) !(tiny)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast iqatom')

!broadcast ljcod
call mpi_bcast(ljcod, size(ljcod), mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast ljcod')

!broadcast qconn(nstates,nat_solute, nqat)
if (nodeid .ne. 0) then
allocate(qconn(nstates,nat_solute, nqat),stat=alloc_status)
call check_alloc('qconn')
end if
call mpi_bcast(qconn, size(qconn), mpi_integer2, 0, mpi_comm_world, ierr) ! (tiny)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qconn')


! --- periodic boundary condition data ---




! --- shake data ---

!if (shake_solute .or. shake_solvent .or. shake_hydrogens) then
! shake stuff

!if (nodeid .eq. 0) write (*,'(80a)') 'shake data'

!if (nodeid .ne. 0) then
! allocate shake arrays
! add code here to allocate shake array!
!end if

! bake all shake data into a big packet & bcast
! add code here to broadcast shake data
!end if

! --- lrf data ---

if (use_lrf) then
! lrf stuff

if (nodeid .eq. 0) write (*,'(80a)') 'lrf data'

! allocate arrays
if (nodeid .ne. 0) call allocate_lrf_arrays

!mpi_integer4 is used instead of mpi_ai_integer
!change to mpi_type_create, see note above or note2 in sizes.f90
! iwhich_cgp
call mpi_bcast(iwhich_cgp, natom, mpi_integer4, 0, mpi_comm_world, ierr) !(ai)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast lrf parameters')

! lrf
ftype(:) = mpi_real8
blockcnt(1) = 3                                 ! real(8) cgp_cent(3)
fdisp(1) = 0
blockcnt(2) = 1                                 ! real(8) phi0
fdisp(2) = 3*8
blockcnt(3) = 3                                 ! real(8) phi1(3)
fdisp(3) = 3*8 + 8
blockcnt(4) = 9                                 ! real(8) phi2(9)
fdisp(4) = 3*8 + 8 + 3*8
blockcnt(5) = 27                                ! real(8) phi3(27)
fdisp(5) = 3*8 + 8 + 3*8 + 9*8
call mpi_type_create_struct(5, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_create_struct')
call mpi_type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_commit')
call mpi_bcast(lrf, ncgp, mpitype_batch, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast shake parameters')
call mpi_type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_free')
end if !(use_lrf)

! --- data from the topo module ---

if (nodeid .eq. 0) write (*,'(80a)') 'topo data'

! allocate topology arrays
if (nodeid .ne. 0) then
! don't allocate memory for stuff we don't need
! these array size variables are actually used
max_cgp=ncgp
max_atyps = natyps
max_14long = n14long
max_exlong = nexlong
max_atom = natom

call topo_allocate_atom(alloc_status)
call check_alloc('topology arrays')
call topo_allocate_potential(alloc_status)
call check_alloc('topology arrays')
allocate(istart_mol(nmol+1), &
stat=alloc_status)
call check_alloc('topology arrays')
end if

! broadcast excl
call mpi_bcast(excl, natom, mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast excl')
! broadcast istart_mol
call mpi_bcast(istart_mol, nmol+1, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast istart_mol')

! bcast iac, crg and cgpatom 
call mpi_bcast(iac, natom, mpi_integer2, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast iac')
call mpi_bcast(crg, natom, mpi_real, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast crg')
call mpi_bcast(cgpatom, natom, mpi_integer4, 0, mpi_comm_world, ierr) !(ai)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast cgpatom')

! cgp
!use mpi_type_create_struct here too
ftype(:) = mpi_integer4 !(ai)
blockcnt(:) = 1
fdisp(1) = 0                            ! integer(ai) iswitch
fdisp(2) = ai                           ! integer(ai) first
fdisp(3) = ai + ai                      ! integer(ai) last
call mpi_type_create_struct(3, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_create_struct')
call mpi_type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_commit')
call mpi_bcast(cgp, ncgp, mpitype_batch, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast cgp')
call mpi_type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_free')

! iaclib
ftype(:) = mpi_real8
blockcnt(1) = 1                                 ! real(8) mass
fdisp(1) = 0
blockcnt(2) = nljtyp                            ! real(8) avdw(nljtyp)
fdisp(2) = 8
blockcnt(3) = nljtyp                            ! real(8) bvdw(nljtyp)
fdisp(3) = 8 + 8*nljtyp
call mpi_type_create_struct(3, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_create_struct')
call mpi_type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_commit')
call mpi_bcast(iaclib, max_atyps, mpitype_batch, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast iaclib')
call mpi_type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_type_free')

! list14 and listex share the same format: logical listxx(max_nbr_range,max_atom)
call mpi_bcast(list14, size(list14), mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast list14')
call mpi_bcast(listex, size(listex), mpi_logical, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast listex')

! list14long and listexlong share the same format: integer(ai) listxxlong(2,max_nxxlong)
call mpi_bcast(list14long, 2*n14long, mpi_integer4, 0, mpi_comm_world, ierr) !(ai)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast list14long')
call mpi_bcast(listexlong, 2*nexlong, mpi_integer4, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast listexlong')

! --- data from the qatom module ---

if (nodeid .eq. 0) write (*,'(80a)') 'qatom data'

! allocate memory
if (nodeid .ne. 0) then
allocate(iqseq(nqat), &
qiac(nqat,nstates), &
qcrg(nqat,nstates), &
qavdw(nqlib,nljtyp), &
qbvdw(nqlib,nljtyp), &
eq(nstates), &
sc_lookup(nqat,natyps+nqat,nstates), &
stat=alloc_status)
call check_alloc('q-atom arrays')
end if

!broadcast sc_lookup(nqat,natyps+nqat,nstates)
call mpi_bcast(sc_lookup, size(sc_lookup), mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast sc_lookup')

! integer(ai) ::  iqseq(nqat)
!change to mpi_type_create  (ai)
call mpi_bcast(iqseq, nqat, mpi_integer4, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast iqseq')

!  integer ::  qiac(nqat,nstates)
call mpi_bcast(qiac, size(qiac), mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qiac')

! real(4) ::  qcrg(nqat,nstates)
call mpi_bcast(qcrg, size(qcrg), mpi_real4, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qcrg')

if(qvdw_flag) then
!mn20030409-> havn't tried with qvdw_flag == .true.
! qavdw and qbvdw share the same format: real(8) qxvdw(nqlib,nljtyp)
call mpi_bcast(qavdw, size(qavdw), mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qavdw')
call mpi_bcast(qbvdw, size(qbvdw), mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast qbvdw')
end if

if (nstates .gt. 0) then
! broadcast eq(:)%lambda
allocate(temp_lambda(1:nstates), stat=alloc_status)
call check_alloc('q-atom energy array')
if (nodeid .eq. 0) temp_lambda(1:nstates) = eq(1:nstates)%lambda
call mpi_bcast(temp_lambda, nstates, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast eq%lambda')
if (nodeid .ne. 0) eq(1:nstates)%lambda = temp_lambda(1:nstates)
deallocate(temp_lambda)
end if

if (nodeid .eq. 0) then 
call centered_heading('end of initiation', '-')
print *
end if

!finally allocate for  slaves:e_send, eq_send
!for master :e_recv,d_recv
call allocate_mpi  

end subroutine init_nodes
#endif


!-----------------------------------------------------------------------

subroutine init_shake
!
! initialize shake constraints
!
!locals
integer                                         ::      mol, b, ia, ja, constr, angle
real(8)                                         :: exclshk
integer                                         ::      src, trg
integer                                         ::      solute_shake_constraints

!allocate molecule list
allocate(shake_mol(nmol), stat=alloc_status)
call check_alloc('shake molecule array')

shake_mol(:)%nconstraints = 0
mol = 0
exclshk = 0.

!count bonds to be constrained in each molecule
!also count shake constraints involving excluded atoms
do b=1,nbonds
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1))
                !new molecule
                mol = mol +1
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
        if((shake_hydrogens .and. (.not. heavy(ia) .or. .not. heavy(ja))) .or. &
           (shake_solute .and. ia <= nat_solute) .or. &
           (shake_solvent .and. ia > nat_solute)) then
           shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1

           if( .not. use_pbc ) then
                if(excl(ia)) exclshk = exclshk + 0.5
                if(excl(ja)) exclshk = exclshk + 0.5
           end if

        end if

end do
!count extra shake constraints from fep file in appropriate molecule
do b = 1, nqshake
    ia=iqshake(b)
        mol = 1
        do while(ia >= istart_mol(mol+1))
                mol = mol + 1
        end do
        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
end do

!allocate bond lists for each molecule
do mol = 1, nmol
        !allocate(sbtemp(nconstr(mol), status = alloc_status)
        allocate(shake_mol(mol)%bond(shake_mol(mol)%nconstraints), stat = alloc_status)
        call check_alloc('shake bond array')
        !shake_mol(mol)%bonds => sbtemp
end do

mol = 0
!add the constraint
do b=1,nbonds
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1)) 
                !new molecule
                mol = mol +1
                shake_mol(mol)%nconstraints = 0
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
if((shake_hydrogens .and. (.not. heavy(ia) .or. .not. heavy(ja))) .or.&
           (shake_solute .and. ia <= nat_solute) .or. &
           (shake_solvent .and. ia > nat_solute)) then
                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%i = ia
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%j = ja
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%dist2 = &
                        bondlib(bnd(b)%cod)%bnd0**2
                !set the bond code to -1 for shaken bonds
                !bnd(b) will be deleted by shrink_topology
                bnd(b)%cod = -1
        end if
end do

!add extra shake constraints from fep file to appropriate molecule
do b = 1, nqshake
    ia=iqshake(b)
        ja=jqshake(b)
        mol = 1
        do while(ia >= istart_mol(mol+1))
                mol = mol + 1
        end do
        !see if already shaken
        do constr = 1, shake_mol(mol)%nconstraints
                if((ia == shake_mol(mol)%bond(constr)%i .and. &
                        ja == shake_mol(mol)%bond(constr)%j) .or. &
                   (ja == shake_mol(mol)%bond(constr)%i .and. &
                        ia == shake_mol(mol)%bond(constr)%j)) then
                        !found it: will overwrite 
                        !also decrement number of constraints
                        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints - 1
                        exit
                end if
        end do
        !constr now contains the right index
        shake_mol(mol)%bond(constr)%i = ia
        shake_mol(mol)%bond(constr)%j = ja
        shake_mol(mol)%bond(constr)%dist2 = &
                dot_product(eq(1:nstates)%lambda,qshake_dist(b,1:nstates))**2
        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
end do

!get total number of shake constraints in solute (used for separate scaling of temperatures)
solute_shake_constraints = sum(shake_mol(1:nmol-nwat)%nconstraints)


!remove molecules with zero constraints from list
trg = 1
src = 2
do while(src <= nmol)
        if(shake_mol(trg)%nconstraints == 0) then
                shake_mol(trg) = shake_mol(src)
                !clear source
                shake_mol(src)%nconstraints = 0
                nullify(shake_mol(src)%bond) 
                src = src + 1
        else
                trg = trg + 1
                if(trg == src) src = src + 1
        end if
end do
shake_molecules = trg

!total number of constraints
shake_constraints = sum(shake_mol(1:shake_molecules)%nconstraints)
write(*,100) shake_constraints
write(*,101) shake_molecules
100     format(/,'number of shake constraints             = ',i10)
101     format('no. molecules with shake constraints    = ',i10)
! calculate #degrees of freedom
ndegf=3*natom-shake_constraints    !changed from ndegf=3*natom-3-shake_constraints, center of mass position is not constrained in the simulation, but is constrained for initial temperatures....
ndegfree=ndegf-3*nexats+exclshk

ndegf_solvent = ndegf - 3*nat_solute + solute_shake_constraints
ndegf_solute = ndegf - ndegf_solvent

ndegfree_solvent = 3*(natom - nat_solute) - (shake_constraints - solute_shake_constraints)
ndegfree_solute = ndegfree - ndegfree_solvent

if (ndegfree_solvent*ndegfree_solute .eq. 0) then    ! if either solvent or solute have 0 degrees of freedom, turn off separate scaling (in case it's on) and do not print detailed temperatures
        detail_temps = .false.
        separate_scaling = .false.
else
        detail_temps = .true.
end if



!clear angles which are shaken (i and k atoms shaken)
do mol=1, shake_molecules
        do constr = 1, shake_mol(mol)%nconstraints
        ia = shake_mol(mol)%bond(constr)%i
        ja = shake_mol(mol)%bond(constr)%j
                do angle = 1, nangles
                        if((ang(angle)%i == ia .and. ang(angle)%k == ja) .or. &
                           (ang(angle)%i == ja .and. ang(angle)%k == ia)) then
                                ang(angle)%cod = 0
                                exit
                        end if
                end do
        end do
end do
end subroutine init_shake

!-----------------------------------------------------------------------

subroutine initial_shaking
!
! initial shaking
!
integer                                         :: niter


xx(:)=x(:)
niter=shake(xx, x)      
write(*,100) 'x', niter
100     format('initial ',a,'-shaking required',i4,&
' interations per molecule on average.')

xx(:)=x(:)-dt*v(:)
niter=shake(x, xx)      
write(*,100) 'v', niter

v(:)=(x(:)-xx(:))/dt


end subroutine initial_shaking

!-----------------------------------------------------------------------
logical function initialize()                  
! local variables
character                                       :: text*80
integer                                         :: i,j,length
real(8)                                         :: stepsize
real(8)                                         :: lamda_tmp(max_states)
integer                                         :: fu, fstat
real(8)                                         ::      rjunk
integer                                         ::      ijunk

! local parameters
integer                                         :: num_args
character(200)                          :: infilename
logical                                         ::      yes
logical                                         ::      need_restart
character(len=80)                       ::      instring
logical                                         ::      inlog
integer                                         ::      mask_rows

! this subroutine will init:
!  nsteps, stepsize, dt
!  temp0, tau_t, iseed, tmaxw
!  use_lrf, nbcycle, rcpp, rcww, rcpw, rcq
!  shake_solute, shake_solvent, shake_hydrogens
! fk_pshell
!  fk_wsphere=-1, wpol_restr, wpol_born
!  fkwpol=-1, dwmz=-1 (values  ized to -1 will
!    be set in water_sphere, once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle, [trj_file], [ene_file]
!  fep_file
!  nstates, eq (allocating memory for eq)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)

! read name of input file from the command line
num_args = command_argument_count()
if (num_args .lt. 1) call die('no input file specified on the command line')
#if defined(cray)
call pxfgetarg(num_args, infilename, 200, i)
#elif defined(mpich)
call getarg(1, infilename)
#else
call getarg(num_args, infilename)
#endif
text = 'reading input from '//infilename
call centered_heading(trim(text), '-')

initialize = .true. 

if(.not. prm_open_section('pbc', infilename)) then
        box = .false.
        write(*,'(a)') 'boundary: sphere'
else
        box = .true.
        write(*,'(a)') 'boundary: periodic box'
        if( .not. prm_get_logical_by_key('rigid_box_centre', rigid_box_centre, .false. ) ) then
                write(*,'(a)') '>>> error: rigid_box_centre must be on or off'
                initialize = .false.
        end if
        write(*,'(a,a3)') 'rigid box centre ', onoff(rigid_box_centre)
        if( .not. prm_get_logical_by_key('constant_pressure', constant_pressure, .false.) ) then
                write(*,'(a)') '>>> error: constant_pressure must be on or off'
                initialize = .false.
        end if

        if( constant_pressure ) then
                write(*,'(a)') 'npt-ensemble'
                volume_try = 0
                volume_acc = 0
                if( .not. prm_get_real8_by_key('max_volume_displ', max_vol_displ) ) then
                        initialize = .false.
                        write(*,'(a)') '>>> error: maximum volume displacement not specified (section pbc)'
                else
                        write(*,5) max_vol_displ
                end if
5       format ('maximum volume displacemet = ', f10.3)

                if( .not. prm_get_integer_by_key('pressure_seed', pressure_seed)) then
                        pressure_seed = 3781
                end if

                                write(*, '(a, i4 )' ) 'pressure seed: ', pressure_seed

                if( .not. prm_get_real8_by_key('pressure', pressure) ) then
                        pressure = 1.0
                end if
                write(*,9) pressure
9       format ('pressure = ',f10.3,'  bar')
                !convert pressure to strange internal unit
                pressure = pressure * 1.43836e-5
        else
                write(*,'(a)') 'nvt-ensemble'
                if( prm_get_line_by_key('control_box', instring) ) then
                        read(instring, *) new_boxl(:)
                        control_box = .true.
                        write(*,'(a, 3f10.3)')'boxsize will be changed to: ', new_boxl
                else
                        control_box = .false.
                end if
        end if !section constant_pressure

        yes = prm_get_logical_by_key('put_solvent_back_in_box', put_solvent_back_in_box)

        yes = prm_get_logical_by_key('put_solute_back_in_box', put_solute_back_in_box)


        if(put_solute_back_in_box .and. put_solvent_back_in_box) then
           write(*,'(a)') 'solute and solvent molecules will be put back in box.'
        else
                if (put_solute_back_in_box) then
                write(*,'(a)') 'only solute molecules will be put back in box.'
                else
                        if (put_solvent_back_in_box) then
                                write(*,'(a)') 'only solvent molecules will be put back in box.'
                        else
                                write(*,'(a)') 'no molecules will be put back in box.'                          
                        end if
                end if
        end if



end if !section pbc


if(.not. prm_open_section('md')) then
        call prm_close
        ! open input file
        fu = freefile()
        open(unit=fu, file=infilename, action='read', form='formatted', status='old', iostat=fstat)
        if (fstat .ne. 0) call die('error opening input file '//infilename)
                initialize = old_initialize(fu)
                close(fu)
        return
end if

need_restart = .false. !flag for restart file required
if(.not. prm_get_integer_by_key('steps', nsteps)) then
        write(*,*) '>>> error: steps not specified (section md)'
        initialize = .false.
end if
if(.not. prm_get_real8_by_key('stepsize', stepsize)) then
        write(*,*) '>>> error: stepsize not specified (section md)'
        initialize = .false.
end if
write (*,10) nsteps, stepsize
10      format ('number of md steps =',i10,'  stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462*stepsize

! --- temperature etc.
if(.not. prm_get_real8_by_key('temperature', temp0)) then
        write(*,*) '>>> error: temperature not specified (section md)'
        initialize = .false.
end if
if(.not. prm_get_real8_by_key('bath_coupling', tau_t)) then
        write(*,*) 'temperature bath relaxation time tau_t set to default'
        tau_t = tau_t_default
end if
write (*,15) temp0,tau_t
tau_t=0.020462*tau_t
if(temp0 <= 0) then
        write(*,'(a)') &
                '>>> error: no dynamics at zero temperature!'
        initialize = .false.
end if
if(tau_t < dt) then
        write(*,'(a)') '>>> error: tau_t must be >= stepsize.'
        initialize = .false.
end if

yes = prm_get_logical_by_key('separate_scaling', separate_scaling, .false.)
if(separate_scaling) then
           write(*,'(a)') 'solute and solvent atoms coupled separately to heat bath.'
else 
           write(*,'(a)') 'solute and solvent atoms coupled together to heat bath.'
end if

15      format ('target temperature =',f10.2,'  t-relax time     =',f10.2)

yes = prm_get_integer_by_key('random_seed', iseed, 1) 
if(.not. prm_get_real8_by_key('initial_temperature', tmaxw)) then
        iseed = 0 !set iseed = 0 if no initial temp
        need_restart = .true.
end if

if (iseed > 0) write (*,16) tmaxw, iseed
16      format ('initial velocities will be generated from maxwell distribution:',&
                /,'maxwell temperature=',f10.2,' random number seed=',i10)

! --- shake, lrf
if(.not. prm_get_logical_by_key('shake_solvent', shake_solvent, .true.)) then
        write(*,'(a)') '>>> error: shake_solvent must be on or off.'
        initialize = .false.
end if
write(*,17) 'all solvent bonds', onoff(shake_solvent)
17      format('shake constaints on ',a,t42,': ',a3)

if(.not. prm_get_logical_by_key('shake_solute', shake_solute, .false.)) then
        write(*,'(a)') '>>> error: shake_solute must be on or off.'
        initialize = .false.
end if 
write(*,17) 'all solute bonds', onoff(shake_solute)

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .false.)) then
        write(*,'(a)') '>>> error: shake_hydrogens must be on or off.'
        initialize = .false.
end if 
write(*,17) 'all bonds to hydrogen', onoff(shake_hydrogens)


       yes = prm_get_logical_by_key('lrf', use_lrf, .false.)
       if(use_lrf) then
               write(*,20) 'lrf taylor expansion outside cut-off'
       else 
               write(*,20) 'standard cut-off'
       end if

20      format ('nonbonded method   = ',a)

yes = prm_get_logical_by_key('force_rms', force_rms, .false.)
if(force_rms) then
        write(*,22) 
end if
22      format ('r.m.s. force will be calculated.')


! --- rcpp, rcww, rcpw, rcq, rclrf
if(.not. prm_open_section('cut-offs')) then
        write(*,'(a)') 'no cut-offs section, default cut-offs used'
        rcpp = rcpp_default
        rcww = rcww_default
        rcpw = rcpw_default
        rcq = rcq_default
        rclrf = rclrf_default
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
        if(use_lrf) then
                if(.not. prm_get_real8_by_key('lrf', rclrf, rclrf_default)) then
                        write(*,'(a)') 'lrf cut-off set to default'
                end if
                if(rclrf < rcpp .or. rclrf < rcpw .or. rclrf < rcww) then
                        write(*,'(a)') &
                                '>>> error; lrf cut-off must not be smaller than solute or solvent cut-offs!'
                        initialize = .false.
                end if
        end if
end if

write (*,25) rcpp,rcww,rcpw,rcq
if(use_lrf) write(*,26) rclrf
25      format ('cut-off radii for non-bonded interactions:',/, &
                'solute-solute:    ',f6.2,/,&
                'solvent-solvent:  ',f6.2,/,&
                'solute-solvent:   ',f6.2,/,&
                'q-atom-non-q-atom:',f6.2)
26      format ('lrf:              ',f6.2)

30      format ('>>> warning: ingnoring obsolete keyword ',a,'.')
! --- simulation sphere

if( .not. box ) then
        if(.not. prm_open_section('sphere')) then
                        fk_pshell = fk_pshell_default
                        print*,'radius of inner restrained shell set to 85% of exclusion shell radius.'
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
                    print*,'radius of inner restrained shell set to 85% of exclusion shell radius.'
                                        rexcl_i = shell_default
                    write(*,50) rexcl_i
                end if
50 format('radius of inner restrained shell       =    ',f8.3) 
                if(.not. prm_get_real8_by_key('shell_force', fk_pshell)) then
                        write(*,'(a)') 'shell force constant set to default'
                        fk_pshell = fk_pshell_default
                end if
                if(fk_pshell > 0) then
                        write(*,47) fk_pshell
                end if
47              format('shell restraint force constant         =',f8.2)

                yes = prm_get_logical_by_key('exclude_bonded', exclude_bonded, .false.)
                if(exclude_bonded) then
                        write(*,'(a)') &
                                'bonded interactions outside the sphere will be eliminated'
                end if
       end if

        ! --- solvent 
        inlog = prm_open_section('solvent')
        if(.not. inlog) inlog = prm_open_section('water') !try also the old name
        if(.not. inlog) then       !defaults
                fk_wsphere = -1
                dwmz = -1
                awmz = -1
                wpol_restr = wpol_restr_default
                wpol_born = wpol_restr_default
                fkwpol = -1 
        else
                if(prm_get_real8_by_key('radius', rwat_in)) then
                        write(*,'(a,f8.2)') 'target solvent radius =',rwat_in
                end if
                if(prm_get_line_by_key('centre', instring)) then
                        write(*,30) 'centre'
                end if
                if(prm_get_real8_by_key('pack', rjunk)) then
                        write(*,30) 'pack'
                end if


          if(.not. prm_get_real8_by_key('radial_force', fk_wsphere)) then
                write(*,'(a)') 'solvent radial restraint force constant set to default'
                fk_wsphere = -1 ! this will be set in water_sphere, once target radius is known
          end if
          yes=prm_get_logical_by_key('polarisation', wpol_restr, wpol_restr_default)
          !default is on when pol. restr is on, otherwise off
          yes=prm_get_logical_by_key('charge_correction', wpol_born, wpol_restr)
          if(wpol_born .and. .not. wpol_restr) then
                write(*,'(a)') '>>> error: charge_correction on requires polarisation on (section solvent)'
                initialize = .false.
          end if
          if(.not. prm_get_real8_by_key('polarisation_force', fkwpol)) then
                write(*,'(a)') 'solvent polarisation force constant set to default'
                fkwpol = -1 ! this will be set in water_sphere, once target radius is known
          end if
          yes = prm_get_real8_by_key('morse_depth', dwmz, -1._8)                        
          yes = prm_get_real8_by_key('morse_width', awmz, -1._8)                        
          if(prm_get_string_by_key('model', instring)) then
                write(*,30) 'model'
          end if
        end if !if (.not. inlog)
end if !if( .not. box )


if(.not. prm_open_section('intervals')) then
        write(*,'(a)') 'non-bond list update interval set to default.'
        nbcycle = nb_cycle_default
        write(*,'(a)') 'energy summary interval set to default.'
        iout_cycle = iout_cycle_default
        itemp_cycle = iout_cycle_default
        iene_cycle = 0 !no energy
        itrj_cycle = 0 !no trajectory

        ivolume_cycle = ivolume_cycle_default


else
        if(.not. prm_get_integer_by_key('non_bond', nbcycle)) then
                write(*,'(a)') 'non-bond list update interval set to default.'
                nbcycle = nb_cycle_default
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

write(*,84) nbcycle
84      format('non-bonded pair list update interval   =',i8)
86      format('energy summary print-out interval      =',i8)
87      format('temperature print-out interval         =',i8)
88      format('trajectory write interval              =',i8)
89      format('energy file write interval             =',i8)
83  format('volume change interval                 =',i8)

if(iout_cycle > 0) then
        write (*,86) iout_cycle
else
        write(*,'(a)') 'no energy summaries written.'
        iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if
if(itemp_cycle > 0) then
        write (*,87) itemp_cycle
else
        write(*,'(a)') 'no temperatures written.'
        itemp_cycle = -999999999 ! make sure mod(istep, itemp_cycle) never = 0
end if
if(itrj_cycle > 0) then
        write (*,88) itrj_cycle
else
        itrj_cycle = -999999999 !no energy
        write(*,'(a)') 'no trajectory written.'
end if
if(iene_cycle > 0) then
        write (*,89) iene_cycle
else
        iene_cycle = -999999999 !no energy
        write(*,'(a)') 'no energy file written.'
end if
if( constant_pressure ) then
        write(*,83) ivolume_cycle
end if

!read trajectory atom mask
mask_rows = prm_count('trajectory_atoms')
if(itrj_cycle > 0) then
        if(mask_rows == 0) then
                write(*,'(a)') 'all atoms will be included in the trajectory.'
                yes = trj_store_mask('all')
        else
                do i=1,mask_rows
                        yes = prm_get_line(text)
                        yes = trj_store_mask(text)
                end do
        end if
elseif(mask_rows == 0) then
        write(*,'(a)') 'ignoring section trajectory_atoms.'
end if

if(.not. prm_open_section('files')) then
        write(*,'(a)') '>>> error: files section not found.'
        initialize = .false.
else
        if(.not. prm_get_string_by_key('topology', top_file)) then
                write(*,'(a)') '>>> error: topology not specified (section files)'
                initialize = .false.
        end if
        write (*,60) trim(top_file)
60              format ('topology file      = ',a)

        if(.not. prm_get_string_by_key('restart', restart_file)) then
                restart = .false.
                if(need_restart) then
                        write(*,'(a)') '>>> error: restart file required when initial temp. not given.'
                        initialize = .false.
                end if
        else
                restart = .true.
        end if

        if(restart) then
                write (*,65) trim(restart_file)
        else
                write (*,'(a)') 'initial coordinates taken from topology.'
                if(iseed == 0) then
                        write(*,'(a)') &
                                '>>> error: need a random number seed to generate initial velocities, aborting.'
                        initialize = .false.
                end if
        end if
65              format ('initial coord. file= ',a)

        if(.not. prm_get_string_by_key('final', xfin_file)) then
                write(*,'(a)') '>>> error: final co-ordinate file not specified (section files, keyword final)'
                initialize = .false.
        end if
        write (*,80) trim(xfin_file)
80              format ('final coord. file  = ',a)

        if(.not. prm_get_string_by_key('trajectory', trj_file)) then
                if(itrj_cycle > 0) then
                        write(*,'(a)') '>>> error: trajectory file name required to write trajectory!'
                        initialize = .false.
                end if
        else
                if(itrj_cycle < 0) then
                        write(*,*) '>>> error: trajectory file given but no output interval'
                        initialize = .false.
                end if
                if(itrj_cycle > 0) write (*,90) trim(trj_file)
        end if
90              format ('trajectory file    = ',a)

        if(.not. prm_get_string_by_key('energy', ene_file)) then
                if(iene_cycle > 0) then
                        write(*,'(a)') '>>> error: energy file name required to write energies!'
                        initialize = .false.
                end if
        else

                if(iene_cycle < 0) then

                        write(*,'(a)') '>>> error: energy file given but no energy interval'

                        initialize=.false.

                end if
                if(iene_cycle > 0) write (*,94) trim(ene_file)
        end if
94              format ('energy output file = ',a)

        if(.not. prm_get_string_by_key('fep', fep_file)) then
                write(*,'(a)') 'no fep file.'
                !initialize = .false. !this condition is ok.
                fep_file = ''
        else
                write (*,95) trim(fep_file)
95                      format ('fep input file     = ',a,/)
        end if
        if(.not. prm_get_string_by_key('restraint', exrstr_file)) then
                implicit_rstr_from_file = 0
        else
                implicit_rstr_from_file = 1
                write (*,104) trim(exrstr_file)
104                     format ('external rstr file = ',a,/)
        end if
        if(prm_get_string_by_key('water', instring)) then
                write(*,30) 'water'
        end if
end if                  

! --- states, eq
nstates = 0
if(prm_open_section('lambdas')) then
        do while(prm_get_field(instring))
                nstates = nstates + 1
                read(instring, *, iostat=fstat) lamda_tmp(nstates)
                if(fstat /= 0) then
                        write(*,'(a)') '>>> error: invalid lambda value.'
                        initialize = .false.
                        exit
                end if
        end do
end if
if(nstates == 0 .and. fep_file /= '') then
        if(fep_file /= '') then
                write(*,'(a)') 'defaulting to single fep state.'
                nstates = 1
                lamda_tmp(1) = 1.
        end if
end if
if(nstates > 0 ) then
        if(fep_file == '') then
                write(*,'(a)') '>>> error: fep file required to use lambdas!'
                initialize = .false.
        else
                ! allocate memory for eq
                allocate(eq(nstates), stat=alloc_status)
                call check_alloc('q-atom energy array')

                ! init eq%lambda
                eq(1:nstates)%lambda = lamda_tmp(1:nstates)
                write (*,98) (eq(i)%lambda,i=1,nstates)
98                      format ('lambda-values      = ',10f8.5)
        end if
end if

!       --- restraints:
write (*,'(/,a)') 'listing of restraining data:'

! --- nrstr_seq, [rstseq]
nrstr_seq = prm_count('sequence_restraints')
109 format (/,'no. of sequence restraints =',i10)
if ( nrstr_seq .gt. 0 ) then
        ! allocate memory for rstseq
        write (*,109) nrstr_seq
        allocate(rstseq(nrstr_seq), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,110)
110             format ('  atom_i  atom_j      fc  h-flag to_centre')
        do i=1,nrstr_seq
                ! read rstseq(i)
                yes = prm_get_line(text)
                rstseq(i)%to_centre = 0 
                read(text,*, end=111, err=111) rstseq(i)
111                     write(*,112) rstseq(i)
112                     format (2i8,f8.2,i8,i10)
  end do
end if

! --- nrstr_pos, [rstpos]
nrstr_pos = prm_count('atom_restraints')
115 format (/,'no. of position restratints =',i10)
if ( nrstr_pos .gt. 0 ) then
        write (*,115) nrstr_pos
        ! allocate memory for rstpos
        allocate(rstpos(nrstr_pos), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,120)
120             format ('atom_i      x0      y0      z0     fcx     fcy     fcz   state')
        do i=1,nrstr_pos ! read rstpos(i)
                yes = prm_get_line(text)
                read(text,*, iostat=fstat) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
                        (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
                if(fstat /= 0) then
                        write(*,'(a)') '>>> error: invalid atom restraint data.'
                        initialize = .false.
                        exit
                end if
                write (*,122) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
                        (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
        end do
122             format (i6,6f8.2,i8)
end if

! --- nrstr_dist, [rstdis]
nrstr_dist = prm_count('distance_restraints')
125     format (/,'no. of distance restraints =',i10)
if ( nrstr_dist .gt. 0 ) then
        write (*,125) nrstr_dist
        ! allocate memory for rstdis
        allocate(rstdis(nrstr_dist), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,130)
130             format ('atom_i atom_j   dist1   dist2   fc        state') 
        do i=1,nrstr_dist
                yes=prm_get_line(text)
                ! read rstdis(i)
                if(scan(text, ':') > 0) then !got res:atnr
                  !store in i&j as res:atnr and assign atom nr after topology is read (prep_coord)
                  read(text,*, iostat=fstat) rstdis(i)%itext,rstdis(i)%jtext,rstdis(i)%d1,& 
                     rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
                else !plain numbers
                  read(text,*, iostat=fstat) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,&
                        rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
                  rstdis(i)%itext = 'nil'
                  rstdis(i)%jtext = 'nil'
                end if
                if(fstat /= 0) then
                  write(*,'(a)') '>>> error: invalid distance restraint data.'
                  initialize = .false.
                  exit
                end if
                write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%d2,rstdis(i)%fk, &
                        rstdis(i)%ipsi
        end do
132             format (i6,1x,i6,3f8.2,i8)
end if


if (.not. box )then
! --- nrstr_wall, [rstwal]
nrstr_wall = prm_count('wall_restraints')
135 format (/,'no. of wall sequence restraints=',i10)
if ( nrstr_wall .gt. 0) then
        write (*,135) nrstr_wall
        ! allocate memory for rstwal
        allocate(rstwal(nrstr_wall), stat=alloc_status)
        call check_alloc('restraint list')
        write (*,140)
140             format ('atom_i atom_j   dist.      fc  amorse  dmorse  h-flag')
        do i=1,nrstr_wall
                ! read rstwal(:)
                yes = prm_get_line(text)
                read(text,*, iostat=fstat) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
                        rstwal(i)%amorse, rstwal(i)%dmorse, rstwal(i)%ih
                if(fstat /= 0) then
                        write(*,'(a)') '>>> error: invalid wall restraint data.'
                        initialize = .false.
                        exit
                end if
                write (*,142) rstwal(i)     
        end do
142             format (i6,1x,i6,4f8.2,i8)
end if
end if

call prm_close
end function initialize


!-----------------------------------------------------------------------------------


logical function old_initialize(fu)
!arguments
integer                                         ::      fu
! local variables
integer                                         ::      iuse_indip, shake_flag
character                                       :: text*80, watmodel*80
integer                                         :: i,j,length
integer                                         :: irestart
real(8)                                         :: stepsize
real(8)                                         :: lamda_tmp(max_states)
integer                                         :: fstat
integer                                         ::      nbmethod
integer                                         ::      iwpol_restr
real(8)                                         ::      rjunk

!this is called by initialize to read old-style input file which is 
!alreadu open as unit fu

! this subroutine will init:
!  nsteps, stepsize, dt
!  temp0, tau_t, iseed, tmaxw
!  usr_lrf, nbcycle, rcpp, rcww, rcpw, rcq
!  shake_solvent, shake_solute, shake_hydrogens
!  fk_pshell


!  fk_wsphere=-1, wpol_restr, wpol_born fkwpol=-1, dwmz=-1, awmz=-1
!    (values initialized to -1 will be set in water_sphere, 
!    once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle [trj_file], [ene_file]
!  fep_file
!  nstates, eq (allocating memory for eq)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)


write(*,1)
1       format('>>> warning: entering unsupported compatibility mode',/ &
           '             to read version 2 input file.',/,&
           '             new features unavailable.')

old_initialize = .true. 

!use default values for new features not in old kind of input.
rclrf = 999.
exclude_bonded = .false.
force_rms = .false.
shake_hydrogens = .false.
itemp_cycle = iout_cycle_default
awmz = -1

! --- nsteps, stepsize
!read (fu,*, iostat=stat) nsteps,stepsize
if (.not. prm_get_int_real8(nsteps,stepsize)) then
        old_initialize = .false.
                call die("wrong input format.")
end if

write (*,10) nsteps, stepsize
10      format ('number of md steps =',i10,'  stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462*stepsize

! --- temp0, tau_t, iseed, tmaxw
read(fu,'(a80)') text !read line into buffer
!now read buffer (avoid reading more lines from inpunt in search for more values)
read(text,*, err=17, end=17) temp0,tau_t, iseed,tmaxw 
17      write (*,15) temp0,tau_t
if(temp0 <= 0) then
        write(*,'(a)') &
                '>>> error: no dynamics at zero temperature! aborting.'
        old_initialize = .false.
end if

15      format ('target temperature =',f10.2,'  t-relax time     =',f10.2)
if (iseed > 0) write (*,16) tmaxw, iseed
16      format ('initial velocities will be generated from maxwell distribution:',&
                /,'maxwell temperature=',f10.2,' random number seed=',i10)
tau_t=0.020462*tau_t
if(tau_t < dt) then
        write(*,'(a)') '>>> error: tau_t must be >= stepsize.'
        old_initialize = .false.
end if

! --- nbmethod, nbcycle, rcpp, rcww, rcpw, rcq
read (fu,*) nbmethod,nbcycle,rcpp,rcww,rcpw,rcq
if(nbmethod == 2) then
        use_lrf = .true.
else 
        use_lrf = .false.
end if
write (*,20) nbmethod,nbcycle
20      format ('nonbonded method   =',i10,'  nb update cycle  =',i10,/)
write (*,25) rcpp,rcww,rcpw,rcq
25      format ('cutoffs are: rcpp  =',f6.2,'  rcww =',f6.2,'  rcpw =', &
f6.2,'  rcqp =',f6.2,/)

! --- shake_flag
read (fu,*) shake_flag
shake_solvent = .false.
shake_solute = .false.
if(shake_flag >= 1) then
        shake_solvent = .true.
end if
if(shake_flag == 2) then
        shake_solute = .true.
end if

write (*,30) shake_flag
30      format ('shake method       =',i10)

! --- iuse_indip
read (fu,*) 
write (*,35) 
35      format ('ignoring induced dipole flag.')

! --- protein center: xpcent(:)
read (fu,*) 
write (*,40)
40      format ('ignoring solute centre.')

! --- rexcl_o, rexcl_i, fk_pshell
read (fu,*) rjunk, rjunk, fk_pshell
write(*,44)
write (*,45) fk_pshell
44      format ('ignoring exclusion and shell radii.')
45      format ('restrained shell force const.          =',f8.2)

! --- water center: xwcent(:)
read (fu,*) 
write (*,50) 
50      format ('ignoring solvent centre.')

! set default values before reading
! done this way because the sgi compiler initialises values to be read to zero

read(fu,'(a80)') text ! read line into buffer
! now read buffer (avoid reading more lines from input in search for more values)
read(text, fmt=*, err=58, end=58) rjunk, rjunk, fk_wsphere, iwpol_restr, fkwpol, dwmz
goto 59

! set default values:
58      dwmz = -1
if (fkwpol .eq. 0) then
  fkwpol = -1
        if (fk_wsphere .eq. 0) then
          fk_wsphere = -1
  end if
end if
59      if(iwpol_restr == 0) then
        wpol_restr = .false.
        wpol_born = .false.
elseif(iwpol_restr == 1) then
        wpol_restr = .true.
        wpol_born = .true.
elseif(iwpol_restr == 2) then
        wpol_restr = .true.
        wpol_born = .false.
else
        call die('unknown water polarisation restraining mode')
end if
write(*,57)
57      format('ignoring solvent radius and min. packing distance.')
! --- top_file
read (fu,'(a80)') text
call get_fname (text,length,top_file)
write (*,60) top_file(1:length)
60      format ('topology file      = ',a)

! --- restart, [restart_file]
read (fu,*) irestart
if ( irestart .eq. 1 ) then
        restart = .true.
        read (fu,'(a80)') text
        call get_fname (text,length,restart_file)
        write (*,65) restart_file(1:length)
else
        restart = .false.
        write (*,'(a)') 'initial coordinates taken from topology.'
        if(iseed == 0) then
                write(*,'(a)') &
                        'error: need a random number seed to generate initial velocities, aborting.'
                call die('invalid data in input')
        end if
end if
65 format ('initial coord. file= ',a)

! --- xfin_file
read (fu,'(a80)') text
call get_fname (text,length,xfin_file)
write (*,80) xfin_file(1:length)
80 format ('final coord. file  = ',a,/)

! --- itrj_cycle, iene_cycle, iout_cycle, [trj_file], [ene_file]
read (fu,*) itrj_cycle, iene_cycle, iout_cycle
write (*,85) itrj_cycle, iene_cycle, iout_cycle
85 format ('trajectory, energy and output cycles   =',3i8,/)

if ( itrj_cycle .gt. 0 ) then
read (fu,'(a80)') text
call get_fname (text,length,trj_file)
write (*,90) trj_file(1:length)
else
write (*,'(a)') 'no trajectory written.'
 itrj_cycle = -999999999 !make sure mod(istep, itrj_cycle) never = 0
end if
90 format ('trajectory file    = ',a)

if ( iene_cycle .gt. 0 ) then
read (fu,'(a80)') text
call get_fname (text,length,ene_file)
write (*,94) ene_file(1:length)
else
write (*,'(a)') 'no energy file written'
 iene_cycle = -999999999 ! make sure mod(istep, iene_cycle) never = 0
end if
94 format ('energy output file = ',a)

if(iout_cycle == 0) then
write(*,'(a)') 'no energy summaries written.'
iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if

! --- fep_file
read (fu,'(a80)') text
call get_fname (text,length,fep_file)
write (*,95) fep_file(1:length)
95 format ('fep input file     = ',a,/)

! --- nstates, eq
read (fu,*) nstates, (lamda_tmp(i),i=1,nstates)
if ( nstates .gt. 0 ) then
! allocate memory for eq
allocate(eq(nstates), stat=alloc_status)
call check_alloc('q-atom energy array')

! init eq%lambda
eq(1:nstates)%lambda = lamda_tmp(1:nstates)
write (*,98) (eq(i)%lambda,i=1,nstates)
98  format ('lambda-values      = ',10f8.5)
end if

!       --- restraints:
write (*,'(/,a)') 'listing of restraining data:'

! --- implicit_rstr_from_file, [exrstr_file]
read (fu,*) implicit_rstr_from_file
write (*,101) implicit_rstr_from_file
101 format ('read rstr file     =',i10)
if ( implicit_rstr_from_file .eq. 1 ) then
read (fu,'(a80)') text
call get_fname (text,length,exrstr_file)
write (*,104) exrstr_file(1:length)
else
write (*,105)
end if
104 format ('external rstr file = ',a,/)
105 format ('implicit positional restraints from topology.',/)

! --- nrstr_seq, [rstseq]
read (fu,*) nrstr_seq
write (*,109) nrstr_seq
109 format (/,'no. sequence rstrs =',i10)
if ( nrstr_seq .gt. 0 ) then
! allocate memory for rstseq
allocate(rstseq(nrstr_seq), stat=alloc_status)
call check_alloc('restraint list')
write (*,110)
110 format (1x,'  atom_i  atom_j      fc  h-flag to_centre')
end if
do i=1,nrstr_seq
! read rstseq(i)
read (fu,'(a80)') text
rstseq(i)%to_centre = 0 
read(text,*, end=111, err=111) rstseq(i)
111     write(*,112) rstseq(i)
112 format (2i8,f8.2,i8,i10)
end do

! --- nrstr_pos, [rstpos]
read (fu,*) nrstr_pos
write (*,115) nrstr_pos
115 format (/,'no. position rstrs =',i10)
if ( nrstr_pos .gt. 0 ) then
! allocate memory for rstpos
allocate(rstpos(nrstr_pos), stat=alloc_status)
call check_alloc('restraint list')
write (*,120)
end if
120 format ('atom_i      x0      y0      z0     fcx     fcy     fcz  istate')
do i=1,nrstr_pos
! read rstpos(i)
read (fu,*) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
  (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
write (*,122) rstpos(i)%i,(rstpos(i)%x(j),j=1,3), &
  (rstpos(i)%fk(j),j=1,3), rstpos(i)%ipsi
end do
122 format (i6,6f8.2,i8)

! --- nrstr_dist, [rstdis]
read (fu,*) nrstr_dist
write (*,125) nrstr_dist
125 format ('no. distance rstrs =',i10)
if ( nrstr_dist .gt. 0 ) then
! allocate memory for rstdis
allocate(rstdis(nrstr_dist), stat=alloc_status)
call check_alloc('restraint list')
write (*,130)
end if
130 format ('atom_i atom_j   dist.      fc  istate') 
do i=1,nrstr_dist
! read rstdis(i)
read (fu,*) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%fk, &
  rstdis(i)%ipsi
        rstdis(i)%d2 = rstdis(i)%d1 !no flat-bottom
write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%fk, &
  rstdis(i)%ipsi
end do
132 format (i6,1x,i6,2f8.2,i8)

! --- nrstr_wall, [rstwal]
read (fu,*) nrstr_wall
write (*,135) nrstr_wall
135 format ('no. wall seq. rstrs=',i10)
if ( nrstr_wall .gt. 0) then


! allocate memory for rstwal
allocate(rstwal(nrstr_wall), stat=alloc_status)
call check_alloc('restraint list')
write (*,140)
end if
140 format ('atom_i atom_j   dist.      fc  h-flag')
do i=1,nrstr_wall
! read rstwal(:)
read (fu,*) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
  rstwal(i)%ih
write (*,142) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
  rstwal(i)%ih     
end do
142 format (i6,1x,i6,2f8.2,i8)

read (fu,'(a80)') text
write (*,157) 
157 format ('ignoring water file.')

! --- determine water model
read (fu,'(a80)') text
write (*,160) 
160 format ('ignoring water model.')

end function old_initialize

!-----------------------------------------------------------------------

subroutine lrf_taylor
! *** local variables
integer                                         ::      i,i3,ic
real(8)                                         ::      vij, q
real(8)                                         ::      dr(3),df(3)

! global variables used:
!  e%lrf, natom, excl, iqatom, iwhich_cgp, lrf, x, crg, d

do i = 1, natom
! for every atom:

if ( ( use_pbc .and. (iqatom(i)==0) ) .or. ( (.not. excl(i) ) .and. (iqatom(i)==0) ) ) then
  ! unless excluded atom or q-atom:

  ! find the displacement dr from the center of the charge group
i3 = i*3-3
ic = iwhich_cgp(i)
dr(1) = lrf(ic)%cgp_cent(1) - x(i3+1)
dr(2) = lrf(ic)%cgp_cent(2) - x(i3+2)
dr(3) = lrf(ic)%cgp_cent(3) - x(i3+3)

! --- electric potential
vij=lrf(ic)%phi0 &
     +lrf(ic)%phi1(1)*dr(1)+lrf(ic)%phi1(2)*dr(2)+lrf(ic)%phi1(3)*dr(3) &
     +0.5*(lrf(ic)%phi2(1)*dr(1)+lrf(ic)%phi2(2)*dr(2) &
     +lrf(ic)%phi2(3)*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi2(4)*dr(1)+lrf(ic)%phi2(5)*dr(2) &
     +lrf(ic)%phi2(6)*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi2(7)*dr(1)+lrf(ic)%phi2(8)*dr(2) &
     +lrf(ic)%phi2(9)*dr(3))*dr(3)

e%lrf = e%lrf + .5 * crg(i) * vij

! --- electric field
df(1)=lrf(ic)%phi1(1) &
     +lrf(ic)%phi2(1)*dr(1)+lrf(ic)%phi2(2)*dr(2)+lrf(ic)%phi2(3)*dr(3) &
     +0.5*(lrf(ic)%phi3(1 )*dr(1)+lrf(ic)%phi3(2 )*dr(2) &
     +lrf(ic)%phi3(3 )*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi3(4 )*dr(1)+lrf(ic)%phi3(5 )*dr(2) &
     +lrf(ic)%phi3(6 )*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi3(7 )*dr(1)+lrf(ic)%phi3(8 )*dr(2) &
     +lrf(ic)%phi3(9 )*dr(3))*dr(3)
df(2)=lrf(ic)%phi1(2) &
     +lrf(ic)%phi2(4)*dr(1)+lrf(ic)%phi2(5)*dr(2)+lrf(ic)%phi2(6)*dr(3) &
     +0.5*(lrf(ic)%phi3(10)*dr(1)+lrf(ic)%phi3(11)*dr(2) &
     +lrf(ic)%phi3(12)*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi3(13)*dr(1)+lrf(ic)%phi3(14)*dr(2) &
     +lrf(ic)%phi3(15)*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi3(16)*dr(1)+lrf(ic)%phi3(17)*dr(2) &
     +lrf(ic)%phi3(18)*dr(3))*dr(3)
df(3)=lrf(ic)%phi1(3) &
     +lrf(ic)%phi2(7)*dr(1)+lrf(ic)%phi2(8)*dr(2)+lrf(ic)%phi2(9)*dr(3) &
     +0.5*(lrf(ic)%phi3(19)*dr(1)+lrf(ic)%phi3(20)*dr(2) &
     +lrf(ic)%phi3(21)*dr(3))*dr(1) &
     +0.5*(lrf(ic)%phi3(22)*dr(1)+lrf(ic)%phi3(23)*dr(2) &
     +lrf(ic)%phi3(24)*dr(3))*dr(2) &
     +0.5*(lrf(ic)%phi3(25)*dr(1)+lrf(ic)%phi3(26)*dr(2) &
     +lrf(ic)%phi3(27)*dr(3))*dr(3) 

  ! update d
d(i3+1)=d(i3+1)-crg(i)*df(1)
d(i3+2)=d(i3+2)-crg(i)*df(2)
d(i3+3)=d(i3+3)-crg(i)*df(3)
end if
end do
end subroutine lrf_taylor


!-----------------------------------------------------------------------

subroutine make_pair_lists
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

if( use_pbc ) then
        if (.not. use_lrf)then
                !cutoff
           if(iuse_switch_atom == 1) then
                call nbpplist_box
                call nbpwlist_box
                call nbqplist_box
            else
                call nbpplis2_box
                call nbpwlis2_box
                call nbqplis2_box
            end if
            call nbwwlist_box
        else
            call cgp_centers
                              if ( iuse_switch_atom == 1 ) then 
                                call nbpplist_box_lrf
                                      call nbpwlist_box_lrf
                                      call nbqplist_box
                              else
                                call nbpplis2_box_lrf
                                      call nbpwlis2_box_lrf
                                      call nbqplis2_box
            endif
            call nbwwlist_box_lrf
        endif
                call nbqwlist_box
else !spherical case
        if(.not. use_lrf) then
                ! cutoff
                if( iuse_switch_atom .eq. 1 ) then
                        call nbpplist
                        call nbpwlist
                        call nbqplist
                else
                        call nbpplis2
                        call nbpwlis2
                        call nbqplis2
                end if
                call nbwwlist
        else 
                ! cutoff with lrf
                call cgp_centers ! *** m�ste anropas av alla noder (nollst�ller lrf)
                if( iuse_switch_atom .eq. 1 ) then
                        call nbpplist_lrf   
                        call nbpwlist_lrf   
                        call nbqplist       
                else
                        call nbpplis2_lrf
                        call nbpwlis2_lrf
                        call nbqplis2
                end if
                call nbwwlist_lrf   
        end if

        call nbqwlist   
end if

#if defined (profiling)
profile(1)%time = profile(1)%time + rtime() - start_loop_time
#endif

end subroutine make_pair_lists

!-----------------------------------------------------------------------

subroutine maxwell
! *** local variables
integer                                         :: i,j,k
real(8)                                         :: zero,sd,vg,kt

!       generate maxwellian velocities 
zero  = 0.0
kt = boltz*tmaxw

do i=1,natom
        sd = sqrt (kt/iaclib(iac(i))%mass)
        do j=1,3
                call gauss (zero,sd,vg,iseed)
                k=(i-1)*3+j
                v(k)=vg
        end do
end do

end subroutine maxwell

!-----------------------------------------------------------------------



subroutine temperature(temp,tscale_solute,tscale_solvent,ekinmax)
! calculate the temperature
!arguments
real(8)                                         :: temp,tscale_solute,tscale_solvent,ekinmax

!locals
integer                                         ::      i, i3
real(8)                                         ::      ekin

temp = 0.
temp_solute = 0.
tfree_solute = 0.
texcl_solute = 0.

!get kinetic energies for solute atoms
do i=1,nat_solute
        i3=i*3-3
        ekin = 0.5*iaclib(iac(i))%mass*(v(i3+1)**2+v(i3+2)**2+v(i3+3)**2)
        temp_solute = temp_solute + ekin

        !******pwadded if
        if( use_pbc .or. ( (.not. use_pbc) .and. (.not. excl(i)) ) ) then
                tfree_solute = tfree_solute +ekin
        else
                texcl_solute = texcl_solute +ekin
        end if
        !if ( .not. excl(i)) tfree = tfree + ekin
        if ( ekin .gt. ekinmax ) then
                ! hot atom warning
                write (*,180) i,2.*ekin/boltz/3.
        end if
end do

tfree_solvent = 0.
temp_solvent = 0.
texcl_solvent = 0.
ekin = 0


!get kinetic energies for solvent atoms
do i=nat_solute+1,natom
        i3=i*3-3
        ekin = 0.5*iaclib(iac(i))%mass*(v(i3+1)**2+v(i3+2)**2+v(i3+3)**2)
        temp_solvent = temp_solvent + ekin

        !******pwadded if
        if( use_pbc .or. ( (.not. use_pbc) .and. (.not. excl(i)) ) ) then
                tfree_solvent = tfree_solvent +ekin
        else
                texcl_solvent = texcl_solvent +ekin
        end if
        !if ( .not. excl(i)) tfree = tfree + ekin
        if ( ekin .gt. ekinmax ) then
                ! hot atom warning
                write (*,180) i,2.*ekin/boltz/3.
        end if
end do

tfree = tfree_solvent + tfree_solute
temp = temp_solute + temp_solvent

e%kinetic = temp

temp  = 2.0*temp/boltz/real(ndegf)
tfree = 2.0*tfree/boltz/real(ndegfree)

if (detail_temps) then
        temp_solute  = 2.0*temp_solute /boltz/real(ndegf_solute)
        tfree_solute = 2.0*tfree_solute/boltz/real(ndegfree_solute)
        if ( ndegf_solute .ne. ndegfree_solute) texcl_solute = 2.0*texcl_solute/boltz/real(ndegf_solute - ndegfree_solute)

        temp_solvent  = 2.0*temp_solvent /boltz/real(ndegf_solvent)
        tfree_solvent = 2.0*tfree_solvent/boltz/real(ndegfree_solvent)
        if ( ndegf_solvent .ne. ndegfree_solvent) texcl_solvent = 2.0*texcl_solvent/boltz/real(ndegf_solvent - ndegfree_solvent)
end if


if (separate_scaling) then
        if ( tfree_solvent .ne. 0 ) tscale_solvent = temp0/tfree_solvent - 1.0
        tscale_solvent = sqrt ( 1 + dt/tau_t * tscale_solvent )
        if ( tfree_solute .ne. 0 ) tscale_solute = temp0/tfree_solute - 1.0
        tscale_solute = sqrt ( 1 + dt/tau_t * tscale_solute )
else
        if ( tfree .ne. 0 ) tscale_solvent = temp0/tfree - 1.0
        tscale_solvent = sqrt ( 1 + dt/tau_t * tscale_solvent )
        tscale_solute = tscale_solvent
end if


180 format ('>>> warning: hot atom, i =',i10,' temp(i)=',f10.2)

end subroutine temperature

!-----------------------------------------------------------------------
!******pwchanged 2002-10-01
subroutine md_run

! local variables
integer                         :: i,j,k,niter

integer                         :: i3
real(8)                         :: temp,tlast
real(8)                         ::ekinmax
real(8)                         ::tscale_solute,tscale_solvent

real(8)                         :: time0, time1, time_per_step, startloop
integer(4)                              :: time_completion

#if defined(profiling)
real(8)                                         :: start_loop_time1, start_loop_time2

profile(1)%name = 'nb_update'
profile(2)%name = '   nbwwlist_time'
profile(3)%name = '   nbpplist_time'
profile(4)%name = '   nbpwlist_time'
profile(5)%name = '   nbqplist_time'
profile(6)%name = '   nbqwlist_time'
profile(7)%name = 'shake'
profile(8)%name = 'bonded terms'
profile(9)%name = 'restraints'
profile(10)%name = 'nonbonded terms'
profile(11)%name = 'update vel. & coords.'


#endif

#if defined(profiling)
#if defined(use_mpi)
if (nodeid .eq. 0) then
        allocate(all_node_times(num_profiling_times*numnodes), stat=alloc_status) !vector for storing all node's node_times, used by mpi_gather at end of md_run
        call check_alloc('mpi profiling')
end if
allocate(node_times(num_profiling_times), stat=alloc_status) !each node's profiling times, used at end of md_run by mpi_gather
call check_alloc('mpi profiling')

all_node_times(:) = 0.0
node_times(:) = 0.0

#endif
#endif


!define number of coord to send/recieve
nat3=natom*3

! calculate maximum temperature
!**mn-> only master calc. temp for now.
if (nodeid .eq. 0) then
ekinmax = 1000.0*ndegf*boltz*temp0/2.0/real(natom)

call temperature(temp,tscale_solute,tscale_solvent,ekinmax)
!store old temp
tlast = temp
end if

if (nodeid .eq. 0) then
        ! master node only: print initial temperatures
        write (*,*)
        write (*,120) 'initial', temp, tfree
                if ( detail_temps ) then
                        write (*,120) 'solvent', temp_solvent, tfree_solvent
                        write (*,120) 'solute', temp_solute, tfree_solute
!                       write (*,120) 'excl solute, solvent', texcl_solute, texcl_solvent
                end if
120             format(a7,' temperatures are : ttot =',f10.2,' tfree =',f10.2)
        write (*,*)

        ! init timer
        time0 = rtime()

                ! init timer of total loop time
                startloop = rtime()
end if


!***********************************************************************
!       begin main dynamics loop (verlet leap-frog algorithm)
!***********************************************************************

! no loop (only calc. energies) if compiling with the dum flag
#ifndef dum
do istep = 0, nsteps-1
#endif
        if ( mod(istep, nbcycle) .eq. 0 ) then


                ! every nbcycle steps:

                                !put molecules back in box for nice visualisation, needs to be here to prevent problems with lrf
                                !update cgp_centers for lrf
                                !only call put_back_in_box if using pbc and either solute or solvent should be put back in box
                                if( use_pbc .and. (put_solute_back_in_box .or. put_solvent_back_in_box) ) then 
                                  call put_back_in_box() 
                                end if

                if ((nodeid .eq. 0) .and. (istep > 0)) then
                        ! print timing info
                        call centered_heading('timing', '-')
                        time1 = rtime()
                        time_per_step = (time1-time0)/nbcycle
                        time_completion = int(time_per_step*(nsteps-istep)/60)
                        time0 = time1
                        write(*,222) time_per_step, time_completion
222                             format('seconds per step (wall-clock): ',f5.2,&
                                ' estimated completion in',i6,' minutes')
                end if

                ! update lists of nonbonded interaction pairs
                if (nodeid .eq. 0) then
                   call centered_heading('nonbonded pair list generation', '-')
                end if
                call make_pair_lists
#if defined(dump)
write(*,332) 'solute-solute', 'solute-water', 'water-water', 'q-solute', 'q-water'
write(*,333) nodeid, 'count', nbpp_pair, nbpw_pair, &
     &  nbww_true_pair, nbqp_pair, 3*nqat*nbqw_pair

#if defined(use_mpi)
!reduce totnxx, i.e. collect # pairs found by slave nodes
nbxx(1)=nbpp_pair
nbxx(2)=nbpw_pair
nbxx(3)=nbww_true_pair
nbxx(4)=nbqp_pair
nbxx(5)=3*nqat*nbqw_pair

call mpi_reduce(nbxx,nbxx_tot,5,mpi_integer,mpi_sum,0,mpi_comm_world,ierr) 
if (ierr .ne. 0) call die('run/reduce')

if (nodeid .eq. 0) then
totnbpp=nbxx_tot(1)
totnbpw=nbxx_tot(2)
totnbww=nbxx_tot(3)
totnbqp=nbxx_tot(4)
totnbqw=nbxx_tot(5)
write(*,99) 'total', totnbpp,totnbpw,totnbww,totnbqp,totnbqw
end if
99  format(a10,1x,5(1x,i12))
#endif
332     format('node value ',5a13)
333     format(i4,1x,a5,1x,5(1x,i12))
#endif

        end if ! every nbcycle steps

        ! --- start of time step ---

        ! get potential energy and derivatives from ff
        call pot_energy

if(nodeid .eq. 0) then
        if ( mod(istep,iout_cycle) == 0 .and. monitor_group_pairs > 0) then
           call nonbond_monitor  
        end if

        ! off-diagonals
        if ( noffd .gt. 0 ) call offdiag

#ifndef dum
        ! update velocities from accelerations,
        ! scale velocities & update positions from velocities
#if defined (profiling)
start_loop_time1 = rtime()
#endif

                !solute atoms first
        do i=1,nat_solute
                i3=i*3-3
                v(i3+1)= ( v(i3+1)-d(i3+1)*winv(i)*dt ) * tscale_solute
                xx(i3+1) = x(i3+1)
                x(i3+1) = x(i3+1) + v(i3+1)*dt

                v(i3+2)= ( v(i3+2)-d(i3+2)*winv(i)*dt ) * tscale_solute
                xx(i3+2) = x(i3+2)
                x(i3+2) = x(i3+2) + v(i3+2)*dt

                v(i3+3)= ( v(i3+3)-d(i3+3)*winv(i)*dt ) * tscale_solute
                xx(i3+3) = x(i3+3)
                x(i3+3) = x(i3+3) + v(i3+3)*dt
        end do

                !now solvent atoms
        do i=nat_solute+1,natom
                i3=i*3-3
                v(i3+1)= ( v(i3+1)-d(i3+1)*winv(i)*dt ) * tscale_solvent
                xx(i3+1) = x(i3+1)
                x(i3+1) = x(i3+1) + v(i3+1)*dt

                v(i3+2)= ( v(i3+2)-d(i3+2)*winv(i)*dt ) * tscale_solvent
                xx(i3+2) = x(i3+2)
                x(i3+2) = x(i3+2) + v(i3+2)*dt

                v(i3+3)= ( v(i3+3)-d(i3+3)*winv(i)*dt ) * tscale_solvent
                xx(i3+3) = x(i3+3)
                x(i3+3) = x(i3+3) + v(i3+3)*dt
        end do
#if defined (profiling)
profile(11)%time = profile(11)%time + rtime() - start_loop_time1
#endif

        ! shake if necessary
        if(shake_constraints > 0) then
                niter=shake(xx, x)
                v(:) = (x(:) - xx(:)) / dt
        end if

        ! --- end of time step ---
#if defined (profiling)
start_loop_time2 = rtime()
#endif

        ! calculate temperature and scaling factor
        call temperature(temp,tscale_solute,tscale_solvent,ekinmax)
#if defined (profiling)
profile(12)%time = profile(12)%time + rtime() - start_loop_time2
#endif

end if !if(nodeid .eq. 0)


!change volume
if( use_pbc .and. constant_pressure) then
   if( mod(istep, ivolume_cycle)==0 .and. istep>0 ) then
      call mc_volume
   end if
end if

#if defined(use_mpi)
call mpi_bcast(x, nat3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast x')
#endif

        ! print [intermediate] results (master node only)
        if (nodeid .eq. 0) then
                ! trajectory, energy data, output and backup restart file
                if ( mod(istep,itrj_cycle) == 0 .and. istep > 0) then
                ! write_trj: write x to the trajectory file
                        call write_trj
                end if
                ! energies
                if ( mod(istep, iene_cycle) == 0 .and. istep > 0) then
                ! nrgy_put_ene(unit, e2, offd): print 'e2'=eq and offd to unit 'unit'=11
                        call put_ene(11, eq, offd)
                end if
                ! end-of-line, then call write_out, which will print a report on e and eq
                if ( mod(istep,iout_cycle) == 0 ) then
                        call write_out
                end if
                ! backup file of coordinates and velocities
                if ( mod(istep,1000) .eq. 0 ) then
                        call write_xfin
                end if
                if ( abs(temp-tlast)/temp > temp_print_threshold .or. &
                        (mod(istep, itemp_cycle) == 0 .and. istep > 0)) then
                        ! temperatures
                        tlast = temp
                        write(*,201) istep, temp, tfree
                                                if (detail_temps) then
                                                        write(*,2020) tfree_solute, tfree_solvent
!                                                       write(*,2030) texcl_solute, texcl_solvent

                                                end if

                end if


        end if ! print results


end do ! time step
201      format('temperature at step',i8,':         t_tot=',f10.1,'         t_free=',f10.1)
2020 format('                             t_free_solute=',f10.1,' t_free_solvent=',f10.1)
2030 format('                             t_excl_solute=',f10.1,' t_excl_solvent=',f10.1)


!***********************************************************************
!       end main dynamics loop
!***********************************************************************

! end of qdum exclusion
#else
end if !(nodeid .eq. 0) from far above
#endif

        ! write final trajectory image when istep = nsteps
#ifndef dum
if (nodeid .eq. 0) then
    if ( mod(istep,itrj_cycle) == 0) call write_trj
end if
#endif

        ! write output for final step and final coords
call make_pair_lists
call pot_energy
if (nodeid .eq. 0) then
        write(*,*)
        call write_out
        call write_xfin
end if



if (nodeid .eq. 0) then
        time1 = rtime()
        write (*,202) time1 - startloop
        202 format('total time of main loop:                     ', f15.1,'(s)')
end if
#if defined(profiling)
!print more profiling info

#if defined(use_mpi)
do i=1,num_profiling_times
        node_times(i) = profile(i)%time
end do
call mpi_gather(node_times,num_profiling_times,mpi_real8,all_node_times,num_profiling_times,mpi_real8,0,mpi_comm_world,ierr)
if (ierr .ne. 0) call die('md_run/mpi_gather profiling times')

if (nodeid .eq. 0) then
write (*,210,advance='no')
do j=0,numnodes-1
        write (*,209,advance='no') j
end do
write(*,*)

do i=1,num_profiling_times
        write (*,207,advance='no') profile(i)%name
        do j=0,numnodes-1
                write (*,208,advance='no') all_node_times(i+j*num_profiling_times)
        end do
        write (*,*) ' (s)'
end do

207     format('total time of ',a25,t40,': ')
208     format(f10.1,' ')
209     format(i11)
210 format(t30,'node:     ')

write (*,*)

end if

#else

do i=1,num_profiling_times
        write (*,207) profile(i)%name,profile(i)%time
end do
207     format('total time of ',a25,t40,': ',f15.1,' (s).')
#endif
#endif


end subroutine md_run

!-----------------------------------------------------------------------
subroutine nbpp_count(npp, nppcgp)
! arguments
integer                                         :: npp
integer                                         :: nppcgp(:)

! local variables
integer                                         :: i,j,ig,jg,ia,ja,i3,j3,nl
real(8)                                         :: rcut2,r2
integer                                         :: lj_code

real(8)                                         :: dx, dy, dz
! this routine counts non-bonded solute-solute atom pairs 
! excluding any q-atoms.

! uses the global variables:
!  rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, 
!  iaclib, max_nbr_range, listex, nexlong, listexlong


npp = 0
rcut2 = rcpp*rcpp

igloop: do ig = 1, ncgp_solute
nppcgp(ig) = 0

ia = cgp(ig)%iswitch
 
! skip if excluded group
if ( .not. use_pbc .and. excl(ia) ) cycle igloop

i3 = 3*ia-3

jgloop: do jg = 1, ncgp_solute
  ja = cgp(jg)%iswitch
 
   ! skip if excluded group
  if ( .not. use_pbc .and. excl(ja) ) cycle jgloop

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

  j3 = 3*ja-3

  !******pwadded if-statement 2001-10-01

  if( .not. use_pbc ) then
  r2 = ( x(i3+1) -x(j3+1) )**2 &
        +( x(i3+2) -x(j3+2) )**2 &
        +( x(i3+3) -x(j3+3) )**2
  else
        dx = x(i3+1) -x(j3+1)
        dy = x(i3+2) -x(j3+2)
        dz = x(i3+3) -x(j3+3)
        dx = dx - boxlength(1)*nint(dx*inv_boxl(1))
        dy = dy - boxlength(2)*nint(dx*inv_boxl(2))
        dz = dz - boxlength(3)*nint(dx*inv_boxl(3))
        r2 = dx**2 + dy**2 + dz**2
  end if


  ! skip if outside cutoff
  if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
        i = cgpatom(ia)

        ! skip if q-atom
        if ( iqatom(i)/=0 ) cycle ialoop

jaloop: do ja = cgp(jg)%first, cgp(jg)%last
          j = cgpatom(ja)

          ! skip if q-atom
          if ( iqatom(j)/=0 ) cycle jaloop

          if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

          ! skip if all interactions zero
          lj_code = ljcod(iac(i),iac(j))
          if((crg(i) * crg(j) == 0.) &
                .and. &
                (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(j))%avdw(lj_code) == 0.) &
                .and. &
                (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(j))%bvdw(lj_code) == 0.)) &
                cycle jaloop

          ! check bonded exclusions and 1-4 nbors
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( listex(j-i,i) ) cycle jaloop
                else
                  if ( listex(i-j,j) ) cycle jaloop
                end if
          else
                do nl = 1, nexlong
                  if ( (listexlong(1,nl) .eq. i .and. &
                        listexlong(2,nl) .eq. j      ) .or. &
                        (listexlong(1,nl) .eq. j .and. &
                        listexlong(2,nl) .eq. i      ) ) cycle jaloop
                end do
          end if

          ! passed all tests -- count the pair
          npp = npp + 1
          nppcgp(ig) = nppcgp(ig) + 1

        end do jaloop
  end do ialoop
end do jgloop
end do igloop

end subroutine nbpp_count

!-----------------------------------------------------------------------


subroutine nbpplis2
! local variables
integer                                         :: i,j,ig,jg,ia,ja,i3,j3,nl,inside
real(8)                                         :: rcut2,r2


! for spherical boundary 
!       this routine makes a list of non-bonded solute-solute atom pairs 
!       excluding any q-atoms.

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpp_pair = 0
rcut2 = rcpp*rcpp


igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end

ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

jgloop:     do jg = 1, ncgp_solute

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

ja = cgp(jg)%iswitch
         if ( excl(ja) ) cycle jgloop

!             --- outside cutoff ? ---
        inside = 0
        ia = cgp(ig)%first
        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
   i = cgpatom(ia)
   i3 = 3*i-3

           ja = cgp(jg)%first
           do while ((ja .le. cgp(jg)%last) .and. (inside .eq. 0))
      j = cgpatom(ja)
      j3 = 3*j-3


                  r2 = ( x(i3+1) -x(j3+1) )**2 &
                          +( x(i3+2) -x(j3+2) )**2 &
                          +( x(i3+3) -x(j3+3) )**2

      if ( r2 .le. rcut2 ) then
                    ! one atom pair is within cutoff: set inside
                    inside = 1
                  end if

                  ja = ja + 1
   end do

           ia = ia + 1
end do
        if (inside .eq. 0) cycle jgloop

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last
   i = cgpatom(ia)

   !             --- q-atom ? ---
   if ( iqatom(i)/=0 ) cycle ialoop

jaloop:       do ja = cgp(jg)%first, cgp(jg)%last
      j = cgpatom(ja)

      !             --- q-atom ? ---
      if ( iqatom(j)/=0 ) cycle jaloop

                        ! count once
      if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

      !             --- check bonded exclusions and 1-4 nbors ---

      if ( abs(j-i) .le. max_nbr_range ) then
         if ( i .lt. j ) then
            if ( listex(j-i,i) ) cycle jaloop
         else
            if ( listex(i-j,j) ) cycle jaloop
         end if
      else
         do nl = 1, nexlong
            if ( (listexlong(1,nl) .eq. i .and. &
                 listexlong(2,nl) .eq. j      ) .or. &
                 (listexlong(1,nl) .eq. j .and. &
                 listexlong(2,nl) .eq. i      ) ) cycle jaloop
         end do
      end if

                  ! if out of space then make more space
      if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

      nbpp_pair = nbpp_pair + 1
      nbpp(nbpp_pair)%i = i
      nbpp(nbpp_pair)%j = j 
      nbpp(nbpp_pair)%ljcod = ljcod(iac(i),iac(j))

      if ( abs(j-i) .le. max_nbr_range ) then
         if ( i .lt. j ) then
            if ( list14(j-i,i) ) nbpp(nbpp_pair)%ljcod = 3
         else
            if ( list14(i-j,j) ) nbpp(nbpp_pair)%ljcod = 3
         end if
      else
         do nl = 1, n14long
            if ( (list14long(1,nl) .eq. i .and. &
                 list14long(2,nl) .eq. j      ) .or. &
                 (list14long(1,nl) .eq. j .and. &
                 list14long(2,nl) .eq. i      ) ) then
                 nbpp(nbpp_pair)%ljcod = 3
                                end if
         end do
      end if

                  end do jaloop
           end do ialoop
        end do jgloop
end do igloop
#if defined (profiling) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2


!--------------------------------------------------------------------

subroutine nbpplis2_box
  ! local variables
  integer                                               :: i,j,ig,jg,ia,ja,i3,j3,nl,inside
  real(8)                                               :: rcut2,r2
  real(8)                                               :: dx, dy, dz

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif
 
  ! for periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solute atom pairs 
  !     excluding any q-atoms.

  nbpp_pair = 0
  nbpp_cgp_pair = 0
  rcut2 = rcpp*rcpp

igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
jgloop:     do jg = 1, ncgp_solute

          ! count each charge group pair once only
          if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                   cycle jgloop
                
        !             --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
           i = cgpatom(ia)
           i3 = 3*i-3

                   ja = cgp(jg)%first
                   do while ((ja .le. cgp(jg)%last) .and. (inside .eq. 0))
              j = cgpatom(ja)
              j3 = 3*j-3
                          dx = x(i3+1) - x(j3+1)
                          dy = x(i3+2) - x(j3+2)
                          dz = x(i3+3) - x(j3+3)
                          dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                          dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                          dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                          r2 = dx**2 + dy**2 + dz**2

              if ( r2 .le. rcut2 ) then
                            ! one atom pair is within cutoff: set inside
                            inside = 1

                                if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp

                                nbpp_cgp_pair = nbpp_cgp_pair + 1
                                nbpp_cgp(nbpp_cgp_pair)%i = i
                                nbpp_cgp(nbpp_cgp_pair)%j = j

                          end if

                          ja = ja + 1
           end do

                   ia = ia + 1
        end do
                if (inside .eq. 0) cycle jgloop

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last
           i = cgpatom(ia)

           !             --- q-atom ? ---
           if ( iqatom(i)/=0 ) cycle ialoop

jaloop:       do ja = cgp(jg)%first, cgp(jg)%last
              j = cgpatom(ja)

              !             --- q-atom ? ---
              if ( iqatom(j)/=0 ) cycle jaloop

                                ! count once
              if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

              !             --- check bonded exclusions and 1-4 nbors ---

              if ( abs(j-i) .le. max_nbr_range ) then
                 if ( i .lt. j ) then
                    if ( listex(j-i,i) ) cycle jaloop
                 else
                    if ( listex(i-j,j) ) cycle jaloop
                 end if
              else
                 do nl = 1, nexlong
                    if ( (listexlong(1,nl) .eq. i .and. &
                         listexlong(2,nl) .eq. j      ) .or. &
                         (listexlong(1,nl) .eq. j .and. &
                         listexlong(2,nl) .eq. i      ) ) cycle jaloop
                 end do
              end if

                          ! if out of space then make more space
              if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

              nbpp_pair = nbpp_pair + 1
              nbpp(nbpp_pair)%i = i
              nbpp(nbpp_pair)%j = j 
              nbpp(nbpp_pair)%ljcod = ljcod(iac(i),iac(j))
                          nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair

              if ( abs(j-i) .le. max_nbr_range ) then
                 if ( i .lt. j ) then
                    if ( list14(j-i,i) ) nbpp(nbpp_pair)%ljcod = 3
                 else
                    if ( list14(i-j,j) ) nbpp(nbpp_pair)%ljcod = 3
                 end if
              else
                 do nl = 1, n14long
                    if ( (list14long(1,nl) .eq. i .and. &
                         list14long(2,nl) .eq. j      ) .or. &
                         (list14long(1,nl) .eq. j .and. &
                         list14long(2,nl) .eq. i      ) ) then
                         nbpp(nbpp_pair)%ljcod = 3
                                        end if
                 end do
              end if

                          end do jaloop
                   end do ialoop
                end do jgloop
        end do igloop
#if defined (profiling) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box
!-----------------------------------------------------------------------------
subroutine nbpplis2_box_lrf
  ! local variables
  integer                                               :: i,j,ig,jg,ia,ja,i3,j3,nl,inside
  real(8)                                               :: rcut2,r2
  real(8)                                               :: dx, dy, dz
 
  real(8)                                               ::rclrf2,field0, field1, field2
  real(8)                                               ::dr(3)
  real(8)                                               ::boxshiftx, boxshifty, boxshiftz       
  integer                                               ::inside_lrf, is3
  ! for periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solute atom pairs 
  !     excluding any q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  nbpp_pair = 0
  nbpp_cgp_pair = 0
  rcut2 = rcpp*rcpp
  rclrf2 = rclrf*rclrf

igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
jgloop:     do jg = 1, ncgp_solute

          ! count each charge group pair once only
          if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                   cycle jgloop
                
        !             --- outside cutoff ? ---
                inside = 0
                inside_lrf = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
           i = cgpatom(ia)
           i3 = 3*i-3

                   ja = cgp(jg)%first
                   do while ((ja .le. cgp(jg)%last) .and. (inside .eq. 0))
              j = cgpatom(ja)
              j3 = 3*j-3
                          dx = x(i3+1) - x(j3+1)
                          dy = x(i3+2) - x(j3+2)
                          dz = x(i3+3) - x(j3+3)
                          dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                          dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                          dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                          r2 = dx**2 + dy**2 + dz**2

              if ( r2 .le. rcut2 ) then
                            ! one atom pair is within cutoff: set inside
                            inside = 1

                                if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp

                                nbpp_cgp_pair = nbpp_cgp_pair + 1
                                nbpp_cgp(nbpp_cgp_pair)%i = i
                                nbpp_cgp(nbpp_cgp_pair)%j = j
              elseif (r2 <= rclrf2) then
                            inside_lrf = 1
              end if
                          ja = ja + 1
           end do
                   ia = ia + 1
        end do
                
                if (inside .eq. 1) then
ialoop:    do ia = cgp(ig)%first, cgp(ig)%last
           i = cgpatom(ia)

           !             --- q-atom ? ---
           if ( iqatom(i)/=0 ) cycle ialoop

jaloop:       do ja = cgp(jg)%first, cgp(jg)%last
              j = cgpatom(ja)

              !             --- q-atom ? ---
              if ( iqatom(j)/=0 ) cycle jaloop

                                ! count once
              if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

              !             --- check bonded exclusions and 1-4 nbors ---

              if ( abs(j-i) .le. max_nbr_range ) then
                 if ( i .lt. j ) then
                    if ( listex(j-i,i) ) cycle jaloop
                 else
                    if ( listex(i-j,j) ) cycle jaloop
                 end if
              else
                 do nl = 1, nexlong
                    if ( (listexlong(1,nl) .eq. i .and. &
                         listexlong(2,nl) .eq. j      ) .or. &
                         (listexlong(1,nl) .eq. j .and. &
                         listexlong(2,nl) .eq. i      ) ) cycle jaloop
                 end do
              end if

                          ! if out of space then make more space
              if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

              nbpp_pair = nbpp_pair + 1
              nbpp(nbpp_pair)%i = i
              nbpp(nbpp_pair)%j = j 
              nbpp(nbpp_pair)%ljcod = ljcod(iac(i),iac(j))
                          nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair

              if ( abs(j-i) .le. max_nbr_range ) then
                 if ( i .lt. j ) then
                    if ( list14(j-i,i) ) nbpp(nbpp_pair)%ljcod = 3
                 else
                    if ( list14(i-j,j) ) nbpp(nbpp_pair)%ljcod = 3
                 end if
              else
                 do nl = 1, n14long
                    if ( (list14long(1,nl) .eq. i .and. &
                         list14long(2,nl) .eq. j      ) .or. &
                         (list14long(1,nl) .eq. j .and. &
                         list14long(2,nl) .eq. i      ) ) then
                         nbpp(nbpp_pair)%ljcod = 3
                                        end if
                 end do
              end if

                          end do jaloop
                   end do ialoop
        elseif((inside_lrf ==1) .and. (inside == 0)) then
        ! outside pp-cutoff but inside lrf cut-off use lrf
                
                !ig : jg calculation
                boxshiftx = x(i3+1) - lrf(jg)%cgp_cent(1)
                boxshifty = x(i3+2) - lrf(jg)%cgp_cent(2)
                boxshiftz = x(i3+3) - lrf(jg)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

        do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle

          is3 = i*3-3

          dr(1) = x(is3+1) - lrf(jg)%cgp_cent(1) - boxshiftx
          dr(2) = x(is3+2) - lrf(jg)%cgp_cent(2) - boxshifty
          dr(3) = x(is3+3) - lrf(jg)%cgp_cent(3) - boxshiftz
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        end do

                !jg : ig calculations
                boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
                boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
                boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)
        
                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

            do ja = cgp(jg)%first, cgp(jg)%last

          ! skip if q-atom
          j = cgpatom(ja)
          if ( iqatom(j)/=0 ) cycle

          j3 = j*3-3

          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

      end if ! outside cutoff

        end do jgloop
  end do igloop
#if defined (profiling) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpplis2_lrf
! local variables
integer                                         :: i,j,ig,jg,ia,ja,i3,j3,nl,is 
logical                                         ::      inside
real(8)                                         :: rcut2,r2,field0,field1,field2
real(8)                                         :: dr(3)
real(8)                                         ::      rclrf2

!       this routine makes a list of non-bonded solute-solute atom pairs 
!       excluding any q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpp_pair = 0
rcut2 = rcpp*rcpp
rclrf2 = rclrf*rclrf


igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end

        ! skip if excluded group
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop

jgloop: do jg = 1, ncgp_solute

                ! count each charge group pair once only
                if( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                        ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                        cycle jgloop

                !             --- excluded group ? ---
                ja = cgp(jg)%iswitch
                if ( excl(ja) ) cycle jgloop

                !             --- outside cutoff ? ---
                inside = .false.
pairloop:       do ia=cgp(ig)%first, cgp(ig)%last
                        i = cgpatom(ia)
                        i3 = 3*i-3

                        do ja=cgp(jg)%first, cgp(jg)%last
                                j = cgpatom(ja)
                                j3 = 3*j-3

                                r2 = ( x(i3+1) -x(j3+1) )**2 &
                                   +( x(i3+2) -x(j3+2) )**2 &
                                   +( x(i3+3) -x(j3+3) )**2

                                if ( r2 <= rcut2 ) then
                                        ! one atom pair is within cutoff: set inside
                                        inside = .true.
                                        exit pairloop
                                end if
                        end do
            end do pairloop

                !             --- inside cutoff ? ---
                if (inside) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)

                                !                    --- q-atom ? ---
                                if ( iqatom(i)/=0 ) cycle ialoop
jaloop:                         do ja = cgp(jg)%first, cgp(jg)%last
                                        j = cgpatom(ja)
                                        !                      --- q-atom ? ---
                                        if ( iqatom(j)/=0 ) cycle jaloop
                                        ! count once
                                        if ( ig .eq. jg .and. i .ge. j ) cycle jaloop
                                        !                      --- check bonded exclusions and 1-4 nbors ---
                                        if ( abs(j-i) .le. max_nbr_range ) then
                                                if ( i .lt. j ) then
                                                        if ( listex(j-i,i) ) cycle jaloop
                                                else
                                                        if ( listex(i-j,j) ) cycle jaloop
                                                end if
                                        else
                                                do nl = 1, nexlong
                                                        if ( (listexlong(1,nl) .eq. i .and. &
                                                                listexlong(2,nl) .eq. j      ) .or. &
                                                                (listexlong(1,nl) .eq. j .and. &
                                                                listexlong(2,nl) .eq. i      ) ) cycle jaloop
                                                end do
                                        end if

                                        ! if out of space then make more space
                                        if (nbpp_pair == calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j 
                                        nbpp(nbpp_pair)%ljcod = ljcod(iac(i),iac(j))

                                        if ( abs(j-i) .le. max_nbr_range ) then
                                                if ( i .lt. j ) then
                                                        if ( list14(j-i,i) ) nbpp(nbpp_pair)%ljcod = 3
                                                else
                                                        if ( list14(i-j,j) ) nbpp(nbpp_pair)%ljcod = 3
                                                end if
                                        else
                                                do nl = 1, n14long
                                                        if ( (list14long(1,nl) .eq. i .and. &
                                                                list14long(2,nl) .eq. j      ) .or. &
                                                                (list14long(1,nl) .eq. j .and. &
                                                                list14long(2,nl) .eq. i      ) ) &
                                                                nbpp(nbpp_pair)%ljcod = 3
                                                end do
                                        end if




                        end do jaloop
                        end do ialoop
                elseif(r2 <= rclrf2) then
                        ! outside pp-cutoff but inside lrf cut-off: use lrf

ialoop2:                do ia = cgp(ig)%first, cgp(ig)%last

                                !                    --- q-atom ? ---
                                i = cgpatom(ia)
                                if ( iqatom(i)/=0 ) cycle ialoop2

                                i3 = i*3-3

                                dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
                                dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
                                dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
                                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                                field0=crg(i)/(r2*sqrt(r2))
                                lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
                                lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
                                lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
                                lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
                                field1=3.*field0/r2
                                lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
                                lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
                                lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
                                lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
                                lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
                                lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
                                lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
                                lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
                                lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
                                field2=-field1/r2
                                lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                                lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                                lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                                lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                                lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                                lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
                                lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                                lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
                                lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                                lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                                lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                                lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
                                lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                                lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                                lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                                lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
                                lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                                lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                                lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                                lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
                                lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                                lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
                                lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                                lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                                lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                                lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                                lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
                        end do ialoop2

jaloop2:                do ja = cgp(jg)%first, cgp(jg)%last

                                !                    --- q-atom ? ---
                                j = cgpatom(ja)
                                if ( iqatom(j)/=0 ) cycle jaloop2

                                j3 = j*3-3

                                dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
                                dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
                                dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
                                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                                field0=crg(j)/(r2*sqrt(r2))
                                lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
                                lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
                                lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
                                lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
                                field1=3.*field0/r2
                                lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
                                lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)


                        lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
                                lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
                                lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
                                lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
                                lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
                                lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
                                lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
                                field2=-field1/r2
                                lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                                lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                                lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                                lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                                lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                                lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
                                lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                                lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
                                lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                                lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                                lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                                lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
                                lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                                lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                                lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                                lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
                                lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                                lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                                lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                                lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
                                lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                                lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
                                lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                                lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                                lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                                lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                                lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))

                        end do jaloop2
                end if
        end do jgloop
end do igloop
#if defined (profiling) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_lrf

!-----------------------------------------------------------------------

subroutine nbpplist
! local variables
integer                                         :: i,j,ig,jg,ia,ja,i3,j3,nl
real(8)                                         :: rcut2,r2
integer                                         :: lj_code

! for use with spherical boundary   
!       this routine makes a list of non-bonded solute-solute atom pairs 
!       excluding any q-atoms.

! uses the global variables:
!  rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
!  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

! reset #pairs

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpp_pair = 0
rcut2 = rcpp*rcpp

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
! for every assigned charge group:

! skip if excluded group
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

i3 = 3*ia-3

jgloop: do jg = 1, ncgp_solute
  ! for every charge group:

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

  ! skip if excluded group
  ja = cgp(jg)%iswitch
  if ( excl(ja) ) cycle jgloop

  j3 = 3*ja-3
  r2 = ( x(i3+1) - x(j3+1) )**2 &
          +( x(i3+2) - x(j3+2) )**2 &
          +( x(i3+3) - x(j3+3) )**2

  ! skip if outside cutoff
  if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
        ! for every atom in the charge group (of the outermost loop):
        i = cgpatom(ia)

        ! skip if q-atom
        if ( iqatom(i)/=0 ) cycle ialoop

jaloop: do ja = cgp(jg)%first, cgp(jg)%last
          ! for every atom in the charge group (innermost loop)
          j = cgpatom(ja)

          ! make sure each pair is only counted once
          if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

          ! skip if q-atom
          if ( iqatom(j)/=0 ) cycle jaloop

          lj_code = ljcod(iac(i),iac(j))

          ! skip if all interactions zero
          if((crg(i) * crg(j) == 0.) &
                .and. &
                (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(j))%avdw(lj_code) == 0.) &
                .and. &
                (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(j))%bvdw(lj_code) == 0.)) &
                cycle jaloop

          ! skip bonded exclusions and 1-4 nbors
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( listex(j-i, i) ) cycle jaloop
                else
                  if ( listex(i-j, j) ) cycle jaloop
                end if
          else
                do nl = 1, nexlong
                  if ((listexlong(1, nl) .eq. i .and. listexlong(2, nl) .eq. j) .or. &
                        (listexlong(1, nl) .eq. j .and. listexlong(2, nl) .eq. i) ) &
                        cycle jaloop
                end do
          end if

          ! if out of space then make more space
          if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

          ! all tests passed, add the pair
          nbpp_pair = nbpp_pair + 1
          nbpp(nbpp_pair)%i = i
          nbpp(nbpp_pair)%j = j 
          nbpp(nbpp_pair)%ljcod = lj_code

          ! set ljcod of the pair to 3 if the atoms have bonded interactions
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( list14(j-i, i) ) nbpp(nbpp_pair)%ljcod = 3
                else
                  if ( list14(i-j, j) ) nbpp(nbpp_pair)%ljcod = 3
                end if
          else


        do nl = 1, n14long
                  if ((list14long(1, nl) .eq. i .and. list14long(2, nl) .eq. j) .or. &
                        (list14long(1, nl) .eq. j .and. list14long(2, nl) .eq. i)) &
                        nbpp(nbpp_pair)%ljcod = 3
                end do
          end if
        end do jaloop
  end do ialoop
end do jgloop
end do igloop
#if defined (profiling) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist

!-----------------------------------------------------------------------

subroutine nbpplist_box
  ! local variables
  integer                                               :: i,j,ig,jg,ia,ja,i3,j3,nl,ig_sw, jg_sw
  real(8)                                               :: rcut2,r2
  integer                                               :: lj_code
  real(8)                                               :: dx, dy, dz
  
  ! for use with periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solute atom pairs 
  !     excluding any q-atoms.

  ! uses the global variables:
  !  rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
  rcut2 = rcpp*rcpp

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
        ! for every assigned charge group:
        ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
        i3 = 3*ig_sw-3

jgloop: do jg = 1, ncgp_solute
          ! for every charge group:

          ! count each charge group pair once only
          if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                   cycle jgloop

          jg_sw = cgp(jg)%iswitch !switching atom in charge group jg
          j3 = 3*jg_sw-3
        
          dx = x(i3+1) - x(j3+1)
          dy = x(i3+2) - x(j3+2)
          dz = x(i3+3) - x(j3+3)
          dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
          dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
          dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
          r2 = dx**2 + dy**2 + dz**2
        
          ! skip if outside cutoff
          if ( r2 .gt. rcut2 ) cycle jgloop

          !inside cutoff

          !check if more memory is needed
          if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
          !add the charge group pair
          nbpp_cgp_pair = nbpp_cgp_pair + 1
          nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
          nbpp_cgp(nbpp_cgp_pair)%j = jg_sw

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
                ! for every atom in the charge group ig (of the outermost loop):
                i = cgpatom(ia)

                ! skip if q-atom
                if ( iqatom(i)/=0 ) cycle ialoop

jaloop: do ja = cgp(jg)%first, cgp(jg)%last
                  ! for every atom in the charge group jg (innermost loop)
                  j = cgpatom(ja)

                  ! make sure each pair is only counted once
                  if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

                  ! skip if q-atom
                  if ( iqatom(j)/=0 ) cycle jaloop

                  lj_code = ljcod(iac(i),iac(j))

                  ! skip if all interactions zero
                  if((crg(i) * crg(j) == 0.) &
                        .and. &
                        (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(j))%avdw(lj_code) == 0.) &
                        .and. &
                        (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(j))%bvdw(lj_code) == 0.)) &
                        cycle jaloop

                  ! skip bonded exclusions and 1-4 nbors
                  if ( abs(j-i) .le. max_nbr_range ) then
                        if ( i .lt. j ) then
                          if ( listex(j-i, i) ) cycle jaloop
                        else
                          if ( listex(i-j, j) ) cycle jaloop
                        end if
                  else
                        do nl = 1, nexlong
                          if ((listexlong(1, nl) .eq. i .and. listexlong(2, nl) .eq. j) .or. &
                                (listexlong(1, nl) .eq. j .and. listexlong(2, nl) .eq. i) ) &
                                cycle jaloop
                        end do
                  end if

                  ! if out of space then make more space
                  if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                  ! all tests passed, add the pair
                  nbpp_pair = nbpp_pair + 1
                  nbpp(nbpp_pair)%i = i
                  nbpp(nbpp_pair)%j = j 
                  nbpp(nbpp_pair)%ljcod = lj_code
                  nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
                
                  ! set ljcod of the pair to 3 if the atoms have bonded interactions
                  if ( abs(j-i) .le. max_nbr_range ) then
                        if ( i .lt. j ) then
                          if ( list14(j-i, i) ) nbpp(nbpp_pair)%ljcod = 3
                        else
                          if ( list14(i-j, j) ) nbpp(nbpp_pair)%ljcod = 3
                        end if
                  else

                do nl = 1, n14long
                          if ((list14long(1, nl) .eq. i .and. list14long(2, nl) .eq. j) .or. &
                                (list14long(1, nl) .eq. j .and. list14long(2, nl) .eq. i)) &
                                nbpp(nbpp_pair)%ljcod = 3
                        end do
                  end if
                end do jaloop
          end do ialoop
        end do jgloop
  end do igloop
#if defined (profiling) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist_box

!--------------------------------------------------------------------------------------

subroutine nbpplist_lrf
! local variables
integer                                         :: i,j,ig,jg,ia,ja,i3,j3,nl,is,is3
real(8)                                         :: rcut2,r2,field0,field1,field2
real(8)                                         :: dr(3)
integer                                         :: lj_code
real(8)                                         ::      rclrf2

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!       this routine makes a list of non-bonded solute-solute atom pairs 
!       excluding any q-atoms.

! uses the global variables:
!  nbpp_pair, rcpp, rclrf, cgp, excl, ncgp, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range,
!  listex, nexlong, listexlong, nbpp, alloc_status, list14, n14long, list14long, lrf


nbpp_pair = 0
rcut2 = rcpp*rcpp
rclrf2 = rclrf*rclrf
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
! for every assigned charge group:

! skip if excluded group
is = cgp(ig)%iswitch
if ( excl(is) ) cycle igloop

is3 = 3*is-3

jgloop: do jg = 1, ncgp_solute
  ! for every charge group:


  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop

  ! skip if excluded group
  ja = cgp(jg)%iswitch
  if ( excl(ja) ) cycle jgloop

  j3 = 3*ja-3
  r2 = ( x(is3+1) -x(j3+1) )**2 &
     +( x(is3+2) -x(j3+2) )**2 &
     +( x(is3+3) -x(j3+3) )**2

  ! inside cutoff?
  if ( r2 .le. rcut2 ) then
ialoop: do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle ialoop

jaloop:   do ja = cgp(jg)%first, cgp(jg)%last
                j = cgpatom(ja)

                ! make sure each pair is only counted once
                if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

                ! skip if q-atom
                if ( iqatom(j)/=0 ) cycle jaloop

                lj_code = ljcod(iac(i),iac(j))

                ! skip if all interactions zero
                if ((crg(i) * crg(j) == 0.) &
                  .and. &
                  (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(j))%avdw(lj_code) == 0.) &
                  .and. &
                  (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(j))%bvdw(lj_code) == 0.)) &
                  cycle jaloop

                ! check bonded exclusions and 1-4 nbors
                if ( abs(j-i) .le. max_nbr_range ) then
                  if ( i .lt. j ) then
                        if ( listex(j-i,i) ) cycle jaloop
                  else
                        if ( listex(i-j,j) ) cycle jaloop
                  end if
                else
                  do nl = 1, nexlong
                        if ( (listexlong(1,nl) .eq. i .and. &
                          listexlong(2,nl) .eq. j      ) .or. &
                          (listexlong(1,nl) .eq. j .and. &
                          listexlong(2,nl) .eq. i ) ) &
                          cycle jaloop
                  end do
                end if

                ! if out of space then make more space
                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                ! all tests passed, add the pair
                nbpp_pair = nbpp_pair + 1
                nbpp(nbpp_pair)%i = i
                nbpp(nbpp_pair)%j = j 
                nbpp(nbpp_pair)%ljcod = ljcod(iac(i),iac(j))
!tmp                nbpp_per_cgp(ig)=nbpp_per_cgp(ig)+1

                if ( abs(j-i) .le. max_nbr_range ) then
                  if ( i .lt. j ) then
                        if ( list14(j-i,i) ) nbpp(nbpp_pair)%ljcod = 3
                  else
                        if ( list14(i-j,j) ) nbpp(nbpp_pair)%ljcod = 3
                  end if
                else
                  do nl = 1, n14long
                        if ( (list14long(1,nl) .eq. i .and. &
                          list14long(2,nl) .eq. j      ) .or. &
                          (list14long(1,nl) .eq. j .and. &
                          list14long(2,nl) .eq. i      ) ) &
                          nbpp(nbpp_pair)%ljcod = 3
                  end do
                end if
          end do jaloop
        end do ialoop

  elseif(r2 <= rclrf2) then
        ! outside pp-cutoff but inside lrf cut-off use lrf

        do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle

          i3 = i*3-3

          dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
          dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
          dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        end do

        do ja = cgp(jg)%first, cgp(jg)%last

          ! skip if q-atom
          j = cgpatom(ja)
          if ( iqatom(j)/=0 ) cycle

          j3 = j*3-3

          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

  end if ! outside cutoff

end do jgloop
end do igloop

#if defined (profiling)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_lrf
!----------------lrf version of pw pbc-----------------------
subroutine nbpplist_box_lrf
  ! local variables
  integer                                               :: i,j,ig,jg,ia,ja,i3,j3,nl,ig_sw, jg_sw, is3
  real(8)                                               :: rcut2,r2
  integer                                               :: lj_code
  real(8)                                               :: dx, dy, dz
  
  real(8)                                               ::rclrf2,field0, field1, field2
  real(8)                                               ::dr(3)
  real(8)                                               ::boxshiftx, boxshifty, boxshiftz       

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! for use with periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solute atom pairs 
  !     excluding any q-atoms.

  ! uses the global variables:
  !  rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs
  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
  rcut2 = rcpp*rcpp
  rclrf2 = rclrf*rclrf

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
        ! for every assigned charge group:
        ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
        i3 = 3*ig_sw-3

jgloop: do jg = 1, ncgp_solute
          ! for every charge group:

          ! count each charge group pair once only
          if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                   ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                   cycle jgloop

          jg_sw = cgp(jg)%iswitch !switching atom in charge group jg
          j3 = 3*jg_sw-3
        
          dx = x(i3+1) - x(j3+1)
          dy = x(i3+2) - x(j3+2)
          dz = x(i3+3) - x(j3+3)
          dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
          dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
          dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
          r2 = dx**2 + dy**2 + dz**2
        
          ! skip if outside cutoff
          if ( r2 .le. rcut2 ) then 

        !innside cutoff

            !check if more memory is needed
            if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
            !add the charge group pair
            nbpp_cgp_pair = nbpp_cgp_pair + 1
            nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
            nbpp_cgp(nbpp_cgp_pair)%j = jg_sw

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
            ! for every atom in the charge group ig (of the outermost loop):
                i = cgpatom(ia)

                ! skip if q-atom
                if ( iqatom(i)/=0 ) cycle ialoop

jaloop:   do ja = cgp(jg)%first, cgp(jg)%last
                    ! for every atom in the charge group jg (innermost loop)
                    j = cgpatom(ja)

                    ! make sure each pair is only counted once
                    if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

                    ! skip if q-atom
                    if ( iqatom(j)/=0 ) cycle jaloop

                    lj_code = ljcod(iac(i),iac(j))

                    ! skip if all interactions zero
                    if((crg(i) * crg(j) == 0.) &
                          .and. &
                          (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(j))%avdw(lj_code) == 0.) &
                          .and. &
                          (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(j))%bvdw(lj_code) == 0.)) &
                          cycle jaloop

                    ! skip bonded exclusions and 1-4 nbors
                    if ( abs(j-i) .le. max_nbr_range ) then
                          if ( i .lt. j ) then
                            if ( listex(j-i, i) ) cycle jaloop
                          else
                            if ( listex(i-j, j) ) cycle jaloop
                          end if
                    else
                          do nl = 1, nexlong
                            if ((listexlong(1, nl) .eq. i .and. listexlong(2, nl) .eq. j) .or. &
                                  (listexlong(1, nl) .eq. j .and. listexlong(2, nl) .eq. i) ) &
                                  cycle jaloop
                          end do
                    end if

                    ! if out of space then make more space
                    if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                    ! all tests passed, add the pair
                    nbpp_pair = nbpp_pair + 1
                    nbpp(nbpp_pair)%i = i
                    nbpp(nbpp_pair)%j = j 
                    nbpp(nbpp_pair)%ljcod = lj_code
                    nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
                
                    ! set ljcod of the pair to 3 if the atoms have bonded interactions
                    if ( abs(j-i) .le. max_nbr_range ) then
                          if ( i .lt. j ) then
                            if ( list14(j-i, i) ) nbpp(nbpp_pair)%ljcod = 3
                          else
                            if ( list14(i-j, j) ) nbpp(nbpp_pair)%ljcod = 3
                          end if
                    else

                      do nl = 1, n14long
                       if ((list14long(1, nl) .eq. i .and. list14long(2, nl) .eq. j) .or. &
                              (list14long(1, nl) .eq. j .and. list14long(2, nl) .eq. i)) &
                                  nbpp(nbpp_pair)%ljcod = 3
                          end do
                    end if
        
                  end do jaloop
            end do ialoop
      
          elseif(r2 <= rclrf2) then
        ! outside pp-cutoff but inside lrf cut-off use lrf
                
                !ig : jg calculation
                boxshiftx = x(i3+1) - lrf(jg)%cgp_cent(1)
                boxshifty = x(i3+2) - lrf(jg)%cgp_cent(2)
                boxshiftz = x(i3+3) - lrf(jg)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

        do ia = cgp(ig)%first, cgp(ig)%last

          ! skip if q-atom
          i = cgpatom(ia)
          if ( iqatom(i)/=0 ) cycle

          is3 = i*3-3

          dr(1) = x(is3+1) - lrf(jg)%cgp_cent(1) - boxshiftx
          dr(2) = x(is3+2) - lrf(jg)%cgp_cent(2) - boxshifty
          dr(3) = x(is3+3) - lrf(jg)%cgp_cent(3) - boxshiftz
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        end do

                !jg : ig calculations
                boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
                boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
                boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)
        
                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

            do ja = cgp(jg)%first, cgp(jg)%last

          ! skip if q-atom
          j = cgpatom(ja)
          if ( iqatom(j)/=0 ) cycle

          j3 = j*3-3

          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

      end if ! outside cutoff

        end do jgloop
  end do igloop
#if defined (profiling)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_box_lrf
!-----------------------------------------------------------------------
!******pwchanged 2001-10-01
subroutine nbpw_count(npw, npwcgp)
! arguments
integer                                         :: npw
integer                                         :: npwcgp(:)

! local variables
integer                                         :: i,ig,jg,ia,ja,i3,j3
real(8)                                         :: rcut2,r2
integer                                         :: lj_code

!******pwadded variables 2001-10-01

real(8)                                         :: dx, dy, dz

!       this routine makes a list of non-bonded solute-solvent atom pairs
!       excluding q-atoms.

! uses the global variables:
!  rcpw, ncgp, cgp, excl, nwat, nat_solute, x, cgpatom, iqatom, ljcod, crg, iaclib

npw = 0
rcut2 = rcpw*rcpw

igloop: do ig = 1, ncgp_solute
! for each charge group of the protein:
npwcgp(ig) = 0

! skip if excluded charge group
ia = cgp(ig)%iswitch
if ( .not.use_pbc .and. excl(ia) ) cycle igloop

i3 = 3*ia-3

jgloop: do jg = 1, nwat
  ! for each water molecule:
ja = nat_solute + 3*jg-2
if(.not. use_pbc .and. excl(ja) ) cycle jgloop ! skip excluded waters

j3 = 3*ja-3

  if( .not. use_pbc ) then
        r2 = ( x(i3+1) -x(j3+1) )**2 &
                +( x(i3+2) -x(j3+2) )**2 &
                +( x(i3+3) -x(j3+3) )**2
  else
        dx = x(i3+1) - x(j3+1)
        dy = x(i3+2) - x(j3+2)
        dz = x(i3+3) - x(j3+3)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
  end if

! skip if outside cutoff
if ( r2 .gt. rcut2 ) cycle jgloop

ialoop:  do ia = cgp(ig)%first, cgp(ig)%last
! for each atom in the protein charge group:
i = cgpatom(ia)

! skip if q-atom
if ( iqatom(i)/=0 ) cycle ialoop

jaloop: do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
          ! for every atom of the water molecule:

          ! calculate lj_code for the pair
          lj_code = ljcod(iac(i),iac(ja))

          ! skip pairs with zero interaction
          if((crg(i) * crg(ja) == 0.) &
                  .and. &
                  (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(ja))%avdw(lj_code) == 0.) &
                  .and. &
                  (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(ja))%bvdw(lj_code) == 0.)) &
                  cycle jaloop

          ! count the pair
          npw = npw + 1
          npwcgp(ig) = npwcgp(ig) + 1

        end do jaloop
end do ialoop
end do jgloop
end do igloop

end subroutine nbpw_count

!-----------------------------------------------------------------------

subroutine nbpwlis2
! local variables
integer                                         :: i,ig,jg,ia,ja,i3,j3, inside
real(8)                                         :: rcut2,r2
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

!       this routine makes a list of non-bonded solute-solvent atom pairs
!       excluding q-atoms.


nbpw_pair = 0
rcut2 = rcpw*rcpw


igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end

!          --- excluded group ? ---
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

jgloop:     do jg = 1, nwat
ja = nat_solute + 3*jg-2
if ( excl(ja) ) cycle jgloop ! skip excluded waters
j3 = 3*ja-3

!             --- outside cutoff ? ---
        inside = 0
        ia = cgp(ig)%first
        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
   i = cgpatom(ia)
   i3 = 3*i-3

           r2 = ( x(i3+1) -x(j3+1) )**2 &
                  +( x(i3+2) -x(j3+2) )**2 &
                  +( x(i3+3) -x(j3+3) )**2

   if ( r2 .le. rcut2 ) then
             ! inside cutoff, raise the flag
                 inside = 1
           end if

          ia = ia + 1
end do
        if (inside .eq. 0) cycle jgloop

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last

   !             --- q-atom ? ---
   i = cgpatom(ia)
   if ( iqatom(i)/=0 ) cycle ialoop

           ! if out of space then make more space
           if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

   nbpw_pair = nbpw_pair + 3
   nbpw(nbpw_pair-2)%i = i
   nbpw(nbpw_pair-2)%j = ja 
   nbpw(nbpw_pair-2)%ljcod = ljcod(iac(i),iac(ja))
   nbpw(nbpw_pair-1)%i = i
   nbpw(nbpw_pair-1)%j = ja+1
   nbpw(nbpw_pair-1)%ljcod = ljcod(iac(i),iac(ja+1))
   nbpw(nbpw_pair  )%i = i
   nbpw(nbpw_pair  )%j = ja+2
   nbpw(nbpw_pair  )%ljcod = ljcod(iac(i),iac(ja+2))

end do ialoop
end do jgloop
end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2

!------------------------------------------------------------------------------

!******pwadded 2001-10-18
subroutine nbpwlis2_box
  ! local variables
  integer                                               :: i,ig,jg,ia,ja,i3,j3,ig_atom, inside
  real(8)                                               :: rcut2,r2
  real(8)                                               :: dx, dy,dz
  
  ! for periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solvent atom pairs
  !     excluding q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = rcpw*rcpw

igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end

jgloop:     do jg = 1, nwat

                ja = nat_solute + 3*jg-2
                j3 = 3*ja-3

        !             --- outside cutoff ? ---
                inside = 0
                ig_atom = cgp(ig)%first
                do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                i3 = 3*i-3

                dx = x(i3+1) - x(j3+1)
                dy = x(i3+2) - x(j3+2)
                dz = x(i3+3) - x(j3+3)
                dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                r2 = dx**2 + dy**2 + dz**2
                   
                 if ( r2 .le. rcut2 ) then
                     ! inside cutoff, raise the flag
                         inside = 1
                         if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                         nbpw_cgp_pair = nbpw_cgp_pair + 1
                         nbpw_cgp(nbpw_cgp_pair)%i = i
                         nbpw_cgp(nbpw_cgp_pair)%j = ja
                   end if

                  ig_atom = ig_atom + 1 !ia = ia + 1
        end do
                if (inside .eq. 0) cycle jgloop

ialoop:    do ia = cgp(ig)%first, cgp(ig)%last

           !             --- q-atom ? ---
           i = cgpatom(ia)
           if ( iqatom(i)/=0 ) cycle ialoop

                   ! if out of space then make more space
                   if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

           nbpw_pair = nbpw_pair + 3
           nbpw(nbpw_pair-2)%i = i
           nbpw(nbpw_pair-2)%j = ja 
           nbpw(nbpw_pair-2)%ljcod = ljcod(iac(i),iac(ja))
                   nbpw(nbpw_pair-2)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair-1)%i = i
           nbpw(nbpw_pair-1)%j = ja+1
           nbpw(nbpw_pair-1)%ljcod = ljcod(iac(i),iac(ja+1))
                   nbpw(nbpw_pair-1)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair  )%i = i
           nbpw(nbpw_pair  )%j = ja+2
           nbpw(nbpw_pair  )%ljcod = ljcod(iac(i),iac(ja+2))
                   nbpw(nbpw_pair  )%cgp_pair = nbpw_cgp_pair

        end do ialoop
     end do jgloop
  end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box
!-----------------------------------------------------------------------
subroutine nbpwlis2_box_lrf
  ! local variables
  integer                                               :: i,ig,jg,ia,ja,i3,j3,ig_atom, inside
  real(8)                                               :: rcut2,r2
  real(8)                                               :: dx, dy,dz
  !lrf
  real(8)                                               :: rclrf2, field0, field1, field2
  real(8)                                               :: dr(3)
  integer                                               :: jg_cgp, j, inside_lrf, is3
  real(8)                                               :: boxshiftx, boxshifty, boxshiftz
  
  ! for periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solvent atom pairs
  !     excluding q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = rcpw*rcpw
  rclrf2 = rclrf*rclrf

igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
jgloop:    do jg = 1, nwat
                        ja = nat_solute + 3*jg-2
                                j3 = 3*ja-3

                                !             --- outside cutoff ? ---
                                inside = 0
                                inside_lrf = 0
                                ig_atom = cgp(ig)%first
                         do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                                 i = cgpatom(ig_atom)
                                 i3 = 3*i-3

                                         dx = x(i3+1) - x(j3+1)
                                         dy = x(i3+2) - x(j3+2)
                                         dz = x(i3+3) - x(j3+3)
                                         dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                                         dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                                         dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                                 r2 = dx**2 + dy**2 + dz**2
                   
                                 if ( r2 .le. rcut2 ) then
                                           ! inside cutoff, raise the flag
                                           inside = 1
                                           if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                                           nbpw_cgp_pair = nbpw_cgp_pair + 1
                                           nbpw_cgp(nbpw_cgp_pair)%i = i
                                           nbpw_cgp(nbpw_cgp_pair)%j = ja
                                         elseif (r2 <= rclrf2) then
                             inside_lrf = 1
                                 end if

                             ig_atom = ig_atom + 1 !ia = ia + 1
                     end do

          if (inside .eq. 1) then

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

           !             --- q-atom ? ---
           i = cgpatom(ia)
           if ( iqatom(i)/=0 ) cycle ialoop

                   ! if out of space then make more space
                   if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

           nbpw_pair = nbpw_pair + 3
           nbpw(nbpw_pair-2)%i = i
           nbpw(nbpw_pair-2)%j = ja 
           nbpw(nbpw_pair-2)%ljcod = ljcod(iac(i),iac(ja))
                   nbpw(nbpw_pair-2)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair-1)%i = i
           nbpw(nbpw_pair-1)%j = ja+1
           nbpw(nbpw_pair-1)%ljcod = ljcod(iac(i),iac(ja+1))
                   nbpw(nbpw_pair-1)%cgp_pair = nbpw_cgp_pair
           nbpw(nbpw_pair  )%i = i
           nbpw(nbpw_pair  )%j = ja+2
           nbpw(nbpw_pair  )%ljcod = ljcod(iac(i),iac(ja+2))
                   nbpw(nbpw_pair  )%cgp_pair = nbpw_cgp_pair

        end do ialoop
      elseif((inside_lrf ==1) .and. (inside == 0)) then   
                ! outside pw-cutoff but inside lrf cut-off: use lrf
                !solut : solvent
                jg_cgp = iwhich_cgp(ja)
                
                boxshiftx = x(i3+1) - lrf(jg_cgp)%cgp_cent(1)
                boxshifty = x(i3+2) - lrf(jg_cgp)%cgp_cent(2)
                boxshiftz = x(i3+3) - lrf(jg_cgp)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

ialoop2: do ia = cgp(ig)%first, cgp(ig)%last
          i = cgpatom(ia)

          ! skip q-atoms
          if ( iqatom(i)/=0 ) cycle
                  is3 = i*3-3
                  
                  !jg = ncgp + iw
                  ! calculate dr and (d)r2
                  dr(1) = x(is3+1) - lrf(jg_cgp)%cgp_cent(1) - boxshiftx
                  dr(2) = x(is3+2) - lrf(jg_cgp)%cgp_cent(2) - boxshifty
                  dr(3) = x(is3+3) - lrf(jg_cgp)%cgp_cent(3) - boxshiftz
                                                
                  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf parameters for the charge group
                  field0=crg(i)/(r2*sqrt(r2))
                  lrf(jg_cgp)%phi0=lrf(jg_cgp)%phi0+field0*r2
                  lrf(jg_cgp)%phi1(1)=lrf(jg_cgp)%phi1(1)-field0*dr(1)
                  lrf(jg_cgp)%phi1(2)=lrf(jg_cgp)%phi1(2)-field0*dr(2)
          lrf(jg_cgp)%phi1(3)=lrf(jg_cgp)%phi1(3)-field0*dr(3)
                  field1=3.*field0/r2
                  lrf(jg_cgp)%phi2(1)=lrf(jg_cgp)%phi2(1)+field1*dr(1)*dr(1)-field0
                  lrf(jg_cgp)%phi2(2)=lrf(jg_cgp)%phi2(2)+field1*dr(1)*dr(2)
                  lrf(jg_cgp)%phi2(3)=lrf(jg_cgp)%phi2(3)+field1*dr(1)*dr(3)
                  lrf(jg_cgp)%phi2(4)=lrf(jg_cgp)%phi2(4)+field1*dr(2)*dr(1)
                  lrf(jg_cgp)%phi2(5)=lrf(jg_cgp)%phi2(5)+field1*dr(2)*dr(2)-field0
                  lrf(jg_cgp)%phi2(6)=lrf(jg_cgp)%phi2(6)+field1*dr(2)*dr(3)
                  lrf(jg_cgp)%phi2(7)=lrf(jg_cgp)%phi2(7)+field1*dr(3)*dr(1)
                  lrf(jg_cgp)%phi2(8)=lrf(jg_cgp)%phi2(8)+field1*dr(3)*dr(2)
                  lrf(jg_cgp)%phi2(9)=lrf(jg_cgp)%phi2(9)+field1*dr(3)*dr(3)-field0
                  field2=-field1/r2
                  lrf(jg_cgp)%phi3(1 )=lrf(jg_cgp)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                  lrf(jg_cgp)%phi3(2 )=lrf(jg_cgp)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                  lrf(jg_cgp)%phi3(3 )=lrf(jg_cgp)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                  lrf(jg_cgp)%phi3(4 )=lrf(jg_cgp)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                  lrf(jg_cgp)%phi3(5 )=lrf(jg_cgp)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                  lrf(jg_cgp)%phi3(6 )=lrf(jg_cgp)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
                  lrf(jg_cgp)%phi3(7 )=lrf(jg_cgp)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                  lrf(jg_cgp)%phi3(8 )=lrf(jg_cgp)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
                  lrf(jg_cgp)%phi3(9 )=lrf(jg_cgp)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                  lrf(jg_cgp)%phi3(10)=lrf(jg_cgp)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                  lrf(jg_cgp)%phi3(11)=lrf(jg_cgp)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                  lrf(jg_cgp)%phi3(12)=lrf(jg_cgp)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
                  lrf(jg_cgp)%phi3(13)=lrf(jg_cgp)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                  lrf(jg_cgp)%phi3(14)=lrf(jg_cgp)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                  lrf(jg_cgp)%phi3(15)=lrf(jg_cgp)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                  lrf(jg_cgp)%phi3(16)=lrf(jg_cgp)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
                  lrf(jg_cgp)%phi3(17)=lrf(jg_cgp)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                  lrf(jg_cgp)%phi3(18)=lrf(jg_cgp)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                  lrf(jg_cgp)%phi3(19)=lrf(jg_cgp)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                  lrf(jg_cgp)%phi3(20)=lrf(jg_cgp)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
                  lrf(jg_cgp)%phi3(21)=lrf(jg_cgp)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                  lrf(jg_cgp)%phi3(22)=lrf(jg_cgp)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
                  lrf(jg_cgp)%phi3(23)=lrf(jg_cgp)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                  lrf(jg_cgp)%phi3(24)=lrf(jg_cgp)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                  lrf(jg_cgp)%phi3(25)=lrf(jg_cgp)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                  lrf(jg_cgp)%phi3(26)=lrf(jg_cgp)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                  lrf(jg_cgp)%phi3(27)=lrf(jg_cgp)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
                 enddo ialoop2

                !solvent : solut        
                boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
                boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
                boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

jaloop2: do ja = 1, 3
                  j = nat_solute + 3*jg-3 + ja
                  j3 = j*3-3

                  ! calculate dr and (d)r2
                  dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
                  dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
                  dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
                  
              r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                  ! calculate lrf for the water molecule
                  field0=crg(j)/(r2*sqrt(r2))
                  lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
                  lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
                  lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
                  lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
                  field1=3.*field0/r2
                  lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
                  lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
                  lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
                  lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
                  lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
                  lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
                  lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
                  lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
                  lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
                  field2=-field1/r2
                  lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                  lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                  lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                  lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                  lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                  lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
                  lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                  lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
                  lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                  lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                  lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                  lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
                  lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                  lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                  lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                  lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
                  lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                  lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                  lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                  lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
                  lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                  lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
                  lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                  lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                  lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                  lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                  lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
         enddo jaloop2
      end if

     end do jgloop
  end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpwlis2_lrf
! local variables
integer                                         :: i,j,ig,iw,jg,ia,ja,i3,j3,inside,is,is3
real(8)                                         :: rcut2,r2,field0,field1,field2
real(8)                                         :: dr(3)
real(8)                                         ::      rclrf2
!       this routine makes a list of non-bonded solute-solvent atom pairs
!       excluding q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

nbpw_pair = 0
rcut2 = rcpw*rcpw
rclrf2 = rclrf*rclrf

igloop:  do ig = calculation_assignment%pw%start, calculation_assignment%pw%end

!          --- excluded group ? ---
is = cgp(ig)%iswitch
if ( excl(is) ) cycle igloop
is3 = 3*is-3

iwloop:     do iw = 1, nwat

ja = nat_solute + 3*iw-2
if(excl(ja)) cycle iwloop ! skip excluded waters
j3 = 3*ja-3

!             --- outside cutoff ? ---
        inside = 0
        ia = cgp(ig)%first
        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
   i = cgpatom(ia)
   i3 = 3*i-3

   r2 = ( x(i3+1) -x(j3+1) )**2 &
        +( x(i3+2) -x(j3+2) )**2 &
        +( x(i3+3) -x(j3+3) )**2

   if ( r2 .le. rcut2 ) then
             ! inside cutoff, raise the flag
                 inside = 1
           end if

          ia = ia + 1
end do

if ( inside .eq. 1 ) then

ialoop:       do ia = cgp(ig)%first, cgp(ig)%last

      !             --- q-atom ? ---
      i = cgpatom(ia)
      if ( iqatom(i)/=0 ) cycle ialoop

                  ! if out of space then make more space
                  if (nbpw_pair > calculation_assignment%pw%max-3) call reallocate_nonbondlist_pw

      nbpw_pair = nbpw_pair + 3
      nbpw(nbpw_pair-2)%i = i
      nbpw(nbpw_pair-2)%j = ja 
      nbpw(nbpw_pair-2)%ljcod = ljcod(iac(i),iac(ja))
      nbpw(nbpw_pair-1)%i = i
      nbpw(nbpw_pair-1)%j = ja+1
      nbpw(nbpw_pair-1)%ljcod = ljcod(iac(i),iac(ja+1))
      nbpw(nbpw_pair  )%i = i
      nbpw(nbpw_pair  )%j = ja+2
      nbpw(nbpw_pair  )%ljcod = ljcod(iac(i),iac(ja+2))

   end do ialoop

elseif(r2 <= rclrf2) then   
                ! outside pw-cutoff but inside lrf cut-off: use lrf

ialoop2:      do ia = cgp(ig)%first, cgp(ig)%last

      !              --- q-atom ? ---
      i = cgpatom(ia)
      if ( iqatom(i)/=0 ) cycle ialoop2

      i3 = i*3-3
      !jg = ncgp + iw
                  jg = iwhich_cgp(ja)

      dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
      dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
      dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
      r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

      field0=crg(i)/(r2*sqrt(r2))
      lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
      lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
      lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
      lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
      field1=3.*field0/r2
      lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
      lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
      lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
      lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
      lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
      lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
      lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
      lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
      lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
      field2=-field1/r2
      lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
      lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
      lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
      lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
      lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
      lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
      lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
      lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
      lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
      lrf(jg)%phi3(10)=lrf(jg)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
      lrf(jg)%phi3(11)=lrf(jg)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
      lrf(jg)%phi3(12)=lrf(jg)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
      lrf(jg)%phi3(13)=lrf(jg)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
      lrf(jg)%phi3(14)=lrf(jg)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
      lrf(jg)%phi3(15)=lrf(jg)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
      lrf(jg)%phi3(16)=lrf(jg)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
      lrf(jg)%phi3(17)=lrf(jg)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
      lrf(jg)%phi3(18)=lrf(jg)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
      lrf(jg)%phi3(19)=lrf(jg)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
      lrf(jg)%phi3(20)=lrf(jg)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
      lrf(jg)%phi3(21)=lrf(jg)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
      lrf(jg)%phi3(22)=lrf(jg)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
      lrf(jg)%phi3(23)=lrf(jg)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
      lrf(jg)%phi3(24)=lrf(jg)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
      lrf(jg)%phi3(25)=lrf(jg)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
      lrf(jg)%phi3(26)=lrf(jg)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
      lrf(jg)%phi3(27)=lrf(jg)%phi3(27)         +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))

   end do ialoop2

jaloop2:   do ja = 1, 3

      j = nat_solute + 3*iw-3 + ja

      j3 = j*3-3

      dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
      dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
      dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
      r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

      field0=crg(j)/(r2*sqrt(r2))
      lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
      lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
      lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
      lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
      field1=3.*field0/r2
      lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
      lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
      lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
      lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
      lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
      lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
      lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
      lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
      lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
      field2=-field1/r2
      lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 ) +field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
      lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 ) +field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
      lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 ) +field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
      lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 ) +field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
      lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 ) +field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
      lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 ) +field2*(5.*dr(1)*dr(2)*dr(3))
      lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 ) +field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
      lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 ) +field2*(5.*dr(1)*dr(3)*dr(2))
      lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 ) +field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
      lrf(ig)%phi3(10)=lrf(ig)%phi3(10) +field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
      lrf(ig)%phi3(11)=lrf(ig)%phi3(11) +field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
      lrf(ig)%phi3(12)=lrf(ig)%phi3(12) +field2*(5.*dr(2)*dr(1)*dr(3))
      lrf(ig)%phi3(13)=lrf(ig)%phi3(13) +field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
      lrf(ig)%phi3(14)=lrf(ig)%phi3(14) +field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
      lrf(ig)%phi3(15)=lrf(ig)%phi3(15) +field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
      lrf(ig)%phi3(16)=lrf(ig)%phi3(16) +field2*(5.*dr(2)*dr(3)*dr(1))
      lrf(ig)%phi3(17)=lrf(ig)%phi3(17) +field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
      lrf(ig)%phi3(18)=lrf(ig)%phi3(18) +field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
      lrf(ig)%phi3(19)=lrf(ig)%phi3(19) +field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
      lrf(ig)%phi3(20)=lrf(ig)%phi3(20) +field2*(5.*dr(3)*dr(1)*dr(2))
      lrf(ig)%phi3(21)=lrf(ig)%phi3(21) +field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
      lrf(ig)%phi3(22)=lrf(ig)%phi3(22) +field2*(5.*dr(3)*dr(2)*dr(1))
      lrf(ig)%phi3(23)=lrf(ig)%phi3(23) +field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
      lrf(ig)%phi3(24)=lrf(ig)%phi3(24) +field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
      lrf(ig)%phi3(25)=lrf(ig)%phi3(25) +field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
      lrf(ig)%phi3(26)=lrf(ig)%phi3(26) +field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
      lrf(ig)%phi3(27)=lrf(ig)%phi3(27) +field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))

   end do jaloop2

end if

end do iwloop

end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_lrf

!-----------------------------------------------------------------------

subroutine nbpwlist
! local variables
integer                                         :: i,ig,jg,ia,ja,i3,j3
real(8)                                         :: rcut2,r2
integer                                         :: lj_code


! for use with spherical boundary

!       this routine makes a list of non-bonded solute-solvent atom pairs
!       excluding q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

! reset nbpw_pair
nbpw_pair = 0
rcut2 = rcpw*rcpw

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
! for every charge group:

! skip excluded groups
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop
i3 = 3*ia-3

jgloop: do jg = 1, nwat
  ! for every water molecule:

ja = nat_solute + 3*jg-2
if(excl(ja)) cycle jgloop ! skip excluded waters

  j3 = 3*ja-3
  r2 = ( x(i3+1) -x(j3+1) )**2 &
  +( x(i3+2) -x(j3+2) )**2 &
  +( x(i3+3) -x(j3+3) )**2

  ! skip water outside cutoff
if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
    ! for every atom in the charge group:

        ! find the atom index of the atom in the charge group
i = cgpatom(ia)

!       skip q-atoms
        if ( iqatom(i)/=0 ) cycle ialoop

        ! if out of space then make more space
        if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw

jaloop: do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
          ! for every atom of the water molecule:

          ! calculate lj_code for the pair
          lj_code = ljcod(iac(i),iac(ja))

          ! skip pairs with zero interaction
          if((crg(i) * crg(ja) == 0.) &
                  .and. &
                  (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(ja))%avdw(lj_code) == 0.) &
                  .and. &
                  (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(ja))%bvdw(lj_code) == 0.)) &
                  cycle jaloop

          ! add the pair
          nbpw_pair = nbpw_pair + 1
          nbpw(nbpw_pair)%i = i
          nbpw(nbpw_pair)%j = ja 
          nbpw(nbpw_pair)%ljcod = lj_code
        end do jaloop
end do ialoop
end do jgloop
end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist

!-----------------------------------------------------------------------
!******pwadded 2001-10-18

subroutine nbpwlist_box
  ! local variables
  integer                                               :: i,ig,jg,ia,ja,i3,j3,ig_sw,jg_sw
  real(8)                                               :: rcut2,r2
  integer                                               :: lj_code
  real(8)                                               :: dx, dy, dz

  ! for use with periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solvent atom pairs
  !     excluding q-atoms.
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! reset nbpw_pair
  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = rcpw*rcpw

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
        ! for every charge group:
    ig_sw = cgp(ig)%iswitch
    i3 = 3*ig_sw-3

jgloop: do jg = 1, nwat
          ! for every water molecule:

      jg_sw = nat_solute + 3*jg-2

          j3 = 3*jg_sw-3
          dx = x(i3+1) - x(j3+1)
          dy = x(i3+2) - x(j3+2)
          dz = x(i3+3) - x(j3+3)
          dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
          dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
          dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
          r2 = dx**2 + dy**2 + dz**2

      ! skip water outside cutoff
      if ( r2 .gt. rcut2 ) cycle jgloop

          !inside cut-off
          !check if charge group pair list is big enough
          if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp

          nbpw_cgp_pair = nbpw_cgp_pair + 1
          nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
          nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
            ! for every atom in the charge group:

                ! find the atom index of the atom in the charge group
        i = cgpatom(ia)

       !        skip q-atoms
                if ( iqatom(i)/=0 ) cycle ialoop

                ! if out of space then make more space
                if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw
                   
jaloop: do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
                  ! for every atom of the water molecule:

                  ! calculate lj_code for the pair
                  lj_code = ljcod(iac(i),iac(ja))

                  ! skip pairs with zero interaction
                  if((crg(i) * crg(ja) == 0.) &
                          .and. &
                          (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(ja))%avdw(lj_code) == 0.) &
                          .and. &
                          (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(ja))%bvdw(lj_code) == 0.)) &
                          cycle jaloop

                  ! add the pair
                  nbpw_pair = nbpw_pair + 1
                  nbpw(nbpw_pair)%i = i
                  nbpw(nbpw_pair)%j = ja 
                  nbpw(nbpw_pair)%ljcod = lj_code
                  nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                  !write(*,'(a, i8, a, i7)') 'pw-pair ', nbpw_pair, ' belongs to charge group ' , nbpw(nbpw_pair)%cgp_pair
                end do jaloop
      end do ialoop
    end do jgloop
  end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_box
!-------------------------------------------------------------------------------------

subroutine nbpwlist_lrf
! local variables
integer                                         :: i,j,ig,iw,jg,ia,ja,i3,j3,is,is3
real(8)                                         :: rcut2,r2,field0,field1,field2
real(8)                                         :: dr(3)
integer                                         :: lj_code
real(8)                                         ::      rclrf2

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

!       this routine makes a list of non-bonded solute-solvent atom pairs
!       excluding q-atoms.

! reset nbpw_pair
nbpw_pair = 0
rcut2 = rcpw*rcpw
rclrf2 = rclrf*rclrf

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
! for every charge group:

! skip excluded groups
is = cgp(ig)%iswitch
if ( excl(is) ) cycle igloop
is3 = 3*is-3

iwloop: do iw = 1, nwat
  ! for every water molecule:

ja = nat_solute + 3*iw-2
if(excl(ja)) cycle iwloop ! skip excluded waters
j3 = 3*ja-3

r2 = ( x(is3+1) -x(j3+1) )**2 &
     +( x(is3+2) -x(j3+2) )**2 &
     +( x(is3+3) -x(j3+3) )**2

if ( r2 .le. rcut2 ) then
    ! within the cutoff radix:

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  ! skip q-atoms
  if ( iqatom(i)/=0 ) cycle ialoop

          ! if out of space then make more space
          if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw

jaloop:   do ja = nat_solute + 3*iw-2, nat_solute + 3*iw
                ! calculate the lj_code of the pair
                lj_code = ljcod(iac(i),iac(ja))

                ! skip pairs with zero interaction
                if((crg(i) * crg(ja) == 0.) &
                        .and. &
                        (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(ja))%avdw(lj_code) == 0.) &
                        .and. &
                        (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(ja))%bvdw(lj_code) == 0.)) &
                        cycle jaloop

                ! add the pair
                nbpw_pair = nbpw_pair + 1
                nbpw(nbpw_pair)%i = i
                nbpw(nbpw_pair)%j = ja 
                nbpw(nbpw_pair)%ljcod = lj_code
!tmp                nbpw_per_cgp(ig) = nbpw_per_cgp(ig) + 1
          end do jaloop
end do ialoop

elseif(r2 <= rclrf2) then   
                ! outside pw-cutoff but inside lrf cut-off: use lrf

ialoop2: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  ! skip q-atoms
  if ( iqatom(i)/=0 ) cycle

  i3 = i*3-3
  !jg = ncgp + iw
          jg = iwhich_cgp(ja)

          ! calculate dr and (d)r2
  dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
  dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
  dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf parameters for the charge group
  field0=crg(i)/(r2*sqrt(r2))
  lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
  lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
  lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
  lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
  field1=3.*field0/r2
  lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
  lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
  lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
  lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
  lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
  lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
  lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
  lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
  lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
  field2=-field1/r2
  lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
  lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
  lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
  lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
  lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
  lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
  lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
  lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
  lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
  lrf(jg)%phi3(10)=lrf(jg)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
  lrf(jg)%phi3(11)=lrf(jg)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
  lrf(jg)%phi3(12)=lrf(jg)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
  lrf(jg)%phi3(13)=lrf(jg)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
  lrf(jg)%phi3(14)=lrf(jg)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
  lrf(jg)%phi3(15)=lrf(jg)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
  lrf(jg)%phi3(16)=lrf(jg)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
  lrf(jg)%phi3(17)=lrf(jg)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
  lrf(jg)%phi3(18)=lrf(jg)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
  lrf(jg)%phi3(19)=lrf(jg)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
  lrf(jg)%phi3(20)=lrf(jg)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
  lrf(jg)%phi3(21)=lrf(jg)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
  lrf(jg)%phi3(22)=lrf(jg)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
  lrf(jg)%phi3(23)=lrf(jg)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
  lrf(jg)%phi3(24)=lrf(jg)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
  lrf(jg)%phi3(25)=lrf(jg)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
  lrf(jg)%phi3(26)=lrf(jg)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
  lrf(jg)%phi3(27)=lrf(jg)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
enddo ialoop2

jaloop2: do ja = 1, 3
  j = nat_solute + 3*iw-3 + ja
  j3 = j*3-3

          ! calculate dr and (d)r2
  dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
  dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
  dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
  
  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf for the water molecule
  field0=crg(j)/(r2*sqrt(r2))
  lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
  lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
  lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
  lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
  field1=3.*field0/r2
  lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
  lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
  lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
  lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
  lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
  lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
  lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
  lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
  lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
  field2=-field1/r2
  lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
  lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
  lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
  lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
  lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
  lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
  lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
  lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
  lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
  lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
  lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
  lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
  lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
  lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
  lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
  lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
  lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
  lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
  lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
  lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
  lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
  lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
  lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
  lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
  lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
  lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
  lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
enddo jaloop2

end if

end do iwloop 
end do igloop

#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_lrf
!---------------lrf version of pw pbc-----------------------
subroutine nbpwlist_box_lrf
  ! local variables
  integer                                               :: i,ig,jg,ia,ja,i3,j3,ig_sw,jg_sw
  real(8)                                               :: rcut2,r2
  integer                                               :: lj_code
  real(8)                                               :: dx, dy, dz
  ! lrf
  real(8)                                               :: rclrf2, field0, field1, field2
  real(8)                                               :: dr(3)
  integer                                               :: jg_cgp, j, is3
  real(8)                                               :: boxshiftx, boxshifty, boxshiftz

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! for use with periodic boundary conditions
  !     this routine makes a list of non-bonded solute-solvent atom pairs
  !     excluding q-atoms.

  ! reset nbpw_pair
  nbpw_pair = 0
  nbpw_cgp_pair = 0
  rcut2 = rcpw*rcpw
  rclrf2 = rclrf*rclrf

igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
        ! for every charge group:
    ig_sw = cgp(ig)%iswitch
    i3 = 3*ig_sw-3

jgloop: do jg = 1, nwat
          ! for every water molecule:

      jg_sw = nat_solute + 3*jg-2

          j3 = 3*jg_sw-3
          dx = x(i3+1) - x(j3+1)
          dy = x(i3+2) - x(j3+2)
          dz = x(i3+3) - x(j3+3)
          dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
          dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
          dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
          r2 = dx**2 + dy**2 + dz**2

      ! skip water outside cutoff
      if ( r2 .le. rcut2 ) then

          !inside cut-off
          !check if charge group pair list is big enough
          if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp

          nbpw_cgp_pair = nbpw_cgp_pair + 1
          nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
          nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
              ! for every atom in the charge group:

                  ! find the atom index of the atom in the charge group
          i = cgpatom(ia)

          !     skip q-atoms
                  if ( iqatom(i)/=0 ) cycle ialoop

                  ! if out of space then make more space
                  if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw
                   
jaloop:   do ja = nat_solute + 3*jg-2, nat_solute + 3*jg
                    ! for every atom of the water molecule:

                    ! calculate lj_code for the pair
                    lj_code = ljcod(iac(i),iac(ja))

                    ! skip pairs with zero interaction
                    if((crg(i) * crg(ja) == 0.) &
                            .and. &
                            (iaclib(iac(i))%avdw(lj_code)*iaclib(iac(ja))%avdw(lj_code) == 0.) &
                            .and. &
                            (iaclib(iac(i))%bvdw(lj_code)*iaclib(iac(ja))%bvdw(lj_code) == 0.)) &
                            cycle jaloop

                    ! add the pair
                    nbpw_pair = nbpw_pair + 1
                    nbpw(nbpw_pair)%i = i
                    nbpw(nbpw_pair)%j = ja 
                    nbpw(nbpw_pair)%ljcod = lj_code
                    nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                    !write(*,'(a, i8, a, i7)') 'pw-pair ', nbpw_pair, ' belongs to charge group ' , nbpw(nbpw_pair)%cgp_pair
                  end do jaloop
        end do ialoop

      elseif(r2 <= rclrf2) then   
                ! outside pw-cutoff but inside lrf cut-off: use lrf
                !solut : solvent
                jg_cgp = iwhich_cgp(jg_sw)
                
                boxshiftx = x(i3+1) - lrf(jg_cgp)%cgp_cent(1)
                boxshifty = x(i3+2) - lrf(jg_cgp)%cgp_cent(2)
                boxshiftz = x(i3+3) - lrf(jg_cgp)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

ialoop2: do ia = cgp(ig)%first, cgp(ig)%last
          i = cgpatom(ia)

          ! skip q-atoms
          if ( iqatom(i)/=0 ) cycle
                  is3 = i*3-3
                  
                  !jg = ncgp + iw
          
                  ! calculate dr and (d)r2
                  dr(1) = x(is3+1) - lrf(jg_cgp)%cgp_cent(1) - boxshiftx
                  dr(2) = x(is3+2) - lrf(jg_cgp)%cgp_cent(2) - boxshifty
                  dr(3) = x(is3+3) - lrf(jg_cgp)%cgp_cent(3) - boxshiftz
                                                
                  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          ! calculate lrf parameters for the charge group
                  field0=crg(i)/(r2*sqrt(r2))
                  lrf(jg_cgp)%phi0=lrf(jg_cgp)%phi0+field0*r2
                  lrf(jg_cgp)%phi1(1)=lrf(jg_cgp)%phi1(1)-field0*dr(1)
                  lrf(jg_cgp)%phi1(2)=lrf(jg_cgp)%phi1(2)-field0*dr(2)
          lrf(jg_cgp)%phi1(3)=lrf(jg_cgp)%phi1(3)-field0*dr(3)
                  field1=3.*field0/r2
                  lrf(jg_cgp)%phi2(1)=lrf(jg_cgp)%phi2(1)+field1*dr(1)*dr(1)-field0
                  lrf(jg_cgp)%phi2(2)=lrf(jg_cgp)%phi2(2)+field1*dr(1)*dr(2)
                  lrf(jg_cgp)%phi2(3)=lrf(jg_cgp)%phi2(3)+field1*dr(1)*dr(3)
                  lrf(jg_cgp)%phi2(4)=lrf(jg_cgp)%phi2(4)+field1*dr(2)*dr(1)
                  lrf(jg_cgp)%phi2(5)=lrf(jg_cgp)%phi2(5)+field1*dr(2)*dr(2)-field0
                  lrf(jg_cgp)%phi2(6)=lrf(jg_cgp)%phi2(6)+field1*dr(2)*dr(3)
                  lrf(jg_cgp)%phi2(7)=lrf(jg_cgp)%phi2(7)+field1*dr(3)*dr(1)
                  lrf(jg_cgp)%phi2(8)=lrf(jg_cgp)%phi2(8)+field1*dr(3)*dr(2)
                  lrf(jg_cgp)%phi2(9)=lrf(jg_cgp)%phi2(9)+field1*dr(3)*dr(3)-field0
                  field2=-field1/r2
                  lrf(jg_cgp)%phi3(1 )=lrf(jg_cgp)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                  lrf(jg_cgp)%phi3(2 )=lrf(jg_cgp)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                  lrf(jg_cgp)%phi3(3 )=lrf(jg_cgp)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                  lrf(jg_cgp)%phi3(4 )=lrf(jg_cgp)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                  lrf(jg_cgp)%phi3(5 )=lrf(jg_cgp)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                  lrf(jg_cgp)%phi3(6 )=lrf(jg_cgp)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
                  lrf(jg_cgp)%phi3(7 )=lrf(jg_cgp)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                  lrf(jg_cgp)%phi3(8 )=lrf(jg_cgp)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
                  lrf(jg_cgp)%phi3(9 )=lrf(jg_cgp)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                  lrf(jg_cgp)%phi3(10)=lrf(jg_cgp)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                  lrf(jg_cgp)%phi3(11)=lrf(jg_cgp)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                  lrf(jg_cgp)%phi3(12)=lrf(jg_cgp)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
                  lrf(jg_cgp)%phi3(13)=lrf(jg_cgp)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                  lrf(jg_cgp)%phi3(14)=lrf(jg_cgp)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                  lrf(jg_cgp)%phi3(15)=lrf(jg_cgp)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                  lrf(jg_cgp)%phi3(16)=lrf(jg_cgp)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
                  lrf(jg_cgp)%phi3(17)=lrf(jg_cgp)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                  lrf(jg_cgp)%phi3(18)=lrf(jg_cgp)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                  lrf(jg_cgp)%phi3(19)=lrf(jg_cgp)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                  lrf(jg_cgp)%phi3(20)=lrf(jg_cgp)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
                  lrf(jg_cgp)%phi3(21)=lrf(jg_cgp)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                  lrf(jg_cgp)%phi3(22)=lrf(jg_cgp)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
                  lrf(jg_cgp)%phi3(23)=lrf(jg_cgp)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                  lrf(jg_cgp)%phi3(24)=lrf(jg_cgp)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                  lrf(jg_cgp)%phi3(25)=lrf(jg_cgp)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                  lrf(jg_cgp)%phi3(26)=lrf(jg_cgp)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                  lrf(jg_cgp)%phi3(27)=lrf(jg_cgp)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
                 enddo ialoop2

                !solvent : solut        
                boxshiftx = x(j3+1) - lrf(ig)%cgp_cent(1)
                boxshifty = x(j3+2) - lrf(ig)%cgp_cent(2)
                boxshiftz = x(j3+3) - lrf(ig)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

jaloop2: do ja = 1, 3
                  j = nat_solute + 3*jg-3 + ja
                  j3 = j*3-3

                                  ! calculate dr and (d)r2
                  dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
                  dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
                  dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz
                  
              r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                                  ! calculate lrf for the water molecule
                  field0=crg(j)/(r2*sqrt(r2))
                  lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
                  lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
                  lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
                  lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
                  field1=3.*field0/r2
                  lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
                  lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
                  lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
                  lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
                  lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
                  lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
                  lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
                  lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
                  lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
                  field2=-field1/r2
                  lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
                  lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
                  lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
                  lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
                  lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
                  lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
                  lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
                  lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
                  lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
                  lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
                  lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
                  lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
                  lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
                  lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
                  lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
                  lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
                  lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
                  lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
                  lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
                  lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
                  lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
                  lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
                  lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
                  lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
                  lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
                  lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
                  lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
         enddo jaloop2
      end if
    end do jgloop
  end do igloop
#if defined (profiling)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif    

end subroutine nbpwlist_box_lrf
!-------------------------------------------------------------------------------------
subroutine make_qconn

integer                                         ::       i, iq, is, ia

allocate(qconn(nstates,nat_solute, nqat))
qconn(:,:,:) = 9

do iq = 1, nqat
        qconn(:, iqseq(iq), iq) = 1
end do

do iq = 1, nqat
        do is = 1, nstates
                i = iqseq(iq)
                call find_bonded(origin=i, current=i, level=1, state=is)
        end do
end do

!modify matrix to take special exclusions into account
do i = 1, nexspec
        iq = iqatom(exspec(i)%i)
        if(iq > 0) then
                do is = 1, nstates
                        if(exspec(i)%flag(is)) then
                                qconn(is, exspec(i)%j, iq) = 0 !exclude by setting to 0
                        end if
                end do
        end if
        iq = iqatom(exspec(i)%j)
        if(iq > 0) then
                do is = 1, nstates
                        if(exspec(i)%flag(is)) then
                                qconn(is, exspec(i)%i, iq) = 0 !exclude by setting to 0
                        end if
                end do
        end if
end do

end subroutine make_qconn


!------------------------------------------------------------------------------


recursive subroutine find_bonded(origin, current, level, state)
!args
integer, intent(in)                     ::      origin, current, level, state
!locals
integer                                         ::      b, newcurrent, newlevel

!find q-atom connectivity using the bond list and the q-bond list
!shaken bonds (code -1) must be taken into account, but not 
!redefined bonds in the topology
do b = 1, nbonds_solute
        if(bnd(b)%cod == 0) cycle !skip redefined (but not shaken)
        if(bnd(b)%i  == current) then
                newlevel = level + 1 
                newcurrent = bnd(b)%j
                if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                        qconn(state, newcurrent, iqatom(origin)) = newlevel
                        if(newlevel < 4) then
                                call find_bonded(origin, newcurrent, newlevel, state)
                        end if
                end if
        elseif(bnd(b)%j  == current) then
                newlevel = level + 1 
                newcurrent = bnd(b)%i
                if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                        qconn(state, newcurrent, iqatom(origin)) = newlevel
                        if(newlevel < 4) then
                                call find_bonded(origin, newcurrent, newlevel, state)
                        end if
                end if
        end if
end do
do b = 1, nqbond
        if(qbnd(b)%cod(state) > 0) then
                if(qbnd(b)%i  == current) then
                        newlevel = level + 1 
                        newcurrent = qbnd(b)%j
                        if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                                qconn(state, newcurrent, iqatom(origin)) = newlevel
                                if(newlevel < 4) then
                                        call find_bonded(origin, newcurrent, newlevel, state)
                                end if
                        end if
                elseif(qbnd(b)%j  == current) then
                        newlevel = level + 1 
                        newcurrent = qbnd(b)%i
                        if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                                qconn(state, newcurrent, iqatom(origin)) = newlevel
                                if(newlevel < 4) then
                                        call find_bonded(origin, newcurrent, newlevel, state)
                                end if
                        end if
                end if
        end if
end do

end subroutine find_bonded


!------------------------------------------------------------------------------


integer function nbqq_count()

integer                                         ::      iq, j, jq, is

nbqq_pair(:) = 0

!count q-q
do iq = 1, nqat - 1
        do jq = iq + 1, nqat
                do is = 1, nstates
                        if(qconn(is, iqseq(jq), iq) > 3) then
                                nbqq_pair(is) = nbqq_pair(is)+1
                        end if
                end do
        end do
end do
!count q-non-q
do j = 1, nat_solute
        if(iqatom(j) > 0) cycle
        if(any(qconn(:,j,:) <= 3)) then 
                !bonded or angled to at least one q-atom in any state
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, j, iq) >= 4) then
                                        nbqq_pair(is) = nbqq_pair(is)+1
                                end if
                        end do
                end do
        end if
end do

nbqq_count = maxval(nbqq_pair(:))
end function nbqq_count

!---------------------------------------------------------------------------------

subroutine nbqqlist
integer                                         ::      iq, j, jq, is, i, k,l
real(8)                     :: el_scale
logical                     :: set

nbqq_pair(:) = 0

!list q-q
do iq = 1, nqat - 1
   do jq = iq + 1, nqat
       j = iqseq(jq)
       do is = 1, nstates
           if(qconn(is, j, iq) > 3) then
                   nbqq_pair(is) = nbqq_pair(is)+1
                   nbqq(nbqq_pair(is),is)%iq = iq
                   nbqq(nbqq_pair(is),is)%j = j
                   nbqq(nbqq_pair(is),is)%jq = jq

                   if(qconn(is, j, iq) == 4) then
                         nbqq(nbqq_pair(is),is)%ljcod = 3
                   elseif(.not. qvdw_flag) then
                         nbqq(nbqq_pair(is),is)%ljcod = ljcod(iac(iqseq(iq)),iac(j))
                   else
                         nbqq(nbqq_pair(is),is)%ljcod = 1
                         do i = 1, nqexpnb
                           if ((iq == iqexpnb(i) .and.  &
                             jq == jqexpnb(i))  .or. &
                             ( jq == iqexpnb(i) .and. &
                             iq == jqexpnb(i)))  then
                                  nbqq(nbqq_pair(is),is)%ljcod = 2
                                  exit
                           end if
                         end do
                   end if !if (qconn = 4)
                                   if (nel_scale .eq. 0) then
                     nbqq(nbqq_pair(is),is)%el_scale = 1.0
                                   else
                                     set=.false.
                                     do i=1,nel_scale
                       k=qq_el_scale(i)%iqat
                       l=qq_el_scale(i)%jqat
                       if ((iq == k .and. jq == l) .or. &
                         (iq == l .and. jq == k)) then
                           nbqq(nbqq_pair(is),is)%el_scale = qq_el_scale(i)%el_scale
                                                   set=.true.
                                                   exit
                       end if
                       if (.not. set) nbqq(nbqq_pair(is),is)%el_scale = 1.0
                                     end do
                                   end if
           end if !if (qconn > 3)
       end do
   end do
end do

!list q-non-q
do j = 1, nat_solute
        if(iqatom(j) > 0) cycle
        if(any(qconn(:,j,:) <= 3)) then 
                !bonded or angled to at least one q-atom
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, j, iq) >= 4) then
                                        nbqq_pair(is) = nbqq_pair(is)+1
                                        nbqq(nbqq_pair(is),is)%iq = iq
                                        nbqq(nbqq_pair(is),is)%j = j
                                        nbqq(nbqq_pair(is),is)%jq = 0
                                                                                nbqq(nbqq_pair(is),is)%el_scale=1.0
                                        if(qconn(is, j, iq) == 4) then
                                                nbqq(nbqq_pair(is),is)%ljcod = 3
                                        elseif(qvdw_flag) then
                                                nbqq(nbqq_pair(is),is)%ljcod = 1
                                        else
                                                nbqq(nbqq_pair(is),is)%ljcod = ljcod(iac(iqseq(iq)),iac(j))
                                        end if
                                end if
                        end do
                end do
        end if
end do

end subroutine nbqqlist

!-----------------------------------------------------------------------
!******pwchanged 2001-10-01
subroutine nbqp_count(nqp, nqpcgp)
! arguments
integer                                         :: nqp
integer                                         :: nqpcgp(:)

! local variables
integer                                         :: ig,ia,i,j,iq,i3
real(8)                                         :: rcut2,r2

!******pwadded variables 2001-10-01

real(8)                                         :: dx, dy, dz

!       this routine counts non-bonded atom pairs involving
!       *one* q-atom and *one* non-q-atom, where the latter is *not connected*
!       (meaning not bonded or angled) to any q-atom.
!
!       ( i , j ) pairs correspond to
!       ( iq, j ) with first index being a q-atom with the *q-atom numbering*, &
!                 and the second index is the non-q-atom.


nqp = 0
rcut2 = rcq*rcq



if(nqat==0) return

! --- solute - q-atoms

igloop: do ig = 1, ncgp_solute

nqpcgp(ig) = 0

! skip if excluded group
ia = cgp(ig)%iswitch
if ( .not. use_pbc .and. excl(ia) ) cycle igloop
i3 = 3*ia-3

!******pwadded if 2001-10-01

if( .not. use_pbc ) then
r2 = ( x(i3+1) - xpcent(1) )**2 &
        +( x(i3+2) - xpcent(2) )**2 &
        +( x(i3+3) - xpcent(3) )**2
else
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
end if

! skip if outside cutoff
if ( r2 .gt. rcut2 ) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  ! check if already on qq list
  if(any(qconn(:,i,:) <= 3)) cycle ialoop

  ! count the pairs
  nqp = nqp + nqat
  nqpcgp(ig) = nqpcgp(ig) + nqat

end do ialoop
end do igloop

end subroutine nbqp_count

!-----------------------------------------------------------------------
!******pwchanged 2001-10-01
subroutine nbqw_count(nqw, nqwmol)
! arguments
integer                                         :: nqw
integer                                         :: nqwmol(:)

! local variables
integer                                         :: ig,ia,i,j,iq,i3
real(8)                                         :: rcut2,r2

!******pwadded variables

real(8)                                         :: dx, dy, dz

!       this routine counts water molecules that interact with q-atoms
nqw = 0
rcut2 = rcq*rcq


if(nqat==0) return

! --- solvent - q-atoms

iwloop: do ig = 1, nwat
nqwmol(ig) = 0
ia = nat_solute + 3*ig-2
if(.not. use_pbc .and. excl(ia)) cycle iwloop ! skip excluded waters
i3 = 3*ia-3

!******pwadded if-statement 2001-10-01
if( .not. use_pbc ) then
        r2 = ( x(i3+1) - xpcent(1) )**2 &
                +( x(i3+2) - xpcent(2) )**2 &
                +( x(i3+3) - xpcent(3) )**2
else
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
end if

! skip if outside cutoff
if ( r2 <= rcut2 ) then
        nqw = nqw + 1
        nqwmol(ig) = 3*nqat
end if
end do iwloop

end subroutine nbqw_count

!-----------------------------------------------------------------------

subroutine nbqplis2
! local variables
integer                                         :: ig,ia,i,j,iq,i3,nl,inside
real(8)                                         :: rcut2,r2
integer                                         :: xspec
logical, save                                   :: list_done


!       this routine makes a list of non-bonded atom pairs involving
!       *one* q-atom and *one* non-q-atom, where the latter is *not connected*
!       (meaning not bonded or angled) to any q-atom.
!
!       ( i , j ) pairs correspond to
!       ( iq, j ) with first index being a q-atom with the *q-atom numbering*, &
!                 and the second index is the non-q-atom.

! uses the global variables:
!  nbqs_pair, rcq, cgp, excl, cgpatom, x, xpcent, nqat, iqseq, 
!  qconn, calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

!don't remake pair list if q-atoms interact with all solute atoms (rcq>rexcl_o)
!and list already made
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif




if(list_done .and. rcq > rexcl_o) return

nbqp_pair = 0
rcut2 = rcq*rcq


if(nqat == 0) return


! --- solute - q-atoms

igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
! for every assigned charge group:

! skip if excluded group
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

! check cutoff
inside = 0
ia = cgp(ig)%first
do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
  i = cgpatom(ia)
  i3 = 3*i-3


  r2 = ( x(i3+1) - xpcent(1) )**2 &
        +( x(i3+2) - xpcent(2) )**2 &
        +( x(i3+3) - xpcent(3) )**2

  if ( r2 .le. rcut2 ) then
        inside = 1
  end if

  ia = ia + 1
end do
if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  ! check if already on qq list
  if(any(qconn(:,i,:) <= 3)) cycle ialoop

        ! if out of space then make more space
        if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

        !check special exclusions
        !is already done in make_qconn

        ! store the pair
        nbqp_pair = nbqp_pair + 1
        nbqp(nbqp_pair)%i = iq
        nbqp(nbqp_pair)%j = i
        nbqp(nbqp_pair)%ljcod = ljcod(iac(i),iac(iqseq(iq)))

        ! adjust ljcod for neighbors
        ! check only first state, should be same for all in this list
        if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%ljcod = 3

        !store special q-lj code = 1,1,3 for normal code = 1,2,3 resp.
        nbqp(nbqp_pair)%qljcod = nbqp(nbqp_pair)%ljcod
        if(nbqp(nbqp_pair)%qljcod == 2) nbqp(nbqp_pair)%qljcod = 1

  end do qaloop
end do ialoop
end do igloop

list_done = .true. !save this value
#if defined (profiling)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplis2



!----------------------------------------------------------------------------

subroutine nbqplis2_box
!  ! local variables
  integer                                               :: ig,ia,i,j,iq,i3,nl,inside,ig_atom
  real(8)                                               :: rcut2,r2
  integer                                               :: xspec
  real(8)                                               :: dx, dy, dz
 
  ! for periodic boundary conditions
  !     this routine makes a list of non-bonded atom pairs involving
  !     *one* q-atom and *one* non-q-atom, where the latter is *not connected*
  !     (meaning not bonded or angled) to any q-atom.
  !
  !     ( i , j ) pairs correspond to
  !     ( iq, j ) with first index being a q-atom with the *q-atom numbering*, &
  !               and the second index is the non-q-atom.

  ! uses the global variables:
  !  nbqs_pair, rcq, cgp, excl, cgpatom, x, xpcent, nqat, iqseq, 
  !  qconn, calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


  nbqp_pair = 0
  nbqp_cgp_pair = 0
  rcut2 = rcq*rcq
        
  if(nqat == 0) return

 ! --- solute - q-atoms

igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
    ! for every assigned charge group:

        ! check cutoff
        inside = 0
        ig_atom = cgp(ig)%first
        do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
          i = cgpatom(ig_atom)
          i3 = 3*i-3

        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2

          if ( r2 .le. rcut2 ) then
                inside = 1

                if(nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
                nbqp_cgp_pair = nbqp_cgp_pair + 1
                nbqp_cgp(nbqp_cgp_pair)%i = i !leave %j empty, equals qswitch
          end if

          ig_atom = ig_atom + 1 !ia = ia + 1
        end do
        if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
          i = cgpatom(ia)

          ! check if already on qq list
          if(any(qconn(:,i,:) <= 3)) cycle ialoop

                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

                !check special exclusions
                !is already done in make_qconn

                ! store the pair
                nbqp_pair = nbqp_pair + 1
                nbqp(nbqp_pair)%i = iq
                nbqp(nbqp_pair)%j = i
                nbqp(nbqp_pair)%ljcod = ljcod(iac(i),iac(iqseq(iq)))
                nbqp(nbqp_pair)%cgp_pair = nbqp_cgp_pair

                ! adjust ljcod for neighbors
                ! check only first state, should be same for all in this list
                if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%ljcod = 3
                
                !store special q-lj code = 1,1,3 for normal code = 1,2,3 resp.
                nbqp(nbqp_pair)%qljcod = nbqp(nbqp_pair)%ljcod
                if(nbqp(nbqp_pair)%qljcod == 2) nbqp(nbqp_pair)%qljcod = 1

          end do qaloop
        end do ialoop
  end do igloop
#if defined (profiling)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplis2_box
!-----------------------------------------------------------------------

subroutine nbqplist
! local variables
integer                                         :: ig,ia,i,j,iq,i3,nl
real(8)                                         :: rcut2,r2
integer                                         ::      xspec
logical, save                           ::      list_done = .false.

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!for spherical boundary 
!       this routine makes a list of non-bonded atom pairs involving
!       *one* q-atom and *one* non-q-atom, where the latter is *not connected*
!       (meaning not bonded or angled) to any q-atom.
!
!       ( i , j ) pairs correspond to
!       ( iq, j ) with first index being a q-atom with the *q-atom numbering*, &
!                 and the second index is the non-q-atom.

! global variables used:
!  nbqs_pair, rcq, cgp, excl, x, xpcent, cgpatom, nqat, iqseq,
!  calculation_assignment%qp%max, nbqs, ljcod, nwat, nat_solute

!don't remake pair list if q-atoms interact with all solute atoms (rcq>rexcl_o)
!and list already made
if(list_done .and. rcq > rexcl_o) return

nbqp_pair = 0
rcut2 = rcq*rcq

if(nqat==0) return

! --- solute - q-atoms

igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end

! skip if excluded group
ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop

i3 = 3*ia-3
r2 = ( x(i3+1) - xpcent(1) )**2 &
        +( x(i3+2) - xpcent(2) )**2 &
        +( x(i3+3) - xpcent(3) )**2

! skip if outside cutoff
if ( r2 .gt. rcut2 ) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

  i = cgpatom(ia)

  ! check if already on qq list
  if(any(qconn(:,i,:) <= 3)) cycle ialoop

        ! if out of space then make more space
        if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

        !check special exclusions
        !is already done in make_qconn

        ! add the pair
        nbqp_pair = nbqp_pair + 1
        nbqp(nbqp_pair)%i = iq
        nbqp(nbqp_pair)%j = i
        nbqp(nbqp_pair)%ljcod = ljcod(iac(i),iac(iqseq(iq)))

        ! adjust ljcod for neighbors
        ! check only first state, should be same for all in this list
        if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%ljcod = 3

        !store special q-lj code = 1,1,3 for normal code = 1,2,3 resp.
        nbqp(nbqp_pair)%qljcod = nbqp(nbqp_pair)%ljcod
        if(nbqp(nbqp_pair)%qljcod == 2) nbqp(nbqp_pair)%qljcod = 1
  end do qaloop
!tmp        nbqp_per_cgp(ig) = nbqp_per_cgp(ig) + nqat
end do ialoop
end do igloop

list_done = .true. !save this value

#if defined (profiling)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplist

!-----------------------------------------------------------------------

!******pwadded 2001-10-18
subroutine nbqplist_box
  ! local variables
  integer                                               :: ig,ia,i,j,iq,i3,nl,ig_sw
  real(8)                                               :: rcut2,r2
  integer                                               ::      xspec
  real(8)                                               :: dx, dy, dz
  integer                                               :: ga, gb

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  !for periodic boundary conditions     
  !     this routine makes a list of non-bonded atom pairs involving
  !     *one* q-atom and *one* non-q-atom, where the latter is *not connected*
  !     (meaning not bonded or angled) to any q-atom.
  !
  !     ( i , j ) pairs correspond to
  !     ( iq, j ) with first index being a q-atom with the *q-atom numbering*, &
  !               and the second index is the non-q-atom.

  ! global variables used:
  !  nbqs_pair, rcq, cgp, x, xpcent, cgpatom, nqat, iqseq,
  !  calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute


  nbqp_pair = 0
  nbqp_cgp_pair = 0
  rcut2 = rcq*rcq

  if(nqat==0) return

  ! --- solute - q-atoms
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end

        ig_sw = cgp(ig)%iswitch
        
        i3 = 3*ig_sw-3
        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2

        ! skip if outside cutoff
        if ( r2 .gt. rcut2 ) cycle igloop

        if( nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
        
        nbqp_cgp_pair = nbqp_cgp_pair + 1
        nbqp_cgp(nbqp_cgp_pair)%i = ig_sw
        

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

          i = cgpatom(ia)

          ! check if already on qq list
          if(any(qconn(:,i,:) <= 3)) cycle ialoop

          ! if out of space then make more space
          if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop: do iq = 1, nqat

                
                !check special exclusions
                !is already done in make_qconn

                ! add the pair
                nbqp_pair = nbqp_pair + 1
                nbqp(nbqp_pair)%i = iq
                nbqp(nbqp_pair)%j = i
                nbqp(nbqp_pair)%ljcod = ljcod(iac(i),iac(iqseq(iq)))
                nbqp(nbqp_pair)%cgp_pair = nbqp_cgp_pair
                
                ! adjust ljcod for neighbors
                ! check only first state, should be same for all in this list
                if(qconn(1, i, iq) == 4) nbqp(nbqp_pair)%ljcod = 3
                
                !store special q-lj code = 1,1,3 for normal code = 1,2,3 resp.
                nbqp(nbqp_pair)%qljcod = nbqp(nbqp_pair)%ljcod
                if(nbqp(nbqp_pair)%qljcod == 2) nbqp(nbqp_pair)%qljcod = 1

          end do qaloop
        end do ialoop
   end do igloop

#if defined (profiling)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplist_box
!------------------------------------------------------------------------------------

subroutine nbqwlist

! local variables
integer                                         :: ig,ia,i,i3
real(8)                                         :: rcut2,r2


!       this routine makes a list of water molecules within rcq from xpcent,
! i.e. the q-atom - water non-bond lists which implicitly includes all
! q-atoms with all atoms of the listed water
! waters may not have bonded interactions with q-atoms !
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif



!we don't have to remake the list if q-atoms interact with all waters
!(rcq > rexcl_o) and list already made (nbwq_pair >0)
if(rcq > rexcl_o .and. nbqw_pair > 0) return


nbqw_pair = 0
rcut2 = rcq*rcq


if(nqat==0) return


iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
ia = nat_solute + 3*ig-2
if( excl(ia) ) cycle iwloop ! skip excluded waters

i3 = 3*ia-3


r2 = ( x(i3+1) - xpcent(1) )**2 &
        +( x(i3+2) - xpcent(2) )**2 &
        +( x(i3+3) - xpcent(3) )**2


! store if inside cutoff
if ( r2 <= rcut2 ) then
        nbqw_pair = nbqw_pair + 1
        nbqw(nbqw_pair) = nat_solute + 3*ig-2
end if
end do iwloop
#if defined (profiling)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif

end subroutine nbqwlist

!-----------------------------------------------------------------------

!******pwadded 2001-10-18
subroutine nbqwlist_box

  ! local variables
  integer                                               :: ig,ia,i,i3
  real(8)                                               :: rcut2,r2
  real(8)                                               :: dx, dy, dz

  !     this routine makes a list of water molecules within rcq from xswitch,
  ! i.e. the q-atom - water non-bond lists which implicitly includes all
  ! q-atoms with all atoms of the listed water
  ! waters may not have bonded interactions with q-atoms !
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


        nbqw_pair = 0
        rcut2 = rcq*rcq

        if(nqat==0) return

iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
        ia = nat_solute + 3*ig-2
        i3 = 3*ia-3

        dx = x(i3+1) - x(3*qswitch-2)
        dy = x(i3+2) - x(3*qswitch-1)
        dz = x(i3+3) - x(3*qswitch)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2

        ! store if inside cutoff
        if ( r2 <= rcut2 ) then
                nbqw_pair = nbqw_pair + 1
                nbqw(nbqw_pair) = nat_solute + 3*ig-2
        end if
  end do iwloop
#if defined (profiling)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif

end subroutine nbqwlist_box
!-------------------------------------------------------------------------------------
!******pwchanged 2001-10-01
subroutine nbww_count(nww, nwwmol)
! arguments
integer                                         :: nww
integer                                         :: nwwmol(:)

! local variables
integer                                         :: iw,jw,ia,ja,i3,j3
real(8)                                         :: rcut2,r2

!******pwadded variables                

real(8)                                         :: dx, dy, dz

! this routine counts non-bonded solvent-solvent atom pairs.

nww = 0
rcut2 = rcww*rcww

iwloop: do iw = 1, nwat
nwwmol(iw) = 0

ia = nat_solute + 3*iw-2
if(.not. use_pbc .and. excl(ia)) cycle iwloop ! skip excluded waters

i3 = 3*ia-3

jwloop: do jw = 1, nwat
  ja = nat_solute + 3*jw-2
  if(.not. use_pbc .and. excl(ja)) cycle jwloop ! skip excluded waters

  j3 = 3*ja-3

  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw)) &
           cycle jwloop
  !******pwadded if 2001-10-01

  if( use_pbc ) then
        dx = x(i3+1) - x(j3+1)
        dy = x(i3+2) - x(j3+2)
        dz = x(i3+3) - x(j3+3)
        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
        r2 = dx**2 + dy**2 + dz**2
  else
        r2 = ( x(i3+1) - x(j3+1) )**2 &
           + ( x(i3+2) - x(j3+2) )**2 &
           + ( x(i3+3) - x(j3+3) )**2
  end if
  ! count the pair if inside cutoff
  if ( r2 .le. rcut2 ) then
        nww = nww + 1
        nwwmol(iw) = nwwmol(iw) + 9
  end if

end do jwloop
end do iwloop
end subroutine nbww_count

!-----------------------------------------------------------------------
subroutine nbwwlist
! local variables
integer                                         :: iw,jw,ia,ja,i3,j3
real(8)                                         :: rcut2,r2

! this routine makes a list of non-bonded solvent-solvent atom pairs

! uses the global variables:
!  nbww_pair, rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max,nbww_true_pair
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


nbww_pair = 0
nbww_true_pair=0
rcut2 = rcww*rcww

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end

ia = nat_solute + 3*iw-2

        if (.not.excl(ia)) then
        i3 = 3*ia-3

jwloop:         do jw=1, nwat

                  ja = nat_solute + 3*jw-2
                  if (excl(ja)) cycle jwloop ! skip excluded waters
                  j3 = 3*ja-3


                 ! count each w-w pair once only
                if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
                     ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. (iw .eq. jw))  cycle jwloop

                r2    = ( x(i3+1) -x(j3+1) )**2 &
                      + ( x(i3+2) -x(j3+2) )**2 &
                      + ( x(i3+3) -x(j3+3) )**2

                if ( r2 .le. rcut2 ) then
                                ! inside cutoff: add the pair
                                nbww_true_pair = nbww_true_pair + 9
                                nbww_pair = nbww_pair + 1
                                nbww(nbww_pair) = ja
                end if
                        ! if out of space then make more space  
                        if (nbww_pair .ge. calculation_assignment%ww%max) then
                                call reallocate_nonbondlist_ww
                        end if
                end do jwloop

        end if !if ia not excluded
        ! now mark the end of the list of molecules interacting with molecule iw
        ! by means of a zero element in the list
        nbww_pair = nbww_pair + 1
        nbww(nbww_pair) = 0     
end do iwloop
#if defined (profiling)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist
!------------pwadded 2001-10-18------------------------------------------
subroutine nbwwlist_box

! local variables
integer                                         :: iw,jw,ia,ja,i3,j3
real(8)                                         :: rcut2,r2
real(8)                                         :: dx, dy, dz

! for periodic boundary conditions
! this routine makes a list of non-bonded solvent-solvent atom pairs
! uses the global variables:
!  nbww_pair, rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max,nbww_true_pair
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


nbww_true_pair = 0
nbww_pair = 0
rcut2 = rcww*rcww

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end

        ia = nat_solute + 3*iw-2
        i3 = 3*ia-3

jwloop:         do jw=1, nwat

                        ja = nat_solute + 3*jw-2
                        j3 = 3*ja-3

                        ! count each w-w pair once only
                        if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
                           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
                           (iw .eq. jw)) &
                           cycle jwloop

                        dx = x(i3+1) - x(j3+1)
                        dy = x(i3+2) - x(j3+2)
                        dz = x(i3+3) - x(j3+3)
                        dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
                        dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
                        dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
                        r2 = dx**2 + dy**2 + dz**2

                        if ( r2 .le. rcut2 ) then
                                ! inside cutoff: add the pair
                                nbww_true_pair = nbww_true_pair + 9
                                nbww_pair = nbww_pair + 1
                                nbww(nbww_pair) = ja
                        end if

                        ! if out of space then make more space  
                        if (nbww_pair .ge. calculation_assignment%ww%max) then
                                call reallocate_nonbondlist_ww
                        end if

                end do jwloop

        ! now mark the end of the list of molecules interacting with molecule iw
        ! by means of a zero element in the list
        nbww_pair = nbww_pair + 1
        nbww(nbww_pair) = 0     
end do iwloop
#if defined (profiling)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_box

!---------------------------------------------------------------------------

subroutine nbwwlist_lrf
! local variables
integer                                         :: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3
real(8)                                         :: rcut2,r2,field0,field1,field2
real(8)                                         :: dr(3)
real(8)                                         ::      rclrf2

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

!       this routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max,nbww_true_pair






nbww_true_pair=0
nbww_pair = 0
rcut2 = rcww*rcww
rclrf2 = rclrf*rclrf

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
is  = nat_solute + 3*iw-2
if(.not. excl(is)) then
is3 = 3*is-3

jwloop: do jw = 1, nwat
  ja = nat_solute + 3*jw-2
  if(excl(ja)) cycle jwloop ! skip excluded waters
  j3 = 3*ja-3

  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw) ) &
           cycle jwloop

r2 = ( x(is3+1) -x(j3+1) )**2 &
        +( x(is3+2) -x(j3+2) )**2 &
        +( x(is3+3) -x(j3+3) )**2

  if ( r2 .le. rcut2 ) then
        ! inside cutoff: add the pair
        nbww_pair = nbww_pair + 1
        nbww_true_pair = nbww_true_pair + 9  !to get explicit no. of interactions
        nbww(nbww_pair) = ja 
  elseif(r2 <= rclrf2) then
        ! outside ww-cutoff but inside lrf cut-off: use lrf
        do ia=1,3
          i = nat_solute+iw*3-3+ia
          i3 = i*3-3
          !jg = ncgp + jw
          jg = iwhich_cgp(ja)

          dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1)
          dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2)
          dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)



          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo

        do ja = 1, 3
          j = nat_solute + 3*jw-3 + ja
          j3 = j*3-3
          !ig = ncgp + iw
          ig = iwhich_cgp(is)


          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1)
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2)
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3)
          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo
  end if

  ! if out of space then make more space
  if (nbww_pair .ge. calculation_assignment%ww%max) call reallocate_nonbondlist_ww

end do jwloop
end if !if ia not excluded
 !now mark the end of the list of molecules interacting with molecule iw
 !by means of a zero element in the list
 nbww_pair = nbww_pair + 1
 nbww(nbww_pair) = 0    

end do iwloop

#if defined (profiling)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_lrf
!--------------lrf version of pw pbc-----------------------------------
subroutine nbwwlist_box_lrf
! local variables
integer                                         :: i,j,ig,jg,iw,jw,ia,ja,i3,j3,is,is3
real(8)                                         :: rcut2,r2,field0,field1,field2
real(8)                                         :: dr(3)
real(8)                                         ::      rclrf2
real(8)                                         ::  dx, dy, dz
real(8)                     :: boxshiftx, boxshifty, boxshiftz  

#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!       this routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max,nbww_true_pair

nbww_true_pair=0
nbww_pair = 0
rcut2 = rcww*rcww
rclrf2 = rclrf*rclrf

iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
is  = nat_solute + 3*iw-2
is3 = 3*is-3

jwloop: do jw = 1, nwat
  ja = nat_solute + 3*jw-2
  j3 = 3*ja-3

  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw) ) &
           cycle jwloop

dx = x(is3+1) -x(j3+1)
dy = x(is3+2) -x(j3+2)
dz = x(is3+3) -x(j3+3)

dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
 
r2 = dx*dx + dy*dy + dz*dz

  if ( r2 .le. rcut2 ) then
        ! inside cutoff: add the pair
        nbww_pair = nbww_pair + 1
        nbww_true_pair = nbww_true_pair + 9  !to get explicit no. of interactions
        nbww(nbww_pair) = ja 
  elseif(r2 <= rclrf2) then
        ! outside ww-cutoff but inside lrf cut-off: use lrf
        
                !iw interaction     
            jg = iwhich_cgp(ja)
                
                boxshiftx = x(is3+1) -lrf(jg)%cgp_cent(1)
                boxshifty = x(is3+2) -lrf(jg)%cgp_cent(2)
                boxshiftz = x(is3+3) -lrf(jg)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

                do ia=1,3
          i = nat_solute+iw*3-3+ia
          i3 = i*3-3
          !jg = ncgp + jw
          

          dr(1) = x(i3+1) - lrf(jg)%cgp_cent(1) - boxshiftx
          dr(2) = x(i3+2) - lrf(jg)%cgp_cent(2) - boxshifty
          dr(3) = x(i3+3) - lrf(jg)%cgp_cent(3) - boxshiftz
          
                  r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(i)/(r2*sqrt(r2))
          lrf(jg)%phi0=lrf(jg)%phi0+field0*r2
          lrf(jg)%phi1(1)=lrf(jg)%phi1(1)-field0*dr(1)
          lrf(jg)%phi1(2)=lrf(jg)%phi1(2)-field0*dr(2)
          lrf(jg)%phi1(3)=lrf(jg)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(jg)%phi2(1)=lrf(jg)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(jg)%phi2(2)=lrf(jg)%phi2(2)+field1*dr(1)*dr(2)
          lrf(jg)%phi2(3)=lrf(jg)%phi2(3)+field1*dr(1)*dr(3)
          lrf(jg)%phi2(4)=lrf(jg)%phi2(4)+field1*dr(2)*dr(1)
          lrf(jg)%phi2(5)=lrf(jg)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(jg)%phi2(6)=lrf(jg)%phi2(6)+field1*dr(2)*dr(3)
          lrf(jg)%phi2(7)=lrf(jg)%phi2(7)+field1*dr(3)*dr(1)
          lrf(jg)%phi2(8)=lrf(jg)%phi2(8)+field1*dr(3)*dr(2)
          lrf(jg)%phi2(9)=lrf(jg)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(jg)%phi3(1 )=lrf(jg)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(jg)%phi3(2 )=lrf(jg)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(3 )=lrf(jg)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))


  lrf(jg)%phi3(4 )=lrf(jg)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(5 )=lrf(jg)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(6 )=lrf(jg)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(jg)%phi3(7 )=lrf(jg)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(8 )=lrf(jg)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(jg)%phi3(9 )=lrf(jg)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(10)=lrf(jg)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(jg)%phi3(11)=lrf(jg)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(jg)%phi3(12)=lrf(jg)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(jg)%phi3(13)=lrf(jg)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(14)=lrf(jg)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(jg)%phi3(15)=lrf(jg)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(jg)%phi3(16)=lrf(jg)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(jg)%phi3(17)=lrf(jg)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(18)=lrf(jg)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(19)=lrf(jg)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(jg)%phi3(20)=lrf(jg)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(jg)%phi3(21)=lrf(jg)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(jg)%phi3(22)=lrf(jg)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(jg)%phi3(23)=lrf(jg)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(jg)%phi3(24)=lrf(jg)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(jg)%phi3(25)=lrf(jg)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(jg)%phi3(26)=lrf(jg)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(jg)%phi3(27)=lrf(jg)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo
                
                !jw interaction
                ig = iwhich_cgp(is)
                
                boxshiftx = x(j3+1) -lrf(ig)%cgp_cent(1)
                boxshifty = x(j3+2) -lrf(ig)%cgp_cent(2)
                boxshiftz = x(j3+3) -lrf(ig)%cgp_cent(3)

                boxshiftx = boxlength(1)*nint( boxshiftx*inv_boxl(1) )
                boxshifty = boxlength(2)*nint( boxshifty*inv_boxl(2) )
                boxshiftz = boxlength(3)*nint( boxshiftz*inv_boxl(3) )

        do ja = 1, 3
          j = nat_solute + 3*jw-3 + ja
          j3 = j*3-3
          !ig = ncgp + iw
          
          dr(1) = x(j3+1) - lrf(ig)%cgp_cent(1) - boxshiftx
          dr(2) = x(j3+2) - lrf(ig)%cgp_cent(2) - boxshifty
          dr(3) = x(j3+3) - lrf(ig)%cgp_cent(3) - boxshiftz

          r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

          field0=crg(j)/(r2*sqrt(r2))
          lrf(ig)%phi0=lrf(ig)%phi0+field0*r2
          lrf(ig)%phi1(1)=lrf(ig)%phi1(1)-field0*dr(1)
          lrf(ig)%phi1(2)=lrf(ig)%phi1(2)-field0*dr(2)
          lrf(ig)%phi1(3)=lrf(ig)%phi1(3)-field0*dr(3)
          field1=3.*field0/r2
          lrf(ig)%phi2(1)=lrf(ig)%phi2(1)+field1*dr(1)*dr(1)-field0
          lrf(ig)%phi2(2)=lrf(ig)%phi2(2)+field1*dr(1)*dr(2)
          lrf(ig)%phi2(3)=lrf(ig)%phi2(3)+field1*dr(1)*dr(3)
          lrf(ig)%phi2(4)=lrf(ig)%phi2(4)+field1*dr(2)*dr(1)
          lrf(ig)%phi2(5)=lrf(ig)%phi2(5)+field1*dr(2)*dr(2)-field0
          lrf(ig)%phi2(6)=lrf(ig)%phi2(6)+field1*dr(2)*dr(3)
          lrf(ig)%phi2(7)=lrf(ig)%phi2(7)+field1*dr(3)*dr(1)
          lrf(ig)%phi2(8)=lrf(ig)%phi2(8)+field1*dr(3)*dr(2)
          lrf(ig)%phi2(9)=lrf(ig)%phi2(9)+field1*dr(3)*dr(3)-field0
          field2=-field1/r2
          lrf(ig)%phi3(1 )=lrf(ig)%phi3(1 )+field2*(5.*dr(1)*dr(1)*dr(1)-r2*3.*dr(1))
          lrf(ig)%phi3(2 )=lrf(ig)%phi3(2 )+field2*(5.*dr(1)*dr(1)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(3 )=lrf(ig)%phi3(3 )+field2*(5.*dr(1)*dr(1)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(4 )=lrf(ig)%phi3(4 )+field2*(5.*dr(1)*dr(2)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(5 )=lrf(ig)%phi3(5 )+field2*(5.*dr(1)*dr(2)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(6 )=lrf(ig)%phi3(6 )+field2*(5.*dr(1)*dr(2)*dr(3))
          lrf(ig)%phi3(7 )=lrf(ig)%phi3(7 )+field2*(5.*dr(1)*dr(3)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(8 )=lrf(ig)%phi3(8 )+field2*(5.*dr(1)*dr(3)*dr(2))
          lrf(ig)%phi3(9 )=lrf(ig)%phi3(9 )+field2*(5.*dr(1)*dr(3)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(10)=lrf(ig)%phi3(10)+field2*(5.*dr(2)*dr(1)*dr(1)-r2*dr(2))
          lrf(ig)%phi3(11)=lrf(ig)%phi3(11)+field2*(5.*dr(2)*dr(1)*dr(2)-r2*dr(1))
          lrf(ig)%phi3(12)=lrf(ig)%phi3(12)+field2*(5.*dr(2)*dr(1)*dr(3))
          lrf(ig)%phi3(13)=lrf(ig)%phi3(13)+field2*(5.*dr(2)*dr(2)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(14)=lrf(ig)%phi3(14)+field2*(5.*dr(2)*dr(2)*dr(2)-r2*3.*dr(2))
          lrf(ig)%phi3(15)=lrf(ig)%phi3(15)+field2*(5.*dr(2)*dr(2)*dr(3)-r2*dr(3))
          lrf(ig)%phi3(16)=lrf(ig)%phi3(16)+field2*(5.*dr(2)*dr(3)*dr(1))
          lrf(ig)%phi3(17)=lrf(ig)%phi3(17)+field2*(5.*dr(2)*dr(3)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(18)=lrf(ig)%phi3(18)+field2*(5.*dr(2)*dr(3)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(19)=lrf(ig)%phi3(19)+field2*(5.*dr(3)*dr(1)*dr(1)-r2*dr(3))
          lrf(ig)%phi3(20)=lrf(ig)%phi3(20)+field2*(5.*dr(3)*dr(1)*dr(2))
          lrf(ig)%phi3(21)=lrf(ig)%phi3(21)+field2*(5.*dr(3)*dr(1)*dr(3)-r2*dr(1))
          lrf(ig)%phi3(22)=lrf(ig)%phi3(22)+field2*(5.*dr(3)*dr(2)*dr(1))
          lrf(ig)%phi3(23)=lrf(ig)%phi3(23)+field2*(5.*dr(3)*dr(2)*dr(2)-r2*dr(3))
          lrf(ig)%phi3(24)=lrf(ig)%phi3(24)+field2*(5.*dr(3)*dr(2)*dr(3)-r2*dr(2))
          lrf(ig)%phi3(25)=lrf(ig)%phi3(25)+field2*(5.*dr(3)*dr(3)*dr(1)-r2*dr(1))
          lrf(ig)%phi3(26)=lrf(ig)%phi3(26)+field2*(5.*dr(3)*dr(3)*dr(2)-r2*dr(2))
          lrf(ig)%phi3(27)=lrf(ig)%phi3(27)+field2*(5.*dr(3)*dr(3)*dr(3)-r2*3.*dr(3))
        enddo
  end if

  ! if out of space then make more space
  if (nbww_pair .ge. calculation_assignment%ww%max) call reallocate_nonbondlist_ww

end do jwloop
 !now mark the end of the list of molecules interacting with molecule iw
 !by means of a zero element in the list
 nbww_pair = nbww_pair + 1
 nbww(nbww_pair) = 0    

end do iwloop
#if defined (profiling)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_box_lrf
!---------------------------------------------------------------------------
subroutine nbmonitorlist
!set lj code for the atom pairs in the selected atom groups to be monitored
! local variables
integer         :: i,j,ig,jg,ia,ja,i3,j3,nl,istate,lj_code,maxingroups,par, atomnri
integer         :: grpi,grpj,atomi,atomj,qq_pair,alj,blj


if (monitor_group_pairs == 0) return

!check the size of the largest group
maxingroups=maxval(monitor_atom_group(:)%n)
allocate (special_ljcod(nstates,maxingroups,maxingroups,monitor_group_pairs))

do par=1,monitor_group_pairs
    grpi=monitor_group_pair(par)%i      
    grpj=monitor_group_pair(par)%j 
    do i=1,monitor_atom_group(grpi)%n
      atomi=monitor_atom_group(grpi)%atom(i)
      do j=1,monitor_atom_group(grpj)%n
              atomj=monitor_atom_group(grpj)%atom(j)    
              !assert that atoms are different
              if(atomi == atomj) then 
                      call die('two paired monitor atom groups contain the same atom')
              end if
              do istate=1,nstates
                      !starting guess = use lj_code matrix for topology atom types
                      lj_code = ljcod(iac(atomi),iac(atomj))
                      if(qvdw_flag .and. lj_code == 2) then
                         if((iqatom(atomi) /= 0 .and. iqatom(atomj) == 0) .or.  &
                            (iqatom(atomj)/= 0 .and. iqatom(atomi) == 0)) then
                            !can't use code 2 between q and non-q when q-atom
                            !types are used. q-type 2 params are for exp. 
                            !repulsion, not lj!
                            lj_code = 1
                          end if
                      end if
                      !are atoms of pair in 1-4 position?
                      if(iqatom(atomi) == 0 .and. iqatom(atomj) == 0) then !neither atom is q_atom
                        if (abs(atomj-atomi) .le. max_nbr_range ) then
                          if (atomi .gt. atomj ) then
                            if ( list14(atomi-atomj, atomi) ) lj_code = 3 !3 means 1-4
                          else
                            if ( list14(atomj-atomi, atomj) ) lj_code = 3
                          end if
                        else
                          do nl = 1, n14long
                            if ((list14long(1, nl) .eq. atomi &
                              .and. list14long(2, nl) .eq. atomj) .or. &
                               (list14long(1, nl) .eq. atomj &
                              .and. list14long(2, nl) .eq. atomi)) then
                                  lj_code = 3
                            endif
                          end do
                        endif   !kolla 1-4 interaktioner
                      else      !at least one is q-atom
                            !check q-q pairlist to find lj-code
                        do qq_pair = 1, nbqq_pair(istate)
                          atomnri=iqseq(nbqq(qq_pair,istate)%iq)         !find atom number from qatom number
                          if((atomnri == atomi .and. &     
                              nbqq(qq_pair,istate)%j == atomj ) .or. &
                             (nbqq(qq_pair,istate)%j == atomi .and. &
                              atomnri == atomj)) then   
                              lj_code = nbqq(qq_pair,istate)%ljcod
                                 exit   
                          end if
                        end do
                           !if not found here then the first guess should be used
                      end if
                                special_ljcod(istate,i,j,par)=lj_code
              end do !nstates
      end do   ! monitor_atom_group j
    end do  !monitor_atom_group i
end do  !par
end subroutine nbmonitorlist

!---------------------------------------------------------------------------------------------------------                 

subroutine nonbond_monitor
!monitor nonbonded energies between selected groups of atoms

real(8)  :: dx1,dx2,dx3,r,r2,r6
integer  :: i,j,istate,lj_code,par
integer  :: grpi,grpj,atomi,atomj,qatomi,qatomj,qq_pair, iaci,iacj
real(8)  :: alji,blji,aljj,bljj,qi,qj, vel,vvdw,vwel,vwvdw,vwsum
real(8)  :: r6_hc    !  softcore variables
integer  :: sc_1,sc_2     !  softcore variables, sc_1 is the first index in sc_lookup (the qatom)
logical  :: do_sc    !  softcore variables,   do_sc is a boolean to determine if softcore should be done
                                        ! do_sc is true when atom i or j is a qatom  (and qvdw is true)


do par=1,monitor_group_pairs

        grpi=monitor_group_pair(par)%i      
        grpj=monitor_group_pair(par)%j 
        monitor_group_pair(par)%vel(:)=0
        monitor_group_pair(par)%vlj(:)=0
        monitor_group_pair(par)%vwel = 0
        monitor_group_pair(par)%vwlj = 0
        monitor_group_pair(par)%vwsum= 0     

        do i=1,monitor_atom_group(grpi)%n
            
                        atomi=monitor_atom_group(grpi)%atom(i)
            qatomi = iqatom(atomi)
                        qi   = crg(atomi)
                        iaci = iac(atomi)
            
                            do j=1,monitor_atom_group(grpj)%n
            
                                    atomj=monitor_atom_group(grpj)%atom(j)  
                        qatomj = iqatom(atomj)
                        iacj = iac(atomj)
                        qj   = crg(atomj)
                        dx1  = x(3*atomi-2)-x(3*atomj-2)      ! calculate the distance
                        dx2  = x(3*atomi-1)-x(3*atomj-1)
                        dx3  = x(3*atomi)-x(3*atomj)
                                                
                                                
                                                if (use_pbc) then
                                                        dx1 = dx1 - boxlength(1)*nint( dx1*inv_boxl(1) )
                                                        dx2 = dx2 - boxlength(2)*nint( dx2*inv_boxl(2) )
                                                        dx3 = dx3 - boxlength(3)*nint( dx3*inv_boxl(3) )
                                                end if

                        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
                        r    = sqrt(1/r2) 
                        r6   = r2*r2*r2
                                                r6_hc= r6   !needed for softcore
                                                r6   = 1._8/r6
            
                                    do istate=1,nstates                 
                                                                do_sc = .false.  !default is no softcore

                                                                lj_code  = special_ljcod(istate,i,j,par)
                                alji=iaclib(iaci)%avdw(lj_code)
                                blji=iaclib(iaci)%bvdw(lj_code)
                                aljj=iaclib(iacj)%avdw(lj_code)
                                bljj=iaclib(iacj)%bvdw(lj_code)   
                                if (qatomi /= 0) qi = qcrg(qatomi,istate) 
                                if (qatomj /= 0) qj = qcrg(qatomj,istate) 
                                if (qvdw_flag) then
                                        if (qatomi/=0) then
                                                iaci = qiac(qatomi,istate)
                                                alji = qavdw(iaci,lj_code)
                                                blji = qbvdw(iaci,lj_code)


                                                                                                do_sc = .true.   ! atom i is a q-atom, softcore on
                                                                                                sc_1 = qatomi    ! the first index in sc_lookup should be a qatom number

                                                                                else
                                                                                                
                                                                                                sc_2 = iaci      ! atom i was not a qatom, put atom code i in the second sc_lookup index

                                        endif   


                                        if (qatomj/=0) then                                     
                                                iacj  = qiac(qatomj,istate)                     
                                                aljj  = qavdw(iacj,lj_code)                     
                                                bljj  = qbvdw(iacj,lj_code)

                                                                                                if (do_sc) then   ! do_sc is true if atom i is a qatom
                                                                                                        sc_2 = qatomj + natyps  ! qatomi is sc_1
                                                                                                else
                                                                                                        do_sc = .true.  ! atom i was not a qatom but j is, softcore on
                                                                                                        sc_1 = qatomj   ! qatom j should be the first index
                                                                                                end if

                                                                                else
                                                                                        sc_2 = iacj    ! atom j is not a qatom, should be index 2
                                                                                
                                                                                                
                                        endif
                                                                                
                                                                                if (do_sc) then  ! calculate softcore r6
                                                                                        r6 = r6_hc + sc_lookup (sc_1,sc_2,istate)
                                                                                        r6 = 1._8/r6
                                                                                end if   
                                endif                                                                           
                                vel  = qi*qj*r  
                                if(ivdw_rule==1) then !geometric comb. rule
                                        vvdw = alji*aljj*r6*r6-blji*bljj*r6
                                else !arithmetic
                                        vvdw = blji * bljj * (alji+aljj)**6 * r6 * &
                                                ((alji+aljj)**6 * r6 - 2.0)
                                endif
                                !add up for this pair of atom groups
                                monitor_group_pair(par)%vel(istate)= monitor_group_pair(par)%vel(istate)+vel
                                monitor_group_pair(par)%vlj(istate)= monitor_group_pair(par)%vlj(istate)+vvdw
                        end do  !nstates
                end do   ! monitor_atom_group j
        end do  !monitor_atom_group i
        !calc lambda-weighted sum
        monitor_group_pair(par)%vwel=dot_product(monitor_group_pair(par)%vel(1:nstates),eq(1:nstates)%lambda)
        monitor_group_pair(par)%vwlj=dot_product(monitor_group_pair(par)%vlj(1:nstates),eq(1:nstates)%lambda)
        monitor_group_pair(par)%vwsum= monitor_group_pair(par)%vwlj+monitor_group_pair(par)%vwel
end do !par

end subroutine nonbond_monitor

!-------------------------------------------------------------------------

subroutine nonbon2_pp
! local variables
integer                                         :: ip
real(8)                                         :: alja,blja,dx1a,dx2a,dx3a,r2a,ra,r6a
real(8)                                         :: aljb,bljb,dx1b,dx2b,dx3b,r2b,rb,r6b
real(8)                                         :: vela,v_aa,v_ba,dva
real(8)                                         :: velb,v_ab,v_bb,dvb
type(nb_type), pointer          :: pa
type(nb_type), pointer          :: pb

! global variables used:
!  iaclib, x, crg, el14_scale, d, e

do ip = 1, nbpp_pair - 1, 2
! for every second pair:

! init pointers
pa => nbpp(ip)
pb => nbpp(ip+1)

! calculate alj and blj
alja  = iaclib(iac(pa%i))%avdw(pa%ljcod)+iaclib(iac(pa%j))%avdw(pa%ljcod)
aljb  = iaclib(iac(pb%i))%avdw(pb%ljcod)+iaclib(iac(pb%j))%avdw(pb%ljcod)
blja  = iaclib(iac(pa%i))%bvdw(pa%ljcod)*iaclib(iac(pa%j))%bvdw(pa%ljcod)
bljb  = iaclib(iac(pb%i))%bvdw(pb%ljcod)*iaclib(iac(pb%j))%bvdw(pb%ljcod)
alja  = alja*alja
aljb  = aljb*aljb
alja  = alja*alja*alja
aljb  = aljb*aljb*aljb

! calculate dx, r and r2
dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
dx3b  = x(pb%j*3-0) - x(pb%i*3-0)

r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
ra = sqrt(r2a)
rb = sqrt(r2b) 

r6a   = r2a*r2a*r2a
r6b   = r2b*r2b*r2b

! calculate vel and dv
vela  = crg(pa%i)*crg(pa%j)*ra  
velb  = crg(pb%i)*crg(pb%j)*rb  
if ( pa%ljcod .eq. 3 ) then
  vela = vela*el14_scale
end if
if ( pb%ljcod .eq. 3 ) then
  velb = velb*el14_scale
end if
v_aa  = blja*alja*alja*r6a*r6a 
v_ab  = bljb*aljb*aljb*r6b*r6b 
v_ba  = 2.0*blja*alja*r6a
v_bb  = 2.0*bljb*aljb*r6b
dva   = r2a*( -vela -12.*v_aa +6.*v_ba )
dvb   = r2b*( -velb -12.*v_ab +6.*v_bb )

! update d
d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

! update energies
e%pp%el  = e%pp%el + vela + velb
e%pp%vdw = e%pp%vdw + v_aa - v_ba + v_ab - v_bb
end do

if (ip .eq. nbpp_pair) then
! the last pair:

pa => nbpp(ip)
alja  = iaclib(iac(pa%i))%avdw(pa%ljcod)+iaclib(iac(pa%j))%avdw(pa%ljcod)
blja  = iaclib(iac(pa%i))%bvdw(pa%ljcod)*iaclib(iac(pa%j))%bvdw(pa%ljcod)
alja  = alja*alja
alja  = alja*alja*alja

dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)

r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
ra = sqrt(r2a)

r6a   = r2a*r2a*r2a

vela  = crg(pa%i)*crg(pa%j)*ra
if ( pa%ljcod .eq. 3 ) vela = vela*el14_scale

v_aa  = blja*alja*alja*r6a*r6a 
v_ba  = 2.0*blja*alja*r6a
dva   = r2a*( -vela -12.*v_aa +6.*v_ba )

d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

e%pp%el  = e%pp%el + vela 
e%pp%vdw = e%pp%vdw + v_aa - v_ba 
end if

end subroutine nonbon2_pp

!------------------------------------------------------------------------
subroutine nonbon2_pp_box
  ! local variables
  integer                                               :: ip, ga, gb, group
  real(8)                                               :: alja,blja,dx1a,dx2a,dx3a,r2a,ra,r6a
  real(8)                                               :: aljb,bljb,dx1b,dx2b,dx3b,r2b,rb,r6b
  real(8)                                               :: vela,v_aa,v_ba,dva
  real(8)                                               :: velb,v_ab,v_bb,dvb
  type(nb_type), pointer                :: pa
  type(nb_type), pointer                :: pb

 ! global variables used:
  !  iaclib, x, crg, el14_scale, d, e


        do group = 1, nbpp_cgp_pair
                ga = nbpp_cgp(group)%i !atom index for the two switching atoms
                gb = nbpp_cgp(group)%j

                !the distance between the two switching atoms
                dx1a = x(3*gb-2) - x(3*ga-2)
                dx2a = x(3*gb-1) - x(3*ga-1)
                dx3a = x(3*gb  ) - x(3*ga  )

                nbpp_cgp(group)%x = boxlength(1)*nint( dx1a*inv_boxl(1) )
                nbpp_cgp(group)%y = boxlength(2)*nint( dx2a*inv_boxl(2) )       
                nbpp_cgp(group)%z = boxlength(3)*nint( dx3a*inv_boxl(3) )

        end do

  do ip = 1, nbpp_pair - 1, 2
    ! for every second pair:

        ! init pointers
        pa => nbpp(ip)
        pb => nbpp(ip+1)
        ga = pa%cgp_pair
        gb = pb%cgp_pair

        ! calculate alj and blj
    alja  = iaclib(iac(pa%i))%avdw(pa%ljcod)+iaclib(iac(pa%j))%avdw(pa%ljcod)
    aljb  = iaclib(iac(pb%i))%avdw(pb%ljcod)+iaclib(iac(pb%j))%avdw(pb%ljcod)
    blja  = iaclib(iac(pa%i))%bvdw(pa%ljcod)*iaclib(iac(pa%j))%bvdw(pa%ljcod)
    bljb  = iaclib(iac(pb%i))%bvdw(pb%ljcod)*iaclib(iac(pb%j))%bvdw(pb%ljcod)
    alja  = alja*alja
    aljb  = aljb*aljb
    alja  = alja*alja*alja
    aljb  = aljb*aljb*aljb

        ! calculate dx, r and r2
    dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
    dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
    dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
    dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
    dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
    dx3b  = x(pb%j*3-0) - x(pb%i*3-0)
        dx1a = dx1a - nbpp_cgp(ga)%x         
        dx1b = dx1b - nbpp_cgp(gb)%x 
        dx2a = dx2a - nbpp_cgp(ga)%y    
        dx2b = dx2b - nbpp_cgp(gb)%y                
        dx3a = dx3a - nbpp_cgp(ga)%z  
        dx3b = dx3b - nbpp_cgp(gb)%z   
        


    r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
    r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
        ra = sqrt(r2a)
        rb = sqrt(r2b) 

    r6a   = r2a*r2a*r2a
    r6b   = r2b*r2b*r2b

        ! calculate vel and dv
    vela  = crg(pa%i)*crg(pa%j)*ra  
    velb  = crg(pb%i)*crg(pb%j)*rb  
    if ( pa%ljcod .eq. 3 ) then
          vela = vela*el14_scale
        end if
    if ( pb%ljcod .eq. 3 ) then
          velb = velb*el14_scale
        end if
    v_aa  = blja*alja*alja*r6a*r6a 
    v_ab  = bljb*aljb*aljb*r6b*r6b 
    v_ba  = 2.0*blja*alja*r6a
    v_bb  = 2.0*bljb*aljb*r6b
    dva   = r2a*( -vela -12.*v_aa +6.*v_ba )
    dvb   = r2b*( -velb -12.*v_ab +6.*v_bb )

        ! update d
    d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
    d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
    d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
    d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
    d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
    d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
    d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
    d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
    d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
    d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
    d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
    d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

        ! update energies
    e%pp%el  = e%pp%el + vela + velb
    e%pp%vdw = e%pp%vdw + v_aa - v_ba + v_ab - v_bb
  end do

  if (ip .eq. nbpp_pair) then
    ! the last pair:

        pa => nbpp(ip)
        ga = pa%cgp_pair

    alja  = iaclib(iac(pa%i))%avdw(pa%ljcod)+iaclib(iac(pa%j))%avdw(pa%ljcod)
    blja  = iaclib(iac(pa%i))%bvdw(pa%ljcod)*iaclib(iac(pa%j))%bvdw(pa%ljcod)
    alja  = alja*alja
    alja  = alja*alja*alja

    dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
    dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
    dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
        dx1a = dx1a - nbpp_cgp(ga)%x  
        dx2a = dx2a - nbpp_cgp(ga)%y  
        dx3a = dx3a - nbpp_cgp(ga)%z  


    r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
        ra = sqrt(r2a)
    r6a   = r2a*r2a*r2a

   vela  = crg(pa%i)*crg(pa%j)*ra
    if ( pa%ljcod .eq. 3 ) vela = vela*el14_scale
    v_aa  = blja*alja*alja*r6a*r6a 
    v_ba  = 2.0*blja*alja*r6a
    dva   = r2a*( -vela -12.*v_aa +6.*v_ba )

    d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
    d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
    d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
    d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
    d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
    d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

    e%pp%el  = e%pp%el + vela 
    e%pp%vdw = e%pp%vdw + v_aa - v_ba 
  end if

end subroutine nonbon2_pp_box
!----------------------------------------------------------------------

subroutine nonbon2_pw
! local variables
integer                                         :: ip,i,j,i3,j3,iaci,iacj,ilj
real(8)                                         :: alj,blj,dx1,dx2,dx3,r2,r,r6,r12
real(8)                                         :: vel,v_a,v_b,dv

! global variables used:
!  iac, crg, iaclib, x, d, e

do ip = 1, nbpw_pair
! for every assigned pair:
i    = nbpw(ip)%i
j    = nbpw(ip)%j
i3   = i*3-3
j3   = j*3-3
iaci = iac(i)
iacj = iac(j)
ilj  = nbpw(ip)%ljcod
crg(i)   = crg(i)
crg(j)   = crg(j)
alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
alj  = alj*alj
alj  = alj*alj*alj

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)

r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = 1./r2
r    = sqrt ( r2 ) 
r6   = r2*r2*r2
r12  = r6*r6

! calculate vel and dv
vel  = crg(i)*crg(j)*r
v_a  = blj*alj*alj*r12 
v_b  = 2.0*blj*alj*r6
dv   = r2*( -vel -12.*v_a +6.*v_b )

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
e%pw%el  = e%pw%el + vel       
e%pw%vdw = e%pw%vdw + v_a - v_b
end do

end subroutine nonbon2_pw

!-----------------------------------------------------------------------

subroutine nonbon2_pw_box
  ! local variables
  integer                                               :: ip,i,j,i3,j3,iaci,iacj,ilj
  real(8)                                               :: alj,blj,dx1,dx2,dx3,r2,r,r6,r12
  real(8)                                               :: vel,v_a,v_b,dv
  integer                                               :: group, ga, gb

  ! global variables used:
  !  iac, crg, iaclib, x, d, e


   !compute the peridocal shift for every charge group pair
        do group = 1, nbpw_cgp_pair
                ga = nbpw_cgp(group)%i !atom index for solute switching atom
                gb = nbpw_cgp(group)%j  !atom index for the solvent switching atom

                !the distance between the two switching atoms
                dx1 = x(3*gb-2) - x(3*ga-2)
                dx2 = x(3*gb-1) - x(3*ga-1)
                dx3 = x(3*gb  ) - x(3*ga  )

                nbpw_cgp(group)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
                nbpw_cgp(group)%y = boxlength(2)*nint( dx2*inv_boxl(2) )        
                nbpw_cgp(group)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

        end do


  do ip = 1, nbpw_pair
        ! for every assigned pair:

        i    = nbpw(ip)%i !solute atom 
        j    = nbpw(ip)%j !solvent atom
        group = nbpw(ip)%cgp_pair
        i3   = i*3-3
        j3   = j*3-3
        iaci = iac(i)
        iacj = iac(j)
        ilj  = nbpw(ip)%ljcod
        crg(i)   = crg(i)
        crg(j)   = crg(j)
        alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        alj  = alj*alj
        alj  = alj*alj*alj

        ! calculate dx and r
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        dx1 = dx1 - nbpw_cgp(group)%x  
        dx2 = dx2 - nbpw_cgp(group)%y  
        dx3 = dx3 - nbpw_cgp(group)%z  
        

        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6

        ! calculate vel and dv
        vel  = crg(i)*crg(j)*r
        v_a  = blj*alj*alj*r12 
        v_b  = 2.0*blj*alj*r6
        dv   = r2*( -vel -12.*v_a +6.*v_b )

        ! update forces
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3

        ! update energies
        e%pw%el  = e%pw%el + vel       
        e%pw%vdw = e%pw%vdw + v_a - v_b
  end do

end subroutine nonbon2_pw_box

!-----------------------------------------------------------------------

subroutine nonbon2_qq
! local variables
integer                                         :: istate
integer                                         :: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,ilj
real(8)                                         :: qi,qj,alj,blj,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)                                         :: vel,v_a,v_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq
  i    = iqseq(iq)
  j    = nbqq(ip,istate)%j
  jq   = nbqq(ip,istate)%jq
  i3   = i*3-3
  j3   = j*3-3
  ilj  = nbqq(ip,istate)%ljcod
  qi   = qcrg(iq,istate)
  el_scale = nbqq(ip,istate)%el_scale


  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        if ( jq /= 0) then
          qj = qcrg(jq,istate)
        else
          qj = crg(j)
        end if
  else
        iaci = qiac(iq,istate)
        alj  = qavdw(iaci,ilj)
        blj  = qbvdw(iaci,ilj)
        if ( jq /= 0) then
          iacj = qiac(jq,istate)
          qj   = qcrg(jq,istate)
          if ( ilj .eq. 2 ) then
                alj = alj*qavdw(iacj,ilj)
          else
                alj = alj+qavdw(iacj,ilj)
          end if
          blj = blj*qbvdw(iacj,ilj)
        else
          if ( ilj .eq. 2 ) then
                alj  = qavdw(iaci,1)
                blj  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          alj = alj+iaclib(iacj)%avdw(ilj)
          blj = blj*iaclib(iacj)%bvdw(ilj)
          qj = crg(j)
        end if
  end if

  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)

  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
  r6_hc = r2*r2*r2  !for softcore
  r6   = r6_hc + sc_lookup(iq,jq+natyps,istate)  !softcore
  r6   = 1._8/r6
  r2   = 1./r2
  r    = sqrt ( r2 ) 
  r12  = r6*r6

  ! calculate vel, v_a, v_b and dv
  vel  = qi*qj*r*el_scale
  if ( ilj .eq. 3 ) vel = vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. ilj .eq. 2 ) then
        v_a = alj*exp(-blj/r)
        v_b = 0.0
        dv  = r2*( -vel -blj*v_a/r )*eq(istate)%lambda
  else
        alj  = alj*alj
        alj  = alj*alj*alj
        v_a  = blj*alj*alj*r12 
        v_b  = 2.0*blj*alj*r6
        dv  = r2*( -vel -(12.*v_a -6.*v_b)*r6*r6_hc )*eq(istate)%lambda
  endif

  ! update forces
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        eq(istate)%qq%el  = eq(istate)%qq%el + vel
        eq(istate)%qq%vdw = eq(istate)%qq%vdw + v_a - v_b
  else
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
  end if
end do ! ip

end do ! istate
end subroutine nonbon2_qq

!-----------------------------------------------------------------------
subroutine nonbon2_qq_lib_charges
! local variables
integer                                         :: istate
integer                                         :: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,ilj
real(8)                                         :: qi,qj,alj,blj,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)                                         :: vel,v_a,v_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq
  i    = iqseq(iq)
  j    = nbqq(ip,istate)%j
  jq   = nbqq(ip,istate)%jq
  i3   = i*3-3
  j3   = j*3-3
  ilj  = nbqq(ip,istate)%ljcod
  qi   = crg(i)
  el_scale = nbqq(ip,istate)%el_scale

  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        qj = crg(j)
  else
        iaci = qiac(iq,istate)
        alj  = qavdw(iaci,ilj)
        blj  = qbvdw(iaci,ilj)
        if ( jq /= 0) then
          iacj = qiac(jq,istate)
          if ( ilj .eq. 2 ) then
                alj = alj*qavdw(iacj,ilj)
          else
                alj = alj+qavdw(iacj,ilj)
          end if
          blj = blj*qbvdw(iacj,ilj)
        else
          if ( ilj .eq. 2 ) then
                alj  = qavdw(iaci,1)
                blj  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          alj = alj+iaclib(iacj)%avdw(ilj)
          blj = blj*iaclib(iacj)%bvdw(ilj)
        end if
        qj = crg(j)
  end if

  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)

  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
  r6_hc = r2*r2*r2   !needed for softcore
  r6   = r6_hc + sc_lookup(iq,jq+natyps,istate)  !softcore
  r6   = 1._8/r6
  r12  = r6*r6
  r2   = 1./r2
  r    = sqrt ( r2 ) 

  ! calculate vel, v_a, v_b and dv
  vel  = qi*qj*r*el_scale
  if ( ilj .eq. 3 ) vel = vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. ilj .eq. 2 ) then
        v_a = alj*exp(-blj/r)
        v_b = qbvdw(iaci,1)*qbvdw(iacj,1)*r6
        dv  = r2*( -vel -blj*v_a/r +6.*v_b )*eq(istate)%lambda
        !
        ! ---       change here to exclude 1/r6 attraction
        !
        !              v_b = 0.0
        !              dv  = r2*( -vel -blj*v_a/r )*eq(istate)%lambda
  else
        alj  = alj*alj
        alj  = alj*alj*alj
        v_a  = blj*alj*alj*r12 
        v_b  = 2.0*blj*alj*r6
        dv  = r2*( -vel -(12.*v_a -6.*v_b)*r6*r6_hc )*eq(istate)%lambda
  endif

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        eq(istate)%qq%el  = eq(istate)%qq%el + vel
        eq(istate)%qq%vdw = eq(istate)%qq%vdw + v_a - v_b
  else
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
  end if
end do ! ip

end do ! istate
end subroutine nonbon2_qq_lib_charges

!-----------------------------------------------------------------------

subroutine nonbon2_qp
! local variables
integer                                         :: ip,iq,i,j,i3,j3,iaci,iacj,ilj
integer                                         :: istate
real(8)                                         :: alj,blj,dx1,dx2,dx3,r2,r,r6,r6_hc
real(8)                                         :: vel,v_a,v_b,dv

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute

do ip = 1, nbqp_pair
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip)%i
i    = iqseq(iq)
j    = nbqp(ip)%j
i3   = i*3-3
j3   = j*3-3
iacj = iac(j)
ilj  = nbqp(ip)%ljcod


dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)


r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

r6   = r2*r2*r2
r6_hc = r6     !for softcore

r2   = 1._8/r2
r    = sqrt(r2)



do istate = 1, nstates
  ! for every state:

  ! calculate iaci, alj and blj
  if (.not. qvdw_flag) then
        iaci = iac(i)
        alj  = iaclib(iaci)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)
                r6 = 1._8/r6_hc
  else
        iaci = qiac(iq,istate)
                if ( ilj .eq. 2 ) then
          alj  = qavdw(iaci,1)
          blj  = qbvdw(iaci,1)
        else
          alj  = qavdw(iaci,ilj)
          blj  = qbvdw(iaci,ilj)
                end if

                r6 = r6_hc + sc_lookup(iq,iacj,istate) !this is softcore
                r6 = 1._8/r6
  end if
  alj = alj+iaclib(iacj)%avdw(ilj)
  blj = blj*iaclib(iacj)%bvdw(ilj)
  alj  = alj*alj
  alj  = alj*alj*alj

  ! calculate qi, vel, v_a, v_b and dv
  vel  = qcrg(iq,istate)*crg(j)*r
  if ( ilj .eq. 3 ) vel = vel*el14_scale
  v_a  = blj*alj*alj*r6*r6
  v_b  = 2.0*blj*alj*r6
  dv   = r2*( -vel -(12.*v_a -6.*v_b)*r6*r6_hc )*eq(istate)%lambda   !softcore r6*r6_hc is (r^6/(r^6+alpha))

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein or q-water energies
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
end do ! istate

end do
end subroutine nonbon2_qp
!-----------------------------------------------------------------------

!******pwadded 2001-10-23
subroutine nonbon2_qp_box
  ! local variables
  integer                                               :: ip,iq,i,j,i3,j3,iaci,iacj,ilj
  integer                                               :: istate
  real(8)                                               :: alj,blj,dx1,dx2,dx3,r2,r,r6,r6_hc
  real(8)                                               :: vel,v_a,v_b,dv
  integer                                               :: group, gr, ia

  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute


        !compute the peridocal shift for every charge group pair
        do gr = 1, nbqp_cgp_pair
                ia = nbqp_cgp(gr)%i !atom index for the atom
                
                !the distance between the two switching atoms
                dx1 = x(3*ia-2) - x(3*qswitch-2)
                dx2 = x(3*ia-1) - x(3*qswitch-1)
                dx3 = x(3*ia  ) - x(3*qswitch  )

                nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
                nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )   
                nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

        end do

  do ip = 1, nbqp_pair
        ! for every assigned q-s pair:

        ! init state-invariant variables:
        iq   = nbqp(ip)%i
        i    = iqseq(iq)
        j    = nbqp(ip)%j
        i3   = i*3-3
        j3   = j*3-3
        iacj = iac(j)
        ilj  = nbqp(ip)%ljcod
        group = nbqp(ip)%cgp_pair

        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        dx1 = dx1 - nbqp_cgp(group)%x  
        dx2 = dx2 - nbqp_cgp(group)%y  
        dx3 = dx3 - nbqp_cgp(group)%z  
        

        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r6_hc = r2*r2*r2  !for softcore

        r2   = 1._8/r2
        r    = sqrt(r2)

        do istate = 1, nstates
          ! for every state:

          ! calculate iaci, alj and blj
          if (.not. qvdw_flag) then
                iaci = iac(i)
                alj  = iaclib(iaci)%avdw(ilj)
                blj  = iaclib(iaci)%bvdw(ilj)
                r6 = 1._8/r6_hc
          else
                iaci = qiac(iq,istate)
        if ( ilj .eq. 2 ) then
                  alj  = qavdw(iaci,1)
                  blj  = qbvdw(iaci,1)
                else
                  alj  = qavdw(iaci,ilj)
                  blj  = qbvdw(iaci,ilj)
        end if
                r6 = r6_hc + sc_lookup(iq,iacj,istate)
                r6 = 1._8/r6
          end if
          alj = alj+iaclib(iacj)%avdw(ilj)
          blj = blj*iaclib(iacj)%bvdw(ilj)
          alj  = alj*alj
          alj  = alj*alj*alj

          ! calculate qi, vel, v_a, v_b and dv
          vel  = qcrg(iq,istate)*crg(j)*r
          if ( ilj .eq. 3 ) vel = vel*el14_scale
          v_a  = blj*alj*alj*r6*r6
          v_b  = 2.0*blj*alj*r6
          dv   = r2*( -vel -(12.*v_a -6.*v_b)*r6*r6_hc )*eq(istate)%lambda

          ! update forces
          d(i3+1) = d(i3+1) - dv*dx1
          d(i3+2) = d(i3+2) - dv*dx2
          d(i3+3) = d(i3+3) - dv*dx3
          d(j3+1) = d(j3+1) + dv*dx1
          d(j3+2) = d(j3+2) + dv*dx2
          d(j3+3) = d(j3+3) + dv*dx3
        
          ! update q-protein or q-water energies
                eq(istate)%qp%el  = eq(istate)%qp%el + vel
                eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
        end do ! istate

  end do
end subroutine nonbon2_qp_box

!-----------------------------------------------------------------------

subroutine nonbon2_qw

! local variables
integer                                         :: jw,iq,i,j,iljo, iljh, iaci
integer                                         :: istate
real(8)                                         ::      aljo, bljo, aljh, bljh
real(8)                                         ::      dxo, dyo, dzo, dxh1, dyh1, dzh1, dxh2, dyh2, dzh2
real(8)                                         ::      ro, r2o, r6o, rh1, r2h1, r6h1, rh2, r2h2, r6h2,r6o_hc,r6h1_hc,r6h2_hc
real(8)                                         ::      velo, velh1, velh2, dvo, dvh1, dvh2
real(8)                                         :: v_ao, v_bo, v_ah1, v_bh1, v_ah2, v_bh2 
real(8), save                           ::      ao(2), bo(2), ah(2), bh(2)
integer, save                           ::      iac_ow, iac_hw
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute

if(iac_ow == 0) then !set first time
        iac_ow = iac(nat_solute + 1)
        iac_hw = iac(nat_solute + 2)
        ao(1:2) = iaclib(iac_ow)%avdw(1:2)
        bo(1:2) = iaclib(iac_ow)%bvdw(1:2)
        ah(1:2) = iaclib(iac_hw)%avdw(1:2)
        bh(1:2) = iaclib(iac_hw)%bvdw(1:2)
end if

!loop over listed waters
do jw = 1, nbqw_pair
        j = nbqw(jw) !top. # of o in water iw
        !loop over all q-atoms
        do iq = 1, nqat
                i = iqseq(iq)
                dxo  = x(3*j-2) - x(3*i-2)
                dyo  = x(3*j-1) - x(3*i-1)
                dzo  = x(3*j  ) - x(3*i  )
                dxh1 = x(3*j+1) - x(3*i-2)
                dyh1 = x(3*j+2) - x(3*i-1)
                dzh1 = x(3*j+3) - x(3*i  )
                dxh2 = x(3*j+4) - x(3*i-2)
                dyh2 = x(3*j+5) - x(3*i-1)
                dzh2 = x(3*j+6) - x(3*i  )
                r2o  = 1._8/(dxo*dxo + dyo*dyo + dzo*dzo)
                r2h1 = 1._8/(dxh1*dxh1 + dyh1*dyh1 + dzh1*dzh1)
                r2h2 = 1._8/(dxh2*dxh2 + dyh2*dyh2 + dzh2*dzh2)
                ro = sqrt(r2o)
                r6o = r2o*r2o*r2o
                                r6o_hc = r6o    !softcore
                rh1 = sqrt(r2h1)
                r6h1 = r2h1*r2h1*r2h1
                                r6h1_hc = r6h1  !softcore
                rh2 = sqrt(r2h2)
                r6h2 = r2h2*r2h2*r2h2
                                r6h2_hc = r6h2  !softcore

                iaci = iac(i)

                !reset potential
                dvo = 0
                dvh1 = 0
                dvh2 = 0
                iljo = ljcod(iac_ow, iaci)
                iljh = ljcod(iac_hw, iaci)
                if(.not. qvdw_flag) then
                        !use same lj params for all states
                        aljo  = iaclib(iaci)%avdw(iljo)+ao(iljo)
                        bljo  = iaclib(iaci)%bvdw(iljo)*bo(iljo)
                        aljh  = iaclib(iaci)%avdw(iljh)+ah(iljh)
                        bljh  = iaclib(iaci)%bvdw(iljh)*bh(iljh) 
                        aljo = aljo * aljo
                        aljo = aljo * aljo * aljo
                        aljh  = aljh  * aljh 
                        aljh  = aljh  * aljh  * aljh 
                        v_ao = bljo*aljo*aljo*r6o*r6o 
                        v_bo = 2.0*bljo*aljo*r6o
                        v_ah1= bljh*aljh*aljh*r6h1*r6h1
                        v_bh1= 2.0*bljh*aljh*r6h1
                        v_ah2= bljh*aljh*aljh*r6h2*r6h2
                        v_bh2= 2.0*bljh*aljh*r6h2
                end if
                do istate = 1, nstates ! for every state:

                        ! set new lj params if q-atom types are used
                        if (qvdw_flag) then

                                                                r6o = 1._8/r6o_hc
                                                                r6o = r6o + sc_lookup(iq,iac_ow,istate)   !softcore
                                                                r6o = 1._8/r6o
                                                                r6h1 = 1._8/r6h1_hc
                                                                r6h1 = r6h1 + sc_lookup(iq,iac_hw,istate)   !softcore
                                                                r6h1 = 1._8/r6h1
                                                                r6h2 = 1._8/r6h2_hc
                                                                r6h2 = r6h2 + sc_lookup(iq,iac_hw,istate)   !softcore
                                                                r6h2 = 1._8/r6h2
                                aljo  = qavdw(qiac(iq,istate),1)+ao(iljo)
                                bljo  = qbvdw(qiac(iq,istate),1)*bo(iljo)
                                aljh  = qavdw(qiac(iq,istate),1)+ah(iljh)
                                bljh  = qbvdw(qiac(iq,istate),1)*bh(iljh)
                                aljo = aljo * aljo
                                aljo = aljo * aljo * aljo
                                aljh  = aljh  * aljh 
                                aljh  = aljh  * aljh  * aljh 
                                v_ao = bljo*aljo*aljo*r6o*r6o 
                                v_bo = 2.0*bljo*aljo*r6o
                                v_ah1= bljh*aljh*aljh*r6h1*r6h1
                                v_bh1= 2.0*bljh*aljh*r6h1
                                v_ah2= bljh*aljh*aljh*r6h2*r6h2
                                v_bh2= 2.0*bljh*aljh*r6h2
                        end if


                        ! calculate qi, vel, v_a, v_b and dv
                        velo = crg_ow*qcrg(iq,istate)*ro
                        velh1 = crg_hw*qcrg(iq,istate)*rh1
                        velh2 = crg_hw*qcrg(iq,istate)*rh2
                        dvo  = dvo  + r2o *( -velo  -(12.*v_ao  -6.*v_bo )*r6o/r6o_hc)*eq(istate)%lambda
                        dvh1 = dvh1 + r2h1*( -velh1 -(12.*v_ah1 -6.*v_bh1)*r6h1/r6h1_hc)*eq(istate)%lambda
                        dvh2 = dvh2 + r2h2*( -velh2 -(12.*v_ah2 -6.*v_bh2)*r6h2/r6h2_hc)*eq(istate)%lambda
                        ! update q-water energies
                        eq(istate)%qw%el  = eq(istate)%qw%el + velo + velh1 + velh2
                        eq(istate)%qw%vdw = eq(istate)%qw%vdw + v_ao + v_ah1 + v_ah2 - v_bo - v_bh1 - v_bh2 
                end do !istate

                ! update forces on q-atom
                d(3*i-2) = d(3*i-2) - dvo*dxo - dvh1*dxh1 - dvh2*dxh2
                d(3*i-1) = d(3*i-1) - dvo*dyo - dvh1*dyh1 - dvh2*dyh2
                d(3*i  ) = d(3*i  ) - dvo*dzo - dvh1*dzh1 - dvh2*dzh2

                ! update forces on water
                d(3*j-2) = d(3*j-2) + dvo*dxo
                d(3*j-1) = d(3*j-1) + dvo*dyo
                d(3*j  ) = d(3*j  ) + dvo*dzo
                d(3*j+1) = d(3*j+1) + dvh1*dxh1
                d(3*j+2) = d(3*j+2) + dvh1*dyh1
                d(3*j+3) = d(3*j+3) + dvh1*dzh1
                d(3*j+4) = d(3*j+4) + dvh2*dxh2
                d(3*j+5) = d(3*j+5) + dvh2*dyh2
                d(3*j+6) = d(3*j+6) + dvh2*dzh2

        end do !iq
end do !jw
end subroutine nonbon2_qw

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbon2_qw_box
        ! local variables
        integer                                         :: jw,iq,i,j,iljo, iljh, iaci
        integer                                         :: istate
        real(8)                                         ::      aljo, bljo, aljh, bljh
        real(8)                                         ::      dxo, dyo, dzo, dxh1, dyh1, dzh1, dxh2, dyh2, dzh2
        real(8)                                         ::      ro, r2o, r6o, rh1, r2h1, r6h1, rh2, r2h2, r6h2,r6o_hc,r6h1_hc,r6h2_hc
        real(8)                                         ::      velo, velh1, velh2, dvo, dvh1, dvh2
        real(8)                                         :: v_ao, v_bo, v_ah1, v_bh1, v_ah2, v_bh2
        real(8)                                         :: boxshiftx, boxshifty, boxshiftz, dx, dy, dz 
        real(8), save                           ::      ao(2), bo(2), ah(2), bh(2)
        integer, save                           ::      iac_ow, iac_hw
        ! global variables used:
        !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute
  
        if(iac_ow == 0) then !set first time
                iac_ow = iac(nat_solute + 1)
                iac_hw = iac(nat_solute + 2)
                ao(1:2) = iaclib(iac_ow)%avdw(1:2)
                bo(1:2) = iaclib(iac_ow)%bvdw(1:2)
                ah(1:2) = iaclib(iac_hw)%avdw(1:2)
                bh(1:2) = iaclib(iac_hw)%bvdw(1:2)
        end if

        !loop over listed waters
        do jw = 1, nbqw_pair
                j = nbqw(jw) !top. # of o in water iw

                !compute the periodical shift
                dx = x(3*j-2) - x(3*qswitch-2)
                dy = x(3*j-1) - x(3*qswitch-1)
                dz = x(3*j  ) - x(3*qswitch  )
                boxshiftx = boxlength(1)*nint(dx*inv_boxl(1))
                boxshifty = boxlength(2)*nint(dy*inv_boxl(2))
                boxshiftz = boxlength(3)*nint(dz*inv_boxl(3))

                !loop over all q-atoms
                do iq = 1, nqat
                        i = iqseq(iq)
                        dxo  = x(3*j-2) - x(3*i-2)
                        dyo  = x(3*j-1) - x(3*i-1)
                        dzo  = x(3*j  ) - x(3*i  )
                        dxh1 = x(3*j+1) - x(3*i-2)
                        dyh1 = x(3*j+2) - x(3*i-1)
                        dzh1 = x(3*j+3) - x(3*i  )
                        dxh2 = x(3*j+4) - x(3*i-2)
                        dyh2 = x(3*j+5) - x(3*i-1)
                        dzh2 = x(3*j+6) - x(3*i  )
                        dxo = dxo - boxshiftx
                        dyo = dyo - boxshifty  
                        dzo = dzo - boxshiftz  
                        dxh1 = dxh1 - boxshiftx  
                        dyh1 = dyh1 - boxshifty  
                        dzh1 = dzh1 - boxshiftz  
                        dxh2 = dxh2 - boxshiftx 
                        dyh2 = dyh2 - boxshifty  
                        dzh2 = dzh2 - boxshiftz  
                        r2o  = 1._8/(dxo*dxo + dyo*dyo + dzo*dzo)
                        r2h1 = 1._8/(dxh1*dxh1 + dyh1*dyh1 + dzh1*dzh1)
                        r2h2 = 1._8/(dxh2*dxh2 + dyh2*dyh2 + dzh2*dzh2)
                        ro = sqrt(r2o)
                        r6o = r2o*r2o*r2o
                        r6o_hc = r6o    !needed for softcore
                        rh1 = sqrt(r2h1)
                        r6h1 = r2h1*r2h1*r2h1
                        r6h1_hc = r6h1    !needed for softcore
                        rh2 = sqrt(r2h2)
                        r6h2 = r2h2*r2h2*r2h2
                        r6h2_hc = r6h2    !needed for softcore
                        iaci = iac(i)
                        
                        !reset potential
                        dvo = 0
                        dvh1 = 0
                        dvh2 = 0
                        iljo = ljcod(iac_ow, iaci)
                        iljh = ljcod(iac_hw, iaci)
                        if(.not. qvdw_flag) then
                                !use same lj params for all states
                                aljo  = iaclib(iaci)%avdw(iljo)+ao(iljo)
                                bljo  = iaclib(iaci)%bvdw(iljo)*bo(iljo)
                                aljh  = iaclib(iaci)%avdw(iljh)+ah(iljh)
                                bljh  = iaclib(iaci)%bvdw(iljh)*bh(iljh) 
                                aljo = aljo * aljo
                                aljo = aljo * aljo * aljo
                                aljh  = aljh  * aljh 
                                aljh  = aljh  * aljh  * aljh 
                                v_ao = bljo*aljo*aljo*r6o*r6o 
                                v_bo = 2.0*bljo*aljo*r6o
                                v_ah1= bljh*aljh*aljh*r6h1*r6h1
                                v_bh1= 2.0*bljh*aljh*r6h1
                                v_ah2= bljh*aljh*aljh*r6h2*r6h2
                                v_bh2= 2.0*bljh*aljh*r6h2
                        end if
                        do istate = 1, nstates ! for every state:
                        
                                ! set new lj params if q-atom types are used
                                if (qvdw_flag) then
                                        r6o = 1._8/r6o_hc
                                        r6o = r6o + sc_lookup(iq,iac_ow,istate)   !softcore
                                        r6o = 1._8/r6o
                                        r6h1 = 1._8/r6h1_hc
                                        r6h1 = r6h1 + sc_lookup(iq,iac_hw,istate)   !softcore
                                        r6h1 = 1._8/r6h1
                                        r6h2 = 1._8/r6h2_hc
                                        r6h2 = r6h2 + sc_lookup(iq,iac_hw,istate)   !softcore
                                        r6h2 = 1._8/r6h2
                                        aljo  = qavdw(qiac(iq,istate),1)+ao(iljo)
                                        bljo  = qbvdw(qiac(iq,istate),1)*bo(iljo)
                                        aljh  = qavdw(qiac(iq,istate),1)+ah(iljh)
                                        bljh  = qbvdw(qiac(iq,istate),1)*bh(iljh)
                                        aljo = aljo * aljo
                                        aljo = aljo * aljo * aljo
                                        aljh  = aljh  * aljh 
                                        aljh  = aljh  * aljh  * aljh 
                                        v_ao = bljo*aljo*aljo*r6o*r6o 
                                        v_bo = 2.0*bljo*aljo*r6o
                                        v_ah1= bljh*aljh*aljh*r6h1*r6h1
                                        v_bh1= 2.0*bljh*aljh*r6h1
                                        v_ah2= bljh*aljh*aljh*r6h2*r6h2
                                        v_bh2= 2.0*bljh*aljh*r6h2
                                end if


                                ! calculate qi, vel, v_a, v_b and dv
                                velo = crg_ow*qcrg(iq,istate)*ro
                                velh1 = crg_hw*qcrg(iq,istate)*rh1
                                velh2 = crg_hw*qcrg(iq,istate)*rh2

                                dvo  = dvo  + r2o *( -velo  -(12.*v_ao  -6.*v_bo )*r6o/r6o_hc)*eq(istate)%lambda    !r6o/r6o_hc softcore
                                dvh1 = dvh1 + r2h1*( -velh1 -(12.*v_ah1 -6.*v_bh1)*r6h1/r6h1_hc)*eq(istate)%lambda  !r6h1/r6h1_hc softcore
                                dvh2 = dvh2 + r2h2*( -velh2 -(12.*v_ah2 -6.*v_bh2)*r6h2/r6h2_hc)*eq(istate)%lambda  !r6h2/r6h2_hc softcore
                                ! update q-water energies
                                eq(istate)%qw%el  = eq(istate)%qw%el + velo + velh1 + velh2
                                eq(istate)%qw%vdw = eq(istate)%qw%vdw + v_ao + v_ah1 + v_ah2 - v_bo - v_bh1 - v_bh2 
                        end do !istate
                        
                        ! update forces on q-atom
                        d(3*i-2) = d(3*i-2) - dvo*dxo - dvh1*dxh1 - dvh2*dxh2
                        d(3*i-1) = d(3*i-1) - dvo*dyo - dvh1*dyh1 - dvh2*dyh2
                        d(3*i  ) = d(3*i  ) - dvo*dzo - dvh1*dzh1 - dvh2*dzh2

                        ! update forces on water
                        d(3*j-2) = d(3*j-2) + dvo*dxo
                        d(3*j-1) = d(3*j-1) + dvo*dyo
                        d(3*j  ) = d(3*j  ) + dvo*dzo
                        d(3*j+1) = d(3*j+1) + dvh1*dxh1
                        d(3*j+2) = d(3*j+2) + dvh1*dyh1
                        d(3*j+3) = d(3*j+3) + dvh1*dzh1
                        d(3*j+4) = d(3*j+4) + dvh2*dxh2
                        d(3*j+5) = d(3*j+5) + dvh2*dyh2
                        d(3*j+6) = d(3*j+6) + dvh2*dzh2

                end do !iq
        end do !jw
end subroutine nonbon2_qw_box

!----------------------------------------------------------------------------------------------------

subroutine nonbon2_ww
! local variables
integer                                         :: iw,ip,i,j,i3,j3,ia
integer                                         :: iaci,iacj,ilj,ja
real(8)                                         :: alj,blj
integer                                         :: ipstart
real(8)                                         :: dx1,dx2,dx3,r2,r,r6,r12
real(8)                                         :: vel,v_a,v_b,dv

! global variables used:
!  nat_solute, iac, crg, ljcod, iaclib, x, d, e

ipstart = 1

do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
! for every assigned water molecule:

do ia = 1, 3
  ! for every atom of the current water molecule:
  i    = nat_solute+3*(iw-1)+ia
  i3   = i*3-3
  iaci = iac(i)
  crg(i)   = crg(i)

  ip = ipstart
  do while (nbww(ip) .ne. 0)
        ! loop over the interactions with other water molecules

        ! x-o
        j    = nbww(ip)
        j3   = j*3-3
        iacj = iac(j)
        ilj  = ljcod(iac(i),iac(j))
        crg(j)   = crg(j)
        alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        alj  = alj*alj
        alj  = alj*alj*alj
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6
        vel  = crg(i)*crg(j)*r
        v_a  = blj*alj*alj*r12 
        v_b  = 2.0*blj*alj*r6
        dv   = r2*( -vel -12.*v_a +6.*v_b )
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3
        e%ww%el  = e%ww%el + vel       
        e%ww%vdw = e%ww%vdw + v_a - v_b

        ! x-h1
        j    = j + 1
        j3   = j3 + 3
        iacj = iac(j)
        ilj  = ljcod(iac(i),iac(j))
        crg(j)   = crg(j)
        alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        alj  = alj*alj
        alj  = alj*alj*alj
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6
        vel  = crg(i)*crg(j)*r
        v_a  = blj*alj*alj*r12 
        v_b  = 2.0*blj*alj*r6
        dv   = r2*( -vel -12.*v_a +6.*v_b )
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3
        e%ww%el  = e%ww%el + vel       
        e%ww%vdw = e%ww%vdw + v_a - v_b

        ! x-h2
        j    = j + 1
        j3   = j3 + 3
        iacj = iac(j)
        ilj  = ljcod(iac(i),iac(j))
        crg(j)   = crg(j)
        alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        alj  = alj*alj
        alj  = alj*alj*alj
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 ) 
        r6   = r2*r2*r2
        r12  = r6*r6
        vel  = crg(i)*crg(j)*r
        v_a  = blj*alj*alj*r12 
        v_b  = 2.0*blj*alj*r6
        dv   = r2*( -vel -12.*v_a +6.*v_b )
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3
        e%ww%el  = e%ww%el + vel       
        e%ww%vdw = e%ww%vdw + v_a - v_b

        ip = ip + 1
  end do
end do

ipstart = ip +1                                         ! skip over the 0
end do

end subroutine nonbon2_ww

!----------------------------------------------------------------------------------------------------
subroutine nonbon2_ww_box
  ! local variables
  integer                                               :: iw,ip,i,j,i3,j3,ia
  integer                                               :: iaci,iacj,ilj,ja
  real(8)                                               :: alj,blj
  integer                                               :: ipstart
  real(8)                                               :: dx1,dx2,dx3,r2,r,r6,r12
  real(8)                                               :: vel,v_a,v_b,dv
  real(8)                                               :: ds1, ds2, ds3

  ! global variables used:
  !  nat_solute, iac, crg, ljcod, iaclib, x, d, e

  ipstart = 1

  do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
        ! for every assigned water molecule:

        do ia = 1, 3
          ! for every atom of the current water molecule:
          i    = nat_solute+3*(iw-1)+ia
          i3   = i*3-3
          iaci = iac(i)
          crg(i)   = crg(i)

          ip = ipstart
          do while (nbww(ip) .ne. 0)
                ! loop over the interactions with other water molecules

                ! x-o
                j    = nbww(ip)
                j3   = j*3-3
                iacj = iac(j)
                ilj  = ljcod(iac(i),iac(j))
                crg(j)   = crg(j)
                alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
                blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
                alj  = alj*alj
                alj  = alj*alj*alj

                !distance between this oxygen atom and the oxygen atom of the above watermolecule, iw
                ds1 = x(j3+1) - x( 3*(nat_solute+3*(iw-1)+1) - 2 )
                ds2 = x(j3+2) - x( 3*(nat_solute+3*(iw-1)+1) - 1 )
                ds3 = x(j3+3) - x( 3*(nat_solute+3*(iw-1)+1)     )
                !the peridic shift
                ds1 = boxlength(1)*nint( ds1*inv_boxl(1) )
                ds2 = boxlength(2)*nint( ds2*inv_boxl(2) )
                ds3 = boxlength(3)*nint( ds3*inv_boxl(3) )


                dx1  = x(j3+1) - x(i3+1)
                dx2  = x(j3+2) - x(i3+2)
                dx3  = x(j3+3) - x(i3+3)
                dx1 = dx1 - ds1
                dx2 = dx2 - ds2
                dx3 = dx3 - ds3

                r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
                r2   = 1./r2
                r    = sqrt ( r2 ) 
                r6   = r2*r2*r2
                r12  = r6*r6
                vel  = crg(i)*crg(j)*r
                v_a  = blj*alj*alj*r12 
                v_b  = 2.0*blj*alj*r6
                dv   = r2*( -vel -12.*v_a +6.*v_b )
                d(i3+1) = d(i3+1) - dv*dx1
                d(i3+2) = d(i3+2) - dv*dx2
                d(i3+3) = d(i3+3) - dv*dx3
                d(j3+1) = d(j3+1) + dv*dx1
                d(j3+2) = d(j3+2) + dv*dx2
                d(j3+3) = d(j3+3) + dv*dx3
                e%ww%el  = e%ww%el + vel       
                e%ww%vdw = e%ww%vdw + v_a - v_b

                ! x-h1
                j    = j + 1
                j3   = j3 + 3
                iacj = iac(j)
                ilj  = ljcod(iac(i),iac(j))
                crg(j)   = crg(j)
                alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
                blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
                alj  = alj*alj
                alj  = alj*alj*alj
                dx1  = x(j3+1) - x(i3+1)
                dx2  = x(j3+2) - x(i3+2)
                dx3  = x(j3+3) - x(i3+3)
                dx1 = dx1 - ds1
                dx2 = dx2 - ds2
                dx3 = dx3 - ds3
                
                r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
                r2   = 1./r2
                r    = sqrt ( r2 ) 
                r6   = r2*r2*r2
                r12  = r6*r6
                vel  = crg(i)*crg(j)*r
                v_a  = blj*alj*alj*r12 
                v_b  = 2.0*blj*alj*r6
                dv   = r2*( -vel -12.*v_a +6.*v_b )
                d(i3+1) = d(i3+1) - dv*dx1
                d(i3+2) = d(i3+2) - dv*dx2
                d(i3+3) = d(i3+3) - dv*dx3
                d(j3+1) = d(j3+1) + dv*dx1
                d(j3+2) = d(j3+2) + dv*dx2
                d(j3+3) = d(j3+3) + dv*dx3
                e%ww%el  = e%ww%el + vel       
                e%ww%vdw = e%ww%vdw + v_a - v_b
                
                ! x-h2
                j    = j + 1
                j3   = j3 + 3
                iacj = iac(j)
                ilj  = ljcod(iac(i),iac(j))
                crg(j)   = crg(j)
                alj  = iaclib(iaci)%avdw(ilj)+iaclib(iacj)%avdw(ilj)
                blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
                alj  = alj*alj
                alj  = alj*alj*alj
                dx1  = x(j3+1) - x(i3+1)
                dx2  = x(j3+2) - x(i3+2)
                dx3  = x(j3+3) - x(i3+3)
                dx1 = dx1 - ds1
                dx2 = dx2 - ds2
                dx3 = dx3 - ds3
        
                r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
                r2   = 1./r2
                r    = sqrt ( r2 ) 
                r6   = r2*r2*r2
                r12  = r6*r6
                vel  = crg(i)*crg(j)*r
                v_a  = blj*alj*alj*r12 
                v_b  = 2.0*blj*alj*r6
                dv   = r2*( -vel -12.*v_a +6.*v_b )
                d(i3+1) = d(i3+1) - dv*dx1
                d(i3+2) = d(i3+2) - dv*dx2
                d(i3+3) = d(i3+3) - dv*dx3
                d(j3+1) = d(j3+1) + dv*dx1
                d(j3+2) = d(j3+2) + dv*dx2
                d(j3+3) = d(j3+3) + dv*dx3
                e%ww%el  = e%ww%el + vel       
                e%ww%vdw = e%ww%vdw + v_a - v_b

                ip = ip + 1
          end do
        end do

        ipstart = ip +1                                         ! skip over the 0
  end do

end subroutine nonbon2_ww_box

!-----------------------------------------------------------------------

subroutine nonbond_pp
! local variables
integer                                         :: ip
type(nb_type), pointer          :: pa, pb
real(8)                                         :: dx1a,dx2a,dx3a,r2a,ra,r6a
real(8)                                         :: vela,v_aa,v_ba,dva
real(8)                                         :: dx1b,dx2b,dx3b,r2b,rb,r6b
real(8)                                         :: velb,v_ab,v_bb,dvb

! global variables used:

!  x, crg, el14_scale, iaclib, d, e, 


do ip = 1, nbpp_pair - 1, 2
! for every second pair (two parallel runs to improve performance):


! set up pointers
pa => nbpp(ip)
pb => nbpp(ip+1)

! calculate the distance r
dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
dx3b  = x(pb%j*3-0) - x(pb%i*3-0)
r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
ra    = sqrt ( r2a ) 
rb    = sqrt ( r2b ) 
r6a   = r2a*r2a*r2a
r6b   = r2b*r2b*r2b

! calculate vel and dv
vela  = crg(pa%i)*crg(pa%j)*ra
velb  = crg(pb%i)*crg(pb%j)*rb
if ( pa%ljcod .eq. 3 ) then 
  vela = vela*el14_scale
end if
if ( pb%ljcod .eq. 3 ) then
 velb = velb*el14_scale
end if
v_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%ljcod) &
        *iaclib(iac(pa%j))%avdw(pa%ljcod)
v_ab = r6b*r6b*iaclib(iac(pb%i))%avdw(pb%ljcod) &
        *iaclib(iac(pb%j))%avdw(pb%ljcod)
v_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%ljcod) &
        *iaclib(iac(pa%j))%bvdw(pa%ljcod)
v_bb = r6b*iaclib(iac(pb%i))%bvdw(pb%ljcod) &
        *iaclib(iac(pb%j))%bvdw(pb%ljcod)

dva   = r2a*( -vela -12.*v_aa +6.*v_ba )
dvb   = r2b*( -velb -12.*v_ab +6.*v_bb )

! update d
d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

! update energies
e%pp%el  = e%pp%el + vela + velb
e%pp%vdw = e%pp%vdw + v_aa + v_ab - v_ba - v_bb
end do

if (ip .eq. nbpp_pair) then
! odd #pairs, handle the last pair
pa => nbpp(ip)

dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
ra   = sqrt(r2a)
r6a   = r2a*r2a*r2a

vela  = crg(pa%i)*crg(pa%j)*ra
if ( pa%ljcod .eq. 3 ) vela = vela*el14_scale
v_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%ljcod)*iaclib(iac(pa%j))%avdw(pa%ljcod)
v_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%ljcod)*iaclib(iac(pa%j))%bvdw(pa%ljcod)
dva   = r2a*( -vela -12.*v_aa +6.*v_ba )

d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

e%pp%el  = e%pp%el + vela 
e%pp%vdw = e%pp%vdw + v_aa - v_ba 
end if

end subroutine nonbond_pp

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_pp_box
  ! local variables
  integer                                               :: ip, group, ga, gb
  type(nb_type), pointer                :: pa, pb
  real(8)                                               :: dx1a,dx2a,dx3a,r2a,ra,r6a
  real(8)                                               :: vela,v_aa,v_ba,dva
  real(8)                                               :: dx1b,dx2b,dx3b,r2b,rb,r6b
  real(8)                                               :: velb,v_ab,v_bb,dvb

  ! global variables used:
  !  x, crg, el14_scale, iaclib, d, e, 

  !compute the peridocal shift for every charge group pair
        do group = 1, nbpp_cgp_pair
                ga = nbpp_cgp(group)%i !atom index for the two switching atoms
                gb = nbpp_cgp(group)%j

                !the distance between the two switching atoms
                dx1a = x(3*gb-2) - x(3*ga-2)
                dx2a = x(3*gb-1) - x(3*ga-1)
                dx3a = x(3*gb  ) - x(3*ga  )

                nbpp_cgp(group)%x = boxlength(1)*nint( dx1a*inv_boxl(1) )
                nbpp_cgp(group)%y = boxlength(2)*nint( dx2a*inv_boxl(2) )       
                nbpp_cgp(group)%z = boxlength(3)*nint( dx3a*inv_boxl(3) )

        end do

   do ip = 1, nbpp_pair - 1, 2
    ! for every second pair (two parallel runs to improve performance):

        ! set up pointers
        pa => nbpp(ip)
        pb => nbpp(ip+1)
        ga = pa%cgp_pair
        gb = pb%cgp_pair

        ! calculate the distance r
    dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
    dx1b  = x(pb%j*3-2) - x(pb%i*3-2)
    dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
    dx2b  = x(pb%j*3-1) - x(pb%i*3-1)
    dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
    dx3b  = x(pb%j*3-0) - x(pb%i*3-0)
        dx1a = dx1a - nbpp_cgp(ga)%x         
        dx1b = dx1b - nbpp_cgp(gb)%x 
        dx2a = dx2a - nbpp_cgp(ga)%y    
        dx2b = dx2b - nbpp_cgp(gb)%y                
        dx3a = dx3a - nbpp_cgp(ga)%z  
        dx3b = dx3b - nbpp_cgp(gb)%z   
        
    r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
    r2b   = 1./(dx1b*dx1b + dx2b*dx2b + dx3b*dx3b)
    ra    = sqrt ( r2a ) 
    rb    = sqrt ( r2b ) 
    r6a   = r2a*r2a*r2a
    r6b   = r2b*r2b*r2b

        ! calculate vel and dv
    vela  = crg(pa%i)*crg(pa%j)*ra
    velb  = crg(pb%i)*crg(pb%j)*rb
    if ( pa%ljcod .eq. 3 ) then 
          vela = vela*el14_scale
        end if
        if ( pb%ljcod .eq. 3 ) then
          velb = velb*el14_scale
        end if
        v_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%ljcod) &
                *iaclib(iac(pa%j))%avdw(pa%ljcod)
        v_ab = r6b*r6b*iaclib(iac(pb%i))%avdw(pb%ljcod) &
                *iaclib(iac(pb%j))%avdw(pb%ljcod)
        v_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%ljcod) &
                *iaclib(iac(pa%j))%bvdw(pa%ljcod)
        v_bb = r6b*iaclib(iac(pb%i))%bvdw(pb%ljcod) &
                *iaclib(iac(pb%j))%bvdw(pb%ljcod)
    dva   = r2a*( -vela -12.*v_aa +6.*v_ba )
    dvb   = r2b*( -velb -12.*v_ab +6.*v_bb )

        ! update d
    d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
    d(pb%i*3-2) = d(pb%i*3-2) - dvb*dx1b
    d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
    d(pb%i*3-1) = d(pb%i*3-1) - dvb*dx2b
    d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
    d(pb%i*3-0) = d(pb%i*3-0) - dvb*dx3b
    d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
    d(pb%j*3-2) = d(pb%j*3-2) + dvb*dx1b
    d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
    d(pb%j*3-1) = d(pb%j*3-1) + dvb*dx2b
    d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a
    d(pb%j*3-0) = d(pb%j*3-0) + dvb*dx3b

        ! update energies
    e%pp%el  = e%pp%el + vela + velb
    e%pp%vdw = e%pp%vdw + v_aa + v_ab - v_ba - v_bb
  end do

  if (ip .eq. nbpp_pair) then
    ! odd #pairs, handle the last pair
        pa => nbpp(ip)
        ga = pa%cgp_pair

        dx1a  = x(pa%j*3-2) - x(pa%i*3-2)
        dx2a  = x(pa%j*3-1) - x(pa%i*3-1)
        dx3a  = x(pa%j*3-0) - x(pa%i*3-0)
        dx1a = dx1a - nbpp_cgp(ga)%x  
        dx2a = dx2a - nbpp_cgp(ga)%y  
        dx3a = dx3a - nbpp_cgp(ga)%z  
        
        r2a   = 1./(dx1a*dx1a + dx2a*dx2a + dx3a*dx3a)
        ra   = sqrt(r2a)


        r6a   = r2a*r2a*r2a

        vela  = crg(pa%i)*crg(pa%j)*ra
        if ( pa%ljcod .eq. 3 ) vela = vela*el14_scale
        v_aa = r6a*r6a*iaclib(iac(pa%i))%avdw(pa%ljcod)*iaclib(iac(pa%j))%avdw(pa%ljcod)
        v_ba = r6a*iaclib(iac(pa%i))%bvdw(pa%ljcod)*iaclib(iac(pa%j))%bvdw(pa%ljcod)
        dva   = r2a*( -vela -12.*v_aa +6.*v_ba )

        d(pa%i*3-2) = d(pa%i*3-2) - dva*dx1a
        d(pa%i*3-1) = d(pa%i*3-1) - dva*dx2a
        d(pa%i*3-0) = d(pa%i*3-0) - dva*dx3a
        d(pa%j*3-2) = d(pa%j*3-2) + dva*dx1a
        d(pa%j*3-1) = d(pa%j*3-1) + dva*dx2a
        d(pa%j*3-0) = d(pa%j*3-0) + dva*dx3a

        e%pp%el  = e%pp%el + vela 
        e%pp%vdw = e%pp%vdw + v_aa - v_ba 
  end if

end subroutine nonbond_pp_box

!-----------------------------------------------------------------------

subroutine nonbond_pw
! local variables
integer                                         :: ip,i,j,i3,j3,iaci,iacj,ilj
real(8)                                         :: alj,blj,dx1,dx2,dx3,r2,r,r6,r12
real(8)                                         :: vel,v_a,v_b,dv

! global variables used:
!  iac, crg, iaclib, x, d, e

do ip = 1, nbpw_pair
! for every assigned pair:

i    = nbpw(ip)%i
j    = nbpw(ip)%j
i3   = i*3-3
j3   = j*3-3
iaci = iac(i)
iacj = iac(j)

ilj  = nbpw(ip)%ljcod
alj  = iaclib(iaci)%avdw(ilj)*iaclib(iacj)%avdw(ilj)
blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)

! calculate dx and r
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
r2   = 1./r2
r    = sqrt ( r2 )
r6   = r2*r2*r2
r12  = r6*r6

! calculate vel and dv
vel  = crg(i)*crg(j)*r
v_a  = alj*r12 
v_b  = blj*r6
dv   = r2*( -vel -12.*v_a +6.*v_b )

! update forces
d(i3+1) = d(i3+1) - dv*dx1
d(i3+2) = d(i3+2) - dv*dx2
d(i3+3) = d(i3+3) - dv*dx3
d(j3+1) = d(j3+1) + dv*dx1
d(j3+2) = d(j3+2) + dv*dx2
d(j3+3) = d(j3+3) + dv*dx3

! update energies
e%pw%el  = e%pw%el + vel       
e%pw%vdw = e%pw%vdw + v_a - v_b
end do
end subroutine nonbond_pw

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_pw_box
  ! local variables
  integer                                               :: ip,i,j,i3,j3,iaci,iacj,ilj
  real(8)                                               :: alj,blj,dx1,dx2,dx3,r2,r,r6,r12
  real(8)                                               :: vel,v_a,v_b,dv
  integer                                               :: group, ga, gb

  ! global variables used:
  !  iac, crg, iaclib, x, d, e

  !compute the peridocal shift for every charge group pair
        do group = 1, nbpw_cgp_pair
                ga = nbpw_cgp(group)%i !atom index for the solute switching atoms
                gb = nbpw_cgp(group)%j !atom index for the solvent switching atom 

                !the distance between the two switching atoms
                dx1 = x(3*gb-2) - x(3*ga-2)
                dx2 = x(3*gb-1) - x(3*ga-1)
                dx3 = x(3*gb  ) - x(3*ga  )

                nbpw_cgp(group)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
                nbpw_cgp(group)%y = boxlength(2)*nint( dx2*inv_boxl(2) )        
                nbpw_cgp(group)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

        end do


  do ip = 1, nbpw_pair
        ! for every assigned pair:

        i    = nbpw(ip)%i  ! solute atom 
        j    = nbpw(ip)%j       ! solvent atom 
        group = nbpw(ip)%cgp_pair 
        i3   = i*3-3
        j3   = j*3-3
        iaci = iac(i)
        iacj = iac(j)
        ilj  = nbpw(ip)%ljcod
        alj  = iaclib(iaci)%avdw(ilj)*iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)

        ! calculate dx and r
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        dx1 = dx1 - nbpw_cgp(group)%x  
        dx2 = dx2 - nbpw_cgp(group)%y  
        dx3 = dx3 - nbpw_cgp(group)%z 
        
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r    = sqrt ( r2 )
        r6   = r2*r2*r2
        r12  = r6*r6

        ! calculate vel and dv
        vel  = crg(i)*crg(j)*r
        v_a  = alj*r12 
        v_b  = blj*r6
        dv   = r2*( -vel -12.*v_a +6.*v_b )

        ! update forces
        d(i3+1) = d(i3+1) - dv*dx1
        d(i3+2) = d(i3+2) - dv*dx2
        d(i3+3) = d(i3+3) - dv*dx3
        d(j3+1) = d(j3+1) + dv*dx1
        d(j3+2) = d(j3+2) + dv*dx2
        d(j3+3) = d(j3+3) + dv*dx3

        ! update energies
        e%pw%el  = e%pw%el + vel       
        e%pw%vdw = e%pw%vdw + v_a - v_b
  end do
end subroutine nonbond_pw_box

!-----------------------------------------------------------------------

subroutine nonbond_qq
! local variables
integer                                 :: istate, is
integer                                 :: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,ilj
real(8)                                 :: qi,qj,alj,blj,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)                                 :: vel,v_a,v_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq   !q-atom number
  i    = iqseq(iq)           !atom number
  j    = nbqq(ip,istate)%j   !atom number
  jq   = nbqq(ip,istate)%jq  !q-atom number (if any)
  i3   = i*3-3
  j3   = j*3-3
  ilj  = nbqq(ip,istate)%ljcod
  qi   = qcrg(iq,istate)
  el_scale = nbqq(ip,istate)%el_scale

  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        alj  = iaclib(iaci)%avdw(ilj)*iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)

        if (jq /= 0) then
          qj = qcrg(jq,istate)
        else
          qj = crg(j)
        end if
  else
        iaci = qiac(iq,istate)
        alj  = qavdw(iaci,ilj)
        blj  = qbvdw(iaci,ilj)
        if (jq /= 0) then
          iacj = qiac(jq,istate)
          alj = alj*qavdw(iacj,ilj)
          blj = blj*qbvdw(iacj,ilj)
          qj = qcrg(jq,istate)
        else
          if ( ilj .eq. 2 ) then
                alj  = qavdw(iaci,1)
                blj  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          alj = alj*iaclib(iacj)%avdw(ilj)
          blj = blj*iaclib(iacj)%bvdw(ilj)
          qj = crg(j)
        end if


  end if




  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)
  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

  r6_hc   = r2*r2*r2   !hardcore
  r6   = r2*r2*r2+sc_lookup(iq,natyps+jq,istate)   !use softcore instead. sc is 0 for hardcore mpa
  r6   = 1./r6
  r12  = r6*r6

  r2   = 1./r2
  r    = sqrt ( r2 ) 

  ! calculate vel, v_a, v_b and dv
  vel  = qi*qj*r*el_scale
  if ( ilj .eq. 3 ) vel = vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. ilj .eq. 2 ) then
        v_a = alj*exp(-blj/r)
        v_b = 0.0
        dv  = r2*( -vel -blj*v_a/r )*eq(istate)%lambda
  else
        v_a = alj*r12 
        v_b = blj*r6
        dv  = r2*( -vel - (12.*v_a - 6.*v_b)*r6_hc*r6 )*eq(istate)%lambda
  endif

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        eq(istate)%qq%el  = eq(istate)%qq%el + vel
        eq(istate)%qq%vdw = eq(istate)%qq%vdw + v_a - v_b
  else
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
  end if
end do
end do
end subroutine nonbond_qq

!-----------------------------------------------------------------------

subroutine nonbond_qq_lib_charges
!special version that uses library charges - for transformation of
!solute-surrounding only

! local variables
integer                                 :: istate
integer                                 :: ip,iq,jq,i,j,k,i3,j3,iaci,iacj,ilj
real(8)                                 :: qi,qj,alj,blj,dx1,dx2,dx3,r2,r,r6,r12,r6_hc
real(8)                                 :: vel,v_a,v_b,dv,el_scale

do istate = 1, nstates
! for every state:

do ip = 1, nbqq_pair(istate)
  ! for every pair:

  iq   = nbqq(ip,istate)%iq
  i    = iqseq(iq)
  j    = nbqq(ip,istate)%j
  jq   = nbqq(ip,istate)%jq
  i3   = i*3-3
  j3   = j*3-3
  ilj  = nbqq(ip,istate)%ljcod

  if (jq /= 0) then
          qi   = crg(i)   !use library charges for i
  else
          qi   = qcrg(iq,istate)   !since j is not a qatom we need the fep file charges for i
  end if

  el_scale = nbqq(ip,istate)%el_scale

  if (.not. qvdw_flag) then
        iaci = iac(i)
        iacj = iac(j)
        alj  = iaclib(iaci)%avdw(ilj)*iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)
        qj = crg(j)
  else
        iaci = qiac(iq,istate)
        alj  = qavdw(iaci,ilj)
        blj  = qbvdw(iaci,ilj)
        if (jq /= 0) then
          iacj = qiac(jq,istate)
          alj = alj*qavdw(iacj,ilj)
          blj = blj*qbvdw(iacj,ilj)
        else
          if ( ilj .eq. 2 ) then
                alj  = qavdw(iaci,1)
                blj  = qbvdw(iaci,1)
          end if
          iacj = iac(j)
          alj = alj*iaclib(iacj)%avdw(ilj)
          blj = blj*iaclib(iacj)%bvdw(ilj)
        end if
    qj = crg(j)
  end if

  ! calculate dx and r
  dx1  = x(j3+1) - x(i3+1)
  dx2  = x(j3+2) - x(i3+2)
  dx3  = x(j3+3) - x(i3+3)
  r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

  r6_hc = r2*r2*r2   !hardcore
  r6   = r2*r2*r2+sc_lookup(iq,natyps+jq,istate)   !sc_lookup is softcore fix mpa
  r6   = 1./r6
  r12  = r6*r6

  r2   = 1./r2
  r    = sqrt ( r2 ) 

  ! calculate vel, v_a, v_b and dv
  vel  = qi*qj*r*el_scale
  if ( ilj .eq. 3 ) vel = vel*el14_scale
  if (qvdw_flag .and. jq /= 0 .and. ilj .eq. 2 ) then
        v_a = alj*exp(-blj/r)
        v_b = 0.0
        dv  = r2*( -vel -blj*v_a/r )*eq(istate)%lambda
  else
        v_a = alj*r12 
        v_b = blj*r6
        dv  = r2*( -vel -(12.*v_a -6.*v_b)*r6_hc*r6 )*eq(istate)%lambda
  endif

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update energies
  if ( jq /= 0 ) then
        eq(istate)%qq%el  = eq(istate)%qq%el + vel
        eq(istate)%qq%vdw = eq(istate)%qq%vdw + v_a - v_b
  else  ! j is not a qatom
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
  end if
end do
end do
end subroutine nonbond_qq_lib_charges


!-----------------------------------------------------------------------

subroutine nonbond_qp

!calculate non-bonded interactions between q-atom i and  non-q-atoms j
!using standard atom types

! local variables
integer                                         :: ip,iq,i,j,i3,j3,iaci,iacj,ilj
integer                                         :: istate
real(8)                                         :: alj,blj,dx1,dx2,dx3,r2,r,r6
real(8)                                         :: vel,v_a,v_b,dv

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, 
! qcrg, el14_scale, eq, d, nat_solute

do ip = 1, nbqp_pair
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip)%i
i    = iqseq(iq)
j    = nbqp(ip)%j
i3   = i*3-3
j3   = j*3-3
iacj = iac(j)
ilj  = nbqp(ip)%ljcod
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
!softcore not needed here since nonbond_qp_qvdw is called instead
r2   = 1./r2
r6   = r2*r2*r2
r    = sqrt(r2)
do istate = 1, nstates
  ! for every state:

  ! calculate iaci, alj and blj
        iaci = iac(i)
        alj  = iaclib(iaci)%avdw(ilj)*iaclib(iacj)%avdw(ilj)
        blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)

  ! calculate qi, vel, v_a, v_b and dv
  vel  = qcrg(iq,istate)*crg(j)*r
  if ( ilj .eq. 3 ) vel = vel*el14_scale


  v_a  = alj*r6*r6 
  v_b  = blj*r6        
  dv   = r2*( -vel -12.*v_a +6.*v_b )*eq(istate)%lambda

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein energies
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
end do ! istate

end do
end subroutine nonbond_qp

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_qp_box

        !calculate non-bonded interactions between q-atom i and  non-q-atoms j
        !using standard atom types

  ! local variables
  integer                                               :: ip,iq,i,j,i3,j3,iaci,iacj,ilj
  integer                                               :: istate
  real(8)                                               :: alj,blj,dx1,dx2,dx3,r2,r,r6
  real(8)                                               :: vel,v_a,v_b,dv
  integer                                               :: group, gr, ia
  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, 
  ! qcrg, el14_scale, eq, d, nat_solute


        !compute the peridocal shift for every charge group pair
        do gr = 1, nbqp_cgp_pair
                ia = nbqp_cgp(gr)%i !atom index for the switching atom
                
                !the distance between the two switching atoms
                dx1 = x(3*ia-2) - x(3*qswitch-2)
                dx2 = x(3*ia-1) - x(3*qswitch-1)
                dx3 = x(3*ia  ) - x(3*qswitch  )

                nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
                nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )   
                nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

        end do


  do ip = 1, nbqp_pair
        ! for every assigned q-s pair:
        ! init state-invariant variables:
        iq   = nbqp(ip)%i
        i    = iqseq(iq)
        j    = nbqp(ip)%j
        i3   = i*3-3
        j3   = j*3-3
        iacj = iac(j)
        ilj  = nbqp(ip)%ljcod
        group = nbqp(ip)%cgp_pair
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        dx1 = dx1 - nbqp_cgp(group)%x  
        dx2 = dx2 - nbqp_cgp(group)%y  
        dx3 = dx3 - nbqp_cgp(group)%z  
        
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3
        r2   = 1./r2
        r6   = r2*r2*r2   !softcore not needed here. taken care of in nonbond_qp_qvdw_box
        r    = sqrt(r2)

        do istate = 1, nstates
          ! for every state:

          ! calculate iaci, alj and blj
                iaci = iac(i)
                alj  = iaclib(iaci)%avdw(ilj)*iaclib(iacj)%avdw(ilj)
                blj  = iaclib(iaci)%bvdw(ilj)*iaclib(iacj)%bvdw(ilj)

          ! calculate qi, vel, v_a, v_b and dv
          vel  = qcrg(iq,istate)*crg(j)*r
          if ( ilj .eq. 3 ) vel = vel*el14_scale

          v_a  = alj*r6*r6 
          v_b  = blj*r6
          dv   = r2*( -vel -12.*v_a +6.*v_b )*eq(istate)%lambda

          ! update forces
          d(i3+1) = d(i3+1) - dv*dx1
          d(i3+2) = d(i3+2) - dv*dx2
          d(i3+3) = d(i3+3) - dv*dx3
          d(j3+1) = d(j3+1) + dv*dx1
          d(j3+2) = d(j3+2) + dv*dx2
          d(j3+3) = d(j3+3) + dv*dx3

          ! update q-protein energies
                eq(istate)%qp%el  = eq(istate)%qp%el + vel
                eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
        end do ! istate

  end do
end subroutine nonbond_qp_box

!----------------------------------------------------------------------

subroutine nonbond_qp_qvdw
!calculate nonbonded interactions between q-atom i and non-q-atom j
!using q-atom types

! local variables
integer                                         :: ip,iq,i,j,i3,j3,iaci,iacj,ilj,qlj
integer                                         :: istate
real(8)                                         :: alj,blj,dx1,dx2,dx3,r2,r,r6
real(8)                                         :: vel,v_a,v_b,dv,r6_sc

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute


do ip = 1, nbqp_pair
! for every assigned q-s pair:

! init state-invariant variables:
iq   = nbqp(ip)%i
i    = iqseq(iq)
j    = nbqp(ip)%j
i3   = i*3-3
j3   = j*3-3
iacj = iac(j)
ilj  = nbqp(ip)%ljcod
qlj  = nbqp(ip)%qljcod
dx1  = x(j3+1) - x(i3+1)
dx2  = x(j3+2) - x(i3+2)
dx3  = x(j3+3) - x(i3+3)
r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

r6   = r2*r2*r2
r2   = 1./r2
r    = sqrt(r2)

do istate = 1, nstates
        ! for every state:

        ! calculate iaci, alj and blj
        iaci = qiac(iq,istate)
        alj  = qavdw(iaci,qlj)*iaclib(iacj)%avdw(ilj)
        blj  = qbvdw(iaci,qlj)*iaclib(iacj)%bvdw(ilj)

  ! calculate qi, vel, v_a, v_b and dv
  vel  = qcrg(iq,istate)*crg(j)*r
  if ( ilj .eq. 3 ) vel = vel*el14_scale

  r6_sc = r6 + sc_lookup(iq,iacj,istate) !sc_lookup is softcore fix mpa
  v_a  = alj/(r6_sc*r6_sc)
  v_b  = blj/(r6_sc)
  dv   = r2*( -vel -(12.*v_a -6.*v_b)*(r6/r6_sc) )*eq(istate)%lambda  !r6 is r^6 not 1/r^6, r6_sc is r^6+sc not 1/(r^6+sc)

  ! update forces
  d(i3+1) = d(i3+1) - dv*dx1
  d(i3+2) = d(i3+2) - dv*dx2
  d(i3+3) = d(i3+3) - dv*dx3
  d(j3+1) = d(j3+1) + dv*dx1
  d(j3+2) = d(j3+2) + dv*dx2
  d(j3+3) = d(j3+3) + dv*dx3

  ! update q-protein energies
        eq(istate)%qp%el  = eq(istate)%qp%el + vel
        eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
end do ! istate

end do
end subroutine nonbond_qp_qvdw
!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_qp_qvdw_box
        !calculate nonbonded interactions between q-atom i and non-q-atom j
        !using q-atom types

  ! local variables
  integer                                               :: ip,iq,i,j,i3,j3,iaci,iacj,ilj,qlj
  integer                                               :: istate
  real(8)                                               :: alj,blj,dx1,dx2,dx3,r2,r,r6
  real(8)                                               :: vel,v_a,v_b,dv,r6_sc
  integer                                               :: group, gr, ia


  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute


        !compute the peridocal shift for every charge group pair
        do gr = 1, nbqp_cgp_pair
                ia = nbqp_cgp(gr)%i !atom index for the atom
                
                !the distance between the two switching atoms
                dx1 = x(3*ia-2) - x(3*qswitch-2)
                dx2 = x(3*ia-1) - x(3*qswitch-1)
                dx3 = x(3*ia  ) - x(3*qswitch  )

                nbqp_cgp(gr)%x = boxlength(1)*nint( dx1*inv_boxl(1) )
                nbqp_cgp(gr)%y = boxlength(2)*nint( dx2*inv_boxl(2) )   
                nbqp_cgp(gr)%z = boxlength(3)*nint( dx3*inv_boxl(3) )

        end do


  do ip = 1, nbqp_pair
        ! for every assigned q-s pair:

        ! init state-invariant variables:
        iq   = nbqp(ip)%i
        i    = iqseq(iq)
        j    = nbqp(ip)%j
        i3   = i*3-3
        j3   = j*3-3
        iacj = iac(j)
        ilj  = nbqp(ip)%ljcod
        qlj  = nbqp(ip)%qljcod
        group = nbqp(ip)%cgp_pair
        
        dx1  = x(j3+1) - x(i3+1)
        dx2  = x(j3+2) - x(i3+2)
        dx3  = x(j3+3) - x(i3+3)
        dx1 = dx1 - nbqp_cgp(group)%x  
        dx2 = dx2 - nbqp_cgp(group)%y  
        dx3 = dx3 - nbqp_cgp(group)%z  
        
        r2   = dx1*dx1 + dx2*dx2 + dx3*dx3

        r6   = r2*r2*r2
        r2   = 1./r2
        r    = sqrt(r2)

        do istate = 1, nstates
                ! for every state:

                ! calculate iaci, alj and blj
                iaci = qiac(iq,istate)
                alj  = qavdw(iaci,qlj)*iaclib(iacj)%avdw(ilj)
                blj  = qbvdw(iaci,qlj)*iaclib(iacj)%bvdw(ilj)

          ! calculate qi, vel, v_a, v_b and dv
          vel  = qcrg(iq,istate)*crg(j)*r
          if ( ilj .eq. 3 ) vel = vel*el14_scale


          r6_sc = r6 + sc_lookup(iq,iacj,istate)        !sc_lookup is softcore fix mpa
          v_a  = alj/(r6_sc*r6_sc)
          v_b  = blj/r6_sc         !sc_lookup is softcore fix mpa
          dv   = r2*( -vel - ( (12.*v_a - 6.*v_b)*(r6/r6_sc) ) )*eq(istate)%lambda  !r6 is r^6 not 1/r^6, r6_sc is r^6+sc not 1/(r^6+sc)

          ! update forces
          d(i3+1) = d(i3+1) - dv*dx1
          d(i3+2) = d(i3+2) - dv*dx2
          d(i3+3) = d(i3+3) - dv*dx3
          d(j3+1) = d(j3+1) + dv*dx1
          d(j3+2) = d(j3+2) + dv*dx2
          d(j3+3) = d(j3+3) + dv*dx3

          ! update q-protein energies
                eq(istate)%qp%el  = eq(istate)%qp%el + vel
                eq(istate)%qp%vdw = eq(istate)%qp%vdw + v_a - v_b
        end do ! istate

  end do
end subroutine nonbond_qp_qvdw_box

!-----------------------------------------------------------------------

subroutine nonbond_qw_spc
!calculate non-bonded interactions between q-atoms and spc water molecules
!(optimisations rely on lj params = 0 for water h) using geometric comb. rule

! local variables
integer                                         :: jw,iq,i,j,ilj
integer                                         :: istate
real(8)                                         ::      alj, blj
real(8)                                         ::      dxo, dyo, dzo, dxh1, dyh1, dzh1, dxh2, dyh2, dzh2
real(8)                                         ::      ro, r2o, r6o, rh1, r2h1, r6h1, rh2, r2h2, r6h2
real(8)                                         ::      velo, velh1, velh2, dvo, dvh1, dvh2
real(8)                                         ::  v_a, v_b, r6o_sc, r6o_hc
real(8), save                           ::      ao(2), bo(2)
integer, save                           ::      iac_ow = 0, iac_hw = 0

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute

if(iac_ow == 0) then !set first time
        iac_ow = iac(nat_solute + 1)
        iac_hw = iac(nat_solute + 2)
        ao(1:2) = iaclib(iac(nat_solute + 1))%avdw(1:2)
        bo(1:2) = iaclib(iac(nat_solute + 1))%bvdw(1:2)
end if

!loop over listed waters
do jw = 1, nbqw_pair
        j = nbqw(jw) !top. # of o in water iw
        !loop over all q-atoms
        do iq = 1, nqat
                i = iqseq(iq)
                dxo  = x(3*j-2) - x(3*i-2)
                dyo  = x(3*j-1) - x(3*i-1)
                dzo  = x(3*j  ) - x(3*i  )
                dxh1 = x(3*j+1) - x(3*i-2)
                dyh1 = x(3*j+2) - x(3*i-1)
                dzh1 = x(3*j+3) - x(3*i  )
                dxh2 = x(3*j+4) - x(3*i-2)
                dyh2 = x(3*j+5) - x(3*i-1)
                dzh2 = x(3*j+6) - x(3*i  )
                r2o  = dxo*dxo + dyo*dyo + dzo*dzo
                rh1  = sqrt(1._8/(dxh1*dxh1 + dyh1*dyh1 + dzh1*dzh1))
                rh2  = sqrt(1._8/(dxh2*dxh2 + dyh2*dyh2 + dzh2*dzh2))
                r6o_hc  = r2o*r2o*r2o   !will set r6o to 1/r6o later, needed for softcore
                r2o  = 1._8/r2o
                ro   = sqrt(r2o)
                r2h1 = rh1*rh1
                r6h1 = r2h1*r2h1*r2h1
                r2h2 = rh2*rh2
                r6h2 = r2h2*r2h2*r2h2

                        r6o_sc = r6o_hc   !default is hardcore (i.e. not softcore)  mpa

                !reset potential
                dvo = 0
                dvh1 = 0
                dvh2 = 0
                ilj = ljcod(iac_ow, iac(i))
                if(.not. qvdw_flag) then
                        !use same lj params for all states
                                                r6o  = 1./r6o_hc     !softcore hack, see comment 15 lines up
                        alj  = iaclib(iac(i))%avdw(ilj)
                        blj  = iaclib(iac(i))%bvdw(ilj)
                        v_a  = alj*ao(ilj)*r6o*r6o 
                        v_b  = blj*bo(ilj)*r6o
                end if
                do istate = 1, nstates
                        ! for every state:
                        ! calculate iaci, alj and blj
                        if (qvdw_flag) then
                                alj  = qavdw(qiac(iq,istate),1)
                                blj  = qbvdw(qiac(iq,istate),1)
                                                r6o_sc = r6o_hc + sc_lookup(iq,iac_ow,istate)   !softcore  mpa

                                v_a  = alj*ao(ilj)/(r6o_sc*r6o_sc)
                                v_b  = blj*bo(ilj)/(r6o_sc)
                        end if
                        ! calculate qi, vel, v_a, v_b and dv
                        velo = crg_ow*qcrg(iq,istate)*ro
                        velh1 = crg_hw*qcrg(iq,istate)*rh1
                        velh2 = crg_hw*qcrg(iq,istate)*rh2
                        dvo  = dvo  + r2o*( -velo -( (12.*v_a - 6.*v_b)*(r6o_hc/r6o_sc) ))*eq(istate)%lambda
                        dvh1 = dvh1 - r2h1*velh1*eq(istate)%lambda
                        dvh2 = dvh2 - r2h2*velh2*eq(istate)%lambda
                        ! update q-water energies
                        eq(istate)%qw%el  = eq(istate)%qw%el + velo + velh1 + velh2
                        eq(istate)%qw%vdw = eq(istate)%qw%vdw + v_a - v_b
                end do !istate
                                
                                ! if qvdw_flag is true, then r6o is not the usual 1/ro^6, but rather ro^6. be careful!!! mpa
                                
                                                        
                ! update forces on q-atom
                d(3*i-2) = d(3*i-2) - dvo*dxo - dvh1*dxh1 - dvh2*dxh2
                d(3*i-1) = d(3*i-1) - dvo*dyo - dvh1*dyh1 - dvh2*dyh2
                d(3*i  ) = d(3*i  ) - dvo*dzo - dvh1*dzh1 - dvh2*dzh2
                ! update forces on water
                d(3*j-2) = d(3*j-2) + dvo*dxo
                d(3*j-1) = d(3*j-1) + dvo*dyo
                d(3*j  ) = d(3*j  ) + dvo*dzo
                d(3*j+1) = d(3*j+1) + dvh1*dxh1
                d(3*j+2) = d(3*j+2) + dvh1*dyh1
                d(3*j+3) = d(3*j+3) + dvh1*dzh1
                d(3*j+4) = d(3*j+4) + dvh2*dxh2
                d(3*j+5) = d(3*j+5) + dvh2*dyh2
                d(3*j+6) = d(3*j+6) + dvh2*dzh2

        end do !iq
end do !jw
end subroutine nonbond_qw_spc

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_qw_spc_box
        !calculate non-bonded interactions between q-atoms and spc water molecules
        !(optimisations rely on lj params = 0 for water h) using geometric comb. rule

        ! local variables
        integer                                         :: jw,iq,i,j,ilj
        integer                                         :: istate
        real(8)                                         ::      alj, blj
        real(8)                                         ::      dxo, dyo, dzo, dxh1, dyh1, dzh1, dxh2, dyh2, dzh2
        real(8)                                         ::      ro, r2o, r6o, rh1, r2h1, r6h1, rh2, r2h2, r6h2
        real(8)                                         ::      velo, velh1, velh2, dvo, dvh1, dvh2
        real(8)                                         :: v_a, v_b, r6o_sc, r6o_hc
        real(8)                                         :: dx, dy, dz, boxshiftx, boxshifty, boxshiftz
        real(8), save                           ::      ao(2), bo(2)
        integer, save                           ::      iac_ow = 0, iac_hw = 0

        ! global variables used:
        !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute

        if(iac_ow == 0) then !set first time   
                iac_ow = iac(nat_solute + 1)
                iac_hw = iac(nat_solute + 2)    ! not used !?
                ao(1:2) = iaclib(iac(nat_solute + 1))%avdw(1:2)
                bo(1:2) = iaclib(iac(nat_solute + 1))%bvdw(1:2)
        end if

        !loop over listed waters
        do jw = 1, nbqw_pair
                j = nbqw(jw) !top. # of o in water iw
        
                !compute the periodical shift
                dx = x(3*j-2) - x(3*qswitch-2)
                dy = x(3*j-1) - x(3*qswitch-1)
                dz = x(3*j  ) - x(3*qswitch  )
                boxshiftx = boxlength(1)*nint(dx*inv_boxl(1))
                boxshifty = boxlength(2)*nint(dy*inv_boxl(2))
                boxshiftz = boxlength(3)*nint(dz*inv_boxl(3))   
        
                !loop over all q-atoms
                do iq = 1, nqat
                        i = iqseq(iq)
                        dxo  = x(3*j-2) - x(3*i-2)
                        dyo  = x(3*j-1) - x(3*i-1)
                        dzo  = x(3*j  ) - x(3*i  )
                        dxh1 = x(3*j+1) - x(3*i-2)
                        dyh1 = x(3*j+2) - x(3*i-1)
                        dzh1 = x(3*j+3) - x(3*i  )
                        dxh2 = x(3*j+4) - x(3*i-2)
                        dyh2 = x(3*j+5) - x(3*i-1)
                        dzh2 = x(3*j+6) - x(3*i  )
                        dxo = dxo - boxshiftx  
                        dyo = dyo - boxshifty  
                        dzo = dzo - boxshiftz  
                        dxh1 = dxh1 - boxshiftx  
                        dyh1 = dyh1 - boxshifty  
                        dzh1 = dzh1 - boxshiftz  
                        dxh2 = dxh2 - boxshiftx  
                        dyh2 = dyh2 - boxshifty  
                        dzh2 = dzh2 - boxshiftz  

                        r2o = dxo*dxo + dyo*dyo + dzo*dzo
                        rh1 = sqrt(1._8/(dxh1*dxh1 + dyh1*dyh1 + dzh1*dzh1))
                        rh2 = sqrt(1._8/(dxh2*dxh2 + dyh2*dyh2 + dzh2*dzh2))
                        r6o_hc = r2o*r2o*r2o   !will set r6o = 1/r6o later, need for softcore mpa
                        r2o = 1._8/r2o
                        ro  = sqrt(r2o)
                        r2h1 = rh1*rh1
                        r6h1 = r2h1*r2h1*r2h1
                        r2h2 = rh2*rh2
                        r6h2 = r2h2*r2h2*r2h2

                        r6o_sc = r6o_hc   !default is hardcore (i.e. not softcore)  mpa

                        !reset potential
                        dvo = 0
                        dvh1 = 0
                        dvh2 = 0
                        ilj = ljcod(iac_ow, iac(i))
                        if(.not. qvdw_flag) then
                                !use same lj params for all states
                                alj  = iaclib(iac(i))%avdw(ilj)
                                blj  = iaclib(iac(i))%bvdw(ilj)
                                r6o  = 1._8/r6o_hc   !softcore hack  mpa
                                v_a  = alj*ao(ilj)*r6o*r6o 
                                v_b  = blj*bo(ilj)*r6o
                        end if
                        do istate = 1, nstates
                                ! for every state:

                                ! calculate iaci, alj and blj
                                if (qvdw_flag) then
                                        alj  = qavdw(qiac(iq,istate),1)
                                        blj  = qbvdw(qiac(iq,istate),1)
                                        r6o_sc = r6o_hc + sc_lookup(iq,iac_ow,istate)   !softcore  mpa

                                        v_a  = alj*ao(ilj)/(r6o_sc*r6o_sc)
                                        v_b  = blj*bo(ilj)/r6o_sc
                                end if
                                ! calculate qi, vel, v_a, v_b and dv
                                velo = crg_ow*qcrg(iq,istate)*ro
                                velh1 = crg_hw*qcrg(iq,istate)*rh1
                                velh2 = crg_hw*qcrg(iq,istate)*rh2
                                dvo  = dvo  + r2o*( -velo -( (12.*v_a - 6.*v_b)*(r6o_hc/r6o_sc) ))*eq(istate)%lambda
                                dvh1 = dvh1 - r2h1*velh1*eq(istate)%lambda
                                dvh2 = dvh2 - r2h2*velh2*eq(istate)%lambda
                                ! update q-water energies
                                eq(istate)%qw%el  = eq(istate)%qw%el + velo + velh1 + velh2
                                eq(istate)%qw%vdw = eq(istate)%qw%vdw + v_a - v_b
                        end do !istate                  
                        ! update forces on q-atom
                        d(3*i-2) = d(3*i-2) - dvo*dxo - dvh1*dxh1 - dvh2*dxh2
                        d(3*i-1) = d(3*i-1) - dvo*dyo - dvh1*dyh1 - dvh2*dyh2
                        d(3*i  ) = d(3*i  ) - dvo*dzo - dvh1*dzh1 - dvh2*dzh2
                        ! update forces on water
                        d(3*j-2) = d(3*j-2) + dvo*dxo
                        d(3*j-1) = d(3*j-1) + dvo*dyo
                        d(3*j  ) = d(3*j  ) + dvo*dzo
                        d(3*j+1) = d(3*j+1) + dvh1*dxh1
                        d(3*j+2) = d(3*j+2) + dvh1*dyh1
                        d(3*j+3) = d(3*j+3) + dvh1*dzh1
                        d(3*j+4) = d(3*j+4) + dvh2*dxh2
                        d(3*j+5) = d(3*j+5) + dvh2*dyh2
                        d(3*j+6) = d(3*j+6) + dvh2*dzh2

                end do !iq
        end do !jw
end subroutine nonbond_qw_spc_box

!-----------------------------------------------------------------------

subroutine nonbond_qw_3atom
!calculate non-bonded interactions between q-atoms and 3-atom solvent molecules
!using geometric comb. rule

! local variables
integer                                         :: jw,iq,i,j,ilj1, ilj2, ilj3, iaci
integer                                         :: istate
real(8)                                         ::      dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3
real(8)                                         ::      r_1, r2_1, r6_1, r_2, r2_2, r6_2, r_3, r2_3, r6_3
real(8)                                         ::      vel1, vel2, vel3, dv1, dv2, dv3
real(8)                                         :: v_a1,v_b1, v_a2, v_b2, v_a3, v_b3,r6_1_sc,r6_2_sc,r6_3_sc,r6_1_hc,r6_2_hc,r6_3_hc
real(8), save                           ::      a1(2), b1(2), a2(2), b2(2), a3(2), b3(2)
integer, save                           ::      iac1, iac2, iac3
real, save                                      ::      crg1, crg2, crg3
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute

if(a1(1) == 0.) then !set first time
        iac1 = iac(nat_solute+1)
        iac2 = iac(nat_solute+2)
        iac3 = iac(nat_solute+3)
        crg1 = crg(nat_solute+1)
        crg2 = crg(nat_solute+2)
        crg3 = crg(nat_solute+3)
        a1(1:2) = iaclib(iac1)%avdw(1:2)
        b1(1:2) = iaclib(iac1)%bvdw(1:2)
        a2(1:2) = iaclib(iac2)%avdw(1:2)
        b2(1:2) = iaclib(iac2)%bvdw(1:2)
        a3(1:2) = iaclib(iac3)%avdw(1:2)
        b3(1:2) = iaclib(iac3)%bvdw(1:2)
end if

!loop over listed waters
do jw = 1, nbqw_pair
        j = nbqw(jw) !top. # of o in water iw
        !loop over all q-atoms
        do iq = 1, nqat
                i = iqseq(iq)
                dx1 = x(3*j-2) - x(3*i-2)
                dy1 = x(3*j-1) - x(3*i-1)
                dz1 = x(3*j  ) - x(3*i  )
                dx2 = x(3*j+1) - x(3*i-2)
                dy2 = x(3*j+2) - x(3*i-1)
                dz2 = x(3*j+3) - x(3*i  )
                dx3 = x(3*j+4) - x(3*i-2)
                dy3 = x(3*j+5) - x(3*i-1)
                dz3 = x(3*j+6) - x(3*i  )
                r2_1 = dx1*dx1 + dy1*dy1 + dz1*dz1
                r2_2 = dx2*dx2 + dy2*dy2 + dz2*dz2
                r2_3 = dx3*dx3 + dy3*dy3 + dz3*dz3
                r6_1_hc = r2_1*r2_1*r2_1   !will set r6 to 1/r6 later, needed for softcore
                r2_1 = 1._8/r2_1
                                r_1  = sqrt(r2_1)
                r6_2_hc = r2_2*r2_2*r2_2   !will set r6 to 1/r6 later, needed for softcore
                r2_2 = 1._8/r2_2
                                r_2  = sqrt(r2_2)
                r6_3_hc = r2_3*r2_3*r2_3   !will set r6 to 1/r6 later, needed for softcore
                r2_3 = 1._8/r2_3
                                r_3  = sqrt(r2_3)

                !reset potential
                dv1 = 0
                dv2 = 0
                dv3 = 0
                iaci = iac(i)
                ilj1 = ljcod(iac1, iaci)
                ilj2 = ljcod(iac2, iaci)
                ilj3 = ljcod(iac3, iaci)
                if(.not. qvdw_flag) then
                        !use same lj params for all states
                                                r6_1 = 1._8/r6_1_hc     !softcore hack
                                                r6_2 = 1._8/r6_2_hc     !softcore hack
                                                r6_3 = 1._8/r6_3_hc     !softcore hack

                        v_a1  = iaclib(iaci)%avdw(ilj1)*a1(ilj1)*r6_1*r6_1 
                        v_b1  = iaclib(iaci)%bvdw(ilj1)*b1(ilj1)*r6_1
                        v_a2 = iaclib(iaci)%avdw(ilj2)*a2(ilj2)*r6_2*r6_2
                        v_b2 = iaclib(iaci)%bvdw(ilj2)*b2(ilj2)*r6_2
                        v_a3 = iaclib(iaci)%avdw(ilj3)*a3(ilj3)*r6_3*r6_3
                        v_b3 = iaclib(iaci)%bvdw(ilj3)*b3(ilj3)*r6_3
                end if
                do istate = 1, nstates
                        ! for every state:

                        ! calculate v_a:s and v_b:s for each state
                        if (qvdw_flag) then
                                r6_1_sc = r6_1_hc + sc_lookup(iq,iac1,istate)
                                r6_2_sc = r6_2_hc + sc_lookup(iq,iac2,istate)
                                r6_3_sc = r6_3_hc + sc_lookup(iq,iac3,istate)
                                v_a1 = qavdw(qiac(iq,istate),1)*a1(ilj1)/(r6_1_sc*r6_1_sc)
                                v_b1 = qbvdw(qiac(iq,istate),1)*b1(ilj1)/(r6_1_sc)
                                v_a2 = qavdw(qiac(iq,istate),1)*a2(ilj2)/(r6_2_sc*r6_2_sc)
                                v_b2 = qbvdw(qiac(iq,istate),1)*b2(ilj2)/(r6_2_sc)
                                v_a3 = qavdw(qiac(iq,istate),1)*a3(ilj3)/(r6_3_sc*r6_3_sc)
                                v_b3 = qbvdw(qiac(iq,istate),1)*b3(ilj3)/(r6_3_sc)
                        end if
                        ! calculate  vel, v_a, v_b and dv
                        vel1 = crg1*qcrg(iq,istate)*r_1
                        vel2 = crg2*qcrg(iq,istate)*r_2
                        vel3 = crg3*qcrg(iq,istate)*r_3
                        dv1 = dv1 + r2_1*(-vel1- ((12.*v_a1-6.*v_b1)*(r6_1_hc/r6_1_sc)) )*eq(istate)%lambda
                        dv2 = dv2 + r2_2*(-vel2- ((12.*v_a2-6.*v_b2)*(r6_2_hc/r6_2_sc)) )*eq(istate)%lambda
                        dv3 = dv3 + r2_3*(-vel3- ((12.*v_a3-6.*v_b3)*(r6_3_hc/r6_3_sc)) )*eq(istate)%lambda
                        ! update q-water energies
                        eq(istate)%qw%el  = eq(istate)%qw%el + vel1 + vel2 + vel3
                        eq(istate)%qw%vdw = eq(istate)%qw%vdw + v_a1 - v_b1 &
                                + v_a2 - v_b2 + v_a3 - v_b3
                end do !istate                  
                ! update forces on q-atom
                d(3*i-2) = d(3*i-2) - dv1*dx1 - dv2*dx2 - dv3*dx3
                d(3*i-1) = d(3*i-1) - dv1*dy1 - dv2*dy2 - dv3*dy3
                d(3*i  ) = d(3*i  ) - dv1*dz1 - dv2*dz2 - dv3*dz3
                ! update forces on water
                d(3*j-2) = d(3*j-2) + dv1*dx1
                d(3*j-1) = d(3*j-1) + dv1*dy1
                d(3*j  ) = d(3*j  ) + dv1*dz1
                d(3*j+1) = d(3*j+1) + dv2*dx2
                d(3*j+2) = d(3*j+2) + dv2*dy2
                d(3*j+3) = d(3*j+3) + dv2*dz2
                d(3*j+4) = d(3*j+4) + dv3*dx3
                d(3*j+5) = d(3*j+5) + dv3*dy3
                d(3*j+6) = d(3*j+6) + dv3*dz3

        end do !iq
end do !jw
end subroutine nonbond_qw_3atom

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_qw_3atom_box
        !calculate non-bonded interactions between q-atoms and 3-atom solvent molecules
        !using geometric comb. rule
        ! local variables
        integer                                         :: jw,iq,i,j,ilj1, ilj2, ilj3, iaci
        integer                                         :: istate
        real(8)                                         ::      dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3
        real(8)                                         ::      r_1, r2_1, r6_1, r_2, r2_2, r6_2, r_3, r2_3, r6_3
        real(8)                                         ::      vel1, vel2, vel3, dv1, dv2, dv3
        real(8)                                         :: v_a1,v_b1, v_a2, v_b2, v_a3, v_b3,r6_1_sc,r6_2_sc,r6_3_sc,r6_1_hc,r6_2_hc,r6_3_hc
        real(8)                                         :: boxshiftx, boxshifty, boxshiftz, dx, dy, dz
        real(8), save                           ::      a1(2), b1(2), a2(2), b2(2), a3(2), b3(2)
        integer, save                           ::      iac1, iac2, iac3
        real, save                                      ::      crg1, crg2, crg3
        ! global variables used:
        !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, eq, d, nat_solute
  
        if(a1(1) == 0.) then !set first time
                iac1 = iac(nat_solute+1)
                iac2 = iac(nat_solute+2)
                iac3 = iac(nat_solute+3)
                crg1 = crg(nat_solute+1)
                crg2 = crg(nat_solute+2)
                crg3 = crg(nat_solute+3)
                a1(1:2) = iaclib(iac1)%avdw(1:2)
                b1(1:2) = iaclib(iac1)%bvdw(1:2)
                a2(1:2) = iaclib(iac2)%avdw(1:2)
                b2(1:2) = iaclib(iac2)%bvdw(1:2)
                a3(1:2) = iaclib(iac3)%avdw(1:2)
                b3(1:2) = iaclib(iac3)%bvdw(1:2)
        end if

        !loop over listed waters
        do jw = 1, nbqw_pair
                j = nbqw(jw) !top. # of o in water iw
                
                !compute the periodical shift
                dx = x(3*j-2) - x(3*qswitch-2)
                dy = x(3*j-1) - x(3*qswitch-1)
                dz = x(3*j  ) - x(3*qswitch  )
                boxshiftx = boxlength(1)*nint(dx*inv_boxl(1))
                boxshifty = boxlength(2)*nint(dy*inv_boxl(2))
                boxshiftz = boxlength(3)*nint(dz*inv_boxl(3))
                                
                !loop over all q-atoms
                do iq = 1, nqat
                        i = iqseq(iq)
                        dx1 = x(3*j-2) - x(3*i-2)
                        dy1 = x(3*j-1) - x(3*i-1)
                        dz1 = x(3*j  ) - x(3*i  )
                        dx2 = x(3*j+1) - x(3*i-2)
                        dy2 = x(3*j+2) - x(3*i-1)
                        dz2 = x(3*j+3) - x(3*i  )
                        dx3 = x(3*j+4) - x(3*i-2)
                        dy3 = x(3*j+5) - x(3*i-1)
                        dz3 = x(3*j+6) - x(3*i  )
                        dx1 = dx1 - boxshiftx  
                        dy1 = dy1 - boxshifty  
                        dz1 = dz1 - boxshiftz  
                        dx2 = dx2 - boxshiftx 
                        dy2 = dy2 - boxshifty 
                        dz2 = dz2 - boxshiftz           
                        dx3 = dx3 - boxshiftx 
                        dy3 = dy3 - boxshifty 
                        dz3 = dz3 - boxshiftz  
                        
            r2_1 = dx1*dx1 + dy1*dy1 + dz1*dz1
            r2_2 = dx2*dx2 + dy2*dy2 + dz2*dz2
            r2_3 = dx3*dx3 + dy3*dy3 + dz3*dz3
            r6_1_hc = r2_1*r2_1*r2_1   !will set r6 to 1/r6 later, needed for softcore
            r2_1 = 1._8/r2_1
                        r_1  = sqrt(r2_1)
                        r6_2_hc = r2_2*r2_2*r2_2   !will set r6 to 1/r6 later, needed for softcore
                        r2_2 = 1._8/r2_2
                        r_2  = sqrt(r2_2)
                        r6_3_hc = r2_3*r2_3*r2_3   !will set r6 to 1/r6 later, needed for softcore
                        r2_3 = 1._8/r2_3
                        r_3  = sqrt(r2_3)

                        !reset potential
                        dv1 = 0
                        dv2 = 0
                        dv3 = 0
                        iaci = iac(i)
                        ilj1 = ljcod(iac1, iaci)
                        ilj2 = ljcod(iac2, iaci)
                        ilj3 = ljcod(iac3, iaci)
                        if(.not. qvdw_flag) then
                                !use same lj params for all states
                                r6_1 = 1._8/r6_1_hc     !softcore hack
                                r6_2 = 1._8/r6_2_hc     !softcore hack
                                r6_3 = 1._8/r6_3_hc     !softcore hack

                                v_a1  = iaclib(iaci)%avdw(ilj1)*a1(ilj1)*r6_1*r6_1 
                                v_b1  = iaclib(iaci)%bvdw(ilj1)*b1(ilj1)*r6_1
                                v_a2 = iaclib(iaci)%avdw(ilj2)*a2(ilj2)*r6_2*r6_2
                                v_b2 = iaclib(iaci)%bvdw(ilj2)*b2(ilj2)*r6_2
                                v_a3 = iaclib(iaci)%avdw(ilj3)*a3(ilj3)*r6_3*r6_3
                                v_b3 = iaclib(iaci)%bvdw(ilj3)*b3(ilj3)*r6_3
                        end if
                        do istate = 1, nstates
                                ! for every state:

                                ! calculate v_a:s and v_b:s for each state
                                if (qvdw_flag) then
                                        r6_1_sc = r6_1_hc + sc_lookup(iq,iac1,istate)
                                        r6_2_sc = r6_2_hc + sc_lookup(iq,iac2,istate)
                                        r6_3_sc = r6_3_hc + sc_lookup(iq,iac3,istate)
                                        v_a1 = qavdw(qiac(iq,istate),1)*a1(ilj1)/(r6_1_sc*r6_1_sc)
                                        v_b1 = qbvdw(qiac(iq,istate),1)*b1(ilj1)/(r6_1_sc)
                                        v_a2 = qavdw(qiac(iq,istate),1)*a2(ilj2)/(r6_2_sc*r6_2_sc)
                                        v_b2 = qbvdw(qiac(iq,istate),1)*b2(ilj2)/(r6_2_sc)
                                        v_a3 = qavdw(qiac(iq,istate),1)*a3(ilj3)/(r6_3_sc*r6_3_sc)
                                        v_b3 = qbvdw(qiac(iq,istate),1)*b3(ilj3)/(r6_3_sc)

                                end if
                                ! calculate  vel, v_a, v_b and dv
                                vel1 = crg1*qcrg(iq,istate)*r_1
                                vel2 = crg2*qcrg(iq,istate)*r_2
                                vel3 = crg3*qcrg(iq,istate)*r_3
                                dv1 = dv1 + r2_1*(-vel1- ((12.*v_a1-6.*v_b1)*(r6_1_hc/r6_1_sc)) )*eq(istate)%lambda
                                dv2 = dv2 + r2_2*(-vel2- ((12.*v_a2-6.*v_b2)*(r6_2_hc/r6_2_sc)) )*eq(istate)%lambda
                                dv3 = dv3 + r2_3*(-vel3- ((12.*v_a3-6.*v_b3)*(r6_3_hc/r6_3_sc)) )*eq(istate)%lambda
                                ! update q-water energies
                                eq(istate)%qw%el  = eq(istate)%qw%el + vel1 + vel2 + vel3
                                eq(istate)%qw%vdw = eq(istate)%qw%vdw + v_a1 - v_b1 &
                                        + v_a2 - v_b2 + v_a3 - v_b3
                        end do !istate                  
                        ! update forces on q-atom
                        d(3*i-2) = d(3*i-2) - dv1*dx1 - dv2*dx2 - dv3*dx3
                        d(3*i-1) = d(3*i-1) - dv1*dy1 - dv2*dy2 - dv3*dy3
                        d(3*i  ) = d(3*i  ) - dv1*dz1 - dv2*dz2 - dv3*dz3
                        ! update forces on water
                        d(3*j-2) = d(3*j-2) + dv1*dx1
                        d(3*j-1) = d(3*j-1) + dv1*dy1
                        d(3*j  ) = d(3*j  ) + dv1*dz1
                        d(3*j+1) = d(3*j+1) + dv2*dx2
                        d(3*j+2) = d(3*j+2) + dv2*dy2
                        d(3*j+3) = d(3*j+3) + dv2*dz2
                        d(3*j+4) = d(3*j+4) + dv3*dx3
                        d(3*j+5) = d(3*j+5) + dv3*dy3
                        d(3*j+6) = d(3*j+6) + dv3*dz3

                end do !iq
        end do !jw
end subroutine nonbond_qw_3atom_box

!-----------------------------------------------------------------------

subroutine nonbond_ww_spc
! local variables
integer                                         :: iw,ip,i,j,i3,j3, ia
integer                                         :: ipstart
real(8)                                         :: rox, rh1x, rh2x, r2
real(8)                                         :: dxox,  dyox,  dzox
real(8)                                         :: dxh1x, dyh1x, dzh1x
real(8)                                         :: dxh2x, dyh2x, dzh2x
real(8)                                         :: vel,v_a,v_b,dv
real(8), save                                   :: a_oo, b_oo
integer                                         ::      iac_ow, iac_hw

! global variables used:
!  iaclib, nat_solute, x, crg_ow, e, d, crg_hw

if(a_oo == 0) then !initialize 'em!
iac_ow = iac(nat_solute + 1)
iac_hw = iac(nat_solute + 2)
a_oo = iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow)) &
        *iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow))
b_oo = iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow)) &
        *iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow))
end if

ipstart = 1

do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
ip = ipstart

do while (nbww(ip) .ne. 0)
  ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

  ! --- o - (o,h1,h2) ---
  i3 = (nat_solute+3*(iw-1))*3+1 !point to x for o in mol. iw
  j3 = nbww(ip)*3-2 !point to o in interacting mol.
  ! o - o (x=o)
  dxox = x((j3  ))-x((i3  ))
  dyox = x((j3+1))-x((i3+1))
  dzox = x((j3+2))-x((i3+2))
  rox = dxox*dxox+dyox*dyox+dzox*dzox
  ! o-h1 (x=h1)
  j3 = j3 + 3

  dxh1x = x((j3  ))-x((i3  ))
  dyh1x = x((j3+1))-x((i3+1))
  dzh1x = x((j3+2))-x((i3+2))
  rh1x = dxh1x*dxh1x+dyh1x*dyh1x+dzh1x*dzh1x
  ! o-h2 (x=h2)
  j3 = j3 + 3
  dxh2x = x((j3  ))-x((i3  ))
  dyh2x = x((j3+1))-x((i3+1))
  dzh2x = x((j3+2))-x((i3+2))
  rh2x = dxh2x*dxh2x+dyh2x*dyh2x+dzh2x*dzh2x
  rox=  sqrt(1._8/rox  )
  rh1x= sqrt(1._8/rh1x )
  rh2x= sqrt(1._8/rh2x )
  ! o-o 
  ! lj only for o-o
  r2 = rox * rox
  vel = crg_ow*crg_ow*rox
  v_a  = a_oo*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_oo*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 - 6 !move pointer back to o in interacting mol.
  d((i3  )) = d((i3  )) -(dv*dxox)
  d((j3  )) = d((j3  )) +(dv*dxox)
  d((i3+1)) = d((i3+1)) -(dv*dyox)
  d((j3+1)) = d((j3+1)) +(dv*dyox)
  d((i3+2)) = d((i3+2)) -(dv*dzox)
  d((j3+2)) = d((j3+2)) +(dv*dzox)
  ! o-h1 
  r2 = rh1x * rh1x
  vel = crg_ow*crg_hw*rh1x
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 +3 !point to h1 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxh1x)
  d((j3  )) = d((j3  )) +(dv*dxh1x)
  d((i3+1)) = d((i3+1)) -(dv*dyh1x)
  d((j3+1)) = d((j3+1)) +(dv*dyh1x)
  d((i3+2)) = d((i3+2)) -(dv*dzh1x)
  d((j3+2)) = d((j3+2)) +(dv*dzh1x)
  ! o-h2
  r2 = rh2x * rh2x
  vel = crg_ow*crg_hw*rh2x
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 +3 !point to h2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxh2x)
  d((j3  )) = d((j3  )) +(dv*dxh2x)
  d((i3+1)) = d((i3+1)) -(dv*dyh2x)
  d((j3+1)) = d((j3+1)) +(dv*dyh2x)
  d((i3+2)) = d((i3+2)) -(dv*dzh2x)
  d((j3+2)) = d((j3+2)) +(dv*dzh2x)

  ! --- h1 - (o,h1,h2) ---
  i3 = i3 + 3 !point to x for h1 in mol. iw
  j3 = nbww(ip)*3-2 !point to o in j-mol.
  ! h1 - o (x=o)
  dxox = x((j3  ))-x((i3  ))
  dyox = x((j3+1))-x((i3+1))
  dzox = x((j3+2))-x((i3+2))
  rox = dxox*dxox+dyox*dyox+dzox*dzox
  ! h1-h1 (x=h1)
  j3 = j3 + 3
  dxh1x = x((j3  ))-x((i3  ))
  dyh1x = x((j3+1))-x((i3+1))
  dzh1x = x((j3+2))-x((i3+2))
  rh1x = dxh1x*dxh1x+dyh1x*dyh1x+dzh1x*dzh1x
  ! h1-h2 (x=h2)
  j3 = j3 + 3
  dxh2x = x((j3  ))-x((i3  ))
  dyh2x = x((j3+1))-x((i3+1))
  dzh2x = x((j3+2))-x((i3+2))
  rh2x = dxh2x*dxh2x+dyh2x*dyh2x+dzh2x*dzh2x
  rox=  sqrt(1._8/rox  )
  rh1x= sqrt(1._8/rh1x )
  rh2x= sqrt(1._8/rh2x )
  ! h1-o 
  r2 = rox * rox
  vel = crg_hw*crg_ow*rox
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 - 6 !move pointer back to o in j-mol.
  d((i3  )) = d((i3  )) -(dv*dxox)
  d((j3  )) = d((j3  )) +(dv*dxox)
  d((i3+1)) = d((i3+1)) -(dv*dyox)
  d((j3+1)) = d((j3+1)) +(dv*dyox)
  d((i3+2)) = d((i3+2)) -(dv*dzox)
  d((j3+2)) = d((j3+2)) +(dv*dzox)
  ! h1-h1 
  r2 = rh1x * rh1x
  vel = crg_hw*crg_hw*rh1x
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 +3 !point to h1 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxh1x)
  d((j3  )) = d((j3  )) +(dv*dxh1x)
  d((i3+1)) = d((i3+1)) -(dv*dyh1x)
  d((j3+1)) = d((j3+1)) +(dv*dyh1x)
  d((i3+2)) = d((i3+2)) -(dv*dzh1x)
  d((j3+2)) = d((j3+2)) +(dv*dzh1x)
  ! h1-h2
  r2 = rh2x * rh2x
  vel = crg_hw*crg_hw*rh2x
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 +3 !point to h2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxh2x)
  d((j3  )) = d((j3  )) +(dv*dxh2x)
  d((i3+1)) = d((i3+1)) -(dv*dyh2x)
  d((j3+1)) = d((j3+1)) +(dv*dyh2x)
  d((i3+2)) = d((i3+2)) -(dv*dzh2x)
  d((j3+2)) = d((j3+2)) +(dv*dzh2x)

  ! --- h2 - (o,h1,h2) ---
  i3 = i3 + 3 !point to x for h2 in mol. iw
  j3 = nbww(ip)*3-2 !point to o in j-mol.
  ! h2 - o (x=o)
  dxox = x((j3  ))-x((i3  ))
  dyox = x((j3+1))-x((i3+1))
  dzox = x((j3+2))-x((i3+2))
  rox = dxox*dxox+dyox*dyox+dzox*dzox
  ! h2-h1 (x=h1)
  j3 = j3 + 3
  dxh1x = x((j3  ))-x((i3  ))
  dyh1x = x((j3+1))-x((i3+1))
  dzh1x = x((j3+2))-x((i3+2))
  rh1x = dxh1x*dxh1x+dyh1x*dyh1x+dzh1x*dzh1x
  ! h2-h2 (x=h2)
  j3 = j3 + 3
  dxh2x = x((j3  ))-x((i3  ))
  dyh2x = x((j3+1))-x((i3+1))
  dzh2x = x((j3+2))-x((i3+2))
  rh2x = dxh2x*dxh2x+dyh2x*dyh2x+dzh2x*dzh2x
  rox=  sqrt(1._8/rox  )
  rh1x= sqrt(1._8/rh1x )
  rh2x= sqrt(1._8/rh2x )
  ! h2-o 
  r2 = rox * rox
  vel = crg_hw*crg_ow*rox
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 - 6 !move pointer back to o in j-mol.
  d((i3  )) = d((i3  )) -(dv*dxox)
  d((j3  )) = d((j3  )) +(dv*dxox)
  d((i3+1)) = d((i3+1)) -(dv*dyox)
  d((j3+1)) = d((j3+1)) +(dv*dyox)
  d((i3+2)) = d((i3+2)) -(dv*dzox)
  d((j3+2)) = d((j3+2)) +(dv*dzox)
  ! h2-h1 
  r2 = rh1x * rh1x
  vel = crg_hw*crg_hw*rh1x
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 +3 !point to h1 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxh1x)
  d((j3  )) = d((j3  )) +(dv*dxh1x)
  d((i3+1)) = d((i3+1)) -(dv*dyh1x)
  d((j3+1)) = d((j3+1)) +(dv*dyh1x)
  d((i3+2)) = d((i3+2)) -(dv*dzh1x)
  d((j3+2)) = d((j3+2)) +(dv*dzh1x)
  ! h2-h2
  r2 = rh2x * rh2x
  vel = crg_hw*crg_hw*rh2x
  e%ww%el = e%ww%el + vel
  dv   = r2*( -vel)
  j3 = j3 +3 !point to h2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dxh2x)
  d((j3  )) = d((j3  )) +(dv*dxh2x)
  d((i3+1)) = d((i3+1)) -(dv*dyh2x)
  d((j3+1)) = d((j3+1)) +(dv*dyh2x)
  d((i3+2)) = d((i3+2)) -(dv*dzh2x)
  d((j3+2)) = d((j3+2)) +(dv*dzh2x)

  ip = ip + 1
end do ! while ip

! skip the gap
ipstart = ip + 1                        
end do ! iw

end subroutine nonbond_ww_spc

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_ww_spc_box
  ! local variables
  integer                                               :: iw,ip,i,j,i3,j3, ia
  integer                                               :: ipstart
  real(8)                                               :: rox, rh1x, rh2x, r2
  real(8)                                               :: dxox,  dyox,  dzox
  real(8)                                               :: dxh1x, dyh1x, dzh1x
  real(8)                                               :: dxh2x, dyh2x, dzh2x
  real(8)                                               :: vel,v_a,v_b,dv
  real(8), save                                 :: a_oo, b_oo
  integer                                               ::      iac_ow, iac_hw
  real(8)                                               :: boxshiftx, boxshifty, boxshiftz

  ! global variables used:
  !  iaclib, nat_solute, x, crg_ow, e, d, crg_hw

  if(a_oo == 0) then !initialize 'em!
        iac_ow = iac(nat_solute + 1)
        iac_hw = iac(nat_solute + 2)
        a_oo = iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow)) &
                *iaclib(iac_ow)%avdw(ljcod(iac_ow, iac_ow))
        b_oo = iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow)) &
                *iaclib(iac_ow)%bvdw(ljcod(iac_ow, iac_ow))
  end if

  ipstart = 1

  do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
        ip = ipstart

        do while (nbww(ip) .ne. 0)
          ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

          ! --- o - (o,h1,h2) ---
          i3 = (nat_solute+3*(iw-1))*3+1 !point to x for o in mol. iw
          j3 = nbww(ip)*3-2 !point to o in interacting mol.
          ! o - o (x=o)
          dxox = x((j3  ))-x((i3  ))
          dyox = x((j3+1))-x((i3+1))
          dzox = x((j3+2))-x((i3+2))
          !the periodical shift
          boxshiftx = boxlength(1)*nint( dxox*inv_boxl(1) )
          boxshifty = boxlength(2)*nint( dyox*inv_boxl(2) )
          boxshiftz = boxlength(3)*nint( dzox*inv_boxl(3) )
          ! o-h1 (x=h1)
          j3 = j3 + 3
          dxh1x = x((j3  ))-x((i3  ))
          dyh1x = x((j3+1))-x((i3+1))
          dzh1x = x((j3+2))-x((i3+2))
          ! o-h2 (x=h2)
          j3 = j3 + 3
          dxh2x = x((j3  ))-x((i3  ))
          dyh2x = x((j3+1))-x((i3+1))
          dzh2x = x((j3+2))-x((i3+2))
          dxox = dxox - boxshiftx
          dyox = dyox - boxshifty
          dzox = dzox - boxshiftz
          dxh1x = dxh1x - boxshiftx
          dyh1x = dyh1x - boxshifty
          dzh1x = dzh1x - boxshiftz
          dxh2x = dxh2x - boxshiftx
          dyh2x = dyh2x - boxshifty
          dzh2x = dzh2x - boxshiftz

          rox = dxox*dxox + dyox*dyox + dzox*dzox
          rh1x = dxh1x*dxh1x + dyh1x*dyh1x + dzh1x*dzh1x
          rh2x = dxh2x*dxh2x + dyh2x*dyh2x + dzh2x*dzh2x

          rox=  sqrt(1._8/rox  )
          rh1x= sqrt(1._8/rh1x )
          rh2x= sqrt(1._8/rh2x )
          ! o-o 
          ! lj only for o-o
          r2 = rox * rox
          vel = crg_ow*crg_ow*rox
          v_a  = a_oo*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_oo*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 - 6 !move pointer back to o in interacting mol.
          d((i3  )) = d((i3  )) -(dv*dxox)
          d((j3  )) = d((j3  )) +(dv*dxox)
          d((i3+1)) = d((i3+1)) -(dv*dyox)
          d((j3+1)) = d((j3+1)) +(dv*dyox)
          d((i3+2)) = d((i3+2)) -(dv*dzox)
          d((j3+2)) = d((j3+2)) +(dv*dzox)
          ! o-h1 
          r2 = rh1x * rh1x
          vel = crg_ow*crg_hw*rh1x
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 +3 !point to h1 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dxh1x)
          d((j3  )) = d((j3  )) +(dv*dxh1x)
          d((i3+1)) = d((i3+1)) -(dv*dyh1x)
          d((j3+1)) = d((j3+1)) +(dv*dyh1x)
          d((i3+2)) = d((i3+2)) -(dv*dzh1x)
          d((j3+2)) = d((j3+2)) +(dv*dzh1x)
          ! o-h2
          r2 = rh2x * rh2x
          vel = crg_ow*crg_hw*rh2x
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 +3 !point to h2 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dxh2x)
          d((j3  )) = d((j3  )) +(dv*dxh2x)
          d((i3+1)) = d((i3+1)) -(dv*dyh2x)
          d((j3+1)) = d((j3+1)) +(dv*dyh2x)
          d((i3+2)) = d((i3+2)) -(dv*dzh2x)
          d((j3+2)) = d((j3+2)) +(dv*dzh2x)
                        
          ! --- h1 - (o,h1,h2) ---
          i3 = i3 + 3 !point to x for h1 in mol. iw
          j3 = nbww(ip)*3-2 !point to o in j-mol.
          ! h1 - o (x=o)
          dxox = x((j3  ))-x((i3  ))
          dyox = x((j3+1))-x((i3+1))
          dzox = x((j3+2))-x((i3+2))
          ! h1-h1 (x=h1)
          j3 = j3 + 3
          dxh1x = x((j3  ))-x((i3  ))
          dyh1x = x((j3+1))-x((i3+1))
          dzh1x = x((j3+2))-x((i3+2))
          ! h1-h2 (x=h2)
          j3 = j3 + 3
          dxh2x = x((j3  ))-x((i3  ))
          dyh2x = x((j3+1))-x((i3+1))
          dzh2x = x((j3+2))-x((i3+2))
          dxox = dxox - boxshiftx
          dyox = dyox - boxshifty
          dzox = dzox - boxshiftz
          dxh1x = dxh1x - boxshiftx
          dyh1x = dyh1x - boxshifty
          dzh1x = dzh1x - boxshiftz
          dxh2x = dxh2x - boxshiftx
          dyh2x = dyh2x - boxshifty
          dzh2x = dzh2x - boxshiftz

          rox = dxox*dxox + dyox*dyox + dzox*dzox
          rh1x = dxh1x*dxh1x + dyh1x*dyh1x + dzh1x*dzh1x
          rh2x = dxh2x*dxh2x + dyh2x*dyh2x + dzh2x*dzh2x
          rox=  sqrt(1._8/rox  )
          rh1x= sqrt(1._8/rh1x )
          rh2x= sqrt(1._8/rh2x )
          ! h1-o 
          r2 = rox * rox
          vel = crg_hw*crg_ow*rox
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 - 6 !move pointer back to o in j-mol.
          d((i3  )) = d((i3  )) -(dv*dxox)
          d((j3  )) = d((j3  )) +(dv*dxox)
          d((i3+1)) = d((i3+1)) -(dv*dyox)
          d((j3+1)) = d((j3+1)) +(dv*dyox)
          d((i3+2)) = d((i3+2)) -(dv*dzox)
          d((j3+2)) = d((j3+2)) +(dv*dzox)
  ! h1-h1 
          r2 = rh1x * rh1x
          vel = crg_hw*crg_hw*rh1x
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 +3 !point to h1 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dxh1x)
          d((j3  )) = d((j3  )) +(dv*dxh1x)
          d((i3+1)) = d((i3+1)) -(dv*dyh1x)
          d((j3+1)) = d((j3+1)) +(dv*dyh1x)
          d((i3+2)) = d((i3+2)) -(dv*dzh1x)
          d((j3+2)) = d((j3+2)) +(dv*dzh1x)
          ! h1-h2
          r2 = rh2x * rh2x
          vel = crg_hw*crg_hw*rh2x
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 +3 !point to h2 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dxh2x)
          d((j3  )) = d((j3  )) +(dv*dxh2x)
          d((i3+1)) = d((i3+1)) -(dv*dyh2x)
          d((j3+1)) = d((j3+1)) +(dv*dyh2x)
          d((i3+2)) = d((i3+2)) -(dv*dzh2x)
          d((j3+2)) = d((j3+2)) +(dv*dzh2x)
                        
          ! --- h2 - (o,h1,h2) ---
          i3 = i3 + 3 !point to x for h2 in mol. iw
          j3 = nbww(ip)*3-2 !point to o in j-mol.
          ! h2 - o (x=o)
          dxox = x((j3  ))-x((i3  ))
          dyox = x((j3+1))-x((i3+1))
          dzox = x((j3+2))-x((i3+2))
          ! h2-h1 (x=h1)
          j3 = j3 + 3
          dxh1x = x((j3  ))-x((i3  ))
          dyh1x = x((j3+1))-x((i3+1))
          dzh1x = x((j3+2))-x((i3+2))
          ! h2-h2 (x=h2)
          j3 = j3 + 3
          dxh2x = x((j3  ))-x((i3  ))
          dyh2x = x((j3+1))-x((i3+1))
          dzh2x = x((j3+2))-x((i3+2))
          dxox = dxox - boxshiftx
          dyox = dyox - boxshifty
          dzox = dzox - boxshiftz
          dxh1x = dxh1x - boxshiftx
          dyh1x = dyh1x - boxshifty
          dzh1x = dzh1x - boxshiftz
          dxh2x = dxh2x - boxshiftx
          dyh2x = dyh2x - boxshifty
          dzh2x = dzh2x - boxshiftz

          rox = dxox*dxox + dyox*dyox + dzox*dzox
          rh1x = dxh1x*dxh1x + dyh1x*dyh1x + dzh1x*dzh1x
          rh2x = dxh2x*dxh2x + dyh2x*dyh2x + dzh2x*dzh2x
          rox=  sqrt(1._8/rox  )
          rh1x= sqrt(1._8/rh1x )
          rh2x= sqrt(1._8/rh2x )
          ! h2-o 
          r2 = rox * rox
          vel = crg_hw*crg_ow*rox
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 - 6 !move pointer back to o in j-mol.
          d((i3  )) = d((i3  )) -(dv*dxox)
          d((j3  )) = d((j3  )) +(dv*dxox)
          d((i3+1)) = d((i3+1)) -(dv*dyox)
          d((j3+1)) = d((j3+1)) +(dv*dyox)
          d((i3+2)) = d((i3+2)) -(dv*dzox)
          d((j3+2)) = d((j3+2)) +(dv*dzox)
          ! h2-h1 
          r2 = rh1x * rh1x
          vel = crg_hw*crg_hw*rh1x
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 +3 !point to h1 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dxh1x)
          d((j3  )) = d((j3  )) +(dv*dxh1x)
          d((i3+1)) = d((i3+1)) -(dv*dyh1x)
          d((j3+1)) = d((j3+1)) +(dv*dyh1x)
          d((i3+2)) = d((i3+2)) -(dv*dzh1x)
          d((j3+2)) = d((j3+2)) +(dv*dzh1x)
          ! h2-h2
          r2 = rh2x * rh2x
          vel = crg_hw*crg_hw*rh2x
          e%ww%el = e%ww%el + vel
          dv   = r2*( -vel)
          j3 = j3 +3 !point to h2 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dxh2x)
          d((j3  )) = d((j3  )) +(dv*dxh2x)
          d((i3+1)) = d((i3+1)) -(dv*dyh2x)
          d((j3+1)) = d((j3+1)) +(dv*dyh2x)
          d((i3+2)) = d((i3+2)) -(dv*dzh2x)
          d((j3+2)) = d((j3+2)) +(dv*dzh2x)
                        
          ip = ip + 1
        end do ! while ip

        ! skip the gap
        ipstart = ip + 1                        
  end do ! iw

end subroutine nonbond_ww_spc_box

!-----------------------------------------------------------------------

subroutine nonbond_3atomsolvent
! local variables
integer                                         :: iw,ip,i,j,i3,j3, ia
integer                                         :: ipstart
real(8)                                         ::      r1x, r2x, r3x, r2
real(8)                                         ::      dx1x,  dy1x,  dz1x
real(8)                                         ::      dx2x, dy2x, dz2x
real(8)                                         ::      dx3x, dy3x, dz3x
real(8)                                         ::      vel,v_a,v_b,dv
real(8), save                                   ::      a_11, b_11, a_12, b_12, a_13, b_13
real(8), save                                   ::      a_22, b_22, a_23, b_23, a_33, b_33
real, save                                      ::      crg1, crg2, crg3
integer                                         ::      iac1, iac2, iac3

! global variables used:
!  iaclib, nat_solute, x, e, d

if(a_11 == 0.) then !initialize static variables
iac1 = iac(nat_solute+1)
iac2 = iac(nat_solute+2)
iac3 = iac(nat_solute+3)
crg1 = crg(nat_solute+1)
crg2 = crg(nat_solute+2)
crg3 = crg(nat_solute+3)

a_11 = iaclib(iac1)%avdw(ljcod(iac1, iac1)) &
        *iaclib(iac1)%avdw(ljcod(iac1, iac1))
b_11 = iaclib(iac1)%bvdw(ljcod(iac1, iac1)) &
        *iaclib(iac1)%bvdw(ljcod(iac1, iac1))
a_12 = iaclib(iac1)%avdw(ljcod(iac1, iac2)) &
        *iaclib(iac2)%avdw(ljcod(iac1, iac2))
b_12 = iaclib(iac1)%bvdw(ljcod(iac1, iac2)) &
        *iaclib(iac2)%bvdw(ljcod(iac1, iac2))
a_13 = iaclib(iac1)%avdw(ljcod(iac1, iac3)) &
        *iaclib(iac3)%avdw(ljcod(iac1, iac3))
b_13 = iaclib(iac1)%bvdw(ljcod(iac1, iac3)) &
        *iaclib(iac3)%bvdw(ljcod(iac1, iac3))
a_22 = iaclib(iac2)%avdw(ljcod(iac2, iac2)) &
        *iaclib(iac2)%avdw(ljcod(iac2, iac2))
b_22 = iaclib(iac2)%bvdw(ljcod(iac2, iac2)) &
        *iaclib(iac2)%bvdw(ljcod(iac2, iac2))
a_23 = iaclib(iac2)%avdw(ljcod(iac2, iac3)) &
        *iaclib(iac3)%avdw(ljcod(iac2, iac3))
b_23 = iaclib(iac2)%bvdw(ljcod(iac2, iac3)) &
        *iaclib(iac3)%bvdw(ljcod(iac2, iac3))
a_33 = iaclib(iac3)%avdw(ljcod(iac3, iac3)) &
        *iaclib(iac3)%avdw(ljcod(iac3, iac3))
b_33 = iaclib(iac3)%bvdw(ljcod(iac3, iac3)) &
        *iaclib(iac3)%bvdw(ljcod(iac3, iac3))
end if

ipstart = 1

do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
ip = ipstart

do while (nbww(ip) /= 0)
  ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

  ! --- 1 - (1,2,3) ---
  i3 = (nat_solute+3*(iw-1))*3+1 !point to x for atom 1 in mol. iw
  j3 = nbww(ip)*3-2 !point to 1 in interacting mol.
  ! 1 - 1
  dx1x = x((j3  ))-x((i3  ))
  dy1x = x((j3+1))-x((i3+1))
  dz1x = x((j3+2))-x((i3+2))
  r1x = dx1x*dx1x+dy1x*dy1x+dz1x*dz1x
  ! 1-2
  j3 = j3 + 3
  dx2x = x((j3  ))-x((i3  ))
  dy2x = x((j3+1))-x((i3+1))
  dz2x = x((j3+2))-x((i3+2))
  r2x = dx2x*dx2x+dy2x*dy2x+dz2x*dz2x
  ! 1-3

  j3 = j3 + 3
  dx3x = x((j3  ))-x((i3  ))
  dy3x = x((j3+1))-x((i3+1))
  dz3x = x((j3+2))-x((i3+2))
  r3x = dx3x*dx3x+dy3x*dy3x+dz3x*dz3x
  r1x=  sqrt(1._8/r1x  )
  r2x=  sqrt(1._8/r2x )
  r3x=  sqrt(1._8/r3x )
  ! 1-1 
  r2 = r1x * r1x
  vel = crg1*crg1*r1x
  v_a  = a_11*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_11*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 - 6 !move pointer back to 1 in interacting mol.
  d((i3  )) = d((i3  )) -(dv*dx1x)
  d((j3  )) = d((j3  )) +(dv*dx1x)
  d((i3+1)) = d((i3+1)) -(dv*dy1x)
  d((j3+1)) = d((j3+1)) +(dv*dy1x)
  d((i3+2)) = d((i3+2)) -(dv*dz1x)
  d((j3+2)) = d((j3+2)) +(dv*dz1x)
  ! 1-2 
  r2 = r2x * r2x
  vel = crg1*crg2*r2x
  v_a  = a_12*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_12*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 +3 !point to 2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx2x)
  d((j3  )) = d((j3  )) +(dv*dx2x)
  d((i3+1)) = d((i3+1)) -(dv*dy2x)
  d((j3+1)) = d((j3+1)) +(dv*dy2x)
  d((i3+2)) = d((i3+2)) -(dv*dz2x)
  d((j3+2)) = d((j3+2)) +(dv*dz2x)
  ! 1-3
  r2 = r3x * r3x
  vel = crg1*crg3*r3x
  v_a  = a_13*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_13*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 +3 !point to 3 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx3x)
  d((j3  )) = d((j3  )) +(dv*dx3x)
  d((i3+1)) = d((i3+1)) -(dv*dy3x)
  d((j3+1)) = d((j3+1)) +(dv*dy3x)
  d((i3+2)) = d((i3+2)) -(dv*dz3x)
  d((j3+2)) = d((j3+2)) +(dv*dz3x)

  ! --- 2 - (1,2,3) ---
  i3 = i3 + 3 !point to x for 2 in mol. iw
  j3 = nbww(ip)*3-2 !point to o in j-mol.
  ! 2 - 1 (x=1)
  dx1x = x((j3  ))-x((i3  ))
  dy1x = x((j3+1))-x((i3+1))
  dz1x = x((j3+2))-x((i3+2))
  r1x = dx1x*dx1x+dy1x*dy1x+dz1x*dz1x
  ! 2-2 (x=2)
  j3 = j3 + 3
  dx2x = x((j3  ))-x((i3  ))
  dy2x = x((j3+1))-x((i3+1))
  dz2x = x((j3+2))-x((i3+2))
  r2x = dx2x*dx2x+dy2x*dy2x+dz2x*dz2x
  ! 2-3 (x=3)
  j3 = j3 + 3
  dx3x = x((j3  ))-x((i3  ))
  dy3x = x((j3+1))-x((i3+1))
  dz3x = x((j3+2))-x((i3+2))
  r3x = dx3x*dx3x+dy3x*dy3x+dz3x*dz3x
  r1x=  sqrt(1._8/r1x  )
  r2x=  sqrt(1._8/r2x )
  r3x=  sqrt(1._8/r3x )
  ! 2-1 
  r2 = r1x * r1x
  vel = crg2*crg1*r1x
  v_a  = a_12*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_12*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 - 6 !move pointer back to o in j-mol.
  d((i3  )) = d((i3  )) -(dv*dx1x)
  d((j3  )) = d((j3  )) +(dv*dx1x)
  d((i3+1)) = d((i3+1)) -(dv*dy1x)
  d((j3+1)) = d((j3+1)) +(dv*dy1x)
  d((i3+2)) = d((i3+2)) -(dv*dz1x)
  d((j3+2)) = d((j3+2)) +(dv*dz1x)
  ! 2-2 
  r2 = r2x * r2x
  vel = crg2*crg2*r2x
  v_a  = a_22*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_22*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 +3 !point to 2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx2x)
  d((j3  )) = d((j3  )) +(dv*dx2x)
  d((i3+1)) = d((i3+1)) -(dv*dy2x)
  d((j3+1)) = d((j3+1)) +(dv*dy2x)
  d((i3+2)) = d((i3+2)) -(dv*dz2x)
  d((j3+2)) = d((j3+2)) +(dv*dz2x)
  ! 2-3
  r2 = r3x * r3x
  vel = crg2*crg3*r3x
  v_a  = a_23*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_23*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 +3 !point to 3 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx3x)
  d((j3  )) = d((j3  )) +(dv*dx3x)
  d((i3+1)) = d((i3+1)) -(dv*dy3x)
  d((j3+1)) = d((j3+1)) +(dv*dy3x)
  d((i3+2)) = d((i3+2)) -(dv*dz3x)
  d((j3+2)) = d((j3+2)) +(dv*dz3x)

  ! --- 3 - (1,2,3) ---
  i3 = i3 + 3 !point to x for 3 in mol. iw
  j3 = nbww(ip)*3-2 !point to o in j-mol.
  ! 3 - 1 (x=1)
  dx1x = x((j3  ))-x((i3  ))
  dy1x = x((j3+1))-x((i3+1))
  dz1x = x((j3+2))-x((i3+2))
  r1x = dx1x*dx1x+dy1x*dy1x+dz1x*dz1x
  ! 3-2 (x=2)
  j3 = j3 + 3
  dx2x = x((j3  ))-x((i3  ))
  dy2x = x((j3+1))-x((i3+1))
  dz2x = x((j3+2))-x((i3+2))
  r2x = dx2x*dx2x+dy2x*dy2x+dz2x*dz2x
  ! 3-3 (x=3)
  j3 = j3 + 3
  dx3x = x((j3  ))-x((i3  ))
  dy3x = x((j3+1))-x((i3+1))
  dz3x = x((j3+2))-x((i3+2))
  r3x = dx3x*dx3x+dy3x*dy3x+dz3x*dz3x
  r1x=  sqrt(1._8/r1x  )
  r2x=  sqrt(1._8/r2x )
  r3x=  sqrt(1._8/r3x )
  ! 3-1 
  r2 = r1x * r1x
  vel = crg3*crg1*r1x
  v_a  = a_13*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_13*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 - 6 !move pointer back to o in j-mol.
  d((i3  )) = d((i3  )) -(dv*dx1x)
  d((j3  )) = d((j3  )) +(dv*dx1x)
  d((i3+1)) = d((i3+1)) -(dv*dy1x)
  d((j3+1)) = d((j3+1)) +(dv*dy1x)
  d((i3+2)) = d((i3+2)) -(dv*dz1x)
  d((j3+2)) = d((j3+2)) +(dv*dz1x)
  ! 3-2 
  r2 = r2x * r2x
  vel = crg3*crg2*r2x
  v_a  = a_23*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_23*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 +3 !point to 2 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx2x)
  d((j3  )) = d((j3  )) +(dv*dx2x)
  d((i3+1)) = d((i3+1)) -(dv*dy2x)
  d((j3+1)) = d((j3+1)) +(dv*dy2x)
  d((i3+2)) = d((i3+2)) -(dv*dz2x)
  d((j3+2)) = d((j3+2)) +(dv*dz2x)
  ! 3-3
  r2 = r3x * r3x
  vel = crg3*crg3*r3x
  v_a  = a_33*(r2*r2*r2)*(r2*r2*r2)
  v_b  = b_33*(r2*r2*r2)
  e%ww%vdw = e%ww%vdw + v_a-v_b
  e%ww%el = e%ww%el + vel
  dv = r2*( -vel -12.*v_a +6.*v_b )
  j3 = j3 +3 !point to 3 in j-molecule
  d((i3  )) = d((i3  )) -(dv*dx3x)
  d((j3  )) = d((j3  )) +(dv*dx3x)
  d((i3+1)) = d((i3+1)) -(dv*dy3x)
  d((j3+1)) = d((j3+1)) +(dv*dy3x)
  d((i3+2)) = d((i3+2)) -(dv*dz3x)
  d((j3+2)) = d((j3+2)) +(dv*dz3x)

  ip = ip + 1
end do ! while ip

! skip the gap
ipstart = ip + 1                        
end do ! iw

end subroutine nonbond_3atomsolvent

!-----------------------------------------------------------------------
!******pwadded 2001-10-23
subroutine nonbond_3atomsolvent_box
  ! local variables
  integer                                               :: iw,ip,i,j,i3,j3, ia
  integer                                               :: ipstart
  real(8)                                               ::      r1x, r2x, r3x, r2
  real(8)                                               ::      dx1x,  dy1x,  dz1x
  real(8)                                               ::      dx2x, dy2x, dz2x
  real(8)                                               ::      dx3x, dy3x, dz3x
  real(8)                                               ::      vel,v_a,v_b,dv
  real(8), save                                 ::      a_11, b_11, a_12, b_12, a_13, b_13
  real(8), save                                 ::      a_22, b_22, a_23, b_23, a_33, b_33
  real, save                                    ::      crg1, crg2, crg3
  integer                                               ::      iac1, iac2, iac3
  real(8)                                               ::  boxshiftx, boxshifty, boxshiftz

  ! global variables used:
  !  iaclib, nat_solute, x, e, d

  if(a_11 == 0.) then !initialize static variables
        iac1 = iac(nat_solute+1)
        iac2 = iac(nat_solute+2)
        iac3 = iac(nat_solute+3)
        crg1 = crg(nat_solute+1)
        crg2 = crg(nat_solute+2)
        crg3 = crg(nat_solute+3)

        a_11 = iaclib(iac1)%avdw(ljcod(iac1, iac1)) &
                *iaclib(iac1)%avdw(ljcod(iac1, iac1))
        b_11 = iaclib(iac1)%bvdw(ljcod(iac1, iac1)) &
                *iaclib(iac1)%bvdw(ljcod(iac1, iac1))
        a_12 = iaclib(iac1)%avdw(ljcod(iac1, iac2)) &
                *iaclib(iac2)%avdw(ljcod(iac1, iac2))
        b_12 = iaclib(iac1)%bvdw(ljcod(iac1, iac2)) &
                *iaclib(iac2)%bvdw(ljcod(iac1, iac2))
        a_13 = iaclib(iac1)%avdw(ljcod(iac1, iac3)) &
                *iaclib(iac3)%avdw(ljcod(iac1, iac3))
        b_13 = iaclib(iac1)%bvdw(ljcod(iac1, iac3)) &
                *iaclib(iac3)%bvdw(ljcod(iac1, iac3))
        a_22 = iaclib(iac2)%avdw(ljcod(iac2, iac2)) &
                *iaclib(iac2)%avdw(ljcod(iac2, iac2))
        b_22 = iaclib(iac2)%bvdw(ljcod(iac2, iac2)) &
                *iaclib(iac2)%bvdw(ljcod(iac2, iac2))
        a_23 = iaclib(iac2)%avdw(ljcod(iac2, iac3)) &
                *iaclib(iac3)%avdw(ljcod(iac2, iac3))
        b_23 = iaclib(iac2)%bvdw(ljcod(iac2, iac3)) &
                *iaclib(iac3)%bvdw(ljcod(iac2, iac3))
        a_33 = iaclib(iac3)%avdw(ljcod(iac3, iac3)) &
                *iaclib(iac3)%avdw(ljcod(iac3, iac3))
        b_33 = iaclib(iac3)%bvdw(ljcod(iac3, iac3)) &
                *iaclib(iac3)%bvdw(ljcod(iac3, iac3))
  end if

  ipstart = 1

  do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
        ip = ipstart

        do while (nbww(ip) /= 0)
          ! consider the pair (nat_solute+3*(iw-1)) - nbww(ip)

          ! --- 1 - (1,2,3) ---
          i3 = (nat_solute+3*(iw-1))*3+1 !point to x for atom 1 in mol. iw
                                                                                !corresponds to oxygen in water
          j3 = nbww(ip)*3-2 !point to 1 in interacting mol.
                                                !corresponds to oxygen in water 
          ! 1 - 1
          dx1x = x((j3  ))-x((i3  ))
          dy1x = x((j3+1))-x((i3+1))
          dz1x = x((j3+2))-x((i3+2))
          !the periodical shift
          boxshiftx = boxlength(1)*nint( dx1x*inv_boxl(1) )
          boxshifty = boxlength(2)*nint( dy1x*inv_boxl(2) )
          boxshiftz = boxlength(3)*nint( dz1x*inv_boxl(3) )
          ! 1-2
          j3 = j3 + 3
          dx2x = x((j3  ))-x((i3  ))
          dy2x = x((j3+1))-x((i3+1))
          dz2x = x((j3+2))-x((i3+2))
          ! 1-3
          j3 = j3 + 3
          dx3x = x((j3  ))-x((i3  ))
          dy3x = x((j3+1))-x((i3+1))
          dz3x = x((j3+2))-x((i3+2))
          dx1x = dx1x - boxshiftx
          dy1x = dy1x - boxshifty
          dz1x = dz1x - boxshiftz
          dx2x = dx2x - boxshiftx
          dy2x = dy2x - boxshifty
          dz2x = dz2x - boxshiftz
          dx3x = dx3x - boxshiftx
          dy3x = dy3x - boxshifty
          dz3x = dz3x - boxshiftz

          r1x = dx1x*dx1x + dy1x*dy1x + dz1x*dz1x
          r2x = dx2x*dx2x + dy2x*dy2x + dz2x*dz2x
          r3x = dx3x*dx3x + dy3x*dy3x + dz3x*dz3x

          r1x=  sqrt(1._8/r1x  )
          r2x=  sqrt(1._8/r2x )
          r3x=  sqrt(1._8/r3x )
          ! 1-1 
          r2 = r1x * r1x
          vel = crg1*crg1*r1x
          v_a  = a_11*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_11*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 - 6 !move pointer back to 1 in interacting mol.
          d((i3  )) = d((i3  )) -(dv*dx1x)
          d((j3  )) = d((j3  )) +(dv*dx1x)
          d((i3+1)) = d((i3+1)) -(dv*dy1x)
          d((j3+1)) = d((j3+1)) +(dv*dy1x)
          d((i3+2)) = d((i3+2)) -(dv*dz1x)
          d((j3+2)) = d((j3+2)) +(dv*dz1x)
          ! 1-2 
          r2 = r2x * r2x
          vel = crg1*crg2*r2x
          v_a  = a_12*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_12*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 +3 !point to 2 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dx2x)
          d((j3  )) = d((j3  )) +(dv*dx2x)
          d((i3+1)) = d((i3+1)) -(dv*dy2x)
          d((j3+1)) = d((j3+1)) +(dv*dy2x)
          d((i3+2)) = d((i3+2)) -(dv*dz2x)
          d((j3+2)) = d((j3+2)) +(dv*dz2x)
          ! 1-3
          r2 = r3x * r3x
          vel = crg1*crg3*r3x
          v_a  = a_13*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_13*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 +3 !point to 3 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dx3x)
          d((j3  )) = d((j3  )) +(dv*dx3x)
          d((i3+1)) = d((i3+1)) -(dv*dy3x)
          d((j3+1)) = d((j3+1)) +(dv*dy3x)
          d((i3+2)) = d((i3+2)) -(dv*dz3x)
          d((j3+2)) = d((j3+2)) +(dv*dz3x)
                        
          ! --- 2 - (1,2,3) ---
          i3 = i3 + 3 !point to x for 2 in mol. iw
          j3 = nbww(ip)*3-2 !point to o in j-mol.
          ! 2 - 1 (x=1)
          dx1x = x((j3  ))-x((i3  ))
          dy1x = x((j3+1))-x((i3+1))
          dz1x = x((j3+2))-x((i3+2))
          ! 2-2 (x=2)
          j3 = j3 + 3
          dx2x = x((j3  ))-x((i3  ))
          dy2x = x((j3+1))-x((i3+1))
          dz2x = x((j3+2))-x((i3+2))
          ! 2-3 (x=3)
          j3 = j3 + 3
          dx3x = x((j3  ))-x((i3  ))
          dy3x = x((j3+1))-x((i3+1))
          dz3x = x((j3+2))-x((i3+2))
          dx1x = dx1x - boxshiftx
          dy1x = dy1x - boxshifty
          dz1x = dz1x - boxshiftz
          dx2x = dx2x - boxshiftx
          dy2x = dy2x - boxshifty
          dz2x = dz2x - boxshiftz
          dx3x = dx3x - boxshiftx
          dy3x = dy3x - boxshifty
          dz3x = dz3x - boxshiftz

          r1x = dx1x*dx1x + dy1x*dy1x + dz1x*dz1x
          r2x = dx2x*dx2x + dy2x*dy2x + dz2x*dz2x
          r3x = dx3x*dx3x + dy3x*dy3x + dz3x*dz3x
          r1x=  sqrt(1._8/r1x  )
          r2x=  sqrt(1._8/r2x )
          r3x=  sqrt(1._8/r3x )
          ! 2-1 
          r2 = r1x * r1x
          vel = crg2*crg1*r1x
          v_a  = a_12*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_12*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 - 6 !move pointer back to o in j-mol.
          d((i3  )) = d((i3  )) -(dv*dx1x)
          d((j3  )) = d((j3  )) +(dv*dx1x)
          d((i3+1)) = d((i3+1)) -(dv*dy1x)
          d((j3+1)) = d((j3+1)) +(dv*dy1x)
          d((i3+2)) = d((i3+2)) -(dv*dz1x)
          d((j3+2)) = d((j3+2)) +(dv*dz1x)
          ! 2-2 
          r2 = r2x * r2x
          vel = crg2*crg2*r2x
          v_a  = a_22*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_22*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 +3 !point to 2 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dx2x)
          d((j3  )) = d((j3  )) +(dv*dx2x)
          d((i3+1)) = d((i3+1)) -(dv*dy2x)
          d((j3+1)) = d((j3+1)) +(dv*dy2x)
          d((i3+2)) = d((i3+2)) -(dv*dz2x)
          d((j3+2)) = d((j3+2)) +(dv*dz2x)
          ! 2-3
          r2 = r3x * r3x
          vel = crg2*crg3*r3x
          v_a  = a_23*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_23*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 +3 !point to 3 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dx3x)
          d((j3  )) = d((j3  )) +(dv*dx3x)
          d((i3+1)) = d((i3+1)) -(dv*dy3x)
          d((j3+1)) = d((j3+1)) +(dv*dy3x)
          d((i3+2)) = d((i3+2)) -(dv*dz3x)
          d((j3+2)) = d((j3+2)) +(dv*dz3x)
                        
          ! --- 3 - (1,2,3) ---
          i3 = i3 + 3 !point to x for 3 in mol. iw
          j3 = nbww(ip)*3-2 !point to o in j-mol.
          ! 3 - 1 (x=1)
          dx1x = x((j3  ))-x((i3  ))
          dy1x = x((j3+1))-x((i3+1))
          dz1x = x((j3+2))-x((i3+2))
          ! 3-2 (x=2)
          j3 = j3 + 3
          dx2x = x((j3  ))-x((i3  ))
          dy2x = x((j3+1))-x((i3+1))
          dz2x = x((j3+2))-x((i3+2))
          ! 3-3 (x=3)
          j3 = j3 + 3
          dx3x = x((j3  ))-x((i3  ))
          dy3x = x((j3+1))-x((i3+1))
          dz3x = x((j3+2))-x((i3+2))
          dx1x = dx1x - boxshiftx
          dy1x = dy1x - boxshifty
          dz1x = dz1x - boxshiftz
          dx2x = dx2x - boxshiftx
          dy2x = dy2x - boxshifty
          dz2x = dz2x - boxshiftz
          dx3x = dx3x - boxshiftx
          dy3x = dy3x - boxshifty
          dz3x = dz3x - boxshiftz

          r1x = dx1x*dx1x + dy1x*dy1x + dz1x*dz1x
          r2x = dx2x*dx2x + dy2x*dy2x + dz2x*dz2x
          r3x = dx3x*dx3x + dy3x*dy3x + dz3x*dz3x
          r1x=  sqrt(1._8/r1x  )
          r2x=  sqrt(1._8/r2x )
          r3x=  sqrt(1._8/r3x )
          ! 3-1 
          r2 = r1x * r1x
          vel = crg3*crg1*r1x
          v_a  = a_13*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_13*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 - 6 !move pointer back to o in j-mol.
          d((i3  )) = d((i3  )) -(dv*dx1x)
          d((j3  )) = d((j3  )) +(dv*dx1x)
          d((i3+1)) = d((i3+1)) -(dv*dy1x)
          d((j3+1)) = d((j3+1)) +(dv*dy1x)
          d((i3+2)) = d((i3+2)) -(dv*dz1x)
          d((j3+2)) = d((j3+2)) +(dv*dz1x)
          ! 3-2 
          r2 = r2x * r2x
          vel = crg3*crg2*r2x
          v_a  = a_23*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_23*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 +3 !point to 2 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dx2x)
          d((j3  )) = d((j3  )) +(dv*dx2x)
          d((i3+1)) = d((i3+1)) -(dv*dy2x)
          d((j3+1)) = d((j3+1)) +(dv*dy2x)
          d((i3+2)) = d((i3+2)) -(dv*dz2x)
          d((j3+2)) = d((j3+2)) +(dv*dz2x)
          ! 3-3
          r2 = r3x * r3x
          vel = crg3*crg3*r3x
          v_a  = a_33*(r2*r2*r2)*(r2*r2*r2)
          v_b  = b_33*(r2*r2*r2)
          e%ww%vdw = e%ww%vdw + v_a-v_b
          e%ww%el = e%ww%el + vel
          dv = r2*( -vel -12.*v_a +6.*v_b )
          j3 = j3 +3 !point to 3 in j-molecule
          d((i3  )) = d((i3  )) -(dv*dx3x)
          d((j3  )) = d((j3  )) +(dv*dx3x)
          d((i3+1)) = d((i3+1)) -(dv*dy3x)
          d((j3+1)) = d((j3+1)) +(dv*dy3x)
          d((i3+2)) = d((i3+2)) -(dv*dz3x)
          d((j3+2)) = d((j3+2)) +(dv*dz3x)

          ip = ip + 1
        end do ! while ip

        ! skip the gap
        ipstart = ip + 1                        
  end do ! iw

end subroutine nonbond_3atomsolvent_box

!-----------------------------------------------------------------------

subroutine offdiag
! local variables
integer                                         :: io,i,j,k,l,k3,l3
real(8)                                         :: r

! global variables used:
!  offd, noffd, iqseq, x, hij, offd2

do io = 1, noffd
! for every offd:

i  = offd(io)%i
j  = offd(io)%j
k  = iqseq(offd2(io)%k)
l  = iqseq(offd2(io)%l)
k3 = k*3-3
l3 = l*3-3

r = sqrt ( (x(l3+1)-x(k3+1))**2 + &
  (x(l3+2)-x(k3+2))**2 + &
  (x(l3+3)-x(k3+3))**2 )

hij(i,j) = offd2(io)%a * exp(-offd2(io)%mu*r)
offd(io)%hij = hij(i,j) ! store for save
offd(io)%rkl = r
end do
end subroutine offdiag

!-----------------------------------------------------------------------
subroutine p_restrain
! *** local variables
integer                                         ::      ir,i,j,k,i3,j3,istate, n_ctr
real(8)                                         ::      fk,r2,erst,edum,x2,y2,z2,wgt,b,db,dv, totmass
real(8)                                         ::      dr(3), ctr(3)
real(8)                                         ::      fexp

! global variables used:
!  e, nstates, eq, nrstr_seq, rstseq, heavy, x, xtop, d, nrstr_pos, rstpos, nrstr_dist, 
!  rstdis, nrstr_wall, rstwal, xwcent

! sequence restraints (independent of q-state)
do ir = 1, nrstr_seq
fk = rstseq(ir)%fk

if(rstseq(ir)%to_centre == 1) then     ! put == 1, then equal to 2
  ! restrain to geometrical centre

  ! reset dr & atom counter
  dr(:) = 0.
  n_ctr = 0

  ! calculate deviation from center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
          n_ctr = n_ctr + 1
          dr(:) = dr(:) + x(i*3-2:i*3) - xtop(i*3-2:i*3)
        end if
  end do

  if(n_ctr > 0) then 
    ! only if atoms were found:

        ! form average
        dr(:) = dr(:) / n_ctr 
        r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
        erst    = 0.5*fk*r2
        e%restraint%protein  = e%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
          if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                d(i*3-2:i*3) = d(i*3-2:i*3) + fk*dr(:)*iaclib(iac(i))%mass/12.010
          end if
        end do
  end if
 
else if(rstseq(ir)%to_centre == 2) then     ! put == 1, then equal to 2
  ! restrain to mass centre
  ! reset dr & variable to put masses
  dr(:) = 0.
  totmass = 0.
  
! calculate deviation from mass center
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
          totmass = totmass + iaclib(iac(i))%mass                              ! add masses
          dr(:) = dr(:) + (x(i*3-2:i*3) - xtop(i*3-2:i*3))*iaclib(iac(i))%mass ! massweight distances
                end if
  end do

 if(totmass > 0) then 
    ! only if atoms were found: (i.e has a total mass)

        ! form average
        dr(:) = dr(:)/totmass                                  ! divide by total mass
        r2      = dr(1)**2 + dr(2)**2 + dr(3)**2
        erst    = 0.5*fk*r2
        e%restraint%protein  = e%restraint%protein + erst

        ! apply same force to all atoms
        do i = rstseq(ir)%i, rstseq(ir)%j
          if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
                d(i*3-2:i*3) = d(i*3-2:i*3) + fk*dr(:)
          end if
        end do
  end if

else 
  ! restrain each atom to its topology co-ordinate
  do i = rstseq(ir)%i, rstseq(ir)%j
        if ( heavy(i) .or. rstseq(ir)%ih .eq. 1 ) then
          i3 = i*3-3

          dr(1)   = x(i3+1) - xtop(i3+1)
          dr(2)   = x(i3+2) - xtop(i3+2)
          dr(3)   = x(i3+3) - xtop(i3+3)
       !use the periodically minimal distance:
       if( use_pbc ) then
               dr(1) = dr(1) - boxlength(1)*nint( dr(1)*inv_boxl(1) )
               dr(2) = dr(2) - boxlength(2)*nint( dr(2)*inv_boxl(2) )
               dr(3) = dr(3) - boxlength(3)*nint( dr(3)*inv_boxl(3) )
       end if
          r2      = dr(1)**2 + dr(2)**2 + dr(3)**2

          erst    = 0.5*fk*r2
          e%restraint%protein  = e%restraint%protein + erst

          d(i3+1) = d(i3+1) + fk*dr(1)
          d(i3+2) = d(i3+2) + fk*dr(2)
          d(i3+3) = d(i3+3) + fk*dr(3)
        end if
  end do
end if
end do

! extra positional restraints (q-state dependent)
do ir = 1, nrstr_pos
istate = rstpos(ir)%ipsi
i      = rstpos(ir)%i
i3     = i*3-3

dr(1)  = x(i3+1) - rstpos(ir)%x(1)
dr(2)  = x(i3+2) - rstpos(ir)%x(2)
dr(3)  = x(i3+3) - rstpos(ir)%x(3)

if ( istate .ne. 0 ) then
wgt = eq(istate)%lambda
else
wgt = 1.0
end if

x2 = dr(1)**2
y2 = dr(2)**2
z2 = dr(3)**2

edum = 0.5*rstpos(ir)%fk(1)*x2 + &
   0.5*rstpos(ir)%fk(2)*y2 + &
   0.5*rstpos(ir)%fk(3)*z2

d(i3+1) = d(i3+1) + rstpos(ir)%fk(1)*dr(1)*wgt
d(i3+2) = d(i3+2) + rstpos(ir)%fk(2)*dr(2)*wgt
d(i3+3) = d(i3+3) + rstpos(ir)%fk(3)*dr(3)*wgt

if ( istate .eq. 0 ) then
do k = 1, nstates
eq(k)%restraint = eq(k)%restraint + edum
end do
if ( nstates .eq. 0 ) e%restraint%protein = e%restraint%protein + edum
else
eq(istate)%restraint = eq(istate)%restraint + edum
end if
end do

! atom-atom distance restraints (q-state dependent)
do ir = 1, nrstr_dist
istate = rstdis(ir)%ipsi
i      = rstdis(ir)%i
j      = rstdis(ir)%j
i3     = i*3-3
j3     = j*3-3

dr(1)  = x(j3+1) - x(i3+1)
dr(2)  = x(j3+2) - x(i3+2)
dr(3)  = x(j3+3) - x(i3+3)

! if pbc then adjust lengths according to periodicity - ma
if( use_pbc ) then
        dr(1) = dr(1) - boxlength(1)*nint( dr(1)*inv_boxl(1) )
        dr(2) = dr(2) - boxlength(2)*nint( dr(2)*inv_boxl(2) )
        dr(3) = dr(3) - boxlength(3)*nint( dr(3)*inv_boxl(3) )
end if

if ( istate .ne. 0 ) then
wgt = eq(istate)%lambda
else
wgt = 1.0
end if

b      = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
if(b < rstdis(ir)%d1) then !shorter than d1
        db     = b - rstdis(ir)%d1
elseif(b > rstdis(ir)%d2) then !longer than d2
        db     = b - rstdis(ir)%d2
else
        db = 0
        cycle !skip zero force calculation
endif

edum   = 0.5*rstdis(ir)%fk*db**2
dv     = wgt*rstdis(ir)%fk*db/b

d(j3+1) = d(j3+1) + dr(1)*dv
d(j3+2) = d(j3+2) + dr(2)*dv
d(j3+3) = d(j3+3) + dr(3)*dv
d(i3+1) = d(i3+1) - dr(1)*dv
d(i3+2) = d(i3+2) - dr(2)*dv
d(i3+3) = d(i3+3) - dr(3)*dv

if ( istate .eq. 0 ) then
do k = 1, nstates
eq(k)%restraint = eq(k)%restraint + edum
end do
if ( nstates .eq. 0 ) e%restraint%protein = e%restraint%protein + edum
else
eq(istate)%restraint = eq(istate)%restraint + edum
end if
end do

if( .not. use_pbc ) then
! extra half-harmonic wall restraints
do ir = 1, nrstr_wall
        fk = rstwal(ir)%fk
        do i = rstwal(ir)%i, rstwal(ir)%j
                if ( heavy(i) .or. rstwal(ir)%ih .eq. 1 ) then
                        i3 = i*3-3

                        dr(1) = x(i3+1) - xwcent(1)
                        dr(2) = x(i3+2) - xwcent(2)
                        dr(3) = x(i3+3) - xwcent(3)

                        b     = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
                        db = b - rstwal(ir)%d

                        if(db > 0.) then
                                erst =  0.5 * fk * db**2 - rstwal(ir)%dmorse
                                dv = fk*db/b
                        else
                                fexp = exp(rstwal(ir)%amorse*db)
                                erst = rstwal(ir)%dmorse*(fexp*fexp-2.*fexp)
                                dv=-2.*rstwal(ir)%dmorse*rstwal(ir)%amorse*(fexp-fexp*fexp)/b
                        end if
                        e%restraint%protein = e%restraint%protein + erst

                        d(i3+1) = d(i3+1) + dv*dr(1)
                        d(i3+2) = d(i3+2) + dv*dr(2)
                        d(i3+3) = d(i3+3) + dv*dr(3)
                end if
        end do
end do

end if

end subroutine p_restrain

!-----------------------------------------------------------------------

subroutine pot_energy
! local variables


integer                                 :: istate, i, nat3
integer                                 :: is, j
#if defined (profiling)
real(8)                                 :: start_loop_time1
real(8)                                 :: start_loop_time2
real(8)                                 :: start_loop_time3
#endif

! --- reset all energies

e%potential = 0.0
!e%kinetic = 0.0        ! no need to reset because it will be assigned its final value at once
e%lrf = 0.0
e%p%bond  = 0.0
e%p%angle = 0.0
e%p%torsion = 0.0
e%p%improper = 0.0
e%w%bond  = 0.0
e%w%angle = 0.0
e%w%torsion = 0.0
e%w%improper = 0.0
e%q%bond  = 0.0
e%q%angle = 0.0
e%q%torsion   = 0.0
e%q%improper   = 0.0
e%pp%el  = 0.0
e%pp%vdw = 0.0
e%pw%el  = 0.0
e%pw%vdw = 0.0
e%ww%el  = 0.0
e%ww%vdw = 0.0
e%qx%el    = 0.0
e%qx%vdw   = 0.0
!e%restraint%total = 0.0        ! will be assigned its final value later
e%restraint%fix = 0.0
e%restraint%shell = 0.0
e%restraint%protein = 0.0
e%restraint%solvent_radial = 0.0
e%restraint%water_pol = 0.0
do istate = 1, nstates
!eq(istate)%lambda set by initialize
!eq(istate)%total assigned its final value later
eq(istate)%q%bond = 0.0
eq(istate)%q%angle = 0.0
eq(istate)%q%torsion = 0.0
eq(istate)%q%improper = 0.0
!eq(istate)%qx%el = 0.0 ! assigned its final value later
!eq(istate)%qx%vdw = 0.0        ! assigned its final value later
eq(istate)%qq%el = 0.0
eq(istate)%qq%vdw = 0.0
eq(istate)%qp%el = 0.0
eq(istate)%qp%vdw = 0.0
eq(istate)%qw%el = 0.0
eq(istate)%qw%vdw = 0.0
eq(istate)%restraint = 0.0
end do

!reset derivatives ---
d(:) = 0.

! --- calculate the potential energy and derivatives ---
! *** nonbonds distribueras

#if defined (use_mpi)
if (nodeid .eq. 0) then
!first post recieves for gathering data from slaves
call gather_nonbond
end if
#endif

#if defined (profiling)
start_loop_time2 = rtime()
#endif

! classical nonbonds
call pot_energy_nonbonds
#if defined (profiling)
profile(10)%time = profile(10)%time + rtime() - start_loop_time2
#endif

if (nodeid .eq. 0) then


! classical bond interactions (master node only)
#if defined (profiling)
start_loop_time1 = rtime()
#endif
call pot_energy_bonds
#if defined (profiling)
profile(8)%time = profile(8)%time + rtime() - start_loop_time1
#endif

! various restraints
if( .not. use_pbc ) then
   call fix_shell     !restrain all excluded atoms plus heavy solute atoms in the inner shell.
end if

call p_restrain       !seq. restraints, dist. restaints, etc

if( .not. use_pbc ) then
        if(nwat > 0) then
          call restrain_solvent 
        if (wpol_restr) call watpol
        end if
end if

! q-q nonbonded interactions
if(.not. qq_use_library_charges) then
        if(ivdw_rule .eq. 1 ) then
                call nonbond_qq
        elseif ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq
        end if
else
        if ( ivdw_rule .eq. 1 ) then
                call nonbond_qq_lib_charges
        else if ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq_lib_charges
        end if
end if

! q-atom bonded interactions: loop over q-atom states
do istate = 1, nstates
  ! bonds, angles, torsions and impropers
  call qbond (istate)
  call qangle (istate)
  if(ff_type == ff_charmm) call qurey_bradley(istate)
  call qtorsion (istate)
  call qimproper (istate)
end do
#if defined (profiling)
profile(9)%time = profile(9)%time + rtime() - start_loop_time1 - profile(8)%time
#endif
#if defined(use_mpi)
else  !slave nodes
call gather_nonbond
#endif
end if

if (nodeid .eq. 0) then 
#if (use_mpi)
do i = 1, 3
    call mpi_waitall(numnodes-1,request_recv(1,i),mpi_status,ierr)
end do

!forces and energies are summarised
do i=1,numnodes-1
  d = d + d_recv(:,i)
  e%pp%el   = e%pp%el  + e_recv(i)%pp%el
  e%pp%vdw  = e%pp%vdw + e_recv(i)%pp%vdw
  e%pw%el   = e%pw%el  + e_recv(i)%pw%el
  e%pw%vdw  = e%pw%vdw + e_recv(i)%pw%vdw
  e%ww%el   = e%ww%el  + e_recv(i)%ww%el
  e%ww%vdw  = e%ww%vdw + e_recv(i)%ww%vdw
  e%lrf     = e%lrf    + e_recv(i)%lrf
  eq(1:nstates)%qp%el  = eq(1:nstates)%qp%el  + eq_recv(1:nstates,i)%qp%el
  eq(1:nstates)%qp%vdw = eq(1:nstates)%qp%vdw + eq_recv(1:nstates,i)%qp%vdw
  eq(1:nstates)%qw%el  = eq(1:nstates)%qw%el  + eq_recv(1:nstates,i)%qw%el
  eq(1:nstates)%qw%vdw = eq(1:nstates)%qw%vdw + eq_recv(1:nstates,i)%qw%vdw
end do
#endif

! q-atom energy summary
do istate = 1, nstates
! update eq
eq(istate)%qx%el  = eq(istate)%qq%el +eq(istate)%qp%el +eq(istate)%qw%el
eq(istate)%qx%vdw = eq(istate)%qq%vdw+eq(istate)%qp%vdw+eq(istate)%qw%vdw

eq(istate)%total =  eq(istate)%q%bond + eq(istate)%q%angle   &
  + eq(istate)%q%torsion  + eq(istate)%q%improper + eq(istate)%qx%el &
  + eq(istate)%qx%vdw  + eq(istate)%restraint

! update e with an average of all states
e%q%bond  = e%q%bond  + eq(istate)%q%bond *eq(istate)%lambda
e%q%angle = e%q%angle + eq(istate)%q%angle*eq(istate)%lambda
e%q%torsion   = e%q%torsion   + eq(istate)%q%torsion  *eq(istate)%lambda
e%q%improper   = e%q%improper   + eq(istate)%q%improper  *eq(istate)%lambda
e%qx%el    = e%qx%el    + eq(istate)%qx%el   *eq(istate)%lambda
e%qx%vdw   = e%qx%vdw   + eq(istate)%qx%vdw  *eq(istate)%lambda

! update e%restraint%protein with an average of all states
e%restraint%protein = e%restraint%protein + eq(istate)%restraint*eq(istate)%lambda
end do

! total energy summary
e%restraint%total = e%restraint%fix + e%restraint%shell + &
e%restraint%protein + e%restraint%solvent_radial + e%restraint%water_pol

e%potential = e%p%bond + e%w%bond + e%p%angle + e%w%angle + e%p%torsion + &
e%p%improper + e%pp%el + e%pp%vdw + e%pw%el + e%pw%vdw + e%ww%el + &
e%ww%vdw + e%q%bond + e%q%angle + e%q%torsion + &
e%q%improper + e%qx%el + e%qx%vdw + e%restraint%total + e%lrf
end if

end subroutine pot_energy

!-----------------------------------------------------------------------

subroutine pot_energy_bonds
! bond, angle, torsion and improper potential energy

select case(ff_type)
case(ff_gromos)
        e%p%bond = bond(1, nbonds_solute)
        e%w%bond = bond(nbonds_solute+1, nbonds)
        e%p%angle = angle(1, nangles_solute)
        e%w%angle = angle(nangles_solute+1, nangles)
        e%p%torsion = torsion(1, ntors_solute)
        e%w%torsion = torsion(ntors_solute+1, ntors)
        e%p%improper = improper(1, nimps_solute)
        e%w%improper = improper(nimps_solute+1, nimps)
case(ff_amber)
        e%p%bond = bond(1, nbonds_solute)
        e%w%bond = bond(nbonds_solute+1, nbonds)
        e%p%angle = angle(1, nangles_solute)
        e%w%angle = angle(nangles_solute+1, nangles)
        e%p%torsion = torsion(1, ntors_solute)
        e%w%torsion = torsion(ntors_solute+1, ntors)
        e%p%improper = improper2(1, nimps_solute)
        e%w%improper = improper2(nimps_solute+1, nimps)
case(ff_charmm)
        e%p%bond = bond(1, nbonds_solute)
        e%w%bond = bond(nbonds_solute+1, nbonds)
        e%p%angle = angle(1, nangles_solute)
        e%w%angle = angle(nangles_solute+1, nangles)
        e%p%angle = e%p%angle + urey_bradley(1, nangles_solute)
        e%w%angle = e%w%angle + urey_bradley(nangles_solute+1, nangles)
        e%p%torsion = torsion(1, ntors_solute)
        e%w%torsion = torsion(ntors_solute+1, ntors)
        e%p%improper = improper(1, nimps_solute)
        e%w%improper = improper(nimps_solute+1, nimps)
end select

end subroutine pot_energy_bonds

!-----------------------------------------------------------------------
subroutine pot_energy_nonbonds

!nonbonded interactions

if( use_pbc ) then !periodic box
      
                select case(ivdw_rule)
        case(vdw_geometric)
                call nonbond_pp_box
                call nonbond_pw_box
                if(qvdw_flag) then
                        call nonbond_qp_qvdw_box
                else
                        call nonbond_qp_box
                end if
                if(natom > nat_solute) then !if any solvent
                        if(solvent_type == solvent_spc) then
                                !use the optimised spc routine when possible
                                call nonbond_ww_spc_box
                                call nonbond_qw_spc_box
                        elseif(solvent_type == solvent_3atom) then !otherwise calc. lj with all atoms
                                call nonbond_3atomsolvent_box
                                call nonbond_qw_3atom_box
                        end if
                end if
        case(vdw_arithmetic)
                call nonbon2_pp_box
                call nonbon2_qp_box
                if(natom > nat_solute) then !if any solvent
                        call nonbon2_pw_box
                        call nonbon2_qw_box !no spc-specific optimised routines here
                        call nonbon2_ww_box
                end if
        end select
        
                !lrf pbc
                if (use_lrf) then
                call lrf_taylor
        end if

else !simulation sphere

        select case(ivdw_rule)
        case(vdw_geometric)
                call nonbond_pp   
                call nonbond_pw   
                if(qvdw_flag) then
                        call nonbond_qp_qvdw
                else
                        call nonbond_qp
                end if
                if(natom > nat_solute) then !if any solvent
                        if(solvent_type == solvent_spc) then
                                !use the optimised spc routine when possible
                                call nonbond_ww_spc
                                call nonbond_qw_spc
                        elseif(solvent_type == solvent_3atom) then !otherwise calc. lj with all atoms
                                call nonbond_3atomsolvent
                                call nonbond_qw_3atom
                        end if
                end if
        case(vdw_arithmetic)
                call nonbon2_pp
                call nonbon2_qp
                if(natom > nat_solute) then !if any solvent
                        call nonbon2_pw
                        call nonbon2_qw !no spc-specific optimised routines here
                        call nonbon2_ww
                end if  
        end select

        ! on demand: taylor expansion of the electric field from charge groups beyond rcutoff
        if (use_lrf) then
                call lrf_taylor
        end if

end if

end subroutine pot_energy_nonbonds

!-----------------------------------------------------------------------

subroutine prep_coord


! local variables
integer(4)          :: i,nat3

! --- refresh topology coords. if needed (external restraints file)
if ( implicit_rstr_from_file .eq. 1 ) then
write (*,'(/,a,/)') 'refreshing topology coords for restraining...'
read (12) nat3,(xtop(i),i=1,nat_pro*3)
end if

!assign restraints of kind res:atom their numerical atom numbers
do i=1,nrstr_dist
  if(rstdis(i)%itext .ne. 'nil') then
    if (scan(rstdis(i)%itext,':') .ne. 0) then
      rstdis(i)%i=get_atom_from_resnum_atnum(rstdis(i)%itext)
        else
      read(rstdis(i)%itext,*) rstdis(i)%i
    end if
    if (scan(rstdis(i)%jtext,':') .ne. 0) then
      rstdis(i)%j=get_atom_from_resnum_atnum(rstdis(i)%jtext)
    else
      read(rstdis(i)%jtext,*) rstdis(i)%j
    end if
  end if
end do

! --- make spherical restraining shell lists based on
!     the xtop coords.
if (.not. use_pbc) then

        if(rexcl_i > rexcl_o) then
          call die('inner radius of restrained shell must be < exclusion radius')
        end if
        !first find atoms in shell 
        if (rexcl_i >= 0) then      !if rexcl_i is defined...
          if (rexcl_i <= 1.00) then   !if rexcl_i is defined as fraction of rexcl_o
            rexcl_i = rexcl_i * rexcl_o  !translate to angstrom
          end if 
          if(iuse_switch_atom == 1) then
             call make_shell
          else
             call make_shell2
          end if
        else
          write (*,'(/,a,/)') 'restrained shell not defined!'
        end if
else
        shell(:) = .false.
end if ! .not. use_pbc
! --- read restart file

call allocate_natom_arrays
if(restart) then
        ! topology routine has determined nwat, natom and allocated storage
        call centered_heading('reading restart file','-')
        read (2) nat3
        rewind(2)
        if(nat3 /= 3*natom) then
                write(*,100) nat3/3, natom
100                     format('>>>>> error:',i5,' atoms in restart file not equal to',i5,&
                        ' in topology.')
                call die('wrong number of atoms in restart file')
        end if
        read (2,err=112,end=112) nat3, (x(i),i=1,nat3)
        read (2,err=112,end=112) nat3, (v(i),i=1,nat3)
        write (*,'(a30,i8)')   'total number of atoms        =',natom
        write (*,'(a30,i8,/)') 'number of waters encountered =',nwat

        if( use_pbc) then
                read(2,err=112,end=112) boxlength(:)
                read(2,err=112,end=112) boxcentre(:)
                write(*,*)
                write(*,'(a16,3f8.3)') 'boxlength     =', boxlength(:)
                write(*,'(a16,3f8.3)') 'centre of box =', boxcentre(:)
        end if
        !water polarisation data will be read from restart file in wat_shells
else
        x(1:nat_pro*3) = xtop(1:nat_pro*3)
end if

! clear iqatom atom array
iqatom(:) = 0

return
#if defined(use_mpi)
112 call mpi_abort(mpi_comm_world,1,ierr)
#else
112 stop 'aborting due to errors reading restart file.'
#endif

end subroutine prep_coord

!-----------------------------------------------------------------------

!sort out heavy atoms in restrained shell. use protein center to calculate distance.
!uses coordinates from topology unless 'implicit_rstr_from_file' is specified.
subroutine make_shell
! *** local variables
        integer                                         ::      i,ig,i3
        real(8)                                         ::      rin2,r2

        nshellats = 0
        rin2  = rexcl_i**2

        shell(:) = .false.

        do ig=1,ncgp_solute
       if (.not. excl(cgp(ig)%iswitch) .and. heavy(cgp(ig)%iswitch)) then 
                i3 = 3*cgp(ig)%iswitch-3
                r2 = ( xtop(i3+1) - xpcent(1) )**2 &
                        +( xtop(i3+2) - xpcent(2) )**2 &
                        +( xtop(i3+3) - xpcent(3) )**2

                if(r2 > rin2) then
                        do i=cgp(ig)%first, cgp(ig)%last
        nshellats = nshellats + 1 
                                shell(cgpatom(i)) = .true.
                        end do
                end if
   end if
        end do
        write(*,105) nshellats, rexcl_i, rexcl_o
105     format('found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,')')
end subroutine make_shell

!------------------------------------------------------------------------

!sort out heavy atoms in restrained shell. use protein center to calculate distance.
!use coordinates from topology unless 'implicit_rstr_from_file' is specified
subroutine make_shell2
! *** local variables
        integer                                         ::      i,ig,i3,k
        real(8)                                         ::      rout2,rin2,r2
        real(8), allocatable            ::      cgp_cent(:,:)
        nshellats = 0
        rin2  = rexcl_i**2

        shell(:) = .false.

        allocate(cgp_cent(3,ncgp+nwat))

        cgp_cent(:,:) = 0.

        do ig=1,ncgp_solute
    if (.not. excl(cgp(ig)%iswitch) .and. heavy(cgp(ig)%iswitch)) then
                do i = cgp(ig)%first,cgp(ig)%last
                        i3 = cgpatom(i)*3
                        cgp_cent(:,ig) = cgp_cent(:,ig) + xtop(i3-2:i3)
                end do
        cgp_cent(:,ig) = cgp_cent(:,ig)/real(cgp(ig)%last - cgp(ig)%first +1)
                r2 = dot_product(cgp_cent(:,ig)-xpcent(:),cgp_cent(:,ig)-xpcent(:))

                if ( r2 .gt. rin2 ) then
                        do i=cgp(ig)%first, cgp(ig)%last
                                nshellats = nshellats + 1
                                shell(cgpatom(i)) = .true.
                        end do
                end if
   end if
        end do

        deallocate(cgp_cent) 
        write(*,105) nshellats, rexcl_i, rexcl_o
105     format('found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,'�)')
end subroutine make_shell2

!-----------------------------------------------------------------------

subroutine init_trj
!locals
integer                                         ::      trj_atoms

!initialise trajectory atom mask
if(itrj_cycle > 0) then
        call trj_initialize(frames=nsteps/itrj_cycle, steps_before=itrj_cycle,&
                interval=itrj_cycle, steps=nsteps,      degf=ndegfree, &
                topfile=top_file)

        trj_atoms = trj_commit_mask()
        write(*,100) trj_atoms
        if(.not. trj_create(trj_file)) then
                call die('failure to open trajectory file')
        end if
end if

100     format('coordinates for',i6,' atoms will be written to the trajectory.')
end subroutine init_trj

!-----------------------------------------------------------------------
subroutine prep_sim
! local variables
integer                                         :: i, j, ig, istate

if (nodeid .eq. 0) then 
        write(*,*)
        call centered_heading('initialising dynamics', '-')
end if

! set parameters (bonds, angles, charges,...) & restraints for water   
if(nwat > 0) then
        select case (solvent_type)
        case (solvent_spc, solvent_3atom)
                crg_ow = crg(nat_solute+1)
                crg_hw = -crg_ow/2.0
        case(solvent_general)
                !add appropriate code here
                call die('topology contains mixed or non-3-atomic solvent. this feature is not implemented yet.')
        end select

        if( .not. use_pbc ) then
                call wat_sphere
                if (wpol_restr) call wat_shells

        else !compute charges of the system for box case 
        !(done in subroutine wat_sphere for sphere case)
                !calc. total charge of non-q-atoms
                crgtot = 0.0
                do i = 1, nat_solute
                        if ( iqatom(i)==0 ) crgtot = crgtot + crg(i)
                end do
                write (*,60) crgtot
60 format ('total charge of non-q atoms             = ',f10.2)


        !calc effective charge of whole system at this lambda
                crgqtot = 0.0
                do i = 1, nqat
                        do istate = 1, nstates
                        crgtot = crgtot + qcrg(i,istate)*eq(istate)%lambda
                        crgqtot = crgqtot + qcrg(i,istate)*eq(istate)%lambda
                        end do
                end do

                write (*,70) crgtot

70 format ('total charge of system                  = ',f10.2)

        end if
end if
! set the charge group membership for every topology atom only if using lrf or pbc
if(use_lrf .or. use_pbc) then
        call allocate_lrf_arrays

        do ig = 1, ncgp
                do i = cgp(ig)%first, cgp(ig)%last
                        iwhich_cgp(cgpatom(i)) = ig
                end do
        end do
end if

!       prepare an array of inverse masses
winv(:) = 1./iaclib(iac(:))%mass


if(use_pbc .and. control_box) then
        boxlength(:) = new_boxl(:)
        if ( put_solute_back_in_box .or. put_solvent_back_in_box ) then !only put back in box if either solute or solvent should be put back (qdyn input option)
                call put_back_in_box
        end if
        write(*,'(a)') 'boxsize changed. equilibration may be needed'
end if

if( use_pbc ) then 
        !compute masses of all molecules
        allocate( mol_mass(1:nmol) )
        mol_mass(:) = 0.0

        do i = 1,nmol-1 !all molecules but the last
                do j = istart_mol(i), istart_mol(i+1)-1 !all atoms of molecule
                        mol_mass(i) = mol_mass(i) + iaclib(iac(j))%mass
                end do
        end do

        do j = istart_mol(nmol), natom !last molecule
                mol_mass(nmol) = mol_mass(nmol) + iaclib(iac(j))%mass
        end do

        mol_mass(:) = 1./mol_mass(:)

        !prepare array of masses
        allocate( mass(1:natom) )
        mass(:) = 1.0/winv(:)

end if


!scale charges by sqrt(coulomb_constant) 
crg(:) = crg(:) * sqrt(coulomb_constant)
crg_ow = crg_ow * sqrt(coulomb_constant)
crg_hw = crg_hw * sqrt(coulomb_constant)

if(nqat > 0) then
        qcrg(:,:) = qcrg(:,:) * sqrt(coulomb_constant)
end if
end subroutine prep_sim

!-----------------------------------------------------------------------

subroutine qangle (istate)
! arguments
integer                                         :: istate

! local variables
integer                                         :: ia,i,j,k,ic,i3,j3,k3,im,icoupl,ib
real(8)                                         :: bji,bjk,scp,ang,da,ae,dv,gamma
real(8)                                         :: rji(3),rjk(3),f1,di(3),dk(3)


do ia = 1, nqangle


ic = qang(ia)%cod(istate)
 !skip if angle not present (code 0)
if ( ic > 0 ) then

gamma = 1.0
icoupl = 0

do im = 1, nang_coupl
   if ( iang_coupl(1,im) .eq. ia ) then
      icoupl = im
      ib     = iang_coupl(2,im)
      gamma = emorsed(ib)
          !couple improper to bond breaking not making
          if ( iang_coupl(3,im) .eq. 1) gamma = 1.0_8 - gamma
          exit
   end if
end do

i  = qang(ia)%i
j  = qang(ia)%j
k  = qang(ia)%k
i3 = 3*i-3
j3 = 3*j-3
k3 = 3*k-3

rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)
bji = sqrt ( rji(1)**2 + rji(2)**2 + rji(3)**2 )
bjk = sqrt ( rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
scp = scp / (bji*bjk)
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
ang = acos(scp)
da = ang - qanglib(ic)%ang0
ae = 0.5*qanglib(ic)%fk*da**2
eq(istate)%q%angle = eq(istate)%q%angle + ae*gamma

dv = gamma*qanglib(ic)%fk*da*eq(istate)%lambda
f1 = sin ( ang )
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
di(1) = f1 * ( rjk(1)/(bji*bjk) - scp*rji(1)/bji**2 )
di(2) = f1 * ( rjk(2)/(bji*bjk) - scp*rji(2)/bji**2 )
di(3) = f1 * ( rjk(3)/(bji*bjk) - scp*rji(3)/bji**2 )
dk(1) = f1 * ( rji(1)/(bji*bjk) - scp*rjk(1)/bjk**2 )
dk(2) = f1 * ( rji(2)/(bji*bjk) - scp*rjk(2)/bjk**2 )
dk(3) = f1 * ( rji(3)/(bji*bjk) - scp*rjk(3)/bjk**2 )
d(i3+1) = d(i3+1) + dv*di(1)
d(i3+2) = d(i3+2) + dv*di(2)
d(i3+3) = d(i3+3) + dv*di(3)
d(k3+1) = d(k3+1) + dv*dk(1)
d(k3+2) = d(k3+2) + dv*dk(2)
d(k3+3) = d(k3+3) + dv*dk(3)
d(j3+1) = d(j3+1) - dv*( di(1) + dk(1) )
d(j3+2) = d(j3+2) - dv*( di(2) + dk(2) )
d(j3+3) = d(j3+3) - dv*( di(3) + dk(3) )

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3

   d(i3+1) = d(i3+1) + dmorse_i(1,ib)*ae
   d(i3+2) = d(i3+2) + dmorse_i(2,ib)*ae
   d(i3+3) = d(i3+3) + dmorse_i(3,ib)*ae
   d(j3+1) = d(j3+1) + dmorse_j(1,ib)*ae
   d(j3+2) = d(j3+2) + dmorse_j(2,ib)*ae
   d(j3+3) = d(j3+3) + dmorse_j(3,ib)*ae

end if
end if
end do
end subroutine qangle

!-----------------------------------------------------------------------

subroutine qurey_bradley (istate)
! arguments
integer                                         :: istate



! local variables
integer                                         ::      ia,i,j,k,ic,i3,j3,k3,im,icoupl,ib
real(8)                                         ::      gamma
real(8)                                         ::      rik(3), dik, du, ru, eurey

do ia = 1, nqangle
        ic = qang(ia)%cod(istate)
        !skip if angle not present (code 0)
        if ( ic == 0  .or. qanglib(ic)%ureyfk == 0.) cycle

gamma = 1.0
icoupl = 0

do im = 1, nang_coupl
   if ( iang_coupl(1,im) .eq. ia ) then
      icoupl = im
      ib     = iang_coupl(2,im)
      gamma = emorsed(ib)
          !couple improper to bond breaking not making
          if ( iang_coupl(3,im) .eq. 1) gamma = 1 - gamma
   end if
end do

i  = qang(ia)%i
j  = qang(ia)%j
k  = qang(ia)%k
i3 = 3*i-3
j3 = 3*j-3
k3 = 3*k-3
        rik(1) = x(k3+1) - x(i3+1)
        rik(2) = x(k3+2) - x(i3+2)
        rik(3) = x(k3+3) - x(i3+3)
        dik = sqrt(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
        ru = dik - qanglib(ic)%ureyr0
        eurey = qanglib(ic)%ureyfk*ru**2
        eq(istate)%q%angle = eq(istate)%q%angle + eurey*gamma
        du = gamma*2*(qanglib(ic)%ureyfk*ru/dik)*eq(istate)%lambda


if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3

   d(i3+1) = d(i3+1) + dmorse_i(1,ib)*eurey
   d(i3+2) = d(i3+2) + dmorse_i(2,ib)*eurey
   d(i3+3) = d(i3+3) + dmorse_i(3,ib)*eurey


   d(j3+1) = d(j3+1) + dmorse_j(1,ib)*eurey
   d(j3+2) = d(j3+2) + dmorse_j(2,ib)*eurey
   d(j3+3) = d(j3+3) + dmorse_j(3,ib)*eurey


end if
end do
end subroutine qurey_bradley

!-----------------------------------------------------------------------

subroutine qbond (istate)
! arguments
integer                                         ::      istate

! local variables
integer                                         ::      ib,i,j,ic,i3,j3
real(8)                                         ::      b,db,be,dv,fexp
real(8)                                         ::      rij(3)

do ib = 1, nqbond

        ic = qbnd(ib)%cod(istate)
        !code 0 means bond not present
if ( ic > 0 ) then

i  = qbnd(ib)%i
j  = qbnd(ib)%j
i3 = 3*i-3
j3 = 3*j-3

rij(1) = x(j3+1) - x(i3+1)
rij(2) = x(j3+2) - x(i3+2)
rij(3) = x(j3+3) - x(i3+3)

b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
db = b - qbondlib(ic)%r0

        fexp = exp ( -qbondlib(ic)%amz*db ) 
        be = qbondlib(ic)%dmz*(fexp*fexp-2.*fexp) + 0.5*qbondlib(ic)%fk*db**2
        emorsed(ib) = -(fexp*fexp-2.*fexp)
eq(istate)%q%bond = eq(istate)%q%bond + be
dv = (2.*qbondlib(ic)%dmz*qbondlib(ic)%amz*(fexp-fexp*fexp) + qbondlib(ic)%fk*db)*eq(istate)%lambda/b

d(i3+1) = d(i3+1) - dv*rij(1)
d(i3+2) = d(i3+2) - dv*rij(2)
d(i3+3) = d(i3+3) - dv*rij(3)
d(j3+1) = d(j3+1) + dv*rij(1)
d(j3+2) = d(j3+2) + dv*rij(2)
d(j3+3) = d(j3+3) + dv*rij(3)

!force scaling factor to be 1 when distance is smaller than r0
 if ( db > 0 ) then
        emorsed(ib) = -(fexp*fexp-2.*fexp)
        dmorse_i(:,ib) = +2.*qbondlib(ic)%amz*(fexp-fexp*fexp)*eq(istate)%lambda/b*rij(:)
        dmorse_j(:,ib) = -2.*qbondlib(ic)%amz*(fexp-fexp*fexp)*eq(istate)%lambda/b*rij(:)
 else
        emorsed(ib) = 1
        dmorse_i(:,ib) = 0
        dmorse_j(:,ib) = 0
 end if

end if
end do
end subroutine qbond

!-----------------------------------------------------------------------

subroutine qimproper (istate)
! arguments
integer                                         ::      istate

! local variables
integer                                         ::      i,j,k,l,ip,ic,i3,j3,k3,l3
integer                                         ::      icoupl,im,ib
real(8)                                         ::      bj,bk,scp,phi,sgn,pe,dv,arg,f1,gamma
real(8)                                         ::      rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)                                         ::      rki(3),rlj(3),dp(12),di(3),dl(3)

do ip = 1,nqimp

ic = qimpcod(ip,istate)

if ( ic > 0 ) then

gamma = 1.0
icoupl = 0

do im = 1, nimp_coupl
   if ( iimp_coupl(1,im) .eq. ip ) then
      icoupl = im
      ib     = iimp_coupl(2,im)
      gamma = emorsed(ib)
          !couple improper to bond breaking not making
          if ( iimp_coupl(3,im) .eq. 1) gamma = 1 - gamma
   end if 
end do



i  = iqimp(ip)
j  = jqimp(ip)
k  = kqimp(ip)
l  = lqimp(ip)

i3=i*3-3
j3=j*3-3
k3=k*3-3
l3=l*3-3
rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)
rkl(1) = x(l3+1) - x(k3+1)
rkl(2) = x(l3+2) - x(k3+2)
rkl(3) = x(l3+3) - x(k3+3)
rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)


rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)
bj = sqrt ( rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk = sqrt ( rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))/(bj*bk)
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
phi = acos ( scp )
sgn =  rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
     +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1))
if ( sgn .lt. 0 ) phi = -phi

! ---       energy

arg = phi - qimp0(ic)
arg = arg - 2.*pi*nint(arg/(2.*pi))
dv  = qfkimp(ic)*arg
pe  = 0.5*dv*arg
eq(istate)%q%improper = eq(istate)%q%improper + pe*gamma
dv = dv*gamma*eq(istate)%lambda

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
di(1) = f1 * ( rnk(1)/(bj*bk) - scp*rnj(1)/bj**2 )
di(2) = f1 * ( rnk(2)/(bj*bk) - scp*rnj(2)/bj**2 )
di(3) = f1 * ( rnk(3)/(bj*bk) - scp*rnj(3)/bj**2 )
dl(1) = f1 * ( rnj(1)/(bj*bk) - scp*rnk(1)/bk**2 )
dl(2) = f1 * ( rnj(2)/(bj*bk) - scp*rnk(2)/bk**2 )
dl(3) = f1 * ( rnj(3)/(bj*bk) - scp*rnk(3)/bk**2 )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)
rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)


dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(i3+1) = d(i3+1) + dv*dp(1)
d(i3+2) = d(i3+2) + dv*dp(2)
d(i3+3) = d(i3+3) + dv*dp(3)
d(j3+1) = d(j3+1) + dv*dp(4)
d(j3+2) = d(j3+2) + dv*dp(5)
d(j3+3) = d(j3+3) + dv*dp(6)
d(k3+1) = d(k3+1) + dv*dp(7)
d(k3+2) = d(k3+2) + dv*dp(8)
d(k3+3) = d(k3+3) + dv*dp(9)
d(l3+1) = d(l3+1) + dv*dp(10)
d(l3+2) = d(l3+2) + dv*dp(11)
d(l3+3) = d(l3+3) + dv*dp(12)

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3

   d(i3+1) = d(i3+1) + dmorse_i(1,ib)*pe
   d(i3+2) = d(i3+2) + dmorse_i(2,ib)*pe
   d(i3+3) = d(i3+3) + dmorse_i(3,ib)*pe
   d(j3+1) = d(j3+1) + dmorse_j(1,ib)*pe
   d(j3+2) = d(j3+2) + dmorse_j(2,ib)*pe
   d(j3+3) = d(j3+3) + dmorse_j(3,ib)*pe

end if
end if
end do
end subroutine qimproper

!-----------------------------------------------------------------------

subroutine qtorsion (istate)
! arguments
integer                                         ::      istate

! local variables
integer                                         ::      i,j,k,l,ip,ic,i3,j3,k3,l3
integer                                         ::      icoupl,im,ib
real(8)                                         ::      bj,bk,scp,phi,sgn,pe,dv,arg,f1,gamma
real(8)                                         ::      rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8)                                         ::      rki(3),rlj(3),dp(12),di(3),dl(3)

do ip = 1,nqtor

ic = qtorcod(ip,istate)

if ( ic > 0 ) then

gamma = 1.0
icoupl = 0

do im = 1, ntor_coupl
   if ( itor_coupl(1,im) .eq. ip ) then
      icoupl = im
      ib     = itor_coupl(2,im)
      gamma = emorsed(ib)
          !couple improper to bond breaking not making
          if ( itor_coupl(3,im) .eq. 1) gamma = 1 - gamma
   end if
end do

i  = iqtor(ip)
j  = jqtor(ip)
k  = kqtor(ip)
l  = lqtor(ip)

i3=i*3-3
j3=j*3-3
k3=k*3-3
l3=l*3-3
rji(1) = x(i3+1) - x(j3+1)
rji(2) = x(i3+2) - x(j3+2)
rji(3) = x(i3+3) - x(j3+3)
rjk(1) = x(k3+1) - x(j3+1)
rjk(2) = x(k3+2) - x(j3+2)
rjk(3) = x(k3+3) - x(j3+3)
rkl(1) = x(l3+1) - x(k3+1)
rkl(2) = x(l3+2) - x(k3+2)
rkl(3) = x(l3+3) - x(k3+3)
rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)
bj = sqrt ( rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk = sqrt ( rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))/(bj*bk)
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
phi = acos ( scp )
sgn =  rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
     +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &


    +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1))
if ( sgn .lt. 0 ) phi = -phi




! ---       energy

arg = qrmult(ic)*phi-qdeltor(ic)
pe  = qfktor(ic)*(1.0+cos(arg))
eq(istate)%q%torsion = eq(istate)%q%torsion + pe*gamma
dv = -qrmult(ic)*qfktor(ic)*sin(arg)*gamma*eq(istate)%lambda

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
di(1) = f1 * ( rnk(1)/(bj*bk) - scp*rnj(1)/bj**2 )
di(2) = f1 * ( rnk(2)/(bj*bk) - scp*rnj(2)/bj**2 )
di(3) = f1 * ( rnk(3)/(bj*bk) - scp*rnj(3)/bj**2 )
dl(1) = f1 * ( rnj(1)/(bj*bk) - scp*rnk(1)/bk**2 )
dl(2) = f1 * ( rnj(2)/(bj*bk) - scp*rnk(2)/bk**2 )
dl(3) = f1 * ( rnj(3)/(bj*bk) - scp*rnk(3)/bk**2 )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)


rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)
dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(i3+1) = d(i3+1) + dv*dp(1)
d(i3+2) = d(i3+2) + dv*dp(2)
d(i3+3) = d(i3+3) + dv*dp(3)
d(j3+1) = d(j3+1) + dv*dp(4)
d(j3+2) = d(j3+2) + dv*dp(5)
d(j3+3) = d(j3+3) + dv*dp(6)
d(k3+1) = d(k3+1) + dv*dp(7)
d(k3+2) = d(k3+2) + dv*dp(8)
d(k3+3) = d(k3+3) + dv*dp(9)
d(l3+1) = d(l3+1) + dv*dp(10)
d(l3+2) = d(l3+2) + dv*dp(11)
d(l3+3) = d(l3+3) + dv*dp(12)

if ( icoupl .ne. 0 ) then

   i  = qbnd(ib)%i
   j  = qbnd(ib)%j
   i3 = 3*i-3
   j3 = 3*j-3



  d(i3+1) = d(i3+1) + dmorse_i(1,ib)*pe
   d(i3+2) = d(i3+2) + dmorse_i(2,ib)*pe
   d(i3+3) = d(i3+3) + dmorse_i(3,ib)*pe
   d(j3+1) = d(j3+1) + dmorse_j(1,ib)*pe
   d(j3+2) = d(j3+2) + dmorse_j(2,ib)*pe
   d(j3+3) = d(j3+3) + dmorse_j(3,ib)*pe

end if
end if
end do
end subroutine qtorsion

!-----------------------------------------------------------------------


real(8) function randm (ig)
! arguments
integer                                 ::      ig

! local variables
integer, parameter              ::      m = 100000000
integer, parameter              ::  m1 = 10000
integer, parameter              ::      mult=31415821
integer                                 ::      irandh,irandl,multh,multl
real(8)                                 ::      r
integer, save                           ::      irand = 0
integer, save                           ::  new = 0

if (new .eq. 0) then
new = 1
irand = mod (iabs(ig),m)
end if

! --- multiply irand by mult, but take into account that overflow must
! --- be discarded, and do not generate an error.
irandh = irand / m1
irandl = mod(irand, m1)
multh = mult / m1
multl = mod(mult, m1)

irand = mod(irandh*multl + irandl*multh, m1) * m1 + irandl*multl
irand = mod(irand + 1, m)

! --- convert irand to a real random number between 0 and 1.
r = real(irand / 10) * 10 / real(m)
if ((r .le. 0.e0) .or. (r .gt. 1.e0)) r = 0.e0
randm = r
ig = irand

end function randm

!-----------------------------------------------------------------------

integer function shake(xx, x)
!arguments
real(8)                                         ::      xx(:), x(:)
!       returns no. of iterations

! *** local variables
integer                                         ::      i,j,i3,j3,mol,ic,nits
real(8)                                         ::      xij2,diff,corr,scp,xxij2
real(8)                                         ::      xij(3),xxij(3)
#if defined (profiling)
real(8)                                         :: start_loop_time
start_loop_time = rtime()
#endif

! reset niter
shake = 0

do mol=1,shake_molecules
        ! for every molecule:
        ! reset nits (iterations per molecule)
        nits = 0
        ! reset iready for every constraint
        shake_mol(mol)%bond(:)%ready = .false.
        do !iteration loop
                do ic=1,shake_mol(mol)%nconstraints
                        ! for every constraint:

                        if (.not. shake_mol(mol)%bond(ic)%ready) then
                                ! repeat until done:

                                i = shake_mol(mol)%bond(ic)%i
                                j = shake_mol(mol)%bond(ic)%j
                                i3 = i*3-3
                                j3 = j*3-3
                                xij(1)  = x(i3+1) - x(j3+1)
                                xij(2)  = x(i3+2) - x(j3+2)
                                xij(3)  = x(i3+3) - x(j3+3)
                                xij2    = xij(1)**2+xij(2)**2+xij(3)**2
                                diff    = shake_mol(mol)%bond(ic)%dist2 - xij2
                                if(abs(diff) < shake_tol*shake_mol(mol)%bond(ic)%dist2) then
                                        shake_mol(mol)%bond(ic)%ready = .true. ! in range
                                end if
                                xxij(1) = xx(i3+1) - xx(j3+1)
                                xxij(2) = xx(i3+2) - xx(j3+2)
                                xxij(3) = xx(i3+3) - xx(j3+3)
                                scp = xij(1)*xxij(1)+xij(2)*xxij(2)+xij(3)*xxij(3)
                                corr = diff/(2.*scp*(winv(i)+winv(j)))

                                x(i3+1) = x(i3+1)+xxij(1)*corr*winv(i)
                                x(i3+2) = x(i3+2)+xxij(2)*corr*winv(i)
                                x(i3+3) = x(i3+3)+xxij(3)*corr*winv(i)
                                x(j3+1) = x(j3+1)-xxij(1)*corr*winv(j)
                                x(j3+2) = x(j3+2)-xxij(2)*corr*winv(j)
                                x(j3+3) = x(j3+3)-xxij(3)*corr*winv(j)
                        end if
                end do

            nits = nits+1

                ! see if every constraint is met
                if(all(shake_mol(mol)%bond(1:shake_mol(mol)%nconstraints)%ready)) then
                        exit !from iteration loop
                elseif(nits >= shake_max_iter) then
                        ! fail on too many iterations
                        do ic=1,shake_mol(mol)%nconstraints
                                if (.not. shake_mol(mol)%bond(ic)%ready) then
                                        ! repeat until done:

                                        i = shake_mol(mol)%bond(ic)%i
                                        j = shake_mol(mol)%bond(ic)%j
                                        i3 = i*3-3
                                        j3 = j*3-3
                                        xxij(1) = xx(i3+1) - xx(j3+1)
                                        xxij(2) = xx(i3+2) - xx(j3+2)
                                        xxij(3) = xx(i3+3) - xx(j3+3)
                                        xxij2   = xxij(1)**2+xxij(2)**2+xxij(3)**2
                                        write (*,100) i,j,sqrt(xxij2),&
                                                sqrt(shake_mol(mol)%bond(ic)%dist2)
                                end if
                        end do
                        call die('shake failure')
                end if
100         format ('>>> shake failed, i,j,d,d0 = ',2i6,2f10.5)
        end do

        ! update niter
        shake = shake+nits
end do

! set niter to the average number of iterations per molecule
shake=shake/nmol
#if defined (profiling)
profile(7)%time = profile(7)%time + rtime() - start_loop_time
#endif

end function shake

!-----------------------------------------------------------------------

subroutine shrink_topology
!get rid of bonds and angles where all atoms are excluded
!or where the code has been set to 0 due to q-[bonds|angles|...]

!locals
integer                                         ::      i, removed

if(exclude_bonded) then
        call centered_heading &
                ('eliminating torsions & impropers for excluded atoms', '-')
end if

10      format('reduced number of ',a,t31,'from ',i8,' to ')
12      format(i8)

i = 1
removed = 0
do while(i <= nbonds)
        !if all atoms excluded
        if(bnd(i)%cod <= 0) then
                !bond code either 0 (bond redefined in fep file) 
                !or -1 (bond removed by shake)
                if(i <= nbonds_solute) then
                        bnd(i) = bnd(nbonds_solute)
                        bnd(nbonds_solute) = bnd(nbonds)
                        nbonds_solute = nbonds_solute - 1
                else
                        bnd(i) = bnd(nbonds)
                endif
                nbonds = nbonds - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do

i = 1
do while(i <= nangles)
        !if all atoms excluded
        if(ang(i)%cod == 0) then
                !move last angle to current position
                if(i <= nangles_solute) then
                        ang(i) = ang(nangles_solute)
                        ang(nangles_solute) = ang(nangles)
                        nangles_solute = nangles_solute - 1
                else
                        ang(i) = ang(nangles)
                endif
                nangles = nangles - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do

if(exclude_bonded) write(*,10, advance='no') 'torsions', ntors
i = 1
do while(i <= ntors)
        !if all atoms excluded
        if((exclude_bonded .and. excl(tor(i)%i) .and. excl(tor(i)%j) &
                .and. excl(tor(i)%k) .and. excl(tor(i)%l)) &
                .or. tor(i)%cod == 0) then
                !move last bond to current position
                if(i <= ntors_solute) then
                        tor(i) = tor(ntors_solute)
                        tor(ntors_solute) = tor(ntors)
                        ntors_solute = ntors_solute - 1
                else
                        tor(i) = tor(ntors)
                endif
                ntors = ntors - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do
if(exclude_bonded) write(*, 12) ntors

if(exclude_bonded) write(*,10, advance='no') 'impropers', nimps
i = 1
do while(i <= nimps)
        !if all atoms excluded
        if(exclude_bonded .and. excl(imp(i)%i) .and. excl(imp(i)%j) &
                .and. excl(imp(i)%k) .and. excl(imp(i)%l) &
                .or. imp(i)%cod == 0) then
                if(i <= nimps_solute) then
                        imp(i) = imp(nimps_solute)
                        imp(nimps_solute) = imp(nimps)
                        nimps_solute = nimps_solute - 1
                else
                        imp(i) = imp(nimps)
                endif
                nimps = nimps - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do
if(exclude_bonded) write(*, 12) nimps

end subroutine shrink_topology

!--------------------------------------------------------------------


subroutine stop_cm_translation
! local variables
integer                                         ::      i,j,k
real(8)                                         ::      rmass,totmass
real(8)                                         ::      vcm(3)

! calculate totmass and vcm
totmass = 0.0
do i=1,3
vcm(i) = 0.0
end do
do i=1,natom
rmass = iaclib(iac(i))%mass
totmass=totmass+rmass
do j=1,3
k=(i-1)*3+j
vcm(j)=vcm(j)+rmass*v(k)
end do
end do

! scale vcm
do j=1,3
vcm(j)=vcm(j)/totmass
end do

! update v
do i=1,natom
do j=1,3
k=(i-1)*3+j
v(k)=v(k)-vcm(j)
end do
end do
end subroutine stop_cm_translation

!-----------------------------------------------------------------------
subroutine topology
! local variables
integer                                 ::      nat3
integer                                 ::      i
real(8)                                 ::      box_min, vtemp, vtemp1

!
! read topology
!
! will init:
!  natom
!  lots of stuff from topo_load
!  nwat
!  anglib, torlib, implib (conversion)
!  ljcod
!  [iaclib%bvdw] (conversion)

if(.not. topo_load(top_file, require_version=4.15)) then
        call die('failed to load topology.')
end if
natom = nat_pro

nwat = (natom - nat_solute) / 3
!add extra element to molecule start atom array to keep track of last atom
istart_mol(nmol+1) = nat_pro + 1

! abort if no atoms
if (natom .eq. 0) call die('zero particles to simulate')

! convert libraries from degrees to radians
anglib(1:nangcod)%ang0 = deg2rad*anglib(1:nangcod)%ang0
torlib(1:ntorcod)%paths = 1.0/torlib(1:ntorcod)%paths
torlib(1:ntorcod)%deltor = deg2rad*torlib(1:ntorcod)%deltor
implib(1:nimpcod)%imp0 = deg2rad*implib(1:nimpcod)%imp0

ljcod(:,:) = 1
do i=1,nlj2
ljcod(lj2(i)%i, lj2(i)%j) = 2
ljcod(lj2(i)%j, lj2(i)%i) = 2
end do


!       if arithmetic combination rule (ivdw_rule=2) take sqrt(epsilon) now
if ( ivdw_rule .eq. 2 ) then
do i=1,natyps
  iaclib(i)%bvdw(:) = sqrt(abs(iaclib(i)%bvdw(:)))
end do
end if

!check if same boundary in topology and qdyn input file
if(  ((.not. use_pbc) .and.  box) .or. ( use_pbc .and. (.not. box)) ) then
call die('must have same boundary (sphere or box) in topology and input file')
end if

if( use_pbc ) then
if( any(boxlength(:) == 0.0 ) ) then
        inv_boxl(:) = 0.0
else
        inv_boxl(:) = 1.0/boxlength(:)
end if

!check cut-offs if periodic box used
box_min = min( boxlength(1), boxlength(2), boxlength(3) )

!solute-solute cut-off radii
if( .not. (box_min .gt. rcpp*2) ) then
        call die('solute-solute cut-off radii too large')

!solvent-solvent
else if( .not. (box_min .gt. rcww*2) ) then
        call die('solvent-solvent cut-off radii too large')

!solute-solvent
else if( .not. (box_min .gt. rcpw*2) ) then
        call die('solute-solvent cut-off radii too large')

!q-atom
else if( .not. (box_min .gt. rcq*2) ) then
        call die('q-atom cut-off radii too large')
!lrf
else if( .not. (box_min .gt. rclrf*2) ) then
        call die('lrf cut-off radii too large')
end if
end if
end subroutine topology

!-----------------------------------------------------------------------

real(8) function torsion(istart, iend)
!arguments
integer                                         ::      istart, iend

! local variables
integer                                         ::      ip
real(8)                                         ::      scp,phi,dv,arg,f1
real(8)                                         ::      bjinv, bkinv, bj2inv, bk2inv
real(8), save                                   ::      rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
real(8), save                                   ::      rki(3),rlj(3),dp(12),di(3),dl(3)
type(tor_type), pointer         ::  t
type(torlib_type), pointer      ::      lib

! global variables used:
!  tor, torlib, x, d

! calculate the total energy of all torsion angles
! updates d

torsion = 0.

do ip = iend, istart,-1
        t => tor(ip)
        lib => torlib(t%cod)
rji(1) = x(t%i*3-2) - x(t%j*3-2)
rji(2) = x(t%i*3-1) - x(t%j*3-1)
rji(3) = x(t%i*3-0) - x(t%j*3-0)
rjk(1) = x(t%k*3-2) - x(t%j*3-2)
rjk(2) = x(t%k*3-1) - x(t%j*3-1)
rjk(3) = x(t%k*3-0) - x(t%j*3-0)
rkl(1) = x(t%l*3-2) - x(t%k*3-2)
rkl(2) = x(t%l*3-1) - x(t%k*3-1)
rkl(3) = x(t%l*3-0) - x(t%k*3-0)

rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

        bj2inv = 1./(rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
bk2inv = 1./(rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

        ! calculate scp and phi
scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
if ( scp .gt.  1.0 ) then
                scp =  1.0
                phi = acos (1.0) ! const
        else if ( scp .lt. -1.0 ) then
                scp = -1.0
                phi = acos (-1.0) ! const
        else
    phi = acos ( scp )
        end if
if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
                +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
    +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) .lt. 0) then
                phi = -phi
        end if

! ---       energy

arg = lib%rmult*phi-lib%deltor
torsion = torsion + lib%fk*(1.0+cos(arg))*lib%paths   !lib%paths is previously inverted 
dv = -lib%rmult*lib%fk*sin(arg)*lib%paths

! ---       forces

f1 = sin ( phi ) 
if ( abs(f1) .lt. 1.e-12 ) f1 = 1.e-12
f1 =  -1.0 / f1
di(1) = f1 * ( rnk(1)*(bjinv*bkinv) - scp*rnj(1)*bj2inv )
di(2) = f1 * ( rnk(2)*(bjinv*bkinv) - scp*rnj(2)*bj2inv )
di(3) = f1 * ( rnk(3)*(bjinv*bkinv) - scp*rnj(3)*bj2inv )
dl(1) = f1 * ( rnj(1)*(bjinv*bkinv) - scp*rnk(1)*bk2inv )
dl(2) = f1 * ( rnj(2)*(bjinv*bkinv) - scp*rnk(2)*bk2inv )
dl(3) = f1 * ( rnj(3)*(bjinv*bkinv) - scp*rnk(3)*bk2inv )

rki(1) =  rji(1) - rjk(1)
rki(2) =  rji(2) - rjk(2)
rki(3) =  rji(3) - rjk(3)
rlj(1) = -rjk(1) - rkl(1)
rlj(2) = -rjk(2) - rkl(2)
rlj(3) = -rjk(3) - rkl(3)

dp(1)  = rjk(2)*di(3) - rjk(3)*di(2)
dp(2)  = rjk(3)*di(1) - rjk(1)*di(3)
dp(3)  = rjk(1)*di(2) - rjk(2)*di(1)
dp(4)  = rki(2)*di(3)-rki(3)*di(2)+rkl(2)*dl(3)-rkl(3)*dl(2)
dp(5)  = rki(3)*di(1)-rki(1)*di(3)+rkl(3)*dl(1)-rkl(1)*dl(3)
dp(6)  = rki(1)*di(2)-rki(2)*di(1)+rkl(1)*dl(2)-rkl(2)*dl(1)
dp(7)  = rlj(2)*dl(3)-rlj(3)*dl(2)-rji(2)*di(3)+rji(3)*di(2)
dp(8)  = rlj(3)*dl(1)-rlj(1)*dl(3)-rji(3)*di(1)+rji(1)*di(3)
dp(9)  = rlj(1)*dl(2)-rlj(2)*dl(1)-rji(1)*di(2)+rji(2)*di(1)
dp(10) = rjk(2)*dl(3) - rjk(3)*dl(2)
dp(11) = rjk(3)*dl(1) - rjk(1)*dl(3)
dp(12) = rjk(1)*dl(2) - rjk(2)*dl(1)

d(t%i*3-2) = d(t%i*3-2) + dv*dp(1)
d(t%i*3-1) = d(t%i*3-1) + dv*dp(2)
d(t%i*3-0) = d(t%i*3-0) + dv*dp(3)
d(t%j*3-2) = d(t%j*3-2) + dv*dp(4)
d(t%j*3-1) = d(t%j*3-1) + dv*dp(5)
d(t%j*3-0) = d(t%j*3-0) + dv*dp(6)
d(t%k*3-2) = d(t%k*3-2) + dv*dp(7)
d(t%k*3-1) = d(t%k*3-1) + dv*dp(8)
d(t%k*3-0) = d(t%k*3-0) + dv*dp(9)
d(t%l*3-2) = d(t%l*3-2) + dv*dp(10)
d(t%l*3-1) = d(t%l*3-1) + dv*dp(11)
d(t%l*3-0) = d(t%l*3-0) + dv*dp(12)
end do

end function torsion

!-----------------------------------------------------------------------
subroutine restrain_solvent 
! local variables
integer                                         ::      iw,i,i3
real(8)                                         ::      b,db,erst,dv,fexp
real(8), save                                   ::      dr(3)
real(8)                                         ::      shift

! global variables used:
!  e, boltz, tfree, fk_wsphere, nwat, nat_pro, x, xwcent, rwat, dwmz, awmz, d

if(fk_wsphere /= 0.) then
        shift = sqrt (boltz*tfree/fk_wsphere)
else
        shift = 0.
end if
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        i3 = 3*i-3

        dr(1) = x(i3+1) - xwcent(1)
        dr(2) = x(i3+2) - xwcent(2)
        dr(3) = x(i3+3) - xwcent(3)
        b = sqrt ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
        db = b - (rwat - shift)

        ! calculate erst and dv
        if ( db > 0 ) then
                erst = 0.5 * fk_wsphere * db**2 - dwmz
                dv = fk_wsphere*db/b
        else
                if (b > 0.0) then
                  fexp = exp ( awmz*db )
                  erst = dwmz*(fexp*fexp-2.*fexp)
                  dv = -2.*dwmz*awmz*(fexp-fexp*fexp)/b
                else
                  dv = 0
                  erst = 0
                end if
        end if

        ! update energy and forces
        e%restraint%solvent_radial = e%restraint%solvent_radial + erst
        d(i3+1) = d(i3+1) + dv*dr(1)
        d(i3+2) = d(i3+2) + dv*dr(2)
        d(i3+3) = d(i3+3) + dv*dr(3)
end do
end subroutine restrain_solvent

!-----------------------------------------------------------------------

subroutine wat_sphere
! local variables
integer                                 :: i,i3,kr,isort,int_wat,istate
real(8)                                 :: rc,rnwat
real(8), save                           :: dr(3)
real(8)                                 :: crgexcl

!possibly override target sphere radius from topology
if(rwat_in > 0.) rwat = rwat_in


!calc. total charge of non-excluded non-q-atoms and excluded atoms
crgtot = 0.0
crgexcl = 0.0
do i = 1, nat_solute
        if ( .not. excl(i) ) then
                if ( iqatom(i)==0 ) then
                        crgtot = crgtot + crg(i)
                end if
        else
                crgexcl = crgexcl + crg(i)
        end if
end do
write (*,60) 'non-q atoms', crgtot
write (*,60) 'excluded atoms', crgexcl
60 format ('total charge of ',a,t41,'= ',f10.2)

!calc effective charge of simulation sphere at this lambda
crgqtot = 0.0
do i = 1, nqat
        do istate = 1, nstates
                crgtot = crgtot + qcrg(i,istate)*eq(istate)%lambda
                crgqtot = crgqtot + qcrg(i,istate)*eq(istate)%lambda
        end do
end do
write (*,70) crgtot
70 format ('total charge of system                  = ',f10.2)

if (.not. wpol_born) crgqtot = 0.0 !ignore total q if born corr. is off
if ( nwat .eq. 0 ) return


!       set default values for unspecified optional parameters
if(fk_wsphere == -1) then
        !
        ! to be replaced by function of rc giving appropriate default for any sphere
        !
        fk_wsphere = fk_wsphere_default
end if
if(fkwpol == -1) then
        !
        ! to be replaced by function of rc giving appropriate default for any sphere
        !
        fkwpol = fkwpol_default
end if
if(dwmz == -1) then !use magic function to get suitable dwmz
        dwmz = 0.26*exp(-0.19*(rwat-15.))+0.74
end if
if(awmz == -1) then !use magic for the reach of the morse potential
        awmz = 0.2/(1.+exp(0.4*(rwat-25.)))+0.3
end if

write (*,90) rwat, fk_wsphere, dwmz, awmz
90      format ('target water sphere radius              = ',f10.2,/,&
                'surface inward harmonic force constant  = ',f10.2,/,&
                'surface attraction well depth           = ',f10.2,/,&
                'surface attraction well width           = ',f10.2)
92      format ('water polarisation restraints           : ',a)
if(.not. wpol_restr) then
        write(*,92) 'off'
else if(wpol_born) then
        write(*,92) 'on, born correction enabled'
        write(*, 100) fkwpol
else 
        write(*,92) 'on, born correction disabled'
        write(*, 100) fkwpol
end if
100     format('radial polarisation force constant      = ',f10.2)

end subroutine wat_sphere

!-----------------------------------------------------------------------

subroutine wat_shells
! set up the shells for polarisation restraining

! local variables
real(8)                                         ::      rout, dr, ri, vshell, rshell, drs
integer                                         ::      is, n_insh


integer                                         ::      nwpolr_shell_restart, filestat
integer                                         ::      bndcodw, angcodw

!calc mu_w
!look up bond code for water
bndcodw = bnd(nbonds)%cod
angcodw = ang(nangles)%cod
!find charge of water o = charge of 1st solvent atom
crg_ow = crg(nat_solute + 1)
mu_w = -crg_ow*bondlib(bndcodw)%bnd0*cos(anglib(angcodw)%ang0/2)

! shell widths are drout, 2drout, 3drout
drs = wpolr_layer / drout !number of drouts 

! calc number of shells based on arithmetic series sum formula
nwpolr_shell = int(-0.5 + sqrt(2*drs + 0.25)) 
allocate(wshell(nwpolr_shell), stat=alloc_status)
call check_alloc('water polarisation shell array')

write(*, 100) nwpolr_shell
100     format(/,'setting up ', i1, ' water shells for polarisation restraints.')

if(restart) then !try to load theta_corr from restart file
        read(2, iostat=filestat) nwpolr_shell_restart
        if(filestat /= 0 .or. nwpolr_shell_restart /= nwpolr_shell) then
                write(*,102) 
                wshell(:)%theta_corr = 0.
        else
                backspace(2)
                read(2) nwpolr_shell_restart, wshell(:)%theta_corr
                write(*,103)
        end if
else
        wshell(:)%theta_corr = 0.
end if

102     format('>>> warning: failed to read polarisation restraint data from restart file.')
103     format('loaded polarisation restraint data from restart file.')

write(*,'(a)') 'shell #    outer radius    inner radius'
110     format(i7, 2f16.2)

rout = rwat
n_max_insh = 0
do is = 1, nwpolr_shell 
        wshell(is)%avtheta = 0
        wshell(is)%avn_insh = 0
        wshell(is)%rout = rout
    dr = drout*is
        ri = rout - dr
        wshell(is)%dr = dr
        vshell = rout**3 - ri**3
        n_insh = int(4 * pi/3 * vshell * rho_wat)
        if (n_insh > n_max_insh) n_max_insh = n_insh
rshell = (0.5*(rout**3+ri**3))**(1./3.)


        ! --- note below: 0.98750 = (1-1/epsilon) for water
wshell(is)%cstb = crgqtot*0.98750/(rho_wat*mu_w*4.*pi*rshell**2)
        write(*, 110) is, rout, ri
        rout = rout - dr
end do

n_max_insh = n_max_insh * 1.5 !take largest and add some extra
call allocate_watpol_arrays

end subroutine wat_shells

!-----------------------------------------------------------------------

subroutine watpol
! local variables
integer                                         :: iw,is,i,i3,il,jl,jw,imin,jmin
real(8)                                         :: dr,rw,rshell,rm,rc,scp
real(8)                                         :: tmin,arg,avtdum,dv,f0
real(8), save                                   :: f1(9),f2(3)
real(8), save                                   :: rmu(3),rcu(3)

! global variables used:
!  e, wshell, bndw0, deg2rad, angw0, nwat, theta, theta0, nat_pro, x, xwcent,
!  tdum, nwpolr_shell, list_sh, pi, nsort, istep, itdis_update, fkwpol, d

! reset wshell%n_insh
wshell(:)%n_insh = 0

! calculate theta(:), tdum(:), wshell%n_insh
do iw = 1, nwat
theta(iw)  = 0.0
theta0(iw) = 0.0

i  = nat_solute + iw*3-2
if(excl(i)) cycle ! skip excluded topology waters
i3 = i*3-3

rmu(1) = x(i3+4) + x(i3+7) - 2.*x(i3+1)   !water vector
rmu(2) = x(i3+5) + x(i3+8) - 2.*x(i3+2)
rmu(3) = x(i3+6) + x(i3+9) - 2.*x(i3+3)
rm = sqrt ( rmu(1)**2 + rmu(2)**2 + rmu(3)**2 )
rmu(1) = rmu(1)/rm
rmu(2) = rmu(2)/rm
rmu(3) = rmu(3)/rm

rcu(1) = x(i3+1) - xwcent(1)    !radial vector to ow
rcu(2) = x(i3+2) - xwcent(2) 
rcu(3) = x(i3+3) - xwcent(3) 
rc = sqrt ( rcu(1)**2 + rcu(2)**2 + rcu(3)**2 )
rcu(1) = rcu(1)/rc
rcu(2) = rcu(2)/rc
rcu(3) = rcu(3)/rc

scp = rmu(1)*rcu(1)+rmu(2)*rcu(2)+rmu(3)*rcu(3)   !calculate angle between water vector and radial vector
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
theta(iw) = acos( scp )
tdum(iw) = theta(iw)

if ( rc > wshell(nwpolr_shell)%rout-wshell(nwpolr_shell)%dr ) then
  do is = nwpolr_shell, 2, -1
        if(rc <= wshell(is)%rout) exit
  end do
  wshell(is)%n_insh = wshell(is)%n_insh + 1
  list_sh(wshell(is)%n_insh,is) = iw
end if
end do

! sort the waters according to theta
do is = 1, nwpolr_shell
imin = 0
do il = 1, wshell(is)%n_insh
tmin = 2.*pi
do jl = 1, wshell(is)%n_insh
jw = list_sh(jl,is)
if ( tdum(jw) .lt. tmin ) then
  jmin = jw
  tmin = theta(jw)
end if
end do
imin = imin+1
nsort(imin,is) = jmin
  tdum(jmin) = 99999.

end do

end do

! calculate energy and force
if ( istep .ne. 0 .and. mod(istep,itdis_update) .eq. 0) then
call centered_heading('water polarisation restraint data', '-')
write(*,'(a)') 'shell    <n>    <theta>    theta_0 theta_corr'
do is = 1, nwpolr_shell
wshell(is)%avtheta = wshell(is)%avtheta / real (itdis_update)
wshell(is)%avn_insh = wshell(is)%avn_insh / real (itdis_update)
wshell(is)%theta_corr = wshell(is)%theta_corr + wshell(is)%avtheta-acos(wshell(is)%cstb)
write (*,10) is,wshell(is)%avn_insh,wshell(is)%avtheta/deg2rad, &
     acos(wshell(is)%cstb)/deg2rad,wshell(is)%theta_corr/deg2rad
10        format(i5,1x,f6.1,3x,f8.3,3x,f8.3,3x,f8.3)
wshell(is)%avtheta = 0.0
wshell(is)%avn_insh = 0.0
end do
end if

do is = 1, nwpolr_shell
if(wshell(is)%n_insh == 0) cycle !skip empty shell
avtdum = 0.0
do il = 1, wshell(is)%n_insh
iw = nsort(il,is)
arg = 1. + (1. - 2.*real(il))/real(wshell(is)%n_insh)
theta0(il) = acos ( arg )
theta0(il) = theta0(il)-3.*sin(theta0(il))*wshell(is)%cstb/2.
if ( theta0(il) .lt. 0.0 ) theta0(il) = 0.0
if ( theta0(il) .gt. pi)   theta0(il) = pi

avtdum = avtdum + theta(iw)

e%restraint%water_pol = e%restraint%water_pol + 0.5*fkwpol* &
     (theta(iw)-theta0(il)+wshell(is)%theta_corr)**2

dv = fkwpol*(theta(iw)-theta0(il)+wshell(is)%theta_corr)

i  = nat_solute + iw*3-2
i3 = i*3-3

rmu(1) = x(i3+4) + x(i3+7) - 2.*x(i3+1)
rmu(2) = x(i3+5) + x(i3+8) - 2.*x(i3+2)
rmu(3) = x(i3+6) + x(i3+9) - 2.*x(i3+3)
rm = sqrt ( rmu(1)**2 + rmu(2)**2 + rmu(3)**2 )
rmu(1) = rmu(1)/rm
rmu(2) = rmu(2)/rm
rmu(3) = rmu(3)/rm

rcu(1) = x(i3+1) - xwcent(1)
rcu(2) = x(i3+2) - xwcent(2)
rcu(3) = x(i3+3) - xwcent(3)
rc = sqrt ( rcu(1)**2 + rcu(2)**2 + rcu(3)**2 )
rcu(1) = rcu(1)/rc
rcu(2) = rcu(2)/rc
rcu(3) = rcu(3)/rc



scp = rmu(1)*rcu(1)+rmu(2)*rcu(2)+rmu(3)*rcu(3)
if ( scp .gt.  1.0 ) scp =  1.0
if ( scp .lt. -1.0 ) scp = -1.0
f0 = sin ( acos(scp) )
if ( abs(f0) .lt. 1.e-12 ) f0 = 1.e-12
f0 = -1.0 / f0
f0 = dv*f0

f1(1) = -2.*(rcu(1)-rmu(1)*scp)/rm
f1(2) = -2.*(rcu(2)-rmu(2)*scp)/rm
f1(3) = -2.*(rcu(3)-rmu(3)*scp)/rm
f1(4) =     (rcu(1)-rmu(1)*scp)/rm
f1(5) =     (rcu(2)-rmu(2)*scp)/rm
f1(6) =     (rcu(3)-rmu(3)*scp)/rm
f1(7) =     (rcu(1)-rmu(1)*scp)/rm
f1(8) =     (rcu(2)-rmu(2)*scp)/rm
f1(9) =     (rcu(3)-rmu(3)*scp)/rm

f2(1) = ( rmu(1)-rcu(1)*scp)/rc
f2(2) = ( rmu(2)-rcu(2)*scp)/rc
f2(3) = ( rmu(3)-rcu(3)*scp)/rc

d(i3+1) = d(i3+1) + f0 * ( f1(1) + f2(1) )
d(i3+2) = d(i3+2) + f0 * ( f1(2) + f2(2) )
d(i3+3) = d(i3+3) + f0 * ( f1(3) + f2(3) )
d(i3+4) = d(i3+4) + f0 * ( f1(4) )
d(i3+5) = d(i3+5) + f0 * ( f1(5) )
d(i3+6) = d(i3+6) + f0 * ( f1(6) )
d(i3+7) = d(i3+7) + f0 * ( f1(7) )
d(i3+8) = d(i3+8) + f0 * ( f1(8) )
d(i3+9) = d(i3+9) + f0 * ( f1(9) )
end do

wshell(is)%avtheta = wshell(is)%avtheta + avtdum/real(wshell(is)%n_insh)
wshell(is)%avn_insh = wshell(is)%avn_insh + wshell(is)%n_insh
end do
end subroutine watpol

!----------------------------------------------------------------------------

subroutine write_out
! local variables
integer                                 ::      i,istate

! header line
if(istep >= nsteps) then
write(*,3) 'energy summary'
else
write(*,2) 'energy summary', istep
end if
2 format('======================= ',a15,' at step ',i6,' ========================')
3 format('=========================== final ',a15,' =============================')

! legend line
write(*,4) 'el', 'vdw' ,'bond', 'angle', 'torsion', 'improper'
4 format(16x, 6a10)

! row by row: solute, solvent, solute-solvent, lrf, q-atom
write(*,6) 'solute', e%pp%el, e%pp%vdw, e%p%bond, e%p%angle, e%p%torsion, e%p%improper
6 format(a,t17, 6f10.2)

if(nwat > 0) then
write(*,6) 'solvent', e%ww%el, e%ww%vdw, e%w%bond, e%w%angle, e%w%torsion, e%w%improper
end if

write(*,6) 'solute-solvent', e%pw%el, e%pw%vdw

if(use_lrf) then
write(*,6) 'lrf', e%lrf
end if

if(nqat .gt. 0) then
        write(*,6) 'q-atom', e%qx%el, e%qx%vdw, e%q%bond, e%q%angle, e%q%torsion, e%q%improper
end if

! restraints
write(*,*)
write(*,4) 'total', 'fix', 'slvnt_rad', 'slvnt_pol', 'shell', 'solute'
write(*,6) 'restraints', e%restraint%total, e%restraint%fix, &
        e%restraint%solvent_radial, e%restraint%water_pol, e%restraint%shell, &
        e%restraint%protein
write(*,*)

! totals
if(force_rms) then
        grms = sqrt(dot_product(d(:), d(:))/(3*natom))
        write(*,4) 'total', 'potential', 'kinetic', '', 'rms force'
        write(*,14) 'sum', e%potential+e%kinetic, e%potential, e%kinetic, grms
else
        write(*,4) 'total', 'potential', 'kinetic'
        write(*,6) 'sum', e%potential+e%kinetic, e%potential, e%kinetic
end if
14 format(a,t17, 3f10.2, 10x, f10.2)

! q-atom energies
if(nstates > 0) then
if(istep >= nsteps) then
  write(*,3) 'q-atom energies'
else
  write(*,2) 'q-atom energies', istep
end if

write(*,26) 'el', 'vdw' ,'bond', 'angle', 'torsion', 'improper'

do istate =1, nstates
  write (*,32) 'q-q', istate, eq(istate)%lambda, eq(istate)%qq%el, eq(istate)%qq%vdw
end do
write(*,*)
if(nat_solute > nqat) then !only if there is something else than q-atoms in topology
  do istate =1, nstates
        write (*,32) 'q-prot', istate,eq(istate)%lambda, eq(istate)%qp%el, eq(istate)%qp%vdw
  end do
  write(*,*)
end if

if(nwat > 0) then
  do istate =1, nstates
        write (*,32) 'q-wat', istate, eq(istate)%lambda, eq(istate)%qw%el, eq(istate)%qw%vdw
  end do
  write(*,*)
end if

do istate =1, nstates
  write (*,32) 'q-surr.',istate, eq(istate)%lambda, &
                eq(istate)%qp%el + eq(istate)%qw%el, eq(istate)%qp%vdw &
                + eq(istate)%qw%vdw
end do
write(*,*)

do istate = 1, nstates
  write (*,36) 'q-any', istate, eq(istate)%lambda, eq(istate)%qx%el,&
                eq(istate)%qx%vdw, eq(istate)%q%bond, eq(istate)%q%angle,&
                eq(istate)%q%torsion, eq(istate)%q%improper
end do
write(*,*)

write(*,22) 'total', 'restraint'
do istate = 1, nstates
  write (*,32) 'q-sum', istate, eq(istate)%lambda,&
                eq(istate)%total, eq(istate)%restraint
end do
do i=1,noffd
  write (*,360) offd(i)%i, offd(i)%j, hij(offd(i)%i, offd(i)%j), &
                offd2(i)%k, offd2(i)%l, offd(i)%rkl
360       format ('h(',i2,',',i2,') =',f8.2,' dist. between q-atoms',2i4, ' =',f8.2)
end do
end if

if(monitor_group_pairs > 0) then
        call centered_heading('monitoring selected groups of nonbonded interactions','=')
        write (*,37,advance='no')
        write (*,38) (istate,istate, istate=1,nstates)
        do i=1,monitor_group_pairs
                write (*,39,advance='no') i,monitor_group_pair(i)%vwsum, &
                        monitor_group_pair(i)%vwel,monitor_group_pair(i)%vwlj
                write (*,40) (monitor_group_pair(i)%vel(istate), &
                        monitor_group_pair(i)%vlj(istate), istate=1,nstates)
        end do
end if

write(*,'(80a)') '==============================================================================='


22      format('type   st lambda',2a10)
26      format('type   st lambda',6a10)
32      format (a,t8,i2,f7.4,2f10.2)
36      format (a,t8,i2,f7.4,6f10.2)
37  format ('pair   vwsum    vwel    vwvdw')
38  format (3(i4,':vel',i3,':vvdw'))
39  format (i2,f10.2,f8.2,f9.2)
40  format (3(2f8.2))


if(use_pbc .and. constant_pressure .and. istep>=nsteps ) then
        write(*,*)
        write(*,'(a)') '=========================== volume change summary ==========================='
        write(*,45) boxlength(1)*boxlength(2)*boxlength(3)
        write(*,*)
        write(*,46) 'total', 'accepted', 'ratio'
        write(*,47) 'attempts', volume_try, volume_acc, real(volume_acc)/volume_try
write(*,'(80a)') '==============================================================================='
end if
45 format('final volume: ', f10.3)
46 format(16x, 3a10)
47 format(a,t17, 2i10, f10.3)

end subroutine write_out

!-----------------------------------------------------------------------

subroutine write_trj

if(.not. trj_write(x)) then
        call die('failure to write to trajectory file')
end if

end subroutine write_trj

!-----------------------------------------------------------------------

subroutine write_xfin
! local variables
integer                                         ::      i,nat3

nat3 = natom*3

rewind (3)
write (3) nat3, (x(i),i=1,nat3)
write (3) nat3, (v(i),i=1,nat3)
!save dynamic polarisation restraint data
        if(wpol_restr .and. allocated(wshell)) then
        write (3) nwpolr_shell, wshell(:)%theta_corr
end if

if( use_pbc )then
        write(3) boxlength(:)
        write(3) boxcentre(:)
end if
end subroutine write_xfin

!-----------------------------------------------------------------------
!put molecules back in box for nice visualisation.
!change boxcentre if rigid_box_centre is off.
!update cgp_centers for lrf.
!-----------------------------------------------------------------------
subroutine put_back_in_box

real(8)                         ::      boxc(1:3)
integer                         ::      i, j, starten, slutet
!the borders of the periodic box
real(8)                         ::      x_max, x_min, y_max, y_min, z_max, z_min
real(8)                         :: cm(1:3)
integer                         ::  mvd_mol(1:nmol) !moved molecule 1=yes, 0=no
integer                         ::  k, ig
integer                         ::  pbib_start, pbib_stop

if( .not. rigid_box_centre ) then !if the box is allowed to float around, center on solute if present, otherwise center on solvent
        if( nat_solute > 0) then
                slutet = nat_solute
                !starten = ncgp_solute + 1
                
        else !if no solute present, centre box around solvent
                slutet = natom
                !starten = 1
        end if
        !find center
        boxc(:) = 0.0
        do i = 1,slutet
                boxc(1) = boxc(1) + x( 3*i-2 )
                boxc(2) = boxc(2) + x( 3*i-1 )
                boxc(3) = boxc(3) + x( 3*i   )
        end do
        boxc(:) = boxc(:)/slutet
        boxcentre(:) = boxc(:) !store new boxcentre
        ! starten = ncgp_solute + 1 !solute can not move around the corner
else
        !use boxcenter given in topology, ie. 'moving' solute
        boxc(:) = boxcentre(:)
        !starten = 1
end if

if (nodeid .eq. 0) then
!calculate the borders of the periodic box
x_max = boxc(1) + boxlength(1)/2
x_min = boxc(1) - boxlength(1)/2
y_max = boxc(2) + boxlength(2)/2
y_min = boxc(2) - boxlength(2)/2
z_max = boxc(3) + boxlength(3)/2
z_min = boxc(3) - boxlength(3)/2

mvd_mol(:) = 0 

!pbib_start and pbib_stop are the starting and stopping molecule indexes of which molecules to put back in box
pbib_start = 1
pbib_stop  = nmol
if ( .not. put_solute_back_in_box ) then !we're not putting solute back in box
        pbib_start = nmol - (natom - nat_solute)/(istart_mol(nmol) - istart_mol(nmol-1)) + 1    !(number of mol - number of solvent molecules + 1)
end if
if ( .not. put_solvent_back_in_box ) then !we're not putting solvent back in box
        pbib_stop  = nmol - (natom - nat_solute)/(istart_mol(nmol) - istart_mol(nmol-1))        !(number of mol - number of solvent molecules)
end if          


do i=pbib_start,pbib_stop
        cm(:) =0.0
        do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
                cm(1) = cm(1) + x(j*3-2)*mass(j)
                cm(2) = cm(2) + x(j*3-1)*mass(j)
                cm(3) = cm(3) + x(j*3  )*mass(j)
        end do
        cm(:) = cm(:) * mol_mass(i) !centre of mass of molecule i
                
        !x-direction
        if( cm(1) .gt. x_max) then !position of centre of mass

                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                                x(j*3-2) = x(j*3-2) - boxlength(1)
                                        end do
                                                mvd_mol(i) = 1
                else if ( cm(1) .lt. x_min ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                                x(j*3-2) = x(j*3-2) + boxlength(1)
                                        end do
                                                mvd_mol(i) = 1
        end if

       ! y-direction
        if( cm(2) .gt. y_max) then !position of centre of mass
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                                x(j*3-1) = x(j*3-1) - boxlength(2)
                                        end do
                                                mvd_mol(i) = 1
        else if ( cm(2) .lt. y_min ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                                x(j*3-1) = x(j*3-1) + boxlength(2)
                                        end do
                                                mvd_mol(i) = 1
        end if

        !z-direction
        if( cm(3) .gt. z_max) then !position of centre of mass
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                                x(j*3  ) = x(j*3  ) - boxlength(3)
                                        end do
                                                mvd_mol(i) = 1
        else if ( cm(3) .lt. z_min ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                                x(j*3  ) = x(j*3  ) + boxlength(3)
                                        end do
                                                mvd_mol(i) = 1
                end if
end do !over molecules
end if !if(nodeid .eq. 0)

!lrf: if molecule moved update all cgp_centers from first charge group to the last one
if (use_lrf) then

!broadcast mvd_mol(:) & x(:)
#if defined(use_mpi)
call mpi_bcast(mvd_mol, nmol, mpi_integer, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast mvd_mol(k)')
call mpi_bcast(x, nat3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast x')
#endif

do k=pbib_start,pbib_stop
                if (mvd_mol(k) == 1) then  
                                                do ig=iwhich_cgp(istart_mol(k)),iwhich_cgp(istart_mol(k+1)-1)
                                                                lrf(ig)%cgp_cent(:) = 0
                                                                do i  = cgp(ig)%first, cgp(ig)%last
                                                                                lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:) + x(cgpatom(i)*3-2:cgpatom(i)*3)
                                                                end do

                                                lrf(ig)%cgp_cent(:) = lrf(ig)%cgp_cent(:)/real(cgp(ig)%last - cgp(ig)%first +1)
                                                end do
                end if 
end do
end if

end subroutine put_back_in_box
!----------------------------------------------------------------------------
subroutine mc_volume()

real(8)                                                                 :: old_x(1:3*nat_pro)
type(energies)                                                  :: old_e
type(q_energies), dimension(1:nstates)  :: old_eq
real(8)                                                                 :: old_boxl(1:3), old_inv(1:3)
real(8)                                                                 :: old_v, new_v, deltalength
real(8)                                                                 :: deltav, deltae, deltaw
real(8)                                                                 :: box_min
integer                                                                 :: starten, slutet, i, j, sw_no !indeces
real(8)                                                                 :: randomno !random number
real(8)                                                                 :: new_x, new_y, new_z
real(8)                                                                 :: move_x, move_y, move_z
logical                                                                 :: acc
integer                                                                 :: longest 
real(8)                                                                 :: cubr_vol_ratio
real(8)                                                                 :: cm(1:3)

if (nodeid .eq. 0) then
write(*,8) 'volume change', istep
write(*,*)
write(*,'(a)') '---------- before move'
8 format('======================== ',a14,' at step ',i6,' ========================')
4 format(16x, 3a10)
6 format(a,t17, 3f10.3)
end if  !(nodeid .eq. 0)

!save the old energies,coordinates and forces
old_x(:) = x(:)
old_e = e
old_boxl(:) = boxlength(:)
old_inv(:) = inv_boxl(:)
if (use_lrf) then
        old_lrf(:) = lrf(:)
end if

#if defined(use_mpi)
!update modified coordinates  
call mpi_bcast(x, natom*3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast x')
#endif

call new_potential(old_e)   !compute energies from previous md-step

if (nodeid .eq. 0 ) then
old_e = e                !update to fresh e before changing volume
old_eq = eq(1:nstates)
old_v = old_boxl(1) * old_boxl(2) * old_boxl(3)


!new volume randomized
randomno = randm(pressure_seed) ! 0<=randomno<=1
randomno = randomno*2 - 1    !-1 <= randomno <= 1
deltav = randomno * max_vol_displ
new_v = deltav + old_v
cubr_vol_ratio = (new_v/old_v)**(1./3.)
write(*,4) 'old', 'new', 'delta'
write(*,6) 'volume', old_v, new_v, deltav
write(*,*)

!compute new boxlenth and inv_boxl
boxlength(1) = boxlength(1)*cubr_vol_ratio
boxlength(2) = boxlength(2)*cubr_vol_ratio
boxlength(3) = boxlength(3)*cubr_vol_ratio
inv_boxl(:) = 1.0/boxlength(:)
write(*,10) old_boxl
write(*,2) boxlength
write(*,*)
10 format('old boxlength', 3f10.3)
2 format('new boxlength ', 3f10.3)

!compare cut-offs with new boxsize
box_min = min( boxlength(1), boxlength(2), boxlength(3) )
!solute-solute cut-off radii
if( .not. (box_min .gt. rcpp*2) ) then
        write(*,*) 'solute-solute cut-off radii too large', rcpp
        call die('solute-solute cut-off radii too large')
!solvent-solvent
else if( .not. (box_min .gt. rcww*2) ) then
        write(*,*) 'solvent-solvent cut-off radii too large', rcww
        call die('solvent-solvent cut-off radii too large')
!solute-solvent
else if( .not. (box_min .gt. rcpw*2) ) then
        write(*,*) 'solute-solvent cut-off radii too large', rcpw
        call die('solute-solvent cut-off radii too large')
!q-atom
else if( .not. (box_min .gt. rcq*2) ) then
        write(*,*) 'q-atom cut-off radii too large', rcq
        call die('q-atom cut-off radii too large')
!lrf
else if( .not. (box_min .gt. rclrf*2) ) then
        write(*,*) 'lrf cut-off radii too large', rcq
        call die('lrf cut-off radii too large')
end if

!compute new coordinates after molecules and centre of mass
do i=1,nmol-1 !looping over all molecules but the last one
        cm(:) =0.0
        do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
                cm(1) = cm(1) + x(j*3-2)*mass(j)
                cm(2) = cm(2) + x(j*3-1)*mass(j)
                cm(3) = cm(3) + x(j*3  )*mass(j)
        end do
        cm(:) = cm(:) * mol_mass(i) !centre of mass of molecule i

        move_x = ( ( cm(1)-boxcentre(1) )*boxlength(1)/old_boxl(1) + boxcentre(1) ) - cm(1)
        move_y = ( ( cm(2)-boxcentre(2) )*boxlength(2)/old_boxl(2) + boxcentre(2) ) - cm(2)
        move_z = ( ( cm(3)-boxcentre(3) )*boxlength(3)/old_boxl(3) + boxcentre(3) ) - cm(3)

        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                x(j*3-2) = x(j*3-2) + move_x
                x(j*3-1) = x(j*3-1) + move_y
                x(j*3  ) = x(j*3  ) + move_z
        end do

end do !over molecules

! the last molecule
cm(:) = 0.0
do j = istart_mol(nmol),natom
        cm(1) = cm(1) + x(j*3-2)*mass(j)
        cm(2) = cm(2) + x(j*3-1)*mass(j)
        cm(3) = cm(3) + x(j*3  )*mass(j)
end do

cm(:) = cm(:) * mol_mass(nmol)

move_x = ( ( cm(1)-boxcentre(1) )*boxlength(1)/old_boxl(1) + boxcentre(1) ) - cm(1)
move_y = ( ( cm(2)-boxcentre(2) )*boxlength(2)/old_boxl(2) + boxcentre(2) ) - cm(2)
move_z = ( ( cm(3)-boxcentre(3) )*boxlength(3)/old_boxl(3) + boxcentre(3) ) - cm(3)

do j=istart_mol(nmol) , natom
        x(j*3-2) = x(j*3-2) + move_x
        x(j*3-1) = x(j*3-1) + move_y
        x(j*3  ) = x(j*3  ) + move_z            
end do
end if !nodeid .eq. 0

#if defined(use_mpi)
!update modified coordinates and boxlengths 
call mpi_bcast(x, natom*3, mpi_real8, 0, mpi_comm_world, ierr)
if (ierr .ne. 0) call die('init_nodes/mpi_bcast x')
call mpi_bcast(boxlength, 3, mpi_real8, 0, mpi_comm_world, ierr)
call mpi_bcast(inv_boxl, 3, mpi_real8, 0, mpi_comm_world, ierr)
#endif

!need to update entire lrf... sigh
if (use_lrf) then
        call cgp_centers
        if ( iuse_switch_atom == 1 ) then
                call nbpplist_box_lrf
                call nbpwlist_box_lrf
        else
                call nbpplis2_box_lrf
                call nbpwlis2_box_lrf
        endif
        call nbwwlist_box_lrf
end if


!compute the new potential, in parallel if possible
call new_potential( old_e )

if (nodeid .eq. 0) then
!jamfor nya med gamla
deltae = e%potential - old_e%potential
deltaw = deltae + pressure * deltav - nmol*boltz*temp0*log(new_v/old_v)
write(*,4) 'old', 'new', 'delta'
write(*,6) 'potential', old_e%potential, e%potential, deltae
write(*,*)

!accept or reject
if( deltaw<=0.0 ) then
        acc = .true.
else
        !slumpa tal mellan  0 coh 1
        randomno = randm(pressure_seed)
        if( randomno > exp(- deltaw / boltz / temp0) ) then
                acc = .false.
        else
                acc = .true.
        end if
end if

volume_try = volume_try + 1
write(*,'(a)') '---------- after move' 
if( acc ) then
        write(*,'(a)') 'volume change accepted'
        volume_acc = volume_acc + 1
else
        write(*,'(a)') 'volume change rejected'
        !put stuff back to what they were before
        x(:) = old_x(:)
        e = old_e
        eq(1:nstates) = old_eq(1:nstates)
        boxlength(:) = old_boxl(:)
        inv_boxl(:) = old_inv(:)
                if (use_lrf) then
                        lrf(:) = old_lrf(:)
                end if
end if

write(*,11) boxlength(1)*boxlength(2)*boxlength(3)
write(*,12) boxlength
11 format('final volume: ', f10.3)
12 format('final boxlength: ', 3f10.3)
write(*,*)
write(*,4) 'total', 'accepted', 'ratio'
write(*,7) 'attempts', volume_try, volume_acc, real(volume_acc)/volume_try
write(*,*)
7 format(a,t17, 2i10, f10.3)
write(*,'(80a)') '==============================================================================='
end if !(nodeid .eq. 0)

#if defined(use_mpi)
!make slave nodes put things back if rejected
call mpi_bcast(acc, 1, mpi_logical, 0, mpi_comm_world, ierr)
if (.not. acc) then
  if (nodeid .ne. 0) then
    x(:) = old_x(:)
    boxlength(:) = old_boxl(:)
    inv_boxl(:) = old_inv(:)
  end if
end if
#endif

end subroutine mc_volume        

!-----------------------------------------------------------------------------------------------------
subroutine new_potential( old )

type(energies), intent(in)              :: old  
integer                                                 :: istate,i

!zero all energies
e%potential = 0.0
e%pp%el  = 0.0
e%pp%vdw = 0.0
e%pw%el  = 0.0
e%pw%vdw = 0.0
e%ww%el  = 0.0
e%ww%vdw = 0.0
e%qx%el    = 0.0
e%qx%vdw   = 0.0
e%restraint%protein = 0.0
e%lrf      = 0.0 
do istate = 1, nstates
eq(istate)%qq%el = 0.0
eq(istate)%qq%vdw = 0.0
eq(istate)%qp%el = 0.0
eq(istate)%qp%vdw = 0.0
eq(istate)%qw%el = 0.0
eq(istate)%qw%vdw = 0.0
eq(istate)%restraint = 0.0
end do

!reset derivatives ---
d(:) = 0.0

#if defined (use_mpi)
!first post recieves for gathering data from slaves
if (nodeid .eq. 0) then
call gather_nonbond
end if
#endif

!compute the new potential
select case(ivdw_rule)
        case(vdw_geometric)
                call nonbond_pp_box
                call nonbond_pw_box
                if(qvdw_flag) then
                        call nonbond_qp_qvdw_box
                else
                        call nonbond_qp_box
                end if
                if(natom > nat_solute) then !if any solvent
                if(solvent_type == solvent_spc) then
                        !use the optimised spc routine when possible
                        call nonbond_ww_spc_box
                        call nonbond_qw_spc_box
                elseif(solvent_type == solvent_3atom) then !otherwise calc. lj with all atoms
                        call nonbond_3atomsolvent_box
                call nonbond_qw_3atom_box
                end if
                end if
        case(vdw_arithmetic)
                call nonbon2_pp_box
                call nonbon2_qp_box
                if(natom > nat_solute) then !if any solvent
                        call nonbon2_pw_box
                        call nonbon2_qw_box !no spc-specific optimised routines here
                        call nonbon2_ww_box
                end if
end select


 if (use_lrf) then
     call lrf_taylor
 end if

if (nodeid .eq. 0) then
 call p_restrain

 if(.not. qq_use_library_charges) then
        if(ivdw_rule .eq. 1 ) then
                call nonbond_qq
        elseif ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq
        end if
 else
        if ( ivdw_rule .eq. 1 ) then
                call nonbond_qq_lib_charges
        else if ( ivdw_rule .eq. 2 ) then
                call nonbon2_qq_lib_charges
        end if
 end if
end if

#if defined(use_mpi)
if (nodeid .ne. 0) then  !slave nodes
call gather_nonbond
end if
#endif

if (nodeid .eq. 0) then 
#if (use_mpi)
do i = 1, 3
    call mpi_waitall((numnodes-1),request_recv(1,i),mpi_status,ierr)
end do

!forces and energies are summarised
do i=1,numnodes-1
  d = d + d_recv(:,i)
  e%pp%el   = e%pp%el  + e_recv(i)%pp%el
  e%pp%vdw  = e%pp%vdw + e_recv(i)%pp%vdw
  e%pw%el   = e%pw%el  + e_recv(i)%pw%el
  e%pw%vdw  = e%pw%vdw + e_recv(i)%pw%vdw
  e%ww%el   = e%ww%el  + e_recv(i)%ww%el
  e%ww%vdw  = e%ww%vdw + e_recv(i)%ww%vdw
  e%lrf     = e%lrf    + e_recv(i)%lrf
  eq(1:nstates)%qp%el  = eq(1:nstates)%qp%el  + eq_recv(1:nstates,i)%qp%el
  eq(1:nstates)%qp%vdw = eq(1:nstates)%qp%vdw + eq_recv(1:nstates,i)%qp%vdw
  eq(1:nstates)%qw%el  = eq(1:nstates)%qw%el  + eq_recv(1:nstates,i)%qw%el
  eq(1:nstates)%qw%vdw = eq(1:nstates)%qw%vdw + eq_recv(1:nstates,i)%qw%vdw
end do
#endif

!summation of energies
do istate = 1, nstates
        ! update eq
        eq(istate)%qx%el  = eq(istate)%qq%el +eq(istate)%qp%el +eq(istate)%qw%el
        eq(istate)%qx%vdw = eq(istate)%qq%vdw+eq(istate)%qp%vdw+eq(istate)%qw%vdw
        e%qx%el    = e%qx%el    + eq(istate)%qx%el   *eq(istate)%lambda
        e%qx%vdw   = e%qx%vdw   + eq(istate)%qx%vdw  *eq(istate)%lambda

        ! update e%restraint%protein with an average of all states
        e%restraint%protein = e%restraint%protein + eq(istate)%restraint*eq(istate)%lambda
end do

e%potential = old%p%bond + old%w%bond + old%p%angle + old%w%angle + old%p%torsion + &
old%p%improper + e%pp%el + e%pp%vdw + e%pw%el + e%pw%vdw + e%ww%el + &
e%ww%vdw + old%q%bond + old%q%angle + old%q%torsion + &
old%q%improper + e%qx%el + e%qx%vdw + e%restraint%protein + e%lrf
end if !(nodeid .eq. 0)
end subroutine new_potential

!----------------------------------------------------------------------------------------

#if defined (use_mpi)
!***********************
!subroutine handling summation of nonbonded energies from slave nodes.
!***********************
! use the global vars
!  request_recv, e_send,eq_send,e_recv,eq_recv,d_recv
! allocate  - status


subroutine gather_nonbond()

integer,parameter                       :: vars=3
integer,dimension(3,numnodes-1)         :: tag
integer,dimension(vars)                 :: blockcnt, ftype
integer(kind=mpi_address_kind), dimension(vars) :: fdisp, base
integer                                 :: mpitype_package,mpitype_send
integer                                 :: i,istate

do i=1,numnodes-1
tag(1,i)=numnodes*100+i
tag(2,i)=numnodes*200+i
tag(3,i)=numnodes*300+i
end do

if (nodeid .eq. 0) then        !master

! post receives for each of the d/e/eq_recv structures
! e/eq_recv should really be handled with mpi_type_create_struct
! and d_recv's type should be handled correctly (it's kind=wp8)
! should preferably use size(d_recv, 1) for count
do i = 1,numnodes-1
  call mpi_irecv(d_recv(1,i), natom*3, mpi_real8, i, tag(1,i), mpi_comm_world, &
       request_recv(i,1),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/mpi_irecv d_recv')
  call mpi_irecv(e_recv(i), 3*2+1, mpi_real8, i, tag(2,i), mpi_comm_world, &
       request_recv(i,2),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/mpi_irecv e_recv')
  call mpi_irecv(eq_recv(1,i), nstates*2*2, mpi_real8, i, tag(3,i), mpi_comm_world, &
       request_recv(i,3),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/mpi_irecv eq_recv')
end do

else                  !slave nodes
e_send%pp%el  = e%pp%el
e_send%pp%vdw = e%pp%vdw
e_send%pw%el  = e%pw%el
e_send%pw%vdw = e%pw%vdw
e_send%ww%el  = e%ww%el
e_send%ww%vdw = e%ww%vdw
e_send%lrf    = e%lrf
eq_send(1:nstates)%qp%el  = eq(1:nstates)%qp%el
eq_send(1:nstates)%qp%vdw = eq(1:nstates)%qp%vdw
eq_send(1:nstates)%qw%el  = eq(1:nstates)%qw%el
eq_send(1:nstates)%qw%vdw = eq(1:nstates)%qw%vdw

! see comments above on the irecv part
call mpi_send(d, natom*3, mpi_real8, 0, tag(1,nodeid), mpi_comm_world,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/send d')
call mpi_send(e_send, 3*2+1, mpi_real8, 0, tag(2,nodeid), mpi_comm_world,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/send e_send')
call mpi_send(eq_send, nstates*2*2, mpi_real8, 0, tag(3,nodeid), mpi_comm_world,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/send eq_send')

end if
end subroutine gather_nonbond

#endif
!----------------------------------------------------------------------------------------
!*******************************************************
!will find and return the xtop atom number from 
!  residue number and atom number in residue from
!  library sequence.
! uses global variables: xtop,nres,res
!*******************************************************

integer function get_atom_from_resnum_atnum(aid)
!arguments
character(*), intent(in)        ::      aid     !string=residue:atom
        
!locals
integer                                         ::      separator_pos
character(len=20)                       ::      res_str
character(len=5)                        ::      atom_str
integer                                         ::      filestat
integer                                         ::      resnum, atnum

get_atom_from_resnum_atnum = 0

separator_pos = scan(aid, ':')
if(separator_pos < 2 .or. separator_pos == len_trim(aid)) return !no valid colon found
res_str = aid(1:separator_pos-1)
atom_str = aid(separator_pos+1:len_trim(aid))
read(res_str, *, iostat=filestat) resnum
read(atom_str, *, iostat=filestat) atnum
if(filestat > 0) return

!residue must be in topology
if(resnum < 1 .or. resnum > nres) then                     
  return                                                 
end if

if(atnum .le. (res(resnum+1)%start - res(resnum)%start)) then
  get_atom_from_resnum_atnum = res(resnum)%start + atnum - 1
return
end if

!we have an error: 
write(*, 120) atnum, resnum
call die('error in finding atom number from resnum:atnum.')

120     format('>>>>> error: there is no atom number ',i4,' in residue ',i4,'.')
end function get_atom_from_resnum_atnum

!----------------------------------------------------------------------------------------

end module md




