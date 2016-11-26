! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! qcp.f90
! Q is written by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! Quantum classical path module
! written by Paul Bauer
! based on bisection sampling QCP developed by Dan Thomas Major and Jiali Gao
module QCP

use TOPO
use QATOM
use SIZES
use NRGY
use POTENE
use GLOBALS
use QMATH
use MPIGLOB
! this module does the actual path integral calculations and metropolis sampling of the bead configurations
! code is directly based on the paper, with more bugs than I can count

implicit none
! global QCP values
real(kind=prec) :: beta, wl_lam, tmp_wl_lam, cboltz, hbar, convert, pi_fac
! array of wavelengths for every atom 
real(kind=prec),allocatable :: wl_lam0(:), wl_lam1(:), wl_lam2(:)
! array of square root of mass ratios for every atom, for mass perturbed QCP (I guess)
! only used if user has set masses with qcp_mass
real(kind=prec),allocatable :: sqmass(:)
! energy storage for beads and particles
! one Q energy array for every FEP lambda?
! added in the end to global EQ as last entry in ene_header%array
TYPE(OQ_ENERGIES),allocatable   :: EQtmp(:),EQsave(:),EQpot_ave(:),qcp_EQbead_old(:,:),qcp_EQbead_new(:,:)
real(kind=prec),allocatable     :: qcp_Ebeta(:),qcp_EQbeta(:,:),total_EQPI(:)
TYPE(ENERGIES)                          :: qcp_E_tot
TYPE(OQ_ENERGIES),allocatable           :: qcp_EQ_tot(:)
real(kind=prec)                         :: total_EPI,Epot_ave,qcp_Ering_ave,Esave,Etmp,qcp_Ebead
real(kind=prec)                         :: qcp_Ering_new,qcp_Ering_old,qcp_Ering
real(kind=prec),allocatable             :: qcp_Ebead_old(:),qcp_Ebead_new(:)
! for gaussian deviation algorithm
integer                         :: gauss_set = 0
real(kind=prec)                 :: gval = zero
! coordinates for ring
TYPE(qr_vec),allocatable        :: qcp_coord(:,:),qcp_coord_old(:,:)

contains

subroutine qcp_init
! set up constants after reading in information from inputs
! no matter if input is in MD or Postprocessing
! arguments
! TODO set constants in MISC
integer :: i

convert    = 4.184_prec / 6.022E23_prec
cboltz     = 1.38064852E-23_prec
hbar       = 1.054571800E-34_prec

!all temp dependent variables set with current instantanious temperature
!during actual QCP calc

! TODO need allocation check
allocate(wl_lam0(qcp_atnum),wl_lam1(qcp_atnum),wl_lam2(qcp_atnum))
! one array each for fep states and qcp steps to save energy
allocate(qcp_EQ(nstates))
! allocate bead coordinate array
allocate(qcp_coord(qcp_atnum,qcp_size),qcp_coord_old(qcp_atnum,qcp_size))

allocate(qcp_Ebeta(maxval(qcp_steps)))
allocate(qcp_EQbeta(maxval(qcp_steps),nstates),EQsave(nstates),total_EQPI(nstates),&
        EQpot_ave(nstates),EQtmp(nstates))

allocate(qcp_EQbead_old(qcp_size,nstates),qcp_EQbead_new(qcp_size,nstates))

allocate(qcp_EQ_tot(nstates))

allocate(x_save(natom),d_save(natom))

allocate(qcp_Ebead_old(qcp_size),qcp_Ebead_new(qcp_size))

! will need this later
qcp_EQ(1:nstates)%lambda = EQ(1:nstates)%lambda

if (use_qcp_mass) then
	allocate(sqmass(qcp_atnum))
end if

! set wl values for every atom
do i = 1, qcp_atnum
call qcp_init_beads(qcp_atom(i),qcp_coord(i,:))
end do

! what are TIAC, CHIN, CHAC????
end subroutine qcp_init

subroutine qcp_init_beads(atom,qcpvec)
! create initial guess for bead positions, wont even try to read stuff in
! again stolen from original code by D.T. Major
! arguments
integer                         :: atom
TYPE(qr_vec)                    :: qcpvec(:)
! locals
! displ -> displacment, dependent on heavy atom or not
real(kind=prec)                 :: angle,percent,displ_l,displ_h,displ,dangl
integer                         :: i, irand,jrand,krand

qcpvec = qcpvec * zero

end subroutine qcp_init_beads

function qcp_bisect(qcpvec,wavelength,level)
! this one takes care of the actual bisection step to generate the bead positions
! shamelessly stolen from original code by D.T. Major and J. Gao
! based on following reference for the crazy
! Levy, P., Compostio Math. 7, 283, 1939
! Ceperly, F.M.; Pollock, E.L. Monte Carlo methods in theoretical physics, p 35. ETS Editrice, Pisa 1992

! arguments
! all the coordinates for the beads of a given PI atom
TYPE(qr_vec)                    :: qcpvec(:)
TYPE(qr_vec)                    :: qcp_bisect(size(qcpvec))
! wavelength, used for step
real(kind=prec)                 :: wavelength
! bisection level, I guess
integer                         :: level

! locals
real(kind=prec),parameter       :: half=0.5_prec
integer                         :: n,i,j,nbsect,bbsect,nextbead,medbead,thisbead,bead
TYPE(qr_vec)                    :: midpoint
real(kind=prec)                 :: lbsect

! for a given particle, go through all beads
! perform number of steps depending on bisection level
! starting at highest and working down
        do bead = 1, qcp_size
        do i = level, 1, -1
        lbsect   = q_sqrt(real(2**(i-1),kind=prec))
! set number of bisection steps
        nbsect   = 2**(level-i)
! current bisection size
        bbsect   = 2**i
        thisbead = bead
        nextbead = thisbead + bbsect
        medbead  = int(thisbead+nextbead)/2
        do j = 1, nbsect
        if (thisbead.gt.qcp_size) thisbead = thisbead - qcp_size
        if (nextbead.gt.qcp_size) nextbead = nextbead - qcp_size
        if (medbead .gt.qcp_size) medbead  = medbead  - qcp_size
        midpoint = (qcpvec(thisbead) + qcpvec(nextbead)) * half
        qcp_bisect(medbead) = qcp_gauss(wavelength*lbsect,midpoint)
        thisbead = thisbead + bbsect
        nextbead = nextbead + bbsect
        medbead  = medbead  + bbsect
        end do ! nbsect
        end do ! level
        end do ! bead

end function qcp_bisect

function qcp_center(qcpvec)
! centers all beads of a qcp atom on the middle
! arguments
TYPE(qr_vec)                    :: qcpvec(:)
TYPE(qr_vec)                    :: qcp_center(size(qcpvec))
! locals
integer                         :: i
real(kind=prec)                 :: inv_beads
TYPE(qr_vec)                    :: center

inv_beads = one/(real(qcp_size ,kind=prec))
center = center * zero
do i = 1, qcp_size
center = center + qcpvec(i)
end do
center = center * inv_beads
do i = 1, qcp_size
qcp_center(i) = qcpvec(i) - center
end do

end function qcp_center


TYPE(qr_vec) function qcp_gauss(sigma,coord)
! returns the gaussian distribution weighted coordinate
! reimplemented from code provided by D.T. Major
! from Numerical Recipies, Ch. 7.2, according to Box-Muller algorithm
! separate from general Q gaussian algorith to not be influenced by this
! arguments
TYPE(qr_vec)                    :: coord
real(kind=prec)                 :: sigma
! locals
real(kind=prec)                 :: val1,val2,fac,gdev,rval
real(kind=prec),parameter       :: two=2.0_prec
integer                         :: i

do i = 1, 3

if (gauss_set.eq.0) then

        do 
        val1 = 2.0_prec * qcp_randm() - one
        val2 = 2.0_prec * qcp_randm() - one
        rval = val1**2 + val2**2
        if (rval.gt.zero .and. rval .lt.one) exit
        end do

        fac  = q_sqrt(-two*q_logarithm(rval)/rval)
        gval = val1 * fac
        gdev = val2*fac*sigma
        gauss_set = 1
else
        gdev = gval * sigma
        gauss_set = 0
end if
if (i.eq.1) qcp_gauss%x = coord%x + gdev
if (i.eq.2) qcp_gauss%y = coord%y + gdev
if (i.eq.3) qcp_gauss%z = coord%z + gdev
end do

end function qcp_gauss

subroutine qcp_run(temp,E_init,EQ_init)
! main meaty subroutine that performs the actual run to generate coordinates
! and evaluate energies 
! based on original code by D.T. Major and J.Gao
! after ripping out all the stuff I don't need
! only for bisection now, have to implement rest later
! arguments
real(kind=prec)                         :: temp
TYPE(ENERGIES), INTENT(in)              :: E_init
TYPE(OQ_ENERGIES),INTENT(in)            :: EQ_init(:)
! locals
integer                                 :: ncent,nmove,qcp_loop,istate
integer                                 :: i,ia,j,step
real(kind=prec)                         :: bfactor,deltaH,irand,MCaccept,MCreject,MCrepeat
! in case of MPI this is a mess, because only the pot_ene stuff needs it
! everything else needs to be on the master node to keep it clean

! first set of energies will come from independent calculation, in MD this is ther normal EQ,
! in case of post processing this will be EQ done on the first restart
! first thing needed is saving coordinates and forces to restore them later

x_save = x
d_save = d

! set temperature dependent variables
if (nodeid .eq. 0) then
beta       = convert / (temp * cboltz)
tmp_wl_lam = angstrom * cboltz * temp / hbar
pi_fac     = (0.5_prec * real(qcp_size, kind=prec) * tmp_wl_lam**2 * amu) / convert
tmp_wl_lam = hbar / (2.0_prec * tmp_wl_lam *  angstrom * real(qcp_size,kind=prec) * amu ) 

! set wl values for every atom
do i = 1, qcp_atnum
	wl_lam0(i) = tmp_wl_lam * winv(iqseq(qcp_atom(i)))
	if(use_qcp_mass) sqmass = sqrt(one/winv(iqseq(qcp_atom(i)))/QCP_mass(i))
	wl_lam1(i) = sqrt(wl_lam0(i))
	wl_lam2(i) = 2.0_prec * wl_lam1(i)
end do

end if

do i = 1, qcp_atnum
call qcp_init_beads(qcp_atom(i),qcp_coord(i,:))
end do


qcp_EQ_tot = qcp_EQ_tot * zero
qcp_Ebead_old = zero
qcp_Ebead_new = zero
qcp_Ebead     = zero
Etmp  = zero
Esave = zero
EQtmp = EQtmp * zero
do i = 1, maxval(qcp_steps)
qcp_EQbead_old(i,:) = qcp_EQbead_old(i,:) * zero
qcp_EQbead_new(i,:) = qcp_EQbead_new(i,:) * zero
end do
if (nodeid.eq.0) then
! set number of beads to move, 1 bead for bisection, center all
ncent = qcp_size
nmove = 1
end if

! main loop, first equil, then sampling
! set steps according to qcp_steps
do qcp_loop = 1, 2

if (qcp_loop .eq. 2) then
! sampling needs to set some stuff first

! move all beads, one after the other
! first, each bead for each atom
! only on node 0 for mpi
do i = 1, qcp_size
if (nodeid.eq.0) then
do j = 1, qcp_atnum
ia = qcp_atom(j)
x(ia) = x_save(ia) + qcp_coord(j,i)
end do ! qcp_atnum
end if

! MPI BCast new coordinates to nodes
#ifdef USE_MPI
call MPI_BCast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast x')
#endif

! all atoms at first perturbed position -> get energy
! E and EQ and set to zero in pot_ene, don't need to do this manually
call pot_energy(qcp_E,qcp_EQ,.false.)
! rest set to zero
if (nodeid.eq.0) then
Etmp = Etmp + qcp_E%potential
EQtmp(:) = EQtmp(:) + qcp_EQ(:)
qcp_Ebead_old(i) = qcp_E%potential
do  istate = 1, nstates
qcp_EQbead_old(i,istate) = qcp_EQ(istate)
end do
end if
end do ! qcp_size

! what is factor? CHFACN is one for not using strange CHIN method stuff
if (nodeid.eq.0) then
        Esave = Etmp / real(qcp_size,kind=prec)
        Etmp = Esave
end if

if (nodeid.eq.0) then
        EQsave(:) = EQtmp(:) / real(qcp_size,kind=prec)
        EQtmp = EQsave
end if
end if ! qcp_loop .eq.2

if (nodeid .eq.0) then
qcp_Ering = zero
do i = 1, qcp_atnum
qcp_Ering = qcp_Ering + qcp_energy(qcp_atom(i),qcp_coord(i,:),pi_fac)
end do

qcp_Ering_old = qcp_Ering

! set some average values
! why another one of those
total_EPI = zero
qcp_Ering_ave = zero
Epot_ave  = zero
EQpot_ave = EQpot_ave * zero
total_EQPI = zero
! MC stuff
MCrepeat = one
MCaccept = zero
MCreject = zero

end if

do step = 1 , qcp_steps(qcp_loop)
! save bead coords

if (nodeid.eq.0) then

qcp_coord_old = qcp_coord
! call actual bisection step algortihm
! and recenter ring
! for each ring
qcp_Ering = zero
do i = 1, qcp_atnum
qcp_coord(i,:) = qcp_bisect(qcp_coord(i,:),wl_lam1(i),qcp_level)
qcp_coord(i,:) = qcp_center(qcp_coord(i,:))
qcp_Ering = qcp_Ering + qcp_energy(qcp_atom(i),qcp_coord(i,:),pi_fac)
end do
qcp_Ering_new = qcp_Ering

end if

! needed for selection, different for other algorithms
! those would also need Bcasting of energies then
deltaH = zero

bfactor = q_exp(-deltaH*beta)

! get some random number?
irand = qcp_randm()

! this is the monte carlo part
if (bfactor .ge. irand) then
! energy only calculated for sampling, not for equilibration? weird!
        if(qcp_loop .eq. 2) then
                Esave  = Etmp
                EQsave = EQtmp
                if(qcp_veryverbose) then
                write(*,'(a,i4,a,f12.6,a,f12.6)') 'Values at step ',step,' for Esave = ',Esave,' , and Etmp = ',Etmp
                        do istate = 1 ,nstates
                        write(*,'(a,i4,a,i4)') 'Values at step ',step,' state = ',istate
                        write(*,'(a,f12.6,a,f12.6)') 'EQsave = ',EQsave(istate)%total,&
                                ' , EQtmp = ',EQtmp(istate)%total
                        end do
                end if
                do i = 1, qcp_size
                if (nodeid.eq.0) then
                        do j = 1, qcp_atnum
                        ! only atoms with qcp_coords need update, saved in x_save
                        x(qcp_atom(j)) = x_save(qcp_atom(j)) + qcp_coord(j,i)
                        end do ! qcp_atnum
                end if
#ifdef USE_MPI
call MPI_Bcast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast x')
#endif

                call pot_energy(qcp_E,qcp_EQ,.false.)
                if (nodeid.eq.0) then
                        qcp_Ebead_new(i) = qcp_E%potential
                        qcp_EQbead_new(i,:) = qcp_EQ(:)
                        Esave = Esave + ((qcp_Ebead_new(i)-qcp_Ebead_old(i))/real(qcp_size,kind=prec))
                        EQsave(:) = EQsave(:) + ((qcp_EQbead_new(i,:) - qcp_EQbead_old(i,:)) &
                                /real(qcp_size,kind=prec))
                        qcp_Ebead_old(i) = qcp_Ebead_new(i)
                        qcp_EQbead_old(i,:) = qcp_EQbead_new(i,:)
                end if
                end do ! qcp_size


                ! sum up some stuff
                if (nodeid .eq.0) then
                        qcp_Ebeta(step) = q_exp(-beta*(Etmp-E_init%potential))
                        if(qcp_veryverbose) then
                                write(*,'(a,i4,a,f12.6)') 'Current beta for step ',step,' = ',qcp_Ebeta(step)
                                write(*,'(a,f12.6)') 'Etmp = ',Etmp
                                write(*,'(a,f12.6)') 'E_init = ',E_init%potential
                        end if
                        do istate = 1, nstates
                        qcp_EQbeta(step,istate) = q_exp(-beta*(EQtmp(istate)%total-EQ_init(istate)%total))
                        if(qcp_veryverbose) then
                                write(*,'(a,i4)') 'Current step = ',step
                                write(*,'(a,i4)') 'Current state = ',istate
                                write(*,'(a,f12.6)') 'Current beta = ',qcp_EQbeta(step,istate)
                                write(*,'(a,f12.6)') 'EQtmp = ',EQtmp(istate)%total
                                write(*,'(a,f12.6)') 'EQ_init = ',EQ_init(istate)%total
                        end if

                        end do
                        total_EPI = total_EPI + MCrepeat * qcp_Ebeta(step)
                        total_EQPI(:) = total_EQPI(:) + MCrepeat * qcp_EQbeta(step,:)
                        Epot_ave  = Epot_ave  + MCrepeat * Etmp
                        EQpot_ave(:) = EQpot_ave(:) + (EQtmp(:) * MCrepeat)
                        Etmp = Esave
                        EQtmp = EQsave
                end if
        end if ! qcp_loop .eq. 2
        if (nodeid .eq.0) then
                MCaccept = MCaccept + MCrepeat
                qcp_Ering_ave = qcp_Ering_ave + MCrepeat * qcp_Ering_new
                qcp_Ering_old = qcp_Ering_new 
                MCrepeat = one
        end if
else ! failed to accept step
        if (nodeid.eq.0) then
                qcp_coord = qcp_coord_old
                MCreject = MCreject + one
                MCrepeat = MCrepeat + one
        end if
end if ! bfactor .ge. irand

! write info about qcp steps if desired -> flag is qcp_veryverbose
if (nodeid.eq.0 .and. qcp_veryverbose .and. qcp_loop .eq. 2) then
        write(*,'(a)') 'Information for individual step of QCP PI'
        write(*,'(a)') 'Current status is free particle sampling'
        write(*,'(a,i4)') 'Current conformation = ',step
        write(*,'(a,f12.6)') 'Current total particle energy = ',Esave
        write(*,'(a,10f12.6)') 'Current Q particle energy, for each state = ',EQsave(:)%total
        write(*,'(a,f12.6)') 'Current ring polymer energy = ',qcp_Ering_new
end if

end do ! step

if (nodeid.eq.0) then
! write out info about QCP success rate and energies
        if (qcp_loop .eq.1) then
                write(*,'(a)') 'QCP summary for free particle equilibration'
        else
                write(*,'(a)') 'QCP summary for free particle sampling'
        end if
        write(*,'(a,i4)') 'Total number of steps is ',step-1
        write(*,'(a,i4)') 'Number of accepted MC steps = ',int(MCaccept)
        write(*,'(a,i4)') 'Number of rejected MC steps = ',int(MCreject)
end if

! make averages for saving to file
if (nodeid.eq.0 )then
        if (MCaccept .eq. zero) then
                MCaccept = real(qcp_steps(qcp_loop),kind=prec)
                qcp_Ering_ave = qcp_Ering_ave/MCaccept
        end if
end if
end do ! qcp_loop

if (nodeid .eq.0) then

total_EPI = total_EPI / MCaccept
Epot_ave  = Epot_ave  / MCaccept
total_EQPI(:) = total_EQPI(:) / MCaccept
EQpot_ave(:) = EQpot_ave(:) / MCaccept

if(qcp_verbose) then
        write(*,'(a,f12.6)') 'Total PI energy = ',total_EPI
        write(*,'(a,f12.6)') 'Final PI potential energy = ',Epot_ave
        do i = 1 , nstates
        write(*,'(a,i4)') 'Total QPI energy for state ',i
        write(*,'(f12.6)') total_EQPI(i)
        end do
        do i = 1 , nstates
        write(*,'(a,i4)') 'istate = ',i
        write(*,'(a,f12.6)') 'Final PI Q potential energy = ',EQpot_ave(i)%total
        end do
end if

! now we have all the energies for the PI sampling, can print info for user and write to global EQ_save 
! qcp_pos is set during loading, need this later for mass perturbation with two arrays
call qatom_savetowrite(EQpot_ave,qcp_pos)

! set everything back to how we found it 
x = x_save
d = d_save

end if
! Bcast when using MPI
#ifdef USE_MPI
call MPI_Bcast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast x')
call MPI_Bcast(d,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast d')
#endif

! and we are done 

end subroutine qcp_run


real(kind=prec) function qcp_energy(atom,coord,wavelength)
! calculates the ring energies of the particles
! call this for each ring
! arguments
TYPE(qr_vec)                    :: coord(:)
integer                         :: atom
real(kind=prec)                 :: wavelength
! locals
real(kind=prec)                 :: fk
integer                         :: i

qcp_energy = zero
! will need adapting for different QCP types
fk = wavelength/winv(atom)

do i = 1, qcp_size - 1
qcp_energy = qcp_energy + (fk * q_dist4(coord(i),coord(i+1)))
end do

! this only later for closed chain chain polymer
qcp_energy = qcp_energy + (fk * q_dist4(coord(1),coord(qcp_size)))


end function qcp_energy


real(kind=prec) function qcp_randm()
! same random number generator as in charmm code provided by D.T. Major
!
!-----------------------------------------------------------------------
!
!     RANDOM NUMBER GENERATOR
!     SEPTEMBER, 1986: WLJ & J. M. Briggs
!
!     NOTES ON USAGE:
!     THE STATEMENT "COMMON/RANDS/IRN" MUST APPEAR IN EACH
!     PROGRAM,SUBROUTINE... WHERE THE SEED VALUE (IRN) IS
!     BEING WRITTEN OUT OR BEING READ IN. THIS FUNCTION IS
!     USED IN THE SAME MANNER AS THAT USED ON THE GOULD,
!     (I.E. "RANDN = RANU()" ). IRN SHOULD INITIALLY BE
!     AN ODD, I6 INTEGER BUT CAN GET AS BIG AS I7. JMB
!     IMOD-1 UNIQUE RANDOM NUMBERS ARE GENERATED . A SINGLE
!     PRECISION LINEAR CONGRUENTIAL GENERATOR IS USED WITH
!     SHUFFLING OF THE CONSTANT AND OF WHEN THE SHUFFLING IS
!     DONE. CONSEQUENTLY, THE PERIOD IS EXTREMELY LONG - NONE
!     WAS FOUND IN TESTS GENERATING MILLIONS OF RANDOM NUMBERS.
!-----------------------------------------------------------------------
real(kind=prec) :: RNJ,FAC
integer,save :: ICNT=0   ,ICHG=1167     ,ICHG0=1167  ,ICN0=458753759,  &
                IMUL=1173,ICON=458753759,IMOD=1048573, &
                JMUL=1161,JCON=458716759,JMOD=1048573,JRN=124690
integer         :: IRN

  ICNT = ICNT+1
  IF (ICNT.EQ.ICHG) THEN
     !
     !-----------------------------------------------------------------------
     !
     !     CHANGE ICON USING SECONDARY GENERATOR
     !
     !-----------------------------------------------------------------------
     !
     JRN = JRN*JMUL+JCON
     JRN = MOD(JRN,JMOD)
     RNJ = real(JRN,kind=prec)/real(JMOD,kind=prec)
     IF (RNJ.GT.0.5_prec) THEN
        FAC = one + 0.5_prec*RNJ
     ELSE
        FAC = one - 0.5_prec*RNJ
     ENDIF
     FAC = real(ICN0,kind=prec)*FAC
     ICON = INT(FAC)
     !
     !-----------------------------------------------------------------------
     !
     !     CHANGE ICHG USING SECONDARY GENERATOR
     !
     !-----------------------------------------------------------------------
     !
     JRN = JRN*JMUL+JCON
     JRN = MOD(JRN,JMOD)
     RNJ = real(JRN,kind=prec)/real(JMOD,kind=prec)
     IF (RNJ.GT.0.5_prec) THEN
        FAC = one + 0.5_prec*RNJ
     ELSE
        FAC = one - 0.5_prec*RNJ
     ENDIF
     FAC  = real(ICHG0,kind=prec)*FAC
     ICHG = INT(FAC)
     ICNT = 0
  ENDIF
  !
  !-----------------------------------------------------------------------
  !
  !     GENERATE RANDOM NUMBER
  !
  !-----------------------------------------------------------------------
  !
  IRN = IRN*IMUL+ICON
  IRN = MOD(IRN,IMOD)
  qcp_randm = real(IRN,kind=prec)/real(IMOD,kind=prec)

end function qcp_randm


!-----------------------------------------------------------------------
logical function qcp_initialize()                  
! basically a copy of original initialze function in md.f90
! to make sure all variables are properly read in the same way as for normal md
! and that all variables are assigned
! local variables
character					:: text*200
integer						:: i,j,length,ii
real(kind=prec)						:: stepsize
real(kind=prec)						:: lamda_tmp(max_states)
integer						:: fu, fstat
real(kind=prec)						::	rjunk
integer						::	ijunk

! local parameters
integer						:: num_args
character(200)				:: infilename
logical						::	yes
logical						::	need_restart,temp_defined
character(len=200)				::	instring
logical						::	inlog
integer						::	mask_rows, number
integer,parameter				:: maxmaskrows=30
character(len=200)				:: gc_mask_tmp(maxmaskrows),str,str2
logical                                         :: shake_all,shake_all_solvent,shake_all_solute
logical                                         :: shake_all_hydrogens,shake_all_heavy
character(200)                                  :: qcp_select,qcp_size_name

! read name of input file from the command line
num_args = command_argument_count()
if (num_args .lt. 1) call die('no input file specified on the command line')
#if defined(CRAY)
call pxfgetarg(num_args, infilename, 200, i)
#elif defined(MPICH)
call getarg(1, infilename)
#else
call getarg(num_args, infilename)
#endif
text = 'Reading input from '//infilename
call centered_heading(trim(text), '-')

qcp_initialize = .true. 

if(.not. prm_open_section('PBC', infilename)) then
box = .false.
write(*,'(a)') 'Boundary: sphere'
else
box = .true.
write(*,'(a)') 'Boundary: periodic box'

end if !section PBC


if(.not. prm_open_section('general')) then
call prm_close
! open input file
fu = freefile()
open(unit=fu, file=infilename, action='read', form='formatted', status='old', iostat=fstat)
if (fstat .ne. 0) call die('error opening input file '//infilename)
        qcp_initialize = .false.
        close(fu)
return
end if

need_restart = .false. !flag for restart file required
! need restart file to calculate velocities and temperatures
temp_defined = .false.
! default is that we will use temperatures from restart file velocities
! convert to internal time units once and for all.
! don't need to read in stepsize

! --- Temperature, Thermostat etc.
if(.not. prm_get_real8_by_key('temperature', Temp0)) then
write(*,*) 'temperature not specified (section general), will need restart'
need_restart = .true.
temp_defined = .false.
else
write(*,'(a,f8.3)') 'Temperature set to ',Temp0
temp_defined = .true.
end if
! read in temperature here, if not defined will have to use info from restart

! --- shake, LRF
if(.not. prm_get_logical_by_key('shake_solvent', shake_solvent, .true.)) then
write(*,'(a)') '>>> Error: shake_solvent must be on or off.'
qcp_initialize = .false.
end if
!write(*,17) 'all solvent bonds', onoff(shake_solvent)
17	format('Shake constraints on ',a,' to ',a,' : ',a3)

if(.not. prm_get_logical_by_key('shake_solute', shake_solute, .false.)) then
write(*,'(a)') '>>> Error: shake_solute must be on or off.'
qcp_initialize = .false.
end if 
!write(*,17) 'all solute bonds', onoff(shake_solute)

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .true.)) then
write(*,'(a)') '>>> Error: shake_hydrogens must be on or off.'
qcp_initialize = .false.
end if 

if(.not. prm_get_logical_by_key('shake_heavy', shake_heavy, .false.) ) then
write(*,'(a)') '>>> Error: shake_heavy must be on or off.'
qcp_initialize = .false.
end if


if(prm_get_logical_by_key('shake_all',shake_all)) then
if(shake_all) then
       shake_solute    = .true.
       shake_solvent   = .true.
       shake_hydrogens = .true.
       shake_heavy     = .true.
else
       shake_solute    = .false.
       shake_solvent   = .false.
       shake_hydrogens = .false.
       shake_heavy     = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_solvent',shake_all_solvent)) then
if(shake_all_solvent) then
       shake_solvent   = .true.
       shake_hydrogens = .true.
       shake_heavy     = .true.
else
       shake_solvent   = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_solute',shake_all_solute)) then
if(shake_all_solute) then
       shake_solute    = .true.
       shake_hydrogens = .true.
       shake_heavy     = .true.
else
       shake_solute    = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_hydrogens',shake_all_hydrogens)) then
if(shake_all_hydrogens) then
        shake_solute    = .true.
        shake_solvent   = .true.
        shake_hydrogens = .true.
else
        shake_hydrogens = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_heavy',shake_all_heavy)) then
if(shake_all_heavy) then
         shake_solute  = .true.
         shake_solvent = .true.
         shake_heavy   = .true.
else
         shake_heavy   = .false.
end if
end if

if(shake_solvent) then
if(shake_heavy) write(*,17) 'solvent bonds','heavy atoms',onoff(shake_heavy)
if(shake_hydrogens) write(*,17) 'solvent bonds','hydrogens',onoff(shake_hydrogens)
else
write(*,17) 'solvent bonds','any atoms',onoff(shake_solvent)
endif
if(shake_solute) then
if(shake_heavy) write(*,17) 'solute bonds','heavy atoms',onoff(shake_heavy)
if(shake_hydrogens) write(*,17) 'solute bonds','hydrogens',onoff(shake_hydrogens)
else
write(*,17) 'solute bonds','any atoms',onoff(shake_solute)
endif
!write(*,17) 'all bonds to hydrogen', onoff(shake_hydrogens)


yes = prm_get_logical_by_key('lrf', use_LRF, .true.)
if(use_LRF) then
       write(*,20) 'LRF Taylor expansion outside cut-off'
else 
       write(*,20) 'standard cut-off'
end if

20	format ('Nonbonded method   = ',a)

! --- Rcpp, Rcww, Rcpw, Rcq, RcLRF
! change for PBC, q-atoms interact with everything for now, placeholder until we have 
! fixed it to calculate the interaction between the mirror images
if(.not. prm_open_section('cut-offs')) then
write(*,'(a)') 'No cut-offs section, default cut-offs used'
Rcpp = rcpp_default
Rcww = rcww_default
Rcpw = rcpw_default
if (box) then
Rcq = rcq_default_pbc
else
Rcq = rcq_default_sph
end if
RcLRF = rcLRF_default
else
if(.not. prm_get_real8_by_key('solute_solute', Rcpp, rcpp_default)) then
        write(*,'(a)') 'solute-solute cut-off set to default'
end if
if(.not. prm_get_real8_by_key('solvent_solvent', Rcww, rcww_default)) then
        write(*,'(a)') 'solvent-solvent cut-off set to default'
end if
if(.not. prm_get_real8_by_key('solute_solvent', Rcpw, rcpw_default)) then
        write(*,'(a)') 'solute-solvent cut-off set to default'
end if
if (box) then
if(.not. prm_get_real8_by_key('q_atom', Rcq, rcq_default_pbc)) then
        write(*,'(a)') 'q-atom cut-off set to default for PBC'
end if
else
if(.not. prm_get_real8_by_key('q_atom', Rcq, rcq_default_sph)) then
        write(*,'(a)') 'q-atom cut-off set to default for sphere'
end if
end if
if(use_LRF) then
        if(.not. prm_get_real8_by_key('lrf', rcLRF, rcLRF_default)) then
                write(*,'(a)') 'LRF cut-off set to default'
        end if
        if(RcLRF < rcpp .or. RcLRF < rcpw .or. RcLRF < rcww) then
                if (box) then
                      rcLRF = rcLRF_default_pbc
                      write(*,'(a)') 'LRF cut-off set to default for PBC'
                else
                      write(*,'(a)') &
                        '>>> ERROR; LRF cut-off must not be smaller than solute or solvent cut-offs!'
                      qcp_initialize = .false.
                end if
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

30	format ('>>> WARNING: Ignoring obsolete keyword ',a,'.')
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
            if(rexcl_i < zero) then
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
        qcp_initialize = .false.
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


if(.not. prm_open_section('files')) then
write(*,'(a)') '>>> ERROR: files section not found.'
qcp_initialize = .false.
else
if(.not. prm_get_string_by_key('topology', top_file)) then
        write(*,'(a)') '>>> ERROR: topology not specified (section files)'
        qcp_initialize = .false.
end if
write (*,60) trim(top_file)
60		format ('Topology file      = ',a)

if(.not. prm_get_string_by_key('restart', restart_file)) then
        restart = .false.
        if(need_restart) then
                write(*,'(a)') '>>> ERROR: Restart file required when temp. not given.'
                qcp_initialize = .false.
        end if
else
        restart = .true.
end if

if(restart) then
        write (*,65) trim(restart_file)
end if
65		format ('Initial coord. file= ',a)

if(.not. prm_get_string_by_key('traj_input', trj_file)) then
        write(*,'(a)') '>>> ERROR: Need trajectory file with coordinates to work on'
        qcp_initialize = .false.
else
        write (*,90) trim(trj_file)
end if
90		format ('Trajectory file    = ',a)

if(.not. prm_get_string_by_key('energy', ene_file)) then
        write(*,'(a)') '>>> ERROR: Energy file name required to write energies!'
        qcp_initialize = .false.
else
        write (*,94) trim(ene_file)
        iene_cycle = 1
end if
94		format ('Energy output file = ',a)

if(.not. prm_get_string_by_key('fep', fep_file)) then
        qcp_initialize = .false.
        write(*,'(a)') '>>> ERROR: No FEP file given, needed for QCP'
else
        write (*,95) trim(fep_file)
95		format ('FEP input file     = ',a,/)
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
                qcp_initialize = .false.
                exit
        end if
end do
end if
if(nstates == 0 .and. fep_file /= '') then
if(fep_file /= '') then
        write(*,'(a)') 'Defaulting to single FEP state.'
        nstates = 1
        lamda_tmp(1) = one
end if
end if
if(nstates > 0 ) then
if(fep_file == '') then
        write(*,'(a)') '>>> ERROR: FEP file required to use lambdas!'
        qcp_initialize = .false.
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

! no excluded groups in QCP calculation
!Added Paul Bauer 2014
use_excluded_groups = .false.
ngroups_gc = 0

!new section in *inp files to control QCP behaviour
!will only trigger if FEP file is in use -> if more than 0 states
if(nstates > 0 ) then
	if(.not. prm_open_section('QCP')) then
		write(*,'(a)') '>>> ERROR: No QCP section found'
		qcp_initialize = .false.
                use_qcp = .false.
                QCP_N = QCP_OFF
	else
		write(*,'(a)') 'Found QCP section, will use RPMD to describe atoms in Q region.'
	        use_qcp = .true.
                QCP_N = QCP_ON
! decide on printout level
               if(.not. prm_get_logical_by_key('verbose',qcp_verbose)) then
                       qcp_verbose = .false.
               else
                       if(qcp_verbose) then
                               write(*,'(a)') 'Printing more QCP information'
                       end if
               end if
               if(.not.prm_get_logical_by_key('veryverbose',qcp_veryverbose)) then
                       qcp_veryverbose = .false.
               else
                       if(qcp_veryverbose) then
                               qcp_verbose = .true.
                               write(*,'(a)') 'Printing all QCP information I can find'
                       end if
               end if
!chose which atoms should be treated as ring polymers
!important later when setting up NB list, RP will be treated different from classical
!this section can be overwritten in the FEP file
		if(.not. prm_get_string_by_key('selection', qcp_select)) then
			write(*,'(a)') 'Will default to Hydrogen atoms only treated as RP!'
			qcp_enum = QCP_HYDROGEN
		else
			call upcase(qcp_select)
			 if ((qcp_select .eq. 'HYDROGEN') .or. &
				(qcp_select .eq. 'HYD') .or. &
				(qcp_select .eq. 'H')) then
				qcp_enum = QCP_HYDROGEN
				write(*,'(a)') 'Only treat Hydrogen atoms as RP!'
			elseif ((qcp_select .eq. 'ALL') .or. &
					(qcp_select .eq. 'QATOM') .or. &
					(qcp_select .eq. 'FEP')) then
				qcp_enum = QCP_ALLATOM
				write(*,'(a)') 'Treat all Q-atoms as RP!'
			elseif (qcp_select .eq. 'INDIVIDUAL') then
				qcp_enum = QCP_FEPATOM
				write(*,'(a)') 'Will use information in FEP file to select QCP atoms'
			else
				write(*,'(a)') ' >>> ERROR: No such QCP atom selection!'
				qcp_initialize = .false.
			end if
		end if
!how large should the RP be?
!can again be overwritten in FEP file for each atom itself
		if(.not. prm_get_string_by_key('size',qcp_size_name)) then
			write(*,'(a,i6,a)') 'Will use default sizes for RP, ',QCP_SIZE_DEFAULT, ' ring beads per atom!'
			qcp_size = QCP_SIZE_DEFAULT
		else
			call upcase(qcp_size_name)
			if(qcp_size_name .eq. 'DEFAULT') then
				write(*,'(a,i6,a)') 'Will use default sizes for RP, ',QCP_SIZE_DEFAULT,' ring beads per atom!'
				qcp_size = QCP_SIZE_DEFAULT
			else if (qcp_size_name .eq. 'SMALL') then
				write(*,'(a,i6,a)') 'Will use small RP, ',QCP_SIZE_SMALL,' ring beads per atom!'
				qcp_size = QCP_SIZE_SMALL
			else if (qcp_size_name .eq. 'LARGE') then
				write(*,'(a,i6,a)') 'Will use large RP, ',QCP_SIZE_LARGE,' beads per atom!'
				qcp_size = QCP_SIZE_LARGE
                       else if (qcp_size_name .eq. 'DEFINE') then
                               yes = prm_get_integer_by_key('DEFINE',qcp_size)
                               write(*,'(a,i6,a)') 'Will use RP with ',qcp_size,' beads per atom!'
			else
				write(*,'(a)') ' >>> ERROR: No such QCP size selection!'
				qcp_initialize = .false.
			end if
		end if
! now we know how many beads, set default for qcp_level
! for bisection, this is 2**qcp_level .eq. qcp_size
                write(*,'(a,i6,a)') 'Setting bisection level to default, log_2(',qcp_size,')'
                qcp_level = int(q_log2(real(qcp_size,kind=prec)))
!number of PI steps at each calculation
		if(.not. prm_get_integer_by_key('equilibration_steps',qcp_steps(1))) then
			write(*,'(a,i6,a)') 'Will use default number of PI steps, n = ',QCP_steps_default,' for equilibration'
			qcp_steps(1) = QCP_steps_default
		else
			if(qcp_steps(1) .lt. 1) then
				write(*,'(a)') ' >>> ERROR: Can not use less than 1 PI step per classical step!'
				qcp_initialize = .false.
			end if
			write(*,'(a,i4,a)') 'Will use the following number of PI steps, n = ',qcp_steps(1),' for free particle equilibration!'
		end if
                if(.not.prm_get_integer_by_key('sampling_steps',qcp_steps(2))) then
                        write(*,'(a,i6,a)') 'Will use default number of PI steps, n = ',QCP_steps_default,' for sampling'
                        qcp_steps(2) = QCP_steps_default
                else
                        if(qcp_steps(2) .lt.1 ) then
                                write(*,'(a)') ' >>> ERROR: Can not use less than 1 PI step per classical step!'
                                qcp_initialize = .false.
                        end if
                        write(*,'(a,i4,a)') 'Will use the following number of PI steps, n = ',qcp_steps(2),' for free particle sampling!'
                end if
	end if
else
	write(*,'(a)') 'No RPMD in classical MD'
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
        read(text,*, iostat=fstat) rstpos(i)%i,rstpos(i)%x, &
                rstpos(i)%fk, rstpos(i)%ipsi
        if(fstat /= 0) then
                write(*,'(a)') '>>> ERROR: Invalid atom restraint data.'
                qcp_initialize = .false.
                exit
        end if
        write (*,122) rstpos(i)%i,rstpos(i)%x, &
                rstpos(i)%fk, rstpos(i)%ipsi
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
          qcp_initialize = .false.
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
          qcp_initialize = .false.
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
                qcp_initialize = .false.
                exit
        end if
        write (*,152) rstwal(i)     
end do
152  format (i6,1x,i6,4f8.2,i8)
end if
end if

! final thing to get all atoms from trj
yes = trj_store_mask('all')
number = trj_commit_mask()


call prm_close

end function qcp_initialize

!-----------------------------------------------------------------------
subroutine qcp_prep_coord


! local variables
integer(4)          :: i,nat3
logical             :: old_restart = .false.
TYPE(qr_vec),allocatable :: x_old(:),v_old(:)
TYPE(qr_vec)        :: old_boxlength, old_boxcentre
integer             :: headercheck,myprec
!new variables for differen tprecision restart files
TYPE(qr_vecs),allocatable :: x_single(:),v_single(:)
TYPE(qr_vecs)             :: boxl_single,boxc_single
TYPE(qr_vecd),allocatable :: x_double(:),v_double(:)
TYPE(qr_vecd)             :: boxl_double,boxc_double
#ifndef PGI
TYPE(qr_vecq),allocatable   :: x_quad(:),v_quad(:)
TYPE(qr_vecq)               :: boxl_quad,boxc_quad
#endif
if (prec .eq. singleprecision) then
myprec = -137
elseif (prec .eq. doubleprecision) then
myprec = -1337
#ifndef PGI
elseif (prec .eq. quadprecision) then
myprec = -13337
#endif
else
call die('No such precision')
end if



! --- Refresh topology coords. if needed (external restraints file)
if ( implicit_rstr_from_file .eq. 1 ) then
write (*,'(/,a,/)') 'Refreshing topology coords for restraining...'
read(12) headercheck
  if ((headercheck .ne. -137).and.(headercheck.ne.-1337).and.(headercheck.ne.-13337)) then
!old restart file without header canary
      rewind(12)
      old_restart = .true.
  end if
      read (12) nat3
      rewind(12)
        if(nat3 /= 3*natom) then
          write(*,100) nat3/3, natom
          call die('wrong number of atoms in restart file')
       end if
  if (.not.old_restart) then
     read (12) headercheck
     if (myprec .ne. headercheck) then
       write(*,*) '>>> WARNING: Using mismatched precision in restart file'
       if (headercheck .eq. -137) then
         allocate(x_single(natom))
         read (12,err=112,end=112) nat3, (x_single(i),i=1,natom)
         xtop(1:natom)%x = x_single(1:natom)%x
         xtop(1:natom)%y = x_single(1:natom)%y
         xtop(1:natom)%z = x_single(1:natom)%z
        deallocate(x_single)
       else if (headercheck .eq. -1337) then
         allocate(x_double(natom))
         read (12,err=112,end=112) nat3, (x_double(i),i=1,natom)
         xtop(1:natom)%x = x_double(1:natom)%x
         xtop(1:natom)%y = x_double(1:natom)%y
         xtop(1:natom)%z = x_double(1:natom)%z
         deallocate(x_double)
       else if (headercheck .eq. -13337) then
#ifndef PGI
         allocate(x_quad(natom))
         read (12,err=112,end=112) nat3, (x_quad(i),i=1,natom)
         xtop(1:natom)%x = x_quad(1:natom)%x
         xtop(1:natom)%y = x_quad(1:natom)%y
         xtop(1:natom)%z = x_quad(1:natom)%z
         deallocate(x_quad)
#else
         call die('Quadruple precision not supported in PGI')
#endif
       end if
     else
       read (12,err=112,end=112) nat3, (xtop(i),i=1,natom)
     end if
  else
     allocate(x_old(natom))
     read (12,err=112,end=112) nat3, (x_old(i),i=1,natom)
     xtop(1:natom)%x = x_old(1:natom)%x
     xtop(1:natom)%y = x_old(1:natom)%y
     xtop(1:natom)%z = x_old(1:natom)%z
     deallocate(x_old)
  end if

end if
old_restart =.false.
headercheck=0
!Assign restraints of kind res:atom their numerical atom numbers
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

! --- Make spherical restraining shell lists based on
!     the xtop coords.
if (.not. use_PBC) then

	if(rexcl_i > rexcl_o) then
	  call die('inner radius of restrained shell must be < exclusion radius')
	end if
	!first find atoms in shell 
	if (rexcl_i >= zero) then      !if rexcl_i is defined...
	  if (rexcl_i <= one) then   !if rexcl_i is defined as fraction of rexcl_o
	    rexcl_i = rexcl_i * rexcl_o  !translate to Angstrom
	  end if 
	  if(iuse_switch_atom == 1) then
	     call make_shell(nwat)
	  else
	     call make_shell2(nwat)
	  end if
	else
	  write (*,'(/,a,/)') 'Restrained shell not defined!'
	end if
else
	shell(:) = .false.
end if ! .not. use_PBC
! --- read restart file

call allocate_natom_arrays
if(restart) then
        ! topology routine has determined nwat, natom and allocated storage
        call centered_heading('Reading restart file','-')
        read(2) headercheck
        if ((headercheck .ne. -137).and.(headercheck.ne.-1337).and.(headercheck.ne.-13337)) then
!old restart file without header canary
          rewind(2)
          old_restart = .true.
        end if
          read (2) nat3
          rewind(2)
        if(nat3 /= 3*natom) then
                write(*,100) nat3/3, natom
100			format('>>>>> ERROR:',i5,' atoms in restart file not equal to',i5,&
                        ' in topology.')
                call die('wrong number of atoms in restart file')
        end if
        if (.not.old_restart) then
          read(2) headercheck
          if (myprec .ne. headercheck) then
            write(*,*) '>>> WARNING: Using mismatched precision in restart file'
            if (headercheck .eq. -137) then
              allocate(x_single(natom),v_single(natom))
              read (2,err=112,end=112) nat3, (x_single(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_single(i),i=1,natom)
              x(1:natom)%x = x_single(1:natom)%x
              x(1:natom)%y = x_single(1:natom)%y
              x(1:natom)%z = x_single(1:natom)%z
              v(1:natom)%x = v_single(1:natom)%x
              v(1:natom)%y = v_single(1:natom)%y
              v(1:natom)%z = v_single(1:natom)%z
              deallocate(x_single,v_single)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_single
                 read(2,err=112,end=112) boxc_single
                 boxlength%x = boxl_single%x
                 boxlength%y = boxl_single%y
                 boxlength%z = boxl_single%z
                 boxcentre%x = boxc_single%x
                 boxcentre%y = boxc_single%y
                 boxcentre%z = boxc_single%z
              end if
            else if (headercheck .eq. -1337) then
              allocate(x_double(natom),v_double(natom))
              read (2,err=112,end=112) nat3, (x_double(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_double(i),i=1,natom)
              x(1:natom)%x = x_double(1:natom)%x
              x(1:natom)%y = x_double(1:natom)%y
              x(1:natom)%z = x_double(1:natom)%z
              v(1:natom)%x = v_double(1:natom)%x
              v(1:natom)%y = v_double(1:natom)%y
              v(1:natom)%z = v_double(1:natom)%z
              deallocate(x_double,v_double)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_double
                 read(2,err=112,end=112) boxc_double
                 boxlength%x = boxl_double%x
                 boxlength%y = boxl_double%y
                 boxlength%z = boxl_double%z
                 boxcentre%x = boxc_double%x
                 boxcentre%y = boxc_double%y
                 boxcentre%z = boxc_double%z
              end if
            else if (headercheck .eq. -13337) then
#ifndef PGI
              allocate(x_quad(natom),v_quad(natom))
              read (2,err=112,end=112) nat3, (x_quad(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_quad(i),i=1,natom)
              x(1:natom)%x = x_quad(1:natom)%x
              x(1:natom)%y = x_quad(1:natom)%y
              x(1:natom)%z = x_quad(1:natom)%z
              v(1:natom)%x = v_quad(1:natom)%x
              v(1:natom)%y = v_quad(1:natom)%y
              v(1:natom)%z = v_quad(1:natom)%z
              deallocate(x_quad,v_quad)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_quad
                 read(2,err=112,end=112) boxc_quad
                 boxlength%x = boxl_quad%x
                 boxlength%y = boxl_quad%y
                 boxlength%z = boxl_quad%z
                 boxcentre%x = boxc_quad%x
                 boxcentre%y = boxc_quad%y
                 boxcentre%z = boxc_quad%z
              end if
#else
              call die('Quadruple precision not supported in PGI')
#endif
            end if
          else
              read (2,err=112,end=112) nat3, (x(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v(i),i=1,natom)
              if( use_PBC) then
                read(2,err=112,end=112) boxlength
                read(2,err=112,end=112) boxcentre
              end if
          end if
        else
          allocate(x_old(natom),v_old(natom))
          read (2,err=112,end=112) nat3, (x_old(i),i=1,natom)
          read (2,err=112,end=112) nat3, (v_old(i),i=1,natom)
          write(*,*) 'Read coordinates and velocities from previous version of qdyn'
          x(1:natom)%x = x_old(1:natom)%x
          x(1:natom)%y = x_old(1:natom)%y
          x(1:natom)%z = x_old(1:natom)%z
          v(1:natom)%x = v_old(1:natom)%x
          v(1:natom)%y = v_old(1:natom)%y
          v(1:natom)%z = v_old(1:natom)%z
          deallocate(x_old,v_old)
          if( use_PBC) then
             read(2,err=112,end=112) old_boxlength
             read(2,err=112,end=112) old_boxcentre
             write(*,*) 'Read boxlength and center from previous version of qdyn'
             boxlength%x = old_boxlength%x
             boxlength%y = old_boxlength%y
             boxlength%z = old_boxlength%z
             boxcentre%x = old_boxcentre%x
             boxcentre%y = old_boxcentre%y
             boxcentre%z = old_boxcentre%z
          end if
        end if
        write (*,'(a30,i8)')   'Total number of atoms        =',natom
        write (*,'(a30,i8,/)') 'Number of waters encountered =',nwat

        if( use_PBC) then
                write(*,*)
                write(*,'(a16,3f8.3)') 'Boxlength     =', boxlength
                write(*,'(a16,3f8.3)') 'Centre of box =', boxcentre
        end if
        !water polarisation data will be read from restart file in wat_shells
! now get the temperature from the restart file
        call qcp_temperature(Tfree)
else
        x(1:natom) = xtop(1:natom)
        Tfree = Temp0
end if

! clear iqatom atom array
iqatom(:) = 0

return
#if defined(USE_MPI)
112 call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
112 stop 'Aborting due to errors reading restart file.'
#endif

end subroutine qcp_prep_coord

subroutine qcp_temperature(qcp_T)
! gets temperature from restart file for QCP
! simplefied from normal temp as only free atoms count
! arguments
real(kind=prec)         :: qcp_T
! locals
integer                 :: i
Temp = zero
do i = 1, nat_solute
if (use_PBC .or. ( (.not. use_PBC) .and. (.not. excl(i)))) then
qcp_T = qcp_T + 0.5_prec*iaclib(iac(i))%mass*(qvec_square(v(i)))
end if
end do
do i = nat_solute+1,natom
if (use_PBC .or. ( (.not. use_PBC) .and. (.not. excl(i)))) then
qcp_T = qcp_T + 0.5_prec*iaclib(iac(i))%mass*(qvec_square(v(i)))
end if
end do

qcp_T = 2.0_prec*qcp_T/Boltz/real(Ndegfree,kind=prec)

end subroutine qcp_temperature

subroutine qcp_work
! this does the work when doing post processing
! needs to do the following for each coordinate specified in trajectory file
! read in and replace x
! update nb lists
! call qcp_run
! write out energies (and also somehow ring polymer coordinates)
! locals
integer                         :: frame,nframes
logical                         :: yes

frame = 0
nat3 = natom*3
!only do file handling on master node
if (nodeid.eq.0) then
yes =  trj_open(trj_file)
if(yes) then
! only do this if we actually manage to work on the file, otherwise stop here
!call centered_heading('Beginning work on trajectory','-')
do while(trj_read(x))
! first read is to get file length
! put coordinates straight into the x array
frame = frame + 1
end do
! and rewind the whole thing back to the start
yes = trj_seek(1)
nframes = frame
end if
else
! we are dead, can not open dcd file
call die('Could not open trajectory file')
end if

! now we know the length of the file, we can iterate over the individual frames and Bcast stuff to nodes
#ifdef USE_MPI
call MPI_Bcast(yes,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
if(ierr.ne.0) call die('QPI file read Bcase')
call MPI_Bcast(nframes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(ierr.ne.0) call die('QPI Bcast nframes')
#endif

do frame = 1, nframes
! file reading only on node 0
if(nodeid.eq.0) yes = trj_read(x)

#ifdef USE_MPI
call MPI_Bcast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if(ierr.ne.0) call die ('QPI Bcast x')
#endif
call make_pair_lists(Rcq,Rcq**2,RcLRF**2,Rcpp**2,Rcpw**2,Rcww**2)
! first call to pot_energy to get classical energies of the given configuration
call pot_energy(E,EQ,.true.)
! and call qcp_run now for qcp_energies
call qcp_run(Tfree,E,EQ)
! write stuff to nice data structure

if (nodeid.eq.0) then
call qatom_savetowrite(EQ,1)
! sasve to file
call put_ene(11,EQ_save,OFFD,ene_header%arrays,nstates)
end if
end do

if (nodeid.eq.0) call trj_close


end subroutine qcp_work

subroutine qcp_startup
! here we go again

end subroutine qcp_startup

end module QCP