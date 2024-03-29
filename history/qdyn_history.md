# **qdyn** version history


## Release 5.05
- Added feature to select q-atoms in the fep file under [atoms]
  section by defining residue number/s or name/s., e.g.  

    [atoms]
    residue 23
    residue 25
    or
    [atoms]
    all ALA
    all PHE

- Fixed a problem with the parallel version when using periodic
  boundary condition in conjunction with LRF. The problem affected
  water molecules moving across the border.


## Beta version 5.04
- Added LRF to periodic boundary condition.

- Fixed a bug where Q-lib charges were used in certain qp interactions when
  qq_lib_charges was on

- Added softcore.

- Added feature to reference atoms as <residue number:atom number>
  Atom number is the position of the atom in the residue as defined in 
  the lib file.

- Fixed a bug causing certain atoms in the restrained shell not to be
  restrained.
     
- Fixed parallel algorithm for periodic boundary and constant pressure.

## Beta version 5.01

- Parallel algorithm based on MPI implemented.

- Moved definition of solvent restrain shell from topology to Qdyn input file.
  The shell can thus be changed without remaking the topology.  
  Furthermore the default is changed to NO shell. See manual for details.
       
- Periodical boundary condition implemented. LRF and constant preassure options
  included. See manual for details. 
 
- Electrostatic interactions between q-q-pairs can be scaled down. 
  They are defined in the fep-file:  

    [el_scale]
    qatom1 qatom2 scale_factor

## Release 4.26, 01-05-08
- Fixed a bug in shake initialisation which previously caused Qdyn to
  crash when shake was enabled and there were multiple single-atom
  molecules in the topology.
  Thanks to Dr. Igor I. Baskin for pointing out the problem and
  suggesting a fix.

## Release 4.25, 00-11-09
- Fixed a problem with the Q-atom non-bond list, which prevoiusly
  could contain duplicate entries and incorrect 1-4 pair assignments.

## Release 4.24, 00-09-17
- Fixed problem with redefining multiple torsions & impropers which
  led to a run-time error.

## Release 4.23, 00-08-28
- Corrected two new problems:
1. If solvent (water) section was missing in the input file, the
   radial solvent restraining boundary attraction term "morse width"
   was set to 0. This problem came with version 4.18.
2. Q-atom - SPC water LJ interactions were wrong. This problem
   probable occured around 4.18

## Beta version 4.22, 00-08-25
- Will now silently accept input/FEP files not ending in newline.

## Beta version 4.20, 00-08-21
- Some syntax changes to compile with the Portland Group Linux F90 compiler.

## Beta version 4.19, 00-07-14
- Corrected removal of redefined improper torsions from topology with
  special code for the CHARMM force field by P. Vagedes. 

## Beta version 4.18, 00-07-10
    Changed non-bond routines for TIP3P water to handle any 3-atom solvent
    (the same loop unfolding tricks can still be used). 
    
    Elastic wall restraints now have a boundary attraction term, like the
    solvent boundary. The syntax in the section wall_restraints changed to
    [wall_restraints]
    atom_1 atom_2 radius harmonic_fc Morse_width Morse_depth H_flag

    Solvent boundary attraction range can now be set explicitly in the input
    file using the keyword morse_width. If omitted, an empirical function is
    used as before.

    Target solvent sphere radius may be overridden in the input file using the
    keyword radius in the water section.

    Solvent radial restraints are applied to the atom designated as switching 
    atom in each solvent charge group

## Beta version 4.17, 00-06-13
    Corrected IRIX compilation problem caused by improper syntax.

## Beta version 4.16, 00-03-24
    Fixed bug in shake, causing problems when only solvent, not solute shaken.

## Beta version 4.15, 00-03-02
    New more flexible shake code. Shake can be enabled for hydrogens only 
    using the line
        shake_hydrogens on
    in the MD section of the input file. 
    The new shake routines are adapted for shaking any kind of solvent 
    molecule, not just water. The special code for shaking of water 
    has been removed.
    ***IMPORTANT***
    IT IS NECESSARY TO ADD AN EXTRA BOND TO THE WATER LIBRARY ENTRY 
    to shake the angle of the water molecule (by constraining the H-H
    distance). The correct entry for GROMOS87 SPC water is
    {WAT}
    [info]
        SYBYLtype  GROUP
    [atoms]
         1 OW     O.SPC       -0.820 
         2 H1     H            0.410
         3 H2     H            0.410
    [bonds]
         1    2
         1    3
         2    3  !note extra H-H bond for shake
    [charge_groups]
          1    2    3
    
    Add also the bond parameters to the parameter file:
    H       H            0.000   1.63328311  ! H2O shake
    Set the bond force constant to zero to avoid unnecessary O-H-H angles
    being made.

    Working further towards handling any kind of solvent, modifications have 
    been made to the topology: charge groups are generated also for solvent,
    torsions and impropers in solvent can be generated and energies for those
    are book-kept separately.

    The output of temperature can now be limited by adding e.g
    temperature 10
    to the section intervals, to print the temperature every 10 steps.
    
## Beta version 4.14, 00-02-27
    Modified LRF code, but no changes in function. Data for water molecules 
    no longer duplicated.

## Beta version 4.12, 00-01-13
    Changed the mask module: the property keyword all is no longer used but
    is the default.
    The new syntax for the [trajectory_atoms] section is:
    [properties] [residue] [first [last]]
      properties are one or more of:
        solute: only solute atoms
        excluded: atoms outside the simulation sphere
        restrained: atoms outside the simulation sphere and in the restrained shell
        heavy: heavy (non-hydrogen) atoms
      Each property may be negated by prepending not.
      If no property is given, all atoms are selected.
      residue: first and last are residue numbers, not atom numbers
      first: first atom/residue of a sequence.
      last: last atom/residue of a sequence.
      If first and last are omitted, first=1 and last=number of atoms/residues
      If last is omitted, last=first.

## Beta version 4.11, 00-01-12
    Changed the mask module to allow multiple masks (for Qrms - different
    mask in trajectory and RMS calc.)
    Fixed compilation problems on IRIX.
    Added possibility to select residue numbers in mask specification using
    the property keyword residue. Also made the last atom/residue specification optional.
    The new sytax is:
    properties [residue] [first [last]]
      properties are one or more of:
        all: all atoms
        solute: only solute atoms
        excluded: atoms outside the simulation sphere
        restrained: atoms outside the simulation sphere and in the restrained shell
        heavy: heavy (non-hydrogen) atoms
      Each property may be negated by prepending not.
      residue: first and last are residue numbers, not atom numbers
      first: first atom/residue of a sequence.
      last: last atom/residue of a sequence.
      If first and last are omitted, first=1 and last=number of atoms/residues
      If last is omitted, last=first.

## Beta version 4.10, 00-01-11
    DCD trajectory format introduced. This is the charmm trajectory format 
    which is compatible with trajectory visualisation/animation programs
    like VMD. A subset of atoms can be written to the trajectory by using an
    atom mask given in the new section [trajectory_atoms] in the input file.
    The following syntax is used to specify the set of atoms to include:
        properties [first last]
        properties are one or more of:
            all: all atoms
            solute: only solute atoms
            excluded: atoms outside the simulation sphere
            restrained: atoms outside the simulation sphere and in the 
                        restrained shell
            heavy: heavy (non-hydrogen) atoms
        Each property may be negated by prepending not.
        first and last designate the first and last atom of a sequence.
        The default is first=1 and last=number of atoms
    Example:
    [trajectory_atoms]
    heavy solute not restrained
    all 1897 1940

    This will save disk space in both by writing coordinates in single 
    precision and by the option to save only a subset of atoms.

## Beta version 4.08, 00-01-10
    Allow file names in input file to be up to 200 characters. (Only new style
    input files.)

## Beta version 4.03, 99-12-17
    Problem with shakeing of water solved. 
    The heading 'Removing redefined interactions' is not printed if no
    interactions are redefined.

## Beta version 4.00, 99-12-13
    Uses the new topology standard generated by Qprep4 (including solvent, 
    sphere, exclusion, ...) The following keywords in the input file are now 
    obsolete :
    centre, radius, pack, model, water

## Release 3.73, 99-10-08
    Fixed invalid references to unallocated arrays in code for monitor groups 
    and simulation startup that caused Digital UNIX version to crash when running
    with no FEP file.

## Release 3.72, 99-09-15
    Updates to output for better readability & easier understanding.

## Release 3.71, 99-08-31
    Fixed problem with enabling polarisation restraints in old-style input
    file.
    Added proper heading for table with Q atom type parameters for arithmetic
    combination rule force fields.

## Release 3.70, 99-07-05
    Added support for harmonic potentials as Q-bonds. See FEP file example:
    [bond_types]
    !Morse: D_e, alpha, r_0
    !Harmonic: force_k, r_0
    1  85.0 2.00 1.61 !This is a Morse type
    2 600.0      1.61 !This is a harmonic type

    Note: Mixing Morse and harmonic for a single bond is not advisable,
    although possible.
    Note: Harmonic bonds cannot be coupled to angles. 

## Release 3.68, 99-07-02
    Fixed a bug that caused invalid error messages when reading Qangle-Morse 
    couplings.

## Release 3.67, 99-06-11
    It is now possible to override the automatic calibration of the water radius 
    (based on number & positions of solute & solvent atoms) by entering the
    water radius with negative sign in the input file.

## Release 3.66, 99-04-21
    An offset to add to topology atom numbers can now be specified in the FEP
    file:
    [FEP]
    offset  1896        !atom number offset - add this to all topology atom 
                        !numbers
    offset_residue 189  !calculate the offset to that the first Q-atom in the
                        !first atom of fragment 189
    offset_name MTX     !calc. offset so 1st Q-atom is the first atom of the first
                        !fragment named MTX
    Of course only one of these should be used in a FEP file.
    I think the last one is very useful - one can use the same FEP file for 
    both water and protein ligand perturbation simulations.
    
    Fixed problem with too many molecules in outermost water shell which in 
    an unusual case led to a memory fault.

## Release 3.65, 99-03-28
    Fixed bug in special non-bonded exclusions. (The whole qq non-bond list
    was not checked to remove excluded pairs.
    Introduced Isabella Feierberg's code to monitor non-bonded interactions
    between selected groups of atoms. The groups are defined in the FEP file:
    [monitor_atom_groups]
    1  2   3   !group 1 contains atoms 1 2 and 3 in the _TOPOLOGY_
    4  5       !group 2 contains atoms 4 and 5
    7          !group 3 contains atom 7
    [group_pairs]
    1  2        !monitor sum of nonbonded interactions between group 1 and 2
    1  3        !and between group 1 and 3

## Release 3.64, 99-03-24
    Fixed bug in removing bonds/angles/... with code=0 due to Q-bonds/...
    overriding topology. This was only a problem if the last bond/... was
    replaced.

## Release 3.63, 99-03-15
    Urey-Bradley 1-3 potential for complete charmm conformity.
    Water bond and angle parameters and charge now read from topology, hard-
    coded constants for SPC water used only for backward compatibility.
    Cleaned up data types for q-angles.

## Release 3.62, 99-03-08
    Fixed bug in atom position restraint routine.
    Changed initial velocity routine back to setting velocities for excluded atoms.

## Release 3.61, 99-03-02
    Fixed problem with Q-atom - water interactions and arithmetic LJ 
    combination rule. Fixed problems with LJ interactions of TIP3P water 
    hydrogens using geometric comb. rule. There are now separate Q-atom-
    water and water-water nonbond routines for SPC and TIP3P waters.

## Release 3.60, 99-02-25
    Dynamic water polarisation restraint data now written to & read from 
    restart files. 
    Fixed case conversion problem which led to trouble with TIP3P water.

## Release 3.58, 99-02-18
    Fixed the bug Peter Vagedes found in the force calculation in 
    nonbon2_pp (solute-solute nonbond routine for arithmetic combination rule)

## Release 3.57, 99-02-15
    Fixed problem with water exclusion / running with water sphere in topology
    Fixed printing of changing charges when reading old-style FEP file.
    Fixed "sphere too small" problem solvating with solvest sphere.

## Release 3.56, 99-02-03
    Added the keyword exclude_bonded in the section sphere of the input file
    to control the elimination of bonded interaction between excluded atoms.
    Added a cut-off radius for LRF using the keyword LRF in the cut-offs 
    section.

## Interim version 3.55, 99-02-01
    Split nonbond_qp to two simplier routines, avoiding branches inside
    the loops.

## Interim version 3.54, 99-01-30
    Initial velocity is set to zero for excluded atoms. This should keep
    excluded waters from escaping into space.
    Initial velocities for other atoms should be identical - the same
    number of random numbers are drawn.

## Interim version 3.53, 99-01-29
    Bonded interaction between atoms that are all excluded are now eliminated 
    at startup. The number of bonds, angles, torsions and impropers to 
    evaluate is thus reduced significantly for systems where part of the 
    solute is excluded. 3% speedup when 1000 out of 3000 atoms excluded.
    
## Release 3.52, 99-01-28
    Added listing of all effective Q-atom charges when reading (new) FEP file.
    Exclusion of waters now only outside rwat+2� to avoid excluding waters from a 
    previous simulation

## Release 3.51, 99-01-24
    Fixed problems with pure water simulation: 
    1) topology reading. Lists with zero elements are not read.
    2) NaN in water polarisation for water molecule at centre of sphere.
    Atom type names are now stored in the topology.

## Release 3.50, 98-11-27
    Updated water polarisation restraints. Refecence point is now
    oxygen atom, not dipole centre

## Release 3.49, 98-11-10
    Modified FEP file format (new-style FEP files only). Q-atom types defined in the section atom_types now have alphanumeric names (max 8 characters) instead of numbers. These names are used in the section change_atoms. Example:
        [atom_types]
        !Type Name    Ai     Bi      Ci    ai    Ai(1-4) Bi(1-4)  Mass 
        OZ        600.00  23.25    0    0      600.00  23.25   15.999   
        OE1       550.00  23.01   65  1.581    550.00  23.25   15.999   
        [change_atoms]
         9   OE1   OZ    
         
    Fixed problem in prmfile with the use of the null function which is a
    F95 feature and did not work on the IRIX 6.4/MIPSPRO 7.2  compiler.

## Release 3.48b, 98-10-07
    Fixed bug in solvation routine which caused crashes when replicating a periodic water box.

## Release 3.48, 98-09-25
    Complete rewrite of parameter file reading module to improve portability

## Release 3.47, 98-09-06
    Fixed bug caused by pair-counting array overrun.
    Rewrote Q-atom - non-Q-atom nonbond routines, separating Q-protein and
    Q-water, introduced optimisations in new Q-water routine.
    Qdyn3.47i is up to 30% faster than Qdyn2.1!

## Release 3.44, 98-09-01
    
## Release 3.43, 98-08-28
    Special exclusion between q-atom and non-q-atom or between two non-Q-
    atoms are now implemented.

## Release 3.41, 98-08-25
    New input file format with headings. Many entries are optional and have
    default values. See the new manual!

## Release 3.4, 98-07-28
    New FEP file format with headings like in the new libraries and parameter files. Fixed unimportant bug related to exclusion of bonded/angled atoms from non-bonded interactions. Water is now recognised in the topology (only after all solute atoms) and water interaction energies will be book-kept correctly as opposed to previous versions where water in the topology was consideded as solute. It is thus no longer necessary to make an extra topology with xtal waters to use only for the first step. Waters in the topology can now be excluded like protein atoms if outside the sphere. Shell restraints don't apply to waters, though. The earlier version would not exclude water molecules, this could lead to excessive radial restraining forces on molecules far outside the simulation sphere.

## Qdum3 input checker program, 98-06-18
    Compilation directives to omit the whole dynamics loop if DUM is defined.

## Release 3.30, 98-05-19
    Changed format of offdiagonals in energy file to be same as old Qdyn/Qfep.
    Introduced feature to restrain sequence of atom to geometrical centre of the atom group
    instead of only to each atom's topology co-ordinates.

## Release 3.29, 98-04-28
    Write a final trajectory frame after dynamics loop, e.g. writing trajectory every 500
    steps for 25000 steps now gives 50 frames, not 49.
    Removed obsolete open statement for topology file, shortened lines of MD module to 
    below 100 characters.

## Release 3.28, 98-04-24
    Topology module updated, formatting of topology file changed to make sure 
    there's space between values.

## Release 3.27, 98-01-08
    Revised water polarisation restraints to Johan's latest model with shell widths
    dr, 2dr, 3dr

## Release 3.26, 97-12-11
    Decreased the number of protein-protein and protein-water interaction pairs that
    need to be evaluated by ~10% and ~25%, respectively, by checking if all interaction
    parameters are zero when making lists. (Only nbpplist, nbpplist_lrf, nppwlist, 
    nbpwlist_lrf.)

## Release 3.25, 97-12-10
    Separated out all topology-related stuff shared with qprep into a new module
    TOPO, now used by both Qdyn and Qprep.

## Release 3.22, 97-12-01
    Fixed bug which caused wrong number of degrees of freedom and wrong temperature 
    when shakeing the protein.

## Release 3.21, 97-11-28
    Now uses module NRGY for Q-atom energies and writes the energy file using 
    derived types.

## Release 3.20, 97-11-27
    Restuctured the handling of energies by inroducing a module NRGY with derived 
    types.

## Intermediate version 3.19, 97-11-27
    Changed water polarisation restraints according to J�'s new model.
    Intruduced subroutine water_shells which is called during initialisation
    Cleaned up qdyn.f90 by using 'contains' rather than a common block

## Release 3.09, 97-11-24
    Removed silly message when iseed or T_Maxwell = 0, added check if tau_t < dt

## Release 3.08, 97-10-11
    Mostly cosmetic changes to output format

## Release 3.07, 97-09-25
    Sigmoid function giving Water boundary outward Morse potential depth added.
    The function is used when Dwmz is not given in the input file

## Release 3.06[i], 97-09-22
    Water boundary outward morse potential depth can now be set in input file
    but is optional. All water-sphere related paramaters except the radius are now
    optional. Random no. seed and initial temperatures also made optional.
    Corrected problem with allocating 1-4 neighbour and neighbour
    exclusion lists (always nat_pro elements!)

## Release 3.05[i], 97-09-18
    Introduced inv_sqrt and conditional compilation
    use /define:"INVROOT" to get 'i' version with inv_sqrt

## Release 3.03, 97-09-01
- Fixed setting iqatom to 0 for waters.
- Fixed calculating total charge (must be done after get_fep).
- Introduced 'special exclusions'.
- Changed the structure of the water-water non-bond list.

## md.0504.f90
- New energy summary format.
- changed version number to 2.95i (3.0i will be final).
- nonbond_pp loop unfolded (2 parallell parts)
- removed divisions from angle/torsion/improper
- Uses inv_sqrt functions in:  
  nonbond_pp  
  torsion  
  angle  
  improper  
- No speed improvement on Pentium w 256kb cache

## md.0503d.f90
- OK with mdtest.inp and mdtest2.inp
- 20% faster than previous versions 
- 250 step profiling, same speedup on most routines
- all variable kinds have been set to the smallest possible kind
- default kinds can easily be changed
- bonds and angles have data structures and lists are dynamically allocated

## md.0503.f90   
- OK with mdtest.inp and mdtest2.inp
- fixed allocation / reallocation spaghetti
- solvate now always call watparm
- calls prep_coord before get_fep

