A test using a derivative of 2,4-quinazolinediamine.
--------------------------------------------------------------------------------

A derivative of 2,4-quinazolinediamine is the molecule rendered in the
main Q-logo.   The interest in  this molecule  came about due  to it's
relation to small molecules similar to quinine, and their possible use
as anti-malarial drugs.

We are  using here various sets  of examples to try  and simulate this
molecule  in  spherical  boundary  conditions.  We  have  picked  some
selected  conformers  which have  been  generated  using the  Avogadro
molecular visualization package.

##Avogadro can generate sets of conformers easily.

In this example 6-[(N-methylanilino)methyl]quinazoline-2,4-diamine has
been      downloaded       in      the      sdf       format      from
https://pubchem.ncbi.nlm.nih.gov/compound/487413

Then         open         in         avogadro         and         then
extensions/molecular_mechanics/setup_force_field Select GAFF, steepest
descent, 2000 steps.

After the optimization ends you can go to:
extensions/molecular_mechanics/confomer_search

This will generate conformers using default parameters.

After they are generated you can see them by going to:
extensions/animations

From this you can  find the one that looks most similar  to the one in
the old q-manual and documentation.

##smooth.sh

The main script  to run the test  is smooth.sh It's intended  to run a
smooth heating and  release of restraints.  Analysis scripts  in R are
not yet ready but will be available in a corresponding folder.


##Pymol notes
For enhancing sphere quality in pymol.
set sphere_quality, 2


##References
J. Comput. Aided Mol. Des. 1998 12, 119-131.  
Computation    of    affinity     and    selectivity:    binding    of
2,4-diaminopteridine   and    2,4-diaminoquinazoline   inhibitors   to
dihydrofolate reductases.  
Marelius J, Graffner-Nordberg M, Hansson T, Hallberg A, Aqvist J.  

