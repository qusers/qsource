Qsource testing version.
--------------------------------------------------------------------------------

This is yet another testing version.  This one will be kept as a stale
branch with  the final version  of the  NEW MPI implementation  by Ake
Sandgren. Ake also took care of cleaning  parts of the code such as an
odd  call to  acosd which  only  use to  work on  old IBM  proprietary
compilers.

So there are three testing versions now:

    testing/masoud
    testing/esguerra
    testing/akesandgren

New folders  for the testing  branches should be added  under testing,
say, if  I (username:esguerra) wanted  to create  a new dead  branch I
could do:

##Create new testing branch (ideally not being developed) to test against master.  
    git branch testing/esguerra
    git checkout esguerra
    git add .
    git commit -a
    git push --set-upstream origin testing/esguerra

##To go back to main development (master) branch.  
    git checkout master

To get this version to compile intel fortran is needed.

### At CSB
Using the Makefile from Ake  intel fortran is necessary, else gfortran
does not like the declaration of a signal function in qdyn.f90
This is  overcome by  simply changing the  signal function  from where
it's invoked and taking it out of the interface declaration. 

To get intel fortran into your environment in csb do:
    source /home/apps/intel/bin/ifortvars.sh intel64

And to get intel MPI (but not mpiifort) into the environment:
    /home/apps/intel/impi/5.0.0.028/intel64/bin/mpivars.sh  intel64

### At triolith

module load intel/15.0.1.133
module load impi/5.0.2.044

When you want to clone only one branch and not the full repository use:
    git clone https://github.com/qusers/qsource.git --branch testing/akesandgren --single-branch qakesandgren



