Qsource testing version.
=======================

This is yet another testing version.
This one will be kept as a stale branch with the final version
of the NEW MPI implementation by Ake Sandgren. Ake also took care
of cleaning parts of the code such as an odd call to acosd which only
use to work on old IBM proprietary compilers.

So there are three testing versions now:

    testing/masoud
    testing/esguerra
    testing/akesandgren

New folders for the testing branches should be added under testing,
say, if I (username:esguerra) wanted to create a new dead branch I could do:

##Create new testing branch (ideally not being developed) to test against master.  
    git branch testing/esguerra
    git checkout esguerra
    git add .
    git commit -a
    git push --set-upstream origin testing/esguerra

##To go back to main development (master) branch.  
    git checkout master

