Qsource testing version.
=======================

This is a testing version which Masoud has been using for comparing profiles of the
newish code and this older version in this branch of the repository.

This is another dead branch just like the stable branch and it's meant to be
used as a way to have the old code handy if one wants to do back-on-time comparisons.

If more branches from the past appear they will be included in this testing barnch.

For now there is only one folder in testing, that is:

testing/masoud

So, new folders should be added under testing, say, if I wanted to create a new
dead branch I could do:

#To create new dead branch (not being developed) for testing agaist development.
git branch testing/esguerra
git checkout esguerra
git add .
git commit -a
git push

#To go back to main development branch.
git checkout merge


When you want to clone only one branch and not the full repository use:

    git clone https://github.com/qusers/qsource.git --branch testing/masoud --single-branch qmasoud

