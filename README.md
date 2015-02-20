Qsource development.
--------------------------------------------------------------------------------

Creating  a  hierarchy  of  "folders"  corresponding  to  branches  of
development.

The  idea is  to  migrate the  already  existing personal  development
branches to a development hierarchy that will look something like:

	development/alex
	development/beat
    development/esguerra	
	development/irek
	development/paul

New folders for development branches should be added under development.
Say, if  I (username:esguerra) wanted to  create a new branch  for the
development hierarchy I could do:

##Create new development branch which can be latter merge to master  
    git branch development/esguerra
    git checkout development/esguerra
    git add .
    git commit -a
    git push --set-upstream origin development/esguerra

##To go back to main master branch.  
    git checkout master

When you want to clone only one branch and not the full repository use:   
    git clone https://github.com/qusers/qsource.git --branch development/esguerra --single-branch esguerra


