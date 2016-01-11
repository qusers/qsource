Qsource development.
================================================================================

This is a repository to hold clean Q code.
The code is in the process of being reorganized and cleaned-up
according to best coding practices.

Some basic rules taken from the gromacs developers team follow.

Coding Standards
--------------------------------------------------------------------------------

The GROMACS developers have identified the following important main
points for taking into account into organizing a molecular dynamics
code.  

1. **Code formatting** - how to indent code, where to put braces etc.  
2. **Code constructs** - argument order, return values  
3. Encapsulation using abstract data types  
4. **Interfaces** - the Application Programming Interface should say it
   all  
5. **Comments in code** - need we say more?  
6. **Compilation** - now and later  
7. **Allowed Fortran Features**  
8. **Error Handling**  



1. **Code formatting**  
--------------------------------------------------------------------------------

* No tabs, spaces only.  
* Two spaces for indentation of each level.  
* No more than 80 characters to allow for easy code visualization
  across editors and screens, specially important for easy cellphone
  coding, reading.  
* No trailing whitespaces. They are not seen by default in most
  editors but still count as changes in git.  

