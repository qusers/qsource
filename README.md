Qsource development.
================================================================================

This is a repository to hold clean and readable Q code.
The code is in the process of being reorganized and cleaned-up
according to best coding practices.

We are documenting the code using doxygen. To generate the doxygen
documentation do:  

    cd documentation/manuals/developers
	doxygen DoxygenConfigFortran


Some basic development rules taken from the **gromacs** developers team follow.


Coding Standards
--------------------------------------------------------------------------------

The GROMACS developers have identified the following important main
points for taking into account when organizing a molecular dynamics
code.  

1. **Code formatting** - how to indent code, where to put braces etc.  
2. **Code constructs** - argument order, return values, encapsulation
   using abstract data types  
3. **Interfaces** - the Application Programming Interface should say it
   all  
4. **Comments in code** - comments in code that doxygen can use  
5. **Compilation** - using different hardware  
6. **Allowed Fortran Features**  
7. **Error Handling**  



1. **Code formatting**  
--------------------------------------------------------------------------------

* No tabs, spaces only.  
* Two spaces for indentation of each level.  
* No more than 80 characters to allow for easy code visualization
  across editors and screens, specially important for easy cellphone
  coding, reading.  
* No trailing whitespaces. They are not seen by default in most
  editors but still count as changes in git.  


2. **Code constructs**  
--------------------------------------------------------------------------------


