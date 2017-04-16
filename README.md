Qsource development.
================================================================================

This is a repository to hold clean and readable Q code.  
The code is in the process of being reorganized and cleaned-up
according to best coding practices.  

We are documenting the code using doxygen and ford.   
To generate the doxygen documentation do:  

    cd docs/developers  
    doxygen DoxygenConfigFortran  


To generate the ford docs:  
    
    cd docs/developers  
    ford -d ../../src -o devdocs forddocs.md  


Some basic development rules somewhat similar to those of the
**gromacs** developers team follow.


Coding Standards
--------------------------------------------------------------------------------

The **GROMACS** developers have identified the following important main
points for taking into account when organizing a molecular dynamics
code.  

1.  **Code formatting** - how to indent code, how to start and end subroutines
    etc.  
2.  **Code constructs** - argument order, return values, encapsulation
    using abstract data types  
3.  **Interfaces** - the Application Programming Interface should say it
    all  
4.  **Comments in code** - comments in code that ford can use  
5.  **Compilation** - using different hardware  
6.  **Allowed Fortran Features**  
7.  **Error Handling**  
8.  **Benchmarking**  
9.  **Accuracy Testing**  


1. **Code formatting**  
--------------------------------------------------------------------------------

* No tabs, spaces only.  
* Two spaces for indentation of each level.  
* No more than 80 characters to allow for easy code visualization
  across editors and screens. Specially important for easy mobile
  device coding, reading.  
* No trailing whitespaces. They are not seen by default in most
  editors but still count as changes in git.  
* Column 36 for : on variable declaration.  


2. **Code constructs**  
--------------------------------------------------------------------------------


3. **Interfaces**  
--------------------------------------------------------------------------------


4. **Comments in code**  
--------------------------------------------------------------------------------


5. **Compilation**  
--------------------------------------------------------------------------------

5.1  gfortran from gcc
----------------------
At the moment Q doesn't include any of the features of the fortran
2008 standard. The fortran 2008 standard is the latest fortran
standard. Work is being done on the fortran 2015 standard but it is
not yet released. Many features of fortran 2008 are already included
in the GNU-Compiler-Collection gfortran.

- For compiling Q a relatively new version of gfortran is recommended,
at least  gfortran 4.8. The  idea is  to make Q  an example of  a code
which integrates the latest niceties of modern *FORTRAN*.

- **Fortran 2008**(New Features)  
Coarrays for parallel computing.  
Possibility to use submodules, and submodules of submodules. Of help
for very large programs.  
*do concurrent*  


6. **Allowed Fortran Features**  
--------------------------------------------------------------------------------



7. **Error Handling**  
--------------------------------------------------------------------------------


8. **Benchmarking**  
--------------------------------------------------------------------------------


9. **Accuracy Testing**  
--------------------------------------------------------------------------------





