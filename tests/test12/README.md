n-butane
================================================================================

The aim of this test is to compare Figure 8 of Beglov-Roux, JCP,
1994, 100, 9050-9063 with results obtained with Q.


Pymol scripts
--------------------------------------------------------------------------------

And independent pymol folder contains a script to visualize the
results. One needs just to use:  

    pymol heatrelaxequi.pml
	
And this will load the trajectories.


R scripts
--------------------------------------------------------------------------------





Benchmarks
--------------------------------------------------------------------------------


|  Machine     | Compiler    | Comp. time (min) | Sim. time (ns) | Num Proc. |    Date    |
|:-------------|:-----------:|:----------------:|:--------------:|:---------:|:----------:|
| csb          | ifort       | 79.00            |      0.20      |   8       | 2014-05-22 |
| triolith     | ifort       | 00.00            |      0.20      |   8       | 2014-XX-XX |
| tintin       | ifort       | 00.00            |      0.20      |   8       | 2014-XX-XX |
| csb          | gcc         | 00.00            |      0.20      |   8       | 2014-XX-XX |
| triolith     | gcc         | 00.00            |      0.20      |   8       | 2014-XX-XX |
| tintin       | gcc         | 00.00            |      0.20      |   8       | 2014-XX-XX |


TODO
--------------------------------------------------------------------------------

The following to-do list highlights what needs to be done for expanding the benchmark into
something userful.

- [ ] Make a script that will send runs of 1, 2, 3, 4, 5ns runs to different cluster nodes
      at the same time.
- [ ] Make number of processors gradient script.
- [x] Make a markdown document describing each test.
- [ ] Make a general R script for plotting and making statistics with the benchmarks.
