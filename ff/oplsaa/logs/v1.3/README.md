### Modifications

Some angles and torsions in residues GLH and ASH that contain the atom type ODE 
had wrong values (compared to ffld_server output and GROMACS library).
Additionally, some duplicates were removed from the improper section.

Example:  
```
qoplsaa.prm
-------------------------------

O1  C2  ODE    160.000  126.000
CT  C2  ODE    140.000  117.000
...
CT  CT  C2  ODE     0.583  2.000  180.000 1 
```

```
ash.prm
-------------------------------

ash_OD2  ash_CG  ash_OD1    160.00  121.0000 
ash_CB   ash_CG  ash_OD1    140.00  108.0000 
...
ash_CA  ash_CB  ash_CG  ash_OD1     0.5000  -1.0000    0.0000  1.0000 
ash_CA  ash_CB  ash_CG  ash_OD1     0.2730  -2.0000  180.0000  1.0000 
ash_CA  ash_CB  ash_CG  ash_OD1     0.2250   3.0000    0.0000  1.0000 
```

ash.ffld14 was obtained from ash.pdb:  
```ffld_server -print_parameters -version 14 -ipdb ash.pdb > ash.ffld14```

ash.prm and ash.lib were obtained from ash.ffld14:  
```q_ffld2q.py ash.ffld14 ash -s ash.pdb```


