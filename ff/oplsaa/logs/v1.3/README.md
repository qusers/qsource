### Modifications

Some angles and torsions in residues GLH and ASH that contain the atom type ODE 
had values that did not agree with ffld_server output.

Additionally, one charge group in the library in THR was defined with numbers
instead of atom names for some reason.


Examples of wrong parameter values:  
```
qoplsaa.prm
-------------------------------

O1  C2  ODE    160.000  126.000
CT  C2  ODE    140.000  117.000

...

CT  CT  C2  ODE     0.583  2.000  180.000 1 
O1  C2  ODE HO      2.45   2.000  180.000 1.000   !102
CT  C2  ODE HO      1.5   -1.000    0.000 1.000   !103      
CT  C2  ODE HO      2.45   2.000  180.000 1.000   !103 
```

```
ash.prm
-------------------------------

ash_OD1  ash_CG  ash_OD2    160.00  121.0000 
ash_CB   ash_CG  ash_OD2    140.00  108.0000 

...

ash_CA    ash_CB  ash_CG    ash_OD2   0.5000  -1.0000     0.0000   1.0000
ash_CA    ash_CB  ash_CG    ash_OD2   0.2730  -2.0000   180.0000   1.0000
ash_CA    ash_CB  ash_CG    ash_OD2   0.2250   3.0000     0.0000   1.0000
ash_CB    ash_CG  ash_OD2   ash_HD2   0.7500  -1.0000     0.0000   1.0000
ash_CB    ash_CG  ash_OD2   ash_HD2   2.7500   2.0000   180.0000   1.0000
ash_OD1   ash_CG  ash_OD2   ash_HD2   2.7500   2.0000   180.0000   1.0000
```

ash.pdb was created with Maestro (triolith, schrodinger/2014-1-nsc)  

ash.ffld11 was obtained from ash.pdb with ffld_server:  
```$SCHRODINGER/utilities/ffld_server -print_parameters -version 11 -ipdb ash.pdb > ash.ffld14```

ash.prm and ash.lib were obtained from ash.ffld11 with QTools:  
```q_ffld2q.py ash.ffld14 ash -s ash.pdb```


