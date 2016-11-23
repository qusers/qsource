#!/usr/bin/env python
load ../init.pdb
load_traj ../start.dcd, init
load_traj ../heat1.dcd, init
select waters, resn HOH
select protein, not resn HOH
hide everything
show cartoon, protein
spectrum count, rainbow, protein
deselect
