#!/bin/bash
echo "[MD]
steps                     10000
stepsize                  1
temperature               300.0
bath_coupling             10
separate_scaling          on
lrf                       on

[cut_offs]
solute_solute             10
solute_solvent            10
solvent_solvent           10
q_atom                    99
lrf                       99

[sphere]
shell_force               10
shell_radius              0.8

[trajectory_atoms]
not excluded

[intervals]
non_bond                  20
output                    100
energy                    20
trajectory                5000

[files]
topology                  lig_p.top
fep                       lig_1u.fep
restart                   eq5.re
final                     md_1_00.re
energy                    md_1_00.en
trajectory                md_1_00.dcd

[lambdas]
1.00  0.00" > md_1_00.inp

for i in $(seq -w 1 1 50)
do
    j=$(awk 'BEGIN{printf "%02i", ('''$i'''-1)}')
    k=$(awk 'BEGIN{printf "%2.2f", 1.00-('''$i'''*0.02)}')
    l=$(awk 'BEGIN{printf "%2.2f", ('''$i'''*0.02)}')
echo "[MD]
steps                     10000
stepsize                  1
temperature               300.0
bath_coupling             10
separate_scaling          on
lrf                       on

[cut_offs]
solute_solute             10
solute_solvent            10
solvent_solvent           10
q_atom                    99
lrf                       99

[sphere]
shell_force               10
shell_radius              0.8

[trajectory_atoms]
not excluded

[intervals]
non_bond                  20
output                    100
energy                    20
trajectory                5000

[files]
topology                  lig_p.top
fep                       lig_1u.fep
restart                   md_1_$j.re
final                     md_1_$i.re
energy                    md_1_$i.en
trajectory                md_1_$i.dcd

[lambdas]
"$k"  "$l"
" > md_1_$i.inp
done
