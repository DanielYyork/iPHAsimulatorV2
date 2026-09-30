#!/usr/bin/env bash
set -euo pipefail

gmx grompp -f step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p topol.top -n index.ndx -o step6.0_minimization.tpr -maxwarn 1
gmx mdrun -deffnm step6.0_minimization
