# iPHASimulator Package Data

This folder stores templates and parameter files used by workflow helpers.
It stays as package data and is not flattened into root-level Python modules.

PHA molecule names are generated in code with `iphasimulator.naming`:

| Meaning | Examples |
|---|---|
| Monomer / residue code | `3HB`, `3HO`, `3HDD` |
| Polymer code | `P3HB`, `P3HO`, `P3HDD` |
| Single oligomer chain | `P3HB_4`, `P3HO_8`, `P3HDD_4` |
| Multi-chain system | `25_P3HB_3`, `10_P3HO_8` |
| Head/main/tail residue entries | `3HB_H`, `3HB_M`, `3HB_T` |

Files under `gromacs_mdp/` and `gromacs_solvation/` are simulation templates;
do not rename them to match polymer systems. `route_c_gromacs/` holds the Route C
run files (mdp and scripts) that notebook 06D adds to a CHARMM-GUI Solution Builder
package; see its README.
