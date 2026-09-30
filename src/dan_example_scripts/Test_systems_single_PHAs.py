#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:40:49 2026

@author: daniel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import sys



from src.iphasimulator.build_single_PHA_systems import (
    build_dry_PHA,
    build_solvated_PHA,
    build_solvated_PHA_ions,
)


if __name__ == "__main__":

    root_dir = "structure_database"
    polymer_name = "P3HB_10"

    print("=" * 80)
    print("Testing dry PHA build")
    print("=" * 80)

    dry_result = build_dry_PHA(
        polymer_name=polymer_name,
        root_dir=root_dir,
        forcefield="gaff2",
        box_radius=20.0,
    )

    print(dry_result)

    print("=" * 80)
    print("Testing solvated PHA build")
    print("=" * 80)

    solvated_result = build_solvated_PHA(
        polymer_name=polymer_name,
        root_dir=root_dir,
        forcefield="gaff2",
        water_leaprc="water.tip3p",
        water_box="TIP3PBOX",
        box_radius=20.0,
    )

    print(solvated_result)

    print("=" * 80)
    print("Testing solvated + ions PHA build")
    print("=" * 80)

    ionised_result = build_solvated_PHA_ions(
        polymer_name=polymer_name,
        root_dir=root_dir,
        forcefield="gaff2",
        water_leaprc="water.tip3p",
        water_box="TIP3PBOX",
        box_radius=20.0,
        salt="KCl",
        pos_ion="K+",
        neg_ion="Cl-",
        ion_conc=0.15,
    )

    print(ionised_result)

    print("=" * 80)
    print("All tests finished.")
    print("=" * 80)