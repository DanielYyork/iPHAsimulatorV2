#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:36:58 2026

@author: daniel

This is an initial python script guide to be exported into a jupyter notebook.

This notebook serves as an end user guide on how to build polymer melts for MD simulation.
Pictures should be made available in this github also.

There are a few sections

0. Make the code from src/iphasimualtor available in this script
1. Intialise the melt builder
2. Define inputs
3. 3. Run melt builder
4. Run a test simulation

"""

"""
0. Make code inside src folder accesible
"""
from pathlib import Path
import sys

project_root = Path("/Users/daniel/projects/iPHAsimulatorV2")
src_dir = project_root / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

"""
1. Intialise the melt builder
"""

# Import functions
from src.iphasimulator.pha_melt_builder import PHAMeltBuilder

# Pass the structure database to the melt builder and intialise
melt_builder = PHAMeltBuilder(
    root_dir="../structure_database"
)

"""
2. Define inputs

    - polymer names: list; of polymers to be packed
    - number of polymers: list; of the number of each polymers to be packed 
        - NOTE: These 2 lists must be the same length
        
    - density: int; the desired density of the polymer melt - in mg/cm^3
"""
# Inputs
polymer_names = ["P3HB_10"]
number_of_polymers = [25]
density = 750

"""
3. Run melt builder

    Nothing to add here, just using the inputs already defined
"""

result = melt_builder.generate_polymer_melt(
    polymer_names=polymer_names,
    number_of_polymers=number_of_polymers,
    density=density,
)
    
"""
4. Run a test simulation

It is useful to run a short test simulation of these systems. This is to ensure the packing was carried out properly
and that there will be no undesired errors when running openmm simulations in a HPC environment.

If this 'Packing_check' simulation fails, just rerun this script.

Inputs for test simulation:
    - melt name: this is contained in the 'result' dictionary variable of the packing
    - topology file: this is contained in the 'result' dictionary variable of the packing
    - coordinate file: this is contained in the 'result' dictionary variable of the packing
    - run name: the desired name to give the test simulatin, by default this is "Packing_check"
    - test steps: number of desired steps in the packing check
    - temperature: desired temperature of this packing check
    - timestep: desired timestep of this packing checl
"""

test_result = melt_builder.test_polymer_melt_simulation(
    melt_name = result["melt_name"],
    topology_file = result["topology_file"],
    coordinate_file = result["coordinate_file"],
    run_name="Packing_check",
    test_steps=1000,
    temperature=300,
    timestep=1.0,
)