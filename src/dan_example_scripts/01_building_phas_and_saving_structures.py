#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:36:58 2026

@author: daniel

This is an initial python script guide to be exported into a jupyter notebook.
Much of this serves as a general guide for end users - however, where specified they are some
things included that have been deemed to be important for developers understanding of how the code
works.

There are a few sections

0. Make the code from src/iphasimualtor available in this script
1. Visualising and plotting available PHAs  
2. Building PHA polymers
3. Outputs of building PHA polymers (developer focussed)
4. Visualising PHA polymers
5. Validating sterochemistry of PHAs

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
1. Visualising available PHAs 

    - Saving an image of all avaialble PHA monomers
    - Plotting all available PHA monomers
    - Plotting a single PHA monomer
"""
# Import functions
from iphasimulator.visualisation.visualiser import *

# Define residue code csv - hardcoded for ease (this will need chainging as this is an absolute path)
residue_codes_csv = "/Users/daniel/projects/iPHAsimulatorV2/structure_database/residue_codes.csv"

# Save an image of all available PHA monomers
plot_available_PHA_monomers(
    residue_codes_csv=residue_codes_csv,
    output_file="available_PHA_monomers.png",
    mols_per_row=4)

# Output a plot of all available PHA monomers
show_available_PHA_monomers(
    residue_codes_csv=residue_codes_csv)

# Show the structure of a single PHA monomer - can be found from the printed out image
show_PHA_monomer(
    "3HHx",
    residue_codes_csv=residue_codes_csv)

"""
2. Building PHA polymers

    This uses preprameterized units - it does not cover how to parameterise new PHA units
"""
# Import functions
from src.iphasimulator.build_pha import *

# Intialise builder - this a relative path to the structure database (change as you need - although this should be fine for notebooks)
builder = PHAPolymerBuilder("../structure_database")

# Build a PHA - define type and length
output = builder.build_PHA_polymer("3HB", 20)

"""
3. Interrogating output of building a PHA

    This is more on the developer side - users will just need to know the name of the constructed polymer.
    
    This example below shows how the following filepaths can be retrieved:
        - amber parameters (.prmtop)
        - amber coordinates (.rst7)
        - pdb (.pdb)
        
    As well as these useful paths, specific directory paths can also be returned:
        - amber directory: where amber parameters are stored
        - leap directory: where the tleap script for building the polymer is stored
        - gromacs directory: where groamcs parameters are stored (this may be empty at this point if they are not yet required)
        
    The important thing is that these files can be retrieved using only the inputs that were used to build the polymer:
        - PHA type
        - length of polymer
        
    NOTE: the user does need to know how to do this, but this is a good example of how different files
        are retrieved by different builder functions (i.e. build_polymer_melt)
"""
# Import functions
from src.iphasimulator.pha_filepath_manager import PHAFileManager

# Initialise a filepath manager object - this can retrieve paths 
paths = PHAFileManager("../structure_database")

# Redefine PHA type and built polymer length and retrieve polymer name
PHA_type = "3HB"
length = 20
polymer_name = paths.get_built_PHA_name(PHA_type, length)

# Retreive paths for polymer specific directories
amber_dir = paths.get_built_PHA_amber_dir(PHA_type, length)
leap_dir = paths.get_built_PHA_leap_dir(PHA_type, length)
gromacs_dir = paths.get_built_PHA_gromacs_dir(PHA_type, length)

# Show paths for the files that currently exist
prmtop = amber_dir / f"{polymer_name}.prmtop"
rst7 = amber_dir / f"{polymer_name}.rst7"
pdb = amber_dir / f"{polymer_name}.pdb"

# Print out above defined information
print(f"""
      
      Interrogating file locations for polymer with TYPE: {PHA_type} and LENGTH: {str(length)}
      
      The given polymer name is: {polymer_name}
      
      Files for AMBER are saved at: {amber_dir}
      Files used in construction in the TLEAP programme are saved at: {leap_dir}
      A directory has been created for future GROMACS parameter files at: {gromacs_dir}
      
      A pdb file of this polymer is available at: {pdb}
      
      Amber files have been saved in the amber directory, the specific paths are:
          coordinates: {rst7}
          parameters: {prmtop}""")

"""
4. Validating and visualising PHA polymers
"""
# Define polymer smiles code csv - hardcoded for ease (this will need chainging as this is an absolute path)
polymer_smiles_csv = "/Users/daniel/projects/iPHAsimulatorV2/structure_database/polymer_smiles.csv"

# Output a plot of all available PHA polymers - these are built from SMILES so some may look strange
show_available_PHA_polymers(
    polymer_smiles_csv=polymer_smiles_csv)

# Output a plot a single PHA polymer
show_PHA_polymer(
    polymer_name="P3HB_20",
    polymer_smiles_csv=polymer_smiles_csv)
