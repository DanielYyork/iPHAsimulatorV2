#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.iphasimulator.pha_filepath_manager import PHAFileManager
from src.iphasimulator.sw_openmm import GromacsSimulation


if __name__ == "__main__":

    paths = PHAFileManager("../structure_database")

    polymer_names = ["P3HB_10"]
    number_of_polymers = [25]

    melt_name = paths.get_PHA_melt_name(
        polymer_names,
        number_of_polymers,
    )

    melt_dir = paths.get_PHA_melt_dir(
        polymer_names,
        number_of_polymers,
    )

    topology_file = melt_dir / f"{melt_name}.top"
    coordinate_file = melt_dir / f"{melt_name}.gro"
    
    run_dir = paths.create_named_PHA_melt_simulation_run_dir(
        polymer_names,
        number_of_polymers,
        run_name="Test",
    )

    print("Melt name:", melt_name)
    print("Melt dir:", melt_dir)
    print("Topology file:", topology_file)
    print("Coordinate file:", coordinate_file)
    print("OpenMM output directory:", run_dir)

    if not topology_file.exists():
        raise FileNotFoundError(topology_file)

    if not coordinate_file.exists():
        raise FileNotFoundError(coordinate_file)

    sim = GromacsSimulation(
        paths,
        str(topology_file),
        str(coordinate_file),
        output_dir=str(run_dir),
    )

    print(sim)

    minimized_sim = sim.minimize_energy()

    nvt_sim, data_file = sim.basic_NVT(
        minimized_sim,
        total_steps=3000,
        temp=300,
        filename="test_NVT",
    )

    print("Simulation finished.")
    print("Data file:", data_file)
    print("Output directory:", sim.output_dir)
    
    sim.graph_state_data(data_file)