# -*- coding: utf-8 -*-

import os
import time
import shutil
import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import parmed as pmd

import openmm
import openmm.app as app
from openmm import *
from openmm.app import *
from openmm.unit import *


class DcdWriter:
    """Write DCD trajectory files for OpenMM simulations."""

    def __init__(self, prefix, freq):
        self.dcdReporter = app.DCDReporter(
            f"{prefix}.dcd",
            freq,
            enforcePeriodicBox=True,
        )


class DataWriter:
    """Write OpenMM state data to a comma-separated text file."""

    def __init__(self, prefix, freq, steps):
        self.stateDataReporter = app.StateDataReporter(
            f"{prefix}.txt",
            freq,
            totalSteps=steps,
            step=True,
            time=True,
            speed=True,
            progress=True,
            elapsedTime=True,
            totalEnergy=True,
            kineticEnergy=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            separator=",",
        )


class BuildSimulation:

    """

    Parent class for OpenMM simulations.

    Child classes:

        AmberSimulation

        GromacsSimulation

    """

    savepdb_traj = False

    pressure = 1

    temp = 300

    min_temp = 0

    timestep = 2.0

    friction_coeff = 1.0

    total_steps = 1000

    reporter_freq = 1000

    nonbondedcutoff = 1.0

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

    anneal_parameters = [300, 700, 5, 10, 500000]

    minimized_only = None

    restrain_heavys = False

    def __init__(

        self,

        manager,

        filename,

        output_dir=None,

    ):

        self.manager = manager

        self.filename = filename

        if output_dir is None:

            self.output_dir = os.path.join(

                self.manager.systems_dir,

                self.filename,

                self.timestamp,

            )

        
        else:

            self.output_dir = output_dir
        self.run_name = os.path.basename(self.output_dir)       
        os.makedirs(

            self.output_dir,

            exist_ok=True,

        )

    def type_of_simulation(self):

        if isinstance(self, AmberSimulation):

            return "AMB"

        if isinstance(self, GromacsSimulation):

            return "GRO"

        return "UNKNOWN"

    def get_platform(self, platform_name=None):
        """
        Select and configure the OpenMM platform.

        Priority
        --------
        1. Explicit platform_name argument.
        2. IPHA_OPENMM_PLATFORM environment variable.
        3. CUDA.
        4. OpenCL.
        5. CPU.

        Precision
        ---------
        CUDA
            Mixed precision.

        OpenCL
            Single precision.

        CPU
            OpenMM default precision.

        Parameters
        ----------
        platform_name : str, optional
            Explicit OpenMM platform to use.

        Returns
        -------
        openmm.Platform
            Configured OpenMM platform.
        """

        import os
        from openmm import Platform

        requested_platform = (
            platform_name
            or os.environ.get(
                "IPHA_OPENMM_PLATFORM"
            )
        )

        # ---------------------------------------------------------
        # Helper for platform-specific configuration
        # ---------------------------------------------------------

        def configure_platform(
            platform,
            name,
        ):

            if name == "CUDA":

                platform.setPropertyDefaultValue(
                    "Precision",
                    "mixed",
                )

                print(
                    "CUDA precision: mixed"
                )

            elif name == "OpenCL":

                platform.setPropertyDefaultValue(
                    "Precision",
                    "single",
                )

                print(
                    "OpenCL precision: single"
                )

            return platform

        # ---------------------------------------------------------
        # Explicitly requested platform
        # ---------------------------------------------------------

        if requested_platform is not None:

            platform = Platform.getPlatformByName(
                requested_platform
            )

            platform = configure_platform(
                platform,
                requested_platform,
            )

            print(
                f"Using requested platform: "
                f"{requested_platform}"
            )

            return platform

        # ---------------------------------------------------------
        # Automatic platform selection
        # ---------------------------------------------------------

        for candidate in [
            "CUDA",
            "OpenCL",
            "CPU",
        ]:

            try:

                platform = (
                    Platform.getPlatformByName(
                        candidate
                    )
                )

                platform = configure_platform(
                    platform,
                    candidate,
                )

                print(
                    f"Using {candidate} platform."
                )

                return platform

            except Exception:
                continue

        # ---------------------------------------------------------
        # Nothing worked
        # ---------------------------------------------------------

        raise RuntimeError(
            "No usable OpenMM platform was found."
        )

    def create_openmm_system(
        self,
        ensemble="NVT",
        temp=None,
        pressure=None,
    ):
        if temp is None:
            temp = self.temp

        if pressure is None:
            pressure = self.pressure

        sim_type = self.type_of_simulation()

        if sim_type == "AMB":
            system = self.amb_topology.createSystem(
                nonbondedMethod=app.PME,
                nonbondedCutoff=self.nonbondedcutoff * nanometers,
                constraints=app.HBonds,
            )

        elif sim_type == "GRO":
            system = self.gro_topology.createSystem(
                nonbondedMethod=app.PME,
                nonbondedCutoff=self.nonbondedcutoff * nanometers,
                constraints=app.HBonds,
            )

        else:
            raise ValueError(
                "Simulation type must be AmberSimulation or GromacsSimulation."
            )

        if ensemble == "NPT":
            barostat = MonteCarloBarostat(
                pressure * atmosphere,
                temp * kelvin,
            )
            system.addForce(barostat)

        return system

    def create_openmm_simulation(
        self,
        system,
        integrator,
        platform,
    ):
        sim_type = self.type_of_simulation()

        if sim_type == "AMB":
            return app.Simulation(
                self.amb_topology.topology,
                system,
                integrator,
                platform,
            )

        if sim_type == "GRO":
            return app.Simulation(
                self.gro_topology.topology,
                system,
                integrator,
                platform,
            )

        raise ValueError(
            "Simulation type must be AmberSimulation or GromacsSimulation."
        )

    def minimize_energy(self):
        min_start_time = time.time()

        integrator = LangevinIntegrator(
            self.min_temp * kelvin,
            self.friction_coeff / picoseconds,
            self.timestep * femtoseconds,
        )

        platform = self.get_platform()

        system = self.create_openmm_system(
            ensemble="NVT",
            temp=self.min_temp,
        )

        simulation = self.create_openmm_simulation(
            system,
            integrator,
            platform,
        )

        if self.type_of_simulation() == "AMB":
            simulation.context.setPositions(
                self.amb_coordinates.positions
            )

        elif self.type_of_simulation() == "GRO":
            print("Atoms in GRO:", len(self.gro_coordinates.positions))
            print("Atoms in TOP:", self.gro_topology.topology.getNumAtoms())

            simulation.context.setPositions(
                self.gro_coordinates.positions
            )

        simulation.minimizeEnergy()

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        self.min_pdbname = os.path.join(
            self.output_dir,
            f"min_{self.filename}.pdb",
        )

        with open(self.min_pdbname, "w") as output:
            PDBFile.writeFile(
                simulation.topology,
                state.getPositions(),
                output,
            )

        time_taken = time.time() - min_start_time
        print(f"Minimization completed in {time_taken:.2f} seconds.")

        return simulation

    @classmethod
    def minimize_energy_help(cls):
        print(cls.minimize_energy.__doc__)

    def anneal_NVT(
        self,
        simulation,
        start_temp=None,
        max_temp=None,
        cycles=None,
        quench_rate=None,
        steps_per_cycle=None,
        filename=None,
        save_restart=False,
        restart_name=None,
        verbose=True,
    ):
        anneal_start_time = time.time()

        if filename is None:
            filename = "_anneal_"
        else:
            filename = f"_{filename}_"

        if start_temp is None:
            start_temp = self.anneal_parameters[0]

        if max_temp is None:
            max_temp = self.anneal_parameters[1]

        if cycles is None:
            cycles = self.anneal_parameters[2]

        if quench_rate is None:
            quench_rate = self.anneal_parameters[3]

        if steps_per_cycle is None:
            steps_per_cycle = self.anneal_parameters[4]

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            enforcePeriodicBox=True,
        )

        xyz = state.getPositions()
        vx, vy, vz = state.getPeriodicBoxVectors()

        integrator = LangevinIntegrator(
            start_temp * kelvin,
            self.friction_coeff / picoseconds,
            self.timestep * femtoseconds,
        )

        platform = self.get_platform()

        system = self.create_openmm_system(
            ensemble="NVT",
            temp=start_temp,
        )

        simulation = self.create_openmm_simulation(
            system,
            integrator,
            platform,
        )

        simulation.context.setPeriodicBoxVectors(vx, vy, vz)
        simulation.context.setPositions(xyz)

        total_steps = steps_per_cycle * cycles
    
        if self.savepdb_traj is True:
            output_pdbname = os.path.join(
                self.output_dir,
                self.filename + filename + self.run_name + ".pdb",
            )

            simulation.reporters.append(
                app.PDBReporter(
                    output_pdbname,
                    self.reporter_freq,
                )
            )

        output_dcdname = os.path.join(
            self.output_dir,
            self.filename + filename + self.run_name,
        )

        dcdWriter = DcdWriter(
            output_dcdname,
            self.reporter_freq,
        )

        simulation.reporters.append(
            dcdWriter.dcdReporter
        )

        output_dataname = os.path.join(
            self.output_dir,
            self.filename + filename + self.run_name,
        )

        dataWriter = DataWriter(
            output_dataname,
            self.reporter_freq,
            total_steps,
        )

        simulation.reporters.append(
            dataWriter.stateDataReporter
        )

        increments = int(
            (max_temp - start_temp) / quench_rate
        )

        steps_per_slope = int(
            steps_per_cycle * 0.4
        )

        holding_steps = int(
            steps_per_cycle * 0.1
        )

        steps_at_increment = int(
            steps_per_slope / increments
        )

        if verbose is True:
            print(
                f"""Annealing information:
                 - Number of heating/cooling increments: {increments}
                 - Steps per temperature in-/decrease: {steps_per_slope}
                 - Holding steps at {max_temp} K: {holding_steps}
                 - Steps at heating/cooling increment: {steps_at_increment}
                 - Total simulation time {(total_steps * self.timestep):.0f} fs
                """
            )

        def cycle():
            integrator.setTemperature(
                start_temp * kelvin
            )

            simulation.step(
                steps_at_increment
            )

            for i in range(increments):
                integrator.setTemperature(
                    (start_temp + i * quench_rate) * kelvin
                )

                simulation.step(
                    steps_at_increment
                )

            integrator.setTemperature(
                max_temp * kelvin
            )

            simulation.step(
                holding_steps
            )

            for i in range(increments):
                integrator.setTemperature(
                    (max_temp - i * quench_rate) * kelvin
                )

                simulation.step(
                    steps_at_increment
                )

            integrator.setTemperature(
                start_temp * kelvin
            )

            simulation.step(
                holding_steps
            )

        for _ in range(cycles):
            cycle()

        time_taken = time.time() - anneal_start_time

        print(
            f"Annealing completed in {time_taken:.2f} seconds."
        )

        final_state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        self.final_pdbname = os.path.join(
            self.output_dir,
            "final" + filename + self.filename + ".pdb",
        )

        with open(
            self.final_pdbname,
            "w",
        ) as output:
            PDBFile.writeFile(
                simulation.topology,
                final_state.getPositions(),
                output,
            )

        if save_restart is True:
            self.save_rst(
                simulation,
                restart_name=restart_name,
            )

        return simulation, output_dataname + ".txt"

    @classmethod
    def anneal_help(cls):
        print(
            cls.anneal_NVT.__doc__
        )

    def basic_NPT(
        self,
        simulation,
        total_steps=None,
        temp=None,
        pressure=None,
        filename=None,
        save_restart=False,
        restart_name=None,
        verbose=True,
    ):
        equili_start_time = time.time()

        if filename is None:
            filename = "_basic_NPT_"
        else:
            filename = f"_{filename}_"

        if total_steps is None:
            total_steps = self.total_steps

        if temp is None:
            temp = self.temp

        if pressure is None:
            pressure = self.pressure

        if verbose is True:
            print(
                f"""Basic NPT information:
                - Total steps: {total_steps}
                - Total simulation time: {(total_steps * self.timestep):.0f} fs
                - Temperature: {temp} K
                - Pressure: {pressure} atm
                """
            )

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        xyz = state.getPositions()
        vx, vy, vz = state.getPeriodicBoxVectors()

        integrator = LangevinIntegrator(
            temp * kelvin,
            self.friction_coeff / picoseconds,
            self.timestep * femtoseconds,
        )

        platform = self.get_platform()

        system = self.create_openmm_system(
            ensemble="NPT",
            temp=temp,
            pressure=pressure,
        )

        simulation = self.create_openmm_simulation(
            system,
            integrator,
            platform,
        )

        simulation.context.setPeriodicBoxVectors(
            vx,
            vy,
            vz,
        )

        simulation.context.setPositions(
            xyz
        )

        simulation.context.setVelocitiesToTemperature(
            temp * kelvin
        )

        if self.savepdb_traj is True:
            output_pdbname = os.path.join(
                self.output_dir,
                self.filename
                + "_"
                + str(pressure)
                + filename
                + self.run_name
                + ".pdb",
            )

            simulation.reporters.append(
                app.PDBReporter(
                    output_pdbname,
                    self.reporter_freq,
                )
            )

        output_dcdname = os.path.join(
            self.output_dir,
            self.filename
            + "_"
            + str(pressure)
            + filename
            + self.run_name,
        )

        dcdWriter = DcdWriter(
            output_dcdname,
            self.reporter_freq,
        )

        simulation.reporters.append(
            dcdWriter.dcdReporter
        )

        output_dataname = os.path.join(
            self.output_dir,
            self.filename
            + "_"
            + str(pressure)
            + filename
            + self.run_name,
        )

        dataWriter = DataWriter(
            output_dataname,
            self.reporter_freq,
            total_steps,
        )

        simulation.reporters.append(
            dataWriter.stateDataReporter
        )

        simulation.step(
            total_steps
        )

        time_taken = time.time() - equili_start_time

        print(
            f"Basic NPT completed in {time_taken:.2f} seconds."
        )

        final_state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        self.final_pdbname = os.path.join(
            self.output_dir,
            "final_"
            + filename
            + self.filename
            + "_"
            + str(pressure)
            + "_atm.pdb",
        )

        with open(
            self.final_pdbname,
            "w",
        ) as output:
            PDBFile.writeFile(
                simulation.topology,
                final_state.getPositions(),
                output,
            )

        if save_restart is True:
            self.save_rst(
                simulation,
                restart_name=restart_name,
            )

        return simulation, output_dataname + ".txt"

    @classmethod
    def basic_NPT_help(cls):
        print(
            cls.basic_NPT.__doc__
        )

    def basic_NVT(
        self,
        simulation,
        total_steps=None,
        temp=None,
        filename=None,
        save_restart=False,
        restart_name=None,
        verbose=True,
    ):
        prod_start_time = time.time()

        if filename is None:
            filename = "_basic_NVT_"
        else:
            filename = f"_{filename}_"

        if total_steps is None:
            total_steps = self.total_steps

        if temp is None:
            temp = self.temp

        if verbose is True:
            print(
                f"""Basic NVT information:
                - Total steps: {total_steps}
                - Total simulation time: {(total_steps * self.timestep):.0f} fs
                - Temperature: {temp} K
                """
            )

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        xyz = state.getPositions()
        vx, vy, vz = state.getPeriodicBoxVectors()

        integrator = LangevinIntegrator(
            temp * kelvin,
            self.friction_coeff / picoseconds,
            self.timestep * femtoseconds,
        )

        platform = self.get_platform()

        system = self.create_openmm_system(
            ensemble="NVT",
            temp=temp,
        )

        simulation = self.create_openmm_simulation(
            system,
            integrator,
            platform,
        )

        simulation.context.setPositions(
            xyz
        )

        simulation.context.setPeriodicBoxVectors(
            vx,
            vy,
            vz,
        )

        if self.savepdb_traj is True:
            output_pdbname = os.path.join(
                self.output_dir,
                self.filename
                + filename
                + self.run_name
                + ".pdb",
            )

            simulation.reporters.append(
                app.PDBReporter(
                    output_pdbname,
                    self.reporter_freq,
                )
            )

        output_dcdname = os.path.join(
            self.output_dir,
            self.filename
            + filename
            + self.run_name,
        )

        dcdWriter = DcdWriter(
            output_dcdname,
            self.reporter_freq,
        )

        simulation.reporters.append(
            dcdWriter.dcdReporter
        )

        output_dataname = os.path.join(
            self.output_dir,
            self.filename
            + filename
            + self.run_name,
        )

        dataWriter = DataWriter(
            output_dataname,
            self.reporter_freq,
            total_steps,
        )

        simulation.reporters.append(
            dataWriter.stateDataReporter
        )

        simulation.step(
            total_steps
        )

        time_taken = time.time() - prod_start_time

        print(
            f"Basic NVT completed in {time_taken:.2f} seconds."
        )

        final_state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        self.final_pdbname = os.path.join(
            self.output_dir,
            "final_" + filename + self.filename + ".pdb",
        )

        with open(
            self.final_pdbname,
            "w",
        ) as output:
            PDBFile.writeFile(
                simulation.topology,
                final_state.getPositions(),
                output,
            )

        if save_restart is True:
            self.save_rst(
                simulation,
                restart_name=restart_name,
            )

        return simulation, output_dataname + ".txt"
    
    def thermal_ramp(
        self,
        simulation,
        heating=None,
        quench_rate=None,
        ensemble=None,
        start_temp=None,
        max_temp=None,
        total_steps=None,
        pressure=None,
        filename=None,
        save_restart=False,
        restart_name=None,
    ):
        if ensemble not in ["NVT", "NPT"]:
            print("Please specify 'NVT' or 'NPT' ensemble for the thermal ramp")
            return None

        if heating is None:
            print("Please specify True for heating or False for cooling")
            return None

        if quench_rate is None:
            print("Please specify a quench rate as an integer")
            return None

        thermal_ramp_start_time = time.time()

        if filename is None:
            filename = "_thermal_ramp_"
        else:
            filename = f"_{filename}_"

        if start_temp is None:
            start_temp = self.anneal_parameters[0]

        if max_temp is None:
            max_temp = self.anneal_parameters[1]

        if pressure is None:
            pressure = self.pressure

        if total_steps is None:
            total_steps = self.total_steps

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            enforcePeriodicBox=True,
        )

        xyz = state.getPositions()
        vx, vy, vz = state.getPeriodicBoxVectors()

        initial_temp = start_temp

        integrator = LangevinIntegrator(
            initial_temp * kelvin,
            self.friction_coeff / picoseconds,
            self.timestep * femtoseconds,
        )

        platform = self.get_platform()

        system = self.create_openmm_system(
            ensemble=ensemble,
            temp=initial_temp,
            pressure=pressure,
        )

        simulation = self.create_openmm_simulation(
            system,
            integrator,
            platform,
        )

        simulation.context.setPeriodicBoxVectors(
            vx,
            vy,
            vz,
        )

        simulation.context.setPositions(
            xyz
        )

        method = "heat" if heating is True else "cool"

        output_filename = (
            self.filename
            + filename
            + method
            + self.run_name
        )

        if self.savepdb_traj is True:
            output_pdbname = os.path.join(
                self.output_dir,
                output_filename + ".pdb",
            )

            simulation.reporters.append(
                app.PDBReporter(
                    output_pdbname,
                    self.reporter_freq,
                )
            )

        output_dcdname = os.path.join(
            self.output_dir,
            output_filename,
        )

        dcdWriter = DcdWriter(
            output_dcdname,
            self.reporter_freq,
        )

        simulation.reporters.append(
            dcdWriter.dcdReporter
        )

        output_dataname = os.path.join(
            self.output_dir,
            output_filename,
        )

        dataWriter = DataWriter(
            output_dataname,
            self.reporter_freq,
            total_steps,
        )

        simulation.reporters.append(
            dataWriter.stateDataReporter
        )

        if heating is True:
            if max_temp <= start_temp:
                raise ValueError(
                    f"Heating selected, but max_temp ({max_temp}) "
                    f"<= start_temp ({start_temp})"
                )

            incremental_temps = np.arange(
                start_temp,
                max_temp + abs(quench_rate),
                abs(quench_rate),
            ).tolist()

        else:
            if start_temp <= max_temp:
                raise ValueError(
                    f"Cooling selected, but start_temp ({start_temp}) "
                    f"<= max_temp ({max_temp})"
                )

            incremental_temps = np.arange(
                start_temp,
                max_temp - abs(quench_rate),
                -abs(quench_rate),
            ).tolist()

        if len(incremental_temps) == 0:
            raise ValueError(
                f"No temperature increments generated. "
                f"Check start_temp={start_temp}, "
                f"max_temp={max_temp}, "
                f"quench_rate={quench_rate}."
            )

        steps_at_increment = int(
            total_steps / len(incremental_temps)
        )

        for temp_i in incremental_temps:
            integrator.setTemperature(
                temp_i * kelvin
            )

            if ensemble == "NPT":
                simulation.context.setParameter(
                    MonteCarloBarostat.Temperature(),
                    float(temp_i),
                )

            simulation.step(
                steps_at_increment
            )

        time_taken = time.time() - thermal_ramp_start_time

        print(
            f"Thermal ramp completed in {time_taken:.2f} seconds."
        )

        final_state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        self.final_pdbname = os.path.join(
            self.output_dir,
            "final_" + output_filename + ".pdb",
        )

        with open(
            self.final_pdbname,
            "w",
        ) as output:
            PDBFile.writeFile(
                simulation.topology,
                final_state.getPositions(),
                output,
            )

        if save_restart is True:
            self.save_rst(
                simulation,
                restart_name=restart_name,
            )

        return simulation, output_dataname + ".txt"

    def save_rst(
        self,
        simulation,
        restart_name=None,
        overwrite=False,
    ):
        if restart_name is None:
            restart_name = f"final_{self.filename}"
        else:
            restart_name = f"{self.filename}_{restart_name.strip()}"

        parm = pmd.load_file(
            self.topology_file,
            self.coordinates_file,
        )

        final_state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
        )

        parm.coordinates = final_state.getPositions(
            asNumpy=True
        ).value_in_unit(angstroms)

        parm.velocities = final_state.getVelocities(
            asNumpy=True
        ).value_in_unit(angstroms / picoseconds)

        restart_folder = os.path.join(
            self.manager.systems_dir,
            restart_name,
        )

        os.makedirs(
            restart_folder,
            exist_ok=True,
        )

        rst_filename = os.path.join(
            restart_folder,
            f"{restart_name}.rst7",
        )

        parm.save(
            rst_filename,
            format="rst7",
            overwrite=overwrite,
        )

        new_top_name = os.path.join(
            restart_folder,
            f"{restart_name}.prmtop",
        )

        shutil.copy(
            self.topology_file,
            new_top_name,
        )

        return rst_filename

    def restrain_heavy_atoms(
        self,
        system,
        topology,
        positions,
    ):
        force = openmm.CustomExternalForce(
            "1000*(x-x0)^2 + 1000*(y-y0)^2 + 1000*(z-z0)^2"
        )

        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        for atom in topology.atoms():
            if atom.element.symbol != "H":
                pos = positions[atom.index]

                force.addParticle(
                    atom.index,
                    [
                        pos.x,
                        pos.y,
                        pos.z,
                    ],
                )

        system.addForce(force)

        return system

    def __repr__(self):
        print(
            "Simulation parameters: "
            "('{}', '{}', '{}, {}, {}, {}')".format(
                self.pressure,
                self.temp,
                self.timestep,
                self.friction_coeff,
                self.total_steps,
                self.reporter_freq,
            )
        )

        return (
            "Simulation parameters given in the following format: "
            "('{}', '{}', '{}, {}, {}, {}')".format(
                "pressure",
                "temperature",
                "timestep",
                "friction coefficient",
                "total steps",
                "reporter frequency",
            )
        )

    def __str__(self):
        return f"Simulation object of - {self.filename}"

    @classmethod
    def display_start_time(cls):
        print(
            "Simulation initiated at: ",
            cls.timestamp,
        )

    @classmethod
    def savepdb_trajectories(cls, boolean):
        if not isinstance(boolean, bool):
            print("Please pass True or False to this function.")
            return

        cls.savepdb_traj = boolean

        if boolean:
            print(
                "Simulation will save trajectories in both "
                ".pdb and .dcd format."
            )
        else:
            print(
                "Simulation will save trajectories in .dcd format only."
            )

    @classmethod
    def set_temperature(cls, temp):
        cls.temp = temp

        print(
            "Temperature set to: ",
            str(temp),
            "kelvin",
        )

    @classmethod
    def set_pressure(cls, pressure):
        cls.pressure = pressure

        print(
            "Pressure set to: ",
            str(pressure),
            " atmospheres",
        )

    @classmethod
    def set_timestep(cls, timestep):
        cls.timestep = timestep

        print(
            "Timestep set to: ",
            str(timestep),
        )

    @classmethod
    def set_friction_coeff(cls, friction_coeff):
        cls.friction_coeff = friction_coeff

        print(
            "Friction coefficient set to: ",
            str(friction_coeff),
        )

    @classmethod
    def set_total_steps(cls, total_steps):
        cls.total_steps = total_steps

        print(
            "Total steps for simulation set to: ",
            str(total_steps),
        )

    @classmethod
    def set_reporter_freq(cls, reporter_freq):
        cls.reporter_freq = reporter_freq

        print(
            "Reporter frequency set to every: ",
            reporter_freq,
            " steps",
        )

    @classmethod
    def set_nonbondedcutoff(cls, nonbondedcutoff):
        if not isinstance(nonbondedcutoff, float):
            print(
                "Please pass a float to this method to set "
                "a new nonbondedcutoff."
            )
            print(
                "Example: simulation.set_nonbondedcutoff(5.0)"
            )
            return None

        cls.nonbondedcutoff = nonbondedcutoff

    @classmethod
    def set_anneal_parameters(cls, new_anneal_parameters):
        if len(new_anneal_parameters) != len(cls.anneal_parameters):
            format_str = (
                "Expected format: "
                "[start_temp, max_temp, cycles, quench_rate, steps_per_cycle]"
            )

            raise ValueError(
                f"Invalid parameters provided. {format_str}"
            )

        cls.anneal_parameters = new_anneal_parameters

        print("Anneal parameters set.")
        print("Starting temperature is: ", str(new_anneal_parameters[0]))
        print("Target temperature is: ", str(new_anneal_parameters[1]))
        print("Number of annealing cycles is: ", str(new_anneal_parameters[2]))
        print("The quench rate is: ", str(new_anneal_parameters[3]))
        print("The number of steps per cycle is: ", str(new_anneal_parameters[4]))

    @classmethod
    def set_anneal_parameters_help(cls):
        print(cls.set_anneal_parameters.__doc__)

    @staticmethod
    def graph_state_data(data_file):
        png_file_name = data_file.rsplit(
            ".",
            1,
        )[0] + ".png"

        df = pd.read_csv(
            data_file,
            delimiter=",",
        )

        columns_to_plot = df.columns[3:-2]

        num_rows = (
            len(columns_to_plot) + 1
        ) // 2

        num_cols = min(
            2,
            len(columns_to_plot),
        )

        fig, axes = plt.subplots(
            num_rows,
            num_cols,
            figsize=(12, 4 * num_rows),
        )

        if num_rows == 1 and num_cols == 1:
            axes = np.array([[axes]])

        elif num_rows == 1:
            axes = np.array([axes])

        elif num_cols == 1:
            axes = np.array([[ax] for ax in axes])

        for i, column in enumerate(columns_to_plot):
            row = i // num_cols
            col = i % num_cols

            axes[row, col].plot(
                df["Time (ps)"],
                df[column],
            )

            axes[row, col].set_title(column)
            axes[row, col].set_xlabel("Time (ps)")
            axes[row, col].set_ylabel(column)
            axes[row, col].grid(True)

        for i in range(
            len(columns_to_plot),
            num_rows * num_cols,
        ):
            row = i // num_cols
            col = i % num_cols

            fig.delaxes(
                axes[row, col]
            )

        plt.tight_layout(
            rect=[0, 0, 1, 0.96]
        )

        plt.savefig(
            png_file_name,
            dpi=600,
        )

        plt.show()

    @classmethod
    def graph_state_data_help(cls):
        print(cls.graph_state_data.__doc__)


class GromacsSimulation(BuildSimulation):

    """

    OpenMM simulation using GROMACS topology and coordinate files.

    """

    def __new__(cls, *args, **kwargs):

        if len(args) < 3:

            raise TypeError(

                "Usage:\n"

                "sim = GromacsSimulation(\n"

                "    manager,\n"

                "    topology_file,\n"

                "    coordinates_file,\n"

                "    output_dir=None,\n"

                ")"

            )

        return super().__new__(cls)

    def __init__(

        self,

        manager,

        topology_file,

        coordinates_file,

        output_dir=None,

    ):

        self.manager = manager

        self.filename = os.path.basename(

            topology_file

        ).split(".")[0]

        super().__init__(

            manager,

            self.filename,

            output_dir=output_dir,

        )

        self.coordinates_file = coordinates_file

        self.topology_file = topology_file

        self.gro_coordinates = GromacsGroFile(

            coordinates_file

        )

        self.gro_topology = GromacsTopFile(

            topology_file,

            periodicBoxVectors=self.gro_coordinates.getPeriodicBoxVectors(),

        )

    def __str__(self):

        return f"GROMACS simulation object of - {self.filename}"


class AmberSimulation(BuildSimulation):
    """
    OpenMM simulation using AMBER topology and coordinate files.
    """

    def __new__(cls, *args, **kwargs):
        if len(args) < 3:
            raise TypeError(
                "Usage:\n"
                "sim = AmberSimulation(\n"
                "    manager,\n"
                "    topology_file,\n"
                "    coordinates_file,\n"
                "    output_dir=None,\n"
                ")"
            )

        return super().__new__(cls)

    def __init__(
        self,
        manager,
        topology_file,
        coordinates_file,
        output_dir=None,
    ):
        self.manager = manager

        self.filename = os.path.basename(
            topology_file
        ).split(".")[0]

        super().__init__(
            manager,
            self.filename,
            output_dir=output_dir,
        )

        self.coordinates_file = coordinates_file
        self.topology_file = topology_file

        self.amb_coordinates = app.AmberInpcrdFile(
            coordinates_file
        )

        self.amb_topology = app.AmberPrmtopFile(
            topology_file,
            periodicBoxVectors=self.amb_coordinates.boxVectors,
        )

    def __str__(self):
        return f"Amber simulation object of - {self.filename}"
