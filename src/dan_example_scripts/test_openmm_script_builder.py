#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 00:14:58 2026

@author: daniel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.iphasimulator.openmmscript_builder import OpenMMScriptBuilder


builder = OpenMMScriptBuilder(
    root_dir="../structure_database",
    polymer_names=["P3HB_10"],
    number_of_polymers=[25],
    run_name="BuilderTest",
)

builder.add_minimization()

builder.add_basic_NVT(
    total_steps=3000,
    temp=300,
    filename="test_NVT",
)

builder.add_basic_NPT(
    total_steps=3000,
    temp=300,
    pressure=1,
    filename="test_NPT",
)

builder.add_basic_NPT(
    total_steps=3000,
    temp=300,
    pressure=1,
    filename="test_NPT2",
)

output_script = builder.write_script(
    "dan_example_scripts/generated_builder_test.py"
)

print("Generated script:")
print(output_script)