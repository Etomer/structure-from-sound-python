#!/bin/bash

for file in ./old_matlab_solvers/*; do
    ./convert_matlab_solver_to_python.sh "$file"
done

for file in ./generated_solvers/*; do
   python3 setup.py "$file" build_ext --build-lib solver_so
done