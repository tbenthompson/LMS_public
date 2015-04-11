#!/bin/bash

python lms_code/analysis/mesh_interior.py
mpirun -n 12 python lms_code/analysis/eval_interior.py
python lms_code/analysis/coalesce_interior.py
python lms_code/plots/plot_interior.py
