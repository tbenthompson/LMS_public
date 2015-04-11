import lms_code.lib.rep2 as rep2
from lms_code.analysis.run_bem import bemify, boundary_conditions,\
    assemble, constrain, solve, evaluate_surface_disp
from lms_code.analysis.simplified_bem import create_surface_mesh, \
    set_params
from codim1.core import simple_line_mesh, combine_meshes, ray_mesh

def create_fault_mesh(d):
    top_fault_vert = [0, -1e9]
    top = d['intersection_pt']
    joint = [4.20012e5 + 1.6, -2.006e4 - 5]
    bottom = [3.09134e5 + 1.1, -2.3376e4 - 3]
    detach = simple_line_mesh(d['fault_elements'], bottom, joint)
    d['fault_mesh'] = detach

if __name__ == "__main__":
    d = dict()
    set_params(d)
    create_fault_mesh(d)
    create_surface_mesh(d)
    bemify(d)
    boundary_conditions(d)
    assemble(d)
    # constrain(d)
    solve(d)
    evaluate_surface_disp(d)
    rep2.save("bem_just_detach", d)
