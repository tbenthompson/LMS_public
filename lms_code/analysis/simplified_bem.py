import lms_code.lib.rep2 as rep2
from lms_code.analysis.run_bem import bemify, boundary_conditions,\
    assemble, constrain, solve, evaluate_surface_disp
from codim1.core import simple_line_mesh, combine_meshes, ray_mesh

def set_params(d):
    d['intersection_pt'] = [461386, 1590]
    d['degree'] = 3
    d['fault_elements'] = 30
    d['surface_points'] = 50

    d['quad_max'] = 20
    d['quad_logr'] = 20
    d['quad_oneoverr'] = 20

    d['shear_modulus'] = 30e9
    d['poisson_ratio'] = 0.25

    d['slip_magnitude'] = 1.0


    d['far_per_step'] = 8
    d['far_steps'] = 11
    # faster, for debugging
    # d['far_steps'] = 4
    # d['degree'] = 1
    # d['fault_elements'] = 5
    d['far_mult'] = 250.0
    far_ray_lengths = [d['far_mult']] * d['far_per_step']
    for i in range(1, d['far_steps']):
        new_els = [d['far_mult'] * (2.0 ** float(i))] * d['far_per_step']
        far_ray_lengths.extend(new_els)
    d['far_ray_lengths'] = far_ray_lengths


def create_fault_mesh(d):
    top_fault_vert = [0, -1e9]
    top = d['intersection_pt']
    joint = [4.20012e5 + 1.6, -2.006e4 - 5]
    bottom = [3.09134e5 + 1.1, -2.3376e4 - 3]
    beichuan = simple_line_mesh(d['fault_elements'], joint, top)
    detach = simple_line_mesh(d['fault_elements'], bottom, joint)
    fault_mesh = combine_meshes(beichuan, detach, ensure_continuity = True)
    d['fault_mesh'] = fault_mesh

def create_surface_mesh(d):
    top = d['intersection_pt']

    ray_left_dir = (-1.0, 0.0)
    left_side = ray_mesh(top, ray_left_dir, d['far_ray_lengths'], flip = True)

    ray_right_dir = (1.0, 0.0)
    right_side = ray_mesh(top, ray_right_dir, d['far_ray_lengths'])

    surf_mesh = combine_meshes(left_side, right_side, ensure_continuity = True)
    d['surf_mesh'] = surf_mesh

if __name__ == "__main__":
    d = dict()
    set_params(d)
    create_fault_mesh(d)
    create_surface_mesh(d)
    bemify(d)
    boundary_conditions(d)
    assemble(d)
    constrain(d)
    solve(d)
    evaluate_surface_disp(d)
    rep2.save("bem_simple", d)
