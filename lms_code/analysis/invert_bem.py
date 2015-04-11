from lms_code.analysis.run_bem import bemify, boundary_conditions,\
    assemble, constrain, solve, evaluate_surface_disp, set_params,\
    create_fault_mesh, create_surface_mesh, pin_ends_constraint,\
    apply_jump_constraint
import lms_code.lib.rep2 as rep2
from codim1.core.tools import plot_mesh
import sys
import matplotlib.pyplot as plt
sys.setrecursionlimit(50000)


def set_params(d):
    d['intersection_pt'] = (461386, 1590)
    d['degree'] = 3
    d['skip_vertices'] = 1
    d['fault_refine'] = 1
    d['quad_max'] = 20
    d['quad_logr'] = 20
    d['quad_oneoverr'] = 20
    d['surface_points'] = 10
    d['shear_modulus'] = 30e9
    d['poisson_ratio'] = 0.25
    d['slip_magnitude'] = 1.0
    d['far_per_step'] = 5
    d['far_steps'] = 14
    d['far_mult'] = 5000.0
    far_ray_lengths = [1.0] * d['far_per_step']
    for i in range(1, d['far_steps']):
        new_els = [d['far_mult'] * (2.0 ** float(i))] * d['far_per_step']
        far_ray_lengths.extend(new_els)
    d['far_ray_lengths'] = far_ray_lengths


if __name__ == "__main__":

    # Compute the BEM solution
    joint_x = 4.20012e5 + 1.6
    beichuan = lambda x: x > joint_x + 10
    detachment = lambda x: x < joint_x - 10
    faults = dict()
    faults['beichuan'] = beichuan
    faults['detachment'] = detachment


    solns = dict()
    for name in faults.keys():
        print name
        geom = rep2.load('lms_geometry')
        d = dict()
        set_params(d)
        create_fault_mesh(d, geom, faults[name])
        create_surface_mesh(d, geom)
        bemify(d)
        # plot_mesh(d['combined_mesh'])
        # plt.show()
        boundary_conditions(d)
        assemble(d)
        pin_ends_constraint(d['matrix'], d['rhs'], d['surf_mesh'], [0, 0], [0, 0])
        if name == 'beichuan':
            apply_jump_constraint(d['fault_mesh'], d['matrix'], d['rhs'])
        solve(d)
        evaluate_surface_disp(d)
        solns[name] = d
        # plt.plot(d['x'][0, :], d['u_soln'][0, :])
        # plt.show()

    rep2.save("bem_all_details_inverse1", solns)
