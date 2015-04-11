import lms_code.lib.rep2 as rep2
from codim1.core import *
from codim1.fast_lib import *
from lms_code.analysis.run_bem import bemify,\
    assemble, constrain, solve, evaluate_surface_disp,get_slip_magnitude
from lms_code.analysis.simplified_bem import create_surface_mesh, \
    set_params
from codim1.core import simple_line_mesh, combine_meshes, ray_mesh

def boundary_conditions(d):
    apply_to_elements(d['surf_mesh'], "bc",
                    BC("traction", ZeroBasis()), non_gen = True)

    # Slip in the fault tangential direction
    def make_trac_with_magnitude(e):
        trac = get_slip_magnitude(e)
        print trac
        return BC("displacement", ConstantBasis(trac))

    apply_to_elements(d['fault_mesh'], "bc", make_trac_with_magnitude)

def create_fault_mesh(d):
    top_fault_vert = [0, -1e9]
    top = d['intersection_pt']
    joint = [4.20012e5 + 1.6, -2.006e4 - 5]
    bottom = [3.09134e5 + 1.1, -2.3376e4 - 3]
    beichuan = simple_line_mesh(d['fault_elements'], joint, top)
    detach = simple_line_mesh(d['fault_elements'], bottom, joint)
    # d['beichuan_mesh'] = beichuan
    d['detachment_mesh'] = detach
    d['fault_mesh'] = detach
    # d['fault_mesh'] = combine_meshes(beichuan, detach, ensure_continuity = True)

if __name__ == "__main__":
    d = dict()
    set_params(d)
    d['skip_vertices'] = 30
    d['fault_elements'] = 5
    d['far_steps'] = 3
    d['far_mult'] = 5000.0
    far_ray_lengths = [1.0] * d['far_per_step']
    for i in range(1, d['far_steps']):
        new_els = [d['far_mult'] * (2.0 ** float(i))] * d['far_per_step']
        far_ray_lengths.extend(new_els)
    d['far_ray_lengths'] = far_ray_lengths

    create_fault_mesh(d)
    create_surface_mesh(d)
    bemify(d)
    boundary_conditions(d)
    assemble(d)
    import ipdb; ipdb.set_trace()
    constrain(d, jump = False)
    solve(d)
    evaluate_surface_disp(d)
    rep2.save("bem_traction", d)
