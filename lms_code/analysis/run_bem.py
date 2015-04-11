from codim1.core import *
from codim1.fast_lib import *
from codim1.assembly import *
from codim1.assembly.sgbem import sgbem_assemble
from codim1.post import *
from codim1.core.tools import plot_mesh
from math import atan, cos
import numpy as np
import lms_code.lib.rep2 as rep2
import sys
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

def transform_fault_mesh(intersection_pt, vertices, elements, refine,
                         criteria = lambda x: True):
    # Puts the fault mesh into the codim1 required format
    # also, puts the element list in the same format
    # Selects the distance from the NW as the x axis and the depth as Y

    new_v_list = []
    kmtom = 1000
    def is_copy(v):
        for v in new_v_list:
            if np.all(np.abs(v - new_v) < [0.01, 0.01]):
                return True
        return False

    for i, (k, v) in enumerate(vertices.iteritems()):
        in_meters = v[3] * kmtom
        if criteria(in_meters):
            new_v = np.array((in_meters, v[2]))
            if is_copy(new_v):
                continue
            new_v_list.append(new_v)

    new_v_list.sort(key = lambda x: -x[0])
    if criteria(intersection_pt):
        new_v_list.insert(0, intersection_pt)
    new_el_list = [(i, i + 1) for i in range(len(new_v_list) - 1)]

    return np.array(new_v_list), np.array(new_el_list)

def create_fault_mesh(d, geom, criteria = lambda x: True):
    vertices = geom['twod_vertices']
    elements = geom['twod_segments']
    vertices, elements = transform_fault_mesh(d['intersection_pt'],
                                              vertices,
                                              elements,
                                              d['fault_refine'],
                                              criteria = criteria)
    fault_mesh = from_vertices_and_etov(vertices, elements, flip = True)
    d['fault_mesh'] = fault_mesh
    # from matplotlib import pyplot as plt
    # plt.plot(vertices[:, 0], vertices[:, 1], 'o-')
    # plt.show()

def create_surface_mesh(d, geom):
    vx = geom['elevation_dist']
    # convert from km to meters
    vx = 1000 * np.array(vx)
    vy = geom['elevations']

    vertices = np.array(zip(vx, vy))[::d['skip_vertices']]
    fst_index_past = next(i for i in range(len(vertices))
                          if vertices[i][0] > d['intersection_pt'][0])
    vertices = np.insert(vertices, fst_index_past,
                         np.array(d['intersection_pt']), 0)
    elements = np.array([(i, i + 1) for i in range(len(vertices) - 1)])

    surf_mesh = from_vertices_and_etov(vertices, elements)

    # extend the mesh into the far field
    ray_left_dir = (-1.0, 0.0)
    mesh2 = ray_mesh(vertices[0, :], ray_left_dir,
                     d['far_ray_lengths'], flip = True)

    ray_right_dir = (1.0, 0.0)
    mesh3 = ray_mesh(vertices[-1, :], ray_right_dir, d['far_ray_lengths'])
    surf_mesh = combine_meshes(mesh2, combine_meshes(surf_mesh, mesh3),
                               ensure_continuity = True)


    d['surf_mesh'] = surf_mesh

def bemify(d):
    combined_mesh = combine_meshes(d['surf_mesh'],
                                   d['fault_mesh'],
                                   ensure_continuity = True)
    d['bf'] = gll_basis(d['degree'])
    min_pts = d['degree'] + 1

    d['qs'] = GLLQuadStrategy(combined_mesh, min_pts, d['quad_max'],
    # d['qs'] = QuadStrategy(combined_mesh, min_pts, d['quad_max'],
                                 d['quad_logr'], d['quad_oneoverr'])
    apply_to_elements(combined_mesh, "qs", d['qs'], non_gen = True)
    apply_to_elements(combined_mesh, "basis", d['bf'], non_gen = True)
    d['combined_mesh'] = combined_mesh

def get_slip_magnitude(e):
    left_end = e.vertex1.loc
    right_end = e.vertex2.loc
    if right_end[0] < left_end[0]:
        print "Backwards elements!"
        left_end, right_end = right_end, left_end
    fault_vector = left_end - right_end
    fault_tangential = fault_vector / np.linalg.norm(fault_vector)
    theta = atan(fault_tangential[1] / fault_tangential[0])
    # Simple slip magnitude rule
    slip = fault_tangential / cos(theta)
    return slip

def boundary_conditions(d):
    apply_to_elements(d['surf_mesh'], "bc",
                    BC("traction", ZeroBasis()), non_gen = True)

    # Slip in the fault tangential direction
    slips = []
    def make_slip_with_magnitude(e):
        slip = get_slip_magnitude(e)
        slips.append(slip)
        # Weird faster at depth rule.
        return BC("crack_displacement", ConstantBasis(slip))

    apply_to_elements(d['fault_mesh'], "bc", make_slip_with_magnitude)


def apply_jump_constraint(fault_mesh, matrix, rhs):
    sorted_elements = sorted(fault_mesh.elements, key = lambda e: e.vertex2.loc[1])
    slip = sorted_elements[-1].bc.basis.evaluate(0, 0.0)
    lse = sorted_elements[-1].vertex2.connected_to[0]
    rse = sorted_elements[-1].vertex2.connected_to[1]
    constraint_dofx = lse.dofs[0, -1]
    constraint_dofy = lse.dofs[1, -1]
    other_dofx = rse.dofs[0, 0]
    other_dofy = rse.dofs[1, 0]
    matrix[constraint_dofx, :] = 0
    matrix[constraint_dofy, :] = 0
    rhs[constraint_dofx] = slip[0]
    rhs[constraint_dofy] = slip[1]
    matrix[constraint_dofx, constraint_dofx] = -1
    matrix[constraint_dofx, other_dofx] = 1
    matrix[constraint_dofy, constraint_dofy] = -1
    matrix[constraint_dofy, other_dofy] = 1

def assemble(d):
    sgbem_dofs(d['combined_mesh'])
    ek = ElasticKernelSet(d['shear_modulus'], d['poisson_ratio'])
    d['matrix'], d['rhs'] = sgbem_assemble(d['combined_mesh'], ek)

def constrain(d, jump = True):
    pin_ends_constraint(d['matrix'], d['rhs'], d['surf_mesh'], [0, 0], [0, 0])
    if jump:
        apply_jump_constraint(d['fault_mesh'], d['matrix'], d['rhs'])

def solve(d):
    soln_coeffs = np.linalg.solve(d['matrix'], d['rhs'])
    d['soln_coeffs'] = soln_coeffs
    apply_coeffs(d['combined_mesh'], d['soln_coeffs'], "soln")

def evaluate_surface_disp(d):
    d['x'], d['u_soln'], t = \
        evaluate_boundary_solution(d['surf_mesh'], d['soln_coeffs'], d['surface_points'])

if __name__ == "__main__":
    geom = rep2.load('lms_geometry')

    # Compute the BEM solution
    d = dict()
    set_params(d)
    create_fault_mesh(d, geom)
    create_surface_mesh(d, geom)
    bemify(d)
    boundary_conditions(d)
    assemble(d)
    constrain(d)
    solve(d)
    evaluate_surface_disp(d)
    rep2.save("bem_all_details", d)
