from lms_code.analysis.run_bem import bemify, boundary_conditions,\
    assemble, constrain, solve, evaluate_surface_disp, set_params,\
    create_fault_mesh, create_surface_mesh, pin_ends_constraint,\
    apply_jump_constraint
import lms_code.lib.rep2 as rep2
from codim1.core.tools import plot_mesh
import sys
import matplotlib.pyplot as plt
from codim1.assembly import sgbem_dofs
from codim1.assembly.sgbem import sgbem_assemble
from codim1.core import *
from codim1.post import evaluate_solution_on_element
import numpy as np
import copy
import scipy.optimize
sys.setrecursionlimit(50000)

def set_params(d):
    d['intersection_pt'] = (461386, 1590)
    d['degree'] = 1
    d['skip_vertices'] = 5
    d['fault_refine'] = 1
    d['quad_max'] = 20
    d['quad_logr'] = 20
    d['quad_oneoverr'] = 20
    d['surface_points'] = 10
    d['shear_modulus'] = 30e9
    d['poisson_ratio'] = 0.25
    d['slip_magnitude'] = 1.0
    d['far_per_step'] = 5
    d['far_steps'] = 5
    d['far_mult'] = 5000.0
    far_ray_lengths = [1.0] * d['far_per_step']
    for i in range(1, d['far_steps']):
        new_els = [d['far_mult'] * (2.0 ** float(i))] * d['far_per_step']
        far_ray_lengths.extend(new_els)
    d['far_ray_lengths'] = far_ray_lengths

def get_soln_val(qi_x, soln_coeffs):
    for e in d['surf_mesh']:
        x0 = e.vertex1.loc[0]
        x1 = e.vertex2.loc[0]
        if x0 < qi_x < x1 or x1 < qi_x < x0:
            x_hat = (qi_x - x0) / (x1 - x0)
            u, t = evaluate_solution_on_element(e, x_hat, soln_coeffs)
            return u[0]
    raise Exception("Desired location is beyond surface mesh. Something\
                     is very seriously wrong. Flee!")

def build_G(d, qi_gps):
    sgbem_dofs(d['combined_mesh'])
    fault_dof_map = dict()
    next_fault_dof = 0
    for e in d['fault_mesh']:
        fault_dof_map[e] = next_fault_dof
        next_fault_dof += 1
    ek = ElasticKernelSet(d['shear_modulus'], d['poisson_ratio'])
    full_matrix, full_rhs = sgbem_assemble(d['combined_mesh'], ek)

    G = np.zeros((len(qi_gps['dist']), next_fault_dof))
    for e in d['fault_mesh']:
        print e.id
        matrix = full_matrix
        vertices = [e.vertex1, e.vertex2]
        elements = [e]
        one_element = Mesh(vertices, elements)
        combined_mesh = combine_meshes(d['surf_mesh'],
                                       one_element,
                                       ensure_continuity = True)
        combined_mesh.total_dofs = d['combined_mesh'].total_dofs
        local_matrix, local_rhs = sgbem_assemble(combined_mesh, ek, just_rhs = True)
        pin_ends_constraint(matrix, local_rhs, d['surf_mesh'], [0, 0], [0, 0])
        # check for surface intersecting element
        if abs(e.vertex2.loc[0] - 461386) < 0.01:
            matrix = copy.copy(matrix)
            apply_jump_constraint(one_element, matrix, local_rhs)
        soln_coeffs = np.linalg.solve(matrix, local_rhs)
        for i, qi_x in enumerate(qi_gps['dist']):
            value = get_soln_val(qi_x, soln_coeffs)
            G[i, fault_dof_map[e]] = value
    inv_d = dict()
    inv_d['G'] = G
    inv_d['bem_data'] = d
    inv_d['fault_dof_map'] = fault_dof_map
    inv_d['bem_matrix'] = full_matrix
    inv_d['bem_rhs'] = full_rhs

    rep2.save("bem_all_details_inverse_all_segs", inv_d)

def regularized_invert(inv_d, qi_gps, rp, type = 'tikhonov'):
    G = inv_d['G']
    f_dof_map = inv_d['fault_dof_map']
    reg_term = np.zeros((G.shape[1], G.shape[1]))
    # lambda * weight of element
    for e in inv_d['bem_data']['fault_mesh']:
        dof = f_dof_map[e]
        if type == 'tikhonov':
            reg_term[dof, dof] = rp * e.length
        elif type == '1st_deriv':
            try:
                left_e = e.neighbors_left[0]
                right_e = e.neighbors_right[0]
                if left_e.bc.type != "crack_displacement" or \
                    right_e.bc.type != "crack_displacement":
                    continue
                left_dof = f_dof_map[left_e]
                right_dof = f_dof_map[right_e]
                # 1st deriv
                reg_term[dof, dof] = 1 * rp / (e.length ** 1)
                reg_term[dof, right_dof] = -1 * rp / (e.length ** 1)
            except IndexError:
                continue
    # regularized least squares matrix
    reg_G = np.vstack((G, reg_term))
    print("Condition number: " + str(np.linalg.cond(reg_G)))
    reg_d = qi_gps['vel'] + [0] * G.shape[1]
    # soln = np.linalg.lstsq(reg_G, reg_d)
    soln = scipy.optimize.nnls(reg_G, reg_d)
    Lx2 = 0.0
    if type == 'tikhonov':
        Lx2 = np.sum(soln[0] ** 2)
    elif type == '1st_deriv':
        for e in inv_d['bem_data']['fault_mesh']:
            dof = f_dof_map[e]
            try:
                left_e = e.neighbors_left[0]
                right_e = e.neighbors_right[0]
                if left_e.bc.type != "crack_displacement" or \
                    right_e.bc.type != "crack_displacement":
                    continue
                left_dof = f_dof_map[left_e]
                right_dof = f_dof_map[right_e]
                # 1st deriv
                Lx2 += ((soln[0][dof] - soln[0][right_dof]) / e.length) ** 2
            except IndexError:
                continue
    return soln, Lx2

def L_plot(inv_d, qi_gps):
    G = inv_d['G']
    # lambda

    # reg_params = 5.0 ** np.linspace(0, 9, 21)
    reg_params = [0.0]#5.0 ** np.linspace(-7, -5, 21)
    x2 = []
    residual2 = []
    type = 'tikhonov'
    for rp in reg_params:
        soln, Lx2_norm = regularized_invert(inv_d, qi_gps, rp, type = type)
        x2.append(Lx2_norm)
        residual2.append(np.sum((G.dot(soln[0]) - qi_gps['vel']) ** 2))

    plt.loglog(x2, residual2, 'o')
    for i,rp in enumerate(reg_params):
        x = x2[i]
        y = residual2[i]
        plt.annotate(i, (x, y))
    plt.show()
    idx = int(raw_input("Type which regularization level to show: "))
    soln, Lx2 = regularized_invert(inv_d, qi_gps, reg_params[idx], type = type)
    print("Lx2: " + str(Lx2))

    plt.plot(qi_gps['dist'], qi_gps['vel'], '.')
    plt.plot(qi_gps['dist'], G.dot(soln[0]), 'o')
    plt.show()

    x = []
    y = []
    for e in inv_d['bem_data']['fault_mesh']:
        x.append(e.vertex2.loc[0])
        y.append(soln[0][inv_d['fault_dof_map'][e]])
        x.append(e.vertex1.loc[0])
        y.append(soln[0][inv_d['fault_dof_map'][e]])
    plt.plot(x, y, 'o-')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Provide an operation name [build, solve] as the first argument!")
        sys.exit()
    width = 0.5
    qi_gps = rep2.load('qi_gps_near_' + str(width))
    if sys.argv[1] == 'build':
        geom = rep2.load('lms_geometry')

        d = dict()
        set_params(d)
        create_fault_mesh(d, geom)
        create_surface_mesh(d, geom)
        bemify(d)
        boundary_conditions(d)
        build_G(d, qi_gps)
    if sys.argv[1] == 'solve':
        inv_d = rep2.load('bem_all_details_inverse_all_segs')
        L_plot(inv_d, qi_gps)
