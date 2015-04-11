import numpy as np
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
from lms_code.lib.volume_mesh import VolumeMesh
import sys
sys.setrecursionlimit(50000)

def get_leftmost(mesh):
    leftmost_pt = np.array([1e9, 1e9])
    for e in mesh:
        for v in [e.vertex1, e.vertex2]:
            if v.loc[0] < leftmost_pt[0]:
                leftmost_pt = v.loc
    return leftmost_pt

def build_decoupled_boundary(mesh, int_params):
    leftmost_pt = get_leftmost(mesh)
    left_of_that = [int_params['min_x'], leftmost_pt[1]]
    the_edge = [0, 1]
    boundary = [[leftmost_pt, left_of_that], [the_edge]]
    return boundary


def mesh_interior(bem_soln, int_params):
    bem_mesh = bem_soln['combined_mesh']
    region = ((int_params['min_x'], int_params['max_x']),
              (int_params['min_y'], int_params['max_y'] + 1))
    decoupled = build_decoupled_boundary(bem_soln['fault_mesh'], int_params)
    vm = VolumeMesh(bem_mesh,
                    bem_soln['fault_mesh'],
                    region,
                    refine_length = int_params['edge_length_threshold'],
                    refine_area = int_params['tri_area_threshold'],
                    near_edge_factor = int_params['near_edge_refine'],
                    extra_edges = [decoupled]
                    )
    return vm

def main():
    model = 'all_details'
    bem_soln = rep2.load('bem_' + model)

    int_params = lms_plot.interior_params()
    vm = mesh_interior(bem_soln, int_params)
    int_eval = vm.get_evaluator()

    rep2.save('interior_mesh_' + model, int_eval)

if __name__ == "__main__":
    main()
