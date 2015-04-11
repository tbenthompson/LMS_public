import lms_code.lib.rep2 as rep2
import numpy as np
from codim1.core import *
from codim1.core.quad_strategy import AdaptiveInteriorQuad2
import matplotlib.pyplot as plt
from codim1.core.tools import plot_mesh
from lms_code.plots import plot_all as lms_plot
from mpi4py import MPI
import itertools
import sys

def plot_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x, y = zip(*points)
    ax.scatter(x, y, marker = 'o', c = 'b', s = 5)
    plt.show()

def get_field_type(type):
    if type == 'u':
        field_normal = [0.0, 0.0]
        field_type = 'displacement'
    elif type == 'sx':
        field_normal = [1.0, 0.0]
        field_type = 'traction'
    elif type == 'sy':
        field_normal = [0.0, 1.0]
        field_type = 'traction'
    else:
        raise Exception('Unknown field type')
    return field_normal, field_type

def pick_field_and_region(int_eval):
    print('Number of points to evaluate: ' + str(int_eval.meshpy_pts.shape[0]))
    int_eval.viz_mesh()
    region_num = raw_input('Type the region #: ')
    name = raw_input('Type the region name: ')
    type = raw_input('[u, sx, sy]: ')
    plt.close('all')
    return region_num, name, type

def create_quad(q_min):
    return AdaptiveInteriorQuad2(q_min, 1e50)

def one_chunk(model, region_num, subset_num, total_subsets, name, type, which_side):
    bem_soln = rep2.load('bem_' + model)
    int_eval = rep2.load('interior_mesh_' + model)

    field_normal, field_type = get_field_type(type)

    q_min = 4

    ek = ElasticKernelSet(bem_soln['shear_modulus'], bem_soln['poisson_ratio'])
    qs = create_quad(q_min)

    int_data = dict()
    int_data['region_num'] = region_num
    int_data['subset_num'] = subset_num
    int_data['total_subsets'] = total_subsets
    int_data['name'] = name
    int_data['type'] = type

    int_eval.choose_subregion(int(region_num), subset_num, total_subsets)
    int_data['tris'] = int_eval.selected_tris
    print 'Number of triangles: ' + str(int_data['tris'].shape)

    soln = int_eval.eval_soln(qs, ek, bem_soln['soln_coeffs'],
                              which_side, field_normal, field_type)
    int_data['soln'] = soln
    rep2.save(interior_fname(model, name, subset_num, type), int_data)

def interior_fname(model, name, subset_num, type):
    return 'interior_' + model + '_' + name + '/' + str(subset_num) + '_' + type

def chunkify(lst, n):
    return [lst[i::n] for i in xrange(n)]

def main(args):
    model = 'all_details'
    int_eval = rep2.load('interior_mesh_' + model)

    if len(args) == 2 and args[1] == 'choose':
        pick_field_and_region(int_eval)
        sys.exit()

    info = lms_plot.interior_params()
    options = [info['regions'], range(info['subsets']), info['types']]
    jobs = list(itertools.product(*options))


    comm = MPI.COMM_WORLD
    total_procs = comm.size
    chunks = chunkify(jobs, total_procs)
    my_chunk = chunks[comm.rank]
    if comm.rank == 0:
        print('There are ' + str(len(jobs)) + ' jobs to run.')
        print('Each of ' + str(comm.size) + ' processor is running ' + str(len(my_chunk)) + ' of the jobs.')
    for (region_num, name), subset_num, type in my_chunk:
        which_side = 'positive'
        if name == 'sichuan':
            which_side = 'negative'
        one_chunk(model, region_num, subset_num,
                  info['subsets'], name, type, which_side)


if __name__ == '__main__':
    main(sys.argv)
