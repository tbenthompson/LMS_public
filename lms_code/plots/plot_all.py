import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from copy import copy

def get_summery_cmap():
    cmap = copy(plt.cm.summer)
    cmap_start = cmap(0.0)
    cmap_end = cmap(1.0)
    cmap_start = [c_v + 0.03 for c_v in cmap_start]
    cmap_end = [c_v for c_v in cmap_end]
    # remove the blue from the high end of the scale
    cmap_end[2] = 0.0
    cdict = {
        'red': ((0.0, cmap_start[0], cmap_start[0]),
                (1.0, cmap_end[0], cmap_end[0])),
        'green': ((0.0, cmap_start[1], cmap_start[1]),
                  (1.0, cmap_end[1], cmap_end[1])),
        'blue': ((0.0, cmap_start[2], cmap_start[2]),
                 (1.0, cmap_end[2], cmap_end[2]))
    }
    cmap = matplotlib.colors.LinearSegmentedColormap('BenCMap', cdict)
    cmap.set_bad(cmap(0.0))
    return cmap

def params():
    import lms_code.lib.rep2 as rep2
    rep2.cfg['repeat_enabled'] = False
    p = dict()
    p['refine_factor'] = 10
    p['refine_criteria'] = lambda e: e.bc.type == 'crack_displacement'
    p['interior_levels'] = np.linspace(-0.8, 0.8, 20)

    p['view_left'] = 150
    p['view_right'] = 700
    p['fig_scale'] = [6.4, 4.8]
    # Location of the block boundary in meters
    p['boundary'] = lambda bem_soln: bem_soln['intersection_pt'][0] - 1
    p['slip_magnitude'] = 8.0
    p['rigid_offset'] = -0.6
    p['cbar_width'] = '3%'
    return p

def interior_params():
    p = dict()
    p['edge_length_threshold'] = 500000.0

    p['min_x'] = 150000
    p['max_x'] = 700000
    p['min_y'] = -25000
    p['max_y'] = 6000
    p['regions'] = [(413, 'tibet'), (414, 'sichuan')]
    p['tri_area_threshold'] = 40000
    p['near_edge_refine'] = 2.0
    p['types'] = ['u', 'sx', 'sy']
    p['subsets'] = 85

    # p['min_x'] = 440000
    # p['max_x'] = 480000
    # p['min_y'] = -3000
    # p['max_y'] = 3000
    # p['regions'] = [(46, 'tibet'), (41, 'sichuan')]
    # p['tri_area_threshold'] = 100000
    # p['near_edge_refine'] = 9.0

    # p['min_x'] = 458000
    # p['max_x'] = 462000
    # p['min_y'] = 0000
    # p['max_y'] = 2000
    # p['regions'] = [(24, 'tibet'), (21, 'sichuan')]
    # p['tri_area_threshold'] = 500000
    # p['near_edge_refine'] = 100.0
    # p['types'] = ['sy']
    # p['subsets'] = 6
    return p

def setup_poster():
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['font.size'] = 18
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    matplotlib.rcParams['lines.linewidth'] = 3
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.rcParams['savefig.format'] = 'pdf'
    matplotlib.rcParams['pdf.fonttype'] = 42

def setup():
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['font.size'] = 12
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    matplotlib.rcParams['lines.linewidth'] = 3
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.rcParams['savefig.format'] = 'pdf'
    matplotlib.rcParams['pdf.fonttype'] = 42
    # matplotlib.rcParams['svg.fonttype'] = 'path'
    # matplotlib.rcParams['svg.image_inline'] = True
