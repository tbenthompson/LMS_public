import numpy as np
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
from lms_code.plots.plot_interior import plot_surface, plot_fault
from lms_code.plots.plot_comparison import plot_lms_gps
import matplotlib.pyplot as plt

model = 'all_details'
geom = rep2.load('lms_geometry')
bem_soln = rep2.load('bem_' + model)
lms_plot.setup()
int_params = lms_plot.interior_params()
plot_params = lms_plot.params()

fig, axarr = plt.subplots(2, sharex = True)

plot_lms_gps(geom['gps_dist'], geom['gps_parallel_vel'],
             geom['gps_parallel_sig'], plot_params, axarr[0])

def plot_fault():
    angles = []
    for e in bem_soln['fault_mesh']:
        v1 = e.vertex1.loc
        v2 = e.vertex2.loc
        dx = v2[0] - v1[0]
        dy = v2[1] - v1[1]
        angles.append(np.arctan2(dy, dx))
    print np.degrees(np.min(angles))
    vs = [[e.vertex1.loc, e.vertex2.loc] for e in bem_soln['fault_mesh']]
    vs = np.sort(np.array([v for pair in vs for v in pair]), axis = 0)
    x = (vs[:, 0] - int_params['min_x']) / 1000.0
    y = vs[:, 1] / 1000.0
    axarr[1].plot(x, y, 'r-')

def plot_surf():
    x = (bem_soln['x'][0, :] - int_params['min_x']) / 1000.0
    y = bem_soln['x'][1, :] / 1000.0
    axarr[1].plot(x, y, 'k-')
    axarr[1].fill_between(x, [-25] * len(x), y, color = '#D8D8D8')

plot_fault()
plot_surf()

axarr[0].set_ylim([-2, 10])
axarr[0].set_ylabel('$v$ (mm/yr)')
axarr[1].set_xlim([0.0, (int_params['max_x'] - int_params['min_x']) / 1000.0])
axarr[1].set_ylim([-25, 5])
axarr[1].set_ylabel('$d$ (km)')
axarr[1].set_xlabel('$x$ (km)')
fig.set_size_inches(plot_params['fig_scale'])
plt.savefig('only_geometry', bbox_inches = 'tight')
