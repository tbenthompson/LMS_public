import numpy as np
import lms_code.plots.plot_all as lms_plot
from lms_code.plots.plot_comparison import plot_lms_gps
import lms_code.lib.rep2 as rep2
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime

june15 = datetime.datetime(year = 2008, month = 6, day = 15)
def get_early(qi_gps, qi_gps_all, before = june15, inv = False):
    dist = []
    vel = []
    sig = []
    for i, gps_idx in enumerate(qi_gps['selected']):
        if inv is False:
            if qi_gps_all['resurvey'][gps_idx] > before:
                continue
        else:
            if qi_gps_all['resurvey'][gps_idx] < before:
                continue

        dist.append(qi_gps['dist'][i] / 1000.0)
        vel.append(qi_gps['vel'][i])
        sig.append(qi_gps['sig'][i])
    return dist, vel, sig

def main():
    qi_gps_all = rep2.load('qi_gps_imported')
    qi_gps = rep2.load('qi_gps_near')
    plot_params = lms_plot.params()
    plot_params['rigid_offset'] = 0.0

    lms_plot.setup()
    fig, axarr = plt.subplots(1, sharex = True)
    axarr = [axarr]
    # divider = make_axes_locatable(axarr[0])
    # hidden_ax = divider.append_axes("right", plot_params['cbar_width'], pad="3%")
    # hidden_ax.set_visible(False)


    def plot(inv, ax, label):
        dist, vel, sig = get_early(qi_gps, qi_gps_all, inv = inv)
        color = '#FF0000'
        if inv:
            color = '#0000FF'
        ax.plot(np.array(dist) - plot_params['view_left'],
                   np.array(vel) + plot_params['rigid_offset'],
                   'o', color = color, label = label,
                   zorder = 1000)

    plot(False, axarr[0], 'Early')
    plot(True, axarr[0], 'Late')
    axarr[0].set_ylabel('$v_{\\textrm{both}}$ (cm)')
    # axarr[1].set_ylabel('$v_{\\textrm{late}}$ (cm)')
    axarr[0].set_xlabel('$x$ (km)')
    axarr[0].set_xlim([0, 550])
    # axarr[1].set_xlim([0, 550])
    zoomed = False
    if zoomed:
        axarr[0].set_ylim([-20, 20])
        axarr[1].set_ylim([-20, 20])
    plot_params['fig_scale'][1] /= 2.0
    plt.legend()
    fig.set_size_inches(plot_params['fig_scale'])
    filename = 'qi_gps_xsec.pdf'
    if zoomed:
        filename = 'qi_gps_xsec_zoomed.pdf'
    fig.savefig(filename, bbox_inches = 'tight')



if __name__ == '__main__':
    main()
