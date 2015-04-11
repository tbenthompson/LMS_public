import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from lms_code.lib import rep2
import lms_code.plots.plot_all as lms_plot
from lms_code.plots.plot_comparison import do_comparison_plot
from lms_code.plots.plot_interior import do_interior_plot, cntr_opts

def main():
    which_model = 'all_details'
    show_numbers = True
    bem_soln = rep2.load('bem_' + which_model)
    geom = rep2.load('lms_geometry')
    lsqr_est = rep2.load('shortening_estimate_' + which_model)
    plot_params = lms_plot.params()
    plot_params['which_model'] = which_model

    # Do the geometry and the comparison plot
    fig, axarr = plt.subplots(3, sharex = True)
    do_comparison_plot(bem_soln, geom, lsqr_est, plot_params, show_numbers, fig, axarr)

    # Do the interior plot
    cntr_opts[0] = lambda levels: {
        'levels': levels,
        'extend': 'neither',
        'linewidths': 0.7
    }

    def get_f(all):
        return all['interseis_ux']
    levels = np.linspace(-1.0, 7.0, 17)
    mask, colored = do_interior_plot(fig, axarr[2], which_model, get_f, levels,
                                     None, True, False)
    divider = make_axes_locatable(axarr[2])
    colorbar_ax = divider.append_axes("right", plot_params['cbar_width'], pad="3%")
    cbar = plt.colorbar(colored, cax = colorbar_ax)
    cbar.set_label('$\delta v$ (mm/yr)')
    axarr[2].set_ylabel('$d$ (km)')

    for i, label in enumerate(['(a)', '(b)', '(c)']):
        axarr[i].text(0.95, 0.9,
                      label,
                      transform = axarr[i].transAxes,
                      fontweight = 'bold')

    fig.set_size_inches([plot_params['fig_scale'][0], plot_params['fig_scale'][1] * 1.5])
    plt.savefig('stack_figure_' + plot_params['which_model'])

if __name__ == '__main__':
    main()
