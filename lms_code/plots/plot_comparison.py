import lms_code.lib.rep2 as rep2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import lms_code.plots.plot_all as lms_plot
from codim1.core.tools import plot_mesh
from codim1.post import evaluate_solution_on_element
from functools import partial

def block_vels(x, coseis_x, coseis_y, boundary, left_block, right_block, mask_offset = 1):
    block = np.where(x <= boundary, left_block, right_block)
    interseismic = block - coseis_x

    mask = np.ones(x.shape[0], dtype = bool)
    mask_me = np.where(x == [xv for xv in x if xv > boundary][0])[0] - mask_offset
    mask[mask_me] = False
    return x[mask], interseismic[mask], coseis_y[mask]

def bem_to_view_space(plot_params, x):
    return (np.array(x) / 1000.0) - (plot_params['view_left'], 0)

def plot_surface(plot_params, mesh, ax):
    pts_per_element = 5
    color = '#000000'
    x_hat = np.linspace(0.0, 1.0, pts_per_element)
    all_pts = []
    for e in mesh:
        points = map(e.mapping.get_physical_point, x_hat)
        points = map(partial(bem_to_view_space, plot_params), points)
        all_pts.extend(points)
    x = [p[0] for p in all_pts]
    y = [p[1] for p in all_pts]
    ax.plot(x, y, color = color)
    ax.fill_between(x, [-25] * len(x), y, color = '#D8D8D8')

def plot_fault(bem_soln, plot_params, shortening, mesh, fig, ax):
    def displacement_magnitude(e, x_hat):
        u = evaluate_solution_on_element(e, x_hat, bem_soln['soln_coeffs'])[0]
        u *= 2
        u_mag = np.sqrt(u[0] ** 2 + u[1] ** 2)
        return u_mag

    pts_per_element = 10
    x_hat = np.linspace(0.0, 1.0, pts_per_element)
    x_hat_extra = -0.1
    all_lines = []
    disp_mags = []
    for e in mesh:
        points = map(e.mapping.get_physical_point, x_hat)
        points = map(partial(bem_to_view_space, plot_params), points)
        lines = zip(points[:-1], points[1:])
        all_lines.extend(lines)
        # The left hand point of a line segment controls the width
        disp_mags.extend(map(partial(displacement_magnitude, e), x_hat)[:-1])

    linewidths = [(u - 0.8) * 32 for u in disp_mags]

    #http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
    lc = matplotlib.collections.LineCollection(all_lines,
                                               linewidths = linewidths,
                                               cmap = plt.get_cmap('coolwarm'))
    # lc.set_clim(min(disp_mags), max(disp_mags))
    lc.set_array(shortening * np.array(disp_mags))
    ax.add_collection(lc)

    divider = make_axes_locatable(ax)
    colorbar_ax = divider.append_axes("right", plot_params['cbar_width'], pad="3%")
    cbar = plt.colorbar(lc, cax=colorbar_ax)
    cbar.set_label('$s$ (mm/yr)')

def plot_lms_mesh(bem_soln, plot_params, shortening, fig, ax):
    plot_surface(plot_params, bem_soln['surf_mesh'], ax)
    plot_fault(bem_soln, plot_params, shortening, bem_soln['fault_mesh'], fig, ax)
    ax.set_xlim(0, plot_params['view_right'] - plot_params['view_left'])
    ax.grid(False)

def plot_lms_vels(bem_soln, plot_params, ax, slip, tibet,
                  sichuan, idx, label = "", plot_y = True,
                  error = None):
    # Calculate the plotting velocities from the coseismic
    # and the long term block velocities.
    # interseismic = block - coseismic
    x, interseismic, uy = block_vels(bem_soln['x'][0, :],
                                     slip * bem_soln['u_soln'][0, :],
                                     slip * bem_soln['u_soln'][1,:],
                                     plot_params['boundary'](bem_soln),
                                     tibet,
                                     sichuan,
                                     mask_offset = 0)
    x /= 1000.0

    color = ['r-', 'b-', 'g-', 'c-', 'm-', 'g-', 'b-'][idx]
    shifted_x = x - plot_params['view_left']
    # plot error window before the main line so it's behind.
    if error is not None:
        error *= 1.0
        low_error = interseismic - error
        high_error = interseismic + error
        ax.fill_between(shifted_x, low_error, high_error, facecolor = '#AAAAAA')

    ax.plot(shifted_x, interseismic, color, label = label)

    if plot_y:
        ax.plot(shifted_x, uy, 'g-')
    ax.set_xlim([0, plot_params['view_right'] - plot_params['view_left']])
    ax.grid(False)


def plot_lms_gps(gps_dist, gps_vel, gps_sig, plot_params, ax, color = '#000000'):
    conv = matplotlib.colors.ColorConverter()
    # Old version
    # gps_color = conv.to_rgba('#C3C3E5', alpha = 0.4)
    gps_color = conv.to_rgba(color)
    for x, v, sig in zip(gps_dist, gps_vel, gps_sig):
        xs = [x - plot_params['view_left'], x - plot_params['view_left']]
        ys = [v - sig + plot_params['rigid_offset'], v + sig + plot_params['rigid_offset']]
        ax.plot(xs, ys, '-', color = gps_color, zorder = 1000, linewidth = 2.0)
    ax.plot(np.array(gps_dist) - plot_params['view_left'],
                 np.array(gps_vel) + plot_params['rigid_offset'],
                 'o', color = gps_color, mfc = 'white', mew = 1.5,
                 zorder = 1000)

def main(bem_soln, geom, lsqr_est, plot_params, show_numbers):
    fig, axarr = plt.subplots(2, sharex = True)
    do_comparison_plot(bem_soln, geom, lsqr_est, plot_params, show_numbers, fig, axarr)
    fig.set_size_inches(plot_params['fig_scale'])
    plt.savefig('comparison_plot_' + plot_params['which_model'])

def do_comparison_plot(bem_soln, geom, lsqr_est, plot_params, show_numbers, fig, axarr):
    shortening = lsqr_est['lsqr_shortening']
    shortening_error = lsqr_est['lsqr_shortening_error']
    tibet = lsqr_est['lsqr_tibet'] + plot_params['rigid_offset']
    sichuan = tibet - shortening

    lms_plot.setup()
    slip_options = np.linspace(0, 12, 7)
    divider = make_axes_locatable(axarr[0])
    hidden_ax = divider.append_axes("right", plot_params['cbar_width'], pad="3%")
    hidden_ax.set_visible(False)

    plot_lms_mesh(bem_soln, plot_params, shortening, fig, axarr[1])
    plot_lms_vels(bem_soln, plot_params, axarr[0], shortening,
                  tibet, sichuan, 0, plot_y = False,
                  error = shortening_error)
    plot_lms_gps(geom['gps_dist'], geom['gps_parallel_vel'],
                 geom['gps_parallel_sig'], plot_params, axarr[0])

    axarr[0].set_ylim([-2, 10])
    axarr[0].set_ylabel('$v$ (mm/yr)')
    axarr[1].set_ylim([-25, 5])
    axarr[1].set_ylabel('$d$ (km)')
    axarr[-1].set_xlabel('$x$ (km)')
    if show_numbers:
        axarr[1].text(265, -0.5, '9.0', fontsize = 12)
        axarr[1].text(260, -8, '8.0', fontsize = 12)
        axarr[1].text(248, -16, '7.0', fontsize = 12)
        axarr[1].text(202, -20.5, '6.0', fontsize = 12)

        axarr[0].text(0.50, 0.9,
                      '$\\delta v = %.1f \\pm %.1f$ mm/yr'%(shortening, shortening_error),
                      transform = axarr[0].transAxes,
                      fontweight = 'bold')


if __name__ == "__main__":
    import plot_comparison
    # which_model = 'coarse'
    which_model = 'just_thrust'
    # which_model = 'simple'
    # show_numbers = False
    # which_model = 'all_details'
    show_numbers = True
    bem_soln = rep2.load('bem_' + which_model)
    geom = rep2.load('lms_geometry')
    lsqr_est = rep2.load('shortening_estimate_' + which_model)
    plot_params = lms_plot.params()
    plot_params['which_model'] = which_model
    rep2.repeat(plot_comparison, "plot1", bem_soln, geom, lsqr_est, plot_params, show_numbers)
