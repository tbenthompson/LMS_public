import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
# from matplotlib.mlab import griddata
from scipy.interpolate import griddata

cmap = plt.cm.coolwarm
cntr_opts = [lambda levels: {
    'levels': levels,
    'extend': 'both',
    'linewidths': 0.7
}]

def plot_surface(ax, bem_soln, int_params, fill_between = False):
    x = (bem_soln['x'][0, :] - int_params['min_x']) / 1000.0
    y = bem_soln['x'][1, :] / 1000.0
    ax.plot(x, y, 'k-')
    if fill_between:
        ax.fill_between(x, [-25] * len(x), y, color = '#D8D8D8')

def plot_fault(ax, bem_soln, int_params, linestyle = 'k-'):
    vs = [[e.vertex1.loc, e.vertex2.loc] for e in bem_soln['fault_mesh']]
    vs = np.sort(np.array([v for pair in vs for v in pair]), axis = 0)
    x = (vs[:, 0] - int_params['min_x']) / 1000.0
    y = vs[:, 1] / 1000.0
    ax.plot(x, y, linestyle)

def triplotter(ax, field, f_x, f_y, tris, levels):
    triang = mtri.Triangulation(f_x, f_y, tris)

    cntrf = ax.tricontourf(triang, field, cmap = cmap, **cntr_opts[0](levels))
    cntr = ax.tricontour(triang, field,
                         colors = '#333333',
                         linestyles = 'solid',
                         **cntr_opts[0](levels))
    return cntrf

def ptplotter(ax, field, f_x, f_y, levels, int_params, bem_soln, mask):
    n = (600, 600)
    xi = np.linspace(0.0, int_params['max_x'] - int_params['min_x'], n[0])
    xi /= 1000.0
    yi = np.linspace(int_params['min_y'], int_params['max_y'], n[1])
    yi /= 1000.0
    xi, yi = np.meshgrid(xi, yi)

    input_mask = np.logical_or(np.logical_not(np.isinf(field)), np.abs(field) < 1e-15)
    f_x = f_x[input_mask]
    f_y = f_y[input_mask]
    field = field[input_mask]

    zi = griddata((f_x, f_y), field, (xi, yi),
                  method = 'linear', fill_value = 0)

    if mask is None:
        mask = np.zeros(n)
        for e in bem_soln['surf_mesh']:
            x, y = e.vertex1.loc
            x -= int_params['min_x']
            x /= 1000.0
            if x < 0 or x > int_params['max_x'] / 1000.0:
                continue
            y /= 1000.0
            in_range = np.abs(xi - x) < 3
            higher_than_surf = yi > y
            mask = mask + np.all([in_range, higher_than_surf], axis = 0)
        mask = np.where(mask > 0, True, False)

    xi = np.ma.array(xi, mask = mask)
    yi = np.ma.array(yi, mask = mask)
    zi = np.ma.array(zi, mask = mask)

    cntrf = ax.contourf(xi, yi, zi, mask = mask,
                        cmap = cmap,
                        **cntr_opts[0](levels))
    cntr = ax.contour(xi, yi, zi, mask = mask,
                      colors = '#333333',
                      linestyles = 'solid',
                      **cntr_opts[0](levels))

    return mask, cntrf

def do_interior_plot(fig, ax, model, field_grabber, levels, mask, disc, plot_tris):
    all = rep2.load('interior_' + model + '_coalesced')
    bem_soln = rep2.load('bem_' + model)
    int_params = lms_plot.interior_params()

    field = field_grabber(all)
    if disc:
        x = all['disc_x']
        y = all['disc_y']
        tris = all['disc_tris']['tibet'] + all['disc_tris']['sichuan']
    else:
        x = all['x']
        y = all['y']
        tris = all['tris']

    x -= int_params['min_x']
    x /= 1000.0
    y /= 1000.0


    if disc or plot_tris:
        cntrf = triplotter(ax, field, x, y, tris, levels)
    else:
        mask, cntrf = ptplotter(ax, field, x, y, levels, int_params, bem_soln, mask)

    plot_surface(ax, bem_soln, int_params)
    plot_fault(ax, bem_soln, int_params)

    ax.set_ylim([int_params['min_y'] / 1000.0, int_params['max_y'] / 1000.0])
    ax.set_xlim([0.0, (int_params['max_x'] - int_params['min_x']) / 1000.0])
    return mask, cntrf

def get_energy_density(all):
    return all['energy_density']

def post_interseis(fig, ax, cbar, filename):
    cbar.set_label('$\delta_{\\textrm{v}}$ (mm/yr)')
    post_default(fig, ax, cbar, filename)

def post_log_energy_density(fig, ax, cbar, filename):
    cbar.set_label('$\log_{10}(E)$')
    post_default(fig, ax, cbar, filename)

def post_energy_density(fig, ax, cbar, filename):
    cbar.set_label('$E$ (Pa/yr)')
    post_default(fig, ax, cbar, filename)

def post_default(fig, ax, cbar, filename):
    ax.set_ylabel('$d$ (km)')
    ax.set_xlabel('$x$ (km)')
    plot_params = lms_plot.params()
    plot_params['fig_scale'][1] /= 2.0
    fig.set_size_inches(plot_params['fig_scale'])
    plt.savefig(filename, bbox_inches = 'tight')

def main():
    model = 'all_details'
    stress_steps = 21
    stress_levels = np.linspace(-1e6, 1e6, 21)
    stress_levels = np.delete(stress_levels, (stress_steps - 1) / 2)

    strain_steps = 41
    strain_levels = np.linspace(-2e-5, 2e-5, 41)
    strain_levels = np.delete(strain_levels, (strain_steps - 1) / 2)

    levels = dict()
    levels['interseis_ux'] = np.linspace(-1.0, 7.0, 17)
    levels['ux'] = levels['uy'] = np.linspace(-0.4, 1.1, 30)
    levels['sxx'] = levels['syy'] = levels['sxy'] = stress_levels
    levels['exx'] = levels['exy'] = levels['eyy'] = strain_levels
    levels['energy_density'] = np.linspace(0, 0.00050, 21)
    levels['log_energy_density'] = np.linspace(-16, -2, 15)

    disc = dict()
    disc['interseis_ux'] = disc['ux'] = disc['uy'] = True

    log = dict()
    log['log_energy_density'] = True

    field = dict()
    field['log_energy_density'] = 'energy_density'

    postproc = dict()
    postproc['interseis_ux'] = post_interseis
    postproc['log_energy_density'] = post_log_energy_density
    postproc['energy_density'] = post_energy_density

    fields = []
    # fields.extend(['ux', 'uy'])
    # fields.extend(['sxx', 'syy', 'sxy'])
    # fields.extend(['exx', 'exy', 'eyy'])
    fields.append('interseis_ux')
    # fields.append('energy_density')
    fields.append('log_energy_density')
    mask = None
    for f in fields:
        def get_f(all):
            field_name = field.get(f, f)
            if log.get(f, False):
                return np.log(all[field_name])
            return all[field_name]
        print("Making " + f)
        fig, ax = plt.subplots(1)
        mask, colored = do_interior_plot(fig, ax, model, get_f, levels[f],
                    mask, disc.get(f, False), False)
        cbar = plt.colorbar(colored)
        post_fnc = postproc.get(f, post_default)
        post_fnc(fig, ax, cbar, f)


if __name__ == '__main__':
    lms_plot.setup()
    main()
