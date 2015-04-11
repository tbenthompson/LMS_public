# Template taken from http://matplotlib.org/examples/pylab_examples/scatter_hist.html
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
import copy
import math
import numpy as np
import scipy
import random
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.ticker import NullFormatter


lms_plot.setup()
# the random data
model = 'all_details'
bem_soln = rep2.load('bem_' + model)
filename = 'bootstrap_' + model
# filename = 'bootstrap_alternate'
shortenings = rep2.load(filename)
depths = []
short_mults = []
every_distance = 30.0
replications = 15
min_d = 0
max_d = 23000
scale = 'contour'
plot_y = True
colored_1d = False

fig = plt.figure()

for e in bem_soln['fault_mesh']:
    if min_d < e.vertex1.loc[1] < -max_d:
        continue
    l = int(math.ceil(e.length / every_distance))
    x_hat = np.linspace(0.0, 1.0, l)
    normal = e.mapping.get_normal(0.0)
    angle = (math.pi / 2.0) - math.atan(normal[1] / -normal[0])
    shortening_mult = 1.0 / math.cos(angle)
    pts = [e.mapping.get_physical_point(x_h) for x_h in x_hat]
    short_mults.extend([shortening_mult] * len(pts))
    depths.extend([abs(p[1]) for p in pts])

# plt.plot(depths)
# plt.plot(short_mults)
# plt.show()

shortening_from_depth = []
replicated_depths = []
for d, s_mult in zip(depths, short_mults):
    sample = random.sample(shortenings, replications)
    replicated_depths.extend([d / 1000.0] * replications)
    shortening_from_depth.extend([s * s_mult for s in sample])

x = replicated_depths
y = shortening_from_depth

nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.6
bottom, height = 0.1, 0.6
bottom_h = left_h = bottom+height+0.05

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
cbar_zone = [left_h, bottom_h, 0.02, 0.2]

# start with a rectangular Figure
plt.figure(1, figsize=(8,8))

ax2D = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axCbar = plt.axes(cbar_zone)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)

#
x_min = min_d / 1000.0
x_max = max_d / 1000.0
x_steps = int((((max_d - min_d) / 1000.0) + 1) * 2.0)
xbins = np.linspace(x_min, x_max, x_steps)

y_min = 0.0
y_max = 14.0
y_steps = 30
ybins = np.linspace(y_min, y_max, y_steps)

cmap = lms_plot.get_summery_cmap()

if scale == 'log':
    ax2D.hist2d(x, y, bins = [xbins, ybins],
                    cmap = cmap,
                    norm = matplotlib.colors.LogNorm())#vmax = 135, vmin = 0)
    ax2D.set_ylabel('$\log(s)$ (mm/yr)')
elif scale == 'linear':
    ax2D.hist2d(x, y, bins = [xbins, ybins], cmap = cmap, vmax = 135, vmin = 0)
    ax2D.set_ylabel('$s$ (mm/yr)')
elif scale == 'scatter':
    xy = np.vstack([x, y])

    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = np.array(x)[idx], np.array(y)[idx], z[idx]
    ax2D.scatter(x, y, c = np.log(z), s=50, edgecolor='')#, cmap = cmap)
elif scale == 'contour':
    cxbins = np.linspace(x_min, x_max, 500)
    cybins = np.linspace(y_min, y_max, 500)
    min_density = 8e-4

    xy = np.vstack([x, y])

    # KERNEL DENSITY ESTIMATE HELPS FIND
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = np.array(x)[idx], np.array(y)[idx], z[idx]
    Z = scipy.interpolate.griddata((x, y), z,
                                  (cxbins[None, :], cybins[:, None]),
                                  fill_value = min_density)
    low_indices = Z < min_density
    Z[low_indices] = min_density

    cntr_opts = dict(
        norm = matplotlib.colors.LogNorm(vmin = np.min(Z), vmax = (10 * np.max(Z))),
        levels = min_density * (1.35 ** np.arange(0.0, 16.0))
    )
    cntf = ax2D.contourf(cxbins, cybins, Z, cmap = cmap, **cntr_opts)
    ax2D.contour(cxbins, cybins, Z, cmap = plt.cm.binary, **cntr_opts)
    cbar = plt.colorbar(cntf, cax = axCbar,
                 ticks = cntr_opts['levels'][[0, 3, 6, 9, 12, 15]],
                 format = "%.1e")
    cbar.set_label('$\log_{10}(N)$')

ax2D.set_xlabel('$d$ (km)')
ax2D.set_ylabel('$s$ (mm/yr)')

# now determine nice limits by hand:
ax2D.set_xlim( (min_d / 1000.0, max_d / 1000.0) )
ax2D.set_ylim( (0.0, 14.0) )


def color_hist(N, bins, patches):
    if colored_1d:
        #colors for the histogram
        fracs = N.astype(float)/N.max()
        norm = matplotlib.colors.Normalize(fracs.min(), fracs.max())

        for thisfrac, thispatch in zip(fracs, patches):
            color = cmap(norm(thisfrac))
            thispatch.set_facecolor(color)

default_color = '#777777'
N, bins, patches = axHistx.hist(x, bins=xbins, color = default_color)
color_hist(N, bins, patches)

axHistx.set_xlim(ax2D.get_xlim())
axHistx.set_ylabel('$N$')


ticks = [t for i, t in enumerate(axHistx.get_yticks()) if i % 2 == 0]
axHistx.set_yticks(ticks)

if plot_y:
    axHisty = plt.axes(rect_histy)
    axHisty.yaxis.set_major_formatter(nullfmt)
    N, bins, patches = axHisty.hist(y, bins=ybins, color = default_color, orientation='horizontal')
    axHisty.set_ylim(ax2D.get_ylim())
    ticks = [t for i, t in enumerate(axHisty.get_xticks()) if i % 3 == 0]
    axHisty.set_xticks(ticks)
    axHisty.set_xlabel('$N$')

    color_hist(N, bins, patches)


fig.set_size_inches(5.5, 5.5)
plt.savefig('depth_slip_' + scale + '.pdf')
