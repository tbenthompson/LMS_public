import copy
import numpy as np
import datetime
import lms_code.lib.rep2 as rep2
from lms_code.plots.plot_map import base, topo, gps_opts,\
    even_grid, fault_trace, get_topo, labels
from lms_code.plots import plot_all as lms_plot
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde

def throw_out(data, idx, reason):
    for k in data.keys():
        del data[k][idx]

def get_early(qigps):
    may15 = datetime.datetime(year = 2008, month = 5, day = 15)
    junefirst = datetime.datetime(year = 2008, month = 6, day = 1)
    june15 = datetime.datetime(year = 2008, month = 6, day = 15)
    early = [i for i,d in qigps['resurvey'].iteritems() if d < june15]
    return early

def get_list(d):
    return [d[k] for k in d.keys()]

def get_sublist(d, indices, inverse = False):
    if inverse:
        return [d[k] for k in d.keys() if k not in indices]
    return [d[i] for i in indices]

def plot_vector(m, lon, lat, u, v, scale):
    x, y = m(lon, lat)
    gps_opts['scale'] = scale
    m.quiver(x, y, u, v, **gps_opts)

    if scale > 50.0:
        loc = [103.60, 32.8]
    else:
        loc = [103.80, 29.2]
    x, y = m(loc[0], loc[1])
    plt.text(x - 32000, y - 2000, str(scale) + ' cm', fontsize = 10)
    m.quiver([x], [y], [scale], [0.0], color = '#0055CC', **gps_opts)

def plot_kde(m, lon, lat, u, v):
    # Okay for 2D data if value and pts and center are 2-vectors
    def gaussian(locs, values, centers, widths):
        l1 = (((locs - centers) / widths) ** 2)
        l2 = np.sqrt(np.sum(l1, axis = 1))
        return values * np.exp(-l2)[:, np.newaxis]

    def run():
        centers = np.vstack([lon, lat]).T
        data = np.vstack([u, v]).T
        n = 50
        width = 0.2
        eval_lons, eval_lats = even_grid(m, n, n)
        out_data = np.empty((2, n, n))
        for i in range(n):
            print i
            for j in range(n):
                loc = np.array([eval_lons[i, j], eval_lats[i, j]])
                effects = gaussian(loc, data, centers, width)
                influence = gaussian(loc, 1.0, centers, width)
                regress = np.sum(effects, axis = 0) / np.sum(influence, axis = 0)
                out_data[:, i, j] = regress
        return out_data
    out_data = rep2.run_or_load('qi_kde', run)
    m.imshow(out_data[0, :, :], cmap = plt.cm.jet)

def toss_large(q, thresh = 20):
    for i in q['EW'].keys():
        vel = np.sqrt(q['EW'][i] ** 2 + q['SN'][i] ** 2)
        if vel > thresh:
            throw_out(q, i, "Too fast!")

def makefig(settings):
    inv = settings[0]
    big = settings[1]
    lms_plot.setup()
    if big:
        llcrnr = (101.0, 29.0)
        urcrnr = (106.0, 34.5)
    else:
        llcrnr = (102.4, 30.5)
        urcrnr = (105.5, 33.0)

    qigps = rep2.load('qi_gps_imported')
    throw_out(qigps, 310, "It is before the earthquake?!")

    lon = get_list(qigps['lon'])
    lat = get_list(qigps['lat'])
    u = get_list(qigps['EW'])
    v = get_list(qigps['SN'])
    q = copy.copy(qigps)

    if big:
        toss_large(q)

    early = get_early(q)
    early_lon = get_sublist(q['lon'], early, inverse = inv)
    early_lat = get_sublist(q['lat'], early, inverse = inv)
    early_u = get_sublist(q['EW'], early, inverse = inv)
    early_v = get_sublist(q['SN'], early, inverse = inv)

    plt.figure()
    m = base(llcrnr, urcrnr)
    dem = get_topo(m, llcrnr, urcrnr, res = 6)
    topo(m, rgb = dem)
    scale = 200.0
    if big:
        scale = 20.0
    plot_vector(m, early_lon, early_lat, early_u, early_v, scale)
    # plot_kde(m, early_lon, early_lat, early_u, early_v)
    # plt.colorbar()
    fault_trace(m)
    root = 'qietal11'
    if big:
        root += '_big'
    if inv is True:
        root += '_late'
    else:
        root += '_early'
    plt.savefig(root + '.pdf')

def main():
    from multiprocessing import Pool
    p = Pool(4)
    p.map(makefig, [[True, True],
                    [False, True],
                    [True, False],
                    [False, False]])

if __name__ == "__main__":
    main()
