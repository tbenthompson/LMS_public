from scipy.io import loadmat
import numpy as np
import pyproj
import StringIO
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import Beach
from lms_code.lib.enjoy_dem import get_config, dem_interp_value
from lms_code.lib.extract_fault import *
from lms_code.lib.geometry import *
from lms_code.lib.gps import get_stations
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
from lms_code.plots.shading import set_shade, hillshade


rep2.cfg['repeat_enabled'] = False


def inside(lon, lat, llcrnr, urcrnr, looseness = 0.0):
    if llcrnr[0] - looseness <= lon <= urcrnr[0] + looseness and \
        llcrnr[1] - looseness <= lat <= urcrnr[1] + looseness:
        return True
    return False

def base(llcrnr, urcrnr):
    m = Basemap(projection='tmerc',
                lon_0=103.5,
                lat_0=31.5,
                resolution="c",
                llcrnrlon = llcrnr[0],
                llcrnrlat = llcrnr[1],
                urcrnrlon = urcrnr[0],
                urcrnrlat = urcrnr[1])
    return m

def even_grid(m, llcrnr, urcrnr, res_lon, res_lat):
    llcrnr_proj = m(*llcrnr)
    urcrnr_proj = m(*urcrnr)
    xs = np.linspace(llcrnr_proj[0], urcrnr_proj[0], res_lat)
    ys = np.linspace(llcrnr_proj[1], urcrnr_proj[1], res_lon)
    XS, YS = np.meshgrid(xs, ys)
    lons, lats = m(XS, YS, inverse = True)
    return lons, lats

def get_topo(m, llcrnr, urcrnr, res = 15):
    min_topo = -22000
    max_topo = 6000
    res_lat = 100 * res
    res_lon = 100 * res

    cfg = get_config('data/dem/e100n40.Bathymetry.srtm')
    srtm = np.load('data/dem/e100n40.npy')

    lons, lats = even_grid(m, llcrnr, urcrnr, res_lon, res_lat)

    dem = np.empty((res_lat, res_lon))
    for i in range(res_lat):
        for j in range(res_lon):
            dem[i, j] = dem_interp_value(srtm, cfg, (lats[i, j], lons[i, j]))
    rgb = set_shade(dem, cmap = plt.cm.Greys)
    return rgb

def topo(m, rgb = None):
    if rgb is None:
        rgb = get_topo(m, 15)
    m.imshow((1 - rgb) / 2.0 + 0.35)

def get_fault_trace(intersect_depth = 0.0):
    x, y, z, tris = get_fault_data()
    lat, lon = convert_fault_data_to_latlon(x, y)
    threeD_vertices = np.array(zip(lat, lon, z))
    verts, segs = intersection_surface_plane(threeD_vertices,
                                             tris,
                                             np.array([0.0, 0.0, intersect_depth]),
                                             np.array([0.0, 0.0, 1.0]))
    return verts, segs

def depth_contours(m):
    ds = np.linspace(-27000, 0000, 10)
    for d in ds:
        fault_trace(m, 'b', d, width = 0.6)

def fault_trace(m, color, depth = 0.0, width = 2.0):
    fault_verts, fault_segs = get_fault_trace(depth)
    fault_lats = [v[0] for v in fault_verts]
    fault_lons = [v[1] for v in fault_verts]
    fault_x, fault_y = m(fault_lons, fault_lats)
    for seg in fault_segs:
        x_seg = [fault_x[v] for v in seg]
        y_seg = [fault_y[v] for v in seg]
        m.plot(x_seg, y_seg, color = color, linestyle = '-', linewidth = width)

def beachballs(m):
    # mechanism from:
    # http://neic.usgs.gov/neis/eq_depot/2008/eq_080512_ryan/neic_ryan_cmt.html
    lats = [30.969]
    lons = [103.186]
    mech = [2, 47, 45]

    x, y = m(lons, lats)
    ax = plt.gca()
    b = Beach(mech, xy=(x[0], y[0]), width=3e4, linewidth=1, zorder = 101)
    ax.add_collection(b)

gps_opts = {
    'scale_units': 'inches',
    'scale': 24.0,
    'angles': 'xy',
    'width': 0.0075,
    'zorder': 100,
}

def unselected_gps(m, llcrnr, urcrnr):
    # Un-selected stations
    all_gps = loadmat('data/gps/CalaisGanSocquetBanerjeeApel_in_Calais_RF_trim_Clean_Tog.sta.data.mat')
    lats = all_gps['station_lat'][0]
    lons = all_gps['station_lon'][0]
    us = all_gps['station_east_vel'][0]
    vs = all_gps['station_north_vel'][0]
    # make sure every station we plot is inside the viewing region
    # it isn't bad to plot the outside stations except that the ends of
    # their arrows stick out into the plot and look bad.
    inview = filter(lambda (lon,lat,u,v): inside(lon, lat, llcrnr, urcrnr, looseness = 1.0),
                    zip(lons, lats, us, vs))
    lons, lats, us, vs = zip(*inview)
    xs, ys = m(lons, lats)
    Q = m.quiver(xs, ys, us, vs, color = '#222222', **gps_opts)

def selected_gps(m, geom):
    # Selected stations
    stations = get_stations(geom['which_gps_profile'])
    lat = [s['lat'] for s in stations]
    lon = [s['lon'] for s in stations]
    x, y = m(lon, lat)
    u = [s['east_vel'] for s in stations]
    v = [s['north_vel'] for s in stations]
    m.quiver(x, y, u, v, color = '#f03b20', **gps_opts)

def gps_scalebar(m):
    scalebar_length = 10
    loc = [103.60, 34.4]
    x, y = m(loc[0], loc[1])
    plt.text(x - 85000, y - 5000, str(scalebar_length) + ' mm/yr', fontsize = 10)
    m.quiver([x], [y], [scalebar_length], [0.0], color = '#0055CC', **gps_opts)

def labels(m, llcrnr, urcrnr):
    plt.text(310000, 20000, 'Sichuan Basin')
    plt.text(5000, 430000, 'Tibetan Plateau')
    plt.text(50000, 375000, 'Longriba Fault', fontsize = 8, rotation = 45)
    plt.text(322000, 324000, 'Beichuan', fontsize = 8, rotation = 45)
    plt.text(332000, 296000, 'Fault', fontsize = 8, rotation = 45)
    label_params = dict(
        linewidth = 0.0,
        labels = [1, 0, 0, 1]
    )
    m.drawparallels(np.linspace(llcrnr[1] + 1.0, urcrnr[1], 5),
                    fmt = lambda x: r'$%.0f ^{\circ}$'%x + 'N',
                    **label_params)
    m.drawmeridians(np.linspace(llcrnr[0], urcrnr[0], 6),
                    fmt = lambda x: r'$%.0f ^{\circ}$'%x + 'E',
                    **label_params)

def fault_shaded(m):
    x, y, z, tris = get_fault_data()
    lat, lon = convert_fault_data_to_latlon(x, y)
    x_map, y_map = m(lon, lat)
    import matplotlib
    tri = matplotlib.tri.Triangulation(x_map, y_map, tris)
    verts = np.concatenate((tri.x[tri.triangles][..., None],
                         tri.y[tri.triangles][..., None]), axis=2)
    color = (0.0, 0.0, 1.0, 0.15)
    collection = matplotlib.collections.PolyCollection(verts,
        facecolors = [color] * verts.shape[0])
    collection.set_edgecolor('none')
    ax = plt.gca()
    ax.add_collection(collection)
    # m.pcolor(np.array(x_map), np.array(y_map), np.array(tris), tri = True)

def other_faults(m):
    def close_to_beichuan_penguan(s):
        for lon, lat in s:
            d = (lat - 31) ** 2 + (lon - 103.5) ** 2
            if d < 0.1:
                return True
            d = (lat - 32) ** 2 + (lon - 105) ** 2
            if d < 0.1:
                return True
        return False

    with open('data/fault_geometry/taylor_yin_faults_only.xy', 'r') as f:
        text = f.read()
    data = np.genfromtxt(StringIO.StringIO(text), names = 'x, y', dtype = [float, float])
    segments = [[]]
    for lon,lat in zip(data['x'], data['y']):
        if np.isnan(lon) or np.isnan(lat):
            if len(segments[-1]) > 0:
                segments.append([])
            continue
        segments[-1].append((lon,lat))
    for s in segments:
        if close_to_beichuan_penguan(s):
            continue
        lons, lats = zip(*s)
        xs, ys = m(lons, lats)
        m.plot(xs, ys, 'k-', linewidth = 0.7)

def plot_xsec(m, geom):
    plt_params = lms_plot.params()

    g = pyproj.Geod(ellps='WGS84')

    pt1 = np.array(geom['basin_pt'])
    pt2 = np.array(geom['tibet_pt'])
    fwd, back, dist = g.inv(pt1[1], pt1[0], pt2[1], pt2[0])
    Xp = g.fwd(pt1[1], pt1[0], fwd, plt_params['view_left'] * 1000)
    X = g.fwd(pt2[1], pt2[0], back, (800 - plt_params['view_right']) * 1000)

    lon = [X[0], Xp[0]]
    lat = [X[1], Xp[1]]
    x, y = m(lon, lat)
    plt.plot(x, y, color = '#000000', linewidth = 2)
    shift = 3000
    plt.text(x[0] + shift, y[0], '$\\textrm{X}$')
    plt.text(x[1] + shift, y[1], '$\\textrm{X}^{\'}$')

def plot_map(steps = False):
    llcrnr = (101.0, 29.0)
    urcrnr = (106.0, 34.5)
    lms_plot.setup()
    m = base(llcrnr, urcrnr)
    dem = get_topo(m, llcrnr, urcrnr, 15)
    topo(m, dem)
    labels(m, llcrnr, urcrnr)
    if steps:
        plt.savefig('lms_map_only_topo')

    fault_shaded(m)
    other_faults(m)
    depth_contours(m)
    fault_trace(m, '#feb24c')
    beachballs(m)
    if steps:
        plt.savefig('lms_map_with_fault')

    geom = rep2.load('lms_geometry')
    unselected_gps(m, llcrnr, urcrnr)
    selected_gps(m, geom)
    gps_scalebar(m)
    if steps:
        plt.savefig('lms_map_with_gps')

    plot_xsec(m, geom)
    plt.savefig('lms_map_all')

if __name__ == "__main__":
    import plot_map
    rep2.repeat(plot_map, "plot_map", True)
