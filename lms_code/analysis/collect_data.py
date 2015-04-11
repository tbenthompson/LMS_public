from lms_code.lib.enjoy_dem import extract_line, get_config, get_srtm_data
from lms_code.lib.extract_fault import *
from lms_code.lib.geometry import *
from lms_code.lib.gps import get_stations
from lms_code.lib.smooth import smooth
from lms_code.lib.vertex_fixer import vertex_fixer, process_fixes

import lms_code.lib.rep2 as rep2

import os
from pyproj import Geod
import matplotlib
from matplotlib import pyplot as plt

def set_params(d):
    d['which_gps_profile'] = 'properly_aligned3'

def get_plane():
    basin_pt = [29.439518, 106.485528]
    tibet_pt = [34.7327312, 100.8766929]
    # Make the plane
    lat_plane = [tibet_pt[0], tibet_pt[0], basin_pt[0], basin_pt[0]]
    lon_plane = [tibet_pt[1], tibet_pt[1], basin_pt[1], basin_pt[1]]
    z_plane = np.array([8000, -30000, 8000, -30000])
    triangles_plane = [(0, 1, 2), (2, 1, 3)]

    # Get the defining plane data.
    plane_pt1 = (lat_plane[0], lon_plane[0], z_plane[0])
    plane_pt2 = (lat_plane[1], lon_plane[1], z_plane[1])
    plane_pt3 = (lat_plane[2], lon_plane[2], z_plane[2])
    plane_n = plane_normal_from_three_points(plane_pt1, plane_pt2, plane_pt3)
    return tibet_pt, basin_pt, (plane_n, plane_pt1, plane_pt2, plane_pt3)

def set_plane(d):
    d['tibet_pt'], d['basin_pt'], d['xsec_plane'] = get_plane()

def dem_xsec(d):
    cfg = get_config('data/dem/e100n40.Bathymetry.srtm')
    if not os.path.exists('data/dem/e100n40.npy'):
        data = get_srtm_data('data/dem/e100n40.Bathymetry.srtm', cfg)
        np.save('data/dem/e100n40.npy', data)
    dem = np.load('data/dem/e100n40.npy')
    dem_cfg = cfg
    d['elevations'], d['elevation_dist'], d['elevation_latlon_list'] = \
        extract_line(dem, dem_cfg, d['tibet_pt'], d['basin_pt'], 500)

def tris_xsec(d):
    # Pull the fault data out of the tsurfs
    x, y, d['flt_z'], d['flt_triangles'] = get_fault_data()
    d['flt_lat'], d['flt_lon'] = convert_fault_data_to_latlon(x, y)
    threeD_vertices = np.array(zip(d['flt_lat'], d['flt_lon'], d['flt_z']))
    verts, segs = intersection_surface_plane(threeD_vertices,
                                             d['flt_triangles'],
                                             d['xsec_plane'][1],
                                             d['xsec_plane'][0])
    d['twod_vertices'] = verts
    d['twod_segments'] = segs

def project_gps(lon_3d, lat_3d, xsec_plane):
    proj_lons = []
    proj_lats = []
    for lon, lat in zip(lon_3d, lat_3d):
        xyz = (lat, lon, 0.0)
        (proj_lat, proj_lon, z), normal_project =\
            project_onto_plane(xyz, xsec_plane[1], xsec_plane[0])
        proj_lons.append(proj_lon)
        proj_lats.append(proj_lat)
    return proj_lons, proj_lats

def gps_xsec(d):
    # Grab the gps stations
    d['gps_stations'] = get_stations(d['which_gps_profile'])

    gps_lat_3d = [s['lat'] for s in d['gps_stations']]
    gps_lon_3d = [s['lon'] for s in d['gps_stations']]
    d['gps_lon'], d['gps_lat'] = project_gps(gps_lon_3d, gps_lat_3d, d['xsec_plane'])
    d['gps_parallel_vel'] = [s['parallel_vel'] for s in d['gps_stations']]

    # It appears that Brendan's tool output variances instead of std dev. So,
    # I take the sqrt to get 1-sigma
    d['gps_parallel_sig'] = np.sqrt([s['parallel_sigma'] for s in d['gps_stations']])

def smooth_elevations(d):
    half_window_size = 15
    d['elevations'] = smooth(d['elevations'], half_window_size * 2 + 1)
    d['elevations'] = d['elevations'][half_window_size:-half_window_size]

def remove_duplicate_verts(d):
    d['twod_vertices'] = {i:v for i, v in enumerate(d['twod_vertices'])}
    for (i, v1) in d['twod_vertices'].items():
        for (j, v2) in d['twod_vertices'].items():
            if j <= i:
                continue
            if sum(abs(v1 - v2)) > 0.01:
                continue

            del d['twod_vertices'][j]
            for (seg_idx, seg) in enumerate(d['twod_segments']):
                if seg[0] == j:
                    d['twod_segments'][seg_idx] = (i, seg[1])
                if seg[1] == j:
                    d['twod_segments'][seg_idx] = (seg[0], i)

def calculate_distances(d):
    g = Geod(ellps='WGS84')
    far_left_lat = d['elevation_latlon_list'][0][0]
    far_left_lon = d['elevation_latlon_list'][1][0]
    for i, v in d['twod_vertices'].iteritems():
        pt_lat = v[0]
        pt_lon = v[1]
        az1, az2, dist = g.inv(far_left_lon, far_left_lat, pt_lon, pt_lat)
        km_dist = dist / 1000.0
        d['twod_vertices'][i] = [v[0], v[1], v[2], km_dist]

    d['gps_dist'] = []
    for pt_lat, pt_lon in zip(d['gps_lat'], d['gps_lon']):
        az1, az2, dist = g.inv(far_left_lon, far_left_lat, pt_lon, pt_lat)
        km_dist = dist / 1000.0
        d['gps_dist'].append(km_dist)

def plot_2d_geometry(d):
    # Set some plotting style variables
    matplotlib.rcParams['lines.linewidth'] = 3

    # Should we plot latitude or longitude or distance
    latlonidx = 1

    # Put it all together in a nice plot
    fig, ax = plt.subplots(1)
    lines = []
    twod_vertices = d['twod_vertices']
    for seg in d['twod_segments']:
        v1 = (twod_vertices[seg[0]][latlonidx], twod_vertices[seg[0]][2])
        v2 = (twod_vertices[seg[1]][latlonidx], twod_vertices[seg[1]][2])
        lines.append((v1, v2))
    lc = matplotlib.collections.LineCollection(lines, colors=['r'])



    conv = matplotlib.colors.ColorConverter()
    gps_color = conv.to_rgba('#C3C3E5', alpha = 0.4)
    ax2 = ax.twinx()
    ax2.errorbar([d['gps_lat'], d['gps_lon'], 0, d['gps_dist']][latlonidx],
                  d['gps_parallel_vel'], yerr=d['gps_parallel_sig'],
                  fmt='o', color=gps_color)
    ax2.set_ylabel('Parallel GPS velocities (mm/yr)')
    ax2.set_ylim([0.0, 12.0])
    ax2.grid(False)

    ax.add_collection(lc)
    ax.plot(d['elevation_latlon_list'][latlonidx], d['elevations'], 'g')
    ax.set_xlabel(["Latitude", "Longitude", "", "Kilometers"][latlonidx] + " (degrees)")
    ax.set_ylabel("Elevation/Depth (meters)")
    ax.set_title("Longmenshan cross section")
    ax.set_xlim([[30, 35],[101.0, 106]][latlonidx])
    ax.grid(False)
    fig.set_size_inches([13.5, 7.5])
    filename = 'lms_geom'
    plt.savefig(filename + '.pdf')

    plt.figure()
    plt.plot(d['gps_lat'], d['gps_lon'], '.')
    plt.savefig('confirm_gps_projected.pdf')
    # plt.show()

if __name__ == "__main__":
    d = dict()
    set_params(d)
    set_plane(d)
    dem_xsec(d)
    tris_xsec(d)
    gps_xsec(d)
    smooth_elevations(d)
    remove_duplicate_verts(d)
    vertex_fixer(d)
    process_fixes(d)
    calculate_distances(d)
    plot_2d_geometry(d)
    rep2.save("lms_geometry", d)
