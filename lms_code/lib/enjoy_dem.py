from math import floor, ceil, radians, sin, cos, atan2, sqrt
import numpy as np
import re
import struct

def get_config(filename):
    """Parse some of the relevant descriptors of the dem file."""
    config = dict()
    with open(filename + '.ers') as f:
        line = f.readline()
        while line:
            if 'NrOfCellsPerLine' in line:
                m = re.search(r'([0-9]+)', line)
                config['cells_per_row'] = int(m.group())

            if 'NrOfLines' in line:
                m = re.search(r'([0-9]+)', line)
                config['rows'] = int(m.group())

            if 'Xdimension' in line:
                m = re.search(r'0.([0-9]+)', line)
                config['lon_step'] = float(m.group())

            if 'Ydimension' in line:
                m = re.search(r'0.([0-9]+)', line)
                config['lat_step'] = float(m.group())

            # I assume that the upper left coordinates are full degrees
            # no minutes or seconds. Easy to fix if necessary.
            if 'Latitude' in line:
                m = re.search(r'[0-9]+', line)
                config['max_lat'] = float(m.group())
            if 'Longitude' in line:
                m = re.search(r'[0-9]+', line)
                config['min_lon'] = float(m.group())

            line = f.readline()
    config['min_lat'] = config['max_lat'] - \
                        config['lat_step'] * config['rows']
    config['max_lon'] = config['min_lon'] + \
                        config['lon_step'] * config['cells_per_row']
    return config

def get_srtm_data(filename, cfg):
    """
    This function reads in a srtm30 file from Dave Sandwell's website at:
    ftp://topex.ucsd.edu/pub/srtm30_plus/

    An alternative would be to use the gtopo function. However, my bem code
    is written in python, thus a simple python reader is handy.
    """
    cells_per_row = cfg['cells_per_row']
    rows = cfg['rows']
    print cells_per_row, rows
    data = np.empty((rows, cells_per_row))
    idx_row = 0
    idx_col = 0
    with open(filename, 'rb') as f:
        bytes = f.read(2)
        while bytes:
            # bizarrely, the srtm30 data is in big-endian form?!
            val = struct.unpack(">h", bytes)[0]
            # if the value is really big, that means the data is missing
            # I just uniformly code this as -1
            if val > 60000:
                val = -1
            # Just a check to make sure we are reading real earth
            # elevation data. the maximum elevation is less than 10000.
            assert(val < 10000)
            data[idx_row, idx_col] = val
            idx_col += 1
            if idx_col >= cells_per_row:
                idx_row += 1
                idx_col = 0
                if idx_row % 100 == 0:
                    print 'Row: ' + str(idx_row)
            bytes = f.read(2)
    return data


def dem_interp_value(dem, cfg, lat_lon):
    """
    Perform a 2d linear interpolation between adjacent points to
    allow sampling of arbitrary latitude and longitude positions.
    """
    lat_idx, lon_idx = get_idx(cfg, lat_lon)
    if lat_idx < 0 or lon_idx < 0:
        raise Exception("(dem_interp_value)The requested point is outside" +
                        " the provided DEM's coverage.")

    pts = []
    pts.append((floor(lat_idx), floor(lon_idx)))
    pts.append((ceil(lat_idx), floor(lon_idx)))
    pts.append((ceil(lat_idx), ceil(lon_idx)))
    pts.append((floor(lat_idx), ceil(lon_idx)))

    pts_dist = []
    for p in pts:
        d = (abs(lat_idx - p[0]), abs(lon_idx - p[1]))
        pts_dist.append(d)
        # If the distance to a point is effectively zero, then just
        # return the value at that point.
        if abs(d[0]) < 1e-4 and abs(d[1]) < 1e-4:
           return dem[int(p[0]), int(p[1])]

    value = 0
    for p, d in zip(pts, pts_dist):
        value += dem[int(p[0]), int(p[1])] * (1 - d[0]) * (1 - d[1])
    return value

def get_idx(cfg, lat_lon):
    """
    Convert from latitude and longitude into dem index coordinates.
    """
    lat_idx = (-lat_lon[0] + cfg['max_lat']) / cfg['lat_step']
    lon_idx = (lat_lon[1] - cfg['min_lon']) / cfg['lon_step']
    return (lat_idx, lon_idx)

def extract_line(dem, cfg, start_lat_lon, end_lat_lon, num_pts):
    """
    Basically a 2d linspace for a DEM.
    Runs a straight line in latitude, longitude space from
    start_lat_lon to end_lat_lon with num_pts sampled.
    """
    lat_list = np.linspace(start_lat_lon[0], end_lat_lon[0], num_pts)
    lon_list = np.linspace(start_lat_lon[1], end_lat_lon[1], num_pts)
    values = []
    dist = []
    last_latlon = None
    for (lat, lon) in zip(lat_list, lon_list):
        values.append(dem_interp_value(dem, cfg, (lat, lon)))
        if last_latlon is None:
            dist.append(0)
        else:
            dist.append(dist[-1] + get_distance(last_latlon, (lat, lon)))
        last_latlon = (lat, lon)
    latlon_list = (lat_list, lon_list)
    return np.array(values), dist, latlon_list


def get_distance(lat_lon1, lat_lon2):
    """
    Returns the great circle distance between two points in kilometers
    """
    R = 6371;
    phi1 = radians(lat_lon1[0])
    phi2 = radians(lat_lon2[0])
    dphi = radians(lat_lon2[0]-lat_lon1[0])
    dlam = radians(lat_lon2[1]-lat_lon1[1])

    a = sin(dphi/2) * sin(dphi/2) + \
        cos(phi1) * cos(phi2) * sin(dlam/2) * sin(dlam/2)
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    d = R * c;
    return d
