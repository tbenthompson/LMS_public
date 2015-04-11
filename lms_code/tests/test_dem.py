import numpy as np
from lms_code.lib.enjoy_dem import get_config, dem_interp_value, get_idx
from lms_code.lib.enjoy_dem import extract_line, get_distance

cfg = get_config('data/dem/e100n40.Bathymetry.srtm')
dem = np.load('data/dem/e100n40.npy')

def test_get_config():
    assert(cfg['cells_per_row'] == 4800)
    assert(cfg['rows'] == 6000)
    assert(cfg['lon_step'] == 0.008333333333333333)
    assert(cfg['lat_step'] == 0.008333333333333333)
    assert(cfg['min_lat'] == -10)
    assert(cfg['max_lat'] == 40)
    assert(cfg['min_lon'] == 100)
    assert(cfg['max_lon'] == 140)

def test_get_idx():
    lat = cfg['max_lat'] - cfg['lat_step'] * 20.5
    lon = cfg['min_lon'] + cfg['lat_step'] * 20.5
    np.testing.assert_almost_equal(get_idx(cfg, (lat, lon)), (20.5, 20.5))

def test_interp_value():
    a = dem[20, 20]
    b = dem[21, 20]
    c = dem[20, 21]
    d = dem[21, 21]
    average = (a + b + c + d) / 4.0
    lat = cfg['max_lat'] - cfg['lat_step'] * 20.5
    lon = cfg['min_lon'] + cfg['lat_step'] * 20.5
    value = dem_interp_value(dem, cfg, (lat, lon))
    np.testing.assert_almost_equal(value, average)

def test_interp_value2():
    a = dem[20, 20]
    b = dem[21, 20]
    average = (a + b) / 2.0
    lat = cfg['max_lat'] - cfg['lat_step'] * 20.5
    lon = cfg['min_lon'] + cfg['lat_step'] * 20.0
    value = dem_interp_value(dem, cfg, (lat, lon))
    np.testing.assert_almost_equal(value, average)

def test_interp_value_exact_point():
    lat = cfg['max_lat'] - cfg['lat_step'] * 80.0
    lon = cfg['min_lon'] + cfg['lat_step'] * 80.0
    value = dem_interp_value(dem, cfg, (lat, lon))
    assert(value == dem[80, 80])

def test_extract_line():
    start_lat = cfg['max_lat'] - cfg['lat_step'] * 20.0
    start_lon = cfg['min_lon'] + cfg['lat_step'] * 22.0
    end_lat = cfg['max_lat'] - cfg['lat_step'] * 22.0
    end_lon = cfg['min_lon'] + cfg['lat_step'] * 22.0
    values, dist, latlon_list = extract_line(dem, cfg, (start_lat, start_lon),
                                             (end_lat, end_lon), 3)
    np.testing.assert_almost_equal(dist[-1], 1.8532, 3)
    np.testing.assert_almost_equal(values, dem[20:23, 22])

def test_distance():
    start_lat = cfg['max_lat'] - cfg['lat_step'] * 20.0
    start_lon = cfg['min_lon'] + cfg['lat_step'] * 22.0
    end_lat = cfg['max_lat'] - cfg['lat_step'] * 22.0
    end_lon = cfg['min_lon'] + cfg['lat_step'] * 22.0
    d = get_distance((start_lat, start_lon), (end_lat, end_lon))
    np.testing.assert_almost_equal(d, 1.8532, 3)
