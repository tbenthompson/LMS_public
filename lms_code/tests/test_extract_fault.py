from lms_code.lib.extract_fault import *
from lms_code.lib.geometry import *
# from mayavi.mlab import triangular_mesh, show
import pyproj
import numpy as np

def test_extract_fault():
    x, y, z, triangles = get_fault_data()
    lat, lon = convert_fault_data_to_latlon(x, y)
    p = pyproj.Proj(proj='utm', zone = 48, ellps = 'WGS84')
    inv_x, inv_y = zip(*[p(x_v, y_v) for (x_v, y_v) in zip(lon, lat)])
    np.testing.assert_almost_equal(inv_x, x, 5)
    np.testing.assert_almost_equal(inv_y, y, 5)

def test_plot_tris():
    x, y, z, triangles = get_fault_data()
    lat, lon = convert_fault_data_to_latlon(x, y)
    # Plot them.
    # triangular_mesh(np.array(lon), np.array(lat), z / 60000.0, triangles)
    # show()

if __name__ == "__main__":
    test_plane_intersection()
