import pyproj
import numpy as np
import scipy.io
import subprocess

def extract_tsurfs():
    # Run Brendan's tsurf parser
    cmd = ['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-r', 'run(\'ReadAndSaveCfm.m\');quit;']
    p = subprocess.Popen(cmd, cwd = './data/fault_geometry')
    p.wait()

def get_fault_data():
    # Load the .mat file
    matfile = scipy.io.loadmat('data/fault_geometry/Cfm1.2Pre_V7.mat')
    # Put the triangles into a nicer data form
    x = []
    y = []
    # z is depth
    z = []
    triangles = []
    # Walk to .mat file data structure to grab the values and
    # pull them out into a vertices and triangles structure
    cfm_tris = matfile['Cfm'][0][0][0]
    n_cfm_tris = cfm_tris.shape[0]
    for (idx, t) in enumerate(cfm_tris):
        x.extend(t['x'][0][0].tolist())
        y.extend(t['y'][0][0].tolist())
        z.extend(t['z'][0][0].tolist())
        triangles.append([3 * idx, 3 * idx + 1, 3 * idx + 2])
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    return x, y, z, triangles

def convert_fault_data_to_latlon(x, y):
    # Use Proj to convert from UTM48 to Latlon
    p = pyproj.Proj(proj = 'utm', zone = 48, ellps = 'WGS84')
    lon, lat = zip(*[p(x_v, y_v, inverse = True) for (x_v, y_v) in zip(x, y)])
    return lat, lon
