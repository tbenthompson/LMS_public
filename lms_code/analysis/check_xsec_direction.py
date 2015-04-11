from lms_code.lib import rep2
from lms_code.lib.geometry import *
from lms_code.lib.extract_fault import *
from lms_code.plots.plot_map import base
import matplotlib.pyplot as plt
import numpy as np

geom = rep2.load('lms_geometry')
t_pt = geom['tibet_pt']
b_pt = geom['basin_pt']
t_pt[1] -= 0.9

x, y, flt_z, flt_triangles = get_fault_data()
flt_lat, flt_lon = convert_fault_data_to_latlon(x, y)
threeD_vertices = np.array(zip(flt_lat, flt_lon, flt_z))
verts, segs = intersection_surface_plane(threeD_vertices,
                                         flt_triangles,
                                         np.array([0.0, 0.0, 0.0]),
                                         np.array([0.0, 0.0, 1.0]))
start = 300
end = 305
verts = np.array(verts)
# plt.plot(verts[:, 1], verts[:, 0], 'ro')
# plt.plot(verts[[start,end], 1], verts[[start, end], 0], 'o')
# plt.show()


step = 40
start_opts = range(1, verts.shape[0], step)
end_opts = range(0, verts.shape[0], step)

normals = []
for start in start_opts:
    for end in end_opts:
        vec = verts[end, :] - verts[start, :]
        normal = np.cross(vec, [0, 0, 1])
        normal /= np.linalg.norm(normal)
        if normal[0] < 0:
            normal = -normal
        normals.append(normal)

llcrnr = (101.0, 29.0)
urcrnr = (106.0, 34.5)
m = base(llcrnr, urcrnr)
surf_x, surf_y = m(verts[:, 1], verts[:, 0])
m.plot(surf_x, surf_y, 'o')
xsec_x, xsec_y = m([t_pt[1], b_pt[1]], [t_pt[0], b_pt[0]])
for normal in normals:
    xsec2_x, xsec2_y = m([b_pt[1] + normal[1] * 10, b_pt[1]],
                         [b_pt[0] + normal[0] * 10, b_pt[0]])
    m.plot(xsec2_x, xsec2_y, 'r-')
m.plot(xsec_x, xsec_y, 'b-')
plt.show()
