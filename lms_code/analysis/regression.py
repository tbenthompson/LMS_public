from matplotlib import pyplot as plt
import numpy as np
import lms_code.lib.rep2 as rep2

geom = rep2.load('lms_geometry')
gps_x = geom['gps_dist']
gps_v = geom['gps_parallel_vel']

def gaussian(pts, value, center, width):
    return value * np.exp(-(((pts - center) / width) ** 2))

pts = 100
smooth = 50.0

x = np.linspace(0, 800.0, pts)
effects = np.array([gaussian(x, vel, loc, smooth) for (loc, vel) in zip(gps_x, gps_v)])
influence = np.array([gaussian(x, 1.0, loc, smooth) for loc in gps_x])
regression = np.sum(effects, axis = 0) / np.sum(influence, axis = 0)

plt.plot(gps_x, gps_v, '.')
plt.plot(x, regression)
plt.show()
