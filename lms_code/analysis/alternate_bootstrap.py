import numpy as np
import lib.rep2 as rep2
from lib.least_squares import build_least_squares
import lms_code.plot_all as lms_plot
from lms_code.powersets import cull_to_acceptable

bem_soln = dict()
geom = dict()
plot_params = lms_plot.params()

rep2.load('bem_all_details', bem_soln)
bem_soln['boundary'] = bem_soln['intersection_pt'][0] - 1
rep2.load('lms_geometry', geom)

samples = 1000
means = geom['gps_parallel_vel']
stddevs = geom['gps_parallel_sig']
data = np.empty((len(means), samples))
data = np.array([np.random.normal(means[i], stddevs[i], samples)
                 for i in range(len(means))])
print data.shape
results = []
for i in range(samples):
    print i
    A, b = build_least_squares(geom['gps_dist'],
                               data[:, i],
                               bem_soln['x'],
                               bem_soln['u_soln'],
                               plot_params['boundary'](bem_soln))
    shortening = np.linalg.lstsq(A, b)[0][0]
    print shortening
    results.append(shortening)

culled_sorted = cull_to_acceptable(results, [0, 18])
rep2.save('bootstrap_alternate', culled_sorted)
