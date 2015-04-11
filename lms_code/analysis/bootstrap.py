import numpy as np
from itertools import combinations
import lms_code.lib.rep2 as rep2
from lms_code.lib.least_squares import build_least_squares

def bootstrap(model):
    bem_soln = dict()
    geom = dict()

    rep2.load('bem_' + model, bem_soln)
    bem_soln['boundary'] = bem_soln['intersection_pt'][0] - 1
    rep2.load('lms_geometry', geom)
    import ipdb; ipdb.set_trace()

    print "Calculating least squares estimates for powerset."
    A, b = build_least_squares(geom['gps_dist'],
                               geom['gps_parallel_vel'],
                               geom['gps_parallel_sig'],
                               bem_soln['x'],
                               bem_soln['u_soln'],
                               bem_soln['boundary'])
    n_obs = len(b)
    powerset = []
    # We want any set that has at least one element
    for i in range(2, n_obs + 1):
        powerset.extend(combinations(range(n_obs), i))
    powerset = [list(s) for s in powerset]
    results = []
    for i, p in enumerate(powerset):
        cur_A = A[p, :]
        cur_b = b[p]
        shortening = np.linalg.lstsq(cur_A, cur_b)[0][0]
        results.append(shortening)
        if i % 10000 == 0:
            print i
    return results

def cull_to_acceptable(results, acceptable):
    print "Culling Data outside reasonable range (" + \
          str(acceptable[0]) + \
          ", " + str(acceptable[1]) + ") mm/yr."

    return np.array([res for res in sorted(results)
            if acceptable[0] < res < acceptable[1]])

if __name__ == "__main__":
    model = 'all_details'
    results = bootstrap(model)
    culled_sorted = cull_to_acceptable(results, [0, 18])
    rep2.save('bootstrap_' + model, culled_sorted)
