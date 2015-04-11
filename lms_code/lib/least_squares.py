import numpy as np
from copy import copy

def build_least_squares(gps_x, gps_v, gps_err, x_soln, v_soln, block_boundary):
    """
    Build the least squares problem Ax = b where x is a vector consisting of
    three elements --
    x[0] -- the shortening rate
    x[1] -- the tibet profile parallel velocity
    x[2] -- the sichuan profile parallel velocity
    So for a station in Tibet:    x[1] - x[0] * m = G
    and for a station in Sichuan: x[2] - x[0] * m = G
    where m is the boundary element calculated motion for the coseismic
    motion and G is the measured GPS velocity.
    But, we actually want x[1] - x[0] = x[2] for kinematic consistency, so
    condensing that relation out of the equations, we get:
    for a station in Tibet:    x[1] - x[0] * m       = G
    for a station in Sichuan:  x[1] - x[0] * (m + 1) = G

    Solve it later using np.linalg.lstsq(A, b).
    """
    b = copy(gps_v)
    gps_x, b = zip(*sorted(zip(gps_x, b)))
    remove = 0
    b = b[remove:]
    A = np.zeros((len(b), 2))

    for i, pos in enumerate(gps_x[remove:]):
        pos_m = pos * 1000
        x_chosen = 1e18

        # position less than the boundary and we are in Tibet
        # position greater than the boundary and we are in Sichuan
        sign = 1
        if pos_m <= block_boundary:
            A[i, 1] = 1
        else:
            # Sichuan = Tibet - Shortening
            # otherwise the model is kinematically inconsistent
            A[i, 1] = 1
            A[i, 0] = -1
        for j, x in enumerate(x_soln[0, :]):
            if abs(x - pos_m) < 250 and abs(x - x_chosen) > 500:
                A[i, 0] -= v_soln[0, j]
                # d['boundary'] is block boundary defined previously.
                x_chosen = x
                break
    b = np.array(b)
    weights = 1.0 / (gps_err ** 2)
    w_b = b * np.sqrt(weights)
    w_A = A * np.sqrt(weights)[:, None]
    return w_A, w_b
