import lms_code.lib.rep2 as rep2
from lms_code.lib.least_squares import build_least_squares
import lms_code.plots.plot_all as lms_plot
import numpy as np
import matplotlib.pyplot as plt

def main(offset):
    geom = rep2.load('lms_geometry')
    # which_model = 'coarse'
    # which_model = 'just_thrust'
    # which_model = 'simple'
    which_model = 'all_details'
    bem_soln = rep2.load('bem_' + which_model)
    plot_params = lms_plot.params()

    geom['gps_parallel_vel'][0] += offset
    A, b = build_least_squares(geom['gps_dist'],
                               geom['gps_parallel_vel'],
                               geom['gps_parallel_sig'],
                               bem_soln['x'],
                               bem_soln['u_soln'],
                               plot_params['boundary'](bem_soln))

    soln = np.linalg.lstsq(A, b)[0]
    shortening = soln[0]
    tibet = soln[1]
    sichuan = soln[1] - soln[0]

    n = len(b)
    p = 2

    hat_matrix = A.dot(np.linalg.inv(A.T.dot(A))).dot(A.T)
    leverage = np.diag(hat_matrix)

    residual = b - A.dot(soln)
    std_dev_residuals = np.std(residual)


    I = np.identity(n)
    cov_e = (I - hat_matrix).dot((I - hat_matrix).T)
    std_dev_e = np.sqrt(np.diag(cov_e))
    student_resid = residual / std_dev_e

    cooks_d = ((student_resid ** 2) * leverage) / (p * (1 - leverage))

    cross_validation_statistic = (residual / (1 - leverage)) ** 2


    total_sum_sqs = np.sum((b - np.mean(b)) ** 2)
    residual_sum_sqs = np.sum((b - A.dot(soln)) ** 2)
    residual_variance = residual_sum_sqs / n
    FVU = (residual_sum_sqs / total_sum_sqs)
    R2 = 1 - FVU
    R2_bar = 1 - (1 - R2) * ((n - 1) / (n - p -1))
    model_error = np.sqrt(np.linalg.inv(A.T.dot(A)))

    print("Leverage: " + str(leverage))
    print("Residuals: " + str(residual))
    print("Residual_sigma: " + str(std_dev_residuals))
    print("Studentized residuals: " + str(student_resid))
    print("Cook's distance: " + str(cooks_d))
    print("Cross validation statistic: " + str(cross_validation_statistic))
    print("Shortening Rate: " + str(shortening))
    print("Tibet velocity: " + str(tibet))
    print("Sichuan velocity: " + str(sichuan))
    print("Shortening rate std dev: " + str(model_error[0, 0]))
    print("Tibet velocity std dev: " + str(model_error[1, 1]))
    print("Tibet-Shortening corr coeff: " + str(model_error[0, 1]))
    print("Residual Variance: " + str(residual_variance))
    print("FVU: " + str(FVU))
    print("R2: " + str(R2))
    print("R2_bar: " + str(R2_bar))

    d = dict()
    d['lsqr_shortening'] = shortening
    d['lsqr_tibet'] = tibet
    d['lsqr_shortening_error'] = model_error[0, 0]
    rep2.save('shortening_estimate_' + which_model, d)
    return shortening

def debug():
    print(A)
    print(A.dot(soln))
    print(b)
    plt.plot(b, 'r-')
    plt.plot(A.dot(soln), 'b-')
    plt.show()

    from plot_1 import plot_lms_gps, plot_lms_vels
    lms_plot.setup()
    fig, ax = plt.subplots()
    plot_lms_vels(bem_soln, plot_params, ax, shortening, tibet, sichuan, 0, label = "", plot_y = False)
    plot_lms_gps(geom, plot_params, ax)
    ax.set_xlabel("x (kilometers)")
    plt.show()

def sensitivity():
    off = np.linspace(-5, 5, 10)
    shortenings = []
    s2 = []
    for o in off:
        shortenings.append(main(o))
        s2.append(0.379 * o)
    plt.plot(off, shortenings)
    plt.plot(off, s2)
    plt.show()

if __name__ == '__main__':
    main(0.0)
