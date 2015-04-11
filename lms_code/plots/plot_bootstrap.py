import random
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
from lms_code.analysis.run_bem import get_slip_magnitude
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import sys

def main():
    model = 'all_details'
    filename = 'bootstrap_' + model
    # filename = 'bootstrap_alternate'
    shortenings = rep2.load(filename)
    est = rep2.load('shortening_estimate_' + model)
    bem_soln = rep2.load('bem_' + model)
    lms_plot.setup()
    plot_range = (0.0, 10.0, 21)
    fig, ax = plt.subplots()
    slip = False
    if slip:
        new_shortenings = []
        every_distance = 5.0
        replications = 75
        total_length = 0.0
        slip_sum = 0.0
        for e in bem_soln['fault_mesh']:
            total_length += e.length
            slip_mag = np.linalg.norm(get_slip_magnitude(e))
            slip_sum += slip_mag * e.length
            l = int(math.ceil(e.length / every_distance))
            for i in range(l):
                sample = random.sample(shortenings, replications)
                new_shortenings.extend([s * slip_mag for s in sample])
        slip = slip_sum / total_length
        print slip
        shortenings = np.array(new_shortenings)
        est['lsqr_shortening'] *= slip

    boot_mean = np.mean(shortenings)
    boot_var = np.sum((shortenings - boot_mean) ** 2) / shortenings.shape[0]
    boot_std = np.sqrt(boot_var)
    boot_median = np.median(shortenings)
    boot_mad = np.median(np.abs(shortenings - boot_median))
    boot_normal_scaled_mad = 1.4826 * boot_mad

    print "Shortening mean from bootstrap: " + str(boot_mean)
    print "Shortening std dev from bootstrap: " + str(boot_std)

    print "Shortening median from bootstrap: " + str(boot_median)
    print "Shortening MAD from bootstrap: " + str(boot_mad)
    print "Shortening normal-scaled MAD from bootstrap: " + str(boot_normal_scaled_mad)


    print "Binning data"
    xs = np.linspace(*plot_range)
    hist, bin_edges = np.histogram(shortenings, xs)
    cs = []
    cumulative = 0

    n = len(shortenings)
    one_sigma_width = 0.6827
    outside_one_sigma = (1.0 - one_sigma_width) / 2.0
    one_sigma_min = int(round(n * outside_one_sigma))
    one_sigma_max = int(round(n * (outside_one_sigma + one_sigma_width)))
    for bin in hist:
        if one_sigma_min < cumulative < one_sigma_max:
            cs.append('r')
        else:
            cs.append('#AAAAAA')
        cumulative += bin

    print "Plotting Data"
    ax.bar(xs[:-1], hist, color = cs, width = 0.35)
    if slip:
        ax.set_xlabel('$s$ (mm/yr)')
    else:
        ax.set_xlabel('$\delta v$ (mm/yr)')
    ax.set_ylabel('$N$')
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    if 'alternate' not in filename:
        # Least squares estimate
        text_left = 0.75
        top_text_pos = 313000
        ax.text(text_left, top_text_pos, r'$\delta v_{\mathrm{L}}$ = ' +
                            "%.1f"%est['lsqr_shortening'] +
                            r' $\pm$ ' +
                            "%.1f"%est['lsqr_shortening_error'] +
                            r' mm/yr', color = 'b')

        # Bootstrap median estimate
        ax.text(text_left, top_text_pos - 16000, r'$\delta v_{\mathrm{M}}$ = ' +
                            "%.1f"%boot_median +
                            r' $\pm$ ' +
                            "%.1f"%boot_normal_scaled_mad +
                            r' mm/yr', color = 'r')
        ax.plot([est['lsqr_shortening']], [340000], 'b.', markersize=16)
        ax.plot([est['lsqr_shortening'] - est['lsqr_shortening_error'],
                 est['lsqr_shortening'] + est['lsqr_shortening_error']],
                [340000, 340000], 'b-', linewidth=2)
    if slip:
        filename += '_slip'
    plt.tight_layout(pad = 0.1)
    plt.savefig(filename + "_plot")

if __name__ == "__main__":
    main()
