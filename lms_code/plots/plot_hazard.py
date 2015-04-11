import matplotlib.pyplot as plt
import numpy as np
from lms_code.analysis.run_bem import get_slip_magnitude
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot

def main():
    lms_plot.setup()
    fig = plt.figure()

    which_model = 'all_details'
    bem_soln = rep2.load('bem_' + which_model)
    shortening = rep2.load('shortening_estimate_' + which_model)

    est = shortening['lsqr_shortening']
    est_low = est - shortening['lsqr_shortening_error']
    est_high = est + shortening['lsqr_shortening_error']

    total_length = 0.0
    slip = 0.0
    slip_low = 0.0
    slip_high = 0.0
    joint = [4.20012e5 + 1.6, -2.006e4 - 5]
    for e in bem_soln['fault_mesh']:
        if e.vertex1.loc[0] < joint[0] - 10:
            continue
        total_length += e.length
        slip_mag = np.linalg.norm(get_slip_magnitude(e))
        slip += e.length * est * slip_mag
        slip_low += e.length * est_low * slip_mag
        slip_high += e.length * est_high * slip_mag
    s = (slip / total_length) / 1000
    s_low = (slip_low / total_length) / 1000
    s_high = (slip_high / total_length) / 1000
    slip_err = s_high - s
    # s = 6.1 / 1000
    # s_low = 4.6 / 1000
    # s_high = 7.6 / 1000

    T = np.linspace(0, 3000, 100)
    d = T * s
    T_high = d / s_low
    T_low = d / s_high

    wenchuan_d = 4.0
    wenchuan_T_low = wenchuan_d / s_low
    wenchuan_T = wenchuan_d / s
    wenchuan_T_high = wenchuan_d / s_high
    print("Wenchuan recurrence: " + str(wenchuan_T) + " (low: " + str(wenchuan_T_low) + ", high: " + str(wenchuan_T_high) + ")")

    a_wells = 6.93
    b_wells = 0.82
    mag7_ad = np.exp((7.0 - a_wells) / b_wells)
    mag7_T = mag7_ad / s

    paleo_T = 2300
    paleo_ad = paleo_T * s
    paleo_mag = (np.log(paleo_ad) * b_wells) + a_wells

    plt.plot(d, T, 'k-')
    plt.fill_between(d, T_low, T_high, facecolor = '#AAAAAA')
    plt.plot([0, paleo_ad + 100], [paleo_T, paleo_T], 'k--')
    plt.plot([wenchuan_d, mag7_ad, paleo_ad], [wenchuan_T, mag7_T, paleo_T],
             linestyle = 'None',
             marker = 'o',
             markeredgewidth = 4.0,
             markeredgecolor = (0, 0, 0, 1.0),
             markerfacecolor = (1, 1, 1, 1.0),
             markersize = 15)

    # Plot Wenchuan
    text = 'Wenchuan-like $\\textrm{M}_{\\textrm{w}}$ 7.9 (' + '%.0f'%wenchuan_d + ' m, ' +\
            '%.0f'%wenchuan_T + ' years)'
    plt.annotate(text, (wenchuan_d, wenchuan_T),
                xytext = (wenchuan_d + 0.5, wenchuan_T - 50))

    # Plot the Mw 7 pt
    text = 'Typical $\\textrm{M}_{\\textrm{w}}$ 7.0 (' + '%.0f'%mag7_ad + ' m, ' +\
            '%.0f'%mag7_T + ' years)'
    plt.annotate(text, (mag7_ad, mag7_T),
                xytext = (mag7_ad + 0.9, mag7_T - 30))

    # Plot the paleoseismic pt
    text = 'Low paleoseismic estimate'
    plt.text(1.7, 2350, text)
    text =  '($Ran$ $et$ $al.$ 2010)'
    plt.text(1.7, 2200, text)
    text = '$\\textrm{M}_{\\textrm{w}}$ ' + '%0.f'%paleo_mag + ', ' + '%0.f'%paleo_ad + ' m'
    plt.annotate(text, (paleo_ad, paleo_T),
                xytext = (paleo_ad - 3.2, paleo_T + 30))

    plt.text(2.0, 40, '($Wells$ $and$ $Coppersmith$ 1994)')
    plt.text(0.5, 1800, 'average slip rate = ' + '%.1f'%(s * 1000) + ' $\pm$ %.1f'%(slip_err * 1000) + ' mm/yr')
    plt.ylabel('$T$ (years)')
    plt.xlabel('$d$ (meters)')
    plt.ylim([0, 2500])
    plt.xlim([0, 2500 * s])
    width = 7.0
    fig.set_size_inches([width, (6.0 / 8.0) * width])
    plt.savefig('hazard_' + which_model)

if __name__ == '__main__':
    main()

