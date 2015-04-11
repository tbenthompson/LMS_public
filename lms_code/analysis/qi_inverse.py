import numpy as np
import matplotlib.pyplot as plt
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot

solns = rep2.load("bem_all_details_inverse1")
width = 0.2
qi_gps = rep2.load('qi_gps_near_' + str(width))
# plt.plot(solns['detachment']['x'][0, :], solns['detachment']['u_soln'][0, :])
# plt.plot(solns['beichuan']['x'][0, :], solns['beichuan']['u_soln'][0, :])
# plt.show()

gps_x = np.array(qi_gps['dist'])
gps_v = qi_gps['vel']
# plt.plot(gps_x, gps_v, '.')
# plt.plot(solns['beichuan']['x'][0, :], 5000*solns['beichuan']['u_soln'][0, :])
# plt.show()
det_v = []
bei_v = []

def get_closest_soln(gps_x_val, soln_xs, soln_vs):
    for sol_x, sol_v in zip(soln_xs, soln_vs):
        if abs(x - sol_x) < 250:
            return sol_v

for x in gps_x:
    det_v.append(get_closest_soln(x, solns['detachment']['x'][0, :],
                                  solns['detachment']['u_soln'][0, :]))
    bei_v.append(get_closest_soln(x, solns['beichuan']['x'][0, :],
                                  solns['beichuan']['u_soln'][0, :]))

det_v = np.array(det_v)
bei_v = np.array(bei_v)



n = (100, 100)
det_slips = np.linspace(0.0, 500.0, n[0])
bei_slips = np.linspace(0.0, 500.0, n[1])
det_slips_M, bei_slips_M = np.meshgrid(det_slips, bei_slips)

E = np.empty((n[1], n[0]))
for i in range(n[1]):
    for j in range(n[0]):
        v = det_slips_M[i, j] * det_v + bei_slips_M[i, j] * bei_v
        E[i, j] = np.sqrt(np.sum((v - gps_v) ** 2) / len(gps_v))

min_idx = np.unravel_index(np.argmin(E), E.shape)

lms_plot.setup()
cmap = lms_plot.get_summery_cmap()
levels = np.linspace(15, 120.0, 31)
cntrf = plt.contourf(det_slips_M, bei_slips_M, E, cmap = cmap, levels = levels)
plt.contour(det_slips_M, bei_slips_M, E, colors = 'k', linestyles = 'solid',
            levels = levels, linewidths = 0.9)
cmap = plt.colorbar(cntrf)
cmap.set_label('$E_{\\textrm{RMSE}}$')
plt.xlabel('$s_{\\textrm{d}}$')
plt.ylabel('$s_{\\textrm{b}}$')
plt.savefig('beichuan_detach_error.pdf')

min_bei = bei_slips_M[min_idx[0], min_idx[1]]
min_det = det_slips_M[min_idx[0], min_idx[1]]
plt.figure()
plot_params = lms_plot.params()
plt.plot(gps_x / 1000.0 - plot_params['view_left'], gps_v, '.', label = 'Data')
plt.plot(gps_x / 1000.0 - plot_params['view_left'], min_bei * bei_v + min_det * det_v,
         'o', label = 'Best fit model')
plt.xlim([0, plot_params['view_right'] - plot_params['view_left']])
plt.xlabel('$x$ (km)')
plt.ylabel('$v$ (cm)')
plt.legend(loc = 'lower left')
plt.savefig('beichuan_detach_error_best.pdf')
