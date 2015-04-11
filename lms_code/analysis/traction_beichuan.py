import lms_code.lib.rep2 as rep2
import matplotlib.pyplot as plt

trac = rep2.load('bem_traction')
detach = rep2.load('bem_just_detach')
thrust = rep2.load('bem_just_thrust')

import ipdb; ipdb.set_trace()
plt.plot(trac['x'][0, :], trac['u_soln'][0, :])
plt.plot(detach['x'][0, :], detach['u_soln'][0, :])
plt.show()
