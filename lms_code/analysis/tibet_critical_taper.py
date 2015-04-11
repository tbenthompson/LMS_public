import numpy as np
import matplotlib.pyplot as plt

qi_geom = False

g = 9.8
surf_angle = 0.0 #alpha

detach_depth_delta = 4.5 * 1000
detach_length = 150.0 * 1000
if qi_geom:
    detach_angle = np.deg2rad(7.0)#np.arctan(detach_depth_delta / detach_length)
else:
    detach_angle = np.arctan(detach_depth_delta / detach_length)
print detach_angle
density = 2700.0 #rho
basal_pres_ratio = 0.370 #lambda_b
fluid_pres_ratio = 0.370 #lambda
thickness = 20.0 * 1000 #H
internal_friction_set = np.linspace(0.4, 1.4, 4) # mu
internal_fric_angle_set = np.arctan(internal_friction_set) # mu = tan(phi)
basal_cohesion_set = [0.0e6]#[0.0e6, 5.0e6] #S_b
wedge_cohesion_set = [5.0e6]#[0.0e6, 20.0e6] #C
basal_friction = np.linspace(0.0, 0.3, 10000)
min_b_fric = 100
max_b_fric = -100
for basal_cohesion in basal_cohesion_set:
    for wedge_cohesion in wedge_cohesion_set:
        for internal_fric in internal_friction_set:
            internal_fric_angle = np.arctan(internal_fric)
            lhs = surf_angle + detach_angle
            rhs_numer = detach_angle + basal_friction * (1 - basal_pres_ratio) + (basal_cohesion / (density * g * thickness))
            rhs_denom = 1 + (2 * (1 - fluid_pres_ratio) * (np.sin(internal_fric_angle) / (1 - np.sin(internal_fric_angle)))) + \
                        (wedge_cohesion / (density * g * thickness))
            rhs = rhs_numer / rhs_denom
            residual = np.abs(rhs - lhs)
            optimal_fric = basal_friction[np.argmin(residual)]
            print("Basal cohesion(Pa): %.2f \n"%basal_cohesion +\
                  "    Wedge cohesion: %.2f \n"%wedge_cohesion +\
                  "        Internal Friction: %.2f \n"%internal_fric +\
                  "        ---> Optimal Friction constant:" + str(optimal_fric))
            min_b_fric = min(min_b_fric, optimal_fric)
            max_b_fric = max(max_b_fric, optimal_fric)
print "Minimum basal friction: %.3f"%min_b_fric
print "Maximum basal friction: %.3f"%max_b_fric
