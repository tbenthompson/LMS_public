from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import scipy.ndimage
import lms_code.lib.rep2 as rep2
from lms_code.lib.enjoy_dem import extract_line, get_config,\
    get_srtm_data, dem_interp_value

def plot_dem():
    # Load the DEM
    cfg = get_config('./data/dem/e100n40.Bathymetry.srtm')
    dem = np.load('./data/dem/e100n40.npy')

    # Determine the xsection normal and tangential direction
    geom = rep2.load('lms_geometry')
    xsec_normal = np.array(geom['xsec_plane'][0])[:2]
    xsec_tangent = np.array(geom['tibet_pt']) - np.array(geom['basin_pt'])
    xsec_tangent /= np.linalg.norm(xsec_tangent)

    # Extract the elevations right along the xsec
    elevations, _, elevation_latlon_list = extract_line(dem, cfg, geom['tibet_pt'], geom['basin_pt'], 100)

    n_tangential = len(elevation_latlon_list[0])
    n_normal = 101
    tangential_scale = 10.0
    normal_scale = 5.0
    pts = []
    z = np.empty((n_tangential, n_normal))
    tangential_dir = np.linspace(0.0, tangential_scale, n_tangential)
    normal_dir = np.linspace(0.0, normal_scale, n_normal)
    for i, t_dist in enumerate(tangential_dir):
        t_frac = t_dist / tangential_scale
        xsec_pt = [
            elevation_latlon_list[0][0]*(1-t_frac) + elevation_latlon_list[0][-1]*(t_frac),
            elevation_latlon_list[1][0]*(1-t_frac) + elevation_latlon_list[1][-1]*(t_frac)
        ]

        for j, n_dist in enumerate(normal_dir):
            pt_to_add = xsec_pt + n_dist * xsec_normal
            pts.append(pt_to_add)
            z[i, j] = dem_interp_value(dem, cfg, pt_to_add)
    plt_n_scale = 2.0
    plt_t_scale = 1.0
    T_grid, N_grid = np.meshgrid(plt_n_scale * normal_dir, plt_t_scale * tangential_dir)

    #Smooth the z coordinate
    z = scipy.ndimage.gaussian_filter(z, sigma=2, order=0)

    # Plot the 3D field
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.plot_surface(T_grid, N_grid, z, cmap = plt.cm.Spectral)
    return ax

def main():
    # ax = plot_dem()
    interior = rep2.load('interior_all_details_coalesced')
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    x = interior['disc_x']
    y = interior['disc_y']
    tris = interior['disc_tris']['tibet']
    f = interior['interseis_ux']

    triang = mtri.Triangulation(x, y, tris)
    ax.plot_trisurf(triang, f, cmap = plt.cm.Spectral)
    # ax.tricontourf(triang, f, levels = np.linspace(-1.0, 7.5, 18))
    plt.show()


if __name__ == "__main__":
    main()


# fT = np.linspace(tangential_dir[0], tangential_dir[-1], n_tangential)
# fZ = np.linspace(-5000, 7000, 50)
# fTmat, fZmat = np.meshgrid(fT, fZ)
# f = (fTmat - 1) ** 2 + fZmat / 5000.0
# for i in range(n_tangential):
#     f[fZmat[:,i] > z[i,0],i] = np.nan
# cset = ax.contourf(f, fTmat, fZmat, zdir='x',
#                    antialiased = True, offset = -0.01, alpha = 1.0, extend = 'min')
#
ax.set_xlabel('X')
ax.set_xlim(tangential_dir[0], tangential_dir[-1])
ax.set_ylabel('Y')
ax.set_ylim(normal_dir[0], normal_dir[-1])

plt.show()
