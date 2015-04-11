import numpy as np
import matplotlib.pyplot as plt
import lms_code.lib.rep2 as rep2
import lms_code.plots.plot_all as lms_plot
from lms_code.analysis.eval_interior import interior_fname

def calc_stress(poisson, sxx, sxy, syy):
    # In plane strain:
    # e_zz = 0 = (1/E) * (s_zz - mu * (s_yy + s_xx))
    # solve for s_zz
    szz = poisson * (sxx + syy)
    stress = [[sxx, sxy, 0], [sxy, syy, 0], [0, 0, szz]]
    return stress

def calc_strains(youngs, poisson, stress):
    diag = stress[0][0] + stress[1][1] + stress[2][2]
    delta = np.diag([1,1,1])
    strain = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            strain[i][j] = (1 / youngs) * (stress[i][j] - poisson * (diag * delta[i][j] - stress[i][j]))
    return strain

def calc_energy(stress, strain):
    energy = 0
    for i in range(3):
        for j in range(3):
            energy += stress[i][j] * strain[i][j]
    return energy

def fix_tris(tris, pts):
    new_pts = []
    new_tris = []
    mapper = dict()
    for t in tris:
        next = len(new_pts)
        new_tris.append([next, next + 1, next + 2])
        for p in t:
            mapper[p] = len(new_pts)
            new_pts.append(pts[p])
    return np.array(new_pts), np.array(new_tris), mapper

def main():
    model = 'all_details'
    info = lms_plot.interior_params()

    bem_soln = rep2.load('bem_' + model)
    shear = bem_soln['shear_modulus']
    poisson = bem_soln['poisson_ratio']
    youngs = 2 * shear * (1 + poisson)

    shortening_est = rep2.load('shortening_estimate_' + model)
    shortening = shortening_est['lsqr_shortening']
    block = dict()
    block['tibet'] = shortening_est['lsqr_tibet']
    block['sichuan'] = block['tibet'] - shortening

    int_eval = rep2.load('interior_mesh_' + model)

    pts = int_eval.meshpy_pts
    all_x = pts[:, 0]
    all_y = pts[:, 1]
    all_tris = []

    disc_map = dict()
    for (region_num, name) in info['regions']:
        disc_map[name] = dict()

    disc_x = []
    disc_x.extend(all_x)
    disc_x.extend(all_x)

    disc_y = []
    disc_y.extend(all_y)
    disc_y.extend(all_y)

    disc_tris = dict()
    disc_tris['tibet'] = []
    disc_tris['sichuan'] = []

    data = dict()
    data_disc = dict()
    data_disc['interseis_u'] = np.zeros((2, 2 * len(all_x)))
    first_type = True
    for type in info['types']:
        data[type] = np.zeros((2, len(all_x)))
        data_disc[type] = np.zeros((2, 2 * len(all_x)))
        for (region_num, name) in info['regions']:
            print name
            idx_offset = 0
            if name == 'tibet':
                idx_offset = len(all_x)

            for subset in range(info['subsets']):
                filename = interior_fname(model, name, subset, type)
                cur_file = rep2.load(filename)

                # Handle tris
                if first_type:
                    cur_tris = cur_file['tris']
                    all_tris.extend(cur_tris)
                    disc_tris[name].extend(cur_tris + idx_offset)

                # Handle pts
                for k,v in cur_file['soln'].iteritems():
                    data[type][:, k] = v
                    data_disc[type][:, k + idx_offset] = v

                # Handle interseismic velocity field
                if type == 'u':
                    for k,v in cur_file['soln'].iteritems():
                        data_disc['interseis_u'][:, k + idx_offset] = \
                            -shortening * v + block[name]

        first_type = False

    all = dict()
    all['x'] = all_x
    all['y'] = all_y
    all['disc_tris'] = disc_tris
    all['disc_x'] = np.array(disc_x)
    all['disc_y'] = np.array(disc_y)
    all['tris'] = all_tris

    def to_mks(d):
        return 0.001 * -shortening * d
    if 'u' in data_disc:
        all['ux'] = to_mks(data_disc['u'][0, :])
        all['uy'] = to_mks(data_disc['u'][1, :])

        all['interseis_ux'] = data_disc['interseis_u'][0, :]

    if 'sx' in data:
        all['sxx'] = to_mks(data['sx'][0, :])
        all['sxy'] = to_mks(data['sx'][1, :])
    if 'sy' in data:
        all['syy'] = to_mks(data['sy'][1, :])

    if 'sx' in data and 'sy' in data:
        stress = calc_stress(poisson, all['sxx'], all['sxy'], all['syy'])
        strain = calc_strains(youngs, poisson, stress)
        all['energy_density'] = calc_energy(stress, strain)
        all['exx'] = strain[0][0]
        all['exy'] = strain[0][1]
        all['eyy'] = strain[1][1]

    rep2.save('interior_' + model + '_coalesced', all)

if __name__ == "__main__":
    main()
