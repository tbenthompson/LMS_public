from scipy.io import loadmat



# which_profile = "LMS"
# which_profile = "longer_profile2"
# which_profile = "properly_aligned3"
def get_stations(which_profile):
    all_gps = loadmat('data/gps/CalaisGanSocquetBanerjeeApel_in_Calais_RF_trim_Clean_Tog.sta.data.mat')
    all_station_names = all_gps['station_name']

    prof_gps = loadmat('data/gps/' + which_profile + '.profile.mat')
    names = prof_gps['profile_station_names']
    stations = []
    for (idx, n) in enumerate(names):
        s = dict()
        s['name'] = n
        s['distance'] = prof_gps['parallel_distance'][idx][0]
        s['parallel_vel'] = prof_gps['parallel_vel'][idx][0]
        s['normal_vel'] = prof_gps['perpendicular_vel'][idx][0]
        s['parallel_sigma'] = prof_gps['par_sig'][idx][0]
        s['normal_sigma'] = prof_gps['per_sig'][idx][0]
        for (all_idx, all_name) in enumerate(all_station_names):
            if all_name == n:
                s['lat'] = all_gps['station_lat'][0][all_idx]
                s['lon'] = all_gps['station_lon'][0][all_idx]
                s['east_vel'] = all_gps['station_east_vel'][0][all_idx]
                s['north_vel'] = all_gps['station_north_vel'][0][all_idx]
        stations.append(s)
    return stations

if __name__ == "__main__":
    print get_stations(properly_aligned3)
