import pyproj
import datetime
import numpy as np
from StringIO import StringIO
import lms_code.lib.rep2 as rep2
from lms_code.analysis.collect_data import get_plane
from lms_code.lib.geometry import project_onto_plane

def remove_comments(f):
    line = f.readline()
    if line.startswith('#'):
        print line
        return remove_comments(f)
    return line

with open('data/qi_gps/all', 'r') as f:
    first_line = remove_comments(f)
    text = first_line + f.read()

f = StringIO(text)
names = 'lon, lat, EW, SN, Sew, Sns, coff, Name,\
         Up, Sup, resurvey, who1, who2'
data = np.genfromtxt(f, names = names, dtype = None)


save_me = {k: {i:v for i,v in enumerate(data[k])} for k in data.dtype.names}
save_me['resurvey'] = {i:datetime.datetime.strptime(r, "%d-%b-%Y")
                       for i,r in enumerate(data['resurvey'])}


# project onto xsec
tibet_pt, basin_pt, xsec_plane = get_plane()
xsec_vector = np.array(basin_pt) - np.array(tibet_pt)
xsec_vector /= np.linalg.norm(xsec_vector)

g = pyproj.Geod(ellps = 'WGS84')
selected_stations = []
xsec_lon = []
xsec_lat = []
xsec_dist = []
xsec_vel = []
xsec_sig = []
width = 0.2
for i in range(len(data['lon'])):
    lon = data['lon'][i]
    lat = data['lat'][i]
    proj_pt, remainder = project_onto_plane((lat, lon, 0.0),
                                            xsec_plane[1], xsec_plane[0])
    distance = np.linalg.norm(remainder)
    if distance > width:
        continue
    selected_stations.append(i)
    xsec_lon.append(proj_pt[1])
    xsec_lat.append(proj_pt[0])
    az1, az2, dist = g.inv(tibet_pt[1], tibet_pt[0], proj_pt[1], proj_pt[0])
    xsec_dist.append(dist)
    par_vel = xsec_vector.dot([data['SN'][i], data['EW'][i]])
    par_sig = xsec_vector.dot([data['Sns'][i], data['Sew'][i]])
    xsec_vel.append(par_vel)
    xsec_sig.append(par_sig)
near_qi = dict()
near_qi['selected'] = selected_stations
near_qi['lon'] = xsec_lon
near_qi['lat'] = xsec_lat
near_qi['dist'] = xsec_dist
near_qi['vel'] = xsec_vel
near_qi['sig'] = xsec_sig

geom = rep2.load('lms_geometry')
rep2.save('qi_gps_imported', save_me)
rep2.save('qi_gps_near_' + str(width), near_qi)
