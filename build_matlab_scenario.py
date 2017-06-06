import numpy as np 
import pandas as pd 
from scipy.spatial.distance import pdist, squareform
import scipy.io as sio

# time we care about
t0 = pd.to_datetime('2016-01-21 00:00:00')
t1 = pd.to_datetime('2016-01-22 00:00:00')

# read tables
print('reading tables')
#av_veh = pd.read_csv('ignored_assets/av_veh.csv', parse_dates=[0]).set_index('time_bucket')
paxin_table = pd.read_csv('ignored_assets/paxin_table.csv', parse_dates=[0]).set_index('time_bucket')
paxout_table = pd.read_csv('ignored_assets/paxout_table.csv', parse_dates=[0]).set_index('time_bucket')
#rebin_table = pd.read_csv('ignored_assets/rebin_table.csv', parse_dates=[0]).set_index('time_bucket')
#rebout_table = pd.read_csv('ignored_assets/rebout_table.csv', parse_dates=[0]).set_index('time_bucket')
#starts_table = pd.read_csv('ignored_assets/starts_table.csv', parse_dates=[0]).set_index('time_bucket')
#breaks_table = pd.read_csv('ignored_assets/breaks_table.csv', parse_dates=[0]).set_index('time_bucket')

# distance stuff
print('building distance matrix')
station_locations = pd.read_csv('inferred_locations.csv')
D = squareform(pdist(station_locations[['x','y']].as_matrix()))
stations = list(station_locations['start_district_hash'])

# parameters
T = 24 * 12
N = len(stations)

# road graph
print('building road graph')
adjacencies = np.zeros((N,), dtype=np.object)

for i,s in enumerate(stations):
    adjacencies[i] = range(1,len(stations) + 1) # matlab friendly

D2 = np.ceil(D / 5.) # discretize to 5 minute travel times

np.fill_diagonal(D2, 1.) # if you travel within a zone, you must spend one time period

flowsout = paxout_table[stations]
flowsin = paxin_table[stations]

problem  ={
        'RoadGraph':adjacencies,     
        'TravelTimes': D2,
        'T':float(T),
        'FlowsIn': flowsin[t0:t1].as_matrix(),
        'FlowsOut': flowsout[t0:t1].as_matrix()
    }

print('saving mat')
sio.savemat('didi.mat',problem)
print('done!')