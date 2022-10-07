import sys
import yt
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import time as tim
from mpi4py import MPI

from yt.units.unit_object import Unit
from yt.fields.derived_field import ValidateSpatial
from yt.data_objects.time_series import DatasetSeries
from yt.frontends.boxlib.api import AMReXDataset

class AMReXDatasetSeries(DatasetSeries): 
    _dataset_cls = AMReXDataset


# Get the current working directory
basedir = os.getcwd() 
# basedir = '/Users/mikemeehan/Research/PeleData/bplume_Ri4_Re200/'
os.chdir(basedir)

def main():

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "stixgeneral",
        # "font.sans-serif": ["Helvetica"],
        "font.size": 16})

    yt.funcs.mylog.setLevel(50) # eliminate output from yt load
    CPU_time = tim.time() # start timer
    yt.enable_parallelism() # init parallel processes
    comm = MPI.COMM_WORLD # start MPI processes

    if yt.is_root():
        print('---> Starting python script to extract slice data <---')

    # Set time limits for averaging
    starttime = 0.000000001 # start time of analysis
    endtime   = 99999 # end time of analysis
    nskip     = 1    # plt interval (useful for downsampling)

    ### Ray to query (small "pressure tap" region)
    ray_start = [0.60*8.5498, 0.0, 0.0]
    ray_end = [0.66*8.5498, 0.0, 0.0]

    variables = [('boxlib', 'pressure')]

    # Create the array of plt directories and time of plt file
    ts = AMReXDatasetSeries("plt*[0-9]", parallel=True)

    # Initialize storage array for all times and plt files
    time_storage = {}

    # Iterate through plt's to find time indices
    # this method leaves None's and index numbers for plt's not in time range
    for sto, ds in ts.piter(storage=time_storage):
        if 'old' in str(ds):
            continue
        if ds.current_time >= starttime and ds.current_time <= endtime:
            sto.result = float(ds.current_time)
            sto.result_id = str(ds)

    # Convert the storage dictionary values to time and plt file names
    time1   = np.array(list(time_storage.values()))
    time    = np.array([x for x in time1 if x is not None])
    numplt1 = np.array(list(time_storage.keys()))
    numplt  = np.array([x for x in numplt1 if x.startswith("plt")])

    # Sort these
    itime = np.argsort(time)
    numplt = numplt[itime]
    time = time[itime]

    # Print out plt number and exact starting and ending time
    if yt.is_root():
        print('start plt:  %s'    % numplt[0])
        print('end plt:    %s'    % numplt[-1])
        print('start time: %0.6f' % time[0])
        print('end time:   %0.6f' % time[-1])
        print('nPlots:     %i' % len(numplt))

    # Define number of time steps and a string to the time
    nt     = len(numplt)
    ttext  = 't%0.4f-%0.4f' % (time[0], time[-1])

    # Allocate arrays for variables
    pressure_storage = {}

    comm.Barrier()

    # reset amrex series
    ts = AMReXDatasetSeries(numplt, parallel=True)

    # Start loop
    if yt.is_root():
        print('---> Looping Through Time Series <---\n', flush=True)
    
    for sto, ds in ts.piter(storage=pressure_storage):

        count = np.argmin(abs(np.array(time)-float(ds.current_time)))
        print('budgets, time: %0.6f sec' % ds.current_time, flush=True)

        # Probe Point, Line, and Box
        dd = ds.ray(ray_start, ray_end)

        # Read velocities
        time = ds.current_time.value
        xs = dd[('boxlib', 'x')]
        ps = dd[('boxlib', 'pressure')]

        # Convert pressure to atm
        ps /= 1.013e+6

        # compute average pressure
        p_avg = np.mean(ps)

        # print('Pressure = %0.2f atm\n' % p_avg, flush=True)
        sto.result = float(p_avg)
        sto.result_id = float(time)

    # Print storage
    print(pressure_storage)

    # Convert to arrays
    ts = np.array(list(pressure_storage.keys()))
    ps = np.array(list(pressure_storage.values()))

    # Sort
    ps = ps[np.argsort(ts)]
    ts = ts[np.argsort(ts)]

    # Plot these arrays
    fig, ax = plt.subplots(1, 1)
    ax.plot(ts, ps, c='b', linewidth=1)
    ax.set_xlabel('t [s]')
    ax.set_ylabel('Pressure [atm]')
    ax.set_title('Pressure Time Series')
    ax.grid()

    plt.savefig('./pressure_timeseries.png', dpi=200)


if __name__ == "__main__": main()
