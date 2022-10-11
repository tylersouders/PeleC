import sys
#sys.path.append('/p/home/whitmans/yt-conda/lib/python3.7/site-packages')
import mpi4py
import yt
import numpy as np
import os
import glob
# import warnings

# Silence MPL deprecation warnings baked into yt
# warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# Silence output from yt loading
yt.funcs.mylog.setLevel(50) # eliminate output from yt load

def main():
    yt.enable_parallelism()
    if yt.is_root():
        print('Starting python script to make images and movies')

    # Get the current working directory
    basedir = os.getcwd()

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Some User Options
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    make_movie = True # Flag to enable/disable ffmpeg movie generation
    slice_txt_format = 'Slice*00000.png' # format of yt outputs (may change with a lot of images)

    # Times to include in image processing
    start_time = 0.0000
    end_time = 10000
    nskip = 1

    # Where to take slices
    inlet_slice = False
    if inlet_slice:
       slice_center=[0.0001, 0.0, 0.0]
       slice_width=((0.512,'cm'),(0.512,'cm'))
       slice_normal = 'x'
       dir_suffix = '_inlet'
    else:
       slice_center = [6.0, 0.0, 0.00]
       slice_width = ((12.0, 'cm'), (4.0, 'cm'))
       slice_normal = 'z'
       dir_suffix = '_side'

    ## Change region
    xlo = np.array([-3.4822, -0.256])
    xhi = np.array([8.5498, 0.256])
    xlen = xhi - xlo

    # slice_center = [0.0, 0.0, 0.01]
    slice_center = [xlo[0]+xlen[0]/2.0, xlo[1]+xlen[1]/2.0, 0.01]
    # slice_width = ((4.0, 'cm'), (1.6, 'cm'))
    slice_width = ((xlen[0], 'cm'), (xlen[1], 'cm'))
    slice_normal = 'z'

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Find index of the start and end times
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numplt, time = group_timeseries(start_time, end_time, nskip)

    # Define number of time steps and a string to the time
    nt     = len(numplt)
    dt     = (end_time-start_time)/(nt-1)
    ttext  = 't%0.4f-%0.4f' % (time[0], time[-1])
    imagedir = basedir + '/images_' + ttext + dir_suffix

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create images
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Create parallelizable time series, parallel=1 just means 1 ds at a time
    ts = yt.DatasetSeries(numplt, parallel=True)

    # Loop and generate images
    for ds in ts.piter():

        max_level = ds.index.max_level
        ref = int(np.product(ds.ref_factors[0:max_level]))
        dims = ds.domain_dimensions*ref
        nx, ny, nz = dims

        count = np.argmin(abs(np.array(time)-float(ds.current_time)))

        if yt.is_root():
           # print('time: %0.5f sec, ctime: %0.5f, count = %i, plt = %s' % \
           #     (ds.current_time, (ds.current_time-starttime)/dt, count, str(ds)))
           print('time: %0.5f sec' % ds.current_time)

        # Call plotting function (keeps things cleaner)
        drawplot(ds, slice_normal, slice_center, slice_width)

        # Change the names of all the files
        all_plt_imgs = glob.glob(str(ds) + '_Slice*.png')
        # print(all_imgs)
        for img in all_plt_imgs:
            if len(str(ds)) == 8:
                os.system('mv ' + img + ' ' + img[9:-4] + str(count).zfill(5) + '.png')
            elif len(str(ds)) == 9:
                os.system('mv ' + img + ' ' + img[10:-4] + str(count).zfill(5) + '.png')

    # Move images and clean up - make movies if desired
    process_plots(imagedir, slice_txt_format, make_movie)

def process_plots(imagedir, slicetextformat, do_movie):

    # Processing plots
    if yt.is_root():

        # Check for and create directory for images
        if not os.path.exists(imagedir):

            os.mkdir(imagedir)

        os.system('mv Slice_* ' + imagedir)
    
        # Change to new directory with all the images in it
        os.chdir(imagedir)
        os.system('rm -rf *.mov')
        all_imgs = glob.glob(slicetextformat)

        for count, img in enumerate(all_imgs):

            all_imgs[count] = img[0:-9]

        print(all_imgs)

        if do_movie:

            # Iterate throught all slices and save as movies
            for img in all_imgs:

                # print(img)
                os.system('ffmpeg -framerate 15 -i '+img+'%5d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+img+'.mov')                


def group_timeseries(starttime, endtime, nskip):

    if yt.is_root():

        print('Grouping plt files and constructing yt series')

    # Create the array of plt directories and time of plt file
    all_plt_files = glob.glob("plt?????") + glob.glob("plt??????")
    all_ts = yt.DatasetSeries(all_plt_files, parallel=True)

    # Initialize storage array for all times and plt files
    time_storage = {}

    # Iterate through plt's to find time indices
    # this method leaves None's and index numbers for plt's not in time range
    for sto, ds in all_ts.piter(storage=time_storage):
        if ds.current_time >= starttime and ds.current_time <= endtime:
            sto.result = float(ds.current_time)
            sto.result_id = str(ds)

    # Convert the storage dictionary values to time and plt file names
    time1   = np.array(list(time_storage.values()))
    time    = [x for x in time1 if x is not None]
    numplt1 = np.array(list(time_storage.keys()))
    numplt  = [x for x in numplt1 if x.startswith("plt")]

    # Sort these
    numplt.sort()
    time.sort()

    # Reduce these to number we are skipping
    numplt = numplt[0::nskip]
    time   = time[0::nskip]

    # Print out plt number and exact starting time
    startplt  = numplt[0]
    endplt    = numplt[-1]
    starttime = time[0]
    endtime   = time[-1]

    if yt.is_root():
        print('start plt:  %s'    % startplt)
        print('end plt:    %s'    % endplt)
        print('start time: %0.6f' % starttime)
        print('end time:   %0.6f' % endtime)    

    return numplt, time


def drawplot(ds, slice_normal, slice_center, slice_width):
# Put all the plotting settings here so that you don't wanna cry when you look at this script

    #### Velocity Plots
    sx = yt.SlicePlot(ds, slice_normal, 'x_velocity', origin="native", center=slice_center, width=slice_width)
    sx.set_cmap(field="x_velocity", cmap='bwr')
    sx.annotate_timestamp(corner='lower_left', draw_inset_box=True)
    sx.set_log('x_velocity',False)
    sx.set_zlim('x_velocity', -20000.0, 75000.0)
    sx.save()

    ## Mach
    sx = yt.SlicePlot(ds, slice_normal, 'MachNumber', origin="native", center=slice_center, width=slice_width)
    sx.set_cmap(field="MachNumber", cmap='bwr')
    sx.set_log('MachNumber', False)
    sx.set_zlim('MachNumber', 0.0, 2.0)
    sx.save()

    # sx = yt.SlicePlot(ds, slice_normal, 'y_velocity', origin="native", center=slice_center, width=slice_width)
    # sx.set_cmap(field="y_velocity", cmap='bwr')
    # sx.annotate_timestamp(corner='lower_left', draw_inset_box=True)
    # sx.set_log('y_velocity',False)
    # sx.set_zlim('y_velocity', -1350.0, 1350.0)
    # sx.save()

    # sx = yt.SlicePlot(ds, slice_normal, 'z_velocity', origin="native", center=slice_center, width=slice_width)
    # sx.set_cmap(field="z_velocity", cmap='bwr')
    # sx.annotate_timestamp(corner='lower_left', draw_inset_box=True)
    # sx.set_log('z_velocity',False)
    # sx.set_zlim('z_velocity', -300.0, 300.0)
    # sx.save()

    ### Vorticity Plot
    sx = yt.SlicePlot(ds,slice_normal,'magvort',origin="native", center=slice_center, width=slice_width) 
    sx.set_cmap(field="magvort", cmap='Hue Sat Lightness 2')
    sx.set_log('magvort',False)
    # sx.annotate_timestamp(corner='lower_left',draw_inset_box=True)
    # sx.annotate_grids(min_level=3)
    plot = sx.plots['magvort']
    sx.set_zlim('magvort', 5e1, 5e4)
    cbar  = plot.cb
    sx._setup_plots()
    # cbar.set_ticks([1e0,1e2,1e3,1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,1e5])
    # cbar.set_ticklabels(['$10^0$','$10^2$','$10^3$','$10^{4}$'])
    #cbar.set_label('$(\omega\cdot\omega)^{1/2}$ $(\mathrm{1/s})$')
    sx.save()

    ### Pressure Plot
    sx = yt.SlicePlot(ds, slice_normal, ("boxlib", "pressure"), origin="native", center=slice_center, width=slice_width)
    sx.set_cmap(field=("boxlib", "pressure"), cmap='kelp')
    # sx.annotate_timestamp(corner='lower_left', draw_inset_box=True)
    sx.set_log(("boxlib", "pressure"), False)
    sx.set_zlim('pressure', 1000000.0, 20000000.0)
    sx.save()

    ### Divergence Plot
    # sx = yt.SlicePlot(ds, slice_normal, 'divu', origin="native", center=slice_center, width=slice_width)
    # sx.set_cmap(field="divu", cmap='bwr')
    # sx.annotate_timestamp(corner='lower_left', draw_inset_box=True)
    # sx.set_log('divu',False)
    # sx.set_zlim('divu', -5000.0, 5000.0)
    # sx.save()

    ### Density Plot
    sx = yt.SlicePlot(ds, slice_normal, ("gas", "density"), origin="native", center=slice_center, width=slice_width)
    sx.set_cmap(field=("gas", "density"), cmap='jet')
    # sx.annotate_timestamp(corner='lower_left', draw_inset_box=True)
    sx.set_log(("gas", "density"),False)
    sx.set_zlim(("gas", "density"), 0.001, 0.03)
    sx.save()

    ### Temperature plot
    sx = yt.SlicePlot(ds, 
    slice_normal, 
    ('boxlib', 'Temp'), 
    origin="native", 
    center=slice_center, 
    width=slice_width)
    sx.set_cmap(field='Temp', cmap='hot')
    sx.set_log('Temp', False)
    sx.set_zlim('Temp', 100, 500)
    sx.save()

if __name__ == "__main__":

    main()
