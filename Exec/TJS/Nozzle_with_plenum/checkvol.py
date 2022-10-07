import numpy as np
import yt

def main():

    file = 'plt00001'

    xlo = np.array([-3.4822, -0.256])
    xhi = np.array([8.5498, 0.256])
    xlen = xhi - xlo

    ds = yt.load(file)

    slice_center = [0.0, 0.0, 0.01]
    slice_center = [xlo[0]+xlen[0]/2.0, xlo[1]+xlen[1]/2.0, 0.01]
    slice_width = ((4.0, 'cm'), (1.6, 'cm'))
    slice_width = ((xlen[0], 'cm'), (xlen[1], 'cm'))
    slice_normal = 'z'

    sx = yt.SlicePlot(ds, 
        slice_normal, 
        ('boxlib', 'vfrac'),
        origin="native")
    sx.save()

    # sx = yt.SlicePlot(ds, 
    #     slice_normal, 
    #     ('boxlib', 'x_velocity'), 
    #     origin="native", 
    #     center=slice_center, 
    #     width=slice_width)
    # sx.set_cmap(field="x_velocity", cmap='bwr')
    # sx.set_log('x_velocity', False)
    # sx.set_zlim('x_velocity', 0.0, 30000.0)
    # sx.save()

    # sx = yt.SlicePlot(ds, 
    #     slice_normal, 
    #     ('boxlib', 'pressure'), 
    #     origin="native", 
    #     center=slice_center, 
    #     width=slice_width)
    # sx.set_cmap(field="pressure", cmap='hot')
    # sx.set_log('pressure', False)
    # sx.set_zlim('pressure', 1000000.0, 2500000.0)
    # sx.save()

    # sx = yt.SlicePlot(ds, 
    #     slice_normal, 
    #     ('boxlib', 'Temp'), 
    #     origin="native", 
    #     center=slice_center, 
    #     width=slice_width)
    # sx.set_cmap(field='Temp', cmap='hot')
    # sx.set_log('Temp', False)
    # sx.set_zlim('Temp', 250, 400)
    # sx.save()



if __name__ == "__main__": main()