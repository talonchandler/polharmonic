import numpy as np
from polharmonic import util
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mayavi import mlab
import matplotlib.image as mpimg
import subprocess

class IntensityField:
    """An IntensityField represents the data collected from a DistributionField.  
    It consists of an array of intensities g. 
    """
    def __init__(self, file_names=None, x0=0, y0=0, z0=0, xspan=10, yspan=10, zspan=10, cal=None):
        if file_names is not None:
            g = np.zeros((xspan, yspan, zspan, len(file_names)))
            for i, file_name in enumerate(file_names):
                im = util.tiff2array(file_name, x=x0, y=y0, z=z0, width=xspan, height=yspan, slices=zspan)
                if cal is None:
                    g[:,:,:,i] = im
                else:
                    g[:,:,:,i] = im/cal[i]
            g = g/g.max()
            self.g = g
        else:
            self.g = 0
                
        
    def plot_int_field(self, output_file='out.pdf', shape=(2, 4),
                       row_labels=['$x$-illumination\n $z$-detection\n 1.1 NA', '$z$-illumination\n $x$-detection\n 0.71 NA'], col_labels=['$\phi = 0^{\circ}$', '$\phi = 45^{\circ}$', '$\phi = 90^{\circ}$', '$\phi = 135^{\circ}$'],
                       line_start=(0,0), line_end=(0,0),
                       roi_upper_left=(0, 0), roi_wh=(0,0), plot_roi=False, zslice=0, d=50, mag=1):
        # Create axes
        inches = 3
        fig, axs = plt.subplots(shape[0], shape[1], figsize=(shape[1]*inches, shape[0]*inches), gridspec_kw={'height_ratios':[1]*shape[0], 'hspace':0, 'wspace':0})

        # Load data and plot
        for i, ax in enumerate(axs.flatten()):

            # Plot image
            im = self.g[:,:,:,i]
            mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 800))
            mlab.clf()
            mlab.points3d(0,0,0, color=(1, 1, 1))

            for i in 0.1*np.arange(9):
                color = cm.hot(i)
                obj = mlab.contour3d(im, contours=[i], color=color[:-1], opacity=0.8) 
            mlab.gcf().scene.parallel_projection = True
            mlab.outline(extent=[0,im.shape[0],0,im.shape[1],0,im.shape[2]], line_width=2*mag)
            mlab.view(azimuth=45, elevation=45, distance=d, focalpoint=None,
                      roll=None, reset_roll=True, figure=None)
            #mlab.show()
            mlab.savefig('temp.png', magnification=mag)
            subprocess.call(['convert', '-density', str(300), '-units', 'PixelsPerInch', 'temp.png', 'temp2.png'])
            im = mpimg.imread('temp2.png')
            subprocess.call(['rm', 'temp.png', 'temp2.png'])

            ax.imshow(im, interpolation=None, vmin=0, vmax=1)            
            ax.set_axis_off()

        # Label rows and columns
        for i, label in enumerate(col_labels):
            axs[0, i].annotate(label, xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
        for i, label in enumerate(row_labels):
            axs[i, 0].annotate(label, xy=(0,0), xytext=(-0.1, 0.5), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False, rotation=90)

        fig.savefig(output_file, dpi=400)
