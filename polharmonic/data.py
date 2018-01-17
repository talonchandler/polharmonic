import numpy as np
from polharmonic import util
import matplotlib.pyplot as plt
from mayavi import mlab

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
                
        
    def plot_int_field(self, output_file='out.pdf', shape=(2, 4), row_labels=[], col_labels=[],
                       line_start=(0,0), line_end=(0,0),
                       roi_upper_left=(0, 0), roi_wh=(0,0), plot_roi=False, zslice=0):
        # Create axes
        inches = 3
        fig, axs = plt.subplots(2*shape[0], shape[1], figsize=(3.5*inches, 3*inches), gridspec_kw={'height_ratios':[1, 0.2]*shape[0], 'hspace':0.1, 'wspace':0.1})
        main_axs = axs[::2,:]
        line_axs = axs[1::2,:]

        # Load data and plot
        for i, (main_ax, line_ax) in enumerate(zip(main_axs.flatten(), line_axs.flatten())):

            # Plot image
            im = self.g[:,:,:,i]
            mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1000, 1000))
            mlab.clf()
            mlab.points3d(0,0,0, color=(1, 1, 1))
            obj = mlab.contour3d(im, contours=[.2, .4, .6, .8, .9], transparent=True)
            
            mlab.outline(extent=[0,im.shape[0],0,im.shape[1],0,im.shape[2]])
            mlab.view(azimuth=45, elevation=45, distance=50, focalpoint=None,
                      roll=None, reset_roll=True, figure=None)

            # TODO save 3d figure
            mlab.show()

            main_ax.imshow(im[:,:,zslice], interpolation=None, vmin=0, vmax=1)
            main_ax.set_axis_off()
            line_ax.set_axis_off()

            # Plot roi box on image
            if not plot_roi:
                rx0, ry0 = roi_upper_left
                rx1, ry1 = rx0 + roi_wh[0], ry0 + roi_wh[1],
                main_ax.plot([rx0, rx1, rx1, rx0, rx0], [ry0, ry0, ry1, ry1, ry0], 'g-')

            # Plot line on image
            x0, y0 = line_start
            x1, y1 = line_end
            main_ax.plot([x0, x1], [y0, y1], 'r-')
            main_ax.plot([x0], [y0], 'r.', markersize=6)

            # Extract line profile
            num = 1000
            x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)
            zi = []
            for (xi, yi) in zip(x, y):
                zi.append(im[int(yi), int(xi)])

            # Plot profile
            line_ax.plot(zi, '-r')
            line_ax.plot(zi[0], '.r', markersize=6)
            line_ax.set_ylim(0, 1)
            line_ax.set_xlim(-0.05*len(zi), len(zi))

        # Label rows and columns
        for i, label in enumerate(col_labels):
            axs[0, i].annotate(label, xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
        for i, label in enumerate(row_labels):
            axs[2*i, 0].annotate(label, xy=(0,0), xytext=(-0.1, 0.5), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False, rotation=90)

        fig.savefig(output_file, dpi=200)

