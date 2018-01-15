import numpy as np
from polharmonic import util
import matplotlib.pyplot as plt
from mayavi import mlab

class IntensityField:
    """An IntensityField represents the data collected from Distribution Field.  
    It consists of an array of intensities g. 
    """
    def __init__(self, file_names, x0=0, y0=0, z0=0, xspan=10, yspan=10, zspan=10, cal=None):
        g = np.zeros((xspan, yspan, len(file_names)))
        for i, file_name in enumerate(file_names):
            im = util.tiff2array(file_name, x=x0, y=y0, z=z0, width=xspan, height=yspan, slices=zspan)
            if cal is None:
                g[:,:,i] = im
            else:
                g[:,:,i] = im/cal[i]
        g = g/g.max()
        self.g = g

    def plot_int_field(self, output_file='out.pdf', shape=(2, 4), row_labels=[], col_labels=[],
                       line_start=(0,0), line_end=(0,0),
                       roi_upper_left=(0, 0), roi_wh=(0,0), plot_roi=False, zslice=0, dim=2):
        # Create axes
        inches = 3
        fig, axs = plt.subplots(2*shape[0], shape[1], figsize=(3.5*inches, 3*inches), gridspec_kw={'height_ratios':[1, 0.2]*shape[0], 'hspace':0.1, 'wspace':0.1})
        main_axs = axs[::2,:]
        line_axs = axs[1::2,:]

        # Load data and plot
        for i, (main_ax, line_ax) in enumerate(zip(main_axs.flatten(), line_axs.flatten())):

            # Plot image
            im = self.g[:,:,i]
            # if dim == 3:
            #     x, y, z = np.ogrid[-5:5:64j, -5:5:64j, -5:5:64j]
            #     scalars = x * x * 0.5 + y * y + z * z * 2.0
            #     obj = mlab.contour3d(scalars, contours=8, transparent=True)
            #     import pdb; pdb.set_trace() 

            #     # sys.stdout.flush()            
            #     # sys.stdout.write("Plotting: "+ str(j) + '/' + str(N) + '\r')

            #     # mlab.triangular_mesh(p*xyz[:,0] + space*i[0] - shift, p*xyz[:,1] + space*i[1] - shift, p*xyz[:,2], triangles, color=(1, 0, 0))
            #     # mlab.triangular_mesh(n*xyz[:,0] + space*i[0] - shift, n*xyz[:,1] + space*i[1] - shift, n*xyz[:,2], triangles, color=(0, 0, 1))

            #     # View and save
            #     # mlab.view(azimuth=45, elevation=45, distance=d, focalpoint=None,
            #     #           roll=None, reset_roll=True, figure=None)
            #     # mlab.savefig(filename, magnification=mag)
            #     # subprocess.call(['convert', filename, '-transparent', 'white', filename])
            #     mlab.show()

            main_ax.imshow(im, interpolation=None, vmin=0, vmax=1)
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

