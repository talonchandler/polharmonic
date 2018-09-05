from polharmonic import ill, det, micro, util, shcoeffs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import os
import subprocess

class MultiMicroscope:
    """A MultiMicroscope represents an experiment that collects intensity data 
    under several different conditions (different polarization states or 
    illumination schemes).

    A MultiMicroscope mainly consists of a list of Microscopes.
    """
    def __init__(self, ill_optical_axes=[[0,0,1], [1,0,0]],
                 det_optical_axes=[[1,0,0], [0,0,1]],
                 ill_nas=2*[0], det_nas=2*[0.8], n_samp=1.33, sigma_ax=0.33):

        m = [] # List of microscopes

        # Cycle through paths
        for i, det_optical_axis in enumerate(det_optical_axes):
            ill_ = ill.Illuminator(optical_axis=ill_optical_axes[i],
                                   na=ill_nas[i], n=n_samp)
            det_ = det.Detector(optical_axis=det_optical_axes[i],
                                na=det_nas[i], n=n_samp, sigma_ax=sigma_ax)
            m.append(micro.Microscope(ill=ill_, det=det_)) # Add microscope

        self.micros = m
        self.N = len(m)
        self.jmax = m[0].h(0, 0, 0).jmax
        self.nbranch = self.N*self.micros[0].n_len

    def calc_SVD(self, n_px=2**6):
        w = 2.0
        self.xcoords = np.linspace(-w, w, n_px),
        X, Y, Z = np.mgrid[-w:w:n_px*1j, -w:w:n_px*1j, -w:w:n_px*1j]
        
        # For each position calculate K and solve eigenequation
        sigma = np.zeros((n_px, n_px, n_px, self.nbranch))
        for index, x in np.ndenumerate(X):
            u, s, v = self.calc_point_SVD(x, Y[index], Z[index])
            sigma[index[0], index[1], index[2], :] = s

        self.sigma = sigma
        self.sigma_max = np.max(sigma)

    def calc_point_SVD(self, x, y, z):
        # TODO: Generalize to n views
        HH0 = self.micros[0].H(x, y, z).coeffs
        HH1 = self.micros[1].H(x, y, z).coeffs
        HH = np.concatenate((HH0, HH1))
            
        K = np.dot(HH, HH.T)
        mu, v = np.linalg.eigh(K)
        u = np.dot(HH.T, v)

        return u[:,::-1], np.sqrt(mu[::-1]), v[:,::-1] # Transpose

    def plot_SVS_3D(self, filename='svs.pdf', n_px=2**6, marks=np.array([[0,0,0], [0.5,0,0], [1,0,0], [1.5,0,0]])):
        print('Plotting: ' + filename)

        # Layout windows
        inches = 1.5
        rows = self.sigma.shape[-1] + 1
        cols = marks.shape[0] + 1
        f = plt.figure(figsize=(inches*cols, inches*(rows - 0.75)))
        widths = cols*[1]
        heights = [1]*(rows - 1) + [0.05]
        spec = gridspec.GridSpec(ncols=cols, nrows=rows, width_ratios=widths,
                                 height_ratios=heights)

        # For saving singular function pngs
        folder = 'singular'
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        # For each branch
        # self.sigma[0,0,0,:] = 10 # Marker for testing
        # self.sigma[-1,0,0,:] = 10 # Marker for testing
        for row in range(rows):
            for col in range(cols):
                if col == 0 and row != rows - 1:
                    mini_spec = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=spec[row, col])
                    for a in range(2):
                        for b in range(2):
                            # Annotate 
                            ax = f.add_subplot(mini_spec[a, b])
                            ax.axis('off')
                            sigma = self.sigma[:,:,:,row]/self.sigma_max
                            sigma_plot = np.zeros(sigma.shape[0:2])

                            if a == 0 and b == 1:
                                ax.annotate('', xy=(-0.1,-0.1), xytext=(-0.6, -0.1), xycoords='axes fraction', textcoords='axes fraction', va='center', arrowprops=dict(arrowstyle="<-", shrinkB=0, lw=0.5))
                                ax.annotate('', xy=(-0.1,-0.1), xytext=(-0.1, 0.4), xycoords='axes fraction', textcoords='axes fraction', ha='center', arrowprops=dict(arrowstyle="<-", shrinkB=0, lw=0.5))
                                ax.annotate('', xy=(-0.1,-0.1), xytext=(-0.1, -0.6), xycoords='axes fraction', textcoords='axes fraction', ha='center', arrowprops=dict(arrowstyle="<-", shrinkB=0, lw=0.5))
                                ax.annotate('', xy=(-0.1,-0.1), xytext=(0.4, -0.1), xycoords='axes fraction', textcoords='axes fraction', va='center', arrowprops=dict(arrowstyle="<-", shrinkB=0, lw=0.5))
                                ax.annotate('$z$', xy=(0,0), xytext=(-0.7, -0.1), xycoords='axes fraction', textcoords='axes fraction', va='center', ha='center', fontsize=8)
                                ax.annotate('$y$', xy=(0,0), xytext=(-0.1, 0.5), xycoords='axes fraction', textcoords='axes fraction', va='center', ha='center', fontsize=8)
                                ax.annotate('$z$', xy=(0,0), xytext=(-0.1, -0.7), xycoords='axes fraction', textcoords='axes fraction', va='center', ha='center', fontsize=8)
                                ax.annotate('$x$', xy=(0,0), xytext=(0.5, -0.1), xycoords='axes fraction', textcoords='axes fraction', va='center', ha='center', fontsize=8)
                                sigma_plot = np.max(sigma, axis=2).T
                                ax.plot(marks[:,0], marks[:,1], 'xk', ms=1.5, mew=0.2)
                            if a == 1 and b == 1:
                                sigma_plot = np.max(sigma, axis=1)[:,::-1].T
                                ax.plot(marks[:,0], marks[:,2], 'xk', ms=1.5, mew=0.2)
                            if a == 0 and b == 0:
                                sigma_plot = np.max(sigma, axis=0)[:,::-1]
                                ax.plot(marks[:,2], marks[:,1], 'xk', ms=1.5, mew=0.2)                                
                            ax.imshow(sigma_plot, cmap="bwr", vmin=-1, vmax=1, interpolation='none', extent=[-2,2,-2,2], origin='lower')
                            ax.set_xlim([-2.05,2.05])
                            ax.set_ylim([-2.05,2.05])            
                elif row != rows - 1:
                    ax = f.add_subplot(spec[row, col])
                    cl = col - 1
                    print(row, cl)                    
                    u, s, v = self.calc_point_SVD(marks[cl, 0], marks[cl, 1], marks[cl, 2])
                    if np.max(u[:,row]) > 0:
                        sh = shcoeffs.SHCoeffs(u[:,row])
                        shfilename = folder+'/'+str(row)+str(cl)+'.png'
                        sh.plot_dist(filename=shfilename, r=1.1)
                        from PIL import Image
                        im1 = np.asarray(Image.open(shfilename))
                        ax.imshow(im1, interpolation='none')
                        ax.annotate('$\sigma='+'{:.2f}'.format(s[row]/self.sigma_max)+'$', xy=(1, 1), xytext=(0.5, 0),
                                    textcoords='axes fraction', ha='center', va='center')

                    ax.axis('off')
                elif row == rows - 1 and col == 0:
                    ax = f.add_subplot(spec[row, col])
                    # Colorbars
                    X, Y = np.meshgrid(np.linspace(0, 1, 100),
                                       np.linspace(0, 1, 100))
                    ax.imshow(X, cmap="bwr", vmin=-1, vmax=1, interpolation='none', extent=[0,1,0,1], origin='lower', aspect='auto')
                    # ax.contour(X, levels, colors='k',linewidths=0.5, extent=[0,1,0,1], origin='lower',)
                    ax.set_xlim([0,1])
                    ax.set_ylim([0,1])
                    ax.tick_params(direction='out', bottom=True, top=False)
                    ax.xaxis.set_ticks([0, 0.5, 1.0])
                    ax.yaxis.set_ticks([])

        # f, axs = plt.subplots(rows, cols,
        #                       figsize=(inches*cols, inches*(rows - 0.75)),
        #                       gridspec_kw={'hspace':0.05, 'wspace':0.10,
        #                                    'height_ratios':[1]*(rows - 1) + [0.05]})
        

        # # Label top row with arrow
        # axs[0,0].annotate('', xy=(-0.03, 1.1), xytext=(0.55, 1.1),
        #                   textcoords='axes fraction', xycoords='axes fraction', ha='center', va='center',
        #                   arrowprops=dict(arrowstyle='<->, head_width=0.05, head_length=0.1',
        #                                   connectionstyle="arc3", linewidth=0.5),)
        # axs[0,0].annotate(r'$2\textrm{NA}/\lambda$', xy=(0,0), xytext=(0.25, 1.2),
        #                   textcoords='axes fraction', xycoords='axes fraction', ha='center', va='center', fontsize=7)

        # # Top row singular values        
        # for j, ax in enumerate(axs[:-1,0]):
        #     ax.axis('off')
        #     ax.annotate('$j='+str(j)+'$', xy=(1, 1), xytext=(-0.15, 0.5),
        #                 textcoords='axes fraction', ha='center', va='center', rotation=90)
        #     extent= [-2,2,-2,2]
        #     origin = 'lower'
        #     ax.imshow(self.sigma[:,:,j]/self.sigma_max, cmap="bwr", vmin=-1, vmax=1, interpolation='none', extent=extent, origin=origin)
        #     ax.set_xlim([-2.05,2.05])
        #     ax.set_ylim([-2.05,2.05])            
        #     levels = [-0.1, -1e-5, 1e-5, 0.1]
        #     ct = ax.contour(self.sigma[:,:,j]/self.sigma_max, levels, colors='k',linewidths=0.5, extent=extent, origin=origin)
        #     for mark in marks:
        #         ax.plot(mark[0]*np.cos(mark[1]), mark[0]*np.sin(mark[1]), 'xk', ms=2.5, mew=0.5)

        # # Colorbars
        # X, Y = np.meshgrid(np.linspace(0, 1, 100),
        #                      np.linspace(0, 1, 100))
        # axs[-1,0].imshow(X, cmap="bwr", vmin=-1, vmax=1, interpolation='none', extent=[0,1,0,1], origin='lower', aspect='auto')
        # axs[-1,0].contour(X, levels, colors='k',linewidths=0.5, extent=[0,1,0,1], origin='lower',)
        # axs[-1,0].set_xlim([0,1])
        # axs[-1,0].set_ylim([0,1])
        # axs[-1,0].tick_params(direction='out', bottom=True, top=False)
        # axs[-1,0].xaxis.set_ticks([0, 0.5, 1.0])
        # axs[-1,0].yaxis.set_ticks([])
        # for j, ax in enumerate(axs[-1, 1:]):
        #     ax.axis('off')

        # # For saving singular function pngs
        # folder = 'singular'
        # if not os.path.exists(folder):
        #     os.makedirs(folder)

        # # Object space singular functions
        # for j, ax in np.ndenumerate(axs[:-1,1:]):
        #     if self.det.optical_axis == [0,0,1]: # z-detection
        #         u, s, v = self.calc_point_SVD(marks[j[1]][0], marks[j[1]][1], 0)
        #     elif self.det.optical_axis == [1,0,0]: # x-detection
        #         u, s, v = self.calc_point_SVD(0, marks[j[1]][0], marks[j[1]][1])

        #     # Labels
        #     if j[0] == 0:
        #         rho_label = str(marks[j[1]][0])
        #         if marks[j[1]][0] < 1e-3:
        #             rho_label = '0'
        #         ax.annotate(r'$\rho = '+rho_label+'$', xy=(1, 1), xytext=(0.5, 1.15),
        #                     textcoords='axes fraction', ha='center', va='center')
        #         ax.annotate(r'$\phi_{\rho} = '+str(marks[j[1]][1])+'$', xy=(1, 1), xytext=(0.5, 0.95),
        #                     textcoords='axes fraction', ha='center', va='center')
        #     ax.axis('off')
        #     ax.annotate('$\sigma='+'{:.2f}'.format(s[j[0]]/self.sigma_max)+'$', xy=(1, 1), xytext=(0.5, 0),
        #                 textcoords='axes fraction', ha='center', va='center')

        #     # Create singular function plots
        #     sh = shcoeffs.SHCoeffs(u[:,j[0]])
        #     shfilename = folder+'/'+str(j[0])+str(j[1])+'.png'
        #     sh.plot_dist(filename=folder+'/'+str(j[0])+str(j[1])+'.png', r=1.1)

        #     from PIL import Image
        #     im1 = np.asarray(Image.open(shfilename))
        #     ax.imshow(im1, interpolation='none')
        
        f.savefig(filename, bbox_inches='tight')

    # def plot_SVS_3D(self, filename, mag=1, show=False):
    #     from mayavi import mlab
    #     print('Plotting: ' + filename)

    #     inches = 4
    #     f, axs = plt.subplots(1, self.nbranch,
    #                           figsize=(inches*self.N, inches),
    #                           gridspec_kw={'hspace':0.0, 'wspace':0.0})

    #     folder = 'singular3'
    #     if not os.path.exists(folder):
    #         os.makedirs(folder)
        
    #     for j, ax in enumerate(axs):
    #         mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
    #         mlab.clf()

    #         half = int(self.sigma.shape[0]/2) + 1
    #         scalar = self.sigma[:,:,:,j]/self.sigma_max
    #         src = mlab.pipeline.scalar_field(scalar)
    #         mlab.pipeline.image_plane_widget(src,
    #                                          plane_orientation='x_axes',
    #                                          slice_index=half,
    #                                          colormap='bwr',
    #                                          vmin=-1, vmax=1)
    #         mlab.pipeline.image_plane_widget(src,
    #                                          plane_orientation='y_axes',
    #                                          slice_index=half,
    #                                          colormap='bwr',
    #                                          vmin=-1, vmax=1)
    #         mlab.pipeline.image_plane_widget(src,
    #                                          plane_orientation='z_axes',
    #                                          slice_index=half,
    #                                          colormap='bwr',
    #                                          vmin=-1, vmax=1)
            
    #         scalar = self.sigma[:,:half,:,j]/self.sigma_max            
    #         s = mlab.contour3d(scalar, contours=[.9, .7, .5, .3, .1], colormap='bwr', vmax=1, vmin=-1, transparent=True)

    #         s.scene.light_manager.light_mode = "vtk"
    #         mlab.view(azimuth=45, elevation=45, distance=5*half, focalpoint=None,
    #                   roll=None, reset_roll=True, figure=None)

    #         ax.annotate('$j='+str(j)+'$', xy=(1, 1), xytext=(0.5, 1.1),
    #                     textcoords='axes fraction', ha='center', va='center')

    #         brfilename = folder+'/'+str(j)+'.png'
    #         mlab.savefig(brfilename, magnification=mag)
    #         subprocess.call(['convert', filename, '-transparent', 'white', filename])

    #         from PIL import Image
    #         ax.axis('off')
    #         im1 = np.asarray(Image.open(brfilename))
    #         ax.imshow(im1, interpolation='none')
    #         if show:
    #             mlab.show()

    #     f.savefig(filename, bbox_inches='tight')

    def plot_scene(self, filename):
        print('Plotting: ' + filename)
        scene_string = ''
        for m in self.micros:
            scene_string += m.scene_string()
        util.draw_scene(scene_string, filename=filename, save_file=True)
        
    def plot_frames(self, folder='out'):
        if not os.path.exists(folder):
            os.makedirs(folder)
        for i, m in enumerate(self.micros):
            m.plot_scene(filename=folder+'/scene'+str(i)+'.pdf')
            m.plot(m.H, filename=folder+'/otf'+str(i)+'.pdf')

