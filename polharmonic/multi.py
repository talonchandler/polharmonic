from polharmonic import ill, det, micro, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
                 ill_nas=2*[0], det_nas=2*[0.8], n_samp=1.33):

        m = [] # List of microscopes

        # Cycle through paths
        for i, det_optical_axis in enumerate(det_optical_axes):
            ill_ = ill.Illuminator(optical_axis=ill_optical_axes[i],
                                   na=ill_nas[i], n=n_samp)
            det_ = det.Detector(optical_axis=det_optical_axes[i],
                                na=det_nas[i], n=n_samp)
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

    def plot_SVS_3D(self, filename, mag=1, show=False):
        from mayavi import mlab
        print('Plotting: ' + filename)

        inches = 4
        f, axs = plt.subplots(1, self.nbranch,
                              figsize=(inches*self.N, inches),
                              gridspec_kw={'hspace':0.0, 'wspace':0.0})

        folder = 'singular3'
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        for j, ax in enumerate(axs):
            mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
            mlab.clf()

            half = int(self.sigma.shape[0]/2) + 1
            scalar = self.sigma[:,:,:,j]/self.sigma_max
            src = mlab.pipeline.scalar_field(scalar)
            mlab.pipeline.image_plane_widget(src,
                                             plane_orientation='x_axes',
                                             slice_index=half,
                                             colormap='bwr',
                                             vmin=-1, vmax=1)
            mlab.pipeline.image_plane_widget(src,
                                             plane_orientation='y_axes',
                                             slice_index=half,
                                             colormap='bwr',
                                             vmin=-1, vmax=1)
            mlab.pipeline.image_plane_widget(src,
                                             plane_orientation='z_axes',
                                             slice_index=half,
                                             colormap='bwr',
                                             vmin=-1, vmax=1)
            
            scalar = self.sigma[:,:half,:,j]/self.sigma_max            
            s = mlab.contour3d(scalar, contours=[.9, .7, .5, .3, .1], colormap='bwr', vmax=1, vmin=-1, transparent=True)

            s.scene.light_manager.light_mode = "vtk"
            mlab.view(azimuth=45, elevation=45, distance=5*half, focalpoint=None,
                      roll=None, reset_roll=True, figure=None)

            ax.annotate('$j='+str(j)+'$', xy=(1, 1), xytext=(0.5, 1.1),
                        textcoords='axes fraction', ha='center', va='center')

            brfilename = folder+'/'+str(j)+'.png'
            mlab.savefig(brfilename, magnification=mag)
            subprocess.call(['convert', filename, '-transparent', 'white', filename])

            from PIL import Image
            ax.axis('off')
            im1 = np.asarray(Image.open(brfilename))
            ax.imshow(im1, interpolation='none')
            if show:
                mlab.show()

        f.savefig(filename, bbox_inches='tight')

    def plot_SVS(self, filename='svs.pdf'):
        print('Plotting: ' + filename)

        # Layout windows
        inches = 1
        f, axs = plt.subplots(1, self.N,
                              figsize=(inches*self.N, inches),
                              gridspec_kw={'hspace':0.0, 'wspace':0.05})

        for j, ax in enumerate(axs):
            ax.axis('off')
            ax.annotate('$j='+str(j)+'$', xy=(1, 1), xytext=(0.5, 1.15),
                        textcoords='axes fraction', ha='center', va='center')
            
            ax.imshow(self.mu[:,:,j], cmap="bwr", vmin=-1, vmax=1, interpolation='none')
            levels = [-0.1, -1e-5, 1e-5, 0.1]
            ct = ax.contour(self.mu[:,:,j], levels, colors='k',linewidths=0.5)
            
        f.savefig(filename, bbox_inches='tight')

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

