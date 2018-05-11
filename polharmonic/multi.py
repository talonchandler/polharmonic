from polharmonic import ill, det, micro, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import os

class MultiMicroscope:
    """A MultiMicroscope represents an experiment that collects intensity data 
    under several different conditions (different polarization states or 
    illumination schemes).

    A MultiMicroscope mainly consists of a list of Microscopes.
    """
    def __init__(self, ill_optical_axes=[[0,0,1]], det_optical_axes=[[0,0,1]],
                 ill_nas=[0.8], det_nas=[0.8], n_samp=1.33,
                 ill_pols=[[1,0,0], [1/np.sqrt(2), 1/np.sqrt(2), 0], [0,1,0], [-1/np.sqrt(2), 1/np.sqrt(2), 0]],
                 det_pols=4*[None]):

        m = [] # List of microscopes

        # Cycle through paths
        for i, det_optical_axis in enumerate(det_optical_axes):
            # Cycle through polarizations
            for j, pol in enumerate(ill_pols):
                ill_ = ill.Illuminator(optical_axis=ill_optical_axes[i],
                                       na=ill_nas[i], n=n_samp,
                                       polarizer=ill_pols[j])
                det_ = det.Detector(optical_axis=det_optical_axes[i],
                                    na=det_nas[i], n=n_samp,
                                    polarizer=det_pols[j])
                m.append(micro.Microscope(ill=ill_, det=det_)) # Add microscope

        self.micros = m
        self.N = len(m)
        self.jmax = m[0].h(0).jmax

    def calc_SVD(self, n_px=2**6):
        w = 2.05
        [X, Y] = np.meshgrid(np.linspace(-w, w, n_px),
                             np.linspace(-w, w, n_px))
        R = np.sqrt(X**2 + Y**2)
        Phi = np.nan_to_num(np.arctan(Y/X))

        # For each position and frame calculate H
        H = np.zeros((n_px, n_px, self.jmax, self.N))
        for index, r in np.ndenumerate(R):
            for n, m in enumerate(self.micros):
                H[index[0], index[1], :, n] = (m.H(r, Phi[index])).coeffs

        # For each position calculate K and solve eigenequation
        mu = np.zeros((n_px, n_px, self.N))
        for index, r in np.ndenumerate(R):
            K = np.dot(H[index].T, H[index])
            w, v = np.linalg.eigh(K)
            mu[index[0], index[1], :] = w[::-1]

        self.mu = np.sqrt(mu/np.max(mu))

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

