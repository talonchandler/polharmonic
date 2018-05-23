import numpy as np
from polharmonic import util, ill, det
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

class Microscope:
    """
    A Microscope represents an experiment that collects a single frame of 
    intensity data.  

    A Microscope is specified by its illumination path (an Illuminator object),
    and its detection path (a Detector object).
    """
    def __init__(self, ill=ill.Illuminator(), det=det.Detector(),
                 color=(0,1,0.3)):
        self.ill = ill
        self.det = det
        self.h0 = self.ill.h()*self.det.h(0)
        self.n_len = self.h0.coeffs.shape[0]
        self.j_len = self.h0.coeffs.shape[1]
        self.hnorm = (self.ill.h()*self.det.h(0)).coeffs[0, 0]
        self.Hnorm = (self.ill.H()*self.det.H(0)).coeffs[0, 0]
        self.color = color

    def h(self, r, phi=0):
        return self.ill.h()*self.det.h(r, phi=phi)/self.hnorm

    def H(self, nu, phi_nu=0):
        return self.ill.H()*self.det.H(nu, phi_nu=phi_nu)/self.Hnorm

    def plot(self, func=h, filename='micro.pdf', n_px=2**6, plot_m=[-2, 0, 2]):
        print('Plotting: ' + filename)
        
        # Calculate data for plotting
        w = 2.05
        [X, Y] = np.meshgrid(np.linspace(-w, w, n_px),
                             np.linspace(-w, w, n_px))
        R = np.sqrt(X**2 + Y**2)
        Phi = np.nan_to_num(np.arctan(Y/X))

        data = np.zeros((n_px, n_px, self.n_len, self.j_len))
        for index, r in np.ndenumerate(R):
            data[index[0], index[1], :, :] = (func(r, Phi[index])).coeffs
                
        # Layout windows
        if plot_m is None:
            mcols = self.h0.rmax
        else:
            mcols = len(plot_m)

        cols = mcols*self.n_len
        inches = 1
        f, axs = plt.subplots(self.h0.rmax, cols,
                              figsize=(inches*cols, inches*self.h0.rmax),
                              gridspec_kw={'hspace':0.0, 'wspace':0.05})
        
        for index, ax in np.ndenumerate(axs):
            l = 2*index[0]
            if plot_m is None:
                m = index[1] - int(self.h0.mmax/2)
            else:
                plot_n = [-2, 0, 2]
                n_ind = index[1]//3
                n = plot_n[n_ind]
                m = plot_m[index[1]%3]
            j = util.lm2j(l, m)
            ax.axis('off')
            
            if index[0] == 0:
                ax.annotate('$m='+str(m)+'$', xy=(1, 1), xytext=(0.5, 1.15),
                            textcoords='axes fraction', ha='center', va='center')
                if index[1] % 3 == 1:
                    ax.annotate('$n='+str(n)+'$', xy=(1, 1), xytext=(0.5, 1.4),
                                textcoords='axes fraction', ha='center', va='center')
                    
                
            if index[1] == 0:
                ax.annotate('$l='+str(l)+'$', xy=(1, 1), xytext=(-0.15, 0.5),
                            textcoords='axes fraction', ha='center', va='center',
                            rotation=90)
            
            if j is not None:
                ax.imshow(data[:,:,n,j], cmap="bwr", vmin=-1, vmax=1, interpolation='none')
                levels = [-0.1, -1e-5, 1e-5, 0.1]
                ct = ax.contour(data[:,:,n,j], levels, colors='k',linewidths=0.5)
            else:
                ax.imshow(np.zeros(data[:,:,0,0].shape), cmap="bwr", vmin=-1, vmax=1, interpolation='none')
                
        f.savefig(filename, bbox_inches='tight')

    def calc_SVD(self, n_px=2**6):
        w = 2.05
        self.xcoords = np.linspace(-w, w, n_px),
        [X, Y] = np.meshgrid(self.xcoords, self.xcoords)
        R = np.sqrt(X**2 + Y**2)
        Phi = np.nan_to_num(np.arctan(Y/X))

        # For each position and frame calculate K and solve eigenequation
        mu = np.zeros((n_px, n_px, self.n_len))
        for index, r in np.ndenumerate(R):
            HH = self.H(r, Phi[index]).coeffs
            K = np.dot(HH, HH.T)
            w, v = np.linalg.eigh(K)
            mu[index[0], index[1], :] = w[::-1]

        self.mu = np.sqrt(mu/np.max(mu))

    def plot_SVS(self, filename='svs.pdf', n_px=2**6):
        print('Plotting: ' + filename)

        # Layout windows
        inches = 1
        f, axs = plt.subplots(1, self.n_len + 1,
                              figsize=(inches*(self.n_len + 1), inches),
                              gridspec_kw={'hspace':0.0, 'wspace':0.05})

        for j, ax in enumerate(axs[1:]):
            ax.axis('off')
            ax.annotate('$j='+str(j)+'$', xy=(1, 1), xytext=(0.5, 1.15),
                        textcoords='axes fraction', ha='center', va='center')

            axs[0].plot(self.xcoords[0], self.mu[:,int(self.mu.shape[0]/2),j], 'k-', lw=0.5)
            ax.imshow(self.mu[:,:,j], cmap="bwr", vmin=-1, vmax=1, interpolation='none')
            levels = [-0.1, -1e-5, 1e-5, 0.1]
            ct = ax.contour(self.mu[:,:,j], levels, colors='k',linewidths=0.5)

        axs[0].set_xlim([-2.05, 2.05])
        axs[0].set_ylim([0, 1.05])
        f.savefig(filename, bbox_inches='tight')

    def scene_string(self):
        asy_string = ''
        ill_string = "circle(optical_axis, alpha, false, color);\n"
        ill_string = ill_string.replace('optical_axis', str(tuple(self.det.optical_axis)))
        ill_string = ill_string.replace('alpha', str(np.arcsin(self.ill.na/self.ill.n)))
        ill_string = ill_string.replace('color', str(self.color))
        asy_string += ill_string

        if self.ill.polarizer is not None:
            pol_string = "arrow(optical_axis, polarizer, color, false);\n"
            pol_string = pol_string.replace('optical_axis', str(tuple(self.ill.optical_axis)))
            pol_string = pol_string.replace('polarizer', str(tuple(self.ill.polarizer)))
            pol_string = pol_string.replace('color', str(self.color))
            asy_string += pol_string

        if self.det.polarizer is not None:
            pol_string = "arrow(optical_axis, polarizer, color, true);\n"
            pol_string = pol_string.replace('optical_axis', str(tuple(self.det.optical_axis)))
            pol_string = pol_string.replace('polarizer', str(tuple(self.det.polarizer)))
            pol_string = pol_string.replace('color', str(self.color))
            asy_string += pol_string
        
        det_string = "circle(optical_axis, alpha, true, color);\n"
        det_string = det_string.replace('optical_axis', str(tuple(self.det.optical_axis)))
        det_string = det_string.replace('alpha', str(np.arcsin(self.det.na/self.det.n)))
        det_string = det_string.replace('color', str(self.color))
        asy_string += det_string
        
        return asy_string

    def plot_scene(self, filename):
        print('Plotting: ' + filename)        
        util.draw_scene(self.scene_string(), filename=filename, save_file=True)

        
