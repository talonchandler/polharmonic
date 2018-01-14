import numpy as np
from polharmonic import util, sft
from mayavi import mlab
import subprocess
from cvxopt import matrix, solvers

class DistributionField:
    """A DistributionField represents many fluorophore distributions. It 
    consists of an array of Distributions."""
    def __init__(self, sh_arr=None):
        self.sh_arr = sh_arr

    def make_positive(self, B):
        for i in np.ndindex(self.sh_arr.shape[:2]):
            d = Distribution(self.sh_arr[i])
            d.make_positive(B)
            self.sh_arr[i] = d.sh
            
    def plot_dist_field(self, B, xyz, triangles, filename=None, r=1, mag=1, show=False):
        # Calculate radii
        radii = np.einsum('ij,klj->kli', B, self.sh_arr)

        # Create figure
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
        mlab.clf()

        space = 2
        shift = space*np.max(radii.shape[:2])/2
        
        # Plot
        for i in np.ndindex(radii.shape[:2]):
            # Split into positive and negatives
            n = radii[i].clip(max=0) 
            p = radii[i].clip(min=0)*(-1)

            mlab.triangular_mesh(p*xyz[:,0] + space*i[0] - shift, p*xyz[:,1] + space*i[1] - shift, p*xyz[:,2], triangles, color=(1, 0, 0))
            mlab.triangular_mesh(n*xyz[:,0] + space*i[0] - shift, n*xyz[:,1] + space*i[1] - shift, n*xyz[:,2], triangles, color=(0, 0, 1))
        
        # View and save
        mlab.view(azimuth=45, elevation=45, distance=3, focalpoint=None,
                  roll=None, reset_roll=True, figure=None)
        mlab.savefig(filename, magnification=mag)
        subprocess.call(['convert', filename, '-transparent', 'white', filename])
        if show:
            mlab.show()

class Distribution:
    """A Distribution represents a fluorophore distribution. It has redundant
    representations in the angular domain (orientation distribution function)
    and the angular frequency domain (spherical harmonic coefficients).

    """
    def __init__(self, sh=None):
        self.sh = sh

    def make_positive(self, B):
        N = B.shape[1]
        M = B.shape[0]
        P = matrix(2*np.identity(N), tc='d')
        q = matrix(-2*self.sh, tc='d')
        G = matrix(-B, tc='d')
        h = matrix(np.zeros(M), tc='d')
        sol = solvers.qp(P, q, G, h)
        self.sh = np.array(sol['x']).flatten()
        
    def plot_dist(self, B, xyz, triangles, filename=None, r=1, mag=1, show=False):
        # Calculate radii
        radii = np.matmul(B, self.sh)

        # Split into positive and negatives
        n = radii.clip(max=0) 
        p = radii.clip(min=0)*(-1)

        # Create figure
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
        mlab.clf()
        
        # Plot
        mlab.triangular_mesh(p*xyz[:,0], p*xyz[:,1], p*xyz[:,2], triangles, color=(1, 0, 0))
        mlab.triangular_mesh(n*xyz[:,0], n*xyz[:,1], n*xyz[:,2], triangles, color=(0, 0, 1))
        
        # View and save
        mlab.view(azimuth=45, elevation=45, distance=3, focalpoint=None,
                  roll=None, reset_roll=True, figure=None)
        mlab.savefig(filename, magnification=mag)
        subprocess.call(['convert', filename, '-transparent', 'white', filename])
        if show:
            mlab.show()
