import numpy as np
from polharmonic import util, sft
from mayavi import mlab
import subprocess
from cvxopt import matrix, solvers
import sys

class DistributionField:
    """A DistributionField represents many fluorophore distributions. It 
    consists of an array of Distributions."""
    def __init__(self, sh_arr=None):
        self.sh_arr = sh_arr

    def make_positive(self, B, max_l=None):
        for i in np.ndindex(self.sh_arr.shape[:2]):
            d = Distribution(self.sh_arr[i])
            d.make_positive(B, max_l=max_l)
            self.sh_arr[i] = d.sh
            
    def plot_dist_field(self, B, xyz, triangles, filename=None,
                        d=50, r=1, mag=1, show=False, mask_threshold=0.01):
        # Calculate radii
        radii = r*np.einsum('ij,klmj->klmi', B, self.sh_arr)
        mask = np.max(radii, axis=-1) > mask_threshold
        
        # Create figure
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
        mlab.clf()

        space = 1

        N = np.sum(mask)
        j = 1
        
        # Plot
        for i in np.ndindex(radii.shape[:-1]):
            if mask[i]:
                sys.stdout.flush()            
                sys.stdout.write("Plotting: "+ str(j) + '/' + str(N) + '\r')

                # Split into positive and negatives
                n = radii[i].clip(max=0) 
                p = radii[i].clip(min=0)*(-1)
                j += 1

                if i == (2, 8, 0):
                    mlab.triangular_mesh(p*xyz[:,0] + space*i[0], p*xyz[:,1] + space*i[1], p*xyz[:,2] + space*i[2], triangles, color=(0, 1, 0))
                else:
                    mlab.triangular_mesh(p*xyz[:,0] + space*i[0], p*xyz[:,1] + space*i[1], p*xyz[:,2] + space*i[2], triangles, color=(1, 0, 0))

                if np.max(n) > mask_threshold:
                    mlab.triangular_mesh(n*xyz[:,0] + space*i[0], n*xyz[:,1] + space*i[1], n*xyz[:,2] + space*i[2], triangles, color=(0, 0, 1))

        if radii.shape[2] != 1:
            mlab.outline(extent=[0,radii.shape[0],0,radii.shape[1],0,radii.shape[2]])
            mlab.points3d(0,0,0, color=(1, 1, 1))
        
        # View and save
        mlab.view(azimuth=45, elevation=45, distance=d, focalpoint=None,
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

    def make_positive(self, B, max_l=None):
        N = B.shape[1]
        M = B.shape[0]
        P = matrix(2*np.identity(N), tc='d')
        q = matrix(-2*self.sh, tc='d')
        G = matrix(-B, tc='d')
        h = matrix(np.zeros(M), tc='d')
        if max_l is None:
            sol = solvers.qp(P, q, G, h)
        else:
            J = util.maxl2maxj(max_l)
            A = matrix(np.identity(N)[N-J:,:], tc='d')
            b = matrix(np.zeros(J), tc='d')
            sol = solvers.qp(P, q, G, h, A, b)
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
