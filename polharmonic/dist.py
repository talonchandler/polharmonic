import numpy as np
from polharmonic import util, sft
from mayavi import mlab
import subprocess
from cvxopt import matrix, solvers
import sys

class DistributionField:
    """A DistributionField represents many fluorophore distributions. It 
    consists of an array of Distributions."""
    def __init__(self, sh_arr=None, f_arr=None):
        if f_arr is not None and sh_arr is not None:
            print("Warning: sh_arr and f_arr are redundant.")
        elif f_arr is None:
            self.sh_arr = sh_arr
            self.f_arr = None
        elif sh_arr is None:
            self.sh_arr = None
            self.f_arr = f_arr

    def calc_f_arr(self, B):
        self.f_arr = np.einsum('ij,klmj->klmi', B, self.sh_arr)

    def make_positive(self, B, max_l=None):
        for i in np.ndindex(self.sh_arr.shape[:2]):
            d = Distribution(self.sh_arr[i])
            d.make_positive(B, max_l=max_l)
            self.sh_arr[i] = d.sh
            
    def plot_dist_field(self, B, xyz, triangles, filename=None,
                        d=50, r=1, s=1, mag=1, show=False, mask=None):
        
        # Calculate radii
        radii = r*self.f_arr
        
        # Create figure
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 800))
        mlab.clf()

        N = np.sum(mask)
        j = 1

        # If single direction plot quiver, otherwise plot mesh
        # Black magic ahead for faster plotting
        nz_idx = np.nonzero(np.sum(radii, axis=3))
        nz = radii[nz_idx]
        if np.nonzero(nz)[0].shape[0] == nz.shape[0]:
            x, y, z, tp = np.nonzero(radii)
            scale = radii[x, y, z, tp]
            u = xyz[tp, 0]
            v = xyz[tp, 1]
            w = xyz[tp, 2]
            sss = np.random.random(*x.shape)
            mlab.quiver3d(x, y, z, u, v, w, scalars=scale, mode='cylinder',
                          scale_factor=r, scale_mode='scalar')
        else:
            # TODO replace for loop with single tri_mesh
            for i in range(nz.shape[0]):
                sys.stdout.flush()            
                sys.stdout.write("Plotting: "+ str(j) + '/' + str(N) + '\r')
                idx = [x[i] for x in nz_idx]
                r = nz[i,:]
                j += 1
                mlab.triangular_mesh(r*xyz[:,0] + s*idx[0],
                                     r*xyz[:,1] + s*idx[1],
                                     r*xyz[:,2] + s*idx[2],
                                     triangles, color=(1, 0, 0))

        if radii.shape[2] != 1: # If 3D dataset
            extent = [0, radii.shape[0], 0, radii.shape[1], 0, radii.shape[2]]
            mlab.outline(extent=extent, line_width=2*mag)
            mlab.points3d(0,0,0, color=(1, 1, 1))

        # View and save
        mlab.gcf().scene.parallel_projection = True
        mlab.view(azimuth=45, elevation=45, distance=d, focalpoint=None,
                  roll=None, reset_roll=True, figure=None)
        if show:
            mlab.show()
        mlab.savefig(filename, magnification=mag)

class Distribution:
    """A Distribution represents a fluorophore distribution. It has redundant
    representations in the angular domain (orientation distribution function)
    and the angular frequency domain (spherical harmonic coefficients).

    """
    def __init__(self, sh=None, f=None):
        if f is None:
            self.sh = sh
        if sh is None:
            self.f = f

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
