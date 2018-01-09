import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
from scipy.special import sph_harm

def my_Znm(l, m, theta, phi):
    if m > 0:
        return (sph_harm(m, l, phi, theta) + np.conj(sph_harm(m, l, phi, theta)))/(np.sqrt(2))
    elif m == 0:
        return sph_harm(m, l, phi, theta)
    elif m < 0:
        return  (sph_harm(m, l, phi, theta) - np.conj(sph_harm(m, l, phi, theta)))/(np.sqrt(2)*1j)

def listmatrix2arraymatrix(arr, data_index='Hlm', row_index='i', col_index='j', lab_index='lm'):
    out_arr = np.zeros((max(arr[row_index])+1, max(arr[col_index])+1), dtype='complex')
    labels = np.zeros(max(arr[col_index]+1), dtype='U3')    
    for entry in arr:
        out_arr[entry[row_index], entry[col_index]] = entry[data_index]
        labels[entry[col_index]] = entry[lab_index]
    return out_arr, labels

def listmatrix2arraymatrix_real(arr, data_index='Hlm', row_index='i', col_index='j', lab_index='lm'):
    out_arr = np.zeros((max(arr[row_index])+1, max(arr[col_index])+1), dtype='complex')
    labels = np.zeros(max(arr[col_index]+1), dtype='U4')
    for entry in arr:
        out_arr[entry[row_index], entry[col_index]] = entry[data_index]
        labels[entry[col_index]] = str(entry[lab_index[0]]) + ',' + str(entry[lab_index[1]])
    return out_arr, labels


def plot_spherical(coeffs, labels, filename='spherical.png', ax=None, cax=None, r=1, mag=2):
    # Create a sphere
    theta, phi = np.mgrid[0:np.pi:101j, 0:2*np.pi:101j]
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
    mlab.clf()
    
    rc = 0*sph_harm(0, 0, phi, theta)
    for coeff, label in zip(coeffs, labels):
        l = int(label[0])
        m = int(label[2])
        if m == 0:
            rc += coeff*sph_harm(m, l, phi, theta)
        else:
            rc += (1/np.sqrt(2))*coeff*sph_harm(m, l, phi, theta)
            rc += (1/np.sqrt(2))*((-1)**m)*coeff.conjugate()*sph_harm(-m, l, phi, theta)

    # TODO Check that imaginary part is small
    if np.max(np.imag(rc)) > 1e-3:
        import pdb; pdb.set_trace() 
    rc = np.real(rc)
    #s /= s.max()
    mlab.mesh(rc*x, rc*y, rc*z, color=(1, 0, 0), representation='surface')

    mlab.view(azimuth=45, elevation=45, distance=5, focalpoint=None,
              roll=None, reset_roll=True, figure=None)
    mlab.savefig(filename, magnification=mag)
    # mlab.show()

def plot_spherical_real(coeffs, labels, filename='spherical.png', ax=None, cax=None, r=1, mag=2):
    # Create a sphere
    theta, phi = np.mgrid[0:np.pi:101j, 0:2*np.pi:101j]
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
    mlab.clf()
    
    rc = 0
    for coeff, label in zip(coeffs, labels):
        l = int(label.split(',')[0])
        m = int(label.split(',')[1])
        rc += coeff*my_Znm(l, m, theta, phi)
        
    # TODO Check that imaginary part is small
    # if np.max(np.imag(rc)) > 1e-3:
    #     import pdb; pdb.set_trace()
    rc = np.real(rc)
    
    n = rc.clip(max=0)
    p = rc.clip(min=0)*(-1)

    mlab.mesh(p*x, p*y, p*z, color=(1, 0, 0), representation='surface')
    mlab.mesh(n*x, n*y, n*z, color=(0, 0, 1), representation='surface')

    mlab.view(azimuth=45, elevation=45, distance=5, focalpoint=None,
              roll=None, reset_roll=True, figure=None)
    mlab.savefig(filename, magnification=mag)
    # mlab.show()
