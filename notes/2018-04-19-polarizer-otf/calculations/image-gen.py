import numpy as np
from util import *
import matplotlib.pyplot as plt
import os

def psf_gen(n=1.5, NA=1.1, M=20, lamb=0.588, pol=None,
            px_width=7.4, x_px=257, y_px=257, 
            out_path='img-output/'):

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    x_FOV = px_width*x_px/M
    y_FOV = px_width*y_px/M
    [X, Y] = np.meshgrid(np.linspace(-x_FOV/2, x_FOV/2, x_px),
                         np.linspace(-y_FOV/2, y_FOV/2, y_px))
    R = np.sqrt(X**2 + Y**2)
    Phi = np.nan_to_num(np.arctan(Y/X))

    funcs = [h00, h20, h22, H00, H20, H22]
    file_names = ['hh00an', 'hh20an', 'hh22an', 'H00an', 'H20an', 'H22an']
    for func, file_name in zip(funcs, file_names):
        file_path = out_path+file_name+'.tiff'
        print('Generating '+file_path)
        psf_arr = (func(R*NA/lamb, phi=Phi, n=n, NA=NA, phi_p=pol)).astype('float32')
        save_tiff(psf_arr, file_path)

psf_gen(n=1.5, NA=1.3, M=80, pol=None, out_path='img-output/NA13/')
psf_gen(n=1.5, NA=1.1, M=80, pol=None, out_path='img-output/NA11/')
psf_gen(n=1.5, NA=0.8, M=80, pol=None, out_path='img-output/NA08/')
psf_gen(n=1.5, NA=1.3, M=80, pol=0, out_path='img-output/NA13pol/')
psf_gen(n=1.5, NA=1.1, M=80, pol=0, out_path='img-output/NA11pol/')
psf_gen(n=1.5, NA=0.8, M=80, pol=0, out_path='img-output/NA08pol/')
        
