from polharmonic import det, ill, micro, gaunt, multi
import numpy as np

m = multi.MultiMicroscope(ill_optical_axes=[[0,0,1]], det_optical_axes=[[0,0,1]],
                          ill_nas=[0.8], det_nas=[0.8], n_samp=1.33,
                          ill_pols=[[1,0,0], [1/np.sqrt(2), 1/np.sqrt(2), 0], [0,1,0], [-1/np.sqrt(2), 1/np.sqrt(2), 0]],
                          det_pols=4*[None])
m.plot_frames(folder='pol_illum')

m = multi.MultiMicroscope(ill_optical_axes=[[0,0,1]], det_optical_axes=[[0,0,1]],
                          ill_nas=[0.8], det_nas=[0.8], n_samp=1.33,
                          ill_pols=4*[None],
                          det_pols=[[1,0,0], [1/np.sqrt(2), 1/np.sqrt(2), 0], [0,1,0], [-1/np.sqrt(2), 1/np.sqrt(2), 0]])
m.plot_frames(folder='pol_detect')
