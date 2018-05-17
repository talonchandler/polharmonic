from polharmonic import det, ill, micro, gaunt, multi
import numpy as np

m = multi.MultiMicroscope(ill_optical_axes=[[0,0,1]], det_optical_axes=[[0,0,1]],
                          ill_nas=[0.8], det_nas=[0.8], n_samp=1.33,
                          ill_pols=[[1,0,0], [1/np.sqrt(2), 1/np.sqrt(2), 0], [0,1,0], [-1/np.sqrt(2), 1/np.sqrt(2), 0]],
                          det_pols=4*[None])
m.calc_SVD(n_px=2**7)
m.plot_SVS(filename='pol_illum/svs.pdf')
m.plot_scene(filename='pol_illum/scene.pdf')
# m.plot_frames(folder='pol_illum')

m = multi.MultiMicroscope(ill_optical_axes=[[0,0,1]], det_optical_axes=[[0,0,1]],
                          ill_nas=[0.8], det_nas=[0.8], n_samp=1.33,
                          ill_pols=4*[None],
                          det_pols=[[1,0,0], [1/np.sqrt(2), 1/np.sqrt(2), 0], [0,1,0], [-1/np.sqrt(2), 1/np.sqrt(2), 0]])
m.calc_SVD(n_px=2**7)
m.plot_SVS(filename='pol_detect/svs.pdf')
m.plot_scene(filename='pol_detect/scene.pdf')
# m.plot_frames(folder='pol_detect')

# 3-polarizer
m = multi.MultiMicroscope(ill_optical_axes=[[0,0,1]], det_optical_axes=[[0,0,1]],
                          ill_nas=[0.8], det_nas=[0.8], n_samp=1.33,
                          ill_pols=3*[None],
                          det_pols=[[1,0,0], [-0.5, np.sqrt(3)/2, 0], [-0.5, -np.sqrt(3)/2, 0]])
m.calc_SVD(n_px=2**7)
m.plot_SVS(filename='pol_detect3/svs.pdf')
m.plot_scene(filename='pol_detect3/scene.pdf')
m.plot_frames(folder='pol_detect3')

