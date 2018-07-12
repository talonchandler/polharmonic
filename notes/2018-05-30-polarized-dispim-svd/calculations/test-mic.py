from polharmonic import det, ill, micro
import numpy as np

#n_px=2**4 + 1
n_px=2**7 + 1
folder='out2/'

iz = ill.Illuminator(optical_axis=[0,0,1], na=0, polarizer=True, illuminate_all=False)
dx = det.Detector(optical_axis=[1,0,0], polarizer=False)
mx = micro.Microscope(ill=iz, det=dx)
# mx.plot(mx.h, filename=folder+'hhillx.pdf', n_px=n_px, contours=False,  plot_m=[-2, 1, 0, 1, 2])
mx.plot(mx.H, filename=folder+'Hillz.pdf', n_px=n_px, plot_m=[-2, -1, 0, 1, 2])
mx.calc_SVD(n_px=n_px)
mx.plot_SVS(filename=folder+'SVSz.pdf')

ix = ill.Illuminator(optical_axis=[1,0,0], na=0, polarizer=True, illuminate_all=False)
dz = det.Detector(optical_axis=[0,0,1], polarizer=False)
mz = micro.Microscope(ill=ix, det=dz)
# mz.plot(mz.h, filename=folder+'hhillx.pdf', n_px=n_px, contours=False,  plot_m=[-2, -1, 0, 1, 2])
mz.plot(mz.H, filename=folder+'Hillx.pdf', n_px=n_px, plot_m=[-2, -1, 0, 1, 2])
mz.calc_SVD(n_px=n_px)
mz.plot_SVS(filename=folder+'SVSx.pdf')
