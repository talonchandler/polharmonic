from polharmonic import det, ill, micro
import numpy as np

n_px=2**7 + 1
folder='out/'

i1 = ill.Illuminator(polarizer=False)
d1 = det.Detector(polarizer=True) # True
m1 = micro.Microscope(ill=i1, det=d1)
m1.plot(m1.h, filename=folder+'hhdet.pdf', n_px=n_px, contours=False)
# m1.plot(m1.H, filename=folder+'Hdet.pdf', n_px=n_px, contours=True)
# m1.calc_SVD(n_px=n_px)
# m1.plot_SVS(filename=folder+'SVSdet.pdf')


i2 = ill.Illuminator(polarizer=True)
d2 = det.Detector(polarizer=False)
m2 = micro.Microscope(ill=i2, det=d2)
m1.plot(m2.h, filename=folder+'hhill.pdf', n_px=n_px, contours=False)
# m2.plot(m2.H, filename=folder+'Hill.pdf', n_px=n_px)
# m2.calc_SVD(n_px=n_px)
# m2.plot_SVS(filename=folder+'SVSill.pdf')
