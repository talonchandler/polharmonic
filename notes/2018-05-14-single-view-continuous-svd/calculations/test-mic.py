from polharmonic import det, ill, micro
import numpy as np

n_px=2**8

iii = ill.Illuminator(polarizer=False)
ddd = det.Detector(polarizer=True)
mmm = micro.Microscope(ill=iii, det=ddd)
# mmm.plot(mmm.H, filename='Hdet.pdf', n_px=n_px)
mmm.calc_SVD(n_px=n_px)
mmm.plot_SVS(filename='SVSdet.pdf')

iii = ill.Illuminator(polarizer=True)
ddd = det.Detector(polarizer=False)
mmm = micro.Microscope(ill=iii, det=ddd)
# mmm.plot(mmm.H, filename='Hill.pdf', n_px=n_px)
mmm.calc_SVD(n_px=n_px)
mmm.plot_SVS(filename='SVSill.pdf')
