from polharmonic import det, shcoeffs, gaunt
import numpy as np

det = det.Detector()

x = shcoeffs.SHCoeffs([1,1,0,1,0,0.5])
y = shcoeffs.SHCoeffs([1,0,0,1,0,0])
z = x*y
x.plot(ax='test.pdf')
y.plot(ax='test2.pdf')
z.plot(ax='test3.pdf')
print(x*y)
