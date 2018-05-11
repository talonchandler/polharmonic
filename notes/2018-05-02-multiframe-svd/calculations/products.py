from polharmonic import det, shcoeffs
import numpy as np

x = shcoeffs.SHCoeffs([1,0,0,1,0,0])
y = shcoeffs.SHCoeffs([1,1,0,1,0,1])
xy = x*y
x.plot(folder='mult/x')
y.plot(folder='mult/y')
xy.plot(folder='mult/xy')
