from polharmonic import tfcoeffs, shcoeffs
import numpy as np


x = shcoeffs.SHCoeffs([0,1,0,0])
y = shcoeffs.SHCoeffs([0,1,0,1])
z = tfcoeffs.TFCoeffs([x,y,y])
z2 = z+z
import pdb; pdb.set_trace() 
