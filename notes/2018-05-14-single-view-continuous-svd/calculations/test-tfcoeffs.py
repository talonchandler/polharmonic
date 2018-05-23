from polharmonic import tfcoeffs, shcoeffs
import numpy as np


x = [1,0,1,0]
y = [0,1,0,1]
z1 = tfcoeffs.TFCoeffs([x,y,y])
z2 = tfcoeffs.TFCoeffs([x,4*[0],4*[0]])
z1*z2
import pdb; pdb.set_trace() 

