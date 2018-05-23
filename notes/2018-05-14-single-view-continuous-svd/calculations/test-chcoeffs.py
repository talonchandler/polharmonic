from polharmonic import gaunt, chcoeffs
import numpy as np

x = chcoeffs.CHCoeffs([0,1,0])
y = chcoeffs.CHCoeffs([0,1,0])
z = x*y

# z = tfcoeffs.TFCoeffs([x,y,y])
# z2 = z+z
import pdb; pdb.set_trace() 
