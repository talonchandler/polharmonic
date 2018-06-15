from polharmonic import det, ill, micro, shcoeffs

a = shcoeffs.SHCoeffs([0,0,1,0,0,0])
a.plot_dist('normal.png', force_positive=False)
print(a)

b = a.rotate()
b.plot_dist('rotated.png', force_positive=False)
print(b)
