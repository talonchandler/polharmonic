from polharmonic import det, ill, micro, multi
import numpy as np

#n_px=2**4 + 1
n_px=2**4 + 1
folder='multi-out/'

mm = multi.MultiMicroscope(sigma_ax=0.33)
# mm.micros[0].plot(mm.micros[0].H, filename=folder+'H0.pdf', n_px=n_px, plot_m=[-2, -1, 0, 1, 2])
mm.calc_SVD(n_px=n_px)
mm.plot_SVS_3D(filename=folder+'SVS3Dx.pdf')
mm.plot_SVS_3D(filename=folder+'SVS3Dy.pdf', marks=np.array([[0,0,0], [0, 0.5 ,0], [0,1,0], [0,1.5,0]]))
mm.plot_SVS_3D(filename=folder+'SVS3Dz.pdf', marks=np.array([[0,0,0], [0, 0, 0.5], [0,0,1], [0,0,1.5]]))


