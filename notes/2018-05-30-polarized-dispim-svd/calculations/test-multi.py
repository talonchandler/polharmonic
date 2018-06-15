from polharmonic import det, ill, micro, multi
import numpy as np

#n_px=2**4 + 1
n_px=2**3 + 1
folder='multi-out/'

mm = multi.MultiMicroscope()
# mm.micros[0].plot(mm.micros[0].H, filename=folder+'H0.pdf', n_px=n_px, plot_m=[-2, -1, 0, 1, 2])
mm.calc_SVD(n_px=n_px)
mm.plot_SVS_3D(filename=folder+'SVS3D.pdf', show=False)

# mm.micros[0].calc_SVD(n_px=n_px)
# mm.micros[0].plot_SVS(filename=folder+'SVS0.pdf')
# mm.micros[1].calc_SVD(n_px=n_px)
# mm.micros[1].plot_SVS(filename=folder+'SVS1.pdf')

