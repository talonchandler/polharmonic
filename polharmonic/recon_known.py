import util
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec
import h5py
import dill
import sft

# Generate truth orientations
nx, ny = 10, 10
t = np.linspace(0, np.pi/2, nx)
p = np.linspace(0, np.pi, ny)
tv, pv = np.meshgrid(t, p)
tvv = np.expand_dims(tv, 2)
pvv = np.expand_dims(pv, 2)
tp_field = np.concatenate([tvv, pvv], 2)

# Expand onto spherical harmonics
sph_field, labels = sft.field_tp_sft(tp_field, max_l=7)

# Save the input 
# util.plot_spherical_field(sph_field, labels, filename='tp_test.png', ax=None,
#                           r=0.25, mag=10, dist=50, show=False, gridn=75)

# Load asymmetric dispim system matrix
folder = 'dispim_asymmetric'
file_name = 'psi_dispim_asymmetric.dat'
data_file = folder+'/'+file_name
f = open(data_file, 'rb')
sys_matrix = dill.load(f)
labels_sys = sys_matrix[1]

psi = sys_matrix[0]
psi_pinv = np.linalg.pinv(psi)

# Truncate input to get rid of harmonics that don't pass
mask = np.in1d(labels, labels_sys)
labels_trunc = labels[mask]
sph_field_trunc = sph_field[:, :, mask]

# Apply forward model Psi*F = g
g = np.einsum('ij,klj->kli', psi, sph_field_trunc)

# Perform reconstruction Psi^+*g = F (image of coefficients)
F = np.einsum('ij,klj->kli', psi_pinv, g)

# Normalize F
F = F/F.max()

util.plot_spherical_field(F, labels_trunc, filename='recon_test.png', ax=None,
                          r=0.25, mag=10, dist=50, show=False, gridn=75)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
