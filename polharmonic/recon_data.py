import util
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec
import h5py
import dill

# Import data
name_head = '/Users/Talon/GoogleDrive/projects/dispim-data/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif',
         'B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif']
input_files = [name_head + name for name in names]
width = 10
height = 10
data = np.zeros((width, height, 8))
for i, input_file in enumerate(input_files):
    im = util.tiff2array(input_file, x=215, y=585, z=175, width=width, height=height, slices=1)
    data[:,:,i] = im

# Load asymmetric dispim system matrix
folder = 'dispim_asymmetric'
file_name = 'psi_dispim_asymmetric.dat'
data_file = folder+'/'+file_name
f = open(data_file, 'rb')
sys_matrix = dill.load(f)
labels = sys_matrix[1]

psi = sys_matrix[0]
psi_pinv = np.linalg.pinv(psi)

import pdb; pdb.set_trace() 
# Perform reconstruction Psi^+*g = F (image of coefficients)
F = np.einsum('ij,klj->kli', psi_pinv, data)

# Normalize F
F = F/F.max()

util.plot_spherical_field(F, labels, filename='recon_test.png', ax=None, r=1, mag=1, show=True)
#util.plot_spherical(F[20,25,:], labels, filename='recon_test.png', ax=None, r=1, mag=1, show=True)



# Save result
h5f = h5py.File('data.h5', 'w')
h5f.create_dataset('result', data=result)
h5f.close()

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
