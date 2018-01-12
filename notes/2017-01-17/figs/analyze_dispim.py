import numpy as np
import dill
from util import *
from microscopes import *
np.set_printoptions(precision=2, suppress=True)

# Main script
print("Working...")

# Specify saving paths
folder = 'dispim'
file_name = 'psi_dispim.dat'
data_file = folder+'/'+file_name

# Calculate the system matrix
# sys_matrix = dispim_sys_matrix(NAx=0.7, NAz=1.1, n=1.33, phi_pols=[0, np.pi/4, np.pi/2, 3*np.pi/4])

# # Save matrix_file
# dill.dump(sys_matrix, open(data_file,'wb'))

# Load system matrix from file
f = open(data_file, 'rb')
sys_matrix = dill.load(f)

import pdb; pdb.set_trace() 
generate_plots(*sys_matrix, folder)
