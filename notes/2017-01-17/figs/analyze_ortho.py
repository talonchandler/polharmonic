import numpy as np
import dill
from polharmonic.util import *
from polharmonic.microscopes import *
np.set_printoptions(precision=2, suppress=True)

# Main script
print("Working...")

# Specify saving paths
folder = 'ortho'
file_name = 'psi_ortho.dat'
data_file = folder+'/'+file_name

# Calculate the system matrix
sys_matrix = ortho_sys_matrix(NA=0.8, n=1.33, phi_pols=[0, np.pi/4, np.pi/2, 3*np.pi/4])

# Save matrix_file
dill.dump(sys_matrix, open(data_file,'wb'))

# Load 
f = open(data_file, 'rb')
sys_matrix = dill.load(f)

# Generate plots
generate_plots(*sys_matrix, folder)
