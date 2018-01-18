from polharmonic import data, util, multi
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import dill
np.set_printoptions(precision=3, suppress=True)

# Import data
name_head = '/Users/Talon/GoogleDrive/projects/dispim-data/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif',
         'B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif']
input_files = np.array([name_head + name for name in names])

# Calibration data
cal = np.ones(8)
# cal = np.array([1.06, 1.0, 1.0, 1.03, 
#                 1.08, 1.05, 1.0, 1.04])
#cal[4:] = ((1.1/0.71)**2)*cal[4:]

# Redindexing
idx = [0, 1, 2, 3, 4, 5, 6, 7] #[2, 3, 0, 1, 6, 7, 4, 5] 
input_files = input_files[idx]
cal=cal[idx]

# Build dispim microscope
exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                            det_nas=[1.1, 0.8], max_l=4, n_pts=250)

# Load asymmetric dispim system matrix
folder = 'asym_dispim'
file_name = folder+'.dat'
data_file = folder+'/'+file_name

# Large volume (two actin strands)
xs, ys, zs = 50, 50, 12 
x0, y0, z0 = 180, 560, 150
d = 130

# Small volumes (one actin strand)
# xs, ys, zs = 10, 10, 10, 
# x0, y0, z0 = 240, 585, 175
# d = 50

# Load and plot intensity fields
intf = data.IntensityField(file_names=input_files, x0=x0, y0=y0, z0=z0, xspan=xs, yspan=ys, zspan=zs, cal=cal)
#intf.plot_int_field(output_file=folder+'/'+'data.pdf', shape=(2,4), line_start=(0, xs-1), line_end=(ys-1, 0), d=d, mag=1)

# Load system model
f = open(data_file, 'rb')
exp.psi = dill.load(f)

# Calculate B
exp.calc_B_matrix()

# Reconstruct phantom
df = exp.recon_dist_field(intf, mask_threshold=0.6, prior='single')
df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/data_recon_shuffle.png',
                    r=0.4, d=100, mag=1, show=True)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
