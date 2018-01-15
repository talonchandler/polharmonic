from polharmonic import data, util, multi
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import dill

def recon(folder, input_files, exp, shape, cal):
    # Load asymmetric dispim system matrix
    file_name = folder+'.dat'
    data_file = folder+'/'+file_name

    xspan = 15
    yspan = 15
    intf = data.IntensityField(input_files, x0=235, y0=585, z0=175, xspan=xspan, yspan=yspan, zspan=1, cal=cal)
    intf.plot_int_field(output_file=folder+'/'+'data.pdf', shape=shape, line_start=(0, xspan-1), line_end=(yspan-1, 0))
    
    # Load 
    f = open(data_file, 'rb')
    exp.psi = dill.load(f)

    # Calculate B
    exp.calc_B_matrix()

    # Reconstruct phantom
    df = exp.recon_dist_field(intf.g)
    df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/data_recon_shuffle.png',
                        r=0.2, d=100, mag=1, show=True)

    print('Total time: '+str(np.round(time.time() - start, 2)))
    os.system('say "done"')

## A-view only
# original
name_head = '/Users/Talon/GoogleDrive/projects/dispim-data/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif']
cal = np.array([1.06, 1.0, 1.0, 1.03])
idx = [2, 3, 0, 1]# Default [0, 1, 2, 3]
input_files = np.array([name_head + name for name in names])

# Build ortho microscope
exp = multi.MultiMicroscope(ill_thetas=[90], det_thetas=[0],
                            det_nas=[0.8], max_l=4, n_pts=200)
recon(folder='ortho', input_files=input_files[idx], exp=exp, shape=(1,4), cal=cal[idx])

## DISPIM
# Import data
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif',
         'B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif']
input_files = np.array([name_head + name for name in names])
# Calibration data
cal = np.array([1.06, 1.0, 1.0, 1.03,
                1.0, 1.04, 1.08, 1.05])

idx = [2, 3, 0, 1, 6, 7, 4, 5] # Default [0, 1, 2, 3, 4, 5, 6, 7]

# Build dispim microscope
exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                            det_nas=[0.8, 0.8], max_l=4, n_pts=500)

recon(folder='dispim', input_files=input_files[idx], exp=exp, shape=(2,4), cal=cal[idx])



