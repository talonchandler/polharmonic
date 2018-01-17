from polharmonic import data, util, multi
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import dill

def recon(folder, input_files, exp, shape, cal):
    # Load asymmetric dispim system matrix
    file_name = folder+'.dat'
    data_file = folder+'/'+file_name

    # Large (two strand)
    # xs, ys, zs = 50, 50, 12 
    # x0, y0, z0 = 190, 560, 150 

    # Small
    xs, ys, zs = 10, 10, 10, 
    x0, y0, z0 = 240, 585, 175

    intf = data.IntensityField(file_names=input_files, x0=x0, y0=y0, z0=z0, xspan=xs, yspan=ys, zspan=zs, cal=cal)
    intf.plot_int_field(output_file=folder+'/'+'data.pdf', shape=shape, line_start=(0, xs-1), line_end=(ys-1, 0))
    
    # Load 
    f = open(data_file, 'rb')
    exp.psi = dill.load(f)

    # Calculate B
    exp.calc_B_matrix()

    # Reconstruct phantom
    df = exp.recon_dist_field(intf, mask_threshold=0.8)
    df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/data_recon_shuffle.png',
                        r=0.1, d=100, mag=1, show=True)

    import pdb; pdb.set_trace() 
    print('Total time: '+str(np.round(time.time() - start, 2)))
    os.system('say "done"')

# Import data
name_head = '/Users/Talon/GoogleDrive/projects/dispim-data/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif',
         'B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif']
input_files = np.array([name_head + name for name in names])

# Calibration data
cal = np.array([1.06, 1.0, 1.0, 1.03, 
                1.08, 1.05, 1.0, 1.04])

#cal[4:] = ((1.1/0.71)**2)*cal[4:]

idx = [0, 1, 2, 3, 4, 5, 6, 7] #[2, 3, 0, 1, 6, 7, 4, 5]

# Build dispim microscope
exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                            det_nas=[1.1, 0.8], max_l=4, n_pts=250)

recon(folder='asym_dispim', input_files=input_files[idx], exp=exp, shape=(2,4), cal=cal[idx])



