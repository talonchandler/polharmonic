from polharmonic import data, util, multi
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import dill
import matplotlib.image as mpimg
np.set_printoptions(precision=3, suppress=True)

# Import data
name_head = '/Users/Talon/GoogleDrive/projects/dispim-data/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif',
         'B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif']
input_files = np.array([name_head + name for name in names])

# Calibration data
cal = np.array([1.06, 1.0, 1.0, 1.03, 
                1.08, 1.05, 1.0, 1.04])

# HACKING
# cal = np.array([1.06, 1.0, 1.0, 1.03, 
#                 0.97, 1.05, 1.0, 1.15])

cal_ratio = 1.1#((1.1/0.71)**2)
cal[4:] = cal_ratio*cal[4:]

# Redindexing
idx = [0, 3, 2, 1, 4, 5, 6, 7]
input_files = input_files[idx]
cal=cal[idx]

# Build dispim microscope
exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                            det_nas=[1.1, 0.8], max_l=4, n_pts=250)
exp_ortho = multi.MultiMicroscope(ill_thetas=[90], det_thetas=[0],
                                  det_nas=[1.1], max_l=4, n_pts=250)
exp_ortho2 = multi.MultiMicroscope(ill_thetas=[0], det_thetas=[90],
                                  det_nas=[0.8], max_l=4, n_pts=250)

# Load asymmetric dispim system matrix
folder = 'asym_dispim'
file_name = folder+'.dat'
data_file = folder+'/'+file_name

folder_ortho = 'ortho'
file_name_ortho = folder_ortho+'.dat'
data_file_ortho = folder_ortho+'/'+file_name_ortho

# Large volume (two actin strands)
xs, ys, zs = 50, 50, 20 
x0, y0, z0 = 180, 545, 150
d = 130

# X-Large volume (many actin strands)
# xs, ys, zs = 100, 100, 40 
# x0, y0, z0 = 150, 530, 140
# d = 130

# Small volumes (one actin strand)
# xs, ys, zs = 20, 20, 20, 
# x0, y0, z0 = 235, 580, 170
# d = 50

# Load and plot intensity fields
intf = data.IntensityField(file_names=input_files, x0=x0, y0=y0, z0=z0, xspan=xs, yspan=ys, zspan=zs, cal=cal)
intf.plot_int_field(output_file=folder+'/'+'data.pdf', shape=(2,4), line_start=(0, xs-1), line_end=(ys-1, 0), d=d, mag=2)

# Load system model
f = open(data_file, 'rb')
exp.psi = dill.load(f)

f = open(data_file_ortho, 'rb')
exp_ortho.psi = dill.load(f)

exp_ortho2.psi = exp.psi[4:, :]

# Calculate B
exp.calc_B_matrix()
exp_ortho.calc_B_matrix()
exp_ortho2.calc_B_matrix()

# Reconstruct phantom
mask = np.max(intf.g, axis=-1) > 0.6
df = exp.recon_dist_field(intf, mask=mask, prior='single')
df.plot_dist_field(exp.B, exp.xyz, exp.triangles,
                   filename=folder+'/data_both.png', r=1, d=100,
                   mag=4, show=True, mask=mask)

# Split data
intf_ortho = data.IntensityField()
intf_ortho.g = intf.g[:, :, :, :4]
df_ortho = exp_ortho.recon_dist_field(intf_ortho, mask=mask, prior='single')
df_ortho.plot_dist_field(exp_ortho.B, exp_ortho.xyz, exp_ortho.triangles,
                         filename=folder+'/data_ortho1.png', r=1, d=100,
                         mag=4, show=False, mask=mask)

intf_ortho2 = data.IntensityField()
intf_ortho2.g = intf.g[:, :, :, 4:]
df_ortho2 = exp_ortho2.recon_dist_field(intf_ortho2, mask=mask, prior='single')
df_ortho2.plot_dist_field(exp_ortho2.B, exp_ortho2.xyz, exp_ortho2.triangles,
                         filename=folder+'/data_ortho2.png', r=1, d=100,
                          mag=4, show=False, mask=mask)

# Collect results
ims = [ folder+'/data_ortho1.png', folder+'/data_ortho2.png', folder+'/data_both.png',]
fig, axs = plt.subplots(1, 3, figsize=(8, 3), gridspec_kw={'hspace':-0.3, 'wspace':-0.1})
for i, im in enumerate(ims):
    image = mpimg.imread(im)
    axs[i].imshow(image, interpolation=None, vmin=0, vmax=1)
    axs[i].set_axis_off()
axs[0].annotate('1.1 NA data only', xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
axs[1].annotate('0.71 NA data only', xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
axs[2].annotate('All data', xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
fig.savefig(folder+'/recons.pdf', dpi=1600)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
