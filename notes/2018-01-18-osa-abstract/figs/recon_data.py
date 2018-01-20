from polharmonic import data, util, multi
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import dill
import matplotlib.image as mpimg
np.set_printoptions(precision=3, suppress=True)

# High level params
folder = 'asym_dispim'
int_field_mag = 1
recon_mag = 1
recon_mask_threshold = 0.6

# Import data
name_head = '/Users/Talon/GoogleDrive/projects/dispim-data/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif',
         'B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif']
input_files = np.array([name_head + name for name in names])

# Calibration data
cal = np.array([1.06, 1.0, 1.0, 1.03, 
                1.08, 1.05, 1.0, 1.04])

# Reindexing 
idx = [0, 3, 2, 1, 4, 5, 6, 7]
input_files = input_files[idx]
cal = cal[idx]

# Whole volume
# xs, ys, zs = 499, 800, 338
# x0, y0, z0 = 0, 0, 0
# d = 300

# ROI 1
# Large volume (two actin strands)
# xs, ys, zs = 50, 50, 20 
# x0, y0, z0 = 180, 545, 150
# d = 130
# recon_mask_threshold = 0.6
# skip = 2

# X-large volume
xs, ys, zs = 200, 250, 100 
x0, y0, z0 = 100, 500, 100
d = 130
skip = 4
recon_mask_threshold = 0.25

# Small volumes (one actin strand)
# xs, ys, zs = 20, 20, 20, 
# x0, y0, z0 = 235, 580, 170
# d = 50
# skip = 1

# # ROI 2 (same error - ~25 degrees off for ~(1, 0, 1) actin filaments
# xs, ys, zs = 50, 50, 20 
# x0, y0, z0 = 170, 150, 150
# d = 130

# # ROI 3 (slightly smaller ~15 degree error for ~(1, 0, 1)
# xs, ys, zs = 50, 50, 20 
# x0, y0, z0 = 386, 744, 285
# d = 130

# # ROI 4 (larger error ~30 degrees off for ~(1, 0, 1)
# xs, ys, zs = 50, 50, 20 
# x0, y0, z0 = 138, 185, 125
# d = 130
# recon_mask_threshold = 0.5

# ROI 5 
# xs, ys, zs = 50, 50, 20 
# x0, y0, z0 = 357, 243, 265
# d = 130
# recon_mask_threshold = 0.45
# skip = 2

# Load and plot intensity fields
intf = data.IntensityField()
intf.load_from_file(file_names=input_files, x0=x0, y0=y0, z0=z0, xspan=xs, yspan=ys, zspan=zs, cal=cal)
row_labels = ['$x$-illumination\n $z$-detection\n 1.1 NA', '$z$-illumination\n $x$-detection\n 0.71 NA']
col_labels = ['$\phi = 0^{\circ}$', '$\phi = 45^{\circ}$', '$\phi = 90^{\circ}$', '$\phi = 135^{\circ}$']
# intf.plot(output_file=folder+'/'+'data.pdf', shape=(2,4),
#           row_labels=row_labels, col_labels=col_labels,
#           d=d, mag=int_field_mag, dpi=400, show=False)

# Build microscopes
exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                            det_nas=[1.1, 0.8], max_l=4, n_pts=250)
exp_ortho = multi.MultiMicroscope(ill_thetas=[90], det_thetas=[0],
                                  det_nas=[1.1], max_l=4, n_pts=250)
exp_ortho2 = multi.MultiMicroscope(ill_thetas=[0], det_thetas=[90],
                                  det_nas=[0.8], max_l=4, n_pts=250)

# Load asymmetric dispim system matrix
file_name = folder+'.dat'
data_file = folder+'/'+file_name

# Load system models
f = open(data_file, 'rb')
exp.psi = dill.load(f)
exp_ortho.psi = exp.psi[:4, :]
exp_ortho2.psi = exp.psi[4:, :]

# Calculate B
exp.calc_B_matrix()
exp_ortho.calc_B_matrix()
exp_ortho2.calc_B_matrix()

# Reconstruct with all data
threshold_mask = np.max(intf.g, axis=-1) > recon_mask_threshold
sparse_mask = np.zeros(threshold_mask.shape)
sparse_mask[::skip, ::skip, ::skip] = 1
mask = np.logical_and(threshold_mask, sparse_mask)

df = exp.recon_dist_field(intf, mask=mask, prior='single')
df.plot_dist_field(exp.B, exp.xyz, exp.triangles,
                   filename=folder+'/data_both.png', r=skip, d=d,
                   mag=recon_mag, show=True, mask=mask)

# Reconstruct with single view data
intf_ortho = data.IntensityField(g=intf.g[:, :, :, :4])
df_ortho = exp_ortho.recon_dist_field(intf_ortho, mask=mask, prior='single')
df_ortho.plot_dist_field(exp_ortho.B, exp_ortho.xyz, exp_ortho.triangles,
                         filename=folder+'/data_ortho1.png', r=skip, d=d,
                         mag=recon_mag, show=False, mask=mask)

intf_ortho2 = data.IntensityField(g=intf.g[:, :, :, 4:])
df_ortho2 = exp_ortho2.recon_dist_field(intf_ortho2, mask=mask, prior='single')
df_ortho2.plot_dist_field(exp_ortho2.B, exp_ortho2.xyz, exp_ortho2.triangles,
                         filename=folder+'/data_ortho2.png', r=skip, d=d,
                          mag=recon_mag, show=False, mask=mask)

# Plot reconstruction results
ims = [ folder+'/data_ortho1.png', folder+'/data_ortho2.png', folder+'/data_both.png',]
fig, axs = plt.subplots(1, 3, figsize=(12, 4), gridspec_kw={'hspace':0, 'wspace':0})
for i, im in enumerate(ims):
    image = mpimg.imread(im)
    axs[i].imshow(image, interpolation=None, vmin=0, vmax=1)
    axs[i].set_axis_off()
axs[0].annotate('1.1 NA data only', xy=(0,0), xytext=(0.5, 1.07), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
axs[1].annotate('0.71 NA data only', xy=(0,0), xytext=(0.5, 1.07), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
axs[2].annotate('All data', xy=(0,0), xytext=(0.5, 1.07), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
fig.savefig(folder+'/recon.pdf', dpi=800, bbox_inches='tight')

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
