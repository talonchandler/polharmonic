from polharmonic import data, util, multi
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import dill
import matplotlib.image as mpimg
np.set_printoptions(precision=3, suppress=True)
import matplotlib as mpl
import subprocess
mpl.rcParams['font.serif'] = 'Times'

def recon_roi(name, xs, ys, zs, x0, y0, z0, recon_mask_threshold, d, skip, dpi=800, mag=1, note=''):
    # High level params
    mfolder = 'asym_dispim'
    folder = name
    if not os.path.exists(folder):
        os.mkdir(folder)
    int_field_mag = mag
    recon_mag = mag

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

    # Write to logfile
    f = open(folder+'/log.txt', 'w')
    f.write('name: ' + folder + '\n')
    f.write('note: ' + note + '\n')        
    f.write('roi coords (px) = ' + str(xs) + ', ' + str(ys) + ', ' + str(zs) + '\n' )
    f.write('span (px) = ' + str(xs) + ', ' + str(ys) + ', ' + str(zs) + '\n' )
    f.write('span (um) = ' + str(0.135*xs) + ', ' + str(0.135*ys) + ', ' + str(0.135*zs) + '\n' )
    f.write('skip (px) = ' + str(skip) + '\n')
    f.write('threshold (normalized) = ' + str(recon_mask_threshold) + '\n')

    f.close()

    # Load and plot intensity fields
    intf = data.IntensityField()
    intf.load_from_file(file_names=input_files, x0=x0, y0=y0, z0=z0, xspan=xs, yspan=ys, zspan=zs, cal=cal)
    row_labels = ['$x$-illumination\n $z$-detection\n 1.1 NA', '$z$-illumination\n $x$-detection\n 0.71 NA']
    col_labels = ['$\phi_{\mathrm{pol}} = 0^{\circ}$', '$\phi_{\mathrm{pol}} = 45^{\circ}$', '$\phi_{\mathrm{pol}} = 90^{\circ}$', '$\phi_{\mathrm{pol}} = 135^{\circ}$']
    intf.plot(output_file=folder+'/'+'data.pdf', shape=(2,4),
              row_labels=row_labels, col_labels=col_labels,
              d=d, mag=int_field_mag, dpi=400, show=False)
    import pdb; pdb.set_trace() 

    # Build microscopes
    exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                                det_nas=[1.1, 0.8], max_l=4, n_pts=250)
    exp_ortho = multi.MultiMicroscope(ill_thetas=[90], det_thetas=[0],
                                      det_nas=[1.1], max_l=4, n_pts=250)
    exp_ortho2 = multi.MultiMicroscope(ill_thetas=[0], det_thetas=[90],
                                      det_nas=[0.8], max_l=4, n_pts=250)

    # Load asymmetric dispim system matrix
    file_name = mfolder+'.dat'
    data_file = mfolder+'/'+file_name

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
                       mag=recon_mag, show=False, mask=mask)

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

    # Setup viewing window

    from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                      mark_inset)

    # Plot reconstruction results
    ims = [ folder+'/data_ortho1.png', folder+'/data_ortho2.png', folder+'/data_both.png',]

    ncols = 3
    nrows = 1

    inset_hfrac = .15
    inset_vfrac = .15

    inset_hfrac_offset = .77
    inset_vfrac_offset = .32

    top_pad = 0
    bottom_pad = 0
    left_pad = 0
    right_pad = 0

    hspace = .0
    vspace = .0

    ax_width = (1 - left_pad - right_pad - (ncols - 1) * hspace) / ncols
    ax_height = (1 - top_pad - bottom_pad - (nrows - 1) * vspace) / nrows

    fig = plt.figure()

    ax_lst = []
    for j in range(ncols):
        for k in range(nrows):
            a_bottom = bottom_pad + k * ( ax_height + vspace)
            a_left = left_pad + j * (ax_width + hspace)

            inset_bottom = a_bottom + inset_vfrac_offset * ax_height
            inset_left = a_left + inset_hfrac_offset * ax_width

            ax = fig.add_axes([a_left, a_bottom, ax_width, ax_height])
            ax_in = fig.add_axes([inset_left, inset_bottom, ax_width * inset_hfrac, ax_height *  inset_vfrac])
            ax_lst.append((ax,ax_in))

    for i, (ax, ax_in) in enumerate(ax_lst):
        image = mpimg.imread(ims[i])
        ax.imshow(image, interpolation=None, vmin=0, vmax=1)
        ax.set_axis_off()

        xyz, tp = util.fibonacci_sphere(1000, xyz=True)
        sphere_im = util.plot_sphere(filename='scale.png', directions=tp, data=xyz,
                                     show=False, vis_px=500)
        ax_in.imshow(sphere_im, interpolation='none')
        ax_in.set_axis_off()

    ax_lst[0][0].annotate('1.1 NA data only', xy=(0,0), xytext=(0.5, 1.04), textcoords='axes fraction', va='center', ha='center', fontsize=7, annotation_clip=False)
    ax_lst[1][0].annotate('0.71 NA data only', xy=(0,0), xytext=(0.5, 1.04), textcoords='axes fraction', va='center', ha='center', fontsize=7, annotation_clip=False)
    ax_lst[2][0].annotate('All data', xy=(0,0), xytext=(0.5, 1.04), textcoords='axes fraction', va='center', ha='center', fontsize=7, annotation_clip=False)
    fig.savefig(folder+'/recon.pdf', dpi=dpi, bbox_inches='tight')

    subprocess.Popen('rm '+folder+'/*.png', shell=True)

    print('Total time: '+str(np.round(time.time() - start, 2)))
    os.system('say "done"')

# Main calls
# recon_roi(name='whole', xs=499, ys=800, zs=338, x0=0, y0=0, z0=0,
#           recon_mask_threshold=0.25, d=300, skip=8,
#           dpi=800, mag=3)

# recon_roi(name='roi1', xs=50, ys=50, zs=20, x0=180, y0=545, z0=150,
#           recon_mask_threshold=0.6, d=130, skip=2,
#           dpi=800, mag=3, note='roi from previous email')

recon_roi(name='roi2', xs=100, ys=100, zs=40, x0=180, y0=545, z0=150,
          recon_mask_threshold=0.4, d=130, skip=3,
          dpi=1000, mag=5, note='zoomed out version of roi1')

# recon_roi(name='roi3', xs=200, ys=250, zs=100, x0=100, y0=500, z0=100,
#           recon_mask_threshold=0.25, d=130, skip=4,
#           dpi=800, mag=3, note='~half of a cell body')

# recon_roi(name='roi4', xs=100, ys=100, zs=40, x0=325, y0=225, z0=240,
#           recon_mask_threshold=0.3, d=130, skip=3,
#           dpi=800, mag=3, note='edge of cell...attachment to cover slip?')
