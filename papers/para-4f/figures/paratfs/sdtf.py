import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.gridspec as gridspec
import scipy.misc
from scipy import special
import paratf

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times"

inches = 2
rows = 1
cols = 3
widths = [1,1,1]
heights = [1]
hspace = 0.05
wspace = 0.15

f = plt.figure(figsize=(inches*(np.sum(widths) + widths[0]*(wspace)*(cols - 1)),
                        inches*(np.sum(heights))))

outer_grid = gridspec.GridSpec(ncols=cols, nrows=rows,
                               width_ratios=widths, height_ratios=heights,
                               hspace=hspace, wspace=wspace)

alphas = [0.25, 0.5, 0.75]
alpha_strs = ['0.25', '0.5', '0.75']
for col in range(cols):
    alpha = alphas[col]
    ax = f.add_subplot(outer_grid[0, col])
    ax.set_xlim([0,1])
    ax.set_ylim([-0.125,1])

    thetas = [0, np.pi/8, np.pi/4, 3*np.pi/8, np.pi/2]
    theta_strs = ['0', '$\pi/8$', '$\pi/4$', '$3\pi/8$', '$\pi/2$']
    for theta in thetas:
        nucr = np.linspace(1e-3, 1, 100)
        ax.plot(nucr, paratf.sdtf(nucr, theta, NA=alpha)/paratf.sdtf(1e-3, np.pi/2, NA=alpha), '-', c=[1 - theta/(np.pi/2), 0, 0], lw=0.5)

    ax.plot([0,2],[0,0], 'k--', lw=0.25, zorder=1)

    # Monopole
    ax.plot(nucr, paratf.chat(0, nucr)/paratf.chat(0, 0), '-', c=[0, 1, 0], lw=1.5, zorder=0, alpha=1)        

    ax.tick_params(direction='out', left=True, right=False, width=0.75)
    ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    if col == 0:
        ax.set_yticklabels(['0','0.25','0.5','0.75','1'])
    else:
        ax.set_yticklabels([])
    ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0','0.25','0.5','0.75','1'])
    
    ax.text(0.5, 1.08, 'NA$/n_o ='+alpha_strs[col]+'$', ha='center', va='center', rotation=0, transform=ax.transAxes)

    if col == 0:
        for i, theta_str in enumerate(theta_strs):
            if i == 4:
                theta_str = '$\\vartheta =' + theta_str[1:-1] + '$'
            ax.text(0.75, 0.4 + 0.1*i, theta_str, ha='right', va='center', transform=ax.transAxes)
            l1 = lines.Line2D([0.8,0.9], 2*[0.4 + 0.1*i], c=[1 - thetas[i]/(np.pi/2), 0, 0], lw=0.5, transform=ax.transAxes)
            f.lines.extend([l1])

        ax.text(0.75, 0.8 + 0.1*1, 'Monopole', ha='right', va='center', transform=ax.transAxes)
        l1 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*1], c=[0, 1, 0], lw=1.5, transform=ax.transAxes)
        f.lines.extend([l1])

        
        ax.text(-0.4, 0.5, 'Renormalized\n dipole spatial transfer function \n $H(\\nu, \\vartheta)/H(0, \\pi/2)$', ha='center', va='center', rotation=90, transform=ax.transAxes)
        ax.text(0.5, -0.25, 'Scaled spatial frequency $\\nu/\\nu_c$', ha='center', va='center', transform=ax.transAxes)
        # ax.text(0.75, 0.8 + 0.1*1, '$\ell = 0$', ha='right', va='center', transform=ax.transAxes)
        # l1 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*1], c=[1, 0, 0], lw=0.5, transform=ax.transAxes)
        # ax.text(0.75, 0.8 + 0.1*0, '$\ell = 2$', ha='right', va='center', transform=ax.transAxes)
        # l2 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*0], c=[0, 0, 1], lw=0.5, transform=ax.transAxes)
        # f.lines.extend([l1, l2])

    if col != 0:
        ax.text(0.5, -0.25, '$\\nu/\\nu_c$', ha='center', va='center', transform=ax.transAxes)
f.savefig('sdtf.pdf', dpi=500, bbox_inches='tight')
