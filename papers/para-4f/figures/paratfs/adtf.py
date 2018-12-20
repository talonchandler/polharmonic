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
    ax.set_xlim([0,2])
    ax.set_ylim([0,1])

    nucr = np.linspace(1e-3, 2, 100)
    ax.plot(nucr, paratf.adtf(nucr, 0, NA=alpha)/paratf.adtf(1e-3, 0, NA=alpha), '-', c=[1, 0, 0], lw=0.5)
    ax.plot(nucr, paratf.adtf(nucr, 2, NA=alpha)/paratf.adtf(1e-3, 0, NA=alpha), '-', c=[0, 0, 1], lw=0.5)
    ax.plot([0,2],[0,0], 'k--', lw=0.25, zorder=1)

    # Monopole
    ax.plot(nucr, (paratf.jinc(0, nucr)/paratf.jinc(0, 0))**2, '-', c=[0, 1, 0], lw=1.5, zorder=0, alpha=1)        

    ax.tick_params(direction='out', left=True, right=False, width=0.75)
    ax.yaxis.set_ticks([-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
    if col == 0:
        ax.set_yticklabels(['-0.5', '-0.25', '0','0.25','0.5','0.75','1'])
    else:
        ax.set_yticklabels([])
    ax.xaxis.set_ticks([0, 0.5, 1, 1.5, 2])
    ax.set_xticklabels(['0','0.5','1','1.5','2'])
    
    ax.text(0.5, 1.08, 'NA$/n_o ='+alpha_strs[col]+'$', ha='center', va='center', rotation=0, transform=ax.transAxes)
    if col == 0:
        ax.text(-0.4, 0.5, 'Renormalized\n dipole angular transfer function \n $H_{\ell}^m(r)/H_0^0(0)$', ha='center', va='center', rotation=90, transform=ax.transAxes)
        ax.text(0.5, -0.25, 'Scaled distance $\\nu_c r$', ha='center', va='center', transform=ax.transAxes)
        ax.text(0.75, 0.8 + 0.1*0, '$\ell = 0$', ha='right', va='center', transform=ax.transAxes)
        l1 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*0], c=[1, 0, 0], lw=0.5, transform=ax.transAxes)
        ax.text(0.75, 0.8 + 0.1*-1, '$\ell = 2$', ha='right', va='center', transform=ax.transAxes)
        l2 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*-1], c=[0, 0, 1], lw=0.5, transform=ax.transAxes)
        ax.text(0.75, 0.8 + 0.1*1, 'Monopole', ha='right', va='center', transform=ax.transAxes)
        l3 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*1], c=[0, 1, 0], lw=1.5, transform=ax.transAxes)
        f.lines.extend([l1, l2, l3])

    if col != 0:
        ax.text(0.5, -0.25, '$\\nu_c r$', ha='center', va='center', transform=ax.transAxes)
f.savefig('adtf.pdf', dpi=500, bbox_inches='tight')
