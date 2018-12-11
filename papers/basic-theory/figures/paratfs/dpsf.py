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

def an(krda, n):
    return (n/np.pi)*((special.jn(n, krda)/(krda))**2)

def hdet(krd, theta, alpha):
    return an(krd*alpha, 1)*np.sin(theta)**2 + 0.5*an(krd*alpha, 2)*(alpha**2)*(np.cos(theta)**2)

alphas = [0.25, 0.5, 0.75]
alpha_strs = ['0.25', '0.5', '0.75']
for col in range(cols):
    alpha = alphas[col]
    ax = f.add_subplot(outer_grid[0, col])
    ax.set_xlim([0,2])
    ax.set_ylim([0,1])

    thetas = [0, np.pi/8, np.pi/4, 3*np.pi/8, np.pi/2]
    theta_strs = ['0', '$\pi/8$', '$\pi/4$', '$3\pi/8$', '$\pi/2$']
    for theta in thetas:
        nucr = np.linspace(1e-3, 2, 100)
        ax.plot(nucr, paratf.dpsf(nucr, theta, NA=alpha)/paratf.dpsf(1e-3, np.pi/2, NA=alpha), '-', c=[1 - theta/(np.pi/2), 0, 0], lw=0.5)

    # Monopole
    ax.plot(nucr, (paratf.jinc(0, nucr)/paratf.jinc(0, 0))**2, '-', c=[0, 1, 0], lw=1.5, zorder=0, alpha=1)        

    ax.tick_params(direction='out', left=True, right=False, width=0.75)
    ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    if col == 0:
        ax.set_yticklabels(['0','0.25','0.5','0.75','1'])
    else:
        ax.set_yticklabels([])
    ax.xaxis.set_ticks([0, 0.5, 1, 1.5, 2])
    ax.set_xticklabels(['0','0.5','1','1.5','2'])
    
    ax.text(0.5, 1.08, 'NA$/n_o ='+alpha_strs[col]+'$', ha='center', va='center', rotation=0, transform=ax.transAxes)
    if col == 0:
        ax.text(-0.4, 0.5, 'Renormalized\n dipole point spread function \n $h(r, \\vartheta)/h(0, \\pi/2)$', ha='center', va='center', rotation=90, transform=ax.transAxes)
        ax.text(0.5, -0.25, 'Scaled distance $\\nu_c r$', ha='center', va='center', transform=ax.transAxes)
        for i, theta_str in enumerate(theta_strs):
            if i == 4:
                theta_str = '$\\vartheta =' + theta_str[1:-1] + '$'
            ax.text(0.75, 0.4 + 0.1*i, theta_str, ha='right', va='center', transform=ax.transAxes)
            l1 = lines.Line2D([0.8,0.9], 2*[0.4 + 0.1*i], c=[1 - thetas[i]/(np.pi/2), 0, 0], lw=0.5, transform=ax.transAxes)
            f.lines.extend([l1])
        ax.text(0.75, 0.8 + 0.1*1, 'Monopole', ha='right', va='center', transform=ax.transAxes)
        l1 = lines.Line2D([0.8,0.9], 2*[0.8 + 0.1*1], c=[0, 1, 0], lw=1.5, transform=ax.transAxes)
        f.lines.extend([l1])


    if col != 0:
        ax.text(0.5, -0.25, '$\\nu_c r$', ha='center', va='center', transform=ax.transAxes)
f.savefig('dpsf.pdf', dpi=500, bbox_inches='tight')
