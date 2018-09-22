import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.misc
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times"


inches = 2
rows = 1
cols = 1
widths = [1]
heights = [1]
hspace = 0
wspace = 0

f = plt.figure(figsize=(inches*(np.sum(widths) + widths[0]*(wspace)*(cols - 1)),
                        inches*(np.sum(heights))))

outer_grid = gridspec.GridSpec(ncols=cols, nrows=rows,
                               width_ratios=widths, height_ratios=heights,
                               hspace=hspace, wspace=wspace)

def hexc(theta, alpha):
    return np.sin(theta)**2 + 0.5*(alpha**2)*(np.cos(theta)**2)
        
for col in range(cols):
    ax = f.add_subplot(outer_grid[0, col])
    ax.set_xlim([0,np.pi/2])
    ax.set_ylim([0,1])

    for alpha in [0, np.pi/8, np.pi/4, 3*np.pi/8]:
        theta = np.linspace(0, np.pi/2, 100)
        ax.plot(theta, hexc(theta, alpha), '-', lw=0.5, color=[alpha/(3*np.pi/8),0,0])

    ax.text(0.03, 0.85,'$\\alpha = 3\\pi/8$', ha='left', va='center', transform=ax.transData)
    ax.text(0.03, 0.4,'$\\pi/4$', ha='left', va='center', transform=ax.transData)
    ax.text(0.03, 0.18,'$\\pi/8$', ha='left', va='center', transform=ax.transData)
    ax.text(0.03, 0.035,'$0$', ha='left', va='center', transform=ax.transData)
    
    ax.tick_params(direction='out', left=True, right=False, width=0.75)
    ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0','0.25','0.5','0.75','1'])
    ax.xaxis.set_ticks([0, np.pi/8, np.pi/4, 3*np.pi/8, np.pi/2])
    ax.set_xticklabels(['0','$\pi/8$','$\pi/4$','$3\pi/8$','$\pi/2$'])

    ax.text(0.5, -0.25, 'Inclination angle $\\vartheta$', ha='center', va='center', transform=ax.transAxes)
    ax.text(-0.3, 0.5, 'Excitation kernel $h_{\\textrm{exc}}(\\vartheta)$', ha='center', va='center', rotation=90, transform=ax.transAxes)
    

f.savefig('excitation.pdf', dpi=500, bbox_inches='tight')
