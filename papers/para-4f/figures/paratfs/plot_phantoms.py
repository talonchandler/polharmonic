from paratf import *
from phantoms import *
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.gridspec as gridspec
outer_gridspec = gridspec.GridSpec
inner_gridspec = gridspec.GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times"
plt.rc('path',simplify_threshold=.0001)

def plot_phantom(ph_i, phantom_func, x_head, x_labels, y_head, y_labels, y_tail=True):
    # Simulation settings
    lamb = 0.5 # um
    NA = 0.75
    n = 1.0
    nuc = 2*NA/lamb
    oversamp = 20
    dx = 1/(oversamp*2*nuc)#0.05#0.005 # um
    FOVmin = -1
    FOVmax = 4

    line = np.linspace(FOVmin, FOVmax, (FOVmax - FOVmin)/dx)
    rx, ry = np.meshgrid(line, line)

    # Plot settings
    def xbraces(ax):
        ax.text(-0.5, 4.65, x_head, ha='right', va='center', transform=ax.transData)
        for x in range(4):
            w = 0.4
            h = 0.05
            l1 = lines.Line2D([-w+x,w+x], [4.25, 4.25], c='k', lw=0.5, transform=ax.transData)
            l2 = lines.Line2D([w+x,w+x], [4.25, 4.25-h], c='k', lw=0.5, transform=ax.transData)
            l3 = lines.Line2D([-w+x,-w+x], [4.25-h, 4.25], c='k', lw=0.5, transform=ax.transData)
            l4 = lines.Line2D([x,x], [4.25+h, 4.25], c='k', lw=0.5, transform=ax.transData)
            f.lines.extend([l1, l2, l3, l4])
            ax.text(x, 4.65, x_labels[x], ha='center', va='center', transform=ax.transData)

    def ybraces(ax):
        ax.text(4.75, 4.25, y_head, ha='center', va='center', transform=ax.transData, rotation=-90)
        if y_tail:
            ax.text(4.75, -1.5, '($\mu$m)', ha='center', va='center', transform=ax.transData, rotation=-90)
        for x in range(4):
            w = 0.4
            h = 0.05
            l1 = lines.Line2D([4.25, 4.25], [-w+x,w+x], c='k', lw=0.5, transform=ax.transData)
            l2 = lines.Line2D([4.25, 4.25-h], [w+x,w+x], c='k', lw=0.5, transform=ax.transData)
            l3 = lines.Line2D([4.25-h, 4.25], [-w+x,-w+x], c='k', lw=0.5, transform=ax.transData)
            l4 = lines.Line2D([4.25+h, 4.25], [x,x], c='k', lw=0.5, transform=ax.transData)
            f.lines.extend([l1, l2, l3, l4])
            ax.text(4.75, x, y_labels[x], ha='center', va='center', transform=ax.transData, rotation=-90)

    inches = 2
    outer_rows = 1
    outer_cols = 3
    inner_cols = 2
    widths = [1,1,1]
    inner_widths = [1,0.05]
    heights = [1]
    hspace = 0.05
    wspace = 0.15
    wspace_inner = 0.05

    f = plt.figure(figsize=(inches*(np.sum(widths) + widths[0]*(wspace)*(outer_cols - 1) + 0.5*widths[0]*wspace_inner*outer_cols),
                            inches*(np.sum(heights)) + heights[0]*hspace*(outer_rows - 1)))

    outer_grid = gridspec.GridSpec(ncols=outer_cols, nrows=outer_rows,
                                   width_ratios=widths, height_ratios=heights,
                                   hspace=hspace, wspace=wspace)

    for outer_col in range(outer_cols):
        inner_grid = inner_gridspec(1, inner_cols,
                                    subplot_spec=outer_grid[0, outer_col],
                                    width_ratios=inner_widths, height_ratios=[1],
                                    hspace=0, wspace=wspace_inner)
        for inner_col in range(inner_cols):        
            ax = f.add_subplot(inner_grid[0, inner_col])
            if inner_col == 0:
                xbraces(ax)
                if outer_col == 0:
                    ax.text(0.5, 1.3, 'Dipole density: $f_{(ph'+str(ph_i)+')}$', ha='center', va='center', transform=ax.transAxes)
                    if ph_i in [2, 4]:
                        ff = phantom2f(rx, ry)        
                        im = ax.imshow(ff, cmap='gray', origin='lower', interpolation='none', extent=[FOVmin, FOVmax, FOVmin, FOVmax])
                    else:
                        ff = np.zeros_like(rx)
                        im = ax.imshow(ff, cmap='gray', origin='lower', interpolation='none', extent=[FOVmin, FOVmax, FOVmin, FOVmax])
                        xx = np.array([0,1,2,3])
                        for i in range(4):
                            ax.plot(xx, 4*[i], 'wx', markersize=4)
                        
                if outer_col == 1:
                    ax.text(0.5, 1.3, 'Scaled irradiance: $g_{(ph'+str(ph_i)+')}$', ha='center', va='center', transform=ax.transAxes)
                    g = phantom_func(rx, ry, NA=NA, n=n, nuc=nuc, dx=dx)
                    im = ax.imshow(g, cmap='gray', origin='lower', interpolation='none', extent=[FOVmin, FOVmax, FOVmin, FOVmax])
                if outer_col == 2:
                    ax.text(0.5, 1.3, 'Scaled irradiance profiles: $g_{(ph'+str(ph_i)+')}$', ha='center', va='center', transform=ax.transAxes)
                    ybraces(ax)
                    for i in range(4):
                        yind = np.argmin(np.abs(line - i))
                        y = g[yind,:]
                        ax.plot(line, y/np.max(g) + i, '-k', lw=0.5)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.yaxis.set_ticks([])
                ax.xaxis.set_ticks([0, 1, 2, 3])
                ax.set_xticklabels(['0','1', '2', '3'])
                ax.yaxis.set_ticks([0, 1, 2, 3])
                ax.set_xlabel('$x$ position ($\mu$m)')
                if outer_col == 0:
                    ax.set_ylabel('$y$ position ($\mu$m)')
                    ax.set_yticklabels(['0','1', '2', '3'])
                else:
                    ax.set_yticklabels([])
                ax.set_xlim([-1,4])
                ax.set_ylim([-1,4])
                ax.set_aspect('equal')

            if inner_col == 1:
                if outer_col == 2:
                    ax.axis('off')
                else:
                    cc = np.linspace(0,1,100)
                    if ph_i in [1, 3]:
                        if outer_col == 1:
                            im = ax.imshow(np.array([cc,cc]).T, cmap='gray', origin='lower', interpolation='none', aspect=20, extent=[0,1,0,1])
                        else:
                            ax.axis('off')
                    else:
                        im = ax.imshow(np.array([cc,cc]).T, cmap='gray', origin='lower', interpolation='none', aspect=20, extent=[0,1,0,1])
                    ax.set_xlim([0,1])
                    ax.set_ylim([0,1])
                    ax.xaxis.set_ticks([])
                    ax.set_xticklabels([])
                    ax.yaxis.set_ticks([0, 1])
                    ax.set_yticklabels(['0', '1'])
                    ax.yaxis.tick_right()

    f.savefig('ph'+str(ph_i)+'.pdf', dpi=500, bbox_inches='tight')

plot_phantom(1, phantom1g,
             '$\\vartheta = \quad$', ['$0$', '$\\frac{\\pi}{6}$', '$\\frac{\\pi}{3}$', '$\\frac{\\pi}{2}$'],
             '$\\varphi= \quad$', ['$0$', '$\\frac{\\pi}{4}$', '$\\frac{\\pi}{2}$', '$\\frac{3\\pi}{4}$'],
             y_tail=False)
plot_phantom(2, phantom2g,
             '$\\vartheta = \quad$', ['$0$', '$\\frac{\\pi}{6}$', '$\\frac{\\pi}{3}$', '$\\frac{\\pi}{2}$'],
             '$D= \quad$', ['$0.15$', '$0.3$', '$0.45$', '$0.6$'])
plot_phantom(3, phantom3g,
             '$\\vartheta\' = \quad$', ['$0$', '$\\frac{\\pi}{6}$', '$\\frac{\\pi}{3}$', '$\\frac{\\pi}{2}$'],
             '$\\Delta= \quad$', ['$0$', '$\\frac{\\pi}{6}$', '$\\frac{\\pi}{3}$', '$\\frac{\\pi}{2}$'],
             y_tail=False)
plot_phantom(4, phantom4g,
             '$\\Delta = \quad$', ['$0$', '$\\frac{\\pi}{6}$', '$\\frac{\\pi}{3}$', '$\\frac{\\pi}{2}$'],
             '$D= \quad$', ['$0.15$', '$0.3$', '$0.45$', '$0.6$'])
