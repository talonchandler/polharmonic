import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.text as text
import matplotlib.gridspec as gridspec
from scipy import special
outer_gridspec = gridspec.GridSpec
inner_gridspec = gridspec.GridSpecFromSubplotSpec
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times"

# Physical parameters
lamb = 0.3
k = 2*np.pi/lamb
alpha = 0.5

inches = 3.5
outer_rows = 1
outer_cols = 4
inner_rows = 5
inner_cols = 2
widths = [(1/3)*(6/5)]*outer_cols
heights = [1]*outer_rows
wspace = 0.1
wspace_inner = 0.2
hspace = 0
f = plt.figure(figsize=(inches*(np.sum(widths) + widths[0]*(wspace)*(outer_cols - 1) + 0.5*widths[0]*wspace_inner*outer_cols),
                        inches*(np.sum(heights)) + heights[0]*hspace*(outer_rows - 1)))
                                
outer_grid = outer_gridspec(ncols=outer_cols, nrows=outer_rows,
                            width_ratios=widths, height_ratios=heights,
                            hspace=hspace, wspace=wspace)
for outer_row in range(outer_rows):
    for outer_col in range(outer_cols):
        inner_grid = inner_gridspec(inner_rows, inner_cols,
                                    subplot_spec=outer_grid[outer_row, outer_col],
                                    width_ratios=[1]*inner_cols, height_ratios=[1]*inner_rows,
                                    hspace=0, wspace=wspace_inner)
        for inner_row in range(inner_rows):
            for inner_col in range(inner_cols):
                ax = f.add_subplot(inner_grid[inner_row, inner_col])
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.axis('off')

                # Draw lenses
                if inner_col == 0 and (inner_row == 1 or inner_row == 3):
                    # Kludge to smooth sharp edges
                    x1 = np.linspace(0.5,1,50)
                    x2 = np.linspace(1,0.5,50)
                    x3 = np.linspace(0.5,0,50)
                    x4 = np.linspace(0,0.5,50)
                    y1 = -0.4*x1*(x1-1) + 0.5
                    y2 = 0.4*x2*(x2-1) + 0.5
                    y3 = 0.4*x3*(x3-1) + 0.5
                    y4 = -0.4*x4*(x4-1) + 0.5
                    xx = np.concatenate([x1, x2, x3, x4]) 
                    yy = np.concatenate([y1, y2, y3, y4])
                    l1 = lines.Line2D(xx, yy, c='k', lw=1, transform=ax.transAxes)
                    f.lines.extend([l1])

                # Draw labels
                def f_lines(ax):
                    l1 = lines.Line2D([-0.15,-0.15], [0.2,0.5], c='k', lw=1, transform=ax.transAxes)
                    l2 = lines.Line2D([-0.2,-0.1], [0.5,0.5], c='k', lw=1, transform=ax.transAxes)
                    l3 = lines.Line2D([-0.15,-0.15], [-0.2,-0.5], c='k', lw=1, transform=ax.transAxes)
                    l4 = lines.Line2D([-0.2,-0.1], [-0.5,-0.5], c='k', lw=1, transform=ax.transAxes)
                    f.lines.extend([l1, l2, l3, l4])
                
                if outer_col == 0 and inner_col == 0:
                    row_labels = ['Detector\n irradiance',
                                  'Tube\n lens',
                                  'Back focal\n plane\n electric field',
                                  'Objective\n lens',
                                  'Far-field\n radiation\n pattern']
                    ax.text(-0.3,0.5,row_labels[inner_row], ha='right', va='center', transform=ax.transAxes)
                    if inner_row < 4:
                        ax.text(-0.15, 0, '$f$', ha='center', va='center', transform=ax.transAxes)
                        f_lines(ax)

                def brace_lines(ax):
                    wsi = wspace_inner
                    l1 = lines.Line2D([0,2+wsi], [-0.2,-0.2], c='k', lw=1, transform=ax.transAxes)
                    l2 = lines.Line2D([2+wsi,2+wsi], [-0.2,-0.1], c='k', lw=1, transform=ax.transAxes)
                    l3 = lines.Line2D([0,0], [-0.2,-0.1], c='k', lw=1, transform=ax.transAxes)
                    l4 = lines.Line2D([1+wsi/2,1+wsi/2], [-0.2,-0.3], c='k', lw=1, transform=ax.transAxes)
                    f.lines.extend([l1, l2, l3, l4])

                if inner_row == 4:
                    if outer_col == 0 and inner_col == 0:
                        brace_lines(ax)
                        ax.text(1.0+wspace_inner/2,-0.4,'a) Monopole radiator', ha='center', va='top', transform=ax.transAxes)
                    if outer_col == 1 and inner_col == 0:
                        brace_lines(ax)
                        ax.text(1.0+wspace_inner/2,-0.4,'b) $\\vartheta = \pi/2$\n dipole radiator', ha='center', va='top', transform=ax.transAxes)
                    if outer_col == 3 and inner_col == 0:
                        brace_lines(ax)                        
                        ax.text(1.0+wspace_inner/2,-0.4,'d) $\\vartheta = \pi/10$\n dipole radiator', ha='center', va='top', transform=ax.transAxes)
                    if outer_col == 2 and inner_col == 0:
                        brace_lines(ax)                        
                        ax.text(1.0+wspace_inner/2,-0.4,'c) $\\vartheta = 0$\n dipole radiator', ha='center', va='top', transform=ax.transAxes)

                if inner_row == 0:
                    if inner_col == 0:
                        ax.text(0.5,1.2,'Axial\n view', ha='center', va='bottom', transform=ax.transAxes)
                    else:
                        ax.text(0.5,1.2,'Transverse\n view', ha='center', va='bottom', transform=ax.transAxes)

                # Draw detector planes
                if inner_row == 0:
                    ax.set_xlim([-1,1])
                    ax.set_ylim([-1,1])
                    N = 1000
                    x = np.linspace(-1, 1, N)
                    xx, yy = np.meshgrid(x, x)
                    r = np.sqrt(xx**2 + yy**2)
                    # a = 25
                    a = k*alpha
                    if outer_col == 0 or outer_col == 1:
                        z = (special.jn(1,a*r)/(a*r))**2
                    if outer_col == 3:
                        theta = np.pi/10
                        z = (1 - np.cos(theta)**2)*(special.jn(1,a*r)/(a*r))**2 + 0.5*np.cos(theta)**2*(special.jn(2,a*r)/(a*r))**2
                    if outer_col == 2:
                        z = (special.jn(2,a*r)/(a*r))**2
                    if inner_col == 0:
                        l1 = lines.Line2D(x, z[:,N//2]/np.max(z), c='k', lw=1, transform=ax.transData, zorder=0)
                        f.lines.extend([l1])
                    if inner_col == 1:
                        ax.imshow(z/np.max(z), vmin=0, vmax=1, cmap='gray', interpolation='bicubic', origin='lower', extent=[-1,1,-1,1])

                # Draw dots
                if inner_row == 0 or inner_row == 2:
                    ax.scatter([0], [0.5], c='k', marker='o', clip_on=False, s=[8], zorder=10, transform=ax.transAxes)
                    ax.scatter([1], [0.5], c='w', marker='o', edgecolor='k', lw=0.5, clip_on=False, s=[10], zorder=10, transform=ax.transAxes)

                # Draw back focal planes
                if inner_row == 2:
                    ecolor = [1,.3,.3]
                    ax.set_xlim([-1,1])
                    ax.set_ylim([-1,1])                    
                    grid = np.array([-0.66,-0.33,0,0.33,0.66])
                    if inner_col == 1:
                        theta = np.linspace(0, 2*np.pi, 100)
                        x = np.cos(theta)
                        y = np.sin(theta)
                        ax.plot(x, y, '-k', lw=1, clip_on=False, zorder=1)
                    if outer_col == 0 and inner_col == 0:
                        ax.plot([-1,1], [0,0], '-', color=[1,0.5,0.5], zorder=0)
                    if outer_col == 0 and inner_col == 1:
                        theta = np.linspace(0, 2*np.pi, 100)
                        x = np.cos(theta)
                        y = np.sin(theta)
                        ax.fill(x, y, c=ecolor, clip_on=False, zorder=0)
                        ax.fill(x, -y, c=ecolor, clip_on=False, zorder=0)
                    if outer_col == 1 and inner_col == 0:
                        y = np.zeros(grid.shape)
                        u = 0.3*np.ones(grid.shape)
                        v = y
                        ax.quiver(grid, y, u, v, angles='xy', scale_units='xy', scale=1, color=ecolor, pivot='mid', width=0.02)
                    if outer_col == 1 and inner_col == 1:
                        xi, yi = np.meshgrid(grid, grid)
                        x = xi[xi**2 + yi**2 < 0.7]
                        y = yi[xi**2 + yi**2 < 0.7]
                        u = 0.3*np.ones(x.shape)
                        v = np.zeros(y.shape)
                        ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1, color=ecolor, pivot='mid', width=0.02)
                    if outer_col == 3 and inner_col == 0:
                        theta = np.pi/10
                        y = np.zeros(grid.shape)
                        u = -0.1*np.ones(y.shape) + 0.3*grid*np.cos(theta)
                        v = y
                        ax.quiver(grid, y, u, v, angles='xy', scale_units='xy', scale=1, color=ecolor, pivot='mid', width=0.02)
                    if outer_col == 2 and inner_col == 0:
                        y = np.zeros(grid.shape)
                        u = 0.4*grid
                        v = y
                        ax.quiver(grid, y, u, v, angles='xy', scale_units='xy', scale=1, color=ecolor, pivot='mid', width=0.02)
                    if outer_col == 3 and inner_col == 1:
                        xi, yi = np.meshgrid(grid, grid)
                        x = xi[xi**2 + yi**2 < .7]
                        y = yi[xi**2 + yi**2 < .7]
                        theta = np.pi/10
                        u = -0.1*np.ones(x.shape) + 0.3*x*np.cos(theta)
                        v = 0.3*y*np.sin(theta)
                        ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1, color=ecolor, pivot='mid', width=0.02)
                    if outer_col == 2 and inner_col == 1:
                        xi, yi = np.meshgrid(grid, grid)
                        x = xi[xi**2 + yi**2 < .7]
                        y = yi[xi**2 + yi**2 < .7]
                        u = 0.4*x
                        v = 0.4*y
                        ax.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1, color=ecolor, pivot='mid', width=0.02)

                # Draw far-field radiation patterns
                if inner_row == 4 and inner_col == 0:
                    N = 1000
                    x, y = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
                    r = np.sqrt(x**2 + y**2)
                    theta = np.arctan2(y, x)
                    ax.set_xlim([-1,1])
                    ax.set_ylim([-1,1])
                    if outer_col == 0:
                        z = np.where(r < 1, np.cos(2*np.pi*r/lamb)/(0.1 + r), 0)
                        l1 = lines.Line2D([0,0,1,0], [0,2,2,0], c='k', lw=1.0, transform=ax.transData, zorder=3, ls=':')
                        tt = np.linspace(np.arctan(2),np.pi/2,20)
                        l4 = lines.Line2D(0.7*np.cos(tt), 0.7*np.sin(tt), c='k', lw=1.0, transform=ax.transData, zorder=3)
                        ax.text(0.63,1.1,'$\\alpha$', ha='center', va='center', transform=ax.transAxes)
                        ax.text(1.2,1.5,'$\\sin\\alpha \\approx \\alpha$', ha='left', va='center', transform=ax.transAxes)
                        f.lines.extend([l1, l4])
                    if outer_col == 1:
                        z = np.where(r < 1, np.cos(theta - np.pi/2)**2*np.cos(2*np.pi*r/lamb)/(0.1 + r), 0)
                    if outer_col == 3:
                        z = np.where(r < 1, np.cos(theta + np.pi/10)**2*np.cos(2*np.pi*r/lamb)/(0.1 + r), 0)
                    if outer_col == 2:
                        z = np.where(r < 1, np.cos(theta)**2*np.cos(2*np.pi*r/lamb)/(0.1 + r), 0)
                    ax.imshow(z, vmin=-np.max(z), vmax=np.max(z), cmap='bwr', extent=[-1,1,-1,1], interpolation='bicubic', origin='lower')
                    
                    
f.savefig('microscope.pdf', dpi=500, bbox_inches='tight')
