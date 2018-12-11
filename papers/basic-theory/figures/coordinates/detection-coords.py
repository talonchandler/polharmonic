import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.text as text
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
patches = []
from scipy import special
outer_gridspec = gridspec.GridSpec
inner_gridspec = gridspec.GridSpecFromSubplotSpec
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times"
plt.rcParams["font.size"] = 11
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{bm}"]

# Physical parameters
lamb = 0.3
k = 2*np.pi/lamb
alpha = 0.5

inches = 2
outer_rows = 1
outer_cols = 1
inner_rows = 2
inner_cols = 7
widths = [7/2]*outer_cols
heights = [1]*outer_rows
wspace = 0
wspace_inner = 0
hspace = 0
hspace_inner = 0.2
f = plt.figure(figsize=(inches*(np.sum(widths) + widths[0]*(wspace)*(outer_cols - 1) + 0.5*widths[0]*wspace_inner*outer_cols),
                        inches*(np.sum(heights) + heights[0]*(hspace)*(outer_rows - 1) + 0.5*heights[0]*hspace_inner*outer_rows)))

                                
outer_grid = outer_gridspec(ncols=outer_cols, nrows=outer_rows,
                            width_ratios=widths, height_ratios=heights,
                            hspace=hspace, wspace=wspace)
for outer_row in range(outer_rows):
    for outer_col in range(outer_cols):
        inner_grid = inner_gridspec(inner_rows, inner_cols,
                                    subplot_spec=outer_grid[outer_row, outer_col],
                                    width_ratios=[1]*inner_cols, height_ratios=[1]*inner_rows,
                                    hspace=hspace_inner, wspace=wspace_inner)
        for inner_row in range(inner_rows):
            for inner_col in range(inner_cols):
                ax = f.add_subplot(inner_grid[inner_row, inner_col])
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.axis('off')

                # Draw lenses
                if inner_row == 0 and inner_col == 1:
                    theta = np.linspace(np.pi/5,-np.pi/5+0.04,50)
                    x1 = np.cos(theta) - 0.5
                    y1 = np.sin(theta) + 0.5
                    l1 = lines.Line2D(x1, y1, c='k', lw=1, transform=ax.transAxes, ls='--')
                    f.lines.extend([l1])
                if inner_row == 0 and inner_col == 6:
                    theta = np.linspace(-np.pi/11 + np.pi,np.pi/11 + np.pi,50)
                    x1 = 2*np.cos(theta) + 0.5
                    y1 = 2*np.sin(theta) + 0.5
                    l1 = lines.Line2D(x1, y1, c='k', lw=1, transform=ax.transAxes, ls='--')
                    f.lines.extend([l1])

                # Draw labels
                def f_lines(ax):
                    ys = 0.1
                    l1 = lines.Line2D([0.15,0.5], [1.5-ys,1.5-ys], c='k', lw=1, transform=ax.transAxes)
                    l2 = lines.Line2D([0.5,0.5], [1.55-ys,1.45-ys], c='k', lw=1, transform=ax.transAxes)
                    l3 = lines.Line2D([-0.2,-0.5], [1.5-ys,1.5-ys], c='k', lw=1, transform=ax.transAxes)
                    l4 = lines.Line2D([-0.5,-0.5], [1.55-ys,1.45-ys], c='k', lw=1, transform=ax.transAxes)
                    f.lines.extend([l1, l2, l3, l4])
                def f_lines_long(ax):
                    ys = 0.1
                    l1 = lines.Line2D([-0.35,0.5], [1.5-ys,1.5-ys], c='k', lw=1, transform=ax.transAxes)
                    l2 = lines.Line2D([0.5,0.5], [1.55-ys,1.45-ys], c='k', lw=1, transform=ax.transAxes)
                    l3 = lines.Line2D([-0.7,-1.5], [1.5-ys,1.5-ys], c='k', lw=1, transform=ax.transAxes)
                    l4 = lines.Line2D([-1.5,-1.5], [1.55-ys,1.45-ys], c='k', lw=1, transform=ax.transAxes)
                    f.lines.extend([l1, l2, l3, l4])
                
                if outer_row == 0 and inner_row == 0:
                    col_labels = ['Object\n plane',
                                  'Objective\n reference\n sphere',
                                  'Pupil\n plane',
                                  '',
                                  'Tube lens\n reference\n sphere',
                                  '',
                                  'Detector\n plane']
                    ax.text(0.5,1.85,col_labels[inner_col], ha='center', va='center', transform=ax.transAxes)
                    if inner_col > 0 and inner_col < 3:
                        ax.text(0, 1.4, '$f_0$', ha='center', va='center', transform=ax.transAxes)
                        if inner_col == 1:
                            ax.text(0, 1.1, '$n_0$', ha='center', va='center', transform=ax.transAxes)
                        if inner_col == 2:
                            ax.text(0, 1.1, '$n_1$', ha='center', va='center', transform=ax.transAxes)
                        f_lines(ax)
                    if inner_col == 4 or inner_col == 6 :
                        ax.text(-0.5, 1.4, '$f_1$', ha='center', va='center', transform=ax.transAxes)
                        if inner_col == 4:
                            ax.text(-0.5, 1.1, '$n_1$', ha='center', va='center', transform=ax.transAxes)
                        if inner_col == 6:
                            ax.text(-0.5, 1.1, '$n_1$', ha='center', va='center', transform=ax.transAxes)

                        f_lines_long(ax)

                def brace_lines(ax):
                    wsi = wspace_inner
                    l1 = lines.Line2D([0,2+wsi], [-0.2,-0.2], c='k', lw=1, transform=ax.transAxes)
                    l2 = lines.Line2D([2+wsi,2+wsi], [-0.2,-0.1], c='k', lw=1, transform=ax.transAxes)
                    l3 = lines.Line2D([0,0], [-0.2,-0.1], c='k', lw=1, transform=ax.transAxes)
                    l4 = lines.Line2D([1+wsi/2,1+wsi/2], [-0.2,-0.3], c='k', lw=1, transform=ax.transAxes)
                    f.lines.extend([l1, l2, l3, l4])

                def my_axes(x0, y0, theta, strings, ax, dot=True):
                    ax.arrow(x0, y0, 0.4*np.cos(theta), 0.4*np.sin(theta), head_width=0.05, fc='k', transform=ax.transData, clip_on=False)
                    ax.arrow(x0, y0, 0.4*np.cos(theta + np.pi/2), 0.4*np.sin(theta + np.pi/2), head_width=0.05, fc='k', transform=ax.transData, clip_on=False)
                    ax.scatter([x0], [y0], edgecolors='k', facecolors='w', marker='o', lw=0.5, clip_on=False, s=[10], zorder=9, transform=ax.transData)
                    ax.text(x0 + 0.7*np.cos(theta), y0 + 0.7*np.sin(theta),'$'+strings[0]+'$', ha='center', va='center', transform=ax.transData)
                    ax.text(x0 + 0.3*np.cos(theta + np.pi), y0+0.3*np.sin(theta + np.pi),'$'+strings[1]+'$', ha='center', va='center', transform=ax.transData)
                    ax.text(x0 + 0.7*np.cos(theta + np.pi/2), y0 + 0.7*np.sin(theta + np.pi/2),'$'+strings[2]+'$', ha='center', va='center', transform=ax.transData)                                            
                    if dot:
                        ax.scatter([x0], [y0], c='k', marker='.', lw=0.5, clip_on=False, s=[6.5], zorder=10, transform=ax.transData)
                    else:
                        ax.scatter([-2], [0], c='k', marker='x', lw=0.5, clip_on=False, s=[6.5], zorder=10, transform=ax.transData)
                    
                if inner_col == 0:
                    if inner_row == 0:
                        ax.text(-2.75,0,'Axial view', ha='center', va='center', transform=ax.transData, rotation=90)
                        my_axes(-2,0,0,['\mathbf{\hat{z}}','\mathbf{\hat{y}}','\mathbf{\hat{x}}'], ax, dot=True)
                    #     rr = 1.7
                    #     theta = np.arctan(0.5)
                    #     x0 = rr*np.cos(theta)
                    #     y0 = rr*np.sin(theta)                        
                    #     ax.arrow(x0, y0, 0.4*np.cos(theta + np.pi/2), 0.4*np.sin(theta + np.pi/2), head_width=0.05, fc='k', transform=ax.transData, clip_on=False)
                    #     ax.scatter([x0], [y0], edgecolors='k', facecolors='w', marker='o', lw=0.5, clip_on=False, s=[10], zorder=9, transform=ax.transData)
                    #     ax.scatter([x0], [y0], c='k', marker='.', lw=0.5, clip_on=False, s=[6.5], zorder=10, transform=ax.transData)
                    #     ax.text(x0-0.4, y0+0.1, '$\hat{\\boldsymbol{\phi}}$', ha='center', va='center', transform=ax.transData)
                    #     ax.text(x0-0.4, y0+0.6, '$\hat{ \\boldsymbol{\\theta}}$', ha='center', va='center', transform=ax.transData)

                    #     x0 = 2.45
                    #     y0 = 1
                    #     theta = 0
                    #     ax.arrow(x0, y0, 0.4*np.cos(theta + np.pi/2), 0.4*np.sin(theta + np.pi/2), head_width=0.05, fc='k', transform=ax.transData, clip_on=False)
                    #     ax.scatter([x0], [y0], edgecolors='k', facecolors='w', marker='o', lw=0.5, clip_on=False, s=[10], zorder=9, transform=ax.transData)
                    #     ax.scatter([x0], [y0], c='k', marker='.', lw=0.5, clip_on=False, s=[6.5], zorder=10, transform=ax.transData)
                    #     ax.text(x0+0.3, y0, '$\hat{\\boldsymbol{\phi}}$', ha='center', va='center', transform=ax.transData)
                    #     ax.text(x0+0.3, y0+0.5, '$\hat{\\boldsymbol{\\rho}}$', ha='center', va='center', transform=ax.transData)
                        # my_axes(rr*np.cos(theta),rr*np.sin(theta),theta, ['\hat{\\theta}','\hat{\phi}','\mathbf{\hat{x}}'], ax, dot=False)
                    else:
                        ax.text(-2.75,0,'Transverse view', ha='center', va='center', transform=ax.transData, rotation=90)
                        my_axes(-2,0,0,['\mathbf{\hat{y}}','\mathbf{\hat{z}}','\mathbf{\hat{x}}'], ax, dot=False)

                # Draw dots
                if inner_col == 0 or inner_col == 2 or inner_col == 6:
                    ax.scatter([0.5], [0], c='k', marker='o', clip_on=False, s=[8], zorder=10, transform=ax.transAxes)
                    ax.scatter([0.5], [1], c='w', marker='o', edgecolor='k', lw=0.5, clip_on=False, s=[10], zorder=10, transform=ax.transAxes)

                # Draw coordinates
                if inner_col == 0:
                    ax.set_xlim([-1,1])
                    ax.set_ylim([-1,1])
                    if inner_row == 0:
                        ax.arrow(0, 0, 0, 0.4, head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0, 0.7,'$\mathbf{r}_o$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        ax.arrow(0, 0, 0.4*np.cos(-np.pi/4), 0.4*np.sin(-np.pi/4), head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0.75*np.cos(-np.pi/4),0.75*np.sin(-np.pi/4),'$\mathbf{\hat{s}}_o', ha='center', va='center', transform=ax.transData, clip_on=False)

                        # Draw circles
                        for y in range(-3, 3+1):                            
                            circle = mpatches.Circle((0,0.275*y), 0.1, ec="none")
                            patches.append(circle)
                        colors = np.linspace(0.85, 0.85, len(patches))
                        collection = PatchCollection(patches, cmap=plt.cm.gray, alpha=1)
                        collection.set_array(np.array(colors))
                        collection.set_clim([0,1])
                        ax.add_collection(collection)

                    if inner_row == 1:
                        ax.axis('on')
                        r_theta = np.pi/6
                        ax.arrow(0, 0, 0.4*np.cos(r_theta), 0.4*np.sin(r_theta), head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0.75*np.cos(r_theta), 0.75*np.sin(r_theta),'$\mathbf{r}_o$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        # theta = np.linspace(0,2*np.pi,100)
                        # ax.plot(np.cos(theta), np.sin(theta), '-k', lw=1, clip_on=False)

                        # Draw angles
                        l1 = lines.Line2D([0,0], [0,.95], c='k', lw=0.5, transform=ax.transData, zorder=10, ls=':')
                        tt = np.linspace(r_theta, np.pi/2, 20)
                        l3 = lines.Line2D(0.15*np.cos(tt), 0.15*np.sin(tt), c='k', lw=0.5, transform=ax.transData, zorder=3)
                        ax.text(0.5*np.cos(np.pi/3), 0.5*np.sin(np.pi/3),'$\\phi_o$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        f.lines.extend([l1, l3])
                        
                        # r_theta = 5*np.pi/8
                        # ax.arrow(0, 0, 0.4*np.cos(r_theta), 0.4*np.sin(r_theta), head_width=0.05, fc='k', transform=ax.transData)
                        # ax.text(0.75*np.cos(r_theta), 0.75*np.sin(r_theta),'$\mathbf{r}_o$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        s_theta = -np.pi/4
                        ax.arrow(0, 0, 0.4*np.cos(s_theta), 0.4*np.sin(s_theta), head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0.75*np.cos(s_theta),0.75*np.sin(s_theta),'$\mathbf{\hat{s}}_o', ha='center', va='center', transform=ax.transData, clip_on=False)

                        # Draw angles
                        l1 = lines.Line2D([0,0], [0,.95], c='k', lw=0.5, transform=ax.transData, zorder=10, ls=':')
                        tt = np.linspace(np.pi/2,s_theta,20)
                        l3 = lines.Line2D(0.1*np.cos(tt), 0.1*np.sin(tt), c='k', lw=0.5, transform=ax.transData, zorder=3)
                        ax.text(0.4*np.cos(0),0.4*np.sin(0),'$\\varphi$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        f.lines.extend([l3])

                        # Draw circles
                        for x in range(-3, 3+1):
                            for y in range(-3, 3+1):
                                circle = mpatches.Circle((0.275*x,0.275*y), 0.1, ec="none", zorder=10)
                                patches.append(circle)
                        colors = np.linspace(0.85, 0.85, len(patches))
                        collection = PatchCollection(patches, cmap=plt.cm.gray, alpha=1)
                        collection.set_array(np.array(colors))
                        collection.set_clim([0,1])
                        ax.add_collection(collection)
                        
                if inner_col == 2:
                    ax.set_xlim([-1,1])
                    ax.set_ylim([-1,1])
                    if inner_row == 0:
                        ax.arrow(0, 0, 0, 0.4, head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0, 0.7,'$\mathbf{r}_p$', ha='center', va='center', transform=ax.transData, clip_on=False)
                    if inner_row == 1:
                        r_theta = np.pi/6
                        ax.arrow(0, 0, 0.4*np.cos(r_theta), 0.4*np.sin(r_theta), head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0.75*np.cos(r_theta), 0.75*np.sin(r_theta),'$\mathbf{r}_p$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        theta = np.linspace(0,2*np.pi,100)
                        ax.plot(np.cos(theta), np.sin(theta), '-k', lw=1, clip_on=False)

                        # Draw angles
                        l1 = lines.Line2D([0,0], [0,.95], c='k', lw=0.5, transform=ax.transData, zorder=10, ls=':')
                        tt = np.linspace(r_theta, np.pi/2, 20)
                        l3 = lines.Line2D(0.15*np.cos(tt), 0.15*np.sin(tt), c='k', lw=0.5, transform=ax.transData, zorder=3)
                        ax.text(0.5*np.cos(np.pi/3), 0.5*np.sin(np.pi/3),'$\\phi_p$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        f.lines.extend([l1, l3])

                if inner_col == 6:
                    ax.set_xlim([-1,1])
                    ax.set_ylim([-1,1])
                    if inner_row == 0:
                        ax.arrow(0, 0, 0, 0.4, head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0, 0.7,'$\mathbf{r}_d$', ha='center', va='center', transform=ax.transData, clip_on=False)
                    if inner_row == 1:
                        ax.axis('on')

                        ax.arrow(0, 0, 0.4*np.cos(np.pi/6), 0.4*np.sin(np.pi/6), head_width=0.05, fc='k', transform=ax.transData)
                        ax.text(0.75*np.cos(np.pi/6), 0.75*np.sin(np.pi/6),'$\mathbf{r}_d$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        
                        # Draw angles
                        l1 = lines.Line2D([0,0], [0,.95], c='k', lw=0.5, transform=ax.transData, zorder=10, ls=':')
                        tt = np.linspace(r_theta, np.pi/2, 20)
                        l3 = lines.Line2D(0.15*np.cos(tt), 0.15*np.sin(tt), c='k', lw=0.5, transform=ax.transData, zorder=3)
                        ax.text(0.5*np.cos(np.pi/3), 0.5*np.sin(np.pi/3),'$\\phi_d$', ha='center', va='center', transform=ax.transData, clip_on=False)
                        f.lines.extend([l1, l3])


                if inner_row == 0 and inner_col == 0:
                    # Draw alpha
                    # l1 = lines.Line2D([0,2,2,0], [0,0,1,0], c='k', lw=0.5, transform=ax.transData, zorder=3, ls=':')
                    tt = np.linspace(0, np.pi/5)
                    l4 = lines.Line2D(0.5*np.cos(tt), 0.5*np.sin(tt), c='k', lw=0.5, transform=ax.transData, zorder=3)
                    ax.text(0.8*np.cos(np.pi/10), 0.8*np.sin(np.pi/10),'$\\alpha$', ha='center', va='center', transform=ax.transData)
                    # ax.text(2.3,0.5,'$\\sin\\alpha \\approx \\alpha$', ha='left', va='center', transform=ax.transData)
                    f.lines.extend([l4])

                    # Draw vartheta
                    # l1 = lines.Line2D([2,0,2,0], [0,0,1,0], c='k', lw=0.5, transform=ax.transData, zorder=3, ls=':')
                    l1 = lines.Line2D([2,0,2*np.cos(np.pi/5)], [0,0,2*np.sin(np.pi/5)], c='k', lw=0.5, transform=ax.transData, zorder=3, ls=':')                    
                    tt = np.linspace(-np.pi/4,0,20)
                    l4 = lines.Line2D(0.2*np.cos(tt), 0.2*np.sin(tt), c='k', lw=0.5, transform=ax.transData, zorder=3)
                    ax.text(0.7*np.cos(-np.pi/8), 0.7*np.sin(-np.pi/8),'$\\vartheta$', ha='center', va='center', transform=ax.transData)
                    f.lines.extend([l1, l4])
                    
f.savefig('detection-coords.pdf', dpi=500, bbox_inches='tight')
