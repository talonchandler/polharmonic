import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Collect results
folder = 'asym_dispim'
ims = [ folder+'/data_ortho1.png', folder+'/data_ortho2.png', folder+'/data_both.png',]
fig, axs = plt.subplots(1, 3, figsize=(12, 4), gridspec_kw={'hspace':0, 'wspace':0})
for i, im in enumerate(ims):
    image = mpimg.imread(im)
    axs[i].imshow(image, interpolation=None, vmin=0, vmax=1)
    axs[i].set_axis_off()
axs[0].annotate('1.1 NA data only', xy=(0,0), xytext=(0.5, 1.0), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
axs[1].annotate('0.71 NA data only', xy=(0,0), xytext=(0.5, 1.0), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
axs[2].annotate('All data', xy=(0,0), xytext=(0.5, 1.0), textcoords='axes fraction', va='center', ha='center', fontsize=12, annotation_clip=False)
fig.savefig(folder+'/recons.pdf', dpi=800, bbox_inches='tight')
