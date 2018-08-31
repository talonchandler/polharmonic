from polaris import data, spang, phantom
from polaris.micro import multi
import numpy as np

# Make output folder
folder = './helix-noise/'
import os
if not os.path.exists(folder):
    os.makedirs(folder)

# Generate phantom
vox_dim = (130,130,130) # nm
# helix phantom
px = (64,64,64)
phant = phantom.three_helix(vox_dim=vox_dim, px=px) 
# bead phantom - try varying orientation and kappa
# px = (32,32,32)
# phant = phantom.bead(orientation=[1,0,0], kappa=30, vox_dim=vox_dim, px=px)

# Calculate phantom statistics and save
phant.calc_stats()
# masking is necessary for fast rendering
mask = phant.density > 0.1 
# try "interact=True" to interact with the phantom
# phant.visualize(folder+'phantom/', mask=mask, interact=False, video=True,
#                 n_frames=15) 
phant.save_stats(folder+'phantom/')
phant.save_summary(folder+'phantom.pdf', mask=mask)

# Specify microscope
# try "det_nas = [0.8, 0.8]" for symmetric diSPIM
data1 = data.Data(g=np.zeros(phant.f.shape[0:3]+(4,2)), vox_dim=vox_dim,
                  det_nas=[1.1, 0.71])
m = multi.MultiMicroscope(phant, data1, n_samp=1.33, lamb=525,
                          spang_coupling=True)

# Calculate system matrix
m.calc_H()

# Generate data using forward model
# set "snr" to a positive number to simulate Poisson noise
data1.g = m.fwd(phant.f, snr=10) 

# Save data
data1.save_mips(folder+'data.pdf')
# try "diSPIM_format=False" for a 5D ImageJ hyperstack (pol = c, view = t)
data1.save_tiff(folder+'data/', diSPIM_format=True) 

# Calculate pseudoinverse solution
# set "eta" to a positive number for Tikhonov regularization
etas = 10**(np.linspace(-5, 0, 10))
errors = []
for i, eta in enumerate(etas):
    phant1 = spang.Spang(np.zeros(phant.f.shape))
    phant1.f = m.pinv(data1.g, eta=eta)

    # Calculate reconstruction statistics and save
    phant1.calc_stats()
    mask = phant.density > 0.1
    # phant.visualize(folder+'phantom-recon/', mask=mask, interact=False, video=True,
    #                 n_frames=15)
    
    phant1.save_stats(folder+'phantom-recon-'+str(i)+'/')
    phant1.save_summary(folder+'phantom-recon-'+str(i)+'.pdf', mask=mask)

    error = np.linalg.norm((phant1.f - phant.f).flatten())
    print(error)
    errors.append(error)

import matplotlib.pyplot as plt
f, ax = plt.subplots(1, 1, figsize=(4, 4))
ax.plot(etas, errors, '.k-', lw=0.5)
ax.set_xscale('log')
ax.set_xlabel('Tikhonov regularization parameter $\eta$')
ax.set_ylabel(r'Mean square error = $\frac{1}{N}\sum_{\texttt{j},\texttt{r}_o}^N (\hat{\theta}_{\texttt{j},\texttt{r}_o} - \theta_{\texttt{j},\texttt{r}_o})^2$')
f.savefig(folder+'mse.pdf', bbox_inches='tight')
