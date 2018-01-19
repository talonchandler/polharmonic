import numpy as np
import dill
from polharmonic import util, multi, dist, sft
import os; import time; start = time.time(); print('Running...')
np.set_printoptions(precision=3, suppress=True)

# Specify saving path
folder = 'ortho'
if not os.path.exists(folder):
    os.makedirs(folder)
data_file = folder+'/'+folder+'.dat'

# Build microscope
exp = multi.MultiMicroscope(ill_thetas=[90], det_nas=[1.1], max_l=4, n_pts=250)

# Calculate and save (comment this on repeat runs)
exp.calc_sys_matrix()
dill.dump(exp.psi, open(data_file,'wb'))

import pdb; pdb.set_trace() 
# Load 
f = open(data_file, 'rb')
exp.psi = dill.load(f)

# Calculate B
exp.calc_B_matrix()

# Generate phantom of true orientations
nx, ny = 10, 10
t = np.linspace(0, np.pi/2, nx)
p = np.linspace(0, np.pi, ny)
tv, pv = np.meshgrid(t, p)
tvv = np.expand_dims(tv, 2)
pvv = np.expand_dims(pv, 2)
tp_field = np.concatenate([tvv, pvv], 2)

# Expand onto spherical harmonics
sh_field = sft.field_tp_sft(tp_field, max_l=4)
df = dist.DistributionField(sh_arr=sh_field)
df.make_positive(exp.B)
df.sh_arr[:,:,6:] = 0 # Low pass to simplify phantom
df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/phantom.png',
                   mag=5, show=True)

# Reconstruct phantom
g = exp.calc_intensity_field(df)
df2 = exp.recon_dist_field(g)
df2.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/phantom_recon.png',
                    mag=5, show=True)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
