import numpy as np
import dill
from polharmonic import util, multi, dist, sft
import os; import time; start = time.time(); print('Running...')
np.set_printoptions(precision=3, suppress=True)

# Specify saving path
folder = 'recon-pipeline'
if not os.path.exists(folder):
    os.makedirs(folder)
data_file = 'asym_dispim/asym_dispim.dat'

# Build microscope
exp = multi.MultiMicroscope(ill_thetas=[90, 0], det_thetas=[0, 90],
                            det_nas=[1.1, 0.8], max_l=4, n_pts=5000)

# Calculate and save (comment this on repeat runs)
# exp.calc_sys_matrix()
# dill.dump(exp.psi, open(data_file,'wb'))

# Load 
f = open(data_file, 'rb')
exp.psi = dill.load(f)

# Calculate B
exp.calc_B_matrix()

# Generate phantom of true orientations
nx, ny = 12, 12
t = np.linspace(0, np.pi/2, nx)
p = np.linspace(0, 2*np.pi, ny)
tv, pv = np.meshgrid(t, p)
tvv = np.expand_dims(tv, 2)
pvv = np.expand_dims(pv, 2)
tp_field = np.concatenate([tvv, pvv], 2)

# Expand onto spherical harmonics
import pdb; pdb.set_trace()
tp_known_field = util.fibonacci_sphere(exp.n_pts)
sh_field = sft.field_tp_sft(tp_field, max_l=4)
df = dist.DistributionField(sh_arr=sh_field)
df.make_positive(exp.B)
df.sh_arr[:,:,6:] = 0 # Low pass to simplify phantom
df.sh_arr = np.expand_dims(df.sh_arr, 2)
df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/phantom.png',
                   r=0.7, mag=2, show=False)

# Reconstruct phantom
intf = exp.calc_intensity_field(df)
import pdb; pdb.set_trace() 
df2 = exp.recon_dist_field(intf)
df2.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename=folder+'/phantom_recon.png',
                    r=0.7, mag=2, show=True)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
