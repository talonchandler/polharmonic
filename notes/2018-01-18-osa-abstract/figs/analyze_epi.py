import numpy as np
import dill
from polharmonic import util, multi, dist, sft
import os; import time; start = time.time(); print('Running...')
np.set_printoptions(precision=3, suppress=True)

# Specify saving path
folder = 'epi'
if not os.path.exists(folder):
    os.makedirs(folder)
data_file = folder+'/'+folder+'.dat'

# Build microscope
exp = multi.MultiMicroscope(max_l=4, n_pts=5000)
#exp.calc_sys_matrix()

# Save
#dill.dump(exp.psi, open(data_file,'wb'))

# Load 
f = open(data_file, 'rb')
exp.psi = dill.load(f)

# Calculate B
exp.calc_B_matrix()

# # Test the dist.make_positive() function
# d = dist.Distribution(sh=np.array([0,1,1,0,0,0,0,0,0,0,0,0,0,0,0]))
# d.plot_dist(exp.B, exp.xyz, exp.triangles, filename='before.png')
# d.make_positive(exp.B)
# d.plot_dist(exp.B, exp.xyz, exp.triangles, filename='after.png')

# # Find intensities
# g = exp.calc_intensities(d)

# # Try reconstruction
# Fstar = exp.recon_dist(g)
# d2 = dist.Distribution(sh=Fstar)
# d2.plot_dist(exp.B, exp.xyz, exp.triangles, filename='recon.png')

# Test field of distributions
# Generate truth orientations
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
#df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename='field_before.png', show=True)
df.make_positive(exp.B)
#df.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename='field_after.png', show=True)

g = exp.calc_intensity_field(df)
df2 = exp.recon_dist_field(g)
df2.plot_dist_field(exp.B, exp.xyz, exp.triangles, filename='field_recon.png', show=True)

# Plot forward model matrix
#exp.plot_matrix(folder, skip_sch=True)

# Plot svs
# exp.plot_svs(folder)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
