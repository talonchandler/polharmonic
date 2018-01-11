import numpy as np
from util import *
from sympy import *
from sympy.matrices.dense import *
import sympy.functions.special.spherical_harmonics as sh
from sft import *
import dill

# Global initialize symbols
theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)
Theta = Symbol('Theta', real=True)
Phi = Symbol('Phi', real=True)

# Helper functions
def NAn2AB(NA, n):
    alpha = np.arcsin(NA/n)
    return alpha2A(alpha), alpha2B(alpha)

def alpha2A(alpha):
    ca = np.cos(alpha)
    return 0.25 - 0.375*ca + 0.125*(ca**3)

def alpha2B(alpha):
    ca = np.cos(alpha)
    return 0.1875*ca + 0.1875*(ca**3)

# Analytic Hlm functions for specific acquisitions 
# z light-sheet illumination, x widefield detection
def illzdetx(phi_pol, NA, n):
    A, B = NAn2AB(NA, n)    
    prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*(1 - (sin(theta)**2)*(cos(phi)**2)))
    return sft(prf)

# x light-sheet illumination, z widefield detection
def illxdetz(phi_pol, NA, n):
    A, B = NAn2AB(NA, n)    
    prf = ((sin(theta)*sin(phi)*sin(phi_pol) - cos(theta)*cos(phi_pol))**2)*2*(A + B*(sin(theta)**2))
    return sft(prf)

# z light-sheet illumination (small NA), z widefield detection
def illzdetz(phi_pol, NA, n):
    A, B = NAn2AB(NA, n)
    prf = (sin(theta)**2)*(cos(phi - phi_pol)**2)*2*(A + B*sin(theta)**2)
    return sft(prf)

# uniform illumination (all directions), single detector along (Theta, Phi)
def illalldetpixel(Theta, Phi):
    prf = 1 - (sin(Theta)*cos(Phi)*sin(theta)*cos(phi) + sin(Theta)*sin(Phi)*sin(theta)*sin(phi) + cos(Theta)*cos(theta))**2
    return sft(prf)

# Complete system matrix functions
# Use matrix_file to load in a file and skip building the matrix
def dispim_sys_matrix(NAx, NAz, n, phi_pols):
    sys_matrix = []
    phi_pol = Symbol('phi_pol', real=True)
    
    ixdz = illxdetz(phi_pol, NAz, n)
    izdx = illzdetx(phi_pol, NAx, n)

    # Populate system matrix
    for i, pol in enumerate(phi_pols):
        for j, (Hlm, l, m) in enumerate(izdx):
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i, j, l, m, 'izdx,'+str(np.rad2deg(pol))))
    for i, pol in enumerate(phi_pols):            
        for j, (Hlm, l, m) in enumerate(ixdz):            
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i+len(phi_pols), j, l, m, 'ixdz,'+str(np.rad2deg(pol))))

    dtypes = [('Hlm', 'c16'), ('i', 'i4'), ('j', 'i4'), ('l', 'i4'), ('m', 'i4'), ('pol', 'U10')]
    
    arr = np.array(sys_matrix, dtype=dtypes)
    arr, col_labels, row_labels = list2matrix(np.real(arr))
    
    # Remove zero-filled columns
    zero_mask = np.all(np.abs(arr) < 1e-10, axis=0)
    arr = arr[:, ~zero_mask]
    col_labels = col_labels[~zero_mask]

    # Create schematic string
    sch_strings = []
    for i, phi_pol in enumerate(phi_pols):
        asy_string = ''
        illum_string = "mydot(theta, color);\n"
        illum_string = illum_string.replace('theta', str(0))
        illum_string = illum_string.replace('color', str((1,0,0)))
        asy_string += illum_string

        pol_string = "arrow(theta, phi_pol, color, false);\n"
        pol_string = pol_string.replace('theta', str(0))
        pol_string = pol_string.replace('phi_pol', str(np.rad2deg(phi_pol)))
        pol_string = pol_string.replace('color', str((1,0,0)))
        asy_string += pol_string

        detect_string = "circle(theta, alpha, true, color);\n"
        detect_string = detect_string.replace('theta', '1.57')
        detect_string = detect_string.replace('alpha', str(np.arcsin(NAx/n)))
        detect_string = detect_string.replace('color', str((1,0,0)))
        asy_string += detect_string

        sch_strings.append(asy_string)
        
    for i, phi_pol in enumerate(phi_pols):
        asy_string = ''
        illum_string = "mydot(theta, color);\n"
        illum_string = illum_string.replace('theta', '1.57')
        illum_string = illum_string.replace('color', str((0,0,1)))
        asy_string += illum_string

        pol_string = "arrow(theta, phi_pol, color, false);\n"
        pol_string = pol_string.replace('theta', '90')
        pol_string = pol_string.replace('phi_pol', str(np.rad2deg(phi_pol)))
        pol_string = pol_string.replace('color', str((0,0,1)))
        asy_string += pol_string

        detect_string = "circle(theta, alpha, true, color);\n"
        detect_string = detect_string.replace('theta', str(0))
        detect_string = detect_string.replace('alpha', str(np.arcsin(NAz/n)))
        detect_string = detect_string.replace('color', str((0,0,1)))
        asy_string += detect_string

        sch_strings.append(asy_string)        

    return arr, col_labels, row_labels, sch_strings

def epi_sys_matrix(NA, n, phi_pols):
    sys_matrix = []
    phi_pol = Symbol('phi_pol', real=True)

    izdz = illzdetz(phi_pol, NA, n)

    # Populate system matrix
    for i, pol in enumerate(phi_pols):
        for j, (Hlm, l, m) in enumerate(izdz):
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i, j, l, m, str(np.rad2deg(pol))))

    dtypes = [('Hlm', 'c16'), ('i', 'i4'), ('j', 'i4'), ('l', 'i4'), ('m', 'i4'), ('pol', 'U10')]

    arr = np.array(sys_matrix, dtype=dtypes)
    arr, col_labels, row_labels = list2matrix(np.real(arr))

    # Remove zero-filled columns
    zero_mask = np.all(np.abs(arr) < 1e-10, axis=0)
    arr = arr[:, ~zero_mask]
    col_labels = col_labels[~zero_mask]

    # Create schematic string
    sch_strings = []
    for i, phi_pol in enumerate(phi_pols):
        asy_string = ''
        illum_string = "mydot(theta, color);\n"
        illum_string = illum_string.replace('theta', str(0))
        illum_string = illum_string.replace('color', str((1,0,0)))
        asy_string += illum_string

        pol_string = "arrow(theta, phi_pol, color, false);\n"
        pol_string = pol_string.replace('theta', str(0))
        pol_string = pol_string.replace('phi_pol', str(np.rad2deg(phi_pol)))
        pol_string = pol_string.replace('color', str((1,0,0)))
        asy_string += pol_string

        detect_string = "circle(theta, alpha, true, color);\n"
        detect_string = detect_string.replace('theta', str(0))
        detect_string = detect_string.replace('alpha', str(np.arcsin(NA/n)))
        detect_string = detect_string.replace('color', str((1,0,0)))
        asy_string += detect_string

        sch_strings.append(asy_string)        

    return arr, col_labels, row_labels, sch_strings

def ortho_sys_matrix(NA, n, phi_pols):
    sys_matrix = []
    phi_pol = Symbol('phi_pol', real=True)

    izdz = illxdetz(phi_pol, NA, n)

    # Populate system matrix
    for i, pol in enumerate(phi_pols):
        for j, (Hlm, l, m) in enumerate(izdz):
            sys_matrix.append((Hlm.evalf(subs={phi_pol: pol}), i, j, l, m, str(np.rad2deg(pol))))

    dtypes = [('Hlm', 'c16'), ('i', 'i4'), ('j', 'i4'), ('l', 'i4'), ('m', 'i4'), ('pol', 'U10')]

    arr = np.array(sys_matrix, dtype=dtypes)
    arr, col_labels, row_labels = list2matrix(np.real(arr))

    # Remove zero-filled columns
    zero_mask = np.all(np.abs(arr) < 1e-10, axis=0)
    arr = arr[:, ~zero_mask]
    col_labels = col_labels[~zero_mask]

    # Create schematic string
    sch_strings = []
    for i, phi_pol in enumerate(phi_pols):
        asy_string = ''
        illum_string = "mydot(theta, color);\n"
        illum_string = illum_string.replace('theta', '1.57')
        illum_string = illum_string.replace('color', str((1,0,0)))
        asy_string += illum_string

        pol_string = "arrow(theta, phi_pol, color, false);\n"
        pol_string = pol_string.replace('theta', '90')
        pol_string = pol_string.replace('phi_pol', str(np.rad2deg(phi_pol)))
        pol_string = pol_string.replace('color', str((1,0,0)))
        asy_string += pol_string
        
        detect_string = "circle(theta, alpha, true, color);\n"
        detect_string = detect_string.replace('theta', str(0))
        detect_string = detect_string.replace('alpha', str(np.arcsin(NA/n)))
        detect_string = detect_string.replace('color', str((1,0,0)))
        asy_string += detect_string

        sch_strings.append(asy_string)        

    return arr, col_labels, row_labels, sch_strings


# Do all plotting 
def generate_plots(psi, col_labels, row_labels, sch_strings, folder):
    print("Plotting non-zero spherical harmonics...")
    sph_tex_string = plot_multiple_spherical(np.identity(len(col_labels)), col_labels, filename=folder+'/sph')
    
    print("Plotting microscope schematics...") 
    asy_tex_string = ''
    for i, sch_string in enumerate(sch_strings):
        # draw_scene(sch_string, filename=folder+'/sch'+str(i)+'.png', dpi=300,
        #            save_file=True, chop=True)
        tex_string = "\picbig{2.0}{schXXX}\\\\"
        asy_tex_string += tex_string.replace('XXX', str(i))
        
    print("Plotting latex matrix...")
    create_latex_matrix(sph_tex_string, asy_tex_string, psi, folder+'/matrix')
    
    # SVD
    U, s, Vh = np.linalg.svd(psi.real)
    V = Vh.conj().T
    
    # Plot individual singular distributions
    plot_multiple_spherical(V, col_labels, filename=folder+'/sv', r=0.8)

    # Plot singular value spectrum TODO
