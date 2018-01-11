import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mayavi import mlab
from scipy.special import sph_harm
from skimage import io

# Convert a list of coefficients into a system matrix with labels
def list2matrix(arr, data_i='Hlm', row_i='i', col_i='j', lab_col_i='lm', lab_row_i='pol'):
    out_arr = np.zeros((max(arr[row_i])+1, max(arr[col_i])+1), dtype='complex')
    col_labels = np.zeros(max(arr[col_i]+1), dtype='U4')
    row_labels = np.zeros(max(arr[row_i]+1), dtype='U10')
    for entry in arr:
        out_arr[entry[row_i], entry[col_i]] = entry[data_i]
        col_labels[entry[col_i]] = str(entry[lab_col_i[0]]) + ',' + str(entry[lab_col_i[1]])
        row_labels[entry[row_i]] = entry[lab_row_i]
    return np.real(out_arr), col_labels, row_labels

def spherical_mesh(coeffs, labels, r=1, gridn=151):
    # Create sampling grid
    theta, phi = np.mgrid[0:np.pi:(gridn*1j), 0:2*np.pi:(gridn*1j)]
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    # Calculate values
    rc = 0
    for coeff, label in zip(coeffs, labels):
        l = int(label.split(',')[0])
        m = int(label.split(',')[1])
        rc += coeff*spZnm(l, m, theta, phi)
        
    # Check that imaginary part is small before discarding it
    if np.max(np.imag(rc)) > 1e-3:
        import pdb; pdb.set_trace()
    rc = np.real(rc)

    # Separate positive and negative components
    n = rc.clip(max=0) 
    p = rc.clip(min=0)*(-1)

    return x, y, z, p, n

# Plot an orientation distribution from its spherical harmonic coefficients
def plot_spherical(coeffs, labels, filename='spherical.png', 
                   ax=None, r=1, mag=1, show=False):
    # Create figure
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
    mlab.clf()

    x, y, z, p, n = spherical_mesh(coeffs, labels, r=1, gridn=50)
    
    # Plot
    mlab.mesh(p*x, p*y, p*z, color=(1, 0, 0), representation='surface')
    mlab.mesh(n*x, n*y, n*z, color=(0, 0, 1), representation='surface')

    # View and save
    mlab.view(azimuth=45, elevation=45, distance=3, focalpoint=None,
              roll=None, reset_roll=True, figure=None)
    mlab.savefig(filename, magnification=mag)
    subprocess.call(['convert', filename, '-transparent', 'white', filename])
    if show:
        mlab.show()

# Plot a field of orientation distributions
def plot_spherical_field(coeffs_matrix, labels, filename='spherical.png',
                         ax=None, r=1, mag=1, dist=10, show=False, gridn=100):
    # Create figure
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 400))
    mlab.clf()
    
    spacing = 2
    shift = spacing*gridn/2
    for i in np.ndindex(coeffs_matrix.shape[:2]):
        print(i)
        x, y, z, p, n = spherical_mesh(coeffs_matrix[i], labels, r=r, gridn=gridn)
        mlab.mesh(p*x + spacing*i[0] - shift, p*y + spacing*i[1] - shift, p*z, color=(1, 0, 0), representation='surface')
        mlab.mesh(n*x + spacing*i[0] - shift, n*y + spacing*i[1] - shift, n*z, color=(0, 0, 1), representation='surface')

    # View and save
    mlab.view(azimuth=45, elevation=45, distance=dist, focalpoint=None,
              roll=None, reset_roll=True, figure=None)
    mlab.savefig(filename, magnification=mag)
    subprocess.call(['convert', filename, '-transparent', 'white', filename])
    if show:
        mlab.show()
        
# SciPy real spherical harmonics with identical interface to SymPy's Znm
# Useful for faster numerical evaluation of Znm
def spZnm(l, m, theta, phi):
    if m > 0:
        return (sph_harm(m, l, phi, theta) +
                np.conj(sph_harm(m, l, phi, theta)))/(np.sqrt(2))
    elif m == 0:
        return sph_harm(m, l, phi, theta)
    elif m < 0:
        return  (sph_harm(m, l, phi, theta) -
                 np.conj(sph_harm(m, l, phi, theta)))/(np.sqrt(2)*1j)

def draw_scene(scene_string, filename='out.png', my_ax=None, dpi=300,
               save_file=False, chop=True):
    asy_string = """
    import three;
    import graph3;
    settings.outformat = "pdf";
    settings.prc = true;
    settings.embed= true;
    settings.render=16;

    size(6cm,6cm);
    currentprojection = orthographic(1, 1, 1);

    void circle(real Theta, real Alpha, bool dash, triple color) {
      triple normal = expi(Theta, 0);
      real h = 1 - sqrt(2 - 2*cos(Alpha) - sin(Alpha)^2);
      real radius = sin(Alpha);
      path3 mycircle = circle(c=h*normal, r=radius, normal=normal);
      if (dash) {
        draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
      } else {
        draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
      }
    }

    void ellipse(real Theta, real Phi, real a, real b, real theta, bool dash, triple color) {
      triple normal = expi(Theta, Phi);
      real a_scaled = a/max(a, b);
      real b_scaled = b/max(a, b);      
      path3 mycircle = rotate(degrees(Phi), Z)*rotate(degrees(Theta), Y)*shift(Z)*rotate(degrees(theta), Z)*scale(a_scaled, b_scaled, 1)*circle(c=O, r=0.05, normal=Z);
      if (dash) {
        draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
      } else {
        draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
      }
    }

    void mydot(real Theta, triple color) {
      triple normal = expi(Theta, 0);
      dot(normal, p=rgb(xpart(color), ypart(color), zpart(color)));
    }

    void arrow(real Theta, real Phi_Pol, triple color, bool dash) {
      if (dash) {
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      } else {
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      }
    }

    void watson(real Theta, real Phi, real kappa, real x, real y, real z) {
     int n_phi = 10;
     int n_theta = 10;

     real max_radius = 0;
     if(kappa > 0){
       max_radius = exp(kappa);
     }
     else{
       max_radius = 1.0;
     }

     for(int i=0; i <= n_theta; ++i) {
       real Theta_i = pi*i/n_theta;
       real weight = exp(kappa*(cos(Theta_i)**2))/max_radius;     
       path3 mycircle = circle(c=Z*weight*cos(Theta_i), r=weight*sin(Theta_i));
       draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*mycircle);
     }

     triple f(real t) {
       real weight = exp(kappa*(cos(t)**2))/max_radius;
       return (0, weight*sin(t), weight*cos(t));
     }
     path3 phi_path = graph(f, 0, 2pi, operator ..);

     for(int i=0; i <= n_phi; ++i) {
       real Phi_i = 2*pi*i/n_theta;
       draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*rotate(angle=degrees(Phi_i), u=(0,0,0), v=(0,0,1))*phi_path);
     }
    }
    real len = 10;
    draw((-len,-len)--(len,-len)--(len,len)--(-len,len)--(-len,-len), white);

    draw(unitsphere, surfacepen=material(diffusepen=white+opacity(0.1), emissivepen=grey, specularpen=white));

    // Draw points on sphere
    dotfactor = 7;
    dot(X); 
    dot(Y); 

    circle(0, pi/2, false, (0, 0, 0));
    """

    asy_string += scene_string
    asy_string += "dot(Z);shipout(scale(4.0)*currentpicture.fit());"

    text_file = open("temp.asy", "w")
    text_file.write(asy_string)
    text_file.close()

    subprocess.call(['asy', 'temp.asy'])
    subprocess.call(['convert', '-density', str(dpi), '-units', 'PixelsPerInch', 'temp.pdf', 'temp.png'])
    im = mpimg.imread('temp.png')

    # Chop top of im to make it square and fix asy error
    if chop:
        im = im[int(im.shape[1]*0.075):,:,:]
    
    f = plt.figure(figsize=(5, 5), frameon=False)
    local_ax = plt.axes([0, 0, 1, 1]) # x, y, width, height
    if my_ax == None:
        my_ax = local_ax

    for ax in [local_ax, my_ax]:
        #draw_axis(ax)
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

        # Plot
        ax.imshow(im, interpolation='none')

    # Save
    if save_file:
        f.savefig(filename, dpi=dpi)

    subprocess.call(['rm', 'temp.asy', 'temp.pdf', 'temp.png'])
    return ax

# Plots a set of spherical harmonics and returns a string for tex 
def plot_multiple_spherical(coeffs, labels, filename='sph', r=1.0):
    sph_tex_string = ''
    for i in range(len(labels)):
        plot_spherical(coeffs[:,i], labels, filename=filename+str(i)+'.png', r=r)
        tex_string = "\pic{1.0}{sphXXX}\\\\"
        sph_tex_string += tex_string.replace('XXX', str(i))
    return sph_tex_string

# Create matrix in latex and save
def create_latex_matrix(sph_tex_string, asy_tex_string, psi, filename):
    tex_template = r"""
    \documentclass[preview, border={0pt 10pt 0pt 0pt}]{standalone}
    \usepackage{amsmath}
    \usepackage{graphicx}
    \setcounter{MaxMatrixCols}{20}
    \def\xpic#1#2{\includegraphics[trim=0 5cm 0 5cm, width=#1em]{#2}}
    \def\xpicbig#1#2{\includegraphics[trim=0 0 0 0, width=#1em]{#2}}
    \def\pic#1#2{{%
      \mathchoice
        {\xpic{#1}{#2}}%
        {\xpic{#1}{#2}}%
        {\xpic{\defaultscriptratio}{#2}}%
        {\xpic{\defaultscriptscriptratio}{#2}}}}
    \def\picbig#1#2{{% mathord
      \mathchoice
        {\xpicbig{#1}{#2}}%
        {\xpicbig{#1}{#2}}%
        {\xpicbig{\defaultscriptratio}{#2}}%
        {\xpicbig{\defaultscriptscriptratio}{#2}}}}

    \begin{document}
    \begin{align*}
      \begin{bmatrix}
         ASY_STRING
      \end{bmatrix} =
      \begin{bmatrix}
         ARR_STRING
      \end{bmatrix}
      \begin{bmatrix}
        SPH_STRING   
      \end{bmatrix}
    \end{align*}
    \end{document}
    """
    tex_template = tex_template.replace('ASY_STRING', asy_tex_string)
    tex_template = tex_template.replace('SPH_STRING', sph_tex_string)

    # Convert psi to tex
    np.savetxt('psi.csv', psi.real, delimiter=' & ', fmt='%2.2f', newline=' \\\\')
    with open ("psi.csv", "r") as myfile:
        arr_tex_string = myfile.readlines()
    tex_template = tex_template.replace('ARR_STRING', ' '.join(arr_tex_string))
    
    with open(filename+'.tex', 'w') as text_file:
        print(tex_template, file=text_file)
    subprocess.call(['latexmk', '-cd', filename+'.tex'])

def tiff2array(filename, x=0, y=0, z=0, width=None, height=None, slices=None):
    # Opens .tiff and returns array starting at x, y, z with width, height, and
    # slices dimensions. "None" means return the whole dimension.
    im = io.imread(filename)
    shape = im.shape
    x_min = x
    if width is None:
        x_max = shape[2]
    else:
        x_max = x + width
    y_min = y
    if height is None:
        y_max = shape[1]
    else:
        y_max = y + height
    z_min = z
    if slices is None:
        z_max = shape[0]
    else:
        z_max = z + slices
    im = im[z_min:z_max, y_min:y_max, x_min:x_max]
    if slices == 1:
        return im.squeeze()
    else:
        return im

