import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# import matplotlib
# from mayavi import mlab
from scipy.special import sph_harm
# import scipy.misc
# from skimage import io
# import vispy
# from polharmonic import visuals

# SciPy real spherical harmonics with identical interface to SymPy's Znm
# Useful for faster numerical evaluation of Znm
def spZnm(l, m, theta, phi):
    if m > 0:
        return np.real((sph_harm(m, l, phi, theta) +
                np.conj(sph_harm(m, l, phi, theta)))/(np.sqrt(2)))
    elif m == 0:
        return np.real(sph_harm(m, l, phi, theta))
    elif m < 0:
        return  np.real((sph_harm(m, l, phi, theta) -
                 np.conj(sph_harm(m, l, phi, theta)))/(np.sqrt(2)*1j))

# # Create matrix in latex and save 
# def create_latex_matrix(lhs_string, rhs_string, array, filename):
#     tex_template = r"""
#     \documentclass[preview, border={0pt 10pt 0pt 0pt}]{standalone}
#     \usepackage{amsmath}
#     \usepackage{graphicx}
#     \setcounter{MaxMatrixCols}{20}
#     \def\xpic#1#2{\includegraphics[trim=0 5cm 0 5cm, width=#1em]{#2}}
#     \def\xpicbig#1#2{\includegraphics[trim=0 0 0 0, width=#1em]{#2}}
#     \def\pic#1#2{{%
#       \mathchoice
#         {\xpic{#1}{#2}}%
#         {\xpic{#1}{#2}}%
#         {\xpic{\defaultscriptratio}{#2}}%
#         {\xpic{\defaultscriptscriptratio}{#2}}}}
#     \def\picbig#1#2{{% mathord
#       \mathchoice
#         {\xpicbig{#1}{#2}}%
#         {\xpicbig{#1}{#2}}%
#         {\xpicbig{\defaultscriptratio}{#2}}%
#         {\xpicbig{\defaultscriptscriptratio}{#2}}}}

#     \begin{document}
#     \begin{align*}
#       \begin{bmatrix}
#          LHS_STRING
#       \end{bmatrix} =
#       \begin{bmatrix}
#          ARR_STRING
#       \end{bmatrix}
#       \begin{bmatrix}
#         RHS_STRING   
#       \end{bmatrix}
#     \end{align*}
#     \end{document}
#     """
#     tex_template = tex_template.replace('LHS_STRING', lhs_string)
#     tex_template = tex_template.replace('RHS_STRING', rhs_string)

#     # Convert psi to tex
#     np.savetxt('temp.csv', array, delimiter=' & ', fmt='%2.2f', newline=' \\\\')
#     with open ("temp.csv", "r") as myfile:
#         arr_tex_string = myfile.readlines()
#     tex_template = tex_template.replace('ARR_STRING', ' '.join(arr_tex_string))
    
#     with open(filename+'.tex', 'w') as text_file:
#         print(tex_template, file=text_file)
#     subprocess.call(['latexmk', '-cd', filename+'.tex'])
#     subprocess.call(['rm', 'temp.csv'])

# # Opens .tiff and returns array starting at x, y, z with width, height, and
# # slices dimensions. "None" means return the whole dimension.
# def tiff2array(filename, x=0, y=0, z=0, width=None, height=None, slices=None):
#     im = io.imread(filename)
#     shape = im.shape
#     x_min = x
#     if width is None:
#         x_max = shape[2]
#     else:
#         x_max = x + width
#     y_min = y
#     if height is None:
#         y_max = shape[1]
#     else:
#         y_max = y + height
#     z_min = z
#     if slices is None:
#         z_max = shape[0]
#     else:
#         z_max = z + slices
#     im = im[z_min:z_max, y_min:y_max, x_min:x_max]
#     im = np.swapaxes(im,0,2)
#     return im

# Returns "equally" spaced points on a unit sphere in spherical coordinates.
# http://stackoverflow.com/a/26127012/5854689
def fibonacci_sphere(n, xyz=False):
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
    theta = np.arccos(z)
    phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
    if xyz:
        return np.vstack((np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))).T
    else:
        return np.vstack((theta, phi)).T


# def sphere_LUT_points():
#     n = 256
#     z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
#     theta = np.arccos(z)
#     phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
#     xyz = np.abs(np.vstack((np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta))).T)
#     i = np.lexsort((xyz[:,2], xyz[:,1], xyz[:,0]))
#     return xyz[i]

# def sphere_LUT():
#     points = sphere_LUT_points()
#     lut = np.zeros((256, 4))
#     for i in range(256):
#         lut[i, :] = [points[i, 0], points[i, 1], points[i, 2], 1]
#     return (255*lut).astype('uint8')

# def uvw2color(u, v, w):
#     lut = sphere_LUT()
#     color_idx = []
#     for i in range(u.shape[0]):
#         uvw_comp = 255*np.abs(np.array([u[i], v[i], w[i]]))
#         amin = np.argmin(np.linalg.norm(uvw_comp - lut[:,:3], ord=2, axis=1)**2)
#         color_idx.append(amin)
#     return np.array(color_idx)

# # Convert between spherical harmonic indices (l, m) and matrix index (j)
# def j2lm(j):
#     if j < 0:
#         return None
#     l = 0
#     while True:
#         x = 0.5*l*(l+1)
#         if abs(j - x) <= l:
#             return l, int(j-x)
#         else:
#             l = l+2

# def lm2j(l, m):
#     if abs(m) > l or l%2 == 1:
#         return None
#     else:
#         return int(0.5*l*(l+1) + m)

# def maxl2maxj(l):
#     return int(0.5*(l + 1)*(l + 2))

# def tp2xyz(theta, phi):
#     # Convert spherical to cartesian
#     return np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)

# def plot_sphere(filename=None, directions=None, data=None, show=False,
#                 dpi=500, vis_px=500):

#     # Setup viewing window
#     vispy.use('PyQt4')
#     canvas = vispy.scene.SceneCanvas(keys='interactive', bgcolor='white',
#                                      size=(vis_px, vis_px), show=show, dpi=dpi)
#     my_cam = vispy.scene.cameras.turntable.TurntableCamera(fov=0, elevation=40, azimuth=135,
#                                                            scale_factor=2.2)

#     view = canvas.central_widget.add_view(camera=my_cam)

#     # Plot dots
#     dots = vispy.scene.visuals.Markers(parent=view.scene)
#     dots.antialias = False
#     dots.set_data(pos=np.array([[1.01,0,0],[0,1.01,0],[0,0,1.01]]),
#                   edge_color='black', face_color='black', size=vis_px/50)

#     colors = np.abs(data)
    
#     # Plot sphere
#     sphere = visuals.MySphere(parent=view.scene, radius=1.0,
#                               directions=directions, colors=colors)
    
#     # Display or save
#     im = canvas.render()

#     if show:
#         #visuals.MyXYZAxis(parent=view.scene, origin=[0,1.3,-0.3], length=0.2)
#         vispy.app.run()
    
#     return im
#     #scipy.misc.imsave(filename, im)

# Convert between all spherical harmonic indices (l, m) and matrix index (i)
def i2lm(i):
    if i < 0:
        return None
    l = 0
    while True:
        x = l*(l+1)
        if abs(i - x) <= l:
            return l, int(i-x)
        else:
            l = l+1

def lm2i(l, m):
    if abs(m) > l:
        return None
    else:
        return int(l*(l+1) + m)

def maxl2maxi(l):
    return int(l*(l+2))

# Convert between spherical harmonic indices (l, m) and matrix index (j)
def j2lm(j):
    if j < 0:
        return None
    l = 0
    while True:
        x = 0.5*l*(l+1)
        if abs(j - x) <= l:
            return l, int(j-x)
        else:
            l = l+2

def lm2j(l, m):
    if abs(m) > l or l%2 == 1:
        return None
    else:
        return int(0.5*l*(l+1) + m)

def maxl2maxj(l):
    return int(0.5*(l + 1)*(l + 2))

# Draw microscope schematic (move to Microscope class?)
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

    void circle(triple optical_axis, real Alpha, bool dash, triple color) {
      real h = 1 - sqrt(2 - 2*cos(Alpha) - sin(Alpha)^2);
      real radius = sin(Alpha);
      if (dash) {
        path3 mycircle = circle(c=h*optical_axis, r=radius-0.01, normal=optical_axis);
        draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
      } else {
        path3 mycircle = circle(c=h*optical_axis, r=radius, normal=optical_axis);
        draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
      }
    }

    void mydot(triple pos, triple color) {
      dot(pos, p=rgb(xpart(color), ypart(color), zpart(color)));
    }

    void arrow(triple optical_axis, triple polarizer, triple color, bool dash) {
      if (dash) {
        draw(optical_axis--(optical_axis+0.2*polarizer), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        draw(optical_axis--(optical_axis-0.2*polarizer), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      } else {
        draw(optical_axis--(optical_axis+0.2*polarizer), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
        draw(optical_axis--(optical_axis-0.2*polarizer), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
      }
    }

    real len = 10;
    draw((-len,-len)--(len,-len)--(len,len)--(-len,len)--(-len,-len), white);

    draw(unitsphere, surfacepen=material(diffusepen=white+opacity(0.1), emissivepen=grey, specularpen=white));

    // Draw points on sphere
    dotfactor = 7;
    dot(X); 
    dot(Y); 

    circle((0,0,1), pi/2, false, (0, 0, 0));
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

