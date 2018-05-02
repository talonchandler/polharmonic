import numpy as np
from PIL import Image
from scipy import special

# PSF functions
def scalar_a(x):
    if x == 0:
        return 1.0
    else:
        return (special.jn(1,2*np.pi*x)/(np.pi*x))**2
a = np.vectorize(scalar_a)

def s_b(x, NA=0.8, n=1.33):
    if x == 0:
        return 0
    else:
        return (NA/n)**2*(special.jn(2,2*np.pi*x)/(np.pi*x))**2
b = np.vectorize(s_b)

def h00(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    if phi_p==None:
        return a(r_o) + 2*b(r_o, NA, n)
    else:
        return a(r_o) + 4*b(r_o, NA, n)*(np.cos(phi - phi_p)**2)
        
def h20(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    if phi_p==None:
        return (1/np.sqrt(5))*(-a(r_o) + 4*b(r_o, NA, n))
    else:
        return (1/np.sqrt(5))*(-a(r_o) + 8*b(r_o, NA, n)*(np.cos(phi - phi_p)**2))

def h22(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    if phi_p==None:
        return np.zeros(r_o.shape)
    else:
        return np.sqrt(3.0/5.0)*a(r_o)*(np.cos(phi_p)**2 - np.sin(phi_p)**2)

def h2_2(r_o, phi=0, NA=0.8, n=1.33, phi_p=None):
    if phi_p==None:
        return np.zeros(r_o.shape)
    else:
        return -2*np.sqrt(3.0/5.0)*a(r_o)*np.cos(phi_p)*np.sin(phi_p)
    
# OTF functions
def myacos(x):
    return np.nan_to_num(np.arccos(np.abs(x/2)))
    
def mysqrt(x):
    return np.nan_to_num((np.abs(x/2))*np.sqrt(1 - (np.abs(x/2))**2))

def mycubert(x):
    return np.nan_to_num((np.abs(x/2))*((1 - (np.abs(x/2))**2)**(1.5)))

def A(x):
    return (2/np.pi)*(myacos(x) - mysqrt(x))

def B(x, NA=0.8, n=1.33):
    N = (1.0/np.pi)*((NA/n)**2)
    poly = (3.0 - 2.0*(np.abs(x/2)**2))
    return N*(myacos(x) - poly*mysqrt(x))

def C(x, NA=0.8, n=1.33):
    N = (1.0/np.pi)*((NA/n)**2)
    poly  = (4.0/3.0)*(1.0 - 1.0*(np.abs(x/2)**2))
    return -N*poly*mysqrt(x)
    
def H00(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    N = (1 + (NA/n)**2)    
    if phi_p==None:
        return (A(x) + 2*B(x, NA=NA, n=n))/N
    else:
        return (A(x) + 2*B(x, NA=NA, n=n) + 2*C(x, NA=NA, n=n)*(np.cos(2*(phi-phi_p))))/N
    
def H20(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    N = np.sqrt(5)*(1 + (NA/n)**2)
    if phi_p==None:
        return (-A(x) + 4*B(x, NA=NA, n=n))/N
    else:
        return (-A(x) + 4*B(x, NA=NA, n=n) + 4*C(x, NA=NA, n=n)*(np.cos(2*(phi-phi_p))))/N

def H22(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    if phi_p==None:
        return np.zeros(x.shape)
    else:
        return np.sqrt(3.0/5.0)*(A(x)*(np.cos(phi_p)**2 - np.sin(phi_p)**2))/(1 + (NA/n)**2)

def H2_2(x, phi=0, NA=0.8, n=1.33, phi_p=None):
    if phi_p==None:
        return np.zeros(x.shape)
    else:
        return -2*np.sqrt(3.0/5.0)*(A(x)*np.cos(phi_p)*np.sin(phi_p))/(1 + (NA/n)**2)
    
# File I/O
def save_tiff(image, filename):
    im = Image.fromarray(image) # float32
    im.save(filename, "TIFF")

def load_tiff(filename):
    image = Image.open(filename, mode='r')
    return np.array(image, dtype='float32')

def cs(arr, row=False):
    if row:
        return arr[np.int(arr.shape[0]/2), :]
    else:
        return arr[:, np.int(arr.shape[0]/2)]

def ds(arr):
    return np.diagonal(arr)
    
    
# Fourier transform
def myfft(image, pad=1000):
    N = image.shape[0]
    padded_image = np.pad(image, pad_width=pad, mode='constant')
    F = np.fft.fftshift(
        np.fft.fftn(
            np.fft.ifftshift(padded_image)
        ))
    xF = np.fft.fftshift(np.fft.fftfreq(2*pad + N, 4/N))
    return xF, np.abs(F)

# Convert between spherical harmonic indices (l, m) and matrix index (j)
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
