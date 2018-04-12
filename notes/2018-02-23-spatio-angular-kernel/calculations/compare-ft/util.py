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

def h00(x, NA=0.8, n=1.33):
   return a(x) + 2*b(x, NA, n)

def h20(x, NA=0.8, n=1.33):
   return (-a(x) + 4*b(x, NA, n))/np.sqrt(5)

# OTF functions
def myacos(x):
    return np.nan_to_num(np.arccos(np.abs(x/2)))
    
def mysqrt(x):
    return np.nan_to_num((np.abs(x)/2)*np.sqrt(1 - (np.abs(x)/2)**2))
    
def A(x):
    return (2/np.pi)*(myacos(x) - mysqrt(x))

def B(x, NA=0.8, n=1.33):
    return ((NA**2)/(2*np.pi*(n**2)))*((1 - x**2)*myacos(x) + (1 + ((x**2)/2))*mysqrt(x))

def H00(x, NA=0.8, n=1.33):
    return (A(x) + 2*B(x, NA=NA, n=n))/(1 + 0.5*(NA/n)**2)

def H20(x, NA=0.8, n=1.33):
    return (-A(x) + 4*B(x, NA=NA, n=n))/(np.sqrt(5)*(1 + 0.5*(NA/n)**2))

# File I/O
def save_tiff(image, filename):
    im = Image.fromarray(image) # float32
    im.save(filename, "TIFF")

def load_tiff(filename):
    image = Image.open(filename, mode='r')
    return np.array(image, dtype='float32')

def cs(arr):
    return arr[:, np.int(arr.shape[0]/2)]

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
