import matplotlib.pyplot as plt
import numpy as np
from scipy import special
from scipy.fftpack import fft
import matplotlib
matplotlib.rcParams.update({'font.size': 12})

# PSFs
def a(x):
    return (special.jn(1,2*np.pi*x)/(np.pi*x))**2
def b(x, NA=0.8, n=1.33):
    return (NA/n)**2*(special.jn(2,2*np.pi*x)/(np.pi*x))**2
def psf(x, NA=0.8, n=1.33):
   return a(x) + 2*b(x, NA, n)
def psf2(x, NA=0.8, n=1.33):
   return (-a(x) + 4*b(x, NA, n))/np.sqrt(5)

# Helpers
def acos(x):
    return np.arccos(np.abs(0.5*x))
def mysqrt(x):
    return (np.abs(x)/2)*np.sqrt(1 - (np.abs(x)/2)**2)

def A(x):
    return (2/np.pi)*(acos(x) - mysqrt(x))
def B(x, NA=0.8, n=1.33):
    return ((NA**2)/(2*np.pi*(n**2)))*((1 - x**2)*acos(x) + (1 + ((x**2)/2))*mysqrt(x))

def H00(x, NA=0.8, n=1.33):
    return (A(x) + 2*B(x, NA=NA, n=n))/(1 + 0.5*(NA/n)**2)

def H20(x, NA=0.8, n=1.33):
    return (-A(x) + 4*B(x, NA=NA, n=n))/(np.sqrt(5)*(1 + 0.5*(NA/n)**2))

x = np.linspace(-1, 1, 1000)
NA = 0.8
n = 1.33
y = psf(x, NA=NA)
y2 = psf2(x, NA=NA)
y3 = psf2(x, NA=NA)*(-np.sqrt(5))
f, axs = plt.subplots(1, 2, figsize=(10,5))
#axs[0].plot(x, y0, '-', c='g', lw=1, label=r"xxx$h_0^0(r_d')$")
axs[0].plot(x, y, '-', c='k', lw=1, label=r"$h_0^0(r_d')$")
axs[0].plot(x, y2, '--', c='k',lw=1, label=r"$h_0^2(r_d')$")
axs[0].plot(x, y3, ':', c='k',lw=1, label=r"Normalized $h_0^2(\nu)$")
axs[0].set_xlim([-1, 1])
axs[0].set_ylim([-1, 1])
axs[0].set_xlabel(r"$r_d'\, \textrm{NA}/(M\lambda)$")
axs[0].set_title(r'Spatio-angular PSF')
axs[0].get_xaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
axs[0].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
#axs[0].get_yaxis().set_ticks([0.447])
axs[0].xaxis.labelpad = 10
axs[0].legend(frameon=False, loc=8)

x = np.linspace(-2, 2, 1000)
y = H00(x, NA=NA, n=n)
y2 = H20(x, NA=NA, n=n)#*(np.sqrt(5)*(1 + 0.5*(NA/n)**2))/(-1 + ((NA/n)**2))
y3 = H20(x, NA=NA, n=n)*(np.sqrt(5)*(1 + 0.5*(NA/n)**2))/(-1 + ((NA/n)**2))
#axs[1].plot(x, y0, '-', c='g', lw=1, label=r"xxx$h_0^0(r_d')$")
axs[1].plot(x, y, '-', c='k', lw=1, label=r"$H_0^0(\nu)$")
axs[1].plot(x, y2, '--', c='k',lw=1, label=r"$H_0^2(\nu)$")
axs[1].plot(x, y3, ':', c='k',lw=1, label=r"Normalized $H_0^2(\nu)$")
axs[1].set_xlim([-1, 1])
axs[1].set_ylim([-1, 1])
axs[1].set_xlabel(r"$\nu\, \textrm{NA}/(M\lambda)$")
axs[1].set_title(r'Spatio-angular OTF')
axs[1].get_xaxis().set_ticks([-2, -1, 0, 1, 2])
axs[1].get_yaxis().set_ticks([-1, -0.5, 0, 0.5, 1])
#axs[1].get_yaxis().set_ticks([0.447])
axs[1].xaxis.labelpad = 10
axs[1].legend(frameon=False, loc=8)

plt.savefig('psfs.pdf', dpi=800)
