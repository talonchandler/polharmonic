from paratf import *
import numpy as np

from scipy.special import *

def ads_cone(l, theta_p, phi_p, delta):
    if l == 0:
        return Lambda(l)**(-1)
    elif l % 2 == 0:
        Ylm = np.real(sph_harm(0, l, phi_p, theta_p))
        Pl = lpmv(-1, l, np.cos(delta))
        total = np.real(Ylm*Pl/np.tan((delta)/2))
        return total

def phantom_rect(rx, ry, j, k):
    D = 0.15*k+0.15
    return rect((1/D)*np.sqrt((rx - j)**2 + (ry - k)**2))/(D**2)
    
def phantom1g(rx, ry, NA=0.5, n=1.0, nuc=10, dx=0.1):
    out = np.zeros_like(rx)
    for j in range(4):
        for k in range(4):
            out += dpsf(np.sqrt((rx - j)**2 + (ry - k)**2), j*np.pi/6,
                        NA=NA, n=n, nuc=nuc)
    return out

def phantom2g(rx, ry, NA=0.5, n=1.0, nuc=10, dx=0.1):
    freqs = np.fft.fftfreq(rx.shape[-1], d=dx)
    nux, nuy = np.meshgrid(freqs, freqs)
    nu = np.sqrt(nux**2 + nuy**2)

    out = np.zeros_like(rx)
    for j in range(4):
        for k in range(4):
            dd = phantom_rect(rx, ry, j, k)
            sds = np.fft.fft2(dd)
            filt = sdtf(nu, j*np.pi/6, NA=NA, n=n, nuc=nuc)
            out += np.real(np.fft.ifft2(filt*sds))
    return out

def phantom2f(rx, ry):
    out = np.zeros_like(rx)
    for j in range(4):
        for k in range(4):
            out += phantom_rect(rx, ry, j, k)
    return out
    
def phantom3g(rx, ry, NA=0.5, n=1.0, nuc=10, dx=0.1):
    out = np.zeros_like(rx)
    for l in [0, 2]:
        for j in range(4):
            for k in range(4):
                Hlm = adtf(np.sqrt((rx - j)**2 + (ry - k)**2), l,
                           NA=NA, n=n, nuc=nuc)
                Flm = ads_cone(l, theta_p=j*np.pi/6, phi_p=0, delta=k*np.pi/6 + 1e-3)
                out += Hlm*Flm
    return out

def phantom4g(rx, ry, NA=0.5, n=1.0, nuc=10, dx=0.1):
    freqs = np.fft.fftfreq(rx.shape[-1], d=dx)
    nux, nuy = np.meshgrid(freqs, freqs)
    nu = np.sqrt(nux**2 + nuy**2)

    out = np.zeros_like(rx)
    for l in [0, 2]:
        for j in range(4):
            for k in range(4):
                dd = phantom_rect(rx, ry, j, k)
                ads = dd*ads_cone(l, theta_p=0, phi_p=0, delta=j*np.pi/6 + 1e-3)
                sads = np.fft.fft2(ads)
                filt = satf(nu, l, NA=NA, n=n, nuc=nuc)
                out += np.real(np.fft.ifft2(filt*sads))
    return out
