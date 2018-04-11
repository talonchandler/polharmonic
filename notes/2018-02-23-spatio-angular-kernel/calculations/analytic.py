import matplotlib.pyplot as plt
import numpy as np
# from scipy import special
# from scipy.fftpack import fft

import matplotlib
matplotlib.rcParams.update({'font.size': 12})

# def I0(x):
#     return (2*special.jn(1,x)/x)**2

# def I1(x, NA=0.8, n=1.33):
#     return (NA/n)**2*(2*special.jn(2,x)/x)**2

# def psf(x, NA=0.8, n=1.33):
#    return I0(x) + 2*I1(x, NA, n)

# def psf2(x, NA=1.1, n=1.33):
#    return (-I0(x) + 4*I1(x, NA, n))/np.sqrt(5)

def ac(NA=0.8, M=40, lamb=500):
    return NA/(M*lamb)

def w1(nu, NA=0.8, n=1.33, M=40, lamb=500):
    return ((NA**2)/(2*(n**2)))*(1 - ((nu/ac(NA=NA, M=M, lamb=lamb))**2))

def w2(nu, NA=0.8, n=1.33, M=40, lamb=500):
    return ((NA**2)/(2*(n**2)))*(1 + 0.5*((nu/ac(NA=NA, M=M, lamb=lamb))**2))

def H00(nu, NA=0.8, n=1.33, M=40, lamb=500):
    return norm(nu, NA=NA, n=n, M=M, lamb=lamb)*((w1(nu, NA=NA, n=n, M=M, lamb=lamb) + 1)*np.arccos(nu/(2*ac(NA=NA, M=M, lamb=lamb))) + (w2(nu) - 1)*(nu/(2*ac()))*np.sqrt(1 - (nu/(2*ac(NA=NA, M=M, lamb=lamb)))**2))

def H20(nu, NA=0.8, n=1.33, M=40, lamb=500):
    return norm(nu)*(1.0/np.sqrt(5))*((2*w1(nu, NA=NA, n=n, M=M, lamb=lamb) - 1)*np.arccos(nu/(2*ac())) + (2*w2(nu) + 1)*(nu/(2*ac(NA=NA, M=M, lamb=lamb)))*np.sqrt(1 - (nu/(2*ac(NA=NA, M=M, lamb=lamb)))**2))

def norm(nu, NA=0.8, n=1.33, M=40, lamb=500):
    return 2.0/(np.pi*(w1(0, NA=NA, n=n, M=M, lamb=lamb) + 1))

import pdb; pdb.set_trace() 

# Number of sample points
N = 500000
T = 1.0 / 100 # sample spacing
x = np.linspace(-N*T/2, N*T/2, N)
#y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
y = psf(x)
y2 = psf2(x)
xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
yf = fft(y)
yf2 = fft(y2)

f, axs = plt.subplots(1, 2, figsize=(10,5))
axs[0].plot(x, y, '-', c='k', lw=1, label=r'${I_0^*}^2 + 2{I_1^*}^2$')
axs[0].plot(x, y2, '--', c='k',lw=1, label=r'$\frac{1}{\sqrt{5}}(-{I_0^*}^2 + 4{I_1^*}^2)$')
axs[0].set_xlim([-10, 10])
axs[0].set_ylim([-1.1, 1.1])
axs[0].set_xlabel(r'Position $r_d$')
axs[0].set_ylabel(r'Paraxial point spread function $h^*(r_d)$')
axs[0].set_title(r'$h^*(r_d) \propto ({I_0^*}^2 + 2{I_1^*}^2)Y_0^0 + \frac{1}{\sqrt{5}}(-{I_0^*}^2 + 4{I_1^*}^2)Y_2^0$')
axs[0].get_xaxis().set_ticks([])
axs[0].xaxis.labelpad = 15
axs[0].legend(frameon=False, loc=8)
dc = np.max(2.0/N * np.abs(yf[0:N//2]))
axs[1].plot(xf, 2.0/N * np.abs(yf[0:N//2])/dc, '-', c='k', lw=1, label=r'$\mathcal{F}\{{I_0^*}^2 + 2{I_1^*}^2\}$')
axs[1].plot(xf, 2.0/N * np.abs(yf2[0:N//2]/dc), '--', c='k', lw=1, label=r'$\mathcal{F}\{\frac{1}{\sqrt{5}}(-{I_0^*}^2 + 4{I_1^*}^2)\}$')
axs[1].set_xlim([0, 0.4])
axs[1].set_ylim([0, 1])
axs[1].set_xlabel(r'Spatial frequency $\nu_o$')
axs[1].set_ylabel(r'Paraxial optical transfer function $H_l^l_m(\nu_o)$')
axs[1].set_title(r'${H_l^m}^*(\nu_o) = \mathcal{F}\{h^*(r_d)\}$')
axs[1].get_xaxis().set_ticks([])
axs[1].xaxis.labelpad = 15
axs[1].legend(frameon=False)

plt.savefig('ft.pdf', dpi=800)
