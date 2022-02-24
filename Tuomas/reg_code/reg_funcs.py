import numpy as np
import numpy.fft as fft
from math import pi

def linreg(X,y):
    """ Computes and returns the beta-coefficient
        of a linear regression."""
    return np.sum((X-np.mean(X))*(y-np.mean(y)))/np.sum((X-np.mean(X))**2)

def fftmc(X,Y,t_Beta,Num):
    FX = fft.fft(X)
    theta = 2*pi*np.random.rand(int(len(X)/2),Num)
    fac = np.cos(theta)+1j*np.sin(theta)
    fac = np.concatenate([np.fliplr(np.conjugate(fac)),fac])
    sFX = (fac.T*FX).T
    FFX = np.real(fft.ifft(sFX))
    sigcoeff2 = np.zeros(Num)
    for i in range(Num):
        sigcoeff2[i] = linreg(FFX[:,i],Y)
    sigk = sigcoeff2[np.abs(sigcoeff2) < np.abs(t_Beta)]
    return len(sigk)/Num
