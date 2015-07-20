import numpy as np
from numpy import fft, random
from numpy import add, sum, zeros, sqrt, arange 

def ScreenGenerator(nfft, r0, nx, ny):
    """Generate an infinite series of rectangular phase screens
    Uses an FFT screen generator to make a large screen and then
    returns non-overlapping subsections of it"""
    while 1:
        layers = GenerateTwoScreens(nfft, r0)
        for iLayer in range(2):
            for iy in range(int(nfft/ny)):
                for ix in range(int(nfft/nx)):
                    yield layers[iLayer][iy*ny:iy*ny+ny, ix*nx:ix*nx+nx]

def GenerateTwoScreens(nfft, r0):
    """

      Generate phase screens with a Kolmogorov spectrum of atmospheric
      disturbances [c.f. Tatarski 1961,1971], such that the phase structure
      function is given by

           D(r) = <[phi(r')-phi(r'+r)]**2>
                = 6.88*(r/r0)**5/3

       where r0 is the Fried parameter.

       This version returns two screens, because it's easier to do it that
       way.
    """
    C = sqrt(0.0229*(float(nfft)/r0)**(5.0/3.0))
    # Generate a 2-d array populated with rsquared=xsquared+ysquared
    r = arange(nfft)
    r[nfft/2:] = nfft-r[nfft/2:]
    rsq = r**2
    rsq = add.outer(rsq,rsq)
    rsq[0, 0] = 1.0 # To solve pole at origin problem
    sample = random.normal(size=(nfft, nfft))+1j*random.normal(size=(nfft, nfft))
    sample *= C*rsq**(-11.0/12.0)
    result = fft.fft2(sample)
    return(result.real, result.imag)

Noll=[1.0299,0.582,0.134,0.111,0.0880,0.0648,0.0587,0.0525,0.0463,0.0401,0.0377,0.0352,0.0328,0.0304,0.0279,0.0267,0.0255,0.0243,0.0232,0.0220,0.0208]


def NollVariance(screen, zernike):
    """ Return the values for the variance of a phase screen as increasing
    numbers of Zernikes are removed"""
    screen = screen*zernike[0] # Cut out a circle
    normalisation = sum(sum(zernike[0],-1))
    result = zeros(len(zernike))
    for i in range(len(zernike)):
        screen = screen - zernike[i]*sum(sum(zernike[i]*screen,-1))/normalisation
        result[i] = sum(sum(screen**2,-1))/normalisation
    return result

def test(diameter=32, screenSize=256, numIter=100):
    '''Test Noll theory on simulated turbulence'''
    import sys
    from Zernike15 import Zernike15

    zernike = Zernike15(diameter, diameter/2.0)
    screenGenerator = ScreenGenerator(screenSize, float(diameter),
                                      diameter, diameter)
    numScreen = 0
    variance = zeros(len(zernike))
    for i in range(numIter):
        screen = next(screenGenerator)
        variance = variance + NollVariance(screen, zernike)
    variance = variance/numIter
    print('Zernikes removed     Residual variance    Noll Result')
    for i in range(len(variance)):
        print('%16d %16f %16f' % (i, variance[i], Noll[i]))


if __name__ == '__main__':
    test(numIter=1000)
