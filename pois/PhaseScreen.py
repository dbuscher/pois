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
