import numpy as np
from Zernike import ZernikeGrid
from PhaseScreen import ScreenGenerator
import sys

def Atmosphere(numTelescope,r0,gridSize,screenSize=1024):
    """
    Return a generator for atmospheric wavefront across a set of 
    numTelescope telescopes. Use next() to get the next set of screens.
    """
    screenGenerators=[ScreenGenerator(screenSize,r0,gridSize,gridSize)
                      for i in range(numTelescope)]
    while 1:
        yield np.array([next(screen) for screen in screenGenerators])

def RadiusGrid(gridSize):
    """
    Return a square grid with values of the radius from the centre 
    of each gridpoint
    """
    x,y=np.mgrid[0:gridSize,0:gridSize]
    x = x-(gridSize-1.0)/2.0
    y = y-(gridSize-1.0)/2.0
    return np.abs(x+1j*y)

def CircularMaskGrid(gridSize, diameter=None):
    """
    Return a square grid with ones inside and zeros outside a given 
    diameter circle
    """
    if diameter is None: diameter=gridSize
    return np.less_equal(RadiusGrid(gridSize),diameter/2.0)


def ComplexPupil(pupils,diameter=None):
    return exp(1j*pupils)*CircularMaskGrid(pupils.shape[-1],diameter)

def AdaptiveOpticsCorrect(pupils,diameter,maxRadial,numRemove=None):
    """
    Correct a wavefront using Zernike rejection up to some maximal order. 
    Can operate on multiple telescopes in parallel.
    Note that this version removes the piston mode as well
    """
    gridSize=pupils.shape[-1]
    pupilsVector=np.reshape(pupils,(-1,gridSize**2))
    zernikes=np.reshape(ZernikeGrid(gridSize,maxRadial,diameter),(-1,gridSize**2))
    if numRemove is None: numRemove=zernikes.shape[0]
    numScreen=pupilsVector.shape[0]
    normalisation=1.0/np.sum(zernikes[0])
    # Note extra iteration to remove residual piston
    for i in list(range(numRemove))+[0,]:
        amplitudes=np.inner(zernikes[i],pupilsVector)*normalisation
        pupilsVector=pupilsVector-zernikes[i]*amplitudes[:,np.newaxis]
    return np.reshape(pupilsVector,pupils.shape)

def FibreMode(gridSize,rmode):
    """
    Return a pupil-plane Gaussian mode with 1/e radius given by 
    rmode, normalised so that integral power over the mode is unity
    """ 
    return np.exp(-(RadiusGrid(gridSize)/rmode)**2/(np.sqrt(pi/2)*rmode))

def FibreCouple(pupils,rmode):
    """
    Return the complex amplitudes coupled into a set of fibers
    """
    gridSize=pupils.shape[-1]
    pupilsVector=np.reshape(pupils,(-1,gridSize**2))
    mode=reshape(FibreMode(gridSize,rmode),(gridSize**2,))
    return np.inner(pupilsVector,mode)

def MonomodeCombine(pupils,rmode):
    """
    Return the instantaneous coherent fluxes and photometric fluxes for a
    multiway monomode fibre combiner
    """
    amplitudes=FibreCouple(pupils,rmode)
    cc=conj(amplitudes)
    fluxes=amplitudes*cc
    coherentFluxes=[amplitudes[i]*cc[j]
                    for i in range(1,len(amplitudes))
                    for j in range(i)]
    return fluxes,coherentFluxes

def MultimodeCombine(pupils):
    """
    Return the instantaneous coherent fluxes and photometric fluxes for a
    multiway multimode combiner (no spatial filtering)
    """
    gridSize=pupils.shape[-1]
    pupilsVector=np.reshape(pupils,(-1,gridSize**2))
    cc=conj(pupils)
    fluxes=np.inner(pupilsVector,cc)
    coherentFluxes=[np.inner(pupilsVector[i],cc[j])
                    for i in range(1,len(amplitudes))
                    for j in range(i)]

Noll=[1.0299,0.582,0.134,0.111,0.0880,0.0648,0.0587,0.0525,0.0463,0.0401,0.0377,0.0352,0.0328,0.0304,0.0279,0.0267,0.0255,0.0243,0.0232,0.0220,0.0208]

def main(gridSize=32, r0=32.0, screenSize=1024, numIter=10000, numRemove=6,numTelescope=1):
    screenGenerator=Atmosphere(numTelescope,r0,gridSize,screenSize)
    aperture=CircularMaskGrid(gridSize)
    normalisation=np.sum(aperture)
    variance=0.0
    for i in range(numIter):
        pupils = next(screenGenerator)
        screen=AdaptiveOpticsCorrect(pupils,gridSize,
                                     maxRadial=5,numRemove=numRemove)
        screen=screen*aperture
        variance = variance + np.sum(screen**2)
    print(normalisation)
    variance = variance/numIter/normalisation
    print('Residual variance:',variance)
    print('Noll 1976 result:',Noll[numRemove-1])

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        exec(sys.argv[1])
