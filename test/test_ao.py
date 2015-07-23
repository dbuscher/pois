from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import numpy as np
from poisimulator import *

# Results from Noll (1976)
Noll=[1.0299,0.582,0.134,0.111,0.0880,0.0648,0.0587,0.0525,0.0463,0.0401,0.0377,0.0352,0.0328,0.0304,0.0279,0.0267,0.0255,0.0243,0.0232,0.0220,0.0208]

Noll=[1.0299,0.582,0.134,0.111,0.0880,0.0648,0.0587,0.0525,0.0463,0.0401,0.0377,0.0352,0.0328,0.0304,0.0279,0.0267,0.0255,0.0243,0.0232,0.0220,0.0208]

def AoCorrect(gridSize=32, r0=32.0, screenSize=1024, numIter=10000, numRemove=6,numTelescope=1):
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
    return(variance,Noll[numRemove-1])

def test_ao():
    variance,noll=AoCorrect(numIter=1000)
    print('Residual variance:',variance)
    print('Noll 1976 result:',noll)
    assert(np.abs(variance-noll)<0.1*noll)
