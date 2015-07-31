#!/usr/bin/env python
import sys
from pois import *
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import time

def VisibilityStats(numTelescope,pupilSize,r0,radialOrder,
                    screenSize, numIter):
    """Return visibility statistics from a simulation of an interferometer 
    over numIter phase screens"""
    # Set up per-iteration storage
    numCoherent=numTelescope*(numTelescope-1)//2
    cohFluxes=np.zeros((numIter,numCoherent),dtype=np.complex128)
    fluxes=np.zeros((numIter,numTelescope))
    fibreCohFluxes=np.zeros((numIter,numCoherent),dtype=np.complex128)
    fibreFluxes=np.zeros((numIter,numTelescope))
    # Run simulation
    for iter, phaseScreens in enumerate(PhaseScreens(numTelescope,r0,pupilSize,
                                                     screenSize,numIter)):
        pupils=AdaptiveOpticsCorrect(phaseScreens,pupilSize,radialOrder)
        complexPupils=ComplexPupil(pupils)
        fluxes[iter],cohFluxes[iter]=MultimodeCombine(complexPupils)
        fibreFluxes[iter],fibreCohFluxes[iter]=SingleModeCombine(complexPupils)
    # Postprocess results
    rawPowerSpectrum=abs(cohFluxes)**2
    fibrePowerSpectrum=abs(fibreCohFluxes)**2
    return {'Vsq': np.mean(rawPowerSpectrum)/np.mean(fluxes)**2,
            'stdVsq':np.std(rawPowerSpectrum)/np.mean(fluxes)**2,
            'fibrePspec':np.mean(fibrePowerSpectrum)/np.mean(fibreFluxes)**2,
            'stdFibrePspec':np.std(fibrePowerSpectrum)/np.mean(fibreFluxes)**2,
            'couple':np.mean(fibreFluxes)/np.mean(fluxes),
            'stdCouple':np.std(fibreFluxes)/np.mean(fluxes),
            }


def ChoosePupilDiameter(dr0,minPupilDiameter=32,minR0=3.0):
    """
    Select a minimum pupil diameter so that both r0 and the pupil are
    adequately sampled
    """
    return(minPupilDiameter if minPupilDiameter/dr0 > minR0
           else int(np.ceil(dr0*minR0)))

def main(numIter=1000,screenSize=1024,
         numDiameter=20,minDiameter=0.1,maxDiameter=30,
         radialOrders=(1,)):
    """Run a set of simulations for different interferometer parameters"""
    results=[]
    for radialOrder in radialOrders:
        for dr0 in np.logspace(np.log10(minDiameter),np.log10(maxDiameter),
                               numDiameter):
            diameter=ChoosePupilDiameter(dr0,minR0=6.0)
            r=VisibilityStats(numTelescope=2,
                       pupilSize=diameter,
                       r0=float(diameter)/dr0,
                       radialOrder=radialOrder,
                       numIter=numIter,
                       screenSize=screenSize)
            r.update({'d/r0':dr0,
                      'radialOrder':radialOrder,
                      'numIter':numIter,
                      'pupilSize':diameter,
                      'screenSize':screenSize,})
            results.append(r)
    results=Table(results,names=("radialOrder","d/r0","Vsq","stdVsq","fibrePspec","stdFibrePspec","couple","stdCouple","numIter","pupilSize", "screenSize"))
    ascii.write(results,time.strftime("tmp%y%m%d-%H%M.dat"),
                format='fixed_width',
                bookend=False,
                delimiter=None,
                formats={"d/r0":"%8.4f",
                         "Vsq":"%10.4g",
                         "stdVsq":"%10.4g",
                         "fibrePspec":"%10.4g",
                         "stdFibrePspec":"%10.4g",
                         "couple":"%10.4g",
                         "stdCouple":"%10.4g"})

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        exec(sys.argv[1])
