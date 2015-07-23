#!/usr/bin/env python
import sys
from Simulate import ChoosePupilDiameter
from BigScreen import *
from Telemetry import Telemetry
from Interferometer import *
from astropy.io import ascii
import numpy as np

def ChoosePupilDiameter(dr0,minPupilDiameter=32,minR0=3.0):
    return(minPupilDiameter if minPupilDiameter/dr0 > minR0
           else int(np.ceil(dr0*minR0)))

def main(numIter=100000,screenSize=1024,
         numDiameter=50,minDiameter=0.1,maxDiameter=30):
    print "numRemove","d/r0","Vsq","stdVsq","fibrePspec","stdFibrePspec","couple","stdCouple","numIter","pupilSize", "screenSize"
    for numRemove in (2,5,9,14,20):
#    for numRemove in (9,):
        for dr0 in np.logspace(np.log10(minDiameter),np.log10(maxDiameter),
                               numDiameter):
            diameter=ChoosePupilDiameter(dr0,minR0=6.0)
            telemetry=Simulate(diameter=diameter,r0=float(diameter)/dr0,
                               numRemove=numRemove, numIter=numIter,
                               screenSize=screenSize)
            telemetry.Compress()
            v=2.0*telemetry.body["rawFringe"][:,0]
            vsq=v**2
            fiberAmplitudes=telemetry.body["fiberAmplitudes"]
            fiberCohFlux=fiberAmplitudes[:,0]*np.conj(fiberAmplitudes[:,1])
            fiberPowerSpectrum=abs(fiberCohFlux)**2
            couple=2.0*abs(fiberAmplitudes[:,0])**2
            print "%3d %8.4f %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %5d %4d %4d" \
              %(numRemove, dr0,
                np.mean(vsq),
                np.std(vsq,dtype=np.float64),
                np.mean(fiberPowerSpectrum),
                np.std(fiberPowerSpectrum,dtype=np.float64),
                np.mean(couple),
                np.std(couple,dtype=np.float64),
                numIter,diameter,screenSize)

def ScreenSequence(nx,ny,r0,screenSize,
                   timeSampling=None, windDirection=30,):
    if timeSampling:
        return ScrollingScreens(screenSize,r0,nx,ny,
                                dx=timeSampling*r0*sin(windDirection/180.*pi),
                                dy=timeSampling*r0*cos(windDirection/180.*pi))
    else:
        return StaticScreens(screenSize, r0, nx, ny)

def Simulate(diameter, r0, timeSampling=None, windDirection=30,
             screenSize=1024, numIter=10000, numRemove=2,
             makeImages=None):
    screenGenerators=[ScreenSequence(diameter,diameter,r0,
                                     screenSize=screenSize,
                                     timeSampling=timeSampling,
                                     windDirection=windDirection)
                      for i in range(2)]
    telemetry=Telemetry()
    beams=Beams(2,diameter,diameter,numRemove)
    fiberCombiner=FiberCombiner(diameter,diameter)
    pinholeCombiner=PinholeCombiner(diameter,diameter)
    rawCombiner=UnfilteredCombiner(diameter,diameter)
    selfCombiner=SelfCombiner()
    for iter in xrange(numIter):
        beams.Propagate([screenGenerators[0].next(),
                         screenGenerators[1].next(),])
        fiberCombiner.Couple(beams)
        rawCombiner.Combine(beams)
        telemetry.Add({'rawFringe':rawCombiner.fringe,
                       'rawFlux':rawCombiner.flux,
                       'fiberAmplitudes':fiberCombiner.amplitudes,
                       'fiberOpd':fiberCombiner.opd,
                       'aoCoefficients':beams.ao.amplitudes,
        #              'pinholeAmplitudes':pinholeCombiner.amplitudes,
        #               'selfFringe':selfCombiner.fringe,
        #               'selfFlux':selfCombiner.flux,
                       })
    return telemetry

class SelfCombiner:
    def __init__(self):
        pass
    def Combine(self,beams):
        "Use a single phase screen with twice the variance to simulate two"
        pupil=exp(1j*np.sqrt(2)*beams.screens[0])*beams.aperture
        flux = 2*Sum2d(abs(pupil)**2)
        amplitude = Sum2d(pupil*beams.aperture)
        self.fringe=Argand(amplitude)
        self.flux=flux


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        exec(sys.argv[1])
