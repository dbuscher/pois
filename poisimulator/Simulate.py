#!/usr/bin/env python
# Do phase comparison between raw and filtered wavefronts
import sys
from numpy import cos,sin
from BigScreen import *
from Interferometer import *
from Telemetry import Telemetry

def test():
    t=Simulate(simtype='continuous',diameter=32,r0=32,numIter=100)
    print t.header

def ChoosePupilDiameter(dr0,minPupilDiameter=32,minR0=3.0):
    return(minPupilDiameter if minPupilDiameter/dr0 > minR0
           else int(np.ceil(dr0*minR0)))

def Simulate(simtype, diameter, r0, timeSampling=0.1, windDirection=30,
             screenSize=1024, numIter=10000, numRemove=2,
             makeImages=None):
    screenGenerators=[]
    if simtype is 'continuous':
        dx=timeSampling*r0*sin(windDirection/180.*pi)
        dy=timeSampling*r0*cos(windDirection/180.*pi)
        for i in (0,1):
            screenGenerators.append(ScrollingScreens(screenSize, r0,
                                                     diameter, diameter,
                                                     dx=dx,dy=dy))
    elif simtype is 'discrete':
        for i in (0,1):
            screenGenerators.append(StaticScreens(screenSize, r0,
                                                  diameter, diameter))
    else:
        raise 'Programmer error'
    telemetry=Telemetry({'identifier':'Simulate',
                         'simtype':simtype,
                         'revision':(1,3),
                         'diameter':diameter, 'r0':r0, 'screenSize':screenSize,
                         'numRemove':numRemove,'numIter':numIter,
                         'timeSampling':timeSampling,
                         'windDirection':windDirection,
                         'hasImages':makeImages
                         })
    SimulateTwoBeam(screenGenerators,telemetry,diameter,numIter,numRemove,
                    makeImages=makeImages)
    return telemetry

def SimulateTwoBeam(screenGenerators,telemetry,diameter,numIter,numRemove,
                    makeImages=None):
    beams=Beams(2,diameter,diameter,numRemove)
    fiberCombiner=FiberCombiner(diameter,diameter)
    pinholeCombiner=PinholeCombiner(diameter,diameter)
    rawCombiner=UnfilteredCombiner(diameter,diameter)
    for iter in xrange(numIter):
        beams.Propagate([screenGenerators[0].next(),
                         screenGenerators[1].next(),])
        fiberCombiner.Couple(beams)
        pinholeCombiner.Couple(beams)
        rawCombiner.Combine(beams)
        telemetry.Add({'rawFringe':rawCombiner.fringe,
                       'rawFlux':rawCombiner.flux,
                       'fiberAmplitudes':fiberCombiner.amplitudes,
                       'fiberOpd':fiberCombiner.opd,
                       'pinholeAmplitudes':pinholeCombiner.amplitudes,
                       'aoCoefficients':beams.ao.amplitudes,
                       })
        if makeImages:
            beams.MakeImages()
            telemetry.Add({'images':beams.images})

# Boilerplate test execution
if __name__ == '__main__':
    if len(sys.argv) == 1:
        test()
    else:
        exec(sys.argv[1])
