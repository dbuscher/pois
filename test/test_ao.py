from numpy import *
import sys
from BigScreen import *
from Interferometer import *


Noll=[1.0299,0.582,0.134,0.111,0.0880,0.0648,0.0587,0.0525,0.0463,0.0401,0.0377,0.0352,0.0328,0.0304,0.0279,0.0267,0.0255,0.0243,0.0232,0.0220,0.0208]

def main(diameter=32, r0=32.0, screenSize=1000, numIter=10000, numRemove=5):
    screenGenerator = StaticScreens(screenSize, r0,
                                    diameter, diameter)
    ao=AdaptiveOptics(diameter,diameter, numRemove) 
    aperture=CircularMask(diameter/2.0,diameter,diameter)
    normalisation = Sum2d(aperture)
    variance = 0
    for i in xrange(numIter):
        screen = screenGenerator.next()
        screen=ao.CorrectSingle(screen)
        screen=screen*aperture
        variance = variance + sum(sum(screen**2,-1))
    print normalisation
    variance = variance/numIter/normalisation
    print 'Residual variance:',variance
    print 'Noll 1976 result:',Noll[numRemove]

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        exec(sys.argv[1])
