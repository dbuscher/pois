"""
Test that optimum fiber coupling is for a mode size where the 1/e
radius is 0.90 times the aperture radius

Also, the optimum coupled power is 81%, i.e. coupled amplitude is 90%
"""
from numpy import *
import sys
from Interferometer import Fiber
from NumUtil import *

def main(diameter=32):
    screenSize=diameter
    for modeSize in arange(0.5,1.1,0.01):
        f=Fiber(diameter,diameter,modeSize)
        aperture=CircularMask(diameter/2.0, screenSize, screenSize)
        aperture=aperture/sqrt(Sum2d(aperture))
        print modeSize, f.Couple(f.mode), f.Couple(aperture)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        exec(sys.argv[1])
