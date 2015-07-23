"""
Test that optimum fiber coupling is for a mode size where the 1/e
radius is 0.90 times the aperture radius

Also, the optimum coupled power is 81%, i.e. coupled amplitude is 90%
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import numpy as np
from poisimulator import *

def Efficiency(modeSize,gridSize=32):
    diameter=gridSize
    pupil=CircularMaskGrid(gridSize,diameter)
    print(np.sqrt(np.sum(pupil)))
    pupil=pupil/(np.sum(pupil))
    return FibreCouple(pupil,diameter*modeSize)

if __name__ == "__main__":
    for modeSize in np.arange(0.5,1.1,0.01):
        print(modeSize,Efficiency(modeSize))
