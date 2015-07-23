from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import numpy as np
from poisimulator import *

def CouplingEfficiency(modeSize,gridSize=100):
    diameter=gridSize
    pupil=CircularMaskGrid(gridSize,diameter)
    pupil=pupil/np.sqrt(np.sum(pupil))
    return FibreCouple(pupil,diameter*modeSize)[0]**2

def test_fibre():
    """Test that optimum fiber coupling to a uniform circular beam is
    for a mode size where the 1/e radius is 0.90 times the aperture
    radius
    
    Also, the optimum coupled power is about 81%
    """
    assert CouplingEfficiency(0.9)>0.81 and CouplingEfficiency(0.9)<0.82
    assert CouplingEfficiency(0.9)>CouplingEfficiency(0.95)
    assert CouplingEfficiency(0.9)>CouplingEfficiency(0.85)
    
if __name__ == "__main__":
    for modeSize in np.arange(0.5,1.1,0.01):
        print(modeSize,CouplingEfficiency(modeSize))
