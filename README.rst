Python optical interferometry simulation (pois)
===============================================

|results.png| A python package providing the building blocks to simulate the
operation of a ground-based optical interferometer perturbed by
atmospheric seeing perturbations. It is provided as supplementary
material for the book `“Practical Optical Interferometry”`_, and is
derived from the code which was used to provide data for many of the
figures in the book.


The package includes functions to:

-  generate simulated atmospheric turbulent wavefront perturbations
-  correct these perturbations using adaptive optics
-  combine beams from an arbitrary number of telescopes, with or without
   spatial filtering, to provide complex fringe visibility measurements.

The code has been written following a functional-programming style (in
other words minimising "side-effects" in the code where possible) in order
to try and make it modular and extensible.

Requirements
------------

The simulator runs under Python3 and requires ``numpy``. Some of the
test code uses ``astropy`` for manipulating and saving data
tables.

Installation
------------

On unix-like systems do

::

    pip3 install pois

or if that does not work because of file permission errors, then
::

    sudo pip3 install pois

 
Alternatively download and unpack a copy of this repository and then use

::

    python3 setup.py install


This should install the package into Python path.

Basic usage
-----------

An interferometric simulation can be written as a for loop with some custom data processing, for example:

.. code:: python

    from pois import *

    results=[]
    for phaseScreens in PhaseScreens(numTelescope=3,
                                     r0=15,
				     pupilSize=30,
                                     screenSize=1024,
				     numIter=1000)):
        pupils=AdaptiveOpticsCorrect(phaseScreens,pupilSize=30,radialOrder=5)
        complexPupils=ComplexPupil(pupils)
        fluxes,coherentFluxes=SingleModeCombine(complexPupils)
	# process() is a user-written data processing function
	results.append(process(fluxes,coherentFluxes)) 


More example code using the package is in the `tests`_ directory. The file
`tests/test\_visibility.py`_ includes a complete simulation for
determining visibility losses and single-mode fibre coupling losses as a
function of the diameter of the telescopes.

Pronunciation
-------------

The package name should be pronounced as it would be in the phrase
“petits pois”.

Licencing
---------

The code is licenced under a 2-clause BSD licence (see `LICENCE.txt`_).

.. _“Practical Optical Interferometry”: https://dbuscher.github.io/practical-optical-interferometry/
.. _tests: tests
.. _tests/test\_visibility.py: tests/test_visibility.py
.. _LICENCE.txt: LICENCE.txt

.. |results.png| image:: SNR-vs-diameter.png

