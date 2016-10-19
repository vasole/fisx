====
fisx
====

Main development website: https://github.com/vasole/fisx

.. image:: https://travis-ci.org/vasole/fisx.svg?branch=master
    :target: https://travis-ci.org/vasole/fisx

.. image:: https://ci.appveyor.com/api/projects/status/github/vasole/fisx?branch=master&svg=true
    :target: https://ci.appveyor.com/project/vasole/fisx

This software library implements formulas to calculate, given an experimental setup, the expected x-ray fluorescence intensities. The library accounts for secondary and tertiary excitation, K, L and M shell emission lines and de-excitation cascade effects. The basic implementation is written in C++ and a Python binding is provided.

Account for secondary excitation is made via the reference:

D.K.G. de Boer, X-Ray Spectrometry 19 (1990) 145-154

with the correction mentioned in:

D.K.G. de Boer et al, X-Ray Spectrometry 22 (1993) 33-28

Tertiary excitation is accounted for via an appproximation.

The accuracy of the corrections has been tested against experimental data and Monte Carlo simulations.

License
-------

This code is relased under the MIT license as detailed in the LICENSE file.

Installation
------------

To install the library for Python just use ``pip install fisx``. If you want build the library for python use from the code source repository, just use one of the ``pip install .`` or the ``python setup.py install`` approaches. It is convenient (but not mandatory) to have cython >= 0.17 installed for it.

Testing
-------

To run the tests **after installation** run::

    python -m fisx.tests.testAll

Example
-------

There is a `web application <http://fisxserver.esrf.fr>`_ using this library for calculating expected x-ray count rates.

This piece of Python code shows how the library can be used via its python binding.

.. code-block:: python

  from fisx import Elements
  from fisx import Material
  from fisx import Detector
  from fisx import XRF

  elementsInstance = Elements()
  elementsInstance.initializeAsPyMca()
  # After the slow initialization (to be made once), the rest is fairly fast.
  xrf = XRF()
  xrf.setBeam(16.0) # set incident beam as a single photon energy of 16 keV
  xrf.setBeamFilters([["Al1", 2.72, 0.11, 1.0]]) # Incident beam filters
  # Steel composition of Schoonjans et al, 2012 used to generate table I
  steel = {"C":  0.0445, 
           "N":  0.04,
           "Si": 0.5093,
           "P":  0.02,
           "S":  0.0175,
           "V":  0.05,
           "Cr":18.37,
           "Mn": 1.619,
           "Fe":64.314, # calculated by subtracting the sum of all other elements
           "Co": 0.109,
           "Ni":12.35,
           "Cu": 0.175,
           "As": 0.010670,
           "Mo": 2.26,
           "W":  0.11,
           "Pb": 0.001}
  SRM_1155 = Material("SRM_1155", 1.0, 1.0)
  SRM_1155.setComposition(steel)
  elementsInstance.addMaterial(SRM_1155)
  xrf.setSample([["SRM_1155", 1.0, 1.0]]) # Sample, density and thickness
  xrf.setGeometry(45., 45.)               # Incident and fluorescent beam angles
  detector = Detector("Si1", 2.33, 0.035) # Detector Material, density, thickness
  detector.setActiveArea(0.50)            # Area and distance in consistent units
  detector.setDistance(2.1)               # expected cm2 and cm.
  xrf.setDetector(detector)
  Air = Material("Air", 0.0012048, 1.0)
  Air.setCompositionFromLists(["C1", "N1", "O1", "Ar1", "Kr1"],
                              [0.0012048, 0.75527, 0.23178, 0.012827, 3.2e-06])
  elementsInstance.addMaterial(Air)
  xrf.setAttenuators([["Air", 0.0012048, 5.0, 1.0],
                      ["Be1", 1.848, 0.002, 1.0]]) # Attenuators
  fluo = xrf.getMultilayerFluorescence(["Cr K", "Fe K", "Ni K"],
                                       elementsInstance,
                                       secondary=2,
                                       useMassFractions=1)
  print("Element   Peak          Energy       Rate      Secondary  Tertiary")
  for key in fluo:
      for layer in fluo[key]:
          peakList = list(fluo[key][layer].keys())
          peakList.sort()
          for peak in peakList:
              # energy of the peak
              energy = fluo[key][layer][peak]["energy"]
              # expected measured rate
              rate = fluo[key][layer][peak]["rate"]
              # primary photons (no attenuation and no detector considered)
              primary = fluo[key][layer][peak]["primary"]
              # secondary photons (no attenuation and no detector considered)
              secondary = fluo[key][layer][peak]["secondary"]
              # tertiary photons (no attenuation and no detector considered)
              tertiary = fluo[key][layer][peak].get("tertiary", 0.0)
              # correction due to secondary excitation
              enhancement2 = (primary + secondary) / primary
              enhancement3 = (primary + secondary + tertiary) / primary
              print("%s   %s    %.4f     %.3g     %.5g    %.5g" % \
                                 (key, peak + (13 - len(peak)) * " ", energy,
                                 rate, enhancement2, enhancement3))

