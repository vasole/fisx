fisx
====

This software library implements formulas to calculate, given an experimental setup, the expected x-ray fluorescence intensities. The library accounts for secondary and tertiary excitation, K, L and M shell emission lines and de-excitation cascade effects. The basic implementation is written in C++ and a Python binding is provided.

License
-------

This code is relased under the MIT license as detailed in the LICENSE file.

Installation
------------

To build the library for python use, just use the "python setup.py install" approach. It is convenient (but not mandatory) to have cython >= 0.17 installed.


Example
-------

This piece of Python code shows how the library can be used via its python binding.

```python
from fisx import Elements
from fisx import Material
from fisx import Detector
from fisx import XRF
elementsInstance = Elements()
elementsInstance.initializeAsPyMca()
xrf = XRF()
xrf.setBeam(16.0) # set incident beam as a single photon energy of 16 keV
xrf.setBeamFilters([["Al1", 2.72, 0.11, 1.0]]) # Incident beam filters
steel = {"Cr": 0.1837,  "Fe": 0.654, "Ni": 0.1235, "Mo": 0.0226, "Mn": 0.01619}
SRM_1155 = Material("SRM_1155", 1.0, 1.0)
SRM_1155.setComposition(steel)
elementsInstance.addMaterial(SRM_1155)
xrf.setSample([["SRM_1155", 1.0, 1.0]]) # Sample, density and thickness
xrf.setGeometry(45., 45.)               # Incident and fluorescent beam angles
detector = Detector("Si1", 2.33, 0.035) # Detector Material, density, thickness
xrf.setDetector(detector)
Air = Material("Air", 0.0012048, 1.0)
Air.setCompositionFromLists(["C1", "N1", "O1", "Ar1", "Kr1"],
                            [0.0012048, 0.75527, 0.23178, 0.012827, 3.2e-06])
elementsInstance.addMaterial(Air)
xrf.setAttenuators([["Air", 0.0012048, 5.0, 1.0],
                    ["Be1", 1.848, 0.002, 1.0]]) # Attenuators
fluo = xrf.getMultilayerFluorescence(["Cr K", "Fe K"],
                                     elementsInstance,
                                     secondary=1)
print("Element   Peak      Energy      Rate     Correction Factor")
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
            # correction due to secondary excitation
            enhancement = (primary + secondary) / primary
            print("%s   %s    %.4f     %.3g     %.3g" % \
                  (key, peak, energy, rate, enhancement))

```
