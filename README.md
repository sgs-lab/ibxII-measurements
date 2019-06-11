# IBX II Measurements

Information Barrier Experimental II is an information barrier setup for nuclear warhead verification with an Apple IIe ([Details](http://www.vintageverification.org)). An early version was presented at [34C3](https://media.ccc.de/v/34c3-8994-vintage_computing_for_trusted_radiation_measurements_and_a_world_free_of_nuclear_weapons).

This repository holds measurement data and the python processing scripts that can transform the measurements into figures. A [second repository](https://github.com/sgs-lab/ibxII) holds assembler code and hardware design.

## Radioactive Sources 

Four test sources from [Spectrum Techniques](http://www.spectrumtechniques.com) were used.

Isotope | ID | Initial Activity | Production Month | Activity (June 2019)
--------|----|------------------|------------------|---------------------
Na-22 | 14-6760 | 37kBq | 11/2014 | 11kBq
Cs-137 | 14-6759 | 37kBq | 11/2014 | 33kBq
Co-60 | 14-6758 | 37kBq | 10/2014 | 20kBq       
Co-60 | 14-6766 | 37kBq | 10/2014 | 20kBq

In addition, thoriated welding rods are used as a calibration source. Twelve 1-inch pieces of 4mm thick welding rods surround the scintillator crystal in a 3D-printed holder, around the front end of the detector (away from the photomultimplier).

## Measurements to Create a Lookup Table

Energy calibration of the gamma spectra recorded with IBX II is divided in two parts. First, the software on the Apple II identifies a peak at 2.614 MeV from the decay of Ti-208. Ti-208 is part of the decay chain of thorium. This allows the calibration on a linear scale. To reduce the computational load on the MOS6502 processor (1MHz, 8bit, no floating point capabilities) we use a simple lookup table for the correction of a non-linearity. The lookup table is precalculated using Python scripts on a modern computer.

The current lookup table uses measurements of the Na-22 and the Cs-137 source. Both sources were located 35 mm away from the surface of the NaI crystal, measured as the distance between surfaces. Vertical placement was at the center line of the crystal. The Na-22 was placed closer to the NaI crystal. The voltage on the PMT was set to 1000V.

150 individual measurements were carried out on June 05, 2019. Each measurement had 2^18 (262144 or 0x40000) counts, and the measurements were stored on two Apple II disk images (`lookup-35mm-na-cs-01/02.nib` in the repository). The spectra can be found in the `lookup-35mm-na-cs-01/02-data` folders.

All spectra were individually calibrated linearly, so that the Ti-208 peak would lie in channel 226. This was done using a python routine that replicates the Apple II assembler code exactly. Using the average of all measurements, the positions of the five peaks listed in the table below were identified and used to fit a second order polynomial.

Energy (keV) | Source
-------------|-------
238.632 | Pb-212 (thorium decay chain)
511.0   | Na-22
661.64  | Cs-137
1274.5 | Na-22
2614.511 | Ti-208 (thorium decay chain)

The resulting lookup table was used for the following measurements.

Using the same sources, and the Python calibration/lookup routines, 150 measurement (same number of counts, 2^18) were carried out using a Mirion Osprey Base on June 07, 2019. The results are stored in `osprey-35mm-na-cs`.

## Measurements for "reference"/"valid" and "invalid" items

The key idea feature of IBX II is to compare gamma spectra and identify similar distributions. We call one spectrum "valid item", and a second, slightly different spectrum "invalid item". The measurements for valid items are also used to create a "template", a reference item.

The valid item was represented by two Co-60 sources (IDs: 14-6758, 14-6766). They were placed 35 mm away from the surface of the NaI crystal, measured as the distance between surfaces. Vertical placement was at the center line of the crystal. The voltage on the PMT was 1000V. We took 150 individual measurements on June 06, 2019. Each measurement had 2^18 (262144 or 0x40000) counts, and the measurements were stored on two Apple II disk images (`valid-35mm-na-cs-01/02.nib` in the repository). The spectra can be found in the `valid-35mm-na-cs-01/02-data` folders. The measurements include results for which calibration/lookup has been carried out on the Apple II. These are used to create the figures (`M<NNN>.dat` files). The original, uncalibrated data is also stored (`M<NNN>_original.dat` files).

For the invalid item, the Cs-137 source (ID: 14-6759) was added to the two Co-60 sources from the valid item. The Cs-137 was placed in a distance of 140mm away from the surface of the NaI crystal, measured as the distance between surfaces. Vertical placement was at the center line of the crystal. The voltage on the PMT was 1000V. We took 75 individual measurements on June 06, 2019 and 75 on June 07, 2019. Each measurement had 2^18 (262144 or 0x40000) counts, and the measurements were stored on two Apple II disk images (`invalid-35mm-na-cs-01/02.nib` in the repository). The spectra can be found in the `invalid-35mm-na-cs-01/02-data` folders. The measurements include results for which calibration/lookup has been carried out on the Apple II. These are used to create the figures (`M<NNN>.dat` files). The original, uncalibrated data is also stored (`M<NNN>_original.dat` files).

## Python Scripts

* intcal.py is a small module that implements calibration and lookup using only integer multiplications (blueprint for assembler, replicates assembler functions exactly)
* createlookup.py can be used to recreate the lookup table for non-linear corrections in spectrum. It currently uses the two sets of measurements described above
* create-plots.py creates a number of plots
* fullconvert.py can be used to convert the binary files from the Apple IIe into text files with the spectrum
* debughelper.py outputs a single spectrum data set usable for assembler
