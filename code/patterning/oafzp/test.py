'''
Test script for generating a segment of an off-axis Fresnel Zone Plate (FZP) in GDSII 

This script calculates the radii of the Fresnel zones based on the wavelength and focal length,
and generates a GDSII file containing the zones as circular arcs. The script also includes functionality to visualize the generated FZP.

Grant van Riessen, La Trobe University, 2025

'''

import numpy as np
import math 
from zputils import *
from oafzp import OAFZP

wavelength = convertUnits(185, 'nm')  # Wavelength in meters
focalLength = convertUnits(20, 'mm')  # ZP Focal length in meters
diameter = convertUnits(960, 'um')    # Diameter of the ZP in meters

# Initialise the Fresnel Zone Plate
fzp = OAFZP(wavelength, 
            focalLength, 
            diameter,
            apertureDim = (diameter/6.,diameter/6.),     # Dimensions of the photon block (width, height) in meters
            aperturePosition = (3*diameter/6.,0),  # Position of the photon block aperture center in meters, relative to the FZP center.
            apertureOverlap = 0.2,  # Zones will overlap the photon block by this fraction of the aperture dimension in each direction..
            frameScale= 1.50,     # Scale factor for the outer frame in which the FZP and photon block are contained. Relative to the PB aperture.
            mergeInnerZones = 0,  # merge inner zones with photon block 
            gdsTolerance = 1e-5,    
            gdsPrecision = 5e-9,
            )

fzp.build()  # Build the FZP

fzp.describe()  # Print a description of the FZP

fzp.writeGDS("fzp.gds")  # Write the GDSII file

#fzp.writeSVG("fzp.svg")  # Write the SVG file

fzp.display()  # Display the FZP in a viewer
