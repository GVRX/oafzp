# Off-Axis Fresnel Zone Plate (OAFZP) Generator

This repository contains Python code for generating **off-axis Fresnel zone plate (OAFZP)** designs as subsections of conventional FZP and exporting them in **GDSII** format. The code supports flexible definition of the FZP geometry and a photon-block layer, and export areas. The code is intended to be used to develop lithographic fabrication workflows.  

## Usage

1. Install dependencies:

   ```bash
   pip install numpy matplotlib shapely gdspy tqdm
   ```

1. Import and run the generator in Python

```python
from zputils import *
from oafzp import OAFZP

# define zone plate parameters
wavelength  = convertUnits(6.70, 'nm')  # Wavelength in meters
focalLength = convertUnits(21., 'mm')  # ZP Focal length in meters
diameter    = convertUnits(960., 'um')    # Diameter of the ZP in meters
aperture    = convertUnits(300., 'um')  # Aperture side length in meters

# Initialise the Fresnel Zone Plate
fzp = OAFZP(wavelength, 
            focalLength, 
            diameter,
            apertureDim = (aperture,aperture),     # Dimensions of the photon block (width, height) in meters
            aperturePosition = (diameter/2. - aperture/2.,0),  # Position of the photon block aperture center in meters, relative to the FZP center.
            apertureOverlap = 0.2,  # Zones will overlap the photon block by this fraction of the aperture dimension in each direction..
            frameScale= 1.50,     # Scale factor for the outer frame in which the FZP and photon block are contained. Relative to the PB aperture.
            mergeInnerZones = 3,  # number of inner zones to merge with photon block 
            gdsTolerance = 1e-5,    
            gdsPrecision = 5e-9,
            )

fzp.build()  # generate to the OAFZP
```

3. Write to GDSII file and/or display in GDS viewer

```python
fzp.writeGDS("fzp.gds")  
fzp.display()  
```

##License

This project is released under the GNU GPL2 license. See [a relative link](LICENSE)
for details.

Grant van Riessen, 2025