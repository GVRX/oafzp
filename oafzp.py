
'''  Class for generating an an off-axis Fresnel Zone Plate (FZP) in GDSII format.

Grant van Riessen, La Trobe University, 2025

'''

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from shapely.geometry import Polygon, LineString, Point, MultiPoint
import gdspy
import math
from bisect import bisect_right
from zputils import zones, numberOfZones, convertUnits, longestPolyLine


ANG_HALF_PI = (-1.0 * np.pi/2, 1.0 * np.pi/2)
GDS_PRECISION = 1.e-9


class OAFZP:
    """
    Generate an off-axis Fresnel Zone Plate (FZP) in GDSII format.
    
    The FZP is defined by a set of circular arcs, each representing a Fresnel zone. 
    The arcs are defined by their inner and outer radii, which are calculated based on the wavelength and focal length parameters

    The FZP will be surrounded by a photon block, which is a rectangular region that encompasses the outermost zone.
    Overlap of the FZP with the photon block is controlled by the scale_ap parameter, which scales the aperture of the photon block relative to the FZP.
    The outer rectangle in which the FZP and photon block are contained is controlled by the scale_out parameter, which scales the outer 
    rectangle relative to the photon block aperture.

    """

    def __init__(self, 
                 wavelength, # wavelength of the zone plate in units of meters,
                 focalLength, # focal length of the zone plate in units of meters,
                 diameter,  # diameter of the zone plate in units of meters,
                 apertureDim, # Dimensions of the aperture rectangle in which the FZP and photon block are contained, as a tuple (width, height) in meters.
                 aperturePosition=(0, 0), # Position of the aperture rectangle center in meters, relative to the FZP center.
                 apertureOverlap=0.2, # Controls overlap with the FZP. Zones will overlap with the photon block aperture by this fractio nof the aperture dimension in each direction
                 angularRange=ANG_HALF_PI, # Angular range for the arcs: limit speed up the generation of the FZP -- fewer points.
                 frameScale=1.5, # Scale factor for the outer rectangle in which the FZP and photon block are contained. Relative to the PB aperture.
                 mergeInnerZones = 0, # merge inner three zones with the photon block aperture 
                 gdsTolerance=1e-6, # Tolerance for GDSII file 
                 gdsPrecision=1e-9,# Precision for GDSII file (in units of meters)):
                 ):

        self.gdsUnits = GDS_PRECISION # Units for GDSII file (in units of meters)  

        self.gdsPrecision = gdsPrecision
        self.gdsTolerance = gdsTolerance

        self.wavelength = wavelength
        self.focalLength = focalLength
        self.diameter = diameter
        self.apertureDim  = tuple(v / self.gdsUnits for v in apertureDim)
        self.angularRange = angularRange
        self.apertureOverlap = apertureOverlap
        self.aperturePosition = tuple(v / self.gdsUnits for v in aperturePosition)
        self.scaleOut = frameScale
        self.mergeInnerZones = mergeInnerZones

        # Initialize the GDS library 
        self.lib = self.createGdsLibrary()

        #initialize the cell, which is the container for the FZP and photon block
        self.cell = self.lib.new_cell('cellFZP')

        # Define layers for the photon block and Fresnel zones in the GDSII file.
        self.ldPhotonBlock = {"layer": 1, "datatype": 3}
        self.ldZones = {"layer": 2, "datatype": 3}

        # Get radii and widths of the Fresnel zones, and the center coordinates of the FZP.
        self.radii, self.widths = self.zones()

        self.dr = np.min(self.widths)   # outermost zone width in meters
        self.N = np.size(self.widths)  # Number of zones within the aperture 


    def describe(self):
        """Print a description of the Fresnel Zone Plate."""
        print(f"Fresnel zone plate parameters:")
        print(f"Wavelength: {self.wavelength} m")
        print(f"Focal length: {self.focalLength} m")
        print(f"Diameter: {self.diameter} m")
        print(f"Number of zones (full ZP): {self.N}")
        print(f"Outer zone width: {self.dr} nm")
        print(f"Angular range: {self.angularRange}")
        print(f"Photon block dimensions: {self.apertureDim} m")
        print(f"Aperture overlap: {self.apertureOverlap}")
        print(f"Scale frame: {self.scaleOut}")    
        print(f"FZP center: {self.center()} m")
        print(f"\nGDSII Tolerance: {self.gdsTolerance} m")
        print(f"GDSII units: {self.gdsUnits} m")
        print(f"GDSII precision: {self.gdsPrecision} m")


    def createGdsLibrary(self):
        """Create and return a new GDSII library."""
        return gdspy.GdsLibrary(unit=self.gdsUnits, precision=self.gdsPrecision)

    def zones(self):
        """Return the radii and widths of Fresnel zones."""
        return zones(self.wavelength, self.focalLength, numberOfZones(self.wavelength, self.focalLength, self.diameter))
    
    def center(self):
        # Return center coordinates of the Fresnel Zone Plate as a tuple in meters
        return (np.max(self.radii)/2., 
                np.max(self.radii)/2.
                )   
    
    def center_gdsUnits(self):
        """Return the center coordinates of the Fresnel Zone Plate in GDS units."""
        return tuple(v / self.gdsUnits for v in self.center())  

    def roi(self):
        """
        Return the region of interest (ROI) as a rectangle that coincides with the nominal size (before scaling) of the photon block.
        The rectangle is centered at aperturePosition

       """
        center = tuple(map(lambda x, y: x + y, self.center_gdsUnits(), self.aperturePosition))
        return gdspy.Rectangle(
            (center[0] - self.apertureDim[0] / 2, center[1] - self.apertureDim[1] / 2),
            (center[0] + self.apertureDim[0] / 2, center[1] + self.apertureDim[1] / 2)
        )

    def addAperture(self):
        """Create photon block aperture
        
        The photon block is defined as a rectangle that encompasses the outermost zone, offset from outer edge along x-axis
        by the half-width of the photon block aperture.
        The aperture rectangle is scaled by the scaleAp factor, and the outer rectangle is scaled by the scaleOut factor.
        """
        (x0, y0), (x1, y1) = self.roi().polygons[0][0], self.roi().polygons[0][2]
        apertureRect = gdspy.Rectangle((x0, y0), (x1, y1), 
                                       **self.ldPhotonBlock)

        rect = self.dilate(self.roi(), self.scaleOut)  # Dilate the ROI rectangle by the scale factor
        (x0, y0), (x1, y1) = rect.polygons[0][0], rect.polygons[0][2]
        frameRect = gdspy.Rectangle((x0, y0), (x1, y1),
                                         **self.ldPhotonBlock)

        if self.mergeInnerZones == 0:
            aperture = gdspy.boolean(frameRect, apertureRect, "not", **self.ldPhotonBlock)
        else:
            # merge partially covered inner zones into photon block. 
            # we that placement tolerance would not allow outer zones to be handled this way.
            i = self.zoneSubset[self.mergeInnerZones-1][0] # Index of the third zone in the subset
            print(f"Adding inner zone for index {i} with radius {self.radii[i]} and width {self.widths[i]}.")
            innerZone = gdspy.Round(
                    center=self.center_gdsUnits(),
                    radius=self.radii[i]/self.gdsUnits + self.widths[i]/self.gdsUnits,    
                    inner_radius=self.radii[i]/self.gdsUnits - 3. * self.widths[i]/self.gdsUnits,
                    initial_angle=self.angularRange[0],
                    final_angle=self.angularRange[1],
                    tolerance=self.gdsTolerance,
                    **self.ldPhotonBlock
                )
            inner = gdspy.boolean(frameRect, innerZone, "and", **self.ldPhotonBlock)  

            a  = gdspy.boolean(frameRect, inner, "and", **self.ldPhotonBlock)
            b  = gdspy.boolean(frameRect, apertureRect, "not", **self.ldPhotonBlock)
            aperture = gdspy.boolean(a, b, "or", **self.ldPhotonBlock)

        self.cell.add(aperture)

    def dilate(self, rect, scaleFactor):
        """
        Dilates a gdspy.Rectangle about its center by a given scale factor.
        """
        # Extract coordinates of the original rectangle
        (x0, y0), (x1, y1) = rect.polygons[0][0], rect.polygons[0][2]  

        # Compute center and half-widths
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        hw, hh = (x1 - x0) / 2, (y1 - y0) / 2

        return gdspy.Rectangle((cx - hw * scaleFactor, cy - hh * scaleFactor), 
                               (cx + hw * scaleFactor,cy + hh * scaleFactor))
    

    def generateZones(self):
        """Generate Fresnel zones and add valid zones to the cell."""

        roi = self.dilate(self.roi(),1.0+self.apertureOverlap)
        (x0, y0), (x1, y1) = roi.polygons[0][0], roi.polygons[0][2]  
        print(f"ROI coordinates: ({x0:.6f}, {y0:.6f}) to ({x1:.6f}, {y1:.6f})")       

        self.zoneSubset = []
        for i, (r, w) in enumerate(zip(self.radii, self.widths)):
            rNm, wNm = int(r / self.gdsUnits), int(w  / self.gdsUnits)
            arc = gdspy.Round(
                            center=self.center_gdsUnits(),
                            radius=rNm + wNm,
                            inner_radius=rNm,
                            initial_angle=self.angularRange[0],
                            final_angle=self.angularRange[1],
                            tolerance=self.gdsTolerance,
                            **self.ldZones
                            )

            zone = gdspy.boolean(roi, arc, "and", **self.ldZones)

            # Check if the intersects with the aperture rectangle
            if zone is not None:
                    self.cell.add(zone)
                    self.zoneSubset.append((i, rNm, wNm))

        print (f"First zone included in aperture: {self.zoneSubset[0] if self.zoneSubset else 'None'}")   
        print (f"Last zone included in aperture: {self.zoneSubset[-1] if self.zoneSubset else 'None'}")  
        print(f"{len(self.zoneSubset)} of {i} zones included in aperture rectangle.")   

    def writeGDS(self, filename='fzp.gds'):
        """Save the GDS library to a file."""
        self.lib.write_gds(filename)
        print(f"GDS file saved as {filename}")
    
    def writeSVG(self, filename='fzp.svg'):
        """Save the cell as an SVG file."""         
        self.cell.write_svg(filename)
        print(f"SVG file saved as {filename}")

    def display(self):
        """Display the generated FZP in a viewer."""
        gdspy.LayoutViewer()

    def build(self):
        """
        Generate the off-axis FZP and its photon block aperture  
        """      
        self.generateZones()
        self.addAperture()
