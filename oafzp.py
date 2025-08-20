
'''  Class for generating an an off-axis Fresnel Zone Plate (FZP) in GDSII format.

Refactored to use the gdstk library for GDSII file generation instead of gdspy.

Grant van Riessen, La Trobe University, 2025

'''

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from tqdm import tqdm
from shapely.geometry import Polygon, LineString, Point, MultiPoint
import gdstk
import math
import uuid
from bisect import bisect_right
from zputils import zones, numberOfZones, convertUnits, longestPolyLine
from p_tqdm import p_map
from itertools import repeat
from functools import partial
from typing import Tuple


ANG_HALF_PI = (-1.0 * np.pi/2, 1.0 * np.pi/2)
GDS_PRECISION = 1.e-9

def zoneIntersectionPoints(i, r_nm, w_nm,
                               roi_bbox, center, ang0, ang1, 
                               tol, layer_kwargs):
        """
        Return (i, [points...]) where each 'points' is an (N,2) list for a clipped polygon.
        """
        x0, y0, x1, y1 = roi_bbox
        rect = gdstk.rectangle((x0, y0), (x1, y1), **layer_kwargs)
        arc  = gdstk.ellipse(center=center,
                            radius=r_nm + w_nm,
                            inner_radius=r_nm,
                            initial_angle=ang0,
                            final_angle=ang1,
                            tolerance=tol,
                            **layer_kwargs)
        clipped = gdstk.boolean(rect, arc, "and", **layer_kwargs) or []
        # Return plain lists (pickle-friendly)
        return i, [pts.points.tolist() for pts in clipped]
    


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
        self.lib = gdstk.Library(unit=self.gdsUnits, precision=self.gdsPrecision)

        #initialize the cell, which is the container for the FZP and photon block
        uniqueStr = f"cellFZP_{uuid.uuid4().hex[:6]}"
        self.cell = self.lib.new_cell(uniqueStr)

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

    def rectCenterMinMaxDistance(self,
                                    cx: float, cy: float,
                                    x0: float, y0: float, x1: float, y1: float) -> Tuple[float, float]:
        """Min/max distance from (cx,cy) to axis-aligned rectangle [x0,x1]x[y0,y1]."""
        # min distance (point to rect)
        dx = max(x0 - cx, 0.0, cx - x1)
        dy = max(y0 - cy, 0.0, cy - y1)
        dmin = math.hypot(dx, dy)
        # max distance (farthest corner)
        dmax = max(
            math.hypot(x0 - cx, y0 - cy),
            math.hypot(x0 - cx, y1 - cy),
            math.hypot(x1 - cx, y0 - cy),
            math.hypot(x1 - cx, y1 - cy),
        )
        return dmin, dmax

    def zoneRadialOverlap(self, r_inner: float, r_outer: float, dmin: float, dmax: float) -> bool:
        """Annulus [r_inner, r_outer] overlaps radial span [dmin, dmax]?"""
        return not (r_outer < dmin or r_inner > dmax)

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
        return gdstk.rectangle(
            (center[0] - self.apertureDim[0] / 2, center[1] - self.apertureDim[1] / 2),
            (center[0] + self.apertureDim[0] / 2, center[1] + self.apertureDim[1] / 2)
        )

    def addAperture(self):
        """Create photon block aperture
        
        The photon block is defined as a rectangle that encompasses the outermost zone, offset from outer edge along x-axis
        by the half-width of the photon block aperture.
        The aperture rectangle is scaled by the scaleAp factor, and the outer rectangle is scaled by the scaleOut factor.
        """
        roi_poly = self.roi()            # gdstk.Polygon
        pts = roi_poly.points            # numpy array, shape (4, 2)
        (x0, y0) = pts[0]                # opposite corners (CCW order)
        (x1, y1) = pts[2]

        apertureRect = gdstk.rectangle((x0, y0), (x1, y1), **self.ldPhotonBlock)

        rect = self.dilate(self.roi(), self.scaleOut)  # Dilate the ROI rectangle by the scale factor
       
        pts = rect.points            # numpy array, shape (4, 2)
        (x0, y0) = pts[0]            # opposite corners
        (x1, y1) = pts[2]

        frameRect = gdstk.rectangle((x0, y0), (x1, y1), **self.ldPhotonBlock)

        if self.mergeInnerZones == 0:
            aperture = gdstk.boolean(frameRect, apertureRect, "not", **self.ldPhotonBlock)
        else:
            # merge partially covered inner zones into photon block. 
            # we that placement tolerance would not allow outer zones to be handled this way.
            i = self.zoneSubset[self.mergeInnerZones-1][0] # Index of the third zone in the subset
            print(f"Adding inner zone for index {i} with radius {self.radii[i]} and width {self.widths[i]}.")
           
            innerZone = gdstk.ellipse(
                                center=self.center_gdsUnits(),
                                radius=self.radii[i] / self.gdsUnits + self.widths[i] / self.gdsUnits,
                                inner_radius=self.radii[i] / self.gdsUnits - 3.0 * self.widths[i] / self.gdsUnits,
                                initial_angle=self.angularRange[0],
                                final_angle=self.angularRange[1],
                                tolerance=self.gdsTolerance,       # curve approximation
                                **self.ldPhotonBlock               # e.g. layer=..., datatype=...
                                )

  
            # First booleans
            inner = gdstk.boolean(frameRect, innerZone, "and", **self.ldPhotonBlock) or []
            a     = gdstk.boolean(frameRect, inner,      "and", **self.ldPhotonBlock) or []
            b     = gdstk.boolean(frameRect, apertureRect, "not", **self.ldPhotonBlock) or []

            # Final OR
            aperture = gdstk.boolean(a, b, "or", **self.ldPhotonBlock) or []

        # Add polygons to the cell: UNPACK the list
        if aperture:
            self.cell.add(*aperture)

    def dilate(self, rect, scaleFactor):
        """
        Dilates a gdspy.Rectangle about its center by a given scale factor.
        """
        # Extract coordinates of the original rectangle
        points = rect.points        # numpy array, shape (4, 2)
        (x0, y0) = points[0]        # lower-left
        (x1, y1) = points[2]        # upper-right

        # Compute center and half-widths
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        hw, hh = (x1 - x0) / 2, (y1 - y0) / 2

        return gdstk.rectangle((cx - hw * scaleFactor, cy - hh * scaleFactor), 
                               (cx + hw * scaleFactor,cy + hh * scaleFactor))
    

    def intersectingZones(self, roi):
        """
        Prefilter zones cheaply, then compute clipped geometry in parallel once.
        Adds intersecting polygons to self.cell and returns a list of booleans A.
        """

        # Prepare only picklable constants to broadcast to workers...
        # ROI bbox as (x0,y0,x1,y1) from polygon points
        pts = roi.points
        x0, y0 = float(np.min(pts[:, 0])), float(np.min(pts[:, 1]))
        x1, y1 = float(np.max(pts[:, 0])), float(np.max(pts[:, 1]))
        roi_bbox = (x0, y0, x1, y1)

        center = tuple(self.center_gdsUnits())
        ang0, ang1 = self.angularRange
        tol = self.gdsTolerance
        ld  = dict(self.ldZones)  # ensure a plain dict

        #  prefilter zones based on geometry (radial only; conservative — no false negatives) ---
        dmin, dmax = self.rectCenterMinMaxDistance(center[0], center[1], x0, y0, x1, y1)

        survivors = []  # (i, r_nm, w_nm)
        for i, (r, w) in enumerate(zip(self.radii, self.widths)):
            r_nm = float(r / self.gdsUnits)
            w_nm = float(w / self.gdsUnits)
            if self.zoneRadialOverlap(r_nm, r_nm + w_nm, dmin, dmax):
                survivors.append((i, r_nm, w_nm))
        
        # could add angular prefilter here...

        # Early-out if nothing can intersect
        A = [False] * len(self.radii)
        if not survivors:
            return A

        # Parallel compute ONCE and return raw points
        results = p_map(
                    zoneIntersectionPoints,
                    [i for (i, _, _) in survivors],
                    [r_nm for (_, r_nm, _) in survivors],
                    [w_nm for (_, _, w_nm) in survivors],
                    repeat(roi_bbox), repeat(center),
                    repeat(ang0), repeat(ang1), repeat(tol), repeat(ld),
                    #num_cpus=num_cpus
                    )       

        # Rewrap points as gdstk.Polygon (cheap) and add to the cell; build boolean list A
        for i, poly_pts_list in results:
            if poly_pts_list:
                A[i] = True
                for pts_list in poly_pts_list:
                    self.cell.add(gdstk.Polygon(np.asarray(pts_list), **ld))

        return A


    def generateZones(self):
        """Generate Fresnel zones and add valid zones to the cell."""

        roi = self.dilate(self.roi(),1.0+self.apertureOverlap)
        points = roi.points        # numpy array, shape (4, 2)
        (x0, y0) = points[0]        # lower-left
        (x1, y1) = points[2]        # upper-right
        #print(f"ROI coordinates: ({x0:.6f}, {y0:.6f}) to ({x1:.6f}, {y1:.6f})")       

        # Get subset of zones that intersect with the aperture rectangle 
        self.zoneSubset = self.intersectingZones(roi)
        
        #print (f"First zone included in aperture: {self.zoneSubset[0] if self.zoneSubset else 'None'}")   
        #print (f"Last zone included in aperture: {self.zoneSubset[-1] if self.zoneSubset else 'None'}")  
        #print(f"{len(self.zoneSubset)} of {i} zones included in aperture rectangle.")   
    
    def arc_tolerance(r, base=2.0, rel=2e-3):
        # units = your layout units (e.g., nm if that’s what you use)
        # use the larger of an absolute and a relative tolerance
        return max(base, rel * r)

    def writeGDS(self, filename='fzp.gds'):
        """Save the GDS library to a file."""
        self.lib.write_gds(filename,max_points=8190)
        print(f"GDS file saved as {filename}")
    
    def writeSVG(self, filename='fzp.svg'):
        """Save the cell as an SVG file."""         
        self.cell.write_svg(filename)
        print(f"SVG file saved as {filename}")

    def writeOAS(self, filename='fzp.oas',compression=1):
        """Save the cell as an OAS file."""

        # Faster/smaller: OASIS
        # 0 = fastest, 9 = smallest (default 6)
        self.lib.write_oas(
            f"{filename}",
            compression_level=compression,        
            detect_rectangles=True,
            detect_trapezoids=True,
            #circletolerance=0,          # allow auto-detection of circles
        )
        print(f"OAS file saved as {filename}")


    def build(self):
        """
        Generate the off-axis FZP and its photon block aperture  
        """      
        self.generateZones()
        self.addAperture()
