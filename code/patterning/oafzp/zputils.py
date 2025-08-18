'''
File contains utility functions for defining Fresnel Zone Plate (FZP) 

Grant van Riessen, La Trobe University, 2025

'''

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from shapely.geometry import Polygon, LineString, Point, MultiPoint
import gdspy
import math

# OpenCascade imports
from OCC.Core.gp import gp_Pnt, gp_Circ, gp_Ax2, gp_Dir
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX
from OCC.Core.BRep import BRep_Tool


def convertUnits(value, unit):
    if unit == 'm':
        return value
    elif unit == 'cm':
        return value / 100
    elif unit == 'mm':
        return value / 1000
    elif unit == 'um':
        return value / 1e6
    elif unit == 'nm':
        return value / 1e9
    else:
        raise ValueError(f"Unknown unit: {unit}")
    
def outermostZoneWidth(wavelength, focal_length, diameter):
    return (wavelength * focal_length) / diameter

def zones(wavelength, focalLength, numZones):    
    # Calculate the  inner radii and width of each odd zone.  
    radii  = [np.sqrt(n * wavelength * focalLength + (n ** 2) * (wavelength ** 2) / 4) for n in range(0, numZones *2)[::2]]
    widths = [np.sqrt((n + 1) * wavelength * focalLength + ((n + 1) ** 2) * (wavelength ** 2) / 4) - radii[i] for i, n in enumerate(range(0, numZones*2)[::2])]
   
    return radii, widths  

def sizeFZP (radii, widths):
    # Calculate the size of the FZP based on the radii, as a tuple in meters
    return ( 2 * (np.max(radii) + 2.*np.min(widths)), 
             2 * (np.max(radii) + 2.*np.min(widths))
    )


def section_points(edge1, edge2):
    """Return list of (x,y) vertices where edge1 and edge2 intersect."""
    pts = []
    sect = BRepAlgoAPI_Section(edge1, edge2)
    sect.Build()
    if not sect.IsDone():
        return pts
    exp = TopExp_Explorer(sect.Shape(), TopAbs_VERTEX)
    while exp.More():
        p = BRep_Tool.Pnt(exp.Current())
        pts.append((p.X(), p.Y()))
        exp.Next()
    return pts


def intersections(gdspy_arc, outer_radius, inner_radius, arc_center, angularRange, gdspy_rect, DEBUG=False):
    # get exact intersection points for both the outer and inner boundaries of the annular arc, 
    # Using OCCT library because this functionality is not available in gdspy.

    theta0, theta1 = angularRange[0], angularRange[1]  


     # --------------------------------------------------------------------------
    # 2.  Convert rectangle into OCCT edges
    # --------------------------------------------------------------------------
    rect_coords = gdspy_rect.polygons[0]
    rect_pts    = [gp_Pnt(float(p[0]), float(p[1]), 0.0) for p in rect_coords]
    rect_edges  = [
        BRepBuilderAPI_MakeEdge(rect_pts[i], rect_pts[(i + 1) % 4]).Edge()
        for i in range(4)
    ]

    # --------------------------------------------------------------------------
    # 3.  Build trimmed OCCT edges for the outer *and* inner arcs
    # --------------------------------------------------------------------------
    ax2    = gp_Ax2(gp_Pnt(float(arc_center[0]), float(arc_center[1]), 0.0),
                    gp_Dir(0, 0, 1))

    outer_edge = BRepBuilderAPI_MakeEdge(
        gp_Circ(ax2, float(outer_radius)), theta0, theta1
    ).Edge()

    inner_edge = None
    if inner_radius > 0:
        inner_edge = BRepBuilderAPI_MakeEdge(
            gp_Circ(ax2, float(inner_radius)), theta0, theta1
        ).Edge()

    # --------------------------------------------------------------------------
    # 4.  Intersect each rectangle edge with both arc edges
    # --------------------------------------------------------------------------
    crossings_outer = []
    crossings_inner = []

    for re in rect_edges:
        crossings_outer += section_points(re, outer_edge)
        if inner_edge:
            crossings_inner += section_points(re, inner_edge)

    # Deduplicate (optional) — OpenCascade may return the same vertex twice
    crossings_outer = list(dict.fromkeys(crossings_outer))
    crossings_inner = list(dict.fromkeys(crossings_inner))

    # --------------------------------------------------------------------------
    # 5.  Report
    # --------------------------------------------------------------------------
    if DEBUG: 
        print(f"Outer boundary crossings ({len(crossings_outer)}):")
        for i, (x, y) in enumerate(crossings_outer, 1):
            print(f"  {i}: ({x:.6f}, {y:.6f})")

        if inner_edge:
            print(f"\nInner boundary crossings ({len(crossings_inner)}):")
            for i, (x, y) in enumerate(crossings_inner, 1):
                print(f"  {i}: ({x:.6f}, {y:.6f})")

    return crossings_inner + crossings_outer

def longestPolyLine(poly):
    """
    Find the maximum length of any edge (line segment) in a gdspy Polygon or PolygonSet.

    Parameters
    ----------
    zone : gdspy.Polygon or gdspy.PolygonSet
        The polygon(s) defining the clipped arc region.

    Returns
    -------
    float
        Length of the longest edge in the polygon(s).
    """
    if poly is not None and len(poly.polygons) > 0:
        return 0.0

    polygons = zone.polygons if isinstance(zone, gdspy.PolygonSet) else [zone.points]
    max_length = 0.0

    for poly in polygons:
        n = len(poly)
        for i in range(n):
            p1 = poly[i]
            p2 = poly[(i + 1) % n]
            edge_length = np.linalg.norm(p2 - p1)
            if edge_length > max_length:
                max_length = edge_length

    return max_length

def numberOfZones(wavelength, focal_length, diameter):
    """
    Calculate the number of Fresnel zones in a zone plate.

    Parameters
    ----------
    wavelength : float
        Wavelength of light in meters.
    focal_length : float
        Focal length of the zone plate in meters.
    diameter : float
        Diameter of the zone plate in meters.

    Returns
    -------
    int
        Number of Fresnel zones (rounded down to nearest whole number).
    """
    N = (diameter ** 2) / (4 * wavelength * focal_length)
    return int(N)


def zonesInAnnulus(wavelength, focal_length, r_center, annulus_width):
    '''
    Calculate the number of Fresnel zones lying within an annular region
    of specified width and center radius.

    Parameters
    ----------
    wavelength : float
        Wavelength in meters.
    focal_length : float
        Focal length of the Fresnel zone plate in meters.
    r_center : float
        Radial distance from the zone plate center to the center of the aperture (in meters).
    annulus_width : float
        Width of the annulus/aperture (in meters).

    Returns
    -------
    int
        Number of Fresnel zones that intersect the annular region.
    '''
    r_outer = r_center + annulus_width / 2
    r_inner = max(0.0, r_center - annulus_width / 2)

    n_max = math.floor(r_outer**2 / (wavelength * focal_length))
    n_min = math.ceil(r_inner**2 / (wavelength * focal_length))

    return max(0, n_max - n_min + 1)