#!/usr/bin/env python
"""
Classes for sample layout.

This creates a dictionary with three entries:
'extra-cellular', 'cytoplasm' and 'nucleus'.

Each entry contains a list of lists. The lists are ordered
by z plane, then by FOV. The list elements are shapely 
Polygons that describe the space occupied by each type of
entry.

For example:

{'extra-cellular' :
  [[p1, p2, p3],
   [p4, p5, p6]]}

Describes the extra-cellular space for 2 z planes using 3 
polygons per z plane.

Polygons in the same z plane of the same type can overlap.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import scipy.spatial
import shapely
import shapely.geometry
import shapely.ops

import mersim.base as base
import mersim.util as util


class Sample2DTissue(base.SimulationBase):
    """
    A '2D' tissue sample, the cells have the same shape in every
    z plane.

    Note: Only creates 'cytoplasm' and 'nucleus'.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)
        
        umPerPix = simParams.get_microscope().get_microns_per_pixel()
        
        nucleusSize = self._parameters["nucleusSize"]/umPerPix
        cellSize = self._parameters["cellSize"]/umPerPix
        shrink = 0.9
        if "shrink" in self._parameters:
            shrink = self._parameters["shrink"]

        # Figure out total imaging area.
        #
        [allFOV, fovUnion] = util.all_fov(simParams)
        [minx, miny, maxx, maxy] = fovUnion.bounds
        
        bndRect = shapely.geometry.Polygon([[minx - cellSize, miny - cellSize],
                                            [minx - cellSize, maxy + cellSize],
                                            [maxx + cellSize, maxy + cellSize],
                                            [maxx + cellSize, miny - cellSize]])

        # Grid random cell center positions.
        #
        xv = np.arange(minx - cellSize, maxx + cellSize, cellSize)
        yv = np.arange(miny - cellSize, maxy + cellSize, cellSize)

        margin = 1.2
        cellLocs = []
        for i in range(xv.size):
            for j in range(yv.size):
                pnt = shapely.geometry.Point(np.random.uniform(xv[i] + margin*nucleusSize,
                                                               xv[i] + cellSize - margin*nucleusSize),
                                             np.random.uniform(yv[j] + margin*nucleusSize,
                                                               yv[j] + cellSize - margin*nucleusSize))

                if bndRect.contains(pnt):
                    cellLocs.append(pnt)

        # Initial cell outlines using Voronoi diagram.
        #
        cx = [c.x for c in cellLocs]
        cy = [c.y for c in cellLocs]
        points = np.column_stack((cx, cy))
        vor = scipy.spatial.Voronoi(points)
    
        # Make shapely Polygons.
        lines = [
            shapely.geometry.LineString(vor.vertices[line])
            for line in vor.ridge_vertices
            if -1 not in line
        ]
        polygons = list(shapely.ops.polygonize(lines))
    
        # Associate center with polygon
        roughCells = []
        for i in range(points.shape[0]):
            pnt = shapely.geometry.Point(points[i,0], points[i,1])
            if fovUnion.contains(pnt):
                for poly in polygons:
                    if poly.contains(pnt):
                        roughCells.append([poly, pnt])
                        break

        # Scale down slightly so that cells don't touch.
        #
        scaledCells = []
        for cell in roughCells:
            [edge, center] = cell
            coords = edge.exterior.coords.xy
            x = list(map(lambda x: shrink * (x - center.x) + center.x, list(coords[0])))
            y = list(map(lambda x: shrink * (x - center.y) + center.y, list(coords[1])))
            scaledCells.append([shapely.geometry.Polygon(zip(x,y)), center])

        # Add nucleus and create new polygon for cytoplasm.
        #
        cytoNucs = []
        for cell in scaledCells:
            [edge, center] = cell

            # Add nucleus.
            angle = np.linspace(0.0, 2.0*np.pi, 20)
            x = center.x + nucleusSize * np.cos(angle)
            y = center.y + nucleusSize * np.sin(angle)
            nucleus = shapely.geometry.Polygon(zip(x,y))
    
            # Create cytoplasm.
            extL = [pt for pt in edge.exterior.coords]
            intL = [pt for pt in nucleus.exterior.coords][::-1]
            cytoplasm = shapely.geometry.Polygon(extL, [intL])

            cytoNucs.append([cytoplasm, nucleus])

        # Organize by type and save.
        #
        allCyto = []
        allNucleus = []
        for zi in range(simParams.get_number_z()):
            tmpCyto = []
            tmpNucleus = []
            for cell in cytoNucs:
                [cyto, nucleus] = cell            
                tmpCyto.append(cyto)
                tmpNucleus.append(nucleus)

            allCyto.append(tmpCyto)
            allNucleus.append(tmpNucleus)

        self.save_data({"cytoplasm" : allCyto,
                        "nucleus" : allNucleus})

        # Reference images.
        #
        for zi in range(simParams.get_number_z()):
            fig = plt.figure(figsize = (8,8))

            # Draw FOV.
            for elt in allFOV:
                coords = elt.exterior.coords.xy
                plt.plot(list(coords[0]), list(coords[1]), color = 'gray')

            # Draw cytoplasm (exterior).
            for elt in allCyto[zi]:
                coords = elt.exterior.coords.xy
                plt.plot(list(coords[0]), list(coords[1]), color = 'black')

            # Draw nucleus (exterior).
            for elt in allNucleus[zi]:
                coords = elt.exterior.coords.xy
                plt.plot(list(coords[0]), list(coords[1]), color = 'blue')

            ax = plt.gca()
            ax.set_aspect('equal', 'datalim')
            
            plt.title("z plane {0:d}".format(zi))
            plt.xlabel("pixels")
            plt.ylabel("pixels")

            fname = "z_{0:d}.pdf".format(zi)
            fig.savefig(os.path.join(self.get_path(), fname),
                        format='pdf',
                        dpi=100)

            plt.close()


class SampleUniform(base.SimulationBase):
    """
    Uniform sample, i.e. no cells, etc.
    """
    def run_task(self, config, simParams):
        """
        This creates a rectangle that is centered on the FOV. The
        FOV may overlap.

        Output dimensions are pixels.
        """
        super().run_task(config, simParams)

        fovSize = simParams.get_microscope().get_image_dimensions()
        umPerPix = simParams.get_microscope().get_microns_per_pixel()
        
        sx = self._parameters["margin"]/umPerPix
        sy = self._parameters["margin"]/umPerPix
        ex = fovSize[0] - self._parameters["margin"]/umPerPix
        ey = fovSize[1] - self._parameters["margin"]/umPerPix

        # Create save polygons.
        allPoly = []
        for zi in range(simParams.get_number_z()):

            tmp = []
            for fov in range(simParams.get_number_positions()):
                [px, py] = simParams.get_fov_xy(fov)

                exRect = shapely.geometry.Polygon([[px + sx, py + sy],
                                                   [px + sx, py + ey],
                                                   [px + ex, py + ey],
                                                   [px + ex, py + sy]])
                tmp.append(exRect)

            allPoly.append(tmp)

        self.save_data({"extra-cellular" : allPoly})

        # Reference images.
        allFOV = []
        for fov in range(simParams.get_number_positions()):
            allFOV.append(simParams.get_fov_rect(fov))

        for zi in range(simParams.get_number_z()):
            fig = plt.figure(figsize = (8,8))

            # Draw FOV.
            for elt in allFOV:
                coords = elt.exterior.coords.xy
                x = list(coords[0])
                y = list(coords[1])
                plt.plot(x, y, color = 'gray')

            # Draw extra-cellular space.
            for elt in allPoly[zi]:
                coords = elt.exterior.coords.xy
                x = list(coords[0])
                y = list(coords[1])
                plt.plot(x, y, color = 'black')

            ax = plt.gca()
            ax.set_aspect('equal', 'datalim')
            
            plt.title("z plane {0:d}".format(zi))
            plt.xlabel("pixels")
            plt.ylabel("pixels")

            fname = "z_{0:d}.pdf".format(zi)
            fig.savefig(os.path.join(self.get_path(), fname),
                        format='pdf',
                        dpi=100)

            plt.close()
