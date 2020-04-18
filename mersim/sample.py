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
import os
import shapely
import shapely.geometry
import shapely.ops

import mersim.base as base


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
        
        sx = -0.5*fovSize[0] + self._parameters["margin"]/umPerPix
        sy = -0.5*fovSize[1] + self._parameters["margin"]/umPerPix
        ex = sx + fovSize[0] - self._parameters["margin"]/umPerPix
        ey = sy + fovSize[1] - self._parameters["margin"]/umPerPix

        # Create save polygons.
        allPoly = []
        for zi in range(simParams.get_number_z()):

            tmp = []
            for fov in range(simParams.get_number_positions()):
                [px, py] = simParams.get_fov_xy(fov)

                exRect = shapely.geometry.Polygon([[px - sx, py - sy],
                                                   [px - sx, py + sy],
                                                   [px + sx, py + sy],
                                                   [px + sx, py - sy]])
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
