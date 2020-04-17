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
import shapely

import mersim.base as base


class SampleUniform(base.SimulationBase):
    """
    Uniform sample, i.e. no cells, etc.
    """
    def layout(self, config, simParams):
        """
        This creates a rectangle that is centered on the FOV. The
        FOV may overlap.

        All dimensions are microns.
        """
        fov_size = simParams.get_microscope().get_fov_size()
        sx = -0.5*fov_size[0] + self._parameters["margin"]
        sy = -0.5*fov_size[1] + self._parameters["margin"]
        ex = sx + image_size[0] - self._parameters["margin"]
        ey = sy + image_size[1] - self._parameters["margin"]

        all_poly = []
        for zi in range(simParams.get_number_z()):

            tmp = []
            for fov in range(simParams.get_number_positions()):
                [px, py] = simParams.get_positions().getFOVXY(fov)

                fov_rect = shapely.geometry.Polygon([px - sx, py - sy],
                                                    [px - sx, py + sy],
                                                    [px + sx, py + sy],
                                                    [px + sx, py - sy])
                tmp.append(fov_rect)

            all_poly.append(tmp)

        self.save_data({"extra-cellular" : all_poly})
