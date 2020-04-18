#!/usr/bin/env python
"""
Classes for barcode layout, intensity calculation and adding to an image.
"""
import numpy as np
import shapely

import mersim.base as base


class BarcodeIntensityGaussian(base.SimulationBase):
    """
    Barcodes with a Gaussian intensity distribution.
    """
    def intensity(self, config, simParams):

        # Load barcode positions.
        [code_x, code_y, code_z, code_id] = config["layout_barcodes"].load_data()

        # Get barcode information.
        barcodes = simParams.get_codebook().get_barcodes()

        # Random intensity information.
        code_int = numpy.zeros((code_x.size, barcodes.shape[1]))
        for i in range(code_x.size):
            b_int = numpy.random.normal(self.parameters["intensity_mean"],
                                        self.parameters["intensity_sigma"],
                                        barcodes.shape[1])
            d_mask = numpy.random.uniform(size = barcodes.shape[1])
            d_mask[(d_mask < self.parameters["dropout_rate"])] = 0.0
            d_mask[(d_mask >= self.parameters["dropout_rate"])] = 1.0

            code_int[i,:] = b_int * d_mask

        # Save by position, by z for each position. We include the
        # barcode ID to make it easier to compare MERlin results to
        # ground truth.
        #
        for fov in range(simParams.get_number_positions()):
            fovRect = simParams.get_fov_rect(fov, simParams)
            ox, oy = simParams.get_fov_origin(fov, simParams)

            for zi in range(simParams.get_number_z()):
                
                tmp_x = []
                tmp_y = []
                tmp_id = []
                tmp_int = []
                for j in range(code_x.size):
                    if (code_z[j] == zi):
                        pnt = shapely.geometry.point(code_x[j], code_y[j])
                        if fov_rect.contains(pnt):
                            tmp_x.append(code_x[j] - ox)
                            tmp_y.append(code_y[j] - oy)
                            tmp_id.append(code_id[j])
                            tmp_int.append(code_int[j,:])

                tmp_x = numpy.array(tmp_x)
                tmp_y = numpy.array(tmp_y)
                tmp_id = numpy.array(tmp_id)
                tmp_int = numpy.array(tmp_int)

            self.save_data([tmp_x, tmp_y, tmp_id, tmp_int], fov)

    
class BarcodeLocationsUniform(base.SimulationBase):
    """
    Uniform array of barcodes (in extra-cellular space).
    """
    def layout(self, config, simParams):
        """
        This creates random barcodes in extra-cellular space.
        """
        n_barcodes = simParams.get_number_barcodes()
        n_z = simParams.get_number_z()

        # Load polygons describing sample geometry.
        sample_data = config["layout_sample"].load_data()
        ex_polygons = sample_data["extra-cellular"]

        # Create unions for each z plane.
        z_plane_bounds = []
        z_plane_polys = []
        for elt in ex_polygons:
            p_union = shapely.ops.unary_union(elt)
            z_plane_bounds.append(p_union.bounds)
            z_plane_polys.append(p_union)

        # Calculate area.
        total_area = 0.0
        z_plane_areas = numpy.zeros(len(z_plane_polys))
        for elt in z_plane_polys:
            total_area += elt.area
            z_plane_areas.append(elt.area)

        p_z = z_plane_areas/total_area

        # Barcodes randomly positioned in the polygons across z-planes.
        n_pts = int(self.parameters["density"] * total_area)

        code_x = numpy.zeros(n_pts)
        code_y = numpy.zeros(n_pts)
        code_z = numpy.zeros(n_pts, dtype = numpy.int)
        code_id = numpy.zeros(n_pts, dtype = numpy.int)

        cnt = 0
        while(cnt < n_pts):

            # Choose random Z plane.
            zv = np.random.choice(n_z, p = p_z)

            # Choose random XY.
            minx, miny, maxx, maxy = z_plane_bounds[zv]
            poly = z_plane_poly[zv]
            
            pnt = shapely.geometry.point(numpy.random.uniform(bounds[minx, maxx]),
                                         numpy.random.uniform(bounds[miny, maxy]))

            if poly.contains(pnt):
                code_x[cnt] = pnt.x
                code_y[cnt] = pnt.y
                code_z[cnt] = zv
                code_id[cnt] = np.random.choice(n_barcodes)
                
                cnt += 1

        # Save barcodes. The 'barcode_intensity' task will use this to
        # create barcode information for each field.
        self.save_data([code_x, code_y, code_z, code_id])
