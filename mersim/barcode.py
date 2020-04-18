#!/usr/bin/env python
"""
Classes for barcode layout, intensity calculation and adding to an image.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import shapely
import shapely.geometry

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
    def run_task(self, config, simParams):
        """
        This creates random barcodes in extra-cellular space.
        """
        super().run_task(config, simParams)
        
        nBarcodes = simParams.get_number_barcodes()
        nZ = simParams.get_number_z()

        # Load polygons describing sample geometry.
        sampleData = config["layout_sample"].load_data()
        exPolygons = sampleData["extra-cellular"]

        # Create unions for each z plane.
        zPlaneBounds = []
        zPlanePolys = []
        for elt in exPolygons:
            pUnion = shapely.ops.unary_union(elt)
            zPlaneBounds.append(pUnion.bounds)
            zPlanePolys.append(pUnion)

        # Calculate area.
        totalArea = 0.0
        zPlaneAreas = np.zeros(len(zPlanePolys))
        for i, elt in enumerate(zPlanePolys):
            totalArea += elt.area
            zPlaneAreas[i] = elt.area

        pZ = zPlaneAreas/totalArea

        # Barcodes randomly positioned in the polygons across z-planes.
        #
        # FIXME: Should do by z plane, as I don't think it is correctly
        #        Choosing in it's current form.
        #
        nPts = int(self._parameters["density"] * totalArea)

        codeX = np.zeros(nPts)
        codeY = np.zeros(nPts)
        codeZ = np.zeros(nPts, dtype = np.int)
        codeID = np.zeros(nPts, dtype = np.int)

        cnt = 0
        while(cnt < nPts):
            if ((cnt % 10000) == 0):
                print("  ", cnt, nPts)

            # Choose random Z plane.
            zv = np.random.choice(nZ, p = pZ)

            # Choose random XY.
            minx, miny, maxx, maxy = zPlaneBounds[zv]
            poly = zPlanePolys[zv]

            pnt = shapely.geometry.Point(np.random.uniform(minx, maxx),
                                         np.random.uniform(miny, maxy))

            if poly.contains(pnt):
                codeX[cnt] = pnt.x
                codeY[cnt] = pnt.y
                codeZ[cnt] = zv
                codeID[cnt] = np.random.choice(nBarcodes)
                
                cnt += 1

        # Save barcodes. The 'barcode_intensity' task will use this to
        # create barcode information for each field.
        self.save_data([codeX, codeY, codeZ, codeID])

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
            zPoly = zPlanePolys[zi]
            coords = zPoly.exterior.coords.xy
            x = list(coords[0])
            y = list(coords[1])
            plt.plot(x, y, color = 'black')

            # Draw barcode locations.
            mask = (codeZ == zi)
            plt.scatter(codeX[mask], codeY[mask], marker = 'x')
            
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
