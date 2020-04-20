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
    def run_task(self, config, simParams):

        # Load barcode positions.
        [codeX, codeY, codeZ, codeId] = config["layout_barcodes"].load_data()

        # Get barcode information.
        barcodes = simParams.get_codebook().get_barcodes()

        # Random intensity information.
        codeInt = np.zeros((codeX.size, barcodes.shape[1]))
        for i in range(codeX.size):
            bInt = np.random.normal(self._parameters["intensity_mean"],
                                    self._parameters["intensity_sigma"],
                                    barcodes.shape[1])
            dMask = np.random.uniform(size = barcodes.shape[1])
            dMask[(dMask < self._parameters["dropout_rate"])] = 0.0
            dMask[(dMask >= self._parameters["dropout_rate"])] = 1.0

            codeInt[i,:] = bInt * dMask

        # Save by position, by z for each position. We include the
        # barcode ID to make it easier to compare MERlin results to
        # ground truth.
        #
        fovSize = simParams.get_microscope().get_image_dimensions()

        for fov in range(simParams.get_number_positions()):
            fovRect = simParams.get_fov_rect(fov)
            ox, oy = simParams.get_fov_origin(fov)
            print(ox, oy)

            for zi in range(simParams.get_number_z()):
                
                tmpX = []
                tmpY = []
                tmpId = []
                tmpInt = []
                for j in range(codeX.size):
                    if (codeZ[j] == zi):
                        pnt = shapely.geometry.Point(codeX[j], codeY[j])
                        if fovRect.contains(pnt):
                            tmpX.append(codeX[j] - ox)
                            tmpY.append(codeY[j] - oy)
                            tmpId.append(codeId[j])
                            tmpInt.append(codeInt[j,:])

                tmpX = np.array(tmpX)
                tmpY = np.array(tmpY)
                tmpId = np.array(tmpId)
                tmpInt = np.array(tmpInt)

                self.save_data([tmpX, tmpY, tmpId, tmpInt], fov, zi)

                # Make plots.
                fig = plt.figure(figsize = (8,8))

                plt.scatter(tmpX, tmpY, marker = 'x')
                plt.xlim(0, fovSize[0])
                plt.ylim(0, fovSize[1])

                plt.title("fov {0:d}, z {1:d}".format(fov, zi))
                plt.xlabel("pixels")
                plt.ylabel("pixels")

                fname = "fov_{0:d}_{1:d}.pdf".format(fov, zi)
                fig.savefig(os.path.join(self.get_path(), fname),
                            format='pdf',
                            dpi=100)
                
                plt.close()

    
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
