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

        # Calculate number of bar codes in each z plane.
        umPerPix = simParams.get_microscope().get_microns_per_pixel()
        density = self._parameters["density"] * umPerPix * umPerPix

        totalPnts = 0
        zPlaneCounts = np.zeros(len(zPlanePolys), dtype = np.int)
        for i, elt in enumerate(zPlanePolys):
            totalPnts += int(elt.area * density)
            zPlaneCounts[i] = int(elt.area * density)

        # Barcodes randomly positioned in the polygons in each z-plane.
        nPts = int

        codeX = np.zeros(totalPnts)
        codeY = np.zeros(totalPnts)
        codeZ = np.zeros(totalPnts, dtype = np.int)
        codeID = np.zeros(totalPnts, dtype = np.int)

        npts = 0
        for zv in range(nZ):
            cnt = 0
            print("  creating {0:d} localizations for z plane {1:d}".format(zPlaneCounts[zv], zv))

            minx, miny, maxx, maxy = zPlaneBounds[zv]
            poly = zPlanePolys[zv]
            
            while(cnt < zPlaneCounts[zv]):

                # Choose random XY.
                pnt = shapely.geometry.Point(np.random.uniform(minx, maxx),
                                             np.random.uniform(miny, maxy))

                if poly.contains(pnt):
                    codeX[npts] = pnt.x
                    codeY[npts] = pnt.y
                    codeZ[npts] = zv
                    codeID[npts] = np.random.choice(nBarcodes)
                    npts += 1
                    cnt += 1
        print()

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
