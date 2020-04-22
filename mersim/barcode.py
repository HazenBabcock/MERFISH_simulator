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
import mersim.util as util


class BarcodeImage(base.ImageBase):
    """
    Make barcode images.
    """
    def make_image(self, config, simParams, fov, iRound, desc):
        image = super().make_image(config, simParams, fov, iRound, desc)
        psf = config["microscope_psf"]

        # Figure out which bit this is an image of.
        color = str(desc[1])
        bitNum = int(desc[2]) - 1
        zPos = float(desc[4])

        # Initialize PSF.
        psf.initialize(config, simParams, color)

        # Load barcode data.
        [codeX, codeY, codeZ, codeId, codeInt] = config["barcode_intensity"].load_data(fov)

        # Draw PSFs for barcodes that are 'on' for this bit.
        for i in range(codeX.size):
            if (codeInt[i,bitNum] > 0.0):
                [x, y, psfImage] = psf.get_psf(codeX[i], codeY[i], codeZ[i], zPos, color)

                # If the PSF is too dim to be relevant psfImage will be None.
                if psfImage is not None:
                    psfImage = psfImage * codeInt[i, bitNum]
                    util.add_images(image, psfImage, x, y)

        return image


class BarcodeImageUniformBackground(BarcodeImage):
    """
    Make barcode images with a uniform background.
    """
    def make_image(self, config, simParams, fov, iRound, desc):
        image = super().make_image(config, simParams, fov, iRound, desc)
        image += self._parameters["background"]
        return image


class BarcodeIntensityGaussian(base.SimulationBase):
    """
    Barcodes with a Gaussian intensity distribution.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)

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

            codeInt[i,:] = bInt * dMask * barcodes[codeId[i],:]

        # Add intensity information, save by position. We include the
        # barcode ID to make it easier to compare MERlin results to
        # ground truth.
        #
        fovSize = simParams.get_microscope().get_image_dimensions()

        for fov in range(simParams.get_number_positions()):
            fovRect = simParams.get_fov_rect(fov)
            ox, oy = simParams.get_fov_origin(fov)

            tmpX = []
            tmpY = []
            tmpZ = []
            tmpId = []
            tmpInt = []
            for j in range(codeX.size):
                pnt = shapely.geometry.Point(codeX[j], codeY[j])
                if fovRect.contains(pnt):
                    tmpX.append(codeX[j] - ox)
                    tmpY.append(codeY[j] - oy)
                    tmpZ.append(codeZ[j])
                    tmpId.append(codeId[j])
                    tmpInt.append(codeInt[j,:])

            tmpX = np.array(tmpX)
            tmpY = np.array(tmpY)
            tmpZ = np.array(tmpZ)
            tmpId = np.array(tmpId)
            tmpInt = np.array(tmpInt)

            self.save_data([tmpX, tmpY, tmpZ, tmpId, tmpInt], fov)

            # Make plots.
            fig = plt.figure(figsize = (8,8))

            plt.scatter(tmpX, tmpY, marker = 'x')
            plt.xlim(0, fovSize[0])
            plt.ylim(0, fovSize[1])

            plt.title("fov {0:d}".format(fov))
            plt.xlabel("pixels")
            plt.ylabel("pixels")
            
            fname = "fov_{0:d}.pdf".format(fov)
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
        print("  {0:d} non blank barcodes found.".format(nBarcodes))
        nZ = simParams.get_number_z()
        zPos = simParams.get_z_positions()

        # We assume that the z positions are equally spaced.
        deltaZ = 0.0
        if (len(zPos) > 1):
            deltaZ = zPos[1] - zPos[0]

        # Check whether barcodes should always be in focus.
        inFocus = False
        if "in_focus" in self._parameters:
            inFocus = self._parameters["in_focus"]

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
        codeZ = np.zeros(totalPnts)
        codeID = np.zeros(totalPnts, dtype = np.int)

        npts = 0
        for zv in range(nZ):
            cnt = 0
            print("  creating {0:d} barcodes for z plane {1:d}".format(zPlaneCounts[zv], zv))

            minx, miny, maxx, maxy = zPlaneBounds[zv]
            poly = zPlanePolys[zv]
            
            while(cnt < zPlaneCounts[zv]):

                # Choose random XY.
                pnt = shapely.geometry.Point(np.random.uniform(minx, maxx),
                                             np.random.uniform(miny, maxy))

                if poly.contains(pnt):
                    codeX[npts] = pnt.x
                    codeY[npts] = pnt.y
                    if inFocus:
                        codeZ[npts] = zv*deltaZ
                    else:
                        codeZ[npts] = zv*deltaZ + np.random.uniform(0, deltaZ)
                    codeID[npts] = np.random.choice(nBarcodes)
                    npts += 1
                    cnt += 1
        print()

        # Save barcodes. The 'barcode_intensity' task will use this to
        # create barcode information for each field.
        self.save_data([codeX, codeY, codeZ, codeID])

        # Save a text version for easier MERlin comparison.
        #
        cellIndex = np.zeros(codeX.size, dtype = np.int) - 1
        with open(os.path.join(self.get_path(), "barcodes.csv"), "w") as fp:
            fp.write(",".join(["barcode_id", "global_x", "global_y", "cell_index"]) + "\n")
            np.savetxt(fp,
                       np.column_stack([codeID,
                                        codeX * umPerPix,
                                        codeY * umPerPix,
                                        cellIndex]),
                       fmt = "%d,%.2f,%.2f,%d")

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
            if isinstance(zPoly, shapely.geometry.MultiPolygon):
                tmp = zPoly
            else:
                tmp = [zPoly]

            for poly in tmp:
                coords = poly.exterior.coords.xy
                x = list(coords[0])
                y = list(coords[1])
                plt.plot(x, y, color = 'black')

            # Draw barcode locations.
            mask = ((codeZ/deltaZ).astype(np.int) == zi)
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
