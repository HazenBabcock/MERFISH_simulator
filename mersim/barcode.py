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


def random_barcodes(polygons, nBarcodes, density, deltaZ, inFocus):
    """
    Random barcodes in each z plane.
    """
    codeL = [None, None, None, None]
    for zv, zPlane in enumerate(polygons):
        for poly in zPlane:
            [tmpX, tmpY] = util.random_points_in_shape(poly, density)
            tmpZ = zv * deltaZ * np.ones(tmpX.size)
            if not inFocus:
                tmpZ += np.random.uniform(0, deltaZ, tmpX.size)

            tmpID = np.random.choice(nBarcodes, tmpX.size)

            codeL = util.concat(codeL, [tmpX, tmpY, tmpZ, tmpID])

    return codeL


class BarcodeImage(base.ImageBase):
    """
    Make barcode images.
    """
    def foreground(self, config, simParams, fov, iRound, desc):
        image = super().foreground(config, simParams, fov, iRound, desc)
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
    def foreground(self, config, simParams, fov, iRound, desc):
        image = super().foreground(config, simParams, fov, iRound, desc)
        image += self._parameters["background"]
        return image


class BarcodeIntensityGaussian(base.SimulationBase):
    """
    Barcodes with a Gaussian intensity distribution.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)

        # Load barcode positions.
        [codeX, codeY, codeZ, codeId] = config["barcode_layout"].load_data()

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
    
    def run_task(self, config, simParams):
        """
        This creates random barcodes.
        """
        super().run_task(config, simParams)

        umPerPix = simParams.get_microscope().get_microns_per_pixel()
        density = self.get_parameter("density") * umPerPix * umPerPix

        nBarcodes = simParams.get_number_barcodes()
        print("  {0:d} non blank barcodes found.".format(nBarcodes))
        nZ = simParams.get_number_z()
        deltaZ = simParams.get_z_delta()

        # Check whether barcodes should always be in focus.
        inFocus = self.get_parameter("in_focus", False)

        # Load polygons describing sample geometry.
        sampleData = config["sample_layout"].load_data()

        codeL = [None, None, None, None]
        for pType in ['extra-cellular', 'cytoplasm', 'nucleus']:
            if pType in sampleData:
                polygons = sampleData[pType]
        
                assert (nZ == len(polygons))

                # Random barcodes in each z plane.
                tmpL = random_barcodes(polygons, nBarcodes, density, deltaZ, inFocus)
                codeL = util.concat(codeL, tmpL)

        # Save barcodes. The 'barcode_intensity' task will use this to
        # create barcode information for each field.
        #
        # Note: Array is x, y, z, ID.
        #
        [codeX, codeY, codeZ, codeID] = codeL
        self.save_data([codeX, codeY, codeZ, codeID])
        print("  created {0:d} barcode locations".format(codeX.size))

        # Save a text version for easier MERlin comparison.
        #
        cellIndex = np.zeros(codeX.size, dtype = np.int) - 1
        with open(os.path.join(self.get_path(), "barcodes.csv"), "w") as fp:
            fp.write(",".join(["barcode_id", "global_x", "global_y",
                               "cell_index", "z"]) + "\n")
            np.savetxt(fp,
                       np.column_stack([codeID,
                                        codeX * umPerPix,
                                        codeY * umPerPix,
                                        cellIndex,
                                        codeZ]),
                       fmt = "%d,%.2f,%.2f,%d,%.3f")

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

            # Draw sample geometry.
            for pType in ['extra-cellular', 'cytoplasm', 'nucleus']:
                if pType in sampleData:
                    zPolygons = sampleData[pType][zi]
                    for poly in zPolygons:
                        coords = poly.exterior.coords.xy
                        plt.plot(coords[0], coords[1], color = 'black')

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
