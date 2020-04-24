#!/usr/bin/env python
"""
Classes for creating cell stain images.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import shapely
import shapely.geometry

import mersim.base as base
import mersim.util as util


def uniform_spots(polygons, spacing, deltaZ, inFocus, dither):
    """
    Uniformly spaced spots in each Z plane.
    """
    locL = [None, None, None]
    for zv, zPlane in enumerate(polygons):
        for poly in zPlane:
            [tmpX, tmpY] = util.uniform_points_in_shape(poly, spacing, dither)
            tmpZ = zv * deltaZ * np.ones(tmpX.size)
            if not inFocus:
                tmpZ += np.random.uniform(0, deltaZ, tmpX.size)

            locL = util.concat(locL, [tmpX, tmpY, tmpZ])

    return locL


class DAPIImage(base.ImageBase):
    """
    Make DAPI images.
    """
    def foreground(self, config, simParams, fov, iRound, desc):
        image = super().foreground(config, simParams, fov, iRound, desc)
        psf = config["microscope_psf"]

        # Figure out color and z plane.
        color = str(desc[1])
        zPos = float(desc[4])

        # Initialize PSF.
        psf.initialize(config, simParams, color)

        # Load location data.
        [locX, locY, locZ, locInt] = config["dapi_intensity"].load_data(fov)

        # Draw PSFs.
        for i in range(locX.size):
            [x, y, psfImage] = psf.get_psf(locX[i], locY[i], locZ[i], zPos, color)

            # If the PSF is too dim to be relevant psfImage will be None.
            if psfImage is not None:
                psfImage = psfImage * locInt[i]
                util.add_images(image, psfImage, x, y)

        return image


class DAPIIntensityGaussian(base.SimulationBase):
    """
    DAPI dyes with a Gaussian intensity distribution.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)

        # Load barcode positions.
        [locX, locY, locZ] = config["dapi_layout"].load_data()

        # Random intensity information.
        locInt = np.random.normal(self._parameters["intensity_mean"],
                                  self._parameters["intensity_sigma"],
                                  locX.size)
        locInt = np.clip(locInt, 1, None)

        # Add intensity information, save by position.
        #
        fovSize = simParams.get_microscope().get_image_dimensions()

        for fov in range(simParams.get_number_positions()):
            fovRect = simParams.get_fov_rect(fov)
            ox, oy = simParams.get_fov_origin(fov)

            tmpX = []
            tmpY = []
            tmpZ = []
            tmpInt = []
            for j in range(locX.size):
                pnt = shapely.geometry.Point(locX[j], locY[j])
                if fovRect.contains(pnt):
                    tmpX.append(locX[j] - ox)
                    tmpY.append(locY[j] - oy)
                    tmpZ.append(locZ[j])
                    tmpInt.append(locInt[j])

            tmpX = np.array(tmpX)
            tmpY = np.array(tmpY)
            tmpZ = np.array(tmpZ)
            tmpInt = np.array(tmpInt)

            self.save_data([tmpX, tmpY, tmpZ, tmpInt], fov)

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


class DAPIUniform(base.SimulationBase):

    def run_task(self, config, simParams):
        """
        Uniformly spaced dyes in the nucleus.
        """
        super().run_task(config, simParams)

        spacing = self.get_parameter("spacing")
        dither = self.get_parameter("dither", 0.0)
        inFocus = self.get_parameter("in_focus", False)
        
        nZ = simParams.get_number_z()
        deltaZ = simParams.get_z_delta()

        # Load polygons describing sample geometry.
        sampleData = config["sample_layout"].load_data()

        if 'nucleus' in sampleData:
            polygons = sampleData['nucleus']
            assert (nZ == len(polygons))
            locL = uniform_spots(polygons, spacing, deltaZ, inFocus, dither)
        else:
            return

        # Save locations.
        #
        # Note: Array is x, y, z.
        #
        [locX, locY, locZ] = locL
        self.save_data([locX, locY, locZ])
        print("  created {0:d} DAPI locations".format(locX.size))

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

            # Draw locations.
            mask = ((locZ/deltaZ).astype(np.int) == zi)
            plt.scatter(locX[mask], locY[mask], marker = '.')

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


class PolyTImage(base.ImageBase):
    pass

