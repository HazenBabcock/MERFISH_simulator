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


class CellStainImage(base.ImageBase):
    """
    Make cell stain images.
    """
    def foreground(self, config, simParams, fov, iRound, desc):
        image = super().foreground(config, simParams, fov, iRound, desc)
        psf = config["microscope_psf"]

        # Figure out color and z plane.
        color = str(desc[1])
        zPos = float(desc[4])

        # Initialize PSF.
        psf.initialize(config, simParams, color)

        # Load images.
        fovImages = []
        for zi in range(simParams.get_number_z()):
            fovImages.append(config[self.intensity_name].load_data(fov, zi))

        # Convolve images with PSFs.
        zVals = simParams.get_z_positions()

        for zi in range(simParams.get_number_z()):
            [x, y, psfImage] = psf.get_psf(0.0, 0.0, zVals[zi], zPos, color)
            if psfImage is not None:
                image += util.convolve(fovImages[zi], psfImage)

        return image
 

class CellStainIntensityGaussian(base.SimulationBase):
    """
    Cell stain dyes with a Gaussian intensity distribution.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)

        # Load DAPI images.
        locImages = config[self.layout_name].load_data()

        # Load FOV.
        [allFOV, fovUnion] = util.all_fov(simParams)
        minx, miny, maxx, maxy = list(map(int, fovUnion.bounds))
                    
        # Add intensity information, save by position.
        #
        fovSize = simParams.get_microscope().get_image_dimensions()

        for zi in range(simParams.get_number_z()):
            locImage = locImages[zi].astype(np.float32)
            locImage = locImages[zi] * np.random.normal(self._parameters["intensity_mean"],
                                                        self._parameters["intensity_sigma"],
                                                        locImages[zi].shape)
            locImage = np.clip(locImage, 1, None)
            
            for fov in range(simParams.get_number_positions()):
                ox, oy = list(map(int, simParams.get_fov_origin(fov)))
                ox -= minx
                oy -= miny

                fovImage = locImage[ox:ox + fovSize[0],oy:oy + fovSize[1]]

                # Make plots.
                fig = plt.figure(figsize = (8,8))

                plt.imshow(np.transpose(fovImage), cmap = 'gray')
                plt.xlim(0, fovSize[0])
                plt.ylim(0, fovSize[1])
                
                plt.title("fov {0:d}".format(fov))
                plt.xlabel("pixels")
                plt.ylabel("pixels")
            
                fname = "fov_{0:d}_{1:d}.pdf".format(fov, zi)
                fig.savefig(os.path.join(self.get_path(), fname),
                            format='pdf',
                            dpi=100)
                
                plt.close()
                
                self.save_data(fovImage, fov, zi)


class CellStainUniform(base.SimulationBase):

    def run_task(self, config, simParams):
        """
        Uniformly spaced dyes.
        """
        super().run_task(config, simParams)

        nZ = simParams.get_number_z()

        # Load FOV.
        [allFOV, fovUnion] = util.all_fov(simParams)
        
        # Load polygons describing sample geometry.
        sampleData = config["sample_layout"].load_data()

        locImages = []
        if self.region_stained in sampleData:
            polygons = sampleData[self.region_stained]
            assert (nZ == len(polygons))

            for zi, zPlane in enumerate(polygons):
                locImages.append(util.uniform_fill(zPlane, fovUnion.bounds))
        else:
            return

        # Save locations.
        #
        self.save_data(locImages)

        # Reference images.

        for zi in range(simParams.get_number_z()):
            fig = plt.figure(figsize = (8,8))

            # Draw FOV.
            for elt in allFOV:
                coords = elt.exterior.coords.xy
                x = list(coords[0])
                y = list(coords[1])
                plt.plot(x, y, color = 'gray')

            # Draw DAPI array.
            minx, miny, maxx, maxy = fovUnion.bounds
            plt.imshow(1 - np.transpose(locImages[zi]),
                       extent = [minx, maxx, miny, maxy],
                       origin = ('lower', 'lower'),
                       cmap = 'gray')

            # Draw sample geometry.
            for pType in ['extra-cellular', 'cytoplasm', 'nucleus']:
                if pType in sampleData:
                    zPolygons = sampleData[pType][zi]
                    for poly in zPolygons:
                        coords = poly.exterior.coords.xy
                        plt.plot(coords[0], coords[1], color = 'blue')

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



class DAPIImage(CellStainImage):
    """
    Make DAPI images.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.intensity_name = "dapi_intensity"


class DAPIIntensityGaussian(CellStainIntensityGaussian):
    """
    DAPI dyes with a Gaussian intensity distribution.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.layout_name = "dapi_layout"


class DAPIUniform(CellStainUniform):
    """
    Uniformly space DAPI dyes.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.region_stained = "nucleus"


class PolyTImage(CellStainImage):
    """
    Make polyT images.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.intensity_name = "polyt_intensity"
            

class PolyTIntensityGaussian(CellStainIntensityGaussian):
    """
    PolyT dyes with a Gaussian intensity distribution.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.layout_name = "polyt_layout"


class PolyTUniform(CellStainUniform):
    """
    Uniformly spaced polyT dyes.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.region_stained = "cytoplasm"
