#!/usr/bin/env python
"""
Classes for fiducial layout, intensity calculation and adding to an image.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import shapely
import shapely.geometry
import shapely.ops

import mersim.base as base
import mersim.util as util


class FiducialImage(base.ImageBase):
    """
    Make fiducial images.
    """
    def foreground(self, config, simParams, fov, iRound, desc):
        image = super().foreground(config, simParams, fov, iRound, desc)
        psf = config["microscope_psf"]
        color = str(desc[1])

        # Initialize PSF.
        psf.initialize(config, simParams, color)

        # Load barcode data.
        [fidX, fidY, fidInt] = config["fiducial_intensity"].load_data(fov)

        # Draw PSFs for fiducials.
        for i in range(fidX.size):
            [x, y, psfImage] = psf.get_psf(fidX[i], fidY[i], 0.0, 0.0, color)

            # If the PSF is to dim too be relevant psfImage will be None.
            if psfImage is not None:
                psfImage = psfImage * fidInt[i]
                util.add_images(image, psfImage, x, y)

        return image


class FiducialImageUniformBackground(FiducialImage):
    """
    Make fiducial images with a uniform background.
    """
    def foreground(self, config, simParams, fov, iRound, desc):
        image = super().foreground(config, simParams, fov, iRound, desc)    
        image += self._parameters["background"]
        return image

    
class FiducialIntensityGaussian(base.SimulationBase):
    """
    Fiducials with Gaussian intensity distribution.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)

        # Load fiducial positions.
        [fidX, fidY] = config["fiducial_layout"].load_data()

        # Random normal intensities.
        fidInt = np.random.normal(self._parameters["intensity_mean"],
                                  self._parameters["intensity_sigma"],
                                  fidX.size)
        
        # Add intensity information, save by position.
        #
        fovSize = simParams.get_microscope().get_image_dimensions()

        for fov in range(simParams.get_number_positions()):
            fovRect = simParams.get_fov_rect(fov)
            ox, oy = simParams.get_fov_origin(fov)

            tmpX = []
            tmpY = []
            tmpInt = []
            for j in range(fidX.size):
                pnt = shapely.geometry.Point(fidX[j], fidY[j])
                if fovRect.contains(pnt):
                    tmpX.append(fidX[j] - ox)
                    tmpY.append(fidY[j] - oy)
                    tmpInt.append(fidInt[j])

            tmpX = np.array(tmpX)
            tmpY = np.array(tmpY)
            tmpInt = np.array(tmpInt)

            self.save_data([tmpX, tmpY, tmpInt], fov)

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
            

class FiducialLocationsUniform(base.SimulationBase):
    """
    Barcodes uniformly distributed in each FOV.
    """
    def run_task(self, config, simParams):
        super().run_task(config, simParams)

        fovSize = simParams.get_microscope().get_image_dimensions()
        umPerPix = simParams.get_microscope().get_microns_per_pixel()
        
        sx = self._parameters["margin"]/umPerPix
        sy = self._parameters["margin"]/umPerPix
        ex = fovSize[0] - self._parameters["margin"]/umPerPix
        ey = fovSize[1] - self._parameters["margin"]/umPerPix

        # Create a polygon for each FOV, adjusted for margin.
        fovPoly = []
        for fov in range(simParams.get_number_positions()):
            [px, py] = simParams.get_fov_xy(fov)

            fovRect = shapely.geometry.Polygon([[px + sx, py + sy],
                                                [px + sx, py + ey],
                                                [px + ex, py + ey],
                                                [px + ex, py + sy]])
            fovPoly.append(fovRect)

        # Randomly place fiducials across FOV.
        fovUnion = shapely.ops.unary_union(fovPoly)
        density = self._parameters["density"] * umPerPix * umPerPix

        [fidX, fidY] = util.random_points_in_shape(fovUnion, density)
        self.save_data([fidX, fidY])

        # Reference image.
        allFOV = []
        for fov in range(simParams.get_number_positions()):
            allFOV.append(simParams.get_fov_rect(fov))
        
        fig = plt.figure(figsize = (8,8))

        # Draw FOV.
        for elt in allFOV:
            coords = elt.exterior.coords.xy
            x = list(coords[0])
            y = list(coords[1])
            plt.plot(x, y, color = 'gray')

        # Draw fiducial bounding polygon.
        if isinstance(fovUnion, shapely.geometry.MultiPolygon):
            tmp = fovUnion
        else:
            tmp = [fovUnion]

        for poly in tmp:
            coords = poly.exterior.coords.xy
            x = list(coords[0])
            y = list(coords[1])
            plt.plot(x, y, color = 'black')
            
        # Draw fiducials.
        plt.scatter(fidX, fidY, marker = 'x')
            
        ax = plt.gca()
        ax.set_aspect('equal', 'datalim')
            
        plt.title("fiducials")
        plt.xlabel("pixels")
        plt.ylabel("pixels")

        fname = "fiducials.pdf"
        fig.savefig(os.path.join(self.get_path(), fname),
                    format='pdf',
                    dpi=100)
        
        plt.close()

