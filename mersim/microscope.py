#!/usr/bin/env python
"""
Classes that describe the microscope PSF.
"""
import numpy as np
import os
import pickle

import mersim.base as base


class BasePSF(base.Base):
    """
    Base class for PSFs.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)
        
        if not os.path.exists(self.get_path()):
            os.mkdir(self.get_path())

    
class GaussianPSF(BasePSF):
    """
    3D Gaussian PSF.

    This version is discretized in X/Y.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.psfSize = self._parameters["size"]
        self.psfCenter = int(self.psfSize/2)
                    
        # PSF z spacing in microns.
        self.zSpacing = self._parameters["z_spacing"]
        
        # For PSF memoization.
        self.psfs = {}

    def get_psf(self, x, y, z, zPos, color):
        psf = self.psfs[color]
        zv = int(abs(z - zPos)/self.zSpacing)
        if (zv >= psf.shape[0]):
            return [0, 0, None]
        else:
            return [int(np.round(x)) - self.psfCenter,
                    int(np.round(y)) - self.psfCenter,
                    psf[zv,:,:]]

    def initialize(self, config, simParams, color):
        if not color in self.psfs:
            umPerPix = simParams.get_microscope().get_microns_per_pixel()

            # We are assuming color is wavelength in nanometers. Calculate
            # PSF sigmas in microns.
            #
            sigmaX = (0.5e-3 * float(color) / (2.0 * self._parameters["NA"]))
            sigmaZ = (0.5e-3 * float(color) * 2.0/self._parameters["NA"])
            
            xv = (np.arange(self.psfSize) - self.psfCenter) * umPerPix
            exv = np.exp(-xv*xv/(2.0*sigmaX*sigmaX))
            exv = exv/np.sum(exv)
            
            zv = self.zSpacing * np.arange(self.psfSize)
            ezv = np.exp(-zv*zv/(2.0*sigmaZ*sigmaZ))
            
            psf = np.zeros((self.psfSize, self.psfSize, self.psfSize))
            for i in range(self.psfSize):
                psf[i,:,:] = ezv[i] * np.outer(exv, exv)

            self.psfs[color] = psf

            # Save PSF.
            fname = os.path.join(self.get_path(), color + ".bin")
            with open(fname, "wb") as fp:
                pickle.dump(psf, fp)


class IdealCamera(base.Base):

    def camera_image(self, image):        
        return np.random.poisson(image, image.shape) + self._parameters["offset"]


class CameraBasic(base.Base):

    def camera_image(self, image):
        image = np.random.poisson(image, image.shape).astype(np.float)
        image += np.random.normal(0.0, self._parameters["read_noise"], image.shape)
        return image + self._parameters["offset"]
