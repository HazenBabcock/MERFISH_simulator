#!/usr/bin/env python
"""
Classes that describe the microscope PSF.
"""
import numpy as np
import os
import pickle

import microscPSF.microscPSF as ms_psf

import mersim.base as base


#
# PSFs
#

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


class GibsonLanniPSF(BasePSF):
    """
    Gibson-Lanni PSF.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.psfSize = self._parameters["size"]
        self.psfCenter = int(self.psfSize/2)
                    
        # PSF z spacing in microns.
        self.zSpacing = self._parameters["z_spacing"]

        # Set GL-PSF parameters.
        self.mp = ms_psf.m_params
        for elt in self._parameters:
            if elt in self.mp:
                self.mp[elt] = self._parameters[elt]

        # For PSF memoization.
        self.psfs = {}

    def get_psf(self, x, y, z, zPos, color):
        psf = self.psfs[color][str(zPos)]
        zv = int(z/self.zSpacing)
        if (zv >= psf.shape[0]):
            return [0, 0, None]
        else:
            return [int(np.round(x)) - self.psfCenter,
                    int(np.round(y)) - self.psfCenter,
                    psf[zv,:,:]]

    def initialize(self, config, simParams, color):
        if not color in self.psfs:
            umPerPix = simParams.get_microscope().get_microns_per_pixel()

            # Nominal wavelength in microns.
            #
            wvl = float(color)*0.001

            # A radial mask so that we don't need to have the PSF cover a huge
            # area in order not to be visible as a square when far out of
            # focus.
            #
            [mx, my] = np.mgrid[ -self.psfSize/2.0 : self.psfSize/2.0,
                                 -self.psfSize/2.0 : self.psfSize/2.0]
            k = np.sqrt(mx*mx + my*my)
            r_mask = np.ones((self.psfSize, self.psfSize))
            r_mask[(k > self.psfSize/2.0)] = 0.0
            
            # Figure out possible z positions.
            #
            # We assume that if there is only a single Z position it's at Z = 0.
            # Not sure if it is necessary but we'll always generate at least
            # two Z planes.
            #
            zPos = simParams.get_z_positions()
            deltaZ = self.zSpacing
            if (len(zPos) > 1):
                deltaZ = max(zPos[1] - zPos[0], deltaZ)
            maxZ = simParams.get_z_positions()[-1] + 1.1*deltaZ
            pz = np.arange(0.0, maxZ, self.zSpacing)

            # At each Z position compute that a point source X distance above
            # the coverslip would look like.
            #
            zPSFs = {}
            for elt in zPos:
                psfXYZ = ms_psf.gLXYZParticleScan(self.mp,
                                                  umPerPix,
                                                  self.psfSize,
                                                  pz,
                                                  wvl = wvl,
                                                  zv = -elt)
                for i in range(psfXYZ.shape[0]):
                    psfXYZ[i,:,:] = psfXYZ[i,:,:] * r_mask                    
                zPSFs[str(elt)] = psfXYZ
            self.psfs[color] = zPSFs

            # Save PSF.
            fname = os.path.join(self.get_path(), color + ".bin")
            with open(fname, "wb") as fp:
                pickle.dump(zPSFs, fp)


#
# Cameras
#
class IdealCamera(base.Base):

    def camera_image(self, image):        
        return np.random.poisson(image, image.shape) + self._parameters["offset"]


class CameraBasic(base.Base):

    def camera_image(self, image):
        image = np.random.poisson(image, image.shape).astype(float)
        image += np.random.normal(0.0, self._parameters["read_noise"], image.shape)
        return image + self._parameters["offset"]
