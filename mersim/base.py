#!/usr/bin/env python
"""
Simulation base task.
"""
import numpy as np
import os
import pickle


class Base(object):

    def __init__(self,
                 dataPath = None,
                 parameters = None,
                 taskName = None,
                 **kwds):

        super().__init__(**kwds)

        self._parameters = parameters
        self._taskName = taskName

        self._filePath = os.path.join(dataPath, taskName)

    def get_path(self):
        return self._filePath

    def get_parameter(self, pname, default = None):
        if default is None:
            return self._parameters[pname]
        else:
            if pname in self._parameters:
                return self._parameters[pname]
            else:
                return default


class ImageBase(Base):
    """
    Base class for image creation.
    """
    def background(self, config, simParams, fov, iRound, desc):
        """
        Create background for an image.
        """
        return None

    def foreground(self, config, simParams, fov, iRound, desc):
        """
        Create foreground for an image.

        The default behavior is to return image that is all zeros.
        """
        imageSize = simParams.get_microscope().get_image_dimensions()
        return np.zeros(imageSize, dtype = np.float)
        
    def make_image(self, config, simParams, fov, iRound, desc):
        """
        Depending desc[0], which is the 'channelName' field call the
        foreground() or background() method as appropriate.
        """
        # Remove any numbers from 'channelName' field.
        imtype = ''.join([i for i in desc[0] if not i.isdigit()])
        if self._taskName.startswith(imtype):
            return self.foreground(config, simParams, fov, iRound, desc)
        else:
            return self.background(config, simParams, fov, iRound, desc)


class SimulationBase(Base):
    """
    Base class for laying out different aspects of the simulation.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)
        
        if not os.path.exists(self.get_path()):
            os.mkdir(self.get_path())

    def get_fov_name(self, fov, zi):
        if fov is None:
            return "data.bin"
        else:
            if zi is None:
                return "fov_{0:d}.bin".format(fov)
            else:
                return "fov_{0:d}_{1:d}.bin".format(fov, zi)

    def load_data(self, fov = None, zi = None):
        fname = os.path.join(self._filePath, self.get_fov_name(fov, zi))
        if not os.path.exists(fname):
            return None

        with open(fname, "rb") as fp:
            aObject = pickle.load(fp)
        return aObject

    def run_task(self, config, simParams):
        print()
        print("Running:", self._taskName)

    def save_data(self, aObject, fov = None, zi = None):
        fname = os.path.join(self.get_path(), self.get_fov_name(fov, zi))
        with open(fname, "wb") as fp:
            pickle.dump(aObject, fp)
