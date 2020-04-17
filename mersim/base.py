#!/usr/bin/env python
"""
Simulation base task.
"""
import os
import pickle


class SimulationBase(object):

    def __init__(self,
                 dataPath = None,
                 parameters = None,
                 taskname = None,
                 **kwds):

        super().__init__(**kwds)

        self._parameters = parameters

        self._filePath = os.path.join(dataPath, taskname)

    def getFOVName(self, fov = None):
        if fov is None:
            return "data.bin"
        else:
            return "fov_{0:d}.bin".format(fov)
        
    def load_data(self, fov = None):
        fname = self.getFOVName(fov)
        with open(fname, "rb") as fp:
            a_object = pickle.load(fp)
        return a_object
        
    def save_data(self, a_object, fov = None):
        fname = self.getFOVName(fov)
        with open(fname, "wb") as fp:
            pickle.dump(a_object, fp)

