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
                 taskName = None,
                 **kwds):

        super().__init__(**kwds)

        self._parameters = parameters
        self._taskName = taskName

        self._filePath = os.path.join(dataPath, taskName)

        if not os.path.exists(self._filePath):
            os.mkdir(self._filePath)

    def get_fov_name(self, fov = None):
        if fov is None:
            return "data.bin"
        else:
            return "fov_{0:d}.bin".format(fov)

    def get_path(self):
        return self._filePath
        
    def load_data(self, fov = None):
        fname = os.path.join(self._filePath, self.get_fov_name(fov))
        with open(fname, "rb") as fp:
            a_object = pickle.load(fp)
        return a_object

    def run_task(self, config, simParams):
        print("Running:", self._taskName)

    def save_data(self, a_object, fov = None):
        fname = os.path.join(self._filePath, self.get_fov_name(fov))
        with open(fname, "wb") as fp:
            pickle.dump(a_object, fp)
