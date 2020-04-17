#!/usr/bin/env python
"""
MERFISH simulator core.
"""
import importlib
import json
import os

import mersim
import mersim.merlin_interface as mint


class SimulationParameters(object):
    """
    Stores all of the simulation parameters.
    """    
    def __init__(self,
                 codebook = None,
                 dataOrganization = None,
                 microscope = None,
                 positions = None,
                 **kwds):
        super().__init__(**kwds)

        self._codebook = mint.Codebook(codebook)
        self._dataOrganization = mint.DataOrganization(dataOrganization)
        self._microscope = mint.Microscope(microscope)
        self._positions = mint.Positions(positions)

    def get_codebook(self):
        return self._codebook

    def get_data_organization(self):
        return self._dataOrganization

    def get_imaging_rounds(self):
        return self._dataOrganization.get_imaging_rounds()

    def get_microscope(self):
        return self._microscope

    def get_number_positions(self):
        return self._positions.get_number_positions()

    def get_positions(self):
        return self._positions

    
def simulate(config, simParams, dataPath):
    """
    Performs the simulation.
    """
    # Create simulator objects from config file.
    #
    for elt in config:
        a_module = importlib.import_module(config[elt]["module"])
        a_class = getattr(a_module, config[elt]["task"])
        parameters = {}
        if "parameters" in config[elt]:
            parameters = config[elt]["parameters"]
        a_object = a_class(dataPath = dataPath,
                           parameters = parameters,
                           task_name = elt)
        config[elt] = a_object

        
    # Simulation initializations.
    #
    
    # Layout sample.
    config["layout_sample"].layout(config, simParams)

    # Layout barcodes.
    config["layout_barcodes"].layout(config, simParams)

    # Barcode intensities.
    config["barcode_intensity"].intensity(config, simParams)

    # Layout fiducials.
    config["layout_fiducials"].layout(config, simParams)

    # Fiducial intensities.
    config["fiducial_intensity"].intensity(config, simParams)

    
    # Simulated movie creation.
    #
    
    # Iterate over positions.
    for fov in range(simParams.get_number_positions()):

        # Iterate over imaging round.
        for iRound in simParams.get_imaging_rounds():

            createMovie(config, simParams, fov, iRound)

    
