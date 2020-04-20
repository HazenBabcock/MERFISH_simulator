#!/usr/bin/env python
"""
MERFISH simulator core.
"""
import argparse
import importlib
import os
import shapely
import shapely.geometry

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
        self._positions = mint.Positions(positions,
                                         1.0/self._microscope.get_microns_per_pixel())

    def get_codebook(self):
        return self._codebook

    def get_data_organization(self):
        return self._dataOrganization

    def get_fov_origin(self, fov):
        [px, py] = self.get_fov_xy(fov)
        fovSize = self._microscope.get_image_dimensions()
        sx = 0.5*fovSize[0]
        sy = 0.5*fovSize[1]        
        return [px - sx, py - sy]
        
    def get_fov_rect(self, fov):
        [px, py] = self.get_fov_xy(fov)
        fovSize = self._microscope.get_image_dimensions()
        sx = 0.5*fovSize[0]
        sy = 0.5*fovSize[1]

        return shapely.geometry.Polygon([[px - sx, py - sy],
                                         [px - sx, py + sy],
                                         [px + sx, py + sy],
                                         [px + sx, py - sy]])

    def get_fov_xy(self, fov):
        return self._positions.get_fov_xy(fov)

    def get_imaging_rounds(self):
        return self._dataOrganization.get_imaging_rounds()

    def get_microscope(self):
        return self._microscope

    def get_number_barcodes(self):
        return self._codebook.get_barcodes().shape[0]

    def get_number_positions(self):
        return self._positions.get_number_positions()

    def get_number_z(self):
        return len(self._dataOrganization.get_z_positions())
        
    def get_positions(self):
        return self._positions


def build_parser():
    parser = argparse.ArgumentParser(description='Simulate MERFISH data.')

    parser.add_argument('dataset',
                        help='directory where the simulation data is stored')
    parser.add_argument('-a', '--simulation-parameters',
                        help='name of the simulation parameters file to use')
    parser.add_argument('-o', '--data-organization',
                        help='name of the data organization file to use')
    parser.add_argument('-c', '--codebook',
                        help='name of the codebook to use')
    parser.add_argument('-m', '--microscope-parameters',
                        help='name of the microscope parameters to use')
    parser.add_argument('-p', '--positions',
                        help='name of the position file to use')
    return parser

    
def simulate(config, simParams, dataPath):
    """
    Performs the simulation.
    """
    # Create simulator objects from config file.
    #
    for elt in config:
        aModule = importlib.import_module(config[elt]["module"])
        aClass = getattr(aModule, config[elt]["task"])
        parameters = {}
        if "parameters" in config[elt]:
            parameters = config[elt]["parameters"]
        aObject = aClass(dataPath = dataPath,
                         parameters = parameters,
                         taskName = elt)
        config[elt] = aObject

    # Simulation initializations.
    #
    
    # Layout sample.
    config["layout_sample"].run_task(config, simParams)

    # Layout barcodes.
    config["layout_barcodes"].run_task(config, simParams)

    # Barcode intensities.
    config["barcode_intensity"].run_task(config, simParams)

    # Layout fiducials.
#    config["layout_fiducials"].run_task(config, simParams)

    # Fiducial intensities.
#    config["fiducial_intensity"].run_task(config, simParams)

    
    # Simulated movie creation.
    #
    
    # Iterate over positions.
#    for fov in range(simParams.get_number_positions()):

        # Iterate over imaging round.
#        for iRound in simParams.get_imaging_rounds():

#            createMovie(config, simParams, fov, iRound)
   

if (__name__ == "__main__"):

    parser = build_parser()
    args = parser.parse_args()

    config = mint.loadParameters(args.simulation_parameters)
    dataPath = mint.getDataDirectory(args.dataset)
    simParams = SimulationParameters(codebook = args.codebook,
                                     dataOrganization = args.data_organization,
                                     microscope = args.microscope_parameters,
                                     positions = args.positions)

    simulate(config, simParams, dataPath)
