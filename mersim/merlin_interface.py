#!/usr/bin/env python
"""
Classes that use and/or depend on MERlin.
"""
import csv
import json
import os
import numpy as np
import pandas

import merlin
import merlin.data.dataorganization as merlin_do
import merlin.data.codebook as merlin_co

import mersim


def getDataDirectory(filePath):
    """
    Return full path to data.
    """
    if os.path.exists(filePath):
        return filePath

    else:
        filePath = os.sep.join([merlin.DATA_HOME, filePath])
        if not os.path.exists(filePath):
            os.mkdir(filePath)

        return filePath


def loadParameters(filePath):
    """
    Return full path to data.
    """
    if not os.path.exists(filePath):
        filePath = os.path.join(mersim.SIMULATION_PARAMETERS_HOME, filePath)

    with open(filePath) as fp:
        config = json.load(fp)

    return config


class Codebook(merlin_co.Codebook):
    """
    MERlin Codebook like class.
    """
    def __init__(self, filePath,
                 codebookIndex: int = 0,
                 codebookName: str = None):    
        """
        Initialization without MERlin DataSet object.
        """
        if not os.path.exists(filePath):
            filePath = os.sep.join([merlin.CODEBOOK_HOME, filePath])

        newVersion = True
        with open(filePath, 'r') as f:
            if 'version' in f.readline():
                newVersion = False

        if newVersion:
            self._data = pandas.read_csv(filePath)
        else:
            headerLength = 3
            barcodeData = pandas.read_csv(
                filePath, header=headerLength, skipinitialspace=True,
                usecols=['name', 'id', 'barcode'],
                converters={'barcode': merlin_co._parse_barcode_from_string})
            with open(filePath, 'r') as inFile:
                csvReader = csv.reader(inFile, delimiter=',')
                header = [row for i, row in enumerate(csvReader)
                          if i < headerLength]

            bitNames = [x.strip() for x in header[2][1:]]

            self._data = self._generate_codebook_dataframe(
                    barcodeData, bitNames)

        if not codebookName:
            codebookName = os.path.splitext(os.path.basename(filePath))[0]
        self._codebookName = codebookName
        self._codebookIndex = codebookIndex


class DataOrganization(merlin_do.DataOrganization):
    """
    MERlin DataOrganization like class.
    """
    def __init__(self, filePath):
        """
        Initialization without MERlin DataSet object.
        """
        if not os.path.exists(filePath):
            filePath = os.sep.join(
                [merlin.DATA_ORGANIZATION_HOME, filePath])

        self.data = pandas.read_csv(
            filePath,
            converters={'frame': merlin_do._parse_int_list,
                        'zPos': merlin_do._parse_list})
        self.data['readoutName'] = self.data['readoutName'].str.strip()

        stringColumns = ['readoutName', 'channelName', 'imageType',
                         'imageRegExp', 'fiducialImageType', 'fiducialRegExp']
        self.data[stringColumns] = self.data[stringColumns].astype('str')

    def get_frame_description(self, iRound, fi):
        """
        Return description of what a particular frame contains.
        """
        df = self.data.loc[self.data['imagingRound'] == iRound]

        # Check if it's a fiducial.
        if (fi == df['fiducialFrame'].iloc[0]):
            color = df['fiducialColor'].iloc[0]
            return ['fiducial', color]

        # Look for fi in 'frame' field.
        allZ = self.get_z_positions()
        for index, row in df.iterrows():
            if fi in row['frame']:
                zi = np.where(row['frame'] == fi)[0][0]
                zPos = row['zPos'][zi]
                zV = np.where(allZ == zPos)[0][0]
                return [row['channelName'], row['color'], row['bitNumber'], zV, zPos]

        return ["blank"]
                        
    def get_imaging_rounds(self):
        """
        Return the imaging rounds in this simulation.
        """
        return sorted(np.unique(self.data['imagingRound'].to_numpy()))

    def get_image_type(self, iRound):
        """
        Return 'imageType'.
        """
        df = self.data.loc[self.data['imagingRound'] == iRound]
        return df['imageType'].iloc[0]
        
    def get_number_frames(self, iRound):
        """
        Return the number of frames in this imaging round.
        """
        df = self.data.loc[self.data['imagingRound'] == iRound]
        max_frame = np.amax([y for x in df['frame'] for y in x])
        fiducial = df['fiducialFrame'].iloc[0]
        return max(max_frame, fiducial) + 1


class Positions(object):
    """
    Provides positions of simulated images.
    """
    def __init__(self, filePath, umToPix):
        if not os.path.exists(filePath):
            filePath = os.sep.join(
                [merlin.POSITION_HOME, filePath])

        self.umToPix = umToPix

        self.data = pandas.read_csv(filePath,
                                    names = ["x", "y"])

        self.data = self.data * self.umToPix

    def get_fov_xy(self, fov):
        return self.data.loc[fov,:].values.tolist()

    def get_fov_xy_um(self, fov):
        return [x/self.umToPix for x in self.data.loc[fov,:].values.tolist()]

    def get_number_positions(self):
        return self.data.shape[0]
    
class Microscope(object):
    """
    Provides information about the microscope.
    """
    def __init__(self, filePath):
        if not os.path.exists(filePath):
            filePath = os.sep.join(
                [merlin.MICROSCOPE_PARAMETERS_HOME, filePath])

        with open(filePath) as fp:
            self.microscopeParameters = json.load(fp)

    def get_image_dimensions(self):
        return self.microscopeParameters["image_dimensions"]

    def get_microns_per_pixel(self):
        return self.microscopeParameters["microns_per_pixel"]

    
