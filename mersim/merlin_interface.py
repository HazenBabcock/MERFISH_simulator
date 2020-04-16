#!/usr/bin/env python
"""
Classes that use and/or depend on MERlin.

Hazen 04/20
"""
import csv
import json
import os
import pandas

import merlin
import merlin.data.dataorganization as merlin_do
import merlin.data.codebook as merlin_co


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


class Positions(object):
    """
    Provides positions of simulated images.
    """
    def __init__(self, filePath):
        if not os.path.exists(filePath):
            filePath = os.sep.join(
                [merlin.POSITION_HOME, filePath])

        self.data = pandas.read_csv(filePath,
                                    names = ["x", "y"])

    
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

    
