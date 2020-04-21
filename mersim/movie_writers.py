#!/usr/bin/env python
"""
Classes for saving movies.
"""
import numpy as np
import tifffile

from xml.dom import minidom
from xml.etree import ElementTree

import mersim.base as base


class ImageWriter(base.Base):

    def save_stack(self, filePath, stack, pos):
        """
        Create the XML file that MERlin will look at to determine the image position.
        """
        xml_out = ElementTree.Element("settings")
        acquisition_block = ElementTree.SubElement(xml_out, "acquisition")
        stage_position_block = ElementTree.SubElement(acquisition_block, "stage_position")
        stage_position_block.text = "{0:.2f},{1:.2f}".format(*pos)
        stage_position_block.set('type', 'custom')

        rough_string = ElementTree.tostring(xml_out, 'utf-8')
        reparsed = minidom.parseString(rough_string)

        with open(filePath + ".xml", "w") as fp:
            fp.write(reparsed.toprettyxml(indent="  ", encoding = "ISO-8859-1").decode())


class TiffWriter(ImageWriter):
    
    def save_stack(self, filePath, stack, pos):
        super().save_stack(filePath, stack, pos)

        with tifffile.TiffWriter(filePath + ".tif") as tf:
            for elt in stack:
                tmp = np.clip(elt, 0, 2**16-1)
                tf.save(elt.astype(np.uint16))
