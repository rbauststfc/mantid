# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
# std imports
import sys

from mantidqt.utils.qt import import_qt
from matplotlib.image import AxesImage

from .lineplots import CursorTracker, cursor_info

# Constants
DBLMAX = sys.float_info.max

ImageInfoWidget = import_qt('.._common', 'mantidqt.widgets', 'ImageInfoWidget')


class ImageInfoTracker(CursorTracker):
    def __init__(self, image: AxesImage, transpose_xy: False, widget: ImageInfoWidget):
        """
        Update the image that the widget refers too.
        :param: An AxesImage instance to track
        :param: transpose_xy: If true the cursor position should be transposed
                before sending to the table update
        """
        super().__init__(image_axes=image.axes, autoconnect=False)
        self._image = image
        self._transpose_xy = transpose_xy
        self._widget = widget

    def on_cursor_at(self, xdata: float, ydata: float):
        """
        Update the image table for the given coordinates
        :param xdata: X coordinate of cursor in data space
        :param ydata: Y coordinate of cursor in data space
        """
        if self._image is None:
            return

        cinfo = cursor_info(self._image, xdata, ydata)
        if cinfo is not None:
            arr, _, (i, j) = cinfo
            if (0 <= i < arr.shape[0]) and (0 <= j < arr.shape[1]):
                if self._transpose_xy:
                    ydata, xdata = xdata, ydata
                self._widget.updateTable(xdata, ydata, arr[i, j])

    def on_cursor_outside_axes(self):
        """Update the image table given the mouse has moved out of the image axes"""
        self._widget.updateTable(DBLMAX, DBLMAX, DBLMAX)
