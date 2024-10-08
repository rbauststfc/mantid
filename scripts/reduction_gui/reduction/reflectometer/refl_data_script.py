# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
"""
Legacy class that old LR reduction options.
This is still in use for backward compatibility.
"""

import xml.dom.minidom
from reduction_gui.reduction.scripter import BaseScriptElement


class DataSets(BaseScriptElement):
    DataPeakSelectionType = "narrow"
    DataPeakPixels = [120, 130]
    DataPeakDiscreteSelection = "N/A"
    DataBackgroundFlag = False
    DataBackgroundRoi = [115, 137, 123, 137]
    DataTofRange = [9600.0, 21600.0]
    TofRangeFlag = True

    data_x_range_flag = True
    data_x_range = [115, 210]

    tthd_value = "N/A"
    ths_value = "N/A"

    norm_x_range_flag = True
    norm_x_range = [115, 210]

    NormFlag = True
    NormPeakPixels = [120, 130]
    NormBackgroundFlag = False
    NormBackgroundRoi = [115, 137]

    # Data files
    # data_files = [66421]
    # norm_file = 66196
    data_files = [0]
    norm_file = 0

    # Q range
    q_min = 0.001
    q_step = 0.001
    auto_q_binning = False

    # Angle offset
    angle_offset = 0.0
    angle_offset_error = 0.0

    # scaling factor file
    scaling_factor_file = ""
    scaling_factor_file_flag = True
    slits_width_flag = True

    # geometry correction
    geometry_correction_switch = False

    # incident medium list and selected value
    incident_medium_list = ["H2O"]
    incident_medium_index_selected = 0

    # 4th column of ASCII file (precision)
    fourth_column_flag = True
    fourth_column_dq0 = 0.0009
    fourth_column_dq_over_q = 0.045

    # how to treat overlap values
    overlap_lowest_error = True
    overlap_mean_value = False

    # Primary fraction range
    clocking_from = None
    clocking_to = None

    def __init__(self):
        super(DataSets, self).__init__()
        self.reset()

    # pylint: disable = unused-argument, arguments-differ
    def to_script(self, for_automated_reduction=False):
        """
        Generate reduction script
        @param execute: if true, the script will be executed
        """
        raise RuntimeError("refl_data_script.DataSets.to_script is deprecated")

    def update(self):
        """
        Update transmission from reduction output
        """
        pass

    def to_xml(self):
        """
        Create XML from the current data.
        """
        _xml = "<RefLData>\n"
        _xml += "<peak_selection_type>%s</peak_selection_type>\n" % self.DataPeakSelectionType
        _xml += "<from_peak_pixels>%s</from_peak_pixels>\n" % str(self.DataPeakPixels[0])
        _xml += "<to_peak_pixels>%s</to_peak_pixels>\n" % str(self.DataPeakPixels[1])
        _xml += "<peak_discrete_selection>%s</peak_discrete_selection>\n" % self.DataPeakDiscreteSelection
        _xml += "<background_flag>%s</background_flag>\n" % str(self.DataBackgroundFlag)
        _xml += "<back_roi1_from>%s</back_roi1_from>\n" % str(self.DataBackgroundRoi[0])
        _xml += "<back_roi1_to>%s</back_roi1_to>\n" % str(self.DataBackgroundRoi[1])
        _xml += "<back_roi2_from>%s</back_roi2_from>\n" % str(self.DataBackgroundRoi[2])
        _xml += "<back_roi2_to>%s</back_roi2_to>\n" % str(self.DataBackgroundRoi[3])
        _xml += "<tof_range_flag>%s</tof_range_flag>\n" % str(self.TofRangeFlag)
        _xml += "<from_tof_range>%s</from_tof_range>\n" % str(self.DataTofRange[0])
        _xml += "<to_tof_range>%s</to_tof_range>\n" % str(self.DataTofRange[1])
        _xml += "<data_sets>%s</data_sets>\n" % ",".join([str(i) for i in self.data_files])
        _xml += "<x_min_pixel>%s</x_min_pixel>\n" % str(self.data_x_range[0])
        _xml += "<x_max_pixel>%s</x_max_pixel>\n" % str(self.data_x_range[1])
        _xml += "<x_range_flag>%s</x_range_flag>\n" % str(self.data_x_range_flag)

        _xml += "<tthd_value>%s</tthd_value>\n" % str(self.tthd_value)
        _xml += "<ths_value>%s</ths_value>\n" % str(self.ths_value)

        _xml += "<norm_flag>%s</norm_flag>\n" % str(self.NormFlag)
        _xml += "<norm_x_range_flag>%s</norm_x_range_flag>\n" % str(self.norm_x_range_flag)
        _xml += "<norm_x_max>%s</norm_x_max>\n" % str(self.norm_x_range[1])
        _xml += "<norm_x_min>%s</norm_x_min>\n" % str(self.norm_x_range[0])

        _xml += "<norm_from_peak_pixels>%s</norm_from_peak_pixels>\n" % str(self.NormPeakPixels[0])
        _xml += "<norm_to_peak_pixels>%s</norm_to_peak_pixels>\n" % str(self.NormPeakPixels[1])
        _xml += "<norm_background_flag>%s</norm_background_flag>\n" % str(self.NormBackgroundFlag)
        _xml += "<norm_from_back_pixels>%s</norm_from_back_pixels>\n" % str(self.NormBackgroundRoi[0])
        _xml += "<norm_to_back_pixels>%s</norm_to_back_pixels>\n" % str(self.NormBackgroundRoi[1])
        _xml += "<norm_dataset>%s</norm_dataset>\n" % str(self.norm_file)

        # Q cut
        _xml += "<q_min>%s</q_min>\n" % str(self.q_min)
        _xml += "<q_step>%s</q_step>\n" % str(self.q_step)
        _xml += "<auto_q_binning>%s</auto_q_binning>\n" % str(self.auto_q_binning)
        _xml += "<overlap_lowest_error>%s</overlap_lowest_error>\n" % str(self.overlap_lowest_error)
        _xml += "<overlap_mean_value>%s</overlap_mean_value>\n" % str(self.overlap_mean_value)

        # Angle offset
        _xml += "<angle_offset>%s</angle_offset>\n" % str(self.angle_offset)
        _xml += "<angle_offset_error>%s</angle_offset_error>\n" % str(self.angle_offset_error)

        # scaling factor file name
        _xml += "<scaling_factor_flag>%s</scaling_factor_flag>\n" % str(self.scaling_factor_file_flag)
        _xml += "<scaling_factor_file>%s</scaling_factor_file>\n" % str(self.scaling_factor_file)
        _xml += "<slits_width_flag>%s</slits_width_flag>\n" % str(self.slits_width_flag)

        # geometry correction
        _xml += "<geometry_correction_switch>%s</geometry_correction_switch>\n" % str(self.geometry_correction_switch)

        # incident medium
        _xml += "<incident_medium_list>%s</incident_medium_list>\n" % str(self.incident_medium_list[0])
        _xml += "<incident_medium_index_selected>%s</incident_medium_index_selected>\n" % str(self.incident_medium_index_selected)

        # fourth column precision
        _xml += "<fourth_column_flag>%s</fourth_column_flag>\n" % str(self.fourth_column_flag)
        _xml += "<fourth_column_dq0>%s</fourth_column_dq0>\n" % str(self.fourth_column_dq0)
        _xml += "<fourth_column_dq_over_q>%s</fourth_column_dq_over_q>\n" % str(self.fourth_column_dq_over_q)

        # Primary fraction
        if self.clocking_from is not None and self.clocking_to is not None:
            _xml += "<clocking_from>%s</clocking_from>\n" % str(self.clocking_from)
            _xml += "<clocking_to>%s</clocking_to>\n" % str(self.clocking_to)

        _xml += "</RefLData>\n"

        return _xml

    def from_xml(self, xml_str):
        self.reset()
        dom = xml.dom.minidom.parseString(xml_str)
        self.from_xml_element(dom)
        # element_list = dom.getElementsByTagName("RefLData")
        # if len(element_list)>0:
        #    instrument_dom = element_list[0]

    def from_xml_element(self, instrument_dom):
        """
        Read in data from XML
        @param xml_str: text to read the data from
        """
        # Peak selection
        self.DataPeakSelectionType = BaseScriptElement.getStringElement(instrument_dom, "peak_selection_type")

        # Peak from/to pixels
        self.DataPeakPixels = [
            BaseScriptElement.getIntElement(instrument_dom, "from_peak_pixels"),
            BaseScriptElement.getIntElement(instrument_dom, "to_peak_pixels"),
        ]

        # data metadata
        _tthd_value = BaseScriptElement.getStringElement(instrument_dom, "tthd_value")
        if _tthd_value == "":
            _tthd_value = "N/A"
        self.tthd_value = _tthd_value

        _ths_value = BaseScriptElement.getStringElement(instrument_dom, "ths_value")
        if _ths_value == "":
            _ths_value = "N/A"
        self.ths_value = _ths_value

        # low resolution range
        self.data_x_range_flag = BaseScriptElement.getBoolElement(instrument_dom, "x_range_flag", default=DataSets.data_x_range_flag)

        self.data_x_range = [
            BaseScriptElement.getIntElement(instrument_dom, "x_min_pixel"),
            BaseScriptElement.getIntElement(instrument_dom, "x_max_pixel"),
        ]

        self.norm_x_range_flag = BaseScriptElement.getBoolElement(instrument_dom, "norm_x_range_flag", default=DataSets.norm_x_range_flag)

        self.norm_x_range = [
            BaseScriptElement.getIntElement(instrument_dom, "norm_x_min"),
            BaseScriptElement.getIntElement(instrument_dom, "norm_x_max"),
        ]

        # discrete selection string
        self.DataPeakDiscreteSelection = BaseScriptElement.getStringElement(instrument_dom, "peak_discrete_selection")

        # background flag
        self.DataBackgroundFlag = BaseScriptElement.getBoolElement(instrument_dom, "background_flag", default=DataSets.DataBackgroundFlag)

        # background from/to pixels
        self.DataBackgroundRoi = [
            BaseScriptElement.getIntElement(instrument_dom, "back_roi1_from"),
            BaseScriptElement.getIntElement(instrument_dom, "back_roi1_to"),
            BaseScriptElement.getIntElement(instrument_dom, "back_roi2_from"),
            BaseScriptElement.getIntElement(instrument_dom, "back_roi2_to"),
        ]

        # from TOF and to TOF
        self.TofRangeFlag = BaseScriptElement.getBoolElement(instrument_dom, "tof_range_flag", default=DataSets.TofRangeFlag)
        self.DataTofRange = [
            BaseScriptElement.getFloatElement(instrument_dom, "from_tof_range"),
            BaseScriptElement.getFloatElement(instrument_dom, "to_tof_range"),
        ]

        self.data_files = BaseScriptElement.getIntList(instrument_dom, "data_sets")

        # with or without norm
        self.NormFlag = BaseScriptElement.getBoolElement(instrument_dom, "norm_flag", default=DataSets.NormFlag)

        # Peak from/to pixels
        self.NormPeakPixels = [
            BaseScriptElement.getIntElement(instrument_dom, "norm_from_peak_pixels"),
            BaseScriptElement.getIntElement(instrument_dom, "norm_to_peak_pixels"),
        ]

        # background flag
        self.NormBackgroundFlag = BaseScriptElement.getBoolElement(
            instrument_dom, "norm_background_flag", default=DataSets.NormBackgroundFlag
        )

        # background from/to pixels
        self.NormBackgroundRoi = [
            BaseScriptElement.getIntElement(instrument_dom, "norm_from_back_pixels"),
            BaseScriptElement.getIntElement(instrument_dom, "norm_to_back_pixels"),
        ]

        self.norm_file = BaseScriptElement.getIntElement(instrument_dom, "norm_dataset")

        # Q cut
        self.q_min = BaseScriptElement.getFloatElement(instrument_dom, "q_min", default=DataSets.q_min)
        self.q_step = BaseScriptElement.getFloatElement(instrument_dom, "q_step", default=DataSets.q_step)
        self.auto_q_binning = BaseScriptElement.getBoolElement(instrument_dom, "auto_q_binning", default=False)

        # overlap_lowest_error
        self.overlap_lowest_error = BaseScriptElement.getBoolElement(instrument_dom, "overlap_lowest_error", default=True)
        self.overlap_mean_value = BaseScriptElement.getBoolElement(instrument_dom, "overlap_mean_value", default=False)

        # Angle offset
        self.angle_offset = BaseScriptElement.getFloatElement(instrument_dom, "angle_offset", default=DataSets.angle_offset)
        self.angle_offset_error = BaseScriptElement.getFloatElement(
            instrument_dom, "angle_offset_error", default=DataSets.angle_offset_error
        )

        # scaling factor file and options
        self.scaling_factor_file = BaseScriptElement.getStringElement(instrument_dom, "scaling_factor_file")
        self.slits_width_flag = BaseScriptElement.getBoolElement(instrument_dom, "slits_width_flag")
        self.scaling_factor_file_flag = BaseScriptElement.getBoolElement(instrument_dom, "scaling_factor_flag")

        # geometry correction switch
        self.geometry_correction_switch = BaseScriptElement.getBoolElement(instrument_dom, "geometry_correction_switch")

        # incident medium selected
        if BaseScriptElement.getStringList(instrument_dom, "incident_medium_list") != []:
            self.incident_medium_list = BaseScriptElement.getStringList(instrument_dom, "incident_medium_list")
            self.incident_medium_index_selected = BaseScriptElement.getIntElement(instrument_dom, "incident_medium_index_selected")
        else:
            self.incident_medium_list = ["H2O"]
            self.incident_medium_index_selected = 0

        # fourth column (precision)
        self.fourth_column_flag = BaseScriptElement.getBoolElement(instrument_dom, "fourth_column_flag")
        self.fourth_column_dq0 = BaseScriptElement.getFloatElement(instrument_dom, "fourth_column_dq0")
        self.fourth_column_dq_over_q = BaseScriptElement.getFloatElement(instrument_dom, "fourth_column_dq_over_q")

        # Primary fraction
        self.clocking_from = BaseScriptElement.getIntElement(instrument_dom, "clocking_from", default=None)
        self.clocking_to = BaseScriptElement.getIntElement(instrument_dom, "clocking_to", default=None)

    def reset(self):
        """
        Reset state
        """
        self.DataPeakSelectionType = DataSets.DataPeakSelectionType
        self.DataBackgroundFlag = DataSets.DataBackgroundFlag
        self.DataBackgroundRoi = DataSets.DataBackgroundRoi
        self.DataPeakDiscreteSelection = DataSets.DataPeakDiscreteSelection
        self.DataPeakPixels = DataSets.DataPeakPixels
        self.DataTofRange = DataSets.DataTofRange
        self.TofRangeFlag = DataSets.TofRangeFlag
        self.data_files = DataSets.data_files

        self.tthd_value = DataSets.tthd_value
        self.ths_value = DataSets.ths_value

        self.NormFlag = DataSets.NormFlag
        self.NormBackgroundFlag = DataSets.NormBackgroundFlag
        self.NormBackgroundRoi = DataSets.NormBackgroundRoi
        self.NormPeakPixels = DataSets.NormPeakPixels
        self.norm_file = DataSets.norm_file
        self.data_x_range_flag = DataSets.data_x_range_flag
        self.data_x_range = DataSets.data_x_range
        self.norm_x_range_flag = DataSets.norm_x_range_flag
        self.norm_x_range = DataSets.norm_x_range

        # Q range
        self.q_min = DataSets.q_min
        self.q_step = DataSets.q_step
        self.auto_q_binning = DataSets.auto_q_binning

        # overlap_lowest_error
        self.overlap_lowest_error = DataSets.overlap_lowest_error
        self.overlap_mean_value = DataSets.overlap_mean_value

        # Angle offset
        self.angle_offset = DataSets.angle_offset
        self.angle_offset_error = DataSets.angle_offset_error

        # scaling factor file and options
        self.scaling_factor_file = DataSets.scaling_factor_file
        self.slits_width_flag = DataSets.slits_width_flag
        self.scaling_factor_file_flag = DataSets.scaling_factor_file_flag

        # geometry correction
        self.geometry_correction_switch = DataSets.geometry_correction_switch

        # incident medium selected
        self.incident_medium_list = DataSets.incident_medium_list
        self.incident_medium_index_selected = DataSets.incident_medium_index_selected

        # 4th column (precision)
        self.fourth_column_flag = DataSets.fourth_column_flag
        self.fourth_column_dq0 = DataSets.fourth_column_dq0
        self.fourth_column_dq_over_q = DataSets.fourth_column_dq_over_q

        # Primary fraction
        self.clocking_from = DataSets.clocking_from
        self.clocking_to = DataSets.clocking_to
