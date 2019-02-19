# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
# pylint: disable=no-init,too-many-instance-attributes
from __future__ import (absolute_import, division, print_function)
from mantid.api import (AlgorithmFactory, AnalysisDataService, MatrixWorkspace, MatrixWorkspaceProperty,
                        PythonAlgorithm, PropertyMode, WorkspaceGroup)
from mantid.kernel import Direction


def string_ends_with(string, delimiter):
    delimiter_length = len(delimiter)
    return string[-delimiter_length:] == delimiter if len(string) > delimiter_length else False


def exists_in_ads(workspace_name):
    return AnalysisDataService.doesExist(workspace_name)


def get_ads_workspace(workspace_name):
    return AnalysisDataService.retrieve(workspace_name) if exists_in_ads(workspace_name) else None


def contains_workspace(group, workspace):
    return group.contains(workspace) if isinstance(group, WorkspaceGroup) else False


def filter_by_contents(workspace_names, workspace):
    return [name for name in workspace_names if contains_workspace(get_ads_workspace(name), workspace)]


def filter_by_name_end(workspace_names, delimiter):
    return [name for name in workspace_names if string_ends_with(name, delimiter)]


def find_result_group_containing(workspace_name, group_extension):
    workspace_names = AnalysisDataService.Instance().getObjectNames()
    result_groups = filter_by_name_end(workspace_names, group_extension)
    groups = filter_by_contents(result_groups, workspace_name)
    return groups[0] if groups else None


def is_equal(a, b, tolerance):
    return abs(a - b) <= tolerance * max(abs(a), abs(b))


def get_bin_index_of_value(workspace, value):
    x_axis = workspace.getAxis(0)
    for i in range(0, x_axis.length()):
        if is_equal(x_axis.getValue(i), value, 0.000001):
            return i

    raise ValueError('The corresponding bin in the input workspace could not be found.')


def get_x_insertion_index(input_workspace, single_fit_workspace):
    single_fit_x_axis = single_fit_workspace.getAxis(0)
    bin_value = float(single_fit_x_axis.label(0))
    return get_bin_index_of_value(input_workspace, bin_value)


class IndirectReplaceFitResult(PythonAlgorithm):
    _input_workspace = None
    _single_fit_workspace = None
    _output_workspace = None

    _end_row = None
    _column_index = None
    _insertion_x_index = None

    _result_group = None
    _allowed_extension = '_Result'

    def category(self):
        return "Workflow\\DataHandling;Inelastic\\Indirect"

    def summary(self):
        return 'Replaces a fit result within the Input Workspace with the corresponding fit result found in the ' \
               'Single Fit Workspace.'

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty('InputWorkspace', '',
                                                     optional=PropertyMode.Mandatory,
                                                     direction=Direction.Input),
                             doc='The result workspace containing the poor fit value which needs replacing. It\'s name '
                                 'must end with _Result.')

        self.declareProperty(MatrixWorkspaceProperty('SingleFitWorkspace', '',
                                                     optional=PropertyMode.Mandatory,
                                                     direction=Direction.Input),
                             doc='The result workspace containing the result data from a single fit. It\'s name must '
                                 'end with _Result.')

        self.declareProperty(MatrixWorkspaceProperty('OutputWorkspace', '',
                                                     direction=Direction.Output,
                                                     optional=PropertyMode.Optional),
                             doc='The name of the output workspace.')

    def validateInputs(self):
        issues = dict()

        input_name = self.getPropertyValue('InputWorkspace')
        single_fit_name = self.getPropertyValue('SingleFitWorkspace')
        output_name = self.getPropertyValue('OutputWorkspace')

        input_workspace = self.getProperty('InputWorkspace').value
        single_fit_workspace = self.getProperty('SingleFitWorkspace').value

        if not string_ends_with(input_name, self._allowed_extension):
            issues['InputWorkspace'] = 'The input workspace must have a name ending in ' \
                                       '{0}'.format(self._allowed_extension)

        if not string_ends_with(single_fit_name, self._allowed_extension):
            issues['SingleFitWorkspace'] = 'This workspace must have a name ending in ' \
                                           '{0}'.format(self._allowed_extension)

        if not output_name:
            issues['OutputWorkspace'] = 'No OutputWorkspace name was provided.'

        if not isinstance(input_workspace, MatrixWorkspace):
            issues['InputWorkspace'] = 'The input workspace must be a matrix workspace.'
        else:
            input_x_axis = input_workspace.getAxis(0)
            if not input_x_axis.isNumeric():
                issues['InputWorkspace'] = 'The input workspace must have a numeric x axis.'
            if len(input_workspace.readY(0)) < 2:
                issues['InputWorkspace'] = 'The input workspace must contain result data from a fit involving 2 or ' \
                                           'more spectra.'

        if not isinstance(single_fit_workspace, MatrixWorkspace):
            issues['SingleFitWorkspace'] = 'The single fit workspace must be a matrix workspace.'
        else:
            single_fit_x_axis = single_fit_workspace.getAxis(0)
            if not single_fit_x_axis.isNumeric():
                issues['SingleFitWorkspace'] = 'The single fit workspace must have a numeric x axis.'
            if len(single_fit_workspace.readY(0)) > 1:
                issues['SingleFitWorkspace'] = 'The single fit workspace must contain data from a single fit.'

        if isinstance(input_workspace, MatrixWorkspace) and isinstance(single_fit_workspace, MatrixWorkspace):
            if input_workspace.getNumberHistograms() != single_fit_workspace.getNumberHistograms():
                issues['InputWorkspace'] = 'The input workspace and single fit workspace must have the same number ' \
                                           'of histograms.'

        return issues

    def _setup(self):
        self._input_workspace = self.getPropertyValue('InputWorkspace')
        self._single_fit_workspace = self.getPropertyValue('SingleFitWorkspace')
        self._output_workspace = self.getPropertyValue('OutputWorkspace')

        self._end_row = get_ads_workspace(self._single_fit_workspace).getNumberHistograms() - 1
        self._column_index = len(get_ads_workspace(self._single_fit_workspace).readY(0)) - 1
        self._insertion_x_index = get_x_insertion_index(get_ads_workspace(self._input_workspace),
                                                        get_ads_workspace(self._single_fit_workspace))

        self._result_group = find_result_group_containing(self._input_workspace, self._allowed_extension + 's')

    def PyExec(self):
        self._setup()
        self._copy_data()

        self.setProperty('OutputWorkspace', self._output_workspace)

        self._add_workspace_to_group()

    def _copy_data(self):
        copy_algorithm = self.createChildAlgorithm(name='CopyDataRange', startProgress=0.1,
                                                   endProgress=1.0, enableLogging=True)
        copy_algorithm.setAlwaysStoreInADS(True)

        args = {"InputWorkspace": self._single_fit_workspace, "DestWorkspace": self._input_workspace,
                "StartWorkspaceIndex": 0, "EndWorkspaceIndex": self._end_row, "XMin": self._column_index,
                "XMax": self._column_index, "InsertionYIndex": 0, "InsertionXIndex": self._insertion_x_index,
                "OutputWorkspace": self._output_workspace}

        for key, value in args.items():
            copy_algorithm.setProperty(key, value)

        copy_algorithm.execute()

    def _add_workspace_to_group(self):
        group_workspace = get_ads_workspace(self._result_group)
        if self._input_workspace == self._output_workspace:
            group_workspace.remove(self._input_workspace)
        group_workspace.addWorkspace(get_ads_workspace(self._output_workspace))


# Register algorithm with Mantid
AlgorithmFactory.subscribe(IndirectReplaceFitResult)
