# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

import numbers
import numpy as np
from typing import Dict

from mantid.api import AlgorithmFactory, PythonAlgorithm, Progress
from mantid.api import WorkspaceFactory, AnalysisDataService

# noinspection PyProtectedMember
from mantid.simpleapi import GroupWorkspaces
import abins
from abins.abinsalgorithm import AbinsAlgorithm


# noinspection PyPep8Naming,PyMethodMayBeStatic
class Abins2D(PythonAlgorithm, AbinsAlgorithm):
    _ab_initio_program = None
    _vibrational_or_phonon_data_file = None
    _temperature = None
    _bin_width = None
    _sample_form = None
    _instrument_name = None
    _setting = None
    _atoms = None
    _sum_contributions = None
    _save_ascii = None
    _scale_by_cross_section = None
    _calc_partial = None
    _out_ws_name = None
    _num_quantum_order_events = None
    _autoconvolution = None

    def category(self) -> str:
        return "Simulation"

    def summary(self) -> str:
        return "Calculates inelastic neutron scattering over 2D (q,ω) space."

    def PyInit(self) -> None:
        from abins.constants import TWO_DIMENSIONAL_INSTRUMENTS
        # Declare properties for all Abins Algorithms
        self.declare_common_properties()

        # Declare additional properties
        self.declareProperty(name="Autoconvolution", defaultValue=False,
                             doc="Estimate higher quantum orders by convolution with fundamental spectrum.")

        # Declare instrument properties
        self.declare_instrument_properties(
            default="TwoDMap", choices=TWO_DIMENSIONAL_INSTRUMENTS,
            multiple_choice_settings=[('Chopper', 'settings', 'Chopper package')],
            freeform_settings=[('IncidentEnergy', '4100', 'Incident energy in wavenumber'),
                               ('ChopperFrequency', '400', 'Chopper frequency in Hz')])

    def validateInputs(self) -> Dict[str,str]:
        """
        Performs input validation. Use to ensure the user has defined a consistent set of parameters.
        """
        issues = dict()
        issues = self.validate_common_inputs(issues)
        issues.update(self._validate_instrument_settings(name='Chopper', parameter='settings'))

        self._check_advanced_parameter()

        if not isinstance(float(self.getProperty("IncidentEnergy").value), numbers.Real):
            issues["IncidentEnergy"] = "Incident energy must be a real number"

        instrument_name = self.getProperty("Instrument").value
        allowed_frequencies = abins.parameters.instruments[instrument_name]['chopper_allowed_frequencies']
        default_frequency = abins.parameters.instruments[instrument_name].get('chopper_default_frequency', None)

        if (self.getProperty("ChopperFrequency").value) == '' and (default_frequency is None):
            issues["ChopperFrequency"] = "This instrument does not have a default chopper frequency"
        elif int(self.getProperty("ChopperFrequency").value) not in allowed_frequencies:
            issues["ChopperFrequency"] = (f"This chopper frequency is not valid for the instrument {instrument_name}. "
                                          "Valid frequencies: " + ", ".join([str(freq) for freq in allowed_frequencies]))

        return issues

    def PyExec(self):
        # 0) Create reporter to report progress
        # Before calculating S, we use 10% of the bar for two steps
        begin, end, steps = 0, 0.1, 2
        prog_reporter = Progress(self, begin, end, steps)

        # 1) get input parameters from a user
        self._get_properties()
        prog_reporter.report("Input data from the user has been collected.")

        # 2) read ab initio data
        ab_initio_data = abins.AbinsData.from_calculation_data(self._vibrational_or_phonon_data_file,
                                                               self._ab_initio_program)
        prog_reporter.report("Vibrational/phonon data has been read.")

        # 3) calculate S
        # Reset reporter to span range 10%-80%; s_calculator will decide how many steps are appropriate
        # so insert placeholder "1" for now.
        prog_reporter.resetNumSteps(1, 0.1, 0.8)

        s_calculator = abins.SCalculatorFactory.init(filename=self._vibrational_or_phonon_data_file,
                                                     temperature=self._temperature,
                                                     sample_form=self._sample_form,
                                                     abins_data=ab_initio_data,
                                                     autoconvolution=self._autoconvolution,
                                                     instrument=self._instrument,
                                                     quantum_order_num=self._num_quantum_order_events,
                                                     bin_width=self._bin_width)
        s_calculator.progress_reporter = prog_reporter
        s_data = s_calculator.get_formatted_data()
        raise Exception(s_data)

        # Hold reporter at 80% for this message
        prog_reporter.resetNumSteps(1, 0.8, 0.80000001)
        prog_reporter.report("Dynamical structure factors have been determined.")
        # Now determine number of remaining messages and set reporter for rest of run:
        n_messages = 3 + bool(self._sum_contributions) + bool(self._save_ascii)
        prog_reporter.resetNumSteps(n_messages, 0.8, 1)

        # 4) get atoms for which S should be plotted
        atoms_data = ab_initio_data.get_atoms_data()
        atom_numbers, atom_symbols = self.get_atom_selection(atoms_data=atoms_data, selection=self._atoms)
        prog_reporter.report("Atoms, for which dynamical structure factors should be plotted, have been determined.")

        workspaces = [self._create_dummy_workspace('dummy')]
        GroupWorkspaces(InputWorkspaces=workspaces, OutputWorkspace=self._out_ws_name)
        self.setProperty('OutputWorkspace', self._out_ws_name)

        # # 5) create workspaces for atoms in interest
        # workspaces = []
        # workspaces.extend(self.create_workspaces(atoms_symbols=atom_symbols, s_data=s_data, atoms_data=atoms_data,
        #                                          max_quantum_order=self._num_quantum_order_events))
        # workspaces.extend(self.create_workspaces(atom_numbers=atom_numbers, s_data=s_data, atoms_data=atoms_data,
        #                                          max_quantum_order=self._num_quantum_order_events))
        # prog_reporter.report("Workspaces with partial dynamical structure factors have been constructed.")

        # # 6) Create a workspace with sum of all atoms if required
        # if self._sum_contributions:
        #     self.create_total_workspace(workspaces)
        #     prog_reporter.report("Workspace with total S has been constructed.")

        # GroupWorkspaces(InputWorkspaces=workspaces, OutputWorkspace=self._out_ws_name)

        # # 8) save workspaces to ascii_file
        # if self._save_ascii:
        #     self.write_workspaces_to_ascii(ws_name=self._out_ws_name, scale=(1.0 / self._bin_width))
        #     prog_reporter.report("All workspaces have been saved to ASCII files.")

        # # 9) set  OutputWorkspace
        # self.setProperty('OutputWorkspace', self._out_ws_name)
        # prog_reporter.report("Group workspace with all required  dynamical structure factors has been constructed.")

    def _fill_s_workspace(self, s_points=None, workspace=None, protons_number=None, nucleons_number=None):
        """
        Puts S into workspace(s).

        :param s_points: dynamical factor for the given atom
        :param workspace:  workspace to be filled with S
        :param protons_number: number of protons in the given type fo atom
        :param nucleons_number: number of nucleons in the given type of atom
        """
        from abins.constants import FUNDAMENTALS, TWO_DIMENSIONAL_INSTRUMENTS

        if self._instrument.get_name() not in TWO_DIMENSIONAL_INSTRUMENTS:
            raise ValueError("Instrument {self._instrument_name} is not supported by this version of Abins")

        # only FUNDAMENTALS [data is 3d with length 1 in axis 0]
        if s_points.shape[0] == FUNDAMENTALS:
            self._fill_s_2d_workspace(s_points=s_points[0], workspace=workspace, protons_number=protons_number,
                                      nucleons_number=nucleons_number)

        # total workspaces [data is 2d array of S]
        elif s_points.shape[0] == abins.parameters.instruments[self._instrument.get_name()]['q_size']:
            self._fill_s_2d_workspace(s_points=s_points, workspace=workspace, protons_number=protons_number,
                                      nucleons_number=nucleons_number)

        # Multiple quantum order events [data is 3d table of S using axis 0 for quantum orders]
        else:
            dim = s_points.shape[0]
            partial_wrk_names = []

            for n in range(dim):
                seed = "quantum_event_%s" % (n + 1)
                wrk_name = workspace + "_" + seed
                partial_wrk_names.append(wrk_name)

                self._fill_s_2d_workspace(s_points=s_points[n], workspace=wrk_name, protons_number=protons_number,
                                          nucleons_number=nucleons_number)

                GroupWorkspaces(InputWorkspaces=partial_wrk_names, OutputWorkspace=workspace)

    def _create_dummy_workspace(self, name):
        wrk = WorkspaceFactory.create("Workspace2D", NVectors=1, XLength=2, YLength=1)
        wrk.setX(0, [0, 1])
        wrk.setY(0, [0])
        AnalysisDataService.addOrReplace(name, wrk)
        return wrk

    def _fill_s_2d_workspace(self, s_points=None, workspace=None, protons_number=None, nucleons_number=None):
        from abins.constants import Q_BEGIN, Q_END
        from mantid.api import NumericAxis

        if protons_number is not None:
            s_points = s_points *  self.get_cross_section(scattering=self._scale_by_cross_section,
                                                          protons_number=protons_number,
                                                          nucleons_number=nucleons_number)

        n_q_bins, n_freq_bins = s_points.shape

        wrk = WorkspaceFactory.create("Workspace2D", NVectors=n_freq_bins, XLength=n_q_bins + 1, YLength=n_q_bins)

        freq_axis = NumericAxis.create(n_freq_bins)

        q_size = abins.parameters.instruments[self._instrument.get_name()]['q_size']
        q_bins = np.linspace(start=Q_BEGIN, stop=Q_END, num=q_size + 1)

        freq_offset = (self._bins[1] - self._bins[0]) / 2
        for i, freq in enumerate(self._bins[1:]):
            wrk.setX(i, q_bins)
            wrk.setY(i, s_points[:, i].T)
            freq_axis.setValue(i, freq + freq_offset)
        freq_axis.setUnit("Energy_inWavenumber")
        wrk.replaceAxis(1, freq_axis)

        AnalysisDataService.addOrReplace(workspace, wrk)

        self.set_workspace_units(workspace, layout="2D")

    def _check_advanced_parameter(self):
        """
        Checks if parameters from abins.parameters are valid. If any parameter
        is invalid then RuntimeError is thrown with meaningful message.
        """

        message = " in scripts/abins/parameters.py. "

        self._check_common_advanced_parameters(message)
        self._check_2d_parameters(message)

    def _check_2d_parameters(self, message_end=None):
        # check 2D resolution
        resolution_2d = abins.parameters.instruments['TwoDMap']['resolution']
        if not (isinstance(resolution_2d, float) and resolution_2d > 0):
            raise RuntimeError("Invalid value of abins.parameters"
                               ".instruments['TwoDMap']['resolution']"
                               + message_end)

    def _get_properties(self):
        """
        Loads all properties to object's attributes.
        """
        self.get_common_properties()
        self._autoconvolution = self.getProperty("Autoconvolution").value

        self._instrument_kwargs = {"setting": self.getProperty("Chopper").value,
                                   "chopper_frequency": self.getProperty("ChopperFrequency")}
        self.set_instrument()

        instrument_params = abins.parameters.instruments[self._instrument.get_name()]

        # Incident energy currently handled differently, but this should be changed to an
        # instrument parameter as well
        instrument_params['e_init'] = float(self.getProperty("IncidentEnergy").value)

        # Sampling mesh is determined by
        # abins.parameters.sampling['min_wavenumber']
        # abins.parameters.sampling['max_wavenumber'],
        # while abins.parameters.sampling['bin_width'] is set from user input
        step = abins.parameters.sampling['bin_width'] = self._bin_width
        start = abins.parameters.sampling['min_wavenumber']
        stop = abins.parameters.sampling['max_wavenumber'] + step
        self._bins = np.arange(start=start, stop=stop, step=step, dtype=abins.constants.FLOAT_TYPE)


AlgorithmFactory.subscribe(Abins2D)
