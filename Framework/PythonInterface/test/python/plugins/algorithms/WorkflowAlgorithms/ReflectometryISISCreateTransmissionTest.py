# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2024 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
import numpy as np

from mantid import config
from mantid.api import AnalysisDataService, MatrixWorkspace, WorkspaceGroup
from mantid.simpleapi import CreateWorkspace
from testhelpers import create_algorithm

from plugins.algorithms.WorkflowAlgorithms.ReflectometryISISCreateTransmission import Prop


class ReflectometryISISCreateTransmissionTest(unittest.TestCase):
    _CONFIG_KEY_FACILITY = "default.facility"
    _CONFIG_KEY_INST = "default.instrument"

    _OUTPUT_WS_NAME = "out_ws"

    _LOAD_ALG = "LoadNexus"
    _FLOOD_ALG = "ApplyFloodWorkspace"
    _BACK_SUB_ALG = "ReflectometryBackgroundSubtraction"
    _TRANS_WS_ALG = "CreateTransmissionWorkspaceAuto"

    def setUp(self):
        self._oldFacility = config[self._CONFIG_KEY_FACILITY]
        if self._oldFacility.strip() == "":
            self._oldFacility = "TEST_LIVE"
        self._oldInstrument = config[self._CONFIG_KEY_INST]
        config[self._CONFIG_KEY_FACILITY] = "ISIS"
        config[self._CONFIG_KEY_INST] = "POLREF"

    def tearDown(self):
        AnalysisDataService.clear()
        config[self._CONFIG_KEY_FACILITY] = self._oldFacility
        config[self._CONFIG_KEY_INST] = self._oldInstrument

    def test_correct_output(self):
        output_ws = self._run_algorithm(self._create_args("INTER13460", back_sub_roi=""))
        self.assertIsNotNone(output_ws)
        self.assertIsInstance(output_ws, MatrixWorkspace)
        expected_history = [self._LOAD_ALG, self._FLOOD_ALG, self._TRANS_WS_ALG]
        self._check_history(output_ws, expected_history)

    def test_correct_output_for_workspace_groups(self):
        expected_history = [
            self._LOAD_ALG,
            self._FLOOD_ALG,
            self._FLOOD_ALG,
            self._BACK_SUB_ALG,
            self._BACK_SUB_ALG,
            self._TRANS_WS_ALG,
            self._TRANS_WS_ALG,
        ]
        self._run_workspace_group_test(expected_history=expected_history)

    def test_output_workspace_group_returned_when_run_as_child(self):
        self._run_workspace_group_test()

    def test_flood_correction_skipped_if_not_requested(self):
        expected_history = [self._LOAD_ALG, self._BACK_SUB_ALG, self._BACK_SUB_ALG, self._TRANS_WS_ALG, self._TRANS_WS_ALG]
        self._run_workspace_group_test(perform_flood=False, expected_history=expected_history)

    def test_background_subtraction_skipped_if_not_requested(self):
        expected_history = [self._LOAD_ALG, self._FLOOD_ALG, self._FLOOD_ALG, self._TRANS_WS_ALG, self._TRANS_WS_ALG]
        self._run_workspace_group_test(back_sub_roi="", expected_history=expected_history)

    def _create_args(self, input_run, perform_flood=True, back_sub_roi="100-200"):
        args = {
            Prop.INPUT_RUN: input_run,
            Prop.TRANS_ROI: "4",
            Prop.I0_MON_IDX: "2",
            Prop.MON_WAV_MIN: "2.5",
            Prop.MON_WAV_MAX: "10.0",
            Prop.OUTPUT_WS: self._OUTPUT_WS_NAME,
        }

        if perform_flood:
            args[Prop.FLOOD_WS] = self._create_flood_ws()

        if back_sub_roi:
            args[Prop.BACK_SUB_ROI] = back_sub_roi

        return args

    @staticmethod
    def _create_flood_ws():
        """Creates a MatrixWorkspace with a single bin of data. The workspace has 256 spectra with values from
        0.0 to 2.56 in steps of ~0.01"""
        flood_ws = CreateWorkspace(DataX=[0.0, 1.0], DataY=np.linspace(0.0, 2.56, 256), NSpec=256, UnitX="TOF")
        return flood_ws

    def _setup_algorithm(self, args, set_child):
        alg = create_algorithm("ReflectometryISISCreateTransmission", **args)
        alg.setChild(set_child)
        alg.setRethrows(True)
        return alg

    def _run_algorithm(self, args, set_child=False):
        alg = self._setup_algorithm(args, set_child)
        alg.execute()
        if set_child:
            return alg.getProperty(Prop.OUTPUT_WS).value
        else:
            return AnalysisDataService.retrieve(self._OUTPUT_WS_NAME)

    def _check_history(self, ws, expected):
        """Return true if expected algorithm names are found in the workspace history"""
        history = ws.getHistory()
        self.assertFalse(history.empty())
        child_alg_histories = history.getAlgorithmHistory(history.size() - 1).getChildHistories()
        self.assertEqual([alg.name() for alg in child_alg_histories], expected)

    def _run_workspace_group_test(self, perform_flood=True, back_sub_roi="100-200", expected_history=None):
        run_as_child = expected_history is None
        # Omitting the instrument prefix from the run number should use the default instrument
        output_ws = self._run_algorithm(self._create_args("14966", perform_flood, back_sub_roi), run_as_child)

        self.assertIsNotNone(output_ws)
        self.assertIsInstance(output_ws, WorkspaceGroup)

        if not run_as_child:
            self._check_history(output_ws[0], expected_history)


if __name__ == "__main__":
    unittest.main()
