// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2024 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllConfig.h"
#include "IALCDataLoadingView.h"
#include "MantidAPI/IAlgorithm_fwd.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidKernel/System.h"
#include <atomic>

namespace MantidQt::CustomInterfaces {

class MANTIDQT_MUONINTERFACE_DLL ALCDataLoadingModel {

public:
  ALCDataLoadingModel();
  ~ALCDataLoadingModel() = default;
  void load(const IALCDataLoadingView *view);
  void cancelLoading() const;
  Mantid::API::MatrixWorkspace_sptr exportWorkspace();
  bool checkCustomGrouping(const std::string &detGroupingType, const std::string &forwardGrouping,
                           const std::string &backwardGrouping);
  void updateAutoLoadCancelled();
  bool loadFilesFromWatchingDirectory(const std::string &firstFile, const std::vector<std::string> &files,
                                      const std::string &runsText);
  static std::string getPathFromFiles(std::vector<std::string> files);

  // Getters
  bool getLoadingData();
  Mantid::API::MatrixWorkspace_sptr getLoadedData();
  std::vector<std::string> &getLogs();
  std::vector<std::string> &getPeriods();
  Mantid::API::MatrixWorkspace_sptr getWsForMuonInfo();
  double getMinTime() const;
  std::string &getRunsText();

  // Setters
  void setLoadingData(bool isLoading);
  void setLoadedData(const Mantid::API::MatrixWorkspace_sptr &data);
  void setLogs(const Mantid::API::MatrixWorkspace_sptr &ws);
  void setPeriods(const Mantid::API::Workspace_sptr &ws);
  void setWsForMuonInfo(const std::string &filename);
  void setDirectoryChanged(bool hasDirectoryChanged);
  void setFilesToLoad(const std::vector<std::string> &files);

private:
  static std::string isCustomGroupingValid(const std::string &group, bool &isValid);
  static int extractRunNumber(const std::string &file);

  /// Last loaded data workspace
  Mantid::API::MatrixWorkspace_sptr m_loadedData;

  /// Number of detectors for current first run
  size_t m_numDetectors;

  // bool to state whether loading data
  std::atomic_bool m_loadingData;

  // Loading algorithm
  Mantid::API::IAlgorithm_sptr m_LoadingAlg;

  /// Flag for changes in watched directory
  std::atomic_bool m_directoryChanged;

  /// Last run loaded by auto
  int m_lastRunLoadedAuto;

  /// Files that are loaded
  std::vector<std::string> m_filesToLoad;

  /// Last run added by auto was addes as range
  std::atomic_bool m_wasLastAutoRange;

  /// Variables used to update available muon info
  Mantid::API::MatrixWorkspace_sptr m_wsForInfo;
  std::vector<std::string> m_periods;
  std::vector<std::string> m_logs;
  double m_minTime;
  std::string m_runsText;
};
} // namespace MantidQt::CustomInterfaces
