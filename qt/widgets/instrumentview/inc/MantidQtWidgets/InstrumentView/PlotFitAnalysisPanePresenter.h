// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2014 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllOption.h"
#include "MantidQtWidgets/Common/ObserverPattern.h"
#include "MantidQtWidgets/InstrumentView/PlotFitAnalysisPaneModel.h"
#include "MantidQtWidgets/InstrumentView/PlotFitAnalysisPaneView.h"
#include <string>

namespace MantidQt {
namespace MantidWidgets {

class EXPORT_OPT_MANTIDQT_INSTRUMENTVIEW IPlotFitAnalysisPanePresenter {

public:
  IPlotFitAnalysisPanePresenter(){};
  virtual ~IPlotFitAnalysisPanePresenter(){};
  virtual void destructor() = 0;
  virtual IPlotFitAnalysisPaneView *getView() = 0;
  virtual std::string getCurrentWS() = 0;
  virtual void clearCurrentWS() = 0;
  virtual void peakCentreEditingFinished() = 0;
  virtual void fitClicked() = 0;
  virtual void updateEstimateClicked() = 0;
  virtual void addSpectrum(const std::string &wsName) = 0;
};

class EXPORT_OPT_MANTIDQT_INSTRUMENTVIEW PlotFitAnalysisPanePresenter : public QObject,
                                                                        public IPlotFitAnalysisPanePresenter {
  Q_OBJECT

public:
  explicit PlotFitAnalysisPanePresenter(IPlotFitAnalysisPaneView *m_view, IPlotFitAnalysisPaneModel *m_model);
  ~PlotFitAnalysisPanePresenter() {
    delete m_model;
    delete m_fitObserver;
  };
  void destructor() override { this->~PlotFitAnalysisPanePresenter(); };
  IPlotFitAnalysisPaneView *getView() override { return m_view; };
  std::string getCurrentWS() override { return m_currentName; };
  void clearCurrentWS() override { m_currentName = ""; };
  void peakCentreEditingFinished() override;
  void fitClicked() override;
  void updateEstimateClicked() override;
  void addSpectrum(const std::string &wsName) override;

private:
  void updatePeakCentreInViewFromModel();

  VoidObserver *m_peakCentreObserver;
  VoidObserver *m_fitObserver;
  VoidObserver *m_updateEstimateObserver;
  IPlotFitAnalysisPaneView *m_view;
  IPlotFitAnalysisPaneModel *m_model;
  std::string m_currentName;
};
} // namespace MantidWidgets
} // namespace MantidQt
