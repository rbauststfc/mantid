// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>

#include "Common/DataValidationHelper.h"
#include "MantidQtWidgets/Common/AddWorkspaceDialog.h"
#include "Processor/MomentsModel.h"
#include "Processor/MomentsPresenter.h"
#include "Processor/MomentsView.h"

#include "../Common/MockObjects.h"
#include "../QENSFitting/MockObjects.h"

#include "MantidFrameworkTestHelpers/IndirectFitDataCreationHelper.h"
#include "MantidKernel/WarningSuppressions.h"

using namespace Mantid::API;
using namespace Mantid::IndirectFitDataCreationHelper;
using namespace MantidQt::CustomInterfaces;
using namespace MantidQt::CustomInterfaces::Inelastic;
using namespace testing;

class MomentsPresenterTest : public CxxTest::TestSuite {
public:
  static MomentsPresenterTest *createSuite() { return new MomentsPresenterTest(); }

  static void destroySuite(MomentsPresenterTest *suite) { delete suite; }

  void setUp() override {
    m_view = std::make_unique<NiceMock<MockMomentsView>>();
    m_outputPlotView = std::make_unique<NiceMock<MockOutputPlotOptionsView>>();

    auto model = std::make_unique<NiceMock<MockMomentsModel>>();
    m_model = model.get();

    ON_CALL(*m_view, getPlotOptions()).WillByDefault(Return((m_outputPlotView.get())));
    m_presenter = std::make_unique<MomentsPresenter>(nullptr, m_view.get(), std::move(model));

    m_workspace = createWorkspace(5);
    m_ads = std::make_unique<SetUpADSWithWorkspace>("workspace_test", m_workspace);
  }

  void tearDown() override {
    AnalysisDataService::Instance().clear();

    TS_ASSERT(Mock::VerifyAndClearExpectations(m_view.get()));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_model));
    TS_ASSERT(Mock::VerifyAndClearExpectations(m_outputPlotView.get()));

    m_presenter.reset();
    m_view.reset();
    m_outputPlotView.reset();
  }

  ///----------------------------------------------------------------------
  /// Unit Tests that test the signals, methods and slots of the presenter
  ///----------------------------------------------------------------------

  void test_handleScaleChanged_sets_correct_bool_property() {

    EXPECT_CALL(*m_model, setScale(true)).Times((Exactly(1)));
    m_presenter->handleScaleChanged(true);

    EXPECT_CALL(*m_model, setScaleValue(true)).Times((Exactly(1)));
    m_presenter->handleScaleValueChanged(true);
  }

  void test_handleValueChanged_sets_correct_double_property() {
    double value = 0.1;

    EXPECT_CALL(*m_model, setEMin(value)).Times((Exactly(1)));
    m_presenter->handleValueChanged("EMin", value);
    EXPECT_CALL(*m_model, setEMax(value)).Times((Exactly(1)));
    m_presenter->handleValueChanged("EMax", value);
  }

  void test_handleDataReady_does_not_work_with_invalid_data() {
    auto dataSelector = std::make_unique<DataSelector>();
    ON_CALL(*m_view, getDataSelector()).WillByDefault(Return(dataSelector.get()));
    TS_ASSERT(!m_presenter->validate());
    dataSelector.reset();
  }

private:
  NiceMock<MockMomentsModel> *m_model;
  std::unique_ptr<NiceMock<MockOutputPlotOptionsView>> m_outputPlotView;
  std::unique_ptr<NiceMock<MockMomentsView>> m_view;
  std::unique_ptr<MomentsPresenter> m_presenter;

  MatrixWorkspace_sptr m_workspace;
  std::unique_ptr<SetUpADSWithWorkspace> m_ads;
};
