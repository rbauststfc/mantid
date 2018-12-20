// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_CUSTOMINTERFACES_RUNSPRESENTERTEST_H
#define MANTID_CUSTOMINTERFACES_RUNSPRESENTERTEST_H

#include "../../../ISISReflectometry/GUI/Runs/RunsPresenter.h"
#include "../../ReflMockObjects.h"
#include "../RunsTable/MockRunsTablePresenterFactory.h"
#include "../RunsTable/MockRunsTableView.h"
#include "MantidKernel/ConfigService.h"
#include "MantidKernel/make_unique.h"
#include "MantidQtWidgets/Common/Batch/MockJobTreeView.h"
#include "MantidQtWidgets/Common/MockProgressableView.h"
#include "MockRunsView.h"
#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace MantidQt::CustomInterfaces;
using namespace testing;

//=====================================================================================
// Functional tests
//=====================================================================================
class RunsPresenterTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static RunsPresenterTest *createSuite() { return new RunsPresenterTest(); }
  static void destroySuite(RunsPresenterTest *suite) { delete suite; }

  RunsPresenterTest()
      : m_thetaTolerance(0.01), m_instruments{"INTER", "SURF", "CRISP",
                                              "POLREF", "OFFSPEC"},
        m_view(), m_runsTableView(), m_progressView(), m_messageHandler(),
        m_autoreduction(new MockReflAutoreduction),
        m_searcher(new MockReflSearcher) {
    ON_CALL(m_view, table()).WillByDefault(Return(&m_runsTableView));
    ON_CALL(m_runsTableView, jobs()).WillByDefault(ReturnRef(m_jobs));
  }

  void testMakePresenter() { auto presenter = makePresenter(); }

  void testSettingsChanged() {
    // TODO
  }

  void testSearchWithEmptyString() {
    auto presenter = makePresenter();
    auto searchString = std::string("");
    EXPECT_CALL(m_view, getSearchString())
        .Times(1)
        .WillOnce(Return(searchString));
    expectSearchFailed();
    presenter.notify(IRunsPresenter::SearchFlag);
    verifyAndClear();
  }

  void testSearchCatalogLoginFails() {
    auto presenter = makePresenter();
    auto searchString = std::string("test string");
    EXPECT_CALL(m_view, getSearchString())
        .Times(1)
        .WillOnce(Return(searchString));
    // TODO: add expected call to python runner when implemented
    EXPECT_CALL(m_view, noActiveICatSessions()).Times(1);
    expectSearchFailed();
    presenter.notify(IRunsPresenter::SearchFlag);
    verifyAndClear();
  }

  void testSearchSucceeds() {
    // TODO: add this test when python runner implemented
  }

  void testStartNewAutoreduction() {
    auto presenter = makePresenter();
    EXPECT_CALL(m_view, getSearchString()).Times(2);
    EXPECT_CALL(*m_autoreduction, searchStringChanged(_))
        .WillOnce(Return(true));
    EXPECT_CALL(*m_autoreduction, setupNewAutoreduction(_))
        .WillOnce(Return(true));
    expectCheckForNewRuns();
    presenter.notify(IRunsPresenter::StartAutoreductionFlag);
    verifyAndClear();
  }

  void testStartNewAutoreductionWarnsUserIfTableChanged() {
    // TODO
  }

  void testStartRepeatAutoreduction() {
    auto presenter = makePresenter();
    EXPECT_CALL(m_view, getSearchString()).Times(2);
    EXPECT_CALL(*m_autoreduction, searchStringChanged(_))
        .WillOnce(Return(false));
    EXPECT_CALL(*m_autoreduction, setupNewAutoreduction(_))
        .WillOnce(Return(true));
    expectCheckForNewRuns();
    presenter.notify(IRunsPresenter::StartAutoreductionFlag);
    verifyAndClear();
  }

  void testPauseAutoreduction() {
    // TODO
    // auto presenter = makePresenter();
    // presenter.notify(IRunsPresenter::PauseAutoreductionFlag);
    // verifyAndClear();
  }

  void testAutoreductionPollsForNewRunsOnTimerEvent() {
    auto presenter = makePresenter();
    expectCheckForNewRuns();
    presenter.notify(IRunsPresenter::TimerEventFlag);
    verifyAndClear();
  }

  void testICATSearchComplete() {
    // TODO
    // auto presenter = makePresenter();
    // presenter.notify(IRunsPresenter::ICATSearchCompleteFlag);
    // verifyAndClear();
  }

  void testTransferWithNoRowsSelected() {
    auto presenter = makePresenter();
    auto const selectedRows = std::set<int>({});
    EXPECT_CALL(m_view, getSelectedSearchRows())
        .Times(1)
        .WillOnce(Return(selectedRows));
    EXPECT_CALL(m_view, missingRunsToTransfer()).Times(1);
    presenter.notify(IRunsPresenter::TransferFlag);
    verifyAndClear();
  }

  void testTransferWithAutoreductionRunning() {
    auto presenter = makePresenter();
    expectGetValidSearchRowSelection(presenter);
    EXPECT_CALL(*m_autoreduction, running()).Times(1).WillOnce(Return(true));
    expectCreateEndlessProgressIndicator();
    presenter.notify(IRunsPresenter::TransferFlag);
    verifyAndClear();
  }

  void testTransferWithAutoreductionStopped() {
    auto presenter = makePresenter();
    expectGetValidSearchRowSelection(presenter);
    EXPECT_CALL(*m_autoreduction, running()).Times(1).WillOnce(Return(false));
    expectCreatePercentageProgressIndicator();
    presenter.notify(IRunsPresenter::TransferFlag);
    verifyAndClear();
  }

  void testInstrumentChanged() {
    auto presenter = makePresenter();
    auto const instrument = std::string("TEST-instrumnet");
    EXPECT_CALL(m_view, getSearchInstrument())
        .Times(1)
        .WillOnce(Return(instrument));
    EXPECT_CALL(m_mainPresenter, setInstrumentName(instrument)).Times(1);
    presenter.notify(IRunsPresenter::InstrumentChangedFlag);
    verifyAndClear();
  }

  // TODO
  //  void testStartMonitor() {
  //    auto presenter = makePresenter();
  //    EXPECT_CALL(m_view, getMonitorAlgorithmRunner()).Times(1);
  //    EXPECT_CALL(m_view, getSearchInstrument()).Times(1);
  //    expectUpdateViewWhenMonitorStarting();
  //    presenter.notify(IRunsPresenter::StartMonitorFlag);
  //    verifyAndClear();
  //  }
  //
  //  void testStopMonitor() {
  //    auto presenter = makePresenter();
  //    expectUpdateViewWhenMonitorStopped();
  //    presenter.notify(IRunsPresenter::StopMonitorFlag);
  //    verifyAndClear();
  //  }
  //
  //  void testStartMonitorComplete() {
  //    auto presenter = makePresenter();
  //    expectUpdateViewWhenMonitorStarted();
  //    presenter.notify(IRunsPresenter::StartMonitorCompleteFlag);
  //    verifyAndClear();
  //  }

private:
  class RunsPresenterFriend : public RunsPresenter {
    friend class RunsPresenterTest;

  public:
    RunsPresenterFriend(IRunsView *mainView, ProgressableView *progressView,
                        RunsTablePresenterFactory makeRunsTablePresenter,
                        double thetaTolerance,
                        std::vector<std::string> const &instruments,
                        int defaultInstrumentIndex,
                        IReflMessageHandler *messageHandler,
                        boost::shared_ptr<IReflAutoreduction> autoreduction =
                            boost::shared_ptr<IReflAutoreduction>(),
                        boost::shared_ptr<IReflSearcher> searcher =
                            boost::shared_ptr<IReflSearcher>())
        : RunsPresenter(mainView, progressView, makeRunsTablePresenter,
                        thetaTolerance, instruments, defaultInstrumentIndex,
                        messageHandler, autoreduction, searcher) {}
  };

  RunsPresenterFriend makePresenter() {
    auto defaultInstrumentIndex = 0;
    EXPECT_CALL(m_view, subscribe(_)).Times(1);
    EXPECT_CALL(m_runsTableView, subscribe(_)).Times(1);
    EXPECT_CALL(m_view, table()).Times(1).WillOnce(Return(&m_runsTableView));
    EXPECT_CALL(m_runsTableView, jobs()).Times(1).WillOnce(ReturnRef(m_jobs));
    EXPECT_CALL(m_view,
                setInstrumentList(m_instruments, defaultInstrumentIndex))
        .Times(1);
    expectUpdateViewWhenMonitorStopped();

    auto presenter = RunsPresenterFriend(
        &m_view, &m_progressView,
        MockRunsTablePresenterFactory(m_instruments, m_thetaTolerance),
        m_thetaTolerance, m_instruments, defaultInstrumentIndex,
        &m_messageHandler, m_autoreduction, m_searcher);

    presenter.acceptMainPresenter(&m_mainPresenter);
    verifyAndClear();
    return presenter;
  }

  void verifyAndClear() {
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_runsTableView));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_progressView));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_messageHandler));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_autoreduction));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_searcher));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_jobs));
  }

  void expectUpdateViewWhenMonitorStarting() {
    EXPECT_CALL(m_view, setStartMonitorButtonEnabled(false));
    EXPECT_CALL(m_view, setStopMonitorButtonEnabled(false));
  }

  void expectUpdateViewWhenMonitorStarted() {
    EXPECT_CALL(m_view, setStartMonitorButtonEnabled(false));
    EXPECT_CALL(m_view, setStopMonitorButtonEnabled(true));
  }

  void expectUpdateViewWhenMonitorStopped() {
    EXPECT_CALL(m_view, setStartMonitorButtonEnabled(true));
    EXPECT_CALL(m_view, setStopMonitorButtonEnabled(false));
  }

  void expectStopAutoreduction() {
    EXPECT_CALL(m_view, stopTimer()).Times(1);
    EXPECT_CALL(*m_autoreduction, stop()).Times(1);
  }

  void expectSearchFailed() {
    EXPECT_CALL(m_view, getAlgorithmRunner()).Times(0);
    expectStopAutoreduction();
  }

  void expectCheckForNewRuns() {
    EXPECT_CALL(m_view, stopTimer()).Times(1);
    EXPECT_CALL(m_view, startIcatSearch()).Times(1);
  }

  void expectGetValidSearchRowSelection(RunsPresenterFriend &presenter) {
    // Select a couple of rows with random indices
    auto row1Index = 3;
    auto row2Index = 5;
    auto const selectedRows = std::set<int>({row1Index, row2Index});
    EXPECT_CALL(m_view, getSelectedSearchRows())
        .Times(1)
        .WillOnce(Return(selectedRows));
    // Set up a mock search model in the presenter to return something
    // sensible for getRowData
    auto searchModel = boost::make_shared<MockReflSearchModel>(
        "13460", "my title th=0.5", "my location");
    presenter.m_searchModel = searchModel;
  }

  void expectCreateEndlessProgressIndicator() {
    EXPECT_CALL(m_progressView, clearProgress()).Times(1);
    EXPECT_CALL(m_progressView, setProgressRange(_, _)).Times(2);
  }

  void expectCreatePercentageProgressIndicator() {
    EXPECT_CALL(m_progressView, clearProgress()).Times(1);
    EXPECT_CALL(m_progressView, setProgressRange(_, _)).Times(2);
  }

  double m_thetaTolerance;
  std::vector<std::string> m_instruments;
  MockRunsView m_view;
  MockRunsTableView m_runsTableView;
  MockReflBatchPresenter m_mainPresenter;
  MockProgressableView m_progressView;
  MockMessageHandler m_messageHandler;
  boost::shared_ptr<MockReflAutoreduction> m_autoreduction;
  boost::shared_ptr<MockReflSearcher> m_searcher;
  NiceMock<MantidQt::MantidWidgets::Batch::MockJobTreeView> m_jobs;
};

#endif /* MANTID_CUSTOMINTERFACES_RUNSPRESENTERTEST_H */
