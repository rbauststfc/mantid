#ifndef MANTID_ALGORITHMS_ELASTICWINDOWTEST_H_
#define MANTID_ALGORITHMS_ELASTICWINDOWTEST_H_

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <iomanip>

#include "MantidKernel/System.h"

#include "MantidAlgorithms/ConvertUnits.h"
#include "MantidAlgorithms/CreateSampleWorkspace.h"
#include "MantidAlgorithms/ElasticWindow.h"
#include "MantidAlgorithms/Rebin.h"
#include "MantidAlgorithms/SetInstrumentParameter.h"

using namespace Mantid;
using namespace Mantid::Algorithms;
using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::Kernel::Units;

class ElasticWindowTest : public CxxTest::TestSuite
{
public:

  void setUp()
  {
    // Cretae a workspace and format it for the ElasticWindow alorithm

    CreateSampleWorkspace createAlg;
    createAlg.initialize();
    createAlg.setProperty("Function", "User Defined");
    createAlg.setProperty("UserDefinedFunction", "name=Lorentzian,Amplitude=100,PeakCentre=12700,FWHM=20;name=LinearBackground,A0=0.01");
    createAlg.setProperty("XMin", 27000.0);
    createAlg.setProperty("XMax", 28000.0);
    createAlg.setProperty("BinWidth", 10.0);
    createAlg.setProperty("NumBanks", 1);
    createAlg.setProperty("OutputWorkspace", "__ElasticWindowTest_sample");
    createAlg.execute();

    ConvertUnits convertUnitsAlg;
    convertUnitsAlg.initialize();
    convertUnitsAlg.setProperty("InputWorkspace", "__ElasticWindowTest_sample");
    convertUnitsAlg.setProperty("Target", "DeltaE");
    convertUnitsAlg.setProperty("EMode", "Indirect");
    convertUnitsAlg.setProperty("Efixed", 1.555);
    convertUnitsAlg.setProperty("OutputWorkspace", "__ElasticWindowTest_sample");
    convertUnitsAlg.execute();

    Rebin rebinAlg;
    rebinAlg.initialize();
    rebinAlg.setProperty("InputWorkspace", "__ElasticWindowTest_sample");
    rebinAlg.setProperty("Params", "-0.2,0.004,0.2");
    rebinAlg.setProperty("OutputWorkspace", "__ElasticWindowTest_sample");
    rebinAlg.execute();

    SetInstrumentParameter setParamAlg;
    setParamAlg.initialize();
    setParamAlg.setProperty("Workspace", "__ElasticWindowTest_sample");
    setParamAlg.setProperty("ParameterName", "Efixed");
    setParamAlg.setProperty("ParameterType", "Number");
    setParamAlg.setProperty("Value", "1.555");
    setParamAlg.execute();
  }

  /**
   * Test initialization of the algorithm is successful.
   */
  void test_init()
  {
    ElasticWindow alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize());
    TS_ASSERT(alg.isInitialized());
  }

  /**
   * Test running ElasticWindow with just the peak range defined.
   */
  void test_peakOnly()
  {
    ElasticWindow elwinAlg;
    elwinAlg.initialize();

    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("InputWorkspace", "__ElasticWindowTest_sample") );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("Range1Start", -0.1) );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("Range1End", 0.1) );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("OutputInQ", "__ElasticWindowTest_outputQ") );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("OutputInQSquared", "__ElasticWindowTest_outputQsq") );

    TS_ASSERT_THROWS_NOTHING( elwinAlg.execute() );
    TS_ASSERT( elwinAlg.isExecuted() );

    MatrixWorkspace_sptr qWs = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>("__ElasticWindowTest_outputQ");
    verifyQworkspace(qWs);
  }

  /**
   * Test running ElasticWindow with both the peak and bakground ranges defined.
   */
  void test_peakAndBackground()
  {
    ElasticWindow elwinAlg;
    elwinAlg.initialize();

    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("InputWorkspace", "__ElasticWindowTest_sample") );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("Range1Start", -0.04) );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("Range1End", 0.04) );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("Range2Start", 0.05) );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("Range2End", 0.06) );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("OutputInQ", "__ElasticWindowTest_outputQ") );
    TS_ASSERT_THROWS_NOTHING( elwinAlg.setProperty("OutputInQSquared", "__ElasticWindowTest_outputQsq") );

    TS_ASSERT_THROWS_NOTHING( elwinAlg.execute() );
    TS_ASSERT( elwinAlg.isExecuted() );

    MatrixWorkspace_sptr qWs = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>("__ElasticWindowTest_outputQ");
    verifyQworkspace(qWs);

    MatrixWorkspace_sptr q2Ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>("__ElasticWindowTest_outputQsq");
    verifyQ2workspace(q2Ws);
  }

private:

  void verifyQworkspace(MatrixWorkspace_sptr ws)
  {
    std::string unitID = ws->getAxis(0)->unit()->unitID();
    TS_ASSERT_EQUALS(unitID, "MomentumTransfer");
  }

  void verifyQ2workspace(MatrixWorkspace_sptr ws)
  {
    std::string unitID = ws->getAxis(0)->unit()->unitID();
    TS_ASSERT_EQUALS(unitID, "QSquared");
  }

};

#endif /* MANTID_ALGORITHMS_ELASTICWINDOWTEST_H_ */
