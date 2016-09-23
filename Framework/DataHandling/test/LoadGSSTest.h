#ifndef LOADGSSTEST_H_
#define LOADGSSTEST_H_

#include "cxxtest/TestSuite.h"
#include "MantidDataHandling/LoadGSS.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidTestHelpers/ScopedFileHelper.h"
#include "MantidGeometry/Instrument.h"

using namespace Mantid;
using Mantid::DataHandling::LoadGSS;
using ScopedFileHelper::ScopedFile;
using Mantid::Kernel::V3D;

class LoadGSSTest : public CxxTest::TestSuite {
public:
  void test_TheBasics() {
    Mantid::DataHandling::LoadGSS loader;
    TS_ASSERT_THROWS_NOTHING(loader.initialize())
    TS_ASSERT_EQUALS(loader.name(), "LoadGSS")
    TS_ASSERT_EQUALS(loader.version(), 1)
  }

  void test_load_gss_txt() {
    API::IAlgorithm_sptr loader = createAlgorithm();
    loader->setPropertyValue("Filename", "gss.txt");
    TS_ASSERT(loader->execute())
    API::MatrixWorkspace_const_sptr ws = loader->getProperty("OutputWorkspace");
    // Check a few things in the workspace
    checkWorkspace(ws, 8, 816);
    auto x1 = ws->readX(0)[99];
    auto x2 = ws->readX(0)[100];
    auto y = ws->readY(0)[99];
    TS_ASSERT_DELTA((x1 + x2) / 2, 40844.0625, 1e-6);
    TS_ASSERT_DELTA(y, 145304004.625, 1e-6);
  }

  void test_large_x_values() {
    std::string gss =
        "LuBaCo4O7 HR BS P=2GPa Tcryo=130.0K                                   "
        "          \n"
        "# 1 Histograms\n"
        "# File generated by Mantid:\n"
        "# Instrument: WISH\n"
        "# From workspace named : w21552-2foc\n"
        "# with Y multiplied by the bin widths.\n"
        "# Primary flight path 40m \n"
        "# Total flight path 42.2222m, tth 58.308deg, DIFC 10398.8\n"
        "# Data for spectrum :0\n"
        "BANK 1 4399 4399 RALF   166409      124   166409 0.00074 FXYE\n"
        "   115202.20029   123456.00000002        0.00000000\n"
        "   115206.06310  1234567.00000003        0.00000000\n"
        "   115209.92877 12345678.00000004        0.00000000\n"
        "   115213.79731123456789.00000005        0.00000000\n"
        "   115217.66873234567890.00000006        0.00000000\n"
        "\n"
        "# Total flight path 42.060 m, tth 154.257 deg, DIFC 10398.8\n"
        "# Data for spectrum :1\n"
        "BANK 2 4399 4399 RALF   166409      124   166409 0.00074 FXYE\n"
        "   115202.20029        1.00000000        0.00000000\n"
        "   115206.06310        1.00000000        0.00000000\n"
        "   115209.92877        2.00000000        0.00000000\n"
        "   115213.79731        3.00000000        0.00000000\n"
        "   115217.66873        5.00000000        0.00000000\n";
    ScopedFile file(gss, "gss_large_x.txt");
    API::IAlgorithm_sptr loader = createAlgorithm();
    loader->setPropertyValue("Filename", file.getFileName());
	TS_ASSERT(loader->execute());
	TS_ASSERT_EQUALS(loader->isExecuted(), true);
    API::MatrixWorkspace_const_sptr ws = loader->getProperty("OutputWorkspace");
    auto x1 = ws->readX(0)[0];
    auto x2 = ws->readX(0)[1];
    auto dx = x2 - x1;
    auto y = ws->readY(0)[0] * dx;
    TS_ASSERT_DELTA((x1 + x2) / 2, 115202.20029, 1e-6);
    TS_ASSERT_DELTA(y, 123456.00000002, 1e-10);
    x1 = ws->readX(0)[3];
    x2 = ws->readX(0)[4];
    dx = x2 - x1;
    y = ws->readY(0)[3] * dx;
    TS_ASSERT_DELTA(y, 123456789.00000005, 1e-10);

    const auto source = ws->getInstrument()->getSource();
    const auto sample = ws->getInstrument()->getSample();
    TS_ASSERT_DELTA(sample->getDistance(*source), 40., 1e-4);

    const auto det0 = ws->getDetector(0);
    TS_ASSERT_DELTA(det0->getDistance(*sample), 2.2222, 1e-4);
    TS_ASSERT_DELTA(det0->getTwoTheta(V3D(0.0, 0.0, 0.0), V3D(0.0, 0.0, 1.0)) *
                        180. / M_PI,
                    58.308, 1e-4);
    const auto det1 = ws->getDetector(1);
    TS_ASSERT_DELTA(det1->getDistance(*sample), 2.060, 1e-4);
    TS_ASSERT_DELTA(det1->getTwoTheta(V3D(0.0, 0.0, 0.0), V3D(0.0, 0.0, 1.0)) *
                        180. / M_PI,
                    154.257, 1e-4);
  }

  void test_load_gss_ExtendedHeader_gsa() {
    API::IAlgorithm_sptr loader = createAlgorithm();
    loader->setPropertyValue("Filename", "gss-ExtendedHeader.gsa");
    TS_ASSERT(loader->execute())
    // Check a few things in the workspace
    checkWorkspace(loader->getProperty("OutputWorkspace"), 1, 6);
  }

  /** Test LoadGSS with setting spectrum No as bank ID
    */
  void test_load_gss_use_spec() {
    // Set property and execute
    LoadGSS loader;
    loader.initialize();

    loader.setPropertyValue("Filename", "gss1.txt");
    loader.setProperty("OutputWorkspace", "TestWS");
    loader.setProperty("UseBankIDasSpectrumNumber", true);

    TS_ASSERT(loader.execute());

    // Check result
    API::MatrixWorkspace_sptr outws =
        boost::dynamic_pointer_cast<API::MatrixWorkspace>(
            API::AnalysisDataService::Instance().retrieve("TestWS"));
    TS_ASSERT(outws);
    if (!outws)
      return;

    TS_ASSERT_EQUALS(outws->getNumberHistograms(), 3);
    if (outws->getNumberHistograms() != 3)
      return;

    TS_ASSERT_EQUALS(outws->getSpectrum(0).getSpectrumNo(), 1);
    TS_ASSERT_EQUALS(outws->getSpectrum(1).getSpectrumNo(), 3);
    TS_ASSERT_EQUALS(outws->getSpectrum(2).getSpectrumNo(), 5);

    API::AnalysisDataService::Instance().remove("TestWS");
  }

  void test_fails_gracefully_if_passed_wrong_filetype() {
    API::IAlgorithm_sptr loader = createAlgorithm();
    loader->setPropertyValue("Filename", "argus0026287.nxs");
    // Throws different exception type on different platforms!
    TS_ASSERT_THROWS_ANYTHING(loader->execute())

    API::IAlgorithm_sptr loader2 = createAlgorithm();
    loader2->setPropertyValue("Filename", "AsciiExample.txt");
    TS_ASSERT_THROWS(loader2->execute(), std::out_of_range)

    API::IAlgorithm_sptr loader3 = createAlgorithm();
    loader3->setPropertyValue("Filename", "CSP79590.raw");
    TS_ASSERT_THROWS(loader3->execute(), std::out_of_range)

    API::IAlgorithm_sptr loader4 = createAlgorithm();
    loader4->setPropertyValue("Filename", "VULCAN_2916_neutron0_event.dat");
    TS_ASSERT_THROWS(loader4->execute(), std::out_of_range)
  }

private:
  API::IAlgorithm_sptr createAlgorithm() {
    API::IAlgorithm_sptr alg =
        API::AlgorithmManager::Instance().createUnmanaged("LoadGSS");
    alg->initialize();
    alg->setChild(true);
    alg->setPropertyValue("OutputWorkspace", "fakeName");
    return alg;
  }

  void checkWorkspace(const API::MatrixWorkspace_const_sptr &ws, int nHist,
                      int nBins) {
    TS_ASSERT_EQUALS(ws->id(), "Workspace2D")
    TS_ASSERT_EQUALS(ws->getNumberHistograms(), nHist)
    TS_ASSERT_EQUALS(ws->size(), nBins)
    TS_ASSERT_EQUALS(ws->getAxis(0)->unit()->unitID(), "TOF")
  }
};

#endif // LOADGSSTEST_H_
