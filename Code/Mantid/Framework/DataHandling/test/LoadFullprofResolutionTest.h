#ifndef MANTID_DATAHANDLING_LOADFULLPROFRESOLUTIONTEST_H_
#define MANTID_DATAHANDLING_LOADFULLPROFRESOLUTIONTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidDataHandling/LoadFullprofResolution.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidAPI/TableRow.h"
#include "MantidDataHandling/LoadInstrument.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/Component.h"
#include "MantidGeometry/Instrument/FitParameter.h"
#include <fstream>
#include <Poco/File.h>

using Mantid::DataHandling::LoadFullprofResolution;

using namespace Mantid;
using namespace Mantid::DataHandling;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;
using namespace Mantid::API;

using namespace std;

class LoadFullprofResolutionTest : public CxxTest::TestSuite
{
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static LoadFullprofResolutionTest *createSuite() { return new LoadFullprofResolutionTest(); }
  static void destroySuite( LoadFullprofResolutionTest *suite ) { delete suite; }

  //----------------------------------------------------------------------------------------------
  /** Test import from a 1-bank irf file
    */
  void test_1BankCase()
  {
    // 1. Generate file
    string filename("Test1Bank.irf");
    generate1BankIrfFile(filename);

    // 2. Load
    LoadFullprofResolution alg;
    alg.initialize();

    alg.setProperty("Filename", filename);
    alg.setPropertyValue("Banks", "1");
    alg.setProperty("OutputTableWorkspace", "TestBank1Table");

    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    TableWorkspace_sptr outws = boost::dynamic_pointer_cast<TableWorkspace>(
    AnalysisDataService::Instance().retrieve("TestBank1Table"));
    TS_ASSERT(outws);

    TS_ASSERT_EQUALS(outws->columnCount(), 2);
    TS_ASSERT_EQUALS(outws->rowCount(), getExpectedNumberOfRows() );

    // 3. Verify name and value
    map<string, double> parammap;
    parseTableWorkspace(outws, parammap);

    TS_ASSERT_EQUALS(parammap.count("Zero"),1);
    TS_ASSERT_EQUALS(parammap.count("Sig2"),1);
    TS_ASSERT_EQUALS(parammap.count("Beta0t"),1);

    TS_ASSERT_DELTA(parammap["Zero"], -1.00, 0.0001);
    TS_ASSERT_DELTA(parammap["Sig2"], sqrt(514.546), 0.0001);
    TS_ASSERT_DELTA(parammap["Beta0t"], 85.918922, 0.00001);

    // 4. Clean
    AnalysisDataService::Instance().remove("TestBank1Table");
    Poco::File("Test1Bank.irf").remove();

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test import from a 1-bank irf file
    */
  void test_2BankCase()
  {
    // 1. Generate file
    string filename("Test2Bank.irf");
    generate2BankIrfFile(filename);

    // 2. Load
    LoadFullprofResolution alg;
    alg.initialize();

    alg.setProperty("Filename", filename);
    alg.setPropertyValue("Banks", "3");
    alg.setProperty("OutputTableWorkspace", "TestBank3Table");

    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    TableWorkspace_sptr outws = boost::dynamic_pointer_cast<TableWorkspace>(
          AnalysisDataService::Instance().retrieve("TestBank3Table"));
    TS_ASSERT(outws);

    TS_ASSERT_EQUALS(outws->columnCount(), 2);
    TS_ASSERT_EQUALS(outws->rowCount(), getExpectedNumberOfRows());

    // 3. Verify name and value
    map<string, double> parammap;
    parseTableWorkspace(outws, parammap);
    TS_ASSERT_EQUALS(parammap.count("Dtt1"),1);
    TS_ASSERT_EQUALS(parammap.count("Sig1"),1);
    TS_ASSERT_EQUALS(parammap.count("Alph0t"),1);

    TS_ASSERT_DELTA(parammap["Dtt1"], 22586.10156, 0.0001);
    TS_ASSERT_DELTA(parammap["Sig1"], sqrt(10.00), 0.0001);
    TS_ASSERT_DELTA(parammap["Alph0t"], 86.059, 0.00001);

    // 4. Clean
    AnalysisDataService::Instance().remove("TestBank3Table");
    Poco::File("Test2Bank.irf").remove();

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test import all banks from a 2-bank irf file
    */
  void test_LoadAllBankCase()
  {
    // Generate file
    string filename("Test2Bank.irf");
    generate2BankIrfFile(filename);

    // Init LoadFullprofResolution
    LoadFullprofResolution alg;
    alg.initialize();

    // Set up
    alg.setProperty("Filename", filename);
    alg.setProperty("OutputTableWorkspace", "TestBank4Table");

    // Execute
    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    // Check output workspace
    TableWorkspace_sptr outws = boost::dynamic_pointer_cast<TableWorkspace>(
          AnalysisDataService::Instance().retrieve("TestBank4Table"));
    TS_ASSERT(outws);

    // Check table workspace size
    TS_ASSERT_EQUALS(outws->columnCount(), 3);
    TS_ASSERT_EQUALS(outws->rowCount(), getExpectedNumberOfRows());

    // Verify value (including bank number)
    map<string, double> parammap1;
    parseTableWorkspace(outws, parammap1);
    TS_ASSERT_DELTA(parammap1["BANK"], 1.0, 0.0001);
    TS_ASSERT_DELTA(parammap1["Dtt1"], 22580.59157, 0.0001);
    TS_ASSERT_DELTA(parammap1["Sig1"], sqrt(0.00044), 0.0001);
    TS_ASSERT_DELTA(parammap1["Alph0t"], 0.010156, 0.00001);

    map<string, double> parammap2;
    parseTableWorkspace2(outws, parammap2);
    TS_ASSERT_DELTA(parammap2["BANK"], 3.0, 0.0001);
    TS_ASSERT_DELTA(parammap2["Dtt1"], 22586.10156, 0.0001);
    TS_ASSERT_DELTA(parammap2["Sig1"], sqrt(10.00), 0.0001);
    TS_ASSERT_DELTA(parammap2["Alph0t"], 86.059, 0.00001);



    // Clean
    AnalysisDataService::Instance().remove("TestBank4Table");
    Poco::File("Test2Bank.irf").remove();

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test import all banks from a 2-bank irf file
    */
  void test_Load3BankCase()
  {
    // Generate file
    string filename("Test3Bank.irf");
    generate3BankIrfFile(filename);

    // Init LoadFullprofResolution
    LoadFullprofResolution alg;
    alg.initialize();

    // Set up
    alg.setProperty("Filename", filename);
    alg.setProperty("OutputTableWorkspace", "TestBank5Table");
    alg.setPropertyValue("Banks", "2-4");

    // Execute
    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    // Check output workspace
    TableWorkspace_sptr outws = boost::dynamic_pointer_cast<TableWorkspace>(
          AnalysisDataService::Instance().retrieve("TestBank5Table"));
    TS_ASSERT(outws);

    // Check table workspace size
    TS_ASSERT_EQUALS(outws->columnCount(), 4);
    TS_ASSERT_EQUALS(outws->rowCount(), getExpectedNumberOfRows());

    // Verify value
    map<string, double> parammap1;
    parseTableWorkspace(outws, parammap1);
    TS_ASSERT_DELTA(parammap1["Dtt1"], 22580.59157, 0.0001);
    TS_ASSERT_DELTA(parammap1["Sig1"], sqrt(0.00044), 0.0001);
    TS_ASSERT_DELTA(parammap1["Alph0t"], 0.010156, 0.00001);

    // Clean
    AnalysisDataService::Instance().remove("TestBank5Table");
    Poco::File("Test3Bank.irf").remove();

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test import of ALFBE, GAMMA and SIGMA parameters
  *   and check they are given their expected names.
  *   ConvertFullprofToXML relies on features tested here.
  */
  void test_ags_parameters()
  {
    // 1. Generate file
    string filename("TestAGS.irf");
    generate1BankIrfFile(filename);

    // 2. Load
    LoadFullprofResolution alg;
    alg.initialize();

    alg.setProperty("Filename", filename);
    alg.setPropertyValue("Banks", "1");
    alg.setProperty("OutputTableWorkspace", "TestAGSTable");

    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    TableWorkspace_sptr outws = boost::dynamic_pointer_cast<TableWorkspace>(
    AnalysisDataService::Instance().retrieve("TestAGSTable"));
    TS_ASSERT(outws);

    // 3. Verify names 
    map<string, double> parammap;
    parseTableWorkspace(outws, parammap);

    // 3a. ALFBE
    TS_ASSERT_EQUALS(parammap.count("Alph0"),1);
    TS_ASSERT_EQUALS(parammap.count("Beta0"),1);
    TS_ASSERT_EQUALS(parammap.count("Alph1"),1);
    TS_ASSERT_EQUALS(parammap.count("Beta1"),1);
    // 3b. GAMMA
    TS_ASSERT_EQUALS(parammap.count("Gam2"),1);
    TS_ASSERT_EQUALS(parammap.count("Gam1"),1);
    TS_ASSERT_EQUALS(parammap.count("Gam0"),1);
    // 3c. SIGMA
    TS_ASSERT_EQUALS(parammap.count("Sig2"),1);
    TS_ASSERT_EQUALS(parammap.count("Sig1"),1);
    TS_ASSERT_EQUALS(parammap.count("Sig0"),1);

    // 4. Clean
    AnalysisDataService::Instance().remove("TestAGSTable");
    Poco::File("TestAGS.irf").remove();

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test that when the workspace property is used
  **  that parameters are correctly loaded into this workspace
  *   The GEM instrument is used
  */
  void test_workspace()
  {
    // Generate file
    string filename("TestWorskpace.irf");
    generate1BankIrfFile(filename);

    // Load workspace wsName
    load_GEM();

    // Set up algorithm to load into the workspace
    LoadFullprofResolution alg;
    alg.initialize();

    alg.setProperty("Filename", filename);
    alg.setPropertyValue("Banks", "1");
    alg.setProperty("Workspace", wsName);

    // Execute
    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    // Check parameters in workspace T
    // The workspace is a workspace group with one member corresponding to the one bank in the IRF file
    WorkspaceGroup_sptr gws;
    gws = AnalysisDataService::Instance().retrieveWS<WorkspaceGroup>(wsName);
    Workspace_sptr wsi = gws->getItem(0) ;
    auto ws = boost::dynamic_pointer_cast<MatrixWorkspace>(wsi);
    Mantid::Geometry::ParameterMap& paramMap = ws->instrumentParameters();
    boost::shared_ptr<const Mantid::Geometry::Instrument> instr = ws->getInstrument();


    Mantid::Geometry::Parameter_sptr alpha0Param = paramMap.get(&(*instr), "Alpha0", "fitting");
    TS_ASSERT(alpha0Param);
    if(alpha0Param) 
    {
      const Mantid::Geometry::FitParameter& fitParam = alpha0Param->value<Mantid::Geometry::FitParameter>();
      TS_ASSERT_DELTA( boost::lexical_cast<double>(fitParam.getFormula()), 0.000008, 0.0000001);
    }

    Mantid::Geometry::Parameter_sptr beta0Param = paramMap.get(&(*instr), "Beta0", "fitting");
    TS_ASSERT(beta0Param);
    if(beta0Param) 
    {
      const Mantid::Geometry::FitParameter& fitParam = beta0Param->value<Mantid::Geometry::FitParameter>();
      TS_ASSERT_DELTA( boost::lexical_cast<double>(fitParam.getFormula()), 6.251096, 0.0000001);
    }

    Mantid::Geometry::Parameter_sptr alpha1Param = paramMap.get(&(*instr), "Alpha1", "fitting");
    TS_ASSERT(alpha1Param);
    if(alpha1Param) 
    {
      const Mantid::Geometry::FitParameter& fitParam = alpha1Param->value<Mantid::Geometry::FitParameter>();
      TS_ASSERT_DELTA( boost::lexical_cast<double>(fitParam.getFormula()), 0.0, 0.0000001);
    }

    Mantid::Geometry::Parameter_sptr beta1Param = paramMap.get(&(*instr), "Kappa", "fitting");
    TS_ASSERT(beta1Param);
    if(beta0Param) 
    {
      const Mantid::Geometry::FitParameter& fitParam = beta1Param->value<Mantid::Geometry::FitParameter>();
      TS_ASSERT_DELTA( boost::lexical_cast<double>(fitParam.getFormula()), 0.0, 0.0000001);
    }


    boost::shared_ptr<const Mantid::Geometry::IComponent> bank = instr->getComponentByName("bank1");
    Mantid::Geometry::Parameter_sptr sigmaSqParam = paramMap.get(&(*bank), "SigmaSquared", "fitting");
    TS_ASSERT(sigmaSqParam);
    if(sigmaSqParam) 
    {
      const Mantid::Geometry::FitParameter& fitParam = sigmaSqParam->value<Mantid::Geometry::FitParameter>();
      double formulaValueCantreAt0 = fitParam.getValue( 0.0 );    // Value for centre=0.0
      TS_ASSERT_DELTA( formulaValueCantreAt0, 0.355, 0.0000001);
      double formulaValueCantreAt10 = fitParam.getValue( 10.0 );  // Value for centre=10.0
      TS_ASSERT_DELTA( formulaValueCantreAt10, 0.399, 0.0000001);
    }

    Mantid::Geometry::Parameter_sptr gammaParam = paramMap.get(&(*bank), "Gamma", "fitting");
    TS_ASSERT(gammaParam);
    if(gammaParam) 
    {
      const Mantid::Geometry::FitParameter& fitParam = gammaParam->value<Mantid::Geometry::FitParameter>();
      double formulaValueCantreAt0 = fitParam.getValue( 0.0 );    // Value for centre=0.0
      TS_ASSERT_DELTA( formulaValueCantreAt0, 0.0, 0.0000001);
      double formulaValueCantreAt10 = fitParam.getValue( 10.0 );  // Value for centre=10.0
      TS_ASSERT_DELTA( formulaValueCantreAt10, 0.0, 0.0000001);
    }

    // Clean
    Poco::File("TestWorskpace.irf").remove();
  }

  //----------------------------------------------------------------------------------------------
  /** Test that algorithm does not run, 
  *   if neither the OutputTableWorkspace nor Workspace
  **  property is set.
  */
  void test_no_output()
  {
    // Generate file
    string filename("TestNoOutput.irf");
    generate1BankIrfFile(filename);

    // Set up algorithm without specifying OutputTableWorkspace or Workspace
    LoadFullprofResolution alg;
    alg.initialize();

    alg.setProperty("Filename", filename);
    alg.setPropertyValue("Banks", "1");

    // Execute and check that execution failed
    alg.execute();
    TS_ASSERT(!alg.isExecuted());

    // Clean
    Poco::File("TestNoOutput.irf").remove();
  }

  //----------------------------------------------------------------------------------------------
  /** Test that the number from the NPROF is read correctly 
  **  and has correct name in table.
  */
  void test_nprof()
  {
    // Generate file
    string filename("TestNPROF.irf");
    generate3BankIrfFile(filename);

    // Init LoadFullprofResolution
    LoadFullprofResolution alg;
    alg.initialize();

    // Set up
    alg.setProperty("Filename", filename);
    alg.setProperty("OutputTableWorkspace", "TestNPROFTable");

    // Execute
    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());

    // Retrieve output
    TableWorkspace_sptr outws = boost::dynamic_pointer_cast<TableWorkspace>(
    AnalysisDataService::Instance().retrieve("TestNPROFTable"));
    TS_ASSERT(outws);

    // Verify NPROF exists and has a value of 10 for first two banks
    map<string, double> parammap1;
    parseTableWorkspace(outws, parammap1);  // 1st bank
    map<string, double> parammap2;
    parseTableWorkspace2(outws, parammap2); // 2nd bank
    TS_ASSERT_EQUALS(parammap1.count("NPROF"),1);
    TS_ASSERT_DELTA(parammap1["NPROF"], 10.0, 0.0001);
    TS_ASSERT_EQUALS(parammap2.count("NPROF"),1);
    TS_ASSERT_DELTA(parammap2["NPROF"], 10.0, 0.0001);

    // Clean
    AnalysisDataService::Instance().remove("TestNPROFTable");
    Poco::File("TestNPROF.irf").remove();
    
  }

  //----------------------------------------------------------------------------------------------
  /** Test Exception
    */
  void test_WrongInputBankCase()
  {
    // 1. Generate file
    string filename("Test2Bank.irf");
    generate2BankIrfFile(filename);

    // 2. Load
    LoadFullprofResolution alg;
    alg.initialize();

    alg.setProperty("Filename", filename);
    alg.setPropertyValue("Banks", "2");
    alg.setProperty("OutputTableWorkspace", "TestBank3Table");

    alg.execute();

    // 3. Check if failed
    TS_ASSERT(!alg.isExecuted());

    // 4. Clean
    Poco::File("Test2Bank.irf").remove();

    return;
  }


  //----------------------------------------------------------------------------------------------
  /** Parse a TableWorkspace to a map
    */
  void parseTableWorkspace(TableWorkspace_sptr tablews, map<string, double>& parammap)
  {
    parammap.clear();

    size_t numrows = tablews->rowCount();
    for (size_t i = 0; i < numrows; ++i)
    {
      TableRow row = tablews->getRow(i);
      double value;
      string name;
      row >> name >> value;
      parammap.insert(make_pair(name, value));
    }

    return;
  }


  //----------------------------------------------------------------------------------------------
  /** Parse a TableWorkspace's 2nd bank to a map
    */
  void parseTableWorkspace2(TableWorkspace_sptr tablews, map<string, double>& parammap)
  {
    parammap.clear();

    size_t numrows = tablews->rowCount();
    for (size_t i = 0; i < numrows; ++i)
    {
      TableRow row = tablews->getRow(i);
      double value1, value2;
      string name;
      row >> name >> value1 >> value2;
      parammap.insert(make_pair(name, value2));
    }

    return;
  }


  //----------------------------------------------------------------------------------------------
  /** Generate a GEM workspace group with one member
    */
  void load_GEM()
  {
    LoadInstrument loaderGEM;

    TS_ASSERT_THROWS_NOTHING(loaderGEM.initialize());

    //create a workspace with some sample data
    wsName = "LoadFullprofResolutionWorkspace";
    Workspace_sptr ws = WorkspaceFactory::Instance().create("Workspace2D",1,1,1);
    Workspace2D_sptr ws2D = boost::dynamic_pointer_cast<Workspace2D>(ws);
    WorkspaceGroup_sptr gws(new API::WorkspaceGroup);
    gws->addWorkspace( ws2D );

    //put this workspace in the data service
    TS_ASSERT_THROWS_NOTHING(AnalysisDataService::Instance().add(wsName, gws));

    // Path to test input file 
    loaderGEM.setPropertyValue("Filename", "GEM_Definition.xml");
    //inputFile = loaderIDF2.getPropertyValue("Filename");
    loaderGEM.setPropertyValue("Workspace", wsName);
    TS_ASSERT_THROWS_NOTHING(loaderGEM.execute());
    TS_ASSERT( loaderGEM.isExecuted() );

  }

  //----------------------------------------------------------------------------------------------
  /** Generate a 1 bank .irf file
    */
  void generate1BankIrfFile(string filename)
  {
    ofstream ofile;
    ofile.open(filename.c_str());

    if (ofile.is_open())
    {
      ofile << "  Instrumental resolution function for POWGEN/SNS  A Huq  2013-12-03  ireso: 6 \n";
      ofile << "! To be used with function NPROF=10 in FullProf  (Res=6)                       \n";
      ofile << "! ----------------------------------------------  Bank 1  CWL =   0.5330A      \n";
      ofile << "!  Type of profile function: back-to-back exponentials * pseudo-Voigt          \n";
      ofile << "NPROF 10                                                                       \n";
      ofile << "!       Tof-min(us)    step      Tof-max(us)                                   \n";
      ofile << "TOFRG   5000.2300      4.0002  51000.0000                                      \n";
      ofile << "!          Zero    Dtt1                                                        \n";
      ofile << "ZD2TOF     -1.00   22580.59157                                                 \n";
      ofile << "!          Zerot   Dtt1t         Dtt2t    x-cross Width                        \n";
      ofile << "ZD2TOT  933.50214   22275.21084     1.0290  0.0000002  5.0957                  \n";
      ofile << "!     TOF-TWOTH of the bank                                                    \n";
      ofile << "TWOTH     90.00                                                                \n";
      ofile << "!           Sig-2       Sig-1       Sig-0                                      \n";
      ofile << "SIGMA     514.546      0.00044      0.355                                      \n";
      ofile << "!           Gam-2       Gam-1       Gam-0                                      \n";
      ofile << "GAMMA       0.000       0.000       0.000                                      \n";
      ofile << "!         alph0       beta0       alph1       beta1                            \n";
      ofile << "ALFBE    0.000008    6.251096    0.000000    0.000000                          \n";
      ofile << "!         alph0t      beta0t      alph1t      beta1t                           \n";
      ofile << "ALFBT   0.010156   85.918922    0.000000    0.000000                           \n";
      ofile << "END                                                                            \n";

      ofile.close();
    }
    else
    {
      throw runtime_error("Unable to open file to write.");
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Generate a 2 bank .irf file
    */
  void generate2BankIrfFile(string filename)
  {
    ofstream ofile;
    ofile.open(filename.c_str());

    if (ofile.is_open())
    {
      ofile << "  Instrumental resolution function for POWGEN/SNS  A Huq  2013-12-03  ireso: 6 \n";
      ofile << "! To be used with function NPROF=10 in FullProf  (Res=6)                       \n";
      ofile << "! ----------------------------------------------  Bank 1  CWL =   0.5330A      \n";
      ofile << "!  Type of profile function: back-to-back exponentials * pseudo-Voigt          \n";
      ofile << "NPROF 10                                                                       \n";
      ofile << "!       Tof-min(us)    step      Tof-max(us)                                   \n";
      ofile << "TOFRG   5000.2300      4.0002  51000.0000                                      \n";
      ofile << "!          Zero    Dtt1                                                        \n";
      ofile << "ZD2TOF     -1.00   22580.59157                                                 \n";
      ofile << "!          Zerot   Dtt1t         Dtt2t    x-cross Width                        \n";
      ofile << "ZD2TOT  933.50214   22275.21084     1.0290  0.0000002  5.0957                  \n";
      ofile << "!     TOF-TWOTH of the bank                                                    \n";
      ofile << "TWOTH     90.00                                                                \n";
      ofile << "!           Sig-2       Sig-1       Sig-0                                      \n";
      ofile << "SIGMA     514.546      0.00044      0.355                                      \n";
      ofile << "!           Gam-2       Gam-1       Gam-0                                      \n";
      ofile << "GAMMA       0.000       0.000       0.000                                      \n";
      ofile << "!         alph0       beta0       alph1       beta1                            \n";
      ofile << "ALFBE    0.000008    6.251096    0.000000    0.000000                          \n";
      ofile << "!         alph0t      beta0t      alph1t      beta1t                           \n";
      ofile << "ALFBT   0.010156   85.918922    0.000000    0.000000                           \n";
      ofile << "END                                                                            \n";
      ofile << "! ----------------------------------------------  Bank 3  CWL =   1.3330A\n";
      ofile << "!  Type of profile function: back-to-back exponentials * pseudo-Voigt    \n";
      ofile << "NPROF 10                                                                 \n";
      ofile << "!       Tof-min(us)    step      Tof-max(us)                             \n";
      ofile << "TOFRG   9800.0000      5.0000   86000.0000                               \n";
      ofile << "!       Zero   Dtt1                                                      \n";
      ofile << "ZD2TOF     0.00  22586.10156                                             \n";
      ofile << "!       Zerot    Dtt1t       Dtt2t    x-cross    Width                   \n";
      ofile << "ZD2TOT -42.76068   22622.76953    0.30    0.3560    2.4135               \n";
      ofile << "!     TOF-TWOTH of the bank                                              \n";
      ofile << "TWOTH    90.000                                                          \n";
      ofile << "!       Sig-2     Sig-1     Sig-0                                        \n";
      ofile << "SIGMA  72.366    10.000     0.000                                        \n";
      ofile << "!       Gam-2     Gam-1     Gam-0                                        \n";
      ofile << "GAMMA     0.000     2.742      0.000                                     \n";
      ofile << "!          alph0       beta0       alph1       beta1                     \n";
      ofile << "ALFBE        1.500      3.012      5.502      9.639                      \n";
      ofile << "!         alph0t      beta0t      alph1t      beta1t                     \n";
      ofile << "ALFBT       86.059     96.487     13.445      3.435                      \n";

      ofile.close();
    }
    else
    {
      throw runtime_error("Unable to open file to write.");
    }

    return;
  }

  /** Generate a 3 bank .irf file
    */
  void generate3BankIrfFile(string filename)
  {
    ofstream ofile;
    ofile.open(filename.c_str());

    if (ofile.is_open())
    {
      ofile << "  Instrumental resolution function for POWGEN/SNS  A Huq  2013-12-03  ireso: 6 \n";
      ofile << "! To be used with function NPROF=10 in FullProf  (Res=6)                       \n";
      ofile << "! ----------------------------------------------  Bank 2  CWL =   0.5330A      \n";
      ofile << "!  Type of profile function: back-to-back exponentials * pseudo-Voigt          \n";
      ofile << "NPROF 10                                                                       \n";
      ofile << "!       Tof-min(us)    step      Tof-max(us)                                   \n";
      ofile << "TOFRG   5000.2300      4.0002  51000.0000                                      \n";
      ofile << "!          Zero    Dtt1                                                        \n";
      ofile << "ZD2TOF     -1.00   22580.59157                                                 \n";
      ofile << "!          Zerot   Dtt1t         Dtt2t    x-cross Width                        \n";
      ofile << "ZD2TOT  933.50214   22275.21084     1.0290  0.0000002  5.0957                  \n";
      ofile << "!     TOF-TWOTH of the bank                                                    \n";
      ofile << "TWOTH     90.00                                                                \n";
      ofile << "!           Sig-2       Sig-1       Sig-0                                      \n";
      ofile << "SIGMA     514.546      0.00044      0.355                                      \n";
      ofile << "!           Gam-2       Gam-1       Gam-0                                      \n";
      ofile << "GAMMA       0.000       0.000       0.000                                      \n";
      ofile << "!         alph0       beta0       alph1       beta1                            \n";
      ofile << "ALFBE    0.000008    6.251096    0.000000    0.000000                          \n";
      ofile << "!         alph0t      beta0t      alph1t      beta1t                           \n";
      ofile << "ALFBT   0.010156   85.918922    0.000000    0.000000                           \n";
      ofile << "END                                                                            \n";
      ofile << "! ----------------------------------------------  Bank 3  CWL =   0.5339A      \n";
      ofile << "!  Type of profile function: back-to-back exponentials * pseudo-Voigt          \n";
      ofile << "NPROF 10                                                                       \n";
      ofile << "!       Tof-min(us)    step      Tof-max(us)                                   \n";
      ofile << "TOFRG   5000.2300      4.0002  51000.0000                                      \n";
      ofile << "!          Zero    Dtt1                                                        \n";
      ofile << "ZD2TOF     -1.00   22580.59157                                                 \n";
      ofile << "!          Zerot   Dtt1t         Dtt2t    x-cross Width                        \n";
      ofile << "ZD2TOT  933.50214   22275.21084     1.0290  0.0000002  5.0957                  \n";
      ofile << "!     TOF-TWOTH of the bank                                                    \n";
      ofile << "TWOTH     90.00                                                                \n";
      ofile << "!           Sig-2       Sig-1       Sig-0                                      \n";
      ofile << "SIGMA     514.546      0.00044      0.355                                      \n";
      ofile << "!           Gam-2       Gam-1       Gam-0                                      \n";
      ofile << "GAMMA       0.000       0.000       0.000                                      \n";
      ofile << "!         alph0       beta0       alph1       beta1                            \n";
      ofile << "ALFBE    0.000008    6.251096    0.000000    0.000000                          \n";
      ofile << "!         alph0t      beta0t      alph1t      beta1t                           \n";
      ofile << "ALFBT   0.010156   85.918922    0.000000    0.000000                           \n";
      ofile << "END                                                                            \n";
      ofile << "! ----------------------------------------------  Bank 4  CWL =   1.3330A\n";
      ofile << "!  Type of profile function: back-to-back exponentials * pseudo-Voigt    \n";
      ofile << "NPROF 10                                                                 \n";
      ofile << "!       Tof-min(us)    step      Tof-max(us)                             \n";
      ofile << "TOFRG   9800.0000      5.0000   86000.0000                               \n";
      ofile << "!       Zero   Dtt1                                                      \n";
      ofile << "ZD2TOF     0.00  22586.10156                                             \n";
      ofile << "!       Zerot    Dtt1t       Dtt2t    x-cross    Width                   \n";
      ofile << "ZD2TOT -42.76068   22622.76953    0.30    0.3560    2.4135               \n";
      ofile << "!     TOF-TWOTH of the bank                                              \n";
      ofile << "TWOTH    90.000                                                          \n";
      ofile << "!       Sig-2     Sig-1     Sig-0                                        \n";
      ofile << "SIGMA  72.366    10.000     0.000                                        \n";
      ofile << "!       Gam-2     Gam-1     Gam-0                                        \n";
      ofile << "GAMMA     0.000     2.742      0.000                                     \n";
      ofile << "!          alph0       beta0       alph1       beta1                     \n";
      ofile << "ALFBE        1.500      3.012      5.502      9.639                      \n";
      ofile << "!         alph0t      beta0t      alph1t      beta1t                     \n";
      ofile << "ALFBT       86.059     96.487     13.445      3.435                      \n";

      ofile.close();
    }
    else
    {
      throw runtime_error("Unable to open file to write.");
    }

    return;
  }

  /* Return the number of rows the table must have
  */
  int getExpectedNumberOfRows()
  {
    return 29;  // Change this value if you add or remove any rows from the OutputTableWorkspace
  }

  private:
  std::string wsName;  // For Workspace property

};


#endif /* MANTID_DATAHANDLING_LOADFULLPROFRESOLUTIONTEST_H_ */
