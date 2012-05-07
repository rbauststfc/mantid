/*WIKI* 
Calculate the EQSANS detector sensitivity. This workflow algorithm uses the
reduction parameters found in the property manager object passed as the
ReductionProperties parameter to load the given data file, apply all the
necessary corrections to it and compute the sensitivity correction.

Setting the PatchWorkspace property allows you to patch areas of the
detector. All masked pixels in the patch workspace will be patched.
The value assigned to a patched pixel is the average of all unmasked
pixels in this patched pixel's tube.
*WIKI*/
//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidWorkflowAlgorithms/ComputeSensitivity.h"
#include "MantidAPI/FileProperty.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidWorkflowAlgorithms/EQSANSInstrument.h"
#include "MantidAPI/AlgorithmProperty.h"
#include "MantidAPI/PropertyManagerDataService.h"
#include "MantidKernel/PropertyManager.h"

namespace Mantid
{
namespace WorkflowAlgorithms
{

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(ComputeSensitivity)

/// Sets documentation strings for this algorithm
void ComputeSensitivity::initDocs()
{
  this->setWikiSummary("Workflow to calculate EQSANS sensitivity correction.");
  this->setOptionalMessage("Workflow to calculate EQSANS sensitivity correction.");
}

using namespace Kernel;
using namespace API;
using namespace Geometry;
using namespace DataObjects;

void ComputeSensitivity::init()
{
  declareProperty(new API::FileProperty("Filename", "", API::FileProperty::Load, "_event.nxs"),
      "Flood field or sensitivity file.");
  declareProperty(new WorkspaceProperty<>("PatchWorkspace", "", Direction::Input, PropertyMode::Optional),
      "Workspace defining the area of the detector to be patched. All masked pixels in this workspace will be patched.");
  declareProperty("ReductionProperties", "__eqsans_reduction_properties", Direction::Input);
  declareProperty(new WorkspaceProperty<>("OutputWorkspace", "", Direction::Output),
      "Workspace containing the sensitivity correction.");
  declareProperty("OutputMessage", "", Direction::Output);
}

void ComputeSensitivity::exec()
{
  std::string outputMessage = "";
  progress(0.1, "Setting up sensitivity calculation");

  // Reduction property manager
  const std::string reductionManagerName = getProperty("ReductionProperties");
  boost::shared_ptr<PropertyManager> reductionManager;
  if (PropertyManagerDataService::Instance().doesExist(reductionManagerName))
  {
    reductionManager = PropertyManagerDataService::Instance().retrieve(reductionManagerName);
  }
  else
  {
    g_log.notice() << "Could not find property manager" << std::endl;
    reductionManager = boost::make_shared<PropertyManager>();
    PropertyManagerDataService::Instance().addOrReplace(reductionManagerName, reductionManager);
  }

  const std::string outputWS = getPropertyValue("OutputWorkspace");

  // Find beam center
  if (reductionManager->existsProperty("SANSBeamFinderAlgorithm"))
  {
    //const std::string algTxt = reductionManager->getPropertyValue("SANSBeamFinderAlgorithm");

    IAlgorithm_sptr ctrAlg = reductionManager->getProperty("SANSBeamFinderAlgorithm");
    ctrAlg->setPropertyValue("ReductionProperties", reductionManagerName);
    ctrAlg->setChild(true);
    ctrAlg->execute();
    std::string outMsg2 = ctrAlg->getPropertyValue("OutputMessage");
    outputMessage += outMsg2;
  }

  progress(0.2, "Computing sensitivity");

  // Set patch information so that the SANS sensitivity algorithm can
  // patch the sensitivity workspace
  const std::string patchWSName = getPropertyValue("PatchWorkspace");
  if (patchWSName.size()>0)
  {
    IAlgorithm_sptr patchAlg = createSubAlgorithm("EQSANSPatchSensitivity");
    patchAlg->setPropertyValue("PatchWorkspace", patchWSName);
    if (!reductionManager->existsProperty("SensitivityPatchAlgorithm"))
    {
      reductionManager->declareProperty(new AlgorithmProperty("SensitivityPatchAlgorithm"));
    }
    reductionManager->setProperty("SensitivityPatchAlgorithm", patchAlg);
  }

  if (reductionManager->existsProperty("SensitivityAlgorithm"))
  {
    const std::string fileName = getPropertyValue("Filename");
    IAlgorithm_sptr effAlg = reductionManager->getProperty("SensitivityAlgorithm");
    effAlg->setChild(true);
    effAlg->setProperty("Filename", fileName);
    effAlg->setPropertyValue("OutputSensitivityWorkspace", outputWS);
    effAlg->execute();
    MatrixWorkspace_sptr effWS = effAlg->getProperty("OutputSensitivityWorkspace");
    setProperty("OutputWorkspace", effWS);
    std::string outMsg2 = effAlg->getPropertyValue("OutputMessage");
    outputMessage += outMsg2;
    setProperty("OutputMessage", outputMessage);
  } else {
    g_log.error() << "Could not find sensitivity algorithm" << std::endl;
  }
}

} // namespace WorkflowAlgorithms
} // namespace Mantid

