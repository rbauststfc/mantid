// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidAlgorithms/AddAbsorptionWeightedPathLengths.h"
#include "MantidAPI/ExperimentInfo.h"
#include "MantidAPI/Sample.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAlgorithms/SampleCorrections/MCAbsorptionStrategy.h"
#include "MantidAlgorithms/SampleCorrections/RectangularBeamProfile.h"
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidGeometry/Instrument/SampleEnvironment.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/Material.h"
#include "MantidKernel/MersenneTwister.h"

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Geometry;
using namespace Mantid::Kernel;

namespace Mantid {
namespace Algorithms {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(AddAbsorptionWeightedPathLengths)

//----------------------------------------------------------------------------------------------

namespace {

constexpr int DEFAULT_NEVENTS = 1000;
constexpr int DEFAULT_SEED = 123456789;
} // namespace

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void AddAbsorptionWeightedPathLengths::init() {
  declareProperty(
      std::make_unique<WorkspaceProperty<PeaksWorkspace_sptr::element_type>>(
          "InputWorkspace", "", Direction::InOut),
      "An input/output peaks workspace that the path distances will be added "
      "to.");
  auto positiveInt = std::make_shared<Kernel::BoundedValidator<int>>();
  positiveInt->setLower(1);
  declareProperty("EventsPerPoint", DEFAULT_NEVENTS, positiveInt,
                  "The number of \"neutron\" events to generate per peak");
  declareProperty("SeedValue", DEFAULT_SEED, positiveInt,
                  "Seed the random number generator with this value");
  declareProperty("MaxScatterPtAttempts", 5000, positiveInt,
                  "Maximum number of tries made to generate a scattering point "
                  "within the sample. Objects with "
                  "holes in them, e.g. a thin annulus can cause problems "
                  "if this number is too low.\n"
                  "If a scattering point cannot be generated by increasing "
                  "this value then there is most likely a problem with "
                  "the sample geometry.");
}

std::map<std::string, std::string>
AddAbsorptionWeightedPathLengths::validateInputs() {
  PeaksWorkspace_sptr inputWS = getProperty("InputWorkspace");
  std::map<std::string, std::string> issues;
  Geometry::IComponent_const_sptr sample =
      inputWS->getInstrument()->getSample();
  if (!sample) {
    issues["InputWorkspace"] = "Input workspace does not have a Sample";
  } else {
    if (inputWS->sample().hasEnvironment()) {
      issues["InputWorkspace"] = "Sample must not have a sample environment";
    }

    if (inputWS->sample().getMaterial().numberDensity() == 0) {
      issues["InputWorkspace"] =
          "Sample must have a material set up with a non-zero number density";
    }
  }
  return issues;
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void AddAbsorptionWeightedPathLengths::exec() {

  const PeaksWorkspace_sptr inputWS = getProperty("InputWorkspace");
  const int nevents = getProperty("EventsPerPoint");
  const int maxScatterPtAttempts = getProperty("MaxScatterPtAttempts");

  auto instrument = inputWS->getInstrument();
  auto beamProfile = createBeamProfile(*instrument, inputWS->sample());

  const auto npeaks = inputWS->getNumberPeaks();

  // Configure progress
  Progress prog(this, 0.0, 1.0, npeaks);
  prog.setNotifyStep(0.01);
  const std::string reportMsg = "Computing path lengths";

  // Configure strategy
  const int nlambda = 1;
  MCAbsorptionStrategy strategy(*beamProfile, inputWS->sample(),
                                DeltaEMode::Elastic, nevents,
                                maxScatterPtAttempts, true, g_log);

  const int seed = getProperty("SeedValue");
  MersenneTwister rng(seed);

  for (int i = 0; i < npeaks; ++i) {
    IPeak &peak = inputWS->getPeak(i);
    auto peakWavelength = peak.getWavelength();

    std::vector<double> lambdas{peakWavelength}, absFactors(nlambda),
        absFactorErrors(nlambda);

    strategy.calculate(rng, peak.getDetectorPosition(), lambdas, peakWavelength,
                       absFactors, absFactorErrors);
    double mu = inputWS->sample().getMaterial().attenuationCoefficient(
        peakWavelength);                                               // m-1
    double absWeightedPathLength = -log(absFactors[0]) / mu;           // metres
    peak.setAbsorptionWeightedPathLength(absWeightedPathLength * 100); // cm

    prog.report(reportMsg);
  }
}

/**
 * Create the beam profile. Currently only supports Rectangular. The dimensions
 * are either specified by those provided by `SetBeam` algorithm or default
 * to the width and height of the samples bounding box
 * @param instrument A reference to the instrument object
 * @param sample A reference to the sample object
 * @return A new IBeamProfile object
 */
std::unique_ptr<IBeamProfile>
AddAbsorptionWeightedPathLengths::createBeamProfile(
    const Instrument &instrument, const Sample &sample) const {
  const auto frame = instrument.getReferenceFrame();
  const auto source = instrument.getSource();

  auto beamWidthParam = source->getNumberParameter("beam-width");
  auto beamHeightParam = source->getNumberParameter("beam-height");
  double beamWidth(-1.0), beamHeight(-1.0);
  if (beamWidthParam.size() == 1 && beamHeightParam.size() == 1) {
    beamWidth = beamWidthParam[0];
    beamHeight = beamHeightParam[0];
  } else {
    const auto bbox = sample.getShape().getBoundingBox().width();
    beamWidth = bbox[frame->pointingHorizontal()];
    beamHeight = bbox[frame->pointingUp()];
  }
  return std::make_unique<RectangularBeamProfile>(*frame, source->getPos(),
                                                  beamWidth, beamHeight);
}

} // namespace Algorithms
} // namespace Mantid