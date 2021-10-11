// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidCurveFitting/Functions/PowerLaw.h"
#include "MantidAPI/FunctionFactory.h"

namespace Mantid::CurveFitting::Functions {

using namespace CurveFitting;

using namespace Kernel;

using namespace API;

DECLARE_FUNCTION(PowerLaw)

void PowerLaw::init() {
  declareParameter("magnitude", 1.0, "coefficient for linear term");
  declareParameter("exponent", 1.0, "exponent");
  declareParameter("constant", 0.0, "coefficient for constant term");
}

void PowerLaw::function1D(double *out, const double *xValues, const size_t nData) const {
  const double a = getParameter("magnitude");
  const double b = getParameter("exponent");
  const double c = getParameter("constant");

  for (size_t i = 0; i < nData; i++) {
    out[i] = c + a * pow(xValues[i], b);
  }
}

void PowerLaw::functionDeriv1D(Jacobian *out, const double *xValues, const size_t nData) {
  const double a = getParameter("magnitude");
  const double b = getParameter("exponent");

  for (size_t i = 0; i < nData; i++) {
    double diffa = pow(xValues[i], b);
    double diffb = a * b * pow(xValues[i], b) * log(xValues[i]);
    out->set(i, 0, diffa);
    out->set(i, 1, diffb);
    out->set(i, 2, 1);
  }
}

} // namespace Mantid::CurveFitting::Functions