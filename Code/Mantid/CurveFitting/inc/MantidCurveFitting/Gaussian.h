#ifndef MANTID_CURVEFITTING_GAUSSIAN_H_
#define MANTID_CURVEFITTING_GAUSSIAN_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/IPeakFunction.h"

namespace Mantid
{
  namespace CurveFitting
  {
    /**
    Provide gaussian peak shape function interface to IPeakFunction.
    I.e. the function: Height*exp(-0.5*((x-PeakCentre)/Sigma)^2).

    This function actually performs the fitting on 1/Sigma^2 rather than Sigma
    for stability reasons.

    Gauassian parameters:
    <UL>
    <LI> Height - height of peak (default 0.0)</LI>
    <LI> PeakCentre - centre of peak (default 0.0)</LI>
    <LI> Sigma - standard deviation (default 1.0)</LI>
    </UL>

    @author Anders Markvardsen, ISIS, RAL
    @date 19/10/2009

    Copyright &copy; 2007-8 STFC Rutherford Appleton Laboratory

    This file is part of Mantid.

    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>
    Code Documentation is available at: <http://doxygen.mantidproject.org>
    */
    class DLLExport Gaussian : public API::IPeakFunction
    {
    public:
      /// Destructor
      virtual ~Gaussian() {};


      /// overwrite IPeakFunction base class methods
      virtual double centre()const {return getParameter("PeakCentre");};
      virtual double height()const {return getParameter("Height");};
      virtual double width()const {return 2*getParameter("Sigma");};
      virtual void setCentre(const double c) {getParameter("PeakCentre") = c;};
      virtual void setHeight(const double h) {getParameter("Height") = h;};
      virtual void setWidth(const double w) {getParameter("Sigma") = w/2.0;};


      /// overwrite IFunction base class methods
      virtual void init();
      virtual void calJacobianForCovariance(API::Jacobian* out, const double* xValues, const int& nData);
      virtual void setActiveParameter(int i,double value);
      virtual double activeParameter(int i);
      virtual void function(double* out, const double* xValues, const int& nData);
      virtual void functionDeriv(API::Jacobian* out, const double* xValues, const int& nData);


    };

  } // namespace CurveFitting
} // namespace Mantid

#endif /*MANTID_CURVEFITTING_GAUSSIAN_H_*/
