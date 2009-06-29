#include "Instrument3DWidget.h"
#include "InstrumentActor.h"
#include "MantidObject.h"
#include "GLColorMap.h"

#include "MantidAPI/IInstrument.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/SpectraDetectorMap.h"

#include "MantidKernel/Exception.h"

#include "MantidGeometry/Matrix.h"
#include "MantidGeometry/V3D.h"
#include "MantidGeometry/IComponent.h"
#include "MantidGeometry/IObjComponent.h"
#include "MantidGeometry/ICompAssembly.h"
#include "MantidGeometry/Object.h"
#include "MantidGeometry/GeometryHandler.h"

#include "boost/shared_ptr.hpp"

#include <QTimer>
#include <QMessageBox>
#include <QString>

#include <map>
#include <cmath>
#include <cfloat>

using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::Geometry;

static const QRgb BLACK = qRgb(0,0,0);

Instrument3DWidget::Instrument3DWidget(QWidget* parent):
  GL3DWidget(parent),mFastRendering(true), iTimeBin(0), mDataMapping(INTEGRAL),
  mColorMap(GLColorMapQwt()), mInstrumentActor(NULL), mAxisDirection(Mantid::Geometry::V3D(0.0,0.0,1.0)), 
  mAxisUpVector(Mantid::Geometry::V3D(0.0,1.0,0.0)), mDataMinValue(DBL_MAX), mDataMaxValue(-DBL_MAX), 
  mBinMinValue(DBL_MAX), mBinMaxValue(-DBL_MAX), mWkspDataMin(DBL_MAX), mWkspDataMax(-DBL_MAX), 
  mWkspBinMin(DBL_MAX), mWkspBinMax(-DBL_MAX), strWorkspaceName(""), mScaledValues(0)
{
  connect(this, SIGNAL(actorsPicked(const std::set<QRgb>&)), this, SLOT(fireDetectorsPicked(const std::set<QRgb>&)));
  connect(this, SIGNAL(actorHighlighted(QRgb)),this,SLOT(fireDetectorHighligted(QRgb)));
}

Instrument3DWidget::~Instrument3DWidget()
{
	makeCurrent();
}

/**
 * Set the default Axis direction of the model
 */
void Instrument3DWidget::setAxis(const Mantid::Geometry::V3D& direction,const Mantid::Geometry::V3D& up)
{
	mAxisDirection = direction;
	mAxisUpVector = up;
}

/**
 * This method is the slot when the detectors are picked using mouse. This method emits
 * signals the ids of the detector and the spectra index(not spectra number).
 * @param pickedActor the input passed by the the signal.
 */
void Instrument3DWidget::fireDetectorsPicked(const std::set<QRgb>& pickedColors)
{
  std::vector<int> detectorIds;
  for(std::set<QRgb>::const_iterator it = pickedColors.begin(); it!= pickedColors.end(); it++)
  {
    int iDecId = mInstrumentActor->getDetectorIDFromColor(qRed((*it))*65536+qGreen((*it))*256+qBlue((*it)));
    if(iDecId != -1)
    {
      detectorIds.push_back(iDecId);
    }
    
  }
  //convert detector ids to spectra index ids
  std::vector<int> spectraIndices;
  getSpectraIndexList(detectorIds, spectraIndices);
  if( !detectorIds.empty() )
  {
    if( detectorIds.size() == 1)
    {
      //emit the detector id
      emit actionDetectorSelected(detectorIds.front());
      //emit the spectra id
      emit actionSpectraSelected(spectraIndices.front());
    }
    else // If more than one detector selected
    {
      //emit the detector ids
      emit actionDetectorSelectedList(detectorIds);
      //emit the spectra ids
      emit actionSpectraSelectedList(spectraIndices);
    }
    
  }
}

/**
 * This method is the slot when the detector is highlighted using mouse move. This method emits
 * signals the id of the detector and the spectra index(not spectra number).
 * @param pickedActor the input passed by the the signal.
 */
void Instrument3DWidget::fireDetectorHighligted(QRgb pickedColor)
{
  if(pickedColor == BLACK)
  {
    emit actionDetectorHighlighted(-1,-1,-1);
    return;
  }
  int iDecId = mInstrumentActor->getDetectorIDFromColor(qRed(pickedColor)*65536 + qGreen(pickedColor)*256 + qBlue(pickedColor));
  if(iDecId != -1)
  {
    //convert detector id to spectra index id
    std::vector<int> idDecVec(1, iDecId);
    std::vector<int> indexList;
    getSpectraIndexList(idDecVec, indexList);
    MatrixWorkspace_sptr workspace = boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(strWorkspaceName));
    int spectrumNumber(1);
    try
    {
      spectrumNumber = workspace->getAxis(1)->spectraNo(indexList.front());
    }
    catch(Mantid::Kernel::Exception::IndexError&)
    {
      //Not a Workspace2D
    }

    double sum = integrateSingleSpectra(workspace, indexList.front());
    //emit the detector id, spectrum number and count to display in the window
    emit actionDetectorHighlighted(iDecId, spectrumNumber, std::floor(sum));
  }

}

/**
 * This method sets the workspace name input to the widget.
 * @param wsName input workspace name
 */
void Instrument3DWidget::setWorkspace(const std::string& wsName)
{
  if( strWorkspaceName == wsName ) return;

  // Save the workspace name
  strWorkspaceName = wsName;
  MatrixWorkspace_sptr output = boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(strWorkspaceName));
  // Read the instrument geometry
  boost::shared_ptr<Mantid::API::IInstrument> ins = output->getInstrument();
  this->ParseInstrumentGeometry(ins);
  boost::shared_ptr<Mantid::Geometry::IObjComponent> sample = ins->getSample();
  if( sample.get() )
  {
    _trackball->setModelCenter(sample->getPos());
  }
  else
  {
    _trackball->setModelCenter(Mantid::Geometry::V3D(0.0,0.0,0.0));
  }
  defaultProjection(); // Calculate and set projection

  // Calculate bin values, data ranges and integrate data
  calculateBinRange(output);
  calculateColorCounts(output);

}


/**
 * This method parses the instrument information and creates the actors relating to the detectors.
 */
void Instrument3DWidget::ParseInstrumentGeometry(boost::shared_ptr<Mantid::API::IInstrument> ins)
{
	makeCurrent();
	boost::shared_ptr<GLActorCollection> scene = boost::shared_ptr<GLActorCollection>(new GLActorCollection);
	mInstrumentActor = new InstrumentActor(ins, mFastRendering);
	scene->addActor(mInstrumentActor);	 
	this->setActorCollection(scene);
}

/**
 * Calculate the minimum and maximum values of the bins for the set workspace
 */
void Instrument3DWidget::calculateBinRange(Mantid::API::MatrixWorkspace_sptr workspace)
{
  const int nHist = workspace->getNumberHistograms();
  mWkspBinMin = DBL_MAX;
  mWkspBinMax = -DBL_MAX;
  for (int i = 0; i < nHist; ++i)
  {
    const std::vector<double> & values = workspace->readX(i);
    double xtest = values.front();
    if( xtest < mWkspBinMin )
    {
      mWkspBinMin = xtest;
    }
    else if( xtest > mWkspBinMax )
    {
      mWkspBinMax = xtest;
    }
    else {}

    xtest = values.back();
    if( xtest < mWkspBinMin )
    {
      mWkspBinMin = xtest;
    }
    else if( xtest > mWkspBinMax )
    {
      mWkspBinMax = xtest;
    }
    else {}
  }

  // Value has not been preset
  if( (mBinMinValue - DBL_MAX)/DBL_MAX < 1e-08 )
  {
    mBinMinValue = mWkspBinMin;
  }

  // Value has not been preset
  if( (mBinMaxValue + DBL_MAX)/DBL_MAX < 1e-08 )
  {
    mBinMaxValue = mWkspBinMax;
  }

  // Check validity
  if( mBinMinValue < mWkspBinMin || mBinMinValue > mWkspBinMax )
  {
      mBinMinValue = mWkspBinMin;
  }
  
  if( mBinMaxValue > mWkspBinMax || mBinMaxValue < mWkspBinMin )
  {
    mBinMaxValue = mWkspBinMax;
  }
  

}

/**
 * Integrate the workspace
 */
void Instrument3DWidget::calculateColorCounts(boost::shared_ptr<Mantid::API::MatrixWorkspace> workspace)
{
  // This looks like a strange way of doing this but the CompAssemblyActor needs the colours in the same
  // order as it fills its detector lists!
  std::vector<int> detector_list(0);
  mInstrumentActor->getDetectorIDList(detector_list);
  if( detector_list.empty() ) return;
  std::vector<int> index_list;
  getSpectraIndexList(detector_list, index_list);
  
  const int n_spec = index_list.size();
  std::vector<double> integrated_values(n_spec, 0.0);
  mWkspDataMin = DBL_MAX;
  mWkspDataMax = -DBL_MAX;
  for( int i = 0; i < n_spec; ++i )
  {
    int widx = index_list[i];
    if( widx != -1 )
    {
//       std::vector<double>::const_iterator bin_itr = workspace->readX(widx).begin();
//       std::vector<double>::const_iterator bin_end = workspace->readX(widx).end();
//       std::vector<double>::const_iterator data_itr = workspace->readY(widx).begin();
//       std::vector<double>::const_iterator data_end = workspace->readY(widx).end();
//       double sum(0.0);
//       for( ; data_itr != data_end; ++data_itr, ++bin_itr )
//       {
// 	double binvalue = *bin_itr;
// 	if( binvalue >= mBinMinValue && binvalue <= mBinMaxValue )
// 	{
// 	  sum += *data_itr;
// 	}
//       }
//      integrated_values[i] = sum;
      double sum = integrateSingleSpectra(workspace, widx);
      integrated_values[i] = sum;
      if( sum < mWkspDataMin )
      {
	mWkspDataMin = sum;
      }
      else if( sum > mWkspDataMax )
      {
	mWkspDataMax = sum;
      }
      else continue;
      
    }
    else
    {
      integrated_values[i] = -1.0;
    }
  }

  // No preset value
  if( std::abs(mDataMinValue - DBL_MAX)/DBL_MAX < 1e-08 )
  {
    mDataMinValue = mWkspDataMin;
  }

  if( (mDataMaxValue + DBL_MAX)/DBL_MAX < 1e-08 )
  {
    mDataMaxValue = mWkspDataMax;
  }

  // This is the maximum number of colours that we allow for any colour map
  const short max_ncols = mColorMap.getMaxNumberOfColors();
  // This is the number of colours in the current colour map
  const short no_colors = mColorMap.getNumberOfColors();


  //Store values scaled to 1->256
  double wksp_datarange = std::abs(mWkspDataMax - mWkspDataMin);

  if( wksp_datarange < 1e-08 )
  {
    wksp_datarange = 1.0;
  }
  double user_datarange = std::abs(mDataMaxValue - mDataMinValue);
  if( user_datarange < 1e-08 )
  {
    user_datarange = 1.0;
  }

  mScaledValues = std::vector<short>(n_spec, 0);
  std::vector<boost::shared_ptr<GLColor> > colorlist(n_spec);
  std::vector<double>::const_iterator val_end = integrated_values.end();
  int idx(0);
  for( std::vector<double>::const_iterator val_itr = integrated_values.begin(); val_itr != val_end; 
       ++val_itr, ++idx )
  {
    short c_index(no_colors);
    if( (*val_itr) < 0.0 ) 
    {
      mScaledValues[idx] = -1;
    }
    else
    {
      short cache_value = std::floor(1.0 + ((*val_itr - mWkspDataMin)*max_ncols/wksp_datarange) );
      if( cache_value > max_ncols ) cache_value = max_ncols;
      mScaledValues[idx]  = cache_value;

      // Now compute a color index for this color map
      int adjusted_index = std::floor(1.0 + ((*val_itr - mDataMinValue)*no_colors/user_datarange) );
      if( adjusted_index <= (int)no_colors && adjusted_index > 0) 
      {
	c_index = adjusted_index;
      }
      // This will happen for very small ranges
      else if( adjusted_index <= 0 )
      {
	c_index = 1;
      }
      else
      {
	c_index = no_colors;
      }

    }
    colorlist[idx] = mColorMap.getColor(c_index - 1);
  }
  mInstrumentActor->setDetectorColors(colorlist);
}

double Instrument3DWidget::integrateSingleSpectra(Mantid::API::MatrixWorkspace_sptr workspace, const int wksp_index)
{
  std::vector<double>::const_iterator bin_itr = workspace->readX(wksp_index).begin();
  std::vector<double>::const_iterator bin_end = workspace->readX(wksp_index).end();
  std::vector<double>::const_iterator data_end = workspace->readY(wksp_index).end();
  double sum(0.0);
  for( std::vector<double>::const_iterator data_itr = workspace->readY(wksp_index).begin();
       data_itr != data_end; ++data_itr, ++ bin_itr )
  {
    double binvalue = *bin_itr;
    if( binvalue >= mBinMinValue && binvalue <= mBinMaxValue )
    {
      sum += *data_itr;
    }
  }
  return sum;
}

/**
 * For a change in the colour map, just update the color indices
 */
void Instrument3DWidget::updateColorsForNewMap()
{
  double wksp_datarange = std::abs(mWkspDataMax - mWkspDataMin);
  double current_range = std::abs(mDataMaxValue - mDataMinValue);

  const short no_colors = mColorMap.getNumberOfColors();
  const short max_ncols = mColorMap.getMaxNumberOfColors();

  std::vector<boost::shared_ptr<GLColor> > colorlist(mScaledValues.size());
  std::vector<short>::const_iterator val_end = mScaledValues.end();
  int idx(0);
  for( std::vector<short>::const_iterator val_itr = mScaledValues.begin(); val_itr != val_end; 
       ++val_itr, ++idx )
  {
    //Value between 1->256
    short cache_value = *val_itr;
    short c_index(no_colors);
    if( cache_value >= 0 )
    {
      double x_i = ((double)(cache_value - 1) * wksp_datarange / max_ncols ) + mWkspDataMin;
      int adjusted_index = std::floor(1.0 + ((x_i - mDataMinValue)*no_colors/current_range) );
      if( adjusted_index <= (int)no_colors && adjusted_index > 0) 
      {
	c_index = adjusted_index;
      }
      else if( adjusted_index <= 0 )
      {
	c_index = 1;
      }
      else
      {
	c_index = no_colors;
      }

    }
    colorlist[idx] = mColorMap.getColor(c_index - 1);
  }
  mInstrumentActor->setDetectorColors(colorlist);
  mInstrumentActor->refresh();
  update();
}

/**
 * Update the colors based on a change in the maximum data value
 */
void Instrument3DWidget::updateForNewMaxData(const double new_max)
{
  if( std::abs(new_max - mDataMaxValue) / mDataMaxValue < 1e-08 ) return;
  // ratio: old / new
  double range_ratio = std::abs(mWkspDataMax - mWkspDataMin) / std::abs(new_max - mDataMinValue);
  if( std::isnan(range_ratio) || std::isinf(range_ratio) ) range_ratio = 0.0;
  const short no_colors = mColorMap.getNumberOfColors();
  const short max_ncols = mColorMap.getMaxNumberOfColors();

  std::vector<boost::shared_ptr<GLColor> > colorlist(mScaledValues.size());
  std::vector<short>::const_iterator val_end = mScaledValues.end();
  int idx(0);
  for( std::vector<short>::const_iterator val_itr = mScaledValues.begin(); val_itr != val_end; 
       ++val_itr, ++idx )
  {
    short cache_value = *val_itr;
    short c_index(no_colors);
    if( cache_value >= 0 )
    {
      int adjusted_index = std::ceil((double)cache_value * range_ratio * no_colors / max_ncols);
      if( adjusted_index <= (int)no_colors && adjusted_index > 0) 
      {
	c_index = adjusted_index;
      }
      // This will happen for very small ranges
      else if( adjusted_index <= 0 )
      {
	c_index = 1;
      }
      else
      {
	c_index = no_colors;
      }
    }
    colorlist[idx] = mColorMap.getColor(c_index - 1);
  }

  mInstrumentActor->setDetectorColors(colorlist);
  mInstrumentActor->refresh();
  update();

  // Keep the new value
  mDataMaxValue = new_max;
}

/**
 * Update the colors based on a change in the maximum data value
 */
void Instrument3DWidget::updateForNewMinData(const double new_min)
{
  // Don't do anything if no change from the current value
  if( std::abs(new_min - mDataMinValue) < 1e-08 ) return; 

  // ratio: old / new
  double wksp_datarange = std::abs(mWkspDataMax - mWkspDataMin);
  double new_range = std::abs(mDataMaxValue - new_min);
  
  const short no_colors = mColorMap.getNumberOfColors();
  const short max_ncols = mColorMap.getMaxNumberOfColors();

  std::vector<boost::shared_ptr<GLColor> > colorlist(mScaledValues.size());
  std::vector<short>::const_iterator val_end = mScaledValues.end();
  int idx(0);
  for( std::vector<short>::const_iterator val_itr = mScaledValues.begin(); val_itr != val_end; 
       ++val_itr, ++idx )
  {
    short cache_value = *val_itr;
    short c_index(no_colors);
    if( cache_value >= 0 )
    {
      double x_i = ((double)(cache_value - 1) * wksp_datarange / max_ncols ) + mWkspDataMin;
      int adjusted_index = std::floor(1.0 + ((x_i - new_min)*no_colors/new_range) );
      if( adjusted_index <= (int)no_colors && adjusted_index > 0) 
      {
	c_index = adjusted_index;
      }
      else if( adjusted_index <= 0 )
      {
	c_index = 1;
      }
      else
      {
	c_index = no_colors;
      }

    }
    colorlist[idx] = mColorMap.getColor(c_index - 1);
  }
  mInstrumentActor->setDetectorColors(colorlist);
  mInstrumentActor->refresh();
  update();

  // Keep the new value
  mDataMinValue = new_min;

}


/**
 * This method returns the Spectra Index list for the input dectector id list.
 * @param idDecVec is list of detector id's
 */
void Instrument3DWidget::getSpectraIndexList(const std::vector<int>& detIDs, std::vector<int> & wkspIndices) const
{
  wkspIndices.clear();
  if( strWorkspaceName.empty() ) return;
  MatrixWorkspace_sptr workspace = boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(strWorkspaceName));
  // There is no direct way of getting histogram index from the spectra id,
  // get the spectra axis and convert form index to spectra number and create
  // a map.
  const std::vector<int> spectraList = workspace->spectraMap().getSpectra(detIDs);
  Axis* spectraAxis = workspace->getAxis(1);
  std::map<int,int> index_map;
  int n_hist = workspace->getNumberHistograms();
  for (int i = 0; i < n_hist; ++i)
  {
    int current_spectrum = spectraAxis->spectraNo(i);
    index_map[current_spectrum] = i;
  }

  std::vector<int>::const_iterator spec_end = spectraList.end();
  std::vector<int>::const_iterator d_itr = detIDs.begin();
  for( std::vector<int>::const_iterator spec_itr = spectraList.begin(); spec_itr != spec_end;
       ++spec_itr, ++d_itr )
  {
    if( (*d_itr) != -1 )
    {
      wkspIndices.push_back(index_map[*spec_itr]);
    }
    else
    {
      wkspIndices.push_back(-1);
    }
  }
}

/**
 * This method assigns colors to the dectors using the values
 * @param minval input minimum value for scaling of the colormap
 * @param maxval input maximum value for scaling of the colormap
 * @param values input values that coorespond to each detector
 * @param colorMap input color map class which is used for looking up the color to be assigned based on the value of detector.
 */
void Instrument3DWidget::setColorForDetectors(double minval,double maxval,const std::vector<double>& values,const GLColorMap& colorMap)
{
	std::vector<boost::shared_ptr<GLColor> > iColorList;
	int noOfColors=colorMap.getNumberOfColors();
	std::vector<double>::size_type nvals = values.size();
  for(std::vector<double>::size_type i = 0; i < nvals; i++)
	{
		int cIndex;
		if(maxval-minval<0.00000001)
		{
		  cIndex = 0;
		}
		else
		{
		  cIndex=floor(((values[i]-minval)/(maxval-minval))*(noOfColors-1));
		}
		if(cIndex<0)
		{
		  cIndex = noOfColors - 1;
		}
		else if(cIndex>(noOfColors-1))
		{
		  cIndex=(noOfColors-1);
		}
		iColorList.push_back(colorMap.getColor(cIndex));

	} // Looping through the dectors/Actors list
	mInstrumentActor->setDetectorColors(iColorList);
}

/**
 * This method takes the input name as the filename and reads the file for the color index values.
 * NOTE: This method can only read 256 color index with RGB values
 */
void Instrument3DWidget::setColorMapName(const QString & name)
{
  mColorMap.setColorMapFile(name.toStdString());
}

/**
 * This method sets the Time bin values. the value has to be greater than zero
 * @param value input is the time bin value
 */
void Instrument3DWidget::setTimeBin(int value)
{
	if(value>0)
	{
		this->iTimeBin=value;
	}
}


/**
 * Returns workspace name
 */
std::string Instrument3DWidget::getWorkspaceName() const
{
	return strWorkspaceName;
}

/**
 * Returns Colormap
 */
GLColorMapQwt Instrument3DWidget::getColorMap()const
{
	return mColorMap;
}

/**
 * This method takes the input name as the min value Colormap scale.
 */
void Instrument3DWidget::setColorMapMinValue(double minValue)
{
	this->mDataMinValue=minValue;
}

/**
 * This method takes the input name as the min value Colormap scale.
 */
void Instrument3DWidget::setColorMapMaxValue(double maxValue)
{
	this->mDataMaxValue = maxValue;
}

/**
 * This method returns min value. by default will be min value in the current timebin
 */
double Instrument3DWidget::getDataMinValue() const
{
	return this->mDataMinValue;
}

/**
 * This method returns the max value. by default will be max value in the current timebin
 */
double Instrument3DWidget::getDataMaxValue() const
{
	return this->mDataMaxValue;
}

/**
 * Returns the current minimum bin value
 */
double Instrument3DWidget::getBinMinValue() const
{
  return this->mBinMinValue;
}

double Instrument3DWidget::getBinMaxValue() const
{
    return this->mBinMaxValue;
}

/**
 * This method sets the Data mapping type for the color mapping.
 */
void Instrument3DWidget::setDataMappingType(DataMappingType dmType)
{
	mDataMapping=dmType;
}

void Instrument3DWidget::setDataMappingIntegral(double minValue,double maxValue)
{
	this->mBinMinValue = minValue;
	this->mBinMaxValue = maxValue;
	setDataMappingType(INTEGRAL);
	if( this->isVisible() ) 
	{
	  MatrixWorkspace_sptr workspace = 
	    boost::dynamic_pointer_cast<MatrixWorkspace>(AnalysisDataService::Instance().retrieve(strWorkspaceName));
	  calculateColorCounts(workspace);
	  mInstrumentActor->refresh();
	  update();
	}
}

void Instrument3DWidget::setDataMappingSingleBin(int binNumber)
{
	this->iTimeBin=binNumber;
	setDataMappingType(SINGLE_BIN);
}

/**
 * Sets the default view to X direction positive
 */
void Instrument3DWidget::setViewDirectionXPositive()
{
	setViewDirection(XPOSITIVE);
}

/**
 * Sets the default view to Y direction positive
 */
void Instrument3DWidget::setViewDirectionYPositive()
{
	setViewDirection(YPOSITIVE);
}

/**
 * Sets the default view to Z direction positive
 */
void Instrument3DWidget::setViewDirectionZPositive()
{
	setViewDirection(ZPOSITIVE);
}

/**
 * Sets the default view to X direction negative
 */
void Instrument3DWidget::setViewDirectionXNegative()
{
	setViewDirection(XNEGATIVE);
}

/**
 * Sets the default view to Y direction negative
 */
void Instrument3DWidget::setViewDirectionYNegative()
{
	setViewDirection(YNEGATIVE);
}

/**
 * Sets the default view to Z direction negative
 */
void Instrument3DWidget::setViewDirectionZNegative()
{
	setViewDirection(ZNEGATIVE);
}

/**
 * Sets the slow rendering not using display list
 * NOTE: This method will ***NOT*** have any effect after the workspace name is set.
 */
void Instrument3DWidget::setSlowRendering()
{
	mFastRendering=false;
}

/**
 * Sets the fast rendering using display list
 * NOTE: This method will ***NOT*** have any effect after the workspace name is set.
 */
void Instrument3DWidget::setFastRendering()
{
	mFastRendering=true;
}


/**
 * Completely resets the data in the instrument widget. ready for new workspace
 */
void Instrument3DWidget::resetWidget()
{
	iTimeBin = 0;
	strWorkspaceName = "";
	mBinMinValue = DBL_MAX;
	mBinMaxValue = -DBL_MAX;
	mDataMinValue = DBL_MAX;
	mDataMaxValue = -DBL_MAX;
	mDataMapping = INTEGRAL;
	mAxisDirection = Mantid::Geometry::V3D(0.0,0.0,1.0);
	mAxisUpVector = Mantid::Geometry::V3D(0.0,1.0,0.0);
	mScaledValues.clear();
	GL3DWidget::resetWidget();
}

void Instrument3DWidget::setView(const V3D& pos,double xmax,double ymax,double zmax,double xmin,double ymin,double zmin)
{
	////get the centre of the bounding box
	//V3D boundCentre;
	//double xhalf=(xmax-xmin)/2.0;
	//double yhalf=(ymax-ymin)/2.0;
	//double zhalf=(zmax-zmin)/2.0;
	//boundCentre[0]=-1*(xmin+xhalf);
	//boundCentre[1]=-1*(ymin+yhalf);
	//boundCentre[2]=-1*(zmin+zhalf);
	////vector from center to bounding box center
	//V3D vcb=pos-boundCentre;
	//vcb.normalize();
	//V3D zaxis(0,0,1);
	////get the rotation about zaxis
	//Quat rotation(zaxis,vcb);
	//rotation.inverse();
	//_trackball->reset();
	//_trackball->setModelCenter(pos);
	//if(rotation!=Quat(0,0,0,0))
	//	_trackball->setRotation(rotation);
	//_trackball->rotateBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);//rotate the bounding box
	//_viewport->setOrtho(xmin,xmax,ymin,ymax,-zmin,-zmax);
	//_viewport->issueGL();
	//update();

	//Change the View to the axis orientation
	V3D s=mAxisDirection.cross_prod(mAxisUpVector);
	V3D u=s.cross_prod(mAxisDirection);
	double Mat[16];
	Mat[0]=s[0];              Mat[4]=s[1];              Mat[8]=s[2];               Mat[12]=0;
	Mat[1]=u[0];              Mat[5]=u[1];              Mat[9]=u[2];               Mat[13]=0;
	Mat[2]=-mAxisDirection[0];Mat[6]=-mAxisDirection[1];Mat[10]=-mAxisDirection[2];Mat[14]=0;
	Mat[3]=0;                 Mat[7]=0;                 Mat[11]=0;                 Mat[15]=1;
	Quat defaultView;
	defaultView.setQuat(Mat);
	defaultView.normalize();
	//get the rotation to make the center of the bounding box to the view
	V3D boundCentre;
	boundCentre[0]=(xmax+xmin)/2.0;
	boundCentre[1]=(ymax+ymin)/2.0;
	boundCentre[2]=(zmax+zmin)/2.0;
	V3D vcb=boundCentre-pos;
	vcb.normalize();
	V3D zaxis(0,0,1);
	Quat rotation(zaxis,vcb);
	rotation.inverse();
	if(rotation!=Quat(0,0,0,0))
		defaultView=rotation*defaultView;
	_trackball->reset();
	_trackball->setModelCenter(pos);
	if(defaultView!=Quat(0,0,0,0))
		_trackball->setRotation(defaultView);
	_trackball->rotateBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);//rotate the bounding box
	_viewport->setOrtho(xmin,xmax,ymin,ymax,-zmax,-zmin);
	_viewport->issueGL();
	update();
}

/**
 * This method draws the scene using color id. this method is called in pick mode
 */
void Instrument3DWidget::drawSceneUsingColorID()
{
	mInstrumentActor->drawUsingColorID();
}

/**
 * Draws the scene in low resolution. this method is called in interactive mode for faster response
 */
void Instrument3DWidget::setSceneLowResolution()
{
	mInstrumentActor->setObjectResolutionToLow();
}

/**
 * Draws the scene in high resolution. 
 */
void Instrument3DWidget::setSceneHighResolution()
{
	mInstrumentActor->setObjectResolutionToHigh();
}

/**
 * Returns the boundig box of the scene
 * @param minBound :: output min point of the bounding box of scene
 * @param maxBound :: output max point of the bounding box of scene
 */ 
void Instrument3DWidget::getBoundingBox(Mantid::Geometry::V3D& minBound, Mantid::Geometry::V3D& maxBound)
{
	mInstrumentActor->getBoundingBox(minBound,maxBound);
}
