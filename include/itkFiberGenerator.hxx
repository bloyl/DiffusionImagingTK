/**
* @file  itkFiberGenerator.hxx
* @brief This class is used to generate single fiber tracts from
* diffuion images such as DTI or FOD images 
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkFiberGenerator_hxx
#define __itkFiberGenerator_hxx

#include <itkFiberGenerator.h>

#include <itkBinaryThresholdImageFilter.h>
#include <itkNumericTraits.h>

namespace itk
{
template <class TDirPicker, typename TLabelType >
FiberGenerator<TDirPicker,TLabelType>
::FiberGenerator()
{
  this->m_StepLength = -1;
  this->m_CurvatureThreshold = -1;
  this->m_DirectionPicker = NULL;

//  this->m_Points = VtkPointsPointer::New();
  this->m_StoppingMask  = StoppingMaskType::New();
  this->m_ROIImage = NULL;
//  this->m_TraversedROIs = LabelVectorType::New();

  this->m_Initialized = false;
}

template <class TDirPicker, typename TLabelType >
void
FiberGenerator<TDirPicker,TLabelType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Step Length : " << this->m_StepLength << std::endl;
  os << indent << "CurvatureThreshold: " << this->m_CurvatureThreshold << std::endl;
  if (this->m_DirectionPicker)
  {
    os << indent << "DirectionPicker : " << this->m_DirectionPicker << std::endl;
  }
  else
  {
    os << indent << "DirectionPicker : NOT SET" << std::endl;
  }
  //TODO! ADD print support of sometype.
  // os << indent << "Points: " << this->m_Points << std::endl;
  // os << indent << indent << this->m_Points->GetNumberOfPoints() << std::endl;
  // double pt[3];
  // for (int i = 0; i < this->m_Points->GetNumberOfPoints(); ++i)
  // {
  //   this->m_Points->GetPoint(i,pt);
  //   os << indent << indent << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
  // }
  
  os << indent << "Traversed Rois " << std::endl << indent;
  typename LabelSetType::const_iterator it = m_TraversedROIs.begin();
  while(it != m_TraversedROIs.end())
  {
    os << *it << " , " ;
    ++it;
  }
  os << std::endl;
  os << std::endl;
}

template <class TDirPicker, typename TLabelType >
void
FiberGenerator<TDirPicker,TLabelType>
::AddStoppingCriteriaSpatialObject(typename SpatialObject<3>::Pointer mask)
{
  this->m_StoppingMask->AddSpatialObject(mask);
  this->Modified();
}

template <class TDirPicker, typename TLabelType >
void
FiberGenerator<TDirPicker,TLabelType>
::GenerateFiberPoints( const PointType seed, const StartDirectionChoice dirChoice)
{
  if (! this->m_Initialized)
    itkExceptionMacro( "Error fiber Genererator Not initialized");

  //Check the seed point!
  if ( not this->m_DirectionPicker->CanEvaluate(seed) ) 
    itkExceptionMacro( "Error Seed point " << seed << " outside of evaluatable region of Direction Picker!");

  DirectionType startDir;
  try
  {
    startDir = this->m_DirectionPicker->PickStartingDirection(seed);
  }
  catch (itk::ExceptionObject err)
  {
    itkWarningMacro("GenerateFiberPoints error"<<std::endl<<"Picking starting direction at " << seed << std::endl << err.what());
    return;
  }
  GenerateFiberPoints(seed, startDir, dirChoice);
}

template <class TDirPicker, typename TLabelType >
void
FiberGenerator<TDirPicker,TLabelType>
::GenerateFiberPoints( const PointType seed, const DirectionType startDir, const StartDirectionChoice dirChoice)
{
  if (! this->m_Initialized)
    itkExceptionMacro( "Error fiber Genererator Not initialized");

  if ( not this->m_DirectionPicker->CanEvaluate(seed) ) 
    itkExceptionMacro( "Error Seed point " << seed << " outside of evaluatable region of Direction Picker!");

  this->m_Points.clear();
  this->m_TraversedROIs.clear();

  if (dirChoice == ONE_DIRECTION)
  {
    GenerateFiberPointsOneWay( seed, startDir,this->m_Points);
    return;
  }
  else if (dirChoice == BOTH_DIRECTIONS)
  {
    LinePntsContainerType fwdPoints;
    LinePntsContainerType bckPoints;
    GenerateFiberPointsOneWay( seed, startDir, fwdPoints);
    GenerateFiberPointsOneWay( seed, -startDir, bckPoints);

    //Combine the two.
    this->m_Points.reserve( fwdPoints.size() + bckPoints.size()-1 );

    for(typename LinePntsContainerType::reverse_iterator rit = bckPoints.rbegin(); rit != bckPoints.rend(); ++rit)
    {
      this->m_Points.push_back(*rit);
    }

    //Skip the first point...
    for(typename LinePntsContainerType::iterator it = fwdPoints.begin()+1; it != fwdPoints.end(); ++it)
    {
      this->m_Points.push_back(*it);
    }
    return;
  }
  
}

template <class TDirPicker, typename TLabelType >
void
FiberGenerator<TDirPicker,TLabelType>
::GenerateFiberPointsOneWay( const PointType seed, const DirectionType startDir, LinePntsContainerType &locPoints)
{
  //Clear the points
  // locPoints->Reset();
  locPoints.clear();
  
  StoppingMaskConstPointer stopMask = this->GetStoppingMask();
  DirectionPickerConstPointer dirPicker = this->GetDirectionPicker();
  double step    = this->GetStepLength();
  double kThresh = this->GetCurvatureThreshold();
  
  typename LabelImageType::IndexType roiIndex;

  PointType     curPoint  = seed;
  DirectionType curDir    = startDir;
  DirectionType lastDir   = curDir;

  while (true)
  {
    //curvature is too high or point is outside of image
    //   don't add curPoint and STOP.
    if (not dirPicker->CanEvaluate(curPoint)
        || this->ComputeCurvature(lastDir,curDir) > kThresh )
    {
      break;
    }
    
    //add the current point
    // locPoints->InsertNextPoint(curPoint[0],curPoint[1],curPoint[2]);
    LinePointType linePoint;
    linePoint.SetPosition(curPoint);
    locPoints.push_back(linePoint);
    
    //attach labels
    if (this->m_ROIImage)
    {
      this->m_ROIImage->TransformPhysicalPointToIndex(curPoint,roiIndex);
      LabelPixelType roiLabel = this->m_ROIImage->GetPixel(roiIndex);
      if (roiLabel != NumericTraits<LabelPixelType>::Zero)
      {
        if (m_TraversedROIs.find(roiLabel) == m_TraversedROIs.end())
        {
          m_TraversedROIs.insert(roiLabel);
        }
      }
    } //end attachLabel

    //Check if we should stop.
    if ( stopMask->IsInside(curPoint,99999) )
    {
      break;
    }
    //Update the current Point
    lastDir = curDir;
    try
    {
      curDir = dirPicker->PickNextDirection(curPoint,curDir);
    }
    catch (itk::ExceptionObject err)
    {
      itkWarningMacro("GenerateFiberPointsOneWay error"<<std::endl<<"Picking next direction at " << curPoint << std::endl << err.what());
      return;
    }

    curPoint = curPoint + step * curDir;
  }
}

template <class TDirPicker, typename TLabelType >
double
FiberGenerator<TDirPicker,TLabelType>
::ComputeCurvature( const DirectionType inDir, const DirectionType outDir) const
{
  double dotProd = 0;
  for (unsigned int i=0;i<3;++i)
    dotProd += inDir[i]*outDir[i];

  return (vcl_acos(dotProd)/ this->m_StepLength);
}

template <class TDirPicker, typename TLabelType >
typename FiberGenerator<TDirPicker,TLabelType>::LinePointer
FiberGenerator<TDirPicker,TLabelType>
::GetFiber( ) const
{
  LineType::Pointer line = LineType::New();
  std::vector<LineType::LinePointType> points = this->m_Points;
  
  //line->SetPoints(this->m_Points);
  line->SetPoints(points);
  return line;
}

// template <class TDirPicker, typename TLabelType >
// void
// FiberGenerator<TDirPicker,TLabelType>
// ::Modified( )
// {
//   Superclass::Modified();
//   this->m_Initialized = false;
// }


template <class TDirPicker, typename TLabelType >
void
FiberGenerator<TDirPicker,TLabelType>
::Initialize( )
{
  if (this->m_StepLength <=0)
    itkExceptionMacro("Can not initialize FiberGenerator. StepLength must be >0");

  if (this->m_CurvatureThreshold <=0)
    itkExceptionMacro("Can not initialize FiberGenerator. CurvatureThreshold must be >0");

  if (not this->m_DirectionPicker)
    itkExceptionMacro("Can not initialize FiberGenerator. DirectionPicker has not be set");

  if( not this->m_DirectionPicker->IsInitialized() )
    itkExceptionMacro("Can not initialize FiberGenerator. DirectionPicker is not intialized");
  
  this->m_Initialized = true;
}

template <class TDirPicker, typename TLabelType >
double
FiberGenerator<TDirPicker,TLabelType>
::GetFiberLength( ) const
{
  return (this->m_StepLength * this->m_Points.size());
}

}// end namespace itk

#endif
