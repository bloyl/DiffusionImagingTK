/*=========================================================================
 *
 *  Copyright Section of Biomedical Image Analysis
 *            University of Pennsylvania
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/


#ifndef __itkFiberTrackingManger_hxx
#define __itkFiberTrackingManger_hxx

#include <itkFiberTrackingManager.h>

#include <vnl/vnl_sample.h>

namespace itk
{
template <class TDirPicker, typename TLabelType >
FiberTrackingManager<TDirPicker,TLabelType>
::FiberTrackingManager()
{
  this->m_StoppingMasks   = GroupSpatialObjectType::New();
  this->m_SeedMasks       = GroupSpatialObjectType::New();
  this->m_SeedMasks->SetBoundingBoxChildrenDepth(9999);
  this->m_SeedPoints      = PointSetType::New();
  this->m_Fibers          = GroupSpatialObjectType::New();
  this->m_RejectedFibers  = GroupSpatialObjectType::New();
  this->m_FibersROIs      = VectorOfLabelSetsType::New();
  
  //Defaults
  this->m_NumberOfFibersToGen  = 0;
  this->m_NumberOfFibersToKeep = 10;
  this->m_StepLength          = 0.2;
  this->m_CurvatureThreshold  = 3.5355;
  this->m_MinimumFiberLength  = 0;
  
  //Default to all 
  this->m_NumberOfRequiredIncludeROIs = -1;
  
  vnl_sample_reseed();
}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

//  os << indent << "Stopping masks : " << this->m_StoppingMasks << std::endl;
//  os << indent << "Seed masks     : " << this->m_SeedMasks << std::endl;

  const unsigned int numberOfPoints = this->m_SeedPoints->GetNumberOfPoints();
  os << indent << "Seed Points    : " << this->m_SeedPoints << std::endl;
  os << indent << indent << " Num Points    : " << numberOfPoints << std::endl;

}

//TODO expand to multithread...
template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::GenerateFibers( )
{

  PointType seedPoint;

  //Make + initialize fiber generator
  FiberGeneratorPointer ftGen = FiberGeneratorType::New();
  ftGen->SetDirectionPicker(this->GetDirectionPicker());

  ftGen->SetStepLength(this->GetStepLength());
  ftGen->SetCurvatureThreshold(this->GetCurvatureThreshold());
  ftGen->AddStoppingCriteriaSpatialObject(this->GetStoppingMasks());

  ftGen->SetROIImage(this->GetROIImage());

  const double minFiberLength = this->GetMinimumFiberLength();

  bool keepFiber = true;
  
  const LabelSetType incROIS = this->GetIncludeROIs();
  const LabelSetType reqIncROIS = this->GetRequiredIncludeROIs();
  const LabelSetType excROIS = this->GetExcludeROIs();
  const int numReqInc = this->GetNumberOfRequiredIncludeROIs();
  
  //Iterator
  typename LabelSetType::const_iterator fibRoi_it;

  ftGen->Initialize();
  
  unsigned int keptFibers = 0;
  unsigned int keepNFibers = this->GetNumberOfFibersToKeep();
  if (keepNFibers <=0)
  {
    itkExceptionMacro("Error numberOfFibersToKeep ("<<keepNFibers<<") must be positive")
  }

  unsigned int maxNFibers = this->GetNumberOfFibersToGen();
  if (maxNFibers <=keepNFibers)
  {
    maxNFibers = keepNFibers * 100;
  }

  //TODO this for loop should be multithreaded...
  //ftgen isn't thread safe although....
  // but I see no reason not to make seperate ones in each thread.
  for (unsigned int genCounter = 0; genCounter < maxNFibers; genCounter ++)
  {
    seedPoint = GenerateSeedPoint();

    //Track from seed
    ftGen->GenerateFiberPoints(seedPoint,FiberGeneratorType::BOTH_DIRECTIONS);
    const LabelSetType fibROIs = ftGen->GetTraversedROIs();

    //Check if generated track meets out criteria    
    keepFiber = true;
    if (ftGen->GetFiberLength() < minFiberLength)
    {
      keepFiber = false;
    }
    else //Check include/exclude dirs...
    {
      unsigned int incRoiCounter = 0;
      unsigned int reqIncRoiCounter = 0;

      //Any Exclude ROIs is a problem...
      fibRoi_it = fibROIs.begin();
      while (fibRoi_it != fibROIs.end())
      {
        //Any Exlucde Rois are a problem
        if ( excROIS.find(*fibRoi_it) != excROIS.end() )
        {
          keepFiber = false;
          break;
        }
        //count the number of include rois this fibers touches...
        if( incROIS.find(*fibRoi_it) != incROIS.end() )
        {
          incRoiCounter++;
        }
        //count the number of Required rois this fiber touches...
        if( reqIncROIS.find(*fibRoi_it) != reqIncROIS.end() )
        {
          reqIncRoiCounter++;
        }
        fibRoi_it++;
      }

      if ( keepFiber ) //Maybe its already rejected for exclusion ROIs
      { //Fibers must cross all of the required Include ROIs
        keepFiber = ( reqIncRoiCounter == reqIncROIS.size() );
      }

      if ( keepFiber ) //Maybe its already rejected for exclusion ROIs
      {
        if (numReqInc == -1) // Must meed all includsion Rois
        {
          keepFiber = ( incRoiCounter == incROIS.size() );
        }
        else
        {
          keepFiber = ( static_cast<int>(incRoiCounter) >= numReqInc );
        }
      }
    }

    //set seedPoint data accordingly.
    // 0 = unsuccessful seed
    // 1 = successful  seed
    if (keepFiber)
    {
      this->m_SeedPoints->SetPoint(genCounter,seedPoint);
      this->m_SeedPoints->SetPointData(genCounter,1);
      
      //Add fiber to output
      typename FiberType::Pointer fiber = ftGen->GetFiber();
      fiber->SetId(this->m_Fibers->GetNumberOfChildren()+1);
      this->m_Fibers->AddSpatialObject(fiber);
      this->m_FibersROIs->InsertElement(keptFibers,fibROIs);
      keptFibers++;
      if (keptFibers >= keepNFibers)
      {
        break;
      }
    }
    else
    {
      this->m_SeedPoints->SetPoint(genCounter,seedPoint);
      this->m_SeedPoints->SetPointData(genCounter,0);
      this->m_RejectedFibers->AddSpatialObject(ftGen->GetFiber());
    }
  }
}

template <class TDirPicker, typename TLabelType >
typename FiberTrackingManager<TDirPicker,TLabelType>::PointType
FiberTrackingManager<TDirPicker,TLabelType>
::GenerateSeedPoint( ) const
{

  //Get the bounding Box.
  typename GroupSpatialObjectType::BoundingBoxType::Pointer
            boundingBox = this->m_SeedMasks->GetBoundingBox();

  typename GroupSpatialObjectType::PointType
        maxPoint = boundingBox->GetMaximum();
  typename GroupSpatialObjectType::PointType
        minPoint = boundingBox->GetMinimum();

  PointType pt;
  do
  {
    pt[0] = vnl_sample_uniform( minPoint[0], maxPoint[0] );
    pt[1] = vnl_sample_uniform( minPoint[1], maxPoint[1] );
    pt[2] = vnl_sample_uniform( minPoint[2], maxPoint[2] );
  }
  while ( !this->m_SeedMasks->IsInside(pt,9999) );

  return pt;
}

template <class TDirPicker, typename TLabelType >
template <typename  TSeedPixelType>
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddSeedImage(const typename Image<TSeedPixelType,3>::Pointer image, const TSeedPixelType lower_thresh, const TSeedPixelType upper_thresh)
{

  typedef ImageMaskSpatialObject<3>                         MaskType;

  typedef itk::Image< TSeedPixelType , 3 >                  ImageMaskType;
  
  typedef itk::BinaryThresholdImageFilter< ImageMaskType, typename MaskType::ImageType >
                                                            ThresholderType;

  //Set up thresholder
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetOutsideValue(itk::NumericTraits< typename MaskType::ImageType::PixelType>::Zero);
  thresholder->SetInsideValue(itk::NumericTraits< typename MaskType::ImageType::PixelType>::One);
  
  //the threshold is inclusive so is >= 1
  thresholder->SetLowerThreshold(lower_thresh);
  thresholder->SetUpperThreshold(upper_thresh);
  thresholder->InPlaceOn();
  
  thresholder->SetInput( image );
  thresholder->Update();

  typename MaskType::Pointer  spatialObjectMask = MaskType::New();
  spatialObjectMask->SetImage( thresholder->GetOutput() );
  spatialObjectMask->Update();
  spatialObjectMask->DisconnectPipeline();

  this->m_SeedMasks->AddSpatialObject(spatialObjectMask);
  this->m_SeedMasks->ComputeBoundingBox();

}


template <class TDirPicker, typename TLabelType >
template <typename  TMaskPixelType>
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddStoppingCriteria(typename Image<TMaskPixelType,3>::Pointer image, TMaskPixelType lower_thresh, TMaskPixelType upper_thresh)
{

  typedef ImageMaskSpatialObject<3>                         MaskType;

  typedef itk::Image< TMaskPixelType , 3 >                  ImageMaskType;
  
  typedef itk::BinaryThresholdImageFilter< ImageMaskType, typename MaskType::ImageType >
                                                            ThresholderType;

  //Set up thresholder
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetOutsideValue(itk::NumericTraits< typename MaskType::ImageType::PixelType>::Zero);
  thresholder->SetInsideValue(itk::NumericTraits< typename MaskType::ImageType::PixelType>::One);
  
  //the threshold is inclusive so is >= 1
  thresholder->SetLowerThreshold(lower_thresh);
  thresholder->SetUpperThreshold(upper_thresh);
  thresholder->InPlaceOn();
  
  thresholder->SetInput( image );
  thresholder->Update();

  typename MaskType::Pointer  spatialObjectMask = MaskType::New();
  spatialObjectMask->SetImage( thresholder->GetOutput() );
  spatialObjectMask->DisconnectPipeline();
  
  this->m_StoppingMasks->AddSpatialObject(spatialObjectMask);

}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddLabelToIncludeRois( const LabelPixelType roiLabel)
{
  if (m_IncludeROIs.find(roiLabel) == m_IncludeROIs.end())
  {
    m_IncludeROIs.insert(roiLabel);
  }
}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddLabelToIncludeRois( const LabelSetType roiLabels)
{
  typename LabelSetType::const_iterator it = roiLabels.begin();
  while(it != roiLabels.end())
  {
    AddLabelToIncludeRois(*it);
    ++it;
  }
}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddLabelToRequiredIncludeRois( const LabelPixelType roiLabel)
{
  if (m_RequiredIncludeROIs.find(roiLabel) == m_RequiredIncludeROIs.end())
  {
    m_RequiredIncludeROIs.insert(roiLabel);
  }
}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddLabelToRequiredIncludeRois( const LabelSetType roiLabels)
{
  typename LabelSetType::const_iterator it = roiLabels.begin();
  while(it != roiLabels.end())
  {
    AddLabelToRequiredIncludeRois(*it);
    ++it;
  }
}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddLabelToExcludeRois( const LabelPixelType roiLabel)
{
  if (m_ExcludeROIs.find(roiLabel) == m_ExcludeROIs.end())
  {
    m_ExcludeROIs.insert(roiLabel);
  }
}

template <class TDirPicker, typename TLabelType >
void
FiberTrackingManager<TDirPicker,TLabelType>
::AddLabelToExcludeRois( const LabelSetType roiLabels)
{
  typename LabelSetType::const_iterator it = roiLabels.begin();
  while(it != roiLabels.end())
  {
    AddLabelToExcludeRois(*it);
    ++it;
  }
}


}// end namespace itk

#endif
