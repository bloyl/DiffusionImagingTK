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

#ifndef __itkRSHReconImageFilter_hxx_
#define __itkRSHReconImageFilter_hxx_

#include <cmath>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkArray.h>
#include <vnl/vnl_vector.h>
#include <itkProgressReporter.h>

#include <itkRSHReconImageFilter.h>

namespace itk
{

template< class TDiffusionModelCalculator, unsigned int TImageDimension >
RSHReconImageFilter< TDiffusionModelCalculator, TImageDimension >
::RSHReconImageFilter()
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 );
  m_ImageMask = NULL; //Must be suplied by the user
  m_DiffusionModelCalculator = NULL; //Must be suplied by the user
  
  m_ResidualImage = NULL;
  m_CalculateResidualImage = false;

}

template< class TDiffusionModelCalculator, unsigned int TImageDimension >
void
RSHReconImageFilter< TDiffusionModelCalculator, TImageDimension >
::BeforeThreadedGenerateData()
{

  itkDebugMacro( "RSHReconImageFilter::BeforeThreadedGenerateData ")

  if ( m_DiffusionModelCalculator.IsNull() )
  {
    itkExceptionMacro( << "Diffusion Calculator Not Set" );
  }

  //Initialize the m_DiffusionModelCalculator
  // m_DiffusionModelCalculator->InitializeTensorFitting();
  m_DiffusionModelCalculator->InitializeRSHFitting();

  /** Setup both fscores Image and residuals Image */
  //Initialize the residualImage if we are going to calculate it.
  if (m_CalculateResidualImage)
  {
    typename InputImageType::ConstPointer dwiImage =
          static_cast< const InputImageType * >( this->ProcessObject::GetInput(0) );
    m_ResidualImage = ResidualImageType::New();
    m_ResidualImage->CopyInformation(this->ProcessObject::GetInput(0));
    m_ResidualImage->SetRegions(m_ResidualImage->GetLargestPossibleRegion() );
    m_ResidualImage->SetVectorLength( dwiImage->GetNumberOfComponentsPerPixel() );
    m_ResidualImage->Allocate();
  }
  itkDebugMacro( "RSHReconImageFilter::BeforeThreadedGenerateData done")
}

template< class TDiffusionModelCalculator, unsigned int TImageDimension >
void
RSHReconImageFilter< TDiffusionModelCalculator, TImageDimension >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       ThreadIdType threadId)
{

  itkDebugMacro( "RSHReconImageFilter::ThreadedGenerateData Begin")
  
  //Get inputs and outputs
  OutputImagePointer outputImage = static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  InputImageConstPointer dwiImage = static_cast< const InputImageType * >( this->ProcessObject::GetInput(0) );

  //Get the diffusion Calculator
  DiffusionModelCalculatorConstPointer dtCalc = this->GetDiffusionModelCalculator();
  ImageMaskConstPointer mask = this->GetImageMask();

  //Generate iterators 
  ImageRegionConstIteratorWithIndex< InputImageType > git(dwiImage, outputRegionForThread );
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  ImageRegionIterator< ResidualImageType > res_iter;

  git.GoToBegin();
  oit.GoToBegin();

  const unsigned int numberOfGradientImages = dwiImage->GetNumberOfComponentsPerPixel();

  //if we are calculateing the Residual set up an iterator...
  if (m_CalculateResidualImage)
  {
    res_iter = ImageRegionIterator< ResidualImageType >( m_ResidualImage, outputRegionForThread);
    res_iter.GoToBegin();
  }

  //Check is mask had been provided
  bool hasMask = mask.IsNotNull();

  // Support for progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  //Local variables to check supplied mask.
  typename InputImageType::IndexType index;
  typename InputImageType::PointType point;
  
  OutputPixelType       pixValue;
  ResidualPixelType     residualValue( numberOfGradientImages );
  while( !git.IsAtEnd() )
  {
    //Check if the current location is outside of the mask
    if (hasMask)
    {
      index = git.GetIndex();
      dwiImage->TransformIndexToPhysicalPoint(index,point);
      if ( not mask->IsInside(point) )
      {
        //don't compute anything for this pixel
        oit.Set( NumericTraits<OutputPixelType>::Zero );
        ++oit;
        if (m_CalculateResidualImage)
        { 
          residualValue.Fill(NumericTraits<typename ResidualPixelType::ComponentType>::Zero);
          res_iter.Set( residualValue );
          ++res_iter;
        }
        progress.CompletedPixel();
        ++git; // Gradient  image iterator
      }
    }

    //Grab the dwi
    const InputPixelType dwi = static_cast<InputPixelType>(git.Get());
    pixValue = dtCalc->ComputeRSH(dwi);

    ///TODO need to get this from dtCalc...
    residualValue.Fill(NumericTraits<typename ResidualPixelType::ComponentType>::Zero);

    oit.Set( pixValue );
    ++oit;
    if (m_CalculateResidualImage)
    { 
      res_iter.Set( residualValue );
      ++res_iter;
    }
    progress.CompletedPixel();
    ++git; // Gradient  image iterator
  }

  itkDebugMacro( "RSHReconImageFilter::ThreadedGenerateData done")
}

template< class TDiffusionModelCalculator, unsigned int TImageDimension >
void
RSHReconImageFilter< TDiffusionModelCalculator, TImageDimension >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "m_CalculateResidualImage: " << m_CalculateResidualImage << std::endl;
  if (m_CalculateResidualImage)
  {
    os << indent << "m_ResidualImage: " << std::endl << indent << indent << m_ResidualImage << std::endl;
  }
  os << indent << "NEEDS MORE INFO" << std::endl;
}

} //end Itk namespace

#endif
