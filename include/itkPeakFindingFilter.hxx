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

#ifndef __itkPeakFindingFilter_hxx
#define __itkPeakFindingFilter_hxx


#include <itkObjectFactory.h>
#include <itkProgressReporter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <itkPeakFindingFilter.h>

namespace itk
{


template < class TInputImage, class TOutputComponentType >
PeakFindingFilter < TInputImage, TOutputComponentType >
::PeakFindingFilter()
{
  this->numPeaksBeenSet = 0;
  this->numPeaksEachVox = 1;
  this->normalizeVector = 1;
  this->pixelPeakFinder = PeakFindingFilter // a smart pointer
                          < TInputImage, TOutputComponentType >
                          ::PixelPeakFinderType::New();
}


// create an output image object
template < class TInputImage, class TOutputComponentType >
DataObject::Pointer PeakFindingFilter < TInputImage, TOutputComponentType >
::CreateOutputImage()
{
  DataObject::Pointer output = ( TOutputImage::New() ).GetPointer();
  return output.GetPointer();
}


// get a specified output image object
template < class TInputImage, class TOutputComponentType >
typename PeakFindingFilter < TInputImage, TOutputComponentType >::
TOutputImage * PeakFindingFilter < TInputImage, TOutputComponentType >
::GetOutputImage( unsigned int index )
{
  return dynamic_cast < TOutputImage * >
                      (  this->ProcessObject::GetOutput( index )  );
}


// set the number of peaks expected of each voxel
// (zero-vector padding upon any peak unavailability)
template < class TInputImage, class TOutputComponentType >
void PeakFindingFilter < TInputImage, TOutputComponentType >
::SetNumberOfPeaksPerVoxel( int numPeaks )
{  
  if ( this->numPeaksBeenSet ) return;

  this->numPeaksEachVox = numPeaks;
  this->SetNumberOfRequiredOutputs( this->numPeaksEachVox );
  for ( int i = 0; i < this->numPeaksEachVox; i ++ )
  {
    this->SetNthOutput( i, this->CreateOutputImage() );
  }

  this->numPeaksBeenSet = 1;
}


// the engine of the filter
template < class TInputImage, class TOutputComponentType >
void PeakFindingFilter < TInputImage, TOutputComponentType >
::ThreadedGenerateData ( const OutputImageRegionType & outputRegionForThread,
                         int threadId )
{
  itkDebugMacro( "PeakFindingFilter::ThreadedGenerateData()" );

  if ( this->numPeaksBeenSet )
       this->ExplicitGenerateData( outputRegionForThread, threadId );
  else this->ImplicitGenerateData( outputRegionForThread, threadId );
}


// invoked if the number of peaks per voxel is determined without calling
// SetNumberOfPeaksPerVoxel( ... )
template < class TInputImage, class TOutputComponentType >
void PeakFindingFilter < TInputImage, TOutputComponentType >
::ImplicitGenerateData ( const OutputImageRegionType &,
                         int threadId )
{
  itkDebugMacro( "PeakFindingFilter::ImplicitGenerateData()" );

  TOutputComponentType    odValue,  normVec[3];

  InputImageConstPointer  inputPtr  = this->GetInput ( 0 );
  OutputImagePointer      outputPtr = this->GetOutput( 0 );

  ProgressReporter        progress( this, threadId, 
                          outputPtr->GetRequestedRegion().GetNumberOfPixels() );

  typedef ImageRegionConstIterator < TInputImage  >    InputIterator;
  typedef ImageRegionIterator      < TOutputImage >    OutputIterator;

  InputIterator  inputIt  = InputIterator ( inputPtr,  inputPtr->GetRequestedRegion () );
  OutputIterator outputIt = OutputIterator( outputPtr, outputPtr->GetRequestedRegion() );
 
  inputIt.GoToBegin ();
  outputIt.GoToBegin();

  while ( !outputIt.IsAtEnd() )
  {
    InputPixelType  pix_val0 = inputIt.Get();
    OutputPixelType pix_val1;
    pix_val1.Fill(0);

    for ( int i = 0; i < InputPixelTupleSize; i ++ )
    {
      this->pixelPeakFinder->SetCoefficient(  i,  pix_val0.GetNthComponent( i )  );
    }
    
    this->pixelPeakFinder->GetPeak( 0, normVec, odValue, 0 );
    odValue     = this->normalizeVector ? 1.0 : odValue;
    pix_val1[0] = normVec[0] * odValue;
    pix_val1[1] = normVec[1] * odValue;
    pix_val1[2] = normVec[2] * odValue;
    outputIt.Set( pix_val1 );

    ++ inputIt;
    ++ outputIt;
  }

  itkDebugMacro( << "SymRSHPowerImageFilter::ImplicitGenerateData() Complete" );
}


// invoked, with a special treatment, if the number of peaks per voxel is 
// specified by explicitly calling SetNumberOfPeaksPerVoxel( ... )
template < class TInputImage, class TOutputComponentType >
void PeakFindingFilter < TInputImage, TOutputComponentType >
::ExplicitGenerateData ( const OutputImageRegionType &,
                         int threadId )
{ 
  itkDebugMacro( "PeakFindingFilter::ExplicitGenerateData()" );

  int     i;
  TOutputComponentType   odValue, normVec[3];

  typedef ImageRegionConstIterator < TInputImage  >  InputIterator;
  typedef ImageRegionIterator      < TOutputImage >  OutputIterator;

  InputImageConstPointer inputPtr = this->GetInput( 0 );
  InputIterator          inputIt  = InputIterator
                                    ( inputPtr, inputPtr->GetRequestedRegion() );
  inputIt.GoToBegin ();
  
  OutputImagePointer   * outputImages = new OutputImagePointer [ this->numPeaksEachVox ];
  OutputIterator       * outputIts    = new OutputIterator     [ this->numPeaksEachVox ];
  for ( i = 0; i < this->numPeaksEachVox; i ++ )
  {
    outputImages[i] = SmartPointer < TOutputImage >
                                   (  this->GetOutputImage( i )  );
    outputIts[i]    = OutputIterator
                      ( outputImages[i], outputImages[i]->GetRequestedRegion() );
    outputIts[i].GoToBegin();
  }

  ProgressReporter    progress( this, threadId, 
                      outputImages[0]->GetRequestedRegion().GetNumberOfPixels() );

  while ( !outputIts[0].IsAtEnd() )
  {
    InputPixelType  pix_val0 = inputIt.Get();
    OutputPixelType pix_val1;
    pix_val1.Fill(0);

    for ( i = 0; i < InputPixelTupleSize; i ++ )
    {
      this->pixelPeakFinder->SetCoefficient(  i,  pix_val0.GetNthComponent( i )  );
    }
    
    for ( i = 0; i < this->numPeaksEachVox; i ++ )
    {
      this->pixelPeakFinder->GetPeak( i, normVec, odValue, 0 );
      odValue     = this->normalizeVector ? 1.0 : odValue;
      pix_val1[0] = normVec[0] * odValue;
      pix_val1[1] = normVec[1] * odValue;
      pix_val1[2] = normVec[2] * odValue;
      
      outputIts[i].Set( pix_val1 );
      ++ outputIts[i];
    }
 
    ++ inputIt;
  }

  for ( i = 0; i < this->numPeaksEachVox; i ++ ) outputImages[i] = NULL;

  delete [] outputIts;
  delete [] outputImages;
  outputIts    = NULL;
  outputImages = NULL;
  
  itkDebugMacro( << "SymRSHPowerImageFilter::ExplicitGenerateData() Complete" );
}


} // end namespace itk


#endif
