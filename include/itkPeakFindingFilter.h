/**
 * @file  itkPeakFindingFilter.h
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPeakFindingFilter_h
#define __itkPeakFindingFilter_h


#include <itkSmartPointer.h>
#include <itkImageToImageFilter.h>

#include <itkPeakFindingCalculator.h>


namespace itk
{


template < class TInputImage, class TOutputComponentType = double >
class ITK_EXPORT PeakFindingFilter: 
          public ImageToImageFilter
                 <    TInputImage,
                      Image  <  Vector < TOutputComponentType, 3 >,
                                TInputImage::ImageDimension
                             >
                 >
{
public:
  

  typedef Image
          <  Vector < TOutputComponentType, 3 >,            // 3 components per direction
             TInputImage::ImageDimension
          >                                                 TOutputImage;

  // Standard class typedefs
  typedef PeakFindingFilter                                 Self;
  typedef ImageToImageFilter < TInputImage, TOutputImage >  Superclass;
  typedef SmartPointer < Self >                             Pointer;
  typedef SmartPointer < const Self >                       ConstPointer;

  // Image typedefs
  typedef typename Superclass::InputImageType               InputImageType;
  typedef typename InputImageType::PixelType                InputPixelType;
  typedef typename InputPixelType::ComponentType            InputPixelComponentType; // type of SH coefficient
  typedef typename InputImageType::ConstPointer             InputImageConstPointer;
  typedef typename Superclass::OutputImageType              OutputImageType;
  typedef typename OutputImageType::RegionType              OutputImageRegionType;
  typedef typename OutputImageType::PixelType               OutputPixelType;
  typedef typename OutputImageType::Pointer                 OutputImagePointer;

  typedef PeakFindingCalculator 
          < InputPixelType, 
            TOutputComponentType, TOutputComponentType >    PixelPeakFinderType;

  typedef typename PixelPeakFinderType::Pointer             PixelPeakFinderPointer; 

  itkStaticConstMacro ( InputPixelTupleSize, int,           // number of the SH coefficients for each pixel
                        InputPixelType::Dimension );  

  // Method for creation through the object factory
  itkNewMacro ( Self );

  // Run-time type information (and related methods)
  itkTypeMacro ( PeakFindingFilter, ImageToImageFilter );

  // Set the number of peaks expected of each voxel
  // (zero-vector padding upon any peak unavailability)
  void SetNumberOfPeaksPerVoxel( int numPeaks = 1 );

  // Set whether the vector of peak direction is normalized (1) or not (0)
  void SetVectorNormalization  ( int normaliz = 1 )  
       {  this->normalizeVector = normaliz;  }


protected:


  PeakFindingFilter();
 ~PeakFindingFilter() { }

  // the engine of the filter
  void ThreadedGenerateData 
       ( const OutputImageRegionType & outputRegionForThread, int threadId );

  // invoked, with a special treatment, if the number of peaks per voxel is 
  // specified by explicitly calling SetNumberOfPeaksPerVoxel( ... )
  void ExplicitGenerateData 
       ( const OutputImageRegionType & outputRegionForThread, int threadId );

  // invoked if the number of peaks per voxel is determined without calling
  // SetNumberOfPeaksPerVoxel( ... )
  void ImplicitGenerateData
       ( const OutputImageRegionType & outputRegionForThread, int threadId );

  // create an  output image object
  DataObject::Pointer CreateOutputImage();

  // get the specified image object
  TOutputImage   *    GetOutputImage( unsigned int index );


private:


  PeakFindingFilter ( const Self & ); // purposely not implemented
  void operator  =  ( const Self & ); // purposely not implemented

  int                       numPeaksBeenSet;  // number of peaks has been EXPLICITLY set ?
  int                       numPeaksEachVox;
  int                       normalizeVector;
  PixelPeakFinderPointer    pixelPeakFinder;


};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPeakFindingFilter.hxx"
#endif


#endif
