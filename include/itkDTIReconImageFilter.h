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

#ifndef __itkDTIReconImageFilter_h_
#define __itkDTIReconImageFilter_h_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_svd.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>
#include <itkSpatialObject.h>
#include <itkImageToImageFilter.h>

namespace itk
{

template< class TDiffusionModelCalculator,
          unsigned int TImageDimension=3
          >
class ITK_EXPORT DTIReconImageFilter :
  public ImageToImageFilter< VectorImage< typename TDiffusionModelCalculator::DWIPixelType::ComponentType, TImageDimension >,
                             Image< typename TDiffusionModelCalculator::DtType, TImageDimension > >
{

public:
  typedef DTIReconImageFilter<TDiffusionModelCalculator,TImageDimension >
                                                            Self;

  typedef ImageToImageFilter
    < VectorImage< typename TDiffusionModelCalculator::DWIPixelType::ComponentType, TImageDimension >,
    Image< typename TDiffusionModelCalculator::DtType, TImageDimension > >
                                                            Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

   /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DTIReconImageFilter, ImageToImageFilter);

  /** Typedefs for Pixel and Image Types. */
  typedef typename Superclass::OutputImageType              OutputImageType;
  typedef typename OutputImageType::PixelType               OutputPixelType;
  typedef typename OutputImageType::Pointer                 OutputImagePointer;

  typedef typename Superclass::InputImageType               InputImageType;
  typedef typename InputImageType::ConstPointer             InputImageConstPointer;
  typedef typename InputImageType::PixelType                InputPixelType;

  /** Typedefs from diffusion Model calculator */
  typedef TDiffusionModelCalculator                             DiffusionModelCalculatorType;
  typedef typename DiffusionModelCalculatorType::Pointer        DiffusionModelCalculatorPointer;
  typedef typename DiffusionModelCalculatorType::ConstPointer   DiffusionModelCalculatorConstPointer;

  /** Typedefs for computing residual  */
  typedef typename DiffusionModelCalculatorType::PrecisionType  ResidualPrecisionType;
  typedef VectorImage< ResidualPrecisionType, TImageDimension > ResidualImageType;
  typedef typename ResidualImageType::PixelType                 ResidualPixelType;

  /** Random Typedefs  */
  typedef typename Superclass::OutputImageRegionType        OutputImageRegionType;

  /**  Typedefs for the mask */
  typedef SpatialObject< TImageDimension >                  ImageMaskType;
  typedef typename ImageMaskType::Pointer                   ImageMaskPointer;
  typedef typename ImageMaskType::ConstPointer              ImageMaskConstPointer;

  /** Set/Get the image mask.
   *  Only Pixels Inside of the mask will be considered
   */
  itkSetObjectMacro( ImageMask, ImageMaskType );
  itkGetConstObjectMacro( ImageMask, ImageMaskType );

  /** Set/Get the Model Calculator. */
  itkSetObjectMacro( DiffusionModelCalculator, DiffusionModelCalculatorType );
  itkGetConstObjectMacro( DiffusionModelCalculator, DiffusionModelCalculatorType );

  /** Get/set compute residuals flag. */
  itkSetMacro( CalculateResidualImage, bool );
  itkGetConstMacro( CalculateResidualImage, bool );
  itkBooleanMacro( CalculateResidualImage );

  /** Get/set flag for forcing returned tensors to SPD
   *    Multiplies any negative eigenvalues by -1
   */
  itkSetMacro( ForceTensorsToSPD, bool );
  itkGetConstMacro( ForceTensorsToSPD, bool );
  itkBooleanMacro( ForceTensorsToSPD );

  /** Get the Residual image . */
  itkGetConstObjectMacro( ResidualImage, ResidualImageType );

protected:
  DTIReconImageFilter();
  ~DTIReconImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void BeforeThreadedGenerateData();
#if ITK_VERSION_MAJOR == 4
#else
  typedef int                ThreadIdType;
#endif

  void ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread, ThreadIdType);
  
private:
  typename DiffusionModelCalculatorType::Pointer  m_DiffusionModelCalculator;

  /** Image Mask */
  ImageMaskPointer                                m_ImageMask;

  /**Variables for computing residue */
  typename ResidualImageType::Pointer             m_ResidualImage;

  bool                                            m_CalculateResidualImage;
  bool                                            m_ForceTensorsToSPD;

};

}

#include "itkDTIReconImageFilter.hxx"

#endif
