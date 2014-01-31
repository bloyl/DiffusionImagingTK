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


#ifndef __itkVariableLengthVectorCastImageFilter_h
#define __itkVariableLengthVectorCastImageFilter_h

#include <itkUnaryFunctorImageFilter.h>
#include <itkNumericTraitsFixedArrayPixel.h>

namespace itk
{
  
/** \class VariableLengthVectorCastImageFilter
 *
 * @brief Casts input vector pixels to output vector pixel type.
 *
 * This filter needn't exist nor should the VectorCastImageFilter.
 * 
 * The CastImageFilter should work fine provided that vector classes
 * provide templated copy constructors to the other types of vector
 * classes.
 * 
 * This assumes we are converting from a vector to a variable length vector.
 * 
 * This filter is templated over the input image type and 
 * output image type.
 * 
 * The filter expect both images to have the same number of dimensions,
 * and that both the input and output have itk::Vector pixel types
 * of the same VectorDimension.
 *
 * \sa Vector
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  
  
template< class TInput, class TOutput>
class VariableLengthVectorCast
{
public:
  VariableLengthVectorCast() {}
  ~VariableLengthVectorCast() {}
  bool operator!=( const VariableLengthVectorCast & ) const
    {
    return false;
    }
  bool operator==( const VariableLengthVectorCast & other ) const
    {
    return !(*this != other);
    }
  inline TOutput operator()( const TInput & A ) const
    {
    typedef typename TOutput::ValueType OutputValueType;

    TOutput value(TInput::Dimension);
    for( unsigned int k = 0; k < TInput::Dimension; k++ )
      {
      value[k] = static_cast<OutputValueType>( A[k] );
//      std::cerr << k <<  " : " << value[k] << " : " << A[k] << std::endl;
      }
    return value;
    }
}; 
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT VariableLengthVectorCastImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::VariableLengthVectorCast< typename TInputImage::PixelType, 
                                             typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef VariableLengthVectorCastImageFilter                               Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Functor::VariableLengthVectorCast< typename TInputImage::PixelType, 
                         typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(VariableLengthVectorCastImageFilter, 
               UnaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TInputImage::PixelType::ValueType>));
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TOutputImage::PixelType::ValueType>));
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<typename TInputImage::PixelType::ValueType,
                          typename TOutputImage::PixelType::ValueType>));
  /** End concept checking */
#endif

protected:
  VariableLengthVectorCastImageFilter() {}
  virtual ~VariableLengthVectorCastImageFilter() {}
  void GenerateOutputInformation()
  {
    Superclass::GenerateOutputInformation();
    // get pointers to the input and output
    typename Superclass::OutputImagePointer      outputPtr = this->GetOutput();

    // propagate vector length info
    outputPtr->SetVectorLength(Superclass::InputImageType::PixelType::Dimension);

  }

private:
  VariableLengthVectorCastImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
