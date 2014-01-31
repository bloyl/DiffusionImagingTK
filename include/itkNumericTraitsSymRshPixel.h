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


#ifndef __itkNumericTraitsSymRshPixel_h
#define __itkNumericTraitsSymRshPixel_h

#include <itkNumericTraits.h>

#include <itkSymRealSphericalHarmonicRep.h>

// This file is meant to define numeric traits for tensor pixels types in itk
///TODO this should be reformated to match the ITK 4 > numericTraits files...
namespace itk
{

template< typename T, unsigned int D >
class NumericTraits< SymRealSphericalHarmonicRep<T, D> >
{
public:
  typedef T                                                 ValueType;
  typedef SymRealSphericalHarmonicRep<typename NumericTraits<T>::PrintType,D>
                                                            PrintType;
  typedef SymRealSphericalHarmonicRep<typename NumericTraits<T>::AbsType,D>
                                                            AbsType;
  typedef SymRealSphericalHarmonicRep<typename NumericTraits<T>::AccumulateType,D>
                                                            AccumulateType;
  typedef SymRealSphericalHarmonicRep<typename NumericTraits<T>::RealType,D>
                                                            RealType;
  typedef typename NumericTraits<T>::RealType                        ScalarRealType;
  typedef SymRealSphericalHarmonicRep<typename NumericTraits<T>::FloatType,D>
                                                            FloatType;
  static SymRealSphericalHarmonicRep<T,D>  ITKCommon_EXPORT ZeroValue()
  {
    return SymRealSphericalHarmonicRep<T,D>( NumericTraits<T>::Zero );
  }
  static SymRealSphericalHarmonicRep<T,D>  ITKCommon_EXPORT OneValue()
  {
    return SymRealSphericalHarmonicRep<T,D>( NumericTraits<T>::One );
  }

  static const SymRealSphericalHarmonicRep<T,D> ITKCommon_EXPORT
                                                            Zero;
  static const SymRealSphericalHarmonicRep<T,D> ITKCommon_EXPORT
                                                            One;

  /** Fixed length vectors cannot be resized, so an exception will
   *  be thrown if the input size is not valid.  In this case, the
   *  only valid size is 6. If the size is valid the tensor will be
   *  filled with zeros. */
  static void SetLength(SymRealSphericalHarmonicRep<T,D> & m, const unsigned int s)
  {
    if ( s != SymRealSphericalHarmonicRep<T,D>::Dimension )
      {
      itkGenericExceptionMacro(<< "Cannot set the size of a SymRealSphericalHarmonic "
                               << "to " << s << " not the correct length for an order "<< D
                               << "RSH pixel");
      }
    m.Fill(NumericTraits< T >::Zero);
  }

  /** Return the size of the underlying array. */
  static unsigned int GetLength(const SymRealSphericalHarmonicRep<T,D> &)
  {
    return SymRealSphericalHarmonicRep<T,D>::Dimension;
  }


};

template< typename T, unsigned int D >
const SymRealSphericalHarmonicRep<T,D>
NumericTraits< SymRealSphericalHarmonicRep<T, D> >
::Zero = SymRealSphericalHarmonicRep<T,D>( NumericTraits<T>::Zero );

template< typename T, unsigned int D >
const SymRealSphericalHarmonicRep<T,D>
NumericTraits< SymRealSphericalHarmonicRep<T, D> >
::One  = SymRealSphericalHarmonicRep<T,D>( NumericTraits<T>::One );

} // end namespace itk


#endif // __itkNumericTraitsSymRshPixel_h  
