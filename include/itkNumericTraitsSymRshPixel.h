/**
 * @file  itkNumericTraitsSymRshPixel.h
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

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
