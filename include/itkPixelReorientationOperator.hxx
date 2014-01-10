/**
 * @file  itkPixelReorientationOperator.hxx
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPixelReorientationOperator_hxx
#define __itkPixelReorientationOperator_hxx

#include <itkPixelReorientationOperator.h>

namespace itk
{

/**
 * Initialize new instance
 */
template <class TPixelType >
PixelReorientationOperator< TPixelType >
::PixelReorientationOperator()
{

}

template <class TPixelType >
template <typename TMatrixValueType>
typename PixelReorientationOperator< TPixelType >::PixelType
PixelReorientationOperator< TPixelType >
::Reorient(const PixelType p, const Matrix<TMatrixValueType> &) const
{
  return p;
}


} // end namespace itk

#endif
