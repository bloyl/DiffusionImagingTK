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
