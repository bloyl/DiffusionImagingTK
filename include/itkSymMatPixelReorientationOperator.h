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

#ifndef __itkSymMatPixelReorientationOperator_h
#define __itkSymMatPixelReorientationOperator_h

#include <itkObject.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkPixelReorientationOperator.h>

namespace itk
{

/** @class PixelReorientationOperator
 * @brief The Purpose of this class is to reorient symmetric tensors
 * based supplied affine Transforms.
 *
 * The finite strain algorithm is used to compute the rotational
 * component of the affine transform which is used to rotate the tensor.
 *
 * 2 Methods are supplied
 *
 * reorient computes    ---  R * D * Transpose(R)
 * reorientInv computes ---  Transpose(R) * D * R
 *
 */
template <class TCompType >
class ITK_EXPORT PixelReorientationOperator
  < SymmetricSecondRankTensor< TCompType, 3 > >:
    public Object
{
public:

  /** Standard class typedefs. */
  typedef PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >
                                                        Self;

  typedef SymmetricSecondRankTensor< TCompType, 3 >    PixelType;

  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PixelReorientationOperator, Object);

  /**
   * Reorient by the provided matrix
   *
   * this returns R * D * Transpose(R)
   **/
  template<typename TMatrixValueType>
  PixelType Reorient
    (PixelType p, const Matrix<TMatrixValueType, 3u, 3u> & m) const;

  template<typename TMatrixValueType>
  PixelType Reorient
    (PixelType p, const vnl_matrix_fixed<TMatrixValueType, 3u, 3u> & m) const
  {
    return this->Reorient(p, static_cast<Matrix<TMatrixValueType, 3u, 3u> >(m) );
  }

  template<typename TMatrixValueType>
  PixelType Reorient
    (PixelType p, const vnl_matrix<TMatrixValueType> & m) const
  {
    return this->Reorient(p, static_cast<Matrix<TMatrixValueType> >(m) );
  }

  /**
   * Reorient by the inverse of the provided matrix
   *   This is useful when the jacobian being passed in
   *   describes the fixed to moving transform and the moving image is being reoriented
   *
   *    this returns Transpose(R) * D * R
   **/
  template<typename TMatrixValueType>
  PixelType ReorientInv
    (PixelType p, const Matrix<TMatrixValueType, 3u, 3u> & m) const;

  template<typename TMatrixValueType>
  PixelType ReorientInv
    (PixelType p, const vnl_matrix_fixed<TMatrixValueType, 3u, 3u> & m) const
  {
    return this->ReorientInv(p, static_cast<Matrix<TMatrixValueType, 3u, 3u> >(m) );
  }

  template<typename TMatrixValueType>
  PixelType ReorientInv
    (PixelType p, const vnl_matrix<TMatrixValueType> & m) const
  {
    return this->ReorientInv(p, static_cast<Matrix<TMatrixValueType> >(m) );
  }

protected:
  PixelReorientationOperator();
  ~PixelReorientationOperator() {};

private:
  PixelReorientationOperator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //Computes the rotational componant of an affineMatrix
  template <typename TMatrixValueType>
  Matrix<TMatrixValueType>
  RotationFromAffine(const Matrix<TMatrixValueType> f) const;


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSymMatPixelReorientationOperator.hxx"
#endif

#endif
