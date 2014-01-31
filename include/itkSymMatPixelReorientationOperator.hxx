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
 
#ifndef __itkSymMatPixelReorientationOperator_hxx
#define __itkSymMatPixelReorientationOperator_hxx

#include <itkSymMatPixelReorientationOperator.h>

namespace itk
{

/**
 * Initialize new instance
 */
template <class TCompType >
PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >
::PixelReorientationOperator()
{

}

template <class TCompType >
template <typename TMatrixValueType>
typename PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >::PixelType
PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >
::Reorient(PixelType p, const Matrix<TMatrixValueType> & j) const
{
  //Check size of J... should be 3x3!

  PixelType result;
  const unsigned int dimension = PixelType::Dimension;

  typedef Matrix<TMatrixValueType, dimension, dimension>              MatrixType;

  //Determine a rotation matrix from the jacobian
  MatrixType rotMat = this->RotationFromAffine(j);

  MatrixType SCT = rotMat; //p * Transpose(rotMat)
  for(unsigned int r=0; r<dimension; r++)
  {
    for(unsigned int c=0; c<dimension; c++)
    {
      double sum = 0.0;
      for(unsigned int t=0; t<dimension; t++)
      {
        sum += p(r,t) * rotMat(c,t);
      }
      SCT(r,c) = sum;
    }
  }

  //result = rotMat * sct;
  for(unsigned int r=0; r<dimension; r++)
  {
    for(unsigned int c=0; c<dimension; c++)
    {
      double sum = 0.0;
      for(unsigned int t=0; t<dimension; t++)
      {
        sum += rotMat(r,t) * SCT(t,c);
      }
      (result)(r,c) = static_cast<TCompType>( sum );
    }
  }

  return result;
}


template <class TCompType >
template <typename TMatrixValueType>
typename PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >::PixelType
PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >
::ReorientInv(PixelType p, const Matrix<TMatrixValueType> & j) const
{

  //Check size of J... should be 3x3!

  PixelType result;
  const unsigned int dimension = PixelType::Dimension;

  typedef Matrix<TMatrixValueType, dimension, dimension>              MatrixType;

  //Determine a rotation matrix from the jacobian
  MatrixType rotMat = this->RotationFromAffine(j);

  MatrixType SCT = rotMat; //p * rotMat
  for(unsigned int r=0; r<dimension; r++)
  {
    for(unsigned int c=0; c<dimension; c++)
    {
      double sum = 0.0;
      for(unsigned int t=0; t<dimension; t++)
      {
        sum += p(r,t) * rotMat(t,c);
      }
      SCT(r,c) = sum;
    }
  }

  //result = transpose(rotMat) * sct;
  for(unsigned int r=0; r<dimension; r++)
  {
    for(unsigned int c=0; c<dimension; c++)
    {
      double sum = 0.0;
      for(unsigned int t=0; t<dimension; t++)
      {
        sum += rotMat(t,r) * SCT(t,c);
      }
      (result)(r,c) = static_cast<TCompType>( sum );
    }
  }

  return result;
}

//Computes the rotational componant of an affineMatrix
template <class TCompType >
template <typename TMatrixValueType>
Matrix<TMatrixValueType>
PixelReorientationOperator< SymmetricSecondRankTensor< TCompType, 3 > >
::RotationFromAffine(Matrix<TMatrixValueType> f) const
{
  typedef vnl_svd< TMatrixValueType > svdType;
  svdType svd( f.GetVnlMatrix() );
  itkDebugMacro( "In Warper svd:\n"<< svd);
  return svd.U() * svd.V().transpose();
}


} // end namespace itk

#endif
