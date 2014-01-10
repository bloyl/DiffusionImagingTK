/**
 * @file  itkRSHPixelReorientationOperator.h
 * @brief reorient pixels of real spherical harmonic functions based supplied affine Transforms.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkRSHPixelReorientationOperator_h
#define __itkRSHPixelReorientationOperator_h

//#define USE_SPARSE_MATRIX

#include <itkObject.h>
#ifdef USE_SPARSE_MATRIX
#  include <vnl/vnl_sparse_matrix.h>
#else
#  include <itkMatrix.h>
#endif

#include <itkSymRealSphericalHarmonicRep.h>
#include <itkPixelReorientationOperator.h>
#include <itkReplaceSpecialFunctions.h>

namespace itk
{

/** @class PixelReorientationOperator
 * @brief The Purpose of this class is to reorient RSH pixeltypes based supplied affine Transforms.
 *
 * This template works for Symetric spherical harmonics coeffs (supplied by itkSymRealSphericalHarmonicRep.h)
 * Note that this class assumes the basis used in itkSymRealSphericalHarmonicRep
 * and will not yield the correct results for other basis sets
 *
 * The finite strain algorithm is used to compute the rotational
 * component (R) of the affine transform which is used to rotate the tensor.
 *
 * 2 Methods are supplied
 *
 * reorient computes    f' --- f'(g) = f( Transpose(R) * g )
 * reorientInv computes f' --- f'(g) = f( R * g )
 *
 * where f is the original function g is a unit vecotor.
 */
template <class TCompType, unsigned int TOrder>
class ITK_EXPORT PixelReorientationOperator
  <SymRealSphericalHarmonicRep<TCompType,TOrder> >:
    public Object
{
public:

  /** Standard class typedefs. */
  typedef SymRealSphericalHarmonicRep<TCompType,TOrder> PixelType;

  typedef PixelReorientationOperator<PixelType>         Self;

  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  static const int MaxOrder       = TOrder;
  static const int Dimensions     = (MaxOrder+1)*(MaxOrder+2)/2;
  static const int NumberOfOrders = 1+(MaxOrder/2);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PixelReorientationOperator, Object);

  /**
   * Reorient by the provided matrix
   * compute    f' --- f'(g) = f( Transpose(R) * g )
   **/
  template<typename TMatrixValueType>
  PixelType Reorient
    (const PixelType p, const Matrix<TMatrixValueType, 3u, 3u> & m) const;
  template<typename TMatrixValueType>
  PixelType Reorient
    (const PixelType p, const vnl_matrix_fixed<TMatrixValueType, 3u, 3u> & m) const
  {
    return this->Reorient(p, static_cast<Matrix<TMatrixValueType, 3u, 3u> >(m) );
  }
  template<typename TMatrixValueType>
  PixelType Reorient
    (const PixelType p, const vnl_matrix<TMatrixValueType> & m) const
  {
    return this->Reorient(p, static_cast<Matrix<TMatrixValueType> >(m) );
  }

  /**
   * Reorient by the Inverse of the provided matrix
   * reorientInv computes f' --- f'(g) = f( R * g )
   **/
  template<typename TMatrixValueType>
  PixelType ReorientInv
    (const PixelType p, const Matrix<TMatrixValueType, 3u, 3u> & m) const;
  template<typename TMatrixValueType>
  PixelType ReorientInv
    (const PixelType p, const vnl_matrix_fixed<TMatrixValueType, 3u, 3u> & m) const
  {
    return this->ReorientInv(p, static_cast<Matrix<TMatrixValueType, 3u, 3u> >(m) );
  }
  template<typename TMatrixValueType>
  PixelType ReorientInv
    (const PixelType p, const vnl_matrix<TMatrixValueType> & m) const
  {
    return this->ReorientInv(p, static_cast<Matrix<TMatrixValueType> >(m) );
  }

protected:
  //Typedefs for Internal rotaion matrix Types

#ifdef USE_SPARSE_MATRIX
  typedef vnl_sparse_matrix<double>                   RshRotationMatixType;
#else
  typedef Matrix<double, Dimensions, Dimensions>      RshRotationMatixType;
#endif

  PixelReorientationOperator();
  ~PixelReorientationOperator() {};

private:
  PixelReorientationOperator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  template<typename TMatrixValueType>
  const RshRotationMatixType ComputeRotationMatrix(const Matrix<TMatrixValueType, 3u, 3u> & m) const;

  static const RshRotationMatixType ComputeZRotationMatrix(double);

  /**
   * Computes the wigner litte d matrix element(in sp harms) for pi over 2
   */
  static double computeWignerLittleDPlus(int, int, int);

  /**Member variables*/
  // Rz(pi/2)  Ry(pi/2)
  RshRotationMatixType                                m_DMatPiOver2;
  //  Ry(-pi/2) Rz(-pi/2)
  RshRotationMatixType                                m_DMatMinusPiOver2;

  //Computes the rotational componant of an affineMatrix
  template <typename TMatrixValueType>
  Matrix<TMatrixValueType>
  RotationFromAffine(const Matrix<TMatrixValueType> f) const;


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRSHPixelReorientationOperator.hxx"
#endif

#endif
