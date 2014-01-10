/**
 * @file  itkPixelReorientationOperator.h
 * @brief The Purpose of this class is to reorient pixels based supplied affine Transforms.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPixelReorientationOperator_h
#define __itkPixelReorientationOperator_h

#include <itkProcessObject.h>
#include <itkMatrix.h>

namespace itk
{

/** @class PixelReorientationOperator
 * @brief The Purpose of this class is to reorient pixels based supplied affine Transforms.
 *
 * Since only specialized PixelTypes have directional/orientational information the default behavior
 * is not to perform rotation.
 *
 */
template <class TPixelType>
class ITK_EXPORT PixelReorientationOperator:
    public Object
{
public:

  /** Standard class typedefs. */
  typedef PixelReorientationOperator<TPixelType>
                                                        Self;
  typedef TPixelType                                    PixelType;


  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PixelReorientationOperator, Object);

  /**Reorient by the provided matrix
   * Default is to not perform Reorientation
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

  template<typename TMatrixValueType>
  PixelType ReorientInv
    (const PixelType p, const vnl_matrix_fixed<TMatrixValueType, 3u, 3u> & m) const
  {
    return this->Reorient(p, static_cast<Matrix<TMatrixValueType, 3u, 3u> >(m) );
  }

  template<typename TMatrixValueType>
  PixelType ReorientInv
    (const PixelType p, const vnl_matrix<TMatrixValueType> & m) const
  {
    return this->Reorient(p, static_cast<Matrix<TMatrixValueType> >(m) );
  }


protected:
  PixelReorientationOperator();
  ~PixelReorientationOperator() {};

private:
  PixelReorientationOperator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPixelReorientationOperator.hxx"
#endif

//Include specialized versions of the PixelOrientationOperator
#include <itkDTIPixelReorientationOperator.h>
#include <itkSymMatPixelReorientationOperator.h>
#include <itkRSHPixelReorientationOperator.h>

#endif
