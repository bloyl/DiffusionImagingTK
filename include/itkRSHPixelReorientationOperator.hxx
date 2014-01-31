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

#ifndef __itkRSHPixelReorientationOperator_hxx
#define __itkRSHPixelReorientationOperator_hxx

#include <vnl/vnl_det.h>

#include <itkRSHPixelReorientationOperator.h>

#define USE_WIGNER_ROTATION 2

namespace itk
{

/**
 * Initialize new instance
 */
template <class TCompType, unsigned int TOrder >
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::PixelReorientationOperator()
{
  RshRotationMatixType ryPlus,ryMinus;
#ifdef USE_SPARSE_MATRIX
  ryPlus.set_size(Dimensions,Dimensions);
  ryMinus.set_size(Dimensions,Dimensions);
#else
  ryPlus.Fill(0.0);
  ryMinus.Fill(0.0);
#endif


  /** First compute Ry(pi/2) and Ry(-pi/2) */
  int lMax = MaxOrder;
  int bigI = 0;

  for(int l=0; l<lMax+1; l+=2)
  {
    for(int mp=-l; mp<l+1; mp++)
    {
      for(int m=-l; m<l+1; m++)
      {
        //computeSomeSigns
        int sign_m =    ( vcl_abs(m)     % 2 == 1) ? -1 :1;     // (-1)^m
        int sign_mp =   ( vcl_abs(mp)    % 2 == 1) ? -1 :1;     // (-1)^mp
        int sign_mpm =  ( vcl_abs(mp+m)  % 2 == 1) ? -1 :1;     // (-1)^(m+mp)
        int sign_m1 =   ( vcl_abs((m+1)) % 2 == 1) ? -1 :1;     // (-1)^(m+1)

        //Matrix indices
        int mpi = bigI + mp + l;
        int mi = bigI + m + l;
        double plusVal;
        if (m == 0)
        {
          if (mp==0)
          {
            plusVal = computeWignerLittleDPlus(l, mp, m);
          }
          else if (mp<0)
          {
            plusVal = vnl_math::sqrt2 * sign_mp * computeWignerLittleDPlus(l, -mp, m);
          }
          else //mp>0
          {
            plusVal = 0.0;
          }
        }
        else if (m < 0)
        {
          if (mp==0)
          {
            plusVal =vnl_math::sqrt2 * sign_m * computeWignerLittleDPlus(l, mp, -m);
          }
          else if (mp<0)
          {
            plusVal = sign_mpm * computeWignerLittleDPlus(l, -mp, -m)
                              + sign_mp * computeWignerLittleDPlus(l, -mp, m);
          }
          else //mp>0
          {
            plusVal = 0.0;
          }
        }
        else //m >0
        {
        if (mp==0)
          {
            plusVal = 0.0;
          }
          else if (mp<0)
          {
            plusVal = 0.0;
          }
          else //mp>0
          {
            plusVal = computeWignerLittleDPlus(l, mp, m) + sign_m1 * computeWignerLittleDPlus(l,mp,-m);
          }
        }

        ryPlus(mpi,mi) = plusVal;

        //Determine (-1)^(m-mp)
        int sign = 1;
        if (vcl_abs(m-mp) % 2 == 1) sign = -1;
        ryMinus(mpi,mi) = sign * plusVal;

      }
    }
    bigI+= 2 * l+1;

  }

  const RshRotationMatixType zPlus       = ComputeZRotationMatrix(vnl_math::pi_over_2);
  const RshRotationMatixType zMinus      = ComputeZRotationMatrix(-vnl_math::pi_over_2);
  m_DMatPiOver2       = ryPlus * zPlus;
  m_DMatMinusPiOver2  = zMinus * ryMinus;

}


template <class TCompType, unsigned int TOrder >
template<typename TMatrixValueType>
typename PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >::PixelType
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::Reorient(const PixelType p, const Matrix<TMatrixValueType, 3u, 3u> & j) const
{
  typedef Matrix<TMatrixValueType, 3u, 3u>  RotMatrixType;
  PixelType result;

  const unsigned int Dimension =
        SymRealSphericalHarmonicRep<TCompType,TOrder>::Dimension;

  //Determine a rotation matrix from the jacobian
  const RotMatrixType m = this->RotationFromAffine(j);

  const RshRotationMatixType rotMat = ComputeRotationMatrix(m);

  if (rotMat.GetVnlMatrix().has_nans())
  {
    itkExceptionMacro(<<"RSH reorientaion matrix has NANS");
  }

  typedef typename NumericTraits<TCompType>::AccumulateType  AccumulateType;
  for(unsigned int r=0; r< Dimension; r++)
  {
    AccumulateType sum = NumericTraits<AccumulateType>::ZeroValue();
    for(unsigned int c=0; c<Dimension; c++)
    {
      sum += rotMat(r,c) * p[c];
    }
    result[r] = static_cast<TCompType>( sum );
  }

  return result;
}



template <class TCompType, unsigned int TOrder >
template<typename TMatrixValueType>
typename PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >::PixelType
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::ReorientInv(const PixelType p, const Matrix<TMatrixValueType, 3u, 3u> & j) const
{
  typedef Matrix<TMatrixValueType, 3u, 3u>  RotMatrixType;
  PixelType result;

  const unsigned int Dimension =
        SymRealSphericalHarmonicRep<TCompType,TOrder>::Dimension;

  //Determin a rotation matrix from the jacobian
  const RotMatrixType m = this->RotationFromAffine(j);

  const RshRotationMatixType rotMat = ComputeRotationMatrix(m);

  if (rotMat.GetVnlMatrix().has_nans())
  {
    itkExceptionMacro(<<"RSH reorientaion matrix has NANS");
  }

  typedef typename NumericTraits<TCompType>::AccumulateType  AccumulateType;
  for(unsigned int r=0; r< Dimension; r++)
  {
    AccumulateType sum = NumericTraits<AccumulateType>::ZeroValue();
    for(unsigned int c=0; c<Dimension; c++)
    {
      sum += rotMat(c,r) * p[c];
    }
    result[r] = static_cast<TCompType>( sum );
  }

  return result;
}

template <class TCompType, unsigned int TOrder >
const
typename PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::RshRotationMatixType
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::ComputeZRotationMatrix(double alpha)
{
  typedef typename RshRotationMatixType::ValueType            MatrixValueType;
  RshRotationMatixType results;

  if (alpha ==0)
  {
    results.SetIdentity();
  }
  else
  {
    unsigned blockIndex = 0;

    for (int l=0;l<=static_cast<int>(MaxOrder);l+=2)
    {
      for (int m = -l;m<=l;++m)
      {
        if (m == 0)
        {
          results(blockIndex+l,blockIndex+l) = static_cast<MatrixValueType>(1);
        }
        else
        {
          results(blockIndex+l+m,blockIndex+l+m) =
            static_cast<MatrixValueType>(cos(m*alpha));

          int sign = 1;
          if (vcl_abs(m+1) % 2 == 1) sign = -1;
          results(blockIndex+l+m,blockIndex+l-m) =
            static_cast<MatrixValueType>(sign * sin(m*alpha));

        }
      }
      //Incerment the block index.
      blockIndex += 2*l + 1;
    }
  }
  return results;
}

/****EDMONDS implementaiton */
template <class TCompType, unsigned int TOrder >
double
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::computeWignerLittleDPlus(int l, int mp, int m)
{
  //Determine (-1)^(l-mp)
  int sign = 1;
  if ( vcl_abs(l-mp) % 2 == 1) sign = -1;

  //Compute the normalization factor
  double f;

  f = 1; //if m==mp
  if (mp>m)
  {
    for(int i=l+m+1; i<=l+mp; i++)
    {
      f *= i;
    }
    for(int i=l-mp+1; i<=l-m; i++)
    {
      f /= i;
    }
  }
  else if (m>mp)
  {
    for(int i=l-m+1; i<=l-mp; i++)
    {
      f *= i;
    }
    for(int i=l+mp+1; i<=l+m; i++)
    {
      f /= i;
    }
  }

  //Compute bounds on the summation
  int sMin = ( 0 > -(mp+m) ) ? 0 : -(mp+m); // Max[0, -(mp + m)]
  int sMax = ( l-mp < l-m ) ? l-mp : l-m;   // Min[l - mp, l - m];

  double sum = 0;
  for(int s= sMin; s<= sMax; s++) //Not that s > 0 so no need for vcl_abs(s) in (-1)^s
  {
    if (s % 2 == 1) // (-1)^s == -1
    {
      sum -= binomialCoeff(l + m, l - mp - s)*binomialCoeff(l - m, s);
    }
    else // (-1)^s == +1;
    {
      sum += binomialCoeff(l + m, l - mp - s)*binomialCoeff(l - m, s);
    }
  }

  return sign * vcl_sqrt(f) / vcl_pow(2.0,l) * sum;
}

template <class TCompType, unsigned int TOrder >
template < typename TMatrixValueType >
const
typename PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::RshRotationMatixType
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::ComputeRotationMatrix( const Matrix<TMatrixValueType, 3u, 3u> & in ) const
{

  const typename Matrix<TMatrixValueType, 3u, 3u>::InternalMatrixType m = in.GetTranspose();

  //Convert rotation matrix m in zyz euler angles..
  double alpha;
  double beta;
  double gamma;

  if ( vcl_abs(m(2,2)-1.0) <= 1.0E-6 ) //beta is 0 so alpha and gamma are indistiguishable
  {
    beta = 0;
    alpha = 0;
    gamma = vcl_atan2(-m(0,1),m(0,0));
  }
  else if ( vcl_abs(m(2,2) - (-1.0)) <= 1.0E-6  ) //beta is Pi
  {
    alpha = 0;
    beta = vnl_math::pi;
    gamma = vcl_acos(-m(0,0));
  }
  else
  {
    beta  = vcl_acos(m(2,2));
    alpha = vcl_atan2(m(2,1),-m(2,0));
    gamma = vcl_atan2(m(1,2),m(0,2));
  }

  //TODO Should we check alpha, beta and gamma are in bounds?

  #if (USE_WIGNER_ROTATION==1)

  const RshRotationMatixType zRotGamma      = ComputeZRotationMatrix(gamma);

  RshRotationMatixType results = zRotGamma;
  if ( gamma == 0) results.SetIdentity();

  if ( beta != 0)
  {
    results = GetDMatMinusPiOver2()
      * ComputeZRotationMatrix(beta) * GetDMatPiOver2() * results;
  }

  if ( alpha != 0)
  {
    results = ComputeZRotationMatrix(alpha) * results;
  }

#elif (USE_WIGNER_ROTATION==2)
  typedef typename RshRotationMatixType::ValueType            MatrixValueType;

  RshRotationMatixType M3; // M3 = zRotAlpha * GetDMatMinusPiOver2() * zRotBeta
  RshRotationMatixType M1; // M1 = GetDMatPiOver2() * zRotGamma

  const RshRotationMatixType P2 = m_DMatPiOver2;//GetDMatPiOver2();
  const RshRotationMatixType P1 = m_DMatMinusPiOver2;//GetDMatMinusPiOver2();

  //Fill in M1 and M3
  unsigned blockIndex = 0;

  int ip1Sign; // (-1)^(i+1)
  int jp1Sign; // (-1)^(j+1)

  for (int l=0;l<=static_cast<int>(MaxOrder);l+=2)
  {
    ip1Sign = -1; // start at (-1)^(0+1)
    for (int i = 0; i< 2*l+1; ++i)
    {
      jp1Sign = -1; // (-1)^(j+1)
      for (int j = 0; j< 2*l+1; ++j)
      {
        if ( (i==l) )
        {
          if (j==l)
          {
            M1(blockIndex+i,blockIndex+j) =
                    static_cast<MatrixValueType>(P2(blockIndex+i,blockIndex+j));

            M3(blockIndex+i,blockIndex+j) =
                    static_cast<MatrixValueType>(P1(blockIndex+i,blockIndex+j));

          }
          else
          {
            M1(blockIndex+i,blockIndex+j) = static_cast<MatrixValueType>
              ( vcl_cos( (j-l) * gamma) * P2(blockIndex+i,blockIndex+j) +
                jp1Sign * vcl_sin( (l-j) * gamma) * P2(blockIndex+i,blockIndex+2*l-j)
              );

            M3(blockIndex+i,blockIndex+j) = static_cast<MatrixValueType>
              ( vcl_cos( (j-l) * beta) * P1(blockIndex+i,blockIndex+j) +
                jp1Sign * vcl_sin( (l-j) * beta) * P1(blockIndex+i,blockIndex+2*l-j)
              );

          }
        }
        else //i != l
        {
          M1(blockIndex+i,blockIndex+j) = static_cast<MatrixValueType>
              ( vcl_cos( (j-l) * gamma) * P2(blockIndex+i,blockIndex+j) +
                jp1Sign * vcl_sin( (l-j) * gamma) * P2(blockIndex+i,blockIndex+2*l-j)
              );

          M3(blockIndex+i,blockIndex+j) = static_cast<MatrixValueType>(
            vcl_cos( (i-l) * alpha) * ( vcl_cos( (j-l) * beta) * P1(blockIndex+i,blockIndex+j)
                     + jp1Sign * vcl_sin( (l-j) * beta) * P1(blockIndex+i,blockIndex+2*l-j) )
            + ip1Sign * vcl_sin( (i-l) * alpha) * ( vcl_cos( (j-l) * beta) * P1(blockIndex+2*l-i,blockIndex+j)
                     + jp1Sign * vcl_sin( (l-j) * beta) * P1(blockIndex+2*l-i,blockIndex+2*l-j) )
            );

        }
        jp1Sign *= -1;
      }
      ip1Sign *= -1;
    }
    //Incerment the block index.
    blockIndex += 2*l + 1;
  }

  RshRotationMatixType results = M3 * M1;

#endif

  return results;
}

//Computes the rotational componant of an affineMatrix
template <class TCompType, unsigned int TOrder >
template <typename TMatrixValueType>
Matrix<TMatrixValueType>
PixelReorientationOperator< SymRealSphericalHarmonicRep<TCompType,TOrder> >
::RotationFromAffine(Matrix<TMatrixValueType> f) const
{
  typedef vnl_svd< TMatrixValueType > svdType;
  svdType svd( f.GetVnlMatrix() );
  itkDebugMacro( "In Warper svd:\n"<< svd);
  Matrix<TMatrixValueType> retMat = svd.U() * svd.V().transpose();

  //Reorient retMat if it contains an inversion.
  //we can do this since rsh are antipodally symmetric.
  retMat /= vnl_det(retMat.GetVnlMatrix()) ;

  if (retMat.GetVnlMatrix().has_nans())
  {
    itkExceptionMacro(<<"RotationMacroFromAffine Matrix contains nans!");
  }

  return retMat;
}

} // end namespace itk

#endif
