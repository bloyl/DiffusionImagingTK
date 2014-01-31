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
 
#ifndef __itkSymRealSphericalHarmonicRep_hxx
#define __itkSymRealSphericalHarmonicRep_hxx

#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>

#include <itkSymRealSphericalHarmonicRep.h>
#include <itkRSHPixelReorientationOperator.h>

namespace itk
{

/**
 * Assignment Operator
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>&
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator= (const Self& r)
{
  BaseArray::operator=(r);
  return *this;
}


/**
 * Assignment Operator from a scalar constant
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>&
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator= (const ComponentType & r)
{
  BaseArray::operator=(&r);
  return *this;
}

/**
 * Assigment from a plain array
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>&
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator= (const ComponentArrayType r )
{
  BaseArray::operator=(r);
  return *this;
}


/**
 * Returns a temporary copy of a vector
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator+(const Self & r) const
{
  Self result;
  for( unsigned int i=0; i<Dimension; i++)
    {
    result[i] = (*this)[i] + r[i];
    }
  return result;
}

/**
 * Returns a temporary copy of a vector
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator-(const Self & r) const
{
  Self result;
  for( unsigned int i=0; i<Dimension; i++)
    {
    result[i] = (*this)[i] - r[i];
    }
  return result;
}

/**
 * Performs addition in place
 */
template<class T,unsigned int TMaxOrder>
const SymRealSphericalHarmonicRep<T,TMaxOrder> &
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator+=(const Self & r)
{
  for( unsigned int i=0; i<Dimension; i++)
    {
    (*this)[i] += r[i];
    }
  return *this;
}

/**
 * Performs subtraction in place
 */
template<class T,unsigned int TMaxOrder>
const SymRealSphericalHarmonicRep<T,TMaxOrder> &
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator-=(const Self & r)
{
  for( unsigned int i=0; i<Dimension; i++)
    {
    (*this)[i] -= r[i];
    }
  return *this;
}

/**
 * Performs multiplication by a scalar, in place
 */
template<class T,unsigned int TMaxOrder>
const SymRealSphericalHarmonicRep<T,TMaxOrder> &
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator*=(const RealValueType & r)
{
  for( unsigned int i=0; i<Dimension; i++)
    {
    (*this)[i] *= r;
    }
  return *this;
}

/**
 * Performs division by a scalar, in place
 */
template<class T,unsigned int TMaxOrder>
const SymRealSphericalHarmonicRep<T,TMaxOrder> &
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator/=(const RealValueType & r)
{
  for( unsigned int i=0; i<Dimension; i++)
    {
    (*this)[i] /= r;
    }
  return *this;
}

/**
 * Performs multiplication with a scalar
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator*(const RealValueType & r) const
{
  Self result;
  for( unsigned int i=0; i<Dimension; i++)
    {
    result[i] = (*this)[i] * r;
    }
  return result;
}


/**
 * Performs division by a scalar
 */
template<class T,unsigned int TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
SymRealSphericalHarmonicRep<T,TMaxOrder>
::operator/(const RealValueType & r) const
{
  Self result;
  for( unsigned int i=0; i<Dimension; i++)
    {
    result[i] = (*this)[i] / r;
    }
  return result;
}

/**
 * Print content to an ostream
 */
template<class T,unsigned int TMaxOrder>
std::ostream &
operator<<(std::ostream& os,const SymRealSphericalHarmonicRep<T,TMaxOrder> & c )
{
  for(unsigned int i=0; i<c.GetNumberOfComponents(); i++)
    {
    os <<  static_cast<typename NumericTraits<T>::PrintType>(c[i]) << "  ";
    }
  return os;
}

/**
 * Read content from an istream
 */
template<class T,unsigned int TMaxOrder>
std::istream &
operator>>(std::istream& is, SymRealSphericalHarmonicRep<T,TMaxOrder> & dt )
{
  for(unsigned int i=0; i < dt.GetNumberOfComponents(); i++)
    {
    is >> dt[i];
    }
  return is;
}

/**
 * Evaluate the RealSphericalHarmonicRep as a function of theta and phi.
 */
template<class T,unsigned int TMaxOrder>
const typename SymRealSphericalHarmonicRep<T,TMaxOrder>::RealValueType
SymRealSphericalHarmonicRep<T,TMaxOrder>
  ::Evaluate( RealValueType theta, RealValueType phi ) const
{
  RealValueType result = 0;
  for( unsigned int i=0; i<Dimension; i++)
  {
    result += (*this)[i] * Y(i+1,theta,phi);
  }
  return result;
}

template<class T,unsigned int TMaxOrder>
const typename SymRealSphericalHarmonicRep<T,TMaxOrder>::RealValueType
SymRealSphericalHarmonicRep<T,TMaxOrder>
  ::Evaluate(GradientDirectionType Gradient) const
{
  double theta = acos(Gradient[2]);
  double phi   = atan2(Gradient[1],Gradient[0]); // atan2(y,x) = atan(y/x);

  return (*this).Evaluate(theta,phi);
}

template<class T,unsigned int TMaxOrder>
unsigned int 
SymRealSphericalHarmonicRep<T,TMaxOrder>
::GetJ(int l,int m)
{
  return 1 + m + l * (l+1) / 2;
}

template<class T,unsigned int TMaxOrder>
const typename SymRealSphericalHarmonicRep<T,TMaxOrder>::LmVector
SymRealSphericalHarmonicRep<T,TMaxOrder>
::GetLM(unsigned int j)
{
  const int l = 2 * (int) ( ((1 + vcl_sqrt(8 * j - 7)) / 2) / 2);
  const int m = j - 1 - l * (l+1) / 2;
  LmVector retVal;

  retVal[0] = l;
  retVal[1] = m;

  return retVal;
}

template<class T,unsigned int TMaxOrder>
const 
typename SymRealSphericalHarmonicRep<T,TMaxOrder>::RshBasisMatrixType
SymRealSphericalHarmonicRep<T,TMaxOrder>
::ComputeRshBasis( const GradientDirectionContainerType *gradContainer ) 
{

  RshBasisMatrixType basis;
  std::vector<unsigned int> gradientind;
  for(GradientDirectionContainerType::ConstIterator gdcit = gradContainer->Begin();
      gdcit != gradContainer->End(); ++gdcit)
  {
    if(gdcit.Value().one_norm() > 0.0)
    {
      gradientind.push_back(gdcit.Index());
    }
  }

  unsigned int NumberOfCoefficients = Dimension;
  unsigned int numGrads = gradientind.size();

  basis.set_size( numGrads, NumberOfCoefficients );

  for (unsigned int m = 0; m < numGrads; m++)
  {
    /*** Grad directions relation to theta phi
    * x = sin(theta) * cos(phi)
    * y = sin(theta) * sin(phi)
    * z = cos(theta)
    */
    
    GradientDirectionType dir = gradContainer->ElementAt(gradientind[m]);
    dir = dir / dir.two_norm();
    
    double theta, phi;

    if ( dir[2] == 1) // z = 1
    {
      theta =0.0;
      phi   =0.0;
    }
    else
    {
      theta = acos(dir[2]);
      phi   = atan2(dir[1], dir[0]); // atan2(y,x) = atan(y/x);
    }

    for (unsigned int c = 0; c < NumberOfCoefficients; c++)
      basis[m][c]  = Y(c+1,theta,phi);

  }

  return basis;
}

/**
 * Compute the real spherical harmonics the traditional way.
 * Not sure how this stacks up agains other computation stratagies.
 */
template<class T,unsigned int TMaxOrder>
double
SymRealSphericalHarmonicRep<T,TMaxOrder>
::Y( int j, double theta, double phi ) 
{

  LmVector vec = GetLM(j);
  const int l = vec[0];
  const int m = vec[1];
  
  if( m == 0 ) /// Y_l^0
    return K(l,0) * LegendreP(l,m,vcl_cos(theta));
  else if( m < 0 ) /// sqrt2 re(y_l^m)
    return vnl_math::sqrt2 * K(l,m) * vcl_cos(m*phi) * LegendreP(l,m,vcl_cos(theta));
  else ///(m > 0) sqrt2 im(y_l^m)
    return vnl_math::sqrt2* K(l,m) * vcl_sin(m*phi) * LegendreP(l,m,vcl_cos(theta));
}

/**
 * Nomalization factor for the spherical harmonics...
 * vcl_sqrt( ( (2*l+1) * factorial(l-m) ) / ( 4*( vnl_math::pi ) * factorial(l+m) ) );
 * 
 * Use a speed up to compute factorial(l-m) / factorial(l+m)
 * 
 * No Overflow is checked but this should be more robust then using the
 * Factorial method.
 */
template<class T,unsigned int TMaxOrder>
double
SymRealSphericalHarmonicRep<T,TMaxOrder>
::K( int l, int m )
{
  double f = 1; //if m=0
  if (m > 0)
  {
    for(int i=l-m+1; i<l+m+1; i++)
    {
      f /= i; 
    }
  }
  else
  {
    for(int i=l+m+1; i<l-m+1; i++)
      {
        f *= i; 
      }
  }
  return vcl_sqrt( ( (2*l+1) / ( 4*( vnl_math::pi ) ) * f ) );
}

/**
 * Normalize the function to unit strength (integral over the sphere = 1)
 */  
template<class T,unsigned int TMaxOrder>
void
SymRealSphericalHarmonicRep<T,TMaxOrder>
::Normalize( ) 
{
  //only normalize if 0th order is non zero...
  if ( (*this)[0] != 0 &&
       (*this)[0] != 1.0 / ( 2 * vcl_sqrt( vnl_math::pi )) )
  {
    double nFactor = 1.0 / ( 2 * vcl_sqrt( vnl_math::pi ) * (*this)[0]);
    for( unsigned int i=0; i< Dimension; i++)
      (*this)[i] = static_cast<T>(nFactor * (*this)[i]);
  }
}

/**
 * Compute the Gaussian and Mean Curvatures as a function of theta and phi.
 */
template<class T,unsigned int TMaxOrder>
const typename SymRealSphericalHarmonicRep<T,TMaxOrder>::RealValueType
SymRealSphericalHarmonicRep<T,TMaxOrder>
  ::ComputeCurvatures(GradientDirectionType Gradient,
        RealValueType& gaussCurvature, RealValueType& meanCurvature, RealValueType& k1, RealValueType& k2 ) const
{
  RealValueType theta = acos(Gradient[2]);
  RealValueType phi   = atan2(Gradient[1],Gradient[0]); // atan2(y,x) = atan(y/x);

  return (*this).ComputeCurvatures(theta,phi,gaussCurvature,meanCurvature,k1,k2);
}

template<class T,unsigned int TMaxOrder>
const typename SymRealSphericalHarmonicRep<T,TMaxOrder>::RealValueType
SymRealSphericalHarmonicRep<T,TMaxOrder>
  ::ComputeCurvatures( RealValueType theta, RealValueType phi,
        RealValueType& gaussCurvature, RealValueType& meanCurvature, RealValueType& k1, RealValueType& k2 ) const
{

  //Deal with the singlularity at theta == 0
  if (theta < 1e-4 || theta > vnl_math::pi-1e-4) // This threshold was emprically choosen so the isotropic curvature computation would pass.
  {
    //Rotate this pixelType by pi/2 about the Y axis
    typedef PixelReorientationOperator<Self>   ReorienterType;
    typename ReorienterType::Pointer orienter  = ReorienterType::New();
    vnl_matrix<double> rotMat(3,3,0);
    rotMat(0,2) =  1;
    rotMat(1,1) =  1;
    rotMat(2,0) = -1;
    
    //Compute new theta and phi from rotMat * grad
    //TODO With some thought this should be computable directly!
    GradientDirectionType tmpGrad;
    tmpGrad[0] =  cos(theta);
    tmpGrad[1] =  sin(theta) * sin(phi);
    tmpGrad[2] = -sin(theta) * cos(phi);
    RealValueType theta2 = acos(tmpGrad[2]);
    RealValueType phi2   = atan2(tmpGrad[1],tmpGrad[0]); // atan2(y,x) = atan(y/x);

    Self tmp = orienter->ReorientInv(*this,rotMat);
    return tmp.ComputeCurvatures(theta2,phi2,gaussCurvature,meanCurvature,k1,k2);
  } 

  RealValueType result = 0;

  RealValueType lambda = this->Evaluate(theta, phi);
  const RealValueType E = lambda * lambda;
  const RealValueType G = lambda * lambda * vcl_sin(theta) * vcl_sin(theta);

  //Initialize 
  RealValueType e=lambda;                                   // lambda - Psi_{theta,theta}
  RealValueType g=lambda * vcl_sin(theta) * vcl_sin(theta); // lambda (sin(theta)^2) - Psi_{phi,phi}
  RealValueType f=0;                                        // - Psi_{theta,phi}

  // std::cout << theta<< "/" << phi << std::endl;
  // std::cout << e << "/" << E << "/" << g << "/" << G << std::endl;

  RealValueType z            =  vcl_cos(theta);
  RealValueType dz_dTheta    = -vcl_sin(theta);
  RealValueType d2z_dTheta2  = -vcl_cos(theta);

  //First and Second dirvatives of associated legendre polynomials.
  RealValueType diffLegFirst;
  RealValueType diffLegSecond;

  //Placeholders 
  RealValueType tmpFirst;
  RealValueType tmpSecond;

  RealValueType Plmz;
  RealValueType Plplus1mz;
  RealValueType Plplus2mz;

  //Skip i == 0 since it is constant and does not contribute to and partial derivatives.
  for( unsigned int i=0; i<Dimension; i++)
  {
    LmVector vec = GetLM(i+1);
    const int l = vec[0];
    const int m = vec[1];

    //Differentiate LegendreP(l,m,vcl_cos(theta)) with respeact to cos(theta)
    //Taken from LegendreP_derivatives.  LegendreP_derivatives( l, m, z, tmpFirst, tmpSecond );

    Plmz      = LegendreP(l,m,z);
    Plplus1mz = LegendreP(l+1,m,z);
    Plplus2mz = LegendreP(l+2,m,z);

    tmpFirst  = (-((1 + l)*z*Plmz) + (1 - m + l)*Plplus1mz)/(-1 + z*z);
    tmpSecond = ((1 + l + (1 + l)*(2 + l)*z*z)*Plmz + (-1 + m - l)*((5 + 2*l)*z*Plplus1mz + (-2 + m - l)*Plplus2mz))/(-1 + z*z)/(-1 + z*z);
    
    // d_Plm(cos(t))/dt = dcos(t)/dt * d_Plm(z)/dz
    diffLegFirst = dz_dTheta * tmpFirst;
    
    // // d2_Plm(cos(t))/dt2 = (dcos(t)/dt)^2 * d2_Plm(z)/dz2 + d2cos(t)/dt2 * d_Plm(z)/dz
    diffLegSecond = dz_dTheta * dz_dTheta * tmpSecond + d2z_dTheta2 * tmpFirst;
    
    if( m == 0 )      /// Y = K(l,0) * LegendreP(l,m,vcl_cos(theta));
    {
      e -= (*this)[i] * ( vnl_math::sqrt2 * K(l,m) * vcl_cos(m*phi)  * diffLegSecond );
      //No Phi dependence so don't need g and f
    }
    else if( m < 0 )  /// Y = vnl_math::sqrt2 * K(l,m) * vcl_cos(m*phi) * LegendreP(l,m,vcl_cos(theta));
    {
      e -= (*this)[i] * ( K(l,0) * diffLegSecond );
      g -= (*this)[i] * (-m * m * vnl_math::sqrt2 * K(l,m) * vcl_cos(m*phi) * Plmz );
      f -= (*this)[i] * (-m *     vnl_math::sqrt2 * K(l,m) * vcl_sin(m*phi) * diffLegFirst );
    }
    else              ///Y = vnl_math::sqrt2* K(l,m) * vcl_sin(m*phi) * LegendreP(l,m,vcl_cos(theta));
    {
      e -= (*this)[i] * ( vnl_math::sqrt2* K(l,m) * vcl_sin(m*phi) * diffLegSecond );
      g -= (*this)[i] * (-m * m * vnl_math::sqrt2 * K(l,m) * vcl_sin(m*phi) * Plmz );
      f -= (*this)[i] * ( m *     vnl_math::sqrt2 * K(l,m) * vcl_cos(m*phi) * diffLegFirst );
    }
  }

  //  K = (eg -f^2) / (E*G - F^2) *** NOTE THAT F = 0
  gaussCurvature = (e*g - f*f) / ( E*G );
  // H = 1/2 * ( eG - 2fF + gE) / (EG-F^2)  *** NOTE THAT F = 0
  meanCurvature  = 0.5 * ( e*G + g*E) / (E*G);

  //Fix for isotropic 
  double diff = meanCurvature * meanCurvature - gaussCurvature;
  if (diff < -1e-4)
    itkExceptionMacro("Error computing curvature. GuassCurvature > meanCurvature^2")
  else if (diff < 0)
    k1 = k2 = meanCurvature;
  else
  {
    k1 = meanCurvature + vcl_sqrt( meanCurvature * meanCurvature - gaussCurvature);
    k2 = meanCurvature - vcl_sqrt( meanCurvature * meanCurvature - gaussCurvature);    
  }
  return result;
}




} // end namespace itk
 
#endif
