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

#ifndef __itkTrackerDirectionPickerDTI_hxx
#define __itkTrackerDirectionPickerDTI_hxx

#include <itkTrackerDirectionPickerDTI.h>

#include <itkVector.h>

#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_det.h>

namespace itk
{

template <class TInterpType >
TrackerDirectionPickerDTI<TInterpType>
::TrackerDirectionPickerDTI()
{
  m_Method = STT;
  m_NumGenerator = GeneratorType::New();
  m_NumGenerator->Initialize();
  m_TensorLineF = 0.5;
  m_TensorLineG = 0.5;

}

template <class TInterpType >
void
TrackerDirectionPickerDTI<TInterpType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TInterpType >
void
TrackerDirectionPickerDTI<TInterpType>
::UseStreamLineTracking()
{ 
  this->SetMethod(STT);
}

template <class TInterpType >
void
TrackerDirectionPickerDTI<TInterpType>
::UseTensorDeflectionTracking()
{ 
  this->SetMethod(TEND);
}

template <class TInterpType >
void
TrackerDirectionPickerDTI<TInterpType>
::UseTensorLineTracking()
{ 
  this->SetMethod(TLINE);
}

template <class TInterpType >
void
TrackerDirectionPickerDTI<TInterpType>
::UseProbabilisticTracking()
{ 
  this->SetMethod(PROB);
}

template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickStartingDirection(const PointType pt) const
{
  DirectionType dir;
  switch (m_Method)
  {
    case PROB:
      dir = PickStartingDirectionProb(pt);
      break;
    case STT:
    case TEND:
    case TLINE:
    {
      InterpConstPointer interp = this->GetInterpolator();
      //space for eigenSysetm
      typename PixelType::EigenValuesArrayType eigVals;
      typename PixelType::EigenVectorsMatrixType eigVecs;
      
      PixelType dti = interp->Evaluate(pt);
      dti.ComputeEigenAnalysis(eigVals,eigVecs);
      
      dir[0] = eigVecs(2,0);
      dir[1] = eigVecs(2,1);
      dir[2] = eigVecs(2,2);
      
      //Pick a random sign.
      double dot = m_NumGenerator->GetUniformVariate(-0.5, 0.5);
      int signDot = (dot >= 0) - (dot < 0); 
      dir[0] = signDot * dir[0];
      dir[1] = signDot * dir[1];
      dir[2] = signDot * dir[2];

      //Apply the measurementFrame
      dir = interp->GetInputImage()->GetDirection() * dir;
      break;
    }
    default:
      itkExceptionMacro("Invalid method option : "<<m_Method)
  }
  return dir;
}

template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickNextDirection(PointType pt, DirectionType inDir) const
{
  DirectionType dir;
  switch (m_Method)
  {
    case STT:
      dir = PickNextDirectionSST(pt,inDir);
      break;
    case TEND:
      dir = PickNextDirectionTEND(pt,inDir);
      break;
    case TLINE:
      dir = PickNextDirectionTLINE(pt,inDir);
      break;
    case PROB:
      dir = PickNextDirectionProb(pt,inDir);
      break;
  }
  // std::cerr << inDir <<  " : " << dir << std::endl;
  return dir;
}

/**  Specific methods for picking peaks...*/

template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickNextDirectionSST(PointType pt, DirectionType inDir) const
{
  InterpConstPointer interp = this->GetInterpolator();
  
  //space for eigenSysetm
  typename PixelType::EigenValuesArrayType eigVals;
  typename PixelType::EigenVectorsMatrixType eigVecs;
  
  PixelType dti = interp->Evaluate(pt);
  dti.ComputeEigenAnalysis(eigVals,eigVecs);

  DirectionType dir;
  dir[0] = eigVecs(2,0);
  dir[1] = eigVecs(2,1);
  dir[2] = eigVecs(2,2);

  //Apply the measurementFrame
  dir = interp->GetInputImage()->GetDirection() * dir;
 
  double dot = dir[0] * inDir[0] + dir[1] * inDir[1] + dir[2] * inDir[2];
  int signDot = (dot >= 0) - (dot < 0); 

  //dir = signOut * dir;
  dir[0] = signDot * dir[0];
  dir[1] = signDot * dir[1];
  dir[2] = signDot * dir[2];
  
  return dir;
}

/** Returns the outgoing direction based on the diffusion model at 
 *  the supplied point and the incoming direction
 *  Uses the Tensor deflection algorithm ( v_{out} = D * v_{in} )
 */ 
template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickNextDirectionTEND(PointType pt, DirectionType inDir) const
{
  InterpConstPointer interp = this->GetInterpolator();

  //inDir is in LPS coords while dti is in (voxel) coords
  // so some measurement frame business is needed.
  typedef typename ImageType::DirectionType   DirCosMatType;
  DirCosMatType dirCos = interp->GetInputImage()->GetDirection();
  DirCosMatType dirCosInv = DirCosMatType( dirCos.GetInverse() );
  
  PixelType dti = interp->Evaluate(pt);
  
  DirCosMatType tmpMat = dirCos * dti.PostMultiply(dirCosInv);
  
  DirectionType dir = tmpMat * inDir;

  dir.Normalize();
  return dir;
}

/** Returns the outgoing direction based on the diffusion model at 
 *  the supplied point and the incoming direction
 *  Uses the Tensor stream line algorithm
 *    required scalar parameters f and g to be set!
 *    v_{out} = f e1 + ͑ (1 - f͒͑͑) ( (1-g) ͒v_{in} + g D v_{in})
 */ 
template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickNextDirectionTLINE(PointType pt, DirectionType inDir) const
{
  //TODO this can be optimized to cut down on calls to Interpolator
  // I'm not sure if this is correct but we renormalize each of the component vectors...
  DirectionType e1    = PickNextDirectionSST(pt, inDir);
  DirectionType eDefl = PickNextDirectionTEND(pt, inDir);
  DirectionType gProd = (1-m_TensorLineG)*inDir + m_TensorLineG * eDefl;
  gProd.Normalize();

  DirectionType dir = m_TensorLineF * e1 + (1 - m_TensorLineF) * ( gProd );
  dir.Normalize();
  return dir;
}


// Probbaliisitic methods and helper functions
template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickStartingDirectionProb(const PointType pt) const
{ 
  DirectionType dir = PickRandomFromODF(pt);

  //Apply the measurementFrame
  InterpConstPointer interp = this->GetInterpolator();
  dir = interp->GetInputImage()->GetDirection() * dir;

  return dir;
}

template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickNextDirectionProb(PointType pt, DirectionType inDir) const
{
  DirectionType dir = PickRandomFromODF(pt);

  //Apply the measurementFrame
  InterpConstPointer interp = this->GetInterpolator();
  dir = interp->GetInputImage()->GetDirection() * dir;
 
  double dot = dir[0] * inDir[0] + dir[1] * inDir[1] + dir[2] * inDir[2];
  int signDot = (dot >= 0) - (dot < 0); 

  dir[0] = signDot * dir[0];
  dir[1] = signDot * dir[1];
  dir[2] = signDot * dir[2];

  return dir;
}

template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickRandomDirectionOnSphere() const
{
  //Pick Random angles
  double theta; //angle off of z-axis in [0, pi]
  double phi;   //angle off of x-axis in [0,2 pi]

  //Use method from http://mathworld.wolfram.com/SpherePointPicking.html
  double u = m_NumGenerator->GetUniformVariate(0, 1);
  while (u == 0 && u == 1 )
  {
    u = m_NumGenerator->GetUniformVariate(0, 1);
  }
  double v = m_NumGenerator->GetUniformVariate(0, 1);
  while (v == 0 && v == 1 )
  {
    v = m_NumGenerator->GetUniformVariate(0, 1);
  }
  
  theta = vcl_acos(2*v -1);
  phi   = 2 * vnl_math::pi * m_NumGenerator->GetUniformVariate(0, 1);
  
  DirectionType dir;
  
  dir[0] = sin(theta) * cos(phi); 
  dir[1] = sin(theta) * sin(phi);
  dir[2] = cos(theta);

  return dir;
}

template <class TInterpType >
typename TrackerDirectionPickerDTI<TInterpType>::DirectionType
TrackerDirectionPickerDTI<TInterpType>
::PickRandomFromODF(const PointType pt) const
{

  //This all takes place within the internal frame of the tensor...
  InterpConstPointer interp = this->GetInterpolator();
  PixelType dti = interp->Evaluate(pt);

  typename PixelType::EigenValuesArrayType eigVals;
  typename PixelType::EigenVectorsMatrixType eigVecs;
  
  dti.ComputeEigenAnalysis(eigVals,eigVecs);
  
  DirectionType peakDir;
  peakDir[0] = eigVecs(2,0);
  peakDir[1] = eigVecs(2,1);
  peakDir[2] = eigVecs(2,2);

  double us = 1;
  const double max_psi = EvaluateOdf(dti,peakDir);
  double xs = 0;
  DirectionType dir;
  while (us > (xs / max_psi))
  {
    dir = PickRandomDirectionOnSphere();
    us = m_NumGenerator->GetUniformVariate(0,1);
    xs = EvaluateOdf(dti,dir);
  }
  return dir;
}

template <class TInterpType >
double
TrackerDirectionPickerDTI<TInterpType>
::EvaluateOdf(const PixelType dti, const DirectionType dir) const
{
  vnl_matrix_fixed<double,3,3> D;
  D(0,0)          = dti[0];
  D(0,1) = D(1,0) = dti[1];
  D(0,2) = D(2,0) = dti[2];
  D(1,1)          = dti[3];
  D(1,2) = D(2,1) = dti[4];
  D(2,2)          = dti[5];

  vnl_matrix_fixed<double,3,3> Dinv = vnl_inverse(D);
  const vnl_vector<double> g = dir.GetVnlVector();

  // val = 1 / (4 pi) * 1 / (sqrt(det(D)) * 1 / (g' D^(-1) * g)^(3/2)
  double val;
  val = 1. / 4. / vnl_math::pi;
  val /= vcl_sqrt( vnl_det(D) );
  val /= vcl_pow( inner_product(g,Dinv*g), 1.5);

  return val;
}

}// end namespace itk

#endif
