/**
* @file  itkTrackerDirectionPicker.h
* @brief This class is an abstract class that serves as the base for classes used
* to choose 3 dimensional directions from diffusion model images (ie DTI or FOD images)
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkTrackerDirectionPickerRSH_hxx
#define __itkTrackerDirectionPickerRSH_hxx

#include <itkTrackerDirectionPickerRSH.h>

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
TrackerDirectionPickerRSH<TInterpType>
::TrackerDirectionPickerRSH()
{
  m_Method = STT;
  m_NumGenerator = GeneratorType::New();
  m_NumGenerator->Initialize();
  m_PeakCalc = PeakFinderType::New();
  m_PeakCalc->Initialize(3,2);

}

template <class TInterpType >
void
TrackerDirectionPickerRSH<TInterpType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TInterpType >
void
TrackerDirectionPickerRSH<TInterpType>
::UseStreamLineTracking()
{ 
  this->SetMethod(STT);
}

template <class TInterpType >
void
TrackerDirectionPickerRSH<TInterpType>
::UseProbabilisticTracking()
{ 
  this->SetMethod(PROB);
}

template <class TInterpType >
typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
TrackerDirectionPickerRSH<TInterpType>
::PickStartingDirection(const PointType pt) const
{
  DirectionType dir;
  switch (m_Method)
  {
    case PROB:
//      dir = PickStartingDirectionProb(pt);
      break;
    case STT:
    {
      InterpConstPointer interp = this->GetInterpolator();
      
      PixelType rsh = interp->Evaluate(pt);
      PeakDirectionContainerType rawPeaks = m_PeakCalc->ComputePeaks(rsh);

      //Choose a randomn peak
      int peakIndex = m_NumGenerator->GetIntegerVariate(rawPeaks.size()-1);

      dir[0] = rawPeaks[peakIndex][0];
      dir[1] = rawPeaks[peakIndex][1];
      dir[2] = rawPeaks[peakIndex][2];

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
typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
TrackerDirectionPickerRSH<TInterpType>
::PickNextDirection(PointType pt, DirectionType inDir) const
{
  DirectionType dir;
  switch (m_Method)
  {
    case STT:
      dir = PickNextDirectionSST(pt,inDir);
      break;
    case PROB:
//      dir = PickNextDirectionProb(pt,inDir);
      break;
  }
  // std::cerr << inDir <<  " : " << dir << std::endl;
  return dir;
}

/**  Specific methods for picking peaks...*/

// template <class TInterpType >
// typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
// TrackerDirectionPickerRSH<TInterpType>
// ::PickNextDirectionSST(PointType pt, DirectionType inDir) const
// {
//   InterpConstPointer interp = this->GetInterpolator();
  
//   PixelType rsh = interp->Evaluate(pt);
//   PeakDirectionContainerType rawPeaks = m_PeakCalc->ComputePeaks(rsh);

//   //map inDir into model Frame.
//   typedef typename ImageType::DirectionType   DirCosMatType;
//   DirCosMatType dirCos = interp->GetInputImage()->GetDirection();
//   DirCosMatType dirCosInv = DirCosMatType( dirCos.GetInverse() );
//   DirectionType inDirMeas = dirCosInv * inDir;

//   //find the best Peak...
//   DirectionType dir;
//   double maxDot = -1;

//   typename PeakDirectionContainerType::const_iterator iter;
//   for (iter = rawPeaks.begin(); iter != rawPeaks.end(); ++iter )
//   {
//     double dot = vcl_abs( (*iter)[0]*inDirMeas[0] + (*iter)[1]*inDirMeas[1] + (*iter)[2]*inDirMeas[2] );
//     if ( dot > maxDot )
//     {
//       maxDot = dot;
//       dir[0] = (*iter)[0];
//       dir[1] = (*iter)[1];
//       dir[2] = (*iter)[2];
//     }
//   }

//   //Apply the measurementFrame
//   dir = dirCos * dir;
 
//   double dot = dir[0] * inDir[0] + dir[1] * inDir[1] + dir[2] * inDir[2];
//   int signDot = (dot >= 0) - (dot < 0); 

//   //dir = signOut * dir;
//   dir[0] = signDot * dir[0];
//   dir[1] = signDot * dir[1];
//   dir[2] = signDot * dir[2];
  
//   return dir;
// }

template <class TInterpType >
typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
TrackerDirectionPickerRSH<TInterpType>
::PickNextDirectionSST(PointType pt, DirectionType inDir) const
{
  InterpConstPointer interp = this->GetInterpolator();
  
  PixelType rsh = interp->Evaluate(pt);
  
  //map inDir into model Frame.
  typedef typename ImageType::DirectionType   DirCosMatType;
  DirCosMatType dirCos = interp->GetInputImage()->GetDirection();
  DirCosMatType dirCosInv = DirCosMatType( dirCos.GetInverse() );
  DirectionType tmp = dirCosInv * inDir;

  PeakDirectionType inDirMeas;
  inDirMeas[0] = tmp[0]; inDirMeas[1] = tmp[1]; inDirMeas[2] = tmp[2];

  //find the best Peak...
  PeakDirectionType peakDir = m_PeakCalc->FindPeak(rsh, inDirMeas);

  DirectionType dir;
  dir[0] = peakDir[0]; 
  dir[1] = peakDir[1]; 
  dir[2] = peakDir[2]; 
  //Apply the measurementFrame
  dir = dirCos * dir;

  return dir;
}

// Probbaliisitic methods and helper functions
// template <class TInterpType >
// typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
// TrackerDirectionPickerRSH<TInterpType>
// ::PickStartingDirectionProb(const PointType pt) const
// { 
//   DirectionType dir = PickRandomFromODF(pt);

//   // //Apply the measurementFrame
//   // InterpConstPointer interp = this->GetInterpolator();
//   // dir = interp->GetInputImage()->GetDirection() * dir;

//   return dir;
// }

// template <class TInterpType >
// typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
// TrackerDirectionPickerRSH<TInterpType>
// ::PickNextDirectionProb(PointType pt, DirectionType inDir) const
// {
  // DirectionType dir = PickRandomFromODF(pt);

  // //Apply the measurementFrame
  // InterpConstPointer interp = this->GetInterpolator();
  // dir = interp->GetInputImage()->GetDirection() * dir;
 
  // double dot = dir[0] * inDir[0] + dir[1] * inDir[1] + dir[2] * inDir[2];
  // int signDot = (dot >= 0) - (dot < 0); 

  // dir[0] = signDot * dir[0];
  // dir[1] = signDot * dir[1];
  // dir[2] = signDot * dir[2];
  // DirectionType dir;
  // return dir;
// }

// template <class TInterpType >
// typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
// TrackerDirectionPickerRSH<TInterpType>
// ::PickRandomDirectionOnSphere() const
// {
//   //Pick Random angles
//   double theta; //angle off of z-axis in [0, pi]
//   double phi;   //angle off of x-axis in [0,2 pi]

//   //Use method from http://mathworld.wolfram.com/SpherePointPicking.html
//   double u = m_NumGenerator->GetUniformVariate(0, 1);
//   while (u == 0 && u == 1 )
//   {
//     u = m_NumGenerator->GetUniformVariate(0, 1);
//   }
//   double v = m_NumGenerator->GetUniformVariate(0, 1);
//   while (v == 0 && v == 1 )
//   {
//     v = m_NumGenerator->GetUniformVariate(0, 1);
//   }
  
//   theta = vcl_acos(2*v -1);
//   phi   = 2 * vnl_math::pi * m_NumGenerator->GetUniformVariate(0, 1);
  
//   DirectionType dir;
  
//   dir[0] = sin(theta) * cos(phi); 
//   dir[1] = sin(theta) * sin(phi);
//   dir[2] = cos(theta);

//   return dir;
// }

// template <class TInterpType >
// typename TrackerDirectionPickerRSH<TInterpType>::DirectionType
// TrackerDirectionPickerRSH<TInterpType>
// ::PickRandomFromODF(const PointType pt) const
// {

//   //This all takes place within the internal frame of the tensor...
//   InterpConstPointer interp = this->GetInterpolator();
//   PixelType dti = interp->Evaluate(pt);

//   typename PixelType::EigenValuesArrayType eigVals;
//   typename PixelType::EigenVectorsMatrixType eigVecs;
  
//   dti.ComputeEigenAnalysis(eigVals,eigVecs);
  
//   DirectionType peakDir;
//   peakDir[0] = eigVecs(2,0);
//   peakDir[1] = eigVecs(2,1);
//   peakDir[2] = eigVecs(2,2);

//   double us = 1;
//   const double max_psi = EvaluateOdf(dti,peakDir);
//   double xs = 0;
//   DirectionType dir;
//   while (us > (xs / max_psi))
//   {
//     dir = PickRandomDirectionOnSphere();
//     us = m_NumGenerator->GetUniformVariate(0,1);
//     xs = EvaluateOdf(dti,dir);
//   }
//   return dir;
// }

// template <class TInterpType >
// double
// TrackerDirectionPickerRSH<TInterpType>
// ::EvaluateOdf(const PixelType dti, const DirectionType dir) const
// {
//   vnl_matrix_fixed<double,3,3> D;
//   D(0,0)          = dti[0];
//   D(0,1) = D(1,0) = dti[1];
//   D(0,2) = D(2,0) = dti[2];
//   D(1,1)          = dti[3];
//   D(1,2) = D(2,1) = dti[4];
//   D(2,2)          = dti[5];

//   vnl_matrix_fixed<double,3,3> Dinv = vnl_inverse(D);
//   const vnl_vector<double> g = dir.GetVnlVector();

//   // val = 1 / (4 pi) * 1 / (sqrt(det(D)) * 1 / (g' D^(-1) * g)^(3/2)
//   double val;
//   val = 1. / 4. / vnl_math::pi;
//   val /= vcl_sqrt( vnl_det(D) );
//   val /= vcl_pow( inner_product(g,Dinv*g), 1.5);

// }

}// end namespace itk

#endif
