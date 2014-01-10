/**
* @file  itkTrackerDirectionPickerDwi.hxx
* @brief This class is an abstract class that serves as the base for classes used
* to choose 3 dimensional directions from diffusion model images (ie DTI or FOD images)
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkTrackerDirectionPickerDWI_hxx
#define __itkTrackerDirectionPickerDWI_hxx

#include <itkTrackerDirectionPickerDWI.h>

#include <itkVector.h>

#include "vnl/algo/vnl_svd.h"

//TODO remove vnl_det.j
#include "vnl/vnl_rank.h"

namespace itk
{

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::TrackerDirectionPickerDWI()
{
  m_NumGenerator = GeneratorType::New();
  m_NumGenerator->Initialize();

  m_ModelCalculator = DiffModCalculatorType::New();

  m_Method = DTI_STT;
  m_TensorLineF = 0.5;
  m_TensorLineG = 0.5;
  
  m_isInitialized = false;
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::UseDTIStreamLineTracking()
{ 
  this->SetMethod(DTI_STT);
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::UseDTITensorDeflectionTracking()
{ 
  this->SetMethod(DTI_TEND);
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::UseDTITensorLineTracking()
{ 
  this->SetMethod(DTI_TLINE);
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::UseCSAODFStreamLineTracking()
{ 
  this->SetMethod(CSAODF_SST);
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::UseODFStreamLineTracking()
{ 
  this->SetMethod(ODF_SST);
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
typename TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>::DirectionType
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PickStartingDirection(const PointType pt) const
{
  if (! m_isInitialized)
    itkExceptionMacro("itkTrackerDirectionPickerDWI is not initialized")

  DirectionType dir;
  switch (m_Method)
  {
    case DTI_STT:
    case DTI_TEND:
    case DTI_TLINE:
    {
      //space for eigenSysetm
      typename DtType::EigenValuesArrayType eigVals;
      typename DtType::EigenVectorsMatrixType eigVecs;
      
      InterpConstPointer interp = this->GetInterpolator();
      DtType dti = m_ModelCalculator->ComputeTensorWLS_residualBoot( interp->Evaluate(pt), false );
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
    case CSAODF_SST:
    case ODF_SST:
    {
      InterpConstPointer interp = this->GetInterpolator();
      
      RshType rsh;
      if ( m_Method == CSAODF_SST)
        rsh = m_ModelCalculator->ComputeRSH_CSAODF( interp->Evaluate(pt), 0.1, 0.1) ;
      else if ( m_Method == ODF_SST) 
        rsh = m_ModelCalculator->ComputeRSH_ODF( interp->Evaluate(pt) );

      PeakValueContainerType rshPeakVals;
      PeakDirectionContainerType peaks = m_PeakCalc->ComputePeaks(rsh, rshPeakVals);

      if (peaks.size() < 1)
      {
        std::cerr << rsh << std::endl;
        itkExceptionMacro("ERROR no peaks returned at point "<<pt);
      }

      double totalPeakVals = 0;
      for(typename PeakValueContainerType::const_iterator it=rshPeakVals.begin();
              it != rshPeakVals.end(); ++it)
      {
        totalPeakVals += *it;
      }

      if (totalPeakVals == 0)
      {
        std::cerr << rsh << std::endl;
        itkExceptionMacro("ERROR no peakValues returned at point "<<pt);
      }

      //Choose a randomn peak
      // Probability of each peak should be related to it intensity 
      double tmpVal = m_NumGenerator->GetUniformVariate(0,totalPeakVals);
      unsigned int peakIndex;
      double critcalVal=0;
      for (peakIndex = 0; peakIndex < rshPeakVals.size(); peakIndex++)
      {
        critcalVal +=  rshPeakVals[peakIndex];
        if (tmpVal < critcalVal)
          break;
      }

      dir[0] = peaks[peakIndex][0];
      dir[1] = peaks[peakIndex][1];
      dir[2] = peaks[peakIndex][2];

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

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
typename TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>::DirectionType
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PickNextDirection(PointType pt, DirectionType inDir) const
{

  if (! m_isInitialized)
    itkExceptionMacro("itkTrackerDirectionPickerDWI is not initialized")

  DirectionType dir;
  switch (m_Method)
  {
    case DTI_STT:
      dir = PickNextDirectionDTI_SST(pt,inDir);
      break;
    case DTI_TEND:
      dir = PickNextDirectionDTI_TEND(pt,inDir);
      break;
    case DTI_TLINE:
      dir = PickNextDirectionDTI_TLINE(pt,inDir);
      break;
    case CSAODF_SST:
    case ODF_SST:
      dir = PickNextDirectionRSH_SST(pt,inDir);
      break;
    default:
      itkExceptionMacro("Invalid method option : "<<m_Method)
  }
  return dir;
}


/**  Specific methods for picking peaks...*/

template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
typename TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>::DirectionType
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PickNextDirectionDTI_SST(PointType pt, DirectionType inDir) const
{
  InterpConstPointer interp = this->GetInterpolator();

  DtType dti = m_ModelCalculator->ComputeTensorWLS_residualBoot( interp->Evaluate(pt), false );
  
  //space for eigenSysetm
  typename DtType::EigenValuesArrayType eigVals;
  typename DtType::EigenVectorsMatrixType eigVecs;
  
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
template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
typename TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>::DirectionType
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PickNextDirectionDTI_TEND(PointType pt, DirectionType inDir) const
{
  InterpConstPointer interp = this->GetInterpolator();
  DtType dti = m_ModelCalculator->ComputeTensorWLS_residualBoot( interp->Evaluate(pt), false );

  //inDir is in LPS coords while dti is in (voxel) coords
  // so some measurement frame business is needed.
  typedef typename ImageType::DirectionType   DirCosMatType;
  DirCosMatType dirCos = interp->GetInputImage()->GetDirection();
  DirCosMatType dirCosInv = DirCosMatType( dirCos.GetInverse() );
  
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
template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
typename TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>::DirectionType
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PickNextDirectionDTI_TLINE(PointType pt, DirectionType inDir) const
{
  //TODO this can be optimized to cut down on calls to Interpolator
  // I'm not sure if this is correct but we renormalize each of the component vectors...
  DirectionType e1    = PickNextDirectionDTI_SST(pt, inDir);
  DirectionType eDefl = PickNextDirectionDTI_TEND(pt, inDir);
  DirectionType gProd = (1-m_TensorLineG)*inDir + m_TensorLineG * eDefl;
  gProd.Normalize();

  DirectionType dir = m_TensorLineF * e1 + (1 - m_TensorLineF) * ( gProd );
  dir.Normalize();
  return dir;
}

///////TODO THIS Method and the findPeak method needes to be figured out as far as casting direction types around.
template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
typename TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>::DirectionType
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::PickNextDirectionRSH_SST(PointType pt, DirectionType inDir) const
{
  InterpConstPointer interp = this->GetInterpolator();
  
  RshType rsh;
  if ( m_Method == CSAODF_SST)
    rsh = m_ModelCalculator->ComputeRSH_CSAODF( interp->Evaluate(pt), 0.1, 0.1) ;
  else if ( m_Method == ODF_SST) 
    rsh = m_ModelCalculator->ComputeRSH_ODF( interp->Evaluate(pt) );
  else
    itkExceptionMacro("Invalid method option : "<<m_Method)

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




template <class TInterpType, class TPrecisionType, unsigned int TOutputOrder>
void 
TrackerDirectionPickerDWI<TInterpType,TPrecisionType,TOutputOrder>
::Initialize()
{
  //Initialize both fittings...
  m_ModelCalculator->InitializeTensorFitting();
  m_ModelCalculator->InitializeRSHFitting();

  m_PeakCalc = PeakFinderType::New();
  m_PeakCalc->Initialize(3,2);

  m_isInitialized = true;
}

}// end namespace itk

#endif
