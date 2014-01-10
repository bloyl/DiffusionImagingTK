/**
 * @file  itkPeakFindingCalculatorGrid.hxx
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPeakFindingCalculatorGridGrid_hxx
#define __itkPeakFindingCalculatorGridGrid_hxx

#include <itkPeakFindingCalculatorGrid.h>

namespace itk
{

template < typename TPixelType >
PeakFindingCalculatorGrid< TPixelType >
::PeakFindingCalculatorGrid()
{
  m_PeakFinder = PeakFindingCalculatorType::New();
  m_PixelReorienter = ReOrienterType::New();
}

template < typename TPixelType >
PeakFindingCalculatorGrid< TPixelType >
::~PeakFindingCalculatorGrid()
{
  m_ListOfNeighbors.clear();
  m_SampleDirections.clear();
}

template < typename TPixelType >
void
PeakFindingCalculatorGrid< TPixelType >
::Initialize(unsigned int SphereResolution, unsigned int Rigidity )
{

  // declare triangle mesh source
  typedef itk::RegularSphereMeshSource<TriangleMeshType>  SphereMeshSourceType;
  typedef SphereMeshSourceType::PointType PointType;
  typedef SphereMeshSourceType::VectorType VectorType;

  // Declare the type of the gradient image
  typedef itk::TriangleMeshToSimplexMeshFilter<TriangleMeshType, SimplexMeshType>  SimplexFilterType;

  SphereMeshSourceType::Pointer  mySphereMeshSource = SphereMeshSourceType::New();
  PointType center; center.Fill(0);
  PointType::ValueType scaleInit[3] = {1,1,1};
  VectorType scale = scaleInit;

  mySphereMeshSource->SetCenter(center);
  mySphereMeshSource->SetResolution(SphereResolution); 
  mySphereMeshSource->SetScale(scale);

  SimplexFilterType::Pointer simplexFilter = SimplexFilterType::New();
  simplexFilter->SetInput( mySphereMeshSource->GetOutput() );
  simplexFilter->Update();
  
  m_SimplexMesh = simplexFilter->GetOutput();
  m_SimplexMesh->DisconnectPipeline();

  //Now precompute the gradient directions...
  unsigned int numPoints = m_SimplexMesh->GetPoints()->Size();

  PeakDirectionType           gradDir;
  SimplexMeshType::PointType  point;

  typedef  SimplexMeshType::NeighborListType              NeighborsListType;
  NeighborsListType* neighbors = 0;

  m_ListOfNeighbors.clear();
  m_SampleDirections.clear();

  for (unsigned int pointIndex = 0; pointIndex < numPoints; pointIndex++)
  {
    m_SimplexMesh->GetPoint(pointIndex,&point);
    gradDir[0] = point[0]; gradDir[1] = point[1]; gradDir[2] = point[2];
    gradDir.normalize();
    m_SampleDirections.push_back(gradDir);
    m_ListOfNeighbors.push_back( m_SimplexMesh->GetNeighbors( pointIndex, Rigidity) );
  }
}

template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionContainerType
PeakFindingCalculatorGrid< TPixelType >
::ComputePeaks(const PixelType rsh) const
{
  PeakValueContainerType peakVals;
  return ComputePeaks(rsh,peakVals);
}

template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionContainerType
PeakFindingCalculatorGrid< TPixelType >
::ComputePeaksRaw(const PixelType rsh) const
{
  PeakValueContainerType peakVals;
  return ComputePeaksRaw(rsh,peakVals);
}


template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionType
PeakFindingCalculatorGrid< TPixelType >
::FindCriticalDir(const PixelType rsh, const PeakDirectionType inDir) const
{
  PeakValueType peakValue;
  return FindCriticalDir(rsh,inDir,peakValue);
}

template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionType
PeakFindingCalculatorGrid< TPixelType >
::FindPeak(const PixelType rsh, const PeakDirectionType inDir) const
{
  PeakValueType peakValue;
  return FindPeak(rsh,inDir,peakValue);
}

template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionContainerType
PeakFindingCalculatorGrid< TPixelType >
::ComputePeaks(const PixelType rsh, PeakValueContainerType& peakVals) const
{
  PeakDirectionContainerType rawPeaks = ComputePeaksRaw(rsh);
  PeakDirectionContainerType peaks;
  peakVals.clear();
  PeakDirectionType tmpDir;
  
  typename PeakDirectionContainerType::const_iterator iter;
  for (iter = rawPeaks.begin(); iter != rawPeaks.end(); ++iter )
  {
    PeakValueType tmpVal;// = rsh.Evaluate(*iter);
    tmpDir = FindCriticalDir( rsh, (*iter),  tmpVal);
      
    if (addPeak(peaks,tmpDir) && ConfirmPeakUsingCurvature(rsh,tmpDir))
    {
      peaks.push_back(tmpDir);
      peakVals.push_back(tmpVal);
    }
  }
  return peaks;
}


template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionContainerType
PeakFindingCalculatorGrid< TPixelType >
::ComputePeaksRaw(const PixelType rsh, PeakValueContainerType& peakVals) const
{

  typename PeakDirectionContainerType::const_iterator iter;
  typename std::vector<double>::const_iterator peakValIter;

  PeakDirectionContainerType rawPeaks;
  PeakValueContainerType allPeakVals;
  double maxPeakVal = -1000;
  unsigned int numPoints = m_SampleDirections.size();

  //evaluate the pixel at each direction..
  vnl_vector<double> rshValues(numPoints);

  for (unsigned int pointIndex = 0; pointIndex < numPoints; pointIndex++)
  {
    rshValues[pointIndex] = rsh.Evaluate( m_SampleDirections[pointIndex] );
  }

  bool isPeak;
  NeighborsListType* neighbors = 0;

  for (unsigned int pointIndex = 0; pointIndex < numPoints; pointIndex++)
  {
    isPeak = true;
    neighbors = m_ListOfNeighbors[pointIndex];
    for (unsigned int neighborIndex = 0; neighborIndex < neighbors->size(); ++neighborIndex)
    {
      if (rshValues[ (*neighbors)[neighborIndex] ] >= rshValues[ pointIndex ] )
      {
        isPeak = false;
        break;
      }
    }
    if (isPeak)
    {
      rawPeaks.push_back(m_SampleDirections[pointIndex]);
      allPeakVals.push_back( rshValues[pointIndex] );
      if (rshValues[pointIndex] > maxPeakVal)
      {
        maxPeakVal = rshValues[pointIndex];
      }
    }
  }  

  PeakDirectionContainerType peaks;
  peakVals.clear();
  peakValIter = allPeakVals.begin();
  for (iter = rawPeaks.begin(); iter != rawPeaks.end(); ++iter )
  {
    if ( addPeak(peaks,*iter) )
    {
      peaks.push_back(*iter);
      peakVals.push_back(*peakValIter);
    }

    ++peakValIter;
  }
  return peaks;
}

template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionType
PeakFindingCalculatorGrid< TPixelType >
::FindCriticalDir(const PixelType rsh, const PeakDirectionType inDir, PeakValueType& peakValue) const
{
  for ( unsigned int i = 0; i < PixelType::Dimension ; i ++ )
  {
    m_PeakFinder->SetCoefficient(  i,  rsh.GetNthComponent( i )  );
  }
  double inpVec[3];
  inpVec[0] = inDir[0]/inDir.two_norm();
  inpVec[1] = inDir[1]/inDir.two_norm();
  inpVec[2] = inDir[2]/inDir.two_norm(); 
  double outVec[3];
  double val=0;
  int stat = m_PeakFinder->FindPeak( inpVec, outVec, &val, 100 );
  if (stat != 1)
    itkExceptionMacro("FindPeak Failed" << std::endl << "\t RSH : " << rsh << "\t Direction : "<<inDir);

  //find the best Peak...
  PeakDirectionType dir;
  dir[0] = outVec[0]; 
  dir[1] = outVec[1]; 
  dir[2] = outVec[2]; 

  dir.normalize();
  peakValue = static_cast<PeakValueType>(val);
  return dir;
}

template < typename TPixelType >
bool
PeakFindingCalculatorGrid< TPixelType >
::addPeak(const PeakDirectionContainerType peakSet, const PeakDirectionType peak ) const
{
  typename PeakDirectionContainerType::const_iterator iter;
  bool newDir = true;
  for (iter = peakSet.begin(); iter != peakSet.end(); ++iter )
  {
    double dot = vcl_abs(dot_product(peak,*iter));
    if (dot>0.999) //remove directions that differ by less then ~2.5 degrees
    {
      newDir = false;
      break;
    }
  }
  if (newDir)
    return true;

  return false;
}


template < typename TPixelType >
bool
PeakFindingCalculatorGrid< TPixelType >
::ConfirmPeakUsingCurvature(const PixelType rsh, const PeakDirectionType peak) const
{
  //TODO look at this!!!
  return true;
  // double k1,k2,K,H;
  // //Check the curvature
  // rsh.ComputeCurvatures(peak, H, K, k1, k2);
  // std::cout << peak << " : " << k1 <<  " : " << k2 << " : " << (k2 * rsh[0] * vcl_sqrt( 1. / 4. / vnl_math::pi)) << std::endl;
  // if (  (k2 * rsh[0] * vcl_sqrt( 1. / 4. / vnl_math::pi)  > 1)  )
  // {
  //   return true;
  // }
  // else
  // {
  //   return false;
  // }

}

template < typename TPixelType >
typename PeakFindingCalculatorGrid< TPixelType >::PeakDirectionType
PeakFindingCalculatorGrid< TPixelType >
::FindPeak(const PixelType rsh, const PeakDirectionType inDir, PeakValueType& peakValue) const
{
  PeakDirectionType dir;
  //First find the nearest critical value...
  dir = FindCriticalDir(rsh,inDir,peakValue);

  if ( ConfirmPeakUsingCurvature(rsh,dir) )
  {
    //Note that peakValue is already correctly set
    return dir;
  }

  //We need to pick the closest actual peak!
  PeakDirectionContainerType allPeaks = ComputePeaks(rsh);

  double bestDot = 0;
  typename PeakDirectionContainerType::const_iterator iter;

  for (iter = allPeaks.begin(); iter != allPeaks.end(); ++iter )
  {
    double dot = vcl_abs(dot_product(inDir,*iter));
    if (dot>bestDot)
    {
      bestDot = dot;
      dir = *iter;
    }
  }
  peakValue = rsh.Evaluate(dir);
  return dir;
}


} // end namespace


#endif //__itkPeakFindingCalculatorGrid_hxx
