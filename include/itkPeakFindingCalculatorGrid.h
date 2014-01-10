/**
 * @file  itkPeakFindingCalculatorGrid.h
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPeakFindingCalculatorGrid_h
#define __itkPeakFindingCalculatorGrid_h

#include <math.h>
#include <iostream>

#include <itkSymRealSphericalHarmonicRep.h>
#include <itkPixelReorientationOperator.h>

//TODO Possibly remove this!
#include <itkPeakFindingCalculator.h>

#include "itkMesh.h"
#include "itkSimplexMesh.h"
#include "itkSimplexMeshGeometry.h"
#include "itkRegularSphereMeshSource.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkDefaultDynamicMeshTraits.h"

#include <vector>

namespace itk
{

template < typename TPixelType >
class ITK_EXPORT PeakFindingCalculatorGrid: public Object
{
public:

  // Standard class typedefs.
  typedef PeakFindingCalculatorGrid< TPixelType >
                                                        Self;
  typedef SmartPointer < Self >                         Pointer;
  typedef SmartPointer < const Self >                   ConstPointer;

  // Input and output types
  // 
  // Input:  a total of ( MaxOrder + 1 ) * ( MaxOrder + 2 ) / 2 spherical harmonic
  //         coefficients as a representation of an OD function (either FOD or ODF).
  //
  // Output: the peaks (local extrema) of the OD function in the form of normalized
  //         3D vectors/directions and the associated magnitudes (OD values).
  //
  // NOTE:   Each fiber results in a pair of mutually-opposite directions, both with
  //         positive OD values. This rule is used to remove fake/artifact peaks.
  //
  typedef TPixelType                                          PixelType;
  typedef typename PixelType::ComponentType                   InputComponentType;
  typedef typename PixelType::GradientDirectionType           PeakDirectionType;
  typedef typename PixelType::RealValueType                   PeakValueType;
  
  typedef typename std::vector< PeakValueType >               PeakValueContainerType;
  typedef typename std::vector< PeakDirectionType >           PeakDirectionContainerType;
  
  itkStaticConstMacro( MaxOrder,  unsigned int, PixelType::MaxOrder  );
  itkStaticConstMacro( Dimension, unsigned int, PixelType::Dimension );
  
  // Method for creation through the object factory.
  itkNewMacro(Self);  

  // Run-time type information (and related methods).
  itkTypeMacro( PeakFindingCalculatorGrid, Object );

  //Initilize the calculator...
  void Initialize(unsigned int, unsigned int);

  //Compute Peaks ( removing v or -v) using Grid
  PeakDirectionContainerType ComputePeaksRaw(const PixelType, PeakValueContainerType&) const;
  PeakDirectionContainerType ComputePeaksRaw(const PixelType) const;

  //Compute Peaks ( removing v or -v) also do additional gradient ascent
  PeakDirectionContainerType ComputePeaks(const PixelType, PeakValueContainerType&) const;
  PeakDirectionContainerType ComputePeaks(const PixelType) const;

  PeakDirectionType FindPeak(const PixelType, const PeakDirectionType, PeakValueType&) const;
  PeakDirectionType FindPeak(const PixelType, const PeakDirectionType) const;

protected: 

  PeakFindingCalculatorGrid();
  ~PeakFindingCalculatorGrid();

  PeakDirectionType FindCriticalDir(const PixelType, const PeakDirectionType, PeakValueType&) const;
  PeakDirectionType FindCriticalDir(const PixelType, const PeakDirectionType) const;

  bool ConfirmPeakUsingCurvature(const PixelType, const PeakDirectionType) const;

private:

  PeakFindingCalculatorGrid( const Self & ); //purposely not implemented
  void operator = ( const Self & );      //purposely not implemented

  typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double, double>
                                                              MeshTraits;
  typedef itk::Mesh<double,3,MeshTraits>                      TriangleMeshType;
  typedef itk::SimplexMesh<double,3,MeshTraits>               SimplexMeshType;
  typedef SimplexMeshType::Pointer                            SimplexMeshPointer;
  typedef typename SimplexMeshType::NeighborListType          NeighborsListType;

  typedef std::vector< NeighborsListType* >                   ListOfNeighborsListType;

  typedef itk::PeakFindingCalculator < PixelType, double, double >
                                                              PeakFindingCalculatorType;

  typedef itk::PixelReorientationOperator<PixelType>          ReOrienterType;

  //Help funtion for building peak GradientDirs...
  bool addPeak(const PeakDirectionContainerType peakSet, const PeakDirectionType peak ) const;

  //Variables...
  PeakDirectionContainerType                                  m_SampleDirections;
  ListOfNeighborsListType                                     m_ListOfNeighbors;  
  
  SimplexMeshPointer                                          m_SimplexMesh;                                             

  typename PeakFindingCalculatorType::Pointer                 m_PeakFinder;
  typename ReOrienterType::Pointer                            m_PixelReorienter;

}; //Close class
} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPeakFindingCalculatorGrid.hxx"
#endif

  
#endif
