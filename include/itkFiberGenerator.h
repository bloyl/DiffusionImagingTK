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
 
#ifndef __itkFiberGenerator_h
#define __itkFiberGenerator_h

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkImage.h>
#include <itkImageMaskSpatialObject.h>
#include <itkGroupSpatialObject.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLineSpatialObject.h>

#include <itksys/hash_set.hxx>

///TODO this should overload the modified() method to unset m_intialized...

namespace itk
{

/**
 * @class itkFiberGenerator
* @brief This class is used to generate single fiber tracts from
* diffuion images such as DTI or FOD images 
 * 
 **/
template < class TDirPicker, typename TLabelType = unsigned int >
class ITK_EXPORT FiberGenerator:public Object
{
public:

  /** Standard class typedefs. */
  typedef FiberGenerator<TDirPicker, TLabelType>            Self;
  typedef Object                                            Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  typedef TDirPicker                                        DirectionPickerType;
  typedef typename DirectionPickerType::PointType           PointType;
  typedef typename DirectionPickerType::DirectionType       DirectionType;

  typedef typename DirectionPickerType::ConstPointer        DirectionPickerConstPointer;

  enum StartDirectionChoice {ONE_DIRECTION, BOTH_DIRECTIONS};

  /** Typedefs for storing ROI masks (Could probaly be protected) */
  typedef GroupSpatialObject<3>                             StoppingMaskType;
  typedef StoppingMaskType::Pointer                         StoppingMaskPointer;
  typedef StoppingMaskType::ConstPointer                    StoppingMaskConstPointer;

  typedef TLabelType                                        LabelPixelType;
  typedef Image<LabelPixelType,3>                           LabelImageType;
  typedef typename LabelImageType::ConstPointer             LabelImageConstPointer;

  typedef itksys::hash_set< TLabelType >                    LabelSetType;

  typedef itk::LineSpatialObject< 3 >                       LineType;
  typedef typename LineType::Pointer                        LinePointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(FiberGenerator, Object);

  //****************************************************************************
  //Public METHODS  
  //****************************************************************************

  /** get StepLength */
  itkGetConstMacro(StepLength, double);
  
  /** set StepLength */
  itkSetMacro(StepLength, double);
  
  /** get CurvatureThreshold */
  itkGetConstMacro(CurvatureThreshold, double);
  
  /** set CurvatureThreshold
   *  Fibers are stoped when the locale curvature exceeds this value
   *  Curvature is computed as :
   *    acos(theta) / stepLength
   *  where theta is the deflection angle
   */
  itkSetMacro(CurvatureThreshold, double);

  /** Compute the Fiber length. */
  double GetFiberLength( ) const;

  /** Generate a fiber from a seed Pt in either one or both directions. The
   *  initial direction is choosen using the direction picker. */
  void GenerateFiberPoints( const PointType, const StartDirectionChoice);
  
  /** Generate a fiber from a seed Pt and initial direction in either one 
   *  or both directions. */
  void GenerateFiberPoints( const PointType, const DirectionType, const StartDirectionChoice);
  
  /** Set the Direction Picker. */
  itkSetConstObjectMacro(DirectionPicker, DirectionPickerType);
  itkGetConstObjectMacro(DirectionPicker, DirectionPickerType);
  
  /** Get the stopping masks (maybe move to protected)*/
  itkGetConstObjectMacro(StoppingMask, StoppingMaskType);

  /** Add a Stopping criteria from an itkSpatial object */
  void AddStoppingCriteriaSpatialObject(typename SpatialObject<3>::Pointer);

  /** Add an image of rois */
  itkSetConstObjectMacro(ROIImage, LabelImageType);
  itkGetConstObjectMacro(ROIImage, LabelImageType);  

  /** Return a vtkPolyData Structure with this fiber in it*/
  LinePointer GetFiber( ) const;

  void Initialize();
  // void Modified( );

  /** Get list of traversed rois */
  itkGetConstReferenceMacro(TraversedROIs, LabelSetType);
  
protected:

  //Protected Typedefs
  typedef typename LineType::LinePointType                  LinePointType;
  typedef std::vector<LinePointType>                        LinePntsContainerType;

  //Protected methods
  
  /** Compute the Curvature from the incoming and outGoing Directions */
  double ComputeCurvature( const DirectionType , const DirectionType) const;

  /** generate a set of fiber points by moving in a single direction */
  void GenerateFiberPointsOneWay(PointType, DirectionType, LinePntsContainerType& );

  FiberGenerator();
  virtual ~FiberGenerator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  double                                            m_StepLength;
  double                                            m_CurvatureThreshold;

  //Hold the points of a streamline
  LinePntsContainerType                             m_Points;

  //Hold the directionPicker Object
  DirectionPickerConstPointer                       m_DirectionPicker;

  //Hold the spatial masks for both exclusion masks and ROI Labels
  LabelImageConstPointer                            m_ROIImage;
  StoppingMaskPointer                               m_StoppingMask;

  //Store rois that have been passed thorugh
  LabelSetType                                      m_TraversedROIs;

  bool                                              m_Initialized;


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberGenerator.hxx"
#endif

#endif
