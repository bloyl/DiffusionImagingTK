/**
* @file  itkFiberTrackingManager.h
* @brief This class is used to generate a set of fibers using supplied seed masks etc
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkFiberTrackingManager_h
#define __itkFiberTrackingManager_h


#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkVectorContainer.h>
#include <itkGroupSpatialObject.h>
#include <itkPointSet.h>

// #include <vtkSmartPointer.h>
// #include <vtkPolyLine.h>
// #include <vtkPolyData.h>

#include <itkFiberGenerator.h>

namespace itk
{

/**
 * @class itk::FiberTrackingManager
* @brief This class is used to generate a set of fibers using supplied seed masks etc
 * 
 **/
template < class TDirPicker, typename TLabelType = unsigned int >
class ITK_EXPORT FiberTrackingManager:public Object
{
public:

  /** Standard class typedefs. */
  typedef FiberTrackingManager<TDirPicker, TLabelType>      Self;
  typedef Object                                            Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Typedef to hold SpatialObjects
   * such as seed masks and halting criteria masks
   */
  typedef GroupSpatialObject<3>                             GroupSpatialObjectType;
  typedef GroupSpatialObjectType::Pointer                   GroupSpatialObjectPointer;
  typedef GroupSpatialObjectType::ConstPointer              GroupSpatialObjectConstPointer;

  typedef TDirPicker                                        DirectionPickerType;
  typedef typename DirectionPickerType::ConstPointer        DirectionPickerConstPointer;
  typedef typename DirectionPickerType::PointType           PointType;

  typedef FiberGenerator<DirectionPickerType, TLabelType>   FiberGeneratorType;
  typedef typename FiberGeneratorType::Pointer              FiberGeneratorPointer;

  typedef typename FiberGeneratorType::LabelSetType         LabelSetType;
  typedef itk::VectorContainer<unsigned int,LabelSetType>   VectorOfLabelSetsType;
  typedef typename VectorOfLabelSetsType::Pointer           VectorOfLabelSetsPointer;
  typedef typename VectorOfLabelSetsType::ConstPointer      VectorOfLabelSetsConstPointer;
  
  typedef typename FiberGeneratorType::LabelPixelType       LabelPixelType;
  typedef typename FiberGeneratorType::LabelImageType       LabelImageType;
  typedef typename FiberGeneratorType::LabelImageConstPointer
                                                            LabelImageConstPointer;

  typedef itk::PointSet< unsigned short, 3 >                PointSetType;
  typedef typename PointSetType::Pointer                    PointSetPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(FiberGenerator, Object);

  //****************************************************************************
  //Public METHODS  
  //****************************************************************************

  /** Generate a fiber from a seed Pt in either one or both directions. The
   *  initial direction is choosen using the direction picker. */
  void GenerateFibers( );

  /** Add Seed Image criteria via threshholder input imgages */
  template <typename  TSeedPixelType>
  void AddSeedImage(const typename Image<TSeedPixelType,3>::Pointer, const TSeedPixelType, const TSeedPixelType);

  /** Add Stopping criteria via threshholder input imgages*/
  template <typename  TMaskPixelType>
  void AddStoppingCriteria(typename Image<TMaskPixelType,3>::Pointer, TMaskPixelType, TMaskPixelType);

  /** get NumberOfINcludeROis a fiber must touch inorder to be retained*/
  itkGetConstMacro(NumberOfRequiredIncludeROIs, int);
  
  /** set NumberOfINcludeROis a fiber must touch inorder to be retained*/
  itkSetMacro(NumberOfRequiredIncludeROIs, int);

  /** get NumberOfFibers */
  itkGetConstMacro(NumberOfFibersToKeep, unsigned int);
  
  /** set NumberOfFibers */
  itkSetMacro(NumberOfFibersToKeep, unsigned int);

  /** get NumberOfFibers */
  itkGetConstMacro(NumberOfFibersToGen, unsigned int);
  
  /** set NumberOfFibers */
  itkSetMacro(NumberOfFibersToGen, unsigned int);

  /** get StepLength */
  itkGetConstMacro(StepLength, double);
  
  /** set StepLength */
  itkSetMacro(StepLength, double);
  
  /** get Minimum fiber Length */
  itkGetConstMacro(MinimumFiberLength, double);
  
  /** set Minimum fiber Length */
  itkSetMacro(MinimumFiberLength, double);

  /** get CurvatureThreshold */
  itkGetConstMacro(CurvatureThreshold, double);
  
  /** set CurvatureThreshold
   *  Fibers are stoped when the locale curvature exceeds this value
   *  Curvature is computed as :
   *    acos(theta) / stepLength
   *  where theta is the deflection angle */
  itkSetMacro(CurvatureThreshold, double);

  /** Set the Direction Picker. */
  itkSetConstObjectMacro(DirectionPicker, DirectionPickerType);
  itkGetConstObjectMacro(DirectionPicker, DirectionPickerType);

  /** Return the Fibers */
  // itkGetConstReferenceMacro(Fibers,vtkSmartPointer<vtkPolyData>);
  // itkGetConstObjectMacro(Fibers,GroupSpatialObjectType);
  itkGetObjectMacro(Fibers,GroupSpatialObjectType);
  
  /** Return the Rejected Fibers */
  // itkGetConstReferenceMacro(RejectedFibers,vtkSmartPointer<vtkPolyData>);
  // itkGetConstObjectMacro(RejectedFibers,GroupSpatialObjectType);
  itkGetObjectMacro(RejectedFibers,GroupSpatialObjectType);

  /** Get the stopping masks (maybe move to protected)*/
  itkGetObjectMacro(StoppingMasks, GroupSpatialObjectType);

  /** Get the seed masks (maybe move to protected)*/
  itkGetConstObjectMacro(SeedMasks, GroupSpatialObjectType);

  /** Add an image of rois */
  itkSetConstObjectMacro(ROIImage, LabelImageType);
  itkGetConstObjectMacro(ROIImage, LabelImageType);  

  itkGetConstReferenceMacro(IncludeROIs, LabelSetType);  
  itkGetConstReferenceMacro(RequiredIncludeROIs, LabelSetType);  
  itkGetConstReferenceMacro(ExcludeROIs, LabelSetType);  

  void AddLabelToIncludeRois( const LabelPixelType );
  void AddLabelToExcludeRois( const LabelPixelType );
  void AddLabelToRequiredIncludeRois( const LabelPixelType );

  void AddLabelToIncludeRois( const LabelSetType );
  void AddLabelToExcludeRois( const LabelSetType );
  void AddLabelToRequiredIncludeRois( const LabelSetType );

  /** Get the seedPoints object */
  itkGetConstObjectMacro(SeedPoints, PointSetType);  

  /** Get the touched regions of each keptFiber */
  itkGetConstObjectMacro(FibersROIs, VectorOfLabelSetsType);  

protected:
  FiberTrackingManager();
  virtual ~FiberTrackingManager() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  PointType GenerateSeedPoint() const;

private:
  GroupSpatialObjectPointer                           m_StoppingMasks;
  GroupSpatialObjectPointer                           m_SeedMasks;

  //Hold the seed Points that we generated
  PointSetPointer                                     m_SeedPoints;

  unsigned int                                        m_NumberOfFibersToKeep;
  unsigned int                                        m_NumberOfFibersToGen;
  
  double                                              m_StepLength;
  double                                              m_CurvatureThreshold;
  double                                              m_MinimumFiberLength;

  //Hold the Fiber bundles as vtkPolyData
  // vtkSmartPointer<vtkPolyData>                        m_Fibers;
  GroupSpatialObjectPointer                           m_Fibers;

  //Hold the Rejected Fiber bundles as vtkPolyData
  // vtkSmartPointer<vtkPolyData>                        m_RejectedFibers;
  GroupSpatialObjectPointer                           m_RejectedFibers;

  //Hold the directionPicker Object
  DirectionPickerConstPointer                         m_DirectionPicker;

  //Hold the spatial masks for both exclusion masks and ROI Labels
  LabelImageConstPointer                              m_ROIImage;

  //Store 
  LabelSetType                                        m_RequiredIncludeROIs;
  LabelSetType                                        m_IncludeROIs;
  LabelSetType                                        m_ExcludeROIs;
  VectorOfLabelSetsPointer                            m_FibersROIs;
  
  int                                                 m_NumberOfRequiredIncludeROIs;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberTrackingManager.hxx"
#endif

#endif
