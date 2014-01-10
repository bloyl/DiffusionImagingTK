/**
* @file  itkTrackerDirectionPickerImageBase.h
* @brief This class is an abstract class that serves as the base for classes used
* to choose 3 dimensional directions from diffusion model images (ie DTI or FOD images)
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkTrackerDirectionPickerImageBase_h
#define __itkTrackerDirectionPickerImageBase_h

#include <itkObject.h>
#include <itkObjectFactory.h>

namespace itk
{

/**
 * @class itkTrackerDirectionPickerImageBase
 * @brief This class is an abstract class that serves as the base for classes used
 * to choose 3 dimensional directions from diffusion model images (ie DTI or FOD images)
 * 
 **/
template <class TInterpType >
class ITK_EXPORT TrackerDirectionPickerImageBase:public Object
{
public:

  //Public typedefs
  /** Standard class typedefs. */
  typedef TrackerDirectionPickerImageBase < TInterpType >   Self;
  typedef Object                                            Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(TrackerDirectionPickerImageBase, Object);

  typedef TInterpType                                       InterpolatorType;
  typedef typename InterpolatorType::InputImageType         ImageType;
  
  typedef typename ImageType::PixelType                     PixelType;
  typedef typename ImageType::PointType                     PointType;
  typedef typename PointType::VectorType                    DirectionType;  

  /** Pointer types for the image. */
  typedef typename ImageType::Pointer                       ImagePointer;
  typedef typename ImageType::ConstPointer                  ImageConstPointer;

  /** Pointer types for the interpolator. */
  typedef typename InterpolatorType::Pointer                InterpPointer;
  typedef typename InterpolatorType::ConstPointer           InterpConstPointer;

  /** Typedefs for measurement frame. */
  typedef typename ImageType::DirectionType                 MeasurementFrameType;
  
  //****************************************************************************
  //Public METHODS  
  //****************************************************************************

  /** Returns the direction based on the diffusion model at supplied point */ 
  virtual DirectionType PickStartingDirection
            (PointType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction */ 
  virtual DirectionType PickNextDirection
            (PointType, DirectionType) const;

  /** Check the Initialization the TrackerDirectionPickerImageBase */ 
  virtual bool IsInitialized( ) const;

  /** Initialize the TrackerDirectionPickerImageBase */ 
  virtual bool CanEvaluate( PointType ) const;

  /** Set the input interpolator. */
//  virtual void SetInterpolator (const InterpolatorType* _arg);
  itkSetConstObjectMacro(Interpolator, InterpolatorType);
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

protected:

  //Protected Typedefs

  //Protected methods
  TrackerDirectionPickerImageBase();
  virtual ~TrackerDirectionPickerImageBase() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  InterpConstPointer                                        m_Interpolator;
  MeasurementFrameType                                      m_MeasurementFrame;
  bool                                                      m_UserSetMeasurmentFrame;
  bool                                                      m_Initialized;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrackerDirectionPickerImageBase.hxx"
#endif

#endif
