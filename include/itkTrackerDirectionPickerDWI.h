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

#ifndef __itkTrackerDirectionPickerDWI_h
#define __itkTrackerDirectionPickerDWI_h

#include <itkTrackerDirectionPickerImageBase.h>
#include <itkDiffusionModelCalculator.h>
#include <itkPeakFindingCalculatorGrid.h>
#include <itkPeakFindingCalculator.h>

#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkVectorContainer.h>

#include <itkDiffusionTensor3D.h>

#include <vector>

namespace itk
{

/**
 * @class itkTrackerDirectionPickerDWI
 * @brief This class is serves as the base for choosing 3 dimensional directions
 * from DWI images for use in fiber tracking...
 * 
 * this class is inherently stocastic *using bootstrapping*
 **/
template <class TInterpType, class TPrecisionType=double, unsigned int TOutputOrder=4 >
class ITK_EXPORT TrackerDirectionPickerDWI:public TrackerDirectionPickerImageBase<TInterpType>
{
public:

  /** Standard class typedefs. */
  typedef TrackerDirectionPickerDWI <TInterpType,TPrecisionType, TOutputOrder>
                                                            Self;
                                                            
  typedef TrackerDirectionPickerImageBase<TInterpType>      Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(TrackerDirectionPickerDWI, TrackerDirectionPickerImageBase);

  /** typedef alias for the source data container */
  typedef typename Superclass::ImageType                    ImageType;
  typedef typename Superclass::InterpolatorType             InterpolatorType;
  typedef typename Superclass::InterpConstPointer           InterpConstPointer;

  typedef typename Superclass::PointType                    PointType;
  typedef typename Superclass::DirectionType                DirectionType;  
  typedef typename Superclass::PixelType                    PixelType;

  typedef DiffusionModelCalculator<PixelType, TPrecisionType, TOutputOrder>
                                                            DiffModCalculatorType;
  typedef typename DiffModCalculatorType::GradientDirectionContainerType   
                                                            GradientDirectionContainerType;
  typedef typename DiffModCalculatorType::DtType            DtType;
  typedef typename DiffModCalculatorType::RshType           RshType;

  /** enum for support picking methods */
  typedef  enum { DTI_STT, DTI_TEND, DTI_TLINE, CSAODF_SST, ODF_SST } PickingMethodEnumType;

  //****************************************************************************
  //Public METHODS  
  //****************************************************************************

  /** Returns the direction based on the diffusion model at supplied point
    *  All determinisitic methods use the principle eigenvector (e_1) as
    *  the starting direction. Probablisitic method uses its own function.
    */ 
  DirectionType PickStartingDirection (PointType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction */ 
  DirectionType PickNextDirection (PointType, DirectionType) const;

  itkSetMacro(Method, PickingMethodEnumType);
  itkGetConstMacro(Method, PickingMethodEnumType);

  itkSetMacro(TensorLineF, double);
  itkGetConstMacro(TensorLineF, double);
  
  itkSetMacro(TensorLineG, double);
  itkGetConstMacro(TensorLineG, double);

  virtual void SetGradientDirectionContainer(GradientDirectionContainerType * _arg)
  {
    this->m_ModelCalculator->SetGradientDirectionContainer(_arg);
  }
  
  virtual const GradientDirectionContainerType * GetGradientDirectionContainer() const
  {
    return this->m_ModelCalculator->GetGradientDirectionContainer();
  }

  void UseDTIStreamLineTracking();
  void UseDTITensorDeflectionTracking();
  void UseDTITensorLineTracking();
  void UseCSAODFStreamLineTracking();
  void UseODFStreamLineTracking();

  void Initialize();

protected:

  //Protected Typedefs
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator    GeneratorType;
  typedef GeneratorType::Pointer                                    GeneratorPointer; 

  typedef itk::PeakFindingCalculatorGrid<RshType>                   PeakFinderType;
  typedef typename PeakFinderType::PeakDirectionType                PeakDirectionType;
  typedef typename PeakFinderType::PeakDirectionContainerType       PeakDirectionContainerType;
  typedef typename PeakFinderType::PeakValueContainerType           PeakValueContainerType;


  //Protected methods
  TrackerDirectionPickerDWI();
  virtual ~TrackerDirectionPickerDWI() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the SST algorithm (which just returns e_1)
   */ 
  DirectionType PickNextDirectionDTI_SST (PointType, DirectionType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the Tensor deflection algorithm ( v_{out} = D * v_{in} )
   */ 
  DirectionType PickNextDirectionDTI_TEND (PointType, DirectionType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the Tensor stream line algorithm
   *    required scalar parameters f and g to be set!
   *    v_{out} = f e1 + ͑ (1 - f͒͑͑) ( (1-g) ͒v_{in} + g D v_{in})
   */ 
  DirectionType PickNextDirectionDTI_TLINE (PointType, DirectionType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the SST algorithm (which just returns the closes peak to
   *  incoming direction
   */ 
  DirectionType PickNextDirectionRSH_SST (PointType, DirectionType) const;

private:
  PickingMethodEnumType                             m_Method;
  GeneratorPointer                                  m_NumGenerator;

  /** container to hold gradient directions */
  typename DiffModCalculatorType::Pointer           m_ModelCalculator;


  /** Tensor Lines parameters               */
  double                                            m_TensorLineF;
  double                                            m_TensorLineG;

  bool                                              m_isInitialized;

  typename PeakFinderType::Pointer                  m_PeakCalc;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrackerDirectionPickerDWI.hxx"
#endif

#endif
