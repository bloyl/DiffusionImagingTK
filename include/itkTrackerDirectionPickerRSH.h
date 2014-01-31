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

#ifndef __itkTrackerDirectionPickerRSH_h
#define __itkTrackerDirectionPickerRSH_h

#include <itkTrackerDirectionPickerImageBase.h>
#include <itkPeakFindingCalculatorGrid.h>
#include <itkPeakFindingCalculator.h>

#include <itkMersenneTwisterRandomVariateGenerator.h>


namespace itk
{

/**
 * @class itkTrackerDirectionPicker
 * @brief This class is serves as the base for choosing 3 dimensional directions
 * from RSH images for us in fiber tracking...
 * 
 * Terminology and determinisitic methods are implemented from
 * Lazar, M., Weinstein, D. M., Tsuruda, J. S., Hasan, K. M.,
 * Arfanakis, K., Meyerand, M. E., Badie, B., Rowley, H. A.,
 * Haughton, V., Field, A. and Alexander, A. L. (2003),
 * White matter tractography using diffusion tensor deflection.
 * Human Brain Mapping, 18: 306–321. doi: 10.1002/hbm.10102
 *
 * There are options to choose from...
 * 1) STT - Streamlines tracking (STT) techniques model prop-
 * agation in the major eigenvector field of the brain.
 * The major eigenvector direction is assumed to be tangent to the tract pathway.
 *     v_{out} = e_1
 * 
 * 4) PROB - This is a probablisitic framework that has not been thoughrough thought out
 * but uses rejection sampling to choose a direciton from the odf coresponding to the tensor
 *  
 **/
template <class TInterpType >
class ITK_EXPORT TrackerDirectionPickerRSH:public TrackerDirectionPickerImageBase<TInterpType>
{
public:

  /** Standard class typedefs. */
  typedef TrackerDirectionPickerRSH <TInterpType>           Self;
  typedef TrackerDirectionPickerImageBase<TInterpType>      Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(TrackerDirectionPickerRSH, TrackerDirectionPickerImageBase);

  /** typedef alias for the source data container */
  typedef typename Superclass::ImageType                    ImageType;
  typedef typename Superclass::InterpolatorType             InterpolatorType;
  typedef typename Superclass::InterpConstPointer           InterpConstPointer;

  typedef typename Superclass::PointType                    PointType;
  typedef typename Superclass::DirectionType                DirectionType;  
  typedef typename Superclass::PixelType                    PixelType;

  /** enum for support picking methods */
  typedef  enum { STT, PROB } PickingMethodEnumType;

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

  void UseStreamLineTracking();
  void UseProbabilisticTracking();

protected:

  //Protected Typedefs
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator    GeneratorType;
  typedef GeneratorType::Pointer                                    GeneratorPointer; 

  typedef itk::PeakFindingCalculatorGrid<PixelType>                 PeakFinderType;
  typedef typename PeakFinderType::PeakDirectionType                PeakDirectionType;
  typedef typename PeakFinderType::PeakDirectionContainerType       PeakDirectionContainerType;

  //Protected methods
  TrackerDirectionPickerRSH();
  virtual ~TrackerDirectionPickerRSH() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the SST algorithm (which just returns e_1)
   */ 
  DirectionType PickNextDirectionSST (PointType, DirectionType) const;

  /** Returns the direction based on the diffusion model at supplied point
    *  uses a probablisic algorithm (rejection Sampling)
    *  and a marginal odf computed from 
    *   Citation----
    *    Reconstruction of the orientation distribution function in single- and multiple-shell q-ball imaging within constant solid angle.
    *    Aganj I, Lenglet C, Sapiro G, Yacoub E, Ugurbil K, Harel N.
    *    Magn Reson Med. 2010 Aug;64(2):554-66.
    *    PMID: 20535807
    */ 
  // DirectionType PickStartingDirectionProb (PointType) const;

  /** Returns the outgoing direction based on the diffusion model at 
    *  the supplied point and the incoming direction
    *  uses a probablisic algorithm (rejection Sampling)
    *  and a marginal odf computed from 
    *   Citation----
    *    Reconstruction of the orientation distribution function in single- and multiple-shell q-ball imaging within constant solid angle.
    *    Aganj I, Lenglet C, Sapiro G, Yacoub E, Ugurbil K, Harel N.
    *    Magn Reson Med. 2010 Aug;64(2):554-66.
    *    PMID: 20535807
    */ 
  // DirectionType PickNextDirectionProb (PointType, DirectionType) const;

  //Conveince methods for probabilistic tracking
  // DirectionType PickRandomDirectionOnSphere () const;
  // DirectionType PickRandomFromODF (PointType) const;
  // double EvaluateOdf (PixelType, DirectionType) const;

private:

  PickingMethodEnumType                                     m_Method;
  GeneratorPointer                                          m_NumGenerator;
  typename PeakFinderType::Pointer                          m_PeakCalc;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrackerDirectionPickerRSH.hxx"
#endif

#endif
