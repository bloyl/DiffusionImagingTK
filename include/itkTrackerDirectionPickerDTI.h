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

#ifndef __itkTrackerDirectionPickerDTI_h
#define __itkTrackerDirectionPickerDTI_h

#include <itkTrackerDirectionPickerImageBase.h>

#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace itk
{

/**
 * @class itkTrackerDirectionPickerDTI
 * @brief This class is serves as the base for choosing 3 dimensional directions
 * from DTI images for us in fiber tracking...
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
 * 2) TDEF - Use the entire DT to deflect the incoming
 * vector (vin) direction [Weinstein et al., 1999]:
 *     v_{out} = D v_{in}
 *
 * 3) TEND - The original tensor-lines algorithm, described by Weinstein et al. [1999]
 * dynamically modulates the STT and TEND contributions to tract steering:
 *     v_{out} =  f e1 + ͑ (1 - f͒͑͑)( (1-g) ͒v_{in} + g D v_{in})
 *              where f and g are in [0,1] and chosen by the user...
 *
 * 4) PROB - This is a probablisitic framework that has not been thoughrough thought out
 * but uses rejection sampling to choose a direciton from the odf coresponding to the tensor
 *  
 **/
template <class TInterpType >
class ITK_EXPORT TrackerDirectionPickerDTI:public TrackerDirectionPickerImageBase<TInterpType>
{
public:

  /** Standard class typedefs. */
  typedef TrackerDirectionPickerDTI <TInterpType>           Self;
  typedef TrackerDirectionPickerImageBase<TInterpType>      Superclass;

  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(TrackerDirectionPickerDTI, TrackerDirectionPickerImageBase);

  /** typedef alias for the source data container */
  typedef typename Superclass::ImageType                    ImageType;
  typedef typename Superclass::InterpolatorType             InterpolatorType;
  typedef typename Superclass::InterpConstPointer           InterpConstPointer;

  typedef typename Superclass::PointType                    PointType;
  typedef typename Superclass::DirectionType                DirectionType;  
  typedef typename Superclass::PixelType                    PixelType;

  /** enum for support picking methods */
  typedef  enum { STT, TEND, TLINE, PROB } PickingMethodEnumType;

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

  void UseStreamLineTracking();
  void UseTensorDeflectionTracking();
  void UseTensorLineTracking();
  void UseProbabilisticTracking();

protected:

  //Protected Typedefs
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  typedef GeneratorType::Pointer                                 GeneratorPointer; 

  //Protected methods
  TrackerDirectionPickerDTI();
  virtual ~TrackerDirectionPickerDTI() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the SST algorithm (which just returns e_1)
   */ 
  DirectionType PickNextDirectionSST (PointType, DirectionType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the Tensor deflection algorithm ( v_{out} = D * v_{in} )
   */ 
  DirectionType PickNextDirectionTEND (PointType, DirectionType) const;

  /** Returns the outgoing direction based on the diffusion model at 
   *  the supplied point and the incoming direction
   *  Uses the Tensor stream line algorithm
   *    required scalar parameters f and g to be set!
   *    v_{out} = f e1 + ͑ (1 - f͒͑͑) ( (1-g) ͒v_{in} + g D v_{in})
   */ 
  DirectionType PickNextDirectionTLINE (PointType, DirectionType) const;

  /** Returns the direction based on the diffusion model at supplied point
    *  uses a probablisic algorithm (rejection Sampling)
    *  and a marginal odf computed from 
    *   Citation----
    *    Reconstruction of the orientation distribution function in single- and multiple-shell q-ball imaging within constant solid angle.
    *    Aganj I, Lenglet C, Sapiro G, Yacoub E, Ugurbil K, Harel N.
    *    Magn Reson Med. 2010 Aug;64(2):554-66.
    *    PMID: 20535807
    */ 
  DirectionType PickStartingDirectionProb (PointType) const;

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
  DirectionType PickNextDirectionProb (PointType, DirectionType) const;

  //Conveince methods for probabilistic tracking
  DirectionType PickRandomDirectionOnSphere () const;
  DirectionType PickRandomFromODF (PointType) const;
  double EvaluateOdf (PixelType, DirectionType) const;

private:

  PickingMethodEnumType m_Method;
  GeneratorPointer      m_NumGenerator;
  double                m_TensorLineF;
  double                m_TensorLineG;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrackerDirectionPickerDTI.hxx"
#endif

#endif
