/**
* @file  itkDiffusionModelCalculator.h
* @brief This class is designed to faciliated computing diffusion models
* From diffusion weighted MRI
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkDiffusionModelCalculator_h
#define __itkDiffusionModelCalculator_h

#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkVectorContainer.h>
#include <itkDiffusionTensor3D.h>

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkSymRealSphericalHarmonicRep.h>

namespace itk
{

/**
 * @class itkDiffusionModelCalculator
 * @brief 
 * @TODO use a measurment frame variable!
 **/
template < class TDWIPixelType, class TPrecisionType=double, unsigned int TOutputOrder=4 >
class ITK_EXPORT DiffusionModelCalculator:public Object
{
public:

  /** Standard class typedefs. */
  typedef DiffusionModelCalculator <TDWIPixelType, TPrecisionType, TOutputOrder>
                                                            Self;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(DiffusionModelCalculator, Object);

  typedef TDWIPixelType                                     DWIPixelType;
  typedef TPrecisionType                                    PrecisionType;
  
  //OUTPUT Types
  typedef DiffusionTensor3D<PrecisionType>                  DtType;
  typedef SymRealSphericalHarmonicRep< PrecisionType, TOutputOrder >
                                                            RshType;

  /** Holds each magnetic field gradient used to acquire one DWImage */
  typedef vnl_vector_fixed< PrecisionType, 3 >              GradientDirectionType;

  /** Container to hold gradient directions of the 'n' DW measurements */
  typedef VectorContainer< unsigned int,
                                GradientDirectionType >     GradientDirectionContainerType;

  /** Enums for different reconstructions we know */
  typedef  enum { DT_OLS, DT_WLS, DT_WLSBOOT }              TensorReconEnumType;
  typedef  enum { RSH_ADC, RSH_ODF, RSH_CSAODF }            RSHReconEnumType;

  //****************************************************************************
  //Public METHODS  
  //****************************************************************************

  /** Set the Gradient Directions.
   *  The length of each direction should be the sqrt of the Bvalue
   *  These directions are assumed to be expressed in the LPS world coordinate frame
   *  unless the measurement frame is specified.
   */
  itkSetObjectMacro(GradientDirectionContainer, GradientDirectionContainerType);

  /** Get the Gradient Direction Container */
  itkGetConstObjectMacro(GradientDirectionContainer, GradientDirectionContainerType);

  /** Initialize the internal variables for different model estimation */
  void InitializeTensorFitting();
  void InitializeRSHFitting(PrecisionType);

  /** Set/Get the TensorReconMethod
   *   Controls method used by ComputeTensor
   */
  itkSetMacro(TensorReconMethod, TensorReconEnumType);
  itkGetConstMacro(TensorReconMethod, TensorReconEnumType);

  /** Set/Get the RSHReconMethod
   *   Controls method used by ComputeRSH
   */
  itkSetMacro(RSHReconMethod, RSHReconEnumType);
  itkGetConstMacro(RSHReconMethod, RSHReconEnumType);

  //Figure out how to return residuals...
  DtType ComputeTensor(DWIPixelType, bool) const;
  DtType ComputeTensorOLS(DWIPixelType, bool) const;
  DtType ComputeTensorWLS(DWIPixelType, bool) const;
  DtType ComputeTensorWLS_residualBoot(DWIPixelType, bool) const;

  RshType ComputeRSH(DWIPixelType ) const;
  RshType ComputeRSH_ADC(DWIPixelType ) const;
  RshType ComputeRSH_ODF(DWIPixelType ) const;
  RshType ComputeRSH_CSAODF(DWIPixelType, PrecisionType, PrecisionType ) const;

protected:

  //Needed for generating randomn numbers...
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  typedef GeneratorType::Pointer                                 GeneratorPointer; 

  /** Holds the tensor/RSH basis coefficients and inverses*/
  typedef vnl_matrix< PrecisionType >                            CoefficientMatrixType;

  //Protected methods
  DiffusionModelCalculator();
  virtual ~DiffusionModelCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  DtType ForceTensorToSPD(DtType) const;

  //Convience methods for RSH fitting
  itkStaticConstMacro(RSHNumberOfCoefficients,unsigned int, RshType::Dimension);
  itkStaticConstMacro(RSHNumberOfOrders,unsigned int, RshType::NumberOfOrders);
  itkStaticConstMacro(RSHMaxOrder,unsigned int, RshType::MaxOrder);

private:
  /* Randomn Number Generators */
  GeneratorPointer                                  m_NumGenerator;

  /* Tensor basis coeffs */
  CoefficientMatrixType                             m_BMatrix;
  CoefficientMatrixType                             m_BMatrixInverse;

  /* RSH fitting stuff */
  CoefficientMatrixType                             m_RshBasis;
  CoefficientMatrixType                             m_RshBasisPseudoInverse;

  // Compute the indicies of the B0 images and diff weighted images
  std::vector<unsigned int>                         m_B0Indices; 
  std::vector<unsigned int>                         m_DiffWeightedIndices;
  
  /** Lambda for Controlling Beltrami Regularization */
  double                                            m_BeltramiLambda;

  /* Container to hold gradient directions */
  typename GradientDirectionContainerType::Pointer m_GradientDirectionContainer;

  TensorReconEnumType                               m_TensorReconMethod;
  RSHReconEnumType                                  m_RSHReconMethod;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionModelCalculator.hxx"
#endif

#endif
