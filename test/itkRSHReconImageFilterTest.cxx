
#include <itkTestingMacros.h>

#include <itkRSHReconImageFilter.h>
#include <itkDiffusionModelCalculator.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkVariableLengthVector.h>

#include <itkVariableLengthVectorCastImageFilter.h>
#include <itkCastImageFilter.h>

#include <vnl/vnl_random.h>

#include "testingUtils.h"

using namespace DiffusionImagingTK_testing;

int itkRSHReconImageFilterTest( int , char * argv[] )
{

  std::string dwiFile  = argv[1];
  std::string bvecFile = argv[2];
  std::string bvalFile = argv[3];
  std::string maskFile = argv[4];

  //Typedefs for dwi images and conversion!
  typedef short                                                     DwiComponentType;
  typedef itk::VectorImage<DwiComponentType, 3>                     DwiImageType;
  typedef itk::Image<DwiComponentType, 4>                           DwiImageType4D;

  //Typedefs for Diffusion model Calculator
  typedef DwiImageType::PixelType                                   DwiPixelType;
  typedef itk::DiffusionModelCalculator<DwiPixelType, double>     	DTCalculatorType;
  typedef DTCalculatorType::GradientDirectionContainerType          GradientDirectionContainerType;
  
  typedef itk::RSHReconImageFilter<DTCalculatorType,3>			        RSHReconFilterType;
  typedef DwiImageType::Pointer                                     DwiImagePointer;

  RSHReconFilterType::Pointer rshFilter = RSHReconFilterType::New();
  EXERCISE_BASIC_OBJECT_METHODS( rshFilter, RSHReconFilterType );
  TRY_EXPECT_EXCEPTION( rshFilter->Update() );

  // //Load in dwi File
  typedef itk::ImageFileReader< DwiImageType4D >                    DwiReaderType;
  DwiReaderType::Pointer dwiReader = DwiReaderType::New();
  dwiReader->SetFileName( dwiFile );
  dwiReader->Update();

  DwiImagePointer dwiImage = DwiImageType::New();
  convertDWIimage<DwiImageType4D,DwiImageType>(dwiReader->GetOutput(), dwiImage);

  rshFilter->SetInput(dwiImage);
  TRY_EXPECT_EXCEPTION( rshFilter->Update() );

  //Setup diffusion model calculator
  DTCalculatorType::Pointer dtCalc = DTCalculatorType::New();
  rshFilter->SetDiffusionModelCalculator(dtCalc);
  TRY_EXPECT_EXCEPTION( rshFilter->Update() );
  
  GradientDirectionContainerType::Pointer gradDirs = GradientDirectionContainerType::New();
  loadGradDirs<GradientDirectionContainerType>( gradDirs, bvecFile, bvalFile );
  dtCalc->SetGradientDirectionContainer(gradDirs);

  TRY_EXPECT_NO_EXCEPTION( rshFilter->Update() );

  typedef RSHReconFilterType::OutputImageType                       OutputImageType;


  typedef itk::VectorImage
    <OutputImageType::PixelType::ComponentType, 3>
                                                            VectorImageType;

  typedef itk::VariableLengthVectorCastImageFilter<OutputImageType,VectorImageType> CasterType;
  // typedef itk::CastImageFilter <OutputImageType,VectorImageType> CasterType;

  typedef itk::ImageFileWriter< VectorImageType > WriterType;

  //We need to cast the RSH output to vectorImage.
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(rshFilter->GetOutput() );
  caster->Update();
  
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[5]);
  writer->SetInput( caster->GetOutput() );
  writer->Update();


  return EXIT_SUCCESS;
}