
#include <itkTestingMacros.h>

#include <itkDTIReconImageFilter.h>
#include <itkDiffusionModelCalculator.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkVariableLengthVector.h>
#include <vnl/vnl_random.h>

#include "testingUtils.h"

using namespace DiffusionImagingTK_testing;

int itkDTIReconImageFilterTest( int , char * argv[] )
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
  
  typedef itk::DTIReconImageFilter<DTCalculatorType,3>			        DTIReconFilterType;
  typedef DwiImageType::Pointer                                     DwiImagePointer;

  DTIReconFilterType::Pointer dtiFilter = DTIReconFilterType::New();
  EXERCISE_BASIC_OBJECT_METHODS( dtiFilter, DTIReconFilterType );
  TRY_EXPECT_EXCEPTION( dtiFilter->Update() );

  // //Load in dwi File
  typedef itk::ImageFileReader< DwiImageType4D >                    DwiReaderType;
  DwiReaderType::Pointer dwiReader = DwiReaderType::New();
  dwiReader->SetFileName( dwiFile );
  dwiReader->Update();

  DwiImagePointer dwiImage = DwiImageType::New();
  convertDWIimage<DwiImageType4D,DwiImageType>(dwiReader->GetOutput(), dwiImage);

  dtiFilter->SetInput(dwiImage);
  TRY_EXPECT_EXCEPTION( dtiFilter->Update() );

  //Setup diffusion model calculator
  DTCalculatorType::Pointer dtCalc = DTCalculatorType::New();
  dtiFilter->SetDiffusionModelCalculator(dtCalc);
  TRY_EXPECT_EXCEPTION( dtiFilter->Update() );
  
  GradientDirectionContainerType::Pointer gradDirs = GradientDirectionContainerType::New();
  loadGradDirs<GradientDirectionContainerType>( gradDirs, bvecFile, bvalFile );
  dtCalc->SetGradientDirectionContainer(gradDirs);

  TRY_EXPECT_NO_EXCEPTION( dtiFilter->Update() );

  typedef DTIReconFilterType::OutputImageType                       OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType >                   DTIWriterType;

  DTIWriterType::Pointer writer = DTIWriterType::New();
  writer->SetFileName(argv[5]);
  writer->SetInput( dtiFilter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}