/**
 * @file  itkFiberGeneratorTest.cxx
 * @brief Test itkFiberGenerator class.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkFiberGenerator.h>

#include <itkTrackerDirectionPickerDTI.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itksys/SystemTools.hxx>
#include <itkSpatialObjectWriter.h>

#include <itkTestingMacros.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>


namespace{

typedef itk::DiffusionTensor3D<double>                    DtiPixelType;
typedef itk::Image<DtiPixelType,3>                        DtiImageType;
typedef itk::NearestNeighborInterpolateImageFunction<DtiImageType,double>
                                                          DtiNNInterpType;
typedef itk::TrackerDirectionPickerDTI<DtiNNInterpType>   DtiDirPickerType;

typedef itk::ImageMaskSpatialObject<3>                    MaskSpatialObjectType;

int testInitialize()
{

  typedef itk::DiffusionTensor3D<double>   PixelType;
  typedef itk::Image<PixelType,3>          ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                           InterpType;
  
  typedef itk::TrackerDirectionPickerDTI<InterpType>
                                           PickerType;

  typedef  itk::FiberGenerator<PickerType> FiberGeneratorType;
  
  PickerType::Pointer dirPicker = PickerType::New();
  
  FiberGeneratorType::Pointer ftGen = FiberGeneratorType::New();
  std::cout << ftGen << std::endl;
  
  ftGen->SetStepLength(10.2);
  ftGen->SetDirectionPicker(dirPicker);
  
  std::cout << ftGen << std::endl;
  
  return EXIT_FAILURE;
}

MaskSpatialObjectType::Pointer
loadStoppingMask(std::string maskFile, double lT, double hT)
{

  typedef itk::Image< double , 3 >                            ImageMaskType;
  typedef itk::ImageFileReader< ImageMaskType >               ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( maskFile );
  reader->Update();
  
  
  typedef itk::BinaryThresholdImageFilter< ImageMaskType,
       MaskSpatialObjectType::ImageType >           ThresholderType;
  ThresholderType::Pointer thresholder = ThresholderType::New();
  
  thresholder->SetOutsideValue(itk::NumericTraits< MaskSpatialObjectType::ImageType::PixelType>::Zero);
  thresholder->SetInsideValue(itk::NumericTraits< MaskSpatialObjectType::ImageType::PixelType>::One);
  
  //the threshold is inclusive so is >= 1
  thresholder->SetLowerThreshold(lT);
  thresholder->SetUpperThreshold(hT);
  thresholder->InPlaceOn();
  
  thresholder->SetInput( reader->GetOutput() );
  thresholder->Update();

  MaskSpatialObjectType::Pointer spatialObjectMask = MaskSpatialObjectType::New();
  spatialObjectMask->SetImage( thresholder->GetOutput() );
  
  spatialObjectMask->DisconnectPipeline();
  return spatialObjectMask;

}

//void printRois(itk::VectorContainer<unsigned int, unsigned int>::ConstPointer labels)
void printRois(itksys::hash_set< unsigned int > labels)
{
  itksys::hash_set< unsigned int >::const_iterator it;
  it = labels.begin();
  while(it != labels.end())
  {
    std::cout << *it;
    ++it;
    if (it != labels.end()) 
      std::cout << ", ";
  }
  std::cout << std::endl;

}
 
template< class TPickerType >
int test1(std::string dataFile)
{
  typedef TPickerType                                 PickerType;
  typedef typename PickerType::InterpolatorType       InterpolatorType;
  typedef typename PickerType::PointType              PointType;
  typedef typename InterpolatorType::InputImageType   ImageType;
  typedef itk::ImageFileReader< ImageType >           ReaderType;
  typedef itk::FiberGenerator<PickerType>             FiberGeneratorType;
  typedef typename FiberGeneratorType::LinePointer    FiberPointerType;
  typedef itk::SpatialObjectWriter<3>                 WriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( dataFile );
  reader->Update();
  
  typename InterpolatorType::Pointer interp = InterpolatorType::New();
  interp->SetInputImage(reader->GetOutput());

  typename FiberGeneratorType::Pointer ftGen = FiberGeneratorType::New();
  EXERCISE_BASIC_OBJECT_METHODS(ftGen, FiberGeneratorType);

  typename PickerType::Pointer dirPicker = PickerType::New();
  dirPicker->SetInterpolator(interp);
  // dirPicker->UseProbabilisticMethod();

  //initialize the generator;
  ftGen->SetDirectionPicker(dirPicker);
  ftGen->SetStepLength(0.2);
  ftGen->SetCurvatureThreshold(3.44);
  ftGen->Initialize();
  
  PointType seed;
  seed[0] = 14; seed[1] = 6;  seed[2] = 0;

  ftGen->GenerateFiberPoints( seed, FiberGeneratorType::ONE_DIRECTION);
  // vtkSmartPointer<vtkPolyData> polyData = ftGen->GetVTKPolyData();
  FiberPointerType fiber = ftGen->GetFiber();
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(fiber);
  writer->SetFileName("fiber1Dir.meta");
  writer->Update();

  //TODO Check output
  // viewPolyData(polyData);

  ftGen->GenerateFiberPoints( seed, FiberGeneratorType::BOTH_DIRECTIONS);
  // polyData = ftGen->GetVTKPolyData();
  
  //TODO Check output
  // viewPolyData(polyData);

  seed[0] = 10; seed[1] = 0;  seed[2] = 0;
  ftGen->GenerateFiberPoints( seed, FiberGeneratorType::BOTH_DIRECTIONS);
  // ftGen->AddToVTKPolyData( polyData );

  //TODO Check output
  // viewPolyData(polyData);

  ftGen->SetCurvatureThreshold(100000000000000000);
  ftGen->GenerateFiberPoints( seed, FiberGeneratorType::BOTH_DIRECTIONS);
  // ftGen->AddToVTKPolyData( polyData );
  //TODO Check output
  //viewPolyData(polyData);
  

  return EXIT_FAILURE;
}

// template< class TPickerType >
// int test2(std::string dataFile, std::string stopFile, std::string stopFile2, std::string roiFile)
// {
//   typedef TPickerType                                 PickerType;
//   typedef typename PickerType::InterpolatorType       InterpolatorType;
//   typedef typename PickerType::PointType              PointType;
//   typedef typename InterpolatorType::InputImageType   ImageType;
//   typedef itk::ImageFileReader< ImageType >           ReaderType;

//   typedef unsigned int                                LabelPixelType;
//   typedef itk::FiberGenerator<PickerType,LabelPixelType>
//                                                       FiberGeneratorType;

// //  typedef typename FiberGeneratorType::LabelVectorType   LabelVectorType;
//   typedef typename FiberGeneratorType::LabelSetType   LabelSetType;
  
//   typedef itk::Image<LabelPixelType,3>                RoiImageType;
//   typedef itk::ImageFileReader< RoiImageType >        RoiReaderType;

//   typename ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( dataFile );
//   reader->Update();
  
//   typename RoiReaderType::Pointer roiReader = RoiReaderType::New();
//   roiReader->SetFileName( roiFile );
//   roiReader->Update();
  
//   typename InterpolatorType::Pointer interp = InterpolatorType::New();
//   interp->SetInputImage(reader->GetOutput());

//   //Load the mask File...
//   MaskSpatialObjectType::Pointer mask = loadStoppingMask(stopFile,0.1,1.1);
//   MaskSpatialObjectType::Pointer mask2 = loadStoppingMask(stopFile2,0.1,1.1);

//   typename FiberGeneratorType::Pointer ftGen = FiberGeneratorType::New();
  
//   typename PickerType::Pointer dirPicker = PickerType::New();
//   dirPicker->SetInterpolator(interp);
//   // dirPicker->UseProbabilisticMethod();

//   //initialize the generator;
//   ftGen->SetDirectionPicker(dirPicker);
//   ftGen->SetStepLength(0.2);
//   ftGen->SetCurvatureThreshold(3.44);

//   RoiImageType::Pointer roiIm = roiReader->GetOutput();

//   ftGen->SetROIImage(roiIm);
//   TEST_SET_GET(roiIm, ftGen->GetROIImage());

//   ftGen->Initialize();
  
//   PointType seed;
//   seed[0] = 14; seed[1] = 6;  seed[2] = 0;

//   ftGen->GenerateFiberPoints( seed, FiberGeneratorType::BOTH_DIRECTIONS);
//   vtkSmartPointer<vtkPolyData> polyData = ftGen->GetVTKPolyData();
// //  typename LabelVectorType::ConstPointer traversedRois = ftGen->GetTraversedROIs();
//   LabelSetType traversedRois = ftGen->GetTraversedROIs();
  
//   //TODO Check output
//   std::cout << ftGen->GetFiberLength() << std::endl;
//   printRois(traversedRois);
//   //viewPolyData(polyData);

//   ftGen->AddStoppingCriteriaSpatialObject(mask.GetPointer());
//   ftGen->GenerateFiberPoints( seed, FiberGeneratorType::BOTH_DIRECTIONS);
//   polyData = ftGen->GetVTKPolyData();
//   traversedRois = ftGen->GetTraversedROIs();
  
//   //TODO Check output
//   std::cout << ftGen->GetFiberLength() << std::endl;
//   printRois(traversedRois);
//   //viewPolyData(polyData);

//   ftGen->AddStoppingCriteriaSpatialObject(mask2.GetPointer());
//   ftGen->GenerateFiberPoints( seed, FiberGeneratorType::BOTH_DIRECTIONS);
//   polyData = ftGen->GetVTKPolyData();
//   traversedRois = ftGen->GetTraversedROIs();
  
//   //TODO Check output
//   std::cout << ftGen->GetFiberLength() << std::endl;
//   printRois(traversedRois);
//   //viewPolyData(polyData);

//   return EXIT_FAILURE;
// }


} //End empty Namespace

using namespace itk;
int itkFiberGeneratorTest( int , char * argv[] )
{
  std::string dtiFile   =  std::string(argv[1]);
  // std::string maskFile  =  std::string(argv[1])+"sim_dti_orig-stop1.nii.gz";
  // std::string maskFile2 =  std::string(argv[1])+"sim_dti_orig-stop2.nii.gz";
  // std::string roiFile   =  std::string(argv[1])+"sim_dti_orig-ROIs.nii.gz";

  test1<DtiDirPickerType>(dtiFile);
//  test2<DtiDirPickerType>(dtiFile,maskFile,maskFile2,roiFile);
  

  std::cout << "*************************************************************************************" << std::endl
            << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;


  return EXIT_FAILURE;
}
