/**
 * @file  itkSymRealSphericalHarmonicRepTest.cxx
 * @brief Test itkSymRealSphericalHarmonicRep module.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkTrackerDirectionPickerDTI.h>
#include <itkFiberTrackingManager.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>


#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>


namespace{

typedef itk::DiffusionTensor3D<double>                    DtiPixelType;
typedef itk::Image<DtiPixelType,3>                        DtiImageType;
typedef itk::NearestNeighborInterpolateImageFunction<DtiImageType,double>
                                                          DtiNNInterpType;
typedef itk::TrackerDirectionPickerDTI<DtiNNInterpType>   DtiDirPickerType;


void viewPolyData( vtkSmartPointer<vtkPolyData> polyData )
{

  // Setup actor and mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput(polyData);
#else
  mapper->SetInputData(polyData);
#endif
  
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  
  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
                            vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor);

  renderWindow->Render();
  renderWindowInteractor->Start();

}

//test 1 roiSeeding... with a stopMask
template< class TPickerType >
int test( std::string dataFile, std::string roiFile, std::string stopFile,
          std::string outFiberFile, bool useProb)
{

  //Basic Typedefs
  typedef TPickerType                                 PickerType;
  typedef unsigned int                                LabelPixelType;
  typedef itk::FiberTrackingManager<PickerType,LabelPixelType>
                                                      FTManagerType;

  typedef typename FTManagerType::VectorOfLabelSetsConstPointer
                                                      VectorOfLabelSetsConstPointer;
                                                                    
  typedef typename PickerType::InterpolatorType       InterpolatorType;
  typedef typename PickerType::PointType              PointType;
  typedef typename InterpolatorType::InputImageType   ImageType;
  typedef itk::ImageFileReader< ImageType >           ReaderType;

  typedef typename FTManagerType::LabelImageType      RoiImageType;
  typedef itk::ImageFileReader< RoiImageType >        RoiReaderType;

  typedef double                                      StopPixelType;
  typedef itk::Image<StopPixelType,3>                 StopImageType;
  typedef itk::ImageFileReader< StopImageType >       StopReaderType;

  //Read in dataFile
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( dataFile );
  reader->Update();
  
  //Read in Rois
  typename RoiReaderType::Pointer roiReader = RoiReaderType::New();
  roiReader->SetFileName( roiFile );
  roiReader->Update();
  
  //Read in StopFile
  typename StopReaderType::Pointer stopReader = StopReaderType::New();
  stopReader->SetFileName( stopFile );
  stopReader->Update();
  
  //Setup interpolator
  typename InterpolatorType::Pointer interp = InterpolatorType::New();
  interp->SetInputImage(reader->GetOutput());

  //Setup the direction Picker
  typename PickerType::Pointer dirPicker = PickerType::New();
  dirPicker->SetInterpolator(interp);

  // if (useProb)
  //   dirPicker->UseProbabilisticMethod();

  //Make the FiberTracking Manager
  typename FTManagerType::Pointer ftManager = FTManagerType::New();

  try
  {
    reader->Update();
    roiReader->Update();

    ftManager->AddSeedImage(roiReader->GetOutput(),1U,100U);

    //We want to stop when the stop Image drops below 1.
    ftManager->AddStoppingCriteria(stopReader->GetOutput(),0.0,0.0);
    ftManager->SetROIImage(roiReader->GetOutput());

    ftManager->SetNumberOfFibersToKeep(500);
    ftManager->SetMinimumFiberLength(10.0);
    
    ftManager->SetDirectionPicker( dirPicker );

    //Should Error!
    ftManager->GenerateFibers();

    vtkSmartPointer<vtkPolyData> keptFibers = ftManager->GetFibers();
    vtkSmartPointer<vtkPolyData> rejectedFibers = ftManager->GetRejectedFibers();

    std::cout << "Number of keptFibers     : "<<keptFibers->GetNumberOfLines() << std::endl;
    std::cout << "Number of rejectedFibers : "<<rejectedFibers->GetNumberOfLines() << std::endl;
    //viewPolyData(keptFibers);
    //viewPolyData(rejectedFibers);


    //the Fibers generated are stored in LPS coordinate system....
    // SLicer and maybe some others expect vtk points in RAS so lets transform them

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
      transformFilter->SetInput(ftManager->GetFibers());
    #else
      transformFilter->SetInputData(ftManager->GetFibers());
    #endif
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    
    vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    mat->Identity();
    mat->SetElement(0,0,-1);
    mat->SetElement(1,1,-1);
    transform->SetMatrix(mat);
    transformFilter->SetTransform(transform);
    
    transformFilter->Update();

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    writer->SetFileName(outFiberFile.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(transformFilter->GetOutput());
    #else
      writer->SetInputData(transformFilter->GetOutput());
    #endif
    writer->Write();
    
    //viewPolyData(transformFilter->GetOutput());
    
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_FAILURE;
}

} // end empty namespace

using namespace itk;
int itkFiberTrackingManagerTest( int argc, char * argv[] )
{
  std::string dtiFile   =  std::string(argv[1])+"DWIS_dti-scheme_SNR-30_DTI.nii.gz";
  std::string dtiStop1  =  std::string(argv[1])+"DWIS_dti-scheme_SNR-30_DTI_FA_mask.nii.gz";
  std::string roisFile  =  std::string(argv[1])+"seeding_regions.nii.gz";
  std::string seedTestFile  
                        =  std::string(argv[1])+"DWIS_dti-scheme_seedTest_mask.nii.gz";

  std::string singleSeedTestFile  
                        =  std::string(argv[1])+"DWIS_dti-scheme_singleSeedTest_mask.nii.gz";

  std::string outFiberFile   =  std::string(argv[1])+"sim_dti_orig-fibers.vtp";

//  test<DtiDirPickerType>(dtiFile,roisFile,dtiStop1,outFiberFile);
  test<DtiDirPickerType>(dtiFile,seedTestFile,dtiStop1,outFiberFile,false);
  test<DtiDirPickerType>(dtiFile,seedTestFile,dtiStop1,outFiberFile,true);

  test<DtiDirPickerType>(dtiFile,singleSeedTestFile,dtiStop1,outFiberFile,false);
  test<DtiDirPickerType>(dtiFile,singleSeedTestFile,dtiStop1,outFiberFile,true);

  std::cout << "*************************************************************************************" << std::endl
            << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;

  return EXIT_FAILURE;
}
