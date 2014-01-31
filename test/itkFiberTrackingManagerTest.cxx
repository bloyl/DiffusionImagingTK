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

#include <itkTrackerDirectionPickerDTI.h>
#include <itkFiberTrackingManager.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectWriter.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>

namespace{

typedef itk::DiffusionTensor3D<double>                    DtiPixelType;
typedef itk::Image<DtiPixelType,3>                        DtiImageType;
typedef itk::NearestNeighborInterpolateImageFunction<DtiImageType,double>
                                                          DtiNNInterpType;

typedef itk::TrackerDirectionPickerDTI<DtiNNInterpType>   DtiDirPickerType;

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

  typedef typename FTManagerType::GroupSpatialObjectType
                                                      GroupSpatialObjectType;
  typedef typename FTManagerType::GroupSpatialObjectPointer
                                                      GroupSpatialObjectPointer;
  typedef typename FTManagerType::GroupSpatialObjectConstPointer
                                                      GroupSpatialObjectConstPointer;
                                                                    
  typedef typename PickerType::InterpolatorType       InterpolatorType;
  typedef typename PickerType::PointType              PointType;
  typedef typename InterpolatorType::InputImageType   ImageType;
  typedef itk::ImageFileReader< ImageType >           ReaderType;

  typedef typename FTManagerType::LabelImageType      RoiImageType;
  typedef itk::ImageFileReader< RoiImageType >        RoiReaderType;

  typedef double                                      StopPixelType;
  typedef itk::Image<StopPixelType,3>                 StopImageType;
  typedef itk::ImageFileReader< StopImageType >       StopReaderType;

  typedef itk::SpatialObjectWriter<3>                 WriterType;

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

//    ftManager->SetNumberOfFibersToKeep(500);
    ftManager->SetNumberOfFibersToKeep(5);
    ftManager->SetMinimumFiberLength(10.0);
    
    ftManager->SetDirectionPicker( dirPicker );

    //Should Error!
    ftManager->GenerateFibers();

    GroupSpatialObjectPointer keptFibers = ftManager->GetFibers();
    GroupSpatialObjectPointer rejectedFibers = ftManager->GetRejectedFibers();

    // std::cout << "Number of keptFibers     : "<<keptFibers->GetNumberOfLines() << std::endl;
    // std::cout << "Number of rejectedFibers : "<<rejectedFibers->GetNumberOfLines() << std::endl;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(keptFibers);
    writer->SetFileName(outFiberFile);
    writer->Update();

    std::cout << "number of Fibers : " << keptFibers->GetNumberOfChildren() << std::endl;
    // std::cout << keptFibers << std::endl;
  
    typename GroupSpatialObjectType::ChildrenListType * childrenList = keptFibers->GetChildren();
    std::cout << "fibers has " << childrenList->size() << " child" << std::endl;
    typename GroupSpatialObjectType::ChildrenListType::const_iterator it = childrenList->begin();
    while(it != childrenList->end())
    {
      std::cout << "Name of the child of the object 1: ";
      std::cout << (*it)->GetProperty()->GetName() << std::endl;
      std::cout << (*it) << std::endl;
      it++;
    }
    delete childrenList;
    
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
int itkFiberTrackingManagerTest( int, char * argv[] )
{
  std::string dtiFile             = std::string(argv[1]);
  std::string dtiStop1            = std::string(argv[2]);
  std::string roisFile            = std::string(argv[3]);
  std::string seedTestFile        = std::string(argv[4]);
  std::string singleSeedTestFile  = std::string(argv[5]);
  std::string outFiberFile        = std::string(argv[6]);

//  test<DtiDirPickerType>(dtiFile,roisFile,dtiStop1,outFiberFile);
  test<DtiDirPickerType>(dtiFile,seedTestFile,dtiStop1,outFiberFile,false);
  // test<DtiDirPickerType>(dtiFile,seedTestFile,dtiStop1,outFiberFile,true);

  // test<DtiDirPickerType>(dtiFile,singleSeedTestFile,dtiStop1,outFiberFile,false);
  // test<DtiDirPickerType>(dtiFile,singleSeedTestFile,dtiStop1,outFiberFile,true);

  std::cout << "*************************************************************************************" << std::endl
            << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;

  return EXIT_FAILURE;
}
