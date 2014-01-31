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
#include <itkTrackerDirectionPickerDWI.h>
#include <itkFiberTrackingManager.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

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

typedef short                                             DwiPixelType;
typedef itk::VectorImage<DwiPixelType, 3>                 DwiImageType;
typedef itk::Image<DwiPixelType, 4>                       DwiImageType4D;

typedef itk::NearestNeighborInterpolateImageFunction<DtiImageType,double>
                                                          DtiNNInterpType;
typedef itk::NearestNeighborInterpolateImageFunction<DwiImageType,double>
                                                          DwiNNInterpType;

typedef itk::LinearInterpolateImageFunction<DtiImageType,double>
                                                          DtiLinearInterpType;
typedef itk::LinearInterpolateImageFunction<DwiImageType,double>
                                                          DwiLinearInterpType;

void convertDWIimage(DwiImageType4D::Pointer inp4d, DwiImageType::Pointer outputIm )
{


  typedef itk::ImageRegionIteratorWithIndex< DwiImageType >         IteratorType;

  //Set up the gradient image size
  DwiImageType::SizeType  sizeGradImage;
  DwiImageType4D::SizeType size4D = inp4d->GetLargestPossibleRegion().GetSize();
  sizeGradImage[0] = size4D[0];
  sizeGradImage[1] = size4D[1];
  sizeGradImage[2] = size4D[2];
  outputIm->SetVectorLength(size4D[3]);

  DwiImageType::IndexType   indexGradImage = {{ 0, 0, 0 }};
  DwiImageType::RegionType  regionGradImage;
  regionGradImage.SetSize(  sizeGradImage );
  regionGradImage.SetIndex( indexGradImage);
  outputIm->SetRegions( regionGradImage );
  DwiImageType4D::SpacingType img4Dspacing = inp4d->GetSpacing();
  DwiImageType4D::PointType img4Dorigin = inp4d->GetOrigin();
  DwiImageType4D::DirectionType img4Ddir = inp4d->GetDirection();

  DwiImageType::SpacingType gradSpacing;
  DwiImageType::PointType gradOrigin;
  DwiImageType::DirectionType gradDirs;

  gradSpacing[0]  = img4Dspacing[0];  gradSpacing[1] = img4Dspacing[1];   gradSpacing[2] = img4Dspacing[2];
  gradOrigin[0]   = img4Dorigin[0];   gradOrigin[1]  = img4Dorigin[1];    gradOrigin[2] = img4Dorigin[2];

  for (unsigned int i = 0; i<3; ++i)
  {
    for (unsigned int j = 0; j<3; ++j)
    {
      gradDirs[i][j] = img4Ddir[i][j];
    }
  }

  outputIm->SetSpacing( gradSpacing );
  outputIm->SetOrigin( gradOrigin );
  outputIm->SetDirection( gradDirs );

  outputIm->Allocate();

  printf("Done GradIm->Allocate\n");

  //Copy data from img4d to gradim
  //THIS IS SLOW!!!
  IteratorType it( outputIm, outputIm->GetRequestedRegion() );

  //Probably a better way to do this but I don't really know what it is.
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    DwiImageType::IndexType   gradIndex = it.GetIndex();
    DwiImageType::PixelType   gradPix = it.Get();
    DwiImageType4D::IndexType  img4dIndex;
    img4dIndex[0] = gradIndex[0];
    img4dIndex[1] = gradIndex[1];
    img4dIndex[2] = gradIndex[2];

    for ( unsigned int i=0; i<size4D[3]; ++i )
    {
      img4dIndex[3] = i;
      gradPix.SetElement( i, inp4d->GetPixel( img4dIndex ) );
    }
    it.Set( gradPix );
  }

  printf("Done Conversion\n");

}

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

template <typename GradientDirectionContainerType>
void loadGradDirs( typename GradientDirectionContainerType::Pointer gradDirs, std::string bvecFile, std::string bvalFile )
{

  //Parse bvalFile
  std::vector<double> bValues;

  std::ifstream bvalIn(bvalFile.c_str());
  //first count the number of grad directions and make sure the only bvalues in the bval file
  // are 0 and a number...
  std::string line;
  while (! bvalIn.eof() )
  {
    getline(bvalIn,line);
    double val;
    std::stringstream ss(line);
    while (ss >> val)
    {
      bValues.push_back(val);
    }
  }
  
  //ok so now lets process the gradient table...
  //read in each line and put it in a string stream to process...
  std::ifstream bvecIn(bvecFile.c_str());
  getline(bvecIn,line);
  std::stringstream Xss(line);
  getline(bvecIn,line);
  std::stringstream Yss(line);
  getline(bvecIn,line);
  std::stringstream Zss(line);

  typename GradientDirectionContainerType::Element vect3d;

  int counter = 0;
  double x,y,z;
  double scale;
  while (Xss >> x)
  {
    scale = vcl_sqrt(bValues[counter]);
    Yss >> y;
    Zss >> z;
    vect3d[0] = scale * x; vect3d[1] = scale * y; vect3d[2] = scale * z;
    gradDirs->InsertElement( counter, vect3d );
    ++counter;
  }

  if (counter != bValues.size())
  {
    std::cerr << "different number of bvalues and gradients" << std::endl;
    gradDirs->Initialize();
  }

}


void writeVtkPolyDataAsRAS( vtkSmartPointer<vtkPolyData> polyData, std::string outputFile )
{
  //the Fibers generated are stored in LPS coordinate system....
  // SLicer and maybe some others expect vtk points in RAS so lets transform them

  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    transformFilter->SetInput(polyData);
  #else
    transformFilter->SetInputData(polyData);
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

  writer->SetFileName(outputFile.c_str());
  #if VTK_MAJOR_VERSION <= 5
    writer->SetInput(transformFilter->GetOutput());
  #else
    writer->SetInputData(transformFilter->GetOutput());
  #endif
  writer->Write();
}

int test( std::string dwiFile, std::string bvecFile, std::string bvalFile,
          std::string dtiFile, std::string faFile, std::string seedFile,
          std::string outDtiFiberFile,std::string outDwiFiberFile)
{

  typedef double                                      FaPixelType;
  typedef itk::Image<FaPixelType,3>                   FaImageType;
  typedef itk::ImageFileReader< FaImageType >         FaReaderType;

  typedef double                                      SeedPixelType;
  typedef itk::Image<SeedPixelType,3>                 SeedImageType;
  typedef itk::ImageFileReader< SeedImageType >       SeedReaderType;

  typedef unsigned int                                LabelPixelType;


  typedef itk::TrackerDirectionPickerDWI<DwiLinearInterpType>
                                                      DwiLinearDirPickerType;

  typedef itk::TrackerDirectionPickerDTI<DtiLinearInterpType>
                                                      DtiLinearDirPickerType;

  typedef itk::FiberTrackingManager<DtiLinearDirPickerType,LabelPixelType>
                                                      DtiFTManagerType;

  typedef itk::FiberTrackingManager<DwiLinearDirPickerType,LabelPixelType>
                                                      DwiFTManagerType;

  
  typedef itk::ImageFileReader< DtiImageType >        DtiReaderType;
  typedef itk::ImageFileReader< DwiImageType4D >      DwiReaderType;


  try
  {
    //Read in faFile
    FaReaderType::Pointer faReader = FaReaderType::New();
    faReader->SetFileName( faFile );
    faReader->Update();
    
    //Read in seedFile
    SeedReaderType::Pointer seedReader = SeedReaderType::New();
    seedReader->SetFileName( seedFile );
    seedReader->Update();
    
    //Read in dtiFile
    DtiReaderType::Pointer dtiReader = DtiReaderType::New();
    dtiReader->SetFileName( dtiFile );
    dtiReader->Update();

    DwiReaderType::Pointer dwiReader = DwiReaderType::New();
    dwiReader->SetFileName( dwiFile );
    dwiReader->Update();

    DwiImageType::Pointer dwiImage = DwiImageType::New();
    convertDWIimage(dwiReader->GetOutput(), dwiImage);

    //Setup interpolators
    DtiLinearInterpType::Pointer dtiInterp = DtiLinearInterpType::New();
    dtiInterp->SetInputImage(dtiReader->GetOutput());

    DwiLinearInterpType::Pointer dwiInterp = DwiLinearInterpType::New();
    dwiInterp->SetInputImage(dwiImage);


    //Setup the direction Picker
    DtiLinearDirPickerType::Pointer dtiDirPicker = DtiLinearDirPickerType::New();
    dtiDirPicker->SetInterpolator(dtiInterp);
    dtiDirPicker->UseTensorDeflectionTracking();
    

    DwiLinearDirPickerType::Pointer dwiDirPicker = DwiLinearDirPickerType::New();
    //Load the Gradient directions
    typedef DwiLinearDirPickerType::GradientDirectionContainerType  GradientDirectionContainerType;
    GradientDirectionContainerType::Pointer     gradDirs = GradientDirectionContainerType::New();
    loadGradDirs<GradientDirectionContainerType>( gradDirs, bvecFile, bvalFile );

    dwiDirPicker->SetGradientDirectionContainer(gradDirs);
    dwiDirPicker->SetInterpolator(dwiInterp);
    dwiDirPicker->UseDTITensorDeflectionTracking();
    
    std::cerr << dwiImage->GetNumberOfComponentsPerPixel( )  << " : " << gradDirs->Size() << std::endl;
    
    //Make the FiberTracking Manager
    DtiFTManagerType::Pointer dtiFtManager = DtiFTManagerType::New();
    typedef DtiFTManagerType::PointSetType                    DtiPointSetType;
    typedef DtiPointSetType::ConstPointer                     DtiPointSetConstPointer;

    DwiFTManagerType::Pointer dwiFtManager = DwiFTManagerType::New();
    typedef DwiFTManagerType::PointSetType                    DwiPointSetType;
    typedef DwiPointSetType::ConstPointer                     DwiPointSetConstPointer;


    dtiFtManager->AddSeedImage(seedReader->GetOutput(),1.0,100.0);
    dtiFtManager->AddStoppingCriteria(faReader->GetOutput(),0.0,0.15);

    dtiFtManager->SetNumberOfFibersToKeep(100);
    dtiFtManager->SetMinimumFiberLength(10.0);
    
    dtiFtManager->SetDirectionPicker( dtiDirPicker );
    dtiFtManager->GenerateFibers();

    vtkSmartPointer<vtkPolyData> dtiKeptFibers = dtiFtManager->GetFibers();
    vtkSmartPointer<vtkPolyData> dtiRejectedFibers = dtiFtManager->GetRejectedFibers();

    std::cout << "Number of keptFibers     : "<<dtiKeptFibers->GetNumberOfLines() << std::endl;
    std::cout << "Number of rejectedFibers : "<<dtiRejectedFibers->GetNumberOfLines() << std::endl;

    //Try DWI filtering
    dwiFtManager->AddSeedImage(seedReader->GetOutput(),1.0,100.0);
    dwiFtManager->AddStoppingCriteria(faReader->GetOutput(),0.0,0.15);

    dwiFtManager->SetNumberOfFibersToKeep(100);
    dwiFtManager->SetMinimumFiberLength(10.0);
    
    dwiDirPicker->Initialize();

    dwiFtManager->SetDirectionPicker( dwiDirPicker );
    dwiFtManager->GenerateFibers();

    vtkSmartPointer<vtkPolyData> dwiKeptFibers = dwiFtManager->GetFibers();
    vtkSmartPointer<vtkPolyData> dwiRejectedFibers = dwiFtManager->GetRejectedFibers();

    std::cout << "Number of keptFibers     : "<<dwiKeptFibers->GetNumberOfLines() << std::endl;
    std::cout << "Number of rejectedFibers : "<<dwiRejectedFibers->GetNumberOfLines() << std::endl;

    writeVtkPolyDataAsRAS(dtiKeptFibers,outDtiFiberFile);
    writeVtkPolyDataAsRAS(dwiKeptFibers,outDwiFiberFile);

    //viewPolyData(dtiKeptFibers);
    //viewPolyData(dwiKeptFibers);


    // DtiPointSetType::PointType seedPoint;
    // DtiPointSetType::PixelType pointData;

    // DtiPointSetConstPointer dtiSeedPoints = dtiFtManager->GetSeedPoints();
    // for (unsigned int i=0; i<dtiSeedPoints->GetNumberOfPoints(); ++i)
    // {
    //   seedPoint = dtiSeedPoints->GetPoint(i);
    //   dtiSeedPoints->GetPointData(i,&pointData);
    //   if (pointData == 1)
    //   {
    //     std::cout << "keep      : " << seedPoint << std::endl;
    //   }
    //   else
    //   {
    //     std::cout << "rejected  : " << seedPoint << std::endl;
    //   }
    // }

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
int itkFiberTrackingManagerTest2( int argc, char * argv[] )
{
  std::string dwiFile   =  std::string(argv[1])+"100003_t0_low_b.nii.gz";
  std::string bvecFile  =  std::string(argv[1])+"100003_t0_low_b.bvec";
  std::string bvalFile  =  std::string(argv[1])+"100003_t0_low_b.bval";

  std::string dtiFile   =  std::string(argv[1])+"100003_t0_low_b_tensor.nii.gz";

  // std::string seedFile  =  std::string(argv[1])+"100003_t0_low_b_tensor_cc_singSeedMask.nii.gz";
  std::string seedFile  =  std::string(argv[1])+"100003_t0_low_b_tensor_cc_singSeedMaskDil.nii.gz";
  std::string faFile    =  std::string(argv[1])+"100003_t0_low_b_tensor_fa.nii.gz";

  std::string outDtiFiberFile   =  std::string(argv[1])+"100003_t0_low_b_tensor_dtiFibers.vtp";
  std::string outDwiFiberFile   =  std::string(argv[1])+"100003_t0_low_b_tensor_dwiFibers.vtp";

  test( dwiFile, bvecFile, bvalFile, dtiFile, faFile, seedFile, outDtiFiberFile, outDwiFiberFile);

  return EXIT_SUCCESS;
}
